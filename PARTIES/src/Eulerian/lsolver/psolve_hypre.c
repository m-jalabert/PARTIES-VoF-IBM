/*******************************************************************************
 * psolve_hypre.c
 *
 * HYPRE-based pressure solver for the variable-coefficient Poisson equation:
 *
 *     ∇ · ( (1/ρ) ∇φ ) = rhs
 *
 * Uses the HYPRE Struct interface with:
 *   - PCG (conjugate gradient) as the outer solver
 *   - PFMG (parallel semicoarsening multigrid) as the preconditioner
 * 
 * Features:
 *   - "Non-Galerkin" coarse operators to prevent crashes at 1000:1 density jumps.
 *   - Creates and destroys solvers dynamically per-solve to completely 
 *     eliminate HYPRE memory leaks.
 ******************************************************************************/

#ifdef USE_HYPRE

#include "HYPRE_struct_ls.h"
#include "HYPRE.h"

/* -------------------------------------------------------------------------- *
 *  Module-static HYPRE objects (persist across time steps)
 * -------------------------------------------------------------------------- */
static HYPRE_StructGrid    hypre_grid    = NULL;
static HYPRE_StructStencil hypre_stencil = NULL;
static HYPRE_StructMatrix  hypre_A       = NULL;
static HYPRE_StructVector  hypre_b       = NULL;
static HYPRE_StructVector  hypre_x       = NULL;

/* Local owned-cell counts (no ghost) */
static int hypre_Is, hypre_Js, hypre_Ks;
static int hypre_Ie, hypre_Je, hypre_Ke;
static int hypre_nx_local, hypre_ny_local, hypre_nz_local;
static int hypre_n_local;

/* Reusable flat buffers for HYPRE set/get values */
static double *hypre_vals   = NULL;  
static double *hypre_rhs    = NULL;  
static double *hypre_sol    = NULL;  
static int    *hypre_stencil_indices = NULL; 

/* Mapping: stencil entry → (di,dj,dk) offset */
static const int stencil_offsets[7][3] = {
    { 0, 0, 0},   /* center */
    {-1, 0, 0},   /* west   (-x) */
    { 1, 0, 0},   /* east   (+x) */
    { 0,-1, 0},   /* south  (-y) */
    { 0, 1, 0},   /* north  (+y) */
    { 0, 0,-1},   /* back   (-z) */
    { 0, 0, 1}    /* front  (+z) */
};


/*******************************************************************************
 * Pressure_hypre_setup
 ******************************************************************************/
void Pressure_hypre_setup(Cart3d_bag *data_bag)
{
    MAC_grid   *grid   = data_bag->grid;
    Parameters *params = data_bag->params;
    int rank = params->rank;

    HYPRE_Initialize();

    hypre_Is = grid->G_Is;
    hypre_Js = grid->G_Js;
    hypre_Ks = grid->G_Ks;
    hypre_Ie = (grid->G_Ie < grid->NX - 1) ? grid->G_Ie : grid->NX - 1;
    hypre_Je = (grid->G_Je < grid->NY - 1) ? grid->G_Je : grid->NY - 1;
    hypre_Ke = (grid->G_Ke < grid->NZ - 1) ? grid->G_Ke : grid->NZ - 1;

    hypre_nx_local = hypre_Ie - hypre_Is;
    hypre_ny_local = hypre_Je - hypre_Js;
    hypre_nz_local = hypre_Ke - hypre_Ks;
    hypre_n_local  = hypre_nx_local * hypre_ny_local * hypre_nz_local;

    if (rank == 0) {
        printf("[HYPRE] Setting up Struct Matrix and Vectors (Robust PFMG Active)\n");
        printf("[HYPRE] Global grid: %d x %d x %d  (cells)\n", grid->NX - 1, grid->NY - 1, grid->NZ - 1);
    }

    HYPRE_StructGridCreate(PCW, 3, &hypre_grid);

    HYPRE_Int ilower[3] = { hypre_Is, hypre_Js, hypre_Ks };
    HYPRE_Int iupper[3] = { hypre_Ie - 1, hypre_Je - 1, hypre_Ke - 1 };
    HYPRE_StructGridSetExtents(hypre_grid, ilower, iupper);

    HYPRE_Int periodic[3] = {0, 0, 0};
#ifdef XPERIODIC
    periodic[0] = grid->NX - 1;
#endif
#ifdef YPERIODIC
    periodic[1] = grid->NY - 1;
#endif
#ifdef ZPERIODIC
    periodic[2] = grid->NZ - 1;
#endif
    HYPRE_StructGridSetPeriodic(hypre_grid, periodic);
    HYPRE_StructGridAssemble(hypre_grid);

    HYPRE_StructStencilCreate(3, 7, &hypre_stencil);
    for (int s = 0; s < 7; s++) {
        HYPRE_Int offset[3] = { stencil_offsets[s][0],
                                stencil_offsets[s][1],
                                stencil_offsets[s][2] };
        HYPRE_StructStencilSetElement(hypre_stencil, s, offset);
    }

    HYPRE_StructMatrixCreate(PCW, hypre_grid, hypre_stencil, &hypre_A);
    HYPRE_StructMatrixInitialize(hypre_A);

    HYPRE_StructVectorCreate(PCW, hypre_grid, &hypre_b);
    HYPRE_StructVectorInitialize(hypre_b);

    HYPRE_StructVectorCreate(PCW, hypre_grid, &hypre_x);
    HYPRE_StructVectorInitialize(hypre_x);

    hypre_vals = (double *)calloc((size_t)7 * hypre_n_local, sizeof(double));
    hypre_rhs  = (double *)calloc((size_t)hypre_n_local, sizeof(double));
    hypre_sol  = (double *)calloc((size_t)hypre_n_local, sizeof(double));

    hypre_stencil_indices = (int *)malloc(7 * sizeof(int));
    for (int s = 0; s < 7; s++) hypre_stencil_indices[s] = s;
}

/* Is there any Dirichlet BC? */
#if defined(LEFT_INFLOW)   || defined(RIGHT_INFLOW)  || \
    defined(LEFT_DIRICHLET)|| defined(RIGHT_DIRICHLET) || \
    defined(BOTTOM_DIRICHLET) || defined(TOP_DIRICHLET) || \
    defined(BACK_DIRICHLET)   || defined(FRONT_DIRICHLET)
#   define NEED_REFERENCE_PRESSURE 0
#else
#   define NEED_REFERENCE_PRESSURE 1
#endif

/*******************************************************************************
 * Pressure_solve_hypre
 ******************************************************************************/
int Pressure_solve_hypre(Cart3d_bag *data_bag)
{
    Pressure       *p      = data_bag->p;
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double ***rho   = vof->rho;
    double ***phi   = p->deltap;
    double ***rhs   = p->rhs;
    char statement[1024];

    const int NX = grid->NX;
    const int NY = grid->NY;
    const int NZ = grid->NZ;

    const int Is = hypre_Is, Js = hypre_Js, Ks = hypre_Ks;
    const int Ie = hypre_Ie, Je = hypre_Je, Ke = hypre_Ke;

    double t_total0 = MPI_Wtime();
    double t0, t1;
    double cpu_assemble_time = 0.0;
    double hypre_set_values_time = 0.0;
    double hypre_setup_time = 0.0;
    double hypre_solve_time = 0.0;
    double hypre_get_values_time = 0.0;
    double postprocess_time = 0.0;

    /* ================================================================== *
     *  1. Fill matrix coefficients and RHS vector
     * ================================================================== */
    t0 = MPI_Wtime();
    int n = 0;
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++, n++) {

                double beta_xm = 2.0 / (rho[k][j][i] + rho[k][j][i-1]);
                double beta_xp = 2.0 / (rho[k][j][i] + rho[k][j][i+1]);
                double beta_ym = 2.0 / (rho[k][j][i] + rho[k][j-1][i]);
                double beta_yp = 2.0 / (rho[k][j][i] + rho[k][j+1][i]);
                double beta_zm = 2.0 / (rho[k][j][i] + rho[k-1][j][i]);
                double beta_zp = 2.0 / (rho[k][j][i] + rho[k+1][j][i]);
                double cxm, cxp, cym, cyp, czm, czp;
                double row_w = TwodOps_pressure_row_weight(grid, params, i);

                TwodOps_build_scalar_face_coeffs(
                    grid, params, i, j, k,
                    beta_xm, beta_xp, beta_ym, beta_yp, beta_zm, beta_zp,
                    &cxm, &cxp, &cym, &cyp, &czm, &czp);

                cxm *= row_w;
                cxp *= row_w;
                cym *= row_w;
                cyp *= row_w;
                czm *= row_w;
                czp *= row_w;
                double diag = cxm + cxp + cym + cyp + czm + czp;

#ifndef XPERIODIC
                if (i == 0) { diag -= cxm; cxm = 0.0; }
                if (i == NX - 2) { diag -= cxp; cxp = 0.0; }
#endif
#ifndef YPERIODIC
                if (j == 0) { diag -= cym; cym = 0.0; }
                if (j == NY - 2) { diag -= cyp; cyp = 0.0; }
#endif
#ifndef ZPERIODIC
                if (k == 0) { diag -= czm; czm = 0.0; }
                if (k == NZ - 2) { diag -= czp; czp = 0.0; }
#endif

                hypre_vals[0 * hypre_n_local + n] =  diag;
                hypre_vals[1 * hypre_n_local + n] = -cxm;
                hypre_vals[2 * hypre_n_local + n] = -cxp;
                hypre_vals[3 * hypre_n_local + n] = -cym;
                hypre_vals[4 * hypre_n_local + n] = -cyp;
                hypre_vals[5 * hypre_n_local + n] = -czm;
                hypre_vals[6 * hypre_n_local + n] = -czp;

                hypre_rhs[n] = -row_w * rhs[k][j][i];
                hypre_sol[n] = phi[k][j][i];
            }
        }
    }

/* ================================================================== *
 *  2. Remove RHS null-space BEFORE solving pure Neumann pressure
 *
 *  IMPORTANT:
 *  In AXISYM_RZ the matrix rows are multiplied by row_w = r_c[i] to
 *  recover an SPD operator for PCG/PFMG:
 *
 *      A = - row_w * div( beta grad(phi) )
 *      b = - row_w * rhs
 *
 *  Therefore the compatibility correction must subtract a constant
 *  physical divergence from rhs, i.e.
 *
 *      b_i <- b_i - row_w_i * mean_b_per_weight
 *
 *  NOT b_i <- b_i - arithmetic_mean(b).
 *
 *  The latter is equivalent to subtracting a source proportional to 1/r
 *  in axisymmetry and can create a spurious pressure correction near r=0.
 * ================================================================== */
#if NEED_REFERENCE_PRESSURE
{
    double local_b_sum = 0.0, global_b_sum = 0.0;
    double local_w_sum = 0.0, global_w_sum = 0.0;

    int idx = 0;
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++, idx++) {
                const double row_w =
                    TwodOps_pressure_row_weight(grid, params, i);

                local_b_sum += hypre_rhs[idx];
                local_w_sum += row_w;
            }
        }
    }

    MPI_Allreduce(&local_b_sum, &global_b_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);
    MPI_Allreduce(&local_w_sum, &global_w_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);

    if (global_w_sum > 0.0) {
        /*
         * This is algebraic b_sum / sum(row_w).
         * Since hypre_rhs already stores b = -row_w*rhs, this makes
         * sum_i hypre_rhs_i exactly zero while subtracting a constant
         * physical rhs from the projection equation.
         */
        const double mean_b_per_weight = global_b_sum / global_w_sum;

        idx = 0;
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                for (int i = Is; i < Ie; i++, idx++) {
                    const double row_w =
                        TwodOps_pressure_row_weight(grid, params, i);

                    hypre_rhs[idx] -= row_w * mean_b_per_weight;
                }
            }
        }
    }
}
#endif
    t1 = MPI_Wtime();
    cpu_assemble_time += t1 - t0;

    /* ================================================================== *
     *  3. Set matrix values
     * ================================================================== */
    HYPRE_Int ilower[3] = { Is, Js, Ks };
    HYPRE_Int iupper[3] = { Ie - 1, Je - 1, Ke - 1 };

    t0 = MPI_Wtime();
    for (int s = 0; s < 7; s++) {
        HYPRE_StructMatrixSetBoxValues(hypre_A, ilower, iupper, 1, 
                                       &hypre_stencil_indices[s], 
                                       &hypre_vals[s * hypre_n_local]);
    }
    HYPRE_StructMatrixAssemble(hypre_A);

    HYPRE_StructVectorSetBoxValues(hypre_b, ilower, iupper, hypre_rhs);
    HYPRE_StructVectorSetBoxValues(hypre_x, ilower, iupper, hypre_sol);
    HYPRE_StructVectorAssemble(hypre_b);
    HYPRE_StructVectorAssemble(hypre_x);
    t1 = MPI_Wtime();
    hypre_set_values_time += t1 - t0;

    /* ================================================================== *
     *  4. Create fresh solver + preconditioner (Ensures ZERO memory leak)
     * ================================================================== */
    HYPRE_StructSolver solver, precond;

    HYPRE_StructPCGCreate(PCW, &solver);
    HYPRE_StructPCGSetMaxIter(solver, 2000);
    HYPRE_StructPCGSetTol(solver, (double)params->P_CG_ETOL);
    HYPRE_StructPCGSetTwoNorm(solver, 1);
    HYPRE_StructPCGSetRelChange(solver, 0);
    HYPRE_StructPCGSetPrintLevel(solver, 0); 
    HYPRE_StructPCGSetLogging(solver, 1);     

    HYPRE_StructPFMGCreate(PCW, &precond);
    HYPRE_StructPFMGSetMaxIter(precond, 1);
    HYPRE_StructPFMGSetTol(precond, 0.0);
    HYPRE_StructPFMGSetZeroGuess(precond);
    
    /* 2 = Symmetric red-black Gauss-Seidel */
    HYPRE_StructPFMGSetRelaxType(precond, 2);   
    
    /* Slightly increased sweeps to properly damp the 1000:1 interface */
    HYPRE_StructPFMGSetNumPreRelax(precond, 2);  
    HYPRE_StructPFMGSetNumPostRelax(precond, 2);
    
    /* CRITICAL STABILITY FIX: RAPType = 1 enforces strictly non-Galerkin 
     * 7-point coarse operators. This prevents division-by-zero crashes 
     * inside PFMG when the bubble interface distorts heavily. */
    HYPRE_StructPFMGSetRAPType(precond, 1);     
    
    HYPRE_StructPFMGSetSkipRelax(precond, 0);

    HYPRE_StructPCGSetPrecond(solver,
                              HYPRE_StructPFMGSolve,
                              HYPRE_StructPFMGSetup,
                              precond);

    /* ================================================================== *
     *  5. Setup + Solve
     * ================================================================== */
    t0 = MPI_Wtime();
    HYPRE_StructPCGSetup(solver, hypre_A, hypre_b, hypre_x);
    t1 = MPI_Wtime();
    hypre_setup_time += t1 - t0;

    t0 = MPI_Wtime();
    HYPRE_StructPCGSolve(solver, hypre_A, hypre_b, hypre_x);
    t1 = MPI_Wtime();
    hypre_solve_time += t1 - t0;

    HYPRE_Int num_iters;
    double final_res;
    HYPRE_StructPCGGetNumIterations(solver, &num_iters);
    HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res);

    /* ================================================================== *
     *  6. Extract solution 
     * ================================================================== */
    t0 = MPI_Wtime();
    HYPRE_StructVectorGetBoxValues(hypre_x, ilower, iupper, hypre_sol);
    t1 = MPI_Wtime();
    hypre_get_values_time += t1 - t0;

    /* Destroy solver + preconditioner to completely prevent memory leaks */
    HYPRE_StructPCGDestroy(solver);
    HYPRE_StructPFMGDestroy(precond);

    t0 = MPI_Wtime();
    n = 0;
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++, n++) {
                phi[k][j][i] = hypre_sol[n];
            }
        }
    }

/* ================================================================== *
 *  7. Remove null-space from SOLUTION
 *
 *  This is only a gauge choice.  It does not affect velocity because
 *  grad(phi) is unchanged.  Use the same row weight as the pressure
 *  operator for a clean axisymmetric mean.
 * ================================================================== */
#if NEED_REFERENCE_PRESSURE
{
    double local_phi_w_sum = 0.0, global_phi_w_sum = 0.0;
    double local_w_sum     = 0.0, global_w_sum     = 0.0;

    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                const double row_w =
                    TwodOps_pressure_row_weight(grid, params, i);

                local_phi_w_sum += row_w * phi[k][j][i];
                local_w_sum     += row_w;
            }
        }
    }

    MPI_Allreduce(&local_phi_w_sum, &global_phi_w_sum,
                  1, MPI_DOUBLE, MPI_SUM, PCW);
    MPI_Allreduce(&local_w_sum, &global_w_sum,
                  1, MPI_DOUBLE, MPI_SUM, PCW);

    if (global_w_sum > 0.0) {
        const double phi_mean = global_phi_w_sum / global_w_sum;

        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                for (int i = Is; i < Ie; i++) {
                    phi[k][j][i] -= phi_mean;
                }
            }
        }
    }
}
#endif

    /*
     * Mirror the CG path: after the solve, populate the physical pressure
     * ghosts with the same wall BCs used by the face-gradient projection and
     * by the momentum RHS in Velocity.c.
     */
    Pressure_apply_BCs(phi, grid, params);
    Communication_update_ghost_nodes_flow_variable(
        phi, 'h', params->ghost_nodes, data_bag);
    t1 = MPI_Wtime();
    postprocess_time += t1 - t0;

    sprintf(statement,
            "[HYPRE-CPU] PCG+PFMG converged to %g after %d iterations | "
            "wall %.6e s | cpu_assembly %.6e s | h2d_transfer %.6e s | "
            "d2h_transfer %.6e s | set_values %.6e s | setup %.6e s | "
            "solve %.6e s | get_values %.6e s | postprocess %.6e s\n",
            final_res, (int)num_iters, MPI_Wtime() - t_total0,
            cpu_assemble_time, 0.0, 0.0, hypre_set_values_time,
            hypre_setup_time, hypre_solve_time, hypre_get_values_time,
            postprocess_time);
    Display_progress(params, statement);

    return (int)num_iters;
}

/*******************************************************************************
 * Pressure_hypre_destroy
 ******************************************************************************/
void Pressure_hypre_destroy(void)
{
    if (hypre_A)       HYPRE_StructMatrixDestroy(hypre_A);
    if (hypre_b)       HYPRE_StructVectorDestroy(hypre_b);
    if (hypre_x)       HYPRE_StructVectorDestroy(hypre_x);
    if (hypre_stencil) HYPRE_StructStencilDestroy(hypre_stencil);
    if (hypre_grid)    HYPRE_StructGridDestroy(hypre_grid);

    hypre_A       = NULL;
    hypre_b       = NULL;
    hypre_x       = NULL;
    hypre_stencil = NULL;
    hypre_grid    = NULL;

    free(hypre_vals);            hypre_vals   = NULL;
    free(hypre_rhs);             hypre_rhs    = NULL;
    free(hypre_sol);             hypre_sol    = NULL;
    free(hypre_stencil_indices); hypre_stencil_indices = NULL;

    HYPRE_Finalize();
}

#endif /* USE_HYPRE */
