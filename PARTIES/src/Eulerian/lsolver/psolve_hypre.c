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
    char statement[200];

    const int NX = grid->NX;
    const int NY = grid->NY;
    const int NZ = grid->NZ;

    const int Is = hypre_Is, Js = hypre_Js, Ks = hypre_Ks;
    const int Ie = hypre_Ie, Je = hypre_Je, Ke = hypre_Ke;

    const double idx2 = grid->idx_c[1] * grid->idx_c[1];
    const double idy2 = grid->idy_c[1] * grid->idy_c[1];
    const double idz2 = grid->idz_c[1] * grid->idz_c[1];

    /* ================================================================== *
     *  1. Fill matrix coefficients and RHS vector
     * ================================================================== */
    int n = 0;
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++, n++) {

                double cxm = 2.0 / (rho[k][j][i] + rho[k][j][i-1]) * idx2; 
                double cxp = 2.0 / (rho[k][j][i] + rho[k][j][i+1]) * idx2; 
                double cym = 2.0 / (rho[k][j][i] + rho[k][j-1][i]) * idy2; 
                double cyp = 2.0 / (rho[k][j][i] + rho[k][j+1][i]) * idy2; 
                double czm = 2.0 / (rho[k][j][i] + rho[k-1][j][i]) * idz2; 
                double czp = 2.0 / (rho[k][j][i] + rho[k+1][j][i]) * idz2; 

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

                hypre_rhs[n] = -rhs[k][j][i];
                hypre_sol[n] = phi[k][j][i];
            }
        }
    }

    /* ================================================================== *
     *  2. Remove RHS null-space BEFORE solving for pure Neumann setup
     * ================================================================== */
#if NEED_REFERENCE_PRESSURE
    {
        double local_rhs_sum = 0.0, global_rhs_sum = 0.0;
        long   local_cnt = 0,   global_cnt = 0;

        for (int idx = 0; idx < hypre_n_local; idx++) {
            local_rhs_sum += hypre_rhs[idx];
            local_cnt++;
        }

        MPI_Allreduce(&local_rhs_sum, &global_rhs_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);
        MPI_Allreduce(&local_cnt, &global_cnt, 1, MPI_LONG, MPI_SUM, PCW);

        if (global_cnt > 0) {
            double rhs_avg = global_rhs_sum / (double)global_cnt;
            for (int idx = 0; idx < hypre_n_local; idx++) {
                hypre_rhs[idx] -= rhs_avg;
            }
        }
    }
#endif

    /* ================================================================== *
     *  3. Set matrix values
     * ================================================================== */
    HYPRE_Int ilower[3] = { Is, Js, Ks };
    HYPRE_Int iupper[3] = { Ie - 1, Je - 1, Ke - 1 };

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
    HYPRE_StructPCGSetup(solver, hypre_A, hypre_b, hypre_x);
    HYPRE_StructPCGSolve(solver, hypre_A, hypre_b, hypre_x);

    HYPRE_Int num_iters;
    double final_res;
    HYPRE_StructPCGGetNumIterations(solver, &num_iters);
    HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res);

    /* ================================================================== *
     *  6. Extract solution 
     * ================================================================== */
    HYPRE_StructVectorGetBoxValues(hypre_x, ilower, iupper, hypre_sol);

    /* Destroy solver + preconditioner to completely prevent memory leaks */
    HYPRE_StructPCGDestroy(solver);
    HYPRE_StructPFMGDestroy(precond);

    n = 0;
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++, n++) {
                phi[k][j][i] = hypre_sol[n];
            }
        }
    }

    /* ================================================================== *
     *  7. Remove null-space from SOLUTION (pin average to 0)
     * ================================================================== */
#if NEED_REFERENCE_PRESSURE
    {
        double local_sum = 0.0, global_sum = 0.0;
        long   local_cnt = 0,   global_cnt = 0;

        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                for (int i = Is; i < Ie; i++) {
                    local_sum += phi[k][j][i];
                    local_cnt++;
                }
            }
        }

        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);
        MPI_Allreduce(&local_cnt, &global_cnt, 1, MPI_LONG,   MPI_SUM, PCW);

        const double phi_avg = global_sum / (double)global_cnt;

        for (int k = Ks; k < Ke; k++)
            for (int j = Js; j < Je; j++)
                for (int i = Is; i < Ie; i++)
                    phi[k][j][i] -= phi_avg;
    }
#endif

    Communication_update_ghost_nodes_flow_variable(
        phi, 'h', params->ghost_nodes, data_bag);

    sprintf(statement, "HYPRE PCG+PFMG converged to %g after %d iterations\n",
            final_res, (int)num_iters);
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