/*******************************************************************************
 * psolve_hypre_gpu.c
 *
 * CUDA / HYPRE-GPU pressure solver for the variable-coefficient Poisson
 * equation:
 *
 *     div( (1/rho) grad(phi) ) = rhs
 *
 * This GPU implementation intentionally keeps the same 3D Struct grid,
 * 7-point stencil, host-side algebra, row weighting, null-space handling, and
 * pressure post-processing as psolve_hypre.c.  In 2D modes the shared
 * TwodOps_build_scalar_face_coeffs() helper zeros the collapsed z couplings;
 * in 3D the full 7-point stencil remains active.  The HYPRE objects execute
 * under the device policy, and the HYPRE value buffers are unified-memory
 * allocations that are prefetched before HYPRE setup/solve and before host
 * post-processing.
 ******************************************************************************/

#if defined(USE_HYPRE) && defined(GPU_PRESSURE_SOLVER)

#include "HYPRE_struct_ls.h"
#include "HYPRE.h"
#include "HYPRE_utilities.h"
#include <cuda_runtime_api.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

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

/* Unified-memory buffers for HYPRE set/get values */
static double *hypre_vals = NULL;
static double *hypre_rhs  = NULL;
static double *hypre_sol  = NULL;
static int    *hypre_stencil_indices = NULL;

static int hypre_rank = -1;
static int hypre_cuda_device = -1;
static int hypre_device_count = 0;
static int hypre_initialized = 0;
static char hypre_device_name[256] = "unknown";

/* Mapping: stencil entry -> (di,dj,dk) offset */
static const int stencil_offsets[7][3] = {
    { 0, 0, 0},   /* center */
    {-1, 0, 0},   /* west   (-x) */
    { 1, 0, 0},   /* east   (+x) */
    { 0,-1, 0},   /* south  (-y) */
    { 0, 1, 0},   /* north  (+y) */
    { 0, 0,-1},   /* back   (-z) */
    { 0, 0, 1}    /* front  (+z) */
};

static void hypre_gpu_abort(const char *message)
{
    fprintf(stderr, "[HYPRE-GPU] rank %d: %s\n", hypre_rank, message);
    fflush(stderr);
    MPI_Abort(PCW, 1);
}

static void hypre_gpu_check_cuda(cudaError_t err, const char *call,
                                 const char *file, int line)
{
    if (err != cudaSuccess) {
        fprintf(stderr,
                "[HYPRE-GPU] rank %d: CUDA call failed at %s:%d: %s: %s\n",
                hypre_rank, file, line, call, cudaGetErrorString(err));
        fflush(stderr);
        MPI_Abort(PCW, 1);
    }
}

static void hypre_gpu_check_hypre(HYPRE_Int err, const char *call,
                                  const char *file, int line)
{
    if (err != 0) {
        fprintf(stderr,
                "[HYPRE-GPU] rank %d: HYPRE call failed at %s:%d: %s returned %d\n",
                hypre_rank, file, line, call, (int)err);
        fflush(stderr);
        MPI_Abort(PCW, 1);
    }
}

static HYPRE_Int hypre_gpu_check_hypre_allow_conv(HYPRE_Int err,
                                                  const char *call,
                                                  const char *file,
                                                  int line)
{
    if (err == 0)
        return 0;

    if ((err & ~HYPRE_ERROR_CONV) == 0) {
        if (hypre_rank == 0) {
            fprintf(stderr,
                    "[HYPRE-GPU] %s returned HYPRE_ERROR_CONV at %s:%d; "
                    "continuing with the current iterate so the outer "
                    "pressure-correction loop can retry if needed.\n",
                    call, file, line);
            fflush(stderr);
        }
        hypre_gpu_check_hypre(HYPRE_ClearError(HYPRE_ERROR_CONV),
                              "HYPRE_ClearError(HYPRE_ERROR_CONV)",
                              file, line);
        return err;
    }

    hypre_gpu_check_hypre(err, call, file, line);
    return err;
}

#define HYPRE_GPU_CHECK_CUDA(call) \
    hypre_gpu_check_cuda((call), #call, __FILE__, __LINE__)

#define HYPRE_GPU_CHECK_HYPRE(call) \
    hypre_gpu_check_hypre((call), #call, __FILE__, __LINE__)

#define HYPRE_GPU_CHECK_HYPRE_ALLOW_CONV(call) \
    hypre_gpu_check_hypre_allow_conv((call), #call, __FILE__, __LINE__)

static size_t hypre_gpu_vals_bytes(void)
{
    return (size_t)7 * (size_t)hypre_n_local * sizeof(double);
}

static size_t hypre_gpu_vec_bytes(void)
{
    return (size_t)hypre_n_local * sizeof(double);
}

static double hypre_gpu_prefetch(void *ptr, size_t bytes, int destination)
{
    double t0, t1;

    if (ptr == NULL || bytes == 0)
        return 0.0;

    t0 = MPI_Wtime();
    HYPRE_GPU_CHECK_CUDA(cudaMemPrefetchAsync(ptr, bytes, destination, 0));
    HYPRE_GPU_CHECK_CUDA(cudaDeviceSynchronize());
    t1 = MPI_Wtime();

    return t1 - t0;
}

static void hypre_gpu_bind_device(Parameters *params)
{
    MPI_Comm local_comm;
    int local_rank = 0;
    int local_size = 1;
    int proc;
    struct cudaDeviceProp prop;

    hypre_rank = params->rank;

    MPI_Comm_split_type(PCW, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &local_comm);
    MPI_Comm_rank(local_comm, &local_rank);
    MPI_Comm_size(local_comm, &local_size);

    HYPRE_GPU_CHECK_CUDA(cudaGetDeviceCount(&hypre_device_count));
    if (hypre_device_count <= 0)
        hypre_gpu_abort("GPU_PRESSURE_SOLVER was enabled, but CUDA reports no devices.");

    hypre_cuda_device = local_rank % hypre_device_count;
    HYPRE_GPU_CHECK_CUDA(cudaSetDevice(hypre_cuda_device));
    HYPRE_GPU_CHECK_CUDA(cudaGetDeviceProperties(&prop, hypre_cuda_device));
    snprintf(hypre_device_name, sizeof(hypre_device_name), "%s", prop.name);

    MPI_Comm_free(&local_comm);

    for (proc = 0; proc < params->size; proc++) {
        if (params->rank == proc) {
            printf("[HYPRE-GPU] rank %d local_rank %d/%d -> CUDA device %d/%d (%s)\n",
                   params->rank, local_rank, local_size, hypre_cuda_device,
                   hypre_device_count, hypre_device_name);
            fflush(stdout);
        }
        MPI_Barrier(PCW);
    }
}

static void hypre_gpu_print_version(Parameters *params)
{
    char *hypre_version = NULL;

    if (params->rank != 0)
        return;

    if (HYPRE_Version(&hypre_version) == 0 && hypre_version != NULL) {
        printf("[HYPRE-GPU] HYPRE version: %s\n", hypre_version);
        free(hypre_version);
    } else {
        printf("[HYPRE-GPU] HYPRE version: unavailable\n");
    }
}

static void hypre_gpu_malloc_managed(double **ptr, size_t bytes,
                                     const char *name)
{
    if (bytes == 0)
        bytes = sizeof(double);

    HYPRE_GPU_CHECK_CUDA(
        cudaMallocManaged((void **)ptr, bytes, cudaMemAttachGlobal));
    if (*ptr == NULL) {
        fprintf(stderr, "[HYPRE-GPU] rank %d: failed to allocate %s\n",
                hypre_rank, name);
        fflush(stderr);
        MPI_Abort(PCW, 1);
    }
}

/*******************************************************************************
 * Pressure_hypre_setup
 ******************************************************************************/
void Pressure_hypre_setup(Cart3d_bag *data_bag)
{
    MAC_grid   *grid   = data_bag->grid;
    Parameters *params = data_bag->params;
    int rank = params->rank;

    hypre_rank = rank;

    hypre_gpu_bind_device(params);

    HYPRE_GPU_CHECK_HYPRE(HYPRE_Initialize());
    hypre_initialized = 1;
    HYPRE_GPU_CHECK_HYPRE(HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_DeviceInitialize());
    HYPRE_GPU_CHECK_CUDA(cudaDeviceSynchronize());

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
        printf("[HYPRE-GPU] Setting up Struct Matrix and Vectors (CUDA unified memory, device execution)\n");
        printf("[HYPRE-GPU] Global grid: %d x %d x %d  (cells)\n",
               grid->NX - 1, grid->NY - 1, grid->NZ - 1);
        printf("[HYPRE-GPU] Struct layout: legacy 3D grid with 7-point stencil\n");
    }
    hypre_gpu_print_version(params);

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructGridCreate(PCW, 3, &hypre_grid));

    HYPRE_Int ilower[3] = { hypre_Is, hypre_Js, hypre_Ks };
    HYPRE_Int iupper[3] = { hypre_Ie - 1, hypre_Je - 1, hypre_Ke - 1 };
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructGridSetExtents(hypre_grid, ilower, iupper));

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
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructGridSetPeriodic(hypre_grid, periodic));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructGridAssemble(hypre_grid));

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructStencilCreate(3, 7, &hypre_stencil));
    for (int s = 0; s < 7; s++) {
        HYPRE_Int offset[3] = { stencil_offsets[s][0],
                                stencil_offsets[s][1],
                                stencil_offsets[s][2] };
        HYPRE_GPU_CHECK_HYPRE(
            HYPRE_StructStencilSetElement(hypre_stencil, s, offset));
    }

    HYPRE_GPU_CHECK_HYPRE(
        HYPRE_StructMatrixCreate(PCW, hypre_grid, hypre_stencil, &hypre_A));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructMatrixInitialize(hypre_A));

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorCreate(PCW, hypre_grid, &hypre_b));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorInitialize(hypre_b));

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorCreate(PCW, hypre_grid, &hypre_x));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorInitialize(hypre_x));

    hypre_gpu_malloc_managed(&hypre_vals, hypre_gpu_vals_bytes(), "hypre_vals");
    hypre_gpu_malloc_managed(&hypre_rhs, hypre_gpu_vec_bytes(), "hypre_rhs");
    hypre_gpu_malloc_managed(&hypre_sol, hypre_gpu_vec_bytes(), "hypre_sol");

    hypre_stencil_indices = (int *)malloc(7 * sizeof(int));
    if (hypre_stencil_indices == NULL)
        hypre_gpu_abort("failed to allocate hypre_stencil_indices");
    for (int s = 0; s < 7; s++)
        hypre_stencil_indices[s] = s;
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
    double prefetch_h2d_time = 0.0;
    double prefetch_d2h_time = 0.0;
    double cpu_assemble_time = 0.0;
    double hypre_set_values_time = 0.0;
    double hypre_setup_time = 0.0;
    double hypre_solve_time = 0.0;
    double hypre_get_values_time = 0.0;
    double postprocess_time = 0.0;

    if (hypre_cuda_device < 0)
        hypre_gpu_abort("Pressure_hypre_setup must be called before Pressure_solve_hypre.");

    HYPRE_GPU_CHECK_CUDA(cudaSetDevice(hypre_cuda_device));

    prefetch_d2h_time += hypre_gpu_prefetch(hypre_vals, hypre_gpu_vals_bytes(),
                                            cudaCpuDeviceId);
    prefetch_d2h_time += hypre_gpu_prefetch(hypre_rhs, hypre_gpu_vec_bytes(),
                                            cudaCpuDeviceId);
    prefetch_d2h_time += hypre_gpu_prefetch(hypre_sol, hypre_gpu_vec_bytes(),
                                            cudaCpuDeviceId);

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
 *  Keep the same row-weighted compatibility correction as the CPU path.
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

    prefetch_h2d_time += hypre_gpu_prefetch(hypre_vals, hypre_gpu_vals_bytes(),
                                            hypre_cuda_device);
    prefetch_h2d_time += hypre_gpu_prefetch(hypre_rhs, hypre_gpu_vec_bytes(),
                                            hypre_cuda_device);
    prefetch_h2d_time += hypre_gpu_prefetch(hypre_sol, hypre_gpu_vec_bytes(),
                                            hypre_cuda_device);

    /* ================================================================== *
     *  3. Set matrix values
     * ================================================================== */
    HYPRE_Int ilower[3] = { Is, Js, Ks };
    HYPRE_Int iupper[3] = { Ie - 1, Je - 1, Ke - 1 };

    t0 = MPI_Wtime();
    for (int s = 0; s < 7; s++) {
        HYPRE_GPU_CHECK_HYPRE(
            HYPRE_StructMatrixSetBoxValues(hypre_A, ilower, iupper, 1,
                                           &hypre_stencil_indices[s],
                                           &hypre_vals[s * hypre_n_local]));
    }
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructMatrixAssemble(hypre_A));

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorSetBoxValues(hypre_b, ilower, iupper,
                                                        hypre_rhs));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorSetBoxValues(hypre_x, ilower, iupper,
                                                        hypre_sol));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorAssemble(hypre_b));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructVectorAssemble(hypre_x));
    HYPRE_GPU_CHECK_CUDA(cudaDeviceSynchronize());
    t1 = MPI_Wtime();
    hypre_set_values_time += t1 - t0;

    /* ================================================================== *
     *  4. Create fresh solver + preconditioner
     * ================================================================== */
    HYPRE_StructSolver solver, precond;

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPCGCreate(PCW, &solver));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPCGSetMaxIter(solver, 2000));
    HYPRE_GPU_CHECK_HYPRE(
        HYPRE_StructPCGSetTol(solver, (double)params->P_CG_ETOL));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPCGSetTwoNorm(solver, 1));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPCGSetRelChange(solver, 0));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPCGSetPrintLevel(solver, 0));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPCGSetLogging(solver, 1));

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGCreate(PCW, &precond));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetMaxIter(precond, 1));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetTol(precond, 0.0));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetZeroGuess(precond));

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetRelaxType(precond, 2));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetNumPreRelax(precond, 2));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetNumPostRelax(precond, 2));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetRAPType(precond, 1));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGSetSkipRelax(precond, 0));

    HYPRE_GPU_CHECK_HYPRE(
        HYPRE_StructPCGSetPrecond(solver,
                                  HYPRE_StructPFMGSolve,
                                  HYPRE_StructPFMGSetup,
                                  precond));

    /* ================================================================== *
     *  5. Setup + Solve
     * ================================================================== */
    t0 = MPI_Wtime();
    HYPRE_GPU_CHECK_HYPRE(
        HYPRE_StructPCGSetup(solver, hypre_A, hypre_b, hypre_x));
    HYPRE_GPU_CHECK_CUDA(cudaDeviceSynchronize());
    t1 = MPI_Wtime();
    hypre_setup_time += t1 - t0;

    HYPRE_Int solve_err;
    int solve_converged;
    HYPRE_Int num_iters;
    double final_res;

    t0 = MPI_Wtime();
    solve_err = HYPRE_GPU_CHECK_HYPRE_ALLOW_CONV(
        HYPRE_StructPCGSolve(solver, hypre_A, hypre_b, hypre_x));
    HYPRE_GPU_CHECK_CUDA(cudaDeviceSynchronize());
    t1 = MPI_Wtime();
    hypre_solve_time += t1 - t0;

    solve_converged = (solve_err == 0);
    HYPRE_GPU_CHECK_HYPRE(
        HYPRE_StructPCGGetNumIterations(solver, &num_iters));
    HYPRE_GPU_CHECK_HYPRE(
        HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res));

    /* ================================================================== *
     *  6. Extract solution
     * ================================================================== */
    prefetch_h2d_time += hypre_gpu_prefetch(hypre_sol, hypre_gpu_vec_bytes(),
                                            hypre_cuda_device);

    t0 = MPI_Wtime();
    HYPRE_GPU_CHECK_HYPRE(
        HYPRE_StructVectorGetBoxValues(hypre_x, ilower, iupper, hypre_sol));
    HYPRE_GPU_CHECK_CUDA(cudaDeviceSynchronize());
    t1 = MPI_Wtime();
    hypre_get_values_time += t1 - t0;

    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPCGDestroy(solver));
    HYPRE_GPU_CHECK_HYPRE(HYPRE_StructPFMGDestroy(precond));

    prefetch_d2h_time += hypre_gpu_prefetch(hypre_sol, hypre_gpu_vec_bytes(),
                                            cudaCpuDeviceId);

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
 *  7. Remove null-space from SOLUTION with the same row weight.
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

    Pressure_apply_BCs(phi, grid, params);
    Communication_update_ghost_nodes_flow_variable(
        phi, 'h', params->ghost_nodes, data_bag);

    t1 = MPI_Wtime();
    postprocess_time += t1 - t0;

    if (solve_converged) {
        sprintf(statement,
                "[HYPRE-GPU] PCG+PFMG converged to %g after %d iterations | "
                "wall %.6e s | cpu_assembly %.6e s | um_h2d_prefetch %.6e s | "
                "um_d2h_prefetch %.6e s | set_values %.6e s | setup %.6e s | "
                "solve %.6e s | get_values %.6e s | postprocess %.6e s\n",
                final_res, (int)num_iters, MPI_Wtime() - t_total0,
                cpu_assemble_time, prefetch_h2d_time, prefetch_d2h_time,
                hypre_set_values_time, hypre_setup_time, hypre_solve_time,
                hypre_get_values_time, postprocess_time);
    } else {
        sprintf(statement,
                "[HYPRE-GPU] PCG+PFMG returned HYPRE_ERROR_CONV at %g "
                "after %d iterations (tol %g); using current iterate | "
                "wall %.6e s | cpu_assembly %.6e s | um_h2d_prefetch %.6e s | "
                "um_d2h_prefetch %.6e s | set_values %.6e s | setup %.6e s | "
                "solve %.6e s | get_values %.6e s | postprocess %.6e s\n",
                final_res, (int)num_iters, (double)params->P_CG_ETOL,
                MPI_Wtime() - t_total0, cpu_assemble_time,
                prefetch_h2d_time, prefetch_d2h_time, hypre_set_values_time,
                hypre_setup_time, hypre_solve_time, hypre_get_values_time,
                postprocess_time);
    }
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

    if (hypre_vals) HYPRE_GPU_CHECK_CUDA(cudaFree(hypre_vals));
    if (hypre_rhs)  HYPRE_GPU_CHECK_CUDA(cudaFree(hypre_rhs));
    if (hypre_sol)  HYPRE_GPU_CHECK_CUDA(cudaFree(hypre_sol));
    hypre_vals = NULL;
    hypre_rhs  = NULL;
    hypre_sol  = NULL;

    free(hypre_stencil_indices);
    hypre_stencil_indices = NULL;

    if (hypre_initialized) {
        HYPRE_Finalize();
        hypre_initialized = 0;
    }

    hypre_cuda_device = -1;
}



#endif /* USE_HYPRE && GPU_PRESSURE_SOLVER */
