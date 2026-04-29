/******************************************************************************
 * VOF_IBM_wetting.c
 *
 * Characteristic moving-contact-line (MCL) extension for the diffuse-interface
 * VOF + immersed boundary path, following Liu and Ding (2015) section 2.4 and
 * the 2D / axisymmetric usage in Liu et al. (2017).
 *
 * The model extends the liquid color C_L and the liquid-gas chemical potential
 * psi_LG into a thin shell of solid cells along the contact line by tracing
 * two characteristic lines symmetric about the local solid normal.  The
 * production sampler is the quadric interpolation of Liu15 Eq. 20 from three
 * fluid-side stencil points (A, B, C) on the first cell face crossed by the
 * characteristic outside the immersed solid.  The branch is selected by Eq. 19
 * (max for theta <= 90 deg, min otherwise).
 *
 * Only owned local cells are written.  During the 2D MCL pass, scalar reads use
 * an MPI-assembled global plane of the owned C_L/C_G/C_S data, so the
 * characteristic stencil is independent of processor cuts.  There is no linear
 * / nearest-fluid fallback: a cell is updated only when Liu's characteristic
 * projection and Eq. 20 quadratic interpolation succeed.  After the per-cell
 * update, the Liu limiter C_L = min(C_L, 1 - C_S) is applied across the owned
 * domain; mass diagnostics use the centered cell measure from TwodOps.
 *
 * Implementation guards: the MCL machinery activates only when both
 * VOF_DIFFUSE and VOF_IBM are compiled in and the build is one of the two
 * supported 2D modes (TWOD_CARTESIAN or AXISYM_RZ).  The 3D and legacy code
 * paths see no behaviour change.
 ******************************************************************************/

#include "VolumeFraction.h"
#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Conc.h"
#include "Velocity.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Immersed.h"
#include "Display.h"
#include "Array.h"
#include "Cart3d.h"
#include "Lagrangian.h"
#include "VOF_DIFFUSE.h"
#include "TwodOps.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>


#if defined(VOF_DIFFUSE) && defined(VOF_IBM)

void VOF_accumulate_solid_capillary_force(Particle *p, Cart3d_bag *data_bag)
{
    /*
     * Solid capillary / contact-line reaction force is intentionally disabled
     * for the first MCL milestone.  Liu17 sections 6.3-6.4 will need this term
     * before a quantitative trajectory comparison can be claimed; until then
     * keep the accumulators zero so the resolved-particle path links cleanly.
     */
    (void)data_bag;
    DSET_ZERO(p->F_CCF, 3);
    DSET_ZERO(p->T_CCF, 3);
}

#if defined(TWOD_CARTESIAN) || defined(AXISYM_RZ)

/* -------------------------------------------------------------------------- */
/*                              MCL constants                                  */
/* -------------------------------------------------------------------------- */

static const double MCL_CS_SOLID  = 0.05;   /* C_S threshold for the solid    */
static const double MCL_IFACE_MIN = 0.005;  /* minimum interface-band C_L/C_G */
static const double MCL_IFACE_MAX = 0.995;  /* maximum interface-band C_L/C_G */
static const int    MCL_REQUIRED_GHOST_NODES = 3;
static const int    MCL_CONTACT_DEPTH_CELLS  = 1;
static const int    MCL_MAX_ITERS            = 1;
static const int    MCL_MAX_DDA_STEPS        = 32;

#define MCL_GRAD_EPS  1.0e-12
#define MCL_NORM_EPS  1.0e-12
#define MCL_TINY      1.0e-30

/* -------------------------------------------------------------------------- */
/*                          Diagnostics counter struct                         */
/* -------------------------------------------------------------------------- */

typedef enum {
    MCL_TRACE_OK = 0,
    MCL_TRACE_FAIL_NO_FLUID_FACE = 1,
    MCL_TRACE_FAIL_HALO = 2,
    MCL_TRACE_FAIL_QUAD_STENCIL = 3
} MCL_trace_status;

typedef struct {
    long ghost_cells;
    long updated_two_rays;
    long updated_one_ray;
    long fail_zero_normal;
    long fail_no_fluid_face;
    long fail_quad_stencil;
    long fail_halo;
    long fail_both_branches;
    double mass_before;
    double mass_after;
    double mass_removed_by_limiter;
    double mass_restored;
    double mass_restore_residual;
    double max_overflow;
} MCL_diagnostics;

/* -------------------------------------------------------------------------- */
/*                          Coordinate / predicate helpers                     */
/* -------------------------------------------------------------------------- */

typedef struct {
    int active;
    int nx, ny, nz;
    int ng;
    size_t size;
    double *C_S;
    double *C_L;
    double *C_G;
    double *field;
    double ***field_ptr;
} MCL_global_plane_cache;

static MCL_global_plane_cache mcl_global = {
    0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL
};

static inline double mcl_xc(const MAC_grid *g, int i) { return g->xc[i]; }
static inline double mcl_yc(const MAC_grid *g, int j) { return g->yc[j]; }

static inline int mcl_in_halo_x(const MAC_grid *g, int i)
{
    if (mcl_global.active)
        return (i >= -mcl_global.ng) && (i < g->NX + mcl_global.ng);
    return (i >= g->L_Is) && (i < g->L_Ie);
}
static inline int mcl_in_halo_y(const MAC_grid *g, int j)
{
    if (mcl_global.active)
        return (j >= -mcl_global.ng) && (j < g->NY + mcl_global.ng);
    return (j >= g->L_Js) && (j < g->L_Je);
}

static void mcl_global_cache_clear(void)
{
    mcl_global.active = 0;
    mcl_global.field_ptr = NULL;
}

static void mcl_global_cache_free(void)
{
    free(mcl_global.C_S);
    free(mcl_global.C_L);
    free(mcl_global.C_G);
    free(mcl_global.field);
    mcl_global.C_S = NULL;
    mcl_global.C_L = NULL;
    mcl_global.C_G = NULL;
    mcl_global.field = NULL;
    mcl_global.field_ptr = NULL;
    mcl_global.size = 0;
    mcl_global.nx = 0;
    mcl_global.ny = 0;
    mcl_global.nz = 0;
    mcl_global.ng = 0;
    mcl_global.active = 0;
}

static int mcl_global_cache_ensure(const MAC_grid *grid,
                                   const Parameters *params)
{
    const size_t nx = (size_t)grid->NX;
    const size_t ny = (size_t)grid->NY;
    const size_t nz = (size_t)grid->NZ;

    if (nx == 0 || ny == 0 || nz == 0)
        return 0;
    if (nx > SIZE_MAX / ny || nx * ny > SIZE_MAX / nz)
        return 0;

    const size_t n = nx * ny * nz;
    if (n > (size_t)INT_MAX)
        return 0;

    if (mcl_global.C_S != NULL &&
        mcl_global.nx == grid->NX &&
        mcl_global.ny == grid->NY &&
        mcl_global.nz == grid->NZ) {
        mcl_global.ng = params->ghost_nodes;
        return 1;
    }

    mcl_global_cache_free();

    mcl_global.C_S   = (double *)malloc(n * sizeof(double));
    mcl_global.C_L   = (double *)malloc(n * sizeof(double));
    mcl_global.C_G   = (double *)malloc(n * sizeof(double));
    mcl_global.field = (double *)malloc(n * sizeof(double));

    if (mcl_global.C_S == NULL || mcl_global.C_L == NULL ||
        mcl_global.C_G == NULL || mcl_global.field == NULL) {
        mcl_global_cache_free();
        return 0;
    }

    mcl_global.nx = grid->NX;
    mcl_global.ny = grid->NY;
    mcl_global.nz = grid->NZ;
    mcl_global.ng = params->ghost_nodes;
    mcl_global.size = n;
    return 1;
}

static inline size_t mcl_global_index(const MAC_grid *grid, int i, int j, int k)
{
    return ((size_t)k * (size_t)grid->NY + (size_t)j) *
           (size_t)grid->NX + (size_t)i;
}

static int mcl_global_map_x(const MAC_grid *grid, int i, int *mapped)
{
    const int NX = grid->NX;

#ifdef XPERIODIC
    int ii = i % NX;
    if (ii < 0) ii += NX;
    *mapped = ii;
    return 1;
#else
    if (i >= 0 && i < NX) {
        *mapped = i;
        return 1;
    }

    if (i < 0) {
#ifdef AXISYM_RZ
        *mapped = 0;
        return 1;
#elif defined(LEFT_WALL_VELOCITY_FREESLIP)
        int ii = -i - 1;
        if (ii >= 0 && ii < NX) {
            *mapped = ii;
            return 1;
        }
#else
        *mapped = 0;
        return 1;
#endif
    }

    if (i >= NX) {
#ifdef RIGHT_WALL_VELOCITY_FREESLIP
        int ii = 2 * NX - i - 1;
        if (ii >= 0 && ii < NX) {
            *mapped = ii;
            return 1;
        }
#else
        *mapped = NX - 1;
        return 1;
#endif
    }
#endif

    return 0;
}

static int mcl_global_map_y(const MAC_grid *grid, int j, int *mapped)
{
    const int NY = grid->NY;

#ifdef YPERIODIC
    int jj = j % NY;
    if (jj < 0) jj += NY;
    *mapped = jj;
    return 1;
#else
    if (j >= 0 && j < NY) {
        *mapped = j;
        return 1;
    }

    if (j < 0) {
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
        int jj = -j - 1;
        if (jj >= 0 && jj < NY) {
            *mapped = jj;
            return 1;
        }
#else
        *mapped = 0;
        return 1;
#endif
    }

    if (j >= NY) {
#ifdef TOP_WALL_VELOCITY_FREESLIP
        int jj = 2 * NY - j - 1;
        if (jj >= 0 && jj < NY) {
            *mapped = jj;
            return 1;
        }
#else
        *mapped = NY - 1;
        return 1;
#endif
    }
#endif

    return 0;
}

static int mcl_global_map_z(const MAC_grid *grid, int k, int *mapped)
{
    const int NZ = grid->NZ;

#ifdef ZPERIODIC
    int kk = k % NZ;
    if (kk < 0) kk += NZ;
    *mapped = kk;
    return 1;
#else
    if (k < 0) {
        *mapped = 0;
        return 1;
    }
    if (k >= NZ) {
        *mapped = NZ - 1;
        return 1;
    }
    *mapped = k;
    return 1;
#endif
}

static int mcl_global_value(const MAC_grid *grid,
                            const double *buf,
                            int i, int j, int k,
                            double *value)
{
    int ii, jj, kk;

    if (buf == NULL || !mcl_global.active)
        return 0;
    if (!mcl_in_halo_x(grid, i) || !mcl_in_halo_y(grid, j))
        return 0;
    if (!mcl_global_map_x(grid, i, &ii) ||
        !mcl_global_map_y(grid, j, &jj) ||
        !mcl_global_map_z(grid, k, &kk))
        return 0;

    *value = buf[mcl_global_index(grid, ii, jj, kk)];
    return 1;
}

static int mcl_global_refresh(Cart3d_bag *db, double ***field)
{
    MAC_grid       *grid = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof  = db->vof;

    mcl_global_cache_clear();

    if (!mcl_global_cache_ensure(grid, params))
        return 0;

    memset(mcl_global.C_S,   0, mcl_global.size * sizeof(double));
    memset(mcl_global.C_L,   0, mcl_global.size * sizeof(double));
    memset(mcl_global.C_G,   0, mcl_global.size * sizeof(double));
    memset(mcl_global.field, 0, mcl_global.size * sizeof(double));

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
        for (int j = grid->G_Js; j < grid->G_Je; ++j)
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                const size_t idx = mcl_global_index(grid, i, j, k);
                mcl_global.C_S[idx]   = vof->C_S[k][j][i];
                mcl_global.C_L[idx]   = vof->C_L[k][j][i];
                mcl_global.C_G[idx]   = vof->C_G[k][j][i];
                mcl_global.field[idx] = field[k][j][i];
            }

    const int count = (int)mcl_global.size;
    MPI_Allreduce(MPI_IN_PLACE, mcl_global.C_S,   count, MPI_DOUBLE, MPI_SUM, PCW);
    MPI_Allreduce(MPI_IN_PLACE, mcl_global.C_L,   count, MPI_DOUBLE, MPI_SUM, PCW);
    MPI_Allreduce(MPI_IN_PLACE, mcl_global.C_G,   count, MPI_DOUBLE, MPI_SUM, PCW);
    MPI_Allreduce(MPI_IN_PLACE, mcl_global.field, count, MPI_DOUBLE, MPI_SUM, PCW);

    mcl_global.field_ptr = field;
    mcl_global.active = 1;
    return 1;
}

static int mcl_read_C_S(Cart3d_bag *db, int i, int j, int k, double *value)
{
    if (mcl_global.active)
        return mcl_global_value(db->grid, mcl_global.C_S, i, j, k, value);

    if (!mcl_in_halo_x(db->grid, i) || !mcl_in_halo_y(db->grid, j))
        return 0;
    *value = db->vof->C_S[k][j][i];
    return 1;
}

static int mcl_read_C_L(Cart3d_bag *db, int i, int j, int k, double *value)
{
    if (mcl_global.active)
        return mcl_global_value(db->grid, mcl_global.C_L, i, j, k, value);

    if (!mcl_in_halo_x(db->grid, i) || !mcl_in_halo_y(db->grid, j))
        return 0;
    *value = db->vof->C_L[k][j][i];
    return 1;
}

static int mcl_read_C_G(Cart3d_bag *db, int i, int j, int k, double *value)
{
    if (mcl_global.active)
        return mcl_global_value(db->grid, mcl_global.C_G, i, j, k, value);

    if (!mcl_in_halo_x(db->grid, i) || !mcl_in_halo_y(db->grid, j))
        return 0;
    *value = db->vof->C_G[k][j][i];
    return 1;
}

static int mcl_read_field(Cart3d_bag *db,
                          double ***field,
                          int i, int j, int k,
                          double *value)
{
    if (mcl_global.active) {
        if (field == db->vof->C_S)
            return mcl_global_value(db->grid, mcl_global.C_S, i, j, k, value);
        if (field == db->vof->C_L)
            return mcl_global_value(db->grid, mcl_global.C_L, i, j, k, value);
        if (field == db->vof->C_G)
            return mcl_global_value(db->grid, mcl_global.C_G, i, j, k, value);
        if (field == mcl_global.field_ptr)
            return mcl_global_value(db->grid, mcl_global.field, i, j, k, value);
    }

    if (!mcl_in_halo_x(db->grid, i) || !mcl_in_halo_y(db->grid, j))
        return 0;
    *value = field[k][j][i];
    return 1;
}

/* -------------------------------------------------------------------------- */
/*   Per-rank branch cache: stores Liu Eq. 19 branch selected for C_L update.  */
/*   +1 means d1, -1 means d2, 0 means no valid characteristic.               */
/* -------------------------------------------------------------------------- */

static int *mcl_branch_cache = NULL;
static int  mcl_branch_i0 = 0, mcl_branch_j0 = 0, mcl_branch_k0 = 0;
static int  mcl_branch_ni = 0, mcl_branch_nj = 0, mcl_branch_nk = 0;
static size_t mcl_branch_size = 0;

static int mcl_branch_prepare(Cart3d_bag *db)
{
    MAC_grid *grid = db->grid;

    int ni = grid->L_Ie - grid->L_Is;
    int nj = grid->L_Je - grid->L_Js;
    int nk = grid->L_Ke - grid->L_Ks;

    if (ni <= 0 || nj <= 0 || nk <= 0)
        return 0;

    size_t n = (size_t)ni * (size_t)nj * (size_t)nk;

    if (mcl_branch_cache == NULL ||
        mcl_branch_size != n ||
        mcl_branch_i0 != grid->L_Is ||
        mcl_branch_j0 != grid->L_Js ||
        mcl_branch_k0 != grid->L_Ks ||
        mcl_branch_ni != ni ||
        mcl_branch_nj != nj ||
        mcl_branch_nk != nk) {

        int *tmp = (int *)realloc(mcl_branch_cache, n * sizeof(int));
        if (tmp == NULL) {
            free(mcl_branch_cache);
            mcl_branch_cache = NULL;
            mcl_branch_size = 0;
            return 0;
        }

        mcl_branch_cache = tmp;
        mcl_branch_size  = n;
        mcl_branch_i0    = grid->L_Is;
        mcl_branch_j0    = grid->L_Js;
        mcl_branch_k0    = grid->L_Ks;
        mcl_branch_ni    = ni;
        mcl_branch_nj    = nj;
        mcl_branch_nk    = nk;
    }

    memset(mcl_branch_cache, 0, mcl_branch_size * sizeof(int));
    return 1;
}

static inline int mcl_branch_index(Cart3d_bag *db, int i, int j, int k, size_t *idx)
{
    MAC_grid *grid = db->grid;

    if (mcl_branch_cache == NULL)
        return 0;
    if (i < grid->L_Is || i >= grid->L_Ie ||
        j < grid->L_Js || j >= grid->L_Je ||
        k < grid->L_Ks || k >= grid->L_Ke)
        return 0;

    size_t ii = (size_t)(i - mcl_branch_i0);
    size_t jj = (size_t)(j - mcl_branch_j0);
    size_t kk = (size_t)(k - mcl_branch_k0);

    *idx = (kk * (size_t)mcl_branch_nj + jj) * (size_t)mcl_branch_ni + ii;
    return (*idx < mcl_branch_size);
}

static inline void mcl_branch_set(Cart3d_bag *db, int i, int j, int k, int branch)
{
    size_t idx;
    if (mcl_branch_index(db, i, j, k, &idx))
        mcl_branch_cache[idx] = branch;
}

static inline int mcl_branch_get(Cart3d_bag *db, int i, int j, int k)
{
    size_t idx;
    if (mcl_branch_index(db, i, j, k, &idx))
        return mcl_branch_cache[idx];
    return 0;
}

static inline int mcl_is_solid_cell(Cart3d_bag *db, int i, int j, int k)
{
    double cs = 0.0;
    return mcl_read_C_S(db, i, j, k, &cs) && cs >= MCL_CS_SOLID;
}

static inline int mcl_is_fluid_cell(Cart3d_bag *db, int i, int j, int k)
{
    double cs = 0.0;
    return mcl_read_C_S(db, i, j, k, &cs) && cs < MCL_CS_SOLID;
}

static inline int mcl_is_lg_band(double cl, double cg)
{
    return (cl > MCL_IFACE_MIN) && (cl < MCL_IFACE_MAX) &&
           (cg > MCL_IFACE_MIN) && (cg < MCL_IFACE_MAX);
}

static int mcl_is_left_symmetry(const MAC_grid *grid, const Parameters *params)
{
#if defined(TWOD_CARTESIAN) && defined(LEFT_WALL_VELOCITY_FREESLIP)
    (void)params;
    return grid->G_Is == 0;
#else
    (void)grid;
    (void)params;
    return 0;
#endif
}

static int mcl_is_axis_left(const MAC_grid *grid)
{
#ifdef AXISYM_RZ
    return grid->G_Is == 0;
#else
    (void)grid;
    return 0;
#endif
}

/* -------------------------------------------------------------------------- */
/*                MCL-specific normal-vector boundary parity helper            */
/* -------------------------------------------------------------------------- */

static void mcl_apply_normal_parity(Cart3d_bag *db,
                                    double ***nx,
                                    double ***ny,
                                    double ***nz)
{
    MAC_grid   *grid   = db->grid;
    Parameters *params = db->params;

    const int NG  = params->ghost_nodes;
    const int NX  = grid->NX;
    const int NY  = grid->NY;

    const int IsL = grid->L_Is, IeL = grid->L_Ie;
    const int JsL = grid->L_Js, JeL = grid->L_Je;
    const int KsL = grid->L_Ks, KeL = grid->L_Ke;

    /* Z-slab values are bookkeeping; force out-of-plane normal to zero in 2D. */
    if (nz != NULL) {
        for (int k = KsL; k < KeL; ++k)
            for (int j = JsL; j < JeL; ++j)
                for (int i = IsL; i < IeL; ++i)
                    nz[k][j][i] = 0.0;
    }

#ifdef AXISYM_RZ
    if (grid->G_Is == 0) {
        int i0 = 0;
        for (int k = KsL; k < KeL; ++k)
            for (int j = JsL; j < JeL; ++j) {
                /* Force radial component to zero on the axis. */
                nx[k][j][i0] = 0.0;
                /* Odd radial parity, even axial parity across the axis. */
                for (int g = 1; g <= NG; ++g) {
                    int src = i0 + g - 1;
                    int dst = i0 - g;
                    nx[k][j][dst] = -nx[k][j][src];
                    ny[k][j][dst] =  ny[k][j][src];
                }
            }
    }
#else
    if (mcl_is_left_symmetry(grid, params)) {
        int i0 = 0;
        for (int k = KsL; k < KeL; ++k)
            for (int j = JsL; j < JeL; ++j)
                for (int g = 1; g <= NG; ++g) {
                    int src = i0 + g - 1;
                    int dst = i0 - g;
                    nx[k][j][dst] = -nx[k][j][src];
                    ny[k][j][dst] =  ny[k][j][src];
                }
    }
#endif

#if defined(RIGHT_WALL_VELOCITY_FREESLIP) && !defined(XPERIODIC)
    if (grid->G_Ie == NX) {
        int i1 = NX - 1;
        for (int k = KsL; k < KeL; ++k)
            for (int j = JsL; j < JeL; ++j)
                for (int g = 0; g < NG; ++g) {
                    int src = i1 - 1 - g;
                    int dst = i1 + g;
                    nx[k][j][dst] = -nx[k][j][src];
                    ny[k][j][dst] =  ny[k][j][src];
                }
    }
#endif

#if defined(BOTTOM_WALL_VELOCITY_FREESLIP) && !defined(YPERIODIC)
    if (grid->G_Js == 0) {
        int j0 = 0;
        for (int k = KsL; k < KeL; ++k)
            for (int i = IsL; i < IeL; ++i)
                for (int g = 1; g <= NG; ++g) {
                    int src = j0 + g - 1;
                    int dst = j0 - g;
                    nx[k][dst][i] =  nx[k][src][i];
                    ny[k][dst][i] = -ny[k][src][i];
                }
    }
#endif

#if defined(TOP_WALL_VELOCITY_FREESLIP) && !defined(YPERIODIC)
    if (grid->G_Je == NY) {
        int j1 = NY - 1;
        for (int k = KsL; k < KeL; ++k)
            for (int i = IsL; i < IeL; ++i)
                for (int g = 0; g < NG; ++g) {
                    int src = j1 - 1 - g;
                    int dst = j1 + g;
                    nx[k][dst][i] =  nx[k][src][i];
                    ny[k][dst][i] = -ny[k][src][i];
                }
    }
#endif

    /* Defer to the scalar / volume-fraction halo for non-symmetry sides. */
    Communication_update_ghost_nodes_flow_variable(nx, VOLUME_FRACTION, NG, db);
    Communication_update_ghost_nodes_flow_variable(ny, VOLUME_FRACTION, NG, db);
    if (nz != NULL)
        Communication_update_ghost_nodes_flow_variable(nz, VOLUME_FRACTION, NG, db);
}

/* -------------------------------------------------------------------------- */
/*                Solid normals from the diffuse C_S field, Liu Eq. 17         */
/* -------------------------------------------------------------------------- */

void VOF_DIFFUSE_compute_solid_normals_MCL(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    VOF_DIFFUSE_set_boundary_values(vof->C_S, db);

    const int Is = grid->L_Is, Ie = grid->L_Ie;
    const int Js = grid->L_Js, Je = grid->L_Je;
    const int Ks = grid->L_Ks, Ke = grid->L_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            for (int i = Is; i < Ie; ++i) {
                double gx, gy;

                if (i > grid->L_Is && i < grid->L_Ie - 1) {
                    gx = (vof->C_S[k][j][i+1] - vof->C_S[k][j][i-1]) *
                         grid->i2dx_c[i];
                } else if (i == grid->L_Is) {
                    gx = (vof->C_S[k][j][i+1] - vof->C_S[k][j][i]) *
                         grid->idx_c[i];
                } else {
                    gx = (vof->C_S[k][j][i] - vof->C_S[k][j][i-1]) *
                         grid->idx_c[i-1];
                }

                if (j > grid->L_Js && j < grid->L_Je - 1) {
                    gy = (vof->C_S[k][j+1][i] - vof->C_S[k][j-1][i]) *
                         grid->i2dy_c[j];
                } else if (j == grid->L_Js) {
                    gy = (vof->C_S[k][j+1][i] - vof->C_S[k][j][i]) *
                         grid->idy_c[j];
                } else {
                    gy = (vof->C_S[k][j][i] - vof->C_S[k][j-1][i]) *
                         grid->idy_c[j-1];
                }

                double mag = sqrt(gx * gx + gy * gy);
                if (mag < MCL_GRAD_EPS) {
                    vof->nx_IBM[k][j][i] = 0.0;
                    vof->ny_IBM[k][j][i] = 0.0;
                    vof->nz_IBM[k][j][i] = 0.0;
                } else {
                    /* n_s points from solid into fluid: -grad(C_S)/|grad(C_S)|. */
                    vof->nx_IBM[k][j][i] = -gx / mag;
                    vof->ny_IBM[k][j][i] = -gy / mag;
                    vof->nz_IBM[k][j][i] = 0.0;
                }
            }
        }
    }

#ifdef AXISYM_RZ
    if (grid->G_Is == 0) {
        int i0 = 0;
        for (int k = Ks; k < Ke; ++k)
            for (int j = Js; j < Je; ++j)
                vof->nx_IBM[k][j][i0] = 0.0;
    }
#else
    (void)params;
#endif

    mcl_apply_normal_parity(db, vof->nx_IBM, vof->ny_IBM, vof->nz_IBM);
}

/* -------------------------------------------------------------------------- */
/*                       Characteristic direction helper                       */
/* -------------------------------------------------------------------------- */

static int mcl_characteristic_dirs(double nx, double ny, double theta,
                                   double d1[2], double d2[2])
{
    double mag = sqrt(nx * nx + ny * ny);
    if (mag < MCL_NORM_EPS)
        return 0;

    double n0 = nx / mag;
    double n1 = ny / mag;
    double t0 = -n1;       /* tau = (-ny, nx) */
    double t1 =  n0;

    double ct = cos(theta);
    double st = sin(theta);

    d1[0] =  ct * t0 + st * n0;
    d1[1] =  ct * t1 + st * n1;
    d2[0] = -ct * t0 + st * n0;
    d2[1] = -ct * t1 + st * n1;

    double l1 = sqrt(d1[0] * d1[0] + d1[1] * d1[1]);
    double l2 = sqrt(d2[0] * d2[0] + d2[1] * d2[1]);
    if (l1 < MCL_NORM_EPS || l2 < MCL_NORM_EPS)
        return 0;

    d1[0] /= l1; d1[1] /= l1;
    d2[0] /= l2; d2[1] /= l2;
    return 1;
}

/* -------------------------------------------------------------------------- */
/*                       Quadratic interpolation, Liu Eq. 20                   */
/* -------------------------------------------------------------------------- */

static double mcl_quad_interp_eq20(double fA, double fB, double fC,
                                   double lAD, double h)
{
    if (h <= MCL_TINY)
        return fA;
    double r = lAD / h;
    return r * r * (0.5 * fA - fB + 0.5 * fC)
         + r       * (-1.5 * fA + 2.0 * fB - 0.5 * fC)
         + fA;
}

/* -------------------------------------------------------------------------- */
/*       Grid-DDA ray tracer: locate D as the first cell face whose owning     */
/*       cell is fluid-side (C_S < 0.05).  Returns 1 on success.               */
/* -------------------------------------------------------------------------- */

typedef struct {
    int    cell_i;     /* fluid cell containing the entry face        */
    int    cell_j;
    int    face_axis;  /* 0 -> vertical face (x = const), 1 -> horiz. */
    double xD;
    double yD;
} MCL_hit;

static int mcl_bracket_y(const MAC_grid *grid, double yD, int *jlo, int *jhi)
{
    const double eps = 1.0e-12;
    const int Js = mcl_global.active ? -mcl_global.ng : grid->L_Js;
    const int Je = mcl_global.active ? grid->NY + mcl_global.ng : grid->L_Je;

    for (int j = Js; j < Je - 1; ++j) {
        double y0 = grid->yc[j];
        double y1 = grid->yc[j + 1];
        double ymin = (y0 < y1) ? y0 : y1;
        double ymax = (y0 > y1) ? y0 : y1;
        if (yD >= ymin - eps && yD <= ymax + eps) {
            *jlo = j;
            *jhi = j + 1;
            return 1;
        }
    }

    return 0;
}

static int mcl_bracket_x(const MAC_grid *grid, double xD, int *ilo, int *ihi)
{
    const double eps = 1.0e-12;
    const int Is = mcl_global.active ? -mcl_global.ng : grid->L_Is;
    const int Ie = mcl_global.active ? grid->NX + mcl_global.ng : grid->L_Ie;

    for (int i = Is; i < Ie - 1; ++i) {
        double x0 = grid->xc[i];
        double x1 = grid->xc[i + 1];
        double xmin = (x0 < x1) ? x0 : x1;
        double xmax = (x0 > x1) ? x0 : x1;
        if (xD >= xmin - eps && xD <= xmax + eps) {
            *ilo = i;
            *ihi = i + 1;
            return 1;
        }
    }

    return 0;
}



static int mcl_point_valid_for_quad(Cart3d_bag *db, int i, int j, int k)
{
    MAC_grid       *grid = db->grid;

    if (!mcl_in_halo_x(grid, i) || !mcl_in_halo_y(grid, j))
        return 0;

    if (mcl_is_solid_cell(db, i, j, k))
        return 0;

    return 1;
}


static int mcl_try_quad_stencil(
    Cart3d_bag *db,
    double ***field,
    int iA, int jA,
    int iB, int jB,
    int iC, int jC,
    int k,
    double xD,
    double yD,
    int axis,
    double *value,
    double *solid_score
)
{
    MAC_grid       *grid = db->grid;

    if (!mcl_point_valid_for_quad(db, iA, jA, k)) return 0;
    if (!mcl_point_valid_for_quad(db, iB, jB, k)) return 0;
    if (!mcl_point_valid_for_quad(db, iC, jC, k)) return 0;

    double h, lAD;

    if (axis == 0) {
        h = fabs(grid->yc[jB] - grid->yc[jA]);
        lAD = fabs(yD - grid->yc[jA]);
    } else {
        h = fabs(grid->xc[iB] - grid->xc[iA]);
        lAD = fabs(xD - grid->xc[iA]);
    }

    if (h <= MCL_TINY) return 0;
    if (lAD < -1e-12 || lAD > h + 1e-12) return 0;

    double fA, fB, fC;
    if (!mcl_read_field(db, field, iA, jA, k, &fA)) return 0;
    if (!mcl_read_field(db, field, iB, jB, k, &fB)) return 0;
    if (!mcl_read_field(db, field, iC, jC, k, &fC)) return 0;

    *value = mcl_quad_interp_eq20(fA, fB, fC, lAD, h);

    /*
       Liu: choose C nearer the solid boundary.
       Outside the solid, C_S < 0.05, so larger C_S means closer to C_S=0.05.
    */
    if (!mcl_read_C_S(db, iC, jC, k, solid_score)) return 0;

    return 1;
}

static int mcl_sample_at_hit(
    Cart3d_bag *db,
    double ***field,
    const MCL_hit *hit,
    int k0,
    double *out_value,
    int *out_branch_dir
)
{
    (void)out_branch_dir;

    MAC_grid *grid = db->grid;

    double v1 = 0.0, v2 = 0.0;
    double s1 = -DBL_MAX, s2 = -DBL_MAX;
    int ok1 = 0;
    int ok2 = 0;

    if (hit->face_axis == 0) {
        /*
         * D lies on a vertical scalar mesh line x = xc[i_line].  Liu Eq. 20
         * interpolates along that same mesh line using scalar points A,B,C.
         */
        int i_line = hit->cell_i;
        int jlo, jhi;

        if (!mcl_in_halo_x(grid, i_line))
            return 0;
        if (!mcl_bracket_y(grid, hit->yD, &jlo, &jhi))
            return 0;

        /* C on the upper extension of AB. */
        ok1 = mcl_try_quad_stencil(
            db, field,
            i_line, jlo,
            i_line, jhi,
            i_line, jhi + 1,
            k0,
            hit->xD, hit->yD,
            hit->face_axis,
            &v1, &s1
        );

        /* C on the lower extension of AB. */
        ok2 = mcl_try_quad_stencil(
            db, field,
            i_line, jhi,
            i_line, jlo,
            i_line, jlo - 1,
            k0,
            hit->xD, hit->yD,
            hit->face_axis,
            &v2, &s2
        );
    } else {
        /*
         * D lies on a horizontal scalar mesh line y = yc[j_line].
         * Interpolate along x on the scalar line.
         */
        int j_line = hit->cell_j;
        int ilo, ihi;

        if (!mcl_in_halo_y(grid, j_line))
            return 0;
        if (!mcl_bracket_x(grid, hit->xD, &ilo, &ihi))
            return 0;

        /* C on the right extension of AB. */
        ok1 = mcl_try_quad_stencil(
            db, field,
            ilo, j_line,
            ihi, j_line,
            ihi + 1, j_line,
            k0,
            hit->xD, hit->yD,
            hit->face_axis,
            &v1, &s1
        );

        /* C on the left extension of AB. */
        ok2 = mcl_try_quad_stencil(
            db, field,
            ihi, j_line,
            ilo, j_line,
            ilo - 1, j_line,
            k0,
            hit->xD, hit->yD,
            hit->face_axis,
            &v2, &s2
        );
    }

    if (!ok1 && !ok2)
        return 0;

    if (ok1 && ok2) {
        /* Liu: if two C choices exist, keep the one nearer C_S=0.05. */
        *out_value = (s1 >= s2) ? v1 : v2;
        return 1;
    }

    *out_value = ok1 ? v1 : v2;
    return 1;
}

/* -------------------------------------------------------------------------- */
/*       Sample C_L (and remember the branch) along one characteristic.        */
/* -------------------------------------------------------------------------- */

typedef struct {
    int    valid;
    double sample;
    MCL_hit hit;
} MCL_branch_result;

static int mcl_trace_branch_by_scalar_lines(Cart3d_bag *db,
                                            double ***field,
                                            int i0, int j0, int k0,
                                            const double dir[2],
                                            int max_steps,
                                            MCL_branch_result *out,
                                            MCL_trace_status *status)
{
    MAC_grid *grid = db->grid;

    if (out != NULL) {
        out->valid = 0;
        out->sample = 0.0;
    }
    if (status != NULL)
        *status = MCL_TRACE_FAIL_NO_FLUID_FACE;

    double x = mcl_xc(grid, i0);
    double y = mcl_yc(grid, j0);
    int il = i0;
    int jl = j0;
    int saw_quad_candidate = 0;

    for (int step = 0; step < max_steps; ++step) {
        double tx = 1.0e300;
        double ty = 1.0e300;
        int next_i = il;
        int next_j = jl;

        if (dir[0] > MCL_TINY) {
            next_i = il + 1;
            if (!mcl_in_halo_x(grid, next_i)) {
                if (status != NULL) *status = MCL_TRACE_FAIL_HALO;
                return 0;
            }
            tx = (grid->xc[next_i] - x) / dir[0];
        } else if (dir[0] < -MCL_TINY) {
            next_i = il - 1;
            if (!mcl_in_halo_x(grid, next_i)) {
                if (status != NULL) *status = MCL_TRACE_FAIL_HALO;
                return 0;
            }
            tx = (grid->xc[next_i] - x) / dir[0];
        }

        if (dir[1] > MCL_TINY) {
            next_j = jl + 1;
            if (!mcl_in_halo_y(grid, next_j)) {
                if (status != NULL) *status = MCL_TRACE_FAIL_HALO;
                return 0;
            }
            ty = (grid->yc[next_j] - y) / dir[1];
        } else if (dir[1] < -MCL_TINY) {
            next_j = jl - 1;
            if (!mcl_in_halo_y(grid, next_j)) {
                if (status != NULL) *status = MCL_TRACE_FAIL_HALO;
                return 0;
            }
            ty = (grid->yc[next_j] - y) / dir[1];
        }

        if (tx <= MCL_TINY) tx = 1.0e300;
        if (ty <= MCL_TINY) ty = 1.0e300;

        if (tx >= 1.0e299 && ty >= 1.0e299) {
            if (status != NULL) *status = MCL_TRACE_FAIL_NO_FLUID_FACE;
            return 0;
        }

        MCL_hit hit;
        double tmin;
        if (tx <= ty) {
            tmin = tx;
            hit.cell_i = next_i;
            hit.cell_j = jl;
            hit.face_axis = 0;
            il = next_i;
        } else {
            tmin = ty;
            hit.cell_i = il;
            hit.cell_j = next_j;
            hit.face_axis = 1;
            jl = next_j;
        }

        x = x + tmin * dir[0];
        y = y + tmin * dir[1];
        hit.xD = x;
        hit.yD = y;

        saw_quad_candidate = 1;

        double value = 0.0;
        if (mcl_sample_at_hit(db, field, &hit, k0, &value, NULL)) {
            if (out != NULL) {
                out->valid = 1;
                out->sample = value;
                out->hit = hit;
            }
            if (status != NULL) *status = MCL_TRACE_OK;
            return 1;
        }
    }

    if (status != NULL)
        *status = saw_quad_candidate ? MCL_TRACE_FAIL_QUAD_STENCIL
                                     : MCL_TRACE_FAIL_NO_FLUID_FACE;
    return 0;
}


static int mcl_trace_branch(Cart3d_bag *db, double ***field,
                            int i0, int j0, int k0,
                            const double dir[2],
                            int max_steps,
                            MCL_branch_result *out,
                            MCL_trace_status *status)
{
    return mcl_trace_branch_by_scalar_lines(db, field, i0, j0, k0,
                                            dir, max_steps, out, status);
}


/* -------------------------------------------------------------------------- */
/*                        Eq. 18 ghost contact-line predicate                  */
/* -------------------------------------------------------------------------- */

static int mcl_is_ghost_contact_cell(Cart3d_bag *db, int i, int j, int k)
{
    MAC_grid       *grid = db->grid;

    if (!mcl_is_solid_cell(db, i, j, k))
        return 0;

    /*
     * The cell itself in the Eq. 18 band, OR any of its four planar neighbours
     * is fluid-side and in the band.  We scan both signs of i/j to remove the
     * orientation dependence in the original paper notation.
     */
    double cl = 0.0, cg = 0.0;
    if (mcl_read_C_L(db, i, j, k, &cl) &&
        mcl_read_C_G(db, i, j, k, &cg) &&
        mcl_is_lg_band(cl, cg))
        return 1;

    static const int dij[4][2] = { {1, 0}, {-1, 0}, {0, 1}, {0, -1} };
    for (int s = 0; s < 4; ++s) {
        int ni = i + dij[s][0];
        int nj = j + dij[s][1];

        if (!mcl_in_halo_x(grid, ni) || !mcl_in_halo_y(grid, nj))
            continue;
        if (mcl_is_solid_cell(db, ni, nj, k))
            continue;
        if (mcl_read_C_L(db, ni, nj, k, &cl) &&
            mcl_read_C_G(db, ni, nj, k, &cg) &&
            mcl_is_lg_band(cl, cg))
            return 1;
    }

    return 0;
}

/* -------------------------------------------------------------------------- */
/*                Mass diagnostic: integral of C_L over fluid region           */
/* -------------------------------------------------------------------------- */

static double mcl_mass_integral(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    double local = 0.0;
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
        for (int j = grid->G_Js; j < grid->G_Je; ++j)
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                if (vof->C_S[k][j][i] >= MCL_CS_SOLID) continue;
                double measure = TwodOps_cell_measure_c(grid, params, i, j, k);
                local += vof->C_L[k][j][i] * measure;
            }

    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, PCW);
    return global;
}

static double mcl_redistribute_liquid_mass_delta(Cart3d_bag *db, double delta)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    const double tol          = 1.0e-14;
    const double tiny         = 1.0e-30;
    const double relaxed_min  = 1.0e-4;
    const int    max_sweeps   = 4;

    if (fabs(delta) < tol)
        return 0.0;

    /*
     * Important:
     *   We do NOT redistribute into the whole non-solid domain.
     *   The correction is restricted to the diffuse liquid-gas interface.
     *
     * pass = 0: strict interface band
     * pass = 1: slightly relaxed interface band, but still excludes pure gas
     *           and pure liquid.
     */
    for (int pass = 0; pass < 2 && fabs(delta) > tol; ++pass) {

        for (int sweep = 0; sweep < max_sweeps && fabs(delta) > tol; ++sweep) {

            double local_weight = 0.0;

            for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
            for (int j = grid->G_Js; j < grid->G_Je; ++j)
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

                const double cs = vof->C_S[k][j][i];
                if (cs >= MCL_CS_SOLID)
                    continue;

                const double cl = vof->C_L[k][j][i];

                int eligible = 0;
                if (pass == 0) {
                    eligible = (cl > MCL_IFACE_MIN && cl < 1.0 - MCL_IFACE_MIN);
                } else {
                    eligible = (cl > relaxed_min && cl < 1.0 - relaxed_min);
                }

                if (!eligible)
                    continue;

                const double cap = (delta > 0.0) ? (1.0 - cs - cl) : cl;
                if (cap <= 0.0)
                    continue;

                const double measure = TwodOps_cell_measure_c(grid, params, i, j, k);
                if (measure <= 0.0)
                    continue;

                /*
                 * Interface-localized weight.
                 * Maximum at C_L = 0.5, zero in pure phases.
                 */
                const double w = cl * (1.0 - cl);
                if (w <= 0.0)
                    continue;

                local_weight += w * measure;
            }

            double global_weight = 0.0;
            MPI_Allreduce(&local_weight, &global_weight, 1, MPI_DOUBLE, MPI_SUM, PCW);

            if (global_weight <= tiny)
                break;

            double local_applied = 0.0;

            for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
            for (int j = grid->G_Js; j < grid->G_Je; ++j)
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

                const double cs = vof->C_S[k][j][i];
                if (cs >= MCL_CS_SOLID)
                    continue;

                const double cl = vof->C_L[k][j][i];

                int eligible = 0;
                if (pass == 0) {
                    eligible = (cl > MCL_IFACE_MIN && cl < 1.0 - MCL_IFACE_MIN);
                } else {
                    eligible = (cl > relaxed_min && cl < 1.0 - relaxed_min);
                }

                if (!eligible)
                    continue;

                const double measure = TwodOps_cell_measure_c(grid, params, i, j, k);
                if (measure <= 0.0)
                    continue;

                const double cap = (delta > 0.0) ? (1.0 - cs - cl) : cl;
                if (cap <= 0.0)
                    continue;

                const double w = cl * (1.0 - cl);
                if (w <= 0.0)
                    continue;

                /*
                 * Target mass increment in this cell.
                 */
                double dmass = delta * (w * measure) / global_weight;
                double dcl   = dmass / measure;

                /*
                 * Respect local bounds:
                 *   0 <= C_L
                 *   C_L + C_S <= 1
                 */
                if (delta > 0.0) {
                    if (dcl > cap)
                        dcl = cap;
                    if (dcl < 0.0)
                        dcl = 0.0;
                } else {
                    if (-dcl > cap)
                        dcl = -cap;
                    if (dcl > 0.0)
                        dcl = 0.0;
                }

                vof->C_L[k][j][i] = cl + dcl;
                vof->C_G[k][j][i] = clampDouble(1.0 - vof->C_L[k][j][i] - cs, 0.0, 1.0);

                local_applied += dcl * measure;
            }

            double global_applied = 0.0;
            MPI_Allreduce(&local_applied, &global_applied, 1, MPI_DOUBLE, MPI_SUM, PCW);

            delta -= global_applied;

            VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
            VOF_DIFFUSE_update_phase_cache(db);

            if (fabs(global_applied) < tol)
                break;
        }
    }

    /*
     * Nonzero return means the requested correction could not be applied
     * without polluting bulk gas/liquid or violating C_L + C_S <= 1.
     */
    if (fabs(delta) < tol)
        return 0.0;

    return delta;
}

/* -------------------------------------------------------------------------- */
/*       Apply the C_L + C_S <= 1 limiter over the owned domain (Liu15 5.2)    */
/* -------------------------------------------------------------------------- */

static void mcl_apply_full_domain_limiter(Cart3d_bag *db,
                                          MCL_diagnostics *diag)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    double local_removed = 0.0;
    double local_max_over = 0.0;

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

        double cl = vof->C_L[k][j][i];
        double cs = vof->C_S[k][j][i];

        double over = cl + cs - 1.0;
        if (over > local_max_over)
            local_max_over = over;

        if (over > 0.0) {
            double measure = TwodOps_cell_measure_c(grid, params, i, j, k);
            local_removed += over * measure;
            cl = 1.0 - cs;
        }

        if (cl < 0.0) cl = 0.0;
        if (cl > 1.0) cl = 1.0;

        vof->C_L[k][j][i] = cl;

        double cg = 1.0 - cl - cs;
        if (cg < 0.0) cg = 0.0;
        if (cg > 1.0) cg = 1.0;
        vof->C_G[k][j][i] = cg;
    }

    double global_removed = 0.0;
    double global_max_over = 0.0;

    MPI_Allreduce(&local_removed,  &global_removed,  1, MPI_DOUBLE, MPI_SUM, PCW);
    MPI_Allreduce(&local_max_over, &global_max_over, 1, MPI_DOUBLE, MPI_MAX, PCW);

    if (diag != NULL) {
        diag->mass_removed_by_limiter += global_removed;
        if (global_max_over > diag->max_overflow)
            diag->max_overflow = global_max_over;
    }

#ifdef DEBUG_MCL
	    if (global_removed > 1.0e-12) {
	        char msg[256];
	        snprintf(msg, sizeof(msg),
	                 "MCL limiter clipped %.6e liquid mass; redistribution deferred.\n",
	                 global_removed);
	        Display_progress(params, msg);
	    }
#else
    (void)params;
#endif
}

/* -------------------------------------------------------------------------- */
/*           Apply Eqs. 18-20 to C_L for one MCL fixed-point iteration         */
/* -------------------------------------------------------------------------- */

static void mcl_count_trace_fail(MCL_trace_status s1,
                                 MCL_trace_status s2,
                                 long *fail_no_fluid_face,
                                 long *fail_quad_stencil,
                                 long *fail_halo)
{
    if (s1 == MCL_TRACE_FAIL_HALO || s2 == MCL_TRACE_FAIL_HALO) {
        (*fail_halo)++;
    } else if (s1 == MCL_TRACE_FAIL_QUAD_STENCIL ||
               s2 == MCL_TRACE_FAIL_QUAD_STENCIL) {
        (*fail_quad_stencil)++;
    } else {
        (*fail_no_fluid_face)++;
    }
}

static void mcl_iter_apply_C_L(Cart3d_bag *db,
                               double theta,
                               double ***target_field,
                               MCL_diagnostics *diag)
{
    MAC_grid       *grid = db->grid;
    VolumeFraction *vof  = db->vof;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    long ghost_local = 0;
    long updated_two = 0;
    long updated_one = 0;
    long fail_zero_normal = 0;
    long fail_no_fluid_face = 0;
    long fail_quad_stencil = 0;
    long fail_halo = 0;
    long fail_both = 0;

    int max_steps = MCL_MAX_DDA_STEPS;
    if (max_steps < 8) max_steps = 8;

    const int hydrophilic = (theta <= 0.5 * M_PI);

    for (int k = Ks; k < Ke; ++k)
        for (int j = Js; j < Je; ++j)
            for (int i = Is; i < Ie; ++i) {
                if (!mcl_is_ghost_contact_cell(db, i, j, k))
                    continue;
                ghost_local++;
                mcl_branch_set(db, i, j, k, 0);

                double nx = vof->nx_IBM[k][j][i];
                double ny = vof->ny_IBM[k][j][i];
                double d1[2], d2[2];
                if (!mcl_characteristic_dirs(nx, ny, theta, d1, d2)) {
                    fail_zero_normal++;
                    fail_both++;
                    continue;
                }

                MCL_branch_result r1, r2;
                MCL_trace_status s1 = MCL_TRACE_OK;
                MCL_trace_status s2 = MCL_TRACE_OK;

                int ok1 = mcl_trace_branch(db, target_field, i, j, k,
                                           d1, max_steps, &r1, &s1);
                int ok2 = mcl_trace_branch(db, target_field, i, j, k,
                                           d2, max_steps, &r2, &s2);

                if (!ok1 && !ok2) {
                    fail_both++;
                    mcl_count_trace_fail(s1, s2,
                                         &fail_no_fluid_face,
                                         &fail_quad_stencil,
                                         &fail_halo);
                    continue;
                }

                double new_cl;
                int branch = 0;
                if (ok1 && ok2) {
                    if (hydrophilic) {
                        if (r1.sample >= r2.sample) {
                            new_cl = r1.sample;
                            branch = +1;
                        } else {
                            new_cl = r2.sample;
                            branch = -1;
                        }
                    } else {
                        if (r1.sample <= r2.sample) {
                            new_cl = r1.sample;
                            branch = +1;
                        } else {
                            new_cl = r2.sample;
                            branch = -1;
                        }
                    }
                    updated_two++;
                } else if (ok1) {
                    new_cl = r1.sample;
                    branch = +1;
                    updated_one++;
                } else {
                    new_cl = r2.sample;
                    branch = -1;
                    updated_one++;
                }

                if (new_cl < 0.0) new_cl = 0.0;
                if (new_cl > 1.0) new_cl = 1.0;
                double cs = vof->C_S[k][j][i];
                if (new_cl + cs > 1.0) new_cl = 1.0 - cs;
                if (new_cl < 0.0) new_cl = 0.0;

                vof->C_L[k][j][i] = new_cl;
                mcl_branch_set(db, i, j, k, branch);
            }

    long g[8] = { ghost_local, updated_two, updated_one,
                  fail_zero_normal, fail_no_fluid_face,
                  fail_quad_stencil, fail_halo, fail_both };
    long G[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    MPI_Allreduce(g, G, 8, MPI_LONG, MPI_SUM, PCW);

    if (diag != NULL) {
        diag->ghost_cells        += G[0];
        diag->updated_two_rays   += G[1];
        diag->updated_one_ray    += G[2];
        diag->fail_zero_normal   += G[3];
        diag->fail_no_fluid_face += G[4];
        diag->fail_quad_stencil  += G[5];
        diag->fail_halo          += G[6];
        diag->fail_both_branches += G[7];
    }
}

/* -------------------------------------------------------------------------- */
/*                       Public entry: apply MCL on C_L                        */
/* -------------------------------------------------------------------------- */

void VOF_DIFFUSE_apply_contact_angle(Cart3d_bag *db)
{
    Parameters *params = db->params;

    if (params->ghost_nodes < MCL_REQUIRED_GHOST_NODES) {
        char msg[200];
        snprintf(msg, sizeof(msg),
                 "MCL contact-angle update requires ghost_nodes >= %d (have %d).\n",
                 MCL_REQUIRED_GHOST_NODES, params->ghost_nodes);
        Display_throw_warning(msg, params);
    }

    /* Keep C_S boundary halo and the phase cache (C_G, F) in sync first. */
    VOF_DIFFUSE_set_boundary_values(db->vof->C_S, db);
    VOF_DIFFUSE_set_boundary_values(db->vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);
    VOF_DIFFUSE_set_boundary_values(db->vof->C_G, db);

    VOF_DIFFUSE_compute_solid_normals_MCL(db);

    MCL_diagnostics diag;
    diag.ghost_cells             = 0;
    diag.updated_two_rays        = 0;
    diag.updated_one_ray         = 0;
    diag.fail_zero_normal        = 0;
    diag.fail_no_fluid_face      = 0;
    diag.fail_quad_stencil       = 0;
    diag.fail_halo               = 0;
    diag.fail_both_branches      = 0;
    diag.mass_before             = mcl_mass_integral(db);
    diag.mass_after              = 0.0;
    diag.mass_removed_by_limiter = 0.0;
    diag.mass_restored           = 0.0;
    diag.mass_restore_residual   = 0.0;
    diag.max_overflow            = 0.0;

    if (!mcl_branch_prepare(db)) {
        Display_throw_warning("MCL branch cache allocation failed; psi_LG extension will reselect branches from C_L.\n",
                              params);
    }

    const double theta_deg = params->contact_angle_deg;
    const double theta = theta_deg * M_PI / 180.0;
    int warned_global_cache = 0;

    for (int iter = 0; iter < MCL_MAX_ITERS; ++iter) {
        if (mcl_branch_cache != NULL)
            memset(mcl_branch_cache, 0, mcl_branch_size * sizeof(int));

        if (!mcl_global_refresh(db, db->vof->C_L)) {
            mcl_global_cache_clear();
            if (!warned_global_cache) {
                Display_throw_warning("MCL global 2D field cache allocation failed; falling back to local halo tracing.\n",
                                      params);
                warned_global_cache = 1;
            }
        }
        mcl_iter_apply_C_L(db, theta, db->vof->C_L, &diag);
        mcl_global_cache_clear();

        mcl_apply_full_domain_limiter(db, &diag);
        VOF_DIFFUSE_set_boundary_values(db->vof->C_L, db);
        VOF_DIFFUSE_update_phase_cache(db);
        VOF_DIFFUSE_set_boundary_values(db->vof->C_G, db);
    }
    mcl_global_cache_clear();

    diag.mass_after = mcl_mass_integral(db);
    diag.mass_restored = diag.mass_before - diag.mass_after;
    if (fabs(diag.mass_restored) > 1.0e-14) {
        diag.mass_restore_residual =
            mcl_redistribute_liquid_mass_delta(db, diag.mass_restored);
        VOF_DIFFUSE_set_boundary_values(db->vof->C_L, db);
        VOF_DIFFUSE_update_phase_cache(db);
        VOF_DIFFUSE_set_boundary_values(db->vof->C_G, db);
        diag.mass_after = mcl_mass_integral(db);
    }

    /*
     * Display_progress is collective (MPI_Barrier inside) so all ranks must
     * call it; the rank check happens inside.  Diagnostic counters were
     * already MPI-reduced above, so the message is identical on every rank.
     */
    if (params->ghost_nodes >= MCL_REQUIRED_GHOST_NODES) {
        char msg[520];
        snprintf(msg, sizeof(msg),
                 "MCL: theta=%.2fdeg ghost=%ld two=%ld one=%ld "
                 "fail_both=%ld fail_normal=%ld fail_face=%ld "
                 "fail_quad=%ld fail_halo=%ld fallback=0 mass_before=%.6e "
                 "mass_after=%.6e removed=%.3e restored=%.3e "
                 "restore_residual=%.3e max_over=%.3e\n",
                 theta_deg, diag.ghost_cells, diag.updated_two_rays,
                 diag.updated_one_ray, diag.fail_both_branches,
                 diag.fail_zero_normal, diag.fail_no_fluid_face,
                 diag.fail_quad_stencil, diag.fail_halo,
                 diag.mass_before, diag.mass_after,
                 diag.mass_removed_by_limiter, diag.mass_restored,
                 diag.mass_restore_residual, diag.max_overflow);
        Display_progress(params, msg);
    }
}

/* -------------------------------------------------------------------------- */
/*           Public entry: extend psi_LG along the same characteristic         */
/* -------------------------------------------------------------------------- */

void VOF_DIFFUSE_extend_psi_LG_contact_angle(Cart3d_bag *db)
{
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    if (params->ghost_nodes < MCL_REQUIRED_GHOST_NODES)
        return;

    VOF_DIFFUSE_set_boundary_values(vof->psi_LG, db);
    VOF_DIFFUSE_set_boundary_values(vof->C_S, db);
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);
    VOF_DIFFUSE_set_boundary_values(vof->C_G, db);

    VOF_DIFFUSE_compute_solid_normals_MCL(db);

    if (!mcl_global_refresh(db, vof->psi_LG)) {
        mcl_global_cache_clear();
        Display_throw_warning("MCL psi_LG global 2D field cache allocation failed; falling back to local halo tracing.\n",
                              params);
    }

    MAC_grid *grid = db->grid;
    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    int max_steps = MCL_MAX_DDA_STEPS;
    if (max_steps < 8) max_steps = 8;

    const double theta_deg = params->contact_angle_deg;
    const double theta = theta_deg * M_PI / 180.0;
    const int hydrophilic = (theta <= 0.5 * M_PI);

    long local_updated = 0;
    long local_missing_branch = 0;
    long local_psi_fail = 0;

    for (int k = Ks; k < Ke; ++k)
        for (int j = Js; j < Je; ++j)
            for (int i = Is; i < Ie; ++i) {
                if (!mcl_is_ghost_contact_cell(db, i, j, k))
                    continue;

                double nx = vof->nx_IBM[k][j][i];
                double ny = vof->ny_IBM[k][j][i];
                double d1[2], d2[2];
                if (!mcl_characteristic_dirs(nx, ny, theta, d1, d2))
                    continue;

                int branch = mcl_branch_get(db, i, j, k);

                /*
                 * Preferred path: reuse the branch selected during the C_L MCL
                 * pass.  If this function is called before the C_L pass, recover
                 * the Liu Eq. 19 branch from C_L, but still sample psi only along
                 * that selected branch.  There is no try-other-branch fallback.
                 */
                if (branch == 0) {
                    MCL_branch_result c1, c2;
                    MCL_trace_status s1, s2;
                    int ok1 = mcl_trace_branch(db, vof->C_L, i, j, k,
                                               d1, max_steps, &c1, &s1);
                    int ok2 = mcl_trace_branch(db, vof->C_L, i, j, k,
                                               d2, max_steps, &c2, &s2);

                    if (!ok1 && !ok2) {
                        local_missing_branch++;
                        continue;
                    }

                    if (ok1 && ok2) {
                        if (hydrophilic)
                            branch = (c1.sample >= c2.sample) ? +1 : -1;
                        else
                            branch = (c1.sample <= c2.sample) ? +1 : -1;
                    } else {
                        branch = ok1 ? +1 : -1;
                    }

                    mcl_branch_set(db, i, j, k, branch);
                }

                MCL_branch_result p_sample;
                MCL_trace_status ps;
                int sampled_psi = 0;

                if (branch == +1) {
                    sampled_psi = mcl_trace_branch(db, vof->psi_LG, i, j, k,
                                                   d1, max_steps, &p_sample, &ps);
                } else if (branch == -1) {
                    sampled_psi = mcl_trace_branch(db, vof->psi_LG, i, j, k,
                                                   d2, max_steps, &p_sample, &ps);
                }

                if (sampled_psi) {
                    vof->psi_LG[k][j][i] = p_sample.sample;
                    local_updated++;
                } else {
                    local_psi_fail++;
                }
            }

    long local[3] = { local_updated, local_missing_branch, local_psi_fail };
    long global[3] = { 0, 0, 0 };
    MPI_Allreduce(local, global, 3, MPI_LONG, MPI_SUM, PCW);

#ifndef DEBUG_MCL
    (void)global;
#endif

#ifdef DEBUG_MCL
    char msg[240];
    snprintf(msg, sizeof(msg),
             "MCL psi_LG: updated=%ld missing_branch=%ld psi_fail=%ld\n",
             global[0], global[1], global[2]);
    Display_progress(params, msg);
#endif

    mcl_global_cache_clear();
    VOF_DIFFUSE_set_boundary_values(vof->psi_LG, db);
}

#else /* not (TWOD_CARTESIAN || AXISYM_RZ) */

/*
 * 3D path: the full Liu MCL needs orientation handling that Phase A does not
 * cover.  Keep these symbols as no-ops so the build stays linkable when the
 * caller is wrapped in #ifdef VOF_IBM only.
 */
void VOF_DIFFUSE_apply_contact_angle(Cart3d_bag *db)
{
    VOF_DIFFUSE_set_boundary_values(db->vof->C_S, db);
    VOF_DIFFUSE_update_phase_cache(db);
}

void VOF_DIFFUSE_extend_psi_LG_contact_angle(Cart3d_bag *db)
{
    (void)db;
}

void VOF_DIFFUSE_compute_solid_normals_MCL(Cart3d_bag *db)
{
    (void)db;
}

#endif /* TWOD_CARTESIAN || AXISYM_RZ */

#endif /* VOF_DIFFUSE && VOF_IBM */
