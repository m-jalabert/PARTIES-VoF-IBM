#include "VOF_DIFFUSE.h"

#include "VolumeFraction.h"
#include "definitions.h"
#include "Boundary.h"
#include "Communication.h"
#include "Array.h"
#include "Display.h"
#include "Memory.h"
#include "Lagrangian.h"
#include "Interpolate.h"
#include "Particle.h"
#include "TwodOps.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* DIFFUSE_SOLID_MASS_CUTOFF now lives in Eulerian/Include/VOF_DIFFUSE.h */
#define DIFFUSE_MASS_EPS 1.0e-30
#define DIFFUSE_SOLID_GEOM_CUTOFF 1.0e-12

#ifndef VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
#define VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE 0
#endif

/*
 * Log-only audit of the conservative admissibility clip under PHASE_CHANGE
 * (roadmap A.5.1 item 4).  Define to a positive integer N to print the running
 * clip / melt budget every N timesteps; leave undefined for production, in
 * which case the build is bit-identical to before this diagnostic existed.
 * Enable with e.g.  -DVOF_DIFFUSE_CLIP_AUDIT=100  in make.def, or by
 * uncommenting the line below.
 */
/* #define VOF_DIFFUSE_CLIP_AUDIT 100 */

/*
 * Log-only audit of the SEDIMENT half of the admissibility clip (roadmap
 * B.1.4).  diffuse_bound_liquid_fraction trims C_L -> 1 - C_S in every cell,
 * but the fluid-side accounting deliberately excludes cells that are >= 5%
 * sediment, so liquid trimmed inside a resolved grain is not returned by the
 * conservative redistribution -- it is destroyed.  Define to a positive
 * integer N to print the running sediment-trim budget every N timesteps.
 * Enable with e.g.  -DVOF_DIFFUSE_SEDIMENT_CLIP_AUDIT=100  in make.def.
 * Undefined for production => bit-identical to before this diagnostic existed.
 */
/* #define VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT 100 */

/*
 * The sediment channel of the clip costs one extra MPI_Allreduce, so it is
 * only computed when something actually consumes it.  Never true in a Stage-A
 * build (VOF_IBM undefined) => the pre-B.1.4 collective sequence is preserved
 * byte for byte there.
 */
#if defined(VOF_IBM) && (defined(VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE) || \
                         defined(VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT))
#define DIFFUSE_TRACK_SEDIMENT_CLIP 1
#else
#define DIFFUSE_TRACK_SEDIMENT_CLIP 0
#endif

#ifdef VOF_DIFFUSE

#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
static int diffuse_evolving_ice = 0;
static void diffuse_step_ice(Cart3d_bag *db);
#endif
static double diffuse_liquid_mass_integral(Cart3d_bag *db);
static double diffuse_bound_liquid_fraction_split(Cart3d_bag *db,
                                                  double *sediment_delta);
static double diffuse_bound_liquid_fraction(Cart3d_bag *db);
static void diffuse_restore_liquid_mass(Cart3d_bag *db, double target_mass,
                                        const char *where);
#if defined(VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE) && defined(VOF_IBM)
static double diffuse_redistribute_sediment_clip_delta(Cart3d_bag *db,
                                                       double delta);
#endif
#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)

typedef struct {
    double X[3];
    double R;
} Diffuse_solid_geom;

static double diffuse_solid_support_extra(double cn)
{
    if (cn <= 0.0)
        return 0.0;

    const double shift = sqrt(2.0) * log(19.0) * cn;
    const double denom = 2.0 * sqrt(2.0) * cn;
    const double eps = DIFFUSE_SOLID_GEOM_CUTOFF;
    const double x = 1.0 - 2.0 * eps;
    const double atanh_x = 0.5 * log((1.0 + x) / (1.0 - x));
    const double extra = -shift + denom * atanh_x;

    return (extra > 0.0) ? extra : 0.0;
}

static Diffuse_solid_geom *
diffuse_collect_owned_particle_geometry(Particle_list *p_list,
                                        Cart3d_bag *db,
                                        double extra_range,
                                        int *n_local)
{
    int n_particles = 0;
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
    Particle *particles=Particle_collect_owned_centers(p_list,db,extra_range,&n_particles);
#else
    Particle *particles =
        Particle_collect_owned_overlaps(p_list, db, extra_range, -1.0,
                                        1, &n_particles);
#endif

    if (n_particles == 0) {
        *n_local = 0;
        return NULL;
    }

    Diffuse_solid_geom *geom =
        (Diffuse_solid_geom *)malloc(n_particles * sizeof(Diffuse_solid_geom));
    Memory_check_allocation(geom);

    for (int n = 0; n < n_particles; n++) {
        geom[n].X[0] = particles[n].X[0];
        geom[n].X[1] = particles[n].X[1];
        geom[n].X[2] = particles[n].X[2];
        geom[n].R = particles[n].R;
    }

    free(particles);

    *n_local = n_particles;
    return geom;
}

#endif

static inline double diffuse_clamp01(double value)
{
    return clampDouble(value, 0.0, 1.0);
}

static void diffuse_copy_C_to_F(Cart3d_bag *db)
{
    VolumeFraction *vof = db->vof;
    Array_copy_withghost(vof->C_L, vof->F, db->grid, db->params);
}

/*
 * Public phase-cache helper.  MCL/wetting code in VOF_IBM_wetting.c calls this
 * after rebuilding C_S or after writing into C_L; the C_G derived field and the
 * mirrored F array stay in sync.
 */
void VOF_DIFFUSE_update_phase_cache(Cart3d_bag *db)
{
    MAC_grid       *grid = db->grid;
    VolumeFraction *vof  = db->vof;

    int Is = grid->L_Is, Ie = grid->L_Ie;
    int Js = grid->L_Js, Je = grid->L_Je;
    int Ks = grid->L_Ks, Ke = grid->L_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            for (int i = Is; i < Ie; ++i) {
                vof->C_G[k][j][i] = diffuse_clamp01(1.0 - vof->C_L[k][j][i] - vof->C_S[k][j][i]);
            }
        }
    }

    diffuse_copy_C_to_F(db);
}

static double diffuse_inner_product(double ***a, double ***b, Cart3d_bag *db)
{
    MAC_grid   *grid   = db->grid;
    Parameters *params = db->params;
    double local_sum = 0.0;
    double global_sum = 0.0;

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double measure = TwodOps_cell_measure_c(grid, params, i, j, k);
                local_sum += measure * a[k][j][i] * b[k][j][i];
            }
        }
    }

    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);
    return global_sum;
}

static double diffuse_physical_cell_measure(const MAC_grid *grid,
                                            const Parameters *params,
                                            int i, int j, int k)
{
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
    /* G_*e includes the high boundary plane on the last rank. It is a
     * boundary copy, not a physical finite-volume cell. */
    if (i >= grid->NX-1 || j >= grid->NY-1 || k >= grid->NZ-1) return 0.0;
#endif
    double measure = TwodOps_cell_measure_c(grid, params, i, j, k);

    if (params->twod_mode_enabled) {
        const int nk = grid->G_Ke - grid->G_Ks;
        if (nk > 1)
            measure /= (double)nk;
    }

    return measure;
}

static double diffuse_norm(double ***a, Cart3d_bag *db)
{
    return sqrt(diffuse_inner_product(a, a, db));
}

static inline double diffuse_grad_x(double ***f, MAC_grid *grid, int i, int j, int k)
{
    return (f[k][j][i+1] - f[k][j][i-1]) * grid->i2dx_c[i];
}

static inline double diffuse_grad_y(double ***f, MAC_grid *grid, int i, int j, int k)
{
    return (f[k][j+1][i] - f[k][j-1][i]) * grid->i2dy_c[j];
}

static inline double diffuse_grad_z(double ***f, MAC_grid *grid, int i, int j, int k)
{
    return (f[k+1][j][i] - f[k-1][j][i]) * grid->i2dz_c[k];
}

static inline double diffuse_diag_entry(MAC_grid *grid,
                                        Parameters *params,
                                        int i,
                                        int j,
                                        int k,
                                        double lambda)
{
    double lap_diag = TwodOps_scalar_diag_from_face_betas(
        grid, params, i, j, k, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

    return 1.0 + lambda * lap_diag;
}

static void diffuse_apply_helmholtz(double ***x, double ***Ax, double lambda,
                                    Cart3d_bag *db)
{
    MAC_grid   *grid   = db->grid;
    Parameters *params = db->params;

    VOF_DIFFUSE_set_boundary_values(x, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double lap =
                    TwodOps_scalar_laplacian(grid, params, x, i, j, k);

                Ax[k][j][i] = x[k][j][i] - lambda * lap;
            }
        }
    }
}

static int diffuse_solve_helmholtz_pcg(double ***x, double ***b, double lambda,
                                       double *final_rel_res,
                                       Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    double ***r  = vof->ch_res;
    double ***p  = vof->ch_dir;
    double ***z  = vof->bulk_S;
    double ***Ap = vof->lap_C;

    int maxit = (params->ch_iter_max > 0) ? params->ch_iter_max : params->CG_MAXIT;
    double tol = (params->ch_tol > 0.0) ? params->ch_tol : params->CG_ETOL;

    diffuse_apply_helmholtz(x, Ap, lambda, db);
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                r[k][j][i] = b[k][j][i] - Ap[k][j][i];
                z[k][j][i] =
                    r[k][j][i] /
                    diffuse_diag_entry(grid, params, i, j, k, lambda);
                p[k][j][i] = z[k][j][i];
            }
        }
    }

    double rhs_norm = diffuse_norm(b, db);
    if (rhs_norm < 1e-30) {
        rhs_norm = 1.0;
    }

    double rz = diffuse_inner_product(r, z, db);
    double rel_res = diffuse_norm(r, db) / rhs_norm;
    int iter = 0;

    if (rel_res < tol) {
        if (final_rel_res != NULL) {
            *final_rel_res = rel_res;
        }
        VOF_DIFFUSE_set_boundary_values(x, db);
        return iter;
    }

    while (iter < maxit) {
        diffuse_apply_helmholtz(p, Ap, lambda, db);

        double pAp = diffuse_inner_product(p, Ap, db);
        if (fabs(pAp) < 1e-30) {
            break;
        }

        double alpha = rz / pAp;

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    x[k][j][i] += alpha * p[k][j][i];
                    r[k][j][i] -= alpha * Ap[k][j][i];
                }
            }
        }

        ++iter;
        rel_res = diffuse_norm(r, db) / rhs_norm;
        if (rel_res < tol) {
            break;
        }

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    z[k][j][i] =
                        r[k][j][i] /
                        diffuse_diag_entry(grid, params, i, j, k, lambda);
                }
            }
        }

        double rz_new = diffuse_inner_product(r, z, db);
        double beta = rz_new / (rz + 1e-30);

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    p[k][j][i] = z[k][j][i] + beta * p[k][j][i];
                }
            }
        }

        rz = rz_new;
    }

    if (final_rel_res != NULL) {
        *final_rel_res = rel_res;
    }
    VOF_DIFFUSE_set_boundary_values(x, db);
    return iter;
}

static void diffuse_apply_biharmonic_operator(double ***x, double ***Ax,
                                              double a_dt, double bih_coeff,
                                              Cart3d_bag *db)
{
    VolumeFraction *vof = db->vof;
    double ***lap_x = vof->ch_aux1;
    double ***bih_x = vof->psi_LG;
    MAC_grid *grid = db->grid;

    VOF_DIFFUSE_compute_laplacian(x, lap_x, db);
    VOF_DIFFUSE_compute_laplacian(lap_x, bih_x, db);

#if defined(VOF_DIFFUSE_SEDIMENT_CH_MASK) && defined(VOF_IBM)
    /*
     * Roadmap B.1.4-fix ROUND 3.  Sediment cells get an IDENTITY row, so the
     * Cahn-Hilliard solve carries C_L = 0 inside a resolved grain exactly
     * rather than transporting liquid into it and relying on the admissibility
     * clip to take it back out.
     *
     * Why this and not bookkeeping on the trimmed mass: the round-2 A/B showed
     * both post-hoc options produce artifacts (destroy it and the grain is a
     * permanent liquid sink that keeps its own shell dry and suppresses the
     * release; give it back and the shell saturates and releases spuriously).
     * The flux should not exist, so it is removed at source.
     *
     * Symmetry, i.e. is CG still valid: with x == 0 on the sediment set S, the
     * operator is  a_dt*I + bih*P L L P  on the fluid set (P = projection onto
     * fluid; L symmetric => P L L P = (LP)^T (LP), symmetric PSD) and the
     * identity on S.  Symmetric positive definite overall, so PCG applies
     * unchanged.  The caller zeroes x and b on S once; every CG iterate then
     * stays zero there (r = b - Ax = 0 - 0), so no per-iteration masking of x
     * is needed and the Laplacian stencil reads the correct zero.
     *
     * Inert without VOF_IBM (C_S == 0 => the branch never fires), so a Stage-A
     * build is unaffected.
     */
    VolumeFraction *vof_m = db->vof;
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                Ax[k][j][i] =
                    (vof_m->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF)
                        ? x[k][j][i]
                        : a_dt * x[k][j][i] + bih_coeff * bih_x[k][j][i];
            }
        }
    }
#else
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                Ax[k][j][i] = a_dt * x[k][j][i] + bih_coeff * bih_x[k][j][i];
            }
        }
    }
#endif
}

static int diffuse_solve_biharmonic_pcg(double ***x, double ***b,
                                        double a_dt, double bih_coeff,
                                        double *final_rel_res,
                                        Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    double ***r  = vof->ch_res;
    double ***p  = vof->ch_dir;
    double ***z  = vof->bulk_S;
    double ***Ap = vof->lap_C;

    int maxit = (params->ch_iter_max > 0) ? params->ch_iter_max : params->CG_MAXIT;
    double tol = (params->ch_tol > 0.0) ? params->ch_tol : params->CG_ETOL;

    diffuse_apply_biharmonic_operator(x, Ap, a_dt, bih_coeff, db);
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                r[k][j][i] = b[k][j][i] - Ap[k][j][i];
                z[k][j][i] = r[k][j][i];
                p[k][j][i] = z[k][j][i];
            }
        }
    }

    double rhs_norm = diffuse_norm(b, db);
    if (rhs_norm < 1e-30) {
        rhs_norm = 1.0;
    }

    double rz = diffuse_inner_product(r, z, db);
    double rel_res = diffuse_norm(r, db) / rhs_norm;
    int iter = 0;

    if (rel_res < tol) {
        if (final_rel_res != NULL) {
            *final_rel_res = rel_res;
        }
        VOF_DIFFUSE_set_boundary_values(x, db);
        return iter;
    }

    while (iter < maxit) {
        diffuse_apply_biharmonic_operator(p, Ap, a_dt, bih_coeff, db);

        double pAp = diffuse_inner_product(p, Ap, db);
        if (fabs(pAp) < 1e-30) {
            break;
        }

        double alpha = rz / pAp;

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    x[k][j][i] += alpha * p[k][j][i];
                    r[k][j][i] -= alpha * Ap[k][j][i];
                    z[k][j][i] = r[k][j][i];
                }
            }
        }

        ++iter;
        rel_res = diffuse_norm(r, db) / rhs_norm;
        if (rel_res < tol) {
            break;
        }

        double rz_new = diffuse_inner_product(r, z, db);
        double beta = rz_new / (rz + 1e-30);

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    p[k][j][i] = z[k][j][i] + beta * p[k][j][i];
                }
            }
        }

        rz = rz_new;
    }

    if (final_rel_res != NULL) {
        *final_rel_res = rel_res;
    }
    VOF_DIFFUSE_set_boundary_values(x, db);
    return iter;
}

static inline double diffuse_weno5_left(double v0, double v1, double v2,
                                        double v3, double v4)
{
    const double eps = 1e-6;

    double q0 = (2.0 * v0 - 7.0 * v1 + 11.0 * v2) / 6.0;
    double q1 = (-v1 + 5.0 * v2 + 2.0 * v3) / 6.0;
    double q2 = (2.0 * v2 + 5.0 * v3 - v4) / 6.0;

    double b0 = 13.0 / 12.0 * pow(v0 - 2.0 * v1 + v2, 2.0) +
                0.25 * pow(v0 - 4.0 * v1 + 3.0 * v2, 2.0);
    double b1 = 13.0 / 12.0 * pow(v1 - 2.0 * v2 + v3, 2.0) +
                0.25 * pow(v1 - v3, 2.0);
    double b2 = 13.0 / 12.0 * pow(v2 - 2.0 * v3 + v4, 2.0) +
                0.25 * pow(3.0 * v2 - 4.0 * v3 + v4, 2.0);

    double a0 = 0.1 / ((eps + b0) * (eps + b0));
    double a1 = 0.6 / ((eps + b1) * (eps + b1));
    double a2 = 0.3 / ((eps + b2) * (eps + b2));
    double asum = a0 + a1 + a2;

    return (a0 * q0 + a1 * q1 + a2 * q2) / asum;
}

static inline double diffuse_face_state_x(double ***C, int i, int j, int k,
                                          double vel, int weno_order)
{
    if (weno_order == 5) {
        if (vel >= 0.0) {
            return diffuse_weno5_left(C[k][j][i-3], C[k][j][i-2], C[k][j][i-1],
                                      C[k][j][i], C[k][j][i+1]);
        }
        return diffuse_weno5_left(C[k][j][i+2], C[k][j][i+1], C[k][j][i],
                                  C[k][j][i-1], C[k][j][i-2]);
    }

    return (vel >= 0.0) ? C[k][j][i-1] : C[k][j][i];
}

static inline double diffuse_face_state_y(double ***C, int i, int j, int k,
                                          double vel, int weno_order)
{
    if (weno_order == 5) {
        if (vel >= 0.0) {
            return diffuse_weno5_left(C[k][j-3][i], C[k][j-2][i], C[k][j-1][i],
                                      C[k][j][i], C[k][j+1][i]);
        }
        return diffuse_weno5_left(C[k][j+2][i], C[k][j+1][i], C[k][j][i],
                                  C[k][j-1][i], C[k][j-2][i]);
    }

    return (vel >= 0.0) ? C[k][j-1][i] : C[k][j][i];
}

static inline double diffuse_face_state_z(double ***C, int i, int j, int k,
                                          double vel, int weno_order)
{
    if (weno_order == 5) {
        if (vel >= 0.0) {
            return diffuse_weno5_left(C[k-3][j][i], C[k-2][j][i], C[k-1][j][i],
                                      C[k][j][i], C[k+1][j][i]);
        }
        return diffuse_weno5_left(C[k+2][j][i], C[k+1][j][i], C[k][j][i],
                                  C[k-1][j][i], C[k-2][j][i]);
    }

    return (vel >= 0.0) ? C[k-1][j][i] : C[k][j][i];
}

static bool diffuse_find_nearest_particle(double x, double y, double z,
                                          Particle_list *plist, double *best_d2,
                                          Particle **best_particle)
{
    bool found = false;
    Particle *p = plist->start;

    while (p != NULL) {
        double dx = x - p->X[0];
        double dy = y - p->X[1];
        double dz = z - p->X[2];
        double d2 = dx * dx + dy * dy + dz * dz;

        if (d2 < *best_d2) {
            *best_d2 = d2;
            *best_particle = p;
            found = true;
        }
        p = p->next;
    }

    return found;
}

void VOF_DIFFUSE_init(Cart3d_bag *db)
{
    MAC_grid       *grid = db->grid;
    VolumeFraction *vof  = db->vof;

#ifdef VOF_IBM
    VOF_DIFFUSE_compute_C_S(db);
#endif

    for (int k = grid->L_Ks; k < grid->L_Ke; ++k) {
        for (int j = grid->L_Js; j < grid->L_Je; ++j) {
            for (int i = grid->L_Is; i < grid->L_Ie; ++i) {
                vof->C_L[k][j][i] = diffuse_clamp01(vof->F[k][j][i]);
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);

    #if VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
    double initial_mass = diffuse_liquid_mass_integral(db);
    #endif

    diffuse_bound_liquid_fraction(db);
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);

    #if VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
    diffuse_restore_liquid_mass(db, initial_mass, "initial C_L+CS projection");
    #endif

#ifdef VOF_IBM
    VOF_DIFFUSE_apply_contact_angle(db);

    #if VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
    diffuse_restore_liquid_mass(db, initial_mass, "initial MCL projection");
    #endif

#endif

    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);
    VOF_DIFFUSE_set_boundary_values(vof->C_G, db);
}


void VOF_DIFFUSE_set_boundary_values(double ***f, Cart3d_bag *db)
{
    MAC_grid   *grid   = db->grid;
    Parameters *params = db->params;

    int NX = grid->NX;
    int NY = grid->NY;
    int NZ = grid->NZ;

    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    int IsL = grid->L_Is;
    int JsL = grid->L_Js;
    int KsL = grid->L_Ks;
    int IeL = grid->L_Ie;
    int JeL = grid->L_Je;
    int KeL = grid->L_Ke;

    const int NG = params->ghost_nodes;

#ifndef XPERIODIC
    if (Is == 0) {
        int i0 = 0;
        for (int k = KsL; k < KeL; ++k)
        for (int j = JsL; j < JeL; ++j)
        for (int g = 1; g <= NG; ++g) {
#ifdef AXISYM_RZ
            f[k][j][i0 - g] = f[k][j][i0];
#elif defined LEFT_WALL_VELOCITY_FREESLIP
            f[k][j][i0 - g] = f[k][j][i0 + g - 1];
#else
            f[k][j][i0 - g] = f[k][j][i0];
#endif
        }
    }

    if (Ie == NX) {
        int i1 = NX - 1;
        for (int k = KsL; k < KeL; ++k)
        for (int j = JsL; j < JeL; ++j)
        for (int g = 0; g < NG; ++g) {
#ifdef RIGHT_WALL_VELOCITY_FREESLIP
            f[k][j][i1 + g] = f[k][j][i1 - 1 - g];
#else
            f[k][j][i1 + g] = f[k][j][i1 - 1];
#endif
        }
    }
#endif

#ifndef YPERIODIC
    if (Js == 0) {
        int j0 = 0;
        for (int k = KsL; k < KeL; ++k)
        for (int i = IsL; i < IeL; ++i)
        for (int g = 1; g <= NG; ++g) {
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
            f[k][j0 - g][i] = f[k][j0 + g - 1][i];
#else
            f[k][j0 - g][i] = f[k][j0][i];
#endif
        }
    }

    if (Je == NY) {
        int j1 = NY - 1;
        for (int k = KsL; k < KeL; ++k)
        for (int i = IsL; i < IeL; ++i)
        for (int g = 0; g < NG; ++g) {
#ifdef TOP_WALL_VELOCITY_FREESLIP
            f[k][j1 + g][i] = f[k][j1 - 1 - g][i];
#else
            f[k][j1 + g][i] = f[k][j1 - 1][i];
#endif
        }
    }
#endif

#ifndef ZPERIODIC
    if (Ks == 0) {
        int k0 = 0;
        for (int j = JsL; j < JeL; ++j)
        for (int i = IsL; i < IeL; ++i)
        for (int g = 1; g <= NG; ++g) {
#ifdef BACK_WALL_VELOCITY_FREESLIP
            f[k0 - g][j][i] = f[k0 + g - 1][j][i];
#else
            f[k0 - g][j][i] = f[k0][j][i];
#endif
        }
    }

    if (Ke == NZ) {
        int k1 = NZ - 1;
        for (int j = JsL; j < JeL; ++j)
        for (int i = IsL; i < IeL; ++i)
        for (int g = 0; g < NG; ++g) {
#ifdef FRONT_WALL_VELOCITY_FREESLIP
            f[k1 + g][j][i] = f[k1 - 1 - g][j][i];
#else
            f[k1 + g][j][i] = f[k1 - 1][j][i];
#endif
        }
    }
#endif

    Communication_update_ghost_nodes_flow_variable(f, VOLUME_FRACTION, NG, db);
}

void VOF_DIFFUSE_compute_C_S(Cart3d_bag *db)
{
    MAC_grid       *grid = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof  = db->vof;

    const double cn = params->Cn;
    const double shift = sqrt(2.0) * log(19.0) * cn;
    const double denom = 2.0 * sqrt(2.0) * cn;

    Array_set_withghost(vof->C_S, 0.0, grid, params);

    if (db->lag == NULL) {
        VOF_DIFFUSE_update_phase_cache(db);
        return;
    }

#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)

    const double support_extra = diffuse_solid_support_extra(cn);

    /*
     * C_S is an Eulerian solid field, so reconstruct it from rank-owned
     * particle geometry, not from whatever foreign copies neighbor exchange
     * happened to leave locally.  The exchange is targeted to ranks whose
     * subdomain intersects the diffuse solid support.
     */
    Particle_list *lists[2] = {
        db->lag->p_mobile_list,
        db->lag->p_fixed_list
    };

    for (int list_id = 0; list_id < 2; ++list_id) {
        int n_geom = 0;
        Diffuse_solid_geom *geom =
            diffuse_collect_owned_particle_geometry(lists[list_id],
                                                    db,
                                                    support_extra,
                                                    &n_geom);

        for (int n = 0; n < n_geom; ++n) {
            const double *X = geom[n].X;
            const double R  = geom[n].R;
            const double support = R + support_extra;

            for (int k = grid->L_Ks; k < grid->L_Ke; ++k) {
                double dz = grid->zc[k] - X[2];
#if defined(VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT) && defined(ZPERIODIC)
                        dz -= params->Lz * nearbyint(dz/params->Lz);
#endif

                for (int j = grid->L_Js; j < grid->L_Je; ++j) {
                    double dy = grid->yc[j] - X[1];
#if defined(VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT) && defined(YPERIODIC)
                        dy -= params->Ly * nearbyint(dy/params->Ly);
#endif

                    for (int i = grid->L_Is; i < grid->L_Ie; ++i) {
                        double dx = grid->xc[i] - X[0];
#if defined(VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT) && defined(XPERIODIC)
                        dx -= params->Lx * nearbyint(dx/params->Lx);
#endif

                        double dist = sqrt(dx * dx + dy * dy + dz * dz);

#if defined(TWOD_CARTESIAN) && defined(LAG_PARTICLE_RESOLVED)
                        /*
                         * Planar 2D IBM represents particles as cylinders
                         * extruded through the dummy z slab.
                         */
                        dist = sqrt(dx * dx + dy * dy);
#elif defined(AXISYM_RZ) && defined(LAG_PARTICLE_RESOLVED)
                        /*
                         * Axisymmetric IBM stores the meridional sphere in
                         * code coordinates x-y.
                         */
                        dist = sqrt(dx * dx + dy * dy);
#endif
                        if (dist > support)
                            continue;

                        const double arg =
                            (dist - (R - shift)) / (denom + 1.0e-30);
                        const double cs = 0.5 - 0.5 * tanh(arg);

                        if (cs > vof->C_S[k][j][i])
                            vof->C_S[k][j][i] = cs;
                    }
                }
            }
        }

        free(geom);
    }

#else

    Particle_list *lists[2] = {
        db->lag->p_mobile_list,
        db->lag->p_fixed_list
    };

    for (int list_id = 0; list_id < 2; ++list_id) {
        Particle *p = lists[list_id]->start;

        while (p != NULL) {
            for (int k = grid->L_Ks; k < grid->L_Ke; ++k) {
                double dz = grid->zc[k] - p->X[2];

                for (int j = grid->L_Js; j < grid->L_Je; ++j) {
                    double dy = grid->yc[j] - p->X[1];

                    for (int i = grid->L_Is; i < grid->L_Ie; ++i) {
                        double dx = grid->xc[i] - p->X[0];
                        double dist = sqrt(dx * dx + dy * dy + dz * dz);

                        double arg =
                            (dist - (p->R - shift)) / (denom + 1.0e-30);
                        double cs = 0.5 - 0.5 * tanh(arg);

                        if (cs > vof->C_S[k][j][i])
                            vof->C_S[k][j][i] = cs;
                    }
                }
            }

            p = p->next;
        }
    }

#endif

    VOF_DIFFUSE_set_boundary_values(vof->C_S, db);
    diffuse_bound_liquid_fraction(db);
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);
    VOF_DIFFUSE_set_boundary_values(vof->C_G, db);
}


#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
/* Symmetric face-neighbour sum. Physical walls have no transfer; periodic
 * neighbours are read from the freshly exchanged halo. Uniform Cartesian
 * cells have equal volume, so the paired transfers conserve physical volume. */
static double diffuse_neighbour_sum(double ***a,double ***b,
                                    MAC_grid *g,int i,int j,int k)
{
    double sum=0.0;
    for(int axis=0;axis<3;++axis)
    for(int sign=-1;sign<=1;sign+=2) {
        int ii=i+(axis==0?sign:0),jj=j+(axis==1?sign:0),kk=k+(axis==2?sign:0);
#ifndef XPERIODIC
        if(ii<0 || ii>=g->NX-1) continue;
#endif
#ifndef YPERIODIC
        if(jj<0 || jj>=g->NY-1) continue;
#endif
#ifndef ZPERIODIC
        if(kk<0 || kk>=g->NZ-1) continue;
#endif
        sum+=a[kk][jj][ii]*(b ? b[kk][jj][ii] : 1.0);
    }
    return sum;
}

void VOF_DIFFUSE_move_ice_mask(Cart3d_bag *db)
{
    MAC_grid *g=db->grid;
    Parameters *p=db->params;
    VolumeFraction *v=db->vof;
    const int ie=min(g->G_Ie,g->NX-1),je=min(g->G_Je,g->NY-1),ke=min(g->G_Ke,g->NZ-1);
    /* Scratch fields are rebuilt by the following CH stage. No history is
     * borrowed, and all donor states remain frozen during each transfer. */
    double ***ice=v->ch_aux2,***lost=v->ch_aux1;
    double ***room=v->psi,***ratio=v->bulk_S,***accept=v->psi_LG;
    for(int k=g->G_Ks;k<ke;++k)
    for(int j=g->G_Js;j<je;++j)
    for(int i=g->G_Is;i<ie;++i)
        ice[k][j][i]=v->C_S[k][j][i]>=DIFFUSE_SOLID_MASS_CUTOFF ? 0.0 :
            fmax(0.0,1.0-v->C_S[k][j][i]-v->C_L[k][j][i]);
    VOF_DIFFUSE_compute_C_S(db);
    double local=0.0,remaining=0.0;
    for(int k=g->G_Ks;k<ke;++k)
    for(int j=g->G_Js;j<je;++j)
    for(int i=g->G_Is;i<ie;++i) {
        const double cap=v->C_S[k][j][i]>=DIFFUSE_SOLID_MASS_CUTOFF ? 0.0 : 1.0-v->C_S[k][j][i];
        const double ni=fmin(cap,ice[k][j][i]);
        lost[k][j][i]=ice[k][j][i]-ni;
        ice[k][j][i]=ni;
        local+=lost[k][j][i]*diffuse_physical_cell_measure(g,p,i,j,k);
    }
    MPI_Allreduce(&local,&remaining,1,MPI_DOUBLE,MPI_SUM,PCW);
    static double moved=0.0,unplaced=0.0;
    moved+=remaining;
    /* Geometric displacement is not phase change: move real ice to adjacent
     * available pore space, never silently delete it or inject meltwater. */
    for(int pass=0;pass<4 && remaining>1.e-14;++pass) {
        for(int k=g->G_Ks;k<ke;++k)
        for(int j=g->G_Js;j<je;++j)
        for(int i=g->G_Is;i<ie;++i) {
            const double cap=v->C_S[k][j][i]>=DIFFUSE_SOLID_MASS_CUTOFF ? 0.0 : 1.0-v->C_S[k][j][i];
            room[k][j][i]=fmax(0.0,cap-ice[k][j][i]);
        }
        VOF_DIFFUSE_set_boundary_values(room,db);
        for(int k=g->G_Ks;k<ke;++k)
        for(int j=g->G_Js;j<je;++j)
        for(int i=g->G_Is;i<ie;++i) {
            const double available=lost[k][j][i]>0.0 ? diffuse_neighbour_sum(room,NULL,g,i,j,k) : 0.0;
            ratio[k][j][i]=available>0.0 ? lost[k][j][i]/available : 0.0;
        }
        VOF_DIFFUSE_set_boundary_values(ratio,db);
        for(int k=g->G_Ks;k<ke;++k)
        for(int j=g->G_Js;j<je;++j)
        for(int i=g->G_Is;i<ie;++i) {
            const double rate=room[k][j][i]>0.0 ? diffuse_neighbour_sum(ratio,NULL,g,i,j,k) : 0.0;
            accept[k][j][i]=rate>1.0 ? 1.0/rate : 1.0;
        }
        VOF_DIFFUSE_set_boundary_values(accept,db);
        local=0.0;
        for(int k=g->G_Ks;k<ke;++k)
        for(int j=g->G_Js;j<je;++j)
        for(int i=g->G_Is;i<ie;++i) {
            if(room[k][j][i]>0.0)
                ice[k][j][i]+=room[k][j][i]*accept[k][j][i]*diffuse_neighbour_sum(ratio,NULL,g,i,j,k);
            if(ratio[k][j][i]>0.0)
                lost[k][j][i]=fmax(0.0,lost[k][j][i]-ratio[k][j][i]*diffuse_neighbour_sum(room,accept,g,i,j,k));
            local+=lost[k][j][i]*diffuse_physical_cell_measure(g,p,i,j,k);
        }
        MPI_Allreduce(&local,&remaining,1,MPI_DOUBLE,MPI_SUM,PCW);
    }
    unplaced+=remaining;
    if(remaining>1.e-12) {
        char name[80],msg[180];
        snprintf(msg,sizeof(msg),"ICE_MASK remaining=%.17g time=%.17g stage=%d\n",remaining,p->time,p->which_stage);
        Display_progress(p,msg);
        snprintf(name,sizeof(name),"remap_failure_rank%04d.csv",p->rank);
        FILE *f=fopen(name,"w");
        if(f) {
            fprintf(f,"i,j,k,cs,ice,lost,room6,room26\n");
            for(int k=g->G_Ks;k<ke;++k)
            for(int j=g->G_Js;j<je;++j)
            for(int i=g->G_Is;i<ie;++i) if(lost[k][j][i]>1.e-14) {
                double near=0.0;
                for(int dk=-1;dk<=1;++dk)
                for(int dj=-1;dj<=1;++dj)
                for(int di=-1;di<=1;++di)
                    if(j+dj>=0 && j+dj<g->NY-1) near+=room[k+dk][j+dj][i+di];
                fprintf(f,"%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",i,j,k,v->C_S[k][j][i],ice[k][j][i],lost[k][j][i],diffuse_neighbour_sum(room,NULL,g,i,j,k),near);
            }
            fclose(f);
        }
        MPI_Barrier(PCW);
        Display_progress(p,"ICE_MASK ERROR: adjacent pore space cannot hold displaced ice; reduce motion timestep or review release geometry.\n");
        MPI_Abort(PCW,98);
    }
    for(int k=g->G_Ks;k<ke;++k)
    for(int j=g->G_Js;j<je;++j)
    for(int i=g->G_Is;i<ie;++i) {
        const double cap=v->C_S[k][j][i]>=DIFFUSE_SOLID_MASS_CUTOFF ? 0.0 : 1.0-v->C_S[k][j][i];
        v->C_L[k][j][i]=fmax(0.0,cap-ice[k][j][i]);
    }
    if(p->which_stage==2 && p->ntime%100==0) {
        char msg[200];
        snprintf(msg,sizeof(msg),"ICE_MASK ntime=%d redistributed_ice=%.12e unplaced=%.12e\n",p->ntime,moved,unplaced);
        Display_progress(p,msg);
    }
    VOF_DIFFUSE_set_boundary_values(v->C_L,db);
    VOF_DIFFUSE_update_phase_cache(db);
    VOF_DIFFUSE_set_boundary_values(v->C_G,db);
}
#endif


#ifndef VOF_IBM
/*
 * Without VOF_IBM the contact-angle / chemical-potential extension is a no-op:
 * keep the same name so callers in Temporal_int.c and VOF_DIFFUSE_step compile
 * regardless of the IBM flag.
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
#endif

#if defined(VOF_DIFFUSE_SEDIMENT_CH_MASK) && defined(VOF_IBM)
#if defined(AXISYM_RZ)
#error "VOF_DIFFUSE_SEDIMENT_CH_MASK has no axisymmetric face correction yet"
#endif
/*
 * Roadmap B.1.4-fix ROUND 4 -- ZERO-FLUX, not zero-value, at the grain surface.
 *
 * Round 3 pinned C_L = 0 inside the grain with identity rows.  That killed the
 * leak exactly (audit net went to 0.000000e+00) but imposed the wrong physics:
 * a Dirichlet zero says the rock forces ICE against itself, and the measured
 * consequence was a depletion layer ten cells deep -- C_L 0.033 at the surface
 * recovering to 0.536 only by 8-16 cells out, against a flat ~0.68 in the
 * control.  The release criterion collapsed from 0.668 to 0.082.
 *
 * The physically correct condition is that the rock is INDIFFERENT between
 * water and ice: no phase crosses into it, but neither is preferred.  That is
 * homogeneous Neumann on C_L at the sediment boundary -- equivalently the 90
 * degree contact angle the deck already specifies and B.1.2's "no capillary /
 * contact-angle path active".
 *
 * TwodOps_scalar_laplacian is shared with the pressure solver and must not be
 * touched, so instead of rewriting the stencil we CANCEL the face terms that
 * reach into sediment.  Each face of the 5/7-point stencil enters as
 * (phi[nb] - phi[c]) * idx_c * idx_u, so subtracting exactly that term makes
 * the face contribute zero -- an exact no-flux face, not an approximation.
 *
 * Symmetry is preserved: the fluid-sediment coupling is dropped from the fluid
 * row, and the sediment row is already an identity row, so it never coupled
 * back.  The remaining fluid-fluid operator is the original symmetric one, and
 * the biharmonic L(L(.)) composes two such operators, so PCG stays valid.
 */
static double diffuse_sediment_noflux_correction(const MAC_grid   *grid,
                                                 const Parameters *params,
                                                 VolumeFraction   *vof,
                                                 double ***phi,
                                                 int i, int j, int k)
{
    const double cut = DIFFUSE_SOLID_MASS_CUTOFF;
    double corr = 0.0;

    if (vof->C_S[k][j][i + 1] >= cut)
        corr -= (phi[k][j][i + 1] - phi[k][j][i]) * grid->idx_c[i] * grid->idx_u[i];
    if (vof->C_S[k][j][i - 1] >= cut)
        corr += (phi[k][j][i] - phi[k][j][i - 1]) *
                ((i != 0) ? grid->idx_c[i - 1] : grid->idx_c[i]) * grid->idx_u[i];

    if (vof->C_S[k][j + 1][i] >= cut)
        corr -= (phi[k][j + 1][i] - phi[k][j][i]) * grid->idy_c[j] * grid->idy_v[j];
    if (vof->C_S[k][j - 1][i] >= cut)
        corr += (phi[k][j][i] - phi[k][j - 1][i]) *
                ((j != 0) ? grid->idy_c[j - 1] : grid->idy_c[j]) * grid->idy_v[j];

    if (!params->twod_mode_enabled) {
        if (vof->C_S[k + 1][j][i] >= cut)
            corr -= (phi[k + 1][j][i] - phi[k][j][i]) * grid->idz_c[k] * grid->idz_w[k];
        if (vof->C_S[k - 1][j][i] >= cut)
            corr += (phi[k][j][i] - phi[k - 1][j][i]) *
                    ((k != 0) ? grid->idz_c[k - 1] : grid->idz_c[k]) * grid->idz_w[k];
    }

    return corr;
}
#endif

void VOF_DIFFUSE_compute_laplacian(double ***in, double ***lap, Cart3d_bag *db)
{
    MAC_grid   *grid = db->grid;
    Parameters *params = db->params;

    VOF_DIFFUSE_set_boundary_values(in, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                /*
                 * Route the scalar Laplacian through the shared 2D helper so
                 * the Cahn-Hilliard operator follows the same geometry rules
                 * as the pressure solver:
                 *   - full 3D when TWOD_MODE is off
                 *   - planar x-y in TWOD_CARTESIAN
                 *   - axisymmetric r-z in AXISYM_RZ
                 */
                lap[k][j][i] = TwodOps_scalar_laplacian(
                    grid, params, in, i, j, k);

#if defined(VOF_DIFFUSE_SEDIMENT_CH_MASK) && defined(VOF_IBM)
                /* The grain is a hole in the CH domain with no-flux walls. */
                if (db->vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF)
                    lap[k][j][i] = 0.0;
                else
                    lap[k][j][i] += diffuse_sediment_noflux_correction(
                        grid, params, db->vof, in, i, j, k);
#endif
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(lap, db);
}

void VOF_DIFFUSE_compute_bulk_S(Cart3d_bag *db)
{
    MAC_grid       *grid = db->grid;
    VolumeFraction *vof  = db->vof;

    

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double c = vof->C_L[k][j][i];
                double cs = vof->C_S[k][j][i];
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
                if (diffuse_evolving_ice) cs = 0.0;
#endif
                vof->bulk_S[k][j][i] =
                    0.5 * c - 1.5 * c * c + c * c * c -
                    c * cs * (1.0 - c - cs);
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->bulk_S, db);
}

void VOF_DIFFUSE_compute_psi(Cart3d_bag *db)
{
    Parameters     *params = db->params;
    MAC_grid       *grid   = db->grid;
    VolumeFraction *vof    = db->vof;

    VOF_DIFFUSE_compute_laplacian(vof->C_L, vof->lap_C, db);
    VOF_DIFFUSE_compute_bulk_S(db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                vof->psi[k][j][i] = vof->bulk_S[k][j][i] - params->Cn * params->Cn * vof->lap_C[k][j][i];
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->psi, db);
}

void VOF_DIFFUSE_compute_psi_LG(Cart3d_bag *db)
{
    Parameters     *params = db->params;
    MAC_grid       *grid   = db->grid;
    VolumeFraction *vof    = db->vof;

    VOF_DIFFUSE_compute_laplacian(vof->C_L, vof->lap_C, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double c = vof->C_L[k][j][i];
                vof->psi_LG[k][j][i] =
                    c * c * c - 1.5 * c * c + 0.5 * c -
                    params->Cn * params->Cn * vof->lap_C[k][j][i];
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->psi_LG, db);


}

void VOF_DIFFUSE_advect_WENO5(Cart3d_bag *db, double ***rhs_out)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;
    Velocity       *u      = db->u;
    Velocity       *v      = db->v;
    Velocity       *w      = db->w;

    double ***flux_x = vof->flux_x;
    double ***flux_y = vof->flux_y;
    double ***flux_z = vof->flux_z;
    const int collapsed_z = TwodOps_collapsed_component_is_inactive(params);

    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);

    for (int k = grid->L_Ks; k < grid->L_Ke; ++k)
    for (int j = grid->L_Js; j < grid->L_Je; ++j)
    for (int i = grid->L_Is; i < grid->L_Ie; ++i) {
        flux_x[k][j][i] = 0.0;
        flux_y[k][j][i] = 0.0;
        flux_z[k][j][i] = 0.0;
    }

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double vel = u->data[k][j][i];
                double state = diffuse_face_state_x(vof->C_L, i, j, k, vel, params->weno_order);
                flux_x[k][j][i] = vel * state;
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
                if (diffuse_evolving_ice &&
                    (vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF ||
                     vof->C_S[k][j][i-1] >= DIFFUSE_SOLID_MASS_CUTOFF))
                    flux_x[k][j][i] = 0.0;
#endif

                vel = v->data[k][j][i];
                state = diffuse_face_state_y(vof->C_L, i, j, k, vel, params->weno_order);
                flux_y[k][j][i] = vel * state;
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
                if (diffuse_evolving_ice &&
                    (vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF ||
                     vof->C_S[k][j-1][i] >= DIFFUSE_SOLID_MASS_CUTOFF))
                    flux_y[k][j][i] = 0.0;
#endif

                if (!collapsed_z) {
                    vel = w->data[k][j][i];
                    state = diffuse_face_state_z(vof->C_L, i, j, k, vel, params->weno_order);
                    flux_z[k][j][i] = vel * state;
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
                if (diffuse_evolving_ice &&
                    (vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF ||
                     vof->C_S[k-1][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF))
                    flux_z[k][j][i] = 0.0;
#endif
                }
            }
        }
    }

    Communication_update_ghost_nodes_flow_variable(flux_x, FLUX_X, params->ghost_nodes, db);
    Communication_update_ghost_nodes_flow_variable(flux_y, FLUX_Y, params->ghost_nodes, db);
    Communication_update_ghost_nodes_flow_variable(flux_z, FLUX_Z, params->ghost_nodes, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                /*
                 * In 2D mode this uses the physical in-plane divergence only.
                 * The z slab still exists for storage and halo exchange, but
                 * it must not contribute to the transported phase fraction.
                 */
                double div = TwodOps_scalar_flux_divergence(
                    grid, params, flux_x, flux_y, flux_z, i, j, k);

                rhs_out[k][j][i] = -div;
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(rhs_out, db);
}

void VOF_DIFFUSE_explicit_diffusion(Cart3d_bag *db, double ***rhs_out)
{
    MAC_grid       *grid = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof  = db->vof;

    VOF_DIFFUSE_compute_bulk_S(db);
    VOF_DIFFUSE_compute_laplacian(vof->bulk_S, vof->ch_aux1, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                rhs_out[k][j][i] = vof->ch_aux1[k][j][i] / params->Pe_CH;
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(rhs_out, db);
}

void VOF_DIFFUSE_solve_implicit_biharmonic(Cart3d_bag *db, double ***rhs_explicit)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;
    const double BET[] = { BETA };
    char statement[200];

    /*
     * Conc.c-style RK3/CN form for the diffuse CH step:
     *
     *   [ a_dt I + (Cn^2/Pe) L^2 ] C^k
     * = a_dt C^{k-1}
     *   + (GAMMA/BETA)_k R^k + (ZETA/BETA)_k R^{k-1}
     *   - (Cn^2/Pe) L^2 C^{k-1},
     *
     * where a_dt = 1/(BETA_k dt) and R is the explicit CH RHS
     * (-div(uC) + Pe^{-1} L S(C)).
     *
     * This is algebraically equivalent to the previous BETA2/GAMMA/ZETA form,
     * but it now matches the Concentration solver notation and assembly style
     * exactly.
     */
    double a_dt = 1.0 / (params->dt * BET[params->which_stage]);
    double bih_coeff = params->Cn * params->Cn / params->Pe_CH;
    double tol = (params->ch_tol > 0.0) ? params->ch_tol : params->CG_ETOL;
    double final_res = 0.0;
    int iters;

    Array_copy_withghost(vof->C_L, vof->psi, grid, params);

    VOF_DIFFUSE_compute_laplacian(vof->psi, vof->ch_aux1, db);
    VOF_DIFFUSE_compute_laplacian(vof->ch_aux1, vof->psi_LG, db);
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                rhs_explicit[k][j][i] -= bih_coeff * vof->psi_LG[k][j][i];
            }
        }
    }

#if defined(VOF_DIFFUSE_SEDIMENT_CH_MASK) && defined(VOF_IBM)
    /*
     * Roadmap B.1.4-fix ROUND 3 -- enforce the Dirichlet state the identity
     * rows in diffuse_apply_biharmonic_operator expect.  Done once, before the
     * CG: with x and b both zero on the sediment set the iterates stay zero
     * there for the whole solve.  The halo refresh matters because the
     * Laplacian stencil of a neighbouring FLUID cell reads these values.
     */
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
        if (vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF) {
            vof->C_L[k][j][i]     = 0.0;
            rhs_explicit[k][j][i] = 0.0;
        }
    }
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_set_boundary_values(rhs_explicit, db);
#endif

    iters = diffuse_solve_biharmonic_pcg(vof->C_L, rhs_explicit,
                                         a_dt, bih_coeff, &final_res, db);
    

    sprintf(statement,
            "VOF diffuse biharmonic %s to %g after %d iterations\n",
            (final_res <= tol) ? "converged" : "stopped",
            final_res,
            iters);
    Display_progress(params, statement);
}

static double diffuse_liquid_mass_integral(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    double local = 0.0;

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
        if (vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF)
            continue;
        local += vof->C_L[k][j][i] *
                 diffuse_physical_cell_measure(grid, params, i, j, k);
    }

    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, PCW);
    return global;
}

/*
 * Admissibility clip C_L -> [0, 1 - C_S], with the trimmed mass accounted in
 * TWO separate channels (roadmap B.1.4):
 *
 *   return value      cells with C_S <  DIFFUSE_SOLID_MASS_CUTOFF -- the
 *                     fluid-side clip.  This is the Stage-A quantity; the
 *                     caller returns it to the interfacial band, which is what
 *                     makes the clip conservative (A.3 gate, ~5% front-speed
 *                     bias without it).
 *
 *   *sediment_delta   cells with C_S >= DIFFUSE_SOLID_MASS_CUTOFF -- liquid
 *                     trimmed INSIDE a resolved sediment grain.  Historically
 *                     this was silently discarded, which destroys real liquid
 *                     mass and biases the ice-volume budget every ECCO melt
 *                     number rests on (measured 462 cell-units against a
 *                     171-cell-unit grain at t=30 in the 2-D release deck).
 *                     Pass NULL to ignore it.
 *
 * In a Stage-A build C_S == 0 everywhere, so *sediment_delta is identically
 * zero and this function is bit-identical to its pre-B.1.4 form.
 */
static double diffuse_bound_liquid_fraction_split(Cart3d_bag *db,
                                                  double *sediment_delta)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    double local_delta = 0.0;
    double local_sed   = 0.0;

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
        double old_cl = vof->C_L[k][j][i];
        double cs     = vof->C_S[k][j][i];
        double max_cl = 1.0 - cs;

        if (max_cl < 0.0) max_cl = 0.0;
        if (max_cl > 1.0) max_cl = 1.0;

        double new_cl = old_cl;
        if (new_cl < 0.0) new_cl = 0.0;
        if (new_cl > max_cl) new_cl = max_cl;

        vof->C_L[k][j][i] = new_cl;
        vof->C_G[k][j][i] = diffuse_clamp01(1.0 - new_cl - cs);

        double measure = diffuse_physical_cell_measure(grid, params, i, j, k);

        if (cs < DIFFUSE_SOLID_MASS_CUTOFF) {
            local_delta += (new_cl - old_cl) * measure;
        } else {
            local_sed   += (new_cl - old_cl) * measure;
        }
    }

    /*
     * The fluid-side reduction is left EXACTLY as it was (one MPI_DOUBLE,
     * MPI_SUM) rather than folded into a 2-element reduce, so that a build
     * which never asks for the sediment channel issues a byte-for-byte
     * identical sequence of collectives.  MPI reductions are deterministic per
     * (count, datatype, op, communicator) but not across counts, and the
     * Stage-A no-impact gate is required to reproduce lambda to every printed
     * digit -- not worth spending that guarantee to save one collective on the
     * ECCO path.
     */
    double global_delta = 0.0;
    MPI_Allreduce(&local_delta, &global_delta, 1, MPI_DOUBLE, MPI_SUM, PCW);

    if (sediment_delta != NULL) {
        double global_sed = 0.0;
        MPI_Allreduce(&local_sed, &global_sed, 1, MPI_DOUBLE, MPI_SUM, PCW);
        *sediment_delta = global_sed;
    }

    return global_delta;
}

static double diffuse_bound_liquid_fraction(Cart3d_bag *db)
{
    return diffuse_bound_liquid_fraction_split(db, NULL);
}


static int diffuse_mass_redist_eligible(double cl, double cs, int pass)
{
    if (cs >= DIFFUSE_SOLID_MASS_CUTOFF) return 0;

    /* First pass: true diffuse interface only. */
    if (pass == 0)
        return (cl > 0.005 && cl < 0.995);

    /* Second pass: slightly relaxed, but still not bulk gas or bulk liquid. */
    return (cl > 1.0e-4 && cl < 1.0 - 1.0e-4);
}

static double diffuse_redistribute_liquid_mass_delta(Cart3d_bag *db, double delta)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    if (fabs(delta) < 1.0e-14)
        return 0.0;

    for (int pass = 0; pass < 2 && fabs(delta) > 1.0e-14; ++pass) {

        double local_weight = 0.0;

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
        for (int j = grid->G_Js; j < grid->G_Je; ++j)
        for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

            const double cs = vof->C_S[k][j][i];
            const double cl = vof->C_L[k][j][i];

            if (!diffuse_mass_redist_eligible(cl, cs, pass))
                continue;

            const double cap = (delta > 0.0) ? (1.0 - cs - cl) : cl;
            if (cap <= 0.0)
                continue;

            const double measure = diffuse_physical_cell_measure(grid, params, i, j, k);

            /*
             * Weight concentrated at C_L = 0.5.
             * This moves the interface instead of polluting the bulk gas.
             */
            const double w = cl * (1.0 - cl);
            local_weight += w * measure;
        }

        double global_weight = 0.0;
        MPI_Allreduce(&local_weight, &global_weight, 1, MPI_DOUBLE, MPI_SUM, PCW);

        if (global_weight <= DIFFUSE_MASS_EPS)
            continue;

        double local_applied = 0.0;

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
        for (int j = grid->G_Js; j < grid->G_Je; ++j)
        for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

            const double cs = vof->C_S[k][j][i];
            const double cl = vof->C_L[k][j][i];

            if (!diffuse_mass_redist_eligible(cl, cs, pass))
                continue;

            const double measure = diffuse_physical_cell_measure(grid, params, i, j, k);
            const double cap = (delta > 0.0) ? (1.0 - cs - cl) : cl;

            if (cap <= 0.0 || measure <= 0.0)
                continue;

            const double w = cl * (1.0 - cl);
            double dmass = delta * (w * measure) / (global_weight + DIFFUSE_MASS_EPS);
            double dcl   = dmass / measure;

            if (delta > 0.0) {
                if (dcl > cap) dcl = cap;
                if (dcl < 0.0) dcl = 0.0;
            } else {
                if (-dcl > cap) dcl = -cap;
                if (dcl > 0.0) dcl = 0.0;
            }

            vof->C_L[k][j][i] = cl + dcl;
            vof->C_G[k][j][i] = diffuse_clamp01(1.0 - vof->C_L[k][j][i] - cs);

            local_applied += dcl * measure;
        }

        double global_applied = 0.0;
        MPI_Allreduce(&local_applied, &global_applied, 1, MPI_DOUBLE, MPI_SUM, PCW);

        delta -= global_applied;

        VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
        VOF_DIFFUSE_update_phase_cache(db);
    }

    return delta;
}

#if defined(VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE) && defined(VOF_IBM)
/*
 * Host cells for the sediment-clip restore: the ONE-CELL DILATION of the grain,
 * i.e. non-sediment cells having at least one face neighbour that is sediment.
 *
 * An earlier version of this used the C_S threshold band 0 < C_S < cutoff,
 * on the reasoning that VOF_DIFFUSE_compute_C_S places C_S = cutoff = 0.05
 * exactly on r = R, so the sub-cutoff set "is the shell r > R, a few cells
 * thick".  The first half is right and the second half is WRONG, and the 2-D
 * A/B (jobs 20198503/20198504) caught it: C_S is a tanh whose support runs out
 * to diffuse_solid_support_extra(Cn) = -shift + denom*atanh(1 - 2e-12), which
 * for the release deck is 0.102 against R = 0.04.  The "rim" was therefore a
 * disc of radius 3.5 R covering 5.6% of the domain, and depositing into it
 * with F(1-F) weighting melted the entire ice slab (C_L in the far-field ice
 * went 0.000 -> 1.000 by t = 30).  Never key a locality argument on a tanh
 * tail; use the mesh.
 *
 * Reads of C_S at neighbours are safe without a halo exchange: compute_C_S
 * fills the ghosted range L_Is..L_Ie directly from particle geometry on every
 * rank whose subdomain intersects the solid support, so C_S ghosts are already
 * valid.  Only owned cells (G_*) are ever written.
 */
static int diffuse_sediment_adjacent(VolumeFraction *vof, MAC_grid *grid,
                                     int i, int j, int k)
{
    if (vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF)
        return 0;                                   /* inside the grain */

    if (i > grid->L_Is     && vof->C_S[k][j][i-1] >= DIFFUSE_SOLID_MASS_CUTOFF) return 1;
    if (i < grid->L_Ie - 1 && vof->C_S[k][j][i+1] >= DIFFUSE_SOLID_MASS_CUTOFF) return 1;
    if (j > grid->L_Js     && vof->C_S[k][j-1][i] >= DIFFUSE_SOLID_MASS_CUTOFF) return 1;
    if (j < grid->L_Je - 1 && vof->C_S[k][j+1][i] >= DIFFUSE_SOLID_MASS_CUTOFF) return 1;
    if (k > grid->L_Ks     && vof->C_S[k-1][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF) return 1;
    if (k < grid->L_Ke - 1 && vof->C_S[k+1][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF) return 1;

    return 0;
}

/*
 * Per-cell, per-stage ceiling on the restore.  The leak is a diffusive flux
 * across one cell face in one RK stage, so a host cell can never legitimately
 * be owed an O(1) change in C_L.  Without this ceiling, a host set whose total
 * weight is small gets the whole delta concentrated into a handful of cells and
 * saturates them -- which is the second half of how the first implementation
 * destroyed the ice slab.
 */
#define DIFFUSE_SED_MAX_DCL 0.05

/*
 * Return liquid that the admissibility clip trimmed inside a resolved sediment
 * grain (roadmap B.1.4) to the ring of cells immediately around that grain --
 * the cells the spurious Cahn-Hilliard flux actually took it from.
 *
 * Weighting is by C_L, not F(1-F): the flux can only have come from neighbours
 * that HAD liquid, so a bulk-ice neighbour must not receive any.  Cells with
 * essentially no liquid are excluded outright, and if NO adjacent cell has
 * liquid then there was no leak to undo and the delta is refused rather than
 * deposited -- a non-zero trim in that situation is initialisation transient or
 * round-off, and amplifying it is exactly the failure mode being fixed.
 *
 * Whatever cannot be placed is RETURNED AND DROPPED by the caller, never handed
 * to the global band restore.  Teleporting it across the domain is precisely
 * the hazard the CLIP_AUDIT block warns about, and in the 3-D shakeout
 * (job 20198505) the silent global fallback made a broken restore look healthy.
 * A large dropped fraction is a signal to escalate to masking the CH mobility
 * by C_S, not something to paper over.
 *
 * [MPI] one Allreduce for the weight, one for what was applied; every rank
 * divides the same two numbers, so the result is independent of the pencil
 * decomposition.  Only owned cells are written; neighbour reads of C_S stay
 * inside the ghosted range.
 */
static double diffuse_redistribute_sediment_clip_delta(Cart3d_bag *db,
                                                       double delta)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    if (fabs(delta) < 1.0e-14)
        return 0.0;

    double local_weight = 0.0;

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

        if (!diffuse_sediment_adjacent(vof, grid, i, j, k))
            continue;

        const double cl = vof->C_L[k][j][i];
        if (cl <= 0.005)
            continue;                    /* no liquid here to have supplied it */

        local_weight += cl * diffuse_physical_cell_measure(grid, params, i, j, k);
    }

    double global_weight = 0.0;
    MPI_Allreduce(&local_weight, &global_weight, 1, MPI_DOUBLE, MPI_SUM, PCW);

    if (global_weight <= DIFFUSE_MASS_EPS)
        return delta;                    /* nothing adjacent holds liquid */

    double local_applied = 0.0;

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

        if (!diffuse_sediment_adjacent(vof, grid, i, j, k))
            continue;

        const double cs = vof->C_S[k][j][i];
        const double cl = vof->C_L[k][j][i];
        if (cl <= 0.005)
            continue;

        const double measure = diffuse_physical_cell_measure(grid, params, i, j, k);
        if (measure <= 0.0)
            continue;

        double cap = (delta > 0.0) ? (1.0 - cs - cl) : cl;
        if (cap > DIFFUSE_SED_MAX_DCL) cap = DIFFUSE_SED_MAX_DCL;
        if (cap <= 0.0)
            continue;

        double dcl = delta * (cl * measure) / (global_weight + DIFFUSE_MASS_EPS)
                   / measure;

        if (delta > 0.0) {
            if (dcl > cap) dcl = cap;
            if (dcl < 0.0) dcl = 0.0;
        } else {
            if (-dcl > cap) dcl = -cap;
            if (dcl > 0.0) dcl = 0.0;
        }

        vof->C_L[k][j][i] = cl + dcl;
        vof->C_G[k][j][i] = diffuse_clamp01(1.0 - vof->C_L[k][j][i] - cs);

        local_applied += dcl * measure;
    }

    double global_applied = 0.0;
    MPI_Allreduce(&local_applied, &global_applied, 1, MPI_DOUBLE, MPI_SUM, PCW);

    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);

    return delta - global_applied;
}
#endif /* VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE && VOF_IBM */

static void diffuse_restore_liquid_mass(Cart3d_bag *db, double target_mass,
                                        const char *where)
{
    Parameters *params = db->params;

    double current = diffuse_liquid_mass_integral(db);
    double delta = target_mass - current;

    if (fabs(delta) < 1.0e-14)
        return;

    double residual = diffuse_redistribute_liquid_mass_delta(db, delta);
    diffuse_bound_liquid_fraction(db);
    VOF_DIFFUSE_set_boundary_values(db->vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);

#ifdef DEBUG_VOF_DIFFUSE_MASS
    char msg[300];
    snprintf(msg, sizeof(msg),
             "VOF_DIFFUSE mass restore after %s: target=%.16e before=%.16e delta=%.3e residual=%.3e\n",
             where, target_mass, current, delta, residual);
    Display_progress(params, msg);
#else
    (void)params;
    (void)where;
    (void)residual;
#endif
}

#if defined(PHASE_CHANGE) && defined(CONC)
/*
 * Stefan melt-rate source for the Cahn-Hilliard equation (roadmap G.2-G.3),
 * in the superheat-driven phase-field form of Hester et al. (2020) used by
 * Yang et al. (2023):
 *
 *   m = V_G * w(F) ,
 *   V_G  = (St / Pe_T) * ( theta - theta_L(s) ) / melt_band_eps ,
 *   w(F) = F(1-F) / (sqrt(2)*Cn) ,
 *   theta_L(s) = T_melt - liquidus_slope * s      (freezing-point depression)
 *
 * w(F) equals |grad F| on the equilibrium tanh profile (same localization,
 * same unit integral across the band), but unlike the measured gradient it
 * vanishes quadratically at the band tails.  That is essential: the CH
 * relaxation is far slower than the deposition, so a source that keeps
 * loading nearly-saturated cells (F -> 1) pushes them out of [0,1] and the
 * admissibility clip in diffuse_bound_liquid_fraction then silently destroys
 * the melt mass the temperature field already paid latent heat for (observed
 * as a ~10% enthalpy leak and a 4% slow front in the A.3 gate).  With w(F)
 * the clip stays a no-op and the discrete enthalpy theta + F/St changes only
 * through boundary fluxes.
 *
 * m > 0 (local superheat) converts ice into liquid; m < 0 refreezes.  The
 * same discrete field enters the temperature equation as the latent sink
 * -(1/St)*m (Conc_add_latent_heat_RHS); the negative feedback pins the band
 * temperature to the liquidus, and the interface superheat - and with it the
 * front-speed error - vanishes as O(melt_band_eps).  A one-sided flux-jump
 * evaluation of V_G is NOT usable here: on an unpinned diffuse band the
 * temperature smooths out and the liquid/solid one-sided gradients cancel.
 *
 * Called once per RK stage from F^{k-1} and theta^{k-1}, before the CH
 * transport updates C_L.  The previous-stage source is kept in melt_src_old
 * for the ZETA history term of both equations.
 *
 * [MPI] fully pointwise (no stencil): F, theta, s are read at the local cell
 * only, so no halo refresh is needed and the result is rank-count independent.
 */
void VOF_DIFFUSE_compute_melt_rate(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;
    double       ***theta  = db->c[0]->data;
    double       ***F      = vof->C_L;

    const double coeff = params->stefan
                       / (params->Pe[0] * params->melt_band_eps);
    const double band_norm = 1.0 / (sqrt(2.0) * params->Cn);
    const double T_melt = params->T_melt;
    const double m_liq  = params->liquidus_slope;
    const int has_salt  = (params->NConc > 1);
    const int use_sliq  = params->liquid_referenced_salinity && has_salt;

    Array_copy_withghost(vof->melt_src, vof->melt_src_old, grid, params);

    if (params->stefan == 0.0) {
        Array_set_withghost(vof->melt_src, 0.0, grid, params);
        return;
    }

    Array_set_withghost(vof->melt_src, 0.0, grid, params);

    const int i_end = min(grid->NX - 1, grid->G_Ie);
    const int j_end = min(grid->NY - 1, grid->G_Je);
    const int k_end = min(grid->NZ - 1, grid->G_Ke);

    for (int k = grid->G_Ks; k < k_end; ++k) {
        for (int j = grid->G_Js; j < j_end; ++j) {
            for (int i = grid->G_Is; i < i_end; ++i) {

#ifdef VOF_IBM
                /* No Stefan melting of the resolved (sediment) solid. */
                if (vof->C_S[k][j][i] >= DIFFUSE_SOLID_MASS_CUTOFF)
                    continue;
#endif

                double f = F[k][j][i];
                double w = f * (1.0 - f) * band_norm;
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
                /* Sediment excludes volume but is not a meltable phase. */
                const double pore=1.0-vof->C_S[k][j][i];
                w=f*fmax(0.0,pore-f)/(pore+1.0e-30)*band_norm;
#endif
                if (w <= 0.0)
                    continue;

                double s_loc = has_salt ? db->c[1]->data[k][j][i] : 0.0;
                /*
                 * Roadmap A.5.5 variant 1 (liquid-referenced salinity).
                 *
                 * s is a VOLUME-averaged salinity over a cell that may contain
                 * both phases, but the salt lives only in the liquid, so the
                 * salinity the liquidus actually sees is s/(F+delta).  In the
                 * bulk (F=1) this is identity; it acts only where 0 < F < 1
                 * AND salt is present -- i.e. exactly and only in the salty
                 * interfacial film, with zero effect on the freshwater case.
                 * That selectivity is why it is the leading candidate for the
                 * salty melt-rate gap, which the 2026-07-28 resolution test
                 * showed is NOT a resolution artifact.
                 *
                 * Default off: liquid_referenced_salinity = 0 reproduces the
                 * previous behaviour exactly.
                 */
                if (use_sliq && s_loc != 0.0)
                    s_loc /= (f + params->yang_salt_delta);
                double theta_L = T_melt - m_liq * s_loc;

                vof->melt_src[k][j][i] =
                    coeff * (theta[k][j][i] - theta_L) * w;
            }
        }
    }
}
#endif /* PHASE_CHANGE && CONC */

void VOF_DIFFUSE_step(Cart3d_bag *db)
{
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
    diffuse_step_ice(db);
    return;
#endif

    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    const double BET[]  = { BETA };
    const double GAMB[] = { GAMBETA };
    const double ZETB[] = { ZETBETA };

    const int which_stage = params->which_stage;
    const double a_dt = 1.0 / (params->dt * BET[which_stage]);

#if VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
    double mass_target = diffuse_liquid_mass_integral(db);
#endif

    double ***rhs_cur = vof->ch_aux2;
    double ***rhs_tmp = vof->ch_aux1;
    double ***rhs_vec = vof->ch_rhs_nm1;

    /*
     * Explicit CH RHS:
     *   R(C) = -div(u C) + Pe^{-1} Laplacian(S(C))  [+ V_G |grad F| melting]
     */
#if defined(PHASE_CHANGE) && defined(CONC)
    /* Stefan source from F^{k-1}, theta^{k-1}; must precede the C_L update. */
    VOF_DIFFUSE_compute_melt_rate(db);
#endif

    VOF_DIFFUSE_advect_WENO5(db, rhs_cur);
    VOF_DIFFUSE_explicit_diffusion(db, rhs_tmp);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
        rhs_cur[k][j][i] += rhs_tmp[k][j][i];
#if defined(PHASE_CHANGE) && defined(CONC)
        /* Rides through the same GAMMA/ZETA RK combination and the ch_rhs_n
         * history as the advective and diffusive parts. */
        rhs_cur[k][j][i] += vof->melt_src[k][j][i];
#endif
    }

    VOF_DIFFUSE_set_boundary_values(rhs_cur, db);

    if (which_stage == 0) {
        for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
        for (int j = grid->G_Js; j < grid->G_Je; ++j)
        for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
            rhs_vec[k][j][i] =
                a_dt * vof->C_L[k][j][i] +
                GAMB[0] * rhs_cur[k][j][i];
        }
    } else {
        for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
        for (int j = grid->G_Js; j < grid->G_Je; ++j)
        for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
            rhs_vec[k][j][i] =
                a_dt * vof->C_L[k][j][i] +
                GAMB[which_stage] * rhs_cur[k][j][i] +
                ZETB[which_stage] * vof->ch_rhs_n[k][j][i];
        }
    }

    VOF_DIFFUSE_set_boundary_values(rhs_vec, db);

    /*
     * Implicit biharmonic part.
     */
    VOF_DIFFUSE_solve_implicit_biharmonic(db, rhs_vec);

#if defined(PHASE_CHANGE) && defined(CONC)
    /*
     * Conservative admissibility clip.  The melt mass was already paid for by
     * the latent sink, so whatever the [0,1] clip removes (CH transient
     * overshoot behind the advancing front) or creates must be returned to
     * the interfacial band; a plain clip is a hidden enthalpy leak that
     * biases the Stefan front speed (measured ~5% in the A.3 gate).
     * diffuse_redistribute_liquid_mass_delta is MPI-collective and deposits
     * with F(1-F) weighting into cells with headroom.
     */
    {
        double clip_residual = 0.0;

#if DIFFUSE_TRACK_SEDIMENT_CLIP
        double sed_delta = 0.0;
        double sed_residual = 0.0;
        double clip_delta = diffuse_bound_liquid_fraction_split(db, &sed_delta);
#else
        double clip_delta = diffuse_bound_liquid_fraction(db);
#endif

        if (fabs(clip_delta) > 1.0e-14)
            clip_residual = diffuse_redistribute_liquid_mass_delta(db, -clip_delta);

#if defined(VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE) && defined(VOF_IBM)
        /*
         * Roadmap B.1.4 -- the sediment half of the same conservative clip.
         *
         * Liquid the trim removes inside a resolved grain is real mass whose
         * latent heat was already paid, exactly as on the fluid side; the only
         * difference is that it left the accounting domain of
         * diffuse_liquid_mass_integral (which excludes sediment cells) rather
         * than being clipped within it.  Left alone it is destroyed, at a rate
         * that scales with sediment surface area -- 2.71x the grain volume by
         * t=30 in the 2-D release deck, and worse in 3-D and multi-grain.
         *
         * It is put back into the ring of cells immediately around the grain,
         * i.e. where the spurious CH flux took it from.  What will not fit
         * there is DROPPED, deliberately: handing it to the global band restore
         * teleports it across the domain, and in the 3-D shakeout that silent
         * fallback made a broken restore look healthy.
         */
        if (fabs(sed_delta) > 1.0e-14)
            sed_residual = diffuse_redistribute_sediment_clip_delta(db, -sed_delta);
#endif

#if DIFFUSE_TRACK_SEDIMENT_CLIP && defined(VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT)
        /*
         * Roadmap B.1.4 instrumentation -- log only, no physics.  Measures the
         * leak directly instead of inferring it from the meltwater-tracer
         * budget (which conflates it with tracer transport error).
         *
         * Columns:
         *   net       SIGNED running sum -- the mass actually destroyed (or,
         *             with the restore on, still unaccounted).  This is the
         *             physically meaningful number; |.| sums overstate it
         *             badly because the per-stage trims largely cancel.
         *   cum|sed|  running sum of |trim|, an upper bound / activity measure.
         *   dropped   what the adjacent ring could not absorb and was thrown
         *             away rather than teleported.  If this is not small
         *             against net, option (a) is insufficient and the fallback
         *             is masking the CH mobility by C_S (see B.1.4-fix).
         */
        {
            static double sed_abs_sum = 0.0, sed_resid_abs_sum = 0.0;
            static double sed_net_sum = 0.0, sed_dropped_sum = 0.0;
            static double sed_melt_abs_sum = 0.0;

            double local_m = 0.0;
            for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
            for (int j = grid->G_Js; j < grid->G_Je; ++j)
            for (int i = grid->G_Is; i < grid->G_Ie; ++i)
                local_m += fabs(vof->melt_src[k][j][i])
                         * diffuse_physical_cell_measure(grid, params, i, j, k);

            double global_m = 0.0;
            MPI_Allreduce(&local_m, &global_m, 1, MPI_DOUBLE, MPI_SUM, PCW);

            sed_abs_sum       += fabs(sed_delta);
            sed_net_sum       += sed_delta;
            sed_resid_abs_sum += fabs(sed_residual);
            sed_dropped_sum   += sed_residual;
            sed_melt_abs_sum  += global_m * params->dt;

            if (which_stage == 2 &&
                (params->ntime % VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT == 0)) {
                char msg[420];
                snprintf(msg, sizeof(msg),
                         "SED_CLIP_AUDIT ntime=%d sed=%.3e resid=%.3e "
                         "net=%.6e cum|sed|=%.6e dropped=%.6e "
                         "cum|resid|=%.6e cum_melt=%.6e ratio=%.3e\n",
                         params->ntime, sed_delta, sed_residual,
                         sed_net_sum, sed_abs_sum, sed_dropped_sum,
                         sed_resid_abs_sum, sed_melt_abs_sum,
                         sed_melt_abs_sum > 0.0
                             ? fabs(sed_net_sum) / sed_melt_abs_sum : 0.0);
                Display_progress(params, msg);
            }
        }
#elif DIFFUSE_TRACK_SEDIMENT_CLIP
        (void)sed_residual;
#endif

#ifdef VOF_DIFFUSE_CLIP_AUDIT
        /*
         * Roadmap A.5.1 item 4 -- log only, no physics.
         *
         * diffuse_redistribute_liquid_mass_delta is GLOBALLY collective: mass
         * the [0,1] clip removes here is redeposited anywhere on the interface
         * that has F(1-F) headroom, not where it was taken from.  In 1-D that
         * is provably a no-op.  In the 2-D salty run the row-melt CV is
         * 0.35-0.64, so a non-zero clip_delta would teleport melt from fast
         * scallop rows to slow ones and flatten exactly the interfacial
         * dynamics that set convective heat delivery.
         *
         * What to look for: |clip| / melt should be <= O(1e-6).  If it is not,
         * the fix is a band-local restore (cf. VOF_DIFFUSE_BAND_MASS_RESTORE /
         * diffuse_band_conserve_contact_line in the trapped-cavity work), not
         * a tolerance change.  A non-zero residual means the redistribution
         * could not place all the mass, which is strictly worse.
         */
        {
            static double clip_abs_sum = 0.0, resid_abs_sum = 0.0, melt_abs_sum = 0.0;
            double local_melt = 0.0;
            for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
            for (int j = grid->G_Js; j < grid->G_Je; ++j)
            for (int i = grid->G_Is; i < grid->G_Ie; ++i)
                local_melt += fabs(vof->melt_src[k][j][i])
                            * diffuse_physical_cell_measure(grid, params, i, j, k);

            double global_melt = 0.0;
            MPI_Allreduce(&local_melt, &global_melt, 1, MPI_DOUBLE, MPI_SUM, PCW);

            clip_abs_sum  += fabs(clip_delta);
            resid_abs_sum += fabs(clip_residual);
            melt_abs_sum  += global_melt * params->dt;

            if (which_stage == 2 && (params->ntime % VOF_DIFFUSE_CLIP_AUDIT == 0)) {
                char msg[320];
                snprintf(msg, sizeof(msg),
                         "CLIP_AUDIT ntime=%d clip=%.3e residual=%.3e "
                         "cum|clip|=%.6e cum|resid|=%.6e cum_melt=%.6e ratio=%.3e\n",
                         params->ntime, clip_delta, clip_residual,
                         clip_abs_sum, resid_abs_sum, melt_abs_sum,
                         melt_abs_sum > 0.0 ? clip_abs_sum / melt_abs_sum : 0.0);
                Display_progress(params, msg);
            }
        }
#else
        (void)clip_residual;
#endif
    }
#else
    diffuse_bound_liquid_fraction(db);
#endif
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);


    const int last_stage = (which_stage == 2);   /* LSRK3 has 3 substages, idx 0..2 */


#if VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
    diffuse_restore_liquid_mass(db, mass_target, "RK3 soft clamp");
#endif

    if (last_stage) {


    #ifdef VOF_IBM
        /* Pi_theta: characteristic MCL projection, also once per dt */
        VOF_DIFFUSE_apply_contact_angle(db);
        VOF_DIFFUSE_set_boundary_values(vof->C_L, db);

            /* Pi_adm: C_L + C_S <= 1 */
    #if defined(VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE) && defined(PHASE_CHANGE) && defined(CONC)
        /*
         * Roadmap B.1.4 -- second leak channel.  The contact-angle projection
         * can push C_L back over 1 - C_S, and this trim discarded that mass
         * too.  Same treatment as the main clip: rim-local restore, global
         * fallback.  Only the SEDIMENT half is restored here; the fluid half
         * is deliberately left as-is so the Stage-A/VOF_IBM behaviour of this
         * call site is otherwise unchanged.
         */
        {
            double pi_sed = 0.0;
            diffuse_bound_liquid_fraction_split(db, &pi_sed);
            if (fabs(pi_sed) > 1.0e-14)
                (void)diffuse_redistribute_sediment_clip_delta(db, -pi_sed);
        }
    #else
        diffuse_bound_liquid_fraction(db);
    #endif
        VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
        VOF_DIFFUSE_update_phase_cache(db);
    #endif

    #if VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
        diffuse_restore_liquid_mass(db, mass_target, "RK3 contact-angle projection");
    #endif
    }



    Array_copy_withghost(vof->ch_rhs_n, vof->ch_rhs_nm1, grid, params);
    Array_copy_withghost(rhs_cur, vof->ch_rhs_n, grid, params);
}

#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
/* Evolve real ice outside the resolved sediment mask.  The complementary
 * variable makes zero ice an exact invariant of transport: a moving rock
 * cannot nucleate another phase.  CH uses the same symmetric no-flux operator;
 * its polynomial is the binary ice/water potential.  Fluxes into sediment are
 * zero. Positive melt remains positive in melt_src for heat/tracer, but enters
 * ice with the opposite sign. ch_rhs_n stores the ICE RHS in this build.
 */
static void diffuse_step_ice(Cart3d_bag *db)
{
    MAC_grid *g=db->grid;
    Parameters *p=db->params;
    VolumeFraction *v=db->vof;
    const double BET[]={BETA},GAMB[]={GAMBETA},ZETB[]={ZETBETA};
    const int st=p->which_stage;
    const double a=1.0/(p->dt*BET[st]);
    double ***cur=v->ch_aux2,***tmp=v->ch_aux1,***rhs=v->ch_rhs_nm1;
    VOF_DIFFUSE_compute_melt_rate(db);
    double local_source=0.0,global_source=0.0;
    for (int k=g->G_Ks;k<g->G_Ke;++k)
    for (int j=g->G_Js;j<g->G_Je;++j)
    for (int i=g->G_Is;i<g->G_Ie;++i)
        local_source += p->dt*BET[st]*(GAMB[st]*v->melt_src[k][j][i]
            + ZETB[st]*v->melt_src_old[k][j][i])*diffuse_physical_cell_measure(g,p,i,j,k);
    MPI_Allreduce(&local_source,&global_source,1,MPI_DOUBLE,MPI_SUM,PCW);
    static double melt_integral=0.0;
    melt_integral+=global_source;
    if(st==2 && p->ntime%100==0) {
        char msg[200];
        snprintf(msg,sizeof(msg),"MELT_INTEGRAL ntime=%d time=%.12e source=%.12e\n",p->ntime,p->time,melt_integral);
        Display_progress(p,msg);
    }
    for (int k=g->G_Ks;k<g->G_Ke;++k)
    for (int j=g->G_Js;j<g->G_Je;++j)
    for (int i=g->G_Is;i<g->G_Ie;++i)
        v->C_L[k][j][i]=v->C_S[k][j][i]>=DIFFUSE_SOLID_MASS_CUTOFF ? 0.0 :
            fmax(0.0,1.0-v->C_S[k][j][i]-v->C_L[k][j][i]);
    VOF_DIFFUSE_set_boundary_values(v->C_L,db);
    diffuse_evolving_ice=1;
    VOF_DIFFUSE_advect_WENO5(db,cur);
    VOF_DIFFUSE_explicit_diffusion(db,tmp);
    for (int k=g->G_Ks;k<g->G_Ke;++k)
    for (int j=g->G_Js;j<g->G_Je;++j)
    for (int i=g->G_Is;i<g->G_Ie;++i) {
        cur[k][j][i]+=tmp[k][j][i]-v->melt_src[k][j][i];
        if (v->C_S[k][j][i]>=DIFFUSE_SOLID_MASS_CUTOFF) cur[k][j][i]=0.0;
        rhs[k][j][i]=a*v->C_L[k][j][i]+GAMB[st]*cur[k][j][i];
        if(st) rhs[k][j][i]+=ZETB[st]*v->ch_rhs_n[k][j][i];
    }
    VOF_DIFFUSE_set_boundary_values(cur,db);
    VOF_DIFFUSE_set_boundary_values(rhs,db);
    VOF_DIFFUSE_solve_implicit_biharmonic(db,rhs);
    double delta=diffuse_bound_liquid_fraction(db),residual=0.0;
    if(fabs(delta)>1.e-14) residual=diffuse_redistribute_liquid_mass_delta(db,-delta);
    static double clip=0.0,unplaced=0.0;
    clip+=delta;unplaced+=residual;
    if(st==2 && p->ntime%100==0) {
        char msg[200];
        snprintf(msg,sizeof(msg),"ICE_TRANSPORT ntime=%d clip=%.12e unplaced=%.12e\n",p->ntime,clip,unplaced);
        Display_progress(p,msg);
    }
    for (int k=g->G_Ks;k<g->G_Ke;++k)
    for (int j=g->G_Js;j<g->G_Je;++j)
    for (int i=g->G_Is;i<g->G_Ie;++i)
        v->C_L[k][j][i]=v->C_S[k][j][i]>=DIFFUSE_SOLID_MASS_CUTOFF ? 0.0 :
            1.0-v->C_S[k][j][i]-v->C_L[k][j][i];
    diffuse_evolving_ice=0;
    VOF_DIFFUSE_set_boundary_values(v->C_L,db);
    VOF_DIFFUSE_update_phase_cache(db);
    VOF_DIFFUSE_set_boundary_values(v->C_G,db);
    Array_copy_withghost(v->ch_rhs_n,v->ch_rhs_nm1,g,p);
    Array_copy_withghost(cur,v->ch_rhs_n,g,p);
}
#endif



#ifdef SURFACE_TENSION
void VOF_DIFFUSE_compute_f_sigma(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;
    const int collapsed_z = TwodOps_collapsed_component_is_inactive(params);

    /* Match the code's -2*grad(p) convention used in the momentum RHS. */
    const double scale = 2.0 * 6.0 * sqrt(2.0) / (params->Cn * params->We + 1e-30);

    

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
        vof->f_sigma_new_x[k][j][i] = 0.0;
        vof->f_sigma_new_y[k][j][i] = 0.0;
        vof->f_sigma_new_z[k][j][i] = 0.0;
    }

    int i_start_u = max(1, grid->G_Is),      j_start_u = grid->G_Js,            k_start_u = grid->G_Ks;
    int i_end_u   = min(grid->NX-1, grid->G_Ie), j_end_u = min(grid->NY-1, grid->G_Je), k_end_u = min(grid->NZ-1, grid->G_Ke);
#ifdef XPERIODIC
    if (i_end_u == grid->NX-1) i_end_u = grid->NX;
#endif

    int i_start_v = grid->G_Is,              j_start_v = max(1, grid->G_Js),    k_start_v = grid->G_Ks;
    int i_end_v   = min(grid->NX-1, grid->G_Ie), j_end_v = min(grid->NY-1, grid->G_Je), k_end_v = min(grid->NZ-1, grid->G_Ke);
#ifdef YPERIODIC
    if (j_end_v == grid->NY-1) j_end_v = grid->NY;
#endif

    int i_start_w = grid->G_Is,              j_start_w = grid->G_Js,             k_start_w = max(1, grid->G_Ks);
    int i_end_w   = min(grid->NX-1, grid->G_Ie), j_end_w = min(grid->NY-1, grid->G_Je), k_end_w = min(grid->NZ-1, grid->G_Ke);
#ifdef ZPERIODIC
    if (k_end_w == grid->NZ-1) k_end_w = grid->NZ;
#endif

    for (int k = k_start_u; k < k_end_u; ++k)
    for (int j = j_start_u; j < j_end_u; ++j)
    for (int i = i_start_u; i < i_end_u; ++i) {
        double cl = 0.5 * (vof->C_L[k][j][i] + vof->C_L[k][j][i-1]);
        double cg = 0.5 * (vof->C_G[k][j][i] + vof->C_G[k][j][i-1]);
        double cs = 0.5 * (vof->C_S[k][j][i] + vof->C_S[k][j][i-1]);
        if (cl < 0.005 || cl > 0.995 || cg < 0.005 || cg > 0.995 || cs > 0.05) {
            continue;
        }
        double psi_f = 0.5 * (vof->psi_LG[k][j][i] + vof->psi_LG[k][j][i-1]);
        double gradC = (vof->C_L[k][j][i] - vof->C_L[k][j][i-1]) * grid->idx_c[i-1];
        vof->f_sigma_new_x[k][j][i] = scale * psi_f * gradC;
    }

    for (int k = k_start_v; k < k_end_v; ++k)
    for (int j = j_start_v; j < j_end_v; ++j)
    for (int i = i_start_v; i < i_end_v; ++i) {
        double cl = 0.5 * (vof->C_L[k][j][i] + vof->C_L[k][j-1][i]);
        double cg = 0.5 * (vof->C_G[k][j][i] + vof->C_G[k][j-1][i]);
        double cs = 0.5 * (vof->C_S[k][j][i] + vof->C_S[k][j-1][i]);
        if (cl < 0.005 || cl > 0.995 || cg < 0.005 || cg > 0.995 || cs > 0.05) {
            continue;
        }
        double psi_f = 0.5 * (vof->psi_LG[k][j][i] + vof->psi_LG[k][j-1][i]);
        double gradC = (vof->C_L[k][j][i] - vof->C_L[k][j-1][i]) * grid->idy_c[j-1];
        vof->f_sigma_new_y[k][j][i] = scale * psi_f * gradC;
    }

    if (!collapsed_z) {
        for (int k = k_start_w; k < k_end_w; ++k)
        for (int j = j_start_w; j < j_end_w; ++j)
        for (int i = i_start_w; i < i_end_w; ++i) {
            double cl = 0.5 * (vof->C_L[k][j][i] + vof->C_L[k-1][j][i]);
            double cg = 0.5 * (vof->C_G[k][j][i] + vof->C_G[k-1][j][i]);
            double cs = 0.5 * (vof->C_S[k][j][i] + vof->C_S[k-1][j][i]);
            if (cl < 0.005 || cl > 0.995 || cg < 0.005 || cg > 0.995 || cs > 0.05) {
                continue;
            }
            double psi_f = 0.5 * (vof->psi_LG[k][j][i] + vof->psi_LG[k-1][j][i]);
            double gradC = (vof->C_L[k][j][i] - vof->C_L[k-1][j][i]) * grid->idz_c[k-1];
            vof->f_sigma_new_z[k][j][i] = scale * psi_f * gradC;
        }
    }
}

void VOF_DIFFUSE_apply_f_sigma_old(Cart3d_bag *db)
{
    VOF_apply_f_sigma_old(db);
}
#endif

void VOF_DIFFUSE_update_density_viscosity(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    /*
     * Liu15 Eq. (6)-(7):
     *
     *   rho = C_L + lambda_rho,S C_S + lambda_rho,G C_G
     *   mu  = C_L + lambda_mu,S  C_S + lambda_mu,G  C_G
     *
     * with lambda_rho,S = lambda_mu,S = 1 for the stationary solid carrier.
     */
    const double rho_g = params->rho2 / params->rho1;
    const double mu_g  = params->mu2  / params->mu1;
    const double rho_s = 1.0;
    const double mu_s  = 1.0;

    

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {



                vof->rho[k][j][i] = vof->C_L[k][j][i] + rho_s * vof->C_S[k][j][i] + rho_g * vof->C_G[k][j][i];
                vof->mu[k][j][i]  = vof->C_L[k][j][i] + mu_s  * vof->C_S[k][j][i] + mu_g  * vof->C_G[k][j][i];
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->rho, db);
    VOF_DIFFUSE_set_boundary_values(vof->mu, db);
}



#endif
