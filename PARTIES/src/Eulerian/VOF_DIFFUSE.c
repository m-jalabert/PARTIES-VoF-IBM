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

#define DIFFUSE_SOLID_MASS_CUTOFF 0.05
#define DIFFUSE_MASS_EPS 1.0e-30
#define DIFFUSE_SOLID_GEOM_CUTOFF 1.0e-12

#ifndef VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
#define VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE 1
#endif

#ifdef VOF_DIFFUSE

static double diffuse_liquid_mass_integral(Cart3d_bag *db);
static double diffuse_bound_liquid_fraction(Cart3d_bag *db);
static void diffuse_restore_liquid_mass(Cart3d_bag *db, double target_mass,
                                        const char *where);
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
    Particle *particles =
        Particle_collect_owned_overlaps(p_list, db, extra_range, -1.0,
                                        1, &n_particles);

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

static double diffuse_inner_product(double ***a, double ***b, MAC_grid *grid)
{
    double local_sum = 0.0;
    double global_sum = 0.0;

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                local_sum += a[k][j][i] * b[k][j][i];
            }
        }
    }

    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);
    return global_sum;
}

static double diffuse_norm(double ***a, MAC_grid *grid)
{
    return sqrt(diffuse_inner_product(a, a, grid));
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

static inline double diffuse_diag_entry(MAC_grid *grid, int i, int j, int k, double lambda)
{
    double idx = grid->idx_c[i];
    double idy = grid->idy_c[j];
    double idz = grid->idz_c[k];

    return 1.0 + 2.0 * lambda * (idx * idx + idy * idy + idz * idz);
}

static void diffuse_apply_helmholtz(double ***x, double ***Ax, double lambda,
                                    Cart3d_bag *db)
{
    MAC_grid *grid = db->grid;

    VOF_DIFFUSE_set_boundary_values(x, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double idx = grid->idx_c[i];
                double idy = grid->idy_c[j];
                double idz = grid->idz_c[k];

                double lap =
                    (x[k][j][i+1] - 2.0 * x[k][j][i] + x[k][j][i-1]) * idx * idx +
                    (x[k][j+1][i] - 2.0 * x[k][j][i] + x[k][j-1][i]) * idy * idy +
                    (x[k+1][j][i] - 2.0 * x[k][j][i] + x[k-1][j][i]) * idz * idz;

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
                z[k][j][i] = r[k][j][i] / diffuse_diag_entry(grid, i, j, k, lambda);
                p[k][j][i] = z[k][j][i];
            }
        }
    }

    double rhs_norm = diffuse_norm(b, grid);
    if (rhs_norm < 1e-30) {
        rhs_norm = 1.0;
    }

    double rz = diffuse_inner_product(r, z, grid);
    double rel_res = diffuse_norm(r, grid) / rhs_norm;
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

        double pAp = diffuse_inner_product(p, Ap, grid);
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
        rel_res = diffuse_norm(r, grid) / rhs_norm;
        if (rel_res < tol) {
            break;
        }

        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    z[k][j][i] = r[k][j][i] / diffuse_diag_entry(grid, i, j, k, lambda);
                }
            }
        }

        double rz_new = diffuse_inner_product(r, z, grid);
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

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                Ax[k][j][i] = a_dt * x[k][j][i] + bih_coeff * bih_x[k][j][i];
            }
        }
    }
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

    double rhs_norm = diffuse_norm(b, grid);
    if (rhs_norm < 1e-30) {
        rhs_norm = 1.0;
    }

    double rz = diffuse_inner_product(r, z, grid);
    double rel_res = diffuse_norm(r, grid) / rhs_norm;
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

        double pAp = diffuse_inner_product(p, Ap, grid);
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
        rel_res = diffuse_norm(r, grid) / rhs_norm;
        if (rel_res < tol) {
            break;
        }

        double rz_new = diffuse_inner_product(r, z, grid);
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

    /*
     * Make Data_0 the same projected state that will be advanced by the first
     * CH/MCL step.  Without this, Data_0 is the raw geometric initialization,
     * while the first stage immediately enforces C_L + C_S <= 1 and the Liu MCL
     * ghost-contact-line update.  That creates an artificial jump between the
     * initialization output and the first physical step.
     *
     * The target mass is the raw initialized liquid mass over the non-solid
     * region, matching Liu's S = integral_{C_S<0.05} C_L dOmega definition.
     * The projection is conservative: any mass removed by bounding/limiting is
     * redistributed into admissible non-solid capacity.
     */
    double initial_mass = diffuse_liquid_mass_integral(db);

    diffuse_bound_liquid_fraction(db);
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);
    diffuse_restore_liquid_mass(db, initial_mass, "initial C_L+CS projection");

#ifdef VOF_IBM
    VOF_DIFFUSE_apply_contact_angle(db);
    diffuse_restore_liquid_mass(db, initial_mass, "initial MCL projection");
#endif

    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);
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
        for (int g = 1; g <= NG; ++g) {
#ifdef RIGHT_WALL_VELOCITY_FREESLIP
            f[k][j][i1 + g] = f[k][j][i1 - g + 1];
#else
            f[k][j][i1 + g] = f[k][j][i1];
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
        for (int g = 1; g <= NG; ++g) {
#ifdef TOP_WALL_VELOCITY_FREESLIP
            f[k][j1 + g][i] = f[k][j1 - g + 1][i];
#else
            f[k][j1 + g][i] = f[k][j1][i];
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
        for (int g = 1; g <= NG; ++g) {
#ifdef FRONT_WALL_VELOCITY_FREESLIP
            f[k1 + g][j][i] = f[k1 - g + 1][j][i];
#else
            f[k1 + g][j][i] = f[k1][j][i];
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
                const double dz = grid->zc[k] - X[2];

                for (int j = grid->L_Js; j < grid->L_Je; ++j) {
                    const double dy = grid->yc[j] - X[1];

                    for (int i = grid->L_Is; i < grid->L_Ie; ++i) {
                        const double dx = grid->xc[i] - X[0];

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

    VOF_DIFFUSE_update_phase_cache(db);
}

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
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(lap, db);
}

void VOF_DIFFUSE_compute_bulk_S(Cart3d_bag *db)
{
    MAC_grid       *grid = db->grid;
    VolumeFraction *vof  = db->vof;

    VOF_DIFFUSE_update_phase_cache(db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double c = vof->C_L[k][j][i];
                double cs = vof->C_S[k][j][i];
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

#ifdef VOF_IBM
    VOF_DIFFUSE_extend_psi_LG_contact_angle(db);
#endif
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

                vel = v->data[k][j][i];
                state = diffuse_face_state_y(vof->C_L, i, j, k, vel, params->weno_order);
                flux_y[k][j][i] = vel * state;

                if (!collapsed_z) {
                    vel = w->data[k][j][i];
                    state = diffuse_face_state_z(vof->C_L, i, j, k, vel, params->weno_order);
                    flux_z[k][j][i] = vel * state;
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

    iters = diffuse_solve_biharmonic_pcg(vof->C_L, rhs_explicit,
                                         a_dt, bih_coeff, &final_res, db);
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);

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
                 TwodOps_cell_measure_c(grid, params, i, j, k);
    }

    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, PCW);
    return global;
}

static double diffuse_bound_liquid_fraction(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    double local_delta = 0.0;

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

        if (cs < DIFFUSE_SOLID_MASS_CUTOFF) {
            double measure = TwodOps_cell_measure_c(grid, params, i, j, k);
            local_delta += (new_cl - old_cl) * measure;
        }
    }

    double global_delta = 0.0;
    MPI_Allreduce(&local_delta, &global_delta, 1, MPI_DOUBLE, MPI_SUM, PCW);
    return global_delta;
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

            const double measure = TwodOps_cell_measure_c(grid, params, i, j, k);

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

            const double measure = TwodOps_cell_measure_c(grid, params, i, j, k);
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

void VOF_DIFFUSE_step(Cart3d_bag *db)
{
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
     *   R(C) = -div(u C) + Pe^{-1} Laplacian(S(C))
     */
    VOF_DIFFUSE_advect_WENO5(db, rhs_cur);
    VOF_DIFFUSE_explicit_diffusion(db, rhs_tmp);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
        rhs_cur[k][j][i] += rhs_tmp[k][j][i];
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


    const int last_stage = (which_stage == 2);   /* LSRK3 has 3 substages, idx 0..2 */

    /*
    * Soft per-substage clamp.  Numerically keeps C_L in [0,1] for the next CH
    * stencil; does NOT enforce C_L + C_S <= 1 — that is the once-per-dt job of
    * the contact-angle / admissibility projection below.
    */
    for (int k = grid->G_Ks; k < grid->G_Ke; ++k)
    for (int j = grid->G_Js; j < grid->G_Je; ++j)
    for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
        double c = vof->C_L[k][j][i];
        if (c < 0.0) c = 0.0;
        if (c > 1.0) c = 1.0;
        vof->C_L[k][j][i] = c;
    }
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    VOF_DIFFUSE_update_phase_cache(db);

#if VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE
    diffuse_restore_liquid_mass(db, mass_target, "RK3 soft clamp");
#endif

    if (last_stage) {
        /* Pi_adm: C_L + C_S <= 1, applied once per dt as in Liu-Ding Step 2 */
        diffuse_bound_liquid_fraction(db);
        VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
        VOF_DIFFUSE_update_phase_cache(db);

    #ifdef VOF_IBM
        /* Pi_theta: characteristic MCL projection, also once per dt */
        VOF_DIFFUSE_apply_contact_angle(db);
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



#ifdef SURFACE_TENSION
void VOF_DIFFUSE_compute_f_sigma(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;
    const int collapsed_z = TwodOps_collapsed_component_is_inactive(params);

    /* Match the code's -2*grad(p) convention used in the momentum RHS. */
    const double scale = 2.0 * 6.0 * sqrt(2.0) / (params->Cn * params->We + 1e-30);

    VOF_DIFFUSE_update_phase_cache(db);

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

    VOF_DIFFUSE_update_phase_cache(db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {

                double cl = clampDouble(vof->C_L[k][j][i], 0.0, 1.0);
                double cs = clampDouble(vof->C_S[k][j][i], 0.0, 1.0);
                double cg = 1.0 - cl - cs;

                if (cg < 0.0) cg = 0.0;
                if (cg > 1.0) cg = 1.0;

                vof->rho[k][j][i] = cl + rho_s * cs + rho_g * cg;
                vof->mu[k][j][i]  = cl + mu_s  * cs + mu_g  * cg;
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->rho, db);
    VOF_DIFFUSE_set_boundary_values(vof->mu, db);
}



#endif
