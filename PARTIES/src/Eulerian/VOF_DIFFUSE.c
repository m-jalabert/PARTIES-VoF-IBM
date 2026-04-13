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

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef VOF_DIFFUSE

static inline double diffuse_clamp01(double value)
{
    return clampDouble(value, 0.0, 1.0);
}

static void diffuse_copy_C_to_F(Cart3d_bag *db)
{
    VolumeFraction *vof = db->vof;
    Array_copy_withghost(vof->C_L, vof->F, db->grid, db->params);
}

static void diffuse_update_phase_cache(Cart3d_bag *db)
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
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

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
    diffuse_update_phase_cache(db);
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
#ifdef LEFT_WALL_VELOCITY_FREESLIP
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
        diffuse_update_phase_cache(db);
        return;
    }

    Particle_list  *lists[2] = { db->lag->p_mobile_list, db->lag->p_fixed_list };

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
                        double arg = (dist - (p->R - shift)) / (denom + 1e-30);
                        double cs = 0.5 - 0.5 * tanh(arg);

                        if (cs > vof->C_S[k][j][i]) {
                            vof->C_S[k][j][i] = cs;
                        }
                    }
                }
            }
            p = p->next;
        }
    }

    diffuse_update_phase_cache(db);
}

void VOF_DIFFUSE_compute_laplacian(double ***in, double ***lap, Cart3d_bag *db)
{
    MAC_grid *grid = db->grid;

    VOF_DIFFUSE_set_boundary_values(in, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                /*
                 * Finite-volume Laplacian: compute face gradients first, then
                 * take the divergence of the diffusive face fluxes using the
                 * same staggered metric factors as the scalar solver in Conc.c.
                 * Applying this operator twice gives the desired FV
                 * discretization of the biharmonic CH term.
                 */
                double dfdxE = (in[k][j][i+1] - in[k][j][i]) * grid->idx_c[i];
                double dfdxW = (in[k][j][i] - in[k][j][i-1]) *
                               ((i != 0) ? grid->idx_c[i-1] : grid->idx_c[i]);

                double dfdyN = (in[k][j+1][i] - in[k][j][i]) * grid->idy_c[j];
                double dfdyS = (in[k][j][i] - in[k][j-1][i]) *
                               ((j != 0) ? grid->idy_c[j-1] : grid->idy_c[j]);

                double dfdzF = (in[k+1][j][i] - in[k][j][i]) * grid->idz_c[k];
                double dfdzB = (in[k][j][i] - in[k-1][j][i]) *
                               ((k != 0) ? grid->idz_c[k-1] : grid->idz_c[k]);

                lap[k][j][i] =
                    (dfdxE - dfdxW) * grid->idx_u[i] +
                    (dfdyN - dfdyS) * grid->idy_v[j] +
                    (dfdzF - dfdzB) * grid->idz_w[k];
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(lap, db);
}

void VOF_DIFFUSE_compute_bulk_S(Cart3d_bag *db)
{
    MAC_grid       *grid = db->grid;
    VolumeFraction *vof  = db->vof;

    diffuse_update_phase_cache(db);

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

                vel = w->data[k][j][i];
                state = diffuse_face_state_z(vof->C_L, i, j, k, vel, params->weno_order);
                flux_z[k][j][i] = vel * state;
            }
        }
    }

    Communication_update_ghost_nodes_flow_variable(flux_x, FLUX_X, params->ghost_nodes, db);
    Communication_update_ghost_nodes_flow_variable(flux_y, FLUX_Y, params->ghost_nodes, db);
    Communication_update_ghost_nodes_flow_variable(flux_z, FLUX_Z, params->ghost_nodes, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double div =
                    (flux_x[k][j][i+1] - flux_x[k][j][i]) / grid->dx_c[i] +
                    (flux_y[k][j+1][i] - flux_y[k][j][i]) / grid->dy_c[j] +
                    (flux_z[k+1][j][i] - flux_z[k][j][i]) / grid->dz_c[k];

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

void VOF_DIFFUSE_step(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;
    const double BET[] = { BETA };
    const double GAMB[] = { GAMBETA };
    const double ZETB[] = { ZETBETA };
    int which_stage = params->which_stage;
    double a_dt = 1.0 / (params->dt * BET[which_stage]);

    /*
     * Keep the current explicit CH RHS in a dedicated workspace so it survives
     * the implicit solve. vof->psi is reused inside the biharmonic solver as a
     * copy of the old state and therefore cannot safely store RHS history.
     */
    double ***rhs_cur = vof->ch_aux2;
    double ***rhs_tmp = vof->ch_aux1;
    double ***rhs_vec = vof->ch_rhs_nm1;

    VOF_DIFFUSE_advect_WENO5(db, rhs_cur);
    VOF_DIFFUSE_explicit_diffusion(db, rhs_tmp);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                rhs_cur[k][j][i] += rhs_tmp[k][j][i];
            }
        }
    }
    VOF_DIFFUSE_set_boundary_values(rhs_cur, db);

    if (which_stage == 0) {
        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    rhs_vec[k][j][i] =
                        a_dt * vof->C_L[k][j][i] +
                        GAMB[0] * rhs_cur[k][j][i];
                }
            }
        }
    } else {
        for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
            for (int j = grid->G_Js; j < grid->G_Je; ++j) {
                for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                    rhs_vec[k][j][i] =
                        a_dt * vof->C_L[k][j][i] +
                        GAMB[which_stage] * rhs_cur[k][j][i] +
                        ZETB[which_stage] * vof->ch_rhs_n[k][j][i];
                }
            }
        }
    }
    VOF_DIFFUSE_set_boundary_values(rhs_vec, db);

    VOF_DIFFUSE_solve_implicit_biharmonic(db, rhs_vec);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                vof->C_L[k][j][i] = diffuse_clamp01(vof->C_L[k][j][i]);
                if (vof->C_L[k][j][i] < 1e-10) {
                    vof->C_L[k][j][i] = 0.0;
                } else if (vof->C_L[k][j][i] > 1.0 - 1e-10) {
                    vof->C_L[k][j][i] = 1.0;
                }
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    diffuse_update_phase_cache(db);

    Array_copy_withghost(vof->ch_rhs_n, vof->ch_rhs_nm1, grid, params);
    Array_copy_withghost(rhs_cur, vof->ch_rhs_n, grid, params);
}

void VOF_DIFFUSE_build_mcl_mask(Cart3d_bag *db)
{
    MAC_grid       *grid = db->grid;
    VolumeFraction *vof  = db->vof;

    for (int k = grid->L_Ks; k < grid->L_Ke; ++k)
    for (int j = grid->L_Js; j < grid->L_Je; ++j)
    for (int i = grid->L_Is; i < grid->L_Ie; ++i) {
        vof->mcl_mask[k][j][i] = 0;
    }

    diffuse_update_phase_cache(db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double cs = vof->C_S[k][j][i];
                double cl = vof->C_L[k][j][i];
                double cg = vof->C_G[k][j][i];

                bool interface_here = (cl >= 0.005 && cl <= 0.995 &&
                                       cg >= 0.005 && cg <= 0.995);
                bool in_ghost_band = (cs >= 0.05 && cs <= 1.0 && interface_here);

                if (!in_ghost_band) {
                    static const int di[6] = { 1, -1, 0, 0, 0, 0 };
                    static const int dj[6] = { 0, 0, 1, -1, 0, 0 };
                    static const int dk[6] = { 0, 0, 0, 0, 1, -1 };

                    for (int dir = 0; dir < 6; ++dir) {
                        int ii = i + di[dir];
                        int jj = j + dj[dir];
                        int kk = k + dk[dir];

                        double csn = vof->C_S[kk][jj][ii];
                        double cln = vof->C_L[kk][jj][ii];
                        double cgn = vof->C_G[kk][jj][ii];

                        if (csn >= 0.0 && csn <= 0.05 &&
                            cln >= 0.005 && cln <= 0.995 &&
                            cgn >= 0.005 && cgn <= 0.995) {
                            in_ghost_band = true;
                            break;
                        }
                    }
                }

                vof->mcl_mask[k][j][i] = (char) in_ghost_band;
            }
        }
    }
}

void VOF_DIFFUSE_apply_contact_angle(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    const double theta = params->contact_angle_deg * PI / 180.0;
    const double cos_t = cos(theta);
    const double sin_t = sin(theta);
    const double cn = params->Cn;
    const double shift = sqrt(2.0) * log(19.0) * cn;
    const double denom = 2.0 * sqrt(2.0) * cn;

    if (db->lag == NULL) {
        return;
    }

    Particle_list  *mobile = db->lag->p_mobile_list;
    Particle_list  *fixed  = db->lag->p_fixed_list;

    VOF_DIFFUSE_build_mcl_mask(db);
    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        double z = grid->zc[k];
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            double y = grid->yc[j];
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                if (!vof->mcl_mask[k][j][i]) {
                    continue;
                }

                double best_d2 = 1e300;
                Particle *best = NULL;
                diffuse_find_nearest_particle(grid->xc[i], y, z, mobile, &best_d2, &best);
                diffuse_find_nearest_particle(grid->xc[i], y, z, fixed,  &best_d2, &best);

                if (best == NULL) {
                    continue;
                }

                double dx = grid->xc[i] - best->X[0];
                double dy = y - best->X[1];
                double dz = z - best->X[2];
                double dist = sqrt(dx * dx + dy * dy + dz * dz);
                if (dist < 1e-14) {
                    continue;
                }

                double n[3] = { dx / dist, dy / dist, dz / dist };
                double x_mid[3] = {
                    best->X[0] + (best->R - shift) * n[0],
                    best->X[1] + (best->R - shift) * n[1],
                    best->X[2] + (best->R - shift) * n[2]
                };

                double grad[3] = {
                    diffuse_grad_x(vof->C_L, grid, i, j, k),
                    diffuse_grad_y(vof->C_L, grid, i, j, k),
                    diffuse_grad_z(vof->C_L, grid, i, j, k)
                };

                double proj = grad[0] * n[0] + grad[1] * n[1] + grad[2] * n[2];
                double t[3] = {
                    grad[0] - proj * n[0],
                    grad[1] - proj * n[1],
                    grad[2] - proj * n[2]
                };
                double tmag = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);

                if (tmag < 1e-12) {
                    double ref[3] = { 1.0, 0.0, 0.0 };
                    if (fabs(n[0]) > 0.9) {
                        ref[0] = 0.0;
                        ref[1] = 1.0;
                    }
                    double dotrn = ref[0] * n[0] + ref[1] * n[1] + ref[2] * n[2];
                    t[0] = ref[0] - dotrn * n[0];
                    t[1] = ref[1] - dotrn * n[1];
                    t[2] = ref[2] - dotrn * n[2];
                    tmag = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
                }

                if (tmag < 1e-12) {
                    continue;
                }

                t[0] /= tmag;
                t[1] /= tmag;
                t[2] /= tmag;

                double ni[3] = {
                    cos_t * n[0] + sin_t * t[0],
                    cos_t * n[1] + sin_t * t[1],
                    cos_t * n[2] + sin_t * t[2]
                };

                double rel[3] = { grid->xc[i] - x_mid[0], y - x_mid[1], z - x_mid[2] };
                double s = rel[0] * ni[0] + rel[1] * ni[1] + rel[2] * ni[2];

                vof->C_L[k][j][i] = diffuse_clamp01(0.5 * (1.0 + tanh(s / (denom + 1e-30))));
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->C_L, db);
    diffuse_update_phase_cache(db);
}

void VOF_DIFFUSE_compute_f_sigma(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    /* Match the code's -2*grad(p) convention used in the momentum RHS. */
    const double scale = 2.0 * 6.0 * sqrt(2.0) / (params->Cn * params->We + 1e-30);

    diffuse_update_phase_cache(db);

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

void VOF_DIFFUSE_apply_f_sigma_old(Cart3d_bag *db)
{
    VOF_apply_f_sigma_old(db);
}

void VOF_DIFFUSE_update_density_viscosity(Cart3d_bag *db)
{
    MAC_grid       *grid   = db->grid;
    Parameters     *params = db->params;
    VolumeFraction *vof    = db->vof;

    Array_copy_withghost(vof->rho, vof->rho_old, grid, params);
    Array_copy_withghost(vof->mu,  vof->mu_old,  grid, params);

    diffuse_update_phase_cache(db);

    for (int k = grid->G_Ks; k < grid->G_Ke; ++k) {
        for (int j = grid->G_Js; j < grid->G_Je; ++j) {
            for (int i = grid->G_Is; i < grid->G_Ie; ++i) {
                double c = vof->C_L[k][j][i];
                vof->rho[k][j][i] = c + (1.0 - c) * (params->rho2 / params->rho1);
                vof->mu[k][j][i]  = c + (1.0 - c) * (params->mu2 / params->mu1);
            }
        }
    }

    VOF_DIFFUSE_set_boundary_values(vof->rho, db);
    VOF_DIFFUSE_set_boundary_values(vof->mu, db);
    VOF_DIFFUSE_set_boundary_values(vof->rho_old, db);
    VOF_DIFFUSE_set_boundary_values(vof->mu_old, db);
}



#endif
