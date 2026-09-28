/*******************************************************************************
 * velocity_cg.c
 *
 * Conjugate-gradient solver for the implicit velocity sub-problem:
 *
 *     A x = b
 *
 * where A = ρ/(α·dt)·I  −  (1/Re)·∇·(μ∇)
 *
 * Changes from original:
 *
 *   1. Jacobi (diagonal) preconditioner
 * 
 *
 *   2. Batched MPI_Allreduce
 *      The two per-iteration reductions (rz_new and rr for convergence) are
 *      computed in a single loop and communicated in one MPI_Allreduce call
 *      with count=2, halving synchronization latency in the inner loop.
 *
 *   3. Fixed abs() → fabs() bug
 *      abs() on a double in C calls the integer overload; for |r| < 1 it
 *      always returns 0.  Replaced with fabs() throughout.
 *
 *   4. Removed dead variables
 *      imax, jmax, kmax, rmax, rr0, yproc, NPY were computed but never used.
 *
 ******************************************************************************/

#include <math.h>   /* fabs, sqrt */
#include "VOF_DIFFUSE.h"   /* DIFFUSE_SOLID_MASS_CUTOFF */

/******************************************************************************/
/*  innerProd — unchanged from original                                       */
/******************************************************************************/
double innerProd(double ***vec1, double ***vec2, char component,
                 MAC_grid *grid, Parameters *params)
{
    int i, j, k;
    double partialSum, totalSum;

    int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;

    int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
    int Ie = min(grid->G_Ie, NX-1);
    int Je = min(grid->G_Je, NY-1);
    int Ke = min(grid->G_Ke, NZ-1);

    if (component == 'u') { Is = max(Is, 1);
#ifdef XPERIODIC
        Ie = grid->G_Ie;
#endif
    }
    if (component == 'v') { Js = max(Js, 1);
#ifdef YPERIODIC
        Je = grid->G_Je;
#endif
    }
    if (component == 'w') { Ks = max(Ks, 1);
#ifdef ZPERIODIC
        Ke = grid->G_Ke;
#endif
    }

    partialSum = 0.0;
    for (k = Ks; k < Ke; k++)
        for (j = Js; j < Je; j++)
            for (i = Is; i < Ie; i++)
                partialSum += vec1[k][j][i] * vec2[k][j][i];

    MPI_Allreduce(&partialSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, PCW);
    return totalSum;
}


#ifdef VOF
#define HARMONIC(mu1, mu2) (2.0 * (mu1) * (mu2) / ((mu1) + (mu2) + 1e-12))
#define AVG2(a, b)         (0.5 * ((a) + (b)))
#endif


/*******************************************************************************
 * vel_jacobi_diag
 *
 * Returns the diagonal entry of the operator A at cell (i,j,k) for the given
 * velocity component.  Mirrors the diagonal contribution of matVec exactly.
 *
 * For VOF (variable coefficients):
 *
 *   diag = ρ_face/(α·dt)
 *        + (μ_xp + μ_xm)/(dx²·Re)
 *        + (μ_yp + μ_ym)/(dy²·Re)
 *        + (μ_zp + μ_zm)/(dz²·Re)
 *
 * In 2D mode the stored z slab is non-physical, so the z contribution is
 * dropped from both the diagonal and the matrix-vector product.
 *
 * For constant-coefficient case:
 *
 *   diag = 1/(α·dt) + 2/Re·(1/dx² + 1/dy² + 1/dz²)   [same for all cells]
 *
 * In the constant case Jacobi is trivially a scalar multiple of identity and
 * reduces to unpreconditioned CG — no harm, no benefit.
 ******************************************************************************/
static inline double vel_jacobi_diag(
    int i, int j, int k, char component,
    const MAC_grid *grid,
    const Parameters *params,
    double ***rho, double ***mu, double ***fliq, double ***fsol,
    double Re, double alpha_k, double dt)
{
    const int collapsed_z = TwodOps_collapsed_component_is_inactive(params);
    const double iRe = 1.0 / Re;
#ifdef VOF
    int iL = (component == 'u') ? i - 1 : i;
    int jL = (component == 'v') ? j - 1 : j;
    int kL = (component == 'w') ? k - 1 : k;

    double rho_face = 0.5 * (rho[k][j][i] + rho[kL][jL][iL]);
    double d_time   = rho_face / (alpha_k * dt);

#ifdef ICE_PENALIZATION
    /* Same fully implicit Darcy diagonal as in matVec. */
    {
        /* Ice fraction is 1 - C_L - C_S: the resolved solid is NOT ice (see the
         * note in matVec).  fsol == NULL reproduces the pre-B.2 mask exactly. */
        double cs_face = (fsol != NULL)
                       ? 0.5 * (fsol[k][j][i] + fsol[kL][jL][iL]) : 0.0;
        double phi_face = 1.0 - 0.5 * (fliq[k][j][i] + fliq[kL][jL][iL]) - cs_face;
        if (cs_face >= DIFFUSE_SOLID_MASS_CUTOFF) phi_face = 0.0;
        else if (cs_face > DIFFUSE_SOLID_HALO_CUTOFF &&
                 phi_face < DIFFUSE_SOLID_HALO_ICE_TRUST) phi_face = 0.0;
        if (phi_face < 0.0) phi_face = 0.0;
        if (phi_face > 1.0) phi_face = 1.0;
#ifdef VOF_DIFFUSE_ICE_PENAL_THRESHOLD
        /* Yang-style hard mask: rigid where genuinely ice, free elsewhere. */
        phi_face = (phi_face > DIFFUSE_ICE_PENAL_THRESHOLD_VAL) ? 1.0 : 0.0;
#endif
        d_time += rho_face * 2.0 * phi_face / params->darcy_tau;
    }
#else
    (void)fliq;
    (void)fsol;
#endif

    double diag_x = 0.0;
    double diag_y = 0.0;
    double diag_z = 0.0;
    double diag_axisym = 0.0;
    double mu_xp, mu_xm, mu_yp, mu_ym;
    double mu_zp = 0.0, mu_zm = 0.0;

    if (component == 'u') {
        mu_xp = mu[k][j][i];
        mu_xm = mu[k][j][i-1];

        mu_yp = HARMONIC(AVG2(mu[k][j+1][i-1], mu[k][j+1][i]),
                         AVG2(mu[k][j  ][i-1], mu[k][j  ][i]));
        mu_ym = HARMONIC(AVG2(mu[k][j  ][i-1], mu[k][j  ][i]),
                         AVG2(mu[k][j-1][i-1], mu[k][j-1][i]));

        if (!collapsed_z) {
            mu_zp = HARMONIC(AVG2(mu[k+1][j][i-1], mu[k+1][j][i]),
                             AVG2(mu[k  ][j][i-1], mu[k  ][j][i]));
            mu_zm = HARMONIC(AVG2(mu[k  ][j][i-1], mu[k  ][j][i]),
                             AVG2(mu[k-1][j][i-1], mu[k-1][j][i]));
        }

        diag_x = iRe * TwodOps_u_cv_x_diag_coeff(
            grid, params, i,
            mu_xp * grid->idx_u[i],
            mu_xm * grid->idx_u[i-1]);
        diag_y = iRe * (mu_yp * grid->idy_c[j] +
                        mu_ym * grid->idy_c[j-1]) * grid->idy_v[j];
        if (!collapsed_z) {
            diag_z = iRe * (mu_zp * grid->idz_c[k] +
                            mu_zm * grid->idz_c[k-1]) * grid->idz_w[k];
        }
        diag_axisym = TwodOps_axisym_u_radial_linear_coeff(
            grid, params, i,
            iRe * 0.5 * (mu_xp + mu_xm));

    } else if (component == 'v') {
        mu_yp = mu[k][j  ][i];
        mu_ym = mu[k][j-1][i];

        mu_xp = HARMONIC(AVG2(mu[k][j-1][i+1], mu[k][j  ][i+1]),
                         AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]));
        mu_xm = HARMONIC(AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]),
                         AVG2(mu[k][j-1][i-1], mu[k][j  ][i-1]));

        if (!collapsed_z) {
            mu_zp = HARMONIC(AVG2(mu[k+1][j-1][i], mu[k+1][j  ][i]),
                             AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]));
            mu_zm = HARMONIC(AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]),
                             AVG2(mu[k-1][j-1][i], mu[k-1][j  ][i]));
        }

        diag_x = iRe * TwodOps_v_cv_x_diag_coeff(
            grid, params, i,
            mu_xp * grid->idx_c[i],
            mu_xm * grid->idx_c[i-1]);
        diag_y = iRe * (mu_yp * grid->idy_v[j] +
                        mu_ym * grid->idy_v[j-1]) * grid->idy_c[j-1];
        if (!collapsed_z) {
            diag_z = iRe * (mu_zp * grid->idz_c[k] +
                            mu_zm * grid->idz_c[k-1]) * grid->idz_w[k];
        }

    } else { /* 'w' */
        if (!collapsed_z) {
            mu_zp = mu[k  ][j][i];
            mu_zm = mu[k-1][j][i];
        }

        mu_xp = HARMONIC(AVG2(mu[k-1][j][i+1], mu[k  ][j][i+1]),
                         AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]));
        mu_xm = HARMONIC(AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]),
                         AVG2(mu[k-1][j][i-1], mu[k  ][j][i-1]));

        mu_yp = HARMONIC(AVG2(mu[k-1][j+1][i], mu[k  ][j+1][i]),
                         AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]));
        mu_ym = HARMONIC(AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]),
                         AVG2(mu[k-1][j-1][i], mu[k  ][j-1][i]));

        diag_x = iRe * (mu_xp * grid->idx_c[i] +
                        mu_xm * grid->idx_c[i-1]) * grid->idx_u[i];
        diag_y = iRe * (mu_yp * grid->idy_c[j] +
                        mu_ym * grid->idy_c[j-1]) * grid->idy_v[j];
        if (!collapsed_z) {
            diag_z = iRe * (mu_zp * grid->idz_w[k] +
                            mu_zm * grid->idz_w[k-1]) * grid->idz_c[k-1];
        }
    }

    return d_time + diag_x + diag_y + diag_z + diag_axisym;

#else
    double d_time = 1.0 / (alpha_k * dt);
    double diag_x = 0.0;
    double diag_y = 0.0;
    double diag_z = 0.0;
    double diag_axisym = 0.0;

    (void)rho;
    (void)mu;
    (void)fliq;

    if (component == 'u') {
        diag_x = iRe * TwodOps_u_cv_x_diag_coeff(
            grid, params, i,
            grid->idx_u[i],
            grid->idx_u[i-1]);
        diag_y = iRe * (grid->idy_c[j] + grid->idy_c[j-1]) * grid->idy_v[j];
        if (!collapsed_z)
            diag_z = iRe * (grid->idz_c[k] + grid->idz_c[k-1]) * grid->idz_w[k];
        diag_axisym = TwodOps_axisym_u_radial_linear_coeff(
            grid, params, i, iRe);
    } else if (component == 'v') {
        diag_x = iRe * TwodOps_v_cv_x_diag_coeff(
            grid, params, i,
            grid->idx_c[i],
            grid->idx_c[i-1]);
        diag_y = iRe * (grid->idy_v[j] + grid->idy_v[j-1]) * grid->idy_c[j-1];
        if (!collapsed_z)
            diag_z = iRe * (grid->idz_c[k] + grid->idz_c[k-1]) * grid->idz_w[k];
    } else {
        diag_x = iRe * (grid->idx_c[i] + grid->idx_c[i-1]) * grid->idx_u[i];
        diag_y = iRe * (grid->idy_c[j] + grid->idy_c[j-1]) * grid->idy_v[j];
        if (!collapsed_z)
            diag_z = iRe * (grid->idz_w[k] + grid->idz_w[k-1]) * grid->idz_c[k-1];
    }

    return d_time + diag_x + diag_y + diag_z + diag_axisym;
#endif
}


/******************************************************************************/
/*  matVec — 2D-aware implicit velocity operator                              */
/******************************************************************************/
void matVec(double ***Ax, double ***x, char component, Cart3d_bag *data_bag)
{
    int i, j, k;
    const int collapsed_z = TwodOps_collapsed_component_is_inactive(data_bag->params);

    MAC_grid   *grid   = data_bag->grid;
    Parameters *params = data_bag->params;
    int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;

    int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
    int Ie = min(grid->G_Ie, NX-1);
    int Je = min(grid->G_Je, NY-1);
    int Ke = min(grid->G_Ke, NZ-1);

    if (component == 'u') { Is = max(Is, 1);
#ifdef XPERIODIC
        Ie = grid->G_Ie;
#endif
    }
    if (component == 'v') { Js = max(Js, 1);
#ifdef YPERIODIC
        Je = grid->G_Je;
#endif
    }
    if (component == 'w') { Ks = max(Ks, 1);
#ifdef ZPERIODIC
        Ke = grid->G_Ke;
#endif
    }

    double Re = params->Re;
    double iRe = 1.0 / Re;
    const double BET[] = { BETA };
    int stage     = params->which_stage;
    double alpha_k = BET[stage];
    double dt      = params->dt;

    double *idx_u = grid->idx_u;
    double *idy_v = grid->idy_v;
    double *idz_w = grid->idz_w;
    double *idx_c = grid->idx_c;
    double *idy_c = grid->idy_c;
    double *idz_c = grid->idz_c;

#ifdef VOF
    double ***rho = data_bag->vof->rho;
    double ***mu  = data_bag->vof->mu;

#ifdef ICE_PENALIZATION
    /*
     * Darcy/Brinkman ice damping: -(phi_s/darcy_tau) u enters the stage solve
     * fully implicitly as a positive local diagonal (factor 2 from the code's
     * Crank-Nicolson convention, cf. the 2x on explicit sources).  The
     * operator stays SPD, so CG and the Jacobi preconditioner are unchanged.
     */
    #ifdef VOF_DIFFUSE
    double ***fliq = data_bag->vof->C_L;
    double ***fsol = data_bag->vof->C_S;
    #else
    double ***fliq = data_bag->vof->F;
    double ***fsol = NULL;
    #endif
    const double drag_2tau = 2.0 / params->darcy_tau;
#endif

    for (k = Ks; k < Ke; k++) {
        for (j = Js; j < Je; j++) {
            for (i = Is; i < Ie; i++) {

                int iL = (component == 'u') ? i-1 : i;
                int jL = (component == 'v') ? j-1 : j;
                int kL = (component == 'w') ? k-1 : k;

                double rho_face    = 0.5*(rho[k][j][i] + rho[kL][jL][iL]);
                double factor_time = rho_face / (alpha_k * dt);

#ifdef ICE_PENALIZATION
                /* Roadmap B.2 -- the Brinkman ice mask must exclude the RESOLVED solid.
                 * phi_s = 1 - C_L counts a sediment grain (where C_L = 0) as ice and
                 * Darcy-damps it at tau = 1e-4, on top of the IBM forcing that already
                 * represents it -- the grain is pinned and creeps instead of settling
                 * (measured U_y = 5.7e-3 against an O(1) terminal velocity at Ga = 122,
                 * job 20253930).  The ice fraction is 1 - C_L - C_S.
                 * Identical in Stage A, where C_S == 0. */
                /* No Brinkman damping INSIDE the resolved solid: the IBM already
                 * governs there, and the ternary deficit (C_L + C_S < 1 in the
                 * grain's diffuse rim) would otherwise be charged to the ice
                 * phase.  Measured in an ICE-FREE settling deck: phi_s = 0.77
                 * over 0.5R-R, giving a Darcy rate of 15,467 against an inertial
                 * scale of 227 -- the particle was clamped 200-500x too slow
                 * (jobs 20269208/20269209).  The round-4 CH mask makes this worse
                 * by construction, since it pins C_L = 0 wherever C_S >= cutoff. */
                double cs_face = (fsol != NULL)
                               ? 0.5*(fsol[k][j][i] + fsol[kL][jL][iL]) : 0.0;
                double phi_face = 1.0 - 0.5*(fliq[k][j][i] + fliq[kL][jL][iL]) - cs_face;
                if (cs_face >= DIFFUSE_SOLID_MASS_CUTOFF) phi_face = 0.0;
                else if (cs_face > DIFFUSE_SOLID_HALO_CUTOFF &&
                         phi_face < DIFFUSE_SOLID_HALO_ICE_TRUST) phi_face = 0.0;
                if (phi_face < 0.0) phi_face = 0.0;
                if (phi_face > 1.0) phi_face = 1.0;
#ifdef VOF_DIFFUSE_ICE_PENAL_THRESHOLD
                /* Yang-style hard mask: rigid where genuinely ice, free elsewhere. */
                phi_face = (phi_face > DIFFUSE_ICE_PENAL_THRESHOLD_VAL) ? 1.0 : 0.0;
#endif
                factor_time += rho_face * drag_2tau * phi_face;
#endif

                double Ax_val      = factor_time * x[k][j][i];
                double op_val      = 0.0;

                double mu_xp, mu_xm, mu_yp, mu_ym;
                double mu_zp = 0.0, mu_zm = 0.0;

                if (component == 'u') {
                    mu_xp = mu[k][j][i];
                    mu_xm = mu[k][j][i-1];

                    mu_yp = HARMONIC(AVG2(mu[k][j+1][i-1], mu[k][j+1][i  ]),
                                     AVG2(mu[k][j  ][i-1], mu[k][j  ][i  ]));
                    mu_ym = HARMONIC(AVG2(mu[k][j  ][i-1], mu[k][j  ][i  ]),
                                     AVG2(mu[k][j-1][i-1], mu[k][j-1][i  ]));

                    if (!collapsed_z) {
                        mu_zp = HARMONIC(AVG2(mu[k+1][j][i-1], mu[k+1][j][i  ]),
                                         AVG2(mu[k  ][j][i-1], mu[k  ][j][i  ]));
                        mu_zm = HARMONIC(AVG2(mu[k  ][j][i-1], mu[k  ][j][i  ]),
                                         AVG2(mu[k-1][j][i-1], mu[k-1][j][i  ]));
                    }

                    {
                        double flux_xp = mu_xp * (x[k][j][i+1] - x[k][j][i]) * idx_u[i];
                        double flux_xm = mu_xm * (x[k][j][i]   - x[k][j][i-1]) * idx_u[i-1];
                        double flux_yp = mu_yp * (x[k][j+1][i] - x[k][j][i]) * idy_c[j];
                        double flux_ym = mu_ym * (x[k][j][i]   - x[k][j-1][i]) * idy_c[j-1];

                        op_val += iRe * TwodOps_u_cv_x_flux_divergence(
                            grid, params, i, flux_xp, flux_xm);
                        op_val += iRe * (flux_yp - flux_ym) * idy_v[j];

                        if (!collapsed_z) {
                            double flux_zp = mu_zp * (x[k+1][j][i] - x[k][j][i]) * idz_c[k];
                            double flux_zm = mu_zm * (x[k][j][i]   - x[k-1][j][i]) * idz_c[k-1];
                            op_val += iRe * (flux_zp - flux_zm) * idz_w[k];
                        }

                        op_val += TwodOps_axisym_u_radial_linear_term(
                            grid, params, i,
                            iRe * 0.5 * (mu_xp + mu_xm),
                            x[k][j][i]);
                    }

                } else if (component == 'v') {
                    mu_yp = mu[k][j  ][i];
                    mu_ym = mu[k][j-1][i];

                    mu_xp = HARMONIC(AVG2(mu[k][j-1][i+1], mu[k][j  ][i+1]),
                                     AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]));
                    mu_xm = HARMONIC(AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]),
                                     AVG2(mu[k][j-1][i-1], mu[k][j  ][i-1]));

                    if (!collapsed_z) {
                        mu_zp = HARMONIC(AVG2(mu[k+1][j-1][i], mu[k+1][j  ][i]),
                                         AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]));
                        mu_zm = HARMONIC(AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]),
                                         AVG2(mu[k-1][j-1][i], mu[k-1][j  ][i]));
                    }

                    {
                        double flux_xp = mu_xp * (x[k][j][i+1] - x[k][j][i]) * idx_c[i];
                        double flux_xm = mu_xm * (x[k][j][i]   - x[k][j][i-1]) * idx_c[i-1];
                        double flux_yp = mu_yp * (x[k][j+1][i] - x[k][j][i]) * idy_v[j];
                        double flux_ym = mu_ym * (x[k][j][i]   - x[k][j-1][i]) * idy_v[j-1];

                        op_val += iRe * TwodOps_v_cv_x_flux_divergence(
                            grid, params, i, flux_xp, flux_xm);
                        op_val += iRe * (flux_yp - flux_ym) * idy_c[j-1];

                        if (!collapsed_z) {
                            double flux_zp = mu_zp * (x[k+1][j][i] - x[k][j][i]) * idz_c[k];
                            double flux_zm = mu_zm * (x[k][j][i]   - x[k-1][j][i]) * idz_c[k-1];
                            op_val += iRe * (flux_zp - flux_zm) * idz_w[k];
                        }
                    }

                } else {
                    if (!collapsed_z) {
                        mu_zp = mu[k  ][j][i];
                        mu_zm = mu[k-1][j][i];
                    }

                    mu_xp = HARMONIC(AVG2(mu[k-1][j][i+1], mu[k  ][j][i+1]),
                                     AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]));
                    mu_xm = HARMONIC(AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]),
                                     AVG2(mu[k-1][j][i-1], mu[k  ][j][i-1]));

                    mu_yp = HARMONIC(AVG2(mu[k-1][j+1][i], mu[k  ][j+1][i]),
                                     AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]));
                    mu_ym = HARMONIC(AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]),
                                     AVG2(mu[k-1][j-1][i], mu[k  ][j-1][i]));

                    {
                        double flux_xp = mu_xp * (x[k][j][i+1] - x[k][j][i]) * idx_c[i];
                        double flux_xm = mu_xm * (x[k][j][i]   - x[k][j][i-1]) * idx_c[i-1];
                        double flux_yp = mu_yp * (x[k][j+1][i] - x[k][j][i]) * idy_c[j];
                        double flux_ym = mu_ym * (x[k][j][i]   - x[k][j-1][i]) * idy_c[j-1];

                        op_val += iRe * (flux_xp - flux_xm) * idx_u[i];
                        op_val += iRe * (flux_yp - flux_ym) * idy_v[j];

                        if (!collapsed_z) {
                            double flux_zp = mu_zp * (x[k+1][j][i] - x[k][j][i]) * idz_w[k];
                            double flux_zm = mu_zm * (x[k][j][i]   - x[k-1][j][i]) * idz_w[k-1];
                            op_val += iRe * (flux_zp - flux_zm) * idz_c[k-1];
                        }
                    }
                }

                Ax[k][j][i] = Ax_val - op_val;
            }
        }
    }

#else
    double idtimeb = 1.0 / (BET[params->which_stage] * dt);

    for (k = Ks; k < Ke; k++) {
        for (j = Js; j < Je; j++) {
            for (i = Is; i < Ie; i++) {
                double Ax_val = idtimeb * x[k][j][i];
                double op_val = 0.0;

                if (component == 'u') {
                    double flux_xp = (x[k][j][i+1] - x[k][j][i]) * idx_u[i];
                    double flux_xm = (x[k][j][i]   - x[k][j][i-1]) * idx_u[i-1];
                    double flux_yp = (x[k][j+1][i] - x[k][j][i]) * idy_c[j];
                    double flux_ym = (x[k][j][i]   - x[k][j-1][i]) * idy_c[j-1];

                    op_val += iRe * TwodOps_u_cv_x_flux_divergence(
                        grid, params, i, flux_xp, flux_xm);
                    op_val += iRe * (flux_yp - flux_ym) * idy_v[j];

                    if (!collapsed_z) {
                        double flux_zp = (x[k+1][j][i] - x[k][j][i]) * idz_c[k];
                        double flux_zm = (x[k][j][i]   - x[k-1][j][i]) * idz_c[k-1];
                        op_val += iRe * (flux_zp - flux_zm) * idz_w[k];
                    }

                    op_val += TwodOps_axisym_u_radial_linear_term(
                        grid, params, i, iRe, x[k][j][i]);

                } else if (component == 'v') {
                    double flux_xp = (x[k][j][i+1] - x[k][j][i]) * idx_c[i];
                    double flux_xm = (x[k][j][i]   - x[k][j][i-1]) * idx_c[i-1];
                    double flux_yp = (x[k][j+1][i] - x[k][j][i]) * idy_v[j];
                    double flux_ym = (x[k][j][i]   - x[k][j-1][i]) * idy_v[j-1];

                    op_val += iRe * TwodOps_v_cv_x_flux_divergence(
                        grid, params, i, flux_xp, flux_xm);
                    op_val += iRe * (flux_yp - flux_ym) * idy_c[j-1];

                    if (!collapsed_z) {
                        double flux_zp = (x[k+1][j][i] - x[k][j][i]) * idz_c[k];
                        double flux_zm = (x[k][j][i]   - x[k-1][j][i]) * idz_c[k-1];
                        op_val += iRe * (flux_zp - flux_zm) * idz_w[k];
                    }

                } else {
                    double flux_xp = (x[k][j][i+1] - x[k][j][i]) * idx_c[i];
                    double flux_xm = (x[k][j][i]   - x[k][j][i-1]) * idx_c[i-1];
                    double flux_yp = (x[k][j+1][i] - x[k][j][i]) * idy_c[j];
                    double flux_ym = (x[k][j][i]   - x[k][j-1][i]) * idy_c[j-1];

                    op_val += iRe * (flux_xp - flux_xm) * idx_u[i];
                    op_val += iRe * (flux_yp - flux_ym) * idy_v[j];

                    if (!collapsed_z) {
                        double flux_zp = (x[k+1][j][i] - x[k][j][i]) * idz_w[k];
                        double flux_zm = (x[k][j][i]   - x[k-1][j][i]) * idz_w[k-1];
                        op_val += iRe * (flux_zp - flux_zm) * idz_c[k-1];
                    }
                }

                Ax[k][j][i] = Ax_val - op_val;
            }
        }
    }
#endif
}


/******************************************************************************/
/*
 * Velocity_solve_cg
 *
 * Jacobi-preconditioned conjugate-gradient solver for  A x = b.
 *
 * Algorithm (Jacobi PCG):
 *
 *   r  = b - A*x0
 *   d  = M^{-1} r          (M = diag(A), applied cell-by-cell)
 *   rz = <r, M^{-1}r>      (one Allreduce)
 *
 *   loop:
 *     Ad  = A*d
 *     dAd = <d, Ad>         (Allreduce 1)
 *     α   = rz / dAd
 *     x  += α * d
 *     r  -= α * Ad
 *     rz_new = <r, M^{-1}r> \
 *     rr     = <r, r>        |--- single count=2 Allreduce (Allreduce 2)
 *     RMS = sqrt(rr / N_dof)
 *     if RMS < tol: converged
 *     β  = rz_new / rz
 *     d  = M^{-1}r + β*d
 *     rz = rz_new
 *
 * Convergence is checked on the unpreconditioned residual RMS (same metric
 * as before, printed in the log) for comparability with previous results.
 *****************************************************************************/
int Velocity_solve_cg(Velocity *vel, Cart3d_bag *data_bag)
{
    int iters, i, j, k;

    MAC_grid   *grid   = data_bag->grid;
    Parameters *params = data_bag->params;
    char statement[100];

    const int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;

    double ***data = vel->data;
    double ***rhs  = vel->ng_rhs;
    char component = vel->component;

    const double EPS = 1e-14;

    if (TwodOps_collapsed_component_is_inactive(params) && component == 'w') {
        Memory_reset_flow_variable(grid, params, data);
        Memory_reset_noghost_variable(grid, params, rhs);
        Memory_reset_flow_variable(grid, params, vel->d);
        Memory_reset_noghost_variable(grid, params, vel->ng_r);
        Memory_reset_noghost_variable(grid, params, vel->ng_Ad);
#ifdef VOF
        Memory_reset_flow_variable(grid, params, vel->M_inv);
#endif
        Velocity_update_boundaries(data, component, VEL_TYPE_NORMAL, data_bag);
        sprintf(statement, "w-component of velocity converged to 0 after 0 iterations\n");
        Display_progress(params, statement);
        return 0;
    }

    /* ------------------------------------------------------------------ *
     * Index ranges (identical to original)
     * ------------------------------------------------------------------ */
    int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
    int Ie = min(grid->G_Ie, NX-1);
    int Je = min(grid->G_Je, NY-1);
    int Ke = min(grid->G_Ke, NZ-1);

    if (component == 'u') { Is = max(Is, 1);
#ifdef XPERIODIC
        Ie = grid->G_Ie;
#endif
    }
    if (component == 'v') { Js = max(Js, 1);
#ifdef YPERIODIC
        Je = grid->G_Je;
#endif
    }
    if (component == 'w') { Ks = max(Ks, 1);
#ifdef ZPERIODIC
        Ke = grid->G_Ke;
#endif
    }

    const double Re = params->Re;
    const double BET[] = { BETA };
    const double alpha_k = BET[params->which_stage];
    const double dt      = params->dt;

#ifdef VOF
    double ***rho = data_bag->vof->rho;
    double ***mu  = data_bag->vof->mu;
    #ifdef VOF_DIFFUSE
    double ***fliq = data_bag->vof->C_L;   /* liquid fraction for ICE_PENALIZATION */
    double ***fsol = data_bag->vof->C_S;   /* resolved solid: NOT ice (B.2 fix)   */
    #else
    double ***fliq = data_bag->vof->F;
    double ***fsol = NULL;                 /* no resolved-solid field without VOF_DIFFUSE */
    #endif
#else
    double ***rho = NULL;  /* unused in constant-coefficient path */
    double ***mu  = NULL;
    double ***fliq = NULL;
    double ***fsol = NULL;
#endif

    /* ------------------------------------------------------------------ *
     * DOF normaliser for RMS residual (unchanged from original)
     * ------------------------------------------------------------------ */
    double in_tot;
    if (component == 'u') {
#ifdef XPERIODIC
        in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-1));
#else
        in_tot = 1.0 / (double)((NX-2)*(NY-1)*(NZ-1));
#endif
    } else if (component == 'v') {
#ifdef YPERIODIC
        in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-1));
#else
        in_tot = 1.0 / (double)((NX-1)*(NY-2)*(NZ-1));
#endif
    } else {
#ifdef ZPERIODIC
        in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-1));
#else
        in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-2));
#endif
    }

    /* ------------------------------------------------------------------ *
     * CG work arrays (unchanged from original)
     * ------------------------------------------------------------------ */
    double ***r  = vel->ng_r;
    double ***d  = vel->d;
    double ***Ad = vel->ng_Ad;

    /* ------------------------------------------------------------------ *
     * Initialise: r = rhs - A*x,   d = M^{-1} r
     * ------------------------------------------------------------------ */
    matVec(Ad, data, component, data_bag);

    for (k = Ks; k < Ke; k++) {
        for (j = Js; j < Je; j++) {
            for (i = Is; i < Ie; i++) {
                r[k][j][i] = rhs[k][j][i] - Ad[k][j][i];

                /* d = M^{-1} r: divide initial residual by diagonal */
                double diag = vel_jacobi_diag(i, j, k, component, grid, params,
                                              rho, mu, fliq, fsol,
                                              Re, alpha_k, dt);
                d[k][j][i] = r[k][j][i] / (diag + EPS);
            }
        }
    }

    /* Initial rz = <r, M^{-1}r> = <r, d> */
    double rz = innerProd(r, d, component, grid, params);

    /* ------------------------------------------------------------------ *
     * Jacobi-PCG main loop
     * ------------------------------------------------------------------ */
    double RMS = 0.0;
    iters = 0;

    while (iters < params->CG_MAXIT) {

        iters++;

        /* Update ghost/boundary cells for d */
        Velocity_update_boundaries(d, component, VEL_TYPE_CG, data_bag);

        /* Ad = A * d */
        matVec(Ad, d, component, data_bag);

        /* dAd = <d, Ad>  —  Allreduce 1 */
        double dAd = innerProd(d, Ad, component, grid, params);

        double alpha = rz / (dAd + EPS);

        /* Update solution and residual; accumulate local partial sums for
         * rz_new = <r_new, M^{-1}r_new> and rr = <r_new, r_new> in one pass
         * so both can be reduced in a single MPI call (Allreduce 2).        */
        double local_vals[2] = {0.0, 0.0};  /* [rz_new_partial, rr_partial] */

        for (k = Ks; k < Ke; k++) {
            for (j = Js; j < Je; j++) {
                for (i = Is; i < Ie; i++) {

                    data[k][j][i] += alpha * d[k][j][i];
                    r[k][j][i]    -= alpha * Ad[k][j][i];

                    double ri   = r[k][j][i];
                    double diag = vel_jacobi_diag(i, j, k, component, grid, params,
                                                  rho, mu, fliq, fsol,
                                                  Re, alpha_k, dt);
                    local_vals[0] += ri * ri / (diag + EPS);  /* rz_new */
                    local_vals[1] += ri * ri;                  /* rr     */
                }
            }
        }

        /* Single Allreduce for both rz_new and rr — Allreduce 2 */
        double global_vals[2];
        MPI_Allreduce(local_vals, global_vals, 2, MPI_DOUBLE, MPI_SUM, PCW);

        double rz_new = global_vals[0];
        double rr     = global_vals[1];

        /* Convergence check on unpreconditioned residual RMS */
        RMS = sqrt(in_tot * rr);
        if (RMS < params->CG_ETOL) break;

        double beta = rz_new / (rz + EPS);
        rz = rz_new;

        /* d = M^{-1} r + beta * d */
        for (k = Ks; k < Ke; k++) {
            for (j = Js; j < Je; j++) {
                for (i = Is; i < Ie; i++) {
                    double diag = vel_jacobi_diag(i, j, k, component, grid, params,
                                                  rho, mu, fliq, fsol,
                                                  Re, alpha_k, dt);
                    d[k][j][i] = r[k][j][i] / (diag + EPS) + beta * d[k][j][i];
                }
            }
        }
    }

    sprintf(statement,
            "%c-component of velocity converged to %g after %d iterations\n",
            component, RMS, iters);
    Display_progress(params, statement);

    return iters;
}
