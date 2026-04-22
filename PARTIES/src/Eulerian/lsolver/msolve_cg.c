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
 * For constant-coefficient case:
 *
 *   diag = 1/(α·dt) + 2/Re·(1/dx² + 1/dy² + 1/dz²)   [same for all cells]
 *
 * In the constant case Jacobi is trivially a scalar multiple of identity and
 * reduces to unpreconditioned CG — no harm, no benefit.
 ******************************************************************************/
static inline double vel_jacobi_diag(
    int i, int j, int k, char component,
    double ***rho, double ***mu,
    double inv_dx2, double inv_dy2, double inv_dz2,
    double Re, double alpha_k, double dt)
{
#ifdef VOF
    int iL = (component == 'u') ? i - 1 : i;
    int jL = (component == 'v') ? j - 1 : j;
    int kL = (component == 'w') ? k - 1 : k;

    double rho_face = 0.5 * (rho[k][j][i] + rho[kL][jL][iL]);
    double d_time   = rho_face / (alpha_k * dt);

    double mu_xp, mu_xm, mu_yp, mu_ym, mu_zp, mu_zm;

    if (component == 'u') {
        mu_xp = mu[k][j][i];
        mu_xm = mu[k][j][i-1];

        mu_yp = HARMONIC(AVG2(mu[k][j+1][i-1], mu[k][j+1][i]),
                         AVG2(mu[k][j  ][i-1], mu[k][j  ][i]));
        mu_ym = HARMONIC(AVG2(mu[k][j  ][i-1], mu[k][j  ][i]),
                         AVG2(mu[k][j-1][i-1], mu[k][j-1][i]));

        mu_zp = HARMONIC(AVG2(mu[k+1][j][i-1], mu[k+1][j][i]),
                         AVG2(mu[k  ][j][i-1], mu[k  ][j][i]));
        mu_zm = HARMONIC(AVG2(mu[k  ][j][i-1], mu[k  ][j][i]),
                         AVG2(mu[k-1][j][i-1], mu[k-1][j][i]));

    } else if (component == 'v') {
        mu_yp = mu[k][j  ][i];
        mu_ym = mu[k][j-1][i];

        mu_xp = HARMONIC(AVG2(mu[k][j-1][i+1], mu[k][j  ][i+1]),
                         AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]));
        mu_xm = HARMONIC(AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]),
                         AVG2(mu[k][j-1][i-1], mu[k][j  ][i-1]));

        mu_zp = HARMONIC(AVG2(mu[k+1][j-1][i], mu[k+1][j  ][i]),
                         AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]));
        mu_zm = HARMONIC(AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]),
                         AVG2(mu[k-1][j-1][i], mu[k-1][j  ][i]));

    } else { /* 'w' */
        mu_zp = mu[k  ][j][i];
        mu_zm = mu[k-1][j][i];

        mu_xp = HARMONIC(AVG2(mu[k-1][j][i+1], mu[k  ][j][i+1]),
                         AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]));
        mu_xm = HARMONIC(AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]),
                         AVG2(mu[k-1][j][i-1], mu[k  ][j][i-1]));

        mu_yp = HARMONIC(AVG2(mu[k-1][j+1][i], mu[k  ][j+1][i]),
                         AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]));
        mu_ym = HARMONIC(AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]),
                         AVG2(mu[k-1][j-1][i], mu[k  ][j-1][i]));
    }

    return d_time
         + (mu_xp + mu_xm) * inv_dx2 / Re
         + (mu_yp + mu_ym) * inv_dy2 / Re
         + (mu_zp + mu_zm) * inv_dz2 / Re;

#else
    /* Constant-coefficient: diagonal is uniform, Jacobi = scalar shift */
    (void)i; (void)j; (void)k; (void)component; (void)rho; (void)mu;
    return 1.0 / (alpha_k * dt) + 2.0 / Re * (inv_dx2 + inv_dy2 + inv_dz2);
#endif
}


/******************************************************************************/
/*  matVec — unchanged from original                                          */
/******************************************************************************/
void matVec(double ***Ax, double ***x, char component, Cart3d_bag *data_bag)
{
    int i, j, k;

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

    double dx = 1.0 / grid->idx_c[1];
    double dy = 1.0 / grid->idy_c[1];
    double dz = 1.0 / grid->idz_c[1];

    double Re = params->Re;
    const double BET[] = { BETA };
    int stage     = params->which_stage;
    double alpha_k = BET[stage];
    double dt      = params->dt;

#ifdef VOF
    double ***rho = data_bag->vof->rho;
    double ***mu  = data_bag->vof->mu;

    for (k = Ks; k < Ke; k++) {
        for (j = Js; j < Je; j++) {
            for (i = Is; i < Ie; i++) {

                int iL = (component == 'u') ? i-1 : i;
                int jL = (component == 'v') ? j-1 : j;
                int kL = (component == 'w') ? k-1 : k;

                double rho_face    = 0.5*(rho[k][j][i] + rho[kL][jL][iL]);
                double factor_time = rho_face / (alpha_k * dt);
                double Ax_val      = factor_time * x[k][j][i];

                double mu_xp, mu_xm, mu_yp, mu_ym, mu_zp, mu_zm;

                if (component == 'u') {
                    mu_xp = mu[k][j][i];
                    mu_xm = mu[k][j][i-1];

                    mu_yp = HARMONIC(AVG2(mu[k][j+1][i-1], mu[k][j+1][i  ]),
                                     AVG2(mu[k][j  ][i-1], mu[k][j  ][i  ]));
                    mu_ym = HARMONIC(AVG2(mu[k][j  ][i-1], mu[k][j  ][i  ]),
                                     AVG2(mu[k][j-1][i-1], mu[k][j-1][i  ]));

                    mu_zp = HARMONIC(AVG2(mu[k+1][j][i-1], mu[k+1][j][i  ]),
                                     AVG2(mu[k  ][j][i-1], mu[k  ][j][i  ]));
                    mu_zm = HARMONIC(AVG2(mu[k  ][j][i-1], mu[k  ][j][i  ]),
                                     AVG2(mu[k-1][j][i-1], mu[k-1][j][i  ]));

                } else if (component == 'v') {
                    mu_yp = mu[k][j  ][i];
                    mu_ym = mu[k][j-1][i];

                    mu_xp = HARMONIC(AVG2(mu[k][j-1][i+1], mu[k][j  ][i+1]),
                                     AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]));
                    mu_xm = HARMONIC(AVG2(mu[k][j-1][i  ], mu[k][j  ][i  ]),
                                     AVG2(mu[k][j-1][i-1], mu[k][j  ][i-1]));

                    mu_zp = HARMONIC(AVG2(mu[k+1][j-1][i], mu[k+1][j  ][i]),
                                     AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]));
                    mu_zm = HARMONIC(AVG2(mu[k  ][j-1][i], mu[k  ][j  ][i]),
                                     AVG2(mu[k-1][j-1][i], mu[k-1][j  ][i]));

                } else {
                    mu_zp = mu[k  ][j][i];
                    mu_zm = mu[k-1][j][i];

                    mu_xp = HARMONIC(AVG2(mu[k-1][j][i+1], mu[k  ][j][i+1]),
                                     AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]));
                    mu_xm = HARMONIC(AVG2(mu[k-1][j][i  ], mu[k  ][j][i  ]),
                                     AVG2(mu[k-1][j][i-1], mu[k  ][j][i-1]));

                    mu_yp = HARMONIC(AVG2(mu[k-1][j+1][i], mu[k  ][j+1][i]),
                                     AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]));
                    mu_ym = HARMONIC(AVG2(mu[k-1][j  ][i], mu[k  ][j  ][i]),
                                     AVG2(mu[k-1][j-1][i], mu[k  ][j-1][i]));
                }

                double flux_xp = mu_xp * (x[k][j][i+1] - x[k][j][i  ]);
                double flux_xm = mu_xm * (x[k][j][i  ] - x[k][j][i-1]);
                double flux_yp = mu_yp * (x[k][j+1][i] - x[k][j  ][i]);
                double flux_ym = mu_ym * (x[k][j  ][i] - x[k][j-1][i]);
                double flux_zp = mu_zp * (x[k+1][j][i] - x[k  ][j][i]);
                double flux_zm = mu_zm * (x[k  ][j][i] - x[k-1][j][i]);

                double dFlux_dx = (flux_xp - flux_xm) / (dx * dx);
                double dFlux_dy = (flux_yp - flux_ym) / (dy * dy);
                double dFlux_dz = (flux_zp - flux_zm) / (dz * dz);

                Ax[k][j][i] = Ax_val - (dFlux_dx + dFlux_dy + dFlux_dz) / Re;
            }
        }
    }

#else
    double iddx    = grid->idx_c[1] * grid->idx_c[1];
    double iddy    = grid->idy_c[1] * grid->idy_c[1];
    double iddz    = grid->idz_c[1] * grid->idz_c[1];
    double iRe     = 1.0 / params->Re;
    double idtimeb = 1.0 / (BET[params->which_stage] * dt);

    double ac = idtimeb + 2.0*iRe*(iddx + iddy + iddz);
    double ax = -iRe * iddx;
    double ay = -iRe * iddy;
    double az = -iRe * iddz;

    for (k = Ks; k < Ke; k++)
        for (j = Js; j < Je; j++)
            for (i = Is; i < Ie; i++)
                Ax[k][j][i] = ac * x[k][j][i]
                            + ax * (x[k][j][i-1] + x[k][j][i+1])
                            + ay * (x[k][j-1][i] + x[k][j+1][i])
                            + az * (x[k-1][j][i] + x[k+1][j][i]);
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

    /* ------------------------------------------------------------------ *
     * Grid spacings and solver parameters (needed for Jacobi diagonal)
     * ------------------------------------------------------------------ */
    const double dx = 1.0 / grid->idx_c[1];
    const double dy = 1.0 / grid->idy_c[1];
    const double dz = 1.0 / grid->idz_c[1];
    const double inv_dx2 = 1.0 / (dx * dx);
    const double inv_dy2 = 1.0 / (dy * dy);
    const double inv_dz2 = 1.0 / (dz * dz);

    const double Re = params->Re;
    const double BET[] = { BETA };
    const double alpha_k = BET[params->which_stage];
    const double dt      = params->dt;

#ifdef VOF
    double ***rho = data_bag->vof->rho;
    double ***mu  = data_bag->vof->mu;
#else
    double ***rho = NULL;  /* unused in constant-coefficient path */
    double ***mu  = NULL;
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
                double diag = vel_jacobi_diag(i, j, k, component,
                                              rho, mu,
                                              inv_dx2, inv_dy2, inv_dz2,
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
                    double diag = vel_jacobi_diag(i, j, k, component,
                                                  rho, mu,
                                                  inv_dx2, inv_dy2, inv_dz2,
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
                    double diag = vel_jacobi_diag(i, j, k, component,
                                                  rho, mu,
                                                  inv_dx2, inv_dy2, inv_dz2,
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