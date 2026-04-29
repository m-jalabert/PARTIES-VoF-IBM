#ifndef TWOD_OPS_H
#define TWOD_OPS_H

#include "DataTypes.h"

/*
 * Shared 2D operator helpers for the Eulerian layer.
 *
 * The current code keeps the legacy 3D storage layout even in 2D mode, so
 * these helpers provide one place where the physical operator changes from:
 *   - full 3D Cartesian when TWOD_MODE is off
 *   - planar x-y Cartesian when TWOD_CARTESIAN is on
 *   - r-z-like scalar operators when AXISYM_RZ is on
 *
 * The goal is to keep the mode switch explicit and easy to audit instead of
 * repeating small ad-hoc formula changes in every solver.
 */
int TwodOps_collapsed_component_is_inactive(const Parameters *params);

double TwodOps_cell_center_divergence(const MAC_grid *grid,
                                      const Parameters *params,
                                      double ***u_data,
                                      double ***v_data,
                                      double ***w_data,
                                      int i,
                                      int j,
                                      int k,
                                      double idx_u,
                                      double idy_v,
                                      double idz_w);

void TwodOps_build_scalar_face_coeffs(const MAC_grid *grid,
                                      const Parameters *params,
                                      int i,
                                      int j,
                                      int k,
                                      double beta_xm,
                                      double beta_xp,
                                      double beta_ym,
                                      double beta_yp,
                                      double beta_zm,
                                      double beta_zp,
                                      double *cxm,
                                      double *cxp,
                                      double *cym,
                                      double *cyp,
                                      double *czm,
                                      double *czp);

double TwodOps_scalar_diag_from_face_betas(const MAC_grid *grid,
                                           const Parameters *params,
                                           int i,
                                           int j,
                                           int k,
                                           double beta_xm,
                                           double beta_xp,
                                           double beta_ym,
                                           double beta_yp,
                                           double beta_zm,
                                           double beta_zp);

double TwodOps_apply_scalar_operator(const MAC_grid *grid,
                                     const Parameters *params,
                                     double ***phi,
                                     int i,
                                     int j,
                                     int k,
                                     double beta_xm,
                                     double beta_xp,
                                     double beta_ym,
                                     double beta_yp,
                                     double beta_zm,
                                     double beta_zp);

/*
 * Positive row weight that turns the AXISYM_RZ pressure operator into an SPD
 * system for CG / HYPRE-PCG by multiplying each scalar control-volume equation
 * by r_cell. In Cartesian / legacy 3D mode the weight is 1.
 */
double TwodOps_pressure_row_weight(const MAC_grid *grid,
                                   const Parameters *params,
                                   int i);

/*
 * Geometry-aware scalar helpers used by the VOF_DIFFUSE branch.
 *
 * These routines keep the legacy 3D array layout but apply only the physical
 * in-plane operator in 2D mode:
 *   - x-y Laplacian / flux divergence for TWOD_CARTESIAN
 *   - r-z Laplacian / flux divergence for AXISYM_RZ
 *   - full 3D operator when TWOD_MODE is off
 */
double TwodOps_scalar_laplacian(const MAC_grid *grid,
                                const Parameters *params,
                                double ***phi,
                                int i,
                                int j,
                                int k);

double TwodOps_scalar_flux_divergence(const MAC_grid *grid,
                                      const Parameters *params,
                                      double ***flux_x,
                                      double ***flux_y,
                                      double ***flux_z,
                                      int i,
                                      int j,
                                      int k);

/*
 * Momentum helpers for the x/r direction on staggered velocity control
 * volumes. These keep the metric weighting in one place so the explicit
 * momentum kernels and the implicit CG operator use the exact same geometry.
 *
 * In Cartesian / legacy 3D mode they reduce to the original x-divergence.
 * In AXISYM_RZ they apply:
 *   - u-CV  : 1/r_u d(r_c * flux)/dr
 *   - v-CV  : 1/r_c d(r_u * flux)/dr
 *
 * The radial linear term is the extra cylindrical correction in the radial
 * momentum equation:
 *   -2 * coeff_at_u * u_r / r_u^2
 */
double TwodOps_u_cv_x_flux_divergence(const MAC_grid *grid,
                                      const Parameters *params,
                                      int i,
                                      double flux_e,
                                      double flux_w);

double TwodOps_v_cv_x_flux_divergence(const MAC_grid *grid,
                                      const Parameters *params,
                                      int i,
                                      double flux_e,
                                      double flux_w);

double TwodOps_u_cv_x_diag_coeff(const MAC_grid *grid,
                                 const Parameters *params,
                                 int i,
                                 double coeff_e,
                                 double coeff_w);

double TwodOps_v_cv_x_diag_coeff(const MAC_grid *grid,
                                 const Parameters *params,
                                 int i,
                                 double coeff_e,
                                 double coeff_w);

double TwodOps_axisym_u_radial_linear_term(const MAC_grid *grid,
                                           const Parameters *params,
                                           int i,
                                           double coeff_at_u,
                                           double u_value);

double TwodOps_axisym_u_radial_linear_coeff(const MAC_grid *grid,
                                            const Parameters *params,
                                            int i,
                                            double coeff_at_u);

/*
 * Physical control-volume measures used by diagnostics and conserved
 * quantities. These keep the solver storage 3D-shaped while returning the
 * correct physical measure of a cell/domain in each mode.
 */
double TwodOps_cell_measure_c(const MAC_grid *grid,
                              const Parameters *params,
                              int i,
                              int j,
                              int k);

double TwodOps_domain_measure(const MAC_grid *grid,
                              const Parameters *params);

#endif
