#include "TwodOps.h"

int TwodOps_collapsed_component_is_inactive(const Parameters *params)
{
	/*
	 * In both 2D modes the third storage direction remains allocated for MPI
	 * compatibility, but it carries no physical velocity component.
	 */
	return params->twod_mode_enabled;
}

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
                                      double idz_w)
{
	double dudx = (u_data[k][j][i + 1] - u_data[k][j][i]) * idx_u;
	double dvdy = (v_data[k][j + 1][i] - v_data[k][j][i]) * idy_v;

	if (!params->twod_mode_enabled) {
		double dwdz = (w_data[k + 1][j][i] - w_data[k][j][i]) * idz_w;
		return dudx + dvdy + dwdz;
	}

	if (params->twod_cartesian_enabled)
		/* Planar Cartesian divergence uses only the in-plane velocity. */
		return dudx + dvdy;

	/*
	 * Axisymmetric mode reuses x as the radial coordinate. The scalar control
	 * volume is weighted by r, so the divergence becomes 1/r d(r u_r)/dr.
	 */
	return (grid->r_u[i + 1] * u_data[k][j][i + 1] -
	        grid->r_u[i] * u_data[k][j][i]) * idx_u * grid->inv_r_c[i] +
	       dvdy;
}

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
                                      double *czp)
{
	const double idx2 = grid->idx_c[1] * grid->idx_c[1];
	const double idy2 = grid->idy_c[1] * grid->idy_c[1];
	const double idz2 = grid->idz_c[1] * grid->idz_c[1];

	(void)j;
	(void)k;

	*cxm = beta_xm * idx2;
	*cxp = beta_xp * idx2;
	*cym = beta_ym * idy2;
	*cyp = beta_yp * idy2;
	*czm = beta_zm * idz2;
	*czp = beta_zp * idz2;

	if (!params->twod_mode_enabled)
		return;

	/* The collapsed slab is storage-only in 2D mode, so z-couplings vanish. */
	*czm = 0.0;
	*czp = 0.0;

	if (params->twod_cartesian_enabled)
		return;

	/*
	 * Axisymmetric scalar operators keep the same stencil shape but weight the
	 * radial faces by r_face / r_cell.
	 */
	*cxm = grid->r_u[i] * beta_xm * grid->inv_r_c[i] * idx2;
	*cxp = grid->r_u[i + 1] * beta_xp * grid->inv_r_c[i] * idx2;
}

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
                                           double beta_zp)
{
	double cxm, cxp, cym, cyp, czm, czp;

	TwodOps_build_scalar_face_coeffs(grid, params, i, j, k,
	                                 beta_xm, beta_xp, beta_ym, beta_yp,
	                                 beta_zm, beta_zp,
	                                 &cxm, &cxp, &cym, &cyp, &czm, &czp);
	return cxm + cxp + cym + cyp + czm + czp;
}

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
                                     double beta_zp)
{
	double cxm, cxp, cym, cyp, czm, czp;

	TwodOps_build_scalar_face_coeffs(grid, params, i, j, k,
	                                 beta_xm, beta_xp, beta_ym, beta_yp,
	                                 beta_zm, beta_zp,
	                                 &cxm, &cxp, &cym, &cyp, &czm, &czp);

	return cxp * (phi[k][j][i + 1] - phi[k][j][i]) -
	       cxm * (phi[k][j][i] - phi[k][j][i - 1]) +
	       cyp * (phi[k][j + 1][i] - phi[k][j][i]) -
	       cym * (phi[k][j][i] - phi[k][j - 1][i]) +
	       czp * (phi[k + 1][j][i] - phi[k][j][i]) -
	       czm * (phi[k][j][i] - phi[k - 1][j][i]);
}

double TwodOps_pressure_row_weight(const MAC_grid *grid,
                                   const Parameters *params,
                                   int i)
{
	if (!params->axisym_rz_enabled)
		return 1.0;

	/*
	 * Multiplying the axisymmetric scalar equation by r_cell keeps the same
	 * solution but restores symmetry of the discrete Poisson matrix, which is
	 * what the CG/HYPRE-PCG pressure solvers require.
	 */
	return grid->r_c[i];
}

double TwodOps_scalar_laplacian(const MAC_grid *grid,
                                const Parameters *params,
                                double ***phi,
                                int i,
                                int j,
                                int k)
{
	double dfdxE = (phi[k][j][i + 1] - phi[k][j][i]) * grid->idx_c[i];
	double dfdxW = (phi[k][j][i] - phi[k][j][i - 1]) *
	               ((i != 0) ? grid->idx_c[i - 1] : grid->idx_c[i]);

	double dfdyN = (phi[k][j + 1][i] - phi[k][j][i]) * grid->idy_c[j];
	double dfdyS = (phi[k][j][i] - phi[k][j - 1][i]) *
	               ((j != 0) ? grid->idy_c[j - 1] : grid->idy_c[j]);

	if (!params->twod_mode_enabled) {
		double dfdzF = (phi[k + 1][j][i] - phi[k][j][i]) * grid->idz_c[k];
		double dfdzB = (phi[k][j][i] - phi[k - 1][j][i]) *
		               ((k != 0) ? grid->idz_c[k - 1] : grid->idz_c[k]);

		return (dfdxE - dfdxW) * grid->idx_u[i] +
		       (dfdyN - dfdyS) * grid->idy_v[j] +
		       (dfdzF - dfdzB) * grid->idz_w[k];
	}

	if (params->twod_cartesian_enabled)
		/* Planar Cartesian Laplacian: the stored slab has no z diffusion. */
		return (dfdxE - dfdxW) * grid->idx_u[i] +
		       (dfdyN - dfdyS) * grid->idy_v[j];

	/*
	 * Axisymmetric scalar Laplacian:
	 *   1/r d/dr (r dphi/dr) + d2phi/dz^2
	 *
	 * The west radial face coincides with the axis for the first cell, so
	 * r_u[i] = 0 naturally removes the singular face contribution there.
	 */
	return (grid->r_u[i + 1] * dfdxE - grid->r_u[i] * dfdxW) *
	       grid->idx_u[i] * grid->inv_r_c[i] +
	       (dfdyN - dfdyS) * grid->idy_v[j];
}

double TwodOps_scalar_flux_divergence(const MAC_grid *grid,
                                      const Parameters *params,
                                      double ***flux_x,
                                      double ***flux_y,
                                      double ***flux_z,
                                      int i,
                                      int j,
                                      int k)
{
	double div_x = (flux_x[k][j][i + 1] - flux_x[k][j][i]) / grid->dx_c[i];
	double div_y = (flux_y[k][j + 1][i] - flux_y[k][j][i]) / grid->dy_c[j];

	if (!params->twod_mode_enabled) {
		double div_z = (flux_z[k + 1][j][i] - flux_z[k][j][i]) / grid->dz_c[k];
		return div_x + div_y + div_z;
	}

	if (params->twod_cartesian_enabled)
		/* In planar mode only the x-y transport is physical. */
		return div_x + div_y;

	/*
	 * Axisymmetric advective flux divergence:
	 *   1/r d/dr (r F_r) + dF_z/dz
	 */
	return (grid->r_u[i + 1] * flux_x[k][j][i + 1] -
	        grid->r_u[i] * flux_x[k][j][i]) * grid->inv_r_c[i] / grid->dx_c[i] +
	       div_y;
}

double TwodOps_u_cv_x_flux_divergence(const MAC_grid *grid,
                                      const Parameters *params,
                                      int i,
                                      double flux_e,
                                      double flux_w)
{
	if (!params->axisym_rz_enabled)
		return (flux_e - flux_w) * grid->idx_c[i - 1];

	/*
	 * Radial flux divergence on the u control volume:
	 *   1/r_u * d(r * flux_r)/dr
	 *
	 * The east/west faces of the u CV lie at cell centers r_c[i] and r_c[i-1].
	 */
	return (grid->r_c[i] * flux_e - grid->r_c[i - 1] * flux_w) *
	       grid->idx_c[i - 1] * grid->inv_r_u[i];
}

double TwodOps_v_cv_x_flux_divergence(const MAC_grid *grid,
                                      const Parameters *params,
                                      int i,
                                      double flux_e,
                                      double flux_w)
{
	if (!params->axisym_rz_enabled)
		return (flux_e - flux_w) * grid->idx_u[i];

	/*
	 * Radial flux divergence on the v control volume:
	 *   1/r_c * d(r * flux_r)/dr
	 *
	 * The east/west faces of the v CV lie on u locations r_u[i+1] and r_u[i].
	 */
	return (grid->r_u[i + 1] * flux_e - grid->r_u[i] * flux_w) *
	       grid->idx_u[i] * grid->inv_r_c[i];
}

double TwodOps_u_cv_x_diag_coeff(const MAC_grid *grid,
                                 const Parameters *params,
                                 int i,
                                 double coeff_e,
                                 double coeff_w)
{
	if (!params->axisym_rz_enabled)
		return (coeff_e + coeff_w) * grid->idx_c[i - 1];

	return (grid->r_c[i] * coeff_e + grid->r_c[i - 1] * coeff_w) *
	       grid->idx_c[i - 1] * grid->inv_r_u[i];
}

double TwodOps_v_cv_x_diag_coeff(const MAC_grid *grid,
                                 const Parameters *params,
                                 int i,
                                 double coeff_e,
                                 double coeff_w)
{
	if (!params->axisym_rz_enabled)
		return (coeff_e + coeff_w) * grid->idx_u[i];

	return (grid->r_u[i + 1] * coeff_e + grid->r_u[i] * coeff_w) *
	       grid->idx_u[i] * grid->inv_r_c[i];
}

double TwodOps_axisym_u_radial_linear_term(const MAC_grid *grid,
                                           const Parameters *params,
                                           int i,
                                           double coeff_at_u,
                                           double u_value)
{
	if (!params->axisym_rz_enabled)
		return 0.0;

	/*
	 * Cylindrical radial momentum contains the extra linear term
	 *   -2 * coeff_at_u * u_r / r_u^2
	 * which has no Cartesian analogue.
	 */
	return -2.0 * coeff_at_u * u_value * grid->inv_r_u[i] * grid->inv_r_u[i];
}

double TwodOps_axisym_u_radial_linear_coeff(const MAC_grid *grid,
                                            const Parameters *params,
                                            int i,
                                            double coeff_at_u)
{
	if (!params->axisym_rz_enabled)
		return 0.0;

	return 2.0 * coeff_at_u * grid->inv_r_u[i] * grid->inv_r_u[i];
}

double TwodOps_cell_measure_c(const MAC_grid *grid,
                              const Parameters *params,
                              int i,
                              int j,
                              int k)
{
	if (!params->twod_mode_enabled)
		return grid->dx_u[i] * grid->dy_v[j] * grid->dz_w[k];

	if (params->twod_cartesian_enabled) {
		/*
		 * The collapsed slab thickness is pure bookkeeping in planar 2D. We keep
		 * it explicit here so any normalized diagnostic can divide it back out.
		 */
		return grid->dx_u[i] * grid->dy_v[j] * grid->dummy_z_slab_thickness;
	}

	/*
	 * Axisymmetric cell measure = (meridional area) * (azimuthal arc length)
	 *                           = dr * dz * theta_span * r_cell.
	 */
	return grid->dx_u[i] * grid->dy_v[j] * grid->ring_wt_c[i];
}

double TwodOps_domain_measure(const MAC_grid *grid,
                              const Parameters *params)
{
	(void)grid;

	if (!params->twod_mode_enabled)
		return params->Lx * params->Ly * params->Lz;

	if (params->twod_cartesian_enabled)
		return params->Lx * params->Ly * params->twod_slab_thickness;

	return 0.5 * params->axisym_theta_span *
	       (params->xmax * params->xmax - params->xmin * params->xmin) *
	       params->Ly;
}
