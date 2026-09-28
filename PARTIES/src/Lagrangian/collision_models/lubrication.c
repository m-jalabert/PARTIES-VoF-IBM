
#include <math.h>

#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"

#if defined(TWOD_CARTESIAN) && defined(VOF) && defined(VOF_DIFFUSE)
static double lubrication_clamp(double value, double lower, double upper) {
	if (value < lower) return lower;
	if (value > upper) return upper;
	return value;
}

static double lubrication_delta(double r) {
	r = fabs(r);

	if (r <= 0.5)
		return (1.0 + sqrt(-3.0 * r * r + 1.0)) / 3.0;
	else if (r <= 1.5)
		return (5.0 - 3.0 * r - sqrt(-3.0 * (1.0 - r) * (1.0 - r) + 1.0)) / 6.0;
	else
		return 0.0;
}

static void lubrication_update_twod_multiphase(Collision_bag *bag, double h,
		Parameters *params, MAC_grid *grid, VolumeFraction *vof) {

	if (grid == NULL || vof == NULL || vof->C_L == NULL)
		return;

	double *x = grid->xc;
	double *y = grid->yc;
	const int k = grid->G_Ks;

	if (h <= 0.0)
		h = grid->dx_u[1];
	if (h <= 0.0 || k < grid->L_Ks || k >= grid->L_Ke)
		return;

	const double *X = bag->contact_point;
	int i_start = (int)round((X[0] - x[0]) / h) - 1;
	int j_start = (int)round((X[1] - y[0]) / h) - 1;
	int i_end = i_start + 3;
	int j_end = j_start + 3;

	i_start = max(i_start, grid->L_Is);
	j_start = max(j_start, grid->L_Js);
	i_end = min(i_end, grid->L_Ie);
	j_end = min(j_end, grid->L_Je);

	const double mu_g = (params->mu1 > 0.0) ? params->mu2 / params->mu1 : 1.0;
	double weight_sum = 0.0;
	double liquid_sum = 0.0;
	double mu_sum = 0.0;

	for (int j = j_start; j < j_end; ++j) {
		const double wy = lubrication_delta((X[1] - y[j]) / h);
		if (wy <= 0.0)
			continue;

		for (int i = i_start; i < i_end; ++i) {
			const double wx = lubrication_delta((X[0] - x[i]) / h);
			const double weight = wx * wy;

			if (weight <= 0.0)
				continue;

			const double cl = lubrication_clamp(vof->C_L[k][j][i], 0.0, 1.0);
			const double cg = (vof->C_G != NULL) ?
				lubrication_clamp(vof->C_G[k][j][i], 0.0, 1.0) :
				lubrication_clamp(1.0 - cl, 0.0, 1.0);
			// Exclude diffuse solid tails from the fluid-phase viscosity mix.
			const double fluid = cl + cg;
			const double liquid_fraction = (fluid > 1.0e-14) ?
				lubrication_clamp(cl / fluid, 0.0, 1.0) : cl;
			const double mu_ratio =
				liquid_fraction + (1.0 - liquid_fraction) * mu_g;

			weight_sum += weight;
			liquid_sum += weight * liquid_fraction;
			mu_sum += weight * mu_ratio;
		}
	}

	if (weight_sum > 0.0) {
		bag->lub_liquid_fraction = liquid_sum / weight_sum;
		bag->lub_viscosity_ratio = mu_sum / weight_sum;
	}
}
#else
static void lubrication_update_twod_multiphase(Collision_bag *bag, double h,
		Parameters *params, MAC_grid *grid, VolumeFraction *vof) {
	(void)bag;
	(void)h;
	(void)params;
	(void)grid;
	(void)vof;
}
#endif

/******************************************************************************/
void lubrication(Collision_bag *bag, double h, Parameters *params,
		MAC_grid *grid, VolumeFraction *vof) {

	int i;
	double temp, R_eff, R_mean;
	double Flub[3], Tlub[3];

	Particle *p = bag -> p;
	Particle *p2 = bag -> p2;

	// Effective and mean radii
	if (p2 == NULL) {
		R_eff = p->R;
		R_mean = p->R;
	}
	else {
		R_eff = p->R * p2->R / (p->R + p2->R);
		R_mean = 0.5 * (p->R + p2->R);
	}

	// Limit minimum surface distance used to calculate lubrication force based
	// on surface roughness distance
	double abs_roughness = params->roughness * R_mean;
	double surface_distance = max(bag->surface_distance, abs_roughness);

#if defined(TWOD_CARTESIAN) && \
	(defined(LUBRICATION_NORMAL) || defined(LUBRICATION_TANGENTIAL))
	lubrication_update_twod_multiphase(bag, h, params, grid, vof);
	double mu_ratio = bag->lub_viscosity_ratio;
	double slab_thickness = params->twod_slab_thickness;

	if (grid != NULL && grid->dummy_z_slab_thickness > 0.0)
		slab_thickness = grid->dummy_z_slab_thickness;
	if (slab_thickness <= 0.0)
		slab_thickness = h;
	if (mu_ratio <= 0.0)
		mu_ratio = 1.0;
#else
	(void)h;
	(void)grid;
	(void)vof;
#endif

#ifdef LUBRICATION_NORMAL
	// Relative translational velocity normal to contact surface
	double *gn = bag -> gn;

	// Lubrication force, less velocity component
#ifdef TWOD_CARTESIAN
	// 2D parallel-cylinder squeeze film, integrated over the storage slab.
	temp = -3.0 * sqrt(2.0) * PI * mu_ratio / params->Re
	       * pow(R_eff, 1.5) / pow(surface_distance, 1.5)
	       * slab_thickness;
#else
	temp = -6.0 * PI / (surface_distance * params->Re) * R_eff * R_eff;
#endif

	// Evaluate and store lubrication force
	FORI3 Flub[i] = temp * gn[i];
#ifdef TWOD_CARTESIAN
	Flub[2] = 0.0;
#endif
	FORI3 p->Fc[i]  += Flub[i];
#ifdef POST_PROCESS
	FORI3 p->Fl_norm[i]  += Flub[i];
#endif
	if (p2 != NULL) {
		FORI3 p2->Fc[i] -= Flub[i];
#ifdef POST_PROCESS
		FORI3 p2->Fl_norm[i] -= Flub[i];
#endif
	}
#endif

#ifdef LUBRICATION_TANGENTIAL
#ifdef TWOD_CARTESIAN
	double *n = bag -> n;
	double *gt_cp = bag -> gt_cp;
	double torque_base;

	// Leading tangential Couette shear in the parabolic line-contact gap.
	temp = -sqrt(2.0) * PI * mu_ratio / params->Re
	       * sqrt(R_eff / surface_distance) * slab_thickness;

	FORI3 Flub[i] = temp * gt_cp[i];
	Flub[2] = 0.0;
	FORI3 Tlub[i] = 0.0;
	torque_base = n[0] * Flub[1] - n[1] * Flub[0];
	Tlub[2] = bag->R_cp * torque_base;

	p->Fc[0] += Flub[0];
	p->Fc[1] += Flub[1];
	p->Tc[2] += Tlub[2];
#ifdef POST_PROCESS
	p->Fl_tan[0] += Flub[0];
	p->Fl_tan[1] += Flub[1];
#endif
	if (p2 != NULL) {
		p2->Fc[0] -= Flub[0];
		p2->Fc[1] -= Flub[1];
		p2->Tc[2] += bag->R2_cp * torque_base;
#ifdef POST_PROCESS
		p2->Fl_tan[0] -= Flub[0];
		p2->Fl_tan[1] -= Flub[1];
#endif
	}
#else
	// Lubrication force tangential
	double gt_cp_lub[3], gt_cross_n[3];
	double Ft, Fr, Tt, Tr;
	double nu = 1./params->Re;
	double *n = bag -> n;
	double *gt = bag -> gt;
	double *Om_cross_R = bag -> Om_cross_R;

	// Coefficients according to Goldman, Cox, Brenner (CES, 1967)
	double ln_delta = log(surface_distance / R_eff);
	Ft =  (8./15.) * ln_delta - 0.9588;
	Fr = -(2./15.) * ln_delta - 0.2526;
	Tt = -(1./10.) * ln_delta - 0.1895;
	Tr =  (2./ 5.) * ln_delta - 0.3817;

	// Negative sign different to offset negative sign in 'Tr'
	gt_cp_lub[0] = (Tt * gt[0] - Tr * Om_cross_R[0]);
	gt_cp_lub[1] = (Tt * gt[1] - Tr * Om_cross_R[1]);
	gt_cp_lub[2] = (Tt * gt[2] - Tr * Om_cross_R[2]);

	gt_cross_n[0] = (gt_cp_lub[1] * n[2] - gt_cp_lub[2] * n[1]);
	gt_cross_n[1] = (gt_cp_lub[2] * n[0] - gt_cp_lub[0] * n[2]);
	gt_cross_n[2] = (gt_cp_lub[0] * n[1] - gt_cp_lub[1] * n[0]);

	// Evaluate lubrication force and torque
	FORI3 Flub[i] = blend * 6. * PI * nu * R_eff * (Ft * gt[i] + Fr * Om_cross_R[i]);
	FORI3 Tlub[i] = blend * 8. * PI * nu * R_eff * R_eff * gt_cross_n[i];

	// Store lubrication force and torque
	FORI3 p->Fc[i] += Flub[i];
	FORI3 p->Tc[i] += Tlub[i];
#ifdef POST_PROCESS
	FORI3 p->Fl_tan[i] += Flub[i];
#endif
	if (p2 != NULL) {
		FORI3 p2->Fc[i] -= Flub[i];
		FORI3 p2->Tc[i] += Tlub[i];
#ifdef POST_PROCESS
		FORI3 p2->Fl_tan[i] -= Flub[i];
#endif
	}
#endif
#endif
}
