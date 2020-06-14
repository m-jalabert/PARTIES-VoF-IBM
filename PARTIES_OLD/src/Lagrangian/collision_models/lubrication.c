
#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"

/******************************************************************************/
void lubrication(Collision_bag *bag, double h, Parameters *params) {

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

	// Blend lubrication force from a value near zero at surface_distance=2h to
	// a value near one at surface_distance=h
	double blend = 0.5 * erf(6.0 - 4.0*surface_distance/h) + 0.5;

#ifdef LUBRICATION_NORMAL
	// Relative translational velocity normal to contact surface
	double *gn = bag -> gn;

	// Lubrication force, less velocity component
	temp = -6.0 * blend * PI / (surface_distance * params->Re) * R_eff * R_eff;

	// Evaluate and store lubrication force
	FORI3 Flub[i] = temp * gn[i];
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
}
