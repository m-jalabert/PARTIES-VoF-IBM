
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"


/******************************************************************************/
/*
 Assumes uniform square grid
 */
/******************************************************************************/
void dem(Collision *pc, Collision_bag *bag, Parameters *params) {

	int i;
	double spring_force, normal_force[3];
	Particle *p  = bag -> p;
	Particle *p2 = bag -> p2;
	double dn  = pc -> dn;
	double *n  = bag -> n;
	double *gn = bag -> gn;

	spring_force = -(pc->kn) * bag->surface_distance;

	FORI3 normal_force[i] = -n[i] * spring_force - dn * gn[i];
	FORI3 p->Fc[i] += normal_force[i];
#ifdef POST_PROCESS
	FORI3 p->Fc_norm[i] += normal_force[i];
#endif
	if (p2 != NULL) {
		FORI3 p2->Fc[i] -= normal_force[i];
#ifdef POST_PROCESS
		FORI3 p2->Fc_norm[i] -= normal_force[i];
#endif
	}

	bag->normal_force_norm = fabs(DOT(n, normal_force));
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void dem_new_pp_collision(double *kn_ptr, double *dn_ptr, Collision_bag *bag, Parameters *params) {

	Particle  *p  = bag -> p;
	Particle  *p2 = bag -> p2;

	double Tc = params->Ndt_coll * params->dt;
#ifdef SUBSTEP
	Tc = Tc * 15.0;
#endif

	double ln2e = log(params -> e_dry_particles);
	ln2e = ln2e * ln2e;
	ln2e = ln2e / (PI * PI + ln2e);

	// Reduced mass
	double m;
	if (bag->stage == COLL_STAGE_FIXED) {
		m = p->M;
	}
	else {
		m = (p->M * p2->M) / (p->M + p2->M);
	}

	*kn_ptr = PI * PI * m / ( Tc * Tc * ( 1.0 - ln2e ) );
	*dn_ptr = 2.0 * sqrt( m * (*kn_ptr) * ln2e );
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void dem_new_pw_collision(double *kn_ptr, double *dn_ptr, Collision_bag *bag, Parameters *params) {

	Particle  *p  = bag -> p;

	double Tc = params->Ndt_coll * params->dt;
#ifdef SUBSTEP
	Tc = Tc * 15.0;
#endif

	double ln2e = log(params -> e_dry_wall);
	ln2e = ln2e * ln2e;
	ln2e = ln2e / (PI * PI + ln2e);

	// Reduced mass
	double m = p -> M;

	*kn_ptr = PI * PI * m / ( Tc * Tc * ( 1.0 - ln2e ) );
	*dn_ptr = 2.0 * sqrt( m * (*kn_ptr) * ln2e );
}
