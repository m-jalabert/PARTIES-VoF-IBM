/*
 * actm.c
 *
 *  Created on: Sep 27, 2013
 *      Author: daan
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"
#include "Particle.h"
//#include "../main.h" // For verbosity

#include "./actm.h"

#define verbose 0


/******************************************************************************/
/*
 */
/******************************************************************************/
void actm(Collision *pc, Collision_bag *bag, Parameters *params) {

	int i;
	double spring_force, normal_force[3];
	Particle *p  = bag -> p;
	Particle *p2 = bag -> p2;
	double kn = pc -> kn;
	double dn = pc -> dn;
	double surface_distance = bag -> surface_distance;
	double *n = bag -> n;
	double *gn = bag -> gn;

	// Magnitude of spring force (positive)
	spring_force = kn * pow(-surface_distance, 1.5);

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
 * Creates two ParticleCollisions on two particles and calculates the collision constants
 * Also adds the collisions to the end of the particles' linked list of collision objects
 * Assumes no ParticleCollision between these particles exists already (not really important)
 */
/******************************************************************************/
void actm_new_pp_collision(double *kn_ptr, double *dn_ptr, Collision_bag *bag,
		Parameters *params) {

	Particle  *p  = bag -> p;
	Particle  *p2 = bag -> p2;

	double zeta = bag -> surface_distance;
	double u_in = bag -> g_dot_n;
	double e    = params->e_dry_wall;

	// Calculate reduced mass
	double m;

	// If velocity high enough
	double u_in_crit, grav_mag, kn_grav;

	if (bag->stage == COLL_STAGE_FIXED) {
		m = p->M;
		u_in_crit = 4.5 * ST_CRIT / (params->Re * p->rho_s * p->R);
		
	}
	else {
		m = (p->M * p2->M) / (p->M + p2->M);
		u_in_crit = 4.5 * ST_CRIT / (params->Re * p->rho_s * p->R);
		u_in_crit = max(u_in_crit, 4.5 * ST_CRIT / (params->Re * p2->rho_s * p2->R));
	}

	// Collision time
	double Tc = params->Ndt_coll * params->dt;
#ifdef SUBSTEP
	Tc = Tc * 15.0;
#endif

	if(u_in >= u_in_crit) {

		// Calculate coefficients
		actm_calc_coeff(kn_ptr, dn_ptr, m, u_in, Tc, e);
	}
	else { // prevent 1/sqrt(0) if nu_f == 0

		// Calculate coefficients from critical velocity
		actm_calc_coeff(kn_ptr, dn_ptr, m, u_in_crit, Tc, e);

		// Stiffness required to support particle's mass under gravity
		grav_mag = sqrt( DOT(params->grav, params->grav) );
		kn_grav = fabs(p->M * grav_mag * pow(GRAV_OVERLAP * p->R, -1.5));
		kn_grav = max(kn_grav, fabs(p2->M * grav_mag * pow(GRAV_OVERLAP * p2->R, -1.5)));
		*kn_ptr = max(*kn_ptr, kn_grav);
	}

}




/******************************************************************************/
/*
 * Creates a ParticleCollision with a wall
 */
/******************************************************************************/
void actm_new_pw_collision(double *kn_ptr, double *dn_ptr, Collision_bag *bag,
		Parameters *params) {

	Particle  *p  = bag -> p;

	double zeta = bag -> surface_distance;
	double u_in = bag -> g_dot_n;
	double e    = params->e_dry_wall;

	// If velocity high enough
	double u_in_crit, grav_mag, kn_grav;

	// Collision time
	double Tc = params->Ndt_coll * params->dt;
#ifdef SUBSTEP
	Tc = Tc * 15.0;
#endif

 	u_in_crit  =  4.5  * ST_CRIT / (p->rho_s * p->R * params->Re);

	if (u_in >= u_in_crit) {

		// Calculate coefficients
		actm_calc_coeff(kn_ptr, dn_ptr, p->M, u_in, Tc, e);
	}
	else { // prevent 1/sqrt(0) if nu_f == 0

		// Calculate coefficients from critical velocity
		actm_calc_coeff(kn_ptr, dn_ptr, p->M, u_in_crit, Tc, e);

		// Stiffness required to support particle's mass under gravity
		grav_mag = sqrt( DOT(params->grav, params->grav) );
		kn_grav = fabs(p->M * grav_mag * pow(GRAV_OVERLAP * p->R, -1.5));
		*kn_ptr = max(*kn_ptr, kn_grav);
	}

}




/******************************************************************************/
void actm_calc_coeff(double *kn, double *dn, double M,  double gn,
		double Tc, double e) {

	long double lambda,  t_star;
	long double eta,     tau_c0,  alpha;
	double A,       B,       C;
	long double i_a2d2,  Tc_tau_c0;

 	eta       = pow(log(e),2);

	tau_c0    = 3.218;
	alpha     = 1.111;
	A         = 0.716;
	B         = 0.830;
	C         = 0.744;
	i_a2d2    = 1.0 / (pow(alpha,2)*pow(tau_c0,2));
	lambda    = i_a2d2 * (-0.5*C*eta+sqrt(0.25*pow(C,2)*pow(eta,2) + pow(alpha,2)*pow(tau_c0,2)*eta));
	Tc_tau_c0 = Tc / tau_c0;
	t_star    = Tc_tau_c0 * sqrt(1-A*lambda - B*pow(lambda,2));

  	*dn = 2.0 * lambda * M / t_star;
  	*kn = M / sqrt(gn * pow(t_star,5));

}




/******************************************************************************/


#undef verbose
