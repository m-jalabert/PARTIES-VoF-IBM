/*
 * Cohesion.c
 * Author:  J.Withers
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"
#include "Particle.h"
/******************************************************************************/
void cohesion(double lambda, Collision_bag *bag, Parameters *params) {

	int i;
	double VDW, R_eff, cohesive_force[3], k_coh; 

	Particle *p  = bag -> p;
	Particle *p2 = bag -> p2;
	double zeta = bag->surface_distance;
	double D50 = params -> D50;
	double Co = params -> Co;
	double *n = bag -> n;
        double  M50, grav_norm;
	double rho_s = params->rho_s;
	
	// Effective radius particle-wall or particle mobile-fixed	
	if (bag->stage == COLL_STAGE_WALL || bag->stage == COLL_STAGE_FIXED)  R_eff = p->R;
	// Effective radius particle-particle
	else  R_eff = p->R * p2->R / (p->R + p2->R);	
	
// 	Scaling with inertial force
        grav_norm = sqrt(params->grav[0] * params->grav[0] 
                       + params->grav[1] * params->grav[1]
                       + params->grav[2] * params->grav[2]);
	M50       = (rho_s-1) * PI * D50*D50*D50 / 6.0;

// 	Cohesive stiffness constant
	k_coh     = -Co * 8.0 / (lambda * lambda) ;
	
// 	Compute cohesive forces
        VDW  =  k_coh * M50 * grav_norm * R_eff * (zeta * zeta - lambda * zeta);
	
	//Cohesive_force directed from p to p2 
	FORI3 cohesive_force[i] = n[i] * VDW; 

	//Add cohesive force to contact force
	FORI3 p->Fc[i] += cohesive_force[i];
	FORI3 p2->Fc[i] -= cohesive_force[i]; 

}
