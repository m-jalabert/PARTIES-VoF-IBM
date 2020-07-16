/*
 *  erm.c
 *
 *  Created on: January, 2017
 *  Author: B. Vowinckel
 *  Implements the Electrostatic Repulsion Model as
 *  proposed by Mari et al. (JoR, 2014)
 *  F_R = F_ER (R_eff / R_p) * exp(zeta_n/lambda)
 */
#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"

/******************************************************************************/
/*
 * Implements the Electrostatic Repulsion Model as proposed by Mari et al.
 * (Journal of Rheology, 2014)
 *
 *     F_R = F_ER (R_eff / R_mean) * exp(zeta_n/lambda)
 */
/******************************************************************************/
void erm(Collision_bag *bag, Parameters *params) {

	int i;

	Particle *p = bag -> p;
	Particle *p2 = bag -> p2;

	// Normal vector and distance between surfaces
	double *n   = bag -> n;
	double zeta = bag -> surface_distance;

	// Effective and mean radii
	double R_eff, R_mean;
	if (p2 == NULL) {
		R_eff = p->R;
		R_mean = p->R;
	}
	else {
		R_eff = p->R * p2->R / (p->R + p2->R);
		R_mean = 0.5 * (p->R + p2->R);
	}

	// Debye length of 0.05*R
	double kappa = 1./(0.05 * R_mean);

	// Evaluate electrostatic repulsive force
	double F_erm = -params->F_erm * (R_eff / R_mean) * exp(-kappa * zeta);

	// Store electrostatic repulsive force
	FORI3 p->Fc[i] += F_erm * n[i];
	if (p2 != NULL) {
		FORI3 p2->Fc[i] -= F_erm * n[i];
	}
}
