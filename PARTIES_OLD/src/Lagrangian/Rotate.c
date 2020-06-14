#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Rotate.h"

/******************************************************************************/
/*
 * Reference paper: Integrating rotation fron angular velocity - Eva Zupan,
 * Miran Saje.  Advances in Engineering Software : May,2011
 *
 * The time integration scheme that has been used, is described in section 4.1
 * in the above mentioned paper as MP-R
 *
 * NOTE: Our implementation has changed from the paper slightly since we use an
 * angular velocity in the reference frame of the laboratory rather than the
 * particle.  The implementation below "seems" to work, but it is not a
 * mathematically-rigorous derivation.
 */
/******************************************************************************/
void Rotate_particle(Particle* p, Parameters* params) {

	double *Omega, *Omega_old;

	const double BET[] = {BETA};
	double bet = BET[params -> which_stage];
	double dt  = params -> dt;

	Omega = p -> Omega;
	Omega_old = p -> Omega_old;

	// Identity matrix
	double Id[3][3];
	Id[0][0] = 1;
	Id[1][1] = 1;
	Id[2][2] = 1;
	Id[0][1] = 0;
	Id[0][2] = 0;
	Id[1][0] = 0;
	Id[1][2] = 0;
	Id[2][0] = 0;
	Id[2][1] = 0;

	// omega_mdpt = 0.5*(Omega + Omega_old)
	double omega_mdpt[3];
	omega_mdpt[0] = 0.5*(Omega[0] + Omega_old[0]);
	omega_mdpt[1] = 0.5*(Omega[1] + Omega_old[1]);
	omega_mdpt[2] = 0.5*(Omega[2] + Omega_old[2]);

	// omega_cpt = 0.5*(skew(Omega) + skew(Omega_old))
	double omega_cpt[3][3];
	omega_cpt[0][0] = 0;
	omega_cpt[1][1] = 0;
	omega_cpt[2][2] = 0;
	omega_cpt[0][1] = -0.5*(Omega[2] + Omega_old[2]);
	omega_cpt[0][2] =  0.5*(Omega[1] + Omega_old[1]);
	omega_cpt[1][0] =  0.5*(Omega[2] + Omega_old[2]);
	omega_cpt[1][2] = -0.5*(Omega[0] + Omega_old[0]);
	omega_cpt[2][0] = -0.5*(Omega[1] + Omega_old[1]);
	omega_cpt[2][1] =  0.5*(Omega[0] + Omega_old[0]);

	// omega_m2 = square of magnitude of Omega_mdpt
	double omega_m2;
	omega_m2 = DOT(omega_mdpt, omega_mdpt);

	// product = omega_cpt^2
	double product[3][3];
	int cp, dp, zp;
	double sum2;
	for (cp = 0; cp < 3; cp++) {
		for (dp = 0; dp < 3; dp++) {
			sum2 = 0;
			for (zp = 0; zp < 3; zp++) {
				sum2 = sum2 + omega_cpt[cp][zp] * omega_cpt[zp][dp];
			}
			product[cp][dp] = sum2;
		}
	}

	//--------------------------------------------------------------------------
	// The calculation done in the next section is referred in equation 22 in
	// the above mentioned paper
	//--------------------------------------------------------------------------
	double dtt = 2. * dt * bet;
	double tm1 = 2. / (1. + 0.5*dtt*dtt*omega_m2);
	double tm3[3][3];
	int id, jd;
	for (id = 0; id < 3; id++) {
		for (jd = 0; jd < 3; jd++) {
			tm3[id][jd] = Id[id][jd] + tm1 * (0.5*dtt*omega_cpt[id][jd]
				+ 0.25*dtt*dtt*product[id][jd]);
		}
	}

	double sum3;
	int cr,dr,zr;
	for (cr = 0; cr < 3; cr++) {
		for (dr = 0; dr < 3; dr++) {
			sum3 = 0;
			for (zr = 0; zr < 3; zr++) {
//				sum3 = sum3 + p -> Rotn_old[cr][zr] * tm3[zr][dr];
				// This seems to fix the original problem
				sum3 = sum3 + tm3[cr][zr] * p -> Rotn_old[zr][dr];
			}
			p -> Rotn[cr][dr] = sum3;
		}
	}

//	int c1;
//	int c2;
//	printf("................\n");
//	for (c1=0; c1<3 ;c1++){
//		for (c2=0;c2<3;c2++){
//			printf("%g, ",p -> Rotn[c1][c2]);
//		}
//		printf("\n");
//	}

}
