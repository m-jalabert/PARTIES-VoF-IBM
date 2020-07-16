/*
 * actm.h
 *
 *  Created on: Sep 27, 2013
 *      Author: daan
 */

#ifndef ACTM_H_
#define ACTM_H_
#include "DataTypes.h"

void actm(Collision *pc, Collision_bag *bag, Parameters *params);

void actm_new_pp_collision(double *kn_ptr, double *dn_ptr, Collision_bag *bag,
		Parameters *params);
void actm_new_pw_collision(double *kn_ptr, double *dn_ptr, Collision_bag *bag,
		Parameters *params);
void actm_calc_coeff(double *kn, double *dn, double M,  double gn,
		double Tc, double e);

#define NUM_ACTM_STEPS 10.0 // number of timesteps per collision

#endif /* ACTM_H_ */
