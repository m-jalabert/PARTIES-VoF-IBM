
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "DataTypes.h"

#define EPS 1e-14
#undef VERBOSE

/******************************************************************************/
/*
 */
/******************************************************************************/
void atfm_particle(Collision *pc, Collision *pc2, Collision_bag *bag, Parameters *params) {

	int i;
	double tan_force[3];
	Particle *p  = bag -> p;
	Particle *p2 = bag -> p2;
	double *n    = bag -> n;
	double *t    = bag -> t;
	double R_cp  = bag -> R_cp;
	double R2_cp = bag -> R2_cp;
	double tan_force_stat;

	const double BET[] = {BETA};
	double dt_bet = params -> dt * BET[params -> which_stage];

	//--------------------------------------------------------------------------
	// Sliding friction
	//--------------------------------------------------------------------------
	if (pc->type == SLIDING) {

		double friction_force = -params->mu_k * bag->normal_force_norm;

		tan_force[0] = friction_force * t[0]; // force on the particle
		tan_force[1] = friction_force * t[1]; // force on the particle
		tan_force[2] = friction_force * t[2]; // force on the particle

#ifdef VERBOSE
		printf("----------- %d <-> %d sliding ------------\n", p->ID, p2->ID);
		printf("t = [%e, %e, %e]\n",t[0],t[1],t[2]);
		printf("n = [%e, %e, %e]\n",bag->n[0],bag->n[1],bag->n[2]);
		printf("friction_force = %e\n",friction_force);
#endif

		//----------------------------------------------------------------------
		// Store forcing in collision structure
		//----------------------------------------------------------------------
		FORI3 p ->Fc[i] += tan_force[i];
		FORI3 p2->Fc[i] -= tan_force[i];

		p->Tc[0] += R_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
		p->Tc[1] += R_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
		p->Tc[2] += R_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);

		p2->Tc[0] += R2_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
		p2->Tc[1] += R2_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
		p2->Tc[2] += R2_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);

	}
	//--------------------------------------------------------------------------
	// Rolling friction
	//--------------------------------------------------------------------------
	else {

		double Mq = 1.0 / p->M + R_cp * R_cp / p->I_p;
		if (bag->stage != COLL_STAGE_FIXED) {
			Mq += 1.0 / p2->M + R2_cp * R2_cp / p2->I_p;
		}
		Mq = 1.0 / Mq;

		double F[3], T[3];
		double tan_force_norm, tan_force_stat;

		F[0] = ( 2.0 * p->F[0] + p->Fc_old[0] ) / p->M;
		F[1] = ( 2.0 * p->F[1] + p->Fc_old[1] ) / p->M;
		F[2] = ( 2.0 * p->F[2] + p->Fc_old[2] ) / p->M;

		T[0] = ( 2.0 * p->T[0] + p->Tc_old[0] ) * R_cp / p->I_p;
		T[1] = ( 2.0 * p->T[1] + p->Tc_old[1] ) * R_cp / p->I_p;
		T[2] = ( 2.0 * p->T[2] + p->Tc_old[2] ) * R_cp / p->I_p;

		F[0] += -( 2.0 * p2->F[0] + p2->Fc_old[0] ) / p2->M;
		F[1] += -( 2.0 * p2->F[1] + p2->Fc_old[1] ) / p2->M;
		F[2] += -( 2.0 * p2->F[2] + p2->Fc_old[2] ) / p2->M;

		T[0] += ( 2.0 * p2->T[0] + p2->Tc_old[0] ) * R2_cp / p2->I_p;
		T[1] += ( 2.0 * p2->T[1] + p2->Tc_old[1] ) * R2_cp / p2->I_p;
		T[2] += ( 2.0 * p2->T[2] + p2->Tc_old[2] ) * R2_cp / p2->I_p;

		double F_dot_t        = DOT(F, t);
		double Tn_cross_t     = T[0] * ( n[1] * t[2] - n[2] * t[1] );
		       Tn_cross_t    += T[1] * ( n[2] * t[0] - n[0] * t[2] );
		       Tn_cross_t    += T[2] * ( n[0] * t[1] - n[1] * t[0] );

		tan_force_norm = -Mq * ( bag->gt_cp_norm_old / dt_bet + F_dot_t + Tn_cross_t );

		//======================================================================

		// Limiter for sliding based on static friction
		tan_force_stat = -params->mu_s * bag->normal_force_norm;
		if (fabs(tan_force_norm) > fabs(tan_force_stat)) tan_force_norm = tan_force_stat;

		tan_force[0] = tan_force_norm * t[0];
		tan_force[1] = tan_force_norm * t[1];
		tan_force[2] = tan_force_norm * t[2];

#ifdef VERBOSE
		printf("--------- %d <-> %d rolling ----------\n", p->ID, p2->ID);
		printf("t = [%e, %e, %e]\n",t[0],t[1],t[2]);
		printf("n = [%e, %e, %e]\n",bag->n[0],bag->n[1],bag->n[2]);
		printf("Ft = [%e, %e, %e]\n",tan_force[0],tan_force[1],tan_force[2]);
#endif

	//----------------------------------------------------------------------
	// Apply forcing directly to particle
	//----------------------------------------------------------------------
	// Apply tangential force to particles
	FORI3 p ->Fc[i] += tan_force[i];
	FORI3 p2->Fc[i] -= tan_force[i];

	// Apply torque from tangential force to particles
	p->Tc[0] += R_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
	p->Tc[1] += R_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
	p->Tc[2] += R_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);

	p2->Tc[0] += R2_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
	p2->Tc[1] += R2_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
	p2->Tc[2] += R2_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void atfm_wall(Collision *pc, Collision_bag *bag, Parameters *params) {

	int i;
	double tan_force[3];
	Particle *p  = bag -> p;
	double *n    = bag -> n;
	double *t    = bag -> t;
	double R_cp  = bag -> R_cp;
	const double BET[] = {BETA};
	double dt_bet = params -> dt * BET[params -> which_stage];

	//--------------------------------------------------------------------------
	// Sliding friction
	//--------------------------------------------------------------------------
	if(pc->type == SLIDING) {

		double friction_force = -params->mu_k * bag->normal_force_norm;

		tan_force[0] = friction_force * t[0]; // force on the particle
		tan_force[1] = friction_force * t[1]; // force on the particle
		tan_force[2] = friction_force * t[2]; // force on the particle

#ifdef VERBOSE
		printf("----------- %d <-> wall sliding ------------\n", p->ID);
		printf("t = [%e, %e, %e]\n",t[0],t[1],t[2]);
		printf("n = [%e, %e, %e]\n",bag->n[0],bag->n[1],bag->n[2]);
		printf("friction_force = %e\n",friction_force);
#endif

		//----------------------------------------------------------------------
		// Store forcing in collision structure
		//----------------------------------------------------------------------
		FORI3 p->Fc[i] += tan_force[i];

		p->Tc[0] += R_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
		p->Tc[1] += R_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
		p->Tc[2] += R_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);

	}
	//--------------------------------------------------------------------------
	// Rolling friction
	//--------------------------------------------------------------------------
	else {

		double Mw = 1.0 / p->M + R_cp * R_cp / p->I_p;
		Mw = 1.0 / Mw;

		double F[3], T[3];
		double tan_force_norm, tan_force_stat;

		F[0] = ( 2.0 * p->F[0] + p->Fc_old[0] ) / p->M;
		F[1] = ( 2.0 * p->F[1] + p->Fc_old[1] ) / p->M;
		F[2] = ( 2.0 * p->F[2] + p->Fc_old[2] ) / p->M;

		T[0] = ( 2.0 * p->T[0] + p->Tc_old[0] ) * R_cp / p->I_p;
		T[1] = ( 2.0 * p->T[1] + p->Tc_old[1] ) * R_cp / p->I_p;
		T[2] = ( 2.0 * p->T[2] + p->Tc_old[2] ) * R_cp / p->I_p;

		double F_dot_t        = DOT(F, t);
		double Tn_cross_t     = T[0] * ( n[1] * t[2] - n[2] * t[1] );
		       Tn_cross_t    += T[1] * ( n[2] * t[0] - n[0] * t[2] );
		       Tn_cross_t    += T[2] * ( n[0] * t[1] - n[1] * t[0] );
		double fluid_corr     = -(F_dot_t * p->I_p + p->M * R_cp * Tn_cross_t);
		       fluid_corr     = fluid_corr / (p->I_p + p->M * R_cp * R_cp);

		tan_force_norm = -Mw * ( bag->gt_cp_norm_old / dt_bet + F_dot_t + Tn_cross_t );

		//limiter for sliding based on static friction
		tan_force_stat = -params->mu_s * bag->normal_force_norm;
		if (fabs(tan_force_norm) > fabs(tan_force_stat)) tan_force_norm = tan_force_stat;

		tan_force[0] = tan_force_norm * t[0];
		tan_force[1] = tan_force_norm * t[1];
		tan_force[2] = tan_force_norm * t[2];

#ifdef VERBOSE
		printf("--------- %d <-> wall rolling ----------\n", p->ID);
		printf("gt = [%e, %e, %e]\n", bag->gt[0],bag->gt[1],bag->gt[2]);
		printf("gt_cp = [%e, %e, %e]\n", bag->gt_cp[0],bag->gt_cp[1],bag->gt_cp[2]);
		printf("t  = [%e, %e, %e]\n",t[0],t[1],t[2]);
		printf("n  = [%e, %e, %e]\n",bag->n[0],bag->n[1],bag->n[2]);
		printf("Ft = [%e, %e, %e]\n",tan_force[0],tan_force[1],tan_force[2]);
		printf("fluid_corr %e, F_dot_t %e, Tn_cross_t %e\n", fluid_corr, F_dot_t, Tn_cross_t);
		printf("F  = [%e, %e, %e]\n",p->F[0],p->F[1],p->F[2]);
		printf("T  = [%e, %e, %e]\n",p->T[0],p->T[1],p->T[2]);
#endif

		//----------------------------------------------------------------------
		// Apply forcing directly to particle
		//----------------------------------------------------------------------
		// Apply tangential force to particles
		FORI3 p ->Fc[i] += tan_force[i];

		// Apply torque from tangential force to particles
		p->Tc[0] += R_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
		p->Tc[1] += R_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
		p->Tc[2] += R_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);

	}
}




/******************************************************************************/
/*
 * Determines the tangential collision regime, whether it be sliding or rolling
 */
/******************************************************************************/
int atfm_categorize(Collision_bag *bag, Parameters *params, double u_in) {

	double Psi_in = bag->gt_cp_norm / (fabs(u_in) + EPS);

#ifdef VERBOSE
	printf("Tangential collision psi_in: %e, psi_in_rs: %e\n", Psi_in, params->Psi_crit);
#endif

	// Check which regime we are in
	if (Psi_in >= params->Psi_crit) {
		return SLIDING;
	}
	else {
		return ROLLING;
	}

}
#undef VERBOSE
