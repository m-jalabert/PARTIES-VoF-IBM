#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "DataTypes.h"

#undef  VERBOSE
#define EPS 1e-12

/******************************************************************************/
/*
 * Evaluate tangential force 'Ft' based on linear spring, according to [Luding,
 * Granular Matter 2008]:
 *
 *    Ft = -kt * zeta_t - gt * gt_cp
 *    Ft = min( Ft, mu * Fn )
 *
 * where 'kt' and 'dt' are the tangential stiffness and damping coefficients,
 * 'mu' is the coefficient of friction, and 'Fn' is the normal force.  'zeta_t'
 * is the tangential displacement:
 *
 * zeta_t^n = zeta_t^{n-1} + gt_cp * dt_sub
 *
 * During sliding, however, this displacement must be reset to a distance
 * prescribed by the Coulomb force, as shown in the code below.
 *
 * We employ an additional technique suggested by Luding: rotating the
 * displacement 'zeta_t' to be parallel to the tangential vector 't' to prevent
 * the tangential spring from contributing to the normal force as particles
 * slide around each other.
 */
/******************************************************************************/
void lin_tan(Collision *pc, Collision *pc2, Collision_bag *bag, Parameters *params) {

	int i;

	Particle *p   = bag -> p;
	Particle *p2  = bag -> p2;
	double *n     = bag -> n;
	double  R_cp  = bag -> R_cp;
	double  R2_cp = bag -> R2_cp;
	double *gt_cp = bag -> gt_cp;

	double kt = pc -> kt;
	double dt = pc -> dt;

	double tan_force[3];
	double *zeta_t = pc -> zeta_t;
	double *zeta_t_old = pc -> zeta_t_old;

	const double BET[] = {BETA};
	double dt_sub = 2.0 * params -> dt * BET[params -> which_stage];
#ifdef SUBSTEP_OFF
	dt_sub = 15.0 * dt_sub;
#endif
	double temp;

// 	//--------------------------------------------------------------------------
// 	// The stiffness coefficient 'kt' is chosen to relate the half period of
//     // oscillation to the collision time 'Tc' [Shäfer, Dippel, and Wolf,
// 	// Journal de Physique I 1996]:
// 	//     kt = kappa * m * ( PI / Tc )^2
// 	// where 'm' is an effective mass and 'kappa' is a value based on Poisson's
// 	// ratio 'nu' [Thornton, Cummins, and Cleary, Powder Technology 2011]:
// 	//     kappa = 2 * (1 - nu) / (2 - nu)
// 	//--------------------------------------------------------------------------
// 	double Tc, m, kappa, kt;
// 	Tc = params->Ndt_coll * params->dt;
// #ifdef SUBSTEP
// 	Tc = Tc * 15.0;
// #endif
// 	if (bag->stage == COLL_STAGE_WALL || bag->stage == COLL_STAGE_FIXED)  m = p->M;
// 	else  m = (p->M * p2->M) / (p->M + p2->M);
// 	kappa = 2.0 * (1.0 - params->PoissonsRatio) / (2.0 - params->PoissonsRatio);
// //	kappa = 2.0 / 7.0;
// 	kt = PI * PI * m / ( Tc * Tc ) * kappa;
//
// 	//--------------------------------------------------------------------------
// 	// The damping coefficient 'dt' is chosen according to [Thornton, Cummins,
//     // and Cleary, Powder Technology 2013]:
// 	//     dt = 2 * dn * sqrt( m * kt )
// 	// where 'dn' here is the damping coefficient for the linear-spring normal
// 	// contact model:
// 	//     dn = -log(e) / sqrt( PI^2 + log(e)^2 )
// 	// where 'e' is the normal coefficient of restitution.
// 	//--------------------------------------------------------------------------
// 	double dt, e;
// 	if (bag->stage == COLL_STAGE_WALL) e = params->e_dry_wall;
// 	else e = params->e_dry_particles;
// 	temp = log(e);
// 	temp = -temp / sqrt(PI * PI + temp * temp);
// 	dt = 2.0 * temp * sqrt( m * kt );

	//--------------------------------------------------------------------------
	// Tangential force based on tangential stiffness
	//--------------------------------------------------------------------------

	// Rotate previous displacement onto tangent plane
	double n_dot_Ft = DOT(n, zeta_t_old);
	temp = sqrt(DOT(zeta_t_old, zeta_t_old));
	zeta_t_old[0] = zeta_t_old[0] - n_dot_Ft * n[0];
	zeta_t_old[1] = zeta_t_old[1] - n_dot_Ft * n[1];
	zeta_t_old[2] = zeta_t_old[2] - n_dot_Ft * n[2];
	temp = temp / (sqrt(DOT(zeta_t_old, zeta_t_old)) + EPS);
	zeta_t_old[0] = temp * zeta_t_old[0];
	zeta_t_old[1] = temp * zeta_t_old[1];
	zeta_t_old[2] = temp * zeta_t_old[2];

	// Add new displacement
	zeta_t[0] = zeta_t_old[0] + dt_sub * gt_cp[0];
	zeta_t[1] = zeta_t_old[1] + dt_sub * gt_cp[1];
	zeta_t[2] = zeta_t_old[2] + dt_sub * gt_cp[2];

	// New method
	tan_force[0] = -kt * zeta_t[0] - dt * gt_cp[0];
	tan_force[1] = -kt * zeta_t[1] - dt * gt_cp[1];
	tan_force[2] = -kt * zeta_t[2] - dt * gt_cp[2];

	//--------------------------------------------------------------------------
	// Limit tangential force based on friction
	//--------------------------------------------------------------------------
	// Coefficient of friction
	double mu;
	// Kinetic friction coefficient for sliding
	if (pc->type == SLIDING) {
		mu = params -> mu_k;
	}
	// Static friction coefficent for rolling
	else {
		mu = params -> mu_s;
	}

	// Compare desired tangential force to Coulomb friction force
	double tan_force_norm = sqrt(DOT(tan_force, tan_force));
	double friction_force = mu * bag->normal_force_norm;

	//--------------------------------------------------------------------------
	// Sliding occurs if Coulomb force exceeded
	//--------------------------------------------------------------------------
	if (tan_force_norm > friction_force) {

		pc->type = SLIDING;

		temp = friction_force / tan_force_norm;
		tan_force[0] = temp * tan_force[0];
		tan_force[1] = temp * tan_force[1];
		tan_force[2] = temp * tan_force[2];

		// Update displacement
		zeta_t[0] = -( tan_force[0] + dt * gt_cp[0] ) / kt;
		zeta_t[1] = -( tan_force[1] + dt * gt_cp[1] ) / kt;
		zeta_t[2] = -( tan_force[2] + dt * gt_cp[2] ) / kt;

#ifdef VERBOSE
		if (p2 == NULL)
			printf("----------- %d <-> wall sliding ------------\n", p->ID);
		else
			printf("----------- %d <-> %d sliding ------------\n", p->ID, p2->ID);
#endif
	}  // if sliding
	//--------------------------------------------------------------------------
	// Rolling contact
	//--------------------------------------------------------------------------
	else {

		pc->type = ROLLING;

#ifdef VERBOSE
		if (p2 == NULL)
			printf("----------- %d <-> wall rolling ------------\n", p->ID);
		else
			printf("----------- %d <-> %d rolling ------------\n", p->ID, p2->ID);
#endif
	}
#ifdef VERBOSE
	printf("kt    = %e\n", kt);
	printf("n     = [%e, %e, %e]\n", n[0], n[1], n[2]);
	printf("gt_cp = [%e, %e, %e]\n", gt_cp[0], gt_cp[1], gt_cp[2]);
	printf("Ft    = [%e, %e, %e]\n", tan_force[0], tan_force[1], tan_force[2]);
#endif

	//--------------------------------------------------------------------------
	// Apply forcing directly to particle
	//--------------------------------------------------------------------------
	// Apply tangential force to particles
	FORI3 p->Fc[i] += tan_force[i];
#ifdef POST_PROCESS
	FORI3 p->Fc_tan[i] += tan_force[i];
#endif

	// Apply torque from tangential force to particles
	p->Tc[0] += R_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
	p->Tc[1] += R_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
	p->Tc[2] += R_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);

	// Apply force and torque to other particle (unless it is a wall)
	if (p2 != NULL) {
		FORI3 p2->Fc[i] -= tan_force[i];
#ifdef POST_PROCESS
		FORI3 p2->Fc_tan[i] -= tan_force[i];
#endif

		p2->Tc[0] += R2_cp * (n[1] * tan_force[2] - n[2] * tan_force[1]);
		p2->Tc[1] += R2_cp * (n[2] * tan_force[0] - n[0] * tan_force[2]);
		p2->Tc[2] += R2_cp * (n[0] * tan_force[1] - n[1] * tan_force[0]);
	}

	// Copy values to other collision structure
	if (pc2 != NULL) {
		pc2->type = pc->type;
		// Tangential displacement vector changes sign
		FORI3 pc2->zeta_t[i] = -zeta_t[i];
		FORI3 pc2->zeta_t_old[i] = -zeta_t_old[i];
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void lin_tan_new_collision(double *kt_ptr, double *dt_ptr, Collision_bag *bag,
		Parameters *params) {

	Particle  *p  = bag -> p;
	Particle  *p2 = bag -> p2;

	//--------------------------------------------------------------------------
	// The stiffness coefficient 'kt' is chosen to relate the half period of
	// oscillation to the collision time 'Tc' [Shäfer, Dippel, and Wolf,
	// Journal de Physique I 1996]:
	//     kt = kappa * m * ( PI / Tc )^2
	// where 'm' is an effective mass and 'kappa' is a value based on Poisson's
	// ratio 'nu' [Thornton, Cummins, and Cleary, Powder Technology 2011]:
	//     kappa = 2 * (1 - nu) / (2 - nu)
	//--------------------------------------------------------------------------
	double Tc, m, kappa, kt;
	Tc = params->Ndt_coll * params->dt;
#ifdef SUBSTEP
	Tc = Tc * 15.0;
#endif
	if (bag->stage == COLL_STAGE_WALL || bag->stage == COLL_STAGE_FIXED)  m = p->M;
	else  m = (p->M * p2->M) / (p->M + p2->M);
	kappa = 2.0 * (1.0 - params->PoissonsRatio) / (2.0 - params->PoissonsRatio);
//	kappa = 2.0 / 7.0;
	kt = PI * PI * m / ( Tc * Tc ) * kappa;

	//--------------------------------------------------------------------------
	// The damping coefficient 'dt' is chosen according to [Thornton, Cummins,
	// and Cleary, Powder Technology 2013]:
	//     dt = 2 * dn * sqrt( m * kt )
	// where 'dn' here is the damping coefficient for the linear-spring normal
	// contact model:
	//     dn = -log(e) / sqrt( PI^2 + log(e)^2 )
	// where 'e' is the normal coefficient of restitution.
	//--------------------------------------------------------------------------
	double dt, e, temp;
	if (bag->stage == COLL_STAGE_WALL) e = params->e_dry_wall;
	else e = params->e_dry_particles;
	temp = log(e);
	temp = -temp / sqrt(PI * PI + temp * temp);
	dt = 2.0 * temp * sqrt( m * kt );

	// Return values
	*kt_ptr = kt;
	*dt_ptr = dt;
}

#undef VERBOSE
