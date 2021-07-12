#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Velocity.h"
#include "Grid.h"
#include "Display.h"
#include "Outflow.h"
#include "Inflow.h"
#include "Memory.h"
#include "Input.h"
#include "Communication.h"
#include "MyMath.h"
#include "Cart3d.h"
#include "Dtime.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern MPI_Comm comm3d;

/******************************************************************************/
/*
 This function computes dt from the CFL condition
 */
/******************************************************************************/
double Dtime_cfl(Cart3d_bag *data_bag) {

	int i, j, k;
	double T1, T2;
	double u_cfl, v_cfl, w_cfl, v_particle, u_cell, v_cell, w_cell;
	double nu_cell, convective, viscous;

	double dt, dt_cfl;
#ifdef LAG_PARTICLE_RESOLVED
	// Timestep required to resolve collision stiffness
	double dt_kn;
#endif
	double W_u_max;
	double send_data[2], W_min_data[2];
	char message[500];
	FILE *fid;

	T1 = MPI_Wtime();

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	double dt_old = params -> dt_old;
	double settling_speed_max = 0.0;
	double max_idt = 0.0;

	double u_max = 0.0;
	double v_max = 0.0;
	double w_max = 0.0;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Start index of bottom-left-back corner on current processor
	int Is = grid->G_Is;
	int Js = grid->G_Js;
	int Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid->G_Ie;
	int Je = grid->G_Je;
	int Ke = grid->G_Ke;

	// get the array for Velocity at cell center
	double ***u_data_bc = data_bag -> u -> data_bc;
	double ***v_data_bc = data_bag -> v -> data_bc;
	double ***w_data_bc = data_bag -> w -> data_bc;

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;

	double iRe = 1.0 / params->Re;
//#ifdef LES
//	double ***nut = data_bag->smag->nut;
//#elif defined RANS
//	double ***nut = data_bag->rans->nut;
//#endif
#ifdef VAR_VISC
	double ***nu = data_bag->viscosity->nu;
#endif

	// Find the maximum particle settling speed
#ifdef CONC
	int iconc;
	int NConc = params -> NConc;
	for (iconc = 0; iconc < NConc; iconc++) {
		if ( fabs(params->V_s0[iconc]) > settling_speed_max ) {
			settling_speed_max = fabs(params->V_s0[iconc]);
		} // if
	} // for iconc
#endif

	// Default value for constant viscosity case
	nu_cell = iRe;

	// Default value for fully-implicit case
	viscous = 0.0;

	//Calculate dt based on the CFL number
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {

				u_cfl      = fabs(u_data_bc[k][j][i]);
				w_cfl      = fabs(w_data_bc[k][j][i]);
				v_cell     = fabs(v_data_bc[k][j][i]);
				v_particle = fabs(v_data_bc[k][j][i] - settling_speed_max);
				v_cfl = max(v_cell, v_particle);

				convective = u_cfl * idx_u[i] + v_cfl * idy_v[j] + w_cfl * idz_w[k];

#ifdef VAR_VISC
				nu_cell = nu[k][j][i];
#endif

#ifdef FULLY_EXPLICIT

				viscous = 2.0 * nu_cell * ( idx_u[i] * idx_u[i] +
				                            idy_v[j] * idy_v[j] +
				                            idz_w[k] * idz_w[k] );
#elif !defined FULLY_IMPLICIT
				viscous = 2.0 * nu_cell * ( idx_u[i] * idx_u[i] +
				                            idz_w[k] * idz_w[k] );
#endif

				max_idt = max(max_idt, convective + viscous);

#if defined LEFT_OUTFLOW || defined RIGHT_OUTFLOW
				// To find maximum u-Velocity in the domain
				u_max = max(u_cfl, u_max);
				v_max = max(v_cell, v_max);
				w_max = max(w_cfl, w_max);
#endif

			} // for i
		} // for j
	} // for k

	dt_cfl = params->cfl / max_idt;

	// Now, find minimum dt and send it back to all processors
	send_data[0] = dt_cfl;
	send_data[1] = max(max(u_max,v_max), w_max);
	send_data[1] = 1.0/send_data[1]; // inverse u_max to send it with "dt"

	MPI_Allreduce ( (void *)send_data, (void *)W_min_data, 2, MPI_DOUBLE, MPI_MIN, PCW);

	// Global minimum dt
	dt_cfl = W_min_data[0];

	// Update the velocity for the convective outflow boundary condition
#if defined LEFT_OUTFLOW || defined RIGHT_OUTFLOW
	// Global maximum u-velocity from left to right
	W_u_max = 1.0 / W_min_data[1];
	params->U_conv_outflow = W_u_max;
#endif

#ifdef LAG_PARTICLE_RESOLVED
	// Timestep required to resolve particle collision stiffness
	dt_kn = Dtime_lag_particle(data_bag);
#endif

	// Constant time step
	if (params->constant_dt == 1) {
		if (params->default_dt > dt_cfl) {
			sprintf(message, "WARNING: Unstable timestep dt = %g for CFL = %g\n"
					"Stability condition requires dt = %g\n"
					"Continuing anyway...", params->default_dt, params->cfl, dt_cfl);
			Display_throw_warning(message, params);
		}
#ifdef LAG_PARTICLE_RESOLVED
		else if (params->default_dt > dt_kn) {
			sprintf(message, "WARNING: Unstable timestep dt = %g\n"
					"Stiffness stability condition requires dt = %g\n"
					"Continuing anyway...", params->default_dt, dt_kn);
			Display_throw_warning(message, params);
		}
#endif
		dt = params -> default_dt;
	}
	else {
#ifdef LAG_PARTICLE_RESOLVED
		// Set dt to staisfy CFL condition and collision stability
		dt = min(dt_cfl, dt_kn);
#else
		dt = dt_cfl;
#endif

		// This is done to avoid sharp jump in time steps.

		if ( dt > 1.05 * dt_old )
			dt = 1.05 * dt_old;
		if ( dt < 0.95 * dt_old )
			dt = 0.95 * dt_old;

		else


		// Compare to maximum allowable dt
		dt = min(dt, params->max_dt);

		if (params->rank == 0) {
			if (params->time == 0)
				fid = fopen("dthistory.dat", "w");
			else
				fid = fopen("dthistory.dat", "a");
			fprintf(fid, "%+-20.16e %+-20.16e \n", params->time, dt);
			fclose(fid);
		}
	}

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_cfl += T2 - T1;

	return dt;
}

/******************************************************************************/
/*
 This function computes the timestep required to maintain collision stability.
 This is determined by the [currently] hard-coded value 'N_dt_resolve', which is
 the number of timesteps we believe is required to resolve the collision in a
 stable manner.

 Substeps are taken into account.  That is, if substepping is enabled, it would
 allow the fluid timesteps to be 15 times larger.
 */
/******************************************************************************/
double Dtime_lag_particle(Cart3d_bag *data_bag) {

	char message[500];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Grid spacing (assuming uniform grid)
	double h = grid->dx_u[1];

	double dt, dt_old;
	dt_old = params->dt;

	// Maximum collision stiffness for a particle
	double kn_max;
	Collision *pc;

	// Number of timesteps required for stable resolution of collision forces
	double N_dt_resolve = 5.0;

	// Timestep required to resolve collision stiffness
	double dt_kn, dt_kn_min;
	dt_kn_min = params->max_dt;

	//--------------------------------------------------------------------------
	// The following calculations are based on the collision time for elastic
	// spring collisions (no damping):
	//     DEM:  Tc = PI * sqrt(M / kn)
	//     ACTM: sqrt(Tc^5 * u_in) = 18.578 * M / kn
	// For the ACTM spring, we set:
	//     u_in = CFL * h / dt
	// We then want to resolve the collision time with 'N_dt_resolve' timesteps:
	//     Tc = N_dt_resolve * dt
	// Or, if we are using substeps:
	//     Tc = N_dt_resolve * dt / 15
	// Solving for 'dt', these expressions work out to:
	//     dt = factor * sqrt(M / kn)
	// Where 'factor' depends on the collision model used and substepping
	//--------------------------------------------------------------------------
#ifdef ACTM
	#ifdef SUBSTEP
	double factor = sqrt(18.578 / sqrt(params->cfl * h * pow(N_dt_resolve/15.0, 5.0)));
	#else
	double factor = sqrt(18.578 / sqrt(params->cfl * h * pow(N_dt_resolve, 5.0)));
	#endif
#elif defined DEM
	#ifdef SUBSTEP
	double factor = PI / N_dt_resolve;
	#else
	double factor = 15.0 * PI / N_dt_resolve;
	#endif
#endif

	Particle *p = data_bag->lag->p_mobile_list->start;
	while (p != NULL) {

		// Maximum particle collision stiffness
		kn_max = 0.0;
		pc = p->particle_collision;
		while (pc != NULL) {
			kn_max = max(kn_max, pc->kn);
			pc = pc->next;
		}
		pc = p->wall_collision;
		while (pc != NULL) {
			kn_max = max(kn_max, pc->kn);
			pc = pc->next;
		}
		// Stable timestep based on stiffness
		if (kn_max > 0) {
			dt_kn = factor * sqrt(p->M / kn_max);
			dt_kn_min = min(dt_kn, dt_kn_min);
		}

		p = p -> next;
	}

	MPI_Allreduce(MPI_IN_PLACE, &dt_kn_min, 1, MPI_DOUBLE, MPI_MIN, PCW);

	return dt_kn_min;

}


/******************************************************************************/
