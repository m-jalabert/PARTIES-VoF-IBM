#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"

#include "Array.h"
#include "Conc.h"
#include "Communication.h"
#include "Display.h"
#include "Dtime.h"
#include "Grid.h"
#include "Immersed.h"
#include "Inflow.h"
#include "Interpolate.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "MyMath.h"
#include "Outflow.h"
#include "Output.h"
#include "Particle.h"
#include "ParticleOutput.h"
#include "Pressure.h"
#include "Rans.h"
#include "Subgrid.h"
#include "Statistics2d.h"
#include "Strain.h"
#include "Temporal.h"
#include "Two_eqn_rans.h"
#include "Velocity.h"
#include "post_processing.h"
#include "EPforcing.h"


#define TIMEFILE "timesteps.dat"

extern MPI_Comm comm3d;

/******************************************************************************/
/*
 This function integrates the u,v and conc transport equation in time up to
 time_max
 */
/******************************************************************************/
int Temporal_int_rk3(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double time, max_time, dt;
	double output_time_interval;
	double output_time_interval_2d;
	int ntime;

	double T1, T2;
	double Tstart, Tend;

	int rk;
	const double BET2[] = {BETA2};
	const int tsubsteps = 3;

	int int_stop, status;
	FILE *fid;
	FILE *timeFile;
	char message[500];

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
	int pnodes         = params -> ghost_nodes;
	Timer      *timer  = data_bag -> timer;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
	Pressure *p = data_bag -> p;

#ifdef CONC
	int iconc;
	int NConc = params -> NConc;
	Concentration **c = data_bag -> c;
#endif

#ifdef LES
	Subgrid *smag = data_bag -> smag;
	double ***nut = smag->nut;
#endif
#ifdef RANS
	Rans *rans = data_bag -> rans;
	double ***nut = rans->nut;
#endif

#ifdef OUTPUT2D
	Statistics2d *st2d = data_bag -> st2d;
#endif




	/*------------------------------------------------------------------------*/
	/*
	 Get the time properties
	 */
	/*------------------------------------------------------------------------*/
	max_time         = params -> time_max;
	output_time_interval = params -> output_time_interval;
	output_time_interval_2d = params -> output_time_interval_2d;

	//--------------------------------------------------------------------------
	// Resuming simulation
	//--------------------------------------------------------------------------
	if (params->resume) {


		dt     = params -> dt;

		time   = params -> time;
		ntime  = params -> ntime;

		// Do not re-output resume timestep unless resuming from 'stop'
		if (params->noutput > 0) {
			params -> noutput++;
			params -> noutput_2d++;
			params -> output_time += output_time_interval;
			params -> output_time_2d += output_time_interval_2d;

		}
		// Resuming from stopped simulation
		else {
			params->noutput = -params->noutput;
		}
	}
	//--------------------------------------------------------------------------
	// Start from t=0 and default values
	//--------------------------------------------------------------------------
	else {

		params -> dt_old = params -> default_dt;
		params -> dt     = params -> default_dt;
		dt = params -> dt;

		ntime = 1;
		time  = 0.0;
		params -> output_time = 0.0;
		params -> output_time_2d = 0.0;

		params -> time = time;
		params -> noutput = 0;
		params -> noutput_2d = 0;
		params -> ntime = ntime;

		if (params -> rank == 0) {
			timeFile = fopen(TIMEFILE, "w");
			fclose(timeFile);
		}
	} // else


	// Initial output to file. Corresponding to zero velocity.
	Velocity_cell_center(data_bag);

#ifdef VAR_VISC
	Viscosity_set_cell_edges(data_bag);
#endif
	Communication_update_ghost_nodes_flow_variable(p->p_data, CONCENTRATION_PERTURBATION, 1, data_bag);
#ifdef CONC

	for (iconc=0; iconc<NConc; iconc++) {
		Communication_update_ghost_nodes_flow_variable(c[iconc]->data, CONCENTRATION, pnodes, data_bag);

	}
#endif

	//--------------------------------------------------------------------------
	// Saving the initial conditions
	//--------------------------------------------------------------------------

	if (!params->resume) {
		Output_h5_data(data_bag, DTRACE("Output_h5_data"));
		Output_h5_resume(data_bag, DTRACE("Output_h5_resume"));
#ifdef LAG_PARTICLE_RESOLVED

		ParticleOutput_h5(data_bag, params->noutput, DTRACE("ParticleOutput_h5"));

#endif
		if (params -> rank == 0) {
			timeFile = fopen(TIMEFILE, "a");
			fprintf(timeFile, "%f\n", time);
			fclose(timeFile);
		}

#ifdef OUTPUT2D
		//Output_h5_data2d(data_bag, DTRACE("Output_h5_data2d"));
		Statistics2d_computeStatistics(data_bag, DTRACE("Statistics2d_computeStatistics"));
		params->noutput_2d++;

#endif
		params->output_time += output_time_interval;
		params->output_time_2d += output_time_interval_2d;
		params->noutput++;

		Display_progress(params,"Initial flow properties have been saved successfully...\n");
	} // if not resuming


	/*------------------------------------------------------------------------*/
	/*																		  */
	/*	 	 	 	 	 	 Time-stepping									  */
	/*																		  */
	/*------------------------------------------------------------------------*/

	Display_progress(params,"\n*** Beginning main simulation loop ***\n");
	int_stop = 0;


	while ( int_stop == 0 && time - max_time < dt ) {   ///// Here starts time-stepping

		Tstart = MPI_Wtime();
		params -> dt    = dt;
		params -> ntime = ntime;
		params -> time  = time;


		//----------------------------------------------------------------------
		// Save data to file during output time steps
		//----------------------------------------------------------------------

#if defined LAG_PARTICLE_RESOLVED// && !defined SUBSTEP

	#if defined PARTICLE_TRN
		sprintf(message, "Outputting par_trn %d\n", ntime);
		Display_progress(params, message);
		//test_2d_output(data_bag, ntime, DTRACE("test_2d_output"));
		ParticleOutput_h5(data_bag, ntime, DTRACE("ParticleOutput_h5"));
	#endif

		ParticleOutput_dat(data_bag->lag->p_mobile_list, grid, params, DTRACE("ParticleOutput_dat"));

	#endif


		// Regular monitoring  inside the following block

		if ( (fabs(time - params->output_time) < dt &&
		     fabs(time - params->output_time) < fabs(time + dt - params->output_time)) ||
		     (time > params->output_time) ) {


#ifdef CONC
		//----------------------------------------------------------------------
		// Output concentration data
		//----------------------------------------------------------------------
	//	T1 = MPI_Wtime();
		// Integrate (in time) deposited height (from particles

	//	Post_processing_conc(data_bag);

	//	T2 = MPI_Wtime();
	//	timer->Wtime_output += T2 - T1;
#endif

		Output_h5_data(data_bag, DTRACE("Output_h5_data"));
		Output_h5_resume(data_bag, DTRACE("Output_h5_resume"));


#ifdef LAG_PARTICLE_RESOLVED

				ParticleOutput_h5(data_bag, params->noutput, DTRACE("ParticleOutput_h5"));

#endif
			if (params -> rank == 0) {
				timeFile = fopen(TIMEFILE, "a");
				fprintf(timeFile, "%f\n", time);
				fclose(timeFile);
			}

			params->output_time += output_time_interval;
			params->noutput++;
			Display_throw_warning("Runtime data has been saved successfully", params);

#if defined XPERIODIC && defined OUTPUT2D
			xzperiodic_uvel_ave(data_bag,u, v, w, grid, params);
//			xzperiodic_saltsediment_ave(data_bag, c, grid, params);
	#ifdef THERMAL_KADER
			xzperiodic_thermal_ave(data_bag, c, grid, params);
	#endif
	#ifdef RANS
			xzperiodic_rans_nut_ave(data_bag, data_bag->rans, grid, params);
	#endif
#endif
		} // Output writing


#ifdef OUTPUT2D
		if ( (fabs(time - params->output_time_2d) < dt &&
		     fabs(time - params->output_time_2d) < fabs(time + dt - params->output_time_2d)) ||
		     (time > params->output_time_2d) ) {


				Statistics2d_computeStatistics(data_bag, DTRACE("Statistics2d_computeStatistics"));

				params->output_time_2d += output_time_interval_2d;
				params->noutput_2d++;

		}
#endif
		/*--------------------------------------------------------------------*/
		/*
		 Runge Kutta sub-steps
		 */
		/*--------------------------------------------------------------------*/
#ifdef RANS
		//----------------------------------------------------------------------
		// RANS equation integration
		//----------------------------------------------------------------------
	#ifdef TWO_EQUATION_MODEL11
		int nrsteps = 4;
		for (rk = 0; rk < tsubsteps; rk++){

			params -> which_stage = rk;

			for (rsteps = 0; rsteps < nrsteps; rsteps++) {
				Rans_int_equations(data_bag, dt/nrsteps, DTRACE("Rans_int_equations"));
			}

		} // End RK substepping
	#endif
#endif

		//----------------------------------------------------------------------
		// Navier-Stokes and other equations integration
		//----------------------------------------------------------------------

		T1 = MPI_Wtime();
		for (rk = 0; rk < tsubsteps; rk++){
			if (time != 0) {

#if defined LEFT_INFLOW || defined RIGHT_INFLOW
				Inflow_velocity_profile(data_bag, DTRACE("Inflow_velocity_profile"));
#endif

#if defined LEFT_OUTFLOW || defined RIGHT_OUTFLOW
				//--------------------------------------------------------------
				// Impose convective boundary condition at the outlet
				//--------------------------------------------------------------
				Outflow_impose_convective_boundary(data_bag);
#endif

				if (rk == 0){
					//----------------------------------------------------------
					// Cell centered velocity: Communication is done for local
					// ghost nodes within the function
					//----------------------------------------------------------
					Velocity_cell_center(data_bag); // is this necessary??
					//----------------------------------------------------------
					// Compute dt from the CFL condition; find convective
					// outflow velocity
					//----------------------------------------------------------
					params -> dt_old = dt;
					dt = Dtime_cfl(data_bag);
					params -> dt = dt;
				}
			} // if (time != 0)

			params -> which_stage = rk;

#ifdef RANS

			Rans_int_equations(data_bag, dt, DTRACE("Rans_int_equations"));

#endif


			Temporal_int_all_the_equations(data_bag, DTRACE("Temporal_int_all_the_equations"));


			time += dt * BET2[rk];
			params -> time = time;

		} // End RK substepping
		T2 = MPI_Wtime();
		timer->Wtime_intEOM += T2 - T1;

#if defined CONSTANT_MASSFLUX || defined FLUID_OSCILLATION || defined PARTICLE_OSCILLATION
		Velocity_calculate_dpdx(u, data_bag);
#endif
#ifdef LES
		Strain_rate_magnitude(u, v, w, grid, params, smag->st_rate);
	#ifdef SMAG_DYNAMIC
		Dynamic_smag_coeff(data_bag);
	#endif
#endif

		//----------------------------------------------------------------------
		// Update time
		//----------------------------------------------------------------------
		ntime++;
		params -> ntime  = ntime;


		//----------------------------------------------------------------------
		// Iteration display
		//----------------------------------------------------------------------
		sprintf(message, "Iteration %d, dt = %.4g, time = %.4g", ntime, dt, time);
		Display_throw_warning(message, params);

		//----------------------------------------------------------------------
		// Read stop file
		//----------------------------------------------------------------------
		status = 0;
		if (params->rank == 0) {
			fid = fopen("stop.inp","r");
			if (fid == NULL) {
				status = -1;
			}
			else {
				status = fscanf(fid,"%d",&int_stop);
				if (status < 1) {
					status = -2;
				}
				fclose(fid);
			}
		}
		sprintf(message, "Could not read 'stop.inp'");
		Display_assert_error(status, message, params, DTRACE("Display_assert_error"));
		MPI_Bcast(&int_stop, 1, MPI_INT, 0, PCW);

		//----------------------------------------------------------------------
		// Save data, if stopping
		//----------------------------------------------------------------------
		if (int_stop) {
#ifdef OUTPUT2D
			//Output_h5_data2d(data_bag, DTRACE("Output_h5_data2d"));
			Statistics2d_computeStatistics(data_bag, DTRACE("Statistics2d_computeStatistics"));
			params->output_time_2d += output_time_interval_2d;
#endif

			Output_h5_data(data_bag, DTRACE("Output_h5_data"));
			Output_h5_resume(data_bag, DTRACE("Output_h5_resume"));
#if defined LAG_PARTICLE_RESOLVED
			ParticleOutput_h5(data_bag, params->noutput, DTRACE("ParticleOutput_h5"));
	#ifdef PARTICLE_TRN
			ParticleOutput_h5(data_bag, ntime, DTRACE("ParticleOutput_h5"));
	#endif
#endif


			if (params -> rank == 0) {
				timeFile = fopen(TIMEFILE, "a");
				fprintf(timeFile, "%f\n", time);
				fclose(timeFile);
			}

		}

		Tend = MPI_Wtime();
		timer->Wtime_total += Tend - Tstart;

	} // main time-stepping while


	MPI_Barrier(PCW);
	fflush(stdout);

	return ntime;
}




/******************************************************************************/
/*
 This function integrates u,v and conc transport equation for one time step
 */
/******************************************************************************/
void Temporal_int_all_the_equations(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i;
	double T1, T2, Tstart, Tend;
	double G_div_max;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Timer      *timer  = data_bag -> timer;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
	Pressure *p = data_bag -> p;
#ifdef TURB_FORCING
	Fourier *fourier = data_bag -> fourier;
#endif
#ifdef CONC
	Concentration **c = data_bag -> c;
#endif
#ifdef LES
	Subgrid *smag = data_bag -> smag;
#endif
#ifdef RANS
	Rans *rans = data_bag -> rans;
#endif
#ifdef LAG_PARTICLE_RESOLVED
	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;
	#ifdef PARTICLE_RELEASE
	Particle_list *p_release_list  = data_bag -> lag -> p_release_list;

	//printf("Releasing particles\n");
	Particle_release_to_mobile(p_release_list, p_mobile_list, params, DTRACE("Particle_release_to_mobile"));
	//printf("Turning off particles\n");
	//Particle_mobile_turn_off(p_mobile_list, params);
	#endif
#endif

	double dt = params -> dt;
	int which_stage = params -> which_stage;

	int pnodes = params -> ghost_nodes;


	/*------------------------------------------------------------------------*/
	/*
	 Cell centered velocity: Communication is done for local ghost nodes within
	 the function
	 If which_stage == 0, velocity at cell center is calculated before the
	 to Dtime_cfl
	 */
	/*------------------------------------------------------------------------*/
	if (which_stage > 0) {

		Velocity_cell_center(data_bag);
	}


	/*------------------------------------------------------------------------*/
	/*
	 Communicate cell-centered velocity
	 */
	/*------------------------------------------------------------------------*/
	Communication_update_ghost_nodes_flow_variable(u->data_bc, CONCENTRATION_PERTURBATION, pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data_bc, CONCENTRATION_PERTURBATION, pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data_bc, CONCENTRATION_PERTURBATION, pnodes, data_bag);


#ifdef CONC
	int iconc;
	int NConc = params -> NConc;

	for (iconc=0; iconc<NConc; iconc++) {

		Communication_update_ghost_nodes_flow_variable(c[iconc]->data, CONCENTRATION, pnodes, data_bag);
		Conc_set_boundary_values(c[iconc]->data, iconc ,CENTRAL_FULL ,grid , params);
	}
#endif


	/*------------------------------------------------------------------------*/
	/*
	 LES stuff
	 */
	/*------------------------------------------------------------------------*/
#ifdef LES
	T1 = MPI_Wtime();
	if (which_stage != 0)
		Strain_rate_magnitude(u, v, w, grid, params, smag->st_rate);
	Subgrid_smagorinsky_eddy_viscosity(data_bag);
	#ifdef IMMERSED_BOUNDARY
	Immersed_nut_interpolation(grid, params, smag->nut);
	#endif
	#ifdef SCHUMANN
	Schumann_calc_shear(u->data, v->data, w->data, smag->log_law, grid, params);
	#endif
	T2 = MPI_Wtime();
	timer->Wtime_sgs_total += T2 - T1;
#endif

#ifdef VAR_VISC
	// TODO: update 'nut' boundaries
	Viscosity_set_cell_edges(data_bag);
#endif


	/*------------------------------------------------------------------------*/
	/*
	 Calculate the explicitly and semi-implicitly treated part of the RHS for
	 momentum equation
	 */
	/*------------------------------------------------------------------------*/
	T1 = MPI_Wtime();
	Velocity_u_set_implicit_explicit(data_bag);
	Velocity_v_set_implicit_explicit(data_bag);
	Velocity_w_set_implicit_explicit(data_bag);
	T2 = MPI_Wtime();
	timer->Wtime_vel_convective += T2 - T1;


	/*------------------------------------------------------------------------*/
	/*
	 Solve Concentration equations of motion
	 */
	/*------------------------------------------------------------------------*/



#ifdef LAG_PARTICLE_RESOLVED
	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));
	Lagrangian_flag_points(data_bag, DTRACE("Lagrangian_flag_points"));  // New Langragian points are calculated


	#ifdef IBM_SCALAR
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfc);
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));


	#endif
#endif



#ifdef CONC
	Conc_int_equations(data_bag, DTRACE("Conc_int_equations"));
#endif


	/*------------------------------------------------------------------------*/
	/*
	 Set velocity RHS terms
	 */
	/*------------------------------------------------------------------------*/
	int NY = data_bag->grid->NY;


#ifdef TURB_FORCING
if(which_stage == 0){
	Display_progress(params,"Updating turbulent forcing\n");
	b_and_fourier_update(fourier, params);
	ftx_update(data_bag);
	fty_update(data_bag);
	ftz_update(data_bag);

	// if (params->rank == 0) {
	// 	static int nbfileee = 1;
	// 	char filename[50];
	// 	///////file write
	// 	sprintf(filename, "turbu%d.txt",nbfileee++);
	// 	FILE *fil = fopen(filename, "w");
	// 	if (fil == NULL)
	// 	{
	// 		printf("Error opening file!\n");
	// 		exit(1);
	// 	}
	// 	////////
	//
	// 	// Same for all quantities
	// 	int NX = grid -> NX;
	// 	int NY = grid -> NY;
	// 	int NZ = grid -> NZ;
	// 	// Start index of bottom-left-back corner on current processor
	// 	int Is = grid -> G_Is;
	// 	int Js = grid -> G_Js;
	// 	int Ks = grid -> G_Ks;
	// 	// End index of top-right-front corner on current processor
	// 	int Ie = grid -> G_Ie;
	// 	int Je = grid -> G_Je;
	// 	int Ke = grid -> G_Ke;
	// 	// Indices start and end on current processor
	// 	int i_start = max(1,Is); // i=0 not included
	// 	int j_start = max(1,Js);;
	// 	int k_start = max(1,Ks);;
	// 	// Exclude the half cell added
	// 	int i_end = min(NX-1, Ie);
	// 	int j_end = min(NY-1, Je);
	// 	int k_end = min(NZ-1, Ke);
	//
	// 	double ***fturb0 = data_bag -> u -> fturb;
	// 	double ***fturb1 = data_bag -> v -> fturb;
	// 	double ***fturb2 = data_bag -> w -> fturb;
	// 	//print ft in file
	// 	for (int i = i_start; i < i_end; i++) {
	// 		for (int j = j_start; j < j_end; j++) {
	// 			for (int k = k_start; k < k_end; k++) {
	// 				fprintf(fil, "%f %f %f\n", fturb0[k][j][i], fturb1[k][j][i], fturb2[k][j][i]);
	// 			}
	// 			fprintf(fil, "0  \n");
	// 		}
	// 		fprintf(fil, "0 0  \n");
	// 	}
	//
	// 	fclose(fil);  ////text file
	// }


}

#endif


	T1 = MPI_Wtime();
	Velocity_u_set_RHS(data_bag);
	/*if(params->ntime==1){
	test_2d_output(data_bag-> u -> ng_rhs,'p',data_bag, which_stage, DTRACE("test_2d_output"));
}*/
	Velocity_v_set_RHS(data_bag);
	/*if(params->ntime==1){
	test_2d_output(data_bag-> v -> ng_rhs,'p',data_bag, 10+which_stage, DTRACE("test_2d_output"));
}*/

	Velocity_w_set_RHS(data_bag);



#if defined BOUSSINESQ  && defined LAG_PARTICLE_RESOLVED

	if(  !(int)ceil(params->grav[0]) == 0){
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfu);
		Interpolate_add_to_volume_fraction('u', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('u', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	}
	if(  !(int)ceil(params->grav[1]) == 0){
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfv);
		Interpolate_add_to_volume_fraction('v', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('v', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	}
	if( !(int)ceil(params->grav[2]) == 0){
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfw);
		Interpolate_add_to_volume_fraction('w', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('w', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	}
#endif

#ifdef BOUSSINESQ

	 Velocity_add_buoyancy_2_RHS(data_bag); // we need the face centered volume fraction for the buoyancy in the case of particles

#endif


	T2 = MPI_Wtime();
	timer->Wtime_vel_rhs += T2 - T1;


#ifdef LAG_PARTICLE_RESOLVED
	Tstart = MPI_Wtime();
	//--------------------------------------------------------------------------
	// Find intermediate velocity field
	//--------------------------------------------------------------------------



	T1 = MPI_Wtime();
	Velocity_solve_explicit(u, data_bag);
	Velocity_solve_explicit(v, data_bag);
	Velocity_solve_explicit(w, data_bag);
	T2 = MPI_Wtime();
	timer->Wtime_vel_solve += T2 - T1;


	//--------------------------------------------------------------------------
	// Update boundary conditions
	//--------------------------------------------------------------------------
	Velocity_update_boundaries(u->data, 'u', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(v->data, 'v', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(w->data, 'w', VEL_TYPE_NORMAL, data_bag);

	//--------------------------------------------------------------------------
	// Forcing on RHS from particle immersed boundaries
	//--------------------------------------------------------------------------
	#ifdef POST_PROCESS
	if (params->which_stage == 0) {
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_fx_IBM);
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_fy_IBM);
	}
	#endif




		Lagrangian_force(-1, data_bag, DTRACE("Lagrangian_force"));
		Tend = MPI_Wtime();
		timer->Wtime_particle_total += Tend - Tstart;


#endif

/*
printf("vdata NY-3 is %2.5f\n",v->data[0][NY-3][0]);
printf("vdata NY-2 is %2.5f\n",v->data[0][NY-2][0]);
printf("vdata NY-1 is %2.5f\n",v->data[0][NY-1][0]);
printf("vdata 0 is %2.5f\n",v->data[0][0][0]);
printf("vdata 1 is %2.5f\n",v->data[0][1][0]);
printf("vdata 2 is %2.5f\n",v->data[0][2][0]);
*/
	/*------------------------------------------------------------------------*/
	/*
	 Solve for pre-projected velocity
	 */
	/*------------------------------------------------------------------------*/

/*	printf("vdata NY-2 is %2.5f\n",v->data[0][NY-2][0]);
	printf("vdata NY-1 is %2.5f\n",v->data[0][NY-1][0]);
	printf("vdata 0 is %2.5f\n",v->data[0][0][0]);
	printf("vdata 1 is %2.5f\n",v->data[0][1][0]);
	printf("aaaa vdata 2 is %2.5f\n",v->data[0][2][0]);
*/
/*if(params->ntime==1){
		test_2d_output(data_bag->u->data,'u',data_bag, 1, DTRACE("test_2d_output"));
		test_2d_output(data_bag->v->data,'v',data_bag, 1, DTRACE("test_2d_output"));
}*/

	T1 = MPI_Wtime();
	Velocity_solve(u, data_bag);
	Velocity_solve(v, data_bag);
	Velocity_solve(w, data_bag);
	T2 = MPI_Wtime();
	timer->Wtime_vel_solve += T2 - T1;

	/*if(params->ntime==1){
			test_2d_output(data_bag->u->data,'u',data_bag, 2, DTRACE("test_2d_output"));
			test_2d_output(data_bag->v->data,'v',data_bag, 2, DTRACE("test_2d_output"));
	}*/


	/*
	printf("vdata NY-3 is %2.5f\n",v->data[0][NY-3][0]);
	printf("vdata NY-2 is %2.5f\n",v->data[0][NY-2][0]);
	printf("vdata NY-1 is %2.5f\n",v->data[0][NY-1][0]);
	printf("vdata 0 is %2.5f\n",v->data[0][0][0]);
	printf("vdata 1 is %2.5f\n",v->data[0][1][0]);
	printf("oooo vdata 2 is %2.5f\n",v->data[0][2][0]);
*/

	//--------------------------------------------------------------------------
	// Update boundary conditions
	//--------------------------------------------------------------------------


 	Velocity_update_boundaries(u->data, 'u', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(v->data, 'v', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(w->data, 'w', VEL_TYPE_NORMAL, data_bag);

	/*if(params->ntime==1){
			test_2d_output(data_bag->u->data,'u',data_bag, 3, DTRACE("test_2d_output"));
			test_2d_output(data_bag->v->data,'v',data_bag, 3, DTRACE("test_2d_output"));
	}*/

#ifdef LAG_PARTICLE_RESOLVED
	//--------------------------------------------------------------------------
	// Correction to velocity field from particle immersed boundaries
	//--------------------------------------------------------------------------
	#ifndef ONE_WAY
 	Tstart = MPI_Wtime();
	Lagrangian_force(params->N_forcing_loops, data_bag, DTRACE("Lagrangian_force"));
	Tend = MPI_Wtime();
	timer->Wtime_particle_total += Tend - Tstart;

	//--------------------------------------------------------------------------
	// Update boundary conditions
	//--------------------------------------------------------------------------
 	Velocity_update_boundaries(u->data, 'u', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(v->data, 'v', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(w->data, 'w', VEL_TYPE_NORMAL, data_bag);
	#endif
#endif



	/*------------------------------------------------------------------------*/
	/*
	 Immersed boundary stuff
	 */
	/*------------------------------------------------------------------------*/
#ifdef IMMERSED_BOUNDARY
	Immersed_velocity_interpolation(u, grid, params);
	Immersed_velocity_interpolation(v, grid, params);
	Immersed_velocity_interpolation(w, grid, params);
#endif

	/*------------------------------------------------------------------------*/
	/*
	 Solve for pressure and project onto velocity field (to make divergence-
	 free)
	 */
	/*------------------------------------------------------------------------*/
	int poisson_iters = 0;
	do {

		// Prepare and solve for pressure correction
		T1 = MPI_Wtime();
		Pressure_set_RHS(data_bag);
		T2 = MPI_Wtime();


		/*
		if (params->ntime==1){
		test_2d_output(data_bag->u->data,'u',data_bag, nti, DTRACE("test_2d_output"));
		test_2d_output(data_bag->v->data,'v',data_bag, nti, DTRACE("test_2d_output"));
		test_2d_output(data_bag->w->data,'w',data_bag, nti, DTRACE("test_2d_output"));
		test_2d_output(data_bag->p->rhs,'p',data_bag, nti, DTRACE("test_2d_output"));

		nti++;
		}
		*/
		timer->Wtime_p_rhs += T2 - T1;

		T1 = MPI_Wtime();
		Pressure_solve(p, grid, params);
		//printf(" Iteration of Poisson solver: %d\n", poisson_iters);
		T2 = MPI_Wtime();
		timer->Wtime_p_solve += T2 - T1;

		//printf("vertical velocity at bottom is v[0][0][0]=%2.3f\n",v->data[0][0][0]);
		//printf("vertical velocity at top is v[0][NY-1][0]=%2.3f\n",v->data[0][NY-1][0]);
		// use pressure to project a divergence free velocity field

		T1 = MPI_Wtime();
		Pressure_project_velocity(data_bag);
		T2 = MPI_Wtime();
		timer->Wtime_p_project += T2 - T1;



		Velocity_update_boundaries(u->data, 'u', VEL_TYPE_NORMAL, data_bag);
		Velocity_update_boundaries(v->data, 'v', VEL_TYPE_NORMAL, data_bag);
		Velocity_update_boundaries(w->data, 'w', VEL_TYPE_NORMAL, data_bag);



		T1 = MPI_Wtime();
		G_div_max = Pressure_compute_velocity_divergence(data_bag);
		if (params->rank==0)
			printf(" Maximum Velocity Divergence: %2.2e at poisson iteration iteration %d\n", G_div_max, poisson_iters);
		T2 = MPI_Wtime();
		timer->Wtime_p_divergence += T2 - T1;

		poisson_iters++;
		if (poisson_iters > 10) {
			Display_throw_error("Poisson solver did not converge", params, DTRACE("Display_throw_error"));
		}

	} while (G_div_max > 1e-6);


#ifdef VOF_SCALAR
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfc);
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));


		Memory_reset_flow_variable(grid, params, data_bag->lag->vfu);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfv);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfw);

#ifdef VOF_SMOOTH_VELO
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfu_prime);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfv_prime);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfw_prime);
#endif
		Memory_reset_flow_variable(grid, params, data_bag->u->data_vof);
		Memory_reset_flow_variable(grid, params, data_bag->v->data_vof);
		Memory_reset_flow_variable(grid, params, data_bag->w->data_vof);


		Interpolate_add_to_volume_fraction_vof('u', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction_vof('u', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));

		Interpolate_add_to_volume_fraction_vof('v', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction_vof('v', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));

		Interpolate_add_to_volume_fraction_vof('w', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction_vof('w', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));

		Interpolate_bound_to_one( data_bag->lag->vfu, data_bag);
		Interpolate_bound_to_one( data_bag->lag->vfv,  data_bag);
		Interpolate_bound_to_one( data_bag->lag->vfw,  data_bag);

		Lagrangian_calc_compound_vof_fluid_velo(data_bag);

#ifdef VOF_PROJECT
			Velocity_update_boundaries(u->data_vof, 'u', VEL_TYPE_NORMAL, data_bag);
			Velocity_update_boundaries(v->data_vof, 'v', VEL_TYPE_NORMAL, data_bag);
			Velocity_update_boundaries(w->data_vof, 'w', VEL_TYPE_NORMAL, data_bag);

				Pressure_set_RHS_vof(data_bag);
				Pressure_solve(p, grid, params);
				Pressure_project_velocity_vof(data_bag);
#endif




		Velocity_update_boundaries(u->data_vof, 'u', VEL_TYPE_NORMAL, data_bag);
		Velocity_update_boundaries(v->data_vof, 'v', VEL_TYPE_NORMAL, data_bag);
		Velocity_update_boundaries(w->data_vof, 'w', VEL_TYPE_NORMAL, data_bag);
		Communication_update_ghost_nodes_flow_variable(data_bag->lag->vfu, U_VELOCITY_PERTURBATION, pnodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(data_bag->lag->vfv, V_VELOCITY, pnodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(data_bag->lag->vfw, W_VELOCITY, pnodes, data_bag);



#endif






#ifdef LAG_PARTICLE_RESOLVED
	Tstart = MPI_Wtime();
	//--------------------------------------------------------------------------
	// Integrate translational and angular velocities over particle domains
	//--------------------------------------------------------------------------
	T1 = MPI_Wtime();
	// Reset integrals to zero
	{
		Particle *p = p_mobile_list -> start;
		while (p != NULL) {
			DSET_ZERO(p->Int_U, 3);
			DSET_ZERO(p->Int_Omega, 3);
			p = p -> next;
		}
		p = p_fixed_list -> start;
		while (p != NULL) {
			DSET_ZERO(p->Int_U, 3);
			DSET_ZERO(p->Int_Omega, 3);
			p = p -> next;
		}
	}
//#ifndef ONE_WAY
	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfu);
	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfv);
	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfw);

	Interpolate_integrate_momentum(u, p_mobile_list, data_bag, DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(v, p_mobile_list, data_bag, DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(w, p_mobile_list, data_bag, DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(u, p_fixed_list, data_bag, DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(v, p_fixed_list, data_bag, DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(w, p_fixed_list, data_bag, DTRACE("Interpolate_integrate_momentum"));

	T2 = MPI_Wtime();
	timer->Wtime_particle_int += T2 - T1;
//#endif

	//--------------------------------------------------------------------------
	// Advect particles
	//--------------------------------------------------------------------------
	Lagrangian_advect_particles(data_bag, DTRACE("Lagrangian_advect_particles"));
	Tend = MPI_Wtime();
	timer->Wtime_particle_total += Tend - Tstart;
#endif  // LAG_PARTICLE_RESOLVED

#ifdef IMMERSED_BOUNDARY
//	Communication_update_ghost_nodes_flow_variable(u->data, 'u', pnodes, data_bag);
//	Communication_update_ghost_nodes_flow_variable(v->data, 'v', pnodes, data_bag);
//	Communication_update_ghost_nodes_flow_variable(w->data, 'w', pnodes, data_bag);
	Array_copy_withghost(u->data, u->data_old, grid, params);
	Array_copy_withghost(v->data, v->data_old, grid, params);
	Array_copy_withghost(w->data, w->data_old, grid, params);
#endif

}


/******************************************************************************/
/*
 * Simple function to output a slice of 'u_data'
 */
/******************************************************************************/
void test_2d_output(double ***data3d, char nme, Cart3d_bag *data_bag, int counter, Debug_trace *dtrace) {

	int i, j, k;
	char filename[50];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	int Is = grid->G_Is;
	int Js = grid->G_Js;
	int Ks = grid->G_Ks;

	int Ie = grid->G_Ie;
	int Je = grid->G_Je;

	// Dimensions of slice as part of a 3-D array
	int dim[3];
	dim[0] = grid->NX;
	dim[1] = grid->NY;
	dim[2] = 1;

	// Allocate storage for send an receive buffers (initialized to zero)
	// (in production code, only free at beginning of simulation)
	double **G_data2d = Memory_allocate_2D_double_array(dim[0], dim[1]);
	double **W_data2d = Memory_allocate_2D_double_array(dim[0], dim[1]);


	// Copy slice of 3-D data (from local processor)
	k = 0;
	for (j = Js; j < Je; j++) {
		for (i = Is; i < Ie; i++) {
			G_data2d[j][i] = data3d[k][j][i];
		}
	}

	// Describe dimensions of 2-D data
	Indices G_s, G_e, W_e;
	G_s.x_index = Is;
	G_s.y_index = Js;

	G_e.x_index = Ie;
	G_e.y_index = Je;

	W_e.x_index = dim[0];
	W_e.y_index = dim[1];

	// Reduce 2-D data to master processor
	Communication_reduce_2D_arrays(G_data2d, W_data2d, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);

	// Write 2-D data from master processor

	if (nme=='u'){
	sprintf(filename, "./trn/u_%d.h5", counter);
	}
	if (nme=='v'){
	sprintf(filename, "./trn/v_%d.h5", counter);
	}
	if (nme=='w'){
	sprintf(filename, "./trn/w_%d.h5", counter);
	}
	if (nme=='a'){
	sprintf(filename, "./trn/ubc_%d.h5", counter);
	}
	if (nme=='b'){
	sprintf(filename, "./trn/vbc_%d.h5", counter);
	}
	if (nme=='c'){
	sprintf(filename, "./trn/wbc_%d.h5", counter);
	}
	if (nme=='p'){
	sprintf(filename, "./trn/p_%d.h5", counter);
	}
	Output_2d_data(W_data2d, dim, grid->xu, grid->yc, &(grid->zc[k]), filename,
		data_bag, DTRACE("Output_2d_data"));

	// Free storage (in production code, only free at end of simulation)
	Memory_free_2D_double_array(dim[1], G_data2d);
	Memory_free_2D_double_array(dim[1], W_data2d);
}
