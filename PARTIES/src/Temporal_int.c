#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

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
#include "VolumeFraction.h"
#include "VOF_DIFFUSE.h"
#include <math.h>
#include <float.h>


#define TIMEFILE "timesteps.dat"




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

	//Debugging
	int rank = data_bag->params->rank;

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

	#ifdef SLICE_OUTPUT

		double ***slice_data_1, ***slice_data_2;
		char name_data_1, name_data_2;

		if (params->slice_axis==0){ // slice of X-axis
			slice_data_1 = data_bag->w->data;
			slice_data_2 = data_bag->v->data;
			name_data_1 = 'w';
			name_data_2 = 'v';

		} else if (params->slice_axis==1){ // slice of Y-axis
			slice_data_1 = data_bag->u->data;
			slice_data_2 = data_bag->w->data;
			name_data_1 = 'u';
			name_data_2 = 'w';

		} else if (params->slice_axis==2){ // slice of Z-axis
			slice_data_1 = data_bag->u->data;
			slice_data_2 = data_bag->v->data;
			name_data_1 = 'u';
			name_data_2 = 'v';
		} 
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

		#if defined OUTPUT2D || defined SLICE_OUTPUT || defined PARTICLE_TRN
			if( access( "./trn", F_OK ) == -1) {
				mkdir("./trn", 0777);
			}
		#endif


		#if defined OUTPUT2D || defined SLICE_OUTPUT

			#ifdef OUTPUT2D
				Statistics2d_computeStatistics(data_bag, DTRACE("Statistics2d_computeStatistics"));
			#endif


			#ifdef SLICE_OUTPUT
				slice_2d_output(slice_data_1, name_data_1, data_bag, time, DTRACE("slice_2d_output"));
				slice_2d_output(slice_data_2, name_data_2, data_bag, time, DTRACE("slice_2d_output"));
				if (params->slice_p == 1) {
					slice_2d_output(data_bag->p->p_data_avg, 'p', data_bag, time, DTRACE("slice_2d_output"));
				}
			#endif

			params->output_time_2d += output_time_interval_2d;
			params->noutput_2d++;

		#endif

		params->output_time += output_time_interval;
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
		params->dt = dt;


		//----------------------------------------------------------------------
		// Save data to file during output time steps
		//----------------------------------------------------------------------

		#if defined LAG_PARTICLE_RESOLVED// && !defined SUBSTEP

			#if defined PARTICLE_TRN
				sprintf(message, "Outputting par_trn %d\n", ntime);
				Display_progress(params, message);
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
				//	xzperiodic_saltsediment_ave(data_bag, c, grid, params);
				#ifdef THERMAL_KADER
						xzperiodic_thermal_ave(data_bag, c, grid, params);
				#endif
				#ifdef RANS
						xzperiodic_rans_nut_ave(data_bag, data_bag->rans, grid, params);
				#endif
			#endif
		} // Output writing


		#if defined OUTPUT2D || defined SLICE_OUTPUT
			if ( (fabs(time - params->output_time_2d) < dt &&
		     	fabs(time - params->output_time_2d) < fabs(time + dt - params->output_time_2d)) ||
		     	(time > params->output_time_2d) ) {

				#ifdef OUTPUT2D
					Statistics2d_computeStatistics(data_bag, DTRACE("Statistics2d_computeStatistics"));
				#endif

				#ifdef SLICE_OUTPUT
					slice_2d_output(slice_data_1, name_data_1, data_bag, time, DTRACE("slice_2d_output"));
					slice_2d_output(slice_data_2, name_data_2, data_bag, time, DTRACE("slice_2d_output"));
					if (params->slice_p == 1) {
						slice_2d_output(data_bag->p->p_data_avg, 'p', data_bag, time, DTRACE("slice_2d_output"));
					}
				#endif

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

				// If we're running the advection test, force the analytical velocity field.
				if (params->vel_init_type == VEL_INIT_ADVECTION_TEST) {
					Inflow_velocity_profile(data_bag, DTRACE("Inflow_velocity_profile"));
				}

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
					params->dt = dt;
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

		double dpdx = params -> dp_dx;
		double ubulk_target = params -> ubulk_target;

		Velocity_calculate_dpdx(u, data_bag);

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
			#if defined OUTPUT2D || defined SLICE_OUTPUT
				if ( (fabs(time - params->output_time_2d) < dt &&
					fabs(time - params->output_time_2d) < fabs(time + dt - params->output_time_2d)) ||
					(time > params->output_time_2d) ) {

					#ifdef OUTPUT2D
						Statistics2d_computeStatistics(data_bag, DTRACE("Statistics2d_computeStatistics"));
					#endif

					#ifdef SLICE_OUTPUT
						slice_2d_output(slice_data_1, name_data_1, data_bag, time, DTRACE("slice_2d_output"));
						slice_2d_output(slice_data_2, name_data_2, data_bag, time, DTRACE("slice_2d_output"));
						if (params->slice_p == 1) {
							slice_2d_output(data_bag->p->p_data_avg, 'p', data_bag, time, DTRACE("slice_2d_output"));
						}
					#endif

					params->output_time_2d += output_time_interval_2d;
					params->noutput_2d++;
				}
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

#ifdef VOF
	VolumeFraction *vof = data_bag -> vof;
#endif

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

#ifdef LAG_PARTICLE_RESOLVED
	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));
	Lagrangian_flag_points(data_bag, DTRACE("Lagrangian_flag_points"));  // New Langragian points are calculated

	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfc);
	Interpolate_add_to_volume_fraction('c', data_bag->lag->p_mobile_list,
	data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('c', data_bag->lag->p_fixed_list,
	data_bag, DTRACE("Interpolate_add_to_volume_fraction"));

#endif

#ifdef VOF

    int rank = data_bag->params->rank;

    /* ── VOF advection: method-specific branch ────────────────────────── */
    #ifdef VOF_PLIC

		#ifdef VOF_IBM
            VOF_normals_IBM(data_bag);
            Vfc_smoothing(data_bag);
    	#endif
        
        VOF_reconstruct_interface(data_bag);
        VOF_set_advection(data_bag);

    #elif defined VOF_DIFFUSE

		#ifdef VOF_IBM
		VOF_DIFFUSE_compute_C_S(data_bag);
		#endif
		        
        VOF_DIFFUSE_step(data_bag);

    #endif /* advection method */

    // RK stage update of F
    #ifdef VOF_PLIC
    VOF_update_F(data_bag);
    #endif

    #if defined VOF_IBM && defined SURFACE_TENSION
    #ifdef VOF_PLIC
        VOF_extend(data_bag);
    #endif
    #endif

    #ifdef VOF_PLIC
    // Reconstruct interface after update (ready for surface tension)
    VOF_reconstruct_interface(data_bag);
    #endif

    // Update density/viscosity
    VOF_update_density_viscosity(data_bag);

#endif /* VOF */

	/*------------------------------------------------------------------------*/
	/*
	 For a "pure VoF advection test," we skip the velocity solve & Poisson:
	 */
	/*------------------------------------------------------------------------*/
 		if (params->vel_init_type == VEL_INIT_ADVECTION_TEST) 
  		{
		
			//Prescribed velocity test: skipping momentum solver but advancing time;
			return;
		}



	/*------------------------------------------------------------------------*/
	/*
	 Momentum Solver and Poisson Solver
	 */
	/*------------------------------------------------------------------------*/
		else {			

	/*------------------------------------------------------------------------*/
	/*
	 Calculate the explicitly and semi-implicitly treated part of the RHS for
	 momentum equation
	 */
	/*------------------------------------------------------------------------*/


	#ifdef VOF_PLIC
	/* Compute Bussmann-consistent face momentum-flux density arrays */
	VOF_compute_conservative_momentum_fluxes(data_bag);
	#endif

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
	Velocity_v_set_RHS(data_bag);
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


	#ifdef SURFACE_TENSION

		#ifdef VOF_DIFFUSE
		VOF_DIFFUSE_compute_psi_LG(data_bag);
		#ifdef VOF_IBM
		if (params->which_stage == 2) {
			VOF_DIFFUSE_extend_psi_LG_contact_angle(data_bag);
		}
		#endif
		VOF_DIFFUSE_compute_f_sigma(data_bag);
		VOF_apply_f_sigma_old(data_bag);

		#else

    // /* Smooth VOF*/
     VoF_smoothing(data_bag);

    // /*  Reconstruct interface curvature */
     curvature_patel(data_bag);

    VOF_compute_f_sigma(data_bag);     // compute f_σ^k → f_sigma_new  (from F^k)
    VOF_apply_f_sigma_old(data_bag);   // add   f_σ^{k-1} → RHS        (balances ∇p^{k-1})

	#endif

#endif


#ifdef VOF_GRAVITY
  Velocity_add_gravity_2_RHS(data_bag);
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


	/*------------------------------------------------------------------------*/
	/*
	 Solve for pre-projected velocity
	 */
	/*------------------------------------------------------------------------*/


	T1 = MPI_Wtime();
	Velocity_solve(u, data_bag);
	Velocity_solve(v, data_bag);
	Velocity_solve(w, data_bag);
	T2 = MPI_Wtime();
	timer->Wtime_vel_solve += T2 - T1;




	//--------------------------------------------------------------------------
	// Update boundary conditions
	//--------------------------------------------------------------------------


 	Velocity_update_boundaries(u->data, 'u', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(v->data, 'v', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(w->data, 'w', VEL_TYPE_NORMAL, data_bag);

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
		Memory_reset_flow_variable(grid, params, p->deltap);
		Pressure_set_RHS(data_bag);
		T2 = MPI_Wtime();


		timer->Wtime_p_rhs += T2 - T1;

		T1 = MPI_Wtime();

#ifdef VOF // Pressure solved with CG method for VOF

	#ifdef USE_HYPRE
		Pressure_solve_hypre(data_bag);
	#else
		Pressure_compute_preconditioner(data_bag);
		Pressure_solve_cg(data_bag);
	#endif

#else
		Pressure_solve(p, grid, params);
#endif

		T2 = MPI_Wtime();
		timer->Wtime_p_solve += T2 - T1;


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

	// Advance balanced-force arrays for next substep
    #ifdef SURFACE_TENSION
    VOF_swap_f_sigma(data_bag);
    #endif


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
				#ifdef VOF_IBM
					p->Int_rho_scalar = 0.0;     
					DSET_ZERO(p->F_CCF, 3);
					DSET_ZERO(p->T_CCF, 3);
					DSET_ZERO(p->F_CSF_solid, 3);  
                    DSET_ZERO(p->T_CSF_solid, 3);  
				#endif
				p = p -> next;
			}
			p = p_fixed_list -> start;
			while (p != NULL) {
				DSET_ZERO(p->Int_U, 3);
				DSET_ZERO(p->Int_Omega, 3);
				#ifdef VOF_IBM   
					 p->Int_rho_scalar = 0.0;    
					 DSET_ZERO(p->F_CCF, 3);
					 DSET_ZERO(p->T_CCF, 3);  
					 DSET_ZERO(p->F_CSF_solid, 3);  
                     DSET_ZERO(p->T_CSF_solid, 3); 
				#endif
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

		#ifdef SLICE_OUTPUT
			if ( params->slice_axis == 0 && params->center_two_particles == 1 ) {

				FORI3 params->particle_position[0][i] = 0;
				FORI3 params->particle_position[1][i] = 0;

				Particle *par = p_mobile_list -> start;
				while (par != NULL) {
					if ( par->ID == 0 ) {
						FORI3 params->particle_position[0][i] = par->X[i];
					}
					if ( par->ID == 1 )	{
						FORI3 params->particle_position[1][i] = par->X[i];
					}
					par = par->next;	  
				}

				MPI_Allreduce(MPI_IN_PLACE, &(params->particle_position), 6, MPI_DOUBLE, MPI_SUM, PCW);
			} 
		#endif

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

		} // end of else

}	// end of function
