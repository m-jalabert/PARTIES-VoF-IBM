
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"

#include "Communication.h"
#include "Display.h"
#include "Dtime.h"
#include "Grid.h"
#include "Input.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "Output.h"
#include "PartDry.h"
#include "Particle.h"
#include "ParticleInput.h"
#include "ParticleOutput.h"
#include "Resume.h"
#include "Timer.h"

#define TIMEFILE "timesteps.dat"

#ifndef LAG_PARTICLE_RESOLVED
	#error 'LAG_PARTICLE_RESOLVED' not defined for dry particle simulation
#endif

int main(int argc, char **args) {

	int j;
	int iter;
	char message[100];

	int numprocs, rank;
	int ierr;

	MPI_Init(&argc, &args);

	// Get the number of processors and rank of that
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Initialize top-level debug trace structure
	Debug_trace dtrace_parent = {"", "", 0, NULL, NULL};
	Debug_trace *dtrace = &dtrace_parent;

	//--------------------------------------------------------------------------
	// Cart3d_bag: holds all data
	//--------------------------------------------------------------------------
	Cart3d_bag *data_bag = (Cart3d_bag *)malloc(sizeof(Cart3d_bag));
	Memory_check_allocation(data_bag);

	//--------------------------------------------------------------------------
	// Parameters: contains important scalar data
	//--------------------------------------------------------------------------
	Parameters *params = (Parameters *)malloc(sizeof(Parameters));
	Memory_check_allocation(params);

	// Number of processors
	params -> size = numprocs;
	// Current Processor rank
	params -> rank = rank;

	sprintf(message, "parties is running on %d processor(s)", numprocs);
	Display_throw_warning(message, params);

	// Read input parameters from "input.inp" and set the other parameters
	for (j = 0; j < numprocs; j++) {
		if (j == rank) {
			Input_set_parameters(params, "parties.inp");
			data_bag -> params = params;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (params->rank == 0) {
		ierr = Display_parameters(params);
	}
	MPI_Bcast(&ierr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (ierr < 0) {
		Display_throw_error("Error in parameter input", params, DTRACE("Display_throw_error"));
	}

	//--------------------------------------------------------------------------
	// Grid: contains information about the mesh
	//--------------------------------------------------------------------------
	MAC_grid *grid = Grid_create(params, DTRACE("Grid_create"));
	data_bag -> grid = grid;
	Display_progress(params,"Grid has been created successfully...\n");

	Communication_new_xyz_communicator(grid, params);

	//--------------------------------------------------------------------------
	// Timer: contains timing information for the simulation
	//--------------------------------------------------------------------------
	data_bag -> timer = Timer_create(grid, params);

	//--------------------------------------------------------------------------
	// Create Lagrangian structure
	//--------------------------------------------------------------------------
	data_bag -> lag = Lagrangian_create(grid, params);


	/*------------------------------------------------------------------------*/
	/*
	 If resuming simulation
	 */
	/*------------------------------------------------------------------------*/
	if (params->resume) {
		Resume_h5_resume(data_bag, DTRACE("Resume_h5_resume"));
	} // if resuming


	/*------------------------------------------------------------------------*/
	/*
	 Particle Initialization
	 */
	/*------------------------------------------------------------------------*/
	Particle_initialize(data_bag, DTRACE("Particle_initialize"));


	/*------------------------------------------------------------------------*/
	/*
	 Perform 3rd order  Runge-Kutta timestepping to advance to t_final
	 */
	/*------------------------------------------------------------------------*/
	iter = PartDry_int_rk3(data_bag, DTRACE("PartDry_int_rk3"));

	sprintf(message, "Simulation complete, took %d iterations.", iter);
	Display_throw_warning(message, params);

	Display_time_results(data_bag);


	/*------------------------------------------------------------------------*/
	/*
	 Now, release the allocated memory for the variables
	 */
	/*------------------------------------------------------------------------*/
	Lagrangian_destroy(data_bag->lag, grid, params);
	Timer_destroy(data_bag->timer, grid, params);
	Grid_destroy(grid, params);
	Display_progress(params,"grid freed\n");
//	Input_destroy_parameters(params);
//	Display_progress(params,"params freed\n");
	free(data_bag);

	// Free debug trace
	Debug_trace *dtrace_child;
	dtrace = dtrace->child;
	while (dtrace != NULL) {
		dtrace_child = dtrace->child;
		free(dtrace);
		dtrace = dtrace_child;
	}

	Display_progress(params,"data_bag freed\n");
	free(params);
	Communication_finalize();
	printf("Exit here.\n");
	return 0;

}




/******************************************************************************/
/*
 This function integrates the u,v and conc transport equation in time up to
 time_max
 */
/******************************************************************************/
int PartDry_int_rk3(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double time, max_time, dt;
	double output_time_interval;
	int ntime;

	double T1, T2;
	double Tstart, Tend;

	int rk;
	const double BET2[] = {BETA2};
	const int tsubsteps = 3;

	int status, int_stop;
	FILE *fid;
	FILE *timeFile;
	char message[500];

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
	int pnodes         = params -> ghost_nodes;
	Timer      *timer  = data_bag -> timer;


	/*------------------------------------------------------------------------*/
	/*
	 Get the time properties
	 */
	/*------------------------------------------------------------------------*/
	max_time         = params -> time_max;
	output_time_interval = params -> output_time_interval;


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
			params->output_time += output_time_interval;
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
		dt = params -> default_dt;

		ntime = 1;
		time  = 0.0;
		params -> output_time = 0.0;

		params -> time = time;
		params -> noutput = 0;

		if (params -> rank == 0) {
			timeFile = fopen(TIMEFILE, "w");
			fclose(timeFile);
		}
	} // else


	//--------------------------------------------------------------------------
	// Saving the initial conditions
	//--------------------------------------------------------------------------
	if (!params->resume) {

		ParticleOutput_h5(data_bag, params->noutput, DTRACE("ParticleOutput_h5"));

		if (params -> rank == 0) {
			timeFile = fopen(TIMEFILE, "a");
			fprintf(timeFile, "%f\n", time);
			fclose(timeFile);
		}

		params->output_time += output_time_interval;
		params->noutput++;
		Display_progress(params,"Initial flow properties have been saved successfully...\n");
	} // if not resuming


	/*------------------------------------------------------------------------*/
	/*
	 Time-stepping
	 */
	/*------------------------------------------------------------------------*/
	Display_progress(params,"\n*** Beginning main simulation loop ***\n");
	int_stop = 0;
	while ( int_stop == 0 && time - max_time < dt ) {

		Tstart = MPI_Wtime();

		params -> dt    = dt;
		params -> ntime = ntime;
		params -> time  = time;

		//----------------------------------------------------------------------
		// Save data to file during output time steps
		//----------------------------------------------------------------------
		ParticleOutput_dat(data_bag->lag->p_mobile_list, grid, params, DTRACE("ParticleOutput_dat"));
#if defined PARTICLE_TRN
			ParticleOutput_h5(data_bag, ntime, DTRACE("ParticleOutput_h5"));
#endif

		if ( fabs(time - params->output_time) < dt &&
			fabs(time - params->output_time) < fabs(time + dt - params->output_time) ||
			time > params->output_time ) {

			ParticleOutput_h5(data_bag, params->noutput, DTRACE("ParticleOutput_h5"));

			if (params -> rank == 0) {
				timeFile = fopen(TIMEFILE, "a");
				fprintf(timeFile, "%f\n", time);
				fclose(timeFile);
			}
			Output_h5_resume(data_bag, DTRACE("Output_h5_resume"));
			params->output_time += output_time_interval;
			params->noutput++;
			Display_throw_warning("Runtime data has been saved successfully", params);

		} // Output writing

//#ifdef TEST
		if (time != 0) {
			//------------------------------------------------------------------
			// Find dt according to cfl condition and next output write time
			// step
			//------------------------------------------------------------------
			dt = PartDry_get_dt(data_bag);
			params->dt_old = params->dt;
			params->dt = dt;
		} // if
//#endif

		/*--------------------------------------------------------------------*/
		/*
		 Runge Kutta sub-steps
		 */
		/*--------------------------------------------------------------------*/

		//----------------------------------------------------------------------
		// Navier-Stokes and other equations integration
		//----------------------------------------------------------------------
		for (rk = 0; rk < tsubsteps; rk++){

			params -> which_stage = rk;

			PartDry_int_all_the_equations(data_bag, DTRACE("PartDry_int_all_the_equations"));

			time += dt * BET2[rk];
			params -> time = time;

		} // End RK substepping

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
			params->noutput = -params->noutput;
			ParticleOutput_h5(data_bag, params->noutput, DTRACE("ParticleOutput_h5"));
			Output_h5_resume(data_bag, DTRACE("Output_h5_resume"));
#ifdef PARTICLE_TRN
			ParticleOutput_h5(data_bag, ntime, DTRACE("ParticleOutput_h5"));
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
void PartDry_int_all_the_equations(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2, Tstart, Tend;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Timer      *timer  = data_bag -> timer;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;
	Particle *p;

	double dt = params -> dt;
	int which_stage = params -> which_stage;

	Tstart = MPI_Wtime();

	// Reset gravitational force
	p = p_mobile_list -> start;
	while (p != NULL) {
		DSET_ZERO(p->F, 3);
		DSET_ZERO(p->T, 3);
		p = p -> next;
	}

	//--------------------------------------------------------------------------
	// Advect particles
	//--------------------------------------------------------------------------
	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));
	Lagrangian_advect_particles(data_bag, DTRACE("Lagrangian_advect_particles"));

	Tend = MPI_Wtime();
	timer->Wtime_particle_total += Tend - Tstart;

}




/******************************************************************************/
/*
 This function computes dt from the CFL condition and the grid size
 */
/******************************************************************************/
double PartDry_get_dt(Cart3d_bag *data_bag) {

	char message[500];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Grid spacing (assuming uniform grid)
	double h = grid->dx_u[1];

	double dt, dt_old;
	dt_old = params->dt;

	// Timestep required to resolve particle motion
	double dt_lag, dt_lag_min;
	dt_lag_min = params->max_dt;

	// Timestep required to resolve collision stiffness
	double dt_kn_min;

	Particle *p = data_bag->lag->p_mobile_list->start;
	while (p != NULL) {
		// Timestep required for particle to move a distance of cfl*h
		dt_lag = params->cfl * h / sqrt(DOT(p->U, p->U));
		dt_lag_min = min(dt_lag, dt_lag_min);
		p = p -> next;
	}

	MPI_Allreduce(MPI_IN_PLACE, &dt_lag_min, 1, MPI_DOUBLE, MPI_MIN, PCW);

	// Timestep required to resolve particle collision stiffness
	dt_kn_min = Dtime_lag_particle(data_bag);

	// Constant time step
	if (params->constant_dt == 1) {
		if (params->default_dt > dt_lag_min) {
			sprintf(message, "WARNING: Unstable timestep dt = %g\n"
					"CFL stability condition requires dt = %g\n"
					"Continuing anyway...", params->default_dt, dt_lag_min);
			Display_throw_warning(message, params);
		}
		else if (params->default_dt > dt_kn_min) {
			sprintf(message, "WARNING: Unstable timestep dt = %g\n"
					"Stiffness stability condition requires dt = %g\n"
					"Continuing anyway...", params->default_dt, dt_kn_min);
			Display_throw_warning(message, params);
		}
		dt = params -> default_dt;
	}
	else {
		// Set dt to staisfy CFL condition and collision stability
		dt = min(dt_lag_min, dt_kn_min);

		// This is done to avoid sharp jump in time steps.
		if ( dt > 1.1 * dt_old )
			dt = 1.1 * dt_old;
		if ( dt < 0.9 * dt_old )
			dt = 0.9 * dt_old;

		// Compare to maximum allowable dt
		dt = min(dt, params->max_dt);
	}

	return dt;
}
