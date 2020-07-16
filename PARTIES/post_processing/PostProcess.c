
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"

#include "dsolver.h"
#include "Communication.h"
#include "Conc.h"
#include "Display.h"
#include "Dtime.h"
#include "Grid.h"
#include "Input.h"
#include "Interpolate.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "Particle.h"
#include "ParticleInput.h"
#include "ParticleOutput.h"
#include "PostProcess.h"
#include "PpEulerian.h"
#include "Pressure.h"
#include "Resume.h"
#include "Temporal.h"
#include "Timer.h"
#include "Velocity.h"
#include "Viscosity.h"
#include "Vorticity.h"

void Cart3d_setup_lsys_accounting_geometry(Cart3d_bag *data_bag);

int main(int argc, char **args) {

	int iter = 0;
	int i, j, k;
	char message[100];
	Velocity *u, *v, *w;
	Pressure *p;
	double T1, T2;
	double Tstart, Tend;

	int numprocs, rank;
	int ierr;

	MPI_Init(&argc, &args);

	Tstart = MPI_Wtime();

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

	// Read input parameters from PARTIES_INPUT_FILE (macro definition in
	// 'definitions.h') and set the other parameters
	for (j = 0; j < numprocs; j++) {
		if (j == rank) {
			Input_set_parameters(params, PARTIES_INPUT_FILE);
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

	/*------------------------------------------------------------------------*/
	/*
	 Create primitive variables based on the input parameters
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// Create velocity fields
	//--------------------------------------------------------------------------
	u = Velocity_create(grid, params, 'u', DTRACE("Velocity_create"));
	data_bag -> u = u;
	Display_progress(params,"u-Velocity has been created successfully...\n");

	v = Velocity_create(grid, params, 'v', DTRACE("Velocity_create"));
	data_bag -> v = v;
	Display_progress(params,"v-Velocity has been created successfully...\n");

	w = Velocity_create(grid, params, 'w', DTRACE("Velocity_create"));
	data_bag -> w = w;
	Display_progress(params,"w-Velocity has been created successfully...\n");

	Vorticity_create(data_bag);
	Display_progress(params,"Vorticity has been created successfully...\n");

	Velocity_update_boundaries(u->data, 'u', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(v->data, 'v', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(w->data, 'w', VEL_TYPE_NORMAL, data_bag);

	//--------------------------------------------------------------------------
	// Create pressure field
	//--------------------------------------------------------------------------
	data_bag -> p = Pressure_create(grid, params);
	Display_progress(params,"Pressure has been created successfully...\n");

	//--------------------------------------------------------------------------
	// Create viscosity fields
	//--------------------------------------------------------------------------
	data_bag -> viscosity = Viscosity_create(grid, params);
	Display_progress(params,"Viscosity has been created successfully...\n");

#ifdef CONC
	//--------------------------------------------------------------------------
	// Create concentration field
	//--------------------------------------------------------------------------
	int iconc;
	int NConc = params->NConc;
	Concentration **c = (Concentration **)malloc(NConc * sizeof(Concentration *));
	Memory_check_allocation(c);

	for (iconc=0; iconc<NConc; iconc++) {

		c[iconc] = Conc_create(iconc, grid, params);
		sprintf(message,"Concentration field[%d] has been created successfully...\n",iconc);
		Display_progress(params, message);
	}
	data_bag -> c = c;
#endif

#ifdef LAG_PARTICLE_RESOLVED
	//--------------------------------------------------------------------------
	// Create Lagrangian structure
	//--------------------------------------------------------------------------
	data_bag -> lag = Lagrangian_create(grid, params);
#endif

	// Setup the linear system for primitive variables based on the fluid nodes
	Cart3d_setup_lsys_accounting_geometry(data_bag);

	// Setup linear solver
	lsolver_transpose_setup(grid, params);


	/*------------------------------------------------------------------------*/
	/*
	 Perform post-processing operations
	 */
	/*------------------------------------------------------------------------*/
	PostProcess(argc, args, data_bag, DTRACE("PostProcess"));

	/*------------------------------------------------------------------------*/
	/*
	 Now, release the allocated memory for the variables
	 */
	/*------------------------------------------------------------------------*/
	Velocity_destroy(u, grid, params);
	Display_progress(params,"u-velocity freed\n");
	Velocity_destroy(v, grid, params);
	Display_progress(params,"v-velocity freed\n");
	Velocity_destroy(w, grid, params);
	Display_progress(params,"w-velocity freed\n");
	Pressure_destroy(data_bag->p, grid, params);
	Display_progress(params,"pressure freed\n");
	Viscosity_destroy(data_bag->viscosity, grid, params);
	Display_progress(params,"viscosity freed\n");
	Vorticity_destroy(data_bag);
	Display_progress(params,"vorticity freed\n");

#ifdef LAG_PARTICLE_RESOLVED
	Lagrangian_destroy(data_bag->lag, grid, params);
#endif

	// NOTE: something is going wrong in this loop, PETSc throws errors
	// (segmentation faults)
#ifdef CONC
	for(iconc=0; iconc <NConc; iconc++) {
		Conc_destroy(c[iconc], iconc, grid, params);
	} // for iconc
	free(c);
	Display_progress(params,"concentration freed\n");
#endif

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
	Communication_finalize();
	printf("Exit here.\n");
	return 0;
}




/******************************************************************************/
/*
 * Perform post-processing operations
 */
/******************************************************************************/
void PostProcess(int argc, char **args, Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i, i_start, i_step, i_end;
	char h5_resume_linker[] = "Resume.h5";
	char h5_resume_temp[50], h5_resume_original[50], temp[50];
	ssize_t len;

	int j, rk;
	const double BET2[] = {BETA2};

	Parameters *params = data_bag->params;
	MAC_grid   *grid   = data_bag -> grid;

	params->resume = 1;
	parse_input(argc, args, &i_start, &i_step, &i_end, params, DTRACE("parse_input"));

	for (i = i_start; i <= i_end; i += i_step) {

#ifdef LAG_PARTICLE_RESOLVED
		// Clean up particle lists from previous loop iterations
		if (i > i_start) {
			Particle_list_destroy(data_bag->lag->p_mobile_list);
			Particle_list_destroy(data_bag->lag->p_fixed_list);
			data_bag->lag->p_mobile_list = (Particle_list *)malloc(sizeof(Particle_list));
			data_bag->lag->p_fixed_list  = (Particle_list *)malloc(sizeof(Particle_list));
		}
#endif

		// Point "Resume.h5" linker to the desired HDF5 file
		sprintf(h5_resume_temp, "Resume_%d.h5", i);
		if (params->rank == 0) {
			if ((len = readlink(h5_resume_linker, h5_resume_original, sizeof(h5_resume_original))) != -1)
				h5_resume_original[len] = '\0';
			unlink(h5_resume_linker);
			symlink(h5_resume_temp, h5_resume_linker);
		}

		// Read in data from HDF5 files
		Resume_h5_resume(data_bag, DTRACE("Resume_h5_resume"));
		Resume_h5_data(data_bag, DTRACE("Resume_h5_data"));
#ifdef LAG_PARTICLE_RESOLVED
		Particle_initialize(data_bag, DTRACE("Particle_initialize"));
#endif

		// Return linker to its original state
		if (params->rank == 0) {
			unlink(h5_resume_linker);
			if (len != -1)
				symlink(h5_resume_original, h5_resume_linker);
		}

		for (j = 0; j < 2; j++) {
			for (rk = 0; rk < 3; rk++) {
				if (rk == 0){
					Velocity_cell_center(data_bag);
					params->dt_old = params->dt;
//					params->dt = Dtime_cfl(data_bag);
				}
				params -> which_stage = rk;
				Temporal_int_all_the_equations(data_bag, DTRACE("Temporal_int_all_the_equations"));
				params -> time += params->dt * BET2[rk];
			}
		}
*/

#ifdef LAG_PARTICLE_RESOLVED
		ParticleOutput_h5(data_bag, params->noutput, DTRACE("ParticleOutput_h5"));
#endif

#ifdef LAG_PARTICLE_RESOLVED
		// Calculate volume fractions
		Particle_MPI_update(data_bag->lag->p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
		Particle_MPI_update(data_bag->lag->p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfu);
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfc);
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfz);
		Interpolate_add_to_volume_fraction('u', data_bag->lag->p_mobile_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('u', data_bag->lag->p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_mobile_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('z', data_bag->lag->p_mobile_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('z', data_bag->lag->p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Particle_list_remove(data_bag->lag->p_mobile_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));
		Particle_list_remove(data_bag->lag->p_fixed_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));
#endif

		PpEulerian(data_bag, DTRACE("PpEulerian"));
	}
}




/******************************************************************************/
/*
 * Parse input arguments to determing 'start', 'step', and 'end' integers for
 * iterating through the HDF5 files.  The inputs can be one of two options:
 *
 * ./postProcess start:end
 *     - specify 'start' and 'end', with 'step' = 1
 * ./postProcess start:step:end
 *     - specify 'start', 'step', and 'end'
 *
 * Any other inputs should produce an error message
 */
/******************************************************************************/
void parse_input(int argc, char **args, int *i_start, int *i_step, int *i_end,
		Parameters *params, Debug_trace *dtrace) {

	int i, i_last, count;
	char *input;
	int input_ints[3];

	if (argc != 2) {
		Display_throw_error("Invalid number of input arguments:\nTHERE CAN BE ONLY ONE\n"
			"Examples: './postProcess 0:2:12' or './postProcess 1:5'", params, DTRACE("Display_throw_error"));
	}

	input = (char *)malloc((strlen(args[1]) + 1) * sizeof(char));
	strcpy(input, args[1]);

	i = 0;
	i_last = 0;
	count = 0;
	while (input[i] != '\0') {

		// Check that we don't have too many delimited numbers
		if (count > 2) {
			Display_throw_error("Invalid input argument: should have 2 to 3 delimited integers.\n"
				"Examples: './postProcess 0:2:12' or './postProcess 1:5'", params, DTRACE("Display_throw_error"));
		}

		// Colon delimits numbers
		if (input[i] == ':') {
			input[i] = '\0';
			input_ints[count] = atoi(&input[i_last]);
			count++;
			i_last = i + 1;
		}
		// Check that we are reading numbers
		else if (input[i] < '0' || input[i] > '9'){
			Display_throw_error("Invalid input argument: should be colon-delimited integers.\n"
				"Examples: './postProcess 0:2:12' or './postProcess 1:5'", params, DTRACE("Display_throw_error"));
		}

		i++;
	}

	// Grab last number
	input_ints[count] = atoi(&input[i_last]);
	count++;

	// Check that we don't have too few delimited numbers
	if (count < 2) {
		Display_throw_error("Invalid input argument: should have 2 to 3 delimited integers.\n"
			"Examples: './postProcess 0:2:12' or './postProcess 1:5'", params, DTRACE("Display_throw_error"));
	}

	// Output parsed numbers
	if (count == 2) {
		*i_start = input_ints[0];
		*i_step  = 1;
		*i_end   = input_ints[1];
	}
	else if (count == 3) {
		*i_start = input_ints[0];
		*i_step  = input_ints[1];
		*i_end   = input_ints[2];
	}

	free(input);
}




/******************************************************************************/
/*
 This function sets up the linear system LHS matrices for all primitive
 variables using the defined grid (geometry)
 */
/******************************************************************************/
void Cart3d_setup_lsys_accounting_geometry(Cart3d_bag *data_bag)
{
	char message[100];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Set up linear system solvers, based on which cells contain fluid, and
	// which are solid
	Velocity_setup_lsys_accounting_geometry(data_bag->u, grid, params);
	Display_progress(params, "Velocity u linear system has been set-up successfully...\n");

	Velocity_setup_lsys_accounting_geometry(data_bag->v, grid, params);
	Display_progress(params, "Velocity v linear system has been set-up successfully...\n");

	Velocity_setup_lsys_accounting_geometry(data_bag->w, grid, params);
	Display_progress(params, "Velocity w linear system has been set-up successfully...\n");

	Pressure_setup_lsys_accounting_geometry(data_bag->p, grid, params);
	Display_progress(params, "Pressure linear system has been set-up successfully...\n");

#ifdef CONC
	int iconc;
	int NConc = params->NConc;
	for (iconc=0; iconc<NConc; iconc++) {
		Conc_setup_lsys_accounting_geometry(data_bag->c[iconc], grid, params);
		sprintf(message, "Concentration field %d's linear system has been set-up successfully...\n", iconc);
		Display_progress(params, message);
	}
#endif
}
