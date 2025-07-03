//static char help[] = "3-Dimensional Gravity Current Simulation, with Variable Geometry.\n\n";
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Array.h"
#include "Cart3d.h"
#include "Communication.h"
#include "Conc.h"
#include "Display.h"
#include "dsolver.h"
#include "Dtime.h"
#include "Grid.h"
#include "Immersed.h"
#include "Inflow.h"
#include "Input.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "MyMath.h"
#include "Outflow.h"
#include "Output.h"
#include "Particle.h"
#include "ParticleInput.h"
#include "Pressure.h"
#include "Rans.h"
#include "Resume.h"
#include "Schumann.h"
#include "Statistics2d.h"
#include "Strain.h"
#include "Subgrid.h"
#include "Temporal.h"
#include "Timer.h"
#include "Velocity.h"
#include "Viscosity.h"
#include "Initial_Conditions.h"
#include <omp.h>
#include "Interpolate.h"
#include "EPforcing.h"
#include "VolumeFraction.h"




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


#ifdef VOF_PLIC
    //--------------------------------------------------------------------------
	// Create VOF-PLIC field
	//--------------------------------------------------------------------------
    data_bag->vof = VoF_create(grid, params);
    Display_progress(params, "Volume Fraction initialized successfully...\n");
#endif

#ifdef LAG_PARTICLE_RESOLVED
	//--------------------------------------------------------------------------
	// Create Lagrangian structure
	//--------------------------------------------------------------------------
	data_bag -> lag = Lagrangian_create(grid, params);
#endif

#ifdef OUTPUT2D
	//--------------------------------------------------------------------------
	// Create 2-D statistics structure
	//--------------------------------------------------------------------------
	Statistics2d *st2d = Statistics2d_create(grid, params);
	data_bag -> st2d = st2d;
	Display_progress(params,"data_bag has been created successfully...\n");
#endif

#ifdef LES
	//--------------------------------------------------------------------------
	// Create LES structure
	//--------------------------------------------------------------------------
	Subgrid *smag = Subgrid_create(grid, params);
	data_bag -> smag = smag;
	Display_progress(params,"Subgrid has been created successfully...\n");

	#ifdef SMAG_DYNAMIC
	Lesfilter_init(smag, grid, params);
	Dynamic_smag_init(smag, grid, params);
	Display_progress(params,"Filter has been initialized successfully...\n");
	#endif

	#ifdef SCHUMANN
	smag -> log_law = Schumann_create(grid, params);
	Display_progress(params,"Schumann b/c related variables have been initialized successfully...\n");
	#endif
#endif // LES

#ifdef RANS
	//--------------------------------------------------------------------------
	// Create RANS structure
	//--------------------------------------------------------------------------
	data_bag -> rans = Rans_create(grid, params);

	#ifdef TWO_EQUATION_MODEL
	Rans_initialize(data_bag->rans, u, data_bag);
	#endif

	#ifdef SCHUMANN
	data_bag->rans->log_law = Schumann_create(grid, params);
	Display_progress(params,"Schumann b/c related variables have been initialized successfully...\n");
	#endif

	Rans_eddy_viscosity(data_bag->rans, u, params, grid);
#endif

//--------------------------------------------------------------------------
// Random seed
//--------------------------------------------------------------------------
time_t t;
#ifdef TURB_FORCING
if (params->rank == 0) {
	time_t t;
}
MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
srand((unsigned) time(&t));
#else
	srand(time(NULL) + params->rank*10);
#endif


#ifdef TURB_FORCING
//--------------------------------------------------------------------------
// Create Fourier space variable
//--------------------------------------------------------------------------
double L1 = params->xmax -params->xmin;
double L2 = (params->ymax -params->ymin)/L1;
L1 = (params->zmax -params->zmin)/L1;
if (L1 != floor(L1) || L2 != floor(L2)) {
	Display_throw_error("Ly and Lz must be multiples of Lx for the turbulent isotropic forcing", params, DTRACE("Display_throw_error"));
}
Fourier *fourier = (Fourier *)malloc(sizeof(Fourier));
Memory_check_allocation(fourier);
init_turb_forcing(fourier, params);

data_bag -> fourier = fourier;
Display_progress(params,"Turbulent model has been created successfully...\n");
#endif

	/*------------------------------------------------------------------------*/
	/*
	 If the simulation is starting from t=0
	 */
	/*------------------------------------------------------------------------*/
	if (!params->resume) {

		// Use MAC_grid to identify which cells are fluid, and which are solid
		Cart3d_identify_geometry(data_bag);

		// This function initializes the inflow and initial configuration of the
		// concentration fields

		params -> which_stage =0; // used in Outflow_vel_impose_convective_boundary fix that TODO
		Cart3d_initialize_primitive_data(data_bag, DTRACE("Cart3d_initialize_primitive_data"));
		Display_progress(params,"Primitive data initialized\n");

#ifdef IMMERSED_BOUNDARY
//		int pnodes = params -> ghost_nodes;
//		Communication_update_ghost_nodes_flow_variable(u->data, 'u', pnodes, data_bag);
//		Communication_update_ghost_nodes_flow_variable(v->data, 'v', pnodes, data_bag);
//		Communication_update_ghost_nodes_flow_variable(w->data, 'w', pnodes, data_bag);
		Array_copy_withghost(u->data, u->data_old, grid, params);
		Array_copy_withghost(v->data, v->data_old, grid, params);
		Array_copy_withghost(w->data, w->data_old, grid, params);
#endif
	} // if new simulation
	/*------------------------------------------------------------------------*/
	/*
	 If resuming simulation
	 */
	/*------------------------------------------------------------------------*/
	else {
		Resume_h5_resume(data_bag, DTRACE("Resume_h5_resume"));
		Resume_h5_data(data_bag, DTRACE("Resume_h5_data"));

		// Use MAC_grid to identify which cells are fluid, and which are solid
		Cart3d_identify_geometry(data_bag);
	} // if resuming


	/*------------------------------------------------------------------------*/
	/*
	 Particle Initialization
	 */
	/*------------------------------------------------------------------------*/
#ifdef LAG_PARTICLE_RESOLVED
	Particle_initialize(data_bag, DTRACE("Particle_initialize"));
#endif


	/*------------------------------------------------------------------------*/
	/*
	 Setup immersed boundary nodes and coefficients to impose no-slip B.C.
	 for the velocity field on the solid boundary
	 */
	/*------------------------------------------------------------------------*/
#ifdef CONSTANT_MASSFLUX
	params->ubulk = Velocity_ubulk(u, data_bag);
#endif

#ifdef IMMERSED_BOUNDARY
	Immersed_setup_q_immersed_nodes(grid, params, 'u');
	Display_progress(params, "Immersed boundary has been setup successfully for u\n");
	Immersed_setup_q_immersed_nodes(grid, params, 'v');
	Display_progress(params, "Immersed boundary has been setup successfully for v\n");
	Immersed_setup_q_immersed_nodes(grid, params, 'w');
	Display_progress(params, "Immersed boundary has been setup successfully for w\n");
	Immersed_setup_q_immersed_nodes(grid, params, 'c');
	Display_progress(params, "Immersed boundary has been setup successfully for c\n");

	#if defined LES || defined RANS
	Immersed_copy_to_nut(grid, params);
	#endif

	Display_progress(params, "Immersed boundary has been setup successfully\n");
#endif

	// Setup the linear system for primitive variables based on the fluid nodes
	//Cart3d_setup_lsys_accounting_geometry(data_bag);

	// Setup linear solver
	//lsolver_transpose_setup(grid, params);



#if defined CONC && defined LAG_PARTICLE_RESOLVED // Scalar is set zero inside particle

	Particle_MPI_update(data_bag->lag->p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(data_bag->lag->p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));

	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfc);
	Interpolate_add_to_volume_fraction('c', data_bag->lag->p_mobile_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('c', data_bag->lag->p_fixed_list,data_bag, DTRACE("Interpolate_add_to_volume_fraction"));


#ifdef VOF_SCALAR

		Memory_reset_flow_variable(grid, params, data_bag->lag->vfu);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfv);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfw);

		Memory_reset_flow_variable(grid, params, data_bag->u->data_vof);
		Memory_reset_flow_variable(grid, params, data_bag->v->data_vof);
		Memory_reset_flow_variable(grid, params, data_bag->w->data_vof);

#ifdef VOF_SMOOTH_VELO
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfu_prime);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfv_prime);
		Memory_reset_flow_variable(grid, params, data_bag->lag->vfw_prime);
#endif


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

		Lagrangian_calc_compound_vof_fluid_velo(data_bag);

		Communication_update_ghost_nodes_flow_variable(u->data_vof, 'u', params -> ghost_nodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(v->data_vof, 'v', params -> ghost_nodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(w->data_vof, 'w', params -> ghost_nodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(data_bag->lag->vfu, 'u', params -> ghost_nodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(data_bag->lag->vfv, 'v', params -> ghost_nodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(data_bag->lag->vfw, 'w', params -> ghost_nodes, data_bag);
#endif





	Particle_list_remove(data_bag->lag->p_mobile_list, FOREIGN, grid, params,  dtrace);
	Particle_list_remove(data_bag->lag->p_fixed_list, FOREIGN, grid, params,  dtrace);
#endif


#ifdef LAG_PARTICLE_RESOLVED
	#ifdef SCALAR_INIT_MASK
	if (!params->resume){

		for (iconc=0; iconc<NConc; iconc++) {
		Lagrangian_mask_scalar( iconc, data_bag);
		}
	}
#endif
#endif

	/*------------------------------------------------------------------------*/
	/*
	 Display grid information on each processor
	 */
	/*------------------------------------------------------------------------*/
//	Display_DA_3D_info(grid, params);

	Tend = MPI_Wtime();
	data_bag->timer->Wtime_init += Tend - Tstart;





	/*------------------------------------------------------------------------*/
	/*
	 Perform 3rd order  Runge-Kutta timestepping to advance to t_final
	 */
	/*------------------------------------------------------------------------*/
	iter = Temporal_int_rk3(data_bag, DTRACE("Temporal_int_rk3"));

	sprintf(message, "Simulation complete, took %d iterations.", iter);
	Display_throw_warning(message, params);

	Display_time_results(data_bag);


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

#ifdef OUTPUT2D
	Statistics2d_destroy(st2d, grid, params);
	Display_progress(params,"2d statistics freed\n");
#endif
/*
	Communication_finalize();
	printf("Exit here.\n");
	exit(0);
	return;
*/

#ifdef LAG_PARTICLE_RESOLVED
	Lagrangian_destroy(data_bag->lag, grid, params);
#endif

#ifdef VOF_PLIC
    // Destroy Volume Fraction structure
    VOF_destroy(data_bag->vof, grid, params);
    Display_progress(params, "Volume Fraction freed\n");
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




/******************************************************************************/
/*
 This function initializes the inflow and initial configuration of the
 concentration fields
 */
/******************************************************************************/
void Cart3d_initialize_primitive_data(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	char message[100];
	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;


/// I think the nex lines are unncessary
#if defined LEFT_INFLOW || defined RIGHT_INFLOW
	Inflow_velocity_profile(data_bag, DTRACE("Inflow_velocity_profile"));
	Display_progress(params,"u-Velocity inflow profile has been set successfully...\n");
#endif
#if defined LEFT_OUTFLOW || defined RIGHT_OUTFLOW
	Outflow_impose_convective_boundary(data_bag);
	Display_progress(params,"Outflow profile has been set successfully...\n");
#endif


	/*------------------------------------------------------------------------*/
	/*
	 Impose slip-wall conditions
	 */
	/*------------------------------------------------------------------------*/
 	Velocity_update_boundaries(data_bag->u->data, 'u', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(data_bag->v->data, 'v', VEL_TYPE_NORMAL, data_bag);
 	Velocity_update_boundaries(data_bag->w->data, 'w', VEL_TYPE_NORMAL, data_bag);

#ifdef CONC
	int iconc;
	int NConc;
	NConc = params->NConc;
	Concentration **c = data_bag -> c;



	sprintf(message,"Entering the new init routine");
			Display_progress(params, message);
	for (iconc=0; iconc<NConc; iconc++) {

		switch ( data_bag->params->conc_init_type[iconc])
		{
		case(1):
				Conc_init_lock(c[iconc], data_bag, DTRACE("Conc_initialize"));
		break;


		case(2):
				Conc_init_cos(c[iconc], data_bag);
		break;

		case(3):
				Conc_init_sin(c[iconc], data_bag);
		break;

		case(4):
				Conc_init_fullsin(c[iconc], data_bag);

		break;

		case(5):

				Conc_init_xsin(c[iconc], data_bag);

		break;

		case(6):

				Conc_init_RB(c[iconc], data_bag);
		break;

		case(7):

				Conc_init_Ardekani(c[iconc], data_bag);

		break;

		case(8):

				Conc_init_erf_x(c[iconc], data_bag);

		break;


		case(9):

				Conc_init_erf_y(c[iconc], data_bag);

		break;

		case(10):

				Conc_init_unity(c[iconc], data_bag);

		break;

		case(101):

				Conc_init_layer(iconc ,c[iconc], data_bag);

		break;

		case(102):

				Conc_init_intrusion_DD(iconc ,c[iconc], data_bag);

		break;

		case(103):

				Conc_init_intrusion_DD_partlock(iconc ,c[iconc], data_bag);

		break;

		default:
				Conc_init_zero(c[iconc], data_bag); // conc is set to zero

		}


		sprintf(message,"Concentration field %d has been initialized successfully...\n", iconc);
		Display_progress(params, message);
		Conc_set_boundary_values(c[iconc]->data, iconc ,CENTRAL_FULL, grid , params);

	} // for iconc

#endif

#ifdef VOF_PLIC
    // Initialize VOF field based on init_type
    switch(params->init_type) {
        case 1:  // Advection test case
            VoF_init_bubble(data_bag);
            break;
        
        case 2: // stationary droplet test case
			VoF_init_stationary_droplet(data_bag);
            break;
			
		case 3: // ellipsoid droplet test case
			VoF_init_ellipsoid(data_bag);
            break;

		case 4: // rising bubble test case
			VoF_init_rising_bubble(data_bag);
            break;

		case 5: // rising multiple bubbles 
			VoF_init_two_bubbles_coaxial(data_bag);
            break;

		case 6: 
		    VoF_init_vertical_bilayer_Z(data_bag);
			break;
    }
    
        // Set boundary values and update ghost nodes
        VOF_set_boundary_values(data_bag->vof->F, data_bag);		     
        VoF_smoothing(data_bag);
    	curvature_patel(data_bag);
        
        // Update mixture properties (density, viscosity)
        VOF_update_density_viscosity(data_bag);
        Display_progress(params, "VOF initialized successfully...\n");
    
#endif


}




/******************************************************************************/
/*
 This function defines the geometry and then tags the nodes based on the
 location of the node with respect to the solid interface as:
     SOLID, FLUID, IMMERSED, BOUNDARY
 */
/******************************************************************************/
void Cart3d_identify_geometry(Cart3d_bag *data_bag) {

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	int pnodes = params->ghost_nodes;

	// First, generate the location of the interface at the grid centers
	Grid_describe_interface(grid, params);
	if (params->rank==0)
		printf("Grid has been described successfully\n");

	Immersed_find_surface_normal_distance(grid, params, 'u');
	if (params->rank==0) printf("After u\n");
	Immersed_find_surface_normal_distance(grid, params, 'v');
	if (params->rank==0) printf("After v\n");
	Immersed_find_surface_normal_distance(grid, params, 'w');
	if (params->rank==0) printf("After w\n");
	Immersed_find_surface_normal_distance(grid, params, 'c');
	if (params->rank==0) printf("After c\n");
	Communication_update_ghost_nodes_flow_variable(grid->u_sdf, 'u', pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(grid->v_sdf, 'v', pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(grid->w_sdf, 'w', pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(grid->c_sdf, 'c', pnodes, data_bag);

	// This function tags all the nodes (u, v, w, c) based on the surface
	// position.
	// Tagging is done via a levelset represenation of the interface for the
	// solid inteface.
	Grid_identify_geometry(data_bag);

	if (params->rank==0)
		printf("Geometry has been created successfully \n");

#ifdef IMMERSED_BOUNDARY1
	int i, j, k;
	double ***grid_st;
	grid_st =  Memory_allocate_flow_variable(grid, params);
	for (k=grid->G_Ks; k<grid->G_Ke; k++) {
		for (j=grid->G_Js; j<grid->G_Je; j++) {
			for (i=grid->G_Is; i<grid->G_Ie; i++) {
				grid_st[k][j][i] = 0.0;
				if (grid->c_status[k][j][i] == SOLID) grid_st[k][j][i] = 1.0;
			}
		}
	}
	Output_3d_data(grid_st, grid, params, 'c', "Data_imm_c.h5", DTRACE("Output_3d_data"));
	Memory_free_flow_variable(grid, params, grid_st);
#endif

return;
}
