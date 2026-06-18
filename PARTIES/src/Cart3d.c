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
#include "VOF_DIFFUSE.h"

static void Cart3d_twod_abort(Parameters *params, const char *message,
		const char *function, int line)
{
	Debug_trace trace = {"", "", 0, NULL, NULL};

	snprintf(trace.function, sizeof(trace.function), "%s", function);
	snprintf(trace.file, sizeof(trace.file), "%s", __FILE__);
	trace.line = line;
	Display_throw_error(message, params, &trace);
}

#define CART3D_TWOD_ABORT(params, message) \
	Cart3d_twod_abort((params), (message), __func__, __LINE__)

static void Cart3d_twod_assert_all(int local_ok, Parameters *params,
		const char *message)
{
	int global_ok = 0;

	MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_LAND, PCW);
	if (!global_ok)
		CART3D_TWOD_ABORT(params, message);
}

static double Cart3d_twod_exchange_value(int i, int j, int k)
{
	return 1000000.0 * k + 1000.0 * j + (double)i;
}

static void Cart3d_run_twod_pretests(Cart3d_bag *data_bag)
{
	int i, j, k;
	int local_ok = 1;
	int pnodes = 1;
	MAC_grid *grid = data_bag->grid;
	Parameters *params = data_bag->params;
	double ***probe = NULL;

	Cart3d_twod_assert_all(grid->G_Ks == 0 && grid->G_Ke == grid->NZ &&
	                       grid->NK == 1,
	                       params,
	                       "Phase 0 pre-test failed: the z-direction does not reduce to a single physical slab.");

	Cart3d_twod_assert_all(grid->L_Ks == grid->G_Ks - params->ghost_nodes &&
	                       grid->L_Ke == grid->G_Ke + params->ghost_nodes,
	                       params,
	                       "Phase 0 pre-test failed: z ghost extents are not a pure ghost-layer expansion.");

	Cart3d_twod_assert_all(grid->zw[1] > grid->zw[0] &&
	                       grid->zc[0] > grid->zw[0] &&
	                       grid->dummy_z_slab_thickness > 0.0,
	                       params,
	                       "Phase 0 pre-test failed: the dummy slab coordinates are invalid for NZM = 1.");

	if (params->axisym_rz_enabled) {
		Cart3d_twod_assert_all(grid->r_u != NULL && grid->r_c != NULL &&
		                       grid->inv_r_u != NULL && grid->inv_r_c != NULL &&
		                       grid->ring_wt_u != NULL && grid->ring_wt_c != NULL &&
		                       grid->r_u[0] == 0.0 && grid->r_c[0] >= 0.0 &&
		                       grid->ring_wt_u[0] == 0.0,
		                       params,
		                       "Phase 0 pre-test failed: axisymmetric radial metric arrays were not initialized correctly.");
	}
	else {
		Cart3d_twod_assert_all(grid->r_u == NULL && grid->r_c == NULL &&
		                       grid->inv_r_u == NULL && grid->inv_r_c == NULL &&
		                       grid->ring_wt_u == NULL && grid->ring_wt_c == NULL,
		                       params,
		                       "Phase 0 pre-test failed: Cartesian 2D mode unexpectedly allocated axisymmetric radial metrics.");
	}

	probe = Memory_allocate_flow_variable(grid, params);
	local_ok = (probe != NULL) &&
	           (grid->ng_total_nodes ==
	            (grid->G_Ie - grid->G_Is) *
	            (grid->G_Je - grid->G_Js) *
	            (grid->G_Ke - grid->G_Ks));
	Cart3d_twod_assert_all(local_ok,
	                       params,
	                       "Phase 0 pre-test failed: 3D field allocation or no-ghost sizing is inconsistent for NZM = 1.");
	Memory_free_flow_variable(grid, params, probe);

	probe = Memory_allocate_flow_variable(grid, params);
	for (k = grid->G_Ks; k < grid->G_Ke; k++) {
		for (j = grid->G_Js; j < grid->G_Je; j++) {
			for (i = grid->G_Is; i < grid->G_Ie; i++) {
				probe[k][j][i] = Cart3d_twod_exchange_value(i, j, k);
			}
		}
	}

	Communication_update_ghost_nodes_x(probe, CONCENTRATION_PERTURBATION, pnodes, data_bag);
	local_ok = 1;
	if (params->npxminus != MPI_PROC_NULL) {
		for (j = grid->G_Js; j < grid->G_Je && local_ok; j++) {
			for (k = grid->G_Ks; k < grid->G_Ke && local_ok; k++) {
				local_ok = fabs(probe[k][j][grid->G_Is - 1] -
				                Cart3d_twod_exchange_value(grid->G_Is - 1, j, k)) < 1.0e-12;
			}
		}
	}
	if (params->npxplus != MPI_PROC_NULL) {
		for (j = grid->G_Js; j < grid->G_Je && local_ok; j++) {
			for (k = grid->G_Ks; k < grid->G_Ke && local_ok; k++) {
				local_ok = fabs(probe[k][j][grid->G_Ie] -
				                Cart3d_twod_exchange_value(grid->G_Ie, j, k)) < 1.0e-12;
			}
		}
	}
	Cart3d_twod_assert_all(local_ok,
	                       params,
	                       "Phase 0 pre-test failed: x-direction halo exchange did not preserve the in-plane slab data.");

	Communication_update_ghost_nodes_y(probe, CONCENTRATION_PERTURBATION, pnodes, data_bag);
	local_ok = 1;
	if (params->npyminus != MPI_PROC_NULL) {
		for (i = grid->G_Is; i < grid->G_Ie && local_ok; i++) {
			for (k = grid->G_Ks; k < grid->G_Ke && local_ok; k++) {
				local_ok = fabs(probe[k][grid->G_Js - 1][i] -
				                Cart3d_twod_exchange_value(i, grid->G_Js - 1, k)) < 1.0e-12;
			}
		}
	}
	if (params->npyplus != MPI_PROC_NULL) {
		for (i = grid->G_Is; i < grid->G_Ie && local_ok; i++) {
			for (k = grid->G_Ks; k < grid->G_Ke && local_ok; k++) {
				local_ok = fabs(probe[k][grid->G_Je][i] -
				                Cart3d_twod_exchange_value(i, grid->G_Je, k)) < 1.0e-12;
			}
		}
	}
	Cart3d_twod_assert_all(local_ok,
	                       params,
	                       "Phase 0 pre-test failed: y-direction halo exchange did not preserve the in-plane slab data.");
	Memory_free_flow_variable(grid, params, probe);

	probe = Memory_allocate_flow_variable(grid, params);
	for (j = grid->G_Js; j < grid->G_Je; j++) {
		for (i = grid->G_Is; i < grid->G_Ie; i++) {
			probe[grid->G_Ks][j][i] = Cart3d_twod_exchange_value(i, j, grid->G_Ks);
		}
	}
	Communication_update_ghost_nodes_z(probe, CONCENTRATION_PERTURBATION, pnodes, data_bag);
	local_ok = 1;
	for (j = grid->G_Js; j < grid->G_Je && local_ok; j++) {
		for (i = grid->G_Is; i < grid->G_Ie && local_ok; i++) {
			double expected = Cart3d_twod_exchange_value(i, j, grid->G_Ks);
			local_ok = fabs(probe[grid->G_Ks - 1][j][i] - expected) < 1.0e-12 &&
			           fabs(probe[grid->G_Ke - 1][j][i] - expected) < 1.0e-12;
		}
	}
	Cart3d_twod_assert_all(local_ok,
	                       params,
	                       "Phase 0 pre-test failed: z-direction halo exchange did not collapse to a safe periodic slab fill.");
	Memory_free_flow_variable(grid, params, probe);
}

static void Cart3d_validate_twod_startup(Cart3d_bag *data_bag)
{
#ifdef TWOD_MODE
	Parameters *params = data_bag->params;
	MAC_grid *grid = data_bag->grid;
	char message[256];

	Cart3d_twod_assert_all(params->NZM == 1,
	                       params,
	                       "TWOD_MODE requires NZM == 1.");
	Cart3d_twod_assert_all(params->NPZ == 1,
	                       params,
	                       "TWOD_MODE requires NPZ == 1 after domain decomposition.");

#ifndef ZPERIODIC
	CART3D_TWOD_ABORT(params, "TWOD_MODE requires ZPERIODIC.");
#endif

	snprintf(message, sizeof(message),
	         "TWOD_MODE active: z-direction extent [%.6g, %.6g] is treated as a storage-only slab of thickness %.6g.\n",
	         params->zmin, params->zmax, params->twod_slab_thickness);
	Display_progress(params, message);

#ifdef AXISYM_RZ
	Cart3d_twod_assert_all(fabs(params->xmin) < 1.0e-12,
	                       params,
	                       "AXISYM_RZ requires xmin == 0.0 so the left boundary is the axis.");
	Cart3d_twod_assert_all(params->axisym_no_swirl,
	                       params,
	                       "AXISYM_RZ phase 0 requires AXISYM_NO_SWIRL.");
#if defined(LEFT_WALL_VELOCITY_NOSLIP) || defined(LEFT_WALL_VELOCITY_FREESLIP) || defined(LEFT_INFLOW) || defined(LEFT_OUTFLOW) || defined(XPERIODIC)
	CART3D_TWOD_ABORT(params,
	                  "AXISYM_RZ requires the code x-min boundary to be reserved for the symmetry axis; left wall, inflow/outflow, or x-periodic logic is not allowed.");
#endif
#endif

	Cart3d_twod_assert_all(grid->dummy_z_slab_thickness > 0.0,
	                       params,
	                       "TWOD_MODE requires a positive dummy slab thickness.");
	Cart3d_run_twod_pretests(data_bag);
	Display_progress(params, "Phase 0 2D configuration and pre-tests passed.\n");
#else
	(void)data_bag;
#endif
}




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
	Cart3d_validate_twod_startup(data_bag);

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


#ifdef VOF
    //--------------------------------------------------------------------------
	// Create VOF field
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

		Cart3d_identify_geometry(data_bag);

	#ifdef LAG_PARTICLE_RESOLVED
	#if defined(VOF_DIFFUSE) && defined(VOF_IBM)
		/*
		* Diffuse VOF/IBM needs p_fixed/p_mobile before VOF_DIFFUSE_init(),
		* because C_S is built from the particle lists.
		*/
		Particle_initialize(data_bag, DTRACE("Particle_initialize"));
	#endif
	#endif

		params->which_stage = 0;
		Cart3d_initialize_primitive_data(data_bag,
										DTRACE("Cart3d_initialize_primitive_data"));
		Display_progress(params, "Primitive data initialized\n");

	#ifdef LAG_PARTICLE_RESOLVED
	#if defined(VOF_DIFFUSE) && defined(VOF_IBM)
		/*
		 * Particle lists were read before VOF_DIFFUSE_init() so C_S could be
		 * built from them.  Now C_S and rho are initialized, so the
		 * density-weighted startup velocity paint can run with the correct
		 * solid mask.
		 */
		Particle_initialize_velocities(data_bag,
									DTRACE("Particle_initialize_velocities"));
	#endif
	#endif

	#ifdef IMMERSED_BOUNDARY
		Array_copy_withghost(u->data, u->data_old, grid, params);
		Array_copy_withghost(v->data, v->data_old, grid, params);
		Array_copy_withghost(w->data, w->data_old, grid, params);
	#endif
	}	// if new simulation
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
	#if defined(VOF_DIFFUSE) && defined(VOF_IBM)
		/*
		* For new diffuse VOF/IBM runs, particles were already initialized before
		* VOF_DIFFUSE_init().  For resume runs, still initialize them here.
		*/
		if (params->resume) {
			Particle_initialize(data_bag, DTRACE("Particle_initialize"));
		}
	#else
		Particle_initialize(data_bag, DTRACE("Particle_initialize"));
	#endif
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

#ifndef VOF
	// Setup the linear system for primitive variables based on the fluid nodes
	Cart3d_setup_lsys_accounting_geometry(data_bag);

	// Setup linear solver
	lsolver_transpose_setup(grid, params);
#endif


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

#ifdef USE_HYPRE
	Pressure_hypre_setup(data_bag);
#endif





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

#ifdef VOF
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

#ifdef USE_HYPRE
	Pressure_hypre_destroy();
#endif

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

#ifdef VOF
	int init_hydrostatic_pressure = 0;

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

		case 7: // droplet on flat surface test case
			VoF_droplet_flat_plate(data_bag);
			break;

		case 8: // droplet on sphere test case
			VoF_init_droplet_on_sphere(data_bag);
			break;
		case 9:
			VoF_init_bilayer_at_4D(data_bag);
			break;

		case 10:
			VoF_init_droplet_on_sphere_theta(data_bag);
			break;

		case 11:
			VoF_init_all_heavy(data_bag);
			break;

		case 12:
			VoF_init_stationary_droplet_Francois(data_bag);
			break;

		case 13:
			VoF_init_meniscus_154deg(data_bag);
			break;

		case 14: // quasi-2D Rayleigh-Taylor instability in the x-y plane
			VoF_init_rayleigh_taylor_2d(data_bag);
			break;

		case 15: // roadmap Phase 1 planar 2D rising bubble
			if (!params->twod_cartesian_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 15 requires TWOD_CARTESIAN (planar 2D rising bubble).");
			VoF_init_planar_rising_bubble_2d(data_bag);
			break;

		case 16: // temporary legacy alias kept for existing thin-slab bubble inputs
			if (!params->twod_cartesian_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 16 requires TWOD_CARTESIAN.");
			VoF_init_planar_rising_bubble_2d(data_bag);
			break;

		case 17: // roadmap Phase 1 axisymmetric rising bubble on the axis
			if (!params->axisym_rz_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 17 requires AXISYM_RZ (axisymmetric rising bubble).");
			VoF_init_axisymmetric_rising_bubble_2d(data_bag);
			break;

		case 18: // Liu15 section 4.1 planar droplet on a static cylinder
			if (!params->twod_cartesian_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 18 requires TWOD_CARTESIAN.");
			VoF_init_planar_droplet_on_static_cylinder_theta(data_bag);
			break;

		case 19: // Liu17 section 6.2 axisymmetric droplet on a static sphere
			if (!params->axisym_rz_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 19 requires AXISYM_RZ.");
			VoF_init_axisymmetric_droplet_on_static_sphere_theta(data_bag);
			break;

			case 20: // Liu17 section 6.3 planar 2D sinking cylinder from water surface
				if (!params->twod_cartesian_enabled)
					CART3D_TWOD_ABORT(params,
						"init_type = 20 requires TWOD_CARTESIAN.");
				VoF_init_liu17_sinking_cylinder_2d(data_bag);
				break;

			case 21: // Liu17 section 6.4 axisymmetric sphere impact onto water
				if (!params->axisym_rz_enabled)
					CART3D_TWOD_ABORT(params,
						"init_type = 21 requires AXISYM_RZ.");
				VoF_init_liu17_axisymmetric_sphere_impact(data_bag);
				break;

			case 22: // Rajesh/Sauret 2026 section 3.1 capillary bridge
				if (!params->axisym_rz_enabled)
					CART3D_TWOD_ABORT(params,
						"init_type = 22 requires AXISYM_RZ.");
				VoF_init_axisymmetric_capillary_bridge_two_spheres(data_bag);
				break;

			case 23: // full-domain planar 2D pair of sinking cylinders
			if (!params->twod_cartesian_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 23 requires TWOD_CARTESIAN.");
			VoF_init_liu17_two_sinking_cylinders_2d(data_bag);
			break;

			case 24: // full-domain planar 2D triplet of sinking cylinders
			if (!params->twod_cartesian_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 24 requires TWOD_CARTESIAN.");
			VoF_init_liu17_three_sinking_cylinders_2d(data_bag);
			break;

			case 25: // Nguyen et al. APT 2021 section 5.3 static liquid bridge
				if (!params->axisym_rz_enabled)
					CART3D_TWOD_ABORT(params,
						"init_type = 25 requires AXISYM_RZ.");
				VoF_init_nguyen21_axisymmetric_static_bridge(data_bag);
			break;

			case 26: // Liu17 section 6.6 planar 2D self-assembly of floating cylinders
			if (!params->twod_cartesian_enabled)
				CART3D_TWOD_ABORT(params,
					"init_type = 26 requires TWOD_CARTESIAN.");
			VoF_init_liu17_self_assembly_floating_cylinders_2d(data_bag);
			break;


		    }
        // update ghost nodes / initialize the active interface representation
		#ifdef VOF_DIFFUSE
		VOF_DIFFUSE_set_boundary_values(data_bag->vof->F, data_bag);
		VOF_DIFFUSE_init(data_bag);
		#endif
		
		

		#ifdef VOF_PLIC
		// Reconstruct interface (PLIC-specific)
		VOF_set_boundary_values(data_bag->vof->F, data_bag);
        VOF_reconstruct_interface(data_bag);
		#endif

		VOF_update_density_viscosity(data_bag);


		// Bootstrap the balanced-force arrays from the initialized interface
			#ifdef SURFACE_TENSION
			#ifdef VOF_DIFFUSE
			VOF_DIFFUSE_compute_psi_LG(data_bag);
			#ifdef VOF_IBM
			VOF_DIFFUSE_extend_psi_LG_contact_angle(data_bag);
			#endif
			VOF_DIFFUSE_compute_f_sigma(data_bag);
			#endif

		#ifdef VOF_PLIC
		VoF_smoothing(data_bag);
		curvature_patel(data_bag);
		VOF_compute_f_sigma(data_bag);
		#endif
		// Copy new → old so first substep uses f_σ^0
		Array_copy_noghost(data_bag->vof->f_sigma_new_x, data_bag->vof->f_sigma_old_x, grid, params);
		Array_copy_noghost(data_bag->vof->f_sigma_new_y, data_bag->vof->f_sigma_old_y, grid, params);
		Array_copy_noghost(data_bag->vof->f_sigma_new_z, data_bag->vof->f_sigma_old_z, grid, params);
		#endif

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
