#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Communication.h"
#include "Grid.h"
#include "Display.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/******************************************************************************/
/*
 */
/******************************************************************************/
void Display_DA_3D_info(MAC_grid *grid, Parameters *params) {

	int ierr;
	int Is_g, Js_g, Ks_g;
	int Ie_g, Je_g, Ke_g;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int proc;
	int rank;

	// Start index of bottom-left-back corner on current processor including
	//ghost nodes
	Is_g = grid->L_Is;
	Js_g = grid->L_Js;
	Ks_g = grid->L_Ks;

	// End index of top-right-front corner on current processor including ghost
	// nodes
	Ie_g = grid->L_Ie;
	Je_g = grid->L_Je;
	Ke_g = grid->L_Ke;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	rank = params->rank;

	for (proc=0;proc<params->size;proc++) {
		if (proc==params->rank) {
			printf("\nProcessor Rank:%d Is:%d Js:%d Ks:%d\n", params->rank, Is, Js, Ks);
			printf("Processor Rank:%d Ie:%d Je:%d Ke:%d\n", params->rank, Ie, Je, Ke);
			printf("Processor Rank:%d Is_g:%d Js_g:%d Ks_g:%d\n", params->rank, Is_g, Js_g, Ks_g);
			printf("Processor Rank:%d Ie_g:%d Je_g:%d Ke_g:%d\n\n", params->rank, Ie_g, Je_g, Ke_g);
			fflush(stdout);
		}
		fflush(stdout);
		MPI_Barrier(PCW);
	} // for proc

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Display_2D_outflow(double **outflow, MAC_grid *grid, Parameters *params,
		char *q_name, char which_quantity) {

	int j, k;
	int Js, Ks;
	int Je, Ke;
	int proc;

	if (grid->G_Ie == grid->NX) {

		// start index on current processor
		Js = grid->G_Js;
		Ks = grid->G_Ks;

		// end index on current processor
		Je = grid->G_Je;
		Ke = grid->G_Ke;

		for (proc=0;proc<params->size;proc++) {
			if (proc==params->rank) {
				printf("***********************************\n");
				printf("Rank:%d Printing \"%s\"\n", params->rank, q_name);
				printf("StartIndex(k:%d,j:%d) EndIndex(k:%d,j:%d)\n", Ks, Js, Ke, Je);
				printf("q(j:y,k:z)\n");

				for (k=Ks; k<Ke; k++) {
					for (j=Js; j<Je; j++) {

						printf("q(%d,%d)=%f ", j, k, outflow[k][j]);

					} // for j
					printf("\n");
				} // for k
				fflush(stdout);
			}
			fflush(stdout);
			MPI_Barrier(PCW);
		} // for proc
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Display_2D_inflow(double **inflow, MAC_grid *grid, Parameters *params,
		char *q_name, char which_quantity) {

	int j, k;
	int Js, Ks;
	int Je, Ke;
	int proc;

	if (grid->G_Is == 0) {

		// start index on current processor
		Js = grid->G_Js;
		Ks = grid->G_Ks;

		// end index on current processor
		Je = grid->G_Je;
		Ke = grid->G_Ke;

		for (proc=0;proc<params->size;proc++) {
			if (proc==params->rank) {
				printf("***********************************\n");
				printf("Rank:%d Printing \"%s\"\n", params->rank, q_name);
				printf("StartIndex(k:%d,j:%d) EndIndex(k:%d,j:%d)\n", Ks, Js, Ke, Je);
				printf("q(j:y,k:z)\n");

				for (k=Ks; k<Ke; k++) {
					for (j=Js; j<Je; j++) {

						printf("q(%d,%d)=%f ", j, k, inflow[k][j]);

					} // for j
					printf("\n");
				} // for k
				fflush(stdout);
			}
			fflush(stdout);
			MPI_Barrier(PCW);
		} // for proc
	}
}




/******************************************************************************/
/*
 This function displays the parameters
 */
/******************************************************************************/
int Display_parameters(Parameters *params) {

	int iconc, NConc;
	printf("***************************************************************\n");
	printf("********************* COMPILED PARAMETERS *********************\n");
	printf("***************************************************************\n");
	printf("--------------------- BOUNDARY CONDITIONS Velocity ------------\n");
	//--------------------------------------------------------------------------
	// X-Boundaries
	//--------------------------------------------------------------------------
#ifdef XPERIODIC
	printf("X-Boundaries: ... periodic\n");
#endif
#ifdef LEFT_WALL_VELOCITY_NOSLIP
	printf("Left wall: ...... no-slip\n");
#elif defined LEFT_WALL_VELOCITY_FREESLIP
	printf("Left wall: ...... free-slip\n");
#elif defined LEFT_INFLOW
	printf("Left wall: ...... inflow\n");
#elif defined LEFT_OUTFLOW
	printf("Left wall: ...... outflow\n");
#endif
#ifdef RIGHT_WALL_VELOCITY_NOSLIP
	printf("Right wall: ..... no-slip\n");
#elif defined RIGHT_WALL_VELOCITY_FREESLIP
	printf("Right wall: ..... free-slip\n");
#elif defined RIGHT_INFLOW
	printf("Right wall: ..... inflow\n");
#elif defined RIGHT_OUTFLOW
	printf("Right wall: ..... outflow\n");
#endif
	//--------------------------------------------------------------------------
	// Y-Boundaries
	//--------------------------------------------------------------------------
#if defined TOP_WALL_VELOCITY && defined BOTTOM_WALL_VELOCITY
	double wall_vel = params->ubulk_target;
#else
	double wall_vel = 2.0 * params->ubulk_target;
#endif

#if defined YPERIODIC
	printf("Y-Boundaries: ... periodic\n");
#elif BOTTOM_WALL_VELOCITY_NOSLIP
	printf("Bottom wall: .... no-slip\n");
#elif defined BOTTOM_WALL_VELOCITY_FREESLIP
	printf("Bottom wall: .... free-slip\n");
#elif defined BOTTOM_WALL_VELOCITY
	printf("Bottom wall: .... U = %g\n", -wall_vel);
#else  // SCHUMANN
	printf("Bottom wall: .... Schumann\n");
#endif
#ifdef TOP_WALL_VELOCITY_NOSLIP
	printf("Top wall: ....... no-slip\n");
#elif defined TOP_WALL_VELOCITY_FREESLIP
	printf("Top wall: ....... free-slip\n");
#elif defined TOP_WALL_SCHUMANN
	printf("Top wall: ....... Schumann\n");
#elif defined TOP_WALL_VELOCITY
	printf("Top wall: ....... U = %g\n", wall_vel);
#endif
	//--------------------------------------------------------------------------
	// Z-Boundaries
	//--------------------------------------------------------------------------
#ifdef ZPERIODIC
	printf("Z-Boundaries: ... periodic\n");
#endif
#ifdef BACK_WALL_VELOCITY_NOSLIP
	printf("Back wall: ...... no-slip\n");
#elif defined BACK_WALL_VELOCITY_FREESLIP
	printf("Back wall: ...... free-slip\n");
#endif
#ifdef FRONT_WALL_VELOCITY_NOSLIP
	printf("Front wall: ..... no-slip\n");
#elif defined FRONT_WALL_VELOCITY_FREESLIP
	printf("Front wall: ..... free-slip\n");
#endif
	printf("\n");
	printf("--------------------- BOUNDARY CONDITIONS SCALAR(S) -------------\n");
	printf("----------------------A d/dx c + B c= C ----------------------------\n");
	#ifdef XPERIODIC_CONC
	printf("X-Boundaries: ... periodic\n");
	#elif defined CONC_SINGLE_BC_RB || defined CONC_SINGLE_BC_NODIFF
	printf("X-Boundaries: ... adiabatic A=1, B=0, C=0\n");
	#elif defined CONC_SINGLE_BC_ZERO
	printf("X-Boundaries: ... zeroed A=0, B=1, C=0\n");
	#else
	printf("X-Boundaries: ... set during runtime \n");
	#endif
    #ifdef YPERIODIC_CONC
	printf("Y-Boundaries: ... periodic\n");
	#elif defined CONC_SINGLE_BC_RB
	printf("Y-Boundaries: ... North A=0, B=1, C=0\n");
	printf("Y-Boundaries: ... South A=0, B=1, C=1\n");
	#elif  defined CONC_SINGLE_BC_NODIFF
	printf("Y-Boundaries: ... adiabatic A=1, B=0, C=0\n");
	#elif defined CONC_SINGLE_BC_ZERO
	printf("Y-Boundaries: ... zeroed A=0, B=1, C=0\n");
	#else
	printf("Y-Boundaries: ... set during runtime \n");
	#endif
	#ifdef ZPERIODIC_CONC
	printf("Z-Boundaries: ... periodic\n");
	#elif defined CONC_SINGLE_BC_RB && defined CONC_SINGLE_BC_NODIFF
	printf("Z-Boundaries: ... adiabatic A=1, B=0, C=0\n");
	#elif defined CONC_SINGLE_BC_ZERO
	printf("Z-Boundaries: ... zeroed A=0, B=1, C=0\n");
	#else
	printf("Z-Boundaries: ... set during runtime \n");
	#endif




	printf("------------------------- Fluid solver --------------------------\n");
#ifdef DEBUG
	printf("Debugging: ....... yes\n");
#else
	printf("Debugging: ....... no\n");
#endif
#ifdef CONSTANT_MASSFLUX
	printf("dp_dx: ........... variable to achieve u_bulk_target\n");
#else
	printf("dp_dx: ........... constant\n");
#endif
#ifdef FULLY_EXPLICIT
	printf("Viscous terms: ... fully-explicit\n");
#elif defined CG_SOLVE
	printf("Viscous terms: ... fully-implicit (CG)\n");
#elif defined BICG_SOLVE
	printf("Viscous terms: ... fully-implicit (BiCG)\n");
#else  // SEMI_IMPLICIT
	printf("Viscous terms: ... y-semi-implicit (FFT)\n");
#endif
	printf("\n");

#ifdef CONC
	printf("---------------------------- CONC -----------------------------\n");
	#ifdef CONC_CENTRAL
	printf("Advection stencil: ........ central\n");
	#elif defined CONC_BQUICK
	printf("Advection stencil: ........ bquick\n");
	#elif defined CONC_QUICK
	printf("Advection stencil: ........ quick\n");
	#else  // CONC_FTUPWIND
	printf("Advection stencil: ........ upwind\n");
	#endif
	#ifdef CONC_FULLY_IMPLICIT
	printf("Helmholtz solver: ........ CG\n");
	#endif
	#ifdef BOUSSINESQ
	printf("Boussinesq: ...... yes\n");
	#else
	printf("Boussinesq: ...... no\n");
	#endif
	#ifdef CONC_CLIP
	printf("Clip conc: ....... yes\n");
	#else
	printf("Clip conc: ....... no\n");
	#endif
	printf("\n");
#endif

#ifdef LAG_PARTICLE_RESOLVED
	printf("-------------------------- PARTICLE ---------------------------\n");
	#ifdef SUBSTEP
	printf("Temporal sub-stepping: ......... yes\n");
	#else
	printf("Temporal sub-stepping: ......... no\n");
	#endif
	#ifdef ACTM_TEST
	printf("Normal collision model: ........ ACM with substeps\n");
	#elif defined ACTM
	printf("Normal collision model: ........ ACM\n");
	#elif defined DEM
	printf("Normal collision model: ........ DEM\n");
	#else  // repulsive
	printf("Normal collision model: ........ repulsive potential\n");
	#endif
	#ifdef ATFM
		#ifdef ENABLE_ATFM_ROLLING
	printf("Tangential collision model: .... ATFM with rolling\n");
		#else
	printf("Tangential collision model: .... ATFM without rolling\n");
		#endif
		Display_throw_warning("It is inadvised to use the ATFM model", params);
	#elif defined LIN_TAN
	printf("Tangential collision model: .... linear spring-dashpot\n");
	#else
	printf("Tangential collision model: .... none\n");
	#endif
	#ifdef LUBRICATION_NORMAL
	printf("Lubrication model (normal): .... yes\n");
	#else
	printf("Lubrication model (normal): .... no\n");
	#endif
	#ifdef LUBRICATION_TANGENTIAL
	printf("Lubrication model (tangential):  yes\n");
		#ifndef DRY_PARTICLES
		Display_throw_warning("It is inadvised to use tangential lubrication with the fluid solver", params);
		#endif
	#else
	printf("Lubrication model (tangential):  no\n");
	#endif
	#ifdef ELECTROSTATIC_REPULSION
	printf("Electrostatic Repulsion: ....... yes\n");
	#else
	printf("Electrostatic Repulsion: ....... no\n");
	#endif
	#ifdef DRY_COLLISION
	printf("Fluid forces during contact: ... no\n");
	#else
	printf("Fluid forces during contact: ... yes\n");
	#endif
	#ifdef ROUGH_COLLISION
	printf("Start collision contact: ....... at surface roughness\n");
	#else
	printf("Start collision contact: ....... at particle radius\n");
	#endif
	#ifdef LAG_MARKER_FLAG
	printf("Turn off competing Lag markers: ........ yes\n");
	#else
	printf("Turn off competing Lag markers: ........ no\n");
	#endif
	#ifdef LAG_MARKER_PRIORITY
	printf("Use nearby Lag markers: ........ yes\n");
	#else
	printf("Use nearby Lag markers: ........ no\n");
	#endif
	#ifdef SQUIRMER_SWIMMER
	printf("Squirmer model: ................ yes\n");
	#else
	printf("Squirmer model: ................ no\n");
	#endif
	#ifdef COHESION
	printf("Cohesion model: ................ yes\n");
	printf("Bond number:.................... %g\n", params->Co);
	#else
	printf("Cohesion model: ................ no\n");
	#endif
	#ifdef RETRACTION
	printf("Lag markers retracted RETRACTION*dx: .. yes\n");
	#else
	printf("Lag markers retracted xx*dx: .. no\n");
	#endif
	#ifdef IBM_SCALAR
	printf("IBM for scalar fields .......... yes\n");
	#else
	printf("IBM for scalar fields .......... no\n");
	#endif
	#ifdef VOF_SCALAR
	printf("VOF for scalar fields .......... yes\n");
	#else
	printf("VOF for scalar fields .......... no\n");
	#endif
	printf("\n");
#endif

	// TODO: LES

	// TODO: RANS


	printf("***************************************************************\n");
	printf("********************** INPUT PARAMETERS ***********************\n");
	printf("***************************************************************\n");
	printf("-------------------------- GEOMETRY ---------------------------\n");
	printf("xmin: %8.6f, xmax: %8.6f, Lx: %8.6f\n", params->xmin, params->xmax, params->Lx);
	printf("ymin: %8.6f, ymax: %8.6f, Ly: %8.6f\n", params->ymin, params->ymax, params->Ly);
	printf("zmin: %8.6f, zmax: %8.6f, Lz: %8.6f\n", params->zmin, params->zmax, params->Lz);
	printf("\n");

	printf("---------------------------- GRID -----------------------------\n");
	printf("NX: %d, NY: %d, NZ: %d\n", params->NXM, params->NYM, params->NZM);
	Display_flag(params->ImportGridFromFile, "Import Grid From File");
	printf("\n");

	printf("------------------------- SIMULATION --------------------------\n");
	printf("Time max: ............... %g\n", params->time_max);
	printf("Output time interval: ... %g\n", params->output_time_interval);
	printf("Constant dt: ............ ");
	if (params->constant_dt) printf("yes\n");
	else printf("no\n");
	printf("Default dt: ............. %g\n", params->default_dt);
	printf("Maximum dt: ............. %g\n", params->max_dt);
	printf("CFL number: ............. %g\n", params -> cfl);
	Display_flag(params->resume, "Resume");
	printf("\n");

	printf("---------------------------- FLOW -----------------------------\n");
	printf("Re: ............................ %g\n", params->Re);
#ifdef CONSTANT_MASSFLUX
	printf("Target bulk velocity: .......... %g\n", params -> ubulk_target);
	printf("Initial pressure gradient: ..... %g\n", params -> dp_dx);
#else
	printf("Pressure gradient: ............. %g\n", params -> dp_dx);
	printf("Initial bulk velocity: ......... %g\n", params -> ubulk_target);
#endif
	printf("Initial velocity profile: ...... ");
	if (params->vel_init_type == VEL_INIT_ZERO)
		printf("zero flow\n");
	else if (params->vel_init_type == VEL_INIT_UNIFORM)
		printf("uniform flow\n");
	else if (params->vel_init_type == VEL_INIT_LINEAR)
		printf("linear shear flow\n");
	else if (params->vel_init_type == VEL_INIT_POISEUILLE)
		printf("Poiseuille (parabolic) flow\n");
	else if (params->vel_init_type == VEL_INIT_ROT_SHEAR)
		printf("rotational shear flow\n");
	else if (params->vel_init_type == VEL_INIT_PRECURSOR)
		printf("Precursor simulation inflow\n");
	printf("Initial velocity profile y0: ... %g\n", params -> vel_init_y0);
	printf("\n");

#ifdef CONC
	printf("---------------------------- CONC -----------------------------\n");
	NConc = params->NConc;
	for (iconc = 0; iconc < NConc; iconc++) {
		printf("iconc: %d\n", iconc);
		printf("conc_init_type: %d\n", params->conc_init_type[iconc]);
		printf("Pe: %f, V_s0: %f, richardson: %f\n", params->Pe[iconc],
			   params->V_s0[iconc], params->richardson[iconc]);
	}

	printf("\n");

	printf("---------------------------- LOCK -----------------------------\n");
	printf("Lock front x: %f\n", params->x_fr);
	printf("Lock front y: %f\n", params->y_fr);
	printf("Lock front z: %f\n", params->z_fr);
	printf("\n");
#endif

#ifdef LAG_PARTICLE_RESOLVED
	printf("-------------------------- PARTICLE ---------------------------\n");
	printf("Extra forcing loop iterations: ...... %d\n", params->N_forcing_loops);
	printf("Particle/fluid density ratio: ....... %g\n", params->rho_s);
	printf("Gravity: ............................ {%g, %g, %g}\n", params->grav[0], params->grav[1], params->grav[2]);
	printf("Fluid timesteps per collision: ...... %d\n", params->Ndt_coll);
	printf("Particle restitution coefficient: ... %g\n", params->e_dry_particles);
	printf("Wall restitution coefficient: ....... %g\n", params->e_dry_wall);
	printf("Coefficient of kinetic friction: .... %g\n", params->mu_k);
	printf("Coefficient of static friction: ..... %g\n", params->mu_s);
	printf("Poisson's ratio: .................... %g\n", params->PoissonsRatio);
	printf("Roughness length: ................... %g * R_p\n", params->roughness);
	#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL
	printf("Lubrication range: .................. %g * h\n", params->lub_range);
	#endif
	#ifdef ELECTROSTATIC_REPULSION
	printf("Electrostatic repulsion range: ...... %g * h\n", params->erm_range);
	printf("Electrostatic repulsion coeff.: ..... %g\n", params->F_erm);
	#endif
	#ifdef SQUIRMER_SWIMMER
	printf("Swimmer B1: ......................... %g\n", params->B1);
	printf("Swimmer B2: ......................... %g\n", params->B2);
	if (params->target_on == 0) {
		printf("Swimmer target: ..................... off\n");
	}
	else if (params->target_on == 1) {
		printf("Coordinate of swimmer target: ....... {%g, %g, %g}\n", params->target_coord[0], params->target_coord[1], params->target_coord[2]);
	}
	else {
		printf("Swimmer target: ..................... !INVALID!\n");
		return -1;
	}

	#endif
	printf("\n");
#endif

#ifdef LES
	// TODO: LES
#endif

	printf("--------------------------- LSOLVE ----------------------------\n");
	printf("CG error tolerance: ..... %g\n", params->CG_ETOL);
	printf("CG max iterations: ...... %d\n", params->CG_MAXIT);
	printf("\n");

	printf("--------------------------- OUTPUT ----------------------------\n");
	Display_flag(params->ave_height_output, "conc average height output");
	Display_flag(params->front_location_output, "front location output");
	Display_flag(params->front_speed_output, "front speed output");
	Display_flag(params->susp_mass_output, "suspended mass output");
	Display_flag(params->sedim_rate_output, "sediment rate output");
	Display_flag(params->energies_output, "energies_output");
	Display_flag(params->shear_stress_output, "shear_stress_output");
	//Display_flag(params->conc_output_dump, "dump_conc");
	//Display_flag(params->conc_output_deposit_height, "deposit height output");
	printf("***************************************************************\n");
	printf("\n");

	return 0;
}




/******************************************************************************/
/*
 This function displays the results of the time required by various parts of the
 simulation.
 */
/******************************************************************************/
void Display_time_results(Cart3d_bag *data_bag) {

	Timer *timer = data_bag -> timer;

	fflush(stdout);
	fflush(stderr);
	MPI_Barrier(PCW);

	if (data_bag->params->rank == 0) {
		printf("--------------------------------");
		printf(" Timing results --------------------------------\n");
		printf("\n");

		printf("Simulation:\n");
		printf("\tTotal simulation time: .............. %f\n", timer->Wtime_total);
		printf("\tInitialization: ..................... %f\n", timer->Wtime_init);
		printf("\tIntegrating equations of motion: .... %f\n", timer->Wtime_intEOM);
		printf("\tCommunicating 3-D flow variables: ... %f\n", timer->Wtime_comm_3D);
		printf("\tWriting data to HDF5 files: ......... %f\n", timer->Wtime_output);
		printf("\tCalculating and writing 2D Statistics:%f\n", timer->Wtime_output_2d);
		printf("\tEvaluating CFL condition: ........... %f\n", timer->Wtime_cfl);
		printf("\n");

		printf("Momentum:\n");
		printf("\tSolving convective, explicit viscous terms: ... %f\n", timer->Wtime_vel_convective);
		printf("\tSolving right-hand side: ...................... %f\n", timer->Wtime_vel_rhs);
		printf("\tSolving implicit equation: .................... %f\n", timer->Wtime_vel_solve);
		printf("\tCalculating cell-centered values: ............. %f\n", timer->Wtime_vel_cell_center);
		printf("\tCalculating boundary cells: ................... %f\n", timer->Wtime_vel_boundaries);
		printf("\n");

		printf("Pressure:\n");
		printf("\tSolving right-hand side: .......... %f\n", timer->Wtime_p_rhs);
		printf("\tSolving Poisson equation: ......... %f\n", timer->Wtime_p_solve);
		printf("\tProjecting pressure, velocity: .... %f\n", timer->Wtime_p_project);
		printf("\tEvaluating velocity divergence: ... %f\n", timer->Wtime_p_project);
		printf("\n");

#ifdef CONC
		printf("Concentration:\n");
		printf("\tTotal time: ................................... %f\n", timer->Wtime_c_total);
		printf("\tSolving convective, explicit viscous terms: ... %f\n", timer->Wtime_c_convective);
		printf("\tSolving right-hand side: ...................... %f\n", timer->Wtime_c_rhs);
		printf("\tSolving implicit equation: .................... %f\n", timer->Wtime_c_solve);
		printf("\tEvaluating out-of-bounds terms: ............... %f\n", timer->Wtime_c_outofbounds);
		printf("\n");
#endif

#ifdef LES
		printf("LES:\n");
		printf("\tTotal time: ... %f\n", timer->Wtime_sgs_total);
		printf("\tFilter: ....... %f\n", timer->Wtime_sgs_filter);
		printf("\tDynamic: ...... %f\n", timer->Wtime_sgs_dynamic);
		printf("\n");
#endif

#ifdef RANS
		printf("RANS:\n");
		printf("\tTotal time: ................................... %f\n", timer->Wtime_rans_total);
		printf("\tSolving convective, explicit viscous terms: ... %f\n", timer->Wtime_rans_convective);
		printf("\tSolving right-hand side: ...................... %f\n", timer->Wtime_rans_rhs);
		printf("\tSolving implicit equation: .................... %f\n", timer->Wtime_rans_solve);
		printf("\n");
#endif

#ifdef LAG_PARTICLE_RESOLVED
		printf("Resolved Particles:\n");
		printf("\tTotal time: ................... %f\n", timer->Wtime_particle_total);
		printf("\tCommunication: ................ %f\n", timer->Wtime_particle_comm);
		printf("\tForcing Lagrangian markers: ... %f\n", timer->Wtime_particle_forc);
		printf("\tIntegrating velocities: ....... %f\n", timer->Wtime_particle_int);
		printf("\tCollisions: ................... %f\n", timer->Wtime_particle_coll);
		printf("\n");
#endif

		printf("----------------------------------------");
		printf("----------------------------------------\n");
	}

	fflush(stdout);
	fflush(stderr);
	MPI_Barrier(PCW);
}




/******************************************************************************/
/*
 This function printd the status of a flag
 */
/******************************************************************************/
void Display_flag(int flag, char *name) {

	if (flag) {
		printf("%s: YES\n", name);
	} else {
		printf("%s: NO\n", name);
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Display_point(PointType *p, char *name) {

	printf("Printing point:%s ", name);
	printf("P(x,y,z)=(%f,%f,%f)\n", p->x, p->y, p->z);
}




/******************************************************************************/
/*
 This function displays the information for the immersed node
 */
/******************************************************************************/
void Display_immersed_node(ImmersedNode *ib_node) {

	char s[100];
	int g;
	printf("Printing information for the immersed node\n");
	Display_point(&ib_node->im_point, "immersed");
	Display_point(&ib_node->boundary_point, "boundary point");
	Display_point(&ib_node->intersection_point, "intersection point");

}




/******************************************************************************/
/*
 Displays any grid info
 */
/******************************************************************************/
void Display_grid(double *data, int N, char *name) {

	int i;
	printf("Display.c/ *************************************/ \n");
	for (i=0; i<N; i++) {
		printf("Display.c/ gridinfo-%s(%d)=%f\n", name, i, data[i]);
	} // for i
}



/******************************************************************************/
/*
 Displays inputted message if the program reaches a certain point on all
 processors.
 */
/******************************************************************************/
void Display_progress(Parameters *params, char *message) {

	MPI_Barrier(MPI_COMM_WORLD);
	if (params->rank == 0)
		printf("%s", message);
	fflush(stdout);
}




/******************************************************************************/
/*
 * Asserts that p_list->state matches state.  If not, it prints an error message
 * and terminates the program.
 */
/******************************************************************************/
void Display_assert_list_state(Particle_list *p_list, int state,
		Parameters *params, Debug_trace *dtrace) {

	if (p_list->state != state) {
		char message[100];
		char name_list_type[50], name_list_state[50], name_state[50];

		Display_print_list_type(name_list_type, p_list->type);
		Display_print_list_name(name_state, state);
		Display_print_list_name(name_list_state, p_list->state);

		sprintf(message, "Invalid input linked list state '%s'\n"
		        "'p_%s' should contain '%s' particles",
				name_list_state, name_list_type, name_state);
		Display_throw_error(message, params, dtrace);
	}
}



/******************************************************************************/
/*
 * Prints the Particle linked list state 'state' into the character array 'name'
 */
/******************************************************************************/
void Display_print_list_name(char *name, int state) {

	if (state == LIST_STATE_LOCAL)
		sprintf(name, "local");
	else if (state == LIST_STATE_FOREIGN)
		sprintf(name, "foreign");
	else if (state == LIST_STATE_BOTH)
		sprintf(name, "both local and foreign");
	else if (state == LIST_STATE_EDGE)
		sprintf(name, "edge");
	else
		sprintf(name, "UNSUPPORTED STATE");
}




/******************************************************************************/
/*
 * Prints the Particle linked list type 'type' into the character array 'name'
 */
/******************************************************************************/
void Display_print_list_type(char *name, int type) {

	if (type == MOBILE)
		sprintf(name, "mobile");
	else if (type == FIXED)
		sprintf(name, "fixed");
	else
		sprintf(name, "UNSUPPORTED TYPE");
}




/******************************************************************************/
/*
 * Prints list of traced function calls
 */
/******************************************************************************/
void Display_print_dtrace_list(char *message_out, Debug_trace *dtrace) {

	message_out += sprintf(message_out, "Error called by %s (%s:%d)\n",
	                       dtrace->function, dtrace->file, dtrace->line);
	dtrace = dtrace -> parent;
	while (dtrace !=NULL) {
		if (dtrace->parent != NULL) {
			message_out += sprintf(message_out, "             at %s (%s:%d)\n",
			                       dtrace->function, dtrace->file, dtrace->line);
		}
		dtrace = dtrace -> parent;
	}
	message_out += sprintf(message_out, "\n");
}




/******************************************************************************/
/*
 Prints error that is executed from all processors
 */
/******************************************************************************/
void Display_assert_error(int status, const char *message_in,
		Parameters *params, Debug_trace *dtrace) {

	int *status_vec, status_global, i;
	char line1[500], *line2, temp[50], *message;

	// Add function info to error message
	Display_print_dtrace_list(line1, dtrace);
	message = (char *)malloc( (strlen(message_in) + strlen(line1) +
							   50 * params->size) * sizeof(char) );
	strcpy(message, line1);

	// Gather error statuses to all processors
	status_vec = (int *)malloc(params->size * sizeof(int));
	MPI_Allgather(&status, 1, MPI_INT, status_vec, 1, MPI_INT, PCW);

	// Add error statuses and owner processor info to error message
	status_global = 0;
	for (i = 0; i < params->size; i++) {
		if (status_vec[i] < 0) {

			sprintf(temp, "    Processor %d: ierr = %d\n", i, status_vec[i]);
			strcat(message, temp);
			status_global = -1;
		}
	}

	// Add other error message info to message
	strcat(message, message_in);

	// If error exists, display error message and abort
	if (status_global < 0) {

		fflush(stdout);
		MPI_Barrier(PCW);

		if (params -> rank == 0)
			Display_error(message);

		fflush(stdout);
		MPI_Barrier(PCW);

		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	free(status_vec);
	free(message);
}




/******************************************************************************/
/*
 Prints error that is executed from all processors
 */
/******************************************************************************/
void Display_throw_warning(const char *message, Parameters *params) {

	if (params -> rank == 0)
		Display_error(message);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
}




/******************************************************************************/
/*
 * Prints error that is executed from all processors
 */
/******************************************************************************/
void Display_throw_error(const char *message_in, Parameters *params,
		Debug_trace *dtrace) {

	char line1[500], *message;

	// Print nested list of function calls from 'dtrace'
	Display_print_dtrace_list(line1, dtrace);

	// Add 'message_in'
	message = (char *)malloc( (strlen(message_in) + strlen(line1) + 1)
	                         * sizeof(char) );
	strcpy(message, line1);
	strcat(message, message_in);

	// Format and display error message
	if (params -> rank == 0)
		Display_error(message);
	fflush(stdout);

	free(message);
	MPI_Barrier(PCW);
	MPI_Finalize();
	exit(EXIT_FAILURE);
}




/******************************************************************************/
/*
 Prints error message and stops program
 */
/******************************************************************************/
void Display_error(const char *message_in) {

	int i, j;

//	const int STR_LENGTH = 100;
	const int MAX_LINES  = 50;
//	char line_1[STR_LENGTH], line_2[STR_LENGTH];
	int line_start[MAX_LINES], line_length[MAX_LINES];

	int line_count, max_line_length;
	char *error_str, *message, *line;

//	sprintf(line_1, "!!! Error on Processor %d", params->rank);
//	sprintf(line_2, "!!!       in function %s:", function);
//
//	max_line_length = max(strlen(line_1) + 4, strlen(line_2) + 4);
	max_line_length = 0;

	message = (char *)malloc( (strlen(message_in) + 1) * sizeof(char) );
	strcpy(message, message_in);

	// Find number of lines in message and the maximum line length
	line_start[0] = 0;
	line_count = 0;
	for (i = 0; i <= strlen(message_in); i++) {
		if (message[i] == '\n' || message[i] == '\0') {
			message[i] = '\0';
			line_length[line_count] = i - line_start[line_count];
			max_line_length = max(line_length[line_count] + 8, max_line_length);
			line_count++;
			line_start[line_count] = i + 1;
		}
	}

	error_str = (char *)malloc( ( (max_line_length + 1) * (line_count + 4) + 1 )
	                           * sizeof(char) );

	// Line 0
	error_str[0] = '\0';
	for (i = 0; i < max_line_length; i++) {
		strcat(error_str, "*");
	}
	strcat(error_str, "\n");
/*
	// Line 1
	strcat(error_str, line_1);
	for (i = 0; i < max_line_length - strlen(line_1) - 4; i++) {
		strcat(error_str, " ");
	}
	strcat(error_str, " !!!\n");

	// Line 2
	strcat(error_str, line_2);
	for (i = 0; i < max_line_length - strlen(line_2) - 4; i++) {
		strcat(error_str, " ");
	}
	strcat(error_str, " !!!\n");
*/
	// 'message' lines
	for (j = 0; j < line_count; j++) {

		line = &message[line_start[j]];
		strcat(error_str, "!!! ");
		strcat(error_str, line);

		for (i = 0; i < max_line_length - line_length[j] - 8; i++) {
			strcat(error_str, " ");
		}
		strcat(error_str, " !!!\n");
	}

	// Last line
	for (i = 0; i < max_line_length; i++) {
		strcat(error_str, "*");
	}
	strcat(error_str, "\n");

	// Display error
	printf("%s", error_str);

	free(error_str);
	free(message);
}




/******************************************************************************/
/*
 * Create a child in the debug_trace list 'dtrace'.  This function should be
 * called as an input to a function that is being traced, e.g.:
 *
 *     Temporal_int_rk3(data_bag, Display_create_debug_trace_child(dtrace,
 *         "Temporal_int_rk3", __FILE__, __LINE__));
 *
 * I've also created the macro DTRACE to do most of this for you, e.g.:
 *
 *     Temporal_int_rk3(data_bag, DTRACE("Temporal_int_rk3"));
 *
 * Inputs:
 *     dtrace - pointer to Debug_trace list
 *     function - string containing name of function being called
 *     file - string containing filename (use macro '__FILE__')
 *     line - integer containing line number (use macro '__LINE')
 */
/******************************************************************************/
Debug_trace *Display_create_debug_trace_child(Debug_trace *dtrace,
		const char *function, const char *file, const int line) {

	// Create dtrace child based on input arguments
	Debug_trace *child = (Debug_trace *)malloc(sizeof(Debug_trace));
	sprintf(child->function, "%s", function);
	sprintf(child->file, "%s", file);
	child->line = line;

	// Kill sibling and all his children to prevent the apocalpyse (i.e. memory
	// leaks)
	Debug_trace *sibling, *sibling_child;
	sibling = dtrace -> child;
	while (sibling != NULL) {
		sibling_child = sibling->child;
		free(sibling);
		sibling = sibling_child;
	}

	// Establish inheritance for the new (alive) child
	dtrace -> child = child;
	child->parent = dtrace;
	child->child = NULL;

	return child;
}
