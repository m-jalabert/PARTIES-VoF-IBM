#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Collision.h"
#include "Communication.h"
#include "Display.h"
#include "Interpolate.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "Particle.h"
#include "Rotate.h"
#include "VolumeFraction.h"
#include "VOF-CICSAM.h"

#include "qr_solve.h"
#include "r8lib.h"


//------------------------------------------------------------------------------
// Logic test for a local Lagrangian point.  Points that are within the
// processor boundaries will be added to the total force on the particle.
//------------------------------------------------------------------------------
#define LOCAL_POINT(x, y, z) \
    (x >= xu[Is] && x < xu[Ie] && \
     y >= yv[Js] && y < yv[Je] && \
     z >= zw[Ks] && z < zw[Ke])

#define LOCAL_POINT_CONC(x, y, z) \
    (x >= xc[Is] && x < xc[Ie] && \
     y >= yc[Js] && y < yc[Je] && \
     z >= zc[Ks] && z < zc[Ke])


//------------------------------------------------------------------------------
// Logic test for a nearby Lagrangian point.  Points that are within 'range'
// of the processor boundaries will affect local velocity nodes.
//------------------------------------------------------------------------------
#define NEAR_POINT(x, y, z, range) \
    (xu[Is] - x < range && x - xu[Ie] < range && \
     yv[Js] - y < range && y - yv[Je] < range && \
     zw[Ks] - z < range && z - zw[Ke] < range)



/******************************************************************************/
/*
 * Allocates storage for Lagrangian structure and initializes particle data
 */
/******************************************************************************/
Lagrangian *Lagrangian_create(MAC_grid *grid, Parameters *params) {

	Lagrangian *lag = (Lagrangian *)malloc(sizeof(Lagrangian));

	Particle_list *p_mobile_list = (Particle_list *)malloc(sizeof(Particle_list));
	Particle_list *p_fixed_list  = (Particle_list *)malloc(sizeof(Particle_list));
	lag -> p_mobile_list = p_mobile_list;
	lag -> p_fixed_list  = p_fixed_list;

#ifdef PARTICLE_RELEASE
    Particle_list *p_release_list  = (Particle_list *)malloc(sizeof(Particle_list));
    lag -> p_release_list = p_release_list;
#endif

	lag -> send_part = NULL;
	lag -> recv_part = NULL;
	lag -> send_coll = NULL;
	lag -> recv_coll = NULL;

	lag -> send_part_size = 0;
	lag -> recv_part_size = 0;
	lag -> send_coll_size = 0;
	lag -> recv_coll_size = 0;

	lag -> send_Nc = (int *)malloc(sizeof(int));
	lag -> recv_Nc = (int *)malloc(sizeof(int));
	lag -> send_Nc[0] = 0;
	lag -> recv_Nc[0] = 0;

	lag -> send_Nc_size = 1;
	lag -> recv_Nc_size = 1;

	lag -> temp = Memory_allocate_flow_variable(grid, params);
	lag -> ng_temp = Memory_allocate_noghost_variable(grid, params);

#ifdef POST_PROCESS
	lag -> ng_fx_IBM = Memory_allocate_noghost_variable(grid, params);
	lag -> ng_fy_IBM = Memory_allocate_noghost_variable(grid, params);
#endif

	lag -> ng_vfu = Memory_allocate_noghost_variable(grid, params);
	lag -> ng_vfv = Memory_allocate_noghost_variable(grid, params);
	lag -> ng_vfw = Memory_allocate_noghost_variable(grid, params);
#if defined LAG_PARTICLE_RESOLVED
	lag -> ng_vfc = Memory_allocate_noghost_variable(grid, params);
#endif
	lag -> ng_vfz = Memory_allocate_noghost_variable(grid, params);
#ifdef VOF_SCALAR
	lag -> vfu = Memory_allocate_flow_variable(grid, params);
	lag -> vfv = Memory_allocate_flow_variable(grid, params);
	lag -> vfw = Memory_allocate_flow_variable(grid, params);

    #if defined(VOF_SMOOTH_VELO) || defined(VOF_PP_VFPRIME)
	lag -> vfu_prime = Memory_allocate_flow_variable(grid, params);
	lag -> vfv_prime = Memory_allocate_flow_variable(grid, params);
	lag -> vfw_prime = Memory_allocate_flow_variable(grid, params);

	#endif


#endif

	return lag;
}




/******************************************************************************/
/*
 * Frees storage for Lagrangian structure
 */
/******************************************************************************/
void Lagrangian_destroy(Lagrangian *lag, MAC_grid *grid, Parameters *params) {

	Particle_list_destroy(lag -> p_mobile_list);
	Particle_list_destroy(lag -> p_fixed_list);

#ifdef PARTICLE_RELEASE
    Particle_list_destroy(lag -> p_release_list);
#endif

	free(lag -> send_part);
	free(lag -> recv_part);
	free(lag -> send_coll);
	free(lag -> recv_coll);
	free(lag -> send_Nc);
	free(lag -> recv_Nc);
	Memory_free_flow_variable(grid, params, lag -> temp);

	Memory_free_noghost_variable(grid, params, lag -> ng_temp);
	free(lag -> Temp_L);
	free(lag -> Temp_H);
	#ifdef VOF_IBM
	free(lag -> Temp_L_rho);
	#endif


#ifdef POST_PROCESS
	Memory_free_noghost_variable(grid, params, lag -> ng_fx_IBM);
	Memory_free_noghost_variable(grid, params, lag -> ng_fy_IBM);
#endif

	Memory_free_noghost_variable(grid, params, lag->ng_vfu);
	Memory_free_noghost_variable(grid, params, lag->ng_vfv);
	Memory_free_noghost_variable(grid, params, lag->ng_vfw);

#if defined LAG_PARTICLE_RESOLVED
	Memory_free_noghost_variable(grid, params, lag->ng_vfc);
#endif
	Memory_free_noghost_variable(grid, params, lag->ng_vfz);
#ifdef VOF_SCALAR
	Memory_free_flow_variable(grid, params, lag->vfu);
	Memory_free_flow_variable(grid, params, lag->vfv);
	Memory_free_flow_variable(grid, params, lag->vfw);

    #if defined(VOF_SMOOTH_VELO) || defined(VOF_PP_VFPRIME)
	Memory_free_flow_variable(grid, params, lag->vfu_prime);
	Memory_free_flow_variable(grid, params, lag->vfv_prime);
	Memory_free_flow_variable(grid, params, lag->vfw_prime);

	#endif


#endif

	free(lag);
}


/******************************************************************************/
/*
 * Predicts location of particles for evaluating fluid-particle forces
 * Should have only local particles when called
 * Has only local particles afterwards
 */
/******************************************************************************/
void Lagrangian_advect_particles(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Evaluate hydrodynamic forces acting on particles
	Lagrangian_evaluate_fluid_forces(data_bag, DTRACE("Lagrangian_evaluate_fluid_forces"));

#ifdef SUBSTEP
	// Store RK stage number and timestep of fluid solver
	int rk;
	const double BET2[] = {BETA2};
	int which_stage_top = params->which_stage;
	double time_top = params->time;
	double dt_top = params->dt;

	// Number of substeps to take, depending on RK stage
	int N_substeps[] = {8, 2, 5};
	params -> dt = dt_top / 15.0;

	// Solve particle motion based on collisions using substeps
	for (i = 0; i < N_substeps[which_stage_top]; i++) {
//		ParticleOutput_dat(*p_mobile_ptr, grid, params);
		for (rk = 0; rk < 3; rk++) {
			params -> which_stage = rk;
			Lagrangian_integrate_particle_motion(data_bag, DTRACE("Lagrangian_advect_particles"));
			params->time += params->dt * BET2[rk];
		}
	}

	// Reset RK stage number and timestep
	params->which_stage = which_stage_top;
	params->time = time_top;
	params->dt = dt_top;

#else  // not SUBSTEP
	// Evaluate collisions and solve particle motion
	Lagrangian_integrate_particle_motion(data_bag, DTRACE("Lagrangian_advect_particles"));
#endif


}




/******************************************************************************/
/*
 * Collects and evaluates all forces acting on particle due to the IBM and
 * gravity
 *
 * Input linked list state: local and foreign
 * Output linked list state: local
 */
/******************************************************************************/
void Lagrangian_evaluate_fluid_forces(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i,j;
	int ID;
	Particle *p;
	Particle_list *p_list_foreign;
	double *F, *Fc, *T, *Int_U, *Int_Omega, *Int_U_old, *Int_Omega_old;
	double *F_IBM, *T_IBM, *F_rigid, *T_rigid, *X;
	double ***u   = data_bag -> u -> data;
	double ***v   = data_bag -> v -> data;
	double ***w   = data_bag -> w -> data;
	double Mass=0; 			// To be used to calculate the distributed mass for the pressure imposed couette flow set up
	double F_local_sum=0;	// To be used to calculate the local force (on each processor) for the pressure imposed couette flow set up



	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	const double BET[] = {BETA};
	double bet = BET[params -> which_stage];
	double dt  = params -> dt;
	double idt2beta = 1.0 / ( 2.0 * bet * dt );

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	#ifdef VOF_IBM
	//==========================================================================
    // 0. Compute Capillary Forces (CCF)
    //==========================================================================
    // We do this BEFORE collecting forces so that the CCF force calculated
    // on each processor is added to the local particle copy, then summed globally.

    // Loop over mobile particles
    p = p_mobile_list -> start;
    while (p != NULL) {
        // Calculate CCF for this particle
        VOF_accumulate_solid_capillary_force(p, data_bag);
        p = p -> next;
    }

    // Loop over fixed particles (if they interact with capillary interface)
    p = p_fixed_list -> start;
    while (p != NULL) {
         VOF_accumulate_solid_capillary_force(p, data_bag);
         p = p -> next;
    }



	//     // Compute spurious CSF force inside solid (to be subtracted)
    // {
    //     Particle *p = p_mobile_list->start;
    //     while (p != NULL) {
    //         VOF_integrate_CSF_over_solid(p, data_bag);
    //         p = p->next;
    //     }
    //     p = p_fixed_list->start;
    //     while (p != NULL) {
    //         VOF_integrate_CSF_over_solid(p, data_bag);
    //         p = p->next;
    //     }
    // }
	#endif

	//--------------------------------------------------------------------------
	// Collect hydrodynamic forces
	//--------------------------------------------------------------------------


	//#ifndef ONE_WAY
		p_list_foreign = Particle_list_foreign_create(p_mobile_list, data_bag, DTRACE("Particle_list_foreign_create"));
		Particle_MPI_update(p_list_foreign, data_bag, DTRACE("Particle_MPI_update"));
		Lagrangian_collect_forces(p_mobile_list, p_list_foreign, LAG_COLLECT_HYDRO, params, DTRACE("Lagrangian_collect_forces"));
		Particle_list_destroy(p_list_foreign);

		// Remove foreign particles from p_fixed for advecting particles
		p_list_foreign = Particle_list_foreign_create(p_fixed_list, data_bag, DTRACE("Particle_list_foreign_create"));
		Particle_MPI_update(p_list_foreign, data_bag, DTRACE("Particle_MPI_update"));
		Lagrangian_collect_forces(p_fixed_list, p_list_foreign, LAG_COLLECT_HYDRO, params, DTRACE("Lagrangian_collect_forces"));
		Particle_list_destroy(p_list_foreign);
	//#endif


	//--------------------------------------------------------------------------
	// Evaluate fluid forces for mobile particles
	//--------------------------------------------------------------------------
	p = p_mobile_list -> start;
	while (p != NULL) {

		F = p -> F;
		T = p -> T;
		Int_U     = p -> Int_U;
		Int_U_old = p -> Int_U_old;
		Int_Omega     = p -> Int_Omega;
		Int_Omega_old = p -> Int_Omega_old;
		X = p->X;


		F_IBM = p -> F_IBM;
		T_IBM = p -> T_IBM;
		F_rigid = p -> F_rigid;
		T_rigid = p -> T_rigid;


		/* ======================================================================
		 * SYMMETRY CORRECTION FOR HALF-DOMAIN
		 * ====================================================================== 
		 * By applying this to the raw single-stage variables at the very top,
		 * ALL downstream variables (F_IBM, F_CCF_cum, F_rigid) automatically 
		 * inherit the correct full-domain values without RK compounding errors!
		 * ====================================================================== */
		#ifdef LEFT_WALL_VELOCITY_FREESLIP
			// 1. Hydrodynamic Forces & Torques (from Lagrangian_collect_forces)
			F[0] = 0.0;         // Normal force cancels out exactly
			F[1] *= 2.0;        // Tangential drag is doubled
			F[2] *= 2.0;        // Tangential lateral force is doubled

			T[0] *= 2.0;        // Normal torque (rotation in Y-Z plane) is doubled
			T[1] = 0.0;         // Tangential torque cancels out exactly
			T[2] = 0.0;         // Tangential torque cancels out exactly

			// 2. Rigid Body Velocity Integrals (from fluid velocity inside particle)
			Int_U[0] = 0.0;     
			Int_U[1] *= 2.0;    
			Int_U[2] *= 2.0;    

			Int_Omega[0] *= 2.0; 
			Int_Omega[1] = 0.0;  
			Int_Omega[2] = 0.0;  

			#ifdef VOF_IBM
			// 3. Capillary Forces & Torques (from VOF_accumulate_solid_capillary_force)
			p->F_CCF[0] = 0.0;  
			p->F_CCF[1] *= 2.0; 
			p->F_CCF[2] *= 2.0; 

			p->T_CCF[0] *= 2.0; 
			p->T_CCF[1] = 0.0;  
			p->T_CCF[2] = 0.0;  
			#endif
		#endif
		/* ====================================================================== */

		#ifdef FORCES_DAT_OLD
			// Output fluid forces
			if (params->Np_mobile == 1) {
				double F_ddt = idt2beta * (Int_U[1] - Int_U_old[1]);
				double F_grav = p->M * (1.0 - 1.0 / p->rho_s) * params->grav[1];
				double F_IBM = F[1] - F_ddt;

				FILE *fptr;
				if (params -> time == 0) {
					fptr = fopen("forces.dat", "w");
				}
				else {
					fptr = fopen("forces.dat", "a");
				}

				fprintf(fptr, "%.10g, %.10g, %.10g, %.10g, ", params->time, F_IBM, F_ddt, F_grav);
				fclose(fptr);
			}
		#endif

		// Reset forces measured over entire timestep
		if (params->which_stage == 0) {
			DSET_ZERO(F_IBM, 3);
			DSET_ZERO(T_IBM, 3);
			DSET_ZERO(F_rigid, 3);
			DSET_ZERO(T_rigid, 3);
			DSET_ZERO(p->F_coll, 3);
			DSET_ZERO(p->T_coll, 3);
			#ifdef VOF_IBM
				DSET_ZERO(p->F_CCF_cum, 3); // Reset CCF accumulator
				DSET_ZERO(p->T_CCF_cum, 3);
			#endif
			#ifdef POST_PROCESS
				DSET_ZERO(p->Fc_norm_cum, 3);
				DSET_ZERO(p->Fc_tan_cum, 3);
				DSET_ZERO(p->Fl_norm_cum, 3);
				DSET_ZERO(p->Fl_tan_cum, 3);
			#endif
		}

		// IBM force acting on particle over entire timestep

		#ifndef ONE_WAY
			FORI3 F_IBM[i] += 2.0 * bet * F[i];
			FORI3 T_IBM[i] += 2.0 * bet * T[i];

			#ifdef VOF_IBM
				// CCF force acting on particle over entire timestep
				FORI3 p->F_CCF_cum[i] += 2.0 * bet * p->F_CCF[i];
				FORI3 p->T_CCF_cum[i] += 2.0 * bet * p->T_CCF[i];
			    
				// Add instantaneous CCF force to total fluid force
				F[0] += p->F_CCF[0];
				F[1] += p->F_CCF[1];
				F[2] += p->F_CCF[2];

				T[0] += p->T_CCF[0];
				T[1] += p->T_CCF[1];
				T[2] += p->T_CCF[2];

				//----------------------------------------------------------------------
				// Subtract spurious CSF force inside solid
				//----------------------------------------------------------------------
				
				//  F[0] -= p->F_CSF_solid[0];
				//  F[1] -= p->F_CSF_solid[1];
				//  F[2] -= p->F_CSF_solid[2];

				//  T[0] -= p->T_CSF_solid[0];
				//  T[1] -= p->T_CSF_solid[1];
				//  T[2] -= p->T_CSF_solid[2];

			#endif


			//----------------------------------------------------------------------
			// Rigid body force
			//----------------------------------------------------------------------
			F[0] += idt2beta * (Int_U[0] - Int_U_old[0]);
			F[1] += idt2beta * (Int_U[1] - Int_U_old[1]);
			F[2] += idt2beta * (Int_U[2] - Int_U_old[2]);

			T[0] += idt2beta * (Int_Omega[0] - Int_Omega_old[0]);
			T[1] += idt2beta * (Int_Omega[1] - Int_Omega_old[1]);
			T[2] += idt2beta * (Int_Omega[2] - Int_Omega_old[2]);

			// Rigid body force acting on particle over entire timestep
			FORI3 F_rigid[i] += 1.0 / dt * (Int_U[i] - Int_U_old[i]);
			FORI3 T_rigid[i] += 1.0 / dt * (Int_Omega[i] - Int_Omega_old[i]);
		#endif

		Int_U_old[0] = Int_U[0];
		Int_U_old[1] = Int_U[1];
		Int_U_old[2] = Int_U[2];

		Int_Omega_old[0] = Int_Omega[0];
		Int_Omega_old[1] = Int_Omega[1];
		Int_Omega_old[2] = Int_Omega[2];


		#ifdef DRY_COLLISION
			// Reset collision forces if maximum Stokes number > 5
			if (Collision_above_critical(p)) {
				DSET_ZERO(F, 3);
				DSET_ZERO(T, 3);
			}
		#endif

		//----------------------------------------------------------------------
		// Gravitational force
		//----------------------------------------------------------------------
		#ifdef DRY_PARTICLES
			// Normal gravity for dry particles
			F[0] += p->M * params->grav[0];
			F[1] += p->M * params->grav[1];
			F[2] += p->M * params->grav[2];
			// F[0] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[0];
			// F[1] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[1];
			// F[2] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[2];

			#ifdef STOKES_DRAG
				//double top_wall_vel, top_wall_y;
				// Stokes drag for a sheared suspension
				/*#if defined TOP_WALL_VELOCITY && defined BOTTOM_WALL_VELOCITY
				top_wall_vel = params->ubulk_target;
				#else
				top_wall_vel = 2.0 * params->ubulk_target;
				#endif
				#ifdef DOWNWARD_MOVING_WALL
				top_wall_y = (params->vel_init_y0 - params->ymax) * params->time / params->time_max + params->ymax;
				#else
				top_wall_y = params->ymax;
				#endif
				double shear_rate = 2.0 * params->ubulk_target / (top_wall_y - params->ymin);*/

				double U_inf[3]; // Modeled fluid velocity
				//double L=5;//dimension of each Langmuir Cell
				double U_0 = 1.0; // Mean background flow velocity
				//double W=0.0;
				//W=  (p->M * (1.0 - 1.0 / p->rho_s) * params->grav[1] * params->Re) / (6* PI * p->R);
				//double k_wavenumber= params -> k_wavenumber;
				//U_inf[0] =  top_wall_vel - shear_rate * (top_wall_y - p->X[1]);

				U_inf[0] = (U_0/PI)*sin(p->X[0]*PI)*cos(p->X[1]*PI);
				U_inf[1] =  (-1)*(U_0/PI)*cos(p->X[0]*PI)*sin(p->X[1]*PI);
					U_inf[2] =  0.0;
					F[0] += -6.0 * PI * p->R * (p->U[0]-U_inf[0])/ params->Re;
					F[1] += -6.0 * PI * p->R * (p->U[1]-U_inf[1])/ params->Re;
					F[2] += -6.0 * PI * p->R * (p->U[2]-U_inf[2])/ params->Re;

					T[0] +=  0.0;
					T[1] +=  0.0;
				T[2] +=  0.0;
				//T[2] += -8.0 * PI * p->R * p->R * (p->Omega[2]+shear_rate)/ params->Re;
			#endif  // Stokes drag
		#else
			#ifdef ONE_WAY
				// Require uniform grid
				// where is particle
				int ip = (int) round(p->X[0]/(params->xmax-params->xmin)*params->NXM);
				int jp = (int) round(p->X[1]/(params->ymax-params->ymin)*params->NYM);
				int kp = (int) round(p->X[2]/(params->zmax-params->zmin)*params->NZM);

				F[0] += -6.0 * PI * p->R * (p->U[0]-u[kp][jp][ip])/ params->Re;
				F[1] += -6.0 * PI * p->R * (p->U[1]-v[kp][jp][ip])/ params->Re;
				F[2] += -6.0 * PI * p->R * (p->U[2]-w[kp][jp][ip])/ params->Re;

			#endif

		#ifdef VOF_IBM         /* gravity logic ---------------------- */
		{

				p->F[0] += ( p->M  - p->Int_rho_scalar ) * params->grav[0];
				p->F[1] += ( p->M  - p->Int_rho_scalar ) * params->grav[1];
				p->F[2] += ( p->M  - p->Int_rho_scalar ) * params->grav[2];

		}
			#else                  /* original single-phase expression --------------- */
				F[0] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[0];
				F[1] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[1];
				F[2] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[2];
			#endif


			#ifdef OSCILLATION
				if (params->oscillation_frame == 1){ 		// non-inertial (accelerated) frame
					double oscillation = params -> oscillation;
					F[0] -= p->M * (1.0 - 1.0 / p->rho_s) * oscillation;		// oscillation has to be subtracted due to the non-inertial frame
				} 	
				// else if inertial (fixed) frame, there is no need to modify F[0]!
			#endif // OSCILLATION


			// TO REMOVE, JUST TEST
			#ifdef ONE_WAY
				FORI3 p->F_rigid[0] += -6.0 * PI * p->R * (p->U[0]-u[kp][jp][ip])/ params->Re;
				FORI3 p->F_rigid[1] += -6.0 * PI * p->R * (p->U[1]-v[kp][jp][ip])/ params->Re;
				FORI3 p->F_rigid[2] += -6.0 * PI * p->R * (p->U[2]-w[kp][jp][ip])/ params->Re;
			#endif

		#endif

		p = p -> next;
	}

	//--------------------------------------------------------------------------
	// Evaluate fluid forces for fixed particles
	//--------------------------------------------------------------------------
	p = p_fixed_list -> start;
	while (p != NULL) {
		F = p -> F;
		T = p -> T;
		Fc = p -> Fc; 			// To calculate the contact force for the fixed particles
		Int_U     = p -> Int_U;
		Int_U_old = p -> Int_U_old;
		Int_Omega     = p -> Int_Omega;
		Int_Omega_old = p -> Int_Omega_old;

		ID = p -> ID;		// fixed particle ID used to move the upper and lower wall (fixed particles)

		F_IBM = p -> F_IBM;
		T_IBM = p -> T_IBM;
		F_rigid = p -> F_rigid;
		T_rigid = p -> T_rigid;

		// Reset forces measured over entire timestep
		if (params->which_stage == 0) {
			DSET_ZERO(F_IBM, 3);
			DSET_ZERO(T_IBM, 3);
			DSET_ZERO(F_rigid, 3);
			DSET_ZERO(T_rigid, 3);
			DSET_ZERO(p->F_coll, 3);
			DSET_ZERO(p->T_coll, 3);
			#ifdef VOF_IBM
				DSET_ZERO(p->F_CCF_cum, 3); // Reset CCF accumulator
				DSET_ZERO(p->T_CCF_cum, 3);
			#endif
			#ifdef POST_PROCESS
				DSET_ZERO(p->Fc_norm_cum, 3);
				DSET_ZERO(p->Fc_tan_cum, 3);
				DSET_ZERO(p->Fl_norm_cum, 3);
				DSET_ZERO(p->Fl_tan_cum, 3);
			#endif
		}

		// IBM force acting on particle over entire timestep
		#ifndef ONE_WAY
			FORI3 F_IBM[i] += 2.0 * bet * F[i];
			FORI3 T_IBM[i] += 2.0 * bet * T[i];

			#ifdef VOF_IBM
				// CCF force acting on particle over entire timestep
				FORI3 p->F_CCF_cum[i] += 2.0 * bet * p->F_CCF[i];
				FORI3 p->T_CCF_cum[i] += 2.0 * bet * p->T_CCF[i];
			    
				// Add instantaneous CCF force to total fluid force
				F[0] += p->F_CCF[0];
				F[1] += p->F_CCF[1];
				F[2] += p->F_CCF[2];

				T[0] += p->T_CCF[0];
				T[1] += p->T_CCF[1];
				T[2] += p->T_CCF[2];
			#endif

			// Rigid body force acting on particle over entire timestep
			FORI3 F_rigid[i] += 1.0 / dt * (Int_U[i] - Int_U_old[i]);
			FORI3 T_rigid[i] += 1.0 / dt * (Int_Omega[i] - Int_Omega_old[i]);
		#endif
		Int_U_old[0] = Int_U[0];
		Int_U_old[1] = Int_U[1];
		Int_U_old[2] = Int_U[2];

		Int_Omega_old[0] = Int_Omega[0];
		Int_Omega_old[1] = Int_Omega[1];
		Int_Omega_old[2] = Int_Omega[2];

		#ifdef DRY_COLLISION
		// Reset collision forces if maximum Stokes number > 5
			if (Collision_above_critical(p)) {
				DSET_ZERO(F, 3);
				DSET_ZERO(T, 3);
			}
		#endif
			
		//----------------------------------------------------------------------
		// Gravitational force
		//----------------------------------------------------------------------
		#ifdef DRY_PARTICLES
			// Normal gravity for dry particles
			F[0] += p->M * params->grav[0];
			F[1] += p->M * params->grav[1];
			F[2] += p->M * params->grav[2];

			#ifdef STOKES_DRAG

				double U_inf[3]; // Modeled fluid velocity
				//double L=5;//dimension of each Langmuir Cell
				double U_0 = 1.0; // Mean background flow velocity

				U_inf[0] = (U_0/PI)*sin(p->X[0]*PI)*cos(p->X[1]*PI);
				U_inf[1] =  (-1)*(U_0/PI)*cos(p->X[0]*PI)*sin(p->X[1]*PI);
					U_inf[2] =  0.0;
					F[0] += -6.0 * PI * p->R * (p->U[0]-U_inf[0])/ params->Re;
					F[1] += -6.0 * PI * p->R * (p->U[1]-U_inf[1])/ params->Re;
					F[2] += -6.0 * PI * p->R * (p->U[2]-U_inf[2])/ params->Re;

					T[0] +=  0.0;
					T[1] +=  0.0;
				T[2] +=  0.0;
				//T[2] += -8.0 * PI * p->R * p->R * (p->Omega[2]+shear_rate)/ params->Re;
			#endif  // Stokes drag

		#else

    		 #ifdef VOF_IBM         /* gravity logic ---------------------- */
				{
				p->F[0] += ( p->M  - p->Int_rho_scalar ) * params->grav[0];
				p->F[1] += ( p->M  - p->Int_rho_scalar ) * params->grav[1];
				p->F[2] += ( p->M  - p->Int_rho_scalar ) * params->grav[2];
				}
			#else                  /* original single-phase expression --------------- */
				F[0] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[0];
				F[1] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[1];
				F[2] += p->M * (1.0 - 1.0 / p->rho_s) * params->grav[2];
			#endif


		#endif

		p = p -> next;
	}
}




/******************************************************************************/
/*
 * Predicts location of particles for evaluating fluid-particle forces
 * Should have only local particles when called
 * Has only local particles afterwards
 */
/******************************************************************************/
void Lagrangian_integrate_particle_motion(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i;
	int ID; 		// fixed particle ID used to move the upper and lower wall (fixed particles)
	Particle *p;
	Collision *pc;
	Particle_list *p_mobile_list_foreign, *p_fixed_list_foreign;
	double dt_iM, dt_iI, dt_iM1;
	double *X, *X_old, *U, *U_old, *Omega, *Omega_old;
	double *F, *T, *Fc, *Fc_old, *Tc, *Tc_old;
	double *F_coll, *T_coll;
	double F_local_sum=0, F_global_sum=0;	// To calculate total force on each processor locally and communicate to find the global sum of all forces
	double  Mass=0, total_Mass=0;			// To calculate loacl mass on each processor and calculate the global mass from all processors

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	const double BET[] = {BETA};
	const double GAM[] = {GAMMA};
	const double ZET[] = {ZETA};

	double bet = BET[params -> which_stage];
	double gam = GAM[params -> which_stage];
	double zet = ZET[params -> which_stage];

	double dt  = params -> dt;
	double time = params -> time;
	double startup_time = params -> startup_time;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	//--------------------------------------------------------------------------
	// Predict particle locations and velocities
	//--------------------------------------------------------------------------
	p = p_mobile_list -> start;
	while (p != NULL) {

		X = p -> X;
		U = p -> U;
		Omega = p -> Omega;
		X_old = p -> X_old;
		U_old = p -> U_old;
		Omega_old = p -> Omega_old;

		F = p -> F;
		T = p -> T;
		Fc = p -> Fc;
		Tc = p -> Tc;
		Fc_old = p -> Fc_old;
		Tc_old = p -> Tc_old;

		dt_iM = dt / p -> M;
		dt_iI = dt / p -> I_p;

		X_old[0] = X[0];
		X_old[1] = X[1];
		X_old[2] = X[2];

		FORI3 p -> Rotn_old[0][i] = p -> Rotn[0][i];
		FORI3 p -> Rotn_old[1][i] = p -> Rotn[1][i];
		FORI3 p -> Rotn_old[2][i] = p -> Rotn[2][i];

		pc = p -> particle_collision;
		while (pc != NULL) {
			FORI3 pc->zeta_t_old[i] = pc->zeta_t[i];
			pc = pc -> next;
		}
		pc = p -> wall_collision;
		while (pc != NULL) {
			FORI3 pc->zeta_t_old[i] = pc->zeta_t[i];
			pc = pc -> next;
		}

		U_old[0] = U[0];
		U_old[1] = U[1];
		U_old[2] = U[2];

		#ifdef OSCILLATING_PARTICLE
			double amplitude = params -> disp_amplitude;
			double frequency = params -> frequency;
			double angular_vel = 2 * PI * frequency;	// = omega
			double ps = params -> phase_shift_factor;
			
			U[0] = amplitude * angular_vel * (-cos(time * angular_vel + ps * PI)); 
			U[1] = 0;
			U[2] = 0;
		#else
			U[0] = U_old[0] + dt_iM * ( 2.0 * bet * F[0] + gam * Fc[0] + zet * Fc_old[0]);
			U[1] = U_old[1] + dt_iM * ( 2.0 * bet * F[1] + gam * Fc[1] + zet * Fc_old[1]);
			U[2] = U_old[2] + dt_iM * ( 2.0 * bet * F[2] + gam * Fc[2] + zet * Fc_old[2]);
		#endif

		Fc_old[0] = Fc[0];
		Fc_old[1] = Fc[1];
		Fc_old[2] = Fc[2];

#ifdef POST_PROCESS
		FORI3 p->Fc_norm_old[i] = p->Fc_norm[i];
		FORI3 p->Fc_tan_old[i]  = p->Fc_tan[i];
		FORI3 p->Fl_norm_old[i] = p->Fl_norm[i];
		FORI3 p->Fl_tan_old[i]  = p->Fl_tan[i];
#endif

		Omega_old[0] = Omega[0];
		Omega_old[1] = Omega[1];
		Omega_old[2] = Omega[2];
		Omega[0] = Omega_old[0] + dt_iI * ( 2.0 * bet * T[0] + gam * Tc[0] + zet * Tc_old[0] );
		Omega[1] = Omega_old[1] + dt_iI * ( 2.0 * bet * T[1] + gam * Tc[1] + zet * Tc_old[1] );
		Omega[2] = Omega_old[2] + dt_iI * ( 2.0 * bet * T[2] + gam * Tc[2] + zet * Tc_old[2] );
		Tc_old[0] = Tc[0];
		Tc_old[1] = Tc[1];
		Tc_old[2] = Tc[2];

#ifdef STARTUP
		if (params->startup_init == STUP_INIT_TIME && time >= startup_time) {
			params->startup_flag = 0;
		}
		if ((params->startup_init == STUP_INIT_GOND_ST_27 || params->startup_init == STUP_INIT_GOND_ST_152)
																	 && X[1] <= 2*p->R) {
			params->startup_flag = 0;
		}
		if (params->startup_flag) {
			if(params->startup_init == STUP_INIT_TIME) {
				DSET_ZERO(Omega, 3);
				FORI3 U[i] = params->startup_velocity[i];
			}
			if(params->startup_init == STUP_INIT_GOND_ST_27) {
				DSET_ZERO(Omega, 3);
				U[1] = 0.518 * ( exp(-40 * params->time) - 1 );
				}  // Gondret St=27
			if(params->startup_init == STUP_INIT_GOND_ST_152) {
				DSET_ZERO(Omega, 3);
				U[1] = 0.585 * ( exp(-40 * params->time) - 1 ); 
				} // Gondret St=152
//			U[1]=1;
//			Omega[2]=0;
//			U[1] = 0.585 * ( exp(-40 * params->time) - 1 );  // Gondret2
//			U[1] = 10.518 * ( exp(-40 * params->time) - 1 );  // Gondret10d
//			U[1] = 0.385 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=20
//			U[1] = 0.288 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=15
//			U[1] = 0.192 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=10
//			U[1] = 0.096 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=5
//			U[1] = 0.020 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=1
		}
#endif

		X[0] = X_old[0] + dt * bet * ( U[0] + U_old[0] );
		X[1] = X_old[1] + dt * bet * ( U[1] + U_old[1] );
		X[2] = X_old[2] + dt * bet * ( U[2] + U_old[2] );

		Rotate_particle(p, params);

		/* ── Symmetry-plane constraint: pin sphere to x = 0 ── */
		#ifdef LEFT_WALL_VELOCITY_FREESLIP
		X[0]     = 0.0;
		U[0]     = 0.0;
		U_old[0] = 0.0;
		DSET_ZERO(Omega, 3);
		DSET_ZERO(Omega_old, 3);
		#endif

		p = p -> next;
	}	

	if(params->startup_init == STUP_INIT_PRIMPOSED)
		{
			p = p_fixed_list -> start;
			while (p != NULL) {
				F = p -> F;
				Fc = p -> Fc;
				T = p -> T;
				Fc_old = p -> Fc_old;

				#ifdef STARTUP
					if (params->startup_flag) {
						if(params->startup_init == STUP_INIT_PRIMPOSED && ( (p -> ID ) % 2 == 1 ) )
							{
								#ifdef DRY_PARTICLES
									//F[1] += (p->M) * params->grav[1];
									Mass = Mass + (p->M1) ;
									F[1] = Mass * (-9.81);
								#else
									// Reduced gravity for submerged particles
									F[1] += p->M1 * (1.0 - 1.0 / p->rho_prImp) * params->grav[1];
								#endif
								Mass = Mass + p->M1;
								F_local_sum += ( 2.0 * bet * F[1] + gam * Fc[1] + zet * Fc_old[1]) ;
							}
					}
				#endif
				p = p -> next;
			}

				MPI_Allreduce(&F_local_sum, &F_global_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);
				MPI_Allreduce(&Mass, &total_Mass, 1, MPI_DOUBLE, MPI_SUM, PCW);
		}	

	// Update old forces for fixed particles
	if(params->startup_init == STUP_INIT_PRIMPOSED || params->startup_init == STUP_INIT_SHEARFLOW)
		{
			p = p_fixed_list -> start;
			while (p != NULL) 
				{

					X = p -> X;
					U = p -> U;
					Omega = p -> Omega;
					X_old = p -> X_old;
					U_old = p -> U_old;
					Omega_old = p -> Omega_old;
					ID = p -> ID;
					F = p -> F;
					T = p -> T;
					Fc = p -> Fc;
					Tc = p -> Tc;

					Fc_old = p -> Fc_old;
					Tc_old = p -> Tc_old;
					dt_iM = dt / p -> M;
					dt_iM1 = dt / p -> M1;
					dt_iI = dt / p -> I_p;

					X_old[0] = X[0];
					X_old[1] = X[1];
					X_old[2] = X[2];

					FORI3 p -> Rotn_old[0][i] = p -> Rotn[0][i];
					FORI3 p -> Rotn_old[1][i] = p -> Rotn[1][i];
					FORI3 p -> Rotn_old[2][i] = p -> Rotn[2][i];

					pc = p -> particle_collision;
					while (pc != NULL) 
					{
						FORI3 pc->zeta_t_old[i] = pc->zeta_t[i];
						pc = pc -> next;
					}
					pc = p -> wall_collision;
					while (pc != NULL) 
					{
						FORI3 pc->zeta_t_old[i] = pc->zeta_t[i];
						pc = pc -> next;
					}

					U_old[0] = U[0];
					U_old[1] = U[1];
					U_old[2] = U[2];

					Fc_old[0] = Fc[0];
					Fc_old[1] = Fc[1];
					Fc_old[2] = Fc[2];

					#ifdef POST_PROCESS
						FORI3 p->Fc_norm_old[i] = p->Fc_norm[i];
						FORI3 p->Fc_tan_old[i]  = p->Fc_tan[i];
						FORI3 p->Fl_norm_old[i] = p->Fl_norm[i];
						FORI3 p->Fl_tan_old[i]  = p->Fl_tan[i];
					#endif

					Omega_old[0] = Omega[0];
					Omega_old[1] = Omega[1];
					Omega_old[2] = Omega[2];

					Tc_old[0] = Tc[0];
					Tc_old[1] = Tc[1];
					Tc_old[2] = Tc[2];

					#ifdef STARTUP
						if (params->startup_flag) 
						{
							DSET_ZERO(Omega, 3);
							if(params->startup_init == STUP_INIT_SHEARFLOW && ( (p -> ID ) % 2 == 0 ) )
								{
									U[0] = - params->ubulk_target ;
								}
							if(params->startup_init == STUP_INIT_SHEARFLOW && ( (p -> ID ) % 2 == 1 ) )
								{
									U[0] =  params->ubulk_target ;
								}
						}
					#endif

					#ifdef STARTUP
						if (params->startup_flag) 
						{
							DSET_ZERO(Omega, 3);
							if(params->startup_init == STUP_INIT_PRIMPOSED && ( (p -> ID ) % 2 == 0 ) )
								{
									U[0] = - params->ubulk_target ;
								}
							if(params->startup_init == STUP_INIT_PRIMPOSED && ( (p -> ID ) % 2 == 1 ) )
								{
									U[0] = params->ubulk_target ;
									U[1] = U_old[1] + (dt / total_Mass )* F_global_sum ;
									U[2] = 0;
								}
						}
					#endif

			    X[0] = X_old[0] + dt * bet * ( U[0] + U_old[0] );
				X[1] = X_old[1] + dt * bet * ( U[1] + U_old[1] );
				X[2] = X_old[2] + dt * bet * ( U[2] + U_old[2] );

				p = p -> next;

		       	}
		}					

	else{
		// Update old forces for fixed particles
		p = p_fixed_list -> start;
		while (p != NULL) {
			FORI3 p->Fc_old[i] = p->Fc[i];
			FORI3 p->Tc_old[i] = p->Tc[i];
			#ifdef POST_PROCESS
					FORI3 p->Fc_norm_old[i] = p->Fc_norm[i];
					FORI3 p->Fc_tan_old[i]  = p->Fc_tan[i];
					FORI3 p->Fl_norm_old[i] = p->Fl_norm[i];
					FORI3 p->Fl_tan_old[i]  = p->Fl_tan[i];
			#endif
			p = p -> next;
		}
	}
#ifdef STARTUP
	MPI_Allreduce(MPI_IN_PLACE, &(params->startup_flag), 1, MPI_INT, MPI_MIN, PCW);
#endif

	//--------------------------------------------------------------------------
	// Evaluate collision forces based on predicted position and velocity for
	// p_mobile.  Transfer new collision parameters to fixed particles.
	//--------------------------------------------------------------------------
	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));

	p_mobile_list_foreign = Collision_evaluate(data_bag, DTRACE("Collision_evaluate"));
	p_fixed_list_foreign = Particle_list_foreign_create(p_fixed_list, data_bag, DTRACE("Particle_list_foreign_create"));
	Particle_MPI_update(p_mobile_list_foreign, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list_foreign, data_bag, DTRACE("Particle_MPI_update"));

	Lagrangian_collect_forces(p_mobile_list, p_mobile_list_foreign, LAG_COLLECT_COLL, params, DTRACE("Lagrangian_collect_forces"));
	Lagrangian_collect_forces(p_fixed_list, p_fixed_list_foreign, LAG_COLLECT_COLL, params, DTRACE("Lagrangian_collect_forces"));

	Particle_list_destroy(p_mobile_list_foreign);
	Particle_list_destroy(p_fixed_list_foreign);


	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	//--------------------------------------------------------------------------
	// Correct particle positions and velocities
	//--------------------------------------------------------------------------
	p = p_mobile_list -> start;
	while (p != NULL) {

		X = p -> X;
		U = p -> U;
		Omega = p -> Omega;
		X_old = p -> X_old;
		U_old = p -> U_old;
		Omega_old = p -> Omega_old;

		F = p -> F;
		T = p -> T;
		Fc = p -> Fc;
		Tc = p -> Tc;
		Fc_old = p -> Fc_old;
		Tc_old = p -> Tc_old;

		F_coll = p -> F_coll;
		T_coll = p -> T_coll;

		dt_iM = dt / p -> M;
		dt_iI = dt / p -> I_p;

		#ifdef OSCILLATING_PARTICLE
			double amplitude = params -> disp_amplitude;
			double frequency = params -> frequency;
			double angular_vel = 2 * PI * frequency;	// = omega
			double ps = params -> phase_shift_factor;
			
			U[0] = amplitude * angular_vel * (-cos(time * angular_vel + ps * PI)); 
			U[1] = 0;
			U[2] = 0;
		#else
			U[0] = U_old[0] + dt_iM * bet * ( 2.0 * F[0] + Fc[0] + Fc_old[0]);
			U[1] = U_old[1] + dt_iM * bet * ( 2.0 * F[1] + Fc[1] + Fc_old[1]);
			U[2] = U_old[2] + dt_iM * bet * ( 2.0 * F[2] + Fc[2] + Fc_old[2]);
		#endif

		Omega[0] = Omega_old[0] + dt_iI * bet * ( 2.0 * T[0] + Tc[0] + Tc_old[0] );
		Omega[1] = Omega_old[1] + dt_iI * bet * ( 2.0 * T[1] + Tc[1] + Tc_old[1] );
		Omega[2] = Omega_old[2] + dt_iI * bet * ( 2.0 * T[2] + Tc[2] + Tc_old[2] );

		// Store collision force acting over entire timestep
#ifdef SUBSTEP
		FORI3 F_coll[i] += bet * (Fc[i] + Fc_old[i]) / 15.0;
		FORI3 T_coll[i] += bet * (Tc[i] + Tc_old[i]) / 15.0;
	#ifdef POST_PROCESS
		FORI3 p->Fc_norm_cum[i] += bet * (p->Fc_norm[i] + p->Fc_norm_old[i]) / 15.0;
		FORI3 p->Fc_tan_cum[i]  += bet * (p->Fc_tan[i]  + p->Fc_tan_old[i]) / 15.0;
		FORI3 p->Fl_norm_cum[i] += bet * (p->Fl_norm[i] + p->Fl_norm_old[i]) / 15.0;
		FORI3 p->Fl_tan_cum[i]  += bet * (p->Fl_tan[i]  + p->Fl_tan_old[i]) / 15.0;
	#endif
#else
		FORI3 F_coll[i] += bet * (Fc[i] + Fc_old[i]);
		FORI3 T_coll[i] += bet * (Tc[i] + Tc_old[i]);
	#ifdef POST_PROCESS
		FORI3 p->Fc_norm_cum[i] += bet * (p->Fc_norm[i] + p->Fc_norm_old[i]);
		FORI3 p->Fc_tan_cum[i]  += bet * (p->Fc_tan[i]  + p->Fc_tan_old[i]);
		FORI3 p->Fl_norm_cum[i] += bet * (p->Fl_norm[i] + p->Fl_norm_old[i]);
		FORI3 p->Fl_tan_cum[i]  += bet * (p->Fl_tan[i]  + p->Fl_tan_old[i]);
	#endif
#endif

#ifdef STARTUP
		if (params->startup_init == STUP_INIT_TIME && time >= startup_time) {
			params->startup_flag = 0;
		}
		if ((params->startup_init == STUP_INIT_GOND_ST_27 || params->startup_init == STUP_INIT_GOND_ST_152)
																	 && X[1] <= 2*p->R) {
			params->startup_flag = 0;
		}
		if (params->startup_flag) {
			if(params->startup_init == STUP_INIT_TIME) {
				DSET_ZERO(Omega, 3);
				FORI3 U[i] = params->startup_velocity[i];
			}
			if(params->startup_init == STUP_INIT_GOND_ST_27) {
				DSET_ZERO(Omega, 3);
				U[1] = 0.518 * ( exp(-40 * params->time) - 1 );
				}  // Gondret St=27
			if(params->startup_init == STUP_INIT_GOND_ST_152) {
				DSET_ZERO(Omega, 3);
				U[1] = 0.585 * ( exp(-40 * params->time) - 1 ); 
				} // Gondret St=152
//			U[1]=1;
//			Omega[2]=0;
//			U[1] = 0.585 * ( exp(-40 * params->time) - 1 );  // Gondret2
//			U[1] = 10.518 * ( exp(-40 * params->time) - 1 );  // Gondret10d
//			U[1] = 0.385 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=20
//			U[1] = 0.288 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=15
//			U[1] = 0.192 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=10
//			U[1] = 0.096 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=5
//			U[1] = 0.020 * ( exp(-40 * params->time) - 1 );  // Gondret10d: St=1
		}
#endif

		X[0] = X_old[0] + dt * bet * ( U[0] + U_old[0] );
		X[1] = X_old[1] + dt * bet * ( U[1] + U_old[1] );
		X[2] = X_old[2] + dt * bet * ( U[2] + U_old[2] );

		Rotate_particle(p, params);

		/* ── Symmetry-plane constraint: pin sphere to x = 0 ── */
		#ifdef LEFT_WALL_VELOCITY_FREESLIP
		X[0]     = 0.0;
		U[0]     = 0.0;
		U_old[0] = 0.0;
		DSET_ZERO(Omega, 3);
		DSET_ZERO(Omega_old, 3);
		#endif


		//printf("Particle position x,y,z: %f %f %f \n", X[0], X[1], X[2] );


		p = p -> next;
	}
#ifdef STARTUP
	MPI_Allreduce(MPI_IN_PLACE, &(params->startup_flag), 1, MPI_INT, MPI_MIN, PCW);
#endif

	if(params->startup_init == STUP_INIT_PRIMPOSED || params->startup_init == STUP_INIT_SHEARFLOW)
		{
			MPI_Allreduce(&F_local_sum, &F_global_sum, 1, MPI_DOUBLE, MPI_SUM, PCW);
			//return F_global_sum;
			MPI_Allreduce(&Mass, &total_Mass, 1, MPI_DOUBLE, MPI_SUM, PCW);		

			// Store collision force acting over entire timestep for fixed particles
			p = p_fixed_list -> start;
			while (p != NULL) 
			{
				X = p -> X;
				U = p -> U;
				Omega = p -> Omega;
				X_old = p -> X_old;
				U_old = p -> U_old;
				Omega_old = p -> Omega_old;

				F = p -> F;
				T = p -> T;
				Fc = p -> Fc;
				Tc = p -> Tc;
				Fc_old = p -> Fc_old;
				Tc_old = p -> Tc_old;

				F_coll = p -> F_coll;
				T_coll = p -> T_coll;

				dt_iM = dt / p -> M;
				dt_iM1 = dt / p -> M1;
				dt_iI = dt / p -> I_p;


			// Store collision force acting over entire timestep
			#ifdef SUBSTEP
					FORI3 F_coll[i] += bet * (Fc[i] + Fc_old[i]) / 15.0;
					FORI3 T_coll[i] += bet * (Tc[i] + Tc_old[i]) / 15.0;
				#ifdef POST_PROCESS
					FORI3 p->Fc_norm_cum[i] += bet * (p->Fc_norm[i] + p->Fc_norm_old[i]) / 15.0;
					FORI3 p->Fc_tan_cum[i]  += bet * (p->Fc_tan[i]  + p->Fc_tan_old[i]) / 15.0;
					FORI3 p->Fl_norm_cum[i] += bet * (p->Fl_norm[i] + p->Fl_norm_old[i]) / 15.0;
					FORI3 p->Fl_tan_cum[i]  += bet * (p->Fl_tan[i]  + p->Fl_tan_old[i]) / 15.0;
				#endif
			#else
					FORI3 F_coll[i] += bet * (Fc[i] + Fc_old[i]);
					FORI3 T_coll[i] += bet * (Tc[i] + Tc_old[i]);
				#ifdef POST_PROCESS
					FORI3 p->Fc_norm_cum[i] += bet * (p->Fc_norm[i] + p->Fc_norm_old[i]);
					FORI3 p->Fc_tan_cum[i]  += bet * (p->Fc_tan[i]  + p->Fc_tan_old[i]);
					FORI3 p->Fl_norm_cum[i] += bet * (p->Fl_norm[i] + p->Fl_norm_old[i]);
					FORI3 p->Fl_tan_cum[i]  += bet * (p->Fl_tan[i]  + p->Fl_tan_old[i]);
				#endif
			#endif

			#ifdef STARTUP
			if (params->startup_flag) 
			{
				DSET_ZERO(Omega, 3);
				if(params->startup_init == STUP_INIT_SHEARFLOW && ( (p -> ID) % 2 == 0 ) )
				{
					U[0] = - params->ubulk_target ;
				}
				if(params->startup_init == STUP_INIT_SHEARFLOW && ( (p -> ID) % 2 == 1 ) )
				{
					U[0] = params->ubulk_target ;
				}
			}
			#endif

			#ifdef STARTUP
			if (params->startup_flag) 
			{
				DSET_ZERO(Omega, 3);
				if(params->startup_init == STUP_INIT_PRIMPOSED && ( (p -> ID ) % 2 == 0 ) )
				{
					U[0] = - params->ubulk_target ;
					U[1] = 0;
				}
				if(params->startup_init == STUP_INIT_PRIMPOSED && ( (p -> ID ) % 2 == 1 ) )
				{
					U[0] = params->ubulk_target ;
					U[1] = U_old[1] + (dt  / total_Mass )* F_global_sum ;
					U[2] = 0;
				}
			}
			#endif

			X[0] = X_old[0] + dt * bet * ( U[0] + U_old[0] );
			X[1] = X_old[1] + dt * bet * ( U[1] + U_old[1] );
			X[2] = X_old[2] + dt * bet * ( U[2] + U_old[2] );

			p = p -> next;
			}
		}
	else{
		// Store collision force acting over entire timestep for fixed particles
		p = p_fixed_list -> start;
		while (p != NULL) {
			Fc = p -> Fc;
			Tc = p -> Tc;
			Fc_old = p -> Fc_old;
			Tc_old = p -> Tc_old;
			F_coll = p -> F_coll;
			T_coll = p -> T_coll;
			#ifdef SUBSTEP
					FORI3 F_coll[i] += bet * (Fc[i] + Fc_old[i]) / 15.0;
					FORI3 T_coll[i] += bet * (Tc[i] + Tc_old[i]) / 15.0;
				#ifdef POST_PROCESS
					FORI3 p->Fc_norm_cum[i] += bet * (p->Fc_norm[i] + p->Fc_norm_old[i]) / 15.0;
					FORI3 p->Fc_tan_cum[i]  += bet * (p->Fc_tan[i]  + p->Fc_tan_old[i]) / 15.0;
					FORI3 p->Fl_norm_cum[i] += bet * (p->Fl_norm[i] + p->Fl_norm_old[i]) / 15.0;
					FORI3 p->Fl_tan_cum[i]  += bet * (p->Fl_tan[i]  + p->Fl_tan_old[i]) / 15.0;
				#endif
			#else
					FORI3 F_coll[i] += bet * (Fc[i] + Fc_old[i]);
					FORI3 T_coll[i] += bet * (Tc[i] + Tc_old[i]);
				#ifdef POST_PROCESS
					FORI3 p->Fc_norm_cum[i] += bet * (p->Fc_norm[i] + p->Fc_norm_old[i]);
					FORI3 p->Fc_tan_cum[i]  += bet * (p->Fc_tan[i]  + p->Fc_tan_old[i]);
					FORI3 p->Fl_norm_cum[i] += bet * (p->Fl_norm[i] + p->Fl_norm_old[i]);
					FORI3 p->Fl_tan_cum[i]  += bet * (p->Fl_tan[i]  + p->Fl_tan_old[i]);
				#endif
			#endif
			p = p -> next;
		}
	}	
	//--------------------------------------------------------------------------
	// Evaluate collision forces based on final position and velocity for
	// p_mobile.  Transfer new collision parameters to fixed particles.
	//--------------------------------------------------------------------------
	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));

	p_mobile_list_foreign = Collision_evaluate(data_bag, DTRACE("Collision_evaluate"));
	p_fixed_list_foreign = Particle_list_foreign_create(p_fixed_list, data_bag, DTRACE("Particle_list_foreign_create"));
	Particle_MPI_update(p_mobile_list_foreign, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list_foreign, data_bag, DTRACE("Particle_MPI_update"));

	Lagrangian_collect_forces(p_mobile_list, p_mobile_list_foreign, LAG_COLLECT_COLL, params, DTRACE("Lagrangian_collect_forces"));
	Lagrangian_collect_forces(p_fixed_list, p_fixed_list_foreign, LAG_COLLECT_COLL, params, DTRACE("Lagrangian_collect_forces"));

	Particle_list_destroy(p_mobile_list_foreign);
	Particle_list_destroy(p_fixed_list_foreign);





	// Now clean up any completed collisions
#ifdef STORE_COLLISION
	Particle_collision_list_clean(p_mobile_list->start);
	Particle_collision_list_clean(p_fixed_list->start);
#endif

#ifdef FORCES_DAT
	// Output fluid forces
	if (params->Np_mobile == 1) {

		p = p_mobile_list -> start;
		while (p != NULL) {
			FILE *fptr;
			if (params -> time == 0) {
				fptr = fopen("forces.dat", "w");
			}
			else {
				fptr = fopen("forces.dat", "a");
			}

			fprintf(fptr, "%.10g, %.10g, %.10g\n", params->time, p->F[1], p->Fc[1]);
			fclose(fptr);

			p = p -> next;
		}
	}
#endif

}




/******************************************************************************/
/*
 * Add forces on particle from other processors.  This function should only be
 * called immediately after calling Collision_evaulate().
 */
/******************************************************************************/
void Lagrangian_collect_forces(Particle_list *p_list, Particle_list *p_list_foreign,
		int type, Parameters *params, Debug_trace *dtrace) {

	int i;
	int status = 0;
	Particle *p, *pf;

	// Check state of linked lists
	Display_assert_list_state(p_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_list_foreign, LIST_STATE_EDGE, params, DTRACE("Display_assert_list_state"));

	p = p_list->start;
	while (p != NULL) {
		pf = p_list_foreign->start;
		while (pf != NULL) {
			if (p -> ID == pf -> ID) {

				// Hydrodynamic forces
				if (type == LAG_COLLECT_ALL || type == LAG_COLLECT_HYDRO) {
					FORI3 p -> F[i] += pf -> F[i];
					FORI3 p -> T[i] += pf -> T[i];
					
					#ifdef VOF_IBM
					FORI3 p -> F_CCF[i] += pf -> F_CCF[i]; 
					FORI3 p -> T_CCF[i] += pf -> T_CCF[i];
					FORI3 p -> Int_rho[i] += pf -> Int_rho[i];
					p->Int_rho_scalar += pf->Int_rho_scalar;  
					FORI3 p->F_CSF_solid[i] += pf->F_CSF_solid[i]; 
					FORI3 p->T_CSF_solid[i] += pf->T_CSF_solid[i];  
					#endif

					FORI3 p -> Int_U[i] += pf -> Int_U[i];
					FORI3 p -> Int_Omega[i] += pf -> Int_Omega[i];

				}

				// Collision forces
				if (type == LAG_COLLECT_ALL || type == LAG_COLLECT_COLL) {

					status = Particle_collision_list_xfer(pf, p, 'p');

					FORI3 p -> Fc[i] += pf -> Fc[i];
					FORI3 p -> Tc[i] += pf -> Tc[i];

#ifdef POST_PROCESS
					FORI3 p->Fc_norm[i] += pf->Fc_norm[i];
					FORI3 p->Fc_tan[i]  += pf->Fc_tan[i];
					FORI3 p->Fl_norm[i] += pf->Fl_norm[i];
					FORI3 p->Fl_tan[i]  += pf->Fl_tan[i];
#endif
				}
			}
			pf = pf -> next;
			if (status < 0) break;
		}
		p = p -> next;
		if (status < 0) break;
	}

	Display_assert_error(status,
		"Parties cannot handle more than one collision between two\n"
		"particles. Check if domain is one particle diameter wide.\n",
		params, DTRACE("Display_assert_error"));
}




/******************************************************************************/
/*
 Applies forcing onto the fluid flow field using the immersed boundary method to
 achieve the desired velocity, given by the particle's current position and
 velocity.
 This forcing can be carried out at two stages in the program, both of which
 occur before the pressure projection step:
     1. Forcing before implicit viscous terms
          - Forcing is applied to the velocity RHS term (vel -> ng_rhs)
          - Uses intermediate velocity field found using explicit viscous terms
          - Specify 'force_iter' to be a negative number
     2. Forcing after implicit viscous terms
          - Forcing is applied directly to the velocity field (vel -> data)
          - Uses intermediate velocity field found using implicit viscous terms
          - Specify 'force_iter' to be a positive number (including zero).  This
            is the number of times the forcing will be applied this step.
 */
/******************************************************************************/
void Lagrangian_force(int force_iter, Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	int i, j, k, ii, iters, corrector;
	Particle *p;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;

	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = grid -> G_Ie;
	int Je = grid -> G_Je;
	int Ke = grid -> G_Ke;

	const double BET[] = {BETA};
	double dtbeta = BET[params->which_stage] * params->dt;

	double ***u_data = u -> data;
	double ***v_data = v -> data;
	double ***w_data = w -> data;

	double ***u_rhs = u -> ng_rhs;
	double ***v_rhs = v -> ng_rhs;
	double ***w_rhs = w -> ng_rhs;

	double *u_rhs_1d = &(u -> ng_rhs[Ks][Js][Is]);
	double *v_rhs_1d = &(v -> ng_rhs[Ks][Js][Is]);
	double *w_rhs_1d = &(w -> ng_rhs[Ks][Js][Is]);

	double *u_implicit_1d = &(u -> ng_implicit[Ks][Js][Is]);
	double *v_implicit_1d = &(v -> ng_implicit[Ks][Js][Is]);
	double *w_implicit_1d = &(w -> ng_implicit[Ks][Js][Is]);

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	//--------------------------------------------------------------------------
	// 1. Forcing before implicit viscous terms
	//--------------------------------------------------------------------------
	if (force_iter < 0) {
		corrector = 0;
		iters = 1;
	}
	//--------------------------------------------------------------------------
	// 2. Forcing after implicit viscous terms
	//--------------------------------------------------------------------------
	else {
		corrector = 1;
		iters = force_iter;
	}

	for (ii = 0; ii < iters; ii++) {

		Communication_update_ghost_nodes_flow_variable(u->data, 'u', params->ghost_nodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(v->data, 'v', params->ghost_nodes, data_bag);
		Communication_update_ghost_nodes_flow_variable(w->data, 'w', params->ghost_nodes, data_bag);

		// Reset righthand side
		if (corrector) {
			DSET_ZERO(u_rhs_1d, grid->ng_total_nodes);
			DSET_ZERO(v_rhs_1d, grid->ng_total_nodes);
			DSET_ZERO(w_rhs_1d, grid->ng_total_nodes);
		}

		p = p_mobile_list -> start;
		while (p != NULL) {

			// Reset fluid forces
			if (!corrector) {
				DSET_ZERO(p->F, 3);
				DSET_ZERO(p->T, 3);
			}
#ifndef ONE_WAY
			Lagrangian_force_individual(MOBILE, p, corrector, data_bag);
#endif
			p = p -> next;
		}

		p = p_fixed_list -> start;

		while (p != NULL) {

			// Reset fluid forces
			if (!corrector) {
				DSET_ZERO(p->F, 3);
				DSET_ZERO(p->T, 3);
			}
#ifndef ONE_WAY
			Lagrangian_force_individual(FIXED, p, corrector, data_bag);
#endif
			p = p -> next;
		}

		#ifndef ONE_WAY
			// Update velocities
			if (corrector) {
				for (k = Ks; k < Ke; k++) {
					for (j = Js; j < Je; j++) {
						for (i = Is; i < Ie; i++) {
							#ifdef VOF_IBM
							u_data[k][j][i] += dtbeta * u_rhs[k][j][i] / ( 0.5 * (data_bag->vof->rho[k][j][i] + data_bag->vof->rho[k][j][i-1]) );
							v_data[k][j][i] += dtbeta * v_rhs[k][j][i] / ( 0.5 * (data_bag->vof->rho[k][j][i] + data_bag->vof->rho[k][j-1][i]) );
							w_data[k][j][i] += dtbeta * w_rhs[k][j][i] / ( 0.5 * (data_bag->vof->rho[k][j][i] + data_bag->vof->rho[k-1][j][i]) );

							#else
							u_data[k][j][i] += dtbeta * u_rhs[k][j][i];
							v_data[k][j][i] += dtbeta * v_rhs[k][j][i];
							w_data[k][j][i] += dtbeta * w_rhs[k][j][i];
							#endif //VOF_IBM


						}
					}
				}
			}
		#endif
	}

	#ifndef FULLY_EXPLICIT
		if (!corrector) {
			for (i = 0; i < grid->ng_total_nodes; i++) {
				u_rhs_1d[i] -= u_implicit_1d[i];
			}
			for (i = 0; i < grid->ng_total_nodes; i++) {
				v_rhs_1d[i] -= v_implicit_1d[i];
			}
			for (i = 0; i < grid->ng_total_nodes; i++) {
				w_rhs_1d[i] -= w_implicit_1d[i];
			}
		}
	#endif

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_particle_forc += T2 - T1;
}




/******************************************************************************/
/*
 Calculates forcing between an individual particle and the velocity field:
     1. Interpolate velocity field onto each Lagrangian marker
     2. Calculate acceleration required for marker to achieve desired velocity (of
        rigid body motion)
     3. Add the Lagrangian force to the force on the particle
     4. Spread the Lagrangian acceleration onto the Eulerian fluid velocity field
 */
/******************************************************************************/
void Lagrangian_force_individual(int p_type, Particle *p, int corrector, Cart3d_bag *data_bag) {

	int i, j, k, mv;
	double r[3];
	double U_d, F_L;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;

	const double BET[] = {BETA};
	double dt2beta = 2.0 * BET[params->which_stage] * params->dt;
	double idt2beta = 1.0 / dt2beta;

	// Bounds for manipulating Eulerian grid quantities
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	int i_end = grid->G_Ie;
	int j_end = grid->G_Je;
	int k_end = grid->G_Ke;

	// Bounds for determining local Lagrangian points
	int Is = i_start;
	int Js = j_start;
	int Ks = k_start;

	int Ie = min(i_end, grid->NX-1);
	int Je = min(j_end, grid->NY-1);
	int Ke = min(k_end, grid->NZ-1);

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;

	double *xu = grid -> xu;
	double *yv = grid -> yv;
	double *zw = grid -> zw;

	double ***u_data = u -> data;
	double ***v_data = v -> data;
	double ***w_data = w -> data;

	double ***fx_IBM = data_bag->vof->fx_IBM;;
	double ***fy_IBM = data_bag->vof->fy_IBM;
	double ***fz_IBM = data_bag->vof->fz_IBM;


	double ***u_rhs = u -> ng_rhs;
	double ***v_rhs = v -> ng_rhs;
	double ***w_rhs = w -> ng_rhs;

	double rho_s = p -> rho_s;
	double Vol_L = p -> Vol_L;
	double M_L = Vol_L;

	double *X = p -> X;
	double *U = p -> U;
	double *Omega = p -> Omega;
	double *F = p -> F;
	double *T = p -> T;

	double *X_L = p -> X_L;
	double *Y_L = p -> Y_L;
	double *Z_L = p -> Z_L;
	int *flag_L = p -> flag_L;

	double ***temp_f = lag -> ng_temp;
	double *Temp_U_L = lag -> Temp_L;
	double *Temp_F_L = lag -> Temp_L;
	#ifdef POST_PROCESS
		double ***fx_IBM = lag -> ng_fx_IBM;
		double ***fy_IBM = lag -> ng_fy_IBM;
	#endif

	#ifdef SQUIRMER_SWIMMER
		// Amplitude of first two modes
		double B1 = params->B1;
		double B2 = params->B2;

		// Initialize angle to swimming direction and swimming velocity
		double cos_theta;
		double sin_theta;
		double u_theta;

		// Vector of swimming direction (i.e. particle orientation)
		double es[3];
		double es_norm;

	#ifdef SWIMMER_VERTICAL
		es[0] = 0.0;
		es[1] = 1.0;
		es[2] = 0.0;
	#elif defined SWIMMER_HORIZONTAL
		es[0] = 1.0;
		es[1] = 0.0;
		es[2] = 0.0;
	#else
		FORI3 es[i] = p->Rotn[i][1]; //  Result of Rotn.ey
	#endif

		es_norm = sqrt(DOT(es,es)); // Norm of swimming direction vector

		// Vector of tangential direction in spherical coordinates associated to the
		// swimmer's reference frame.
		double e_theta[3];
		double e_theta_norm = 0.0;

		// Target vector (partile center to target)
		double et[3];
		FORI3 et[i] = -X[i] + params->target_coord[i];
		double et_norm = sqrt(DOT(et,et)); //Norm

		// Norm of r (vector from center to Lagrangian marker)
		double r_norm = p->R; // Norm of vector r is always = radius of the particle

		// Dot products
		double es_dot_r;
		double es_dot_et = DOT(es,et);

		//
		// Asymmetry
		//

	//Angle from swimming direction to target
	double cos_theta_t = es_dot_et/es_norm/et_norm;
	//Angle from marker point to target projected on plan normal to swimming direction
	double cos_phi;
	// Asymmetry amplitude function
	double A = params->target_on * (params->coefA + (1-params->coefA) * erf( 2.0 * (1 - cos_theta_t) ));

		// Projection of et onto the plane normal to the swimming direction
		double et_s[3];
		FORI3 et_s[i] = et[i] - es_dot_et / (es_norm * es_norm) * es[i];
		double et_s_norm = sqrt(DOT(et_s,et_s));

		// Projection of r onto the plane normal to the swimming direction
		double r_s[3];
		double r_s_norm;

		// Threshold to avoid division by zero
		double r_s_norm_threshold = 0.001 * r_norm;
	#endif

	#ifdef VOF_IBM
		/* --- new density interpolation ------------------------------------- */
		double *Temp_Rho_L = lag->Temp_L_rho;  /* per-marker mixture density    */

		/* cell-centred mixture density already computed by VoF code */
		double ***rho_cc = data_bag->vof->rho;   /* Eulerian cell-centred mixture density */

		/* interpolate ρ̃ from Eulerian cells to all *local* Lagrangian markers */
		Interpolate_Eul_to_Lag(rho_cc,      /* source field                    */
							Temp_Rho_L,  /* destination                     */
							'c',         /* dummy flag for cell centres     */
							p, grid);

	#endif


	/*------------------------------------------------------------------------*/
	/*
	 u-momentum
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// 1. Interpolate u-velocity onto Lagrangian markers
	//--------------------------------------------------------------------------
	Interpolate_Eul_to_Lag(u_data, Temp_U_L, 'u', p, grid);

	//--------------------------------------------------------------------------
	// 2. Calculate Lagrangian marker forces
	//--------------------------------------------------------------------------
	for (mv = 0; mv < p->N_L_local; mv++) {

		// Skip the forcing for flagged Lagrangian markers
		if (flag_L[mv] == 0) {
			Temp_F_L[mv] = 0.0;
			continue;
		}

		#ifdef LEFT_WALL_VELOCITY_FREESLIP
		// SKIP ALL MARKERS OUTSIDE THE PHYSICAL DOMAIN!
		// The parity logic in the interpolator handles their effect automatically.
		if (X_L[mv] < 0.0) {
			Temp_F_L[mv] = 0.0;
			continue;
		}
        #endif

		r[0] = X_L[mv] - X[0];
		r[1] = Y_L[mv] - X[1];
		r[2] = Z_L[mv] - X[2];

		// Desired velocity of Lagrangian point
		U_d = U[0] + Omega[1] * r[2] - Omega[2] * r[1];

		#ifdef SQUIRMER_SWIMMER
				// Projection of swimming direction onto marker point position vector
				es_dot_r = DOT(es,r);

				// Angle from marker point to swimming direction
				cos_theta = es_dot_r / (es_norm * r_norm);
				sin_theta = sqrt(es_norm*es_norm*r_norm*r_norm-es_dot_r*es_dot_r)/(es_norm*r_norm);

				// Tangential velocity in Swimmer spherical reference frame
				u_theta = B1*sin_theta + B2*sin_theta*cos_theta;

				//------------
				// non-symmetry
				//------------
				// Projection of r onto the plane normal to the swimming direction
				FORI3 r_s[i] = r[i] - es_dot_r/(es_norm*es_norm)*es[i];
				r_s_norm = sqrt(DOT(r_s,r_s));

				if (r_s_norm < r_s_norm_threshold) {
					cos_phi = 0.0;
				}
				else if(et_s_norm < r_s_norm_threshold) {
					cos_phi = 0.0;
				}
				else {
					// Angle to target
					cos_phi = DOT(r_s,et_s)/(r_s_norm*et_s_norm);
				}

				// Apply function
				u_theta = u_theta*sqrt(1-A*cos_phi);

			//		//Progressive inccrease in swimming velocity
			//		u_theta = u_theta * erf(params->time);

					//-----------
					// e_theta
					//-----------
					// Calculating e_theta in lab reference frame
					FORI3 e_theta[i] = r_norm*r_norm*es[i] - es_dot_r*r[i];
					e_theta_norm = sqrt(DOT(e_theta,e_theta));
					FORI3 e_theta[i] = e_theta[i]/e_theta_norm;

					// Adds the u-momentum component of the squiremer velocity to the desired velocity
					U_d = U_d - u_theta*e_theta[0];
		#endif

		/* NEW: convert to momentum forcing with local mixture density */
		#ifdef VOF_IBM
			double rho_L = Temp_Rho_L[mv];
		#else
			double rho_L = 1.0;
		#endif



		// Acceleration required to produce desired velocity at Lagrangian point
		F_L = idt2beta * (U_d - Temp_U_L[mv]);

		//----------------------------------------------------------------------
		// 3. Add force to particle force if local marker point.  Local marker
		//    points are those located on this processor.  This prevents double
		//    counting when we communicate and add these forces later.
		//----------------------------------------------------------------------
		if (LOCAL_POINT(X_L[mv], Y_L[mv], Z_L[mv])) {
			F[0] -= rho_L * M_L * F_L;  // Density-weighted reaction force
			T[1] -= rho_L * M_L * r[2] * F_L;
			T[2] += rho_L * M_L * r[1] * F_L;
		}

		// Store pure velocity forcing
		Temp_F_L[mv] = F_L;
	}

	//--------------------------------------------------------------------------
	// 4. Spread Lagrangian marker acceleration onto flow field
	//--------------------------------------------------------------------------
	Interpolate_Lag_to_Eul(Temp_F_L, temp_f, 'u', p, grid);
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
				#ifdef VOF_IBM
				double rho_face = 0.5 * (data_bag->vof->rho[k][j][i] + data_bag->vof->rho[k][j][i-1]);
				fx_IBM[k][j][i] = 2.0 * rho_face * temp_f[k][j][i];
				u_rhs[k][j][i] += fx_IBM[k][j][i];
				#else
				u_rhs[k][j][i] += 2.0 * temp_f[k][j][i];
				#endif
			}
		}
	}
#ifdef POST_PROCESS
	#ifdef STARTUP
		// to exclude fibm forces due to the fixed particles, acting as walls. (only in volume imposed and pressure imposed setups)
		if ((p_type == MOBILE) || (params->startup_init != STUP_INIT_PRIMPOSED &&  params->startup_init != STUP_INIT_SHEARFLOW)){
			for (k = k_start; k < k_end; k++) {
				for (j = j_start; j < j_end; j++) {
					for (i = i_start; i < i_end; i++) {
						fx_IBM[k][j][i] += 2.0 * BET[params->which_stage] * temp_f[k][j][i];
					}
				}
			}			
		}
	#else
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					fx_IBM[k][j][i] += 2.0 * BET[params->which_stage] * temp_f[k][j][i];
				}
			}
		}
	#endif
#endif


	/*------------------------------------------------------------------------*/
	/*
	 v-momentum
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// 1. Interpolate v-velocity onto Lagrangian markers
	//--------------------------------------------------------------------------
	Interpolate_Eul_to_Lag(v->data, Temp_U_L, 'v', p, grid);

	//--------------------------------------------------------------------------
	// 2. Calculate Lagrangian marker forces
	//--------------------------------------------------------------------------
	for (mv = 0; mv < p->N_L_local; mv++) {

		// Skip the forcing for flagged Lagrangian markers
		if (flag_L[mv] == 0) {
			Temp_F_L[mv] = 0.0;
			continue;
		}

		#ifdef LEFT_WALL_VELOCITY_FREESLIP
		// SKIP ALL MARKERS OUTSIDE THE PHYSICAL DOMAIN!
		// The parity logic in the interpolator handles their effect automatically.
		if (X_L[mv] < 0.0) {
			Temp_F_L[mv] = 0.0;
			continue;
		}
        #endif

		r[0] = X_L[mv] - X[0];
		r[1] = Y_L[mv] - X[1];
		r[2] = Z_L[mv] - X[2];

		// Desired velocity of Lagrangian point
		U_d = U[1] + Omega[2] * r[0] - Omega[0] * r[2];

		#ifdef SQUIRMER_SWIMMER
					// Projection of swimming direction onto marker point position vector
					es_dot_r = DOT(es,r);

					// Angle from marker point to swimming direction
					cos_theta = es_dot_r / (es_norm * r_norm);
					sin_theta = sqrt(es_norm*es_norm*r_norm*r_norm-es_dot_r*es_dot_r)/(es_norm*r_norm);

					// Tangential velocity in Swimmer spherical reference frame
					u_theta = B1*sin_theta + B2*sin_theta*cos_theta;

					//------------
					// non-symmetry
					//------------
					// Projection of r onto the plane normal to the swimming direction
					FORI3 r_s[i] = r[i] - es_dot_r/(es_norm*es_norm)*es[i];
					r_s_norm = sqrt(DOT(r_s,r_s));

					if (r_s_norm < r_s_norm_threshold) {
						cos_phi = 1.0;
					}
					else if(et_s_norm < r_s_norm_threshold) {
						cos_phi = 1.0;
					}
					else {
					// Angle to target
					cos_phi = DOT(r_s,et_s)/(r_s_norm*et_s_norm);
					}

					// Apply function
					u_theta = u_theta*sqrt(1-A*cos_phi);

			//		//Progressive inccrease in swimming velocity
			//		u_theta = u_theta * erf(params->time);

					//-----------
					// e_theta
					//-----------
					// Calculating e_theta in lab reference frame
					FORI3 e_theta[i] = r_norm*r_norm*es[i] - es_dot_r*r[i];
					e_theta_norm = sqrt(DOT(e_theta,e_theta));
					FORI3 e_theta[i] = e_theta[i]/e_theta_norm;

					// Adds the v-momentum component of the squiremer velocity to the desired velocity
					U_d = U_d - u_theta*e_theta[1];
		#endif

		/* NEW: convert to momentum forcing with local mixture density */
		#ifdef VOF_IBM
			double rho_L = Temp_Rho_L[mv];
		#else
			double rho_L = 1.0;
		#endif

		// Acceleration required to produce desired velocity at Lagrangian point 
		F_L = idt2beta * (U_d - Temp_U_L[mv]);

		//----------------------------------------------------------------------
		// 3. Add force to particle force if local marker point.  Local marker
		//    points are those located on this processor.  This prevents double
		//    counting when we communicate and add these forces later.
		//----------------------------------------------------------------------
		if (LOCAL_POINT(X_L[mv], Y_L[mv], Z_L[mv])) {
			F[1] -= rho_L * M_L * F_L;  // Density-weighted reaction force
			T[0] += rho_L * M_L * r[2] * F_L;
			T[2] -= rho_L * M_L * r[0] * F_L;
		}

		// Store pure velocity forcing
		Temp_F_L[mv] = F_L;
	}

	//--------------------------------------------------------------------------
	// 4. Spread Lagrangian marker acceleration onto flow field
	//--------------------------------------------------------------------------
	Interpolate_Lag_to_Eul(Temp_F_L, temp_f, 'v', p, grid);
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
				#ifdef VOF_IBM
				double rho_face = 0.5 * (data_bag->vof->rho[k][j][i] + data_bag->vof->rho[k][j-1][i]);
				fy_IBM[k][j][i] = 2.0 * rho_face * temp_f[k][j][i];
				v_rhs[k][j][i] += fy_IBM[k][j][i];
				#else
				v_rhs[k][j][i] += 2.0 * temp_f[k][j][i];
				#endif
			}
		}
	}
#ifdef POST_PROCESS
	#ifdef STARTUP
		// to exclude fibm forces due to the fixed particles, acting as walls. (only in volume imposed and pressure imposed setups)
		if ((p_type == MOBILE) || (params->startup_init != STUP_INIT_PRIMPOSED &&  params->startup_init != STUP_INIT_SHEARFLOW)){
			for (k = k_start; k < k_end; k++) {
				for (j = j_start; j < j_end; j++) {
					for (i = i_start; i < i_end; i++) {
						fy_IBM[k][j][i] += 2.0 * BET[params->which_stage] * temp_f[k][j][i];
					}
				}
			}			
		}
	#else
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					fy_IBM[k][j][i] += 2.0 * BET[params->which_stage] * temp_f[k][j][i];
				}
			}
		}
	#endif
#endif


	/*------------------------------------------------------------------------*/
	/*
	 w-momentum
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// 1. Interpolate w-velocity onto Lagrangian markers
	//--------------------------------------------------------------------------
	Interpolate_Eul_to_Lag(w->data, Temp_U_L, 'w', p, grid);

	//--------------------------------------------------------------------------
	// 2. Calculate Lagrangian marker forces
	//--------------------------------------------------------------------------
	for (mv = 0; mv < p->N_L_local; mv++) {

		// Skip the forcing for flagged Lagrangian markers
		if (flag_L[mv] == 0) {
			Temp_F_L[mv] = 0.0;
			continue;
		}

		#ifdef LEFT_WALL_VELOCITY_FREESLIP
		// SKIP ALL MARKERS OUTSIDE THE PHYSICAL DOMAIN!
		// The parity logic in the interpolator handles their effect automatically.
		if (X_L[mv] < 0.0) {
			Temp_F_L[mv] = 0.0;
			continue;
		}
        #endif

		r[0] = X_L[mv] - X[0];
		r[1] = Y_L[mv] - X[1];
		r[2] = Z_L[mv] - X[2];

		// Desired velocity of Lagrangian point
		U_d = U[2] + Omega[0] * r[1] - Omega[1] * r[0];

		#ifdef SQUIRMER_SWIMMER
					// Projection of swimming direction onto marker point position vector
					es_dot_r = DOT(es,r);

					// Angle from marker point to swimming direction
					cos_theta = es_dot_r / (es_norm * r_norm);
					sin_theta = sqrt(es_norm*es_norm*r_norm*r_norm-es_dot_r*es_dot_r)/(es_norm*r_norm);

					// Tangential velocity in Swimmer spherical reference frame
					u_theta = B1*sin_theta + B2*sin_theta*cos_theta;

					//------------
					// non-symmetry
					//------------
					// Projection of r onto the plane normal to the swimming direction
					FORI3 r_s[i] = r[i] - es_dot_r/(es_norm*es_norm)*es[i];
					r_s_norm = sqrt(DOT(r_s,r_s));

					if (r_s_norm < r_s_norm_threshold) {
						cos_phi = 1.0;
					}
					else if(et_s_norm < r_s_norm_threshold) {
						cos_phi = 1.0;
					}
					else {
					// Angle to target
					cos_phi = DOT(r_s,et_s)/(r_s_norm*et_s_norm);
					}

					// Apply function
					u_theta = u_theta*sqrt(1-A*cos_phi);

			//		//Progressive inccrease in swimming velocity
			//		u_theta = u_theta * erf(params->time);

					//-----------
					// e_theta
					//-----------
					// Calculating e_theta in lab reference frame
					FORI3 e_theta[i] = r_norm*r_norm*es[i] - es_dot_r*r[i];
					e_theta_norm = sqrt(DOT(e_theta,e_theta));
					FORI3 e_theta[i] = e_theta[i]/e_theta_norm;

					// Adds the w-momentum component of the squiremer velocity to the desired velocity
					U_d = U_d - u_theta*e_theta[2];
		#endif

		/* NEW: convert to momentum forcing with local mixture density */
		#ifdef VOF_IBM
			double rho_L = Temp_Rho_L[mv];
		#else
			double rho_L = 1.0;
		#endif

		// Acceleration required to produce desired velocity at Lagrangian point 
		F_L = idt2beta * (U_d - Temp_U_L[mv]);

		//----------------------------------------------------------------------
		// 3. Add force to particle force if local marker point.  Local marker
		//    points are those located on this processor.  This prevents double
		//    counting when we communicate and add these forces later.
		//----------------------------------------------------------------------
		if (LOCAL_POINT(X_L[mv], Y_L[mv], Z_L[mv])) {
			F[2] -= rho_L * M_L * F_L;  // Density-weighted reaction force
			T[0] -= rho_L * M_L * r[1] * F_L;
			T[1] += rho_L * M_L * r[0] * F_L;
		}

		// Store pure velocity forcing
		Temp_F_L[mv] = F_L;
	}

	//--------------------------------------------------------------------------
	// 4. Spread Lagrangian marker acceleration onto flow field
	//--------------------------------------------------------------------------
	Interpolate_Lag_to_Eul(Temp_F_L, temp_f, 'w', p, grid);
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
				#ifdef VOF_IBM
				double rho_face = 0.5 * (data_bag->vof->rho[k][j][i] + data_bag->vof->rho[k-1][j][i]);
				fz_IBM[k][j][i] = 2.0 * rho_face * temp_f[k][j][i];
				w_rhs[k][j][i] += fz_IBM[k][j][i];
				#else
				w_rhs[k][j][i] += 2.0 * temp_f[k][j][i];
				#endif
			}
		}
	}

}



/******************************************************************************/
/*
 * Generate Lagrangian marker points for particles and then flag them to decide
 * whether or not they will be used for evaluating and spreading IBM forces.
 */
/******************************************************************************/
void Lagrangian_flag_points(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	Particle *p, *p2;
	int mv;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	/*------------------------------------------------------------------------*/
	/*
	 Flag mobile particle Lagrangian points
	 */
	/*------------------------------------------------------------------------*/
	p = p_mobile_list -> start;
	while (p != NULL) {

		// Generate Lagrangian marker points
		Lagrangian_generate_points(p, grid);

		// Reset Lagrangian marker switch to 'on'
		for (mv = 0; mv < p->N_L_local; mv++) {
			p -> flag_L[mv] = 1;
		}

#if defined LAG_MARKER_FLAG || defined LAG_MARKER_PRIORITY
		//----------------------------------------------------------------------
		// Other mobile particles
		//----------------------------------------------------------------------
		p2 = p_mobile_list -> start;
		while (p2 != NULL) {
			if (p2 -> ID != p -> ID)
				Lagrangian_flag_points_individual(p, p2, grid);
			p2 = p2 -> next;
		}

		//----------------------------------------------------------------------
		// Fixed particles
		//----------------------------------------------------------------------
		p2 = p_fixed_list -> start;
		while (p2 != NULL) {
			Lagrangian_flag_points_individual(p, p2, grid);
			p2 = p2 -> next;
		}

		//----------------------------------------------------------------------
		// Walls
		//----------------------------------------------------------------------
		Lagrangian_flag_points_wall(p, grid);
#endif

		p = p -> next;
	}

	/*------------------------------------------------------------------------*/
	/*
	 Flag fixed particle Lagrangian points
	 */
	/*------------------------------------------------------------------------*/
	p = p_fixed_list -> start;
	while (p != NULL) {

		// Generate Lagrangian marker points
		Lagrangian_generate_points(p, grid);

		// Reset Lagrangian marker switch to 'on'
		for (mv = 0; mv < p->N_L_local; mv++) {
			p -> flag_L[mv] = 1;
		}

#if defined LAG_MARKER_FLAG || defined LAG_MARKER_PRIORITY
		//----------------------------------------------------------------------
		// Mobile particles
		//----------------------------------------------------------------------
		p2 = p_mobile_list -> start;
		while (p2 != NULL) {
			Lagrangian_flag_points_individual(p, p2, grid);
			p2 = p2 -> next;
		}

		//----------------------------------------------------------------------
		// Other fixed particles
		//----------------------------------------------------------------------
		p2 = p_fixed_list -> start;
		while (p2 != NULL) {
			if (p2 -> ID != p -> ID)
				Lagrangian_flag_points_individual(p, p2, grid);
			p2 = p2 -> next;
		}

		//----------------------------------------------------------------------
		// Walls
		//----------------------------------------------------------------------
		Lagrangian_flag_points_wall(p, grid);
#endif

		p = p -> next;
	}
}




/******************************************************************************/
/*
 Assumes uniform grid
 */
/******************************************************************************/
void Lagrangian_flag_points_individual(Particle *p, Particle *p2, MAC_grid *grid) {

	int mv;
	double r[3], dist_p_p2, dist_L_p2;
	double h = grid -> dx_u[1];

	double *X = p -> X;
	double  R = p -> R;

	double *X2 = p2 -> X;
	double  R2 = p2 -> R;

	double *X_L = p -> X_L;
	double *Y_L = p -> Y_L;
	double *Z_L = p -> Z_L;
	int *flag_L = p -> flag_L;
	int N_L_local = p -> N_L_local;

#ifdef LAG_MARKER_PRIORITY
	// Particles with smaller IDs keep their markers in case of overlap
	if (p->ID < p2->ID) return;
#endif

	// Distance threshold from center of 'p' to center of 'p2'
	dist_p_p2 = R + R2 + LAG_FLAG_RANGE * h;
	dist_p_p2 = dist_p_p2 * dist_p_p2;

	// Distance threshold from Lagrangian marker on 'p' to center of 'p2'
	dist_L_p2 = R2 + LAG_FLAG_RANGE * h;
	dist_L_p2 = dist_L_p2 * dist_L_p2;

	// Check if other particle is close enough to turn off flags
	r[0] = X2[0] - X[0];
	r[1] = X2[1] - X[1];
	r[2] = X2[2] - X[2];

	if (DOT(r,r) < dist_p_p2) {

		for (mv = 0; mv < N_L_local; mv++) {

			// Check if Lagrangian marker is close enough to turn off its flag
			r[0] = X2[0] - X_L[mv];
			r[1] = X2[1] - Y_L[mv];
			r[2] = X2[2] - Z_L[mv];

			if (flag_L[mv] == 1 && DOT(r,r) < dist_L_p2) {
				flag_L[mv] = 0;
			}
		}
	}
}




/******************************************************************************/
/*
 Assumes uniform grid
 */
/******************************************************************************/
void Lagrangian_flag_points_wall(Particle *p, MAC_grid *grid) {

	int mv;
	double dist_p_wall, dist_L_wall;
	double h = grid -> dx_u[1];

	double *X = p -> X;
	double  R = p -> R;

	double *X_L = p -> X_L;
	double *Y_L = p -> Y_L;
	double *Z_L = p -> Z_L;
	int *flag_L = p -> flag_L;
	int N_L_local = p -> N_L_local;

	// Distance threshold from center of 'p' to wall
	dist_p_wall = R + 0.5*LAG_FLAG_RANGE * h;

	// Distance threshold from Lagrangian marker on 'p' to wall
	dist_L_wall = 0.5*LAG_FLAG_RANGE * h;

#ifndef XPERIODIC
	//--------------------------------------------------------------------------
	// Lower x wall
	//--------------------------------------------------------------------------
	double xmin = grid -> xu[0];
	if (fabs(X[0] - xmin) < dist_p_wall) {

		for (mv = 0; mv < N_L_local; mv++) {

			if (flag_L[mv] == 1 && fabs(X_L[mv] - xmin) < dist_L_wall) {
				flag_L[mv] = 0;
			}
		}
	}

	//--------------------------------------------------------------------------
	// Upper x wall
	//--------------------------------------------------------------------------
	double xmax = grid -> xu[grid -> NX - 1];
	if (fabs(X[0] - xmax) < dist_p_wall) {

		for (mv = 0; mv < N_L_local; mv++) {

			if (flag_L[mv] == 1 && fabs(X_L[mv] - xmax) < dist_L_wall) {
				flag_L[mv] = 0;
			}
		}
	}
#endif

#ifndef YPERIODIC

	//--------------------------------------------------------------------------
	// Lower y wall
	//--------------------------------------------------------------------------
	double ymin = grid -> yv[0];
	if (fabs(X[1] - ymin) < dist_p_wall) {

		for (mv = 0; mv < N_L_local; mv++) {

			if (flag_L[mv] == 1 && fabs(Y_L[mv] - ymin) < dist_L_wall) {
				flag_L[mv] = 0;
			}
		}
	}

	//--------------------------------------------------------------------------
	// Upper y wall
	//--------------------------------------------------------------------------
	double ymax = grid -> yv[grid -> NY - 1];
	if (fabs(X[1] - ymax) < dist_p_wall) {

		for (mv = 0; mv < N_L_local; mv++) {

			if (flag_L[mv] == 1 && fabs(Y_L[mv] - ymax) < dist_L_wall) {
				flag_L[mv] = 0;
			}
		}
	}
#endif
#ifndef ZPERIODIC
	//--------------------------------------------------------------------------
	// Lower z wall
	//--------------------------------------------------------------------------
	double zmin = grid -> zw[0];
	if (fabs(X[2] - zmin) < dist_p_wall) {

		for (mv = 0; mv < N_L_local; mv++) {

			if (flag_L[mv] == 1 && fabs(Z_L[mv] - zmin) < dist_L_wall) {
				flag_L[mv] = 0;
			}
		}
	}

	//--------------------------------------------------------------------------
	// Upper z wall
	//--------------------------------------------------------------------------
	double zmax = grid -> zw[grid -> NZ - 1];
	if (fabs(X[2] - zmax) < dist_p_wall) {

		for (mv = 0; mv < N_L_local; mv++) {

			if (flag_L[mv] == 1 && fabs(Z_L[mv] - zmax) < dist_L_wall) {
				flag_L[mv] = 0;
			}
		}
	}
#endif
}





void find_dx (double x_goal, double *x, double *dx, int *index, int NX) {

	int low_idx = 0;
	int high_idx = NX;

	int border = -1;
	int border_old = -1;

	int Max_Iter = 500;
	int iter_count = 0;

	while(iter_count < Max_Iter){

		iter_count += 1;
		border_old = border;
		border = round((high_idx + low_idx) / 2);

		if(border_old == border) break;

		if (x_goal <= x[border])
			high_idx = border;

		else
			low_idx = border;

	}

	*index = border;

	if(x[*index] > x_goal)
		index = index - 1;

	*dx = x[*index + 1] - x[*index];


	if(iter_count == Max_Iter) {

		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("!!!!!! ERROR IN void find_dx (BINARY SEARCH) !!!!!!\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}

	return;
}



double *mat_mat(int m, int n, int p, int q, double **A, double **B) { //double A[][n], double B[][q]){

	double *C = calloc(m*q, sizeof(double));
	//
	//	for (int i=0; i<n; i++) {
	//		C[i] = calloc(q, sizeof(double));
	//	}


	double sum = 0;
	int cnt=0;
	int i,j,k;


	for (j = 0; j < q; j++) {			// switch i and j loops for row-major order
		for (i = 0; i < m; i++) {
			for (k = 0; k < p; k++) {
				sum += A[i][k]*B[k][j];
			}

			//C[j*m + i] = sum;
			C[cnt] = sum;
			cnt++;
			sum = 0;
		}
	}

	return C;
}







/******************************************************************************/
/*
 Generates 'N_L' points evenly distributed over the surface of a sphere of
 radius 'R' centered on 'X' = {xc, yc, zc}.
 Based on the algorithm of Paul Leopardi: 'A partition of the unit sphere into
 regions of equal area and small diameter,' Electronic Transactions on Numerical
 Analysis, 2006.
 */
/******************************************************************************/
void Lagrangian_generate_points_Leopardi(Particle *p, MAC_grid *grid) {

	int i, j;
	double A_r, A, AC, A_F, A_F_last;
	double delta_I, delta_F;
	double n_I, m_I, m_diff;
	double theta_c, theta_F, theta, phi, delta_phi;
	double R_sin_theta, R_cos_theta, R_sin_theta_h, R_cos_theta_h;
	double x, y, z, x_h, y_h, z_h;
	int n, m;

	double PI4 = 4.0 * PI;

	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	// Evaluating xu[Ie], yv[Je], zw[Ke], not i < Ie
	int Ie = min(grid -> G_Ie, grid -> NX - 1);
	int Je = min(grid -> G_Je, grid -> NY - 1);
	int Ke = min(grid -> G_Ke, grid -> NZ - 1);

	double *xu = grid -> xu;
	double *yv = grid -> yv;
	double *zw = grid -> zw;

	double hgrid = grid -> dx_u[1];


	double R = p -> R;
#ifdef RETRACTION
	R -= RETRACTION*hgrid;
#endif

	double *X = p -> X;
	double *X_L = p -> X_L;
	double *Y_L = p -> Y_L;
	double *Z_L = p -> Z_L;

#ifdef IBM_SCALAR
	double *X_H = p -> X_H;
	double *Y_H = p -> Y_H;
	double *Z_H = p -> Z_H;
#endif
	// Number of Lagrangian marker points
	int N_L = p -> N_L;
	int N_L_local = 0;

	// Range within which Lagrangian points will influence local processor
	// velocity nodes
	double range = DELTA_FUNC_RADIUS * grid -> dx_u[1];

	// Logic check
	if (N_L < 1)
		fprintf(stderr, "Error: invalid number of points '%d' specified in function "
				"generate_points().  Should have N_L > 0\n", N_L);

	/*------------------------------------------------------------------------*/
	/*
	 This algorithm works by dividing the surface of the sphere into 'N_L'
	 regions of equal area, denoted here as 'A_r'.
	 The first two regions are the north and south polar caps.  The rest of the
	 sphere is divided into 'n' rings, or annuli.  Each annulus has a north and
	 south colatitude (polar angle).
	 Each annulus is furthermore divided into 'm' regions having an area of
	 exactly 'A_r'.
	 The points are placed in the area centroid of the regions (as opposed to
	 the colatitude centroid).
	 We use the physics system for spherical coordinate notation:
	 - theta: polar angle (angle from z-axis)
	 - phi: azimuthal angle (angle in x-y plane)
	 */
	/*------------------------------------------------------------------------*/

	// Area of each region
	A_r = 4.0 * PI / N_L;

	// Colatitude of north polar cap
	theta_c = 2.0 * asin(sqrt(A_r / (4.0 * PI)));

	// Ideal annulus angle
	delta_I = sqrt(A_r);

	// Ideal number of annuli
	n_I = (PI - 2.0 * theta_c) / delta_I;

	// Actual number of annuli
	if (N_L < 3)
		n = 0;
	else
		n = max(1, (int)round(n_I));

	//--------------------------------------------------------------------------
	// Coordinates of points in north polar cap
	//--------------------------------------------------------------------------
	x = X[0];
	y = X[1];
	z = X[2] + R;

	// See if coordinate is within processor range of influence
	if (NEAR_POINT(x,y,z,range)) {
		X_L[N_L_local] = x;
		Y_L[N_L_local] = y;
		Z_L[N_L_local] = z;

#ifdef IBM_SCALAR
		X_H[N_L_local] = x;
		Y_H[N_L_local] = y;
		Z_H[N_L_local] = z+hgrid;
#endif

		N_L_local++;
	}

	//--------------------------------------------------------------------------
	// Coordinates of points in south polar cap
	//--------------------------------------------------------------------------
	if (N_L > 1) {
		x = X[0];
		y = X[1];
		z = X[2] - R;

		// See if coordinate is within processor range of influence
		if (NEAR_POINT(x,y,z,range)) {
			X_L[N_L_local] = x;
			Y_L[N_L_local] = y;
			Z_L[N_L_local] = z;
#ifdef IBM_SCALAR
			X_H[N_L_local] = x;
			Y_H[N_L_local] = y;
			Z_H[N_L_local] = z-hgrid;
#endif
			N_L_local++;
		}
	}

	//--------------------------------------------------------------------------
	// If we only have two points, then we are done!
	//--------------------------------------------------------------------------
	if (N_L < 3) {
		p -> N_L_local = N_L_local;
		return;
	}

	//--------------------------------------------------------------------------
	// If we have more than two points, then we will place the additional ones
	// within the 'n' annuli.
	//--------------------------------------------------------------------------

	// "Fitting" annulus angle that would allow the 'n' annuli to fit between
	// the polar caps
	delta_F = (PI - 2.0 * theta_c) / n;

	// Parameter initialization for annulus loop
	theta_F = theta_c;
	A_F_last = A_r;
	m_diff = 0.0;
	AC = A_r;

	// Initialize phi for first annulus
	phi = 0.0;

	//--------------------------------------------------------------------------
	// Here we loop through each annulus, determining how many regions would fit
	// within the ideal annulus size.
	//--------------------------------------------------------------------------
	for (i = 0; i < n; i++) {

		// Southern colatitude of "fitting" annulus
		theta_F += delta_F;

		// Area from north pole to 'theta_F'
		A_F = sin(0.5 * theta_F);
		A_F = 4.0 * PI * A_F * A_F;

		// Ideal number of regions in this annulus (area of annulus divided by
		// the area of a single region)
		m_I = (A_F - A_F_last) / A_r;
		A_F_last = A_F;

		// Actual number of regions in this annulus
		m = (int)round(m_I + m_diff);
		m_diff += m_I - m;

		// Area of this annulus
		A = A_r * m;

		// Actual area from north pole to the center of this annulus -- center
		// in terms of area, not angle, which would have been
		//     (theta_F + theta_F_last) / 2
		AC += 0.5 * A;

		/*--------------------------------------------------------------------*/
		/*
		 Calculate point coordinates
		 Here we loop through all the regions within this annulus using the
		 polar coordinate 'phi', while 'theta' will be our azimuthal angle.
		 We could calculate the central colatitude in this annulus
		     theta = 2.0 * asin(sqrt(AC / (4.0 * PI)));
		 and then get the values
		     R_sin_theta = R * sin(theta);
		     R_cos_theta = R * cos(theta);
		 Or we could evaluate
		     sin( 2.0 * asin( sqrt(AC / (4.0 * PI)) ) )
		 by hand to get the following:
		 */
		/*--------------------------------------------------------------------*/
		delta_phi = 2.0 * PI / m;
		R_sin_theta = R * 2.0 * sqrt( AC * ( PI4 - AC ) ) / PI4;
		R_cos_theta = R * ( PI4 - 2.0 * AC ) / PI4;

		R_sin_theta_h = (R+hgrid)* 2.0 * sqrt( AC * ( PI4 - AC ) ) / PI4;
		R_cos_theta_h = (R+hgrid) * ( PI4 - 2.0 * AC ) / PI4;


		for (j = 0; j < m; j++) {

			// Coordinates of points within this annulus
			x = X[0] + R_sin_theta * cos(phi);
			y = X[1] + R_sin_theta * sin(phi);
			z = X[2] + R_cos_theta;

#ifdef IBM_SCALAR
			x_h = X[0] + R_sin_theta_h * cos(phi);
			y_h = X[1] + R_sin_theta_h * sin(phi);
			z_h = X[2] + R_cos_theta_h;
#endif
			// See if coordinate is within processor range of influence
			if (NEAR_POINT(x,y,z,range)) {
				X_L[N_L_local] = x;
				Y_L[N_L_local] = y;
				Z_L[N_L_local] = z;
#ifdef IBM_SCALAR
				X_H[N_L_local] = x_h;
				Y_H[N_L_local] = y_h;
				Z_H[N_L_local] = z_h;
#endif
				N_L_local++;
			}

			phi += delta_phi;
		}

		AC += 0.5 * A;

		// Set up phi to be offset for next annulus
		phi -= 0.5 * delta_phi;
	}

	p -> N_L_local = N_L_local;
}




#ifndef GRID_UNIFORM
/******************************************************************************/
/*
 Generates 'N_L' points evenly distributed over the surface of a sphere of
 radius 'R' centered on 'X' = {xc, yc, zc}.
 Based on the algorithm of Paul Leopardi: 'A partition of the unit sphere into
 regions of equal area and small diameter,' Electronic Transactions on Numerical
 Analysis, 2006.
 */
/******************************************************************************/
void Lagrangian_generate_points_Bala(Particle *p, MAC_grid *grid) {


	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	// Evaluating xu[Ie], yv[Je], zw[Ke], not i < Ie
	int Ie = min(grid -> G_Ie, grid -> NX - 1);
	int Je = min(grid -> G_Je, grid -> NY - 1);
	int Ke = min(grid -> G_Ke, grid -> NZ - 1);

	double *xu = grid -> xu;
	double *yv = grid -> yv;
	double *zw = grid -> zw;

	double R = p -> R;
	double *X = p -> X;
	double *X_L = p -> X_L;
	double *Y_L = p -> Y_L;
	double *Z_L = p -> Z_L;


	double *N_l_subsec = p -> N_l_subsec;
	int *idx_subsec_low = p -> idx_subsec_low;

	int *idx_x_markers = p -> idx_x_markers;
	int *idx_y_markers = p -> idx_y_markers;
	int *idx_z_markers = p -> idx_z_markers;

	double **A = p -> A;
	double **B = p -> B;

	double sin_theta, cos_theta, sin_phi, cos_phi;

	// Number of Lagrangian marker points
	int N_L = p -> N_L;
	int N_L_local = 0;
	int mark_count = 0;

	// Range within which Lagrangian points will influence local processor
	// velocity nodes
	double range = DELTA_FUNC_RADIUS * grid -> dx_u[1];  //TODO: check limit

	// Logic check
	if (N_L < 1)
		fprintf(stderr, "Error: invalid number of points specified in function "
				"generate_points().  Should have N_L > 0\n");

	int i, j;

	double y_high = p->X[1] + R;
	double y_low = p->X[1] - R;
	double delta_y;
	int idx;
	find_dx(y_low,yv, &delta_y, &idx, grid -> NY); //TODO: check!

	double y_next = yv[idx+1];


	//--------------------------------------------------------------------------
	// Coordinates of points in south polar cap
	//--------------------------------------------------------------------------

	double theta_k = PI;
	double phi_k = 0;
	double phi_old = 0;
	double x_mark, y_mark, z_mark;

	if (N_L > 1) {
		x_mark = X[0];
		y_mark = X[1] - R;
		z_mark = X[2];


		sin_theta = sin(theta_k);
		cos_theta = cos(theta_k);
		sin_phi = sin(phi_k);
		cos_phi = cos(phi_k);

		A[0][mark_count] = sin_theta * cos_phi;
		A[1][mark_count] = sin_theta * sin_phi;
		A[2][mark_count] = cos_theta;
		A[3][mark_count] = (-sin_theta * cos_phi) * sin_theta * cos_phi;
		A[4][mark_count] = (-sin_theta * cos_phi) * cos_theta;
		A[5][mark_count] = (-sin_theta * sin_phi) * sin_theta * sin_phi;
		A[6][mark_count] = (-cos_theta * cos_phi) * cos_theta * cos_phi -
		sin_phi * sin_phi;
		A[7][mark_count] = (-cos_theta * sin_phi) * cos_theta * cos_phi +
		cos_phi * sin_phi;
		A[8][mark_count] = cos_theta * sin_phi * sin_theta;

		B[mark_count][0] = 1;  // Ac00 * ....
		B[mark_count][1] = cos_theta;  // Ac10 * ....
		B[mark_count][2] = -sin_theta * cos_phi;  // Ac11 * ....
		B[mark_count][3] = -sin_theta * sin_phi;  // As11 * .....
		B[mark_count][4] = 0.5 * (3 * cos_theta * cos_theta - 1);  // Ac20 * ....
		B[mark_count][5] = -3 * sin_theta * cos_theta * cos_phi;  // Ac21 * ....
		B[mark_count][6] = -3 * sin_theta * cos_theta * sin_phi;  // As21 * ....
		B[mark_count][7] = 3 * sin_theta * sin_theta * cos(2 * phi_k);  // Ac22 * ....
		B[mark_count][8] = 3 * sin_theta * sin_theta * sin(2 * phi_k);  // As22 * ....

		mark_count++;

		// See if coordinate is within processor range of influence
		if (NEAR_POINT(x_mark,y_mark,z_mark,range)) {
			X_L[N_L_local] = x_mark;
			Y_L[N_L_local] = y_mark;
			Z_L[N_L_local] = z_mark;
//			phi_markers[N_L_local] = phi_k;
//			theta_markers[N_L_local] = theta_k;
			N_L_local++;
		}
	}

	double eps = 0.00007;
	int k = 1;

	int subsec_count = 0;
	double N_l = PI / 3 * ( 12* R*R / (delta_y * delta_y) +1 );
	N_l_subsec[subsec_count] = N_l;
	idx_subsec_low[0] = 0;


	double h_k = 0;

	while (y_mark < y_high - eps) {

		while (k < N_l && y_mark < y_next) {

			h_k = 2* (k-1) / (N_l-1) - 1;
			theta_k = acos(h_k);

			phi_k = phi_old + 3.809 / sqrt(N_l) / sqrt(1 - h_k * h_k);
			phi_old = phi_k;

			sin_theta = sin(theta_k);
			cos_theta = cos(theta_k);
			sin_phi = sin(phi_k);
			cos_phi = cos(phi_k);

			z_mark = X[2] + R * sin_theta * cos_phi;
			x_mark = X[0] + R * sin_theta * sin_phi;
			y_mark = X[1] + R * cos_theta;

			printf("MARK_COUNT = %d\n",mark_count);

			//------------------------------------------------------------------------
			//---------------------- eqn. (5.19) - (5.27) ----------------------------
			A[0][mark_count] = sin_theta * cos_phi;
			A[1][mark_count] = sin_theta * sin_phi;
			A[2][mark_count] = cos_theta;
			A[3][mark_count] = (-sin_theta * cos_phi) * sin_theta * cos_phi;
			A[4][mark_count] = (-sin_theta * cos_phi) * cos_theta;
			A[5][mark_count] = (-sin_theta * sin_phi) * sin_theta * sin_phi;
			A[6][mark_count] = (-cos_theta * cos_phi) * cos_theta * cos_phi -
			sin_phi * sin_phi;
			A[7][mark_count] = (-cos_theta * sin_phi) * cos_theta * cos_phi +
			cos_phi * sin_phi;
			A[8][mark_count] = cos_theta * sin_phi * sin_theta;


			B[mark_count][0] = 1;  // Ac00 * ....
			B[mark_count][1] = cos_theta;  // Ac10 * ....
			B[mark_count][2] = -sin_theta * cos_phi;  // Ac11 * ....
			B[mark_count][3] = -sin_theta * sin_phi;  // As11 * .....
			B[mark_count][4] = 0.5 * (3 * cos_theta * cos_theta - 1);  // Ac20 * ....
			B[mark_count][5] = -3 * sin_theta * cos_theta * cos_phi;  // Ac21 * ....
			B[mark_count][6] = -3 * sin_theta * cos_theta * sin_phi;  // As21 * ....
			B[mark_count][7] = 3 * sin_theta * sin_theta * cos(2 * phi_k);  // Ac22 * ....
			B[mark_count][8] = 3 * sin_theta * sin_theta * sin(2 * phi_k);  // As22 * ....


			// See if coordinate is within processor range of influence
			if (NEAR_POINT(x_mark,y_mark,z_mark,range)) {
				X_L[N_L_local] = x_mark;
				Y_L[N_L_local] = y_mark;
				Z_L[N_L_local] = z_mark;
				N_L_local++;
			}

			mark_count++;
			k++;
		}

		idx ++;
		y_next = yv[idx+1];
		delta_y = y_next - yv[idx];

		N_l = PI / 3 * ( 12* R*R / (delta_y * delta_y) +1 );
		subsec_count ++;
		N_l_subsec[subsec_count] = N_l;
		idx_subsec_low[subsec_count] = N_L_local;

		k =  ceil((h_k + 1) / 2 * (N_l - 1) + 1);
		eps = 0.07 * delta_y;

	}


	x_mark = X[0];
	y_mark = X[1] + R;
	z_mark = X[2];

	theta_k = 0;
	phi_k = 0;

	sin_theta = sin(theta_k);
	cos_theta = cos(theta_k);
	sin_phi = sin(phi_k);
	cos_phi = cos(phi_k);

	A[0][mark_count] = sin_theta * cos_phi;
	A[1][mark_count] = sin_theta * sin_phi;
	A[2][mark_count] = cos_theta;
	A[3][mark_count] = (-sin_theta * cos_phi) * sin_theta * cos_phi;
	A[4][mark_count] = (-sin_theta * cos_phi) * cos_theta;
	A[5][mark_count] = (-sin_theta * sin_phi) * sin_theta * sin_phi;
	A[6][mark_count] = (-cos_theta * cos_phi) * cos_theta * cos_phi -
	sin_phi * sin_phi;
	A[7][mark_count] = (-cos_theta * sin_phi) * cos_theta * cos_phi +
	cos_phi * sin_phi;
	A[8][mark_count] = cos_theta * sin_phi * sin_theta;

	B[mark_count][0] = 1;  // Ac00 * ....
	B[mark_count][1] = cos_theta;  // Ac10 * ....
	B[mark_count][2] = -sin_theta * cos_phi;  // Ac11 * ....
	B[mark_count][3] = -sin_theta * sin_phi;  // As11 * .....
	B[mark_count][4] = 0.5 * (3 * cos_theta * cos_theta - 1);  // Ac20 * ....
	B[mark_count][5] = -3 * sin_theta * cos_theta * cos_phi;  // Ac21 * ....
	B[mark_count][6] = -3 * sin_theta * cos_theta * sin_phi;  // As21 * ....
	B[mark_count][7] = 3 * sin_theta * sin_theta * cos(2 * phi_k);  // Ac22 * ....
	B[mark_count][8] = 3 * sin_theta * sin_theta * sin(2 * phi_k);  // As22 * ....

	mark_count++;

	// See if coordinate is within processor range of influence
	if (NEAR_POINT(x_mark,y_mark,z_mark,range)) {
		X_L[N_L_local] = x_mark;
		Y_L[N_L_local] = y_mark;
		Z_L[N_L_local] = z_mark;
		N_L_local++;
	}


	double *A_red = mat_mat(9, mark_count, mark_count, 9, A, B);
	double b[9] = {0, 0, 0, -4/3*PI*R*R, 0, -4/3*PI*R*R, -8/3*PI*R*R, 0, 0};
	double *A_coeff = qr_solve (9, 9, A_red, b);

	double *Vl = p -> Vl;
	mat_vec(mark_count, 9, B, A_coeff, Vl);

	for (i=1; i < subsec_count; i++) {
		for (j = idx_subsec_low[i]; j < idx_subsec_low[i+1]; j++) {

			Vl[j] *= sqrt((4*PI*R*R)/N_l_subsec[i]);

		}
	}

	free(A_red);		// TODO: better memory management
	free(A_coeff);

	find_dx(X_L[0],xu, &delta_y, &idx_x_markers[0], grid -> NX);
	find_dx(Y_L[0],yv, &delta_y, &idx_y_markers[0], grid -> NY);
	find_dx(Z_L[0],zw, &delta_y, &idx_z_markers[0], grid -> NZ);

	for (i=1; i<mark_count; i++) {

		idx = idx_x_markers[i-1];
		while(xu[idx] < X_L[i]) idx++;
		while(xu[idx] > X_L[i]) idx--;
		idx_x_markers[i] = idx;

		idx = idx_y_markers[i-1];
		while(yv[idx] < Y_L[i]) idx++;
		while(yv[idx] > Y_L[i]) idx--;
		idx_y_markers[i] = idx;

		idx = idx_z_markers[i-1];
		while(zw[idx] < Z_L[i]) idx++;
		while(zw[idx] > Z_L[i]) idx--;
		idx_z_markers[i] = idx;
	}


	p -> N_L_local = N_L_local;
}
#endif

#ifdef IBM_SCALAR
void Lagrangian_heat(int iconc, int heat_iter, Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	int i, iters, corrector;
	Particle *p;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;
	Concentration *c = data_bag -> c[iconc];



	//double *c_rhs_1d = &(c -> ng_rhs[Ks][Js][Is]);

	//	double *c_implicit_1d = &(c -> ng_implicit[Ks][Js][Is]);

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	//--------------------------------------------------------------------------
	// 1. Forcing before implicit diffusion terms
	//--------------------------------------------------------------------------
	if (heat_iter < 0) {
		corrector = 0;
		iters = 1;
	}
	//--------------------------------------------------------------------------
	// 2. Forcing after implicit viscous terms
	//--------------------------------------------------------------------------
	else {
		corrector = 1;
		iters = heat_iter;
	}

	for (i = 0; i < iters; i++) {

		Communication_update_ghost_nodes_flow_variable(c->data, 'c', params->ghost_nodes, data_bag);
		//Conc_set_boundary_values(c->data, iconc ,CENTRAL_FULL ,grid , params);


		p = p_mobile_list -> start;
		while (p != NULL) {


			Lagrangian_heat_individual(iconc, p, corrector, data_bag);  // changes lag->temp_f if corrector == 0 or direct c->data
			p = p -> next;
		}

		p = p_fixed_list -> start;
		while (p != NULL) {



			Lagrangian_heat_individual(iconc , p, corrector, data_bag);
			p = p -> next;
		}

	}



	T2 = MPI_Wtime();
	data_bag->timer->Wtime_particle_forc += T2 - T1;

return;
}
#endif

#ifdef IBM_SCALAR
void Lagrangian_heat_individual(int iconc, Particle *p, int corrector, Cart3d_bag *data_bag) {



	int i, j, k, mv;
	double r[3];
	double C_d, Q_L;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;
	Concentration *c = data_bag -> c[iconc];


	const double BET[] = {BETA};
	double dtbeta = BET[params->which_stage] * params->dt;
	double idtbeta = 1.0 / dtbeta;
	double *Q_p=p->Q_p;
	double *Q_int=c->Q_int;

	if(c->heating_flag == 0){
		Q_int[params->which_stage]=0; c->heating_flag++;
		Q_p[params->which_stage]=0;
	}



	// Bounds for manipulating Eulerian grid quantities
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	int i_end = grid -> G_Ie;
	int j_end = grid -> G_Je;
	int k_end = grid -> G_Ke;

	// Bounds for determining local Lagrangian points
	int Is = i_start;
	int Js = j_start;
	int Ks = k_start;

	int Ie = min(i_end, grid -> NX - 1);
	int Je = min(j_end, grid -> NY - 1);
	int Ke = min(k_end, grid -> NZ - 1);

	double dV = grid->dx_c[1]*grid->dy_c[1]*grid->dz_c[1];


	double *xc = grid -> xc;
	double *yc = grid -> yc;
	double *zc = grid -> zc;

	double ***c_data = c -> data;
	double ***c_rhs = c -> ng_rhs;
    double ***ng_vfc= lag->ng_vfc;

#ifdef MASKING_OUTSIDE
	double *PVol_L= p->PVol_L;
//	double Vol_L_full = p -> Vol_L_full;
#endif

	//double ***c_ibm = c->ng_ibm;
	double tmp=0, tmp2=0;


	double Vol_L = p -> Vol_L;

	double M_L = Vol_L;


	double R2 = p->R*p->R;

	double *X = p -> X;
	double *F = p -> F;


	double *X_L = p -> X_L;
	double *Y_L = p -> Y_L;
	double *Z_L = p -> Z_L;

	double *X_H = p -> X_H;
	double *Y_H = p -> Y_H;
	double *Z_H = p -> Z_H;

	int *flag_L = p -> flag_L;
#ifdef GRID_UNIFORM
	double h = grid -> dx_u[1];
#else
#error 'Nonuniform grid is not supported for the Lagrangian CONC IBM method'
#endif

	double ***temp_f = lag -> ng_temp;
	double *Temp_U_L = lag -> Temp_L;
	double *Temp_U_H = lag -> Temp_H;



	double A = PARTICLE_BC_A;
	double B = PARTICLE_BC_B;
	double C = PARTICLE_BC_C;


	//--------------------------------------------------------------------------
	// 1. Interpolate conc onto Lagrangian markers
	//--------------------------------------------------------------------------
	Interpolate_Eul_to_Lag(c_data, Temp_U_L, 'c', p, grid);
	Interpolate_Eul_to_Lag(c_data, Temp_U_H, 'H', p, grid);
    //Interpolate_Eul_to_Lag(ng_vfc, PVol_L, 'c', p, grid); // TODO create separate function for reducing overhead

//#ifdef MASKING_OUTSIDE
	Interpolate_Eul_to_Lag(ng_vfc, PVol_L, 'c', p, grid); // TODO create separate function for reducing overhead
//#endif

	//--------------------------------------------------------------------------
	// 2. Calculate Lagrangian marker forces
	//--------------------------------------------------------------------------
	for (mv = 0; mv < p->N_L_local; mv++) {

		// Skip the forcing for flagged Lagrangian markers

		if (flag_L[mv] == 0) {

			Temp_U_L[mv] = 0.0;
			continue;
		}

		// Desired conc of Lagrangian points
		C_d = (C- A*Temp_U_H[mv]/h)/(B-A/h);

		// Heat required to produce desired velocity at Lagrangian point
		Temp_U_L[mv]=  (C_d - Temp_U_L[mv]);


		//----------------------------------------------------------------------
		// 3. Interaction with particle could be added here.
		//----------------------------------------------------------------------


		if (LOCAL_POINT_CONC(X_L[mv], Y_L[mv], Z_L[mv])) {
#ifdef MASKING_OUTSIDE
			Q_p[params->which_stage] += PVol_L[mv]*Temp_U_L[mv]*idtbeta*M_L;
#else
			Q_p[params->which_stage] += Temp_U_L[mv]*idtbeta*M_L;
#endif
		}


#ifdef SCALAR_DEBUG  // only works with a single process
	tmp +=	Temp_U_L[mv];
	tmp2 += Temp_U_L[mv]*Temp_U_L[mv];
#endif

	}
#ifdef SCALAR_DEBUG
	tmp=tmp/p->N_L_local;
	tmp2=sqrt(tmp2/p->N_L_local);
	printf("The mean error between the adjusted temp and the current is %g and the  rms value is %g \n", tmp, tmp2 );
#endif


	//--------------------------------------------------------------------------
	// 4. Spread Lagrangian marker forces onto flow field
	//--------------------------------------------------------------------------


	Interpolate_Lag_to_Eul(Temp_U_L, temp_f, 'c', p, grid);




	if (corrector != 0) {
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
#ifdef MASKING_OUTSIDE
					c_data[k][j][i] += temp_f[k][j][i]*ng_vfc[k][j][i];
#else
					c_data[k][j][i] += temp_f[k][j][i];
#endif

				}
			}
		}
	}
	else {


		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
#ifdef MASKING_OUTSIDE
					c_rhs[k][j][i] +=  temp_f[k][j][i]*ng_vfc[k][j][i]*idtbeta ;
#else
					c_rhs[k][j][i] +=  temp_f[k][j][i]*idtbeta ;
#endif
				}
			}
		}
	}
#ifdef GRID_UNIFORM

tmp=0;
if(params->conc_output_Q_int[iconc]){
	for (k = k_start; k < Ke; k++) {
				for (j = j_start; j < Je; j++) {
					for (i = i_start; i < Ie; i++) {
#ifdef MASKING_OUTSIDE
						tmp += temp_f[k][j][i]*idtbeta*ng_vfc[k][j][i] ;
#else
						tmp += temp_f[k][j][i]*idtbeta ;
#endif
					}
				}
			}
	Q_int[params->which_stage] += tmp*dV;

}
#else
printf("Function __funct__ is ready for a nonuniform grid")
#endif


return;
}
#endif



void Lagrangian_mask_scalar(int iconc, Cart3d_bag *data_bag) {



	int i, j, k;


	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;
	Concentration *c = data_bag -> c[iconc];

	double ***data = c->data;
	double ***ng_vfc =  lag->ng_vfc;

	// Bounds for manipulating Eulerian grid quantities
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	int i_end = grid -> G_Ie;
	int j_end = grid -> G_Je;
	int k_end = grid -> G_Ke;

	// Bounds for determining local Lagrangian points
	int Is = i_start;
	int Js = j_start;
	int Ks = k_start;

	int Ie = min(i_end, grid -> NX - 1);
	int Je = min(j_end, grid -> NY - 1);
	int Ke = min(k_end, grid -> NZ - 1);


		for (k = k_start; k < Ke; k++) {
			for (j = j_start; j < Je; j++) {
				for (i = i_start; i < Ie; i++) {
					if(ng_vfc[k][j][i] > 0.99999 ){
				data[k][j][i]= data[k][j][i]*(1-ng_vfc[k][j][i]);
				}}}}

return;
}

#ifdef VOF_SCALAR

void Lagrangian_calc_compound_vof_fluid_velo(Cart3d_bag *data_bag) {

	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;



#if defined(VOF_SMOOTH_VELO) || defined(VOF_PP_VFPRIME)
	double ***vfu = data_bag->lag->vfu_prime;
	double ***vfv = data_bag->lag->vfv_prime;
	double ***vfw = data_bag->lag->vfw_prime;

#else
	double ***vfu = data_bag->lag->vfu;
	double ***vfv = data_bag->lag->vfv;
	double ***vfw = data_bag->lag->vfw;
#endif

	double ***u_vof= data_bag-> u->data_vof;
	double ***v_vof= data_bag-> v->data_vof;
	double ***w_vof= data_bag-> w->data_vof;


	int i_start = data_bag ->grid -> G_Is;
	int j_start = data_bag ->grid -> G_Js;
	int k_start = data_bag ->grid -> G_Ks;


	int i_end =  data_bag ->grid->G_Ie;
	int j_end =  data_bag ->grid->G_Je;
	int k_end =  data_bag ->grid->G_Ke;

	int i,j,k;


	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
	u_vof[k][j][i]=u_data[k][j][i]*(1-vfu[k][j][i])+u_vof[k][j][i];
	v_vof[k][j][i]=v_data[k][j][i]*(1-vfv[k][j][i])+v_vof[k][j][i];
	w_vof[k][j][i]=w_data[k][j][i]*(1-vfw[k][j][i])+w_vof[k][j][i];
	}}}


	return;
}


/******************************************************************************/
/*
*/
#endif
