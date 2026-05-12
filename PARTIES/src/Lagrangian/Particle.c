#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf5.h"

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Communication.h"
#include "Display.h"
#include "Interpolate.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "Particle.h"
#include "ParticleInput.h"
#include "Velocity.h"

MPI_Datatype MPI_PARTICLE;
MPI_Datatype MPI_COLLISION;

// Effective radius of particle (area of influence)
#define R_EFF (p->R)


/******************************************************************************/
/*
  Checks that 'array' can hold 'newsize' elements of type 'array_type'.  If its
  current size (stored in 'lag->array_size') is not large enough, it will
  reallocate the array to hold 2*'newsize' elements.
 */
/******************************************************************************/
#define REALLOC_BUFFER(array, newsize, array_type) \
	do { if ((newsize) > lag->array ## _size) { \
		free(array); \
		array = (array_type *)malloc( 2 * (newsize) * sizeof(array_type) ); \
		Memory_check_allocation(array); \
		lag->array = array; \
		lag->array ## _size = 2 * (newsize); \
	} } while(0)


/******************************************************************************/
/*
  Effectively does this:

      array[index] = value;

  where 'array' is the name of an array stored in 'lag'.  If 'index' is beyond
  the array's length (which is stored in 'lag->array_size'), 'array' will be
  reallocated to have size '2*index+1'.
 */
/******************************************************************************/
#define ADD_TO_BUFFER(array, index, value, array_type) \
	do { \
		if ((index) >= lag->array ## _size) { \
			int i; \
			array_type *new_array; \
			int new_size = 2 * (index) + 1; \
			new_array = (array_type *)malloc( new_size * sizeof(array_type) ); \
			Memory_check_allocation(new_array); \
			for (i = 0; i < lag->array ## _size; i++) { \
				new_array[i] = array[i]; \
			} \
			free(array); \
			array = new_array; \
			lag->array = array; \
			lag->array ## _size = new_size; \
		} \
		array[index] = (value); \
	} while(0)


/******************************************************************************/
/*
  Fills buffer array 'send_Nc' with particle collisions.  Also calculates and
  populates 'send_Nc'.
 */
/******************************************************************************/
#define FILL_COLLISION_SEND_BUFFER() \
	do { \
		int particle_index = 2 * send_Np + 1; \
		int wall_index     = particle_index + 1; \
		ADD_TO_BUFFER(send_Nc, particle_index, 0, int); \
		send_Nc[wall_index] = 0; \
		pc = p -> particle_collision; \
		while (pc != NULL) { \
			ADD_TO_BUFFER(send_coll, send_Nc[0], *pc, Collision); \
			send_coll[send_Nc[0]].next = NULL; \
			pc = pc -> next; \
			send_Nc[0]++; \
			send_Nc[particle_index]++; \
		} \
		pc = p -> wall_collision; \
		while (pc != NULL) { \
			ADD_TO_BUFFER(send_coll, send_Nc[0], *pc, Collision); \
			send_coll[send_Nc[0]].next = NULL; \
			pc = pc -> next; \
			send_Nc[0]++; \
			send_Nc[wall_index]++; \
		} \
	} while(0)


#define NULLIFY_PARTICLE_PTRS(particle) \
	do { \
		particle.X_L = NULL; \
		particle.Y_L = NULL; \
		particle.Z_L = NULL; \
		particle.Vol_L_marker = NULL; \
		particle.flag_L = NULL; \
		particle.particle_collision = NULL; \
		particle.wall_collision = NULL; \
		particle.next = NULL; \
	} while(0)


/******************************************************************************/
/*
 Macro to be used in Particle_MPI_update only.  Communicates particles and
 collisions to/from the receive/send buffers and checks the MPI statuses.

 Variables this macro depends on (but does not modify):
     - src:  source processor (where receiving from)
     - dest: destination processor (where sending to)
     - send_Np: number of particles to be sent
     - send_part: array of particles to be sent
     - send_Nc: array of number of collisions per particle to be sent
     - send_coll: array of collisions to be sent

 Variables this macro modifies (but does not depend on):
     - recv_Np: number of particles to be received
     - recv_part: array of particles to be received
     - recv_Nc: array of number of collisions per particle to be received
     - recv_coll: array of collisions to be received
 */
/******************************************************************************/
#define SENDRECV_PARTICLES_COLLISIONS() \
	do { \
		recv_Np = 0; \
		MPI_Sendrecv(&send_Np, 1, MPI_INT, dest, tag, \
		             &recv_Np, 1, MPI_INT, src, tag, PCW, &status); \
		\
		REALLOC_BUFFER(recv_part, recv_Np, Particle); \
		\
		MPI_Sendrecv(send_part, send_Np, MPI_PARTICLE, dest, tag, \
		             recv_part, recv_Np, MPI_PARTICLE, src, tag, PCW, &status); \
		\
		REALLOC_BUFFER(recv_Nc, 2*recv_Np+1, int); \
		\
		recv_Nc[0] = 0; \
		MPI_Sendrecv(send_Nc, 2*send_Np+1, MPI_INT, dest, tag, \
		             recv_Nc, 2*recv_Np+1, MPI_INT, src, tag, PCW, &status); \
		\
		REALLOC_BUFFER(recv_coll, recv_Nc[0], Collision); \
		\
		MPI_Sendrecv(send_coll, send_Nc[0], MPI_COLLISION, dest, tag, \
		             recv_coll, recv_Nc[0], MPI_COLLISION, src, tag, PCW, &status); \
	} while(0)


/******************************************************************************/
/*
 * Initializes particle information: reads data from input/resume files and
 * allocates storage for Lagrangian marker temp array.
 */
/******************************************************************************/
void Particle_initialize(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Lagrangian *lag    = data_bag -> lag;

	Particle_list *p_mobile_list = lag -> p_mobile_list;
	Particle_list *p_fixed_list  = lag -> p_fixed_list;

#ifdef PARTICLE_RELEASE
	Particle_list *p_release_list  = lag -> p_release_list;
#endif
	Particle_initialize_MPI_datatype();

	if (!params->resume) {

#ifdef STARTUP
		params -> startup_flag = 1;
#endif

		p_fixed_list->ID_start = 0;
		p_fixed_list->type = FIXED;
		ParticleInput_inp(p_fixed_list, grid, params, DTRACE("ParticleInput_inp"));
		Display_progress(params,"Read "P_FIXED_INPUT_FILE"\n");

		p_mobile_list->ID_start = p_fixed_list->Np;
		p_mobile_list->type = MOBILE;
		ParticleInput_inp(p_mobile_list, grid, params, DTRACE("ParticleInput_inp"));
		Display_progress(params,"Read "P_MOBILE_INPUT_FILE"\n");

#ifdef PARTICLE_RELEASE
		p_release_list->ID_start = p_fixed_list->Np + p_mobile_list->Np;
		p_release_list->type = RELEASE;
		ParticleInput_inp(p_release_list, grid, params, DTRACE("ParticleInput_inp"));
		Display_progress(params,"Read "P_RELEASE_INPUT_FILE"\n");

		p_mobile_list->Np = p_mobile_list->Np + p_release_list->Np;
#endif
		Particle_initialize_velocities(data_bag, DTRACE("Particle_initialize_velocities"));
	}
	else {
		ParticleInput_h5(data_bag, DTRACE("ParticleInput_h5"));
		Particle_initialize_volume_fraction(data_bag, DTRACE("Particle_initialize_volume_fraction"));
	}

	//--------------------------------------------------------------------------
	// Allocate grid storage for Lagrangian markers
	//--------------------------------------------------------------------------
	Particle *p;
	int N_L_max, N_L_max_local;

	// Maximum N_L among local particles
	N_L_max_local = 0;
	p = p_mobile_list -> start;
	while (p != NULL) {
		N_L_max_local = max(N_L_max_local, p -> N_L);
		p = p -> next;
	}
	p = p_fixed_list -> start;
	while (p != NULL) {
		N_L_max_local = max(N_L_max_local, p -> N_L);
		p = p -> next;
	}
#ifdef PARTICLE_RELEASE
	p = p_release_list -> start;
	while (p != NULL) {
		N_L_max_local = max(N_L_max_local, p -> N_L);
		p = p -> next;
	}
#endif

	// Maximum N_L among all particles in this linked list
	MPI_Allreduce(&N_L_max_local, &N_L_max, 1, MPI_INT, MPI_MAX, PCW);
	lag -> Temp_L = Memory_allocate_1D_array(GVG_DOUBLE, N_L_max);
	lag -> Temp_H = Memory_allocate_1D_array(GVG_DOUBLE, N_L_max);
	#ifdef VOF_IBM
	lag -> Temp_L_rho = Memory_allocate_1D_array(GVG_DOUBLE, N_L_max);
	#endif


}



/******************************************************************************/
/*
 * Builds Particle MPI datatype
 */
/******************************************************************************/
void Particle_initialize_MPI_datatype() {

	int          blockcounts[1]; // Number of values in each block
	MPI_Datatype types[1];       // Dataype of block
	MPI_Aint     offsets[1];     // Offset of block from beginning of structure

	offsets[0] = 0;
	types[0] = MPI_BYTE;
	blockcounts[0] = sizeof(Particle);

	// Now define structured type and commit it
	MPI_Type_create_struct(1, blockcounts, offsets, types, &MPI_PARTICLE);
	MPI_Type_commit(&MPI_PARTICLE);

	blockcounts[0] = sizeof(Collision);

	// Now define structured type and commit it
	MPI_Type_create_struct(1, blockcounts, offsets, types, &MPI_COLLISION);
	MPI_Type_commit(&MPI_COLLISION);
}




/******************************************************************************/
/*
 * Set initial particle velocities to inital background flow field for immersed
 * particles, or to explicitly set velocities for dry particles
 *
 * Takes in and returns p_mobile as a local list
 */
/******************************************************************************/
void Particle_initialize_velocities(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i;
	Particle *p;
	Particle_list *p_list_foreign;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

#ifdef DRY_PARTICLES
	// Explicitly set velocities for dry particles
	if (p_mobile_list->Np < 4) {
		p = p_mobile_list -> start;
		while (p != NULL) {

			if (p->ID == 0) {
				p->U[0] = 0.0;
				p->U[1] = 0.0;
			}

			p = p -> next;
		}
	}

#else

	// Integrate fluid velocities on each processor subdomain
	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Interpolate_integrate_momentum(data_bag->u, p_mobile_list, data_bag, DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(data_bag->v, p_mobile_list, data_bag, DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(data_bag->w, p_mobile_list, data_bag, DTRACE("Interpolate_integrate_momentum"));

	// Add volume fractions for fixed particles as well
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));
	Interpolate_add_to_volume_fraction('u', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('v', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('w', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	#ifdef VOF_IBM
	Interpolate_add_to_volume_fraction('c', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	#endif

	Particle_list_remove(p_fixed_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));

	#if defined(LAG_PARTICLE_RESOLVED)
	Particle_reduce_oversized_forces_to_owner(p_mobile_list, data_bag);
	#endif

	// Create linked list of foreign particles
	p_list_foreign = Particle_list_foreign_create(p_mobile_list, data_bag, DTRACE("Particle_list_foreign_create"));

	// Receive and combine portions of particles that were on other subdomains
	Particle_MPI_update(p_list_foreign, data_bag, DTRACE("Particle_MPI_update"));
	Lagrangian_collect_forces(p_mobile_list, p_list_foreign, LAG_COLLECT_HYDRO, params, DTRACE("Lagrangian_collect_forces"));
	Particle_list_destroy(p_list_foreign);

	// Set particle velocities based on volume-averaged fluid velocities
	p = p_mobile_list -> start;
	while (p != NULL) {

	    #ifdef VOF_IBM
			// 1. Calculate equilibrium velocity from fluid (likely ~0.0 if quiescent)
			FORI3 {
				if (p->Int_rho[i] > 1e-12) {
					p->U[i]     = p->Int_U[i] / p->Int_rho[i];
					p->Omega[i] = p->Int_Omega[i] / p->Int_rho[i];
				} else {
					p->U[i]     = 0.0;
					p->Omega[i] = 0.0;
				}
			}
	    #else
	        // Single-phase: use solid density
	        FORI3 p->U[i]     = p->Int_U[i] * p->rho_s / p->M;
	        FORI3 p->Omega[i] = p->Int_Omega[i] * p->rho_s / p->I_p;
	    #endif

		// ---------------------------------------------------------------------
		// NEW FIX: Enforce Startup Velocity BEFORE setting U_old
		// This prevents the "Cold Start Shock" (acceleration from 0 to -1 in 1 step)
		// ---------------------------------------------------------------------
		#ifdef VOF_IBM
		#ifdef STARTUP
			// Check if startup is active and mode is time-based (standard release)
			if (params->startup_flag && params->startup_init == STUP_INIT_TIME) {
				FORI3 p->U[i]     = params->startup_velocity[i];
				FORI3 p->Omega[i] = 0.0; // Assume zero rotation for impact
			}
		#endif
		#endif
		// ---------------------------------------------------------------------

		// Now set U_old. Since p->U is now correctly -1.0 (if startup is on),
		// U_old will also be -1.0.
		// Acceleration = (U - U_old)/dt = 0. No unphysical shock.
		FORI3 p->U_old[i]     = p->U[i];
		FORI3 p->Omega_old[i] = p->Omega[i];

#ifdef ONE_WAY
		FORI3 p->U[i]     = 0.0;
		FORI3 p->Omega[i] = 0.0;

		FORI3 p->U_old[i]     = 0.0;
		FORI3 p->Omega_old[i] = 0.0;
#endif

		/* Int_U_old / Int_Omega_old are deferred until after the rigid-body
		 * paint + re-integration below, so they are stored in the same
		 * post-paint reference frame the first time step will see. */

		p = p -> next;
	}

	/* ====================================================================
	 * Cold-start fix: paint rigid-body translation+rotation onto the
	 * fictitious fluid inside each mobile particle, then re-integrate.
	 *
	 * Without this, an imposed-velocity startup leaves the Eulerian fluid
	 * inside the body at its initial state (typically quiescent) while
	 * p->U has just jumped to startup_velocity.  On the first time step
	 * the IBM forcing must drag M_f * U_p of fictitious-fluid momentum
	 * over a single dt, producing an artificial F_rigid impulse of order
	 * (M_f * U_p) / dt that drives sphere over-penetration.
	 *
	 * Painting the interior to U_p + Omega x r and snapshotting Int_U_old
	 * from a fresh integration of that field puts every particle in a
	 * self-consistent state: the first-step IBM correction is an O(dt)
	 * adjustment, not an O(1/dt) shock.  The whole sequence is
	 * multi-rank-safe — each rank paints only the cells it owns, and the
	 * usual foreign-particle exchange handles particles straddling
	 * sub-domain boundaries.
	 * ==================================================================== */

	/* Re-broadcast the updated p->U to ranks that hold neighbouring slices
	 * of each particle (Lagrangian_collect_forces left p_mobile_list in
	 * the LOCAL-only state). */
	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));

	Interpolate_paint_rigid_body_velocity(data_bag->u, p_mobile_list, data_bag,
			DTRACE("Interpolate_paint_rigid_body_velocity"));
	Interpolate_paint_rigid_body_velocity(data_bag->v, p_mobile_list, data_bag,
			DTRACE("Interpolate_paint_rigid_body_velocity"));
	Interpolate_paint_rigid_body_velocity(data_bag->w, p_mobile_list, data_bag,
			DTRACE("Interpolate_paint_rigid_body_velocity"));

	/* Synchronise halos and re-impose physical wall BCs after painting. */
	Communication_update_ghost_nodes_flow_variable(data_bag->u->data, 'u',
			params->ghost_nodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(data_bag->v->data, 'v',
			params->ghost_nodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(data_bag->w->data, 'w',
			params->ghost_nodes, data_bag);
	Velocity_update_boundaries(data_bag->u->data, 'u', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(data_bag->v->data, 'v', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(data_bag->w->data, 'w', VEL_TYPE_NORMAL, data_bag);

	/* Reset per-particle integrators on every copy (owner + foreign) so the
	 * upcoming integration starts clean. */
	{
		Particle *q = p_mobile_list->start;
		while (q != NULL) {
			DSET_ZERO(q->Int_U, 3);
			DSET_ZERO(q->Int_Omega, 3);
#ifdef VOF_IBM
			DSET_ZERO(q->Int_rho, 3);
			q->Int_rho_scalar = 0.0;
#endif
			q = q->next;
		}
	}

	/* Reset staggered-face volume fractions (mobile + fixed contributions
	 * will be re-added below).  ng_vfc is left intact: it is filled only by
	 * fixed particles and is consumed read-only by Int_rho_scalar. */
	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfu);
	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfv);
	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfw);

	Interpolate_integrate_momentum(data_bag->u, p_mobile_list, data_bag,
			DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(data_bag->v, p_mobile_list, data_bag,
			DTRACE("Interpolate_integrate_momentum"));
	Interpolate_integrate_momentum(data_bag->w, p_mobile_list, data_bag,
			DTRACE("Interpolate_integrate_momentum"));

	/* Re-add fixed-particle staggered VF (we cleared ng_vfu/v/w above). */
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));
	Interpolate_add_to_volume_fraction('u', p_fixed_list, data_bag,
			DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('v', p_fixed_list, data_bag,
			DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('w', p_fixed_list, data_bag,
			DTRACE("Interpolate_add_to_volume_fraction"));
	Particle_list_remove(p_fixed_list, FOREIGN, grid, params,
			DTRACE("Particle_list_remove"));

	#if defined(LAG_PARTICLE_RESOLVED)
	Particle_reduce_oversized_forces_to_owner(p_mobile_list, data_bag);
	#endif

	/* Gather post-paint integrals onto owner copies. */
	p_list_foreign = Particle_list_foreign_create(p_mobile_list, data_bag,
			DTRACE("Particle_list_foreign_create"));
	Particle_MPI_update(p_list_foreign, data_bag, DTRACE("Particle_MPI_update"));
	Lagrangian_collect_forces(p_mobile_list, p_list_foreign, LAG_COLLECT_HYDRO,
			params, DTRACE("Lagrangian_collect_forces"));
	Particle_list_destroy(p_list_foreign);

	/* Now snapshot Int_U_old from the self-consistent post-paint integrals.
	 * The first time step's rigid-body correction (Int_U - Int_U_old)/dt
	 * therefore measures the *change* in fictitious-fluid momentum, not the
	 * full impulse needed to spin it up from rest. */
	p = p_mobile_list->start;
	while (p != NULL) {
		FORI3 p->Int_U_old[i]     = p->Int_U[i];
		FORI3 p->Int_Omega_old[i] = p->Int_Omega[i];
		p = p->next;
	}
#endif
}


/******************************************************************************/
/*
 * Initializes volume fraction fields.
 */
/******************************************************************************/
void Particle_initialize_volume_fraction(Cart3d_bag *data_bag, Debug_trace *dtrace)
{
	MAC_grid *grid = data_bag->grid;
	Parameters *params = data_bag->params;
	Lagrangian *lag = data_bag->lag;

	Particle_list *p_mobile_list = lag->p_mobile_list;
	Particle_list *p_fixed_list = lag->p_fixed_list;

	Particle_MPI_update(p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));

	Memory_reset_noghost_variable(grid, params, lag->ng_vfu);
	Memory_reset_noghost_variable(grid, params, lag->ng_vfv);
	Memory_reset_noghost_variable(grid, params, lag->ng_vfw);

	Interpolate_add_to_volume_fraction('u', p_mobile_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('v', p_mobile_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('w', p_mobile_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('u', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('v', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('w', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));

	#ifdef VOF_IBM
	Interpolate_add_to_volume_fraction('c', p_fixed_list, data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	#endif

	Particle_list_remove(p_mobile_list, FOREIGN, grid, params, DTRACE("Particle_MPI_update"));
	Particle_list_remove(p_fixed_list, FOREIGN, grid, params, DTRACE("Particle_MPI_update"));
}


/******************************************************************************/
/*
 * Initialize other parts of particle data structure that will be overwritten
 * before use.  This is done to prevent some arithmetic operations that may
 * reference uninitialized numbers (some of these will be multiplied by zero).
 */
/******************************************************************************/
void Particle_initialize_nonessential_data(Particle *p) {

	memcpy(p->X_old, p->X, 3 * sizeof(double));
	memcpy(p->U_old, p->U, 3 * sizeof(double));
	memcpy(p->Omega_old, p->Omega, 3 * sizeof(double));

	memcpy(p->Fc_old, p->Fc, 3 * sizeof(double));
	memcpy(p->Tc_old, p->Tc, 3 * sizeof(double));

	DSET_ZERO(p->F, 3);
	DSET_ZERO(p->T, 3);

	DSET_ZERO(p->Int_U, 3);
	DSET_ZERO(p->Int_Omega, 3);

	DSET_ZERO(p->F_IBM, 3);
	DSET_ZERO(p->T_IBM, 3);
	DSET_ZERO(p->F_rigid, 3);
	DSET_ZERO(p->T_rigid, 3);
	DSET_ZERO(p->F_coll, 3);
	DSET_ZERO(p->T_coll, 3);

		#ifdef VOF_IBM
			DSET_ZERO(p->F_CCF, 3);
			DSET_ZERO(p->T_CCF, 3);
			DSET_ZERO(p->F_CCF_cum, 3);
			DSET_ZERO(p->T_CCF_cum, 3);
			DSET_ZERO(p->F_body_solid_cum, 3);
		#endif


#ifdef POST_PROCESS
	DSET_ZERO(p->Fc_norm, 3);
	DSET_ZERO(p->Fc_tan, 3);
	DSET_ZERO(p->Fl_norm, 3);
	DSET_ZERO(p->Fl_tan, 3);
	DSET_ZERO(p->Fc_norm_old, 3);
	DSET_ZERO(p->Fc_tan_old, 3);
	DSET_ZERO(p->Fl_norm_old, 3);
	DSET_ZERO(p->Fl_tan_old, 3);
	DSET_ZERO(p->Fc_norm_cum, 3);
	DSET_ZERO(p->Fc_tan_cum, 3);
	DSET_ZERO(p->Fl_norm_cum, 3);
	DSET_ZERO(p->Fl_tan_cum, 3);
#endif


#ifdef POST_PROCESS
	DSET_ZERO(p->Fc_norm, 3);
	DSET_ZERO(p->Fc_tan, 3);
	DSET_ZERO(p->Fl_norm, 3);
	DSET_ZERO(p->Fl_tan, 3);
	DSET_ZERO(p->Fc_norm_old, 3);
	DSET_ZERO(p->Fc_tan_old, 3);
	DSET_ZERO(p->Fl_norm_old, 3);
	DSET_ZERO(p->Fl_tan_old, 3);
	DSET_ZERO(p->Fc_norm_cum, 3);
	DSET_ZERO(p->Fc_tan_cum, 3);
	DSET_ZERO(p->Fl_norm_cum, 3);
	DSET_ZERO(p->Fl_tan_cum, 3);
#endif

}




/******************************************************************************/
/*
 Calculate particle data that can be "derived" from data read from input file
 (e.g. mass and moment of inertia from radius and density)

 Assumptions:
     - uniform grid
 */
/******************************************************************************/
void Particle_calc_derived_data(Particle *p, MAC_grid *grid, Parameters *params) {

	// External data
	double R = p -> R;
	double rho_s = params -> rho_s;
	double rho_prImp = params -> rho_prImp;// density of fixed particles, used only for pressure imposed Couette flow

	// Sphere derived data
	double R2  = R * R;
	double M   = rho_s * 4.0 / 3.0 * PI * R2 * R ;
	double M1  = rho_prImp * 4.0 / 3.0 * PI * R2 * R ; // Mass of fixed particles, used only for pressure imposed Couette flow
	double I_p = rho_s * 8.0 / 15.0 * PI * R2 * R2 * R;

	// Number of Lagrangian marker points determined by grid-to-particle
	// resolution

	double h = grid -> dy_min;

	int N_L = (int)ceil(PI / 3.0 * (12.0 * R2 / (h * h) + 1.0));

	// Volume controlled by each Lagrangian marker point
	double Vol_L = PI * h / (3.0 * N_L) * (12.0 * R2 + h * h);

#if defined(TWOD_CARTESIAN) && defined(LAG_PARTICLE_RESOLVED)
	/*
	 * Phase-3 planar IBM treats each particle as a cylinder extruded through
	 * the storage slab.  Mass and inertia are therefore per slab thickness,
	 * while Vol_L is the regularized line-marker control volume ds*H*h.
	 */
	double slab = (grid->dummy_z_slab_thickness > 0.0) ?
	              grid->dummy_z_slab_thickness : params->twod_slab_thickness;
	if (slab <= 0.0)
		slab = h;
	N_L = max(8, (int)ceil(2.0 * PI * R / h));
	M   = rho_s    * PI * R2 * slab;
	M1  = rho_prImp * PI * R2 * slab;
	I_p = rho_s    * 0.5 * PI * R2 * R2 * slab;
	Vol_L = (2.0 * PI * R / N_L) * slab * h;
#elif defined(AXISYM_RZ) && defined(LAG_PARTICLE_RESOLVED)
	/*
	 * Axisymmetric IBM still represents a physical sphere, but the markers are
	 * meridional rings.  The per-ring control volume varies with radius and is
	 * filled during marker generation; Vol_L remains a safe average fallback.
	 * Marker arc spacing follows Liu et al. (2017) with ds ~ 1.2 h to keep the
	 * IBM correction matrix well-conditioned.
	 */
	N_L = max(4, (int)ceil(PI * R / (1.2 * h)));
	Vol_L = 4.0 * PI * R2 * h / N_L;
#endif

	p -> rho_s = rho_s;
	p -> rho_prImp = rho_prImp;
	p -> M   = M;
	p -> M1  = M1;
	p -> I_p = I_p;
	p -> N_L   = N_L;
	p -> Vol_L = Vol_L;
	p -> N_L_local = 0;

#ifndef GRID_UNIFORM
	p -> N_l_subsec_max = (int)ceil(2*R/h);
#endif

	p -> particle_collision = NULL;
	p -> wall_collision = NULL;
}



/******************************************************************************/
/*
 Frees all the memory associated with a linked particle.
 */
/******************************************************************************/
void Particle_destroy(Particle *p) {
	Particle_collision_list_destroy(p);
	Particle_destroy_internal_arrays(p);
	free(p);
}

/******************************************************************************/
/*
 Frees all the memory associated with a linked list.
 */
/******************************************************************************/
void Particle_list_destroy(Particle_list *p_list) {

	// Particle parameters
	Particle *p, *p_next;

	p = p_list->start;
	while (p != NULL) {
		p_next = p -> next;
		Particle_destroy(p);
		p = p_next;
	}

	free(p_list);

}

void Particle_mobile_turn_off(Particle_list *p_list, Parameters *params) {

	// Particle parameters
	Particle *p, *p_next;

	p = p_list->start;
	while (p != NULL) {
		p_next = p -> next;
		if (p->X[1] > params->y_cut){
			printf("You ! You gon die\n");
			//Particle_destroy(p);
			p->X[1] = params->ymax + 2 * p->R;
			printf("You dead ?\n");
		}
		p = p_next;
	}
}



/******************************************************************************/
/*
 Splits inputted list into local and foreign lists

 Input:
     - 'p_list' containing both local and foreign particles
 Output:
     - 'p_list' containing local particles
	 - returned list containing foreign particles
 */
/******************************************************************************/
Particle_list *Particle_list_foreign_create(Particle_list *p_list, Cart3d_bag *data_bag,
		Debug_trace *dtrace) {

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	// Check state of linked lists
	Display_assert_list_state(p_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	// Create linked list of foreign particles
	Particle_list *p_list_foreign = (Particle_list *)malloc(sizeof(Particle_list));
	// Copy list information
	*p_list_foreign = *p_list;
	p_list_foreign->start = NULL;
//	Particle *p_foreign = NULL;
//	char p_foreign_state = *list_state_ptr;
	// Copy linked list
	Particle_list_copy(p_list, p_list_foreign);
	Particle_list_remove(p_list_foreign, LOCAL, grid, params, DTRACE("Particle_list_remove"));

	// Create linked list of local particles
	Particle_list_remove(p_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));

	return p_list_foreign;
}



/******************************************************************************/
/*
 Copies contents of Particle linked list into new Particle linked list which has
 not yet been allocated.

 Note:
     - Currently allocates space for, but does not copy the contents of,
       internal arrays.

 Inputs:
     p_list_src - source linked list
     p_list_new - new linked list
 */
/******************************************************************************/
void Particle_list_copy(Particle_list *p_list_src, Particle_list *p_list_new) {

	Particle *p, *p_last, *p_src;
	p_list_new->start = NULL;

	p_src = p_list_src -> start;
	while (p_src != NULL) {

		p = (Particle *)malloc(sizeof(Particle));
		Memory_check_allocation(p);

		if (p_list_new->start == NULL)
			p_list_new->start = p;
		else
			p_last -> next = p;

		*p = *p_src;
		Particle_create_internal_arrays(p);
		Particle_collision_list_copy(p_src, p);

		p_src = p_src -> next;
		p_last = p;
	}

	if (p_list_new->start != NULL)
		p -> next = NULL;

}

void Particle_release_to_mobile(Particle_list *p_list_release, Particle_list *p_list_mobile, Parameters *params, Debug_trace *dtrace) {

	Particle *p, *p_last;
	int i;

	Display_assert_list_state(p_list_release, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_list_mobile, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	p = p_list_release->start;
	p_last = p;

	while (p != NULL) {

		if (p->t_part_release < params->time){
			if (p==p_list_release->start){
				p_list_release->start = p->next;
				p->next = p_list_mobile->start;
				p_list_mobile->start = p;
				p = p_list_release->start;
				p_last = p;
			}
			else{
				p_last->next = p->next;
				p->next = p_list_mobile->start;
				p_list_mobile->start = p;
				p = p_last->next;
			}
			// p_list_mobile->start->U[0] = 1;

		}
		else{
			p_last = p;
			p = p->next;
		}
	}

}

/******************************************************************************/
/*
 Adds array of particles to beginning of linked list.

 Inputs:
     p_list  - particle linked list structure
     p_array - array of particles to add
     size    - length of 'p_array'

 Output:
     Pointer to last element of newly enlarged linked list
 */
/******************************************************************************/
void Particle_list_add_array(Particle_list *p_list, Particle *p_array, int p_size,
		Collision *pc_array, int *Nc, MAC_grid *grid) {

	int i, pc_size;
	Particle *p, *p_start_new;

	if (p_size > 0) {

		// Allocate storage
		p_start_new = (Particle *)malloc(sizeof(Particle));
		Memory_check_allocation(p_start_new);
		p = p_start_new;

		// Copy particle data
		*p = p_array[0];
		Particle_create_internal_arrays(p);

		// Copy particle collision list
		pc_size = Nc[1];
		Particle_collision_list_add_array(&(p->particle_collision), pc_array, pc_size);
		pc_array += pc_size;

		// Copy wall collision list
		pc_size = Nc[2];
		Particle_collision_list_add_array(&(p->wall_collision), pc_array, pc_size);
		pc_array += pc_size;

		for (i = 1; i < p_size; i++) {

			// Allocate storage
			p -> next = (Particle *)malloc(sizeof(Particle));
			Memory_check_allocation(p -> next);
			p = p -> next;

			// Copy particle data
			*p = p_array[i];
			Particle_create_internal_arrays(p);

			// Copy particle collision list
			pc_size = Nc[2 * i + 1];
			Particle_collision_list_add_array(&(p->particle_collision), pc_array, pc_size);
			pc_array += pc_size;

			// Copy wall collision list
			pc_size = Nc[2 * i + 2];
			Particle_collision_list_add_array(&(p->wall_collision), pc_array, pc_size);
			pc_array += pc_size;
		}

		p -> next = p_list->start;

		p_list->start = p_start_new;
	}

}




/******************************************************************************/
/*
 Removes particles residing either on this or other processor subdomains and
 returns the number of particles remaining in the linked list.

 Inputs:
     p_list - particle linked list structure to remove particles from
     rm_type - character indicating what to remove from the linked list:
         LOCAL - remove local particles
         FOREIGN - remove foreign particles
	 line - use '__LINE__' to get line number of function call
	 file - use '__FILE__' to get filename of function call

 Assumptions:
 - Particles do not move more than one radius per time step
 */
/******************************************************************************/
void Particle_list_remove(Particle_list *p_list, int rm_type, MAC_grid *grid,
		Parameters *params, Debug_trace *dtrace) {

	// Grid bounds on local subdomain
	double G_xmin = grid -> xu[grid -> G_Is];
	double G_xmax = grid -> xu[min(grid -> G_Ie, grid -> NX - 1)];
	double G_ymin = grid -> yv[grid -> G_Js];
	double G_ymax = grid -> yv[min(grid -> G_Je, grid -> NY - 1)];
	double G_zmin = grid -> zw[grid -> G_Ks];
	double G_zmax = grid -> zw[min(grid -> G_Ke, grid -> NZ - 1)];

	// Particle parameters
	double *X, R;
	Particle *p, *p_last, *p_next, *p_start;

	// Boolean variables
	int outside;     // Is particle outside subdomain?
	int rm_foreign;  // Remove foreign particles?

	char message[100];

	//--------------------------------------------------------------------------
	// Read removal type
	//--------------------------------------------------------------------------
	if (rm_type == LOCAL)
		rm_foreign = 0;
	else if (rm_type == FOREIGN)
		rm_foreign = 1;
	else {
		sprintf(message, "Incorrect input rm_type = '%d'\n"
				"Use 'LOCAL' or 'FOREIGN'", rm_type);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	//--------------------------------------------------------------------------
	// Check current linked list state and assign its new designation
	//--------------------------------------------------------------------------
	// Initial list state should contain both foreign and local particles
	if (p_list->state == LIST_STATE_BOTH) {
		if (rm_foreign)
			p_list->state = LIST_STATE_LOCAL;
		else
			p_list->state = LIST_STATE_FOREIGN;
	}
	// Any other initial state will produce a useless set of particles
	else {
		Display_assert_list_state(p_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));
	}

	//--------------------------------------------------------------------------
	// Remove particles of correct type
	//--------------------------------------------------------------------------
	p_start = p_list -> start;
	p = p_start;
	p_last = NULL;
	while (p != NULL) {

		X = p -> X;
		p_next = p -> next;

		// Debugging
		//		printf("X[0] = %16.14g, G_xmin = %16.14g, G_xmax = %16.14g, X[0] < G_xmin = %d, X[0] >= G_xmax = %d\n",
		//			   X[0], G_xmin, G_xmax, X[0] < G_xmin, X[0] >= G_xmax);

		outside = (X[0] < G_xmin || X[0] >= G_xmax || X[1] < G_ymin || X[1] >= G_ymax ||
				   X[2] < G_zmin || X[2] >= G_zmax);
#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
		/*
		 * In 2D IBM the third coordinate is a storage slab, not an ownership
		 * direction.  Keeping removal in-plane makes the particle exchange
		 * rank-independent when NPZ is forced to one.
		 */
		outside = (X[0] < G_xmin || X[0] >= G_xmax ||
		           X[1] < G_ymin || X[1] >= G_ymax);
#endif
		if (outside && rm_foreign || !outside && !rm_foreign) {

			if (p_last != NULL)
				p_last -> next = p_next;

			Particle_collision_list_destroy(p);
			Particle_destroy_internal_arrays(p);
			if (p == p_start)
				p_start = p_next;

			free(p);
			p = p_next;
		}
		else {
			p_last = p;
			p = p_next;
		}
	}

	p_list -> start = p_start;
}




/******************************************************************************/
/*
 Frees all the memory associated with a particle's Collision linked lists

 Inputs:
     p - pointer to particle containing Collision linked lists
 */
/******************************************************************************/
void Particle_collision_list_destroy(Particle *p) {

	Collision *pc, *pc_next;

	pc = p -> particle_collision;
	while (pc != NULL) {
		pc_next = pc -> next;
		free(pc);
		pc = pc_next;
	}
	pc = p -> wall_collision;
	while (pc != NULL) {
		pc_next = pc -> next;
		free(pc);
		pc = pc_next;
	}

	p -> particle_collision = NULL;
	p -> wall_collision = NULL;

}




/******************************************************************************/
/*
 Copy entirety of linked lists 'particle_collision' and 'wall_collision' from
 'p_src' to 'p_new'.

 Note:
     - Collisions will be copied over in a reverse order linked list
     - 'p_new' should contain no Collision linked lists when this function is
       called
 */
/******************************************************************************/
void Particle_collision_list_copy(Particle *p_src, Particle *p_new) {

	Collision *pc, *pc_src;

	//--------------------------------------------------------------------------
	// Particle collisions
	//--------------------------------------------------------------------------
	pc_src = p_src -> particle_collision;
	p_new -> particle_collision = NULL;
	while (pc_src != NULL) {

		pc = (Collision *)malloc(sizeof(Collision));
		Memory_check_allocation(pc);

		*pc = *pc_src;

		pc -> next = p_new -> particle_collision;
		p_new -> particle_collision = pc;

		pc_src = pc_src -> next;
	}

	//--------------------------------------------------------------------------
	// Wall collisions
	//--------------------------------------------------------------------------
	pc_src = p_src -> wall_collision;
	p_new -> wall_collision = NULL;
	while (pc_src != NULL) {

		pc = (Collision *)malloc(sizeof(Collision));
		Memory_check_allocation(pc);

		*pc = *pc_src;

		pc -> next = p_new -> wall_collision;
		p_new -> wall_collision = pc;

		pc_src = pc_src -> next;
	}

}




/******************************************************************************/
/*
 Transfer collision linked lists from 'p_src' to 'p_new'.  Collisions will be
 transferred if they do not exist in 'p_new' or if 'p_src' is the true collision
 owner, indicated by state 'COLL_STATE_OWNER'.

 Inputs:
     type - character describing which lists to transfer:
         'p' - transfer only particle collisions
         'w' - transfer only wall collisions
         'b' - transfer both particle and wall collisions

 Returns:
     status - value of 0 if transfer is valid
            - value of -1 if two conflicting collision owners are encountered.
              This would happen if there are multiple simultaneous collisions
              between two particles, which might occur if a periodic domain is
              one particle wide.
            - value of -2 if conflicting collision owners for a wall collision.
              This should only happen if something in Collision.c is messed up.

 Note:
     - Transferred collisions will be added to the beginning of the Collision
       lists
     - 'p_new' should contain Collision lists, which should be 'NULL' if they
        are empty
 */
/******************************************************************************/
int Particle_collision_list_xfer(Particle *p_src, Particle *p_new, char type) {

	int status = 0;
	Collision *pc, *pc_src, *pc_next;

	//--------------------------------------------------------------------------
	// Particle collisions
	//--------------------------------------------------------------------------
	pc_src = p_src -> particle_collision;
	while (pc_src != NULL && (type == 'b' || type == 'p')) {

		// See if collision in 'p_src' is contained in 'p_new'
		pc = p_new -> particle_collision;
		while (pc != NULL) {
			if (pc_src->other_ID == pc->other_ID)
				break;
			pc = pc -> next;
		}

		// If collision does not exist in 'p_new', add it to beginning of list
		if (pc == NULL) {
			pc = (Collision *)malloc(sizeof(Collision));
			Memory_check_allocation(pc);

			*pc = *pc_src;

			pc -> next = p_new -> particle_collision;
			p_new -> particle_collision = pc;
		}
		// Throw error for multiple simultaneous collisions
		else if (pc_src->state == COLL_STATE_OWNER && pc->state == COLL_STATE_OWNER) {
			status = -1;
		}
		// Copy contents over from collision owner
		else if (pc_src->state == COLL_STATE_OWNER) {
			pc_next = pc->next;
			*pc = *pc_src;
			pc->next = pc_next;
		}
		// Copy destroy command unless owned by another collision
		else if (pc_src->state == COLL_STATE_DESTROY && pc->state != COLL_STATE_OWNER) {
			pc->state = COLL_STATE_DESTROY;
		}

		pc_src = pc_src -> next;
	}

	//--------------------------------------------------------------------------
	// Wall collisions
	//--------------------------------------------------------------------------
	pc_src = p_src -> wall_collision;
	while (pc_src != NULL && (type == 'b' || type == 'w')) {

		// See if collision in 'p_src' is contained in 'p_new'
		pc = p_new -> wall_collision;
		while (pc != NULL) {
			if (pc_src->other_ID == pc->other_ID)
				break;
			pc = pc -> next;
		}

		// If collision does not exist in 'p_new', add it to beginning of list
		if (pc == NULL) {
			pc = (Collision *)malloc(sizeof(Collision));
			Memory_check_allocation(pc);

			*pc = *pc_src;

			pc -> next = p_new -> wall_collision;
			p_new -> wall_collision = pc;
		}
		// Throw error for multiple simultaneous collisions
		else if (pc_src->state == COLL_STATE_OWNER && pc->state == COLL_STATE_OWNER) {
			status = -2;
		}
		// Copy contents over from collision owner
		else if (pc_src->state == COLL_STATE_OWNER) {
			pc_next = pc->next;
			*pc = *pc_src;
			pc->next = pc_next;
		}
		// Copy destroy command unless owned by another collision
		else if (pc_src->state == COLL_STATE_DESTROY && pc->state != COLL_STATE_OWNER) {
			pc->state = COLL_STATE_DESTROY;
		}

		pc_src = pc_src -> next;
	}

	return status;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Particle_collision_list_add_array(Collision **pc_start_ptr,
		Collision *pc_array, int pc_size) {

	int i;
	Collision *pc;
//	Particle *p, *p_start_new;

	if (pc_size > 0) {
		pc = (Collision *)malloc(sizeof(Collision));
		Memory_check_allocation(pc);
		*pc_start_ptr = pc;
		*pc = pc_array[0];

		for (i = 1; i < pc_size; i++) {
			pc -> next = (Collision *)malloc(sizeof(Collision));
			Memory_check_allocation(pc -> next);
			pc = pc -> next;
			*pc = pc_array[i];
		}

		pc -> next = NULL;
	}
	else {
		*pc_start_ptr = NULL;
	}

}




/******************************************************************************/
/*
 * Looks through collision lists of all particles in list 'p_start', removing
 * any collisions that have the state 'COLL_STATE_DESTROY'.
 */
/******************************************************************************/
#ifdef STORE_COLLISION
void Particle_collision_list_clean(Particle *p_start) {

	Particle *p;
	Collision *pc, *pc_last;

	p = p_start;
	while (p != NULL) {

		// Particle collisions
		pc = p -> particle_collision;
		while (pc != NULL) {
			if (pc->state == COLL_STATE_DESTROY) {
				if (pc == p->particle_collision) {
					p->particle_collision = pc -> next;
					free(pc);
					pc = p->particle_collision;
				}
				else {
					pc_last -> next = pc -> next;
					free(pc);
					pc = pc_last -> next;
				}
			}
			else {
				pc_last = pc;
				pc = pc -> next;
			}
		}

		// Wall collisions
		pc = p -> wall_collision;
		while (pc != NULL) {
			if (pc->state == COLL_STATE_DESTROY) {
				if (pc == p->wall_collision) {
					p->wall_collision = pc -> next;
					free(pc);
					pc = p->wall_collision;
				}
				else {
					pc_last -> next = pc -> next;
					free(pc);
					pc = pc_last -> next;
				}
			}
			else {
				pc_last = pc;
				pc = pc -> next;
			}
		}

		p = p -> next;
	}
}
#endif




/******************************************************************************/
/*
  Allocate memory for pointer arrays within Particle structure
 */
/******************************************************************************/
void Particle_create_internal_arrays(Particle *p) {

	int N_L = p -> N_L;

	p -> X_L = (double *)malloc(3 * N_L * sizeof(double));
	p -> Y_L = &(p->X_L[N_L]);
	p -> Z_L = &(p->X_L[2 * N_L]);
	p -> Vol_L_marker = (double *)malloc(N_L * sizeof(double));
	Memory_check_allocation(p->Vol_L_marker);
	for (int i = 0; i < N_L; i++)
		p->Vol_L_marker[i] = p->Vol_L;
#ifdef IBM_SCALAR
	p -> X_H = (double *)malloc(3 * N_L * sizeof(double));
	p -> Y_H = &(p->X_H[N_L]);
	p -> Z_H = &(p->X_H[2 * N_L]);
#endif


	p -> flag_L = (int *)malloc(N_L * sizeof(int));
//	p -> U_mv = (double *)malloc(3 * N_mv * sizeof(double));
//	p -> V_mv = &(p->U_mv[N_mv]);
//	p -> W_mv = &(p->U_mv[2 * N_mv]);

#ifndef GRID_UNIFORM
	int N_l_subsec_max = p -> N_l_subsec_max;

	p -> idx_x_markers = (int*) calloc(N_L, sizeof(int));
	p -> idx_y_markers = (int*) calloc(N_L, sizeof(int));
	p -> idx_z_markers = (int*) calloc(N_L, sizeof(int));

	p -> Vl = (double*) calloc(N_L, sizeof(double));

	p -> N_l_subsec = (double*) calloc(N_l_subsec_max, sizeof(double));
	p -> idx_subsec_low = (int*) calloc(N_l_subsec_max, sizeof(double));

	p -> A = Memory_allocate_2D_double_array((int)(1.1*p->N_L), 9);
	p -> B = Memory_allocate_2D_double_array(9, (int)(1.1*p->N_L));
	printf("===== N_L_allocated = %d\n", p->N_L);  // TODO: DELETE

#endif
}




/******************************************************************************/
/*
  Free memory for pointer arrays within Particle structure
 */
/******************************************************************************/
void Particle_destroy_internal_arrays(Particle *p) {

	free(p -> X_L);
	free(p -> Vol_L_marker);
#ifdef IBM_SCLAR
	free(p -> X_H);
#endif
	free(p -> flag_L);
//	free(p -> U_mv);



#ifndef GRID_UNIFORM
	free(p -> idx_x_markers);
	free(p -> idx_y_markers);
	free(p -> idx_z_markers);

	free(p -> Vl);

	free(p -> N_l_subsec);
	free(p -> idx_subsec_low);

	free(p -> A);
	free(p -> B);
#endif
}


#if defined(LAG_PARTICLE_RESOLVED)
#ifdef VOF_IBM
#define PARTICLE_FORCE_RECORD_N 28
#else
#define PARTICLE_FORCE_RECORD_N 12
#endif

typedef struct {
	int ID;
	double value[PARTICLE_FORCE_RECORD_N];
} Particle_force_record;

static void Particle_rank_bounds(MAC_grid *grid, Parameters *params,
                                 int rank, double lower[3], double upper[3],
                                 int coords[3])
{
	MPI_Cart_coords(PCW, rank, 3, coords);

	const int Is = (float)(coords[0]    ) / (float)params->NPX * grid->NX;
	const int Ie = (float)(coords[0] + 1) / (float)params->NPX * grid->NX;
	const int Js = (float)(coords[1]    ) / (float)params->NPY * grid->NY;
	const int Je = (float)(coords[1] + 1) / (float)params->NPY * grid->NY;
	const int Ks = (float)(coords[2]    ) / (float)params->NPZ * grid->NZ;
	const int Ke = (float)(coords[2] + 1) / (float)params->NPZ * grid->NZ;

	lower[0] = grid->xu[Is];
	upper[0] = grid->xu[min(Ie, grid->NX - 1)];
	lower[1] = grid->yv[Js];
	upper[1] = grid->yv[min(Je, grid->NY - 1)];
	lower[2] = grid->zw[Ks];
	upper[2] = grid->zw[min(Ke, grid->NZ - 1)];
}

static int Particle_coordinate_is_owned(double x, double lower, double upper,
                                        int is_last_rank)
{
	return (x >= lower && (x < upper || (is_last_rank && x <= upper)));
}

int Particle_center_is_owned_by_rank(const Particle *p, MAC_grid *grid,
                                     Parameters *params, int rank)
{
	double lower[3], upper[3];
	int coords[3];
	Particle_rank_bounds(grid, params, rank, lower, upper, coords);

	if (!Particle_coordinate_is_owned(p->X[0], lower[0], upper[0],
	                                  coords[0] == params->NPX - 1))
		return 0;
	if (!Particle_coordinate_is_owned(p->X[1], lower[1], upper[1],
	                                  coords[1] == params->NPY - 1))
		return 0;
#if defined(TWOD_MODE)
	return 1;
#else
	return Particle_coordinate_is_owned(p->X[2], lower[2], upper[2],
	                                    coords[2] == params->NPZ - 1);
#endif
}

static int Particle_clamped_proc_coord(double x, double xmin, double xmax, int np)
{
	double xi = 0.0;
	if (xmax > xmin)
		xi = (x - xmin) / (xmax - xmin);

	int coord = (int)floor(xi * (double)np);
	if (coord < 0) coord = 0;
	if (coord >= np) coord = np - 1;
	return coord;
}

static double Particle_wrap_periodic_coordinate(double x, double xmin,
                                                double xmax)
{
	const double length = xmax - xmin;
	if (length <= 0.0)
		return x;

	while (x < xmin) x += length;
	while (x >= xmax) x -= length;
	return x;
}

int Particle_center_owner_rank(const Particle *p, MAC_grid *grid,
                               Parameters *params)
{
	int nproc;
	MPI_Comm_size(PCW, &nproc);

	Particle q = *p;
#ifdef XPERIODIC
	q.X[0] = Particle_wrap_periodic_coordinate(q.X[0], params->xmin,
	                                           params->xmax);
#endif
#ifdef YPERIODIC
	q.X[1] = Particle_wrap_periodic_coordinate(q.X[1], params->ymin,
	                                           params->ymax);
#endif
#if !defined(TWOD_MODE)
#ifdef ZPERIODIC
	q.X[2] = Particle_wrap_periodic_coordinate(q.X[2], params->zmin,
	                                           params->zmax);
#endif
#endif

	for (int rank = 0; rank < nproc; rank++)
		if (Particle_center_is_owned_by_rank(&q, grid, params, rank))
			return rank;

	int coords[3];
	coords[0] = Particle_clamped_proc_coord(q.X[0], params->xmin,
	                                        params->xmax, params->NPX);
	coords[1] = Particle_clamped_proc_coord(q.X[1], params->ymin,
	                                        params->ymax, params->NPY);
#if defined(TWOD_MODE)
	coords[2] = 0;
#else
	coords[2] = Particle_clamped_proc_coord(q.X[2], params->zmin,
	                                        params->zmax, params->NPZ);
#endif

	int rank = -1;
	MPI_Cart_rank(PCW, coords, &rank);
	return rank;
}

static int Particle_interval_overlap_shift(double center, double radius,
                                           double lower, double upper,
                                           double period, int periodic,
                                           double *offset_out)
{
	if (offset_out != NULL)
		*offset_out = 0.0;

	if (!periodic || period <= 0.0)
		return (center + radius >= lower && center - radius < upper);

	const int shifts[3] = {0, -1, 1};
	for (int s = 0; s < 3; s++) {
		const double offset = (double)shifts[s] * period;
		const double c = center + offset;
		if (c + radius >= lower && c - radius < upper)
		{
			if (offset_out != NULL)
				*offset_out = offset;
			return 1;
		}
	}

	return 0;
}

static int Particle_overlaps_rank_shift(const Particle *p, MAC_grid *grid,
                                        Parameters *params, int rank,
                                        double extra_range, double shift[3])
{
	double lower[3], upper[3];
	int coords[3];
	Particle_rank_bounds(grid, params, rank, lower, upper, coords);

	const double radius = p->R + extra_range;
	double local_shift[3] = {0.0, 0.0, 0.0};

#ifdef XPERIODIC
	const int xperiodic = 1;
#else
	const int xperiodic = 0;
#endif
#ifdef YPERIODIC
	const int yperiodic = 1;
#else
	const int yperiodic = 0;
#endif
#ifdef ZPERIODIC
	const int zperiodic = 1;
#else
	const int zperiodic = 0;
#endif

	if (!Particle_interval_overlap_shift(p->X[0], radius, lower[0], upper[0],
	                                     params->xmax - params->xmin,
	                                     xperiodic, &local_shift[0]))
		return 0;
	if (!Particle_interval_overlap_shift(p->X[1], radius, lower[1], upper[1],
	                                     params->ymax - params->ymin,
	                                     yperiodic, &local_shift[1]))
		return 0;
#if defined(TWOD_MODE)
	if (shift != NULL)
		memcpy(shift, local_shift, 3 * sizeof(double));
	return 1;
#else
	if (!Particle_interval_overlap_shift(p->X[2], radius, lower[2], upper[2],
	                                     params->zmax - params->zmin,
	                                     zperiodic, &local_shift[2]))
		return 0;
	if (shift != NULL)
		memcpy(shift, local_shift, 3 * sizeof(double));
	return 1;
#endif
}

static int Particle_overlaps_rank(const Particle *p, MAC_grid *grid,
                                  Parameters *params, int rank,
                                  double extra_range)
{
	return Particle_overlaps_rank_shift(p, grid, params, rank, extra_range,
	                                    NULL);
}

Particle *Particle_collect_owned_overlaps(Particle_list *p_list,
                                          Cart3d_bag *data_bag,
                                          double extra_range,
                                          double min_radius,
                                          int include_self,
                                          int *n_recv)
{
	MAC_grid *grid = data_bag->grid;
	Parameters *params = data_bag->params;
	int nproc;
	MPI_Comm_size(PCW, &nproc);

	int *send_counts = (int *)calloc(nproc, sizeof(int));
	int *recv_counts = (int *)calloc(nproc, sizeof(int));
	int *send_displs = (int *)calloc(nproc, sizeof(int));
	int *recv_displs = (int *)calloc(nproc, sizeof(int));
	Memory_check_allocation(send_counts);
	Memory_check_allocation(recv_counts);
	Memory_check_allocation(send_displs);
	Memory_check_allocation(recv_displs);

	for (Particle *p = p_list->start; p != NULL; p = p->next)
		if (p->R > min_radius &&
		    Particle_center_is_owned_by_rank(p, grid, params, params->rank))
			for (int rank = 0; rank < nproc; rank++)
				if ((include_self || rank != params->rank) &&
				    Particle_overlaps_rank(p, grid, params, rank, extra_range))
					send_counts[rank]++;

	MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, PCW);

	int total_send = 0;
	int total_recv = 0;
	for (int rank = 0; rank < nproc; rank++) {
		send_displs[rank] = total_send;
		recv_displs[rank] = total_recv;
		total_send += send_counts[rank];
		total_recv += recv_counts[rank];
	}

	Particle *send_buf =
	    (Particle *)malloc((total_send > 0 ? total_send : 1) * sizeof(Particle));
	Particle *recv_buf =
	    (Particle *)malloc((total_recv > 0 ? total_recv : 1) * sizeof(Particle));
	Memory_check_allocation(send_buf);
	Memory_check_allocation(recv_buf);

	int *offset = (int *)malloc(nproc * sizeof(int));
	Memory_check_allocation(offset);
	memcpy(offset, send_displs, nproc * sizeof(int));

	for (Particle *p = p_list->start; p != NULL; p = p->next) {
		if (p->R <= min_radius ||
		    !Particle_center_is_owned_by_rank(p, grid, params, params->rank))
			continue;

		for (int rank = 0; rank < nproc; rank++) {
			double shift[3] = {0.0, 0.0, 0.0};
			if ((!include_self && rank == params->rank) ||
			    !Particle_overlaps_rank_shift(p, grid, params, rank,
			                                  extra_range, shift))
				continue;

			send_buf[offset[rank]] = *p;
			for (int d = 0; d < 3; d++) {
				send_buf[offset[rank]].X[d] += shift[d];
				send_buf[offset[rank]].X_old[d] += shift[d];
			}
			NULLIFY_PARTICLE_PTRS(send_buf[offset[rank]]);
			offset[rank]++;
		}
	}

	MPI_Alltoallv(send_buf, send_counts, send_displs, MPI_PARTICLE,
	              recv_buf, recv_counts, recv_displs, MPI_PARTICLE, PCW);

	free(send_buf);
	free(send_counts);
	free(recv_counts);
	free(send_displs);
	free(recv_displs);
	free(offset);

	*n_recv = total_recv;
	if (total_recv == 0) {
		free(recv_buf);
		return NULL;
	}
	return recv_buf;
}

static double Particle_min_width_for_rank(MAC_grid *grid, Parameters *params,
                                          int rank)
{
	double lower[3], upper[3];
	int coords[3];
	Particle_rank_bounds(grid, params, rank, lower, upper, coords);

	double width = upper[0] - lower[0];
	width = fmin(width, upper[1] - lower[1]);
#if !defined(TWOD_MODE)
	if (upper[2] > lower[2])
		width = fmin(width, upper[2] - lower[2]);
#endif
	return width;
}

static int Particle_is_oversized_for_owner(const Particle *p, MAC_grid *grid,
                                           Parameters *params, int owner_rank,
                                           double range)
{
	double threshold = Particle_min_width_for_rank(grid, params, owner_rank) - range;
	if (threshold < 0.0) threshold = 0.0;
	return (p->R > threshold);
}

static int Particle_list_has_id(Particle_list *p_list, int id)
{
	for (Particle *p = p_list->start; p != NULL; p = p->next)
		if (p->ID == id)
			return 1;
	return 0;
}

static void Particle_add_flat_copy(Particle_list *p_list, const Particle *src)
{
	Particle *pnew = (Particle *)malloc(sizeof(Particle));
	Memory_check_allocation(pnew);
	*pnew = *src;
	Particle_create_internal_arrays(pnew);
	pnew->particle_collision = NULL;
	pnew->wall_collision = NULL;
	pnew->next = p_list->start;
	p_list->start = pnew;
}

static void Particle_exchange_oversized_overlaps(Particle_list *p_list,
                                                Cart3d_bag *data_bag,
                                                double sub_min,
                                                double range)
{
	double threshold = sub_min - range;
	if (threshold < 0.0) threshold = 0.0;

	int nrecv = 0;
	Particle *recv =
	    Particle_collect_owned_overlaps(p_list, data_bag, 0.0, threshold,
	                                    0, &nrecv);

	for (int n = 0; n < nrecv; n++) {
		if (!Particle_list_has_id(p_list, recv[n].ID))
			Particle_add_flat_copy(p_list, &recv[n]);
	}
	free(recv);
}

static void Particle_pack_force_record(const Particle *p,
                                       Particle_force_record *record)
{
	int k = 0;
	record->ID = p->ID;
	for (int i = 0; i < 3; i++) record->value[k++] = p->F[i];
	for (int i = 0; i < 3; i++) record->value[k++] = p->T[i];
	for (int i = 0; i < 3; i++) record->value[k++] = p->Int_U[i];
	for (int i = 0; i < 3; i++) record->value[k++] = p->Int_Omega[i];
#ifdef VOF_IBM
	for (int i = 0; i < 3; i++) record->value[k++] = p->F_CCF[i];
	for (int i = 0; i < 3; i++) record->value[k++] = p->T_CCF[i];
	for (int i = 0; i < 3; i++) record->value[k++] = p->Int_rho[i];
	record->value[k++] = p->Int_rho_scalar;
	for (int i = 0; i < 3; i++) record->value[k++] = p->F_CSF_solid[i];
	for (int i = 0; i < 3; i++) record->value[k++] = p->T_CSF_solid[i];
#endif
}

static void Particle_zero_hydro_force_fields(Particle *p)
{
	DSET_ZERO(p->F, 3);
	DSET_ZERO(p->T, 3);
	DSET_ZERO(p->Int_U, 3);
	DSET_ZERO(p->Int_Omega, 3);
#ifdef VOF_IBM
	DSET_ZERO(p->F_CCF, 3);
	DSET_ZERO(p->T_CCF, 3);
	DSET_ZERO(p->Int_rho, 3);
	p->Int_rho_scalar = 0.0;
	DSET_ZERO(p->F_CSF_solid, 3);
	DSET_ZERO(p->T_CSF_solid, 3);
#endif
}

static void Particle_add_force_record(Particle *p,
                                      const Particle_force_record *record)
{
	int k = 0;
	for (int i = 0; i < 3; i++) p->F[i] += record->value[k++];
	for (int i = 0; i < 3; i++) p->T[i] += record->value[k++];
	for (int i = 0; i < 3; i++) p->Int_U[i] += record->value[k++];
	for (int i = 0; i < 3; i++) p->Int_Omega[i] += record->value[k++];
#ifdef VOF_IBM
	for (int i = 0; i < 3; i++) p->F_CCF[i] += record->value[k++];
	for (int i = 0; i < 3; i++) p->T_CCF[i] += record->value[k++];
	for (int i = 0; i < 3; i++) p->Int_rho[i] += record->value[k++];
	p->Int_rho_scalar += record->value[k++];
	for (int i = 0; i < 3; i++) p->F_CSF_solid[i] += record->value[k++];
	for (int i = 0; i < 3; i++) p->T_CSF_solid[i] += record->value[k++];
#endif
}

void Particle_reduce_oversized_forces_to_owner(Particle_list *p_list,
                                               Cart3d_bag *data_bag)
{
	MAC_grid *grid = data_bag->grid;
	Parameters *params = data_bag->params;
	const double range = DELTA_FUNC_RADIUS * grid->dx_u[1];

	int nproc;
	MPI_Comm_size(PCW, &nproc);

	int *send_counts = (int *)calloc(nproc, sizeof(int));
	int *recv_counts = (int *)calloc(nproc, sizeof(int));
	int *send_displs = (int *)calloc(nproc, sizeof(int));
	int *recv_displs = (int *)calloc(nproc, sizeof(int));
	Memory_check_allocation(send_counts);
	Memory_check_allocation(recv_counts);
	Memory_check_allocation(send_displs);
	Memory_check_allocation(recv_displs);

	for (Particle *p = p_list->start; p != NULL; p = p->next) {
		const int owner = Particle_center_owner_rank(p, grid, params);
		if (owner >= 0 && Particle_is_oversized_for_owner(p, grid, params,
		                                                  owner, range))
			send_counts[owner]++;
	}

	MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, PCW);

	int total_send = 0;
	int total_recv = 0;
	for (int rank = 0; rank < nproc; rank++) {
		send_displs[rank] = total_send;
		recv_displs[rank] = total_recv;
		total_send += send_counts[rank];
		total_recv += recv_counts[rank];
	}

	Particle_force_record *send_buf =
	    (Particle_force_record *)malloc((total_send > 0 ? total_send : 1) *
	                                    sizeof(Particle_force_record));
	Particle_force_record *recv_buf =
	    (Particle_force_record *)malloc((total_recv > 0 ? total_recv : 1) *
	                                    sizeof(Particle_force_record));
	Memory_check_allocation(send_buf);
	Memory_check_allocation(recv_buf);

	int *offset = (int *)malloc(nproc * sizeof(int));
	Memory_check_allocation(offset);
	memcpy(offset, send_displs, nproc * sizeof(int));

	for (Particle *p = p_list->start; p != NULL; p = p->next) {
		const int owner = Particle_center_owner_rank(p, grid, params);
		if (owner < 0 || !Particle_is_oversized_for_owner(p, grid, params,
		                                                  owner, range))
			continue;

		Particle_pack_force_record(p, &send_buf[offset[owner]]);
		offset[owner]++;
	}

	int *send_bytes = (int *)malloc(nproc * sizeof(int));
	int *recv_bytes = (int *)malloc(nproc * sizeof(int));
	int *sdispl_bytes = (int *)malloc(nproc * sizeof(int));
	int *rdispl_bytes = (int *)malloc(nproc * sizeof(int));
	Memory_check_allocation(send_bytes);
	Memory_check_allocation(recv_bytes);
	Memory_check_allocation(sdispl_bytes);
	Memory_check_allocation(rdispl_bytes);
	for (int rank = 0; rank < nproc; rank++) {
		send_bytes[rank] = send_counts[rank] * (int)sizeof(Particle_force_record);
		recv_bytes[rank] = recv_counts[rank] * (int)sizeof(Particle_force_record);
		sdispl_bytes[rank] = send_displs[rank] * (int)sizeof(Particle_force_record);
		rdispl_bytes[rank] = recv_displs[rank] * (int)sizeof(Particle_force_record);
	}

	MPI_Alltoallv(send_buf, send_bytes, sdispl_bytes, MPI_BYTE,
	              recv_buf, recv_bytes, rdispl_bytes, MPI_BYTE, PCW);

	for (Particle *p = p_list->start; p != NULL; p = p->next) {
		const int owner = Particle_center_owner_rank(p, grid, params);
		if (owner >= 0 && Particle_is_oversized_for_owner(p, grid, params,
		                                                  owner, range))
			Particle_zero_hydro_force_fields(p);
	}

	for (int n = 0; n < total_recv; n++) {
		for (Particle *p = p_list->start; p != NULL; p = p->next) {
			if (p->ID == recv_buf[n].ID &&
			    Particle_center_is_owned_by_rank(p, grid, params, params->rank)) {
				Particle_add_force_record(p, &recv_buf[n]);
				break;
			}
		}
	}

	free(send_buf);
	free(recv_buf);
	free(send_counts);
	free(recv_counts);
	free(send_displs);
	free(recv_displs);
	free(send_bytes);
	free(recv_bytes);
	free(sdispl_bytes);
	free(rdispl_bytes);
	free(offset);
}
#endif /* LAG_PARTICLE_RESOLVED */

/******************************************************************************/
/*
 Update and transmit particles between processors

 Inputs:
     p_list - particle linked list structure to communicate
	 p_list->state - integer indicating the current state of the linked list:
             LIST_STATE_LOCAL - linked list contains only local particles
             LIST_STATE_FOREIGN - linked list contains only foreign particles
             LIST_STATE_BOTH - linked list contains both local and foreign
			       particles
             LIST_STATE_EDGE - linked list contains both local and foreign
			       particles, but only near the processor's edges
	 line - use macro '__LINE__' to get line number of function call
	 file - use macro '__FILE__' to get filename of function call

 Particles are sent to neighboring processor if they extend beyond that
 neighboring edge.
     - Ghost nodes are not taken into account for these considerations.
     - Particle coordinates ('X') are updated for periodic boundaries.

 Assumptions:
     - No y-periodicity
     - uniform grid

 Potential speedups:
     - Use same send buffer for both lower and upper communications
     - Have Np be an input parameter (then Particle_read() should output Np)
 */
/******************************************************************************/
void Particle_MPI_update(Particle_list *p_list, Cart3d_bag *data_bag,
		Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;

	// Grid spacing
	double h = grid -> dx_u[1];

	// Range outside processor boundaries particle has an influence
	double range = DELTA_FUNC_RADIUS * h;
#if defined LAG_MARKER_FLAG || defined LAG_MARKER_PRIORITY
	range += LAG_FLAG_RANGE * h;
#endif
#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL
	range = max(range, (0.5*params->lub_range*h)); // maybe add h for help points
#endif
#ifdef ELECTROSTATIC_REPULSION
	range = max(range, (0.5*params->erm_range*h));
#endif
#ifdef COHESION
	range = max(range, (0.5*params->coh_range*h));
#endif


	// Grid bounds of local subdomain
	double G_xmin = grid -> xu[grid -> G_Is] + range;
	double G_ymin = grid -> yv[grid -> G_Js] + range;
	double G_zmin = grid -> zw[grid -> G_Ks] + range;
	double G_xmax = grid -> xu[min(grid -> G_Ie, grid -> NX - 1)] - range;
	double G_ymax = grid -> yv[min(grid -> G_Je, grid -> NY - 1)] - range;
	double G_zmax = grid -> zw[min(grid -> G_Ke, grid -> NZ - 1)] - range;

	// Periodic parameters
#ifdef XPERIODIC
	double Lx = params -> Lx;
	int xproccoord = params -> xproccoord;
	int NPX = params -> NPX;
#endif
#ifdef YPERIODIC
	double Ly = params -> Ly;
	int yproccoord = params -> yproccoord;
	int NPY = params -> NPY;
#endif
#ifdef ZPERIODIC
	double Lz = params -> Lz;
	int NPZ = params -> NPZ;
	int zproccoord = params -> zproccoord;
#endif

	// Particle parameters
//	int Np;  // Number of particles on current processor
	Particle *p;

	// Number of particles to send to/receive from other processor
	int send_Np, recv_Np;

	// Buffers for sending/receiving particles
	Particle *send_part = lag -> send_part;
	Particle *recv_part = lag -> recv_part;

	Collision *pc;

	// Number of collisions to send to/receive from other processor
	int *send_Nc = lag -> send_Nc;
	int *recv_Nc = lag -> recv_Nc;

	// Buffers for sending/receiving collisions
	Collision *send_coll = lag -> send_coll;
	Collision *recv_coll = lag -> recv_coll;

	// MPI
	int src, dest;
	int tag = 1;
	MPI_Status status;

	// Error message
	char message[100], name0[50], name1[50], name2[50];

	//--------------------------------------------------------------------------
	// Check current linked list state and assign its new designation
	//--------------------------------------------------------------------------
	// Starting with local particles, list will soon have both local and
	// foreign particles
	if (p_list->state == LIST_STATE_LOCAL) {
		p_list->state = LIST_STATE_BOTH;
	}
	// Starting with foreign particles, list will soon have edge particles,
	// both local and foreign
	else if (p_list->state == LIST_STATE_FOREIGN) {
		p_list->state = LIST_STATE_EDGE;
	}
	// Any other initial state will produce a useless set of particles
	else {
		Display_print_list_name(name0, p_list->state);
		Display_print_list_name(name1, LIST_STATE_LOCAL);
		Display_print_list_name(name2, LIST_STATE_FOREIGN);
		sprintf(message, "Invalid input linked list state '%s'\n"
				"Should be '%s' or '%s'", name0, name1, name2);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}


	/*------------------------------------------------------------------------*/
	/*
	 Lower x communication
	 */
	/*------------------------------------------------------------------------*/
#ifdef XPERIODIC
	//--------------------------------------------------------------------------
	// Ensure we don't communicate foreign particles to ourself, which would
	// effectively duplicate the particles and cause double counting of forces
	//--------------------------------------------------------------------------
	if (params -> NPX != 1 || p_list->state != LIST_STATE_EDGE) {
#endif
		//----------------------------------------------------------------------
		// Fill send buffer
		//----------------------------------------------------------------------
		send_Np = 0;
		send_Nc[0] = 0;
		p = p_list->start;
		while (p != NULL) {

			// Send nearby particles
			if (p->X[0] - R_EFF < G_xmin) {

				ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
				NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef XPERIODIC
				// Update particle position to reflect periodicity
				if (xproccoord == 0) {
					send_part[send_Np].X[0]     += Lx;
					send_part[send_Np].X_old[0] += Lx;
				}
#endif
				FILL_COLLISION_SEND_BUFFER();
				send_Np++;
			}

			p = p -> next;
		}

		//----------------------------------------------------------------------
		// Communicate particles and collisions
		//----------------------------------------------------------------------
		dest = params->npxminus;
		src  = params->npxplus;

		SENDRECV_PARTICLES_COLLISIONS();


		/*--------------------------------------------------------------------*/
		/*
		 Upper x communication
		 */
		/*--------------------------------------------------------------------*/

		//----------------------------------------------------------------------
		// Fill send buffer
		//----------------------------------------------------------------------
		send_Np = 0;
		send_Nc[0] = 0;
		p = p_list->start;
		while (p != NULL) {

			// Send nearby particles
			if (p->X[0] + R_EFF > G_xmax) {

				ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
				NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef XPERIODIC
				// Update particle position to reflect periodicity
				if (xproccoord == NPX - 1) {
					send_part[send_Np].X[0]     -= Lx;
					send_part[send_Np].X_old[0] -= Lx;
				}
#endif
				FILL_COLLISION_SEND_BUFFER();
				send_Np++;
			}

			p = p -> next;
		}

		//----------------------------------------------------------------------
		// Add received particles from previous communication to linked list.
		// This placement (after filling send buffer) ensures that we do not
		// send particles we just received back to the host.
		//----------------------------------------------------------------------
		Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

		//----------------------------------------------------------------------
		// Communicate particles and collisions
		//----------------------------------------------------------------------
		dest = params->npxplus;
		src  = params->npxminus;

		SENDRECV_PARTICLES_COLLISIONS();

		//----------------------------------------------------------------------
		// Add received particles to linked list
		//----------------------------------------------------------------------
		Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

#ifdef XPERIODIC
	} // if not communicating foreign particles to self
#endif
/*------------------------------------------------------------------------*/
/*
 Lower y communication
 */
/*------------------------------------------------------------------------*/
#ifdef YPERIODIC
//--------------------------------------------------------------------------
// Ensure we don't communicate foreign particles to ourself, which would
// effectively duplicate the particles and cause double counting of forces
//--------------------------------------------------------------------------
if (params -> NPY != 1 || p_list->state != LIST_STATE_EDGE) {
#endif
	//----------------------------------------------------------------------
	// Fill send buffer
	//----------------------------------------------------------------------
	send_Np = 0;
	send_Nc[0] = 0;
	p = p_list->start;
	while (p != NULL) {

		// Send nearby particles
		if (p->X[1] - R_EFF < G_ymin) {

			ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
			NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef YPERIODIC
			// Update particle position to reflect periodicity
			if (yproccoord == 0) {
				send_part[send_Np].X[1]     += Ly;
				//double xcord= send_part[send_Np].X[0];
				//printf("%d\n",xcord);
				send_part[send_Np].X_old[1] += Ly;
			}
#endif
			FILL_COLLISION_SEND_BUFFER();
			send_Np++;
		}

		p = p -> next;
	}

	//----------------------------------------------------------------------
	// Communicate particles and collisions
	//----------------------------------------------------------------------
	dest = params->npyminus;
	src  = params->npyplus;

	SENDRECV_PARTICLES_COLLISIONS();


	/*--------------------------------------------------------------------*/
	/*
	 Upper y communication
	 */
	/*--------------------------------------------------------------------*/

	//----------------------------------------------------------------------
	// Fill send buffer
	//----------------------------------------------------------------------
	send_Np = 0;
	send_Nc[0] = 0;
	p = p_list->start;
	while (p != NULL) {

		// Send nearby particles
		if (p->X[1] + R_EFF > G_ymax) {

			ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
			NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef YPERIODIC
			// Update particle position to reflect periodicity
			if (yproccoord == NPY - 1) {
				send_part[send_Np].X[1]     -= Ly;
				//double xcord= send_part[send_Np].X[0];
				//printf("%d\n",xcord);
				send_part[send_Np].X_old[1] -= Ly;
			}
#endif
			FILL_COLLISION_SEND_BUFFER();
			send_Np++;
		}

		p = p -> next;
	}

	//----------------------------------------------------------------------
	// Add received particles from previous communication to linked list.
	// This placement (after filling send buffer) ensures that we do not
	// send particles we just received back to the host.
	//----------------------------------------------------------------------
	Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

	//----------------------------------------------------------------------
	// Communicate particles and collisions
	//----------------------------------------------------------------------
	dest = params->npyplus;
	src  = params->npyminus;

	SENDRECV_PARTICLES_COLLISIONS();

	//----------------------------------------------------------------------
	// Add received particles to linked list
	//----------------------------------------------------------------------
	Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

#ifdef YPERIODIC
} // if not communicating foreign particles to self
#endif

	/*------------------------------------------------------------------------*/
	/*
	 Lower y communication
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// Fill send buffer
	//--------------------------------------------------------------------------
	/*send_Np = 0;
	send_Nc[0] = 0;
	p = p_list->start;
	while (p != NULL) {

		// Send nearby particles
		if (p->X[1] - R_EFF < G_ymin) {

			ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
			NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
			FILL_COLLISION_SEND_BUFFER();
			send_Np++;
		}

		p = p -> next;
	}

	//--------------------------------------------------------------------------
	// Communicate particles and collisions
	//--------------------------------------------------------------------------
	dest = params->npyminus;
	src  = params->npyplus;

	SENDRECV_PARTICLES_COLLISIONS();


	/*------------------------------------------------------------------------*/
	/*
	 Upper y communication
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// Fill send buffer
	//--------------------------------------------------------------------------
	/*send_Np = 0;
	send_Nc[0] = 0;
	p = p_list->start;
	while (p != NULL) {

		// Send nearby particles
		if (p->X[1] + R_EFF > G_ymax) {

			ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
			NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
			FILL_COLLISION_SEND_BUFFER();
			send_Np++;
		}

		p = p -> next;
	}

	//--------------------------------------------------------------------------
	// Add received particles from previous communication to linked list.  This
	// placement (after filling send buffer) ensures that we do not send
	// particles we just received back to the host.
	//--------------------------------------------------------------------------
	Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

	//--------------------------------------------------------------------------
	// Communicate particles and collisions
	//--------------------------------------------------------------------------
	dest = params->npyplus;
	src  = params->npyminus;

	SENDRECV_PARTICLES_COLLISIONS();

	//--------------------------------------------------------------------------
	// Add received particles to linked list
	//--------------------------------------------------------------------------
	Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);*/


	/*------------------------------------------------------------------------*/
	/*
	 Lower z communication
	 */
	/*------------------------------------------------------------------------*/
#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
	/*
	 * Phase-3 2D IBM: particles never migrate in the collapsed storage
	 * direction.  X/Y exchanges above still provide all marker support needed
	 * for multi-rank in-plane runs.
	 */
#else
#ifdef ZPERIODIC
	//--------------------------------------------------------------------------
	// Ensure we don't communicate foreign particles to ourself, which would
	// effectively duplicate the particles and cause double counting of forces
	//--------------------------------------------------------------------------
	if (params -> NPZ != 1 || p_list->state != LIST_STATE_EDGE) {
#endif
		//----------------------------------------------------------------------
		// Fill send buffer
		//----------------------------------------------------------------------
		send_Np = 0;
		send_Nc[0] = 0;
		p = p_list->start;
		while (p != NULL) {

			// Send nearby particles
			if (p->X[2] - R_EFF < G_zmin) {

				ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
				NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef ZPERIODIC
				// Update particle position to reflect periodicity
				if (zproccoord == 0) {
					send_part[send_Np].X[2]     += Lz;
					send_part[send_Np].X_old[2] += Lz;
				}
#endif
				FILL_COLLISION_SEND_BUFFER();
				send_Np++;
			}

			p = p -> next;
		}

		//----------------------------------------------------------------------
		// Communicate particles and collisions
		//----------------------------------------------------------------------
		dest = params->npzminus;
		src  = params->npzplus;

		SENDRECV_PARTICLES_COLLISIONS();


		/*--------------------------------------------------------------------*/
		/*
		 Upper z communication
		 */
		/*--------------------------------------------------------------------*/

		//----------------------------------------------------------------------
		// Fill send buffer
		//----------------------------------------------------------------------
		send_Np = 0;
		send_Nc[0] = 0;
		p = p_list->start;
		while (p != NULL) {

			// Send nearby particles
			if (p->X[2] + R_EFF > G_zmax) {

				ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
				NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef ZPERIODIC
				// Update particle position to reflect periodicity
				if (zproccoord == NPZ - 1) {
					send_part[send_Np].X[2]     -= Lz;
					send_part[send_Np].X_old[2] -= Lz;
				}
#endif
				FILL_COLLISION_SEND_BUFFER();
				send_Np++;
			}

			p = p -> next;
		}

		//----------------------------------------------------------------------
		// Add received particles from previous communication to linked list.
		// This placement (after filling send buffer) ensures that we do not
		// send particles we just received back to the host.
		//----------------------------------------------------------------------
		Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

		//----------------------------------------------------------------------
		// Communicate particles and collisions
		//----------------------------------------------------------------------
		dest = params->npzplus;
		src  = params->npzminus;

		SENDRECV_PARTICLES_COLLISIONS();

		//----------------------------------------------------------------------
		// Add received particles to linked list
		//----------------------------------------------------------------------
		Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

#ifdef ZPERIODIC
	} // if not communicating foreign particles to self
#endif
#endif

#if defined(LAG_PARTICLE_RESOLVED)
	/* Exchange only the oversized particles that geometrically overlap this rank. */
	if (p_list->state == LIST_STATE_BOTH) {
		double _sub_x = grid->xu[min(grid->G_Ie, grid->NX-1)] - grid->xu[grid->G_Is];
		double _sub_y = grid->yv[min(grid->G_Je, grid->NY-1)] - grid->yv[grid->G_Js];
		double _sub_min = fmin(_sub_x, _sub_y);
#if !defined(TWOD_MODE)
		double _sub_z = grid->zw[min(grid->G_Ke, grid->NZ-1)] - grid->zw[grid->G_Ks];
		if (_sub_z > 1.0e-30) _sub_min = fmin(_sub_min, _sub_z);
#endif
		Particle_exchange_oversized_overlaps(p_list, data_bag, _sub_min, range);
	}
#endif

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_particle_comm += T2 - T1;
}
