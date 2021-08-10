#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf5.h"

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Display.h"
#include "Interpolate.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "Particle.h"
#include "ParticleInput.h"

MPI_Datatype MPI_PARTICLE;

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
	Particle_list_remove(p_fixed_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));

	// Create linked list of foreign particles
	p_list_foreign = Particle_list_foreign_create(p_mobile_list, data_bag, DTRACE("Particle_list_foreign_create"));

	// Receive and combine portions of particles that were on other subdomains
	Particle_MPI_update(p_list_foreign, data_bag, DTRACE("Particle_MPI_update"));
	Lagrangian_collect_forces(p_mobile_list, p_list_foreign, LAG_COLLECT_HYDRO, params, DTRACE("Lagrangian_collect_forces"));
	Particle_list_destroy(p_list_foreign);

	// Set particle velocities based on volume-averaged fluid velocities
	p = p_mobile_list -> start;
	while (p != NULL) {

		FORI3 p->U[i]     = p->Int_U[i] * p->rho_s / p->M;
		FORI3 p->Omega[i] = p->Int_Omega[i] * p->rho_s / p->I_p;

		FORI3 p->U_old[i]     = p->U[i];
		FORI3 p->Omega_old[i] = p->Omega[i];

#ifdef ONE_WAY
		FORI3 p->U[i]     = 0.0;
		FORI3 p->Omega[i] = 0.0;

		FORI3 p->U_old[i]     = 0.0;
		FORI3 p->Omega_old[i] = 0.0;
#endif

		FORI3 p->Int_U_old[i]     = p->Int_U[i];
		FORI3 p->Int_Omega_old[i] = p->Int_Omega[i];

		p = p -> next;
	}
#endif
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

	// Sphere derived data
	double R2  = R * R;
	double M   = rho_s * 4.0 / 3.0 * PI * R2 * R;
	double I_p = rho_s * 8.0 / 15.0 * PI * R2 * R2 * R;

	// Number of Lagrangian marker points determined by grid-to-particle
	// resolution

	double h = grid -> dy_min;

	int N_L = (int)ceil(PI / 3.0 * (12.0 * R2 / (h * h) + 1.0));

	// Volume controlled by each Lagrangian marker point
	double Vol_L = PI * h / (3.0 * N_L) * (12.0 * R2 + h * h);

	p -> rho_s = rho_s;
	p -> M   = M;
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

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_particle_comm += T2 - T1;
}
