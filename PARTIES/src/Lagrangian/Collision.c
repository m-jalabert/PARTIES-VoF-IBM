

/******************************************************************************/
/*
 This file contains functions that are common to some or all collision models
 */
/******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Cart3d.h"
#include "Collision.h"
#include "Display.h"
#include "Particle.h"

#define GRAV_OVERLAP 0.001  // Fraction of radius to overlap due to gravity
#define ST_CRIT 5.0         // St number below which spring stiffness is fixed
#define ST_FLUID_CRIT 5.0   // St number below which fluid forces do not act

#ifdef ACTM
	#include "collision_models/actm.c"
#elif defined DEM
	#include "collision_models/dem.c"
#else
	#include "collision_models/repulsive.c"
#endif
#include "collision_models/lubrication.c"
#ifdef ATFM
	#include "collision_models/atfm.c"
#elif defined LIN_TAN
	#include "collision_models/lin_tan.c"
#endif
#ifdef ELECTROSTATIC_REPULSION
	#include "collision_models/erm.c"
#endif

// Distance between two particles, squared
#define PP_DIST_SQ(a,b) ( (a[0] - b[0]) * (a[0] - b[0]) + \
                          (a[1] - b[1]) * (a[1] - b[1]) + \
                          (a[2] - b[2]) * (a[2] - b[2]) )

#define MIN_GT_CP_NORM 1e-14

int Collision_get_wall_ID(int dim, int side);
double Collision_pw_distance(int wall_ID, int dim, double *X, Parameters *params);
int Collision_check_ownership(double *midpoint, int stage, MAC_grid *grid);
Collision *Collision_pp_find(Particle *p, Particle *p2);
Collision *Collision_pw_find(Particle *p, int wall_ID);
void Collision_pp_new(Collision **pc_ptr, Collision **pc2_ptr,
		Collision_bag *bag, Parameters *params);
Collision *Collision_pw_new(int wall_ID, Collision_bag *bag, Parameters *params);
void Collision_fill_bag_particle(Collision_bag *bag, double center_distance,
		double surface_distance, Parameters *params);
void Collision_fill_bag_wall(Collision_bag *bag, int dim, int side,
		double center_distance, double surface_distance, Parameters *params);


/******************************************************************************/
/*
 * Checks all collisions associated with particle.
 *     - returns 0 if no collisions are above St = ST_FLUID_CRIT
 *     - returns 1 if at least one collision is above St = ST_FLUID_CRIT
 */
/******************************************************************************/
int Collision_above_critical(Particle *p) {

	Collision *pc;
	int St_gt_crit = 0;

	pc = p -> particle_collision;
	while (pc != NULL && St_gt_crit == 0) {
		if (pc -> St > ST_FLUID_CRIT)
			St_gt_crit = 1;
		pc = pc -> next;
	}

	pc = p -> wall_collision;
	while (pc != NULL && St_gt_crit == 0) {
		if (pc -> St > ST_FLUID_CRIT)
			St_gt_crit = 1;
		pc = pc -> next;
	}

	return St_gt_crit;
}




/******************************************************************************/
/*
 * Evaluates collision forces for p_mobile linked list
 *     - Sets 'p_mobile_list' list to a local state
 *     - Returns 'p_mobile_list_foreign' -- a list of nearby foreign particles
 * Net collision forces acting on the particle are stored in Fc, including both
 * contact and lubrication forces
 */
/******************************************************************************/
Particle_list *Collision_evaluate(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2;
	Particle *p, *pf;
	Particle_list *p_mobile_list_foreign;
	Collision *pc;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_BOTH, params,
		DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_BOTH, params,
		DTRACE("Display_assert_list_state"));

	// Reset mobile collisions
	p = p_mobile_list -> start;
	while (p != NULL) {
		DSET_ZERO(p->Fc, 3);
		DSET_ZERO(p->Tc, 3);
#ifdef POST_PROCESS
		DSET_ZERO(p->Fc_norm, 3);
		DSET_ZERO(p->Fc_tan, 3);
		DSET_ZERO(p->Fl_norm, 3);
		DSET_ZERO(p->Fl_tan, 3);
#endif

		pc = p -> particle_collision;
		while (pc != NULL) {
			pc->state = COLL_STATE_NEUTRAL;
			pc = pc -> next;
		}
		pc = p -> wall_collision;
		while (pc != NULL) {
			pc->state = COLL_STATE_NEUTRAL;
			pc = pc -> next;
		}

		p = p -> next;
	}
	// Reset fixed collisions
	p = p_fixed_list -> start;
	while (p != NULL) {
		DSET_ZERO(p->Fc, 3);
		DSET_ZERO(p->Tc, 3);
#ifdef POST_PROCESS
		DSET_ZERO(p->Fc_norm, 3);
		DSET_ZERO(p->Fc_tan, 3);
		DSET_ZERO(p->Fl_norm, 3);
		DSET_ZERO(p->Fl_tan, 3);
#endif

		pc = p -> particle_collision;
		while (pc != NULL) {
			pc->state = COLL_STATE_NEUTRAL;
			pc = pc -> next;
		}

		p = p -> next;
	}

	// Create linked list of foreign particles
	p_mobile_list_foreign = Particle_list_foreign_create(p_mobile_list, data_bag,
		DTRACE("Particle_list_foreign_create"));

	T1 = MPI_Wtime();
	//--------------------------------------------------------------------------
	// Collisions with local particles
	//--------------------------------------------------------------------------
	p = p_mobile_list -> start;
	while (p != NULL) {

		Collision_particle(p, p_mobile_list->start, COLL_STAGE_MOBILE, grid, params);
		Collision_particle(p, p_fixed_list->start, COLL_STAGE_FIXED, grid, params);
		Collision_particle(p, p_mobile_list_foreign->start, COLL_STAGE_FOREIGN, grid, params);
		Collision_wall(p, grid, params);
		p = p -> next;
	}
	//--------------------------------------------------------------------------
	// Collisions with foreign particles
	//--------------------------------------------------------------------------
	p = p_mobile_list_foreign->start;
	while (p != NULL) {

		Collision_particle(p, p_fixed_list->start, COLL_STAGE_FIXED, grid, params);
		Collision_particle(p, p_mobile_list_foreign->start, COLL_STAGE_FOREIGN, grid, params);
		p = p -> next;
	}
	T2 = MPI_Wtime();
	data_bag->timer->Wtime_particle_coll += T2 - T1;

	return p_mobile_list_foreign;
}

/******************************************************************************/
/*
 * Calculates collision forces and torques between particle 'p' and all the
 * particles in 'p_start'.
 *
 * If this processor owns the collision (the contact point between particles 'p'
 * and 'p2' exists on this processor), the forces and torques are added to
 * 'p->Fc', 'p2->Fc', 'p->Tc', and 'p2->Tc' (for both particles).
 *
 * Assumes uniform square grid
 */
/******************************************************************************/
void Collision_particle(Particle *p, Particle *p_start, int stage,
		MAC_grid *grid, Parameters *params) {

	int i;
	char message[500];

	// "Other" particle for looping through stages
	Particle *p2;

#ifdef STORE_COLLISION
	Collision *pc, *pc2;
#endif

	Collision_bag *bag = (Collision_bag *)malloc(sizeof(Collision_bag));
	bag -> p = p;
	bag -> stage = stage;

	// Distance between particle centers, and its square
	double center_distance,    center_distance_sq;

	// Distance particle centers have to be within to have a collision
	double center_cutoff_dist, center_cutoff_dist_sq;

	// Distance between particle centers below which the long range interaction
	// forces (e.g. lubrication, electrostatic repulsion) act
	double cutoff_dist, cutoff_dist_sq;

	// Distance between particle surfaces, negative value indicates overlap
	double surface_distance;

	// Surface distance at which particles begin to interact
	double interaction_range = 0.0;

	// Grid spacing (assuming uniform grid)
	double h = grid->dx_u[1];

	// Midpoint of collision
	double midpoint[3];

	// Flag for specifying if collision was just created
	int new_collision;

	double temp;

	// Distance at which collision takes place
	double collision_offset = 0.0;

#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL
	interaction_range = max(interaction_range, params->lub_range * h);
#endif
#ifdef ELECTROSTATIC_REPULSION
	interaction_range = max(interaction_range, params->erm_range * h);
#endif
#ifdef COHESION
	interaction_range = max(interaction_range, params->coh_range * h);
#endif

	/*------------------------------------------------------------------------*/
	/*
	 * Calculate forces due to particle-particle interaction between
	 * mobile and local mobile particles -> stage = 0
	 * mobile and fixed particles -> stage = 1 (foreign fixed particles are also kept locally)
	 * mobile and foreign mobile particles -> stage = 2
	 */
	/*------------------------------------------------------------------------*/
	p2 = p_start; // start at local particles
	while (p2 != NULL) {

		//----------------------------------------------------------------------
		// Do not count self-collisions, so loop over all pairs of mobile
		// particles once, in total
		//----------------------------------------------------------------------
		if(p2 == p) {
			break;
		}

		bag -> p2 = p2;
#ifdef ROUGH_COLLISION
		collision_offset = params->roughness * 0.5 * (p->R + p2->R);
#endif

		// Distance particle centers have to be within to have a collision
		center_cutoff_dist = p->R + p2->R + collision_offset;
		center_cutoff_dist_sq = center_cutoff_dist * center_cutoff_dist;

		center_distance_sq = PP_DIST_SQ(p->X, p2->X);
		if (center_distance_sq <= center_cutoff_dist_sq) { // we're in contact

			center_distance = sqrt(center_distance_sq);

			// Point in the middle between particle surfaces
			temp = 0.5 * (center_distance + p->R - p2->R) / center_distance;
			FORI3 midpoint[i] = p->X[i] + temp * (p2->X[i] - p->X[i]);

			// Find out if this processor is supposed to own this collision
			if (Collision_check_ownership(midpoint, stage, grid)) {

				surface_distance = center_distance - center_cutoff_dist;
				Collision_fill_bag_particle(bag, center_distance, surface_distance, params);

#ifdef STORE_COLLISION
				// Look for existing particle collisions
				pc  = Collision_pp_find(p, p2);
				pc2 = Collision_pp_find(p2, p);

				// If pc does not exist, create, calculate and link it
				if (pc == NULL && pc2 == NULL) {
					// Make sure particles are moving towards each other.  The
					// alternate case can occur in rare cases when gravity is
					// involved
					if (bag->g_dot_n > 0) {
						new_collision = 1;
						Collision_pp_new(&pc, &pc2, bag, params);
					}
					else {
						p2 = p2 -> next;
						continue;
					}
				}
				else if (pc == NULL) {
					sprintf(message, "Warning: pc == NULL for p->ID = %d, p2->ID = %d\n", p->ID, p2->ID);
					Display_error(message);
					fflush(stdout);
				}
				else if (pc2 == NULL) {
					sprintf(message, "Warning: pc2 == NULL for p->ID = %d, p2->ID = %d\n", p->ID, p2->ID);
					Display_error(message);
					fflush(stdout);
				}
				else {
					new_collision = 0;
				}
				if(pc  != NULL) pc->state = COLL_STATE_OWNER;
				if(pc2 != NULL) pc2->state = COLL_STATE_OWNER;
#endif

				//--------------------------------------------------------------
				// Normal collision force
				//--------------------------------------------------------------
#ifdef ACTM
				actm(pc, bag, params);
#else
				dem(pc, bag, params);
#endif

				//--------------------------------------------------------------
				// Tangential collision force
				//--------------------------------------------------------------
#ifdef ATFM
				if (bag->gt_cp_norm >= MIN_GT_CP_NORM) {

					atfm_particle(pc, pc2, bag, params);
				} // if tangential velocity large enough
#elif defined LIN_TAN
				lin_tan(pc, pc2, bag, params);
#endif

			} // if processor owns collision
		} // if particles in contact
		else {

#ifdef STORE_COLLISION
			// Remove any existing particle collisions (no longer in contact)
			pc  = Collision_pp_find(p, p2);
			pc2 = Collision_pp_find(p2, p);
			if(pc  != NULL && pc->state != COLL_STATE_OWNER) {
				pc->state = COLL_STATE_DESTROY;
			}
			if(pc2 != NULL && pc2->state != COLL_STATE_OWNER) {
				pc2->state = COLL_STATE_DESTROY;
			}
#endif

#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL || defined ELECTROSTATIC_REPULSION || defined COHESION
			// Check if particles are close enough to interact
			cutoff_dist = p->R + p2->R + interaction_range;
			cutoff_dist_sq = cutoff_dist * cutoff_dist;
			if(center_distance_sq <= cutoff_dist_sq) {

				center_distance = sqrt(center_distance_sq);

				// Point in the middle between particle surfaces
				temp = 0.5 * (center_distance + p->R - p2->R) / center_distance;
				FORI3 midpoint[i] = p->X[i] + temp * (p2->X[i] - p->X[i]);

				// Find out if this processor is supposed to own this collision
				if (Collision_check_ownership(midpoint, stage, grid)) {

					// Don't adjust surface distance for roughness
					surface_distance = center_distance - (p->R + p2->R);
					Collision_fill_bag_particle(bag, center_distance, surface_distance, params);

	#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL
					//----------------------------------------------------------
					// Lubrication force
					//----------------------------------------------------------
					if (surface_distance <= params->lub_range * h) {
						lubrication(bag, h, params);
					}
	#endif
	#ifdef COHESION
					if (surface_distance <= params->coh_range * h) {
					  	cohesion(params->coh_range * h, bag, params);
					}
	#endif
	#ifdef ELECTROSTATIC_REPULSION
					//----------------------------------------------------------
					// Electrostatic repulsion force
					//----------------------------------------------------------
					if (surface_distance <= params->erm_range * h) {
						erm(bag, params);
					}
	#endif
				}  // if processor owns collision
			}  // if within interaction range
#endif  // Long range interaction models

		} // if particles are not in contact

		p2 = p2 -> next;
	}

	free(bag);
}




/******************************************************************************/
/*
 * Calculates collision forces and torques between particle 'p' and walls and
 * adds them to 'p->Fc' and 'p->Tc', respectively.
 *
 * Assumes uniform square grid
 */
/******************************************************************************/
void Collision_wall(Particle *p, MAC_grid *grid, Parameters *params) {

#ifdef STORE_COLLISION
	Collision *pc;
#endif

	Collision_bag *bag = (Collision_bag *)malloc(sizeof(Collision_bag));
	bag -> p = p;
	bag -> p2 = NULL;
	bag -> stage = COLL_STAGE_WALL;

	// Distance between particle centers
	double center_distance;

	// Distance between particle surfaces, negative value indicates overlap
	double surface_distance;

	// Surface distance at which particles begin to interact
	double interaction_range = 0.0;

	// Dimension and side of wall we are referring to, also identified by wall_ID
	int dim, side, wall_ID;

	// Grid spacing (assuming uniform grid)
	double h = grid->dx_u[1];

	// Flag for specifying if collision was just created
	int new_collision;

	// Distance at which collision takes place
	double collision_offset = 0.0;
#ifdef ROUGH_COLLISION
	collision_offset = params->roughness * p->R;
#endif

#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL
	interaction_range = max(interaction_range, params->lub_range * h);
#endif
#ifdef ELECTROSTATIC_REPULSION
	interaction_range = max(interaction_range, params->erm_range * h);
#endif
#ifdef COHESION
	interaction_range = max(interaction_range, params->coh_range * h);
#endif
	/*------------------------------------------------------------------------*/
	// check walls
	/*------------------------------------------------------------------------*/
#if defined XPERIODIC || defined RIGHT_OUTFLOW
	int dim_start = 1;
#else
	int dim_start = 0;
#endif
#ifdef ZPERIODIC
	int dim_end = 1;
#else
	int dim_end = 2;
#endif
int count =0;
#ifdef YPERIODIC
	count=1;
#endif

	// 'dim': {x, y, z} = {0, 1, 2}
	for (dim = dim_start; dim <= dim_end; dim++) {
		if (dim ==1 && count ==1){
  		  continue;
  	  }
		// 'side': {lower wall, upper wall} = {-1, 1}
		for (side = -1; side <= 1; side += 2) {

			wall_ID = Collision_get_wall_ID(dim, side);
			center_distance = Collision_pw_distance(wall_ID, dim, p->X, params);
			surface_distance = center_distance - p->R;

			if (surface_distance <= collision_offset) {

				// Reorient overlap distance to account for offset
				surface_distance -= collision_offset;
				Collision_fill_bag_wall(bag, dim, side, center_distance, surface_distance, params);

#ifdef STORE_COLLISION
				// Look for existing wall collision
				pc = Collision_pw_find(p, wall_ID);

				// If pc does not exist, create, calculate and link it
				if (pc == NULL) {
					// Make sure particle is moving towards wall.  The alternate
					// case can occur in rare cases when gravity is involved
					if (bag->g_dot_n > 0) {
						new_collision = 1;
						pc = Collision_pw_new(wall_ID, bag, params);
					}
					else {
						continue;
					}
				}
				else {
					new_collision = 0;
				}
				if (pc != NULL) pc->state = COLL_STATE_OWNER;
#endif
				//----------------------------------------------------------
				// Normal collision force
				//----------------------------------------------------------
				// '-side' is wall normal direction
#ifdef ACTM
				actm(pc, bag, params);
#else
				dem(pc, bag, params);
#endif

				//----------------------------------------------------------
				// Tangential collision force
				//----------------------------------------------------------
#ifdef ATFM
				if (bag->gt_cp_norm >= MIN_GT_CP_NORM) {

					atfm_wall(pc, bag, params);
				} // if tangential velocity large enough
#elif defined LIN_TAN
				lin_tan(pc, NULL, bag, params);
#endif

			} // if particles in contact
			else {

#ifdef STORE_COLLISION
				// Remove any existing wall collisions (no longer in contact)
				pc = Collision_pw_find(p, wall_ID);
				if (pc != NULL) pc->state = COLL_STATE_DESTROY;
#endif

#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL || defined ELECTROSTATIC_REPULSION || defined COHESION
				// Check if particles are close enough to interact
				if (surface_distance <= interaction_range) {
					Collision_fill_bag_wall(bag, dim, side, center_distance, surface_distance, params);

	#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL
					//----------------------------------------------------------
					// Lubrication force
					//----------------------------------------------------------
					if (surface_distance <= params->lub_range * h) {
						lubrication(bag, h, params);
					}
	#endif
	#ifdef ELECTROSTATIC_REPULSION
					//----------------------------------------------------------
					// Electrostatic repulsion force
					//----------------------------------------------------------
					if (surface_distance <= params->erm_range * h) {
						erm(bag, params);
					}
	#endif
	#ifdef COHESION
					if(surface_distance <= params->coh_range * h) {
					  	cohesion(params-> coh_range * h, bag, params);
					}
	#endif
				}  // if within interaction range
#endif  // Long range interaction models

			}  // if particles are not in contact
		}  // for side
	}  // for dim

	free(bag);
}




/******************************************************************************/
/*
 * Returns the wall ID for the given wall 'dim' and 'side' values
 */
/******************************************************************************/
int Collision_get_wall_ID(int dim, int side) {

	if (dim==0 && side==-1) return WALL_ID_XMIN;
	if (dim==0 && side== 1) return WALL_ID_XMAX;
	if (dim==1 && side==-1) return WALL_ID_YMIN;
	if (dim==1 && side== 1) return WALL_ID_YMAX;
	if (dim==2 && side==-1) return WALL_ID_ZMIN;
	if (dim==2 && side== 1) return WALL_ID_ZMAX;

}




/******************************************************************************/
/*
 * Finds the minimum particle-wall distance in the dimension specified by 'dim'
 * and the direction specified by 'side'
 *
 * Returns a positive value
 *
 * Assumes particles are never in contact with two opposing walls at the same
 * time
 */
/******************************************************************************/
double Collision_pw_distance(int wall_ID, int dim, double *X, Parameters *params) {
	// find wall coordinate
	// ugly
	double edge;
	if (wall_ID == WALL_ID_XMIN) edge = params->xmin;
	if (wall_ID == WALL_ID_XMAX) edge = params->xmax;
	if (wall_ID == WALL_ID_YMIN) edge = params->ymin;
#ifdef DOWNWARD_MOVING_WALL
	if (wall_ID == WALL_ID_YMAX) {
		edge = (params->vel_init_y0 - params->ymax) * params->time
		       / params->time_max + params->ymax;
	}
#else
	if (wall_ID == WALL_ID_YMAX) edge = params->ymax;
#endif
	if (wall_ID == WALL_ID_YMAX) edge = params->ymax;
	if (wall_ID == WALL_ID_ZMIN) edge = params->zmin;
	if (wall_ID == WALL_ID_ZMAX) edge = params->zmax;

	return fabs(edge - X[dim]);
}




/******************************************************************************/
/*
 * See if current processor owns the collision, for the case of local-foreign
 * particle collisions.  It does this by checking ownership of the location of
 * the midpiont of the collision.
 *
 * Returns '1' if it owns the collision, '0' if it does not.
 */
/******************************************************************************/
int Collision_check_ownership(double *midpoint, int stage, MAC_grid *grid) {

	// We own the collision for mobile-mobile collisions
	if (stage == COLL_STAGE_MOBILE) {
		return 1;
	}

	// Grid bounds of local subdomain
	double G_xmin = grid -> xu[grid -> G_Is];
	double G_xmax = grid -> xu[min(grid -> G_Ie, grid -> NX - 1)];
	double G_ymin = grid -> yv[grid -> G_Js];
	double G_ymax = grid -> yv[min(grid -> G_Je, grid -> NY - 1)];
	double G_zmin = grid -> zw[grid -> G_Ks];
	double G_zmax = grid -> zw[min(grid -> G_Ke, grid -> NZ - 1)];

	// Check if midpoint is out of bounds
	if( midpoint[0] < G_xmin || midpoint[0] >= G_xmax ||
	    midpoint[1] < G_ymin || midpoint[1] >= G_ymax ||
	    midpoint[2] < G_zmin || midpoint[2] >= G_zmax ) {

		return 0;
	}
	else {
		return 1;
	}
}




/******************************************************************************/
/*
 * Find an existing particle collision for particle 'p' that involves 'p2'.
 * Return 'NULL' if no such collision found.
 */
/******************************************************************************/
Collision *Collision_pp_find(Particle *p, Particle *p2) {

	Collision *pc;

	pc = p->particle_collision;
	while (pc != NULL && pc->other_ID != p2->ID) {
		pc = pc->next;  // Keep looking until end or we find it
	}
	// Return either the discovered collision, or the NULL terminating value
	return pc;
}




/******************************************************************************/
/*
 * Find an existing wall collision for particle 'p' that involves the wall
 * characterized by 'wall_ID'.  Return 'NULL' if no such collision found.
 */
/******************************************************************************/
Collision *Collision_pw_find(Particle *p, int wall_ID) {

	Collision *pc;

	pc = p->wall_collision;
	while (pc != NULL && (pc->other_ID != wall_ID)) {
		pc = pc->next;  // Keep looking until end or we find it
	}
	// Return either the discovered collision, or the NULL terminating value
	return pc;
}




/******************************************************************************/
/*
 * Creates a new particle-particle collision to be stored.
 *
 * Collision structures are created for both particles 'p' and 'p2' and added to
 * the beginning of their 'paricle_collision' linked lists.
 *
 * Spring and damper coefficients are determined based on the normal collision
 * model.
 *
 * Tangential collision categorization (sliding vs. rolling) also takes place,
 * depending on the tangential collision model.
 */
/******************************************************************************/
void Collision_pp_new(Collision **pc_ptr, Collision **pc2_ptr,
		Collision_bag *bag, Parameters *params) {
	double u_in, m;

	Particle *p  = bag -> p;
	Particle *p2 = bag -> p2;

	Collision *pc  = (Collision *) calloc(1, sizeof(Collision));
	Collision *pc2 = (Collision *) calloc(1, sizeof(Collision));
	*pc_ptr  = pc;
	*pc2_ptr = pc2;

	// Calculate collision parameters.  These functions should only modify 'pc'
#ifdef ACTM
	actm_new_pp_collision(&pc->kn, &pc->dn, bag, params);
#elif defined DEM
	dem_new_pp_collision(&pc->kn, &pc->dn, bag, params);
#else
	pc -> kn = params -> kn;
	pc -> dn = params -> dn;
#endif

	// Compute Stokes number to check for 'dry' or 'wet' collision
	pc -> St = 2.0 / 9.0 * p->rho_s * bag -> g_dot_n * p->R * params->Re;

#ifdef ENABLE_ATFM_ROLLING
	if (pc-> St > ST_CRIT){
		u_in = bag -> g_dot_n;
	}
	else{
		if (bag->stage == COLL_STAGE_FIXED) {
			u_in = 4.5 * ST_CRIT / (params->Re * p->rho_s * p->R);
		}
		else {
			m = (p->M * p2->M) / (p->M + p2->M);
			u_in = 4.5 * ST_CRIT / (params->Re * p->rho_s * p->R);
			u_in = max(u_in, 4.5 * ST_CRIT / (params->Re * p2->rho_s * p2->R));
		}
	}
	pc -> type  = atfm_categorize(bag, params, u_in);
#else
	pc -> type = SLIDING;
#endif
#ifdef LIN_TAN
	lin_tan_new_collision(&pc->kt, &pc->dt, bag, params);
#else
	pc->kt = 0.0;
	pc->dt = 0.0;
#endif

	// Initialize tangential displacement to zero
	DSET_ZERO(pc ->zeta_t, 3);
	DSET_ZERO(pc2->zeta_t, 3);
	DSET_ZERO(pc ->zeta_t_old, 3);
	DSET_ZERO(pc2->zeta_t_old, 3);

	// Initialize state
	pc -> state = COLL_STATE_NEUTRAL;

	// Copy collision parameters to 'pc2'
	*pc2 = *pc;

	// Set unique particle ID
	pc  -> other_ID = p2 -> ID;
	pc2 -> other_ID = p -> ID;

	// Add collision to beginning of linked list
	pc  -> next = p  -> particle_collision;
	pc2 -> next = p2 -> particle_collision;
	p  -> particle_collision = pc;
	p2 -> particle_collision = pc2;
}



/******************************************************************************/
/*
 * Creates a new particle-wall collision to be stored.
 *
 * Collision structures are created for particle 'p' and added to the beginning
 * of its 'paricle_collision' linked list.
 *
 * Spring and damper coefficients are determined based on the normal collision
 * model.
 *
 * Tangential collision categorization (sliding vs. rolling) also takes place,
 * depending on the tangential collision model.
 */
/******************************************************************************/
Collision *Collision_pw_new(int wall_ID, Collision_bag *bag, Parameters *params) {
        double u_in;
	Particle *p = bag -> p;

	Collision *pc = (Collision *) calloc(1, sizeof(Collision));

	// Calculate collision parameters.  These functions should only modify 'pc'
#ifdef ACTM
	actm_new_pw_collision(&pc->kn, &pc->dn, bag, params);
#elif defined DEM
	dem_new_pw_collision(&pc->kn, &pc->dn, bag, params);
#else
	pc -> kn = params -> kn;
	pc -> dn = params -> dn;
#endif

	// Compute Stokes number to check for 'dry' or 'wet' collision
	pc -> St = 2.0 / 9.0 * p->rho_s * bag -> g_dot_n * p->R * params->Re;

#ifdef ENABLE_ATFM_ROLLING
	if (pc-> St > ST_CRIT){
		u_in = bag -> g_dot_n;
	}
	else{
		u_in = 4.5 * ST_CRIT / (params->Re * p->rho_s * p->R);
	}
	pc -> type  = atfm_categorize(bag, params, u_in);
#else
	pc -> type = SLIDING;
#endif
#ifdef LIN_TAN
	lin_tan_new_collision(&pc->kt, &pc->dt, bag, params);
#else
	pc->kt = 0.0;
	pc->dt = 0.0;
#endif

	// Initialize tangential displacement to zero
	DSET_ZERO(pc->zeta_t, 3);
	DSET_ZERO(pc->zeta_t_old, 3);

	// Initialize state
	pc -> state = COLL_STATE_NEUTRAL;

	// Set unique ID to wall ID
	pc -> other_ID = wall_ID;

	// Add collision to beginning of linked list
	pc -> next = p -> wall_collision;
	p -> wall_collision = pc;

	return pc;
}




/******************************************************************************/
/*
 * Fills collision bag with derived variables for particle collision
 */
/******************************************************************************/
void Collision_fill_bag_particle(Collision_bag *bag, double center_distance,
		double surface_distance, Parameters *params) {

	int i;

	// Normal unit vector from center of 'p' to center of 'p2'
	double *n = bag -> n;

	// Tangential unit vector, relative surface velocity between 'p' and 'p2'
	double *t = bag -> t;

	// Relative translational velocities of two particles
	double *g  = bag -> g;   // Total velocity
	double *gn = bag -> gn;  // Component normal to contact surface
	double *gt = bag -> gt;  // Component tangential to contact surface

	// Relative rotational velocities of two particles
	double *Om_cross_R = bag -> Om_cross_R;

	// Relative velocity of contact point between particles
	double *gt_cp = bag -> gt_cp;

	Particle *p = bag -> p;
	Particle *p2 = bag -> p2;

	double R_cp, R2_cp;
	double *Omega = p -> Omega;
	double *Omega2 = p2 -> Omega;

#ifdef ATFM
	double g_old[3];
	double gt_old[3];
	double gt_cp_old[3];
	double *Omega_old = p -> Omega_old;
	double *Omega2_old = p2 -> Omega_old;
#endif

	bag -> surface_distance = surface_distance;

	// Calculate relative translational velocity
	FORI3 g[i] = p->U[i] - p2->U[i];
#ifdef ATFM
	FORI3 g_old[i] = p->U_old[i] - p2->U_old[i];
#endif

	// Calculate normal vector from center of 'p' to center of 'p2'
	FORI3 n[i] = (p2->X[i] - p->X[i]) / center_distance;

	// Calculate relative normal input velocity in direction 'n'
	bag -> g_dot_n = DOT(n, g);

	// Normal component of relative translational velocity
	FORI3 gn[i] = bag->g_dot_n * n[i];

	// Tangential component of relative translational velocity
	FORI3 gt[i] = g[i] - gn[i];
#ifdef ATFM
	FORI3 gt_old[i] = g_old[i] - DOT(n, g_old);
#endif

	// Distance (radius) of contact point from particles 'p' and 'p2'
//	R_cp  = temp * center_distance;
//	R2_cp = center_distance - R_cp;
	R_cp  = p->R;
	R2_cp = p2->R;
	bag -> R_cp  = R_cp;
	bag -> R2_cp = R2_cp;

	Om_cross_R[0] = R_cp  * (Omega[1]  * n[2] - Omega[2]  * n[1])
	              + R2_cp * (Omega2[1] * n[2] - Omega2[2] * n[1]);
	Om_cross_R[1] = R_cp  * (Omega[2]  * n[0] - Omega[0]  * n[2])
	              + R2_cp * (Omega2[2] * n[0] - Omega2[0] * n[2]);
	Om_cross_R[2] = R_cp  * (Omega[0]  * n[1] - Omega[1]  * n[0])
	              + R2_cp * (Omega2[0] * n[1] - Omega2[1] * n[0]);

	// Tangential component of relative velocity of contact point
	gt_cp[0] =  gt[0] + Om_cross_R[0];
	gt_cp[1] =  gt[1] + Om_cross_R[1];
	gt_cp[2] =  gt[2] + Om_cross_R[2];
#ifdef ATFM
	gt_cp_old[0] =  gt_old[0] + R_cp  * (Omega_old[1]  * n[2] - Omega_old[2]  * n[1])
	                          + R2_cp * (Omega2_old[1] * n[2] - Omega2_old[2] * n[1]);
	gt_cp_old[1] =  gt_old[1] + R_cp  * (Omega_old[2]  * n[0] - Omega_old[0]  * n[2])
	                          + R2_cp * (Omega2_old[2] * n[0] - Omega2_old[0] * n[2]);
	gt_cp_old[2] =  gt_old[2] + R_cp  * (Omega_old[0]  * n[1] - Omega_old[1]  * n[0])
	                          + R2_cp * (Omega2_old[0] * n[1] - Omega2_old[1] * n[0]);
#endif

	// Magnitude of contact point velocity
	bag->gt_cp_norm = sqrt(DOT(gt_cp, gt_cp));
#ifdef ATFM
	bag->gt_cp_norm_old = sqrt(DOT(gt_cp_old, gt_cp_old));
#endif

	// Tangential unit normal
	FORI3 t[i] = gt_cp[i] / bag->gt_cp_norm;
}




/******************************************************************************/
/*
 * Fills collision bag with derived variables for wall collision
 */
/******************************************************************************/
void Collision_fill_bag_wall(Collision_bag *bag, int dim, int side,
		double center_distance, double surface_distance, Parameters *params) {

	int i;

	double *n  = bag -> n;
	double *t  = bag -> t;
	double *g  = bag -> g;
	double *gn = bag -> gn;
	double *gt = bag -> gt;
	double *Om_cross_R = bag -> Om_cross_R;
	double *gt_cp = bag -> gt_cp;

	Particle *p = bag -> p;

	double R_cp;
	double *Omega = p -> Omega;

#ifdef ATFM
	double g_old[3];
	double gt_old[3];
	double gt_cp_old[3];
	double *Omega_old = p -> Omega_old;
#endif

	bag -> surface_distance = surface_distance;
	int wall_ID = Collision_get_wall_ID(dim, side);

	// Wall velocity, if enabled
#if defined BOTTOM_WALL_VELOCITY && defined TOP_WALL_VELOCITY
	double wall_vel = params -> ubulk_target;
#else
	double wall_vel = 2.0 * params -> ubulk_target;
#endif

	// Normal vector from particle to wall
	FORI3 n[i] = 0.0;
	n[dim] = side;

	// Calculate relative translational velocity
	FORI3 g[i] = p->U[i];
#ifdef DOWNWARD_MOVING_WALL
	if (wall_ID == WALL_ID_YMAX) {
		g[1] = p->U[1] - (params->vel_init_y0 - params->ymax) / params->time_max;
	}
#endif
#ifdef ATFM
	FORI3 g_old[i] = p->U_old[i];
#endif

	// Normal component of relative translational velocity
	FORI3 gn[i] = 0.0;
	gn[dim] = g[dim];

	// Calculate relative normal input velocity in direction 'n'
	bag->g_dot_n = n[dim] * gn[dim];

	// Tangential component of relative translational velocity
	FORI3 gt[i] = g[i];
#ifdef ATFM
	FORI3 gt_old[i] = g_old[i];
#endif
	gt[dim] = 0;
#ifdef ATFM
	gt_old[dim] = 0;
#endif
#ifdef BOTTOM_WALL_VELOCITY
	if (wall_ID == WALL_ID_YMIN) {
		gt[0] += wall_vel;
	}
#endif
#ifdef TOP_WALL_VELOCITY
	if (wall_ID == WALL_ID_YMAX) {
		gt[0] -= wall_vel;
	}
#endif

	// Distance (radius) of contact point from particle 'p'
//	R_cp = center_distance;
	R_cp = p->R;
	bag -> R_cp = R_cp;

	// Relative rotational velocity of particle and wall
	Om_cross_R[0] = R_cp * (Omega[1] * n[2] - Omega[2] * n[1]);
	Om_cross_R[1] = R_cp * (Omega[2] * n[0] - Omega[0] * n[2]);
	Om_cross_R[2] = R_cp * (Omega[0] * n[1] - Omega[1] * n[0]);

	// Freeslip boundary condition -> zero tangential motion
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
	if (wall_ID == WALL_ID_YMIN) {
		FORI3 gt[i] = 0.0;
		FORI3 Om_cross_R[i] = 0.0;
	}
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
	if (wall_ID == WALL_ID_YMAX) {
		FORI3 gt[i] = 0.0;
		FORI3 Om_cross_R[i] = 0.0;
	}
#endif

	// Tangential component of relative velocity of contact point
	gt_cp[0] = gt[0] + Om_cross_R[0];
	gt_cp[1] = gt[1] + Om_cross_R[1];
	gt_cp[2] = gt[2] + Om_cross_R[2];
#ifdef ATFM
	gt_cp_old[0] = gt_old[0] + R_cp * (Omega_old[1] * n[2] - Omega_old[2] * n[1]);
	gt_cp_old[1] = gt_old[1] + R_cp * (Omega_old[2] * n[0] - Omega_old[0] * n[2]);
	gt_cp_old[2] = gt_old[2] + R_cp * (Omega_old[0] * n[1] - Omega_old[1] * n[0]);
#endif


	// Magnitude of contact point velocity
	bag->gt_cp_norm = sqrt(DOT(gt_cp, gt_cp));
#ifdef ATFM
	bag->gt_cp_norm_old = sqrt(DOT(gt_cp_old, gt_cp_old));
#endif

	// Tangential unit normal
	FORI3 t[i] = gt_cp[i] / bag->gt_cp_norm;
}
