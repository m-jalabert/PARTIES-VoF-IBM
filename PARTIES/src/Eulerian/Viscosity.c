#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Memory.h"
#include "Communication.h"
#include "Viscosity.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/******************************************************************************/
/*
 */
/******************************************************************************/
Viscosity *Viscosity_create(MAC_grid *grid, Parameters *params) {
	
	Viscosity *new_viscosity = (Viscosity *)malloc(sizeof(Viscosity));
	Memory_check_allocation(new_viscosity);
	
	new_viscosity -> nu  = Memory_allocate_flow_variable(grid, params);
	new_viscosity -> nuX = Memory_allocate_flow_variable(grid, params);
	new_viscosity -> nuY = Memory_allocate_flow_variable(grid, params);
	new_viscosity -> nuZ = Memory_allocate_flow_variable(grid, params);
	
	return new_viscosity;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Viscosity_destroy(Viscosity *viscosity, MAC_grid *grid, Parameters *params) {
	
	Memory_free_flow_variable(grid, params, viscosity -> nu);
	Memory_free_flow_variable(grid, params, viscosity -> nuX);
	Memory_free_flow_variable(grid, params, viscosity -> nuY);
	Memory_free_flow_variable(grid, params, viscosity -> nuZ);
	
	free(viscosity);
}




/******************************************************************************/
/*
 Interpolate cell-edge viscosities 'nuX', 'nuY', 'nuZ' from cell-centered 
 viscosity 'nu'
 
     nu[k][j][i]  = viscosity at (xc[i], yc[j], zc[k])
     nuX[k][j][i] = viscosity at (xc[i], yv[j], zw[k])
     nuY[k][j][i] = viscosity at (xu[i], yc[j], zw[k])
     nuZ[k][j][i] = viscosity at (xu[i], yv[j], zc[k])
 
 If using turbulence modeling, also evaluates viscosity 'nu' from subgrid 
 viscosity 'nut'.
 
     NOTE: 'nut' should have its boundaries and ghost cells updated before 
     calling this function.
 */
/******************************************************************************/
void Viscosity_set_cell_edges(Cart3d_bag *data_bag) {
	
	int i, j, k;
	
	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Viscosity *viscosity = data_bag -> viscosity;
	
	double ***nu  = viscosity -> nu;
	double ***nuX = viscosity -> nuX;
	double ***nuY = viscosity -> nuY;
	double ***nuZ = viscosity -> nuZ;
	
	// Processor boundaries, including ghost nodes
	int Is_g = grid -> L_Is;
	int Js_g = grid -> L_Js;
	int Ks_g = grid -> L_Ks;
	
	int Ie_g = grid -> L_Ie;
	int Je_g = grid -> L_Je;
	int Ke_g = grid -> L_Ke;
	
#if defined LES || defined RANS
	
	double iRe = 1.0 / params -> Re;
	
	#ifdef LES
	double ***nut = data_bag -> smag -> nut;
	#elif defined RANS
	double ***nut = data_bag -> rans -> nut;
	#endif
	
	//--------------------------------------------------------------------------
	// Populate 'nu' with eddy viscosity plus viscosity, including ghost nodes
	//--------------------------------------------------------------------------
	for (k = Ks_g; k < Ke_g; k++) {
		for (j = Js_g; j < Je_g; j++) {
			for (i = Is_g; i < Ie_g; i++){
				nu[k][j][i] = iRe + nut[k][j][i];
			}
		}
	}
	
#endif // LES or RANS
	
	//--------------------------------------------------------------------------
	// Interpolate cell edge viscosities from cell-centered 'nu' using first
	// shell of ghost cells
	//--------------------------------------------------------------------------
	for (k = Ks_g+1; k < Ke_g; k++) {
		for (j = Js_g+1; j < Je_g; j++) {
			for (i = Is_g+1; i < Ie_g; i++){
				
				nuX[k][j][i] = 0.25 * ( nu[k][j][i]   + nu[k-1][j][i]
				                      + nu[k][j-1][i] + nu[k-1][j-1][i] );
				
				nuY[k][j][i] = 0.25 * ( nu[k][j][i]   + nu[k][j][i-1]
				                      + nu[k-1][j][i] + nu[k-1][j][i-1] );
				
				nuZ[k][j][i] = 0.25 * ( nu[k][j][i]   + nu[k][j][i-1]
				                      + nu[k][j-1][i] + nu[k][j-1][i-1] );
				
			}
		}
	}
}




/******************************************************************************/
/*
 Sets value of nut at the domain boundaries according to the velocity boundary
 conditions and then updates all ghost nodes
 */
/******************************************************************************/
void Viscosity_update_boundaries(double ***nut, Cart3d_bag *data_bag) {
	
	int i, j, k;
	
	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;
	
	// Start index of bottom-left-back corner on current processor
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;
	
	// End index of top-right-front corner on current processor
	int Ie = min(NX-1, grid -> G_Ie);
	int Je = min(NY-1, grid -> G_Je);
	int Ke = min(NZ-1, grid -> G_Ke);
	
	/*------------------------------------------------------------------------*/
	/*
	 Update ghost nodes.  This takes care of periodic boundaries
	 */
	/*------------------------------------------------------------------------*/
	
#ifndef XPERIODIC
	//--------------------------------------------------------------------------
	// Left wall
	//--------------------------------------------------------------------------
	if (Is == 0) {
		i = 0;
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
	#ifdef LEFT_WALL_VELOCITY_NOSLIP
				nut[k][j][i-1] = -nut[k][j][i];
	#else
				nut[k][j][i-1] = nut[k][j][i];
	#endif
			}
		}
		
		// Now include first ghost cells of left wall
		Is = -1;
	}
	
	//--------------------------------------------------------------------------
	// Right wall
	//--------------------------------------------------------------------------
	if (Ie == NX-1) {
		i = NX - 1;
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
	#ifdef RIGHT_WALL_VELOCITY_NOSLIP
				nut[k][j][i] = -nut[k][j][i-1];
	#else
				nut[k][j][i] = nut[k][j][i-1];
	#endif
			}
		}
		
		// Now include first ghost cells of right wall
		Ie = NX;
	}
#endif // not XPERIODIC
	
	
	//--------------------------------------------------------------------------
	// Bottom wall
	//--------------------------------------------------------------------------
	if (Js == 0) {
		j = 0;
		for (k = Ks; k < Ke; k++) {
			for (i = Is; i < Ie; i++) {
#if defined BOTTOM_WALL_VELOCITY_NOSLIP || defined BOTTOM_WALL_SCHUMANN
				if (i == -1 || i == NX-1) {
					// Copy values into corner cells
					nut[k][j-1][i] = nut[k][j][i];
				}
				else {
					// Negate values in boundary cells
					nut[k][j-1][i] = -nut[k][j][i];
				}
#else
				nut[k][j-1][i] = nut[k][j][i];
#endif
			}
		}
		
		// Now include first ghost cells of bottom wall
		Js = -1;
	}
	
	//--------------------------------------------------------------------------
	// Top wall
	//--------------------------------------------------------------------------
	if (Je == NY-1) {
		j = NY - 1;
		for (k = Ks; k < Ke; k++) {
			for (i = Is; i < Ie; i++) {
#ifdef TOP_WALL_VELOCITY_NOSLIP
				if (i == -1 || i == NX-1) {
					// Copy values into corner cells
					nut[k][j][i] = nut[k][j-1][i];
				}
				else {
					// Negate values in boundary cells
					nut[k][j][i] = -nut[k][j-1][i];
				}
#else
				nut[k][j][i] = nut[k][j-1][i];
#endif
			}
		}
		
		// Now include first ghost cells of top wall
		Je = NY;
	}
	
	
#ifndef ZPERIODIC
	//--------------------------------------------------------------------------
	// Back wall
	//--------------------------------------------------------------------------
	if (Ks == 0) {
		k = 0;
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
	#ifdef BACK_WALL_VELOCITY_NOSLIP
				if (i == -1 || i == NX-1 || j == -1 || j == NY-1) {
					// Copy values into corner cells
					nut[k-1][j][i] = nut[k][j][i];
				}
				else {
					// Negate values in boundary cells
					nut[k-1][j][i] = -nut[k][j][i];
				}
	#else
				nut[k-1][j][i] = nut[k][j][i];
	#endif
			}
		}
	}
	
	//--------------------------------------------------------------------------
	// Front wall
	//--------------------------------------------------------------------------
	if (Ke == NZ-1) {
		k = NZ - 1;
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
	#ifdef FRONT_WALL_VELOCITY_NOSLIP
				if (i == -1 || i == NX-1 || j == -1 || j == NY-1) {
					// Copy values into corner cells
					nut[k][j][i] = nut[k-1][j][i];
				}
				else {
					// Negate values in boundary cells
					nut[k][j][i] = -nut[k-1][j][i];
				}
	#else
				nut[k][j][i] = nut[k-1][j][i];
	#endif
			}
		}
	}
#endif // not ZPERIODIC
	
	
	/*------------------------------------------------------------------------*/
	/*
	 Update ghost nodes.  This takes care of periodic boundaries
	 */
	/*------------------------------------------------------------------------*/
	Communication_update_ghost_nodes_flow_variable(nut, 'c', params->ghost_nodes, data_bag);
}


