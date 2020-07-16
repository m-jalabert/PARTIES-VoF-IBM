#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Array.h"
#include "Cart3d.h"
#include "Communication.h"
#include "Conc.h"
#include "Display.h"
#include "Grid.h"
#include "Immersed.h"
#include "Memory.h"
#include "MyMath.h"
#include "Vorticity.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void Vorticity_create(Cart3d_bag *data_bag){

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Vorticity *new_vor;
	new_vor = (Vorticity *)malloc(sizeof(Vorticity));
	Memory_check_allocation(new_vor);

	/*------------------------------------------------------------------------*/
	/*
	     - vor_*    : value of the * component of vorticity
	     - vor_*_bc : value of the * component of vorticity at the cell_center
	 */
	/*------------------------------------------------------------------------*/


	new_vor->vor_x = Memory_allocate_flow_variable(grid, params);
	new_vor->vor_y = Memory_allocate_flow_variable(grid, params);
	new_vor->vor_z = Memory_allocate_flow_variable(grid, params);

	new_vor->vor_x_bc = Memory_allocate_flow_variable(grid, params);
	new_vor->vor_y_bc = Memory_allocate_flow_variable(grid, params);
	new_vor->vor_z_bc = Memory_allocate_flow_variable(grid, params);

	data_bag->vor = new_vor;

	return;
}


void Vorticity_compute(Cart3d_bag *data_bag) {

	int i,j,k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	// Grid dimensions
	double *xu = grid->xu;
	double *yv = grid->yv;
	double *zw = grid->zw;

	double *xc = grid->xc;
	double *yc = grid->yc;
	double *zc = grid->zc;

	double ***u = data_bag->u->data;
	double ***v = data_bag->v->data;
	double ***w = data_bag->w->data;

	double ***vor_x = data_bag->vor->vor_x;
	double ***vor_y = data_bag->vor->vor_y;
	double ***vor_z = data_bag->vor->vor_z;

	double ***vor_x_bc = data_bag->vor->vor_x_bc;
	double ***vor_y_bc = data_bag->vor->vor_y_bc;
	double ***vor_z_bc = data_bag->vor->vor_z_bc;

//	Velocity_update_boundaries(u, 'u', VEL_TYPE_NORMAL, data_bag);
//	Velocity_update_boundaries(v, 'v', VEL_TYPE_NORMAL, data_bag);
//	Velocity_update_boundaries(w, 'w', VEL_TYPE_NORMAL, data_bag);

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	//--------------------------------------------------------------------------
	// Values for cell-cetentered data
	//--------------------------------------------------------------------------
	// Indices start and end on current processor
	int i_start = max(1,Is);
	int j_start = Js;
	int k_start = max(1,Ks);

	int i_end = Ie+1;
	int j_end = Je+1;
	int k_end = Ke+1;

#ifdef XPERIODIC
	if (i_end == NX-1)
		i_end = NX;
#endif

#ifdef ZPERIODIC
	if (k_end == NZ-1)
		k_end = NZ;
#endif

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++) {
				vor_x[k][j][i] = (w[k][j][i]- w[k][j-1][i])/(yc[j]-yc[j-1]) - (v[k][j][i]- v[k-1][j][i])/(zc[k]-zc[k-1]);
				vor_y[k][j][i] = (u[k][j][i]- u[k-1][j][i])/(zc[k]-zc[k-1]) - (w[k][j][i]- w[k][j][i-1])/(xc[i]-xc[i-1]);
				vor_z[k][j][i] = (v[k][j][i]- v[k][j][i-1])/(xc[i]-xc[i-1]) - (u[k][j][i]- u[k][j-1][i])/(yc[j]-yc[j-1]);
			}
		}
	}

	int i_end_bc;
	int j_end_bc;
	int k_end_bc;

	i_end_bc = Ie;
	j_end_bc = Je;
	k_end_bc = Ke;
/*
	i_end_bc = min(i_end,NX-1);;
	j_end_bc = min(j_end,NY-1);
	k_end_bc = min(k_end,NZ-1);
*/
	for (k=k_start; k<k_end_bc; k++) {
		for (j=j_start; j<j_end_bc; j++) {
			for (i=i_start; i<i_end_bc; i++) {
				vor_x_bc[k][j][i] = 0.25*(vor_x[k][j][i]+vor_x[k][j+1][i]+vor_x[k+1][j][i]+vor_x[k+1][j+1][i]);
				vor_y_bc[k][j][i] = 0.25*(vor_y[k][j][i]+vor_y[k+1][j][i]+vor_y[k][j][i+1]+vor_y[k+1][j][i+1]);
				vor_z_bc[k][j][i] = 0.25*(vor_z[k][j][i]+vor_z[k][j+1][i]+vor_z[k][j][i+1]+vor_z[k][j+1][i+1]);
			}
		}
	}


}

/******************************************************************************/
/*
 This function releases the allocated memory for velocity structure.
 */
/******************************************************************************/
void Vorticity_destroy(Cart3d_bag *data_bag) {

	 MAC_grid *grid = data_bag->grid;
	 Parameters *params = data_bag->params;
	 Vorticity *vor = data_bag->vor;



	// vorticity data

	Memory_free_flow_variable(grid, params, vor->vor_x);
	Memory_free_flow_variable(grid, params, vor->vor_y);
	Memory_free_flow_variable(grid, params, vor->vor_z);

	Memory_free_flow_variable(grid, params, vor->vor_x_bc);
	Memory_free_flow_variable(grid, params, vor->vor_y_bc);
	Memory_free_flow_variable(grid, params, vor->vor_z_bc);

	free(vor);

}
