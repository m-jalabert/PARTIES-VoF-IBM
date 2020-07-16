#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Cart3d.h"
#include "Communication.h"
#include "Conc.h"
#include "Grid.h"
#include "Immersed.h"
#include "Memory.h"
#include "MyMath.h"
#include "Subgrid.h"
#include "Velocity.h"
#include "Viscosity.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/******************************************************************************/
/*
 This function allocates enough memory for the concentration structure based on
 parameters defined in "*params"
 */
/******************************************************************************/
Subgrid *Subgrid_create(MAC_grid *grid, Parameters *params) {

	Subgrid *new_subgrid;
	int NX, NY, NZ, NT;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;

	double ***delta;
	double *dx, *dy, *dz;
	int Ghost_Nodes;

	FILE *fid;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	// indices start and end on current processor
	i_start = Is; // i=NX-1 are not included
	j_start = Js;
	k_start = Ks;

	i_end = Ie;
	j_end = Je;
	k_end = Ke;

	// exclude the half cell added
#ifndef XPERIODIC
	i_end = min(NX-1, Ie);
#endif
	j_end = min(NY-1, Je);
#ifndef ZPERIODIC
	k_end = min(NZ-1, Ke);
#endif

	new_subgrid = (Subgrid *)malloc(sizeof(Subgrid));
	Memory_check_allocation(new_subgrid);

	// Now, create the vectors that hold the local and global data for all the
	// 3D properties

	new_subgrid->nut = Memory_allocate_flow_variable(grid, params);
	new_subgrid->ng_lengthscale_sq = Memory_allocate_noghost_variable(grid, params);

	new_subgrid->st_rate = (Strain_rate *) malloc(sizeof(Strain_rate));
	new_subgrid->st_rate->strain = Memory_allocate_flow_variable(grid, params);

	new_subgrid->Cs = 0.08;
	new_subgrid->Sct = 0.6;

	dx = grid->dx_u;
	dy = grid->dy_v;
	dz = grid->dz_w;

	delta = new_subgrid->ng_lengthscale_sq;
	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				delta[k][j][i] = pow(dx[i]*dy[j]*dz[k],2./3.);
			}
		}
	}

	return new_subgrid;
}




/******************************************************************************/
/*
 Calculate the derivative of concentration
     idir = 1 x derivative
     idir = 2 y derivative
     idir = 3 z derivative
 */
/******************************************************************************/
void Subgrid_concentration_derivative(Concentration *conc, MAC_grid *grid,
		Parameters *params, Subgrid *smag, int idir) {

	int NX, NY, NZ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	int i, j, k;
	double *i2dx_c, *i2dy_c, *i2dz_c;
	double ***c, ***dcdn;

/*	char bin_filename[50];
	FILE *fid, *fid1;

	sprintf(bin_filename,"straininfo%d.dat",params->rank);
	fid = fopen(bin_filename,"w");
	sprintf(bin_filename,"strainmag%d.dat",params->rank);
	fid1 = fopen(bin_filename,"w");
*/

	// Same for all quantities
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	// indices start and end on current processor
	i_start = max(1,Is); // i=0, i=NX-1 are not included
	j_start = max(1,Js);
	k_start = max(1,Ks);

#ifdef XPERIODIC
	i_start = Is;
#endif
#ifdef ZPERIODIC
	k_start = Ks;
#endif

	// exclude the half cell added
	i_end   = min(NX-1, Ie);
	j_end   = min(NY-1, Je);
	k_end   = min(NZ-1, Ke);

	// Get the local velocities at the location where they are defined
	c = conc->data;
	dcdn = smag->work3;


	i2dx_c = grid->i2dx_c;
	i2dy_c = grid->i2dy_c;
	i2dz_c = grid->i2dz_c;

	if (idir == 1) {
		for (k=Ks; k<k_end; k++) {
			for (j=Js; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					dcdn[k][j][i] = ( c[k][j][i+1]-c[k][j][i-1])*i2dx_c[i];
				}
			}
		}
#ifndef XPERIODIC
		if (Is == 0) {
			i=0;
			for (k=Ks; k<k_end; k++) {
				for (j=Js; j<j_end; j++) {
					dcdn[k][j][i] = ( c[k][j][i+1]-c[k][j][i])*i2dx_c[i];
				}
			}
		}
#endif
	}
	else if (idir == 2) {
		for (k=Ks; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=Is; i<i_end; i++){
					dcdn[k][j][i] = ( c[k][j+1][i]-c[k][j-1][i])*i2dy_c[j];
				}
			}
		}
		if (Js == 0) {
			j=0;
			for (k=Ks; k<k_end; k++) {
				for (i=Is; i<i_end; i++) {
					dcdn[k][j][i] = ( c[k][j+1][i]-c[k][j][i])*i2dy_c[j];
				}
			}
		}
	}
	else if (idir == 3) {
		for (k=k_start; k<k_end; k++) {
			for (j=Js; j<j_end; j++) {
				for (i=Is; i<i_end; i++){
					dcdn[k][j][i] = ( c[k+1][j][i]-c[k-1][j][i])*i2dz_c[k];
				}
			}
		}
#ifndef ZPERIODIC
		if (Ks == 0) {
			k=0;
			for (j=Js; j<j_end; j++) {
				for (i=Is; i<i_end; i++) {
					dcdn[k][j][i] = ( c[k+1][j][i]-c[k][j][i])*i2dz_c[k];
				}
			}
		}
#endif
	}

//	Communication_update_ghost_nodes_flow_variable(smag->work3, 'c', 1, data_bag);

	return;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Subgrid_smagorinsky_eddy_viscosity(Cart3d_bag *data_bag) {

	int i, j, k;
	double ***nut, ***deltasq, ***strain, ***Cev, ***mSct, ***Sct;
	int NConc, iconc;

	FILE *fid, *fid1;
	char bin_filename[50];
	int W_found, found;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Subgrid *smag = data_bag -> smag;

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	// exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie);
	int j_end = min(NY-1, grid->G_Je);
	int k_end = min(NZ-1, grid->G_Ke);

	double Cs = smag -> Cs;
	double Cs_sq = Cs * Cs;

	found = 0;
/*
	sprintf(bin_filename,"schmidt%d.dat",params->rank);
	fid = fopen(bin_filename,"w");
	sprintf(bin_filename,"schmidtnan%d.dat",params->rank);
	fid1 = fopen(bin_filename,"w");
*/

	// Get the local velocities at the location where they are defined
	nut = smag->nut;
	strain = smag->st_rate->strain;
	deltasq = smag->ng_lengthscale_sq;
#ifdef SMAG_DYNAMIC
	Cev = smag->ng_Cev;
#endif
	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
#ifdef SMAG_DYNAMIC
				nut[k][j][i] = Cev[k][j][i]*deltasq[k][j][i]*strain[k][j][i];
#else
				nut[k][j][i] = Cs_sq*deltasq[k][j][i]*strain[k][j][i];
#endif
//				fprintf(fid,"%20.12e %20.12e %20.12e %20.12e %20.12e \n",
//						grid->yc[j],u[k][j][i+1],v[k][j+1][i],w[k+1][j][i],nut[k][j][i]);
			}
		}
	}

//	Subgrid_boundary(nut, grid, params);
//	Communication_update_ghost_nodes_flow_variable(nut, 'c', 1, data_bag);
	Viscosity_update_boundaries(nut, data_bag);

#ifdef CONC_DYNAMIC
	int iconc;
	Cev = smag->ng_Cev;
	int NConc = params->NConc;
	for (iconc=0;iconc<NConc;iconc++) {
		mSct = smag->cdev[iconc]->mSct;
		Sct  = smag->cdev[iconc]->Sct;
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++) {
					mSct[k][j][i] = nut[k][j][i]/Sct[k][j][i];
//					mSct[k][j][i] = nut[k][j][i]/0.6;
				}
			}
		}
//		Subgrid_boundary(mSct, grid, params);
//		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->mSct, 'c', 1, data_bag);
		Viscosity_update_boundaries(smag->cdev[iconc]->mSct, data_bag);
		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->Sct, 'c', 1, data_bag);
	}
//	fclose(fid);
//	fclose(fid1);

/*
	MPI_Allreduce (&found, &W_found, 1, MPI_INT, MPI_MAX, PCW);
	if (W_found == 1) {
		printf("terminating rank %d \n",params->rank);
		Communication_finalize();
		exit(0);
	}
*/


#endif
	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Subgrid_boundary(double ***var, MAC_grid *grid, Parameters *params) {

	int i, j, k;

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Start index of bottom-left-back corner on current processor
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid -> G_Ie;
	int Je = grid -> G_Je;
	int Ke = grid -> G_Ke;

	// indices start and end on current processor
	int i_start = Is;
	int j_start = Js;
	int k_start = Ks;

	// exclude the half cell added
	int i_end = min(NX-1, Ie);
	int j_end = min(NY-1, Je);
	int k_end = min(NZ-1, Ke);

	if (Ie == NX) {
		i = NX-1;
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				var[k][j][i] = var[k][j][i-1];
			}
		}
	}

	if (Je == NY) {
		j = NY-1;
		for (k=k_start; k<k_end; k++) {
			for (i=i_start; i<i_end; i++){
				var[k][j][i] = var[k][j-1][i];
			}
		}
	}

	if (Ke == NZ) {
		k = NZ-1;
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				var[k][j][i] = var[k-1][j][i];
			}
		}
	}

	if ( (Ie == NX) && (Je == NY) ) {
		i = NX-1;
		j = NY-1;
		for (k=k_start; k<k_end; k++) {
			var[k][j][i] = var[k][j-1][i-1];
		}
	}

	if ( (Je == NY) && (Ke == NZ) ) {
		j = NY-1;
		k = NZ-1;
		for (i=i_start; i<i_end; i++){
			var[k][j][i] = var[k-1][j-1][i];
		}
	}

	if ( (Ke == NZ) && (Ie == NX) ) {
		k = NZ-1;
		i = NX-1;
		for (j=j_start; j<j_end; j++) {
			var[k][j][i] = var[k-1][j][i-1];
		}
	}
	if ( (Ke == NZ) && (Ie == NX) && (Je == NY) ) {
		i = NX-1;
		j = NY-1;
		k = NZ-1;
		var[k][j][i] = var[k-1][j-1][i-1];
	}

	return;

}
