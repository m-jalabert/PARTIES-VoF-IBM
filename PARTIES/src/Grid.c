/* Note: Function names are in alphabetical order (Not exactly) */
/* grid.c
		This file defines the functions to create MAC_grid objects
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Communication.h"
#include "Display.h"
#include "Grid.h"
#include "Immersed.h"
#include "Memory.h"


/******************************************************************************/
/*
 This function allocates memory for the grid object.

 It also holds the object holding the information regarding solid interface and
 solid object.
 */
/******************************************************************************/
MAC_grid *Grid_create(Parameters *params, Debug_trace *dtrace) {

	int i, j, k;

	double **interface_position;
	int **interface_y_index;
	MAC_grid *new_grid;
	int ierr;
	int success_flag;

	FILE *fileout;

	new_grid = (MAC_grid *)malloc(sizeof(MAC_grid));
	Memory_check_allocation(new_grid);

	// Global (world) number of grids. A half cell is added in each direction.
	// Same on all processors
	int NX = params -> NXM + 1;
	int NY = params -> NYM + 1;
	int NZ = params -> NZM + 1;

	// Create DA for the ghost nodes associated with the ENO scheme (convective
	// terms)
	int Ghost_Nodes = params->ghost_nodes;

	int NPX, NPY, NPZ;
	int NP[3];
	double err_NP, err_buff;

	// Computes number of processors in each direction (up up)
	Grid_find_best_NP(NP, &err_NP, 1, 1, params);
	NPX = NP[0];
	NPY = NP[1];
	NPZ = NP[2];
	// Computes number of processors in each direction (up down) and compares
	// to first result
	Grid_find_best_NP(NP, &err_buff, 1, -1, params);
	if (err_buff < err_NP ){
		NPX = NP[0];
		NPY = NP[1];
		NPZ = NP[2];
		err_NP = err_buff;
	}
	// Computes number of processors in each direction (down up) and compares
	// to previous best
	Grid_find_best_NP(NP, &err_buff, -1, 1, params);
	if (err_buff < err_NP ){
		NPX = NP[0];
		NPY = NP[1];
		NPZ = NP[2];
		err_NP = err_buff;
	}
	// Computes number of processors in each direction (down down) and compares
	// to previous best
	Grid_find_best_NP(NP, &err_buff, -1, -1, params);
	if (err_buff < err_NP ){
		NPX = NP[0];
		NPY = NP[1];
		NPZ = NP[2];
		err_NP = err_buff;
	}

	int color, key;
	int ndims = 3;
	int reorder = 1;
	int *dim_size = (int *) calloc(ndims, sizeof(int));
	int *periods  = (int *) calloc(ndims, sizeof(int));
	int *coords  = (int *) calloc(ndims, sizeof(int));

	dim_size[0] = NPX;
	dim_size[1] = NPY;
	dim_size[2] = NPZ;

	periods[0] = 0;
#ifdef XPERIODIC
	periods[0] = 1;
#endif

	periods[1] = 0;
#ifdef YPERIODIC
	periods[1] = 1;
#endif

	periods[2] = 0;
#ifdef ZPERIODIC
	periods[2] = 1;
#endif


	ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &comm3d);
	ierr = MPI_Comm_rank(PCW, &params->rank);
	ierr = MPI_Cart_get(comm3d, ndims, dim_size, periods, coords);

	params->xproccoord = coords[0];
	params->yproccoord = coords[1];
	params->zproccoord = coords[2];

	free(dim_size);
	free(periods);
	free(coords);

	ierr = MPI_Cart_shift(comm3d, 0, 1, &params->npxminus, &params->npxplus);
	ierr = MPI_Cart_shift(comm3d, 1, 1, &params->npyminus, &params->npyplus);
	ierr = MPI_Cart_shift(comm3d, 2, 1, &params->npzminus, &params->npzplus);

	// create communicator of all processors with the same x-coordinate, commx
	color = params->xproccoord;
	key = params->rank;
	ierr = MPI_Comm_split(PCW,color,key,&params->commx);

	// create communicator of all processors with the same y-coordinate, commy
	color = params->yproccoord;
	key = params->rank;
	ierr = MPI_Comm_split(PCW,color,key,&params->commy);

	// create communicator of all processors with the same z-coordinate, commz
	color = params->zproccoord;
	key = params->rank;
	ierr = MPI_Comm_split(PCW,color,key,&params->commz);

	// create communicator of all processors with same  xy-coordinate, commxy
	color = params->yproccoord;
	key = params->rank;
	ierr = MPI_Comm_split(params->commx,color,key,&params->commxy);

	FILE *fid;
	for (i=0;i<params->size;i++){
		if (i==params->rank) {
			if (i==0)
				fid = fopen("proc_info.dat","w");
			else
				fid = fopen("proc_info.dat","a");
			fprintf(fid,".........Rank=%d........\n",params->rank);
			fprintf(fid, "npxplus = %d npxminus = %d\n",params->npxplus,params->npxminus);
			fprintf(fid, "npyplus = %d npyminus = %d\n",params->npyplus,params->npyminus);
			fprintf(fid, "npzplus = %d npzminus = %d\n",params->npzplus,params->npzminus);
			fprintf(fid, "xproc = %d yproc = %d zproc = %d\n",params->xproccoord, params->yproccoord, params->zproccoord);
			fclose(fid);
		}
		ierr = MPI_Barrier(PCW);
	}

	new_grid -> G_Is = (float) (params->xproccoord    ) / (float)(NPX) * NX;
	new_grid -> G_Ie = (float) (params->xproccoord + 1) / (float)(NPX) * NX;

	new_grid -> G_Js = (float) (params->yproccoord    ) / (float)(NPY) * NY;
	new_grid -> G_Je = (float) (params->yproccoord + 1) / (float)(NPY) * NY;

	new_grid -> G_Ks = (float) (params->zproccoord    ) / (float)(NPZ) * NZ;
	new_grid -> G_Ke = (float) (params->zproccoord + 1) / (float)(NPZ) * NZ;

	new_grid -> L_Is = new_grid->G_Is - Ghost_Nodes;
	new_grid -> L_Ie = new_grid->G_Ie + Ghost_Nodes;

	new_grid -> L_Js = new_grid->G_Js - Ghost_Nodes;
	new_grid -> L_Je = new_grid->G_Je + Ghost_Nodes;

	new_grid -> L_Ks = new_grid->G_Ks - Ghost_Nodes;
	new_grid -> L_Ke = new_grid->G_Ke + Ghost_Nodes;

	params -> NPX = NPX;
	params -> NPY = NPY;
	params -> NPZ = NPZ;

	new_grid -> c_status =  Memory_allocate_flow_variable_int(new_grid, params);
	new_grid -> u_status =  Memory_allocate_flow_variable_int(new_grid, params);
	new_grid -> v_status =  Memory_allocate_flow_variable_int(new_grid, params);
	new_grid -> w_status =  Memory_allocate_flow_variable_int(new_grid, params);


	/*------------------------------------------------------------------------*/
	/*
	 Coordinates of the grid-variables
	 */
	/*------------------------------------------------------------------------*/
	int NX_ghost = NX + 2 * Ghost_Nodes;
	int NY_ghost = NY + 2 * Ghost_Nodes;
	int NZ_ghost = NZ + 2 * Ghost_Nodes;

	// u-, v-, w-faced coordinates
	double *xu = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *yv = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *zw = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	xu += Ghost_Nodes;
	yv += Ghost_Nodes;
	zw += Ghost_Nodes;
	new_grid -> xu = xu;
	new_grid -> yv = yv;
	new_grid -> zw = zw;

	// Cell-centered coordinates
	double *xc = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *yc = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *zc = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	xc += Ghost_Nodes;
	yc += Ghost_Nodes;
	zc += Ghost_Nodes;
	new_grid -> xc = xc;
	new_grid -> yc = yc;
	new_grid -> zc = zc;

	// The physical position of the bottom boundary
	interface_position = Memory_allocate_2D_double_array(NX, NZ);

	// y-index(j) of the (first fluid node) on the bottom geometry
	interface_y_index  = Memory_allocate_2D_int_array(NX, NZ);

	new_grid -> fine_step = 2 * 10; // Should be an even number
	new_grid -> finer_1d_interface_x = Memory_allocate_1D_array(GVG_DOUBLE, new_grid->fine_step*(NX) + 1);
	new_grid -> finer_1d_interface_y = Memory_allocate_1D_array(GVG_DOUBLE, new_grid->fine_step*(NX) + 1);

	int Lx = params -> Lx;
	int Ly = params -> Ly;
	int Lz = params -> Lz;

	new_grid -> NX = NX;
	new_grid -> NY = NY;
	new_grid -> NZ = NZ;
	new_grid -> NT = NX * NY * NZ;

	// Number of grids excluding the last half cell
	new_grid -> NI = NX-1;
	new_grid -> NJ = NY-1;
	new_grid -> NK = NZ-1;

	new_grid -> interface_position = interface_position;
	new_grid -> interface_y_index  = interface_y_index;

	/*------------------------------------------------------------------------*/
	/*
	 Set up staggered grid coordinates:
	     - xu: u-faced x-coordinate
	     - yv: v-faced y-coordinate
	     - zw: w-faced z-coordinate
	     - xc: cell-centered x-coordinate
	     - yc: cell-centered y-coordinate
	     - zc: cell-centered z-coordinate

	 User can read the grid coordinates from the file: grid.inp. Look at the
	 read_grid_from_file for more info (below)
	 */
	/*------------------------------------------------------------------------*/
	if (params -> ImportGridFromFile == YES ) {

#ifdef GRID_UNIFORM
		if (params -> rank == 0) {
			printf("*********************************************************\n");
			printf("!!! Warning. Both GRID_UNIFORM and ImportGridFromFile !!!\n");
			printf("!!! specified. Now generating uniform mesh based on   !!!\n");
            printf("!!! input parameters.                                 !!!\n");
			printf("*********************************************************\n");
		}
		success_flag = 0;
#else
		success_flag = Grid_import_grid_from_file(new_grid, params);
#endif
		/*--------------------------------------------------------------------*/
		/*
		 Could not find the input files, i.e.
		     Grid_x.inp
		     Grid_y.inp
		     Grid_z.inp

		 Switch to uniform grid formaulation
		 */
		/*--------------------------------------------------------------------*/
		if (!success_flag) {
			Grid_generate_mesh_coordinates(new_grid, params, DTRACE("Grid_generate_mesh_coordinates"));
		}
		else if (params->rank==0)
			printf("Grid has been imported from file successfully. \n");

	}
	else {
		Grid_generate_mesh_coordinates(new_grid, params, DTRACE("Grid_generate_mesh_coordinates"));
	}

	/*------------------------------------------------------------------------*/
	/*
	 Calculate staggered grid mesh widths:
	     dx_u[i]   = xu[i+1] - xu[i]
	     dx_c[i]   = xc[i+1] - xc[i]
	     idx_u[i]  = 1.0 / dx_u[i]
	     idx_c[i]  = 1.0 / dx_c[i]
	     i2dx_c[i] = 1.0 / ( xc[i+1] - xc[i-1] )
	 */
	/*------------------------------------------------------------------------*/

	double *dx_u = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *dx_c = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	dx_u += Ghost_Nodes;
	dx_c += Ghost_Nodes;
	new_grid -> dx_u = dx_u;
	new_grid -> dx_c = dx_c;

	double *dy_v = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *dy_c = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	dy_v += Ghost_Nodes;
	dy_c += Ghost_Nodes;
	new_grid -> dy_v = dy_v;
	new_grid -> dy_c = dy_c;

	double *dz_w = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	double *dz_c = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	dz_w += Ghost_Nodes;
	dz_c += Ghost_Nodes;
	new_grid -> dz_w = dz_w;
	new_grid -> dz_c = dz_c;

	double *idx_u  = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *idx_c  = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *i2dx_c = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	idx_u += Ghost_Nodes;
	idx_c += Ghost_Nodes;
	i2dx_c += Ghost_Nodes;
	new_grid -> idx_u  = idx_u;
	new_grid -> idx_c  = idx_c;
	new_grid -> i2dx_c = i2dx_c;

	double *idy_v  = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *idy_c  = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *i2dy_c = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	idy_v  += Ghost_Nodes;
	idy_c  += Ghost_Nodes;
	i2dy_c += Ghost_Nodes;
	new_grid -> idy_v  = idy_v;
	new_grid -> idy_c  = idy_c;
	new_grid -> i2dy_c = i2dy_c;

	double *idz_w  = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	double *idz_c  = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	double *i2dz_c = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	idz_w  += Ghost_Nodes;
	idz_c  += Ghost_Nodes;
	i2dz_c += Ghost_Nodes;
	new_grid -> idz_w  = idz_w;
	new_grid -> idz_c  = idz_c;
	new_grid -> i2dz_c = i2dz_c;

	int i_start = -Ghost_Nodes;
	int j_start = -Ghost_Nodes;
	int k_start = -Ghost_Nodes;

	int i_end = NX + Ghost_Nodes - 1;
	int j_end = NY + Ghost_Nodes - 1;
	int k_end = NZ + Ghost_Nodes - 1;

	//--------------------------------------------------------------------------
	// dx, idx
	//--------------------------------------------------------------------------
	for (i = i_start; i < i_end; i++) {
		dx_u[i] = xu[i+1] - xu[i];
		dx_c[i] = xc[i+1] - xc[i];

		idx_u[i] = 1.0 / dx_u[i];
		idx_c[i] = 1.0 / dx_c[i];
	}

	//--------------------------------------------------------------------------
	// dy, idy
	//--------------------------------------------------------------------------
	for (j = j_start; j < j_end; j++) {
		dy_v[j] = yv[j+1] - yv[j];
		dy_c[j] = yc[j+1] - yc[j];

		idy_v[j] = 1.0 / dy_v[j];
		idy_c[j] = 1.0 / dy_c[j];
	}

	//--------------------------------------------------------------------------
	// dz, idz
	//--------------------------------------------------------------------------
	for (k = k_start; k < k_end; k++) {
		dz_w[k] = zw[k+1] - zw[k];
		dz_c[k] = zc[k+1] - zc[k];

		idz_w[k] = 1.0 / dz_w[k];
		idz_c[k] = 1.0 / dz_c[k];
	}

	//--------------------------------------------------------------------------
	// i2dx
	//--------------------------------------------------------------------------
	for (i = i_start + 1; i < i_end; i++){
		i2dx_c[i] = xc[i+1] - xc[i-1];
		i2dx_c[i] = 1.0 / i2dx_c[i];
	}

	//--------------------------------------------------------------------------
	// i2dy
	//--------------------------------------------------------------------------
	for (j = j_start + 1; j < j_end; j++) {
		i2dy_c[j] = yc[j+1] - yc[j-1];
		i2dy_c[j] = 1.0 / i2dy_c[j];
	}

	//--------------------------------------------------------------------------
	// i2dz
	//--------------------------------------------------------------------------
	for (k = k_start + 1; k < k_end; k++) {
		i2dz_c[k] = zc[k+1] - zc[k-1];
		i2dz_c[k] = 1.0 / i2dz_c[k];
	}


#ifdef BICG_SOLVE
	/*------------------------------------------------------------------------*/
	/*
	 Calculate staggered grid mesh widths required for BICG
	 */
	/*------------------------------------------------------------------------*/

	double *idx2_e = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *idx2_w = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	idx2_e += Ghost_Nodes;
	idx2_w += Ghost_Nodes;
	new_grid -> idx2_e = idx2_e;
	new_grid -> idx2_w = idx2_w;

	double *idy2_n = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *idy2_s = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	idy2_n += Ghost_Nodes;
	idy2_s += Ghost_Nodes;
	new_grid -> idy2_n = idy2_n;
	new_grid -> idy2_s = idy2_s;

	double *idz2_f = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	double *idz2_b = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	idz2_f += Ghost_Nodes;
	idz2_b += Ghost_Nodes;
	new_grid -> idz2_f = idz2_f;
	new_grid -> idz2_b = idz2_b;

	//--------------------------------------------------------------------------
	// idx2
	//--------------------------------------------------------------------------
	for (i = i_start + 1; i < i_end; i++) {
		idx2_e[i] = idx_c[i]   * idx_u[i];
		idx2_w[i] = idx_c[i-1] * idx_u[i];
	}

	//--------------------------------------------------------------------------
	// idy2
	//--------------------------------------------------------------------------
	for (j = j_start + 1; j < j_end; j++) {
		idy2_n[j] = idy_c[j]   * idy_v[j];
		idy2_s[j] = idy_c[j-1] * idy_v[j];
	}

	//--------------------------------------------------------------------------
	// idz2
	//--------------------------------------------------------------------------
	for (k = k_start + 1; k < k_end; k++) {
		idz2_f[k] = idz_c[k]   * idz_w[k];
		idz2_b[k] = idz_c[k-1] * idz_w[k];
	}
#endif


	/*------------------------------------------------------------------------*/
	/*
	 weights to interpolate from center to face centered
	 */
	/*------------------------------------------------------------------------*/
#ifndef GRID_UNIFORM

	double *wc2uW = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *wc2uE = Memory_allocate_1D_array(GVG_DOUBLE, NX_ghost);
	double *wc2vN = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *wc2vS = Memory_allocate_1D_array(GVG_DOUBLE, NY_ghost);
	double *wc2wF = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);
	double *wc2wB = Memory_allocate_1D_array(GVG_DOUBLE, NZ_ghost);

	wc2uW += Ghost_Nodes;
	wc2uE += Ghost_Nodes;
	wc2vN += Ghost_Nodes;
	wc2vS += Ghost_Nodes;
	wc2wF += Ghost_Nodes;
	wc2wB += Ghost_Nodes;


	new_grid -> wc2uW = wc2uW;
	new_grid -> wc2uE = wc2uE;
	new_grid -> wc2vN = wc2vN;
	new_grid -> wc2vS = wc2vS;
	new_grid -> wc2wF = wc2wF;
	new_grid -> wc2wB = wc2wB;

	for (i = i_start; i < i_end; i++) {
			wc2uW[i] = (xc[i] - xu[i]  ) / (xc[i] - xc[i-1]);
			wc2uE[i] = (xu[i] - xc[i-1]) / (xc[i] - xc[i-1]);
		}
		wc2uE[NX-1] = 0.5;
		wc2uW[NX-1] = 0.5;

	for (j = j_start; j < j_end; j++) {
			wc2vS[j] = (yc[j] - yv[j]  ) / (yc[j] - yc[j-1]);
			wc2vN[j] = (yv[j] - yc[j-1]) / (yc[j] - yc[j-1]);
		}
		wc2vN[NY-1] = 0.5;
		wc2vS[NY-1] = 0.5;
	for (k = k_start; k < k_end; k++) {
			wc2wB[k] = (yc[k] - yv[k]  ) / (yc[k] - yc[k-1]);
			wc2wF[k] = (yv[k] - yc[k-1]) / (yc[k] - yc[k-1]);
		}
		wc2wF[NZ-1] = 0.5;
		wc2wB[NZ-1] = 0.5;



#endif


	/*------------------------------------------------------------------------*/
	/*
	 Write grid info to file and find dy_min;
	 */
	/*------------------------------------------------------------------------*/
	new_grid -> dy_min = dy_v[0];
	for (j = 1; j < NY; j++) {
		new_grid -> dy_min = min(new_grid->dy_min, dy_v[j]);
	}

	if (params->rank == 0) {

		//----------------------------------------------------------------------
		// invgrid_info.dat
		//----------------------------------------------------------------------
		fileout = fopen("invgrid_info.dat","w");
		for (i = 0; i < NX; i++) {
			fprintf(fileout,"du = %16.12f dc = %16.12f 2dc = %16.12f \n",
			        idx_u[i],idx_c[i],i2dx_c[i]);
		}
		fprintf(fileout,"\n \n \n");
		for (j = 0; j < NY; j++) {
			fprintf(fileout,"dv = %16.12f dc = %16.12f 2dc = %16.12f \n",
			        idy_v[j],idy_c[j],i2dy_c[j]);
		}
		fprintf(fileout,"\n \n \n");
		for (k = 0; k < NZ; k++) {
			fprintf(fileout,"dw = %16.12f dc = %16.12f 2dc = %16.12f \n",
			        idz_w[k],idz_c[k],i2dz_c[k]);
		}
		fprintf(fileout,"\n \n \n");
#ifndef GRID_UNIFORM
		for (j = 0; j < NY; j++) {
			fprintf(fileout,"wc2vN = %f wc2vS = %f \n",wc2vN[j],wc2vS[j]);
		}
		fprintf(fileout,"\n \n \n");
#endif
		fclose(fileout);

		//----------------------------------------------------------------------
		// cellgrid.dat
		//----------------------------------------------------------------------
		fileout = fopen("cellgrid.dat","w");
		for (i = -Ghost_Nodes; i < NX + Ghost_Nodes; i++) {
			fprintf(fileout,"xu = %16.12f xc = %16.12f \n",
			        new_grid->xu[i],new_grid->xc[i]);
		}
		fprintf(fileout,"\n \n");

		for (j = -Ghost_Nodes; j < NY + Ghost_Nodes; j++) {
			fprintf(fileout,"yv  = %16.12f yc = %16.12f \n",
			        new_grid->yv[j],new_grid->yc[j]);
		}
		fprintf(fileout,"\n \n");

		for (k = -Ghost_Nodes; k < NZ + Ghost_Nodes; k++) {
			fprintf(fileout,"zw = %16.12f zc = %16.12f \n",
			        new_grid->zw[k],new_grid->zc[k]);
		}
		fprintf(fileout,"\n \n");
		fclose(fileout);
	}

	//--------------------------------------------------------------------------
	// grid_info.dat
	//--------------------------------------------------------------------------
	for (i = 0; i < params->size; i++) {
		if (i == params->rank) {
			if (i == 0)
				fileout = fopen("grid_info.dat","w");
			else
				fileout = fopen("grid_info.dat","a");

			fprintf(fileout,".........Rank=%d........\n",params->rank);
			fprintf(fileout,"Is= %d Js= %d Ks = %d\n",new_grid->G_Is,new_grid->G_Js,new_grid->G_Ks);
			fprintf(fileout,"Ie= %d Je= %d Ke = %d\n",new_grid->G_Ie,new_grid->G_Je,new_grid->G_Ke);
			fprintf(fileout,"Is_g= %d Js_g= %d Ks_g = %d\n",new_grid->L_Is,new_grid->L_Js,new_grid->L_Ks);
			fprintf(fileout,"Ie_g= %d Je_g= %d Ke_g = %d\n",new_grid->L_Ie,new_grid->L_Je,new_grid->L_Ke);
			fclose(fileout);
		}
		ierr = MPI_Barrier(PCW);
	}

	new_grid -> total_nodes    = (new_grid->L_Ie - new_grid->L_Is) * (new_grid->L_Je - new_grid->L_Js) * (new_grid->L_Ke - new_grid->L_Ks);
	new_grid -> ng_total_nodes = (new_grid->G_Ie - new_grid->G_Is) * (new_grid->G_Je - new_grid->G_Js) * (new_grid->G_Ke - new_grid->G_Ks);

	new_grid -> u_sdf =  Memory_allocate_flow_variable(new_grid, params);
	new_grid -> v_sdf =  Memory_allocate_flow_variable(new_grid, params);
	new_grid -> w_sdf =  Memory_allocate_flow_variable(new_grid, params);
	new_grid -> c_sdf =  Memory_allocate_flow_variable(new_grid, params);



#ifdef GRID_UNIFORM
	char message[200];
	double dx = params -> Lx / (double)( NX - 1 );
	double dy = params -> Ly / (double)( NY - 1 );
	double dz = params -> Lz / (double)( NZ - 1 );


	if (fabs(dx-dy)/dx > 1e-12 || fabs(dx-dz)/dx > 1e-12 || fabs(dy-dz)/dx > 1e-12) {
			sprintf(message,"GRID_UNIFORM requires a uniform mesh of\n"
			                "dx = dy = dz. Mesh based on input parameters has\n"
			                "dx = %.4g, dy = %.4g, dz = %.4g", dx, dy, dz);
			Display_throw_error(message, params, DTRACE("Display_throw_error"));
		}
#endif



#ifdef YPERIODIC
	new_grid->exchange_slab = Memory_allocate_2D_double_array(NX, NZ);
#endif



	return new_grid;
}




/******************************************************************************/
/*
 This function releases allocated memory for the variable grid-type
 */
/******************************************************************************/
void Grid_destroy(MAC_grid *grid, Parameters *params) {

	int NZ, NY;
	int Ghost_Nodes = params -> ghost_nodes;

	NZ = grid -> NZ;
	NY = grid -> NY;

	Memory_free_flow_variable_int(grid, params, grid -> c_status);
	Memory_free_flow_variable_int(grid, params, grid -> u_status);
	Memory_free_flow_variable_int(grid, params, grid -> v_status);
	Memory_free_flow_variable_int(grid, params, grid -> w_status);

	free(grid -> xu - Ghost_Nodes);
	free(grid -> yv - Ghost_Nodes);
	free(grid -> zw - Ghost_Nodes);

	free(grid -> xc - Ghost_Nodes);
	free(grid -> yc - Ghost_Nodes);
	free(grid -> zc - Ghost_Nodes);

	Memory_free_2D_double_array(NZ, grid -> interface_position);
	Memory_free_2D_int_array(NZ, grid -> interface_y_index);

	free(grid -> finer_1d_interface_x);
	free(grid -> finer_1d_interface_y);

	Memory_free_flow_variable(grid, params, grid -> u_sdf);
	Memory_free_flow_variable(grid, params, grid -> v_sdf);
	Memory_free_flow_variable(grid, params, grid -> w_sdf);
	Memory_free_flow_variable(grid, params, grid -> c_sdf);

	free(grid -> dx_u - Ghost_Nodes);
	free(grid -> dy_v - Ghost_Nodes);
	free(grid -> dz_w - Ghost_Nodes);

	free(grid -> dx_c - Ghost_Nodes);
	free(grid -> dy_c - Ghost_Nodes);
	free(grid -> dz_c - Ghost_Nodes);

	free(grid -> idx_u - Ghost_Nodes);
	free(grid -> idx_c - Ghost_Nodes);
	free(grid -> i2dx_c - Ghost_Nodes);

	free(grid -> idy_v - Ghost_Nodes);
	free(grid -> idy_c - Ghost_Nodes);
	free(grid -> i2dy_c - Ghost_Nodes);

	free(grid -> idz_w - Ghost_Nodes);
	free(grid -> idz_c - Ghost_Nodes);
	free(grid -> i2dz_c - Ghost_Nodes);

#ifdef YPERIODIC
	Memory_free_2D_double_array(NZ, grid->exchange_slab);
#endif

#ifdef BICG_SOLVE
	free(grid -> idx2_e - Ghost_Nodes);
	free(grid -> idx2_w - Ghost_Nodes);

	free(grid -> idy2_n - Ghost_Nodes);
	free(grid -> idy2_s - Ghost_Nodes);

	free(grid -> idz2_f - Ghost_Nodes);
	free(grid -> idz2_b - Ghost_Nodes);
#endif

#ifndef GRID_UNIFORM
	free(grid -> wc2vN - Ghost_Nodes);
	free(grid -> wc2vS - Ghost_Nodes);
	free(grid -> wc2uW - Ghost_Nodes);
	free(grid -> wc2uE - Ghost_Nodes);
	free(grid -> wc2wB - Ghost_Nodes);
	free(grid -> wc2wF - Ghost_Nodes);


#endif


#ifdef IMMERSED_BOUNDARY
	Immersed_destroy(grid -> u_immersed, params, grid);
	Immersed_destroy(grid -> v_immersed, params, grid);
	Immersed_destroy(grid -> w_immersed, params, grid);
	Immersed_destroy(grid -> c_immersed, params, grid);
#endif

	free(grid);
}




/******************************************************************************/
/*
 This function creates the grid coordinates for Uniform mesh!!
 */
/******************************************************************************/
void Grid_generate_mesh_coordinates(MAC_grid *grid, Parameters *params,
		Debug_trace *dtrace) {

	int i, j, k;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;
	int NX, NY, NZ;
	int Ghost_Nodes;
	double dx, dy, dz;
	double xmin, ymin, zmin;

	char message[200];

	// Half cell is added to the far most place
	NX = grid -> NX;
	NY = grid -> NY;
	NZ = grid -> NZ;

	Ghost_Nodes = params -> ghost_nodes;

	i_start = -Ghost_Nodes;
	j_start = -Ghost_Nodes;
	k_start = -Ghost_Nodes;

	i_end = NX + Ghost_Nodes;
	j_end = NY + Ghost_Nodes;
	k_end = NZ + Ghost_Nodes;

	dx = params -> Lx / (double)( NX - 1 );
	dy = params -> Ly / (double)( NY - 1 );
	dz = params -> Lz / (double)( NZ - 1 );

	xmin = params -> xmin;
	ymin = params -> ymin;
	zmin = params -> zmin;

#ifdef LAG_PARTICLE_RESOLVED
	// LAG_PARTICLE_RESOLVED requires that grid has uniform spacing
	if (fabs(dx-dy)/dx > 1e-12 || fabs(dx-dz)/dx > 1e-12 || fabs(dy-dz)/dx > 1e-12) {
		sprintf(message,"LAG_PARTICLES_RESOLVED requires a uniform mesh of\n"
		                "dx = dy = dz. Mesh based on input parameters has\n"
		                "dx = %.4g, dy = %.4g, dz = %.4g", dx, dy, dz);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}
#endif

	for (i = i_start; i < i_end; i++) {

		grid -> xu[i] = xmin + i * dx;
		grid -> xc[i] = xmin + (0.5 + i) * dx;
	}

	for (j = j_start; j < j_end; j++) {

		grid -> yv[j] = ymin + j * dy;
		grid -> yc[j] = ymin + (0.5 + j) * dy;
	}

	for (k = k_start; k < k_end; k++) {

		grid -> zw[k] = zmin + k * dz;
		grid -> zc[k] = zmin + (0.5 + k) * dz;
	}

	return;
}




/******************************************************************************/
/*
 This function describes the interface !!! I have to check this. MOHAMAD
 */
/******************************************************************************/
void Grid_describe_interface(MAC_grid *grid, Parameters *params) {

	double **interface_position;
	double *xc;
	double *yc;
	double *zc;
	double *interface_xfine, *interface_yfine, dx_fine;
	int NX, NZ;
	int i, k;
	double Lx=params->Lx;
	double Ly=params->Ly;
	double Lz=params->Lz;

	int Bottom_Geometry_Shape_Type;

	/*------------------------------------------------------------------------*/
	/*
		0:   Flat surface on the bottom
		1:   Inclinded surface   y=f(x)
		2:   One Errorfucntion ramp at the beginning y=f(x)
		3:   1D Gaussian bump (in y=f(x)
		4:   2D Gaussian bump     y=f(x,z)
		5:   Errorfunction ramp + 1d bump : y=f_1(x) + f_2(x)
		6:   Errorfunction ramp + 2d bump : y =f_1(x) + f_2(x,z)
	 */
	/*------------------------------------------------------------------------*/
	double m_ramp;
	double ramp_start_x;
	double ramp_end_x;
	double ramp_start_y;
	double ramp_end_y;

	double max_height;
	double ramp_center;

	double bump_center_x;
	double bump_height;
	double bump_width;
	double bump_length;
	double bump_start;
	double bump_end;

	double bump_center_z;
	double ramp, bump;

	double channel_end_x   = 1.5 / 2.5 * Lx;
	double channel_width_z = 0.1 / 0.5 * Lz;

	/* case 7 */
	double obstacle_center_x;
	double obstacle_center_z;
	double obstacle_height;
	double obstacle_radius_top;
	double obstacle_radius_bottom;


	/* case 8 */
	double *xr, *zr;
	double d, d2, d2min;
	int c;

	int success = NO;
	char filename[100];
	int ierr;

	double xstep, step_height;
	double S;
	int fine_grid_step;

	// If the flag is on, then try to import the bottom interface from
	// "Bottom_Topography.inp" file
	if (params->ImportBottomInterface) {

		sprintf(filename, "Bottom_Topography.inp");
		success = Grid_import_bottom_interface(grid, params, filename) ;
	}
	// If it was successful, then return. Otherwise, generate the interface
	// based on the default value
	if (success) {

		if (params->rank==0)
			printf("Bottom interface topography has been imported from file successfully.\n");
		return;
	}

	Lx   = params->Lx;
	Ly   = params->Ly;
	Lz   = params->Lz;

	interface_position = grid->interface_position;

	xc = grid->xc;
	yc = grid->yc;
	zc = grid->zc;

	NX = grid->NX;
	NZ = grid->NZ;

	interface_xfine = grid->finer_1d_interface_x;
	interface_yfine = grid->finer_1d_interface_y;
	fine_grid_step = grid->fine_step;

	for (i=0;i<NX;i++) {
		interface_xfine[i*fine_grid_step] = grid->xu[i];  // Invalid write
	}

	for (i=0;i<NX;i++) {
		dx_fine = grid->dx_u[i]/fine_grid_step;
		for (k=1;k<fine_grid_step;k++){
			interface_xfine[i*fine_grid_step + k] = interface_xfine[i*fine_grid_step] + k*dx_fine;
		}
	}

	Bottom_Geometry_Shape_Type = 0;

	if (Bottom_Geometry_Shape_Type !=0) {
	#ifndef IMMERSED_BOUNDARY
		if (params->rank==0) {
			printf("Define IMMERSED_BOUNDARY in Boundary.h\n");
			printf("Program execution stopped \n");
		}
		Communication_finalize();
		exit(0);
	#endif
	}


	/*------------------------------------------------------------------------*/
	/* Other options are
		0:   Flat surface on the bottom
		1:   Inclinded surface   y=f(x)
		2:   One Errorfucntion ramp at the beginning y=f(x)
		3:   1D Gaussian bump (in y=f(x)
		4:   2D Gaussian bump     y=f(x,z)
		5:   Errorfunction ramp + 1d bump : y=f_1(x) + f_2(x)
		6:   Errorfunction ramp + 2d bump : y =f_1(x) + f_2(x,z)
		10:  Cuboidal obstacle, obstacle_width=Lz
		11:  Cuboidal obstacle, obstacle_width<Lz
	 */
	/*------------------------------------------------------------------------*/
	switch (Bottom_Geometry_Shape_Type) {

		case 0: // Flat surface on the bottom
			break;

		case 1: // Inclinded surface

			ramp_start_x = 0.0;
			ramp_end_x   = Lx;
			ramp_start_y = 0.5*Ly;
			ramp_end_y   = 0.0*Ly;
			m_ramp       = (ramp_start_y - ramp_end_y) / (ramp_start_x - ramp_end_x);

			break;

		case 2: // One Errorfucntion ramp at the beginning

			max_height  = 0.5*Ly;
			ramp_center = 0.25*Lx;

			break;

		case 3: // 1D Gaussian bump

			bump_center_x = 0.5*Lx;
			bump_height   = 0.5*Ly;
			bump_width    = 0.5*bump_height; // not exaclty, standard deviation

			break;

		case 4: // 2D Gaussian bump

			bump_center_x = 8.0+0.5*2.0/3.0;
			bump_center_z = 0.5*Lz;
			bump_height   = 0.20*Ly;
			bump_width    = 0.25*bump_height; // not exaclty, standard deviation

			break;

		case 5: // Errorfunction ramp + 1d bump

			max_height  = 0.5*Ly;
			ramp_center = 0.25*Lx;


			bump_center_x = 0.6*Lx;
			bump_height   = 0.5*Ly;
			bump_width    = 0.5*bump_height; // not exaclty, standard deviation

			break;

		case 6: // Errorfunction ramp + 2d bump

			max_height  = 0.5*Ly;
			ramp_center = 0.25*Lx;

			bump_center_x = 0.6*Lx;
			bump_center_z = 0.5*Lz;
			bump_height   = 0.5*Ly;
			bump_width    = 0.41*bump_height; // not exaclty, standard deviation

			break;

		case 7: // Al'Jaidi 's experiment

			channel_end_x   = 1.5/2.5*Lx;
			channel_width_z = 0.1/0.5*Lz;

			obstacle_center_x      = channel_end_x + 5.0*0.3;
			obstacle_center_z      = 0.0;
			obstacle_height        = 0.048*5.0;
			obstacle_radius_top    = 0.09*5.0;
			obstacle_radius_bottom = 0.15*5.0;

			break;
		case 8: // Flow in the river

			// (xr,zr): coordinates of the centerline of the river in th xz
			// plane
			zr = Memory_allocate_1D_array(GVG_DOUBLE, NX);
			xr = grid->xc;

			for (i=0; i<NX; i++) {

				zr[i] = Lz/2.0 + 1.5*erf( (xc[i]-Lx/2.0)/2.0 );
			} // for

			break;
		case 9: // step in x-direction

			xstep = 0.333*Lx;
			step_height = 0.5*Ly;
		break;

		case 10:
			bump_height = 0.18;
			bump_length = 0.35;
			bump_start = 5.0;
			bump_end = bump_start + bump_length;
		break;

		case 11:
			bump_height = 0.18;
			bump_length = 0.35;
			bump_start = 5.0;
			bump_width = Lz/4.;
			bump_end = bump_start + bump_length;
			bump_center_z=0.5*Lz;

	} // switch

/*
	for (i=0; i<fine_grid_step*(NX-1) + 1; i++) {
		bump = bump_height * exp(-(interface_xfine[i]-bump_center_x)*(interface_xfine[i]-bump_center_x)/(2.0*bump_width*bump_width));

		interface_yfine[i] = bump;
	}
*/

	for (k=0; k<NZ; k++) {
		for (i=0; i<NX; i++) {

			switch (Bottom_Geometry_Shape_Type) {

				case 0: // Flat surface on the bottom

					interface_position[k][i] = 0.0;
					interface_position[k][i] = params->ymin;
					break;

				case 1: // Inclinded surface

					if ( (xc[i] > ramp_start_x) && (xc[i] <= ramp_end_x) ){

						ramp = ramp_start_y + m_ramp*(xc[i] - ramp_start_x);
					}
					else {

						ramp = 0.0;
					}

					interface_position[k][i] = ramp;
					break;

				case 2: // One Errorfucntion ramp at the beginning

					ramp = max_height*0.5*(1.0 + 1.0*erf( (ramp_center-xc[i]) / (0.5*max_height)));

					interface_position[k][i] = ramp;
					break;

				case 3: // 1D Gaussian bump

					bump = bump_height * exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width));

					interface_position[k][i] = bump;
					break;

				case 4: // 2D Gaussian bump

					bump = bump_height *
						exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width)) *
						exp(-(zc[k]-bump_center_z)*(zc[k]-bump_center_z)/(2.0*bump_width*bump_width));

					interface_position[k][i] = bump;
					break;

				case 5: // Errorfunction ramp + 1d bump

					ramp = max_height*0.5*(1.0 + 1.0*erf( (ramp_center-xc[i]) / (0.5*max_height)));
					bump = bump_height * exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width));

					interface_position[k][i] = bump + ramp;
					break;

				case 6: // Errorfunction ramp + 2d bump

					ramp = max_height*0.5*(1.0 + 1.0*erf( (ramp_center-xc[i]) / (0.5*max_height)));
					bump = bump_height *
						exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width)) *
						exp(-(zc[k]-bump_center_z)*(zc[k]-bump_center_z)/(2.0*bump_width*bump_width));

					interface_position[k][i] = bump + ramp;
					break;

				case 7: // Al'Jaidi 's experiment

					if (xc[i] <= channel_end_x) {

						if ( zc[k] <= channel_width_z ) { // inside channel
							interface_position[k][i] = 0.0;
						}
						else { // solid

							interface_position[k][i] = Ly;
						} // if-else
					}
					else {
						double r = sqrt( (xc[i] - obstacle_center_x)*(xc[i] - obstacle_center_x)
						     + (zc[k] - obstacle_center_z)*(zc[k] - obstacle_center_z) ) ;

						if (r <= obstacle_radius_bottom) {
							if (r >= obstacle_radius_top) {
								interface_position[k][i] = obstacle_height - (r - obstacle_radius_top) / (obstacle_radius_bottom - obstacle_radius_top ) * obstacle_height;
							}
							else {
								interface_position[k][i] = obstacle_height;
							}
						}
						else {
							interface_position[k][i] = 0.0;
						} // else
					} // else
					break;

				case 8: // Flow in the river

					d2min = 1.0e10;
					for (c=0; c<NX; c++) {

						// Find the distance of the current node from each node
						// on the center line of river
						d2 = (xc[i] - xr[c])*(xc[i] - xr[c]) + (zc[k] - zr[c])*(zc[k] - zr[c]);
						if (d2 < d2min) {

							d2min = d2;
						} // if
					} // for
					d = sqrt(d2min);

					S = 2.0*PI;
//					interface_position[k][i] = 0.5*Ly*(1.0 - cos(2.0*PI*d/S)*exp(-d*d/(0.5*PI*0.5*PI) ) );

					interface_position[k][i] = 0.5*Ly*(1.0 - cos(2.0*PI*d/S)*exp(-0.5*d*d/(0.25*0.25*S*S) ) );
					break;

				case 9: // step in x-direction

					if (xc[i] < xstep) {

						interface_position[k][i] = step_height;
					} else {

						interface_position[k][i] = 0.0;
					}
					break;

				case 10: // cuboidal bump in x-direction

					if ((xc[i] < bump_start) || (xc[i] > bump_end) ) {

						interface_position[k][i] = 0.0;

					} else {

						interface_position[k][i] = bump_height;
					}
					break;
				case 11: // cuboidal bump in x and z-direction

					if ((xc[i] >= bump_start) && (xc[i] <= bump_end) && (zc[k] >= bump_center_z-bump_width/2) && (zc[k]<=bump_center_z+bump_width/2)) {

						interface_position[k][i] = bump_height;

					} else {
						interface_position[k][i] = 0.0;
					}
					break;

				case 12:
					interface_position[k][i] = 2.*0.1*cos(PI*xc[i]);
				break;
				case 13:
//					if (xc[i] < 0.95*Lx)
//						interface_position[k][i] = 0.95*Ly*(1. - (xc[i]-grid->xu[0])/(0.95*Lx)) ;
//					else
//						interface_position[k][i] = 0.0;

					if (xc[i] < 0.03)
						interface_position[k][i] = 0.97;
					else if (xc[i] < 14.2)
						interface_position[k][i] = 0.97*(1. - (xc[i]-grid->xu[0])/(14.17)) ;
					else
						interface_position[k][i] = 0.0;
				break;
				default:

					interface_position[k][i] = 0.0;
			} // switch
		} // for i
	}// for k

#ifdef IMMERSED_BOUNDARY
	for (i=0; i<fine_grid_step*(NX) + 1; i++) {
		//interface_yfine[i] = 0.1*(cos(PI*interface_xfine[i])+1);
//		if (interface_xfine[i] < 0.95*Lx)
//			interface_yfine[i] = 0.95*Ly*(1.0 - (interface_xfine[i]-grid->xu[0])/(0.95*Lx)) ;
//		else
//			interface_yfine[i] = 0.0;

		if (interface_xfine[i] < 0.03)
			interface_yfine[i] = 0.97;
		else if (interface_xfine[i] < 13.03)
			interface_yfine[i] = 0.97*(1. - (interface_xfine[i]-grid->xu[0])/13.03) ;
		else
			interface_yfine[i] = 0.0;
	}

	FILE *fid;
	if (params->rank==0) {
		fid=fopen("interface.dat","w");
		for (i=0; i<fine_grid_step*(NX) + 1; i++) {
			fprintf(fid, "%e %e \n",interface_xfine[i], interface_yfine[i]);
		}
		fclose(fid);
	}
#endif


}




/******************************************************************************/
/*
 This function reads the coordinates of the grid in each direction from file.
 Note that you need three files:

     1- "Grid_x.inp"
     2- "Grid_y.inp"
     3- "Grid_z.inp"

 In each file (txt files), the data is in the form:
     - First N nodes are the cell center coordinates, e.g. if "NX" is the
     "number of cells", then the first NX coordinates are corresponding to the
     xc.
     - Then, next NX+1 are the xu coordinates.

 The same rule is applied to the other 2 coordinates.  If it can not find any of
 those files, then the grid is assumed to be uniform again
 */
/******************************************************************************/
int Grid_import_grid_from_file(MAC_grid *grid, Parameters *params) {

	int success = YES;
	int status;
	int i, j, k;
	char filename[100];

	double xc_, xu_;
	double yc_, yv_;
	double zc_, zw_;

	double dxu_l, dxu_r;
	double dyv_l, dyv_r;
	double dzw_l, dzw_r;

	FILE *file_in;

	int NX = params -> NXM;
	int NY = params -> NYM;
	int NZ = params -> NZM;

	int Ghost_Nodes = params -> ghost_nodes;

	/*------------------------------------------------------------------------*/
	/*
	 First, read the grid in x-direction
	 */
	/*------------------------------------------------------------------------*/
	sprintf(filename, "Grid_x.inp");
	file_in  = fopen(filename, "r");
	if (file_in == NULL) {
		if (params->rank==0)
			printf("Warning!\nCould not open the grid file \"%s\". Using uniform grid formulation\n", filename);
		success = NO;
	}
	else {

		for (i=0; i<NX; i++) {

			status = fscanf(file_in, "%lf\n", &xc_);
			if (status < 1) {
				success = NO;
				break;
			}
			grid->xc[i] = xc_;
		}

		for (i=NX; i<(2*NX+1); i++) {

			status = fscanf(file_in, "%lf\n", &xu_);
			if (status < 1) {
				success = NO;
				break;
			}
			grid->xu[i-NX] = xu_;
		}
		fclose(file_in);

#ifdef XPERIODIC
		for (i = 1; i <= Ghost_Nodes; i++) {
			grid -> xc[NX+i] = grid -> xc[i];
			grid -> xc[0-i]  = grid -> xc[NX-i];

			grid -> xu[NX+i] = grid -> xu[i];
			grid -> xu[0-i]  = grid -> xu[NX-i];
		}
		grid->xc[NX] = grid->xc[0];
#else
		// Make ghost cells same size as boundary cells
		dxu_l = grid -> xu[1]  - grid -> xu[0];
		dxu_r = grid -> xu[NX] - grid -> xu[NX-1];

		for (i = 1; i <= Ghost_Nodes; i++) {
			grid -> xu[NX+i] = grid -> xu[NX] + i * dxu_r;
			grid -> xu[0-i]  = grid -> xu[0]  - i * dxu_l;

			grid -> xc[NX+i] = grid -> xu[NX+i] + 0.5 * dxu_r;
			grid -> xc[0-i]  = grid -> xu[0-i]  - 0.5 * dxu_l;
		}

		// Make first ghost nodes equidistant from boundary
		grid->xc[NX] = grid->xc[NX-1] + 2.0 * (grid->xu[NX] - grid->xc[NX-1]);
		grid->xc[-1] = grid->xc[0]    - 2.0 * (grid->xc[0]  - grid->xu[0]);
#endif

	}

	/*------------------------------------------------------------------------*/
	/*
	 Read the grid in y-direction
	 */
	/*------------------------------------------------------------------------*/
	sprintf(filename, "Grid_y.inp");
	file_in  = fopen(filename, "r");
	if (file_in == NULL) {
		if (params->rank==0)
			printf("Warning!\nCould not open the grid file \"%s\". Using uniform grid formulation\n", filename);
		success = NO;
	}
	else {

		for (j=0; j<NY; j++) {

			status = fscanf(file_in, "%lf\n", &yc_);
			if (status < 1) {
				success = NO;
				break;
			}
			grid->yc[j] = yc_;
		}

		for (j=NY; j<(2*NY+1); j++) {

			status = fscanf(file_in, "%lf\n", &yv_);
			if (status < 1) {
				success = NO;
				break;
			}
			grid->yv[j-NY] = yv_;
		}

		// Make ghost cells same size as boundary cells
		dyv_l = grid -> yv[1]  - grid -> yv[0];
		dyv_r = grid -> yv[NY] - grid -> yv[NY-1];

		for (i = 1; i <= Ghost_Nodes; i++) {
			grid -> yv[NY+i] = grid -> yv[NY] + i * dyv_r;
			grid -> yv[0-i]  = grid -> yv[0]  - i * dyv_l;

			grid -> yc[NY+i] = grid -> yv[NY+i] + 0.5 * dyv_r;
			grid -> yc[0-i]  = grid -> yv[0-i]  - 0.5 * dyv_l;
		}

		// Make first ghost nodes equidistant from boundary
		grid->yc[NY] = grid->yc[NY-1] + 2.0 * (grid->yv[NY] - grid->yc[NY-1]);
		grid->yc[-1] = grid->yc[0]    - 2.0 * (grid->yc[0]  - grid->yv[0]);

	}
	fclose(file_in);

	/*------------------------------------------------------------------------*/
	/*
	 Read the grid in z-direction
	 */
	/*------------------------------------------------------------------------*/
	sprintf(filename, "Grid_z.inp");
	file_in  = fopen(filename, "r");
	if (file_in == NULL) {
		if (params->rank==0)
			printf("Warning!\nCould not open the grid file \"%s\". Using uniform grid formulation\n", filename);
		success = NO;
	}
	else {

		for (k = 0; k < NZ; k++) {

			status = fscanf(file_in, "%lf\n", &zc_);
			if (status < 1) {
				success = NO;
				break;
			}
			grid->zc[k] = zc_;
		}

		for (k = NZ; k < (2 * NZ + 1); k++) {

			status = fscanf(file_in, "%lf\n", &zw_);
			if (status < 1) {
				success = NO;
				break;
			}
			grid->zw[k-NZ] = zw_;
		}

#ifdef ZPERIODIC
		for (i = 1; i <= Ghost_Nodes; i++) {
			grid -> zc[NZ+i] = grid -> zc[i];
			grid -> zc[0-i]  = grid -> zc[NZ-i];

			grid -> zw[NZ+i] = grid -> zw[i];
			grid -> zw[0-i]  = grid -> zw[NZ-i];
		}
		grid->zc[NZ] = grid->zc[0];
#else
		// Make ghost cells same size as boundary cells
		dzw_l = grid -> zw[1]  - grid -> zw[0];
		dzw_r = grid -> zw[NZ] - grid -> zw[NZ-1];

		for (i = 1; i <= Ghost_Nodes; i++) {
			grid -> zw[NZ+i] = grid -> zw[NZ] + i * dzw_r;
			grid -> zw[0-i]  = grid -> zw[0]  - i * dzw_l;

			grid -> zc[NZ+i] = grid -> zw[NZ+i] + 0.5 * dzw_r;
			grid -> zc[0-i]  = grid -> zw[0-i]  - 0.5 * dzw_l;
		}

		// Make first ghost nodes equidistant from boundary
		grid->zc[NZ] = grid->zc[NZ-1] + 2.0 * (grid->zw[NZ] - grid->zc[NZ-1]);
		grid->zc[-1] = grid->zc[0] - 2.0 * (grid->zc[0] - grid->zw[0]);
#endif

	}
	fclose(file_in);

	return success;
}




/******************************************************************************/
/*
 This function imports the cell-center bottom interface exact location. This is
 for the cases where there is not exact analytical function for describing the
 bottom interface. Format would be:

     x_c z_c y_c

 */
/******************************************************************************/
int Grid_import_bottom_interface(MAC_grid *grid, Parameters *params, char *filename) {


	int success = YES;
	int status;
	FILE *file_in;
	int i, k;
	int NX, NZ;
	double xc_, zc_, y_interface;

	NX = params->NXM;
	NZ = params->NZM;

	file_in = fopen(filename, "r");
	if (file_in == NULL) {

		if (params->rank==0)
			printf("Warning!\nCould not open the bottom interface grid file \"%s\". Using the default interface\n", filename);
		success = NO;
	}
	else {

		for (i=0; i<NX; i++) {

			for (k=0; k<NZ; k++) {

				status = fscanf(file_in, "%lf %lf %lf\n", &xc_, &zc_, &y_interface);
				if (status < 1) {
					success = NO;
					break;
				}
				grid->interface_position[k][i] = 30.0*y_interface;
			} /* for k*/
			if (success == NO) {
				break;
			}
		} /* for i */

		fclose(file_in);
	} /* else */

	return success;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Grid_tag_q_nodes_using_surface(MAC_grid *grid, Parameters *params, char which_quantity) {

	double ***q_sdf=NULL;
	int NI, NJ, NK;
	int i, j, k;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int index ;
	int ***status;
	int si_start, si_end;
	int sj_start, sj_end;
	int sk_start, sk_end;

	// Number of grid points excluding the half cell added
	NI = grid->NI;
	NJ = grid->NJ;
	NK = grid->NK;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	if (which_quantity == 'u') {

		// start and end of the for-loops
		// The other nodes are the box boudarieds
		i_start = max(1, Is);  // i = 0, BOUNDARY node
		i_end   = min(Ie, NI); // i = NI, BOUDARY node

		j_start = Js;          // j=0, could be FLUID, IMMERSED or SOLID node
		j_end   = min(Je, NJ); // j=NJ, BOUNDARY node

		k_start = Ks;          // k=0, could be FLUID, IMMERSED or SOLID node
		k_end   = min(Ke, NK); // k=NK, BOUNDARY node

		// signed distance function: distance from the solid interface
		q_sdf = grid->u_sdf;

		status = grid->u_status;

	}
	else if (which_quantity == 'v') {

		// start and end of the for-loops
		// The other nodes are the box boudarieds
		i_start = Is;          // i = 0, BOUNDARY node
		i_end   = min(Ie, NI); // i = NI, BOUDARY node

		j_start = max(1, Js);  // j=0, BOUNDARY node
		j_end   = min(Je, NJ); // j=NJ, BOUNDARY node

		k_start = Ks;          // k=0, could be FLUID, IMMERSED or SOLID node
		k_end   = min(Ke, NK); // k=NK, BOUNDARY node

		// signed distance function: distance from the solid interface
		q_sdf = grid->v_sdf;

		status = grid->v_status;

	}
	else if (which_quantity == 'w') {

		// start and end of the for-loops
		// The other nodes are the box boudarieds
		i_start = Is;          // i = 0, BOUNDARY node
		i_end   = min(Ie, NI); // i = NI, BOUDARY node

		j_start = Js;          // j=0, could be FLUID, IMMERSED or SOLID node
		j_end   = min(Je, NJ); // j=NJ, BOUNDARY node

		k_start = max(1, Ks);  // k=0,  BOUNDARY node
		k_end   = min(Ke, NK); // k=NK, BOUNDARY node

		// signed distance function: distance from the solid interface
		q_sdf = grid->w_sdf;

		status = grid->w_status;

	}
	else if (which_quantity == 'c') {

		// start and end of the for-loops
		// The other nodes are the box boudarieds
		i_start = Is;          // i = 0, could be FLUID, IMMERSED or SOLID node
		i_end   = min(Ie, NI); // i = NI, BOUDARY node

		j_start = Js;          // j=0, could be FLUID, IMMERSED or SOLID node
		j_end   = min(Je, NJ); // j=NJ, BOUNDARY node

		k_start = Ks;          // k=0,  could be FLUID, IMMERSED or SOLID node
		k_end   = min(Ke, NK); // k=NK, BOUNDARY node

		// signed distance function: distance from the solid interface
		q_sdf = grid->c_sdf;

		status = grid->c_status;

	}
	else {
		if (params->rank==0)
			printf("Grid.c/ Could not start tagging the nodes. Unknown quantity\n");
	}


	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++) {

				if ( (q_sdf[k][j][i] < 0.0)  ) {
					status[k][j][i] = SOLID;
				}
				else {
					status[k][j][i] = FLUID;
				}
			} // for i
		} // for j
	} // for k

	si_start = i_start;
	if (si_start == 0) si_start = 1;
	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=si_start; i<i_end; i++) {

				if (status[k][j][i] == FLUID) {
					if ( q_sdf[k][j][i-1] < 0.0 ) {
						status[k][j][i] = IMMERSED;
					}
				}
			} // for i
		} // for j
	} // for k

	si_end = i_end;
	if (si_end == NI) si_end -= 1;
	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<si_end; i++) {

				if (status[k][j][i] == FLUID) {
					if ( q_sdf[k][j][i+1] < 0.0 ) {
						status[k][j][i] = IMMERSED;
					}
				}
			} // for i
		} // for j
	} // for k

	sj_start = j_start;
	if (sj_start == 0) sj_start = 1;
	for (k=k_start; k<k_end; k++) {
		for (j=sj_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++) {

				if (status[k][j][i] == FLUID) {
					if ( q_sdf[k][j-1][i] < 0.0 ) {
						status[k][j][i] = IMMERSED;
					}
				}
			} // for i
		} // for j
	} // for k

	sj_end = j_end;
	if (sj_end == NJ) sj_end -= 1;
	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<sj_end; j++) {
			for (i=i_start; i<i_end; i++) {

				if (status[k][j][i] == FLUID) {
					if ( q_sdf[k][j+1][i] < 0.0 ) {
						status[k][j][i] = IMMERSED;
					}
				}
			} // for i
		} // for j
	} // for k

	sk_start = k_start;
	if (sk_start == 0) sk_start = 1;
	for (k=sk_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++) {

				if (status[k][j][i] == FLUID) {
					if ( q_sdf[k-1][j][i] < 0.0 ) {
						status[k][j][i] = IMMERSED;
					}
				}
			} // for i
		} // for j
	} // for k

	sk_end = k_end;
	if (sk_end == NK) sk_end -= 1;
	for (k=k_start; k<sk_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++) {

				if (status[k][j][i] == FLUID) {
					if ( q_sdf[k+1][j][i] < 0.0 ) {
						status[k][j][i] = IMMERSED;
					}
				}
			} // for i
		} // for j
	} // for k

}




/******************************************************************************/
/*
 This functions sets the BOUNDARY flag for the box boundaries. The last nodes at
 (x=Lx, y=Ly, z=Lz) (This is due to the fact that we have addded a half cell to
 the end of all the grids so they have equal numbers in all directions. Also,
 depedning on the grid (u, v or w), it would do the same for x=0, y=0, z=0
 */
/******************************************************************************/
void Grid_tag_q_box_boundary_nodes(MAC_grid *grid, char which_quantity, Parameters *params) {

	int i, j, k;
	int NX, NY, NZ;
	int ***status;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	// Total number of grid points in the domain, including the half cell added
	// to the very end
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	switch (which_quantity) {
		// pointer to functions
		case 'u':
			status = grid->u_status;
			break;

		case 'v':
			status = grid->v_status;
			break;

		case 'w':
			status = grid->w_status;
			break;

		case 'c':
			status = grid->c_status;
			break;

		default:
			printf("Grid.c/ Could not tag the box boundary nodes. Unknown quantity\n");
	} // switch

	// yz plane
	if ( (which_quantity == 'u') && (Is==0) ) {
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {

				status[k][j][0] = BOUNDARY;
			} // for j
		}// for k
	} // if

	if (Ie == NX ) {
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {

				if (which_quantity == 'u') {
#ifdef XPERIODIC
					if ( (j==NY-1) || (k==NZ-1) )
						status[k][j][NX-1] = BOUNDARY;
					else
						status[k][j][NX-1] = FLUID;
#else
					status[k][j][NX-1] = BOUNDARY;
#endif
				}
				else {
					status[k][j][NX-1] = BOUNDARY;
				}
			} // for j
		}// for k
	}

	// xz plane
	if ( (which_quantity == 'v')  && (Js==0) ){
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				status[k][0][i] = BOUNDARY;
			} // for i
		} // for k
	} // if

	if ( Je==NY ) {
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				if (which_quantity == 'v') {
#ifdef YPERIODIC
					if ( (i==NX-1) || (k==NZ-1) )
						status[k][NY-1][i] = BOUNDARY;
					else
						status[k][NY-1][i] = FLUID;
#else
					status[k][NY-1][i] = BOUNDARY;
#endif
				}
				else {
					status[k][NY-1][i] = BOUNDARY;
				}
			} // for i
		} // for k
	} // if

	// xy plane
	if ( (Ks==0) && (which_quantity == 'w') ) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				status[0][j][i] = BOUNDARY;
			} // for i
		} // for j
	} // if


	if  (Ke==NZ)  {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

				if (which_quantity == 'w') {
#ifdef ZPERIODIC
					if ( (j==NY-1) || (i==NX-1) )
						status[NZ-1][j][i] = BOUNDARY;
					else
						status[NZ-1][j][i] = FLUID;
#else
					status[NZ-1][j][i] = BOUNDARY;
#endif
				}
				else {
					status[NZ-1][j][i] = BOUNDARY;
				}
			} // for i
		} // for j
	} // if

}




/******************************************************************************/
/*
 This function tags all the nodes in the domain based on the bottom surface
 */
/******************************************************************************/
void Grid_identify_geometry(Cart3d_bag *data_bag) {

	int pnodes;
	double ***status;
	int i, j, k;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	status =  Memory_allocate_flow_variable(grid, params);

	int G_Is = grid->G_Is;
	int G_Js = grid->G_Js;
	int G_Ks = grid->G_Ks;

	int G_Ie = grid->G_Ie;
	int G_Je = grid->G_Je;
	int G_Ke = grid->G_Ke;

	int L_Is = grid->L_Is;
	int L_Js = grid->L_Js;
	int L_Ks = grid->L_Ks;

	int L_Ie = grid->L_Ie;
	int L_Je = grid->L_Je;
	int L_Ke = grid->L_Ke;

/*
	// First, tag the box boundary nodes
	Grid_tag_q_box_boundary_nodes(grid, 'u') ;
	Grid_tag_q_box_boundary_nodes(grid, 'v') ;
	Grid_tag_q_box_boundary_nodes(grid, 'w') ;
	Grid_tag_q_box_boundary_nodes(grid, 'c') ;
*/

	// Tag all the nodes based on the location of the inteface for u, v, w and
	// c grid nodes
	Grid_tag_q_nodes_using_surface(grid, params, 'u');
	Grid_tag_q_nodes_using_surface(grid, params, 'v');
	Grid_tag_q_nodes_using_surface(grid, params, 'w');
	Grid_tag_q_nodes_using_surface(grid, params, 'c');

	pnodes = params->ghost_nodes;
	for (k=G_Ks; k<G_Ke; k++) {
		for (j=G_Js; j<G_Je; j++) {
			for (i=G_Is; i<G_Ie; i++) {
				status[k][j][i] = grid->u_status[k][j][i];
			}
		}
	}
	Communication_update_ghost_nodes_flow_variable(status, 'u', pnodes, data_bag);
	for (k=L_Ks; k<L_Ke; k++) {
		for (j=L_Js; j<L_Je; j++) {
			for (i=L_Is; i<L_Ie; i++) {
				grid->u_status[k][j][i] = status[k][j][i];
			}
		}
	}

	for (k=G_Ks; k<G_Ke; k++) {
		for (j=G_Js; j<G_Je; j++) {
			for (i=G_Is; i<G_Ie; i++) {
				status[k][j][i] = grid->v_status[k][j][i];
			}
		}
	}
	Communication_update_ghost_nodes_flow_variable(status, 'v', pnodes, data_bag);
	for (k=L_Ks; k<L_Ke; k++) {
		for (j=L_Js; j<L_Je; j++) {
			for (i=L_Is; i<L_Ie; i++) {
				grid->v_status[k][j][i] = status[k][j][i];
			}
		}
	}

	for (k=G_Ks; k<G_Ke; k++) {
		for (j=G_Js; j<G_Je; j++) {
			for (i=G_Is; i<G_Ie; i++) {
				status[k][j][i] = grid->w_status[k][j][i];
			}
		}
	}
	Communication_update_ghost_nodes_flow_variable(status, 'w', pnodes, data_bag);
	for (k=L_Ks; k<L_Ke; k++) {
		for (j=L_Js; j<L_Je; j++) {
			for (i=L_Is; i<L_Ie; i++) {
				grid->w_status[k][j][i] = status[k][j][i];
			}
		}
	}

	for (k=G_Ks; k<G_Ke; k++) {
		for (j=G_Js; j<G_Je; j++) {
			for (i=G_Is; i<G_Ie; i++) {
				status[k][j][i] = grid->c_status[k][j][i];
			}
		}
	}
	Communication_update_ghost_nodes_flow_variable(status, 'c', pnodes, data_bag);
	for (k=L_Ks; k<L_Ke; k++) {
		for (j=L_Js; j<L_Je; j++) {
			for (i=L_Is; i<L_Ie; i++) {
				grid->c_status[k][j][i] = status[k][j][i];
			}
		}
	}

	Grid_tag_q_box_boundary_nodes(grid, 'u', params) ;
	Grid_tag_q_box_boundary_nodes(grid, 'v', params) ;
	Grid_tag_q_box_boundary_nodes(grid, 'w', params) ;
	Grid_tag_q_box_boundary_nodes(grid, 'c', params) ;
	Grid_set_interface_y_index(grid, params);

 	Memory_free_flow_variable(grid, params, status);


}




/******************************************************************************/
/*
 This function finds the index of the grid for the given x value.

 Note that the returned index (grid) is less than the given input position
 */
/******************************************************************************/
int Grid_get_x_index(double x, MAC_grid *grid, char which_quantity) {

	int index = -1;
	double *xq, *yq;
	int i, NI;

	NI = grid->NI;

	if (which_quantity == 'u')
		xq = grid->xu;
	else
		xq = grid->xc;

	for (i = 0; i < NI; i++) {

		if ( (x >= xq[i]) && (x <= xq[i+1]) ) {
			index = i;
			break;
		} // if
	}
	return index;
}




/******************************************************************************/
/*
 This function finds the index of the grid for the given y value.

 Note that the returned index (grid) is less than the given input position
 */
/******************************************************************************/
int Grid_get_y_index(double y, MAC_grid *grid, char which_quantity) {

	int index = -1;
	double *yq;
	int j, NJ;

	NJ = grid->NJ;

	if (which_quantity == 'v')
		yq = grid->yv;
	else
		yq = grid->yc;

	for (j = 0; j < NJ; j++) {

		if ( (y >= yq[j]) && (y <= yq[j+1]) ) {
			index = j;
			break;
		} // if
	}
	return index;
}




/******************************************************************************/
/*
 This function finds the index of the grid for the given y value.

 Note that the returned index (grid) is less than the given input position
 */
/******************************************************************************/
int Grid_get_z_index(double z, MAC_grid *grid, char which_quantity) {

	int index = -1;
	double *zq;
	int k, NK;

	NK = grid->NK;

	if (which_quantity == 'w')
		zq = grid->zw;
	else
		zq = grid->zc;

	for (k = 0; k < NK; k++) {

		if ( (z >= zq[k]) && (z <= zq[k+1]) ) {
			index = k;
			break;
		} // if
	}
	return index;
}




/******************************************************************************/
/*
 This function assumes that we only have a bottom boundry and sets the first
 y-index correspoding to the IMMERSED node to it. If there is no surface (first
 interior node is a FLUID node, then it sets j=0. If surface does not lay on the
 current processor, it sets it to -1
 */
/******************************************************************************/
void Grid_set_interface_y_index(MAC_grid *grid, Parameters *params) {

	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i, j, k;


	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k = Ks; k < Ke; k++) {
		for (i = Is; i < Ie; i++) {

			// set the index of first IMMERSED (or FLUID for no solid boundary
			// on the bottom) node
			grid->interface_y_index[k][i] = -1;

			for (j=Js; j<Je; j++) {

				if (grid->c_status[k][j][i]  == IMMERSED) {

					grid->interface_y_index[k][i] = j;
					break;
				} else {
					if ( (j == 0) && (grid->c_status[k][j][i] == FLUID) ) {

						grid->interface_y_index[k][i] = j;
						break;
					} // if
				} // if
			} // for j
		} // for i
	} // for k


}

/******************************************************************************/
/*
 This function finds the optimal number of blocks in each direction
 */
/******************************************************************************/
void Grid_find_best_NP(int *NP, double *err_NP, int first_sign, int second_sign, Parameters *params) {

	int NX = params -> NXM + 1;
	int NY = params -> NYM + 1;
	int NZ = params -> NZM + 1;

	int NT_domain, N_coord_length;
	int NPX,NPY,NPZ;
	int NREMYZ, NREMY;
	double f_coord_length;

	NT_domain = NX * NY * NZ;
	NT_domain = (NT_domain + params -> size - 1) / params -> size;
	f_coord_length = pow(NT_domain * 1.0, 1.0 / 3.0);
	N_coord_length = f_coord_length + 0.999;

	NPX = (NX + (N_coord_length - 1) ) / N_coord_length;
	NPX  = min(NPX, params -> size);
	NREMYZ = params -> size / NPX;
	if (NREMYZ * NPX != params -> size) {
		do {
			NPX = NPX + first_sign;
			NREMYZ = params -> size / NPX;
		} while (NREMYZ * NPX != params -> size);
	}
	NPZ = ( NZ + ( N_coord_length - 1 ) ) / N_coord_length;
	NPZ  = min(NPZ, NREMYZ);
	NREMY = NREMYZ / NPZ;
	if (NREMY * NPZ != NREMYZ) {
		do {
			NPZ = NPZ + second_sign;
			NREMY = NREMYZ / NPZ;
		} while (NREMY * NPZ != NREMYZ);
	}
	NPY = NREMY;

	NP[0] = NPX;
	NP[1] = NPY;
	NP[2] = NPZ;

	*err_NP = sqrt( pow(((double) NX / NPX) - f_coord_length , 2.0) +  pow(((double) NY / NPY) - f_coord_length , 2.0) +  pow(((double) NZ / NPZ) - f_coord_length , 2.0));

}
