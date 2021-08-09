#include "definitions.h"
#include "Communication.h"
#include "DataTypes.h"
#include "Grid.h"
#include "Immersed.h"
#include "MyMath.h"
#include "Memory.h"
#include "Display.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#undef IMMERSED_DIAGNOSTIC

MPI_COMM comm3d;

/******************************************************************************/
/*
 This function allocates memory for the immersed structure based on the given
 number of the grids
 */
/******************************************************************************/
Immersed *Immersed_create(MAC_grid *grid, Parameters *params, char which_quantity) {

	int N;
	Immersed *new_immersed;
	int NZ, NY, NX;
	int ierr;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	new_immersed = (Immersed *)malloc(sizeof(Immersed));
	Memory_check_allocation(new_immersed);

	N = Immersed_count_immersed_nodes(grid, params, which_quantity);
//	if (params->rank==0)
//		printf("Found %d immeresed points in the domain on the current processor for the %c-quantity\n", N, which_quantity);

	// Maximum number of immersed nodes in the whole domain
	new_immersed->N = N;

	// Allocate the memory for each individual immersed node
	new_immersed->ib_nodes = (ImmersedNode *)malloc(N*sizeof(ImmersedNode));
	Memory_check_allocation(new_immersed->ib_nodes);

	// u, v, c (or p)
	new_immersed->quantity = which_quantity;

	return new_immersed;
}




/******************************************************************************/
/*
 This function releases the memory allocated for the immersed structure
 */
/******************************************************************************/
void Immersed_destroy(Immersed *q_immersed, Parameters *params, MAC_grid *grid) {

	int ierr;

	free(q_immersed->ib_nodes);

	free(q_immersed);
}




/******************************************************************************/
/*
 This function allocates memory for the immersed structure based on the given
 number of the grids
 */
/******************************************************************************/
void Immersed_copy_to_nut(MAC_grid *grid, Parameters *params) {

	int N;
	Immersed *nut_immersed;
	Immersed *c_immersed;
	ImmersedNode *ib_cnode;
	ImmersedNode *ib_nutnode;
	int g, ibm;
	int ierr;
	double d1, d2, factor;

	c_immersed = grid->c_immersed;
	nut_immersed = (Immersed *)malloc(sizeof(Immersed));
	Memory_check_allocation(nut_immersed);

	N = grid->c_immersed->N;
	// Maximum number of immersed nodes in the whole domain
	nut_immersed->N = N;

	// Allocate the memory for each individual immersed node
	nut_immersed->ib_nodes = (ImmersedNode *)malloc(N*sizeof(ImmersedNode));
	Memory_check_allocation(nut_immersed->ib_nodes);

	// u, v, c (or p)
	nut_immersed->quantity = 'c';

	nut_immersed->boundary_condition =  DIRICHLET;

	for (g=0; g<N; g++) {

		// Get the current immersed node
		ib_cnode = Immersed_get_ib_node(c_immersed, g);
		ib_nutnode = Immersed_get_ib_node(nut_immersed, g);

		ib_nutnode->im_point.x = ib_cnode->im_point.x;
		ib_nutnode->im_point.y = ib_cnode->im_point.y;
		ib_nutnode->im_point.z = ib_cnode->im_point.z;

		ib_nutnode->boundary_point.x = ib_cnode->boundary_point.x;
		ib_nutnode->boundary_point.y = ib_cnode->boundary_point.y;
		ib_nutnode->boundary_point.z = ib_cnode->boundary_point.z;

		ib_nutnode->intersection_point.x = ib_cnode->intersection_point.x;
		ib_nutnode->intersection_point.y = ib_cnode->intersection_point.y;
		ib_nutnode->intersection_point.z = ib_cnode->intersection_point.z;

		ib_nutnode->n.vx = ib_cnode->n.vx;
		ib_nutnode->n.vy = ib_cnode->n.vy;
		ib_nutnode->n.vz = ib_cnode->n.vz;

		// Index of the current ib_node
		ib_nutnode->im_index.x_index = ib_cnode->im_index.x_index;
		ib_nutnode->im_index.y_index = ib_cnode->im_index.y_index;
		ib_nutnode->im_index.z_index = ib_cnode->im_index.z_index;

		for (ibm=0;ibm<3;ibm++) {
			ib_nutnode->fluid_index[ibm].x_index = ib_cnode->fluid_index[ibm].x_index;
			ib_nutnode->fluid_index[ibm].y_index = ib_cnode->fluid_index[ibm].y_index;
			ib_nutnode->fluid_index[ibm].z_index = ib_cnode->fluid_index[ibm].z_index;

			ib_nutnode->fluid_coef[ibm] = ib_cnode->fluid_coef[ibm];
		}

		d1  = (ib_nutnode->intersection_point.x - ib_nutnode->im_point.x)*(ib_nutnode->intersection_point.x - ib_nutnode->im_point.x);
		d1 += (ib_nutnode->intersection_point.y - ib_nutnode->im_point.y)*(ib_nutnode->intersection_point.y - ib_nutnode->im_point.y);
		d1 += (ib_nutnode->intersection_point.z - ib_nutnode->im_point.z)*(ib_nutnode->intersection_point.z - ib_nutnode->im_point.z);
		d1 = sqrt(d1);

		d2  = (ib_nutnode->boundary_point.x - ib_nutnode->im_point.x)*(ib_nutnode->boundary_point.x - ib_nutnode->im_point.x);
		d2 += (ib_nutnode->boundary_point.y - ib_nutnode->im_point.y)*(ib_nutnode->boundary_point.y - ib_nutnode->im_point.y);
		d2 += (ib_nutnode->boundary_point.z - ib_nutnode->im_point.z)*(ib_nutnode->boundary_point.z - ib_nutnode->im_point.z);
		d2 = sqrt(d2);

		factor = d2/(d1+d2);
		ib_nutnode->fluid_coef[0] *= factor;
		ib_nutnode->fluid_coef[1] *= factor;
		ib_nutnode->fluid_coef[2] *= factor;
	} // for g

	grid->nut_immersed = nut_immersed;
	return;
}




/******************************************************************************/
/*
 This function counts the total number of immeresed nodes for each quantity
 */
/******************************************************************************/
int Immersed_count_immersed_nodes(MAC_grid *grid, Parameters *params, char which_quantity) {

	int i, j, k;
	int ***status;
	int N_immersed = 0;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	switch (which_quantity) {

		case 'u':
			status = grid->u_status; break;
		case 'v':
			status = grid->v_status; break;
		case 'w':
			status = grid->w_status; break;
		case 'c':
			status = grid->c_status; break;
		default:
			if (params->rank==0) printf("Can not count immersed nodes. Invalid quantity\n");
	} // switch

	// Go through local part on each processor. Count number of immersed nodes
	// on the current processor
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

				if (status[k][j][i] == IMMERSED ) {

					N_immersed++;
				} // if
			} // for i
		} // for j
	} // for k
	if (N_immersed > 0) printf("rank = %d Immersed_node = %d quantity = %c\n",params->rank, N_immersed, which_quantity);
	fflush(stdout);

	return N_immersed;
}




/******************************************************************************/
/*
 This function sets both global and xyz indices of the immersed nodes in the
 domain on the current processor

 TBC

 */
/******************************************************************************/
void Immersed_set_q_immersed_indices(MAC_grid *grid, Parameters *params, char which_quantity) {

	int i, j, k;
	int ***status;
	Immersed *q_immersed;
	int g_index = 0;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	switch (which_quantity) {

		case 'u':
			status = grid->u_status;
			q_immersed  = grid->u_immersed;
			break;

		case 'v':
			status = grid->v_status;
			q_immersed  = grid->v_immersed;
			break;

		case 'w':
			status = grid->w_status;
			q_immersed  = grid->w_immersed;
			break;

		case 'c':
			status = grid->c_status;
			q_immersed  = grid->c_immersed;
			break;

		default:
			if (params->rank==0) printf("Can not count immersed nodes. Invalid quantity\n");
	} // switch

	// Go through all the nodes in the physical domain on the current processor
	// and set the physical indices for each node and also the the global index
	// of each immersed node. Global index is simply the number of the ib node
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

				if (status[k][j][i] == IMMERSED) {

					// Set the x and y indices of the current ib node
					(&q_immersed->ib_nodes[g_index])->im_index.x_index = i;
					(&q_immersed->ib_nodes[g_index])->im_index.y_index = j;
					(&q_immersed->ib_nodes[g_index])->im_index.z_index = k;

					g_index++;

				}
			} // for i
		} // for j
	} // for k
}




/******************************************************************************/
/*
 This function sets the coordinates of the immersed nodes
 */
/******************************************************************************/
void Immersed_set_q_immersed_coordinates(MAC_grid *grid, Parameters *params,
		char which_quantity) {

	int g;
	int N;
	int i_im, j_im, k_im;
	Immersed *q_immersed;
	ImmersedNode *ib_node;
	double *xq, *yq, *zq;

	switch (which_quantity) {

		case 'u':

			xq = grid->xu;
			yq = grid->yc;
			zq = grid->zc;
			q_immersed = grid->u_immersed;
			break;

		case 'v':

			xq = grid->xc;
			yq = grid->yv;
			zq = grid->zc;
			q_immersed = grid->v_immersed;
			break;

		case 'w':

			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zw;
			q_immersed = grid->w_immersed;
			break;

		case 'c':

			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zc;
			q_immersed = grid->c_immersed;
			break;

		default:

			if (params->rank==0) printf("Immersed.c/ Error setting coordinates. Unknown quantity\n");
	} // switch

	// Total number of immersed nodes on current processor
	N = q_immersed->N;

	// Go through all the immersed nodes and set the coordinates of the
	// immersed nodes
	for (g=0; g<N; g++) {

		// Get the current immersed node
		ib_node = Immersed_get_ib_node(q_immersed, g);

		// Index of the current ib_node
		i_im = ib_node->im_index.x_index;
		j_im = ib_node->im_index.y_index;
		k_im = ib_node->im_index.z_index;

		// Set the coordinates of the current ib_node based on the grid
		// position of the given quantity
		ib_node->im_point.x = xq[i_im];
		ib_node->im_point.y = yq[j_im];
		ib_node->im_point.z = zq[k_im];

	} // for g
}




/******************************************************************************/
/*
 This function sets up all the necessary functions for the immersed nodes
 */
/******************************************************************************/
void Immersed_setup_q_immersed_nodes(MAC_grid *grid, Parameters *params,
		char which_quantity) {

	Immersed *q_immersed;
	int boundary_condition;

	// First, allocate memory for the q_immersed nodes
	q_immersed = Immersed_create(grid, params, which_quantity);

	// Pass the pointer to the grid structure
	switch (which_quantity) {

		case 'u':

			grid->u_immersed = q_immersed;
			q_immersed->boundary_condition =  DIRICHLET;
			break;

		case 'v':

			grid->v_immersed = q_immersed;
			q_immersed->boundary_condition =  DIRICHLET;
			break;

		case 'w':

			grid->w_immersed = q_immersed;
			q_immersed->boundary_condition =  DIRICHLET;
			break;

		case 'c':

			grid->c_immersed = q_immersed;
			q_immersed->boundary_condition =  NEUMANN;
			break;

		default:
			if (params->rank==0) printf("Immeresed.c /Can not setup immersed nodes. Unknown quantity %c.\n", which_quantity);
	} // switch

	// Set the xyz and global indices of the immersed nodes in the physical
	// domain
	Immersed_set_q_immersed_indices(grid, params, which_quantity);

	// Set the coordinates of the immersed nodes
	Immersed_set_q_immersed_coordinates(grid, params, which_quantity);

	// Find the unit normal vector to the surface and also the boundary point
	Immersed_find_boundary_point_and_normal_vector(grid, params, which_quantity);

	// Find the interpolation stencil and the interpolation coefficients
	Immersed_set_q_immersed_interpolation(grid, params, which_quantity);

}




/******************************************************************************/
/*
 This function finds the fluid neighboring nodes for each immersed node. Stores
 the indices and the state of the node comparing to the current node
 */
/******************************************************************************/
void Immersed_set_q_immersed_interpolation(MAC_grid *grid, Parameters *params,
		char which_quantity) {

	int g, N;
	Immersed *q_immersed;
	ImmersedNode *ib_node;
	int success;
	FILE *fid, *fid1;
	char bin_filename[50];
	int ***status;

	switch (which_quantity) {

		case 'u':
			q_immersed = grid->u_immersed;
			status = grid->u_status;
			break;

		case 'v':
			q_immersed = grid->v_immersed;
			status = grid->v_status;
			break;

		case 'w':
			q_immersed = grid->w_immersed;
			status = grid->w_status;
			break;

		case 'c':
			q_immersed = grid->c_immersed;
			status = grid->c_status;
			break;

		default:

			if (params->rank==0) printf("Immersed.c/ Could not set the neighboring immersed nodes. Unknown quantity\n");
	}  // switch

	// Total number of the immersed nodes
	N = q_immersed->N;

	// Go through all the ib nodes and find the neighboring fluid nodes. Store
	// the indices and the coordinates
	for (g=0; g<N; g++) {

#ifdef IMMERSED_DIAGNOSTIC
	sprintf(bin_filename,"intercoeff%c_%d.dat",which_quantity, params->rank);
	fid1 = fopen(bin_filename,"a");
	sprintf(bin_filename,"allnodes%c_%d.dat",which_quantity, params->rank);
	fid = fopen(bin_filename,"a");
#endif

		// Get the current immersed node
		ib_node = Immersed_get_ib_node(q_immersed, g);

		// Set the fluid nodes neighboring the current ib node
		success = Immersed_compute_nearest_plane_interpolation(ib_node, grid, params, q_immersed->boundary_condition, which_quantity);
//		fprintf(fid, " \n after nearest plane call success = %d  i = %d j = %d k = %d point =(%f %f %f) \n",
//				success, ib_node->im_index.x_index, ib_node->im_index.y_index,
//				ib_node->im_index.z_index, ib_node->im_point.x,
//				ib_node->im_point.y, ib_node->im_point.z);

		if (!success) {

#ifdef IMMERSED_DIAGNOSTIC
			fprintf(fid, " first !success i = %d j = %d k = %d point =(%f %f %f) \n",
					ib_node->im_index.x_index, ib_node->im_index.y_index,
					ib_node->im_index.z_index, ib_node->im_point.x,
					ib_node->im_point.y,ib_node->im_point.z);
#endif
			success = Immersed_compute_intersecting_xyz_plane_interpolation(ib_node,
					grid, params, q_immersed->boundary_condition, which_quantity, 1);

		}
		if (!success) {

			success = Immersed_compute_intersecting_xyz_plane_interpolation(ib_node,
					grid, params, q_immersed->boundary_condition, which_quantity, 2);
#ifdef IMMERSED_DIAGNOSTIC
			fprintf(fid, "second !success i = %d j = %d k = %d point =(%f %f %f) \n",
					ib_node->im_index.x_index, ib_node->im_index.y_index,
					ib_node->im_index.z_index, ib_node->im_point.x,
					ib_node->im_point.y, ib_node->im_point.z);
			if (success) fprintf(fid, "third success \n");
#endif
		}
		if (!success) {

#ifdef IMMERSED_DIAGNOSTIC
			fprintf(fid, "Third !success i = %d j = %d k = %d point =(%f %f %f) \n",
					ib_node->im_index.x_index, ib_node->im_index.y_index,
					ib_node->im_index.z_index, ib_node->im_point.x,
					ib_node->im_point.y, ib_node->im_point.z);
#endif
			printf("Interpolation stencil NOT FOUND for (%f %f %f) at index (%d %d %d) for %c \n",ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z,
			ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index, which_quantity);
		}

#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid1,"inter_node(%4d %3d %3d) node1(%4d %3d %3d) node2(%4d %3d %3d) node3(%4d %3d %3d) coef1= %+-8.5f coef2= %+-8.5f coef3= %+-8.5f tcoef= %+-8.5f node1= %1d node2= %1d node3= %1d \n",
			ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index,
			ib_node->fluid_index[0].x_index, ib_node->fluid_index[0].y_index, ib_node->fluid_index[0].z_index,
			ib_node->fluid_index[1].x_index, ib_node->fluid_index[1].y_index, ib_node->fluid_index[1].z_index,
			ib_node->fluid_index[2].x_index, ib_node->fluid_index[2].y_index, ib_node->fluid_index[2].z_index,
			ib_node->fluid_coef[0],ib_node->fluid_coef[1], ib_node->fluid_coef[2],
			ib_node->fluid_coef[0]+ib_node->fluid_coef[1]+ib_node->fluid_coef[2],
			status[ib_node->fluid_index[0].z_index][ib_node->fluid_index[0].y_index][ib_node->fluid_index[0].x_index],
			status[ib_node->fluid_index[1].z_index][ib_node->fluid_index[1].y_index][ib_node->fluid_index[1].x_index],
			status[ib_node->fluid_index[2].z_index][ib_node->fluid_index[2].y_index][ib_node->fluid_index[2].x_index]);

		fclose(fid);
		fclose(fid1);
#endif
	} // for g
}




/******************************************************************************/
/*
 This function finds the nodes used in interpolation stencil and computes the
 interpolation coefficient. The interpolation stencil consists of fluid points
 in the nearest plane, when the normal from the immersed boundary node
 intersects this plane
 */
/******************************************************************************/
int Immersed_compute_nearest_plane_interpolation(ImmersedNode *ib_node,
		MAC_grid *grid, Parameters *params, int boundary_condition,
		char which_quantity) {

	int n;
	int i, j, k;
	int i_immersed, j_immersed, k_immersed;
	int ***status;
	double *xq, *yq, *zq;
	Indices ind_neighbor[27];
	double dist_neighbor[27];
	double dist;
	int NX, NY, NZ;
	int i_minus, j_minus, k_minus;
	int i_plus, j_plus, k_plus;
	int count, switched, size;
	double hold;
	Indices ind_hold;
	PointType *plane, *line;
	double d1, d2;
	int found, success;
	double t, u, v, factor;
	char bin_filename[50];
	char fname3[50];
	FILE *fid;

#ifdef IMMERSED_DIAGNOSTIC
	sprintf(bin_filename,"interpolation%c_%d.dat",which_quantity, params->rank);
	fid = fopen(bin_filename,"a");
#endif


	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	switch (which_quantity) {

		case 'u':
			status = grid->u_status;
			xq = grid->xu;
			yq = grid->yc;
			zq = grid->zc;
			break;

		case 'v':
			status = grid->v_status;
			xq = grid->xc;
			yq = grid->yv;
			zq = grid->zc;
			break;

		case 'w':
			status = grid->w_status;
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zw;
			break;

		case 'c':
			status = grid->c_status;
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zc;
			break;

		default:

			if (params->rank==0) printf("Immersed.c/ Could not set the neighboring immersed nodes. Unknown quantity\n");
	}  // switch

	plane = (PointType *) malloc(3*sizeof(PointType));
	line  = (PointType *) malloc(2*sizeof(PointType));

	// immersed node index
	i_immersed = ib_node->im_index.x_index;
	j_immersed = ib_node->im_index.y_index;
	k_immersed = ib_node->im_index.z_index;


	i_minus = i_immersed - 1;
	i_plus  = i_immersed + 1;
	j_minus = j_immersed - 1;
	j_plus  = j_immersed + 1;
	k_minus = k_immersed - 1;
	k_plus  = k_immersed + 1;



#ifdef IMMERSED_DIAGNOSTIC
	fprintf(fid, "i = %d j = %d k = %d point =(%f %f %f) \n", i_immersed, j_immersed,
			k_immersed, ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z);
	fprintf(fid, "im = %d ip = %d jm = %d jp = %d km = %d kp = %d\n",i_minus,
			i_plus, j_minus, j_plus, k_minus, k_plus);
#endif

	if (i_minus < 0) {
		i_minus = 0;
		i_plus = 1;
	}
	if (j_minus < 0) {
		j_minus = 0;
		j_plus = 1;
	}
	if (k_minus < 0) {
		k_minus = 0;
		k_plus = 1;
	}

	if (i_plus == NX - 1) {
		i_plus = NX - 2;
	}
	if (j_plus == NY - 1) {
		j_plus = NY - 2;
	}
	if (k_plus == NZ - 1) {
		k_plus = NZ - 2;
	}



	// Compute the distance between the immersed node and the neighboring fluid
	// nodes that could be part of the interpolation stencil
	count = 0;
	for (k=k_minus;k<=k_plus;k++) {
		for (j=j_minus;j<=j_plus;j++) {
			for (i=i_minus;i<=i_plus;i++) {
				if ( ! ( (i==i_immersed) && (j==j_immersed) && (k==k_immersed) ) ){
					if (status[k][j][i] == FLUID) {
						ind_neighbor[count].x_index = i;
						ind_neighbor[count].y_index = j;
						ind_neighbor[count].z_index = k;
						dist  = (ib_node->im_point.x - xq[i]) * (ib_node->im_point.x - xq[i]);
						dist += (ib_node->im_point.y - yq[j]) * (ib_node->im_point.y - yq[j]);
						dist += (ib_node->im_point.z - zq[k]) * (ib_node->im_point.z - zq[k]);
						dist_neighbor[count] = sqrt(dist);
						count++;
#ifdef IMMERSED_DIAGNOSTIC
						fprintf(fid, "fluid i = %d j = %d  k = %d \n",i, j, k);
						fprintf(fid, "fluid x = %f y = %f z=%f dist = %f\n",xq[i], yq[j],zq[k],dist);
#endif
					}
				}
			}
		}
	}
#ifdef IMMERSED_DIAGNOSTIC
	fprintf(fid,"immersed point(%f, %f, %f) fluid neihbors = %d\n",ib_node->im_point.x,ib_node->im_point.y,ib_node->im_point.z,count);
#endif


	if (count < 3) {
#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid,"\n\n");
		fclose(fid);
#endif
		return 0;
	}
	size = count-1;

	// Sort the neighboring fluid nodes based distance from the immersed node
	// (using bubble sort)
	switched = 1;
	for(i = 0; i < size && switched; i++) {
		switched = 0;
		for(j = 0; j < size - i; j++) {
			if (dist_neighbor[j] > dist_neighbor[j+1]) {
				switched = 1;
				hold = dist_neighbor[j];
				dist_neighbor[j] = dist_neighbor[j + 1];
				dist_neighbor[j + 1] = hold;

				ind_hold.x_index = ind_neighbor[j].x_index;
				ind_hold.y_index = ind_neighbor[j].y_index;
				ind_hold.z_index = ind_neighbor[j].z_index;
				ind_neighbor[j].x_index = ind_neighbor[j+1].x_index;
				ind_neighbor[j].y_index = ind_neighbor[j+1].y_index;
				ind_neighbor[j].z_index = ind_neighbor[j+1].z_index;
				ind_neighbor[j+1].x_index = ind_hold.x_index;
				ind_neighbor[j+1].y_index = ind_hold.y_index;
				ind_neighbor[j+1].z_index = ind_hold.z_index;
			}
		}
	}

	line[0].x = ib_node->boundary_point.x;
	line[0].y = ib_node->boundary_point.y;
	line[0].z = ib_node->boundary_point.z;
	line[1].x = ib_node->im_point.x;
	line[1].y = ib_node->im_point.y;
	line[1].z = ib_node->im_point.z;
	plane[0].x = xq[ind_neighbor[0].x_index];
	plane[0].y = yq[ind_neighbor[0].y_index];
	plane[0].z = zq[ind_neighbor[0].z_index];
	plane[1].x = xq[ind_neighbor[1].x_index];
	plane[1].y = yq[ind_neighbor[1].y_index];
	plane[1].z = zq[ind_neighbor[1].z_index];
	plane[2].x = xq[ind_neighbor[2].x_index];
	plane[2].y = yq[ind_neighbor[2].y_index];
	plane[2].z = zq[ind_neighbor[2].z_index];
#ifdef IMMERSED_DIAGNOSTIC
	fprintf(fid, "boundary (%f %f %f) and immersed (%f %f %f) \n",
			line[0].x,line[0].y,line[0].z, line[1].x,line[1].y,line[1].z);
#endif

	found = 0;
	size = count;
	count = 2;
	do {
		success = Immersed_line_plane_intersection(line, plane, &t, &u, &v);
#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid, "plane (%f %f %f), (%f %f %f) and (%f %f %f) \n",
				plane[0].x,plane[0].y,plane[0].z,plane[1].x,plane[1].y,plane[1].z,
				plane[2].x,plane[2].y,plane[2].z);
		fprintf(fid,"t=%e u=%e v=%e count = %d success = %d\n",t,u,v,count,success);
#endif
		found = success;
		if ( count == size)
			found = 1;
		if ( (!success) && (count < size) ) {
			plane[2].x = xq[ind_neighbor[count].x_index];
			plane[2].y = yq[ind_neighbor[count].y_index];
			plane[2].z = zq[ind_neighbor[count].z_index];
			count++;
		}
	} while ( !found );


	if (!success) {
#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid,"Stencil not found\n");
		fprintf(fid,"\n");
		fclose(fid);
#endif
		return 0;
	}


	ib_node->fluid_index[0].x_index = ind_neighbor[0].x_index;
	ib_node->fluid_index[0].y_index = ind_neighbor[0].y_index;
	ib_node->fluid_index[0].z_index = ind_neighbor[0].z_index;

	ib_node->fluid_index[1].x_index = ind_neighbor[1].x_index;
	ib_node->fluid_index[1].y_index = ind_neighbor[1].y_index;
	ib_node->fluid_index[1].z_index = ind_neighbor[1].z_index;

	ib_node->fluid_index[2].x_index = ind_neighbor[count].x_index;
	ib_node->fluid_index[2].y_index = ind_neighbor[count].y_index;
	ib_node->fluid_index[2].z_index = ind_neighbor[count].z_index;

	ib_node->fluid_coef[0] = 1 - u - v;
	ib_node->fluid_coef[1] = u;
	ib_node->fluid_coef[2] = v;

	ib_node->intersection_point.x = line[0].x + (line[1].x-line[0].x)*t;
	ib_node->intersection_point.y = line[0].y + (line[1].y-line[0].y)*t;
	ib_node->intersection_point.z = line[0].z + (line[1].z-line[0].z)*t;
#ifdef IMMERSED_DIAGNOSTIC
	fprintf(fid, "u= %f v= %f \n",u,v);
	fprintf(fid, "coef1 = %f coef2 = %f coef3 = %f \n",ib_node->fluid_coef[0],
			ib_node->fluid_coef[1],ib_node->fluid_coef[2]);
	fprintf(fid,"xinter = %f yinter = %f zinter= %f \n",ib_node->intersection_point.x,
			ib_node->intersection_point.y, ib_node->intersection_point.z);
	fprintf(fid,"\n");
	fclose(fid);
#endif

	if (boundary_condition == DIRICHLET) {


		d1  = (ib_node->intersection_point.x - ib_node->im_point.x)*(ib_node->intersection_point.x - ib_node->im_point.x);
		d1 += (ib_node->intersection_point.y - ib_node->im_point.y)*(ib_node->intersection_point.y - ib_node->im_point.y);
		d1 += (ib_node->intersection_point.z - ib_node->im_point.z)*(ib_node->intersection_point.z - ib_node->im_point.z);
		d1 = sqrt(d1);

		d2  = (ib_node->boundary_point.x - ib_node->im_point.x)*(ib_node->boundary_point.x - ib_node->im_point.x);
		d2 += (ib_node->boundary_point.y - ib_node->im_point.y)*(ib_node->boundary_point.y - ib_node->im_point.y);
		d2 += (ib_node->boundary_point.z - ib_node->im_point.z)*(ib_node->boundary_point.z - ib_node->im_point.z);
		d2 = sqrt(d2);

		factor = d2/(d1+d2);
		ib_node->fluid_coef[0] *= factor;
		ib_node->fluid_coef[1] *= factor;
		ib_node->fluid_coef[2] *= factor;
	}
	return 1;

}




/******************************************************************************/
/*
 This function finds the nodes used in interpolation stencil and computes the
 interpolation coefficient.

 The interpolation stencil consists of fluid points in the nearest x-, y- or
 z-plane, when the normal from the immersed boundary node intersects this plane.

 This function is used when the "nearsest plane interpolation" does not work
 (for eg: They might not have enough fluid nodes)
 */
/******************************************************************************/
int Immersed_compute_intersecting_xyz_plane_interpolation(ImmersedNode *ib_node,
		MAC_grid *grid, Parameters *params, int boundary_condition,
		char which_quantity, int shift) {

	int n;
	int i, j, k;
	int i_immersed, j_immersed, k_immersed;
	int ***status;
	double *xq, *yq, *zq;
	Indices *ind_neighbor;
	PointType *x_plane, *y_plane, *z_plane;
	Indices *ind_x_plane, *ind_y_plane, *ind_z_plane;
	double dist_neighbor[4];
	double dist;
	int NX, NY, NZ;
	int i_target, j_target, k_target;
	double x_target, y_target, z_target;
	int i_minus, j_minus, k_minus;
	int i_plus, j_plus, k_plus;
	int count, switched, size;
	double hold;
	Indices ind_hold;
	PointType *plane, *line, *x_interxn, *y_interxn, *z_interxn;
	double d1, d2;
	int found, success;
	double t, u, v, factor;
	double x_dist, y_dist, z_dist, min_dist;
	int x_success, y_success, z_success, min_index;
	double x_u, x_v, y_u, y_v, z_u, z_v;

	i_immersed = ib_node->im_index.x_index;
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	switch (which_quantity) {

		case 'u':
			status = grid->u_status;
			xq = grid->xu;
			yq = grid->yc;
			zq = grid->zc;
			break;

		case 'v':
			status = grid->v_status;
			xq = grid->xc;
			yq = grid->yv;
			zq = grid->zc;
			break;

		case 'w':
			status = grid->w_status;
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zw;
			break;

		case 'c':
			status = grid->c_status;
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zc;
			break;

		default:

			if (params->rank==0)
				printf("Immersed.c/ Could not set the neighboring immersed nodes. Unknown quantity\n");
	}  // switch

	ib_node->fluid_coef[0] = 0.0;
	ib_node->fluid_coef[1] = 0.0;
	ib_node->fluid_coef[2] = 0.0;

	plane = (PointType *) malloc(3*sizeof(PointType));
	line  = (PointType *) malloc(2*sizeof(PointType));
	x_interxn  = (PointType *) malloc(sizeof(PointType));
	y_interxn  = (PointType *) malloc(sizeof(PointType));
	z_interxn  = (PointType *) malloc(sizeof(PointType));
	ind_neighbor = (Indices *) malloc(4*sizeof(Indices));
	x_plane = (PointType *) malloc(3*sizeof(PointType));
	y_plane = (PointType *) malloc(3*sizeof(PointType));
	z_plane = (PointType *) malloc(3*sizeof(PointType));
	ind_x_plane = (Indices *) malloc(3*sizeof(Indices));
	ind_y_plane = (Indices *) malloc(3*sizeof(Indices));
	ind_z_plane = (Indices *) malloc(3*sizeof(Indices));

	// immersed node index
	i_immersed = ib_node->im_index.x_index;
	j_immersed = ib_node->im_index.y_index;
	k_immersed = ib_node->im_index.z_index;


	if ( (ib_node->n.vx + 1e-8) > 0) {
		i_target = i_immersed + shift;
	}
	else {
		i_target = i_immersed - shift;
	}

	if ( (ib_node->n.vy + 1e-8) > 0) {
		j_target = j_immersed + shift;
	}
	else {
		j_target = j_immersed - shift;
	}

	if ( (ib_node->n.vz + 1e-8) > 0) {
		k_target = k_immersed + shift;
	}
	else {
		k_target = k_immersed - shift;
	}

	if (i_target < 0) {
		i_target = 0;
	}
	if (j_target < 0) {
		j_target = 0;
	}
	if (k_target < 0) {
		k_target = 0;
	}

	x_success = 0;
	y_success = 0;
	z_success = 0;
	x_dist = params->Lx;
	y_dist = params->Ly;
	z_dist = params->Lz;

	line[0].x = ib_node->boundary_point.x;
	line[0].y = ib_node->boundary_point.y;
	line[0].z = ib_node->boundary_point.z;
	line[1].x = ib_node->im_point.x;
	line[1].y = ib_node->im_point.y;
	line[1].z = ib_node->im_point.z;

	if ( fabs(ib_node->n.vx) > 1e-8 )  {

		x_target = xq[i_target];
		Immersed_line_xplane_intersection_point(line, x_target, x_interxn);

		j_minus = Grid_get_y_index(x_interxn[0].y, grid, which_quantity);
		j_plus = j_minus + 1;
		k_minus = Grid_get_z_index(x_interxn[0].z, grid, which_quantity);
		k_plus = k_minus + 1;

		/*--------------------------------------------------------------------*/
		/*
		 DO NOT change the order of the nodes below, as the
		 "find_plane_interpolation_coeff" assumes the following order of the
		 nodes:

		               1---------------2
		               |               |
		               |               |
		               |               |
		               0---------------3

		 */
		/*--------------------------------------------------------------------*/
		ind_neighbor[0].x_index = i_target;
		ind_neighbor[0].y_index = j_minus;
		ind_neighbor[0].z_index = k_minus;

		ind_neighbor[1].x_index = i_target;
		ind_neighbor[1].y_index = j_plus;
		ind_neighbor[1].z_index = k_minus;

		ind_neighbor[2].x_index = i_target;
		ind_neighbor[2].y_index = j_plus;
		ind_neighbor[2].z_index = k_plus;

		ind_neighbor[3].x_index = i_target;
		ind_neighbor[3].y_index = j_minus;
		ind_neighbor[3].z_index = k_plus;

		x_dist  = (ib_node->im_point.x - x_interxn[0].x) * (ib_node->im_point.x - x_interxn[0].x);
		x_dist += (ib_node->im_point.y - x_interxn[0].y) * (ib_node->im_point.y - x_interxn[0].y);
		x_dist += (ib_node->im_point.z - x_interxn[0].z) * (ib_node->im_point.z - x_interxn[0].z);

		x_success = Immersed_find_plane_interpolation_coeff(ind_neighbor,
		            x_interxn, 1, which_quantity, grid, x_plane, ind_x_plane,
		            &x_u, &x_v, i_immersed, j_immersed, k_immersed,  params);

	}

	if ( fabs(ib_node->n.vy) > 1e-8 )  {

		y_target = yq[j_target];
		Immersed_line_yplane_intersection_point(line, y_target, y_interxn);

		i_minus = Grid_get_x_index(y_interxn[0].x, grid, which_quantity);
		i_plus = i_minus + 1;
		k_minus = Grid_get_z_index(y_interxn[0].z, grid, which_quantity);
		k_plus = k_minus + 1;

		/*--------------------------------------------------------------------*/
		/*
		 DO NOT change the order of the nodes below, as the
		 "find_plane_interpolation_coeff" assumes the following order of the
		 nodes:

		 1---------------2
		 |               |
		 |               |
		 |               |
		 0---------------3

		 */
		/*--------------------------------------------------------------------*/
		ind_neighbor[0].x_index = i_minus;
		ind_neighbor[0].y_index = j_target;
		ind_neighbor[0].z_index = k_minus;

		ind_neighbor[1].x_index = i_plus;
		ind_neighbor[1].y_index = j_target;
		ind_neighbor[1].z_index = k_minus;

		ind_neighbor[2].x_index = i_plus;
		ind_neighbor[2].y_index = j_target;
		ind_neighbor[2].z_index = k_plus;

		ind_neighbor[3].x_index = i_minus;
		ind_neighbor[3].y_index = j_target;
		ind_neighbor[3].z_index = k_plus;

		y_dist  = (ib_node->im_point.x - y_interxn[0].x) * (ib_node->im_point.x - y_interxn[0].x);
		y_dist += (ib_node->im_point.y - y_interxn[0].y) * (ib_node->im_point.y - y_interxn[0].y);
		y_dist += (ib_node->im_point.z - y_interxn[0].z) * (ib_node->im_point.z - y_interxn[0].z);

		y_success = Immersed_find_plane_interpolation_coeff(ind_neighbor,
		            y_interxn, 2, which_quantity, grid, y_plane, ind_y_plane,
		            &y_u, &y_v, i_immersed, j_immersed, k_immersed, params);

	}

	if ( fabs(ib_node->n.vz) > 1e-8 )  {

		z_target = zq[k_target];
		Immersed_line_zplane_intersection_point(line, z_target, z_interxn);

		i_minus = Grid_get_x_index(z_interxn[0].x, grid, which_quantity);
		i_plus = i_minus + 1;
		j_minus = Grid_get_y_index(z_interxn[0].y, grid, which_quantity);
		j_plus = j_minus + 1;

		/*--------------------------------------------------------------------*/
		/*
		 DO NOT change the order of the nodes below, as the
		 "find_plane_interpolation_coeff" assumes the following order of the
		 nodes:

		 1---------------2
		 |               |
		 |               |
		 |               |
		 0---------------3

		 */
		/*--------------------------------------------------------------------*/
		ind_neighbor[0].x_index = i_minus;
		ind_neighbor[0].y_index = j_minus;
		ind_neighbor[0].z_index = k_target;

		ind_neighbor[1].x_index = i_minus;
		ind_neighbor[1].y_index = j_plus;
		ind_neighbor[1].z_index = k_target;

		ind_neighbor[2].x_index = i_plus;
		ind_neighbor[2].y_index = j_plus;
		ind_neighbor[2].z_index = k_target;

		ind_neighbor[3].x_index = i_plus;
		ind_neighbor[3].y_index = j_minus;
		ind_neighbor[3].z_index = k_target;

		z_dist  = (ib_node->im_point.x - z_interxn[0].x) * (ib_node->im_point.x - z_interxn[0].x);
		z_dist += (ib_node->im_point.y - z_interxn[0].y) * (ib_node->im_point.y - z_interxn[0].y);
		z_dist += (ib_node->im_point.z - z_interxn[0].z) * (ib_node->im_point.z - z_interxn[0].z);

		z_success = Immersed_find_plane_interpolation_coeff(ind_neighbor,
		          z_interxn, 3, which_quantity, grid, z_plane, ind_z_plane,
		          &z_u, &z_v, i_immersed, j_immersed, k_immersed, params);

	}

	if ( ! (x_success || y_success || z_success) ) {
		return 0;
	}

	min_dist = 10*max(max(x_dist, y_dist), z_dist);
	min_index = 10;

	if ( (x_success) && (x_dist < min_dist) ) {
		min_dist = x_dist;
		min_index = 1;
	}

	if ( (y_success) && (y_dist < min_dist) ) {
		min_dist = y_dist;
		min_index = 2;
	}

	if ( (z_success) && (z_dist < min_dist) ) {
		min_dist = z_dist;
		min_index = 3;
	}

	if (min_index == 1) {

		ib_node->fluid_index[0].x_index = ind_x_plane[0].x_index;
		ib_node->fluid_index[0].y_index = ind_x_plane[0].y_index;
		ib_node->fluid_index[0].z_index = ind_x_plane[0].z_index;

		ib_node->fluid_index[1].x_index = ind_x_plane[1].x_index;
		ib_node->fluid_index[1].y_index = ind_x_plane[1].y_index;
		ib_node->fluid_index[1].z_index = ind_x_plane[1].z_index;

		ib_node->fluid_index[2].x_index = ind_x_plane[2].x_index;
		ib_node->fluid_index[2].y_index = ind_x_plane[2].y_index;
		ib_node->fluid_index[2].z_index = ind_x_plane[2].z_index;


		ib_node->fluid_coef[0] = 1 - x_u - x_v;
		ib_node->fluid_coef[1] = x_u;
		ib_node->fluid_coef[2] = x_v;

		ib_node->intersection_point.x = x_interxn[0].x;
		ib_node->intersection_point.y = x_interxn[0].y;
		ib_node->intersection_point.z = x_interxn[0].z;
	}

	if (min_index == 2) {

		ib_node->fluid_index[0].x_index = ind_y_plane[0].x_index;
		ib_node->fluid_index[0].y_index = ind_y_plane[0].y_index;
		ib_node->fluid_index[0].z_index = ind_y_plane[0].z_index;

		ib_node->fluid_index[1].x_index = ind_y_plane[1].x_index;
		ib_node->fluid_index[1].y_index = ind_y_plane[1].y_index;
		ib_node->fluid_index[1].z_index = ind_y_plane[1].z_index;

		ib_node->fluid_index[2].x_index = ind_y_plane[2].x_index;
		ib_node->fluid_index[2].y_index = ind_y_plane[2].y_index;
		ib_node->fluid_index[2].z_index = ind_y_plane[2].z_index;

		ib_node->fluid_coef[0] = 1 - y_u - y_v;
		ib_node->fluid_coef[1] = y_u;
		ib_node->fluid_coef[2] = y_v;

		ib_node->intersection_point.x = y_interxn[0].x;
		ib_node->intersection_point.y = y_interxn[0].y;
		ib_node->intersection_point.z = y_interxn[0].z;
	}

	if (min_index == 3) {

		ib_node->fluid_index[0].x_index = ind_z_plane[0].x_index;
		ib_node->fluid_index[0].y_index = ind_z_plane[0].y_index;
		ib_node->fluid_index[0].z_index = ind_z_plane[0].z_index;

		ib_node->fluid_index[1].x_index = ind_z_plane[1].x_index;
		ib_node->fluid_index[1].y_index = ind_z_plane[1].y_index;
		ib_node->fluid_index[1].z_index = ind_z_plane[1].z_index;

		ib_node->fluid_index[2].x_index = ind_z_plane[2].x_index;
		ib_node->fluid_index[2].y_index = ind_z_plane[2].y_index;
		ib_node->fluid_index[2].z_index = ind_z_plane[2].z_index;

		ib_node->fluid_coef[0] = 1 - z_u - z_v;
		ib_node->fluid_coef[1] = z_u;
		ib_node->fluid_coef[2] = z_v;

		ib_node->intersection_point.x = z_interxn[0].x;
		ib_node->intersection_point.y = z_interxn[0].y;
		ib_node->intersection_point.z = z_interxn[0].z;
	}


	if (boundary_condition == DIRICHLET) {

		d1  = (ib_node->intersection_point.x - ib_node->im_point.x)*(ib_node->intersection_point.x - ib_node->im_point.x);
		d1 += (ib_node->intersection_point.y - ib_node->im_point.y)*(ib_node->intersection_point.y - ib_node->im_point.y);
		d1 += (ib_node->intersection_point.z - ib_node->im_point.z)*(ib_node->intersection_point.z - ib_node->im_point.z);
		d1 = sqrt(d1);

		d2  = (ib_node->boundary_point.x - ib_node->im_point.x)*(ib_node->boundary_point.x - ib_node->im_point.x);
		d2 += (ib_node->boundary_point.y - ib_node->im_point.y)*(ib_node->boundary_point.y - ib_node->im_point.y);
		d2 += (ib_node->boundary_point.z - ib_node->im_point.z)*(ib_node->boundary_point.z - ib_node->im_point.z);
		d2 = sqrt(d2);

		factor = d2/(d1+d2);
		ib_node->fluid_coef[0] *= factor;
		ib_node->fluid_coef[1] *= factor;
		ib_node->fluid_coef[2] *= factor;
	}
	return 1;

}




/******************************************************************************/
/*
 This function returns the ith immersed node
 */
/******************************************************************************/
ImmersedNode *Immersed_get_ib_node(Immersed *q_immersed, int i) {

	return &(q_immersed->ib_nodes[i]);
}




/******************************************************************************/
/*
 It returns 1 if the intersection point is within the triangle formed by three
 points of the plane. It returns 0 otherwise

 Check http://en.wikipedia.org/wiki/Line-plane_intersection
 */
/******************************************************************************/
int Immersed_line_plane_intersection(PointType *line, PointType *plane,
		double *t, double *u, double *v) {

	double a, b, c, d, e, f, g, h, i;
	double inv_a, inv_b, inv_c, inv_d, inv_e, inv_f, inv_g, inv_h, inv_i;
	double detA, inv_detA;
	double rhs_x, rhs_y, rhs_z;

	a = line[0].x - line[1].x;
	d = line[0].y - line[1].y;
	g = line[0].z - line[1].z;

	b = plane[1].x - plane[0].x;
	e = plane[1].y - plane[0].y;
	h = plane[1].z - plane[0].z;

	c = plane[2].x - plane[0].x;
	f = plane[2].y - plane[0].y;
	i = plane[2].z - plane[0].z;

	rhs_x = line[0].x - plane[0].x;
	rhs_y = line[0].y - plane[0].y;
	rhs_z = line[0].z - plane[0].z;

	detA = a*(e*i - h*f) - b*(d*i - f*g) + c*(d*h - e*g);
	if (fabs(detA) > 1e-32) {
		inv_detA = 1.0/detA;
	}
	else {
		*t = 0.0;
		*u = 0.0;
		*v = 0.0;
		return 0;
	}

	inv_a =  inv_detA*(e*i - f*h);
	inv_b = -inv_detA*(b*i - c*h);
	inv_c =  inv_detA*(b*f - c*e);

	inv_d = -inv_detA*(d*i - f*g);
	inv_e =  inv_detA*(a*i - c*g);
	inv_f = -inv_detA*(a*f - c*d);

	inv_g =  inv_detA*(d*h - e*g);
	inv_h = -inv_detA*(a*h - b*g);
	inv_i =  inv_detA*(a*e - b*d);

	*t = inv_a*rhs_x + inv_b*rhs_y + inv_c*rhs_z;
	*u = inv_d*rhs_x + inv_e*rhs_y + inv_f*rhs_z;
	*v = inv_g*rhs_x + inv_h*rhs_y + inv_i*rhs_z;

	if ( (*u > -1e-4) && (*v > -1e-4) && ( (*u + *v) < 1.+1e-4) ) {
		return 1;
	}
	else {
		return 0;
	}
}




/******************************************************************************/
/*
 This function calculates the intersection point of line on a plane

 It returns 1 if the intersection point is within the triangle formed by three
 points of the plane. It returns 0 otherwise

 Check http://en.wikipedia.org/wiki/Line-plane_intersection
 */
/******************************************************************************/
int Immersed_find_plane_interpolation_coeff(Indices *ind_neighbor,
		PointType *interxn, int idir, char which_quantity, MAC_grid *grid,
		PointType *plane, Indices *ind_plane, double *u, double *v,
		int i_immersed, int j_immersed, int k_immersed, Parameters *params) {

	int i, j, k;
	int ***status;
	int (*Immersed_point_location_on_plane)(PointType *, PointType *, double *, double *);
	double *xq, *yq, *zq;
	PointType *plane_neighbor;
	int success;
	double tol;
	char filename[50];
	FILE *fid;

	tol = 1e-6;
	success = 0;
	switch (which_quantity) {

		case 'u':
			status = grid->u_status;
			xq = grid->xu;
			yq = grid->yc;
			zq = grid->zc;
			break;

		case 'v':
			status = grid->v_status;
			xq = grid->xc;
			yq = grid->yv;
			zq = grid->zc;
			break;

		case 'w':
			status = grid->w_status;
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zw;
			break;

		case 'c':
			status = grid->c_status;
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zc;
			break;
	}  /* switch */

#ifdef IMMERSED_DIAGNOSTIC
	sprintf(filename,"interpolation_xyzcoeff%d.dat",params->rank);
	fid = fopen(filename,"a");
	fprintf(fid, "\ni_immersed = %d j_immersed = %d k_immersed = %d idir= %d\n", i_immersed, j_immersed, k_immersed, idir);
	fprintf(fid, "x_inter = %f y_inter = %f z_inter = %f \n", interxn[0].x, interxn[0].y, interxn[0].z);
	fprintf(fid, "xq = %f yq = %f zq = %f \n", xq[i_immersed], yq[j_immersed], zq[k_immersed]);

	fprintf(fid, "n0 = %d n1 = %d n2 = %d n3 = %d \n",
			status[ind_neighbor[0].z_index][ind_neighbor[0].y_index][ind_neighbor[0].x_index],
			status[ind_neighbor[1].z_index][ind_neighbor[1].y_index][ind_neighbor[1].x_index],
			status[ind_neighbor[2].z_index][ind_neighbor[2].y_index][ind_neighbor[2].x_index],
			status[ind_neighbor[3].z_index][ind_neighbor[3].y_index][ind_neighbor[3].x_index]);
	fprintf(fid,"i0 = %d j0 = %d k0 = %d \n",ind_neighbor[0].x_index, ind_neighbor[0].y_index, ind_neighbor[0].z_index);
	fprintf(fid,"i1 = %d j1 = %d k1 = %d \n",ind_neighbor[1].x_index, ind_neighbor[1].y_index, ind_neighbor[1].z_index);
	fprintf(fid,"i2 = %d j2 = %d k2 = %d \n",ind_neighbor[2].x_index, ind_neighbor[2].y_index, ind_neighbor[2].z_index);
	fprintf(fid,"i3 = %d j3 = %d k3 = %d \n",ind_neighbor[3].x_index, ind_neighbor[3].y_index, ind_neighbor[3].z_index);
	fprintf(fid,"x0 = %f y0 = %f z0 = %f \n",xq[ind_neighbor[0].x_index], yq[ind_neighbor[0].y_index], zq[ind_neighbor[0].z_index]);
	fprintf(fid,"x1 = %f y1 = %f z1 = %f \n",xq[ind_neighbor[1].x_index], yq[ind_neighbor[1].y_index], zq[ind_neighbor[1].z_index]);
	fprintf(fid,"x2 = %f y2 = %f z2 = %f \n",xq[ind_neighbor[2].x_index], yq[ind_neighbor[2].y_index], zq[ind_neighbor[2].z_index]);
	fprintf(fid,"x3 = %f y3 = %f z3 = %f \n",xq[ind_neighbor[3].x_index], yq[ind_neighbor[3].y_index], zq[ind_neighbor[3].z_index]);
#endif


	if ( idir == 1) {
		Immersed_point_location_on_plane = &Immersed_point_location_on_xplane;
	}
	else if ( idir == 2) {
		Immersed_point_location_on_plane = &Immersed_point_location_on_yplane;
	}
	else {
		Immersed_point_location_on_plane = &Immersed_point_location_on_zplane;
	}

	if ( (status[ind_neighbor[0].z_index][ind_neighbor[0].y_index][ind_neighbor[0].x_index] == FLUID) &&
	     (status[ind_neighbor[1].z_index][ind_neighbor[1].y_index][ind_neighbor[1].x_index] == FLUID) &&
	     (status[ind_neighbor[2].z_index][ind_neighbor[2].y_index][ind_neighbor[2].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[0].x_index];
		plane[0].y = yq[ind_neighbor[0].y_index];
		plane[0].z = zq[ind_neighbor[0].z_index];

		plane[1].x = xq[ind_neighbor[1].x_index];
		plane[1].y = yq[ind_neighbor[1].y_index];
		plane[1].z = zq[ind_neighbor[1].z_index];

		plane[2].x = xq[ind_neighbor[2].x_index];
		plane[2].y = yq[ind_neighbor[2].y_index];
		plane[2].z = zq[ind_neighbor[2].z_index];

		success = Immersed_point_location_on_plane(interxn, plane, u, v);
#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid, "x0 = %f y0= %f z0 = %f \n", plane[0].x, plane[0].y, plane[0].z);
		fprintf(fid, "x1 = %f y1= %f z1 = %f \n", plane[1].x, plane[1].y, plane[1].z);
		fprintf(fid, "x2 = %f y2= %f z2 = %f \n", plane[2].x, plane[2].y, plane[2].z);
		fprintf(fid, "success1 = %d\n",success);
#endif

		if (success) {
			ind_plane[0].x_index = ind_neighbor[0].x_index;
			ind_plane[0].y_index = ind_neighbor[0].y_index;
			ind_plane[0].z_index = ind_neighbor[0].z_index;

			ind_plane[1].x_index = ind_neighbor[1].x_index;
			ind_plane[1].y_index = ind_neighbor[1].y_index;
			ind_plane[1].z_index = ind_neighbor[1].z_index;

			ind_plane[2].x_index = ind_neighbor[2].x_index;
			ind_plane[2].y_index = ind_neighbor[2].y_index;
			ind_plane[2].z_index = ind_neighbor[2].z_index;
#ifdef IMMERSED_DIAGNOSTIC
			fclose(fid);
#endif
			return success;
		}

	}

	if ( (status[ind_neighbor[0].z_index][ind_neighbor[0].y_index][ind_neighbor[0].x_index] == FLUID) &&
	     (status[ind_neighbor[2].z_index][ind_neighbor[2].y_index][ind_neighbor[2].x_index] == FLUID) &&
	     (status[ind_neighbor[3].z_index][ind_neighbor[3].y_index][ind_neighbor[3].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[0].x_index];
		plane[0].y = yq[ind_neighbor[0].y_index];
		plane[0].z = zq[ind_neighbor[0].z_index];

		plane[1].x = xq[ind_neighbor[2].x_index];
		plane[1].y = yq[ind_neighbor[2].y_index];
		plane[1].z = zq[ind_neighbor[2].z_index];

		plane[2].x = xq[ind_neighbor[3].x_index];
		plane[2].y = yq[ind_neighbor[3].y_index];
		plane[2].z = zq[ind_neighbor[3].z_index];

		success = Immersed_point_location_on_plane(interxn, plane, u, v);
#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid, "x0 = %f y0= %f z0 = %f \n", plane[0].x, plane[0].y, plane[0].z);
		fprintf(fid, "x2 = %f y2= %f z2 = %f \n", plane[1].x, plane[1].y, plane[1].z);
		fprintf(fid, "x3 = %f y3= %f z3 = %f \n", plane[2].x, plane[2].y, plane[2].z);
		fprintf(fid, "success2 = %d\n",success);
#endif

		if (success) {
			ind_plane[0].x_index = ind_neighbor[0].x_index;
			ind_plane[0].y_index = ind_neighbor[0].y_index;
			ind_plane[0].z_index = ind_neighbor[0].z_index;

			ind_plane[1].x_index = ind_neighbor[2].x_index;
			ind_plane[1].y_index = ind_neighbor[2].y_index;
			ind_plane[1].z_index = ind_neighbor[2].z_index;

			ind_plane[2].x_index = ind_neighbor[3].x_index;
			ind_plane[2].y_index = ind_neighbor[3].y_index;
			ind_plane[2].z_index = ind_neighbor[3].z_index;
#ifdef IMMERSED_DIAGNOSTIC
			fclose(fid);
#endif
			return success;
		}

	}

	if ( (status[ind_neighbor[1].z_index][ind_neighbor[1].y_index][ind_neighbor[1].x_index] == FLUID) &&
	     (status[ind_neighbor[2].z_index][ind_neighbor[2].y_index][ind_neighbor[2].x_index] == FLUID) &&
	     (status[ind_neighbor[3].z_index][ind_neighbor[3].y_index][ind_neighbor[3].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[1].x_index];
		plane[0].y = yq[ind_neighbor[1].y_index];
		plane[0].z = zq[ind_neighbor[1].z_index];

		plane[1].x = xq[ind_neighbor[2].x_index];
		plane[1].y = yq[ind_neighbor[2].y_index];
		plane[1].z = zq[ind_neighbor[2].z_index];

		plane[2].x = xq[ind_neighbor[3].x_index];
		plane[2].y = yq[ind_neighbor[3].y_index];
		plane[2].z = zq[ind_neighbor[3].z_index];

		success = Immersed_point_location_on_plane(interxn, plane, u, v);
#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid, "x1 = %f y1= %f z1 = %f \n", plane[0].x, plane[0].y, plane[0].z);
		fprintf(fid, "x2 = %f y2= %f z2 = %f \n", plane[1].x, plane[1].y, plane[1].z);
		fprintf(fid, "x3 = %f y3= %f z3 = %f \n", plane[2].x, plane[2].y, plane[2].z);
		fprintf(fid, "success2 = %d\n",success);
#endif

		if (success) {
			ind_plane[0].x_index = ind_neighbor[1].x_index;
			ind_plane[0].y_index = ind_neighbor[1].y_index;
			ind_plane[0].z_index = ind_neighbor[1].z_index;

			ind_plane[1].x_index = ind_neighbor[2].x_index;
			ind_plane[1].y_index = ind_neighbor[2].y_index;
			ind_plane[1].z_index = ind_neighbor[2].z_index;

			ind_plane[2].x_index = ind_neighbor[3].x_index;
			ind_plane[2].y_index = ind_neighbor[3].y_index;
			ind_plane[2].z_index = ind_neighbor[3].z_index;
#ifdef IMMERSED_DIAGNOSTIC
			fclose(fid);
#endif
			return success;
		}

	}

	/*------------------------------------------------------------------------*/
	/*
	 If the interpolation is done along a x-, y- or z-direction, then we have
	 just two points in the interpolation stencil

	 So we set the third interpolation coefficient to zero and set the
	 coordinate to "immersed point" itself
	 */
	/*------------------------------------------------------------------------*/
	*v = 0.0;
	ind_plane[2].x_index = i_immersed;
	ind_plane[2].y_index = j_immersed;
	ind_plane[2].z_index = k_immersed;

	if ( (status[ind_neighbor[0].z_index][ind_neighbor[0].y_index][ind_neighbor[0].x_index] == FLUID) &&
	     (status[ind_neighbor[1].z_index][ind_neighbor[1].y_index][ind_neighbor[1].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[0].x_index];
		plane[0].y = yq[ind_neighbor[0].y_index];
		plane[0].z = zq[ind_neighbor[0].z_index];

		plane[1].x = xq[ind_neighbor[1].x_index];
		plane[1].y = yq[ind_neighbor[1].y_index];
		plane[1].z = zq[ind_neighbor[1].z_index];

		if (idir == 1) {
			if ( fabs(interxn[0].z - plane[0].z) < tol)  {
				//yline_intersection
				success = Immersed_point_location_on_yline(interxn, plane, u);
			}
		}
		else if (idir == 2) {
			if ( fabs(interxn[0].z - plane[0].z) < tol)  {
				//xline_intersection
				success = Immersed_point_location_on_xline(interxn, plane, u);
			}
		}
		else if (idir == 3) {
			if ( fabs(interxn[0].x - plane[0].x) < tol)  {
				//yline_intersection
				success = Immersed_point_location_on_yline(interxn, plane, u);
			}
		}
		if (success) {
			ind_plane[0].x_index = ind_neighbor[0].x_index;
			ind_plane[0].y_index = ind_neighbor[0].y_index;
			ind_plane[0].z_index = ind_neighbor[0].z_index;

			ind_plane[1].x_index = ind_neighbor[1].x_index;
			ind_plane[1].y_index = ind_neighbor[1].y_index;
			ind_plane[1].z_index = ind_neighbor[1].z_index;

#ifdef IMMERSED_DIAGNOSTIC
			fclose(fid);
#endif
			return success;
		}



	}

	if ( (status[ind_neighbor[0].z_index][ind_neighbor[0].y_index][ind_neighbor[0].x_index] == FLUID) &&
	     (status[ind_neighbor[2].z_index][ind_neighbor[2].y_index][ind_neighbor[2].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[0].x_index];
		plane[0].y = yq[ind_neighbor[0].y_index];
		plane[0].z = zq[ind_neighbor[0].z_index];

		plane[1].x = xq[ind_neighbor[2].x_index];
		plane[1].y = yq[ind_neighbor[2].y_index];
		plane[1].z = zq[ind_neighbor[2].z_index];

		if (idir == 1) {
			// cross-line
		}
		else if (idir == 2) {
			// cross-line
		}
		else if (idir == 3) {
			// cross-line
		}
	}

	if ( (status[ind_neighbor[0].z_index][ind_neighbor[0].y_index][ind_neighbor[0].x_index] == FLUID) &&
	     (status[ind_neighbor[3].z_index][ind_neighbor[3].y_index][ind_neighbor[3].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[0].x_index];
		plane[0].y = yq[ind_neighbor[0].y_index];
		plane[0].z = zq[ind_neighbor[0].z_index];

		plane[1].x = xq[ind_neighbor[3].x_index];
		plane[1].y = yq[ind_neighbor[3].y_index];
		plane[1].z = zq[ind_neighbor[3].z_index];

		if (idir == 1) {
			if ( fabs(interxn[0].y - plane[0].y) < tol)  {
				// zline_intersection
				success = Immersed_point_location_on_zline(interxn, plane, u);
			}
		}
		else if (idir == 2) {
			if ( fabs(interxn[0].x - plane[0].x) < tol)  {
				// zline_intersection
				success = Immersed_point_location_on_zline(interxn, plane, u);
			}
		}
		else if (idir == 3) {
			if ( fabs(interxn[0].y - plane[0].y) < tol)  {
				// xline_intersection
				success = Immersed_point_location_on_xline(interxn, plane, u);
			}
		}

		if (success) {
			ind_plane[0].x_index = ind_neighbor[0].x_index;
			ind_plane[0].y_index = ind_neighbor[0].y_index;
			ind_plane[0].z_index = ind_neighbor[0].z_index;

			ind_plane[1].x_index = ind_neighbor[3].x_index;
			ind_plane[1].y_index = ind_neighbor[3].y_index;
			ind_plane[1].z_index = ind_neighbor[3].z_index;

#ifdef IMMERSED_DIAGNOSTIC
			fclose(fid);
#endif
			return success;
		}
	}

	if ( (status[ind_neighbor[1].z_index][ind_neighbor[1].y_index][ind_neighbor[1].x_index] == FLUID) &&
	     (status[ind_neighbor[2].z_index][ind_neighbor[2].y_index][ind_neighbor[2].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[1].x_index];
		plane[0].y = yq[ind_neighbor[1].y_index];
		plane[0].z = zq[ind_neighbor[1].z_index];

		plane[1].x = xq[ind_neighbor[2].x_index];
		plane[1].y = yq[ind_neighbor[2].y_index];
		plane[1].z = zq[ind_neighbor[2].z_index];

		if (idir == 1) {
			if ( fabs(interxn[0].y - plane[0].y) < tol)  {
				// zline_intersection
				success = Immersed_point_location_on_zline(interxn, plane, u);
			}
		}
		else if (idir == 2) {
			if ( fabs(interxn[0].x - plane[0].x) < tol)  {
				// zline_intersection
				success = Immersed_point_location_on_zline(interxn, plane, u);
			}
		}
		else if (idir == 3) {
			if ( fabs(interxn[0].y - plane[0].y) < tol)  {
				// xline_intersection
				success = Immersed_point_location_on_xline(interxn, plane, u);
			}
		}
		if (success) {
			ind_plane[0].x_index = ind_neighbor[1].x_index;
			ind_plane[0].y_index = ind_neighbor[1].y_index;
			ind_plane[0].z_index = ind_neighbor[1].z_index;

			ind_plane[1].x_index = ind_neighbor[2].x_index;
			ind_plane[1].y_index = ind_neighbor[2].y_index;
			ind_plane[1].z_index = ind_neighbor[2].z_index;

#ifdef IMMERSED_DIAGNOSTIC
			fclose(fid);
#endif
			return success;
		}

	}

	if ( (status[ind_neighbor[1].z_index][ind_neighbor[1].y_index][ind_neighbor[1].x_index] == FLUID) &&
	     (status[ind_neighbor[3].z_index][ind_neighbor[3].y_index][ind_neighbor[3].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[1].x_index];
		plane[0].y = yq[ind_neighbor[1].y_index];
		plane[0].z = zq[ind_neighbor[1].z_index];

		plane[1].x = xq[ind_neighbor[3].x_index];
		plane[1].y = yq[ind_neighbor[3].y_index];
		plane[1].z = zq[ind_neighbor[3].z_index];

		if (idir == 1) {
			// cross-line
		}
		else if (idir == 2) {
			// cross-line
		}
		else if (idir == 3) {
			// cross-line
		}
	}

	if ( (status[ind_neighbor[2].z_index][ind_neighbor[2].y_index][ind_neighbor[2].x_index] == FLUID) &&
	     (status[ind_neighbor[3].z_index][ind_neighbor[3].y_index][ind_neighbor[3].x_index] == FLUID) ) {

		plane[0].x = xq[ind_neighbor[2].x_index];
		plane[0].y = yq[ind_neighbor[2].y_index];
		plane[0].z = zq[ind_neighbor[2].z_index];

		plane[1].x = xq[ind_neighbor[3].x_index];
		plane[1].y = yq[ind_neighbor[3].y_index];
		plane[1].z = zq[ind_neighbor[3].z_index];

		if (idir == 1) {
			if ( fabs(interxn[0].z - plane[0].z) < tol)  {
				// yline_intersection
				success = Immersed_point_location_on_yline(interxn, plane, u);
			}
		}
		else if (idir == 2) {
			if ( fabs(interxn[0].z - plane[0].z) < tol)  {
				// xline_intersection
				success = Immersed_point_location_on_xline(interxn, plane, u);
			}
		}
		else if (idir == 3) {
			if ( fabs(interxn[0].x - plane[0].x) < tol)  {
				// yline_intersection
				success = Immersed_point_location_on_yline(interxn, plane, u);
			}
		}
		if (success) {
			ind_plane[0].x_index = ind_neighbor[2].x_index;
			ind_plane[0].y_index = ind_neighbor[2].y_index;
			ind_plane[0].z_index = ind_neighbor[2].z_index;

			ind_plane[1].x_index = ind_neighbor[3].x_index;
			ind_plane[1].y_index = ind_neighbor[3].y_index;
			ind_plane[1].z_index = ind_neighbor[3].z_index;

#ifdef IMMERSED_DIAGNOSTIC
			fclose(fid);
#endif
			return success;
		}
	}


#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid,"failure after three tries\n");
		fclose(fid);
#endif

		plane[0].x = xq[ind_neighbor[0].x_index];
		plane[0].y = yq[ind_neighbor[0].y_index];
		plane[0].z = zq[ind_neighbor[0].z_index];

		plane[1].x = xq[ind_neighbor[1].x_index];
		plane[1].y = yq[ind_neighbor[1].y_index];
		plane[1].z = zq[ind_neighbor[1].z_index];

		plane[2].x = xq[ind_neighbor[2].x_index];
		plane[2].y = yq[ind_neighbor[2].y_index];
		plane[2].z = zq[ind_neighbor[2].z_index];
		*u = 0.0;
		*v = 0.0;

	return 0;
}




/******************************************************************************/
/*
 This function calculates the intersection point of line on a x-plane
 */
/******************************************************************************/
void Immersed_line_xplane_intersection_point(PointType *line, double x,
		PointType *intersection_point) {

	double t;

	t = (x - line[0].x)/(line[1].x - line[0].x);

	intersection_point[0].x = x;
	intersection_point[0].y = line[0].y + t*(line[1].y - line[0].y);
	intersection_point[0].z = line[0].z + t*(line[1].z - line[0].z);

	return;
}




/******************************************************************************/
/*
 This function calculates the intersection point of line on a y-plane
 */
/******************************************************************************/
void Immersed_line_yplane_intersection_point(PointType *line, double y,
		PointType *intersection_point) {

	double t;

	t = (y - line[0].y)/(line[1].y - line[0].y);

	intersection_point[0].x = line[0].x + t*(line[1].x - line[0].x);
	intersection_point[0].y = y;
	intersection_point[0].z = line[0].z + t*(line[1].z - line[0].z);

	return;

}




/******************************************************************************/
/*
 This function calculates the intersection point of line on a y-plane
 */
/******************************************************************************/
void Immersed_line_zplane_intersection_point(PointType *line, double z,
		PointType *intersection_point) {

	double t;

	t = (z - line[0].z)/(line[1].z - line[0].z);

	intersection_point[0].x = line[0].x + t*(line[1].x - line[0].x);
	intersection_point[0].y = line[0].y + t*(line[1].y - line[0].y);
	intersection_point[0].z = z;

	return;
}




/******************************************************************************/
/*
 This function calculates the intersection point of line on a x-plane

 It returns 1 if the input point is within the triangle formed by three points
 of the plane. It returns 0 otherwise

     Eqn
     (1-u-v) * point_0 + u * point_1 + v * point_2 = point_p
 */
/******************************************************************************/
int Immersed_point_location_on_xplane(PointType *point, PointType *plane,
		double *u, double *v) {

	double a, b, c, d;
	double rhs_1, rhs_2;
	double inv_denom;

	a = plane[1].y - plane[0].y;
	c = plane[1].z - plane[0].z;

	b = plane[2].y - plane[0].y;
	d = plane[2].z - plane[0].z;

	rhs_1 = point[0].y - plane[0].y;
	rhs_2 = point[0].z - plane[0].z;

	inv_denom = 1.0/(a*d - b*c);

	*u = (rhs_1*d - rhs_2*b) * inv_denom;
	*v = (rhs_2*a - rhs_1*c) * inv_denom;

	if ( (*u > -1e-4) && (*v > -1e-4) && ( (*u + *v) < 1.+1e-4) ) {
		return 1;
	}
	else {
		return 0;
	}
}




/******************************************************************************/
/*
 This function calculates the intersection point of line on a y-plane

 It returns 1 if the input point is within the triangle formed by three points
 of the plane. It returns 0 otherwise
     Eqn
     (1-u-v) * point_0 + u * point_1 + v * point_2 = point_p
 */
/******************************************************************************/
int Immersed_point_location_on_yplane(PointType *point, PointType *plane,
		double *u, double *v) {

	double a, b, c, d;
	double rhs_1, rhs_2;
	double inv_denom;

	a = plane[1].x - plane[0].x;
	c = plane[1].z - plane[0].z;

	b = plane[2].x - plane[0].x;
	d = plane[2].z - plane[0].z;

	rhs_1 = point[0].x - plane[0].x;
	rhs_2 = point[0].z - plane[0].z;

	inv_denom = 1.0/(a*d - b*c);

	*u = (rhs_1*d - rhs_2*b) * inv_denom;
	*v = (rhs_2*a - rhs_1*c) * inv_denom;

	if ( (*u > -1e-4) && (*v > -1e-4) && ( (*u + *v) < 1.+1e-4) ) {
		return 1;
	}
	else {
		return 0;
	}
}




/******************************************************************************/
/*
 This function calculates the intersection point of line on a z-plane

 It returns 1 if the input point is within the triangle formed by three points
 of the plane. It returns 0 otherwise
     Eqn
     (1-u-v) * point_0 + u * point_1 + v * point_2 = point_p
 */
/******************************************************************************/
int Immersed_point_location_on_zplane(PointType *point, PointType *plane,
		double *u, double *v) {

	double a, b, c, d;
	double rhs_1, rhs_2;
	double inv_denom;

	a = plane[1].x - plane[0].x;
	c = plane[1].y - plane[0].y;

	b = plane[2].x - plane[0].x;
	d = plane[2].y - plane[0].y;

	rhs_1 = point[0].x - plane[0].x;
	rhs_2 = point[0].y - plane[0].y;

	inv_denom = 1.0/(a*d - b*c);

	*u = (rhs_1*d - rhs_2*b) * inv_denom;
	*v = (rhs_2*a - rhs_1*c) * inv_denom;

	if ( (*u > -1e-4) && (*v > -1e-4) && ( (*u + *v) < 1.+1e-4) ) {
		return 1;
	}
	else {
		return 0;
	}
}




/******************************************************************************/
/*
 This function calculates the intersection point of point on a x-line

 It returns 1 if the input point is between the two points of the line. It
 returns 0 otherwise
     Eqn
     (1-u) * point_0 + u * point_1  = point_p
     or u * (point_1 - point_0) = point_p - point_0 */
/******************************************************************************/
int Immersed_point_location_on_xline(PointType *point, PointType *line, double *u) {

	double a, b;

	a = point[0].x - line[0].x;
	b =  line[1].x - line[0].x;

	*u = a/b;

	if ( (*u > -1e-4) && ( *u  < 1.+1e-4) ) {
		return 1;
	}
	else {
		return 0;
	}
}




/******************************************************************************/
/*
 This function calculates the intersection point of point on a y-line

 It returns 1 if the input point is between the two points of the line. It
 returns 0 otherwise
     Eqn
     (1-u) * point_0 + u * point_1  = point_p
     or u * (point_1 - point_0) = point_p - point_0 */
/******************************************************************************/
int Immersed_point_location_on_yline(PointType *point, PointType *line, double *u) {

	double a, b;

	a = point[0].y - line[0].y;
	b =  line[1].y - line[0].y;

	*u = a/b;

	if ( (*u > -1e-4) && ( *u  < 1.+1e-4) ) {
		return 1;
	}
	else {
		return 0;
	}
}




/******************************************************************************/
/*
 This function calculates the intersection point of point on a z-line

 It returns 1 if the input point is between the two points of the line. It
 returns 0 otherwise
     Eqn
     (1-u) * point_0 + u * point_1  = point_p
     or u * (point_1 - point_0) = point_p - point_0
 */
/******************************************************************************/
int Immersed_point_location_on_zline(PointType *point, PointType *line, double *u) {

	double a, b;

	a = point[0].z - line[0].z;
	b =  line[1].z - line[0].z;

	*u = a/b;

	if ( (*u > -1e-4) && ( *u  < 1.+1e-4) ) {
		return 1;
	}
	else {
		return 0;
	}
}




/******************************************************************************/
/*
 This function returns "approximate" normal distance from the surface i.e
 sdf = y - y_interface

 Accurate normal distance (sdf) is needed for only cells near the surface. So
 for five cells above and below the interface surface, this function calculates
 better approximation of the normal distance.
 */
/******************************************************************************/
void Immersed_find_surface_normal_distance(MAC_grid *grid, Parameters *params,
		char which_quantity) {


	double ***sdf;
	double *xq, *yq;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int istart_fine;
	int i, j, k;
	double dist, y_limit, dot_product;
	int finer_grid_step;
	int NX, NY, NZ;
	int ifine, ifine_end;
	double Ly;
	double *interface_xfine, *interface_yfine;
	double tx, ty;
	double nx, ny;
	double adist;


	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	finer_grid_step = grid->fine_step;
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Ly = params->Ly;
	interface_xfine = grid->finer_1d_interface_x;
	interface_yfine = grid->finer_1d_interface_y;
	ifine_end = grid->fine_step*(NX-1) + 1;

	switch (which_quantity) {

		case 'u':
			xq = grid->xu;
			yq = grid->yc;
			istart_fine = 0;
			sdf = grid->u_sdf;
			break;

		case 'v':
			xq = grid->xc;
			yq = grid->yv;
			istart_fine = finer_grid_step/2;
			sdf = grid->v_sdf;
			Ie = min(Ie, NX-1);
			break;

		case 'w':
			xq = grid->xc;
			yq = grid->yc;
			istart_fine = finer_grid_step/2;
			sdf = grid->w_sdf;
			Ie = min(Ie, NX-1);
			break;

		case 'c':
			xq = grid->xc;
			yq = grid->yc;
			istart_fine = finer_grid_step/2;
			sdf = grid->c_sdf;
			Ie = min(Ie, NX-1);
			break;
	} /* switch */

	y_limit = 0.0;
	dist = 3.*params->Ly;


	for (j=0; j<NY; j++) {
		if (y_limit < grid->dy_v[j]) y_limit = grid->dy_v[j];
	}
	y_limit = 5.*y_limit;

	k = Ks;
	for (j=Js; j<Je; j++) {
		for (i=Is; i<Ie; i++) {

			sdf[k][j][i] = yq[j] - grid->finer_1d_interface_y[istart_fine + i*finer_grid_step] ;
			if (sdf[k][j][i] < y_limit ) sdf[k][j][i] = dist;
		}
	}

	for (j=Js; j<Je; j++) {
		for (i=Is; i<Ie; i++) {
			if ( (j==0) && (which_quantity=='v')) {
				sdf[k][j][i] = 0.0;
			}
			else if (sdf[k][j][i] > 2.*Ly ) {
				for (ifine=1; ifine<ifine_end-1; ifine++) {

					dist = (interface_xfine[ifine] - xq[i])*(interface_xfine[ifine] - xq[i]) +
					     (interface_yfine[ifine] - yq[j])*(interface_yfine[ifine] - yq[j]);
					dist = pow(dist, 0.5);

					if (dist < fabs(sdf[k][j][i]) ) {
						sdf[k][j][i] = dist;
						/*----------------------------------------------------*/
						/*
						 Calculate the dot product of the normal vector to the
						 surface(interface) and the vector joining the immersed
						 node to the point on the boundary

						 Note the normal vector to the surface is
						 (-t_y i + t_x j) where (t_x i + t_y j) is the
						 tangential vector to the interface
						 */
						/*----------------------------------------------------*/
						dot_product =-(interface_yfine[ifine+1]-interface_yfine[ifine-1]) * (xq[i]-interface_xfine[ifine]) +
						         (interface_xfine[ifine+1]-interface_xfine[ifine-1]) * (yq[j]-interface_yfine[ifine]);
						if (dot_product < 0) sdf[k][j][i] = -dist;
					}
/*
				tx = interface_xfine[ifine+1]-interface_xfine[ifine-1];
				ty = interface_yfine[ifine+1]-interface_yfine[ifine-1];
				adist = tx*tx + ty*ty;
				adist = pow(adist, 0.5);
				tx = tx/adist;
				ty = ty/adist;

				nx = xq[i] - interface_xfine[ifine];
				ny = yq[j] - interface_yfine[ifine];
				adist = nx*nx + ny*ny;
				adist = pow(adist, 0.5);
				nx = nx/adist;
				ny = ny/adist;
*/
				}

			}
		}
	}


	for (k=Ks+1; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				sdf[k][j][i] = sdf[Ks][j][i];
			}
		}
	}

	return;
}




/******************************************************************************/
/*
 This function finds the normal vector from the surface to the immersed point
 and the boundary point (i.e the point of intersection of normal vector on the
 surface
 */
/******************************************************************************/
void Immersed_find_boundary_point_and_normal_vector(MAC_grid *grid,
		Parameters *params, char which_quantity) {

	double ***sdf;
	double *xq, *yq, *zq;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int istart_fine;
	int i, j, k;
	double dist, normal;
	int finer_grid_step;
	int NX, NY, NZ;
	int ifine, ifine_end;
	double Ly;
	double *interface_xfine, *interface_yfine;
	int ib;
	int N;
	ImmersedNode *ib_node;
	ImmersedNode *ib_node1;
	Immersed *q_immersed;
	FILE *fid;
	FILE *fid1;
	char bin_filename[50];
	double x_im, y_im, z_im;
	double min_dist;
	int i_im, j_im, k_im;
	double nx, ny, nz;
	double tx, ty;
	int iloc;
	double x1, x2, x3, y1, y2, y3;
	double ax, bx, cx, ay, by, cy;
	double xloc, yloc, s;


	dist = 2.*params->Ly;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	finer_grid_step = grid->fine_step;
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Ly = params->Ly;
	interface_xfine = grid->finer_1d_interface_x;
	interface_yfine = grid->finer_1d_interface_y;
	ifine_end = grid->fine_step*(NX-1) + 1;

	switch (which_quantity) {

		case 'u':
			xq = grid->xu;
			yq = grid->yc;
			zq = grid->zc;
			istart_fine = 0;
			sdf = grid->u_sdf;
			q_immersed  = grid->u_immersed;
			break;

		case 'v':
			xq = grid->xc;
			yq = grid->yv;
			zq = grid->zc;
			istart_fine = finer_grid_step/2;
			sdf = grid->v_sdf;
			Ie = min(Ie, NX-1);
			q_immersed  = grid->v_immersed;
			break;

		case 'w':
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zw;
			istart_fine = finer_grid_step/2;
			sdf = grid->w_sdf;
			Ie = min(Ie, NX-1);
			q_immersed  = grid->w_immersed;
			break;

		case 'c':
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zc;
			istart_fine = finer_grid_step/2;
			sdf = grid->c_sdf;
			Ie = min(Ie, NX-1);
			q_immersed  = grid->c_immersed;
			break;
	} // switch


	// Total number of the immersed nodes
	N = q_immersed->N;
	if (N == 0)
		return;

#ifdef IMMERSED_DIAGNOSTIC
	sprintf(bin_filename,"boundary%c_%d.dat",which_quantity, params->rank);
	fid = fopen(bin_filename,"a");
	sprintf(bin_filename,"nboundary%c_%d.dat",which_quantity, params->rank);
	fid1 = fopen(bin_filename,"a");
#endif

	for (ib=0; ib<N; ib++) {

		ib_node = Immersed_get_ib_node(q_immersed, ib);

		// coordinates of the immersed node
		x_im = ib_node->im_point.x;
		y_im = ib_node->im_point.y;
		z_im = ib_node->im_point.z;

		// (i,j,k) of the current ib node
		i_im = ib_node->im_index.x_index;
		j_im = ib_node->im_index.y_index;
		k_im = ib_node->im_index.z_index;

		iloc = 1;
		min_dist = params->Lx + params->Ly + params->Lz;

		for (ifine=1; ifine<ifine_end-1; ifine++) {

			dist = (interface_xfine[ifine] - x_im)*(interface_xfine[ifine] - x_im) +
			       (interface_yfine[ifine] - y_im)*(interface_yfine[ifine] - y_im);
			dist = pow(dist, 0.5);

			if ( fabs(dist-fabs(sdf[k_im][j_im][i_im]) ) < min_dist) {
				iloc = ifine;
				min_dist = fabs(dist-fabs(sdf[k_im][j_im][i_im]));
			}
		}

		/*--------------------------------------------------------------------*/
		/*
		 Boundary is represented by

		     xs = ax*s*s + bx*s + cx and ys = ay*s*s + by*s + cy

		 where s lies between 0 and 1.

		 (x1, y1), (x2, y2) and (x3, y3) are the actual points on the boundary.
		 They are used to find the coefficients in the above quadratic
		 representation of the boundary
		 */
		/*------------------------------------------------------------------------*/

		x1 = interface_xfine[iloc-1];
		x2 = interface_xfine[iloc];
		x3 = interface_xfine[iloc+1];
		y1 = interface_yfine[iloc-1];
		y2 = interface_yfine[iloc];
		y3 = interface_yfine[iloc+1];

		cx = x1;
		bx = 4.*x2 - x3 - 3.*x1;
		ax = x3 - x1 - bx;

		cy = y1;
		by = 4.*y2 - y3 - 3.*y1;
		ay = y3 - y1 - by;

		s =  Immersed_scoordinate_of_boundary(x_im, y_im, ax, bx, cx, ay, by, cy);
		xloc = ax*s*s + bx*s + cx;
		yloc = ay*s*s + by*s + cy;
		if ((xloc<interface_xfine[iloc-1]) || (xloc>interface_xfine[iloc+1]) ) {
			printf("Not within limits s = %f x1= %f xloc = %f x3 = %f ax = %f bx = %f cx = %f\n", s, x1, xloc, x3, ax, bx, cx);
		}

		// Set the boundary point and the normal vector from the boundary point
		// to the immersed node
		ib_node->boundary_point.x = xloc;
		ib_node->boundary_point.y = yloc;
		ib_node->boundary_point.z = zq[k_im];
		nx = x_im - xloc;
		ny = y_im - yloc;
		nz = 0.0;

		dist = nx*nx + ny*ny;
		dist = pow(dist, 0.5);
		ib_node->n.vx = nx/dist;
		ib_node->n.vy = ny/dist;
		ib_node->n.vz = 0.0;

		tx =  2*ax*s + bx;
		ty =  2*ay*s + by;
		dist = tx*tx + ty*ty;
		dist = pow(dist, 0.5);
		tx = tx/dist;
		ty = ty/dist;

#ifdef IMMERSED_DIAGNOSTIC
		fprintf(fid1, "ib=%d index = ( %d %d %d) tx = %e ty = %e nx = %e ny = %e dot_prod = %e \n",
		        ib, ib_node->im_index.x_index, ib_node->im_index.y_index,ib_node->im_index.z_index,
		        tx, ty, ib_node->n.vx, ib_node->n.vy, tx*nx + ty*ny);

		fprintf(fid, "ib=%d index = (%d %d %d) immersed = ( %f %f %f) boundary point = ( %f %f %f) s = %f\n",
		        ib, ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index,
		        ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z, xloc, yloc, zq[k_im], s);
#endif
		if (  fabs(tx*nx + ty*ny) > 0.0001) {
			printf("large dot product x1 = %f x2 = %f x3 = %f x_im = %f ib = %d \n",x1,x2,x3,x_im,ib);
		}

	} // for ib

#ifdef IMMERSED_DIAGNOSTIC
	fclose(fid);
	fclose(fid1);
#endif

	return;
}




/******************************************************************************/
/*
 Use Newton-Raphson's method to calculate the s-coordinate which minimizes the
 distance between the surface and the immersed point.

 Boundary is represented by

     xs = ax*s*s + bx*s + cx and ys = ay*s*s + by*s + cy

 where s lies between 0 and 1. x & y are the coordinates of the immersed point
 Minimize D^2 = (xs-x)^2 + (ys-ys)^2
 */
/******************************************************************************/
double Immersed_scoordinate_of_boundary(double x, double y, double ax,
		double bx, double cx, double ay, double by, double cy) {

	double s, snew, diff;
	double xs, ys, fs, dfds;
	double dxsds, dysds;

	s = 0.;
	do {
		xs = ax*s*s + bx*s + cx;
		ys = ay*s*s + by*s + cy;
		dxsds = 2*ax*s + bx;
		dysds = 2*ay*s + by;
		fs = (xs-x)*dxsds + (ys-y)*dysds;
		dfds = dxsds*dxsds + dysds*dysds + (xs-x)*2*ax + (ys-y)*2*ay;
		// Newton-Raphson xnew = xold - f/dfdx|at xold
		snew = s - fs/dfds;
		diff = snew - s;
		s = snew;
	} while (fabs(diff) < 1e-8);

	return s;
}




/******************************************************************************/
/*
 This function interpolates velocity at neighbhoring fluid nodes to find the
 velocity at the immersed node
 */
/******************************************************************************/
void Immersed_velocity_interpolation(Velocity *vel, MAC_grid *grid,
		Parameters *params) {

	int g, N;
	Immersed *q_immersed;
	ImmersedNode *ib_node;
	double ***data, ***data_old;
//	double ***p_data;
	int success;
	int i, j, k;
	int i0, j0, k0, i1, j1, k1, i2, j2, k2;
	double coeff0, coeff1, coeff2;
	double eta, vel_corr;
	double press_corr;
	double dp_dxyz, dp_dxyz0, dp_dxyz1, dp_dxyz2;
	double *idxyz;

	switch (vel -> component) {

		case 'u':
			q_immersed = grid->u_immersed;
			idxyz = grid->idx_c;
			break;

		case 'v':
			q_immersed = grid->v_immersed;
			idxyz = grid->idy_c;
			break;

		case 'w':
			q_immersed = grid->w_immersed;
			idxyz = grid->idz_c;
			break;

	}  // switch

	data = vel->data;
	data_old = vel->data_old;
//	p_data = p->p_data;

	// Total number of the immersed nodes
	N = q_immersed->N;



	// Go through all the ib nodes and find the neighboring fluid nodes. Store
	// the indices and the coordinates
	for (g=0; g<N; g++) {

		// Get the current immersed node
		ib_node = Immersed_get_ib_node(q_immersed, g);

		i = ib_node->im_index.x_index;
		j = ib_node->im_index.y_index;
		k = ib_node->im_index.z_index;

		i0 = ib_node->fluid_index[0].x_index;
		j0 = ib_node->fluid_index[0].y_index;
		k0 = ib_node->fluid_index[0].z_index;

		i1 = ib_node->fluid_index[1].x_index;
		j1 = ib_node->fluid_index[1].y_index;
		k1 = ib_node->fluid_index[1].z_index;

		i2 = ib_node->fluid_index[2].x_index;
		j2 = ib_node->fluid_index[2].y_index;
		k2 = ib_node->fluid_index[2].z_index;

		coeff0 = ib_node->fluid_coef[0];
		coeff1 = ib_node->fluid_coef[1];
		coeff2 = ib_node->fluid_coef[2];
		eta = coeff0 + coeff1 +coeff2;
		eta = pow(eta/1-eta, 0.5);
		eta = max(eta,1.0);
		eta = min(eta, 0.0);

		data[k][j][i] = coeff0*data[k0][j0][i0] + coeff1*data[k1][j1][i1] + coeff2*data[k2][j2][i2];
/*
		vel_corr = coeff0*data_old[k0][j0][i0] + coeff1*data_old[k1][j1][i1] + coeff2*data_old[k2][j2][i2];
		vel_corr = data_old[k][j][i] - vel_corr;
		data[k][j][i] = data[k][j][i] + eta*vel_corr;

		if (which_quantity == 'u') {
			dp_dxyz = (p_data[k][j][i]-p_data[k][j][i-1]) * idxyz[i-1];

			dp_dxyz0 = (p_data[k0][j0][i0]-p_data[k0][j0][i0-1]) * idxyz[i0-1];
			dp_dxyz1 = (p_data[k1][j1][i1]-p_data[k1][j1][i1-1]) * idxyz[i1-1];
			dp_dxyz2 = (p_data[k2][j2][i2]-p_data[k2][j2][i2-1]) * idxyz[i2-1];
		}
		else if (which_quantity == 'v') {
			dp_dxyz = (p_data[k][j][i]-p_data[k][j-1][i]) * idxyz[j-1];

			dp_dxyz0 = (p_data[k0][j0][i0]-p_data[k0][j0-1][i0]) * idxyz[j0-1];
			dp_dxyz1 = (p_data[k1][j1][i1]-p_data[k1][j1-1][i1]) * idxyz[j1-1];
			dp_dxyz2 = (p_data[k2][j2][i2]-p_data[k2][j2-1][i2]) * idxyz[j2-1];
		}
		else if (which_quantity == 'w') {
			dp_dxyz = (p_data[k][j][i]-p_data[k-1][j][i]) * idxyz[k-1];

			dp_dxyz0 = (p_data[k0][j0][i0]-p_data[k0-1][j0][i0]) * idxyz[k0-1];
			dp_dxyz1 = (p_data[k1][j1][i1]-p_data[k1-1][j1][i1]) * idxyz[k1-1];
			dp_dxyz2 = (p_data[k2][j2][i2]-p_data[k2-1][j2][i2]) * idxyz[k2-1];
		}
		press_corr = (coeff0+coeff1+coeff2)*dp_dxyz;
		press_corr -= -coeff0*dp_dxyz0 - coeff1*dp_dxyz1 - coeff2*dp_dxyz2;
		data[k][j][i] += press_corr*dtimeb;
*/

	} // for g
}


/******************************************************************************/
/*
 This function interpolates velocity at neighbhoring fluid nodes to find the
 velocity at the immersed node
 */
/******************************************************************************/
void Immersed_conc_interpolation(MAC_grid *grid, Parameters *params,
		Concentration *c) {

	int g, N;
	Immersed *c_immersed;
	ImmersedNode *ib_node;
	double ***data;
	int success;
	int i, j, k;
	int i0, j0, k0, i1, j1, k1, i2, j2, k2;
	double coeff0, coeff1, coeff2;
//    FILE *fid;
//    char bin_filename[50];


	c_immersed = grid->c_immersed;

#ifdef CONC_BQUICK
	data = c->data_temp;
#else
	data = c->data;
#endif

	// Total number of the immersed nodes
	N = c_immersed->N;


//    sprintf(bin_filename,"cinterimm_%d.dat",params->rank);
//    fid = fopen(bin_filename,"a");

	// Go through all the ib nodes and find the neighboring fluid nodes. Store
	// the indices and the coordinates
	for (g=0; g<N; g++) {

		// Get the current immersed node
		ib_node = Immersed_get_ib_node(c_immersed, g);

		i = ib_node->im_index.x_index;
		j = ib_node->im_index.y_index;
		k = ib_node->im_index.z_index;

		i0 = ib_node->fluid_index[0].x_index;
		j0 = ib_node->fluid_index[0].y_index;
		k0 = ib_node->fluid_index[0].z_index;

		i1 = ib_node->fluid_index[1].x_index;
		j1 = ib_node->fluid_index[1].y_index;
		k1 = ib_node->fluid_index[1].z_index;

		i2 = ib_node->fluid_index[2].x_index;
		j2 = ib_node->fluid_index[2].y_index;
		k2 = ib_node->fluid_index[2].z_index;

		coeff0 = ib_node->fluid_coef[0];
		coeff1 = ib_node->fluid_coef[1];
		coeff2 = ib_node->fluid_coef[2];

		data[k][j][i] = coeff0*data[k0][j0][i0] + coeff1*data[k1][j1][i1] + coeff2*data[k2][j2][i2];
// 		fprintf(fid,"node(%4d %3d %3d) node0(%4d %3d %3d) node1(%4d %3d %3d) node2(%4d %3d %3d) coef0= %+-8.5f coef1= %+-8.5f coef2= %+-8.5f d0= %+-8.5f d1= %+-8.5f d2= %+-8.5f \n",
//		i, j, k, i0, j0, k0, i1, j1, k1, i2, j2, k2, coeff0, coeff1, coeff2, data[k0][j0][i0],data[k1][j1][i1],data[k2][j2][i2]);

	} // for g

//	fclose(fid);
}




/******************************************************************************/
/*
 This function interpolates nut at neighbhoring fluid nodes to find nut at the
 immersed node so that nut at the wall is zero
 */
/******************************************************************************/
void Immersed_nut_interpolation(MAC_grid *grid, Parameters *params, Subgrid *smag) {

	int g, N;
	Immersed *nut_immersed;
	ImmersedNode *ib_node;
	double ***data;
	int success;
	int i, j, k;
	int i0, j0, k0, i1, j1, k1, i2, j2, k2;
	double coeff0, coeff1, coeff2;

	nut_immersed = grid->nut_immersed;

	data = smag->nut;

	// Total number of the immersed nodes
	N = nut_immersed->N;

	// Go through all the ib nodes and find the neighboring fluid nodes. Store
	// the indices and the coordinates
	for (g=0; g<N; g++) {

		// Get the current immersed node
		ib_node = Immersed_get_ib_node(nut_immersed, g);

		i = ib_node->im_index.x_index;
		j = ib_node->im_index.y_index;
		k = ib_node->im_index.z_index;

		i0 = ib_node->fluid_index[0].x_index;
		j0 = ib_node->fluid_index[0].y_index;
		k0 = ib_node->fluid_index[0].z_index;

		i1 = ib_node->fluid_index[1].x_index;
		j1 = ib_node->fluid_index[1].y_index;
		k1 = ib_node->fluid_index[1].z_index;

		i2 = ib_node->fluid_index[2].x_index;
		j2 = ib_node->fluid_index[2].y_index;
		k2 = ib_node->fluid_index[2].z_index;

		coeff0 = ib_node->fluid_coef[0];
		coeff1 = ib_node->fluid_coef[1];
		coeff2 = ib_node->fluid_coef[2];

		data[k][j][i] = coeff0*data[k0][j0][i0] + coeff1*data[k1][j1][i1] + coeff2*data[k2][j2][i2];

	} // for g
}


/******************************************************************************/
/*
 This function sets concentration value in the solid part to the
 value in the immersed nodes. This is done to avoid "plume"
 */
/******************************************************************************/
void Immersed_conc_solid(Concentration *c, Cart3d_bag *data_bag) {

	double ***data;
	Immersed *c_immersed;
	int NX, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	Indices G_s, G_e, W_e;
	int j_index;
	int sj_start, sj_end, sj;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	c_immersed = grid->c_immersed;

#ifdef CONC_BQUICK
	data = c->data_temp;
#else
	data = c->data;
#endif

	NX = grid->NX;
	NZ = grid->NZ;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	// Since the 2D arrays are assumed to have [Y][X] index order, pass always
	// the min, max indices based on this rule
	// Start and End indices of the 2D array on the current processor
	G_s.x_index = Is;
	G_s.y_index = Ks;

	G_e.x_index = Ie;
	G_e.y_index = Ke;

	// Total number of grid point in the W_ array
	W_e.x_index = NX;
	W_e.y_index = NZ;

	for (k=Ks; k<Ke; k++) {

		for (i=Is; i<Ie; i++) {

			c->G_conc_immersed[k][i] = 0.0;

			// Index of first interior fluid node (on the bottom boundary
			j_index = grid->interface_y_index[k][i];

			// Get the value of concentration at the immersed node
			if ( ( j_index >= Js) && (j_index < Je) && (j_index!=0) ){

				c->G_conc_immersed[k][i] = data[k][j_index][i];

			} // if

		} // for
	} /* for k*/

	Communication_reduce_2D_arrays(c->G_conc_immersed, c->W_conc_immersed,
			&G_s, &G_e, &W_e, REDUCE_TO_ALL, data_bag);

	for (k=Ks; k<Ke; k++) {

		for (i=Is; i<Ie; i++) {

			// Index of first interior fluid node (on the bottom boundary
			j_index = grid->interface_y_index[k][i];

			if (j_index ==0) break;
			if (Js >= j_index) break;
			if (Je <= j_index - 10) break;

			sj_start = max(j_index-10, Js);
			sj_end = min(j_index, Je);

			for (sj=sj_start; sj < sj_end; sj++) {

				data[k][sj][i] = c->W_conc_immersed[k][i];
			}

		} // for
	} /* for k*/



}

/******************************************************************************/
