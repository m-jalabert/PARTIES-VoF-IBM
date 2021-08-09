#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"
#include "Memory.h"
#include "Communication.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

MPI_COMM comm3d;

/******************************************************************************/
/*
 This function reduces (sums) all the global arrays to the W_array.
 set status:
 1- status: REDUCE_TO_MASTER : result will be stored only on processor zero
 2- status: REDUCE_TO_ALL    : result will be stores on all processors
 */
/******************************************************************************/
void Communication_reduce_2D_arrays(double **G_send_array, double **W_recv_array,
		Indices *G_s, Indices *G_e, Indices *W_e, int status, Cart3d_bag *data_bag) {

	double T1, T2;
	T1 = MPI_Wtime();

	int i, j;
	int index;
	int x_s, x_e;
	int y_s, y_e;
	int X_e, Y_e;
	int W_nt;
	double *send_buffer;
	double *recv_buffer;

	Parameters *params = data_bag -> params;

	// Start index [y_s,x_s]
	// End  index  [y_e-1, x_e-1]

	x_s = G_s->x_index;
	y_s = G_s->y_index;

	x_e = G_e->x_index;
	y_e = G_e->y_index;

	// Last index (+1) of the 2D world array
	X_e = W_e->x_index;
	Y_e = W_e->y_index;

	// Total number of grid points in the World array
	W_nt = X_e * Y_e;


	send_buffer = Memory_allocate_1D_array(GVG_DOUBLE, W_nt);

	// First pack the data into a 1D array
	for (j=y_s; j<y_e; j++) {

		for (i=x_s; i<x_e; i++) {

			// World index of the point in the array
			index = j*X_e + i;
			send_buffer[index] = G_send_array[j][i];
		} // for i
	} // for

	if (status == REDUCE_TO_ALL) {

		recv_buffer = Memory_allocate_1D_array(GVG_DOUBLE, W_nt);

		// This will send to all processors after computing deposit_height
		MPI_Allreduce( (void *)send_buffer, (void *)recv_buffer, W_nt, MPI_DOUBLE, MPI_SUM, PCW);

	}
	else if (status == REDUCE_TO_MASTER) {

		if (params->rank == MASTER) {
			recv_buffer = Memory_allocate_1D_array(GVG_DOUBLE, W_nt);
		}
		else {

			recv_buffer = NULL;
		}

		// Store the final result on processor zero
		MPI_Reduce( (void *)send_buffer, (void *)recv_buffer, W_nt, MPI_DOUBLE,
		           MPI_SUM, MASTER, PCW);
	}

	// Paste back the reduced data to the W_array in the 2D array format
	if (recv_buffer != NULL ) {

		index = 0;
		// Now, unpack the 1D array back into the 2D array
		for (j=0; j<Y_e; j++) {

			for (i=0; i<X_e; i++) {

				W_recv_array[j][i] = recv_buffer[index];
				index++;
			} // for i
		} // for j
	} // if recv_buffer


	free(send_buffer);
	if (recv_buffer != NULL) {

		free(recv_buffer);
	}

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_comm_2D_reduce += T2 - T1;

}




/******************************************************************************/
/*
 This function, communicate between processors and update the ghost nodes on the
 local vector on current processor
 */
/******************************************************************************/
void Communication_update_ghost_nodes_flow_variable(double ***data,
		char component, int pnodes, Cart3d_bag *data_bag) {

	Communication_update_ghost_nodes_x(data, component, pnodes, data_bag);

	Communication_update_ghost_nodes_y(data, component, pnodes, data_bag);

	Communication_update_ghost_nodes_z(data, component, pnodes, data_bag);
}




/******************************************************************************/
/*
 This functions updates the ghost nodes only along x-direction
 */
/******************************************************************************/
void Communication_update_ghost_nodes_x(double ***data, char component,
		int pnodes, Cart3d_bag *data_bag) {

	double T1, T2;
	T1 = MPI_Wtime();

	int i, j, k;

	int ierr ;
	MPI_Status status_first;
	int index, index1, index2;
	int sendtag, recvtag, sendcount, recvcount;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Allocate memory for send and receive buffers
	double *send_buffer = params -> send_buffer;
	double *recv_buffer = params -> recv_buffer;

	// Start index of bottom-left-back corner on current processor
	int Is = grid -> G_Is;

	// End index of top-right-front corner on current processor
	int Ie = grid -> G_Ie;

	// Extra send/receive layer required for periodic boundaries
	int pnodes_send = pnodes;
	int pnodes_recv = pnodes;

#ifdef XPERIODIC
	int NX = grid -> NX;

	if (component == 'u' || component == U_VELOCITY || component == U_VELOCITY_PERTURBATION ) {
		if (Is == 0) {
			Is = 1;
			pnodes_recv = pnodes + 1;
		}
		if (Ie == NX)
			pnodes_send = pnodes + 1;
	}
	else {
		if (Ie == NX)
			Ie = NX - 1;
	}
#endif

	// Start index of bottom-left-back corner on current processor including
	// ghost nodes
	int Js_g = grid->G_Js - pnodes;
	int Ks_g = grid->G_Ks - pnodes;

	// End index of top-right-front corner on current processor including ghost
	// nodes
	int Je_g = grid->G_Je + pnodes;
	int Ke_g = grid->G_Ke + pnodes;

	//--------------------------------------------------------------------------
	// Send to upper X
	//--------------------------------------------------------------------------
	sendtag = 1;
	recvtag = 1;

	sendcount = (Je_g - Js_g) * (Ke_g - Ks_g) * pnodes_send;
	recvcount = (Je_g - Js_g) * (Ke_g - Ks_g) * pnodes_recv;

	// Copy the data that needs to be sent to send_buffer
	if (params->npxplus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * (Je_g - Js_g) * pnodes_send;
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * pnodes_send;
				for (i = 0; i < pnodes_send; i++) {
					index = index2 + i;
					send_buffer[index] = data[k][j][Ie-i-1];
				}
			}
		}
	}

	ierr = MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE,
	                    params->npxplus, sendtag,
	                    recv_buffer, recvcount, MPI_DOUBLE,
	                    params->npxminus, recvtag, PCW, &status_first);

	// Copy the received data
	if (params->npxminus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * (Je_g - Js_g) * pnodes_recv;
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * pnodes_recv;
				for (i = 0; i < pnodes_recv; i++) {
					index = index2 + i;
					data[k][j][Is-i-1] = recv_buffer[index];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Send to lower X
	//--------------------------------------------------------------------------
	sendtag = 2;
	recvtag = 2;

	sendcount = (Je_g - Js_g) * (Ke_g - Ks_g) * pnodes;
	recvcount = (Je_g - Js_g) * (Ke_g - Ks_g) * pnodes;

	// Copy the data that needs to be sent to send_buffer
	if (params->npxminus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * (Je_g - Js_g) * pnodes;
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * pnodes;
				for (i = 0; i < pnodes; i++) {
					index = index2 + i;
					send_buffer[index] = data[k][j][Is+i];
				}
			}
		}
	}

	ierr = MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE,
	                    params->npxminus, sendtag,
	                    recv_buffer, recvcount, MPI_DOUBLE,
	                    params->npxplus, recvtag, PCW, &status_first);

	// Copy the received data
	if (params->npxplus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * (Je_g - Js_g) * pnodes;
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * pnodes;
				for (i = 0; i < pnodes; i++) {
					index = index2 + i;
					data[k][j][Ie+i] = recv_buffer[index];
				}
			}
		}
	}

	T2 = MPI_Wtime();
	data_bag -> timer -> Wtime_comm_3D += T2 - T1;

	return;
}




/******************************************************************************/
/*
 This functions updates the ghost nodes only along y-direction
 */
/******************************************************************************/
void Communication_update_ghost_nodes_y(double ***data, char component,int pnodes, Cart3d_bag *data_bag) {

	double T1, T2;
	T1 = MPI_Wtime();

	int i, j, k;

	int ierr ;
	MPI_Status status_first;
	int index, index1, index2;
	int sendtag, recvtag, sendcount, recvcount;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Allocate memory for send and receive buffers
	double *send_buffer = params -> send_buffer;
	double *recv_buffer = params -> recv_buffer;

	// Start index of bottom-left-back corner on current processor
	int Js = grid -> G_Js;

	// End index of top-right-front corner on current processor
	int Je = grid -> G_Je;


	// Extra send/receive layer required for periodic boundaries
	int pnodes_send = pnodes;
	int pnodes_recv = pnodes;


	#ifdef YPERIODIC
	int NY = grid -> NY;

	if (component == 'v' || component == V_VELOCITY) {
		if (Js == 0) {
			Js = 1;
			pnodes_recv = pnodes + 1;
		}
		if (Je == NY)
			pnodes_send = pnodes + 1;
	}
	else {
		if (Je == NY)
			Je = NY - 1;
	}
	#endif


	// Start index of bottom-left-back corner on current processor including
	// ghost nodes
	int Is_g = grid->G_Is - pnodes;
	int Ks_g = grid->G_Ks - pnodes;

	// End index of top-right-front corner on current processor including ghost
	// nodes
	int Ie_g = grid->G_Ie + pnodes;
	int Ke_g = grid->G_Ke + pnodes;





	//--------------------------------------------------------------------------
	// Send to upper Y
	//--------------------------------------------------------------------------
	sendtag = 5;
	recvtag = 5;

	sendcount = (Ie_g - Is_g) * (Ke_g - Ks_g) * pnodes_send;
	recvcount = (Ie_g - Is_g) * (Ke_g - Ks_g) * pnodes_recv;

	// Copy the data that needs to be sent to send_buffer
	if (params->npyplus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * pnodes_send * (Ie_g - Is_g);
			for (i = Is_g; i < Ie_g; i++) {
				index2 = index1 + (i - Is_g);
				for (j = 0; j < pnodes_send; j++) {
					index = index2 + j * (Ie_g - Is_g);
					send_buffer[index] = data[k][Je-j-1][i];
				}
			}
		}
	}

	ierr = MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE,
	                    params->npyplus, sendtag,
	                    recv_buffer, recvcount, MPI_DOUBLE,
	                    params->npyminus, recvtag, PCW, &status_first);

	// Copy the received data
	if (params->npyminus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * pnodes_recv * (Ie_g - Is_g);
			for (i = Is_g; i < Ie_g; i++) {
				index2 = index1 + (i - Is_g);
				for (j = 0; j < pnodes_recv; j++) {
					index = index2 + j * (Ie_g - Is_g);
					data[k][Js-j-1][i] = recv_buffer[index];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Send to lower Y
	//--------------------------------------------------------------------------
	sendtag = 6;
	recvtag = 6;

	sendcount = (Ie_g - Is_g) * (Ke_g - Ks_g) * pnodes;
	recvcount = (Ie_g - Is_g) * (Ke_g - Ks_g) * pnodes;

	// Copy the data that needs to be sent to send_buffer
	if (params->npyminus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * pnodes * (Ie_g - Is_g);
			for (i = Is_g; i < Ie_g; i++) {
				index2 = index1 + (i - Is_g);
				for (j = 0; j < pnodes; j++) {
					index = index2 + j * (Ie_g - Is_g);
					send_buffer[index] = data[k][Js+j][i];
				}
			}
		}
	}

	ierr = MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE,
	                    params->npyminus, sendtag,
	                    recv_buffer, recvcount, MPI_DOUBLE,
	                    params->npyplus, recvtag, PCW, &status_first);

	if (params->npyplus != MPI_PROC_NULL) {
		for (k = Ks_g; k < Ke_g; k++) {
			index1 = (k - Ks_g) * pnodes * (Ie_g - Is_g);
			for (i = Is_g; i < Ie_g; i++) {
				index2 = index1 + (i - Is_g);
				for (j = 0; j < pnodes; j++) {
					index = index2 + j * (Ie_g - Is_g);
					data[k][Je+j][i] = recv_buffer[index];
				}
			}
		}
	}

#ifdef SHEARED_PERIODIC

	double shear_sigma_u = params->shear_sigma_u;
	double shear_sigma_conc = params->shear_sigma_conc;

	// if we are at the bottom processor copy data into the slab and distribute it to the others
	double **exchange_slab = data_bag->grid->exchange_slab;
	double *xc=grid->xc;
	double *xu=grid->xu;

	int yproc = params-> yproccoord;
	int NX = grid->NX;
	int NZ = grid->NZ;
	int NPY = params->NPY;

	int Lx = params ->Lx;
	int Ly = params ->Ly;

	double xcorr;
	double time=params->time;
	double xmin = params->xmin;

	double ipos_index;
	int ipos_west;
	int ipos_east;
	// exclude the half cell added
	int i_start = grid -> G_Is;
	int k_start = grid -> G_Ks;
	int i_end =  grid->G_Ie;
	int k_end =  grid->G_Ke;

	double hc = grid->dx_c[1];
	double hu = grid->dx_u[1];

	double dh_west;
	DSET_ZERO(&exchange_slab[0][0], NX*NZ );

	double conc_diff;
	double u_diff;
/*
 	i_start = grid -> G_Is - 1;
	k_start = grid -> G_Ks - 1;
	i_end =  grid->G_Ie + 1;
	k_end =  grid->G_Ke + 1;
	if (i_end==NX+1)
	i_end = NX;
	if (k_end==NZ+1)
	k_end = NZ;

	if (i_start==-1)
	i_start = 0;
	if (k_start==-1)
	k_start = 0;
*/
	if ( (component == CONCENTRATION) || (component == CONCENTRATION_PERTURBATION) ){

		conc_diff = shear_sigma_conc*Ly;
		if(component == CONCENTRATION_PERTURBATION) conc_diff = 0;


		if(yproc == 0) {  // We are at the bottom let's copy the data to the exchange
			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][0][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that





		if(yproc == NPY-1) { // We are at the top process
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr= remainder(xc[i] - shear_sigma_u*Ly*time,Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr = xcorr + Lx;  // I assume that xmin=0!!!

							ipos_index = (xcorr - hc/2)/hc;
							ipos_west = floor(ipos_index); // we are also west of the fist cell center and east of the cell center
							ipos_east = ipos_west+1;
					  		dh_west = xcorr-xc[ipos_west];
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

						 // this is the distance to the western cell center


					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the position on the no_ghost grid


						 data[k][NY-1][i] = exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc + conc_diff ;
					}
			}

		}


		DSET_ZERO(&exchange_slab[0][0], NX*NZ );


		if(yproc == NPY-1) {  // We are at the top let's copy the data to the exchange

			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][NY-2][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that



		if(yproc ==0) { // We are at the bottom processes
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr= remainder(xc[i]+shear_sigma_u*Ly*time, Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr= xcorr+Lx;
						 ipos_index = (xcorr-hc/2)/hc;


						 ipos_west = floor(ipos_index); // we are also west of the fist cell center and east of the cell center
						 ipos_east = ipos_west + 1;
						 dh_west= xcorr-xc[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }
						 //printf("ipos_west, ipos_east = %d, %d\n" ,ipos_west,ipos_east);

						/* if (ipos_west==0){
							 if(k==0){
							 	printf("exchange_slab west = %2.5f,\n" , exchange_slab[k][ipos_west]);
								printf("exchange_slab east= %2.5f,\n" , exchange_slab[k][ipos_east]);

						 }}*/
					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the positon on the no_ghost grid


						 data[k][-1][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc - conc_diff ;
					}
			}

		}

	}

	// Communicates the appropriate BC for the shear velocity

	if ( component == U_VELOCITY || component == U_VELOCITY_PERTURBATION ){

		u_diff = shear_sigma_u*Ly;
		if(component == U_VELOCITY_PERTURBATION) u_diff = 0.0;


		if(yproc == 0) {  // We are at the bottom let's copy the data to the exchange
			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][0][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that
		if(yproc == NPY-1) { // We are at the top process
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr= remainder(xu[i]-shear_sigma_u*Ly*time,Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab

						 if(xcorr<0) xcorr= xcorr + Lx;  // I assume that xmin=0!!!

						 ipos_index = (xcorr)/hu;

						 ipos_west = floor(ipos_index); // we are also west of the fist cell center and east of the cell center
						 ipos_east = ipos_west+1;
						 dh_west = xcorr - xu[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the position on the no_ghost grid


						 data[k][NY-1][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hu + u_diff ;
					}
			}

		}
		//printf("I transfered the u-data to the top and udiff = %2.4f\n",u_diff);

		DSET_ZERO(&exchange_slab[0][0], NX*NZ );
		if(yproc == NPY-1) {  // We are at the top let's copy the data to the exchange

			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][NY-2][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that


		if(yproc ==0) { // We are at the bottom processes
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr= remainder(xu[i]+shear_sigma_u*Ly*time, Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr = xcorr + Lx;
						 ipos_index = (xcorr)/hu;


						 ipos_west=floor(ipos_index ); // we are also west of the fist cell center and east of the cell center
						 ipos_east=ipos_west+1;
						 dh_west= xcorr-xu[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the positon on the no_ghost grid


						 data[k][-1][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hu - u_diff ;
					}
			}

		}

	}

	if ( component == V_VELOCITY){

		if(yproc == 0) {  // We are at the bottom let's copy the data to the exchange
			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][0][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that
		if(yproc == NPY-1) {
// We are at the top process
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr= remainder(xc[i] - shear_sigma_u*Ly*time,Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr = xcorr+Lx;  // I assume that xmin=0!!!

						 ipos_index = (xcorr - hc/2)/hc;



						 ipos_west=floor(ipos_index ); // we are also west of the fist cell center and east of the cell center
						 ipos_east=ipos_west+1;
						 dh_west= xcorr-xc[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the position on the no_ghost grid

						 data[k][NY-1][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc;
					}
			}

		}


		DSET_ZERO(&exchange_slab[0][0], NX*NZ );
		if(yproc == NPY-1) {  // We are at the top let's copy the data to the exchange

			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][NY-2][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that


		if(yproc ==0) { // We are at the bottom processes
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr = remainder(xc[i] + shear_sigma_u*Ly*time, Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr= xcorr+Lx;
						 ipos_index = (xcorr-hc/2)/hc;


						 ipos_west=floor(ipos_index ); // we are also west of the fist cell center and east of the cell center
						 ipos_east=ipos_west+1;
						 dh_west= xcorr-xc[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the positon on the no_ghost grid


						 data[k][-1][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc;
					}
			}

		}

		if(yproc == 0) {  // We are at the bottom let's copy the data to the exchange
			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][1][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that
		if(yproc == NPY-1) { // We are at the top process
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr= remainder(xc[i] - shear_sigma_u*Ly*time,Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr = xcorr+Lx;  // I assume that xmin=0!!!

						 ipos_index = (xcorr - hc/2)/hc;



						 ipos_west=floor(ipos_index ); // we are also west of the fist cell center and east of the cell center
						 ipos_east=ipos_west+1;
						 dh_west= xcorr-xc[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the position on the no_ghost grid


						 data[k][NY][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc;
					}
			}

		}


		DSET_ZERO(&exchange_slab[0][0], NX*NZ );
		if(yproc == NPY-1) {  // We are at the top let's copy the data to the exchange

			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][NY-3][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that


		if(yproc ==0) { // We are at the bottom processes
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr = remainder(xc[i] + shear_sigma_u*Ly*time, Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr= xcorr+Lx;
						 ipos_index = (xcorr-hc/2)/hc;


						 ipos_west=floor(ipos_index ); // we are also west of the fist cell center and east of the cell center
						 ipos_east=ipos_west+1;
						 dh_west= xcorr-xc[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the positon on the no_ghost grid


						 data[k][-2][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc;
					}
			}

		}

	}


	if (component == W_VELOCITY){

		if(yproc == 0) {  // We are at the bottom let's copy the data to the exchange
			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][0][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that
		if(yproc == NPY-1) { // We are at the top process
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr= remainder(xc[i] - shear_sigma_u*Ly*time,Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr = xcorr+Lx;  // I assume that xmin=0!!!

						 ipos_index = (xcorr - hc/2)/hc;



						 ipos_west=floor(ipos_index ); // we are also west of the fist cell center and east of the cell center
						 ipos_east=ipos_west+1;
						 dh_west= xcorr-xc[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the position on the no_ghost grid


						 data[k][NY-1][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc;
					}
			}

		}


		DSET_ZERO(&exchange_slab[0][0], NX*NZ );
		if(yproc == NPY-1) {  // We are at the top let's copy the data to the exchange

			for (k = k_start; k < k_end; k++) {
				for (i = i_start; i < i_end; i++) {

					exchange_slab[k][i] = data[k][NY-2][i];

				}
			}
		}
			MPI_Allreduce(MPI_IN_PLACE, &exchange_slab[0][0], NX*NZ, MPI_DOUBLE, MPI_SUM, PCW);
			// Now everyone has the bottom slab, this should be optimized later by a more specific communicator because only the top need to know that


		if(yproc ==0) { // We are at the bottom processes
			for (k = k_start; k < k_end; k++) {
					for (i = i_start; i < i_end; i++) {
						 xcorr = remainder(xc[i] + shear_sigma_u*Ly*time, Lx);  // this is the corresponding position to xc[i] of data with respect to the exchange slab
						 if(xcorr<0) xcorr= xcorr+Lx;
						 ipos_index = (xcorr-hc/2)/hc;


						 ipos_west=floor(ipos_index ); // we are also west of the fist cell center and east of the cell center
						 ipos_east=ipos_west+1;
						 dh_west= xcorr-xc[ipos_west]; // this is the distance to the western cell center
						 if (ipos_west < 0){
							 ipos_west = NX - 2;
						 }

					//	 if (ipos_west < 0) ipos_west= NX-2;
					//	 if (ipos_east > NX-2) ipos_east=0;
						 // no we have the positon on the no_ghost grid


						 data[k][-1][i] =  exchange_slab[k][ipos_west]+ dh_west*(exchange_slab[k][ipos_east]-exchange_slab[k][ipos_west])/hc;
					}
			}

		}

	}

#endif

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_comm_3D += T2 - T1;

	return;
}




/******************************************************************************/
/*
 This functions updates the ghost nodes only along z-direction
 */
/******************************************************************************/
void Communication_update_ghost_nodes_z(double ***data, char component,
		int pnodes, Cart3d_bag *data_bag) {

	double T1, T2;
	T1 = MPI_Wtime();

	int i, j, k;

	int ierr ;
	MPI_Status status_first;
	int index, index1, index2;
	int sendtag, recvtag, sendcount, recvcount;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Allocate memory for send and receive buffers
	double *send_buffer = params -> send_buffer;
	double *recv_buffer = params -> recv_buffer;

	// Start index of bottom-left-back corner on current processor
	int Ks = grid -> G_Ks;

	// End index of top-right-front corner on current processor
	int Ke = grid -> G_Ke;

	// Extra send/receive layer required for periodic boundaries
	int pnodes_send = pnodes;
	int pnodes_recv = pnodes;

#ifdef ZPERIODIC
	int NZ = grid -> NZ;

	if (component == 'w' || component == W_VELOCITY) {
		if (Ks == 0) {
			Ks = 1;
			pnodes_recv = pnodes + 1;
		}
		if (Ke == NZ)
			pnodes_send = pnodes + 1;
	}
	else {
		if (Ke == NZ)
			Ke = NZ - 1;
	}
#endif

	// Start index of bottom-left-back corner on current processor including
	// ghost nodes
	int Is_g = grid->G_Is - pnodes;
	int Js_g = grid->G_Js - pnodes;

	// End index of top-right-front corner on current processor including ghost
	// nodes
	int Ie_g = grid->G_Ie + pnodes;
	int Je_g = grid->G_Je + pnodes;

	//--------------------------------------------------------------------------
	// Send to upper Z
	//--------------------------------------------------------------------------
	sendtag = 3;
	recvtag = 3;

	sendcount = (Ie_g - Is_g) * (Je_g - Js_g) * pnodes_send;
	recvcount = (Ie_g - Is_g) * (Je_g - Js_g) * pnodes_recv;

	// Copy the data that needs to be sent to send_buffer
	if (params->npzplus != MPI_PROC_NULL) {
		for (k = 0; k < pnodes_send; k++) {
			index1 = k * (Je_g - Js_g) * (Ie_g - Is_g);
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * (Ie_g - Is_g);
				for (i = Is_g; i < Ie_g; i++) {
					index = index2 + i - Is_g;
					send_buffer[index] = data[Ke-k-1][j][i];
				}
			}
		}
	}

	ierr = MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE,
	                    params->npzplus, sendtag,
	                    recv_buffer, recvcount, MPI_DOUBLE,
	                    params->npzminus, recvtag, PCW, &status_first);

	// Copy the received data
	if (params->npzminus != MPI_PROC_NULL) {
		for (k = 0; k < pnodes_recv; k++) {
			index1 = k * (Je_g - Js_g) * (Ie_g - Is_g);
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * (Ie_g - Is_g);
				for (i = Is_g; i < Ie_g; i++) {
					index = index2 + i - Is_g;
					data[Ks-k-1][j][i] = recv_buffer[index];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Send to lower Z
	//--------------------------------------------------------------------------
	sendtag = 4;
	recvtag = 4;

	sendcount = (Ie_g - Is_g) * (Je_g - Js_g) * pnodes;
	recvcount = (Ie_g - Is_g) * (Je_g - Js_g) * pnodes;

	// Copy the data that needs to be sent to send_buffer
	if (params->npzminus != MPI_PROC_NULL) {
		for (k = 0; k < pnodes; k++) {
			index1 = k * (Je_g - Js_g) * (Ie_g - Is_g);
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * (Ie_g - Is_g);
				for (i = Is_g; i < Ie_g; i++) {
					index = index2 + i - Is_g;
					send_buffer[index] = data[Ks+k][j][i];
				}
			}
		}
	}

	ierr = MPI_Sendrecv(send_buffer, sendcount, MPI_DOUBLE,
	                    params->npzminus, sendtag,
	                    recv_buffer, recvcount, MPI_DOUBLE,
	                    params->npzplus, recvtag, PCW, &status_first);

	// Copy the received data
	if (params->npzplus != MPI_PROC_NULL) {
		for (k = 0; k < pnodes; k++) {
			index1 = k * (Je_g - Js_g) * (Ie_g - Is_g);
			for (j = Js_g; j < Je_g; j++) {
				index2 = index1 + (j - Js_g) * (Ie_g - Is_g);
				for (i = Is_g; i < Ie_g; i++) {
					index = index2 + i - Is_g;
					data[Ke+k][j][i] = recv_buffer[index];
				}
			}
		}
	}

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_comm_3D += T2 - T1;

	return;
}




/******************************************************************************/
/*
 This function, creates new communicators for each subset of processors that
 share the x or y or z coordinates
 */
/******************************************************************************/
void Communication_new_xyz_communicator(MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int ierr;
	int *istartg, *jstartg, *kstartg, count;
	int *icoord, *jcoord, *kcoord, newcoord;
	int color, key;
	int inext, iprev, jnext, jprev, knext, kprev;
	int sendcount;

	FILE *fid;
	char bin_filename[50];

	int pnodes = params->ghost_nodes+1;
	int NX = grid->NX;
	int NY = grid->NY;
	int NZ = grid->NZ;

	// Start index of bottom-left-back corner on current processor
	int Is = grid->G_Is;
	int Js = grid->G_Js;
	int Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid->G_Ie;
	int Je = grid->G_Je;
	int Ke = grid->G_Ke;

	// Start index of bottom-left-back corner on current processor including
	// ghost nodes
	int Is_g = grid->L_Is;
	int Js_g = grid->L_Js;
	int Ks_g = grid->L_Ks;

	// End index of top-right-front corner on current processor including ghost
	// nodes
	int Ie_g = grid->L_Ie;
	int Je_g = grid->L_Je;
	int Ke_g = grid->L_Ke;


	sendcount = NY*NZ*(pnodes);
	sendcount = NX*NY*(pnodes);

	sendcount = max(sendcount, pnodes * (Je_g - Js_g) * (Ke_g - Ks_g) );
	sendcount = max(sendcount, pnodes * (Ie_g - Is_g) * (Ke_g - Ks_g) );
	sendcount = max(sendcount, pnodes * (Ie_g - Is_g) * (Je_g - Js_g) );

	params -> send_buffer = (double *)calloc(sendcount, sizeof(double));
	params -> recv_buffer = (double *)calloc(sendcount, sizeof(double));

	return;

}




/******************************************************************************/
/*
 This function is used to find the previous and next processors along the coordinate directions
 */
/******************************************************************************/
int proc_find(int *istartg,int *jstartg,int *kstartg, int size, int ilocal,
		int jlocal, int klocal) {

	int i, proc_id;

	proc_id = MPI_PROC_NULL;

	for (i=0;i<size;i++){
		if ( (istartg[i]==ilocal) && (jstartg[i]==jlocal) && (kstartg[i]==klocal) ) proc_id=i;
	}

	return proc_id;
}




/******************************************************************************/
/*
 This function sorts an integer array using bubble_sort downloaded from
 http://www.c.happycodings.com/Sorting_Searching/code4.html
 */
/******************************************************************************/
void bubble_sort(int a[], int size) {

	int switched = 1;
	int hold = 0;
	int i = 0;
	int j = 0;

	size -= 1;

	for(i = 0; i < size && switched; i++) {
		switched = 0;
		for(j = 0; j < size - i; j++) {
			if(a[j] > a[j+1]) {
				switched = 1;
				hold = a[j];
				a[j] = a[j + 1];
				a[j + 1] = hold;
			}
		}
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void xzperiodic_uvel_ave(Cart3d_bag *data_bag,Velocity *u, Velocity *v, Velocity *w,
		MAC_grid *grid, Parameters *params){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fidu, *fidv, *fidw;

	double ***u_data, ***v_data, ***w_data;

	double *udata_ave, *vdata_ave, *wdata_ave;
	double *udata_aveg, *vdata_aveg, *wdata_aveg;
	double *uudata_ave, *vvdata_ave, *wwdata_ave, *uvdata_ave;
	double *uudata_aveg, *vvdata_aveg, *wwdata_aveg, *uvdata_aveg;

	char bin_filename[50];
	double tetime;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	u_data = u->data;
	v_data = v->data;
	w_data = w->data;

	udata_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	vdata_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	wdata_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	udata_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	vdata_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	wdata_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	uudata_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	vvdata_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	wwdata_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	uvdata_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	uudata_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	vvdata_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	wwdata_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	uvdata_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		udata_ave[j] = 0.;
		vdata_ave[j] = 0.;
		wdata_ave[j] = 0.;
		udata_aveg[j] = 0;
		vdata_aveg[j] = 0.;
		wdata_aveg[j] = 0.;

		uudata_ave[j] = 0.;
		vvdata_ave[j] = 0.;
		wwdata_ave[j] = 0.;
		uvdata_ave[j] = 0;
		uudata_aveg[j] = 0.;
		vvdata_aveg[j] = 0.;
		wwdata_aveg[j] = 0.;
		uvdata_aveg[j] = 0;
	}

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Is = max(Is, 1);
	Ie = min(Ie, NX);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				udata_ave[j] = udata_ave[j] + u_data[k][j][i];
				uudata_ave[j] = uudata_ave[j] + u_data[k][j][i]*u_data[k][j][i];
			}
		}
	}

	ierr = MPI_Reduce (udata_ave,  udata_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (uudata_ave, uudata_aveg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);


	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				vdata_ave[j] = vdata_ave[j] + v_data[k][j][i];
				vvdata_ave[j] = vvdata_ave[j] + v_data[k][j][i]*v_data[k][j][i];
			}
		}
	}

	ierr = MPI_Reduce (vdata_ave,  vdata_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (vvdata_ave, vvdata_aveg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ks = max(Ks, 1);
	Ke = min(Ke, NZ);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				wdata_ave[j] = wdata_ave[j] + w_data[k][j][i];
				wwdata_ave[j] = wwdata_ave[j] + w_data[k][j][i]*w_data[k][j][i];
			}
		}
	}

	ierr = MPI_Reduce (wdata_ave,  wdata_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (wwdata_ave, wwdata_aveg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	if (params->rank == 0) {

		tetime=params->time/params->output_time_interval + 4e-1;
		sprintf(bin_filename,"uvelave_t%d.dat",(int)tetime);
		fidu = fopen(bin_filename,"w");
		sprintf(bin_filename,"vvelave_t%d.dat",(int)tetime);
		fidv = fopen(bin_filename,"w");
		sprintf(bin_filename,"wvelave_t%d.dat",(int)tetime);
		fidw = fopen(bin_filename,"w");

		for (j=0;j<NY;j++){
			udata_aveg[j] = udata_aveg[j]/(NX-1.0)/(NZ-1.0);
			uudata_aveg[j] = uudata_aveg[j]/(NX-1.0)/(NZ-1.0);
			uudata_aveg[j] = uudata_aveg[j]-udata_aveg[j]*udata_aveg[j];

			vdata_aveg[j] = vdata_aveg[j]/(NX-1.0)/(NZ-1.0);
			vvdata_aveg[j] = vvdata_aveg[j]/(NX-1.0)/(NZ-1.0);
			vvdata_aveg[j] = vvdata_aveg[j]-vdata_aveg[j]*vdata_aveg[j];

			wdata_aveg[j] = wdata_aveg[j]/(NX-1.0)/(NZ-1.0);
			wwdata_aveg[j] = wwdata_aveg[j]/(NX-1.0)/(NZ-1.0);
			wwdata_aveg[j] = wwdata_aveg[j]-wdata_aveg[j]*wdata_aveg[j];

			fprintf(fidu,"%20.12e %20.12e %20.12e \n",grid->yc[j],udata_aveg[j],uudata_aveg[j]);
			fprintf(fidv,"%e %e %e \n",grid->yv[j],vdata_aveg[j],vvdata_aveg[j]);
			fprintf(fidw,"%e %e %e \n",grid->yc[j],wdata_aveg[j],wwdata_aveg[j]);
		}
		fclose(fidu);
		fclose(fidv);
		fclose(fidw);

	}

	return;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void xzperiodic_saltsediment_ave(Cart3d_bag *data_bag, Concentration **c,
		MAC_grid *grid, Parameters *params){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fid;

	double ***salt, ***sedmt;

	double *salt_ave, *sedmt_ave;
	double *salt_aveg, *sedmt_aveg;
	double *salt_rms, *sedmt_rms;
	double *salt_rmsg, *sedmt_rmsg;

	char bin_filename[50];
	double tetime;
	double zflow,xflow;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	salt = c[0]->data;
	sedmt = c[0]->data;

	salt_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sedmt_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	salt_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sedmt_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	salt_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sedmt_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	salt_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sedmt_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		salt_ave[j] = 0.;
		sedmt_ave[j] = 0.;
		salt_aveg[j] = 0.;
		sedmt_aveg[j] = 0.;

		salt_rms[j] = 0.;
		sedmt_rms[j] = 0.;
		salt_rmsg[j] = 0.;
		sedmt_rmsg[j] = 0.;
	}

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				salt_ave[j] = salt_ave[j] + salt[k][j][i];
				sedmt_ave[j] = sedmt_ave[j] + sedmt[k][j][i];
				salt_rms[j] = salt_rms[j] + salt[k][j][i]*salt[k][j][i];
				sedmt_rms[j] = sedmt_rms[j] + sedmt[k][j][i]*sedmt[k][j][i];
			}
		}
	}
	ierr = MPI_Reduce (salt_ave,  salt_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sedmt_ave, sedmt_aveg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (salt_rms,  salt_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sedmt_rms, sedmt_rmsg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	if (params->rank == 0) {

		tetime=params->time/params->output_time_interval + 4e-1;
		sprintf(bin_filename,"concave_t%d.dat",(int)tetime);
		fid = fopen(bin_filename,"w");

		for (j=0;j<NY;j++){
			salt_aveg[j] = salt_aveg[j]/(NX-1.0)/(NZ-1.0);
			sedmt_aveg[j] = sedmt_aveg[j]/(NX-1.0)/(NZ-1.0);

			salt_rmsg[j] = salt_rmsg[j]/(NX-1.0)/(NZ-1.0);
			sedmt_rmsg[j] = sedmt_rmsg[j]/(NX-1.0)/(NZ-1.0);

			salt_rmsg[j] = salt_rmsg[j] - salt_aveg[j]*salt_aveg[j];
			sedmt_rmsg[j] = sedmt_rmsg[j] - sedmt_aveg[j]*sedmt_aveg[j];

			fprintf(fid,"%e %e %e %e %e\n",grid->yc[j],salt_aveg[j],sedmt_aveg[j],salt_rmsg[j],sedmt_rmsg[j]);
		}

		fclose(fid);
	}

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void xzperiodic_thermal_ave2(Concentration *c, MAC_grid *grid,
		Parameters *params, int flag){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fid;

	double ***salt;

	double *salt_ave;
	double *salt_aveg;
	double *salt_rms;
	double *salt_rmsg;

	char bin_filename[50];
	double tetime;
	double zflow,xflow;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	salt = c->data;

	salt_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	salt_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	salt_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	salt_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		salt_ave[j] = 0.;
		salt_aveg[j] = 0.;

		salt_rms[j] = 0.;
		salt_rmsg[j] = 0.;
	}

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				salt_ave[j] = salt_ave[j] + salt[k][j][i];
				salt_rms[j] = salt_rms[j] + salt[k][j][i]*salt[k][j][i];
			}
		}
	}

	ierr = MPI_Reduce (salt_ave, salt_aveg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (salt_rms, salt_rmsg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	if (params->rank == 0) {

		for (j=0;j<NY;j++){
			salt_aveg[j] = salt_aveg[j]/(NX-1.0)/(NZ-1.0);
			salt_rmsg[j] = salt_rmsg[j]/(NX-1.0)/(NZ-1.0);
		}

		double cbulk;
		cbulk = 0.0;

		for (j=0; j < NY-1; j++){
			cbulk = cbulk + salt_aveg[j]*(grid->yv[j+1]-grid->yv[j]);
		}
		cbulk = cbulk/(grid->yv[NY-1]-grid->yv[0]);

		fid = fopen("tsim_cbulk.dat","a");
		fprintf(fid, "%d %24.16e \n",flag, cbulk);
		fclose(fid);
	}

	return;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void xzperiodic_thermal_ave(Cart3d_bag *data_bag, Concentration **c, MAC_grid *grid,
		Parameters *params){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fid;

	double ***salt;

	double *salt_ave;
	double *salt_aveg;
	double *salt_rms;
	double *salt_rmsg;

	char bin_filename[50];
	double tetime;
	double zflow,xflow;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	salt = c[0]->data;

	salt_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	salt_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	salt_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	salt_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		salt_ave[j] = 0.;
		salt_aveg[j] = 0.;

		salt_rms[j] = 0.;
		salt_rmsg[j] = 0.;
	}



	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				salt_ave[j] = salt_ave[j] + salt[k][j][i];
				salt_rms[j] = salt_rms[j] + salt[k][j][i]*salt[k][j][i];
			}
		}
	}

	ierr = MPI_Reduce (salt_ave, salt_aveg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (salt_rms, salt_rmsg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	if (params->rank == 0) {

		tetime=params->time/params->output_time_interval + 4e-1;
		sprintf(bin_filename,"concave_t%d.dat",(int)tetime);
		fid = fopen(bin_filename,"w");

		for (j=0;j<NY;j++){
			salt_aveg[j] = salt_aveg[j]/(NX-1.0)/(NZ-1.0);
			salt_rmsg[j] = salt_rmsg[j]/(NX-1.0)/(NZ-1.0);
			salt_rmsg[j] = salt_rmsg[j] - salt_aveg[j]*salt_aveg[j];
			fprintf(fid,"%e %e %e \n",grid->yc[j],salt_aveg[j],salt_rmsg[j]);
		}

		fclose(fid);
	}

	return;


}




/******************************************************************************/
/*
 */
/******************************************************************************/
void xzperiodic_nut_ave(Cart3d_bag *data_bag, Subgrid *smag, MAC_grid *grid, Parameters *params){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fid;

	double ***cev, ***nu;
	double ***sc;

	double *cev_ave, *nu_ave;
	double *cev_aveg, *nu_aveg;
	double *cev_rms, *nu_rms;
	double *cev_rmsg, *nu_rmsg;
	double *sc_ave, *sc_rms;
	double *sc_aveg, *sc_rmsg;

	char bin_filename[50];
	double tetime;
	double zflow,xflow;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	cev = smag->ng_Cev;
	nu = smag->nut;
#ifdef CONC_DYNAMIC
	sc = smag->cdev[0]->Sct;
#endif

	cev_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	cev_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	cev_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	cev_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	sc_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		cev_ave[j] = 0.;
		nu_ave[j] = 0.;
		cev_aveg[j] = 0.;
		nu_aveg[j] = 0.;

		cev_rms[j] = 0.;
		nu_rms[j] = 0.;
		cev_rmsg[j] = 0.;
		nu_rmsg[j] = 0.;

		sc_ave[j] = 0.0;
		sc_rms[j] = 0.0;
		sc_aveg[j] = 0.0;
		sc_rmsg[j] = 0.0;
	}

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				cev_ave[j] = cev_ave[j] + cev[k][j][i];
				nu_ave[j] = nu_ave[j] + nu[k][j][i];
				cev_rms[j] = cev_rms[j] + cev[k][j][i]*cev[k][j][i];
				nu_rms[j] = nu_rms[j] + nu[k][j][i]*nu[k][j][i];
#ifdef CONC_DYNAMIC
				sc_ave[j] = sc_ave[j] + sc[k][j][i];
				sc_rms[j] = sc_rms[j] + sc[k][j][i]*sc[k][j][i];
#endif
			}
		}
	}
	ierr = MPI_Reduce (cev_ave, cev_aveg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (nu_ave,  nu_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (cev_rms, cev_rmsg, NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (nu_rms,  nu_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_ave,  sc_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_rms,  sc_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);


	if (params->rank == 0) {

		tetime=params->time/params->output_time_interval + 4e-1;
		sprintf(bin_filename,"nuave_t%d.dat",(int)tetime);
		fid = fopen(bin_filename,"w");

		for (j=0;j<NY;j++){
			cev_aveg[j] = cev_aveg[j]/(NX-1.0)/(NZ-1.0);
			nu_aveg[j] = nu_aveg[j]/(NX-1.0)/(NZ-1.0);

			cev_rmsg[j] = cev_rmsg[j]/(NX-1.0)/(NZ-1.0);
			nu_rmsg[j] = nu_rmsg[j]/(NX-1.0)/(NZ-1.0);

			cev_rmsg[j] = cev_rmsg[j] - cev_aveg[j]*cev_aveg[j];
			nu_rmsg[j] = nu_rmsg[j] - nu_aveg[j]*nu_aveg[j];

			sc_aveg[j] = sc_aveg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j] - sc_aveg[j]*sc_aveg[j];

			fprintf(fid,"%e %e %e %e %e %e %e\n",
					grid->yc[j],cev_aveg[j],nu_aveg[j],cev_rmsg[j],nu_rmsg[j],
					sc_aveg[j],sc_rmsg[j]);
		}

		fclose(fid);
	}

	return;
}


/******************************************************************************/
/*
 */
/******************************************************************************/
void xzperiodic1_rans_nut_ave(Cart3d_bag *data_bag, Rans *rans, MAC_grid *grid, Parameters *params){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fid;

	double ***nu;
	double ***sc;

	double *nu_ave;
	double *nu_aveg;
	double *nu_rms;
	double *nu_rmsg;
	double *sc_ave, *sc_rms;
	double *sc_aveg, *sc_rmsg;

	char bin_filename[50];
	double tetime;
	double zflow,xflow;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	nu = rans->nut;
#ifdef CONC_DYNAMIC
	sc = smag->cdev[0]->Sct;
#endif

	nu_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	sc_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		nu_ave[j] = 0.;
		nu_aveg[j] = 0.;
		nu_rms[j] = 0.;
		nu_rmsg[j] = 0.;

		sc_ave[j] = 0.0;
		sc_rms[j] = 0.0;
		sc_aveg[j] = 0.0;
		sc_rmsg[j] = 0.0;
	}

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				nu_ave[j] = nu_ave[j] + nu[k][j][i];
				nu_rms[j] = nu_rms[j] + nu[k][j][i]*nu[k][j][i];
#ifdef CONC_DYNAMIC
				sc_ave[j] = sc_ave[j] + sc[k][j][i];
				sc_rms[j] = sc_rms[j] + sc[k][j][i]*sc[k][j][i];
#endif
			}
		}
	}
	ierr = MPI_Reduce (nu_ave,  nu_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (nu_rms,  nu_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_ave,  sc_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_rms,  sc_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);


	if (params->rank == 0) {

		tetime=params->time/params->output_time_interval + 4e-1;
		sprintf(bin_filename,"nuave_t%d.dat",(int)tetime);
		fid = fopen(bin_filename,"w");

		for (j=0;j<NY;j++){
			nu_aveg[j] = nu_aveg[j]/(NX-1.0)/(NZ-1.0);
			nu_rmsg[j] = nu_rmsg[j]/(NX-1.0)/(NZ-1.0);

			nu_rmsg[j] = nu_rmsg[j] - nu_aveg[j]*nu_aveg[j];

/*			sc_aveg[j] = sc_aveg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j] - sc_aveg[j]*sc_aveg[j];
			fprintf(fid,"%e %e %e %e %e %e %e\n",
					grid->yc[j],cev_aveg[j],nu_aveg[j],cev_rmsg[j],nu_rmsg[j],
					sc_aveg[j],sc_rmsg[j]);
*/

			fprintf(fid,"%e %e %e \n",
					grid->yc[j],nu_aveg[j],nu_rmsg[j]);
		}

		fclose(fid);

	}

	return;
}


/******************************************************************************/
/*
 */
/******************************************************************************/
void xzperiodic_rans_nut_ave(Cart3d_bag *data_bag, Rans *rans, MAC_grid *grid, Parameters *params){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fid;

	double ***nu;
	double ***sc;

	double *nu_ave;
	double *nu_aveg;
	double *nu_rms;
	double *nu_rmsg;
	double *sc_ave, *sc_rms;
	double *sc_aveg, *sc_rmsg;
	double *tke_ave, *tke_rms, *tke_aveg, *tke_rmsg;
	double *diss_ave, *diss_rms, *diss_aveg, *diss_rmsg;
	double ***tke, ***diss;

	char bin_filename[50];
	double tetime;
	double zflow,xflow;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	nu = rans->nut;
#ifdef CONC_DYNAMIC
	sc = smag->cdev[0]->Sct;
#endif
	tke = rans->two_eqn_rans[0]->data;
	diss = rans->two_eqn_rans[1]->data;

	nu_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	sc_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	tke_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	tke_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	tke_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	tke_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	diss_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	diss_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	diss_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	diss_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		nu_ave[j] = 0.;
		nu_aveg[j] = 0.;
		nu_rms[j] = 0.;
		nu_rmsg[j] = 0.;

		sc_ave[j] = 0.0;
		sc_rms[j] = 0.0;
		sc_aveg[j] = 0.0;
		sc_rmsg[j] = 0.0;

		tke_ave[j] = 0.0;
		tke_aveg[j] = 0.0;
		tke_rms[j] = 0.0;
		tke_rmsg[j] = 0.0;

		diss_ave[j] = 0.0;
		diss_aveg[j] = 0.0;
		diss_rms[j] = 0.0;
		diss_rmsg[j] = 0.0;
	}

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				nu_ave[j] = nu_ave[j] + nu[k][j][i];
				nu_rms[j] = nu_rms[j] + nu[k][j][i]*nu[k][j][i];

				tke_ave[j] = tke_ave[j] + tke[k][j][i];
				tke_rms[j] = tke_rms[j] + tke[k][j][i]*tke[k][j][i];

				diss_ave[j] = diss_ave[j] + diss[k][j][i];
				diss_rms[j] = diss_rms[j] + diss[k][j][i]*diss[k][j][i];
#ifdef CONC_DYNAMIC
				sc_ave[j] = sc_ave[j] + sc[k][j][i];
				sc_rms[j] = sc_rms[j] + sc[k][j][i]*sc[k][j][i];
#endif
			}
		}
	}
	ierr = MPI_Reduce (nu_ave,  nu_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (nu_rms,  nu_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_ave,  sc_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_rms,  sc_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	ierr = MPI_Reduce (tke_ave,  tke_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (tke_rms,  tke_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (diss_ave,  diss_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (diss_rms,  diss_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	if (params->rank == 0) {

		tetime=params->time/params->output_time_interval + 4e-1;
		sprintf(bin_filename,"nuave_t%d.dat",(int)tetime);
		fid = fopen(bin_filename,"w");

		for (j=0;j<NY;j++){
			nu_aveg[j] = nu_aveg[j]/(NX-1.0)/(NZ-1.0);
			nu_rmsg[j] = nu_rmsg[j]/(NX-1.0)/(NZ-1.0);

			nu_rmsg[j] = nu_rmsg[j] - nu_aveg[j]*nu_aveg[j];

			tke_aveg[j] = tke_aveg[j]/(NX-1.0)/(NZ-1.0);
			tke_rmsg[j] = tke_rmsg[j]/(NX-1.0)/(NZ-1.0);
			diss_aveg[j] = diss_aveg[j]/(NX-1.0)/(NZ-1.0);
			diss_rmsg[j] = diss_rmsg[j]/(NX-1.0)/(NZ-1.0);

			tke_rmsg[j] = tke_rmsg[j] - tke_aveg[j]*tke_aveg[j];
			diss_rmsg[j] = diss_rmsg[j] - diss_aveg[j]*diss_aveg[j];

/*			sc_aveg[j] = sc_aveg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j] - sc_aveg[j]*sc_aveg[j];
			fprintf(fid,"%e %e %e %e %e %e %e\n",
					grid->yc[j],cev_aveg[j],nu_aveg[j],cev_rmsg[j],nu_rmsg[j],
					sc_aveg[j],sc_rmsg[j]);
*/

			fprintf(fid,"%e %e %e %e %e %e %e\n",
					grid->yc[j],nu_aveg[j],nu_rmsg[j],tke_aveg[j],tke_rmsg[j],diss_aveg[j],diss_rmsg[j]);
		}

		fclose(fid);

	}

	return;
}


/******************************************************************************/
/*
 */
/******************************************************************************/
void rans_nut_ave(Rans *rans, MAC_grid *grid, Parameters *params, int iter, int rk){

	int NX, NY, NZ;
	int i, j, k, index;
	int ierr ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	FILE *fid;

	double ***nu;
	double ***sc;

	double *nu_ave;
	double *nu_aveg;
	double *nu_rms;
	double *nu_rmsg;
	double *sc_ave, *sc_rms;
	double *sc_aveg, *sc_rmsg;
	double *tke_ave, *tke_rms, *tke_aveg, *tke_rmsg;
	double *diss_ave, *diss_rms, *diss_aveg, *diss_rmsg;
	double ***tke, ***diss;

	char bin_filename[50];
	double tetime;
	double zflow,xflow;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	nu = rans->nut;
#ifdef CONC_DYNAMIC
	sc = smag->cdev[0]->Sct;
#endif
	tke = rans->two_eqn_rans[0]->data;
	diss = rans->two_eqn_rans[1]->data;

	nu_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	nu_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	sc_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	sc_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	tke_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	tke_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	tke_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	tke_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	diss_ave = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	diss_aveg = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	diss_rms = Memory_allocate_1D_array(GVG_DOUBLE,NY);
	diss_rmsg = Memory_allocate_1D_array(GVG_DOUBLE,NY);

	for (j=0;j<NY;j++){
		nu_ave[j] = 0.;
		nu_aveg[j] = 0.;
		nu_rms[j] = 0.;
		nu_rmsg[j] = 0.;

		sc_ave[j] = 0.0;
		sc_rms[j] = 0.0;
		sc_aveg[j] = 0.0;
		sc_rmsg[j] = 0.0;

		tke_ave[j] = 0.0;
		tke_aveg[j] = 0.0;
		tke_rms[j] = 0.0;
		tke_rmsg[j] = 0.0;

		diss_ave[j] = 0.0;
		diss_aveg[j] = 0.0;
		diss_rms[j] = 0.0;
		diss_rmsg[j] = 0.0;
	}

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Ke = min(Ke, NZ-1);

	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				nu_ave[j] = nu_ave[j] + nu[k][j][i];
				nu_rms[j] = nu_rms[j] + nu[k][j][i]*nu[k][j][i];

				tke_ave[j] = tke_ave[j] + tke[k][j][i];
				tke_rms[j] = tke_rms[j] + tke[k][j][i]*tke[k][j][i];

				diss_ave[j] = diss_ave[j] + diss[k][j][i];
				diss_rms[j] = diss_rms[j] + diss[k][j][i]*diss[k][j][i];
#ifdef CONC_DYNAMIC
				sc_ave[j] = sc_ave[j] + sc[k][j][i];
				sc_rms[j] = sc_rms[j] + sc[k][j][i]*sc[k][j][i];
#endif
			}
		}
	}
	ierr = MPI_Reduce (nu_ave,  nu_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (nu_rms,  nu_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_ave,  sc_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (sc_rms,  sc_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	ierr = MPI_Reduce (tke_ave,  tke_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (tke_rms,  tke_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (diss_ave,  diss_aveg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);
	ierr = MPI_Reduce (diss_rms,  diss_rmsg,  NY, MPI_DOUBLE, MPI_SUM, 0, PCW);

	if (params->rank == 0) {

		sprintf(bin_filename,"nutave_t%d.dat",iter*10+rk);
		fid = fopen(bin_filename,"w");

		for (j=0;j<NY;j++){
			nu_aveg[j] = nu_aveg[j]/(NX-1.0)/(NZ-1.0);
			nu_rmsg[j] = nu_rmsg[j]/(NX-1.0)/(NZ-1.0);

			nu_rmsg[j] = nu_rmsg[j] - nu_aveg[j]*nu_aveg[j];

			tke_aveg[j] = tke_aveg[j]/(NX-1.0)/(NZ-1.0);
			tke_rmsg[j] = tke_rmsg[j]/(NX-1.0)/(NZ-1.0);
			diss_aveg[j] = diss_aveg[j]/(NX-1.0)/(NZ-1.0);
			diss_rmsg[j] = diss_rmsg[j]/(NX-1.0)/(NZ-1.0);

			tke_rmsg[j] = tke_rmsg[j] - tke_aveg[j]*tke_aveg[j];
			diss_rmsg[j] = diss_rmsg[j] - diss_aveg[j]*diss_aveg[j];

/*			sc_aveg[j] = sc_aveg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j]/(NX-1.0)/(NZ-1.0);
			sc_rmsg[j] = sc_rmsg[j] - sc_aveg[j]*sc_aveg[j];
			fprintf(fid,"%e %e %e %e %e %e %e\n",
					grid->yc[j],cev_aveg[j],nu_aveg[j],cev_rmsg[j],nu_rmsg[j],
					sc_aveg[j],sc_rmsg[j]);
*/

			fprintf(fid,"%e %e %e %e %e %e %e\n",
					grid->yc[j],nu_aveg[j],nu_rmsg[j],tke_aveg[j],tke_rmsg[j],diss_aveg[j],diss_rmsg[j]);
		}

		fclose(fid);

	}

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Communication_finalize() {

	int ierr = MPI_Finalize();
}
