#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "hdf5.h"

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Display.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "MergeSort.h"
#include "Output.h"
#include "Particle.h"
#include "ParticleOutput.h"


/******************************************************************************/
/*
 Writes particle information to "mobile.dat"
 */
/******************************************************************************/
void ParticleOutput_dat(Particle_list *p_list, MAC_grid *grid, Parameters *params,
		Debug_trace *dtrace) {

	// Only write this file if there is one mobile particle
	if (p_list->Np > 2)
		return;

	// Check state of linked lists
	Display_assert_list_state(p_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	int i;
	MPI_Status status;
	int nproc = params -> size;
	int rank = params -> rank;
	int *Np, Np_local, *Nc;
	Particle *p, *p_last, *send_part, *recv_part;
	FILE *fptr;

	// Count particles on local processor
	Np_local = 0;
	p = p_list->start;
	p_last = NULL;
	while (p != NULL) {
		Np_local++;
		p_last = p;
		p = p -> next;
	}

	// Broadcast local number of particles to receiving processor (0)
	Np = (int *)malloc(nproc * sizeof(int));
	MPI_Gather(&Np_local, 1, MPI_INT, Np, 1, MPI_INT, 0, PCW);

	// Receive particle data at processor 0
	if (rank == 0) {

		for (i = 1; i < nproc; i++) {
			recv_part = (Particle *)malloc(Np[i] * sizeof(Particle));
			Memory_check_allocation(recv_part);

			MPI_Recv(recv_part, Np[i], MPI_PARTICLE, i, 1, PCW, &status);

			Nc = (int *)malloc((2 * Np[i] + 1) * sizeof(int));
			Memory_check_allocation(Nc);
			memset(Nc, 0, (2 * Np[i] + 1) * sizeof(int));

			Particle_list_add_array(p_list, recv_part, Np[i], NULL, Nc, grid);

			free(recv_part);
			free(Nc);
		}
	}
	// Send particle data to processor 0
	else {
		send_part = (Particle *)malloc(Np_local * sizeof(Particle));
		Memory_check_allocation(send_part);
		p = p_list->start;
		for (i = 0; i < Np_local; i++) {
			send_part[i] = *p;
			p = p -> next;
		}
		MPI_Send(send_part, Np_local, MPI_PARTICLE, 0, 1, PCW);
		free(send_part);
	}

	// Print particle data
	if (rank == 0) {
		if (params -> time == 0)
			fptr = fopen("mobile.dat", "w");
		else
			fptr = fopen("mobile.dat", "a");
		p = p_list->start;
		while (p != NULL) {

				fprintf(fptr,
						"% .10g,%d,"
// 						"% .10g,% .10g,% .10g,"
// 						"% .10g,% .10g,% .10g,"
// 						"% .10g,% .10g,% .10g,"
// 						"% .10g,% .10g,% .10g,"
						"%.10g\n",
						params->time, p->ID,
// 						p -> X[0],      p -> X[1],      p -> X[2],
// 						p -> U[0],      p -> U[1],      p -> U[2],
// 						p -> Omega[0],  p -> Omega[1],  p -> Omega[2],
// 						p -> Fc[0],     p -> Fc[1],     p -> Fc[2],
						p -> X[1] );
			p = p -> next;
		}
		fclose(fptr);
	}

	free(Np);

	p_list->state = LIST_STATE_BOTH;
	Particle_list_remove(p_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));
}




/******************************************************************************/
/*
 * Save particle data to Particle_XX.h5
 * where XX is the input 'counter'.  If 'counter == params->noutput', the output
 * will be in the standard directory with Data_XX.h5.  Otherwise, the output
 * will be written to a subfolder
 *
 * This file contains all the particle information necessary to restart the
 * simulation, as well as the particle information wanted for visualization and
 * post-processing.
 */
/******************************************************************************/
void ParticleOutput_h5(Cart3d_bag *data_bag, int counter, Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	int verbose = 1;

	//--------------------------------------------------------------------------
	// HDF5 Handles
	//--------------------------------------------------------------------------
	hid_t   file_id;   // File handles
	hsize_t dim_1d[1]; // Dimensions of data to be written
	herr_t  status;    // Status variable to check return value of HDF5 routines

	//--------------------------------------------------------------------------
	// Strings
	//--------------------------------------------------------------------------
	char h5filename[50];  // HDF5 filename to store 3D data
	char groupname[100];  // HDF5 group name to group data
	char fieldname[100];  // HDF5 name
	char message[500];    // Message for errors/warnings

	//--------------------------------------------------------------------------
	// Databag parameters
	//--------------------------------------------------------------------------
	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	// Particle data
	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	//--------------------------------------------------------------------------
	// Open HDF5 file handles
	//--------------------------------------------------------------------------
	if (counter == params->noutput){
		sprintf(h5filename, "Particle_%d.h5", abs(counter));
	}
	else {
		sprintf(h5filename, "./trn/Particle_%d.h5", counter);
	}
	file_id = Output_h5_open(h5filename, params, DTRACE("Output_h5_open"));

	//--------------------------------------------------------------------------
	// Write timestamp
	//--------------------------------------------------------------------------
	if (verbose) Display_progress(params,"ParticleIO.c: write time\n");
	sprintf(fieldname, "/time");
	dim_1d[0] = 1;
	Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	//--------------------------------------------------------------------------
	// Write domain geometry
	//--------------------------------------------------------------------------
	sprintf(groupname, "/domain");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	dim_1d[0] = 1;
	sprintf(fieldname, "%s/xmin", groupname);
	Output_h5_dataset(&(params->xmin), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/xmax", groupname);
	Output_h5_dataset(&(params->xmax), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/ymin", groupname);
	Output_h5_dataset(&(params->ymin), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/ymax", groupname);
	Output_h5_dataset(&(params->ymax), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/zmin", groupname);
	Output_h5_dataset(&(params->zmin), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/zmax", groupname);
	Output_h5_dataset(&(params->zmax), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	int periodic[] = {0, 0, 0};
#ifdef XPERIODIC
	periodic[0] = 1;
#endif
#ifdef YPERIODIC
	periodic[1] = 1;
#endif
#ifdef ZPERIODIC
	periodic[2] = 1;
#endif
	dim_1d[0] = 3;
	sprintf(fieldname, "%s/periodic", groupname);
	Output_h5_dataset(periodic, H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	//--------------------------------------------------------------------------
	// Write Lagrangian points (not HDF5, debugging)
	//--------------------------------------------------------------------------
	if (0){ // TODO: Debugging
		int i;
		char filename[50];

		Particle *p = p_mobile_list->start;
		if (p != NULL && params->time != 0.0) {

			sprintf(filename, "Lag_pts0_%d.dat", abs(counter));
			FILE *fid0 = fopen(filename, "w");
			sprintf(filename, "Lag_pts1_%d.dat", abs(counter));
			FILE *fid1 = fopen(filename, "w");

			double *X_L = p -> X_L;
			double *Y_L = p -> Y_L;
			double *Z_L = p -> Z_L;
			int *flag_L = p -> flag_L;

			fprintf(fid0, "%.12g, %.12g, %.12g\n", p->X[0], p->X[1], p->X[2]);
			fprintf(fid1, "%.12g, %.12g, %.12g\n", p->X[0], p->X[1], p->X[2]);

			for (i = 0; i < p->N_L_local; i++) {
				if (flag_L[i] == 0) {
					fprintf(fid0, "%.12g, %.12g, %.12g\n", X_L[i], Y_L[i], Z_L[i]);
				}
				else {
					fprintf(fid1, "%.12g, %.12g, %.12g\n", X_L[i], Y_L[i], Z_L[i]);
				}
			}
			fclose(fid0);
			fclose(fid1);
		}
	}

	//--------------------------------------------------------------------------
	// Write mobile particles
	//--------------------------------------------------------------------------
	if (p_mobile_list->Np > 0) {

		if (verbose) Display_progress(params,"ParticleIO.c: write H5 mobile particles\n");
		sprintf(groupname, "/mobile");
		Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));
		ParticleOutput_h5_data(p_mobile_list, file_id, groupname, params, DTRACE("ParticleOutput_h5_data"));
	}

	//--------------------------------------------------------------------------
	// Write fixed particles
	//--------------------------------------------------------------------------
	if (p_fixed_list->Np > 0) {

		if (verbose) Display_progress(params,"ParticleIO.c: write H5 fixed particles\n");
		sprintf(groupname, "/fixed");
		Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));
		ParticleOutput_h5_data(p_fixed_list, file_id, groupname, params, DTRACE("ParticleOutput_h5_data"));
	}

	//--------------------------------------------------------------------------
	// Close HDF5 file handles
	//--------------------------------------------------------------------------
	assert(H5Fclose(file_id) >= 0);

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_output += T2 - T1;

	return;
}




/******************************************************************************/
/*
 * Writes a specific particle linked list to HDF5 file 'file_id'
 */
/******************************************************************************/
void ParticleOutput_h5_data(Particle_list *p_list, hid_t file_id,
		char *groupname, Parameters *params, Debug_trace *dtrace) {


	Particle *p;
	Collision *pc;
	int Np_local = 0;  // Number of local particles
	int Nc, Nc_max, Nc_max_local = 0;  // Maximum number of collisions for one particle

#ifndef DIRTY_FIX
	//printf("A !!!!!! Here still everything OK for id %d!!!!! \n", params->rank);
	//		fflush(stdout);
	// Sort list by particle identity.  This is important for HDF5 hyperslabs!
	MergeSort(&(p_list->start));
	//printf("B !!!!!! Here still everything OK for id %d!!!!! \n", params->rank);
	//		fflush(stdout);

#endif


	// Count number of particles on local processor
	p = p_list -> start;
	while (p != NULL) {
		Np_local++;

		// Count number of particle and wall collisions for this particle
		Nc = 0;
		pc = p -> particle_collision;
		while (pc != NULL) {
			Nc++;
			pc = pc -> next;
		}
		pc = p -> wall_collision;
		while (pc != NULL) {
			Nc++;
			pc = pc -> next;
		}

		// Find maximum number of collisions per particle on local processor
		Nc_max_local = max(Nc_max_local, Nc);

		p = p -> next;
	}

#define OUTPUT_ELEMENT(str, dataSize) \
    ParticleOutput_h5_data_element(p_list, Np_local, file_id, groupname, str, dataSize, \
                                   params, DTRACE("ParticleOutput_h5_data_element"))
	OUTPUT_ELEMENT("R", 1);
	OUTPUT_ELEMENT("X", 3);
	//OUTPUT_ELEMENT("X_old", 3);
	OUTPUT_ELEMENT("Rotn", 9);
	OUTPUT_ELEMENT("U", 3);
	//OUTPUT_ELEMENT("U_old", 3);
	OUTPUT_ELEMENT("Omega", 3);
//	OUTPUT_ELEMENT("F", 3);
//	OUTPUT_ELEMENT("T", 3);
	OUTPUT_ELEMENT("Fc", 3);
	OUTPUT_ELEMENT("Tc", 3);
	OUTPUT_ELEMENT("F_IBM", 3);
	OUTPUT_ELEMENT("T_IBM", 3);
	OUTPUT_ELEMENT("F_rigid", 3);
	OUTPUT_ELEMENT("T_rigid", 3);
	OUTPUT_ELEMENT("F_coll", 3);
	OUTPUT_ELEMENT("T_coll", 3);
	//OUTPUT_ELEMENT("Int_U_old", 3);
	OUTPUT_ELEMENT("Int_Omega_old", 3);

#ifdef POST_PROCESS
	OUTPUT_ELEMENT("Fc_norm_cum", 3);
	OUTPUT_ELEMENT("Fc_tan_cum", 3);
	OUTPUT_ELEMENT("Fl_norm_cum", 3);
	OUTPUT_ELEMENT("Fl_tan_cum", 3);
#endif

#undef OUTPUT_ELEMENT

	// Communicate maximum number collisions per particle
	MPI_Allreduce(&Nc_max_local, &Nc_max, 1, MPI_INT, MPI_MAX, PCW);

	// Write collisions
	ParticleOutput_h5_collision(p_list, Np_local, Nc_max, file_id, groupname, params, DTRACE("ParticleOutput_h5_collision"));
}




/******************************************************************************/
/*
 * Writes particle information to Particle_*.h5 file
 *
 * The data to be written is determined by the input variable 'element', which
 * stores 'dataSize' doubles per particle.  If more elements are added to the
 * Particle structure, they need to be added in this function below in order to
 * support output
 */
/******************************************************************************/
void ParticleOutput_h5_data_element(Particle_list *p_list, int Np_local,
		hid_t file_id, char *groupname, char *element, int dataSize,
		Parameters *params, Debug_trace *dtrace) {

	Particle *p;
	Collision *pc;
	int i, j;
	int Nc_max, Nc_max_local = 1;  // Maximum number of collisions for one particle

	double **data;     // Buffer for particle data to store
	int ndim     = 2;  // Dimension of array to write

	hsize_t dim[ndim], count[ndim], offset[ndim];
	hid_t   plist_id, dataset, dataspace, memspace;
	herr_t  status;

	dim[0] = p_list->Np; // Number of particles globally
	dim[1] = dataSize;   // Number of doubles to store for each particle

	char message[100], filename[50], coll_fieldname[50];

	char fieldname[100];
	sprintf(fieldname, "%s/%s", groupname, element);

	//--------------------------------------------------------------------------
	// Define dataspace, represents required space in file to store particle
	// information from all processors
	//--------------------------------------------------------------------------
	dataspace = H5Screate_simple(ndim, dim, NULL);
	assert(dataspace >= 0);

	// Create dataset using defined dataspace and fieldname
	dataset = H5Dcreate(file_id, fieldname, H5T_NATIVE_DOUBLE, dataspace,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  // Invalid read
	assert(dataset >= 0);

	// Create property list for collective dataset write
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	assert(plist_id >= 0);
	assert(H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) >= 0);

	//--------------------------------------------------------------------------
	// Define memory space, represents required space to transfer only local
	// particle information
	//--------------------------------------------------------------------------
	count[0] = Np_local;
	count[1] = dataSize;
	memspace = H5Screate_simple(ndim, count, NULL);
	assert(memspace >= 0);

	// Select empty set of spaces if there are no particles to write
	if (Np_local == 0) {
		status = H5Sselect_none(dataspace);
		status = H5Sselect_none(memspace);
	}

	// Allocate space for data buffer (same size as 'memspace')
	data = Memory_allocate_2D_double_array(dataSize, max(Np_local, 1));

	// Hyperslab dimensions and offset, holds information for one particle
	count[0] = 1;
	count[1] = dataSize;
	offset[1] = 0;

	//--------------------------------------------------------------------------
	// Select file hyperslab associated with particle ID (output).  Here we
	// addend a slab of particle data to the hyperslab as we find each local
	// particle (thus the SELECT_SET followed by successive SELECT_OR)
	//--------------------------------------------------------------------------
	p = p_list -> start;
	for (j = 0; j < Np_local; j++) {

		offset[0] = p -> ID - p_list->ID_start;
		if (j == 0)
			status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
		else
			status = H5Sselect_hyperslab(dataspace, H5S_SELECT_OR, offset, NULL, count, NULL);

		p = p -> next;
	}  // for j

	//--------------------------------------------------------------------------
	// Fill data buffer with local particle information
	//--------------------------------------------------------------------------
	p = p_list -> start;
	if (strcmp(element, "R") == 0) {
		for (j = 0; j < Np_local; j++) {
			data[j][0] = p -> R;
			p = p -> next;
		}
	}
	else if (strcmp(element, "X") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> X[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Rotn") == 0) {
		for (j = 0; j < Np_local; j++) {
			for (i = 0; i < 9; i++) {
				data[j][i] = p -> Rotn[i/3][i%3];
			}
			p = p -> next;
		}
	}
	else if (strcmp(element, "U") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> U[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Omega") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Omega[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "F") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> F[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Fc") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Fc[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "T") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> T[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Tc") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Tc[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Int_U_old") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Int_U_old[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Int_Omega_old") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Int_Omega_old[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "F_IBM") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> F_IBM[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "T_IBM") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> T_IBM[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "F_rigid") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> F_rigid[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "T_rigid") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> T_rigid[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "F_coll") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> F_coll[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "T_coll") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> T_coll[i];
			p = p -> next;
		}
	}

#ifdef POST_PROCESS
	else if (strcmp(element, "Fc_norm_cum") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Fc_norm_cum[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Fc_tan_cum") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Fc_tan_cum[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Fl_norm_cum") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Fl_norm_cum[i];
			p = p -> next;
		}
	}
	else if (strcmp(element, "Fl_tan_cum") == 0) {
		for (j = 0; j < Np_local; j++) {
			FORI3 data[j][i] = p -> Fl_tan_cum[i];
			p = p -> next;
		}
	}
#endif

	else {
		sprintf(message, "Invalid element \"%s\"", element);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	//--------------------------------------------------------------------------
	// Write data to *.h5 file
	//--------------------------------------------------------------------------
	assert(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, data[0]) >= 0);

	// Close/release resources
	free(data);
	assert(H5Dclose(dataset) >= 0);
	assert(H5Sclose(dataspace) >= 0);
	assert(H5Sclose(memspace) >= 0);
	assert(H5Pclose(plist_id) >= 0);
}




/******************************************************************************/
/*
 * Writes particle collision information to Particle_*.h5 file
 *
 * Collision information is stored as a 2-D array:
 *
 *     data[k][j]
 *
 *     - k      -> particle number
 *     - j == 0 -> number of active particle collisions
 *     - j == 1 -> number of active wall collisions
 *     - j >= 2 -> blocks of length 'dataSize', each containing all the active
 *                 collision information
 *
 * If adding or removing elements from Collision structure:
 *   1. Change 'dataSize' below
 *   2. Add/remove elements from while loop, e.g.:
 *          data[k][j++] = pc -> dn;
 *
 * This function should mirror the changes made to
 * ParticleInput_copy_collision_data()
 */
/******************************************************************************/
void ParticleOutput_h5_collision(Particle_list *p_list, int Np_local, int Nc_max,
		hid_t file_id, char *groupname, Parameters *params, Debug_trace *dtrace) {

	Particle *p;
	Collision *pc;
	int j, k, stage;

	double **data;   // Buffer for particle data to store
	int dataSize = 10; // Number of elements per collision to store
	int ndim     = 2;  // Dimension of array to write

	hsize_t dim[ndim], count[ndim], offset[ndim];
	hid_t   plist_id, dataset, dataspace, memspace;
	herr_t  status;

	dim[0] = p_list->Np; // Number of particles globally
	dim[1] = 2 + Nc_max * dataSize;  // Maximum number of collisions owned by any particle

	char message[100], filename[50];

	char fieldname[100];
	sprintf(fieldname, "%s/collision_data", groupname);

	//--------------------------------------------------------------------------
	// Define dataspace, represents required space in file to store particle
	// information from all processors
	//--------------------------------------------------------------------------
	dataspace = H5Screate_simple(ndim, dim, NULL);
	assert(dataspace >= 0);

	// Create dataset using defined dataspace and fieldname
	dataset = H5Dcreate(file_id, fieldname, H5T_NATIVE_DOUBLE, dataspace,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  // Invalid read
	assert(dataset >= 0);

	// Create property list for collective dataset write
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	assert(plist_id >= 0);
	assert(H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) >= 0);

	//--------------------------------------------------------------------------
	// Define memory space, represents required space to transfer only local
	// particle information
	//--------------------------------------------------------------------------
	count[0] = Np_local;
	count[1] = dim[1];
	memspace = H5Screate_simple(ndim, count, NULL);
	assert(memspace >= 0);

	// Select empty set of spaces if there are no particles to write
	if (Np_local == 0) {
		status = H5Sselect_none(dataspace);
		status = H5Sselect_none(memspace);
	}

	// Allocate space for data buffer (same size as 'memspace')
	data = Memory_allocate_2D_double_array(dim[1], max(Np_local, 1));

	// Hyperslab dimensions and offset, holds information for one particle
	count[0] = 1;
	count[1] = dim[1];
	offset[1] = 0;

	//--------------------------------------------------------------------------
	// Fill data buffer with local particle collision information
	//--------------------------------------------------------------------------
	p = p_list -> start;
	for (k = 0; k < Np_local; k++) {

		//----------------------------------------------------------------------
		// Select file hyperslab associated with particle ID (output).  Here we
		// addend a slab of particle data to the hyperslab as we find each
		// local particle (thus the SELECT_SET followed by successive SELECT_OR)
		//----------------------------------------------------------------------
		offset[0] = p -> ID - p_list->ID_start;
		if (k == 0)
			status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
		else
			status = H5Sselect_hyperslab(dataspace, H5S_SELECT_OR, offset, NULL, count, NULL);

		// Start placing values after Npc and Nwc
		j = 2;
		stage = 0;
		pc = p -> particle_collision;
		if (pc == NULL) {
			pc = p -> wall_collision;
			stage = 1;
		}
		while (pc != NULL) {

			data[k][j++] = pc -> kn;
			data[k][j++] = pc -> dn;
			data[k][j++] = pc -> kt;
			data[k][j++] = pc -> dt;
			data[k][j++] = pc -> St;
			data[k][j++] = pc -> zeta_t[0];
			data[k][j++] = pc -> zeta_t[1];
			data[k][j++] = pc -> zeta_t[2];
			data[k][j++] = (double)(pc -> other_ID);
			data[k][j++] = (double)(pc -> type);

			// Increment number of particle/wall collisions
			if (stage == 0)
				data[k][0] += 1.0;  // Particle collision
			else
				data[k][1] += 1.0;  // Wall collision

			// Move on to next collision, or to wall collisions if done with
			// particles
			pc = pc -> next;
			if (pc == NULL && stage == 0) {
				pc = p -> wall_collision;
				stage = 1;
			}
		}

		p = p -> next;
	}  // for k

	//--------------------------------------------------------------------------
	// Write data to *.h5 file
	//--------------------------------------------------------------------------
	assert(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, data[0]) >= 0);

	// Close/release resources
	free(data);
	assert(H5Dclose(dataset) >= 0);
	assert(H5Sclose(dataspace) >= 0);
	assert(H5Sclose(memspace) >= 0);
	assert(H5Pclose(plist_id) >= 0);
}
