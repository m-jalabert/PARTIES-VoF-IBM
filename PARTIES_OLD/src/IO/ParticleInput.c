#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf5.h"

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Display.h"
#include "Particle.h"
#include "ParticleInput.h"
#include "Lagrangian.h"
#include "Memory.h"


/******************************************************************************/
/*
 * Read particles from input files P_MOBILE_INPUT_FILE and P_FIXED_INPUT_FILE
 * (macros defined in 'definitions.h')
 */
/******************************************************************************/
void ParticleInput_inp(Particle_list *p_list, MAC_grid *grid, Parameters *params,
		Debug_trace *dtrace) {

	int i;

	// Processor rank
	int rank = params -> rank;

	int Np;        // Number of particles on all processors
	int np;        // Number of particles on this processor
	double **X;    // Locations of particles
	double *R;     // Radii of particles
	double *temp;  // Storage for location/radius arrays
	double *ini_rot_theta; // Initial paricle rotation with respect to z-axis in x-y plane
	double *ini_rot_phi; // Initial paricle rotation with respect to y-axis in z-x plane
	double *t_part_release; // Time at which the particle is released

	// Particle linked lists
	Particle *p, *p_last;

	// Grid bounds of local subdomain
	double G_xmin = grid -> xu[grid -> G_Is];
	double G_xmax = grid -> xu[min(grid -> G_Ie, grid -> NX - 1)];
	double G_ymin = grid -> yv[grid -> G_Js];
	double G_ymax = grid -> yv[min(grid -> G_Je, grid -> NY - 1)];
	double G_zmin = grid -> zw[grid -> G_Ks];
	double G_zmax = grid -> zw[min(grid -> G_Ke, grid -> NZ - 1)];

	// File pointer for particle input file
	FILE *fptr;
	char filename[50];
	int count;

	// Variables for error statements
	char message[500] = "";
	int ierr = 0;

	//--------------------------------------------------------------------------
	// Open input file
	//--------------------------------------------------------------------------
	if (rank == 0) {
		if (p_list->type == MOBILE)
			sprintf(filename, P_MOBILE_INPUT_FILE);
		else if (p_list->type == FIXED)
			sprintf(filename, P_FIXED_INPUT_FILE);
#ifdef PARTICLE_RELEASE
		else if (p_list->type == RELEASE)
			sprintf(filename, P_RELEASE_INPUT_FILE);
#endif
		else {
			sprintf(message, "Incorrect input p_list->type = '%d'\n"
			                 "Use 'MOBILE' or 'FIXED'", p_list->type);
			ierr = -1;
		}

		if (ierr == 0) {
			fptr = fopen(filename, "r");
			if (fptr == NULL) {
				sprintf(message, "Could not open %s", filename);
				ierr = -2;
			}
		}
	}
	Display_assert_error(ierr, message, params, DTRACE("Display_assert_error"));

	//--------------------------------------------------------------------------
	// Read number of particles
	//--------------------------------------------------------------------------
	if (rank == 0) {
		count = fscanf(fptr, "%d\n", &Np);
		if (count < 0) {
			sprintf(message, "Incorrect formatting of %s\n"
			        "First line should specify number of particles", filename);
			ierr = -3;
		}
	}
	Display_assert_error(ierr, message, params, DTRACE("Display_assert_error"));

	// Broadcast number of particles to other processors
	MPI_Bcast(&Np, 1, MPI_INT, 0, PCW);
	p_list -> Np = Np;

int N_read_data = 4;

	// Assign storage
#ifdef SQUIRMER_SWIMMER
	N_read_data += 2;
#endif

#ifdef PARTICLE_RELEASE
	N_read_data += 1;
#endif

	temp = (double *)malloc(N_read_data * Np * sizeof(double));
	Memory_check_allocation(temp);

	// 'X' and 'R' point to position and radius information stored in 'temp'
	X = (double **)malloc(Np * sizeof(double *));
	Memory_check_allocation(X);
	for (i = 0; i < Np; i++) {
		X[i] = &temp[3 * i];
	}

	R = &temp[3 * Np];
#ifdef SQUIRMER_SWIMMER
	ini_rot_theta = &temp[4 * (Np)];
	ini_rot_phi = &temp[5 * (Np)];
	#ifdef PARTICLE_RELEASE
	t_part_release = &temp[6 * (Np)];
	#endif
#elif defined PARTICLE_RELEASE
	t_part_release = &temp[4 * (Np)];
#endif


	//-------------------------------------------------------------------
	// Read particle information
	//--------------------------------------------------------------------------
	if (rank == 0) {
		printf("Reading %s\n", filename);
		for (i = 0; i < Np; i++) {

#ifdef SQUIRMER_SWIMMER
	#ifdef PARTICLE_RELEASE
			count = fscanf(fptr, "%lf %lf %lf %lf %lf %lf %lf", &X[i][0], &X[i][1], &X[i][2], &R[i], &ini_rot_theta[i], &ini_rot_phi[i], &t_part_release[i]);
	#else
			count = fscanf(fptr, "%lf %lf %lf %lf %lf %lf", &X[i][0], &X[i][1], &X[i][2], &R[i], &ini_rot_theta[i], &ini_rot_phi[i]);
	#endif

#elif defined PARTICLE_RELEASE
			count = fscanf(fptr, "%lf %lf %lf %lf %lf", &X[i][0], &X[i][1], &X[i][2], &R[i], &t_part_release[i]);
#else
			count = fscanf(fptr, "%lf %lf %lf %lf", &X[i][0], &X[i][1], &X[i][2], &R[i]);
#endif




			if (count < 0) {
				sprintf(message, "Not enough data in %s", filename);
				ierr = -4;
				break;
			}
			else if (count < N_read_data) {
				sprintf(message, "Incorrect formatting on line %d of %s", i+2, filename);
				ierr = -5;
				break;
			}
		}  // for i

		// Close file
		fclose(fptr);
	}
	Display_assert_error(ierr, message, params, DTRACE("Display_assert_error"));

	// Broadcast particle information
	MPI_Bcast(temp, N_read_data * Np, MPI_DOUBLE, 0, PCW);

	// Create local particle linked lists
	p_list -> start = (Particle *)malloc(sizeof(Particle));
	Memory_check_allocation(p_list->start);

	// Number of local particles on this processor subdomain
	np = 0;

	p = p_list -> start;
	//--------------------------------------------------------------------------
	// Add particles that exist on processor subdomain to linked list
	//--------------------------------------------------------------------------
	for (i = 0; i < Np; i++) {

		// Check that particle fits within subdomain
		if (G_xmax-G_xmin < R[i] || G_ymax-G_ymin < R[i] || G_zmax-G_zmin < R[i]) {
			ierr = -6;
			break;
		}

		if (X[i][0] >= G_xmin && X[i][0] <  G_xmax && X[i][1] >= G_ymin &&
			X[i][1] <  G_ymax && X[i][2] >= G_zmin && X[i][2] <  G_zmax) {

			// Assign values to particle
			p -> X[0] = X[i][0];
			p -> X[1] = X[i][1];
			p -> X[2] = X[i][2];
			p -> R = R[i];

			// Fixed and mobile particles should have unique IDs
			p -> ID = i + p_list->ID_start;

			// Initialize velocities and forces to zero
			DSET_ZERO(p->U, 3);
			DSET_ZERO(p->Omega, 3);

			DSET_ZERO(p->Int_U_old, 3);
			DSET_ZERO(p->Int_Omega_old, 3);

			DSET_ZERO(p->Fc, 3);
			DSET_ZERO(p->Tc, 3);


#ifdef PARTICLE_RELEASE
			p -> t_part_release = t_part_release[i];
#endif

#ifdef SQUIRMER_SWIMMER
			// Initialize rotation matrix individually for each swimmer
			p -> Rotn[0][0] = cos(ini_rot_theta[i])*cos(ini_rot_phi[i]);
			p -> Rotn[0][1] = -sin(ini_rot_theta[i])*cos(ini_rot_phi[i]);
			p -> Rotn[0][2] = -sin(ini_rot_phi[i]);

			p -> Rotn[1][0] = sin(ini_rot_theta[i]);
			p -> Rotn[1][1] = cos(ini_rot_theta[i]);
			p -> Rotn[1][2] = 0;

			p -> Rotn[2][0] = cos(ini_rot_theta[i])*sin(ini_rot_phi[i]);
			p -> Rotn[2][1] = -sin(ini_rot_theta[i])*sin(ini_rot_phi[i]);
			p -> Rotn[2][2] = cos(ini_rot_phi[i]);
	#ifdef PARTICLE_RELEASE
			p -> t_part_release = t_part_release[i];
	#endif

#else
			// Initialize rotation matrix
			p -> Rotn[0][0] = 1;
			p -> Rotn[1][1] = 1;
			p -> Rotn[2][2] = 1;
			p -> Rotn[0][1] = 0;
			p -> Rotn[0][2] = 0;
			p -> Rotn[1][0] = 0;
			p -> Rotn[1][2] = 0;
			p -> Rotn[2][0] = 0;
			p -> Rotn[2][1] = 0;
#endif

			p -> particle_collision = NULL;
			p -> wall_collision = NULL;

			Particle_initialize_nonessential_data(p);
			Particle_calc_derived_data(p, grid, params);
			Particle_create_internal_arrays(p);

			// Allocate next space
			p -> next = (Particle *)malloc(sizeof(Particle));
			Memory_check_allocation(p -> next);
			p_last = p;
			p = p -> next;

			np++;
		}  // if
	}  // for i

	Display_assert_error(ierr, "Particle too large for subdomain!", params, DTRACE("Display_assert_error"));

	// Set up linked list terminal pointers
	if (np == 0)
		p_list -> start = NULL;
	else
		p_last -> next = NULL;

	free(p);
	free(temp);
	free(X);

	// Set linked list state
	p_list -> state = LIST_STATE_LOCAL;
}




/******************************************************************************/
/*
 Reads particles from Particle_XX.h5
 */
/******************************************************************************/
void ParticleInput_h5(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	hid_t file_id; // file handles

	char message[100];
	char h5_resume_filename[50];
	char groupname[50];

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Open HDF5 file
	sprintf(h5_resume_filename, "Particle_%d.h5",abs(params->noutput));
	sprintf(message, "Read from File: %s\n", h5_resume_filename);
	Display_progress(params, message);
	file_id = H5Fopen(h5_resume_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	//--------------------------------------------------------------------------
	// Particle data
	//--------------------------------------------------------------------------
	sprintf(groupname, "/fixed");
	p_fixed_list->ID_start = 0;
	p_fixed_list->type = FIXED;
	ParticleInput_h5_data(p_fixed_list, file_id, groupname, grid, params, DTRACE("ParticleInput_h5_data"));

	sprintf(groupname, "/mobile");
	p_mobile_list->ID_start = p_fixed_list->Np;
	p_mobile_list->type = MOBILE;

	ParticleInput_h5_data(p_mobile_list, file_id, groupname, grid, params, DTRACE("ParticleInput_h5_data"));

	H5Fclose(file_id);
}




/******************************************************************************/
/*
 Reads particle information of a single linked list from restart file

 NOTE: 'p_list->ID_start' should be set before calling this function.
 */
/******************************************************************************/
void ParticleInput_h5_data(Particle_list *p_list, hid_t file_id, char *groupname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace) {

	Particle *p, *p_last;
	int ierr;
	int i, j;
	int Np_local;  // Number of local particles
	int Np_global;
	double *X;

	// Check if group exists
	if (H5Lexists(file_id, groupname, H5P_DEFAULT) == 0) {
		p_list->Np = 0;
		p_list->start = NULL;
		p_list->state = LIST_STATE_LOCAL;
		return;
	}

	// Grid bounds of local subdomain
	double G_xmin = grid -> xu[grid -> G_Is];
	double G_xmax = grid -> xu[min(grid -> G_Ie, grid -> NX - 1)];
	double G_ymin = grid -> yv[grid -> G_Js];
	double G_ymax = grid -> yv[min(grid -> G_Je, grid -> NY - 1)];
	double G_zmin = grid -> zw[grid -> G_Ks];
	double G_zmax = grid -> zw[min(grid -> G_Ke, grid -> NZ - 1)];

	double **R_data, **X_data, **Rotn_data, **U_data, **Omega_data;
	double **Int_U_old_data, **Int_Omega_old_data, **Fc_data, **Tc_data;
#define RESUME_ELEMENT(data, str, dataSize) \
	ParticleInput_h5_data_element(data, dataSize, file_id, groupname, str, \
		params, DTRACE("ParticleInput_h5_data_element"))

	Np_global = RESUME_ELEMENT(&R_data, "R", 1);
	RESUME_ELEMENT(&X_data, "X", 3);
	RESUME_ELEMENT(&Rotn_data, "Rotn", 9);
	RESUME_ELEMENT(&U_data, "U", 3);
	RESUME_ELEMENT(&Omega_data, "Omega", 3);
	RESUME_ELEMENT(&Int_U_old_data, "Int_U_old", 3);
	RESUME_ELEMENT(&Int_Omega_old_data, "Int_Omega_old", 3);
	RESUME_ELEMENT(&Fc_data, "Fc", 3);
	RESUME_ELEMENT(&Tc_data, "Tc", 3);
#undef RESUME_ELEMENT

	p_list->start = (Particle *)malloc(sizeof(Particle));
	Memory_check_allocation(p_list->start);
	ierr = 0;
	Np_local = 0;
	p = p_list->start;
	for (j = 0; j < Np_global; j++) {

		// Check that particle fits within subdomain
		if (G_xmax-G_xmin < R_data[j][0] || G_ymax-G_ymin < R_data[j][0] || G_zmax-G_zmin < R_data[j][0]) {
			ierr = -1;
			break;
		}

		// Check if particle belongs to this processor
		X = X_data[j];
		if (X[0] >= G_xmin && X[0] <  G_xmax && X[1] >= G_ymin &&
			X[1] <  G_ymax && X[2] >= G_zmin && X[2] <  G_zmax) {

			// Copy particle information to Particle structure
			p -> ID = j + p_list->ID_start;
			p -> R = R_data[j][0];
			FORI3 p -> X[i]     = X_data[j][i];
			FORI3 p -> U[i]     = U_data[j][i];
			FORI3 p -> Omega[i] = Omega_data[j][i];
			FORI3 p -> Int_U_old[i]     = Int_U_old_data[j][i];
			FORI3 p -> Int_Omega_old[i] = Int_Omega_old_data[j][i];
			FORI3 p -> Fc[i] = Fc_data[j][i];
			FORI3 p -> Tc[i] = Tc_data[j][i];
			for (i = 0; i < 9; i++) {
				p -> Rotn[i/3][i%3] = Rotn_data[j][i];
			}

			// Initialize non-essential particle data
			Particle_initialize_nonessential_data(p);

			// Initialize derivative particle data
			Particle_calc_derived_data(p, grid, params);
			Particle_create_internal_arrays(p);
			Np_local++;

			// Allocate space for next particle
			p_last = p;
			p = (Particle *)malloc(sizeof(Particle));
			Memory_check_allocation(p);
			p_last -> next = p;
		}  // if
	}  // for j

	Display_assert_error(ierr, "Particle too large for subdomain!", params, DTRACE("Display_assert_error"));

	// Set up terminal linked list pointer
	if (Np_local == 0)
		p_list->start = NULL;
	else
		p_last -> next = NULL;

	// Set up linked list number of particles and state
	p_list->Np = Np_global;
	p_list->state = LIST_STATE_LOCAL;

	free(p);
	free(R_data);
	free(X_data);
	free(Rotn_data);
	free(U_data);
	free(Omega_data);
	free(Int_U_old_data);
	free(Int_Omega_old_data);
	free(Fc_data);
	free(Tc_data);

	// Read collisions
	ParticleInput_h5_collision(p_list, file_id, groupname, params, DTRACE("ParticleInput_h5_collision"));

}




/******************************************************************************/
/*
 * Reads particle information from restart file into array 'data', which will be
 * allocated here.  'dataSize' is the number of doubles per particle that
 * 'element' requires for storage (e.g. dataSize = 3 for "X").
 */
/******************************************************************************/
int ParticleInput_h5_data_element(double ***data, int dataSize, hid_t file_id,
		char *groupname, char *element, Parameters *params, Debug_trace *dtrace) {

	int ndim = 2;   // Dimension of array to read
	int Np_global;  // Number of particles, globally

	hsize_t dim_file[ndim];
	hid_t dataset, dataspace;
	herr_t status;

	char fieldname[100], message[100];
	sprintf(fieldname, "%s/%s", groupname, element);

	dataset = H5Dopen2(file_id, fieldname , H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	status    = H5Sget_simple_extent_dims(dataspace, dim_file, NULL);

	// Allocate space to load particle data
	Np_global = (int)dim_file[0];
	*data = Memory_allocate_2D_double_array(dataSize, Np_global);

#ifdef DEBUG2
	printf("Reading \"%s\" of dimensions %d x %d\n", fieldname, Np_global, dataSize);
#endif

	// Check that dimensions of file are what we expect
	if((int)dim_file[1] != dataSize) {
		sprintf(message, "Expected dimensions of \"%s\"\ndon't match file:\n"
				"dataSize = %d dataSize_file = %d\n", fieldname, dataSize, (int)dim_file[1]);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	// Don't need memspace because we are reading all the data in
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data[0][0]);
	sprintf(message, "Read of \"%s\" from *.h5 file failed", fieldname);
	Display_assert_error(status, message, params, DTRACE("Display_assert_error"));

	H5Dclose(dataset);
	H5Sclose(dataspace);

	return Np_global;
}




/******************************************************************************/
/*
 Reads particle collision information from restart file
 */
/******************************************************************************/
void ParticleInput_h5_collision(Particle_list *p_list, hid_t file_id,
		char *groupname, Parameters *params, Debug_trace *dtrace) {

	Particle *p;
	int j, k;

	double **data;  // Buffer for particle data to store
	int ndim = 2;   // Dimension of array to write

	char message[100];

	hsize_t dim_file[ndim], count[ndim], offset[ndim];
	hid_t dataset, dataspace;
	herr_t status;

	char fieldname[100];
	sprintf(fieldname, "%s/collision_data", groupname);

	dataset   = H5Dopen2(file_id, fieldname , H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	status    = H5Sget_simple_extent_dims(dataspace, dim_file, NULL);

	// Allocate space to load particle data
	data = Memory_allocate_2D_double_array((int)dim_file[1], (int)dim_file[0]);

	// Don't need memspace because we are reading all the data in
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data[0]);
	sprintf(message, "Read of \"%s\" from *.h5 file failed", fieldname);
	Display_assert_error(status, message, params, DTRACE("Display_assert_error"));

	p = p_list -> start;
	while (p != NULL) {

		k = p -> ID - p_list->ID_start;

		// Particle collisions
		j = 2;
		j = ParticleInput_copy_collision_data(&(p->particle_collision), j, (int)data[k][0], data[k]);

		// Wall collisions
		ParticleInput_copy_collision_data(&(p->wall_collision), j, (int)data[k][1], data[k]);

		p = p -> next;
	}

	free(data);
	H5Dclose(dataset);
	H5Sclose(dataspace);

}




/******************************************************************************/
/*
 * Reads particle collision information for an individual particle
 *
 * If adding or removing elements from Collision structure:
 *     1. Add/remove elements from for loop, e.g.:
 *          pc -> dn = data[k][j++];
 * This function should mirror the changes made to ParticleOutput_h5_collision()
 */
/******************************************************************************/
int ParticleInput_copy_collision_data(Collision **pc_start_ptr, int j,
		int Nc, double *data) {

	int i;
	Collision *pc, *pc_last;

	pc = (Collision *)malloc(sizeof(Collision));
	Memory_check_allocation(pc);
	*pc_start_ptr = pc;

	for (i = 0; i < Nc; i++) {

		pc -> kn = data[j++];
		pc -> dn = data[j++];
		pc -> kt = data[j++];
		pc -> dt = data[j++];
		pc -> St = data[j++];
		pc -> zeta_t[0] = data[j++];
		pc -> zeta_t[1] = data[j++];
		pc -> zeta_t[2] = data[j++];
		pc -> other_ID = (int)data[j++];
		pc -> type     = (int)data[j++];
		pc -> state    = COLL_STATE_NEUTRAL;

		// Allocate space for next collision
		pc_last = pc;
		pc = (Collision *)malloc(sizeof(Collision));
		Memory_check_allocation(pc);
		pc_last -> next = pc;
	}

	// Set up terminal linked list pointer
	if (Nc == 0)
		*pc_start_ptr = NULL;
	else
		pc_last -> next = NULL;

	free(pc);
	return j;
}
