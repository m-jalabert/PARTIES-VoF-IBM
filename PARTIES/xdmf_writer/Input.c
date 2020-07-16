#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "Input.h"
#include "ini.h"
#include "DataTypes.h"
#include "Marker.h"
#include "Particle.h"
#include "hdf5.h"

#define TRUE 1
#define FALSE 0
#define CHECK_EXISTENCE(print_flag, fieldname) \
	params->print_flag = min( params->print_flag, \
	                          H5Lexists(file_id, fieldname, H5P_DEFAULT) )

#define INPUT_FILE "xdmfWriter.inp"

/******************************************************************************/
/*
 */
/******************************************************************************/
int Input_get_values(Parameters *params, Hyperslab *hyperslab) {

	int ierr;

	// Set defaults
	#define CFG(s, n, default, reader) s->n = default;
	#define CFG_STRING(s, n, default) strcpy(s->n, default);
	#include "default.inp"
	#undef CFG
	#undef CFG_STRING

	// Read input file
	if (ini_parse(INPUT_FILE, handler_params, params) < 0) {
		fprintf(stderr, "******** Error: Can't load '"
			INPUT_FILE ":[params]' ********\n");
		return EXIT_FAILURE;
	}
	if (ini_parse(INPUT_FILE, handler_hs, hyperslab) < 0) {
		fprintf(stderr, "******** Error: Can't load '"
			INPUT_FILE ":[hyperslab]' ********\n");
		return EXIT_FAILURE;
	}

	// Set other values based on inputs
Input_set_values(params, hyperslab);

	Input_set_values_otherData(params, hyperslab);
	Input_set_values_particle(params, hyperslab);
	Input_set_values_marker(params, hyperslab);

	return EXIT_SUCCESS;
}




/******************************************************************************/
/*
 Set values based on read input.  Here we make sure that files exist and that
 the data ranges specified exist.
 */
/******************************************************************************/
void Input_set_values(Parameters *params, Hyperslab *hyperslab) {

	int i, j;

	// Parameters variables
	int iter_start = params -> iter_start;
	int iter_step  = params -> iter_step;
	int iter_end   = params -> iter_end;
	int iter_count;

	// Hyperslab variables
	int NX, NY, NZ;
	int Is_HS, Js_HS, Ks_HS;
	int dI_HS, dJ_HS, dK_HS;
	int NI_HS, NJ_HS, NK_HS;

	double *timesteps = params -> timesteps;

	char h5_filename[50], fieldname[50];
	int search_h5_files;

	hid_t  file_id;    // file handle
	hid_t  dataset_id; // dataset handle
	herr_t status;     // status variable to check return value of HDF5 routines


	// -------------------------------------------------------------------------
	// Search for existence of HDF5 files
	// -------------------------------------------------------------------------
	i = iter_start;
	search_h5_files = 1;
	while (search_h5_files) {

		sprintf(h5_filename, "Data_%d.h5", i);

		if (access(h5_filename, F_OK) != -1) { // File exists
			i += iter_step;
		}
		else { // File does not exist
			search_h5_files = 0;
		}
	}
	if (i == iter_start) {
		printf("******** Could not locate Data_%d.h5 ********\n", i);
		hyperslab -> NX = 0;
		hyperslab -> NY = 0;
		hyperslab -> NZ = 0;
		return;
	}
	i -= iter_step; // Account for last i++

	// -------------------------------------------------------------------------
	// Update 'iter_end' if it was requested or it it exceeded the number of
	// files
	// -------------------------------------------------------------------------
	if (iter_end < 0 || iter_end > i) {
		iter_end = i;
		params -> iter_end = iter_end;
	}


	// -------------------------------------------------------------------------
	// Set 'iter_count'
	// -------------------------------------------------------------------------
	iter_count = (iter_end - iter_start) / iter_step + 1;
	params -> iter_count = iter_count;


	// -------------------------------------------------------------------------
	// From HDF5 files:
	//  - Read timestamps into 'params -> timesteps'
	//  - Read grid dimensions from first file into 'hyperslab -> NX, NY, NZ'
	// -------------------------------------------------------------------------
	timesteps = (double *)malloc(iter_count * sizeof(double));
	params -> timesteps = timesteps;

	for (i = 0; i < iter_count; i++) {

		sprintf(h5_filename, "Data_%d.h5", FILE_INDEX(i));
		file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

		dataset_id = H5Dopen2(file_id, "/time" , H5P_DEFAULT);
		status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
							 H5P_DEFAULT, &timesteps[i]);
		H5Dclose(dataset_id);

		if (i == 0) {

			dataset_id = H5Dopen2(file_id, "/grid/NX" , H5P_DEFAULT);
			status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
								 H5P_DEFAULT, &(hyperslab->NX));
			H5Dclose(dataset_id);

			dataset_id = H5Dopen2(file_id, "/grid/NY" , H5P_DEFAULT);
			status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
								 H5P_DEFAULT, &(hyperslab->NY));
			H5Dclose(dataset_id);

			dataset_id = H5Dopen2(file_id, "/grid/NZ" , H5P_DEFAULT);
			status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
								 H5P_DEFAULT, &(hyperslab->NZ));
			H5Dclose(dataset_id);

			//------------------------------------------------------------------
			// Basic data
			//------------------------------------------------------------------
			params->print_u = H5Lexists(file_id, "/u", H5P_DEFAULT);
//			CHECK_EXISTENCE(print_u, "u");
			params->print_v = H5Lexists(file_id, "/v", H5P_DEFAULT);
			params->print_w = H5Lexists(file_id, "/w", H5P_DEFAULT);
			params->print_p = H5Lexists(file_id, "/p", H5P_DEFAULT);
			params->print_uc = H5Lexists(file_id, "/uc", H5P_DEFAULT);
			params->print_vc = H5Lexists(file_id, "/vc", H5P_DEFAULT);
			params->print_wc = H5Lexists(file_id, "/wc", H5P_DEFAULT);

			//------------------------------------------------------------------
			// Concentration
			//------------------------------------------------------------------
			if (H5Lexists(file_id, "/Conc", H5P_DEFAULT)) {
				j = 0;
				sprintf(fieldname, "/Conc/%d", j);
				while (H5Lexists(file_id, fieldname, H5P_DEFAULT)) {
					j++;
					sprintf(fieldname, "/Conc/%d", j);
				}
				params->print_conc = 1;
				params->NConc = j;
			}
			else {
				params->print_conc = 0;
				params->NConc = 0;
			}

			//------------------------------------------------------------------
			// Particle volume fractions
			//------------------------------------------------------------------
			params->print_vf = H5Lexists(file_id, "/vfu", H5P_DEFAULT);

			//------------------------------------------------------------------
			// LES
			//------------------------------------------------------------------
			if (H5Lexists(file_id, "/LES", H5P_DEFAULT)) {
				params->print_nut = H5Lexists(file_id, "/LES/nut", H5P_DEFAULT);
				params->print_ilm = H5Lexists(file_id, "/LES/ilm", H5P_DEFAULT);
				params->print_imm = H5Lexists(file_id, "/LES/imm", H5P_DEFAULT);
				params->print_conc_itt = H5Lexists(file_id, "/LES/Conc0_itt", H5P_DEFAULT);
				params->print_conc_ikt = H5Lexists(file_id, "/LES/Conc0_ikt", H5P_DEFAULT);
				params->print_conc_alphat = H5Lexists(file_id, "/LES/Conc0_alphat", H5P_DEFAULT);
			}
			else {
				params->print_nut = 0;
				params->print_ilm = 0;
				params->print_imm = 0;
				params->print_conc_itt = 0;
				params->print_conc_ikt = 0;
				params->print_conc_alphat = 0;
			}
		}

		H5Fclose(file_id);
	}


	// -------------------------------------------------------------------------
	// Check hyperslab dimensions. Make them span the whole domain if -1 was
	// specified for NI_HS, NJ_HS, or NK_NS, or if they go beyond the domain
	// upper bounds.
	// -------------------------------------------------------------------------
	NX = hyperslab -> NX;
	NY = hyperslab -> NY;
	NZ = hyperslab -> NZ;

	Is_HS = hyperslab -> Is_HS;
	Js_HS = hyperslab -> Js_HS;
	Ks_HS = hyperslab -> Ks_HS;

	dI_HS = hyperslab -> dI_HS;
	dJ_HS = hyperslab -> dJ_HS;
	dK_HS = hyperslab -> dK_HS;

	NI_HS = hyperslab -> NI_HS;
	NJ_HS = hyperslab -> NJ_HS;
	NK_HS = hyperslab -> NK_HS;

	if (NI_HS < 0 || Is_HS + (NI_HS - 1) * dI_HS > NX - 2)
		NI_HS = (NX - 2 - Is_HS) / dI_HS + 1;

	if (NJ_HS < 0 || Js_HS + (NJ_HS - 1) * dJ_HS > NY - 2)
		NJ_HS = (NY - 2 - Js_HS) / dJ_HS + 1;

	if (NK_HS < 0 || Ks_HS + (NK_HS - 1) * dK_HS > NZ - 2)
		NK_HS = (NZ - 2 - Ks_HS) / dK_HS + 1;

	hyperslab -> NI_HS = NI_HS;
	hyperslab -> NJ_HS = NJ_HS;
	hyperslab -> NK_HS = NK_HS;
}




/******************************************************************************/
/*
 Set values based on read input.  Here we make sure that files exist and that
 the data ranges specified exist.
 */
/******************************************************************************/
void Input_set_values_otherData(Parameters *params, Hyperslab *hyperslab) {

	int i, j;

	// Parameters variables
	int iter_start = params -> iter_start;
	int iter_step  = params -> iter_step;
	int iter_end   = params -> iter_end;
	int iter_count = params -> iter_count;

	double *timesteps = params -> timesteps;
	int *N_pt_marker, *N_poly_marker;

	char h5_filename[50], fieldname[50];
	int search_h5_files;

	hid_t  file_id;     // file handle
	hid_t  group_id;    // group handle
	hid_t  dataset_id;  // dataset handle
	hid_t dataspace_id; // dataspace handle
	hsize_t ndim, *dim; // dataspace dimensions
	herr_t status;      // status variable to check return value of HDF5 routines

	// Not looking for otherData files
	if (strcmp(params->otherData_name, "") == 0) {
		params -> otherData = NULL;
		return;
	}

	// -------------------------------------------------------------------------
	// Search for existence of HDF5 files
	// -------------------------------------------------------------------------
	i = iter_start;
	search_h5_files = 1;
	while (search_h5_files) {

		sprintf(h5_filename, "%s_%d.h5", params->otherData_name, i);

		if (access(h5_filename, F_OK) != -1) { // File exists
			i += iter_step;
		}
		else { // File does not exist
			search_h5_files = 0;
		}
	}
	if (i == iter_start) {
		printf("******** Could not locate %s_%d.h5 ********\n", params->otherData_name, i);
		params -> otherData = NULL;
		return;
	}
	OtherData *otherData = (OtherData *)malloc(sizeof(OtherData));
	params -> otherData = otherData;
	otherData -> next = NULL;
	strcpy(otherData->name, params->otherData_name);

	i -= iter_step; // Account for last i++

	// -------------------------------------------------------------------------
	// Check for mismatch in number of output files between Particle_*.h5 and
	// Marker_*.h5
	// -------------------------------------------------------------------------
	if (hyperslab->NX > 0) {
		if (i < params->iter_end) {
			fprintf(stderr, "******** Error: Mismatch in number of output files:\n"
				"         Data_%d.h5 vs. %s_%d.h5 ********\n", params->iter_end, otherData->name, i);
			exit(EXIT_FAILURE);
		}
		else if (i > params->iter_end) {
			i = params->iter_end;
		}
	}


	if (hyperslab->NX == 0) {
		// ---------------------------------------------------------------------
		// Update 'iter_end' if it was requested or it it exceeded the number of
		// files
		// ---------------------------------------------------------------------
		if (iter_end < 0 || iter_end > i) {
			iter_end = i;
			params -> iter_end = iter_end;
		}


		// -------------------------------------------------------------------------
		// Set 'iter_count'
		// -------------------------------------------------------------------------
		iter_count = (iter_end - iter_start) / iter_step + 1;
		params -> iter_count = iter_count;


		// ---------------------------------------------------------------------
		// From HDF5 files:
		//  - Read timestamps into 'params -> timesteps'
		// ---------------------------------------------------------------------
		timesteps = (double *)malloc(iter_count * sizeof(double));
		params -> timesteps = timesteps;

		for (i = 0; i < iter_count; i++) {

			sprintf(h5_filename, "%s_%d.h5", otherData->name, FILE_INDEX(i));
			file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

			dataset_id = H5Dopen2(file_id, "/time" , H5P_DEFAULT);
			status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
								 H5P_DEFAULT, &timesteps[i]);
			H5Dclose(dataset_id);

			H5Fclose(file_id);
		}
	}


	//--------------------------------------------------------------------------
	// Read data dimensions
	//--------------------------------------------------------------------------
	sprintf(h5_filename, "%s_%d.h5", otherData->name, iter_start);
	file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	dataset_id = H5Dopen2(file_id, "/grid/NX" , H5P_DEFAULT);
	status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
						 H5P_DEFAULT, &(otherData->NX));
	H5Dclose(dataset_id);
	dataset_id = H5Dopen2(file_id, "/grid/NY" , H5P_DEFAULT);
	status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
						 H5P_DEFAULT, &(otherData->NY));
	H5Dclose(dataset_id);
	dataset_id = H5Dopen2(file_id, "/grid/NZ" , H5P_DEFAULT);
	status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
						 H5P_DEFAULT, &(otherData->NZ));
	H5Dclose(dataset_id);

	otherData->data_element = NULL;
	H5Literate(file_id, H5_INDEX_NAME, H5_ITER_INC, NULL, handler_otherData, &(otherData->data_element));

	H5Fclose(file_id);
}




/******************************************************************************/
/*
 Set values based on read input.  Here we make sure that files exist and that
 the data ranges specified exist.
 */
/******************************************************************************/
void Input_set_values_particle(Parameters *params, Hyperslab *hyperslab) {

	int i, j;

	// Parameters variables
	int iter_start = params -> iter_start;
	int iter_step  = params -> iter_step;
	int iter_end   = params -> iter_end;
	int iter_count = params -> iter_count;

	double *timesteps = params -> timesteps;
	int *N_p_mobile, *N_p_fixed;

	char h5_filename[50], fieldname[50];
	int search_h5_files;

	hid_t  file_id;     // file handle
	hid_t  group_id;    // group handle
	hid_t  dataset_id;  // dataset handle
	hid_t dataspace_id; // dataspace handle
	hsize_t ndim, *dim; // dataspace dimensions
	herr_t status;      // status variable to check return value of HDF5 routines


	// -------------------------------------------------------------------------
	// Search for existence of HDF5 files
	// -------------------------------------------------------------------------
	i = iter_start;
	search_h5_files = 1;
	while (search_h5_files) {

		sprintf(h5_filename, "Particle_%d.h5", i);

		if (access(h5_filename, F_OK) != -1) { // File exists
			i += iter_step;
		}
		else { // File does not exist
			search_h5_files = 0;
		}
	}
	if (i == iter_start) {
		printf("******** Could not locate Particle_%d.h5 ********\n", i);
		params -> mobile_data = NULL;
		params -> fixed_data = NULL;
		return;
	}
	i -= iter_step; // Account for last i++

	// -------------------------------------------------------------------------
	// Check for mismatch in number of output files between Particle_*.h5 and
	// Data_*.h5
	// -------------------------------------------------------------------------
	if (hyperslab->NX > 0 || params->otherData != NULL) {
		if (i < params->iter_end) {
			fprintf(stderr, "******** Error: Mismatch in number of output files:\n"
				"         Data_%d.h5/%s_%d.h5 vs. Particle_%d.h5 ********\n",
				params->iter_end, params->otherData_name, params->iter_end, i);
			exit(EXIT_FAILURE);
		}
		else if (i > params->iter_end) {
			i = params->iter_end;
		}
	}


	if (hyperslab->NX == 0 && params->otherData == NULL) {
		// ---------------------------------------------------------------------
		// Update 'iter_end' if it was requested or it it exceeded the number of
		// files
		// ---------------------------------------------------------------------
		if (iter_end < 0 || iter_end > i) {
			iter_end = i;
			params -> iter_end = iter_end;
		}


		// ---------------------------------------------------------------------
		// Set 'iter_count'
		// ---------------------------------------------------------------------
		iter_count = (iter_end - iter_start) / iter_step + 1;
		params -> iter_count = iter_count;


		// ---------------------------------------------------------------------
		// From HDF5 files:
		//  - Read timestamps into 'params -> timesteps'
		// ---------------------------------------------------------------------
		timesteps = (double *)malloc(iter_count * sizeof(double));
		params -> timesteps = timesteps;

		for (i = 0; i < iter_count; i++) {

			sprintf(h5_filename, "Particle_%d.h5", FILE_INDEX(i));
			file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

			dataset_id = H5Dopen2(file_id, "/time" , H5P_DEFAULT);
			status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
								 H5P_DEFAULT, &timesteps[i]);
			H5Dclose(dataset_id);

			H5Fclose(file_id);
		}
	}
	// -------------------------------------------------------------------------
	// From HDF5 files:
	//  - Read number of particles into 'params -> N_p_*'
	// -------------------------------------------------------------------------
	N_p_mobile = (int *)malloc(iter_count * sizeof(int));
	N_p_fixed  = (int *)malloc(iter_count * sizeof(int));
	params -> N_p_mobile = N_p_mobile;
	params -> N_p_fixed  = N_p_fixed;

	for (i = 0; i < iter_count; i++) {

		sprintf(h5_filename, "Particle_%d.h5", FILE_INDEX(i));
		file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

		// Mobile particles
		if (H5Lexists(file_id, MOBILE_GROUP, H5P_DEFAULT)) {
			if (!H5Lexists(file_id, MOBILE_GROUP"/R", H5P_DEFAULT)) {
				fprintf(stderr, "******** Error: cannot locate group:\n"
					"         Particle_%d.h5"MOBILE_GROUP"/R ********\n", FILE_INDEX(i));
				exit(EXIT_FAILURE);
			}
			dataset_id = H5Dopen2(file_id, MOBILE_GROUP"/R" , H5P_DEFAULT);
			dataspace_id = H5Dget_space(dataset_id);

			// Read dimensions
 			ndim = H5Sget_simple_extent_ndims(dataspace_id);
 			dim = (hsize_t *)malloc((int)ndim * sizeof(hsize_t));
 			H5Sget_simple_extent_dims(dataspace_id, dim, NULL);
			params -> N_p_mobile[i] = dim[0];

			free(dim);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		// Fixed particles
		if (H5Lexists(file_id, FIXED_GROUP, H5P_DEFAULT)) {
			if (!H5Lexists(file_id, FIXED_GROUP"/R", H5P_DEFAULT)) {
				fprintf(stderr, "******** Error: cannot locate group:\n"
					"         Particle_%d.h5"FIXED_GROUP"/R ********\n", FILE_INDEX(i));
				exit(EXIT_FAILURE);
			}
			dataset_id = H5Dopen2(file_id, FIXED_GROUP"/R" , H5P_DEFAULT);
			dataspace_id = H5Dget_space(dataset_id);

			// Read dimensions
 			ndim = H5Sget_simple_extent_ndims(dataspace_id);
 			dim = (hsize_t *)malloc((int)ndim * sizeof(hsize_t));
 			H5Sget_simple_extent_dims(dataspace_id, dim, NULL);
			params -> N_p_fixed[i] = dim[0];

			free(dim);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}

		H5Fclose(file_id);
	}

	//--------------------------------------------------------------------------
	// Check existence of particle types
	//--------------------------------------------------------------------------

	params -> mobile_data = NULL;
	params -> fixed_data = NULL;

	sprintf(h5_filename, "Particle_%d.h5", iter_start);
	file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (H5Lexists(file_id, MOBILE_GROUP, H5P_DEFAULT)) {
		group_id = H5Gopen2(file_id, MOBILE_GROUP, H5P_DEFAULT);
		H5Literate(group_id, H5_INDEX_NAME, H5_ITER_INC, NULL, handler_particle, &params->mobile_data);
		H5Gclose(group_id);
	}

	if (H5Lexists(file_id, FIXED_GROUP, H5P_DEFAULT)) {
		group_id = H5Gopen2(file_id, FIXED_GROUP, H5P_DEFAULT);
		H5Literate(group_id, H5_INDEX_NAME, H5_ITER_INC, NULL, handler_particle, &params->fixed_data);
		H5Gclose(group_id);
	}

	H5Fclose(file_id);
}




/******************************************************************************/
/*
 Set values based on read input.  Here we make sure that files exist and that
 the data ranges specified exist.
 */
/******************************************************************************/
void Input_set_values_marker(Parameters *params, Hyperslab *hyperslab) {

	int i, j;

	// Parameters variables
	int iter_start = params -> iter_start;
	int iter_step  = params -> iter_step;
	int iter_end   = params -> iter_end;
	int iter_count = params -> iter_count;

	double *timesteps = params -> timesteps;
	int *N_pt_marker, *N_poly_marker;

	char h5_filename[50], fieldname[50];
	int search_h5_files;

	hid_t  file_id;     // file handle
	hid_t  group_id;    // group handle
	hid_t  dataset_id;  // dataset handle
	hid_t dataspace_id; // dataspace handle
	hsize_t ndim, *dim; // dataspace dimensions
	herr_t status;      // status variable to check return value of HDF5 routines


	// -------------------------------------------------------------------------
	// Search for existence of HDF5 files
	// -------------------------------------------------------------------------
	i = iter_start;
	search_h5_files = 1;
	while (search_h5_files) {

		sprintf(h5_filename, "Marker_%d.h5", i);

		if (access(h5_filename, F_OK) != -1) { // File exists
			i += iter_step;
		}
		else { // File does not exist
			search_h5_files = 0;
		}
	}
	if (i == iter_start) {
		printf("******** Could not locate Marker_%d.h5 ********\n", i);
		params -> marker_data = NULL;
		return;
	}
	i -= iter_step; // Account for last i++

	// -------------------------------------------------------------------------
	// Check for mismatch in number of output files between Particle_*.h5 and
	// Marker_*.h5
	// -------------------------------------------------------------------------
	if (params->mobile_data != NULL || params->fixed_data != NULL) {
		if (i != params->iter_end) {
			fprintf(stderr, "******** Error: Mismatch in number of output files:\n"
				"         Particle_%d.h5 vs. Marker_%d.h5 ********\n", params->iter_end, i);
			exit(EXIT_FAILURE);
		}
	}


	// -------------------------------------------------------------------------
	// From HDF5 files:
	//  - Read number of marker points into 'params -> N_pt_marker'
	//  - Read number of polygons into 'params -> N_poly_marker'
	// -------------------------------------------------------------------------
	N_pt_marker = (int *)malloc(iter_count * sizeof(int));
	N_poly_marker = (int *)malloc(iter_count * sizeof(int));
	params -> N_pt_marker = N_pt_marker;
	params -> N_poly_marker = N_poly_marker;

	for (i = 0; i < iter_count; i++) {

		sprintf(h5_filename, "Marker_%d.h5", FILE_INDEX(i) );
		file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

		if (H5Lexists(file_id, MARKER_GROUP_0, H5P_DEFAULT)) {
			// -----------------------------------------------------------------
			// Number of marker points
			// -----------------------------------------------------------------
			if (!H5Lexists(file_id, MARKER_GROUP_0"/X", H5P_DEFAULT)) {
				fprintf(stderr, "******** Error: cannot locate group:\n"
					"         Marker_%d.h5"MARKER_GROUP_0"/X ********\n", FILE_INDEX(i));
				exit(EXIT_FAILURE);
			}
			dataset_id = H5Dopen2(file_id, MARKER_GROUP_0"/X" , H5P_DEFAULT);
			dataspace_id = H5Dget_space(dataset_id);

			// Read dimensions
 			ndim = H5Sget_simple_extent_ndims(dataspace_id);
 			dim = (hsize_t *)malloc((int)ndim * sizeof(hsize_t));
 			H5Sget_simple_extent_dims(dataspace_id, dim, NULL);
			params -> N_pt_marker[i] = dim[0];

			free(dim);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);

			// -----------------------------------------------------------------
			// Number of marker polygons
			// -----------------------------------------------------------------
			if (!H5Lexists(file_id, MARKER_GROUP_0"/topology", H5P_DEFAULT)) {
				fprintf(stderr, "******** Error: cannot locate group:\n"
					"         Marker_%d.h5"MARKER_GROUP_0"/topology ********\n", FILE_INDEX(i));
				exit(EXIT_FAILURE);
			}
			dataset_id = H5Dopen2(file_id, MARKER_GROUP_0"/topology" , H5P_DEFAULT);
			dataspace_id = H5Dget_space(dataset_id);

			// Read dimensions
 			ndim = H5Sget_simple_extent_ndims(dataspace_id);
 			dim = (hsize_t *)malloc((int)ndim * sizeof(hsize_t));
 			H5Sget_simple_extent_dims(dataspace_id, dim, NULL);
			params -> N_poly_marker[i] = dim[0];
			if (i > 0) {
				if (dim[1] != params->size_poly_marker) {
					fprintf(stderr, "******** Error: polygon size inconsistent:\n"
						"         Marker_%d.h5"MARKER_GROUP_0"/topology dim = [%d, %d] ********\n"
						"         Marker_%d.h5"MARKER_GROUP_0"/topology dim = [%d, %d] ********\n",
						FILE_INDEX(i), dim[0], dim[1], FILE_INDEX(i-1), params->N_poly_marker[i-1], params->size_poly_marker);
					exit(EXIT_FAILURE);
				}
			}
			else {
				params -> size_poly_marker = dim[1];
			}

			free(dim);
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}
		else {
			fprintf(stderr, "******** Error: cannot locate group:\n"
				"         Marker_%d.h5"MARKER_GROUP_0" ********\n", FILE_INDEX(i));
			exit(EXIT_FAILURE);
		}

		H5Fclose(file_id);
	}

	//--------------------------------------------------------------------------
	// Check existence of marker types
	//--------------------------------------------------------------------------

	params -> marker_data = NULL;

	sprintf(h5_filename, "Marker_%d.h5", iter_start);
	file_id   = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (H5Lexists(file_id, MARKER_GROUP_0, H5P_DEFAULT)) {
		group_id = H5Gopen2(file_id, MARKER_GROUP_0, H5P_DEFAULT);
		H5Literate(group_id, H5_INDEX_NAME, H5_ITER_INC, NULL, handler_particle, &params->marker_data);
		H5Gclose(group_id);
	}

	H5Fclose(file_id);
}




/******************************************************************************/
/*
 */
/******************************************************************************/
int handler_params(void* user, const char* section, const char* name,
		const char* value) {

	Parameters *params = (Parameters *)user;
	Hyperslab *hyperslab = (Hyperslab *)user;

	#define CFG(s, n, default, reader) \
		else if (strcmp(section, #s)==0 && strcmp(name, #n)==0) \
			s->n = reader(value);
	#define CFG_STRING(s, n, default) \
		else if (strcmp(section, #s)==0 && strcmp(name, #n)==0) \
			strcpy(s->n, value);

	if (strcmp(section,"params")) {}
	#include "default.inp"
	#undef CFG

	return EXIT_SUCCESS;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
int handler_hs(void* user, const char* section, const char* name,
		const char* value) {

	Parameters *params = (Parameters *)user;
	Hyperslab *hyperslab = (Hyperslab *)user;

	#define CFG(s, n, default, reader) \
		else if (strcmp(section, #s)==0 && strcmp(name, #n)==0) \
			s->n = reader(value);
	#define CFG_STRING(s, n, default) \
		else if (strcmp(section, #s)==0 && strcmp(name, #n)==0) \
			strcpy(s->n, value);

	if (strcmp(section,"hyperslab")) {}
	#include "default.inp"
	#undef CFG

	return EXIT_SUCCESS;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
herr_t handler_particle(hid_t g_id, const char *name, const H5L_info_t *info,
		void *op_data) {

#ifdef DEBUG
	if (info->type == H5G_GROUP)
		printf(" Object with name %s is a group\n", name);
	else if (info->type == H5G_DATASET)
		printf(" Object with name %s is a dataset\n", name);
	else if (info->type == H5G_TYPE)
		printf(" Object with name %s is a datatype\n", name);
	else
		printf(" Object with name %s is unidentifiable\n", name);
	return 0;
#else

	int i;
	hsize_t ndim, *dim;

	// Do not add these datasets as readable
	if (strcmp(name, "collision_data") == 0 ||
	    strcmp(name, "Rotn") == 0 ||
		strcmp(name, "topology") == 0) {
		return 0;
	}
	else {
		// Open dataspace
		hid_t dataset = H5Dopen2(g_id, name , H5P_DEFAULT);
		hid_t dataspace = H5Dget_space(dataset);

		// Read dimensions
		ndim = H5Sget_simple_extent_ndims(dataspace);
		dim = (hsize_t *)malloc((int)ndim * sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dim, NULL);

		// Add information to linked list
		DataParticle **p_data_start_ptr = (DataParticle **)(op_data);
		DataParticle *p_data = (DataParticle *)malloc(sizeof(DataParticle));
		strcpy(p_data->dataname, name);
		p_data -> length = dim[1];
		p_data -> next = *p_data_start_ptr;
		*p_data_start_ptr = p_data;

		free(dim);
		H5Dclose(dataset);
		H5Sclose(dataspace);
	}



#endif
}




/******************************************************************************/
/*
 */
/******************************************************************************/
herr_t handler_otherData(hid_t g_id, const char *name, const H5L_info_t *info,
		void *op_data) {

#ifdef DEBUG
	if (info->type == H5G_GROUP)
		printf(" Object with name %s is a group\n", name);
	else if (info->type == H5G_DATASET)
		printf(" Object with name %s is a dataset\n", name);
	else if (info->type == H5G_TYPE)
		printf(" Object with name %s is a datatype\n", name);
	else
		printf(" Object with name %s is unidentifiable\n", name);
	return 0;
#else

	// Do not add these datasets as readable
	if (strcmp(name, "grid") == 0 || strcmp(name, "time") == 0) {
		return 0;
	}
	else {
		OtherDataElement **data_start_ptr = (OtherDataElement **)(op_data);
		OtherDataElement *data = (OtherDataElement *)malloc(sizeof(OtherDataElement));
		strcpy(data->dataname, name);
		data -> next = *data_start_ptr;
		*data_start_ptr = data;
	}
#endif
}
