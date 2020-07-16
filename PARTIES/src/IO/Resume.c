#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "hdf5.h"

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Communication.h"
#include "Display.h"
#include "Memory.h"
#include "Output.h"
#include "Resume.h"
#include "Velocity.h"
#include "Conc.h"

/******************************************************************************/
/*
 * Reads data from Resume_XX.h5
 */
/******************************************************************************/
void Resume_h5_resume(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	hid_t file_id; // file handles
	hsize_t dim_1d[1], dim_2d[2]; // dimensions of data to be written

	char message[100];
	char h5_resume_filename[50];
	char groupname[50];
	char fieldname[50];

	//--------------------------------------------------------------------------
	// Databag parameters
	//--------------------------------------------------------------------------
	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

#ifdef CONC
	// Concentration data
	Concentration **c = data_bag -> c;
	int NConc = params -> NConc;
	int iconc;
#endif

#ifdef LES
	// LES data
	Subgrid *smag = data_bag -> smag;
#endif

	int pnodes = params -> ghost_nodes;

	/*------------------------------------------------------------------------*/
	/*
	 Read data from Resume_XX.h5
	 */
	/*------------------------------------------------------------------------*/
	Display_progress(params, "Read in Resume.h5\n");
	sprintf(h5_resume_filename, "Resume.h5");
	file_id = H5Fopen(h5_resume_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	sprintf(groupname, "/Resume");

	//--------------------------------------------------------------------------
	// 1D data
	//--------------------------------------------------------------------------
	dim_1d[0] = 1;
	double data[1];
	int data_int[1];

	Resume_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/time", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->time = data[0];

	Resume_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/dt", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->dt = data[0];

	Resume_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/dt_old", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->dt_old = data[0];

	Resume_h5_dataset(data_int, H5T_NATIVE_INT, 1, dim_1d, file_id, "/ntime", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->ntime = data_int[0];

	Resume_h5_dataset(data_int, H5T_NATIVE_INT, 1, dim_1d, file_id, "/noutput", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->noutput = data_int[0];

	Resume_h5_dataset(data_int, H5T_NATIVE_INT, 1, dim_1d, file_id, "/noutput_2d", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->noutput_2d = data_int[0];

	Resume_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/output_time", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->output_time = data[0];

	Resume_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/output_time_2d", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->output_time_2d = data[0];

#ifdef CONSTANT_MASSFLUX
	Resume_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/dp_dx", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->dp_dx = data[0];
#endif

	// Read particle parameters
#ifdef LAG_PARTICLE_RESOLVED
	#ifdef STARTUP
	Resume_h5_dataset(data_int, H5T_NATIVE_INT, 1, dim_1d, file_id, "/startup_flag", params, DTRACE("Resume_h5_dataset"));
	data_bag->params->startup_flag = data_int[0];
	#endif
#endif

	//--------------------------------------------------------------------------
	// Timer data
	//--------------------------------------------------------------------------
	Resume_h5_timer(data_bag->timer, file_id, "/timer", params, DTRACE("Resume_h5_timer"));

	//--------------------------------------------------------------------------
	// 2D data
	//--------------------------------------------------------------------------
	dim_2d[0] = grid->NZ;
	dim_2d[1] = grid->NX;

#ifdef CONC
	// Read deposit Height
	for (iconc = 0; iconc < NConc; iconc++) {
		if (params->conc_output_deposit_height[iconc] == 1) {



			sprintf(fieldname, "%s/depositHeight_%d", groupname, iconc);
			if (grid->G_Js == 0) {
				Resume_h5_dataset(c[iconc]->G_deposit_height[0], H5T_NATIVE_DOUBLE,
				                  2, dim_2d, file_id, fieldname, params, DTRACE("Resume_h5_dataset"));

			}
		}
	}
#endif

	//--------------------------------------------------------------------------
	// 3D data
	//--------------------------------------------------------------------------
#ifdef LES_LAG_AVE
	sprintf(fieldname, "%s/ilm", groupname);
	Resume_h5_flow_variable(smag->ILM, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
	Communication_update_ghost_nodes_flow_variable(smag->ILM, 'c', pnodes, data_bag);

	sprintf(fieldname, "%s/imm", groupname);
	Resume_h5_flow_variable(smag->ILM, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
	Communication_update_ghost_nodes_flow_variable(smag->IMM, 'c', pnodes, data_bag);
#endif

	H5Fclose(file_id);
}

/******************************************************************************/
/*
 * Reads data from Data_XX.h5
 */
/******************************************************************************/
void Resume_h5_data(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	hid_t file_id; // file handles
	hsize_t dim_1d[1], dim_2d[2]; // dimensions of data to be written

	char message[100];
	char h5_resume_filename[50];
	char groupname[50];
	char fieldname[50];

	//--------------------------------------------------------------------------
	// Databag parameters
	//--------------------------------------------------------------------------
	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	// Flow data
	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
	Pressure *p = data_bag -> p;

#ifdef CONC
	// Concentration data
	Concentration **c = data_bag -> c;
	int NConc = params -> NConc;
	int iconc;
#endif

#ifdef LES
	// LES data
	Subgrid *smag = data_bag -> smag;
#endif

#ifdef RANS
	// RANS data
	Rans *rans = data_bag -> rans;
#endif

	int pnodes = params -> ghost_nodes;

	sprintf(h5_resume_filename, "Data_%d.h5", abs(params->noutput));
	sprintf(message, "Read from File: %s\n", h5_resume_filename);
	Display_progress(params, message);
	file_id = H5Fopen(h5_resume_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	Resume_h5_flow_variable(u->data, file_id, "/u", grid, params, DTRACE("Resume_h5_flow_variable"));
	Resume_h5_flow_variable(v->data, file_id, "/v", grid, params, DTRACE("Resume_h5_flow_variable"));
	Resume_h5_flow_variable(w->data, file_id, "/w", grid, params, DTRACE("Resume_h5_flow_variable"));
	Resume_h5_flow_variable(p->p_data, file_id, "/p", grid, params, DTRACE("Resume_h5_flow_variable"));

	Velocity_update_boundaries(u->data, 'u', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(v->data, 'v', VEL_TYPE_NORMAL, data_bag);
	Velocity_update_boundaries(w->data, 'w', VEL_TYPE_NORMAL, data_bag);
	Communication_update_ghost_nodes_flow_variable(p->p_data, 'c', pnodes, data_bag);

#ifdef CONC
	//--------------------------------------------------------------------------
	// Concentration data
	//--------------------------------------------------------------------------
	for (iconc = 0; iconc < NConc; iconc++) {
		sprintf(fieldname, "/Conc/%d",iconc);
					fflush(stdout);
		Resume_h5_flow_variable(c[iconc]->data, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));

		fflush(stdout);

		Communication_update_ghost_nodes_flow_variable(c[iconc]->data, 'c', pnodes, data_bag);
		Conc_set_boundary_values(c[iconc]->data, iconc ,CENTRAL_FULL ,grid , params);
	}
#endif

#ifdef LES
	//--------------------------------------------------------------------------
	// LES data
	//--------------------------------------------------------------------------
	sprintf(groupname, "/LES");

	sprintf(fieldname, "%s/nut", groupname);
	Resume_h5_flow_variable(smag->nut, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
	Communication_update_ghost_nodes_flow_variable(smag->nut, 'c', pnodes, data_bag);

	#ifdef CONC
	for (iconc = 0; iconc < NConc; iconc++) {
		sprintf(fieldname, "%s/Conc%d_alphat", groupname, iconc);
		Resume_h5_flow_variable(smag->cdev[iconc]->mSct, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->mSct, 'c', pnodes, data_bag);
	}
	#endif

	#if defined LES_LAG_AVE && defined CONC
	for (iconc = 0; iconc < NConc; iconc++) {
		sprintf(fieldname, "%s/Conc%d_itt/data", groupname, iconc);
		Resume_h5_flow_variable(smag->cdev[iconc]->ITT, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->ITT, 'c', pnodes, data_bag);
		sprintf(fieldname, "%s/Conc%d_ikt/data", groupname, iconc);
		Resume_h5_flow_variable(smag->cdev[iconc]->IKT, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->IKT, 'c', pnodes, data_bag);
	}
	#endif // LES_LAG_AVE and CONC
#endif // LES

#ifdef RANS
	//--------------------------------------------------------------------------
	// RANS data
	//--------------------------------------------------------------------------
	sprintf(groupname, "/RANS");

	#ifdef TWO_EQUATION_MODEL
	sprintf(fieldname, "%s/tke", groupname);
	Resume_h5_flow_variable(rans->two_eqn_rans[0]->data, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[0]->data, 'c', pnodes, data_bag);
	sprintf(fieldname, "%s/eps", groupname);
	Resume_h5_flow_variable(rans->two_eqn_rans[1]->data, file_id, fieldname, grid, params, DTRACE("Resume_h5_flow_variable"));
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[1]->data, 'c', pnodes, data_bag);
	#endif

	Rans_eddy_viscosity(rans, u, params, grid);
	Communication_update_ghost_nodes_flow_variable(rans->nut, 'c', pnodes, data_bag);
#endif

	H5Fclose(file_id);

}




/******************************************************************************/
/*
 * Read entire dataset from *.h5 file
 */
/******************************************************************************/
void Resume_h5_dataset(void *data, hid_t mem_type, int ndim, hsize_t *dim,
		hid_t file_id, char *fieldname, Parameters *params, Debug_trace *dtrace) {

	hid_t   dataset;        // handles
	hid_t   dataspace;
	hsize_t dim_file[2];    // dataset dimensions
	herr_t  status;
	int     status_n, ndim_file;
	char message[100], filename[50];
	H5Fget_name(file_id, filename, 50);
	int i;

	dataset = H5Dopen2(file_id, fieldname, H5P_DEFAULT);
	sprintf(message, "Open dataset \"%s\" from %s failed", fieldname, filename);
	Display_assert_error(dataset, message, params, DTRACE("Display_assert_error"));

	dataspace = H5Dget_space(dataset);    // dataspace handle
	ndim_file = H5Sget_simple_extent_ndims(dataspace);
	status_n  = H5Sget_simple_extent_dims(dataspace, dim_file, NULL);
	H5Sclose(dataspace);

	// Check that dimensions of data in file are what we expect
	if(ndim != ndim_file) {
		sprintf(message, "Expected number of dimensions don't match file:\n"
		                 "\tndim = %d, ndim_file = %d\n", ndim, ndim_file);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}
	for (i = 0; i < ndim; i++) {
		if(dim[i] != dim_file[i]) {
			sprintf(message, "Expected dimensions don't match file:\n"
			                 "\tdim[%d] = %d, dim_file[%d] = %d\n",
			                  i, (int)dim[i], i, (int)dim_file[i]);
			Display_throw_error(message, params, DTRACE("Display_throw_error"));
		}
	}

	// Don't need a memspace because 'data' and the dataspace are the same size,
	// use 'H5S_ALL' instead
	status = H5Dread(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	sprintf(message, "Read dataset \"%s\" from %s failed", fieldname, filename);
	Display_assert_error(status, message, params, DTRACE("Display_assert_error"));

	H5Dclose(dataset);
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Resume_h5_timer(Timer *timer, hid_t file_id, char *fieldname,
	Parameters *params, Debug_trace *dtrace) {

	hsize_t dim_1d[1] = {N_TIMER_ELEMENTS};
	double data[N_TIMER_ELEMENTS];
	Resume_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Resume_h5_dataset"));

	int j = 0;
	timer -> Wtime_total  = data[j++];
	timer -> Wtime_init   = data[j++];
	timer -> Wtime_intEOM = data[j++];
	timer -> Wtime_output = data[j++];
	timer -> Wtime_output_2d = data[j++];
	timer -> Wtime_cfl    = data[j++];

	timer -> Wtime_comm_3D        = data[j++];
	timer -> Wtime_comm_2D_reduce = data[j++];

	timer -> Wtime_vel_convective  = data[j++];
	timer -> Wtime_vel_rhs         = data[j++];
	timer -> Wtime_vel_solve       = data[j++];
	timer -> Wtime_vel_boundaries  = data[j++];
	timer -> Wtime_vel_cell_center = data[j++];

	timer -> Wtime_p_rhs        = data[j++];
	timer -> Wtime_p_solve      = data[j++];
	timer -> Wtime_p_project    = data[j++];
	timer -> Wtime_p_divergence = data[j++];

	timer -> Wtime_c_total       = data[j++];
	timer -> Wtime_c_convective  = data[j++];
	timer -> Wtime_c_rhs         = data[j++];
	timer -> Wtime_c_solve       = data[j++];
	timer -> Wtime_c_outofbounds = data[j++];

	timer -> Wtime_sgs_total   = data[j++];
	timer -> Wtime_sgs_filter  = data[j++];
	timer -> Wtime_sgs_dynamic = data[j++];

	timer -> Wtime_rans_total      = data[j++];
	timer -> Wtime_rans_convective = data[j++];
	timer -> Wtime_rans_rhs        = data[j++];
	timer -> Wtime_rans_solve      = data[j++];

	timer -> Wtime_particle_total = data[j++];
	timer -> Wtime_particle_comm  = data[j++];
	timer -> Wtime_particle_coll  = data[j++];
	timer -> Wtime_particle_forc  = data[j++];
	timer -> Wtime_particle_int   = data[j++];

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Resume_h5_flow_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace) {

	hid_t       dataset;     // handles
	hid_t       dataspace;
	hid_t       memspace;
	hsize_t     dim[3];
	hsize_t     dim_file[3];        //   dataset dimensions
	herr_t      status;

	hsize_t      count[3];         //  size of the hyperslab in the file
	hsize_t      offset[3];        //  hyperslab offset in the file
	hsize_t      count_out[3];     //  size of the hyperslab in memory
	hsize_t      offset_out[3];    //  hyperslab offset in memory

	double *data_start;

	int ndim = 3;
	int ndim_file;
	int status_n;

	int i;
	char message[100], filename[50];
	H5Fget_name(file_id, filename, 50);

	dim[0] = grid -> NZ;
	dim[1] = grid -> NY;
	dim[2] = grid -> NX;

	dataset = H5Dopen2(file_id, fieldname , H5P_DEFAULT);
	sprintf(message, "Open dataset \"%s\" from %s failed", fieldname, filename);
	Display_assert_error(dataset, message, params, DTRACE("Display_assert_error"));

	dataspace = H5Dget_space(dataset);  // dataspace handle
	ndim_file = H5Sget_simple_extent_ndims(dataspace);
	status_n  = H5Sget_simple_extent_dims(dataspace, dim_file, NULL);

	// Check that dimensions of data in file are what we expect
	if(ndim != ndim_file) {
		sprintf(message, "Expected number of dimensions don't match file:\n"
				"\tndim = %d, ndim_file = %d\n", ndim, ndim_file);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}
	for (i = 0; i < ndim; i++) {
		if(dim[i] != dim_file[i]) {
			sprintf(message, "Expected dimensions don't match file:\n"
					"\tdim[%d] = %d, dim_file[%d] = %d\n",
					i, (int)dim[i], i, (int)dim_file[i]);
			Display_throw_error(message, params, DTRACE("Display_throw_error"));
		}
	}

	//--------------------------------------------------------------------------
	// Dataspace hyperslab - select portion of data that is contained on local
	// processor
	//--------------------------------------------------------------------------
	count[0] = grid->G_Ke - grid->G_Ks;
	count[1] = grid->G_Je - grid->G_Js;
	count[2] = grid->G_Ie - grid->G_Is;

	offset[0] = grid->G_Ks;
	offset[1] = grid->G_Js;
	offset[2] = grid->G_Is;

	assert(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL) >= 0);

	//--------------------------------------------------------------------------
	// Define memory space, represents how local processor data is stored in
	// memory
	//--------------------------------------------------------------------------
	// We have that memspace[0][0][0] == data[L_Ks][L_Js][L_Is]
	data_start = &(data[grid->L_Ks][grid->L_Js][grid->L_Is]);

	count[0] = grid->L_Ke - grid->L_Ks;
	count[1] = grid->L_Je - grid->L_Js;
	count[2] = grid->L_Ie - grid->L_Is;

	memspace = H5Screate_simple(ndim, count, NULL);
	assert(memspace >= 0);

	count[0] = grid->G_Ke - grid->G_Ks;
	count[1] = grid->G_Je - grid->G_Js;
	count[2] = grid->G_Ie - grid->G_Is;

	// Offset by 'ghost_nodes' because
	// memspace[0][0][0] == data[L_Ks][L_Js][L_Is]
	offset[0] = params -> ghost_nodes;
	offset[1] = params -> ghost_nodes;
	offset[2] = params -> ghost_nodes;

	assert(H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL) >= 0);

	//--------------------------------------------------------------------------
	// Read from *.h5 file
	//--------------------------------------------------------------------------
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, data_start);
	sprintf(message, "Read dataset \"%s\" from %s failed", fieldname, filename);
	Display_assert_error(status, message, params, DTRACE("Display_assert_error"));

	assert(H5Dclose(dataset) >= 0);
	assert(H5Sclose(dataspace) >= 0);
	assert(H5Sclose(memspace) >= 0);
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Resume_h5_noghost_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace) {

	hid_t       dataset;     // handles
	hid_t       dataspace;
	hid_t       memspace;
	hsize_t     dim[3];
	hsize_t     dim_file[3];        //   dataset dimensions
	herr_t      status;

	hsize_t      count[3];         //  size of the hyperslab in the file
	hsize_t      offset[3];        //  hyperslab offset in the file
	hsize_t      count_out[3];     //  size of the hyperslab in memory
	hsize_t      offset_out[3];    //  hyperslab offset in memory

	double *data_start;

	int ndim = 3;
	int ndim_file;
	int status_n;

	int i;
	char message[100], filename[50];
	H5Fget_name(file_id, filename, 50);

	dim[0] = grid -> NZ;
	dim[1] = grid -> NY;
	dim[2] = grid -> NX;

	// Open dataset
	dataset = H5Dopen2(file_id, fieldname , H5P_DEFAULT);
	sprintf(message, "Open dataset \"%s\" from %s failed", fieldname, filename);
	Display_assert_error(dataset, message, params, DTRACE("Display_assert_error"));

	dataspace = H5Dget_space(dataset);  // dataspace handle
	ndim_file = H5Sget_simple_extent_ndims(dataspace);
	status_n  = H5Sget_simple_extent_dims(dataspace, dim_file, NULL);

	// Check that dimensions of data in file are what we expect
	if(ndim != ndim_file) {
		sprintf(message, "Expected number of dimensions don't match file:\n"
				"\tndim = %d, ndim_file = %d\n", ndim, ndim_file);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}
	for (i = 0; i < ndim; i++) {
		if(dim[i] != dim_file[i]) {
			sprintf(message, "Expected dimensions don't match file:\n"
					"tdim[%d] = %d, dim_file[%d] = %d\n",
					i, (int)dim[i], i, (int)dim_file[i]);
			Display_throw_error(message, params, DTRACE("Display_throw_error"));
		}
	}

	//--------------------------------------------------------------------------
	// Dataspace hyperslab - select portion of data that is contained on local
	// processor
	//--------------------------------------------------------------------------
	count[0] = grid->G_Ke - grid->G_Ks;
	count[1] = grid->G_Je - grid->G_Js;
	count[2] = grid->G_Ie - grid->G_Is;

	offset[0] = grid->G_Ks;
	offset[1] = grid->G_Js;
	offset[2] = grid->G_Is;

	assert(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL) >= 0);

	//--------------------------------------------------------------------------
	// Define memory space, represents how local processor data is stored in
	// memory
	//--------------------------------------------------------------------------
	// We have that memspace[0][0][0] == data[G_Ks][G_Js][G_Is]
	data_start = &(data[grid->G_Ks][grid->G_Js][grid->G_Is]);

	count[0] = grid->G_Ke - grid->G_Ks;
	count[1] = grid->G_Je - grid->G_Js;
	count[2] = grid->G_Ie - grid->G_Is;

	memspace = H5Screate_simple(ndim, count, NULL);
	assert(memspace >= 0);

	// Offset by '0' because
	// memspace[0][0][0] == data[G_Ks][G_Js][G_Is]
	offset[0] = 0;
	offset[1] = 0;
	offset[2] = 0;

	assert(H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL) >= 0);

	//--------------------------------------------------------------------------
	// Read from *.h5 file
	//--------------------------------------------------------------------------
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, data_start);
	// We don't need a memspace hyperslab because the written data is continuous
//	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace, H5P_DEFAULT,
//	                 &(data[grid->G_Ks][grid->G_Js][grid->G_Is]));
	sprintf(message, "Read dataset \"%s\" from %s failed", fieldname, filename);
	Display_assert_error(status, message, params, DTRACE("Display_assert_error"));

	assert(H5Dclose(dataset) >= 0);
	assert(H5Sclose(dataspace) >= 0);
	assert(H5Sclose(memspace) >= 0);
}
