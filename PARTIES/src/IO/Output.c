/*
 * Output.c: Create output files using HDF5.
 * Note: use HDF5 version 1.8.8! (version 1.8.7 leads to corrupted output files)
 *
 * Author: Roman Fuchs (roman.fuchs@hotmail.com)
 * Date: December 2011
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <hdf5.h>
#include <hdf5_hl.h>


#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"

#include "Communication.h"
#include "Conc.h"
#include "Display.h"
#include "Grid.h"
#include "Immersed.h"
#include "Memory.h"
#include "Output.h"
#include "Particle.h"
#include "Velocity.h"
#include "post_processing.h"

/******************************************************************************/
/*
 * Writes data to Resume_XX.h5
 * where XX is the output iteration defined by 'params->noutput'.
 *
 * This file contains data that is useful only for resuming the simulation, not
 * for viewing afterwards
 */
/******************************************************************************/
void Output_h5_resume(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int status;
	char message[500];

	hid_t file_id; // file handles
	hsize_t dim_1d[1], dim_2d[2]; // dimensions of data to be written

	char h5_resume_filename[50];
	char h5_resume_filename2[50];
	char fieldname[50];
	char groupname[50];

	//--------------------------------------------------------------------------
	// Databag parameters
	//--------------------------------------------------------------------------
	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	// Concentration data
	Concentration **c = data_bag -> c;
	int NConc = params -> NConc;
	int iconc;

	// LES data
	Subgrid *smag = data_bag -> smag;

	//--------------------------------------------------------------------------
	// Open HDF5 file handle
	//--------------------------------------------------------------------------
	Display_progress(params, "Write Resume H5\n");
	sprintf(h5_resume_filename, "Resume_%d.h5", abs(params->noutput));
	file_id = Output_h5_open(h5_resume_filename, params, DTRACE("Output_h5_open"));

	sprintf(groupname, "/Resume");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	//--------------------------------------------------------------------------
	// 1D data
	//--------------------------------------------------------------------------
	dim_1d[0] = 1;

	Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d,
	                  file_id, "/time", params, DTRACE("Output_h5_dataset"));

	Output_h5_dataset(&(params->dt), H5T_NATIVE_DOUBLE, 1, dim_1d,
	                  file_id, "/dt", params, DTRACE("Output_h5_dataset"));

	Output_h5_dataset(&(params->dt_old), H5T_NATIVE_DOUBLE, 1, dim_1d,
	                  file_id, "/dt_old", params, DTRACE("Output_h5_dataset"));

	Output_h5_dataset(&(params->ntime), H5T_NATIVE_INT, 1, dim_1d,
	                  file_id, "/ntime", params, DTRACE("Output_h5_dataset"));

	Output_h5_dataset(&(params->noutput), H5T_NATIVE_INT, 1, dim_1d,
	                  file_id, "/noutput", params, DTRACE("Output_h5_dataset"));

	Output_h5_dataset(&(params->noutput_2d), H5T_NATIVE_INT, 1, dim_1d,
		                  file_id, "/noutput_2d", params, DTRACE("Output_h5_dataset"));
	Output_h5_dataset(&(params->output_time), H5T_NATIVE_DOUBLE, 1, dim_1d,
	                  file_id, "/output_time", params, DTRACE("Output_h5_dataset"));

	Output_h5_dataset(&(params->output_time_2d), H5T_NATIVE_DOUBLE, 1, dim_1d,
		                 file_id, "/output_time_2d", params, DTRACE("Output_h5_dataset"));

	Output_h5_dataset(&(params->dp_dx), H5T_NATIVE_DOUBLE, 1, dim_1d,
	                  file_id, "/dp_dx", params, DTRACE("Output_h5_dataset"));

#ifdef LAG_PARTICLE_RESOLVED
	#ifdef STARTUP
	Output_h5_dataset(&(params->startup_flag), H5T_NATIVE_INT, 1, dim_1d,
	                  file_id, "/startup_flag", params, DTRACE("Output_h5_dataset"));
	#endif
#endif

	//--------------------------------------------------------------------------
	// Timer data
	//--------------------------------------------------------------------------
	Output_h5_timer(data_bag->timer, file_id, "/timer", params, DTRACE("Output_h5_timer"));

	//--------------------------------------------------------------------------
	// 2D data
	//--------------------------------------------------------------------------
	dim_2d[0] = grid->NZ;
	dim_2d[1] = grid->NX;

#ifdef CONC
	for (iconc = 0; iconc < NConc; iconc++) {

		if ( params->conc_output_deposit_height[iconc]) {

			Conc_update_world_deposited_height(c[iconc], data_bag) ;

				sprintf(fieldname, "%s/depositHeight_%d", groupname, iconc);
				Output_h5_dataset(c[iconc]->W_deposit_height, H5T_NATIVE_DOUBLE,2, dim_2d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		}
	}
#endif

	//--------------------------------------------------------------------------
	// 3D data
	//--------------------------------------------------------------------------
#ifdef LES_LAG_AVE
	sprintf(fieldname, "%s/ilm", groupname);
	Output_h5_flow_variable(smag->ILM, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));
	sprintf(fieldname, "%s/imm", groupname);
	Output_h5_flow_variable(smag->IMM, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));
#endif

	//--------------------------------------------------------------------------
	// Close HDF5 file handles
	//--------------------------------------------------------------------------
	H5Fclose(file_id);

	//--------------------------------------------------------------------------
	// Setting linker From Resume_XX.h5 to Resume.h5 and delete old linker
	//--------------------------------------------------------------------------
	sprintf(h5_resume_filename2, "Resume.h5");
	status = 0;
	if (params->rank==0){
		unlink(h5_resume_filename2);
		assert(symlink(h5_resume_filename, h5_resume_filename2) >= 0);
	}

	return;

}




/******************************************************************************/
/*
 * Writes data to Data_XX.h5
 * where XX is the output iteration defined by 'params->noutput'.
 *
 * This file contains data that is useful for visualization or post-processing
 */
/******************************************************************************/
void Output_h5_data(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	int verbose = 1;

	/*------------------------------------------------------------------------*/
	/*
	 Two different *.h5 files are created:
	     - Data_*.h5 contains 3D data
	     - Data2d_*.h5 contains 2D data

	 Hence, we don't need to move/read a large file when just 'small' 2D data is
	 accessed.
	 */
	/*------------------------------------------------------------------------*/

	// HDF5 Handles
	hid_t   file_id;    // File handle
	hsize_t dim_1d[1];  // Dimensions of data to be written
	herr_t  status;     // Status variable to check return value of HDF5 routines

	//--------------------------------------------------------------------------
	// Strings
	//--------------------------------------------------------------------------
	char h5filename[50];    // HDF5 filename to store 3D data
	char groupname[100];    // HDF5 group name to group data
	char subgroupname[100]; // HDF5 subgroup name to group data
	char fieldname[100];    // HDF5 name
	char message[500];      // Message for errors/warnings

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

#ifdef LAG_PARTICLE_RESOLVED
	Lagrangian *lag = data_bag -> lag;
#endif

#ifdef LES
	// LES data
	Subgrid *smag = data_bag -> smag;
#endif

#ifdef RANS
	// RANS data
	Rans *rans = data_bag -> rans;
#endif

	int NX = grid->NX;  // number of values in x-direction
	int NY = grid->NY;  // number of values in x-direction
	int NZ = grid->NZ;  // number of values in z-direction

	/*------------------------------------------------------------------------*/
	/*
	 * Open HDF5 file handle
	 */
	/*------------------------------------------------------------------------*/
	sprintf(h5filename, "Data_%d.h5", abs(params->noutput));
	file_id = Output_h5_open(h5filename, params, DTRACE("Output_h5_open"));


	/*------------------------------------------------------------------------*/
	/*
	 * Write timestamp to *.h5 file
	 */
	/*------------------------------------------------------------------------*/
	sprintf(fieldname, "/time");
	dim_1d[0] = 1;
	Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));


	/*------------------------------------------------------------------------*/
	/*
	 * Write spatial dimensions to *.h5 file
	 */
	/*------------------------------------------------------------------------*/
	sprintf(groupname, "/grid");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	dim_1d[0] = 1;
	sprintf(fieldname, "%s/NX", groupname);
	Output_h5_dataset(&(grid->NX), H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/NY", groupname);
	Output_h5_dataset(&(grid->NY), H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/NZ", groupname);
	Output_h5_dataset(&(grid->NZ), H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	dim_1d[0] = grid->NX;
	sprintf(fieldname, "%s/xu", groupname);
	Output_h5_dataset(grid->xu, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/xc", groupname);
	Output_h5_dataset(grid->xc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	dim_1d[0] = grid->NY;
	sprintf(fieldname, "%s/yv", groupname);
	Output_h5_dataset(grid->yv, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/yc", groupname);
	Output_h5_dataset(grid->yc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	dim_1d[0] = grid->NZ;
	sprintf(fieldname, "%s/zw", groupname);
	Output_h5_dataset(grid->zw, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/zc", groupname);
	Output_h5_dataset(grid->zc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));


	/*------------------------------------------------------------------------*/
	/*
	 * Write three-dimensional flow field variables
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// Flow data
	//--------------------------------------------------------------------------
	// Write u-velocity
	if (verbose) Display_progress(params,"Output.c: write u-velocity\n");
	Output_h5_flow_variable(u->data, file_id, "/u", grid, params, DTRACE("Output_h5_flow_variable"));

	// Write v-velocity
	if (verbose) Display_progress(params,"Output.c: write v-velocity\n");
	Output_h5_flow_variable(v->data, file_id, "/v", grid, params, DTRACE("Output_h5_flow_variable"));

	// Write w-velocity
	if (verbose) Display_progress(params,"Output.c: write w-velocity\n");
	Output_h5_flow_variable(w->data, file_id, "/w", grid, params, DTRACE("Output_h5_flow_variable"));

	// Write pressure
	if (verbose) Display_progress(params,"Output.c: write pressure\n");
	Output_h5_flow_variable(p->p_data, file_id, "/p", grid, params, DTRACE("Output_h5_flow_variable"));


	#ifdef TURB_FORCING
	// Write u-turbulent focing
	if (verbose) Display_progress(params,"Output.c: write u-forcing\n");
	Output_h5_flow_variable(u->fturb, file_id, "/u_forcing", grid, params, DTRACE("Output_h5_flow_variable"));

	// Write v-turbulent focing
	if (verbose) Display_progress(params,"Output.c: write v-forcing\n");
	Output_h5_flow_variable(v->fturb, file_id, "/v_forcing", grid, params, DTRACE("Output_h5_flow_variable"));

	// Write w-turbulent focing
	if (verbose) Display_progress(params,"Output.c: write w-forcing\n");
	Output_h5_flow_variable(w->fturb, file_id, "/w_forcing", grid, params, DTRACE("Output_h5_flow_variable"));
	#endif


#ifdef CONC
	//--------------------------------------------------------------------------
	// Concentration data
	//--------------------------------------------------------------------------
	sprintf(groupname, "/Conc");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	if (verbose) Display_progress(params,"Output.c: write concentration\n");


	for (iconc = 0; iconc < NConc; iconc++) {

		sprintf(fieldname, "%s/%d", groupname, iconc);
		Output_h5_flow_variable(c[iconc]->data, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));

	//	if (params->conc_output_integral[iconc]){

	//	dim_1d[0] = 1;
	//	Output_h5_dataset(&(c[iconc]->integral_conc_fluid), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/integral_conc_fluid", params, DTRACE("Output_h5_dataset"));
	//	Output_h5_dataset(&(c[iconc]->integral_conc_full), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/integral_conc_full", params, DTRACE("Output_h5_dataset"));
	//	 }
	//	if(params->conc_output_Q_int){


			// First, allocate required memory

	//		double recv_buffer[3];
	//		double nt=3;

	//		MPI_Allreduce( (void *) &(c[iconc]->Q_int[0]), (void *)recv_buffer, nt, MPI_DOUBLE, MPI_SUM, PCW);
	//		dim_1d[0] = 3;
	 //	 Output_h5_dataset(recv_buffer, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, "/Q_int", params, DTRACE("Output_h5_dataset"));
 	 //	 }


	} // end for iconc
#endif


#if defined LAG_PARTICLE_RESOLVED

	if (verbose) Display_progress(params,"Output.c: write volume fraction\n");


#ifdef VOF_SCALAR
		if(params->output_vfu == 1){	Output_h5_flow_variable(lag->vfu, file_id, "/vfu", grid, params, DTRACE("Output_h5_flow_variable")); }
		if(params->output_vfv == 1){	Output_h5_flow_variable(lag->vfv, file_id, "/vfv", grid, params, DTRACE("Output_h5_flow_variable")); }
		if(params->output_vfw == 1){	Output_h5_flow_variable(lag->vfw, file_id, "/vfw", grid, params, DTRACE("Output_h5_flow_variable")); }
#else
	if(params->output_vfu == 1){	Output_h5_noghost_variable(lag->ng_vfu, file_id, "/vfu", grid, params, DTRACE("Output_h5_noghost_variable")); }
	if(params->output_vfv == 1){	Output_h5_noghost_variable(lag->ng_vfv, file_id, "/vfv", grid, params, DTRACE("Output_h5_noghost_variable")); }
	if(params->output_vfw == 1){	Output_h5_noghost_variable(lag->ng_vfw, file_id, "/vfw", grid, params, DTRACE("Output_h5_noghost_variable")); }
#endif

if(params->output_vfc == 1){	Output_h5_noghost_variable(lag->ng_vfc, file_id, "/vfc", grid, params, DTRACE("Output_h5_noghost_variable")); }





#ifdef VOF_SCALAR_DEBUG
		if (verbose) Display_progress(params,"Output.c: write u_vof\n");
		Output_h5_flow_variable(u->data_vof, file_id, "/u_vof", grid, params, DTRACE("Output_h5_flow_variable"));
		if (verbose) Display_progress(params,"Output.c: write v_vof\n");
		Output_h5_flow_variable(v->data_vof, file_id, "/v_vof", grid, params, DTRACE("Output_h5_flow_variable"));
		if (verbose) Display_progress(params,"Output.c: write w_vof\n");
		Output_h5_flow_variable(w->data_vof, file_id, "/w_vof", grid, params, DTRACE("Output_h5_flow_variable"));
	#endif
#endif

#ifdef LES
	//--------------------------------------------------------------------------
	// LES data
	//--------------------------------------------------------------------------
	sprintf(groupname, "/LES");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	if (verbose) Display_progress(params,"Output.c: write nut\n");
	sprintf(fieldname, "%s/nut", groupname);
	Output_h5_flow_variable(smag->nut, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));

	#ifdef CONC
	for (iconc = 0; iconc < NConc; iconc++) {
		sprintf(fieldname, "%s/Conc%d_alphat", groupname, iconc);
		Output_h5_flow_variable(smag->cdev[iconc]->mSct, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));
	}
	#endif

	#if defined LES_LAG_AVE && defined CONC
	for (iconc = 0; iconc < NConc; iconc++) {
		sprintf(fieldname, "%s/Conc%d_itt", groupname, iconc);
		Output_h5_flow_variable(smag->cdev[iconc]->ITT, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));
		sprintf(fieldname, "%s/Conc%d_ikt", groupname, iconc);
		Output_h5_flow_variable(smag->cdev[iconc]->IKT, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));
	}
	#endif // LES_LAG_AVE and CONC
#endif // LES

#ifdef RANS
	//--------------------------------------------------------------------------
	// RANS data
	//--------------------------------------------------------------------------
	if (verbose) Display_progress(params,"Output.c: write rans\n");
	sprintf(groupname, "/RANS");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	#ifdef TWO_EQUATION_MODEL
	sprintf(fieldname, "%s/tke", groupname);
	Output_h5_flow_variable(rans->two_eqn_rans[0]->data, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));
	sprintf(fieldname, "%s/eps", groupname);
	Output_h5_flow_variable(rans->two_eqn_rans[1]->data, file_id, fieldname, grid, params, DTRACE("Output_h5_flow_variable"));
	#endif // TWO_EQUATION_MODEL
#endif // RANS


	/*------------------------------------------------------------------------*/
	/*
	 * Close HDF5 file handle
	 */
	/*------------------------------------------------------------------------*/
	assert(H5Fclose(file_id) >= 0);

	if (verbose) Display_progress(params,"Output.c: finished writing data\n");

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_output += T2 - T1;

	return;
}




/******************************************************************************/
/*
 * Writes data to Data2d_XX.h5
 * where XX is the output iteration defined by 'params->noutput'.
 *
 * This file contains 2-D data, typically from flow properties averaged over one
 * of the dimensions.
 */
/******************************************************************************/
void Output_h5_data2d(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	int verbose = 1;

	// HDF5 Handles
	hid_t   file_id;  // File handles
	hsize_t dim_1d[1], dim_2d[2]; // Dimensions of data to be written
	herr_t  status;               // Status variable to check return value of HDF5 routines

	//--------------------------------------------------------------------------
	// Strings
	//--------------------------------------------------------------------------
	char h5filename[50];    // HDF5 filename to store 3D data
	char groupname[100];    // HDF5 group name to group data
	char subgroupname[100]; // HDF5 subgroup name to group data
	char fieldname[100];    // HDF5 name
	char message[500];      // Message for errors/warnings

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

	int NX = grid->NX;  // number of values in x-direction
	int NY = grid->NY;  // number of values in x-direction
	int NZ = grid->NZ;  // number of values in z-direction

//	double totalSuspMass[1];
//	double frontLocation[1];
//	double front_limit;
//	double *averageHeight;
//	double **depositHeight;

	/*------------------------------------------------------------------------*/
	/*
	 * Open HDF5 file handle
	 */
	/*------------------------------------------------------------------------*/
	sprintf(h5filename, "Data2d_%d.h5", abs(params->noutput));
	file_id = Output_h5_open(h5filename, params, DTRACE("Output_h5_open"));


	/*------------------------------------------------------------------------*/
	/*
	 * Write timestamp to *.h5 file
	 */
	/*------------------------------------------------------------------------*/
	sprintf(fieldname, "/time");
	dim_1d[0] = 1;
	Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));


	/*------------------------------------------------------------------------*/
	/*
	 * Write spatial dimensions to *.h5 file
	 */
	/*------------------------------------------------------------------------*/
	sprintf(groupname, "/grid");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	dim_1d[0] = 1;
	sprintf(fieldname, "%s/NX", groupname);
	Output_h5_dataset(&(grid->NX), H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/NY", groupname);
	Output_h5_dataset(&(grid->NY), H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/NZ", groupname);
	Output_h5_dataset(&(grid->NZ), H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	dim_1d[0] = grid->NX;
	sprintf(fieldname, "%s/xc", groupname);
	Output_h5_dataset(grid->xc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	dim_1d[0] = grid->NY;
	sprintf(fieldname, "%s/yc", groupname);
	Output_h5_dataset(grid->yc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	dim_1d[0] = grid->NZ;
	sprintf(fieldname, "%s/zc", groupname);
	Output_h5_dataset(grid->zc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));


	/*------------------------------------------------------------------------*/
	/*
	 * Write two-dimensional flow field variables
	 */
	/*------------------------------------------------------------------------*/

	int i, k;
	double front_limit;

	//--------------------------------------------------------------------------
	// CONC
	//--------------------------------------------------------------------------
#ifdef CONC
	for (iconc = 0; iconc < NConc; iconc++) {

		sprintf(groupname, "/Conc%d", iconc);
		Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

		if (params->susp_mass_output) {
			if (verbose) Display_progress(params,"Output.c: write suspended mass\n");
			Conc_compute_total_suspended_mass(c[iconc], grid);

			dim_1d[0] = 1;
			sprintf(fieldname, "%s/suspendedMass", groupname);

			Output_h5_dataset(&(c[iconc]->W_total_susp_mass), H5T_NATIVE_DOUBLE,
			                 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		}

		if (params->ave_height_output) {
			if (verbose) Display_progress(params,"Output.c: write average height\n");
			Conc_compute_ave_height_x(c, grid, params);

			dim_1d[0] = NX;
			sprintf(fieldname, "%s/averageHeight", groupname);

			Output_h5_dataset(c[iconc]->W_ave_height_x, H5T_NATIVE_DOUBLE,1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		}

		if (params->front_location_output) {

			if (verbose) Display_progress(params,"Output.c: write front location\n");
			if (!params->ave_height_output) {
				Conc_compute_ave_height_x(c, grid, params);
			} // else, already was calculated before, avoid repeated work!

			// Limit to capture the front location
			front_limit = 0.01;
			Conc_find_front_location(c, grid, params, front_limit) ;
//			frontLocation[0] = c[iconc]->W_front_location;
			dim_1d[0] = 1;
			sprintf(fieldname, "%s/frontLocation", groupname);
//			write_1d_to_h5(frontLocation, dim_1d, file_id, fieldname, params);
			Output_h5_dataset(&(c[iconc]->W_front_location), H5T_NATIVE_DOUBLE,
			                 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		}

		if ( params->conc_output_deposit_height[iconc]) {

			if (verbose) Display_progress(params,"Output.c: write deposit height\n");
			Conc_update_world_deposited_height(c[iconc], data_bag);




				dim_2d[0] = grid->NZ;
				dim_2d[1] = grid->NX;
				sprintf(fieldname, "%s/depositHeight", groupname);
				Output_h5_dataset(c[iconc]->W_deposit_height[0], H5T_NATIVE_DOUBLE,2, dim_2d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		}

		MPI_Barrier(PCW);
	}
#endif // CONC

	//--------------------------------------------------------------------------
	// Write energies
	//--------------------------------------------------------------------------
	if (params->energies_output) {
		if (verbose) Display_progress(params,"Output.c: write energies\n");

		Velocity_update_world_energies(u, v, w, grid, params);
#ifdef CONC
		Conc_update_world_potential_energies(c, grid, params);
		Conc_update_world_stokes_dissipation_rate(c, grid, params);
#endif

		dim_1d[0] = 1;

		sprintf(groupname, "/Energies");
		Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

		sprintf(fieldname, "%s/kineticEnergy", groupname);
//		write_1d_to_h5(&(u->W_kinetic_energy), dim_1d, file_id, fieldname, params);
		Output_h5_dataset(&(u->W_kinetic_energy), H5T_NATIVE_DOUBLE, 1, dim_1d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		sprintf(fieldname, "%s/viscDissRate", groupname);
//		write_1d_to_h5(&(u->W_dissipation_rate), dim_1d, file_id, fieldname, params);
		Output_h5_dataset(&(u->W_dissipation_rate), H5T_NATIVE_DOUBLE, 1, dim_1d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

#ifdef CONC
		sprintf(fieldname, "%s/activePotentialEnergy", groupname);
//		write_1d_to_h5(&(c[0]->W_Ep_active), dim_1d, file_id, fieldname, params);
		Output_h5_dataset(&(c[0]->W_Ep_active), H5T_NATIVE_DOUBLE, 1, dim_1d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		sprintf(fieldname, "%s/passivePotentialEnergy", groupname);
//		write_1d_to_h5(&(c[0]->W_Ep_passive), dim_1d, file_id, fieldname, params);
		Output_h5_dataset(&(c[0]->W_Ep_passive), H5T_NATIVE_DOUBLE, 1, dim_1d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));
#endif
	}

	//--------------------------------------------------------------------------
	// Write shear stresses
	//--------------------------------------------------------------------------
	if (params->shear_stress_output) {
		if (verbose) Display_progress(params,"Output.c: write shear stresses\n");

//		Output_write_ascii_shear_stress_bottom(u, v, w, grid, params, noutput);
//		printf("Output.c/ after writing.c ascii\n");

		dim_2d[0] = grid->NZ;
		dim_2d[1] = grid->NX;
		sprintf(groupname, "/Shear");
		Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

		for (k = 0; k < NZ; k++) {
			for (i = 0; i < NX; i++) {
				u->G_u_streak[k][i] = grid->xc[i];
				u->W_u_streak[k][i] = grid->zc[k];
			}
		}
		sprintf(fieldname, "%s/xc_wall", groupname);
//		write_2d_to_h5(u->G_u_streak, dim_2d, group_id, fieldname, grid, params);
		Output_h5_dataset(u->G_u_streak[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "%s/zc_wall", groupname);
//		write_2d_to_h5(u->W_u_streak, dim_2d, group_id, fieldname, grid, params);
		Output_h5_dataset(u->W_u_streak[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		for (k = 0; k < NZ; k++) {
			for (i = 0; i < NX; i++) {
				u->G_u_streak[k][i] = 0.0;
				u->W_u_streak[k][i] = 0.0;
			}
		}

		Velocity_cell_center(data_bag);
		Velocity_wall_shear(data_bag);

		sprintf(fieldname, "%s/BottomShearStress", groupname);
//		write_2d_to_h5(u->W_shear_stress_bottom, dim_2d, group_id, fieldname, grid, params);
		Output_h5_dataset(u->W_shear_stress_bottom[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		sprintf(fieldname, "%s/uShear", groupname);
//		write_2d_to_h5(u->W_u_shear, dim_2d, group_id, fieldname, grid, params);
		Output_h5_dataset(u->W_u_shear[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		sprintf(fieldname, "%s/vShear", groupname);
//		write_2d_to_h5(u->W_v_shear, dim_2d, group_id, fieldname, grid, params);
		Output_h5_dataset(u->W_v_shear[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		sprintf(fieldname, "%s/wShear", groupname);
//		write_2d_to_h5(u->W_w_shear, dim_2d, group_id, fieldname, grid, params);
		Output_h5_dataset(u->W_w_shear[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		// Write fluctuating velocity
		Velocity_cell_center(data_bag);
		Velocity_u_streak(data_bag);
		sprintf(fieldname, "%s/ufluct", groupname);
//		write_2d_to_h5(u->W_u_streak, dim_2d, group_id, fieldname, grid, params);
		Output_h5_dataset(u->W_u_streak[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                 file_id, fieldname, params, DTRACE("Output_h5_dataset"));
	}

	/*------------------------------------------------------------------------*/
	/*
	 * Close HDF5 file handle
	 */
	/*------------------------------------------------------------------------*/
	assert(H5Fclose(file_id) >= 0);

	if (verbose) Display_progress(params,"Output.c: finished writing data2d\n");

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_output += T2 - T1;

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Output_h5_timer(Timer *timer, hid_t file_id, char *fieldname,
		Parameters *params, Debug_trace *dtrace) {

	hsize_t dim_1d[1] = {N_TIMER_ELEMENTS};
	double data[] = {
		timer -> Wtime_total,
		timer -> Wtime_init,
		timer -> Wtime_intEOM,
		timer -> Wtime_output,
		timer -> Wtime_output_2d,
		timer -> Wtime_cfl,

		timer -> Wtime_comm_3D,
		timer -> Wtime_comm_2D_reduce,

		timer -> Wtime_vel_convective,
		timer -> Wtime_vel_rhs,
		timer -> Wtime_vel_solve,
		timer -> Wtime_vel_boundaries,
		timer -> Wtime_vel_cell_center,

		timer -> Wtime_p_rhs,
		timer -> Wtime_p_solve,
		timer -> Wtime_p_project,
		timer -> Wtime_p_divergence,

		timer -> Wtime_c_total,
		timer -> Wtime_c_convective,
		timer -> Wtime_c_rhs,
		timer -> Wtime_c_solve,
		timer -> Wtime_c_outofbounds,

		timer -> Wtime_sgs_total,
		timer -> Wtime_sgs_filter,
		timer -> Wtime_sgs_dynamic,

		timer -> Wtime_rans_total,
		timer -> Wtime_rans_convective,
		timer -> Wtime_rans_rhs,
		timer -> Wtime_rans_solve,

		timer -> Wtime_particle_total,
		timer -> Wtime_particle_comm,
		timer -> Wtime_particle_coll,
		timer -> Wtime_particle_forc,
		timer -> Wtime_particle_int
	};

	Output_h5_dataset(data, H5T_NATIVE_DOUBLE, 1, dim_1d,
	                  file_id, fieldname, params, DTRACE("Output_h5_dataset"));
}




/******************************************************************************/
/*
 Write a 3D array ('data') in HDF5-format. It is written to a *.h5 file given by
 the file handle 'file_id', into a group ('groupname'). The character specifies
 which coordinates are used (cell centered, u velocity, v velocity, w velocity).
 */
/******************************************************************************/
void Output_h5_flow_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace) {

	int ndim = 3;
	double *data_start;
	hsize_t dim[ndim];
	hsize_t data_count[ndim], data_offset[ndim];
	hsize_t mem_count[ndim], mem_offset[ndim];

	//--------------------------------------------------------------------------
	// Define dataspace, represents required space in  *.h5 file to store 3-D
	// data from all processors
	//--------------------------------------------------------------------------
	// Dimensions of entire dataset
	dim[0] = grid -> NZ;
	dim[1] = grid -> NY;
	dim[2] = grid -> NX;

	// Dataspace hyperslab - select portion of data that is contained on local
	// processor
	data_count[0] = grid->G_Ke - grid->G_Ks;
	data_count[1] = grid->G_Je - grid->G_Js;
	data_count[2] = grid->G_Ie - grid->G_Is;

	data_offset[0] = grid->G_Ks;
	data_offset[1] = grid->G_Js;
	data_offset[2] = grid->G_Is;

	//--------------------------------------------------------------------------
	// Define memory space, represents how local processor data is stored in
	// memory
	//--------------------------------------------------------------------------
	// We have that memspace[0][0][0] == data[L_Ks][L_Js][L_Is]
	data_start = &(data[grid->L_Ks][grid->L_Js][grid->L_Is]);

	mem_count[0] = grid->L_Ke - grid->L_Ks;
	mem_count[1] = grid->L_Je - grid->L_Js;
	mem_count[2] = grid->L_Ie - grid->L_Is;

	// Offset by 'ghost_nodes' because
	// memspace[0][0][0] == data[L_Ks][L_Js][L_Is]
	mem_offset[0] = params -> ghost_nodes;
	mem_offset[1] = params -> ghost_nodes;
	mem_offset[2] = params -> ghost_nodes;

	Output_h5_3d_variable(data_start, ndim, dim, data_count, data_offset,
	                      mem_count,  mem_offset, file_id, fieldname, params, DTRACE("Output_h5_3d_variable"));
}




/******************************************************************************/
/*
 Write a 3D array ('data') in HDF5-format. It is written to a *.h5 file given by
 the file handle 'file_id', into a group ('groupname'). The character specifies
 which coordinates are used (cell centered, u velocity, v velocity, w velocity).
 */
/******************************************************************************/
void Output_h5_noghost_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace) {

	int ndim = 3;
	double *data_start;
	hsize_t dim[ndim];
	hsize_t data_count[ndim], data_offset[ndim];
	hsize_t mem_count[ndim], mem_offset[ndim];

	//--------------------------------------------------------------------------
	// Define dataspace, represents required space in  *.h5 file to store 3-D
	// data from all processors
	//--------------------------------------------------------------------------
	// Dimensions of entire dataset
	dim[0] = grid -> NZ;
	dim[1] = grid -> NY;
	dim[2] = grid -> NX;

	// Dataspace hyperslab - select portion of data that is contained on local
	// processor
	data_count[0] = grid->G_Ke - grid->G_Ks;
	data_count[1] = grid->G_Je - grid->G_Js;
	data_count[2] = grid->G_Ie - grid->G_Is;

	data_offset[0] = grid->G_Ks;
	data_offset[1] = grid->G_Js;
	data_offset[2] = grid->G_Is;

	//--------------------------------------------------------------------------
	// Define memory space, represents how local processor data is stored in
	// memory
	//--------------------------------------------------------------------------
	// We have that memspace[0][0][0] == data[G_Ks][G_Js][G_Is]
	data_start = &(data[grid->G_Ks][grid->G_Js][grid->G_Is]);

	mem_count[0] = grid->G_Ke - grid->G_Ks;
	mem_count[1] = grid->G_Je - grid->G_Js;
	mem_count[2] = grid->G_Ie - grid->G_Is;

	// Offset by '0' because
	// memspace[0][0][0] == data[G_Ks][G_Js][G_Is]
	mem_offset[0] = 0;
	mem_offset[1] = 0;
	mem_offset[2] = 0;

	Output_h5_3d_variable(data_start, ndim, dim, data_count, data_offset,
	                      mem_count,  mem_offset, file_id, fieldname, params, DTRACE("Output_h5_3d_variable"));
}




/******************************************************************************/
/*
 Write a 3D array ('data') in HDF5-format. It is written to a *.h5 file given by
 the file handle 'file_id', into a group ('groupname'). The character specifies
 which coordinates are used (cell centered, u velocity, v velocity, w velocity).
 */
/******************************************************************************/
void Output_h5_3d_variable(double *data_start, int ndim, hsize_t *dim,
		hsize_t *data_count, hsize_t *data_offset,
		hsize_t *mem_count,  hsize_t *mem_offset,
		hid_t file_id, char *fieldname, Parameters *params, Debug_trace *dtrace) {

	hid_t   plist_id, dcpl_id, dataset, dataspace, memspace;
	herr_t  status;
	char message[100], filename[50];

	//--------------------------------------------------------------------------
	// Define dataspace, represents required space in  *.h5 file to store 3-D
	// data from all processors
	//--------------------------------------------------------------------------
	dataspace = H5Screate_simple(ndim, dim, NULL);
	assert(dataspace >= 0);

	// Do not initialize dataset when creating it in HDF5
	dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
	assert(dcpl_id >= 0);
	assert(H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER) >= 0);

	// Create dataset using defined dataspace and fieldname
	dataset = H5Dcreate(file_id, fieldname, H5T_NATIVE_DOUBLE, dataspace,
	                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset >= 0);

	// Create property list for collective dataset write
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	assert(plist_id >= 0);
	assert(H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) >= 0);

	assert(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, data_offset, NULL, data_count, NULL) >= 0);

	//--------------------------------------------------------------------------
	// Define memory space, represents how local processor data is stored in
	// memory
	//--------------------------------------------------------------------------
	memspace = H5Screate_simple(ndim, mem_count, NULL);
	assert(memspace >= 0);

	assert(H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offset, NULL, data_count, NULL) >= 0);

	//--------------------------------------------------------------------------
	// Write to *.h5 file
	//--------------------------------------------------------------------------
	assert(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, data_start) >= 0);

	// Close handles
	assert(H5Dclose(dataset) >= 0);
	assert(H5Sclose(dataspace) >= 0);
	assert(H5Sclose(memspace) >= 0);
	assert(H5Pclose(dcpl_id) >= 0);
	assert(H5Pclose(plist_id) >= 0);
}




/******************************************************************************/
/*
 Write entire contents of 'data' to dataset

 Need to pass a pointer to the start of the data
 */
/******************************************************************************/
void Output_h5_dataset(void *data_start, hid_t datatype, int ndim, hsize_t *dim,
		hid_t file_id, char *fieldname, Parameters *params, Debug_trace *dtrace) {

	hid_t  dataspace, dataset;
	herr_t status;

	char message[100], filename[50];

	dataspace = H5Screate_simple(ndim, dim, NULL);
	assert(dataspace >= 0);
	dataset = H5Dcreate(file_id, fieldname, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset >= 0);
	assert(H5Sclose(dataspace) >= 0);

	if(params->rank == 0) {
		assert(H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_start) >= 0);
	}

	assert(H5Dclose(dataset) >= 0);
}




/******************************************************************************/
/*
 Opens an HDF5 file named 'filename' with proper MPI parallel writing access
 */
/******************************************************************************/
hid_t Output_h5_open(char *filename, Parameters *params, Debug_trace *dtrace) {

	// HDF5 handles:
	hid_t fapl_id; // file access property list handle
	hid_t file_id; // file handles
	herr_t status; // status variable to check return value of HDF5 routines

	char groupname[50];
	char fieldname[50];
	char message[100];

	// see http://mail.hdfgroup.org/pipermail/hdf-forum_hdfgroup.org/2011-February/004255.html
	MPI_Info info;

	// Create info to be attached to HDF5 file
	MPI_Info_create(&info);

	// Disables ROMIO's data-sieving
	MPI_Info_set(info, "romio_ds_read", "disable");
	MPI_Info_set(info, "romio_ds_write", "disable");

	// Enable ROMIO's collective buffering
	MPI_Info_set(info, "romio_cb_read", "enable");
	MPI_Info_set(info, "romio_cb_write", "enable");
//	MPI_Info info = MPI_INFO_NULL;

	/*------------------------------------------------------------------------*/
	/*
	 Set file access property list to MPI I/O using PETSc communication world.
	 Then, set file handles.

	 Note, that any old file with the same file name is overwritten
	 (H5F_ACC_TRUNC).  Therefore, the file handle is passed to all subroutines
	 and close at the end of this output routine.
	 */
	/*------------------------------------------------------------------------*/
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	assert(fapl_id >= 0);
//	assert(H5Pset_fapl_mpio(fapl_id, PCW, info) >= 0);
	assert(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info) >= 0);

	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	assert(file_id >= 0);

	assert(H5Pclose(fapl_id) >= 0);

	MPI_Info_free(&info);
	return file_id;
}



/******************************************************************************/
/*
 Creates the group 'groupname' in the HDF5 file 'file_id' and displays the
 error message 'err_text' if the group creation fails.
 */
/******************************************************************************/
void Output_h5_create_group(hid_t file_id, char *groupname, Parameters *params,
		Debug_trace *dtrace) {

	hid_t group_id;
	herr_t status;
	char message[100];

	group_id = H5Gcreate(file_id, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(group_id >= 0);

	assert(H5Gclose(group_id) >= 0);
}




/******************************************************************************/
/*
 This function writes the information for the immersed nodes to the files called
 from gvg.c
 */
/******************************************************************************/
void Output_immersed_info(MAC_grid *grid, Parameters *params, char which_quantity) {

	FILE *file;
	int g, gg;
	int rank;
	char filename[100];
	ImmersedNode *ib_node;
	Immersed *q_immersed;

	/* rank of the current processor */
	rank = params->rank;

	/* filename */
	sprintf(filename, "%c_immersed_P%d.dat", which_quantity, rank);
	file = fopen(filename, "w");
	if (file == NULL) {

		printf("Output.c/ Error opening %s file\n", filename);
	} /* if */

	switch (which_quantity) {

		case 'u':
			q_immersed = grid->u_immersed;
			break;
		case 'v':
			q_immersed = grid->v_immersed;
			break;
		case 'w':
			q_immersed = grid->w_immersed;
			break;
		case 'c':
			q_immersed = grid->c_immersed;
			break;
	} /* switch */

	for (g=0; g<q_immersed->N; g++) {

		ib_node = Immersed_get_ib_node(q_immersed, g);
		fprintf(file, "*********************************************************************\n");
		fprintf(file, "im(%f,%f,%f) at (%d,%d,%d) boundary=(%f,%f,%f) intersection:(%f,%f,%f)\n",
				ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z, ib_node->im_index.x_index,
				ib_node->im_index.y_index, ib_node->im_index.z_index, ib_node->boundary_point.x,
				ib_node->boundary_point.y, ib_node->boundary_point.z, ib_node->intersection_point.x,
				ib_node->intersection_point.y, ib_node->intersection_point.z);
		fprintf(file, "norm(%f,%f,%f) rhs_coef:%f\n", ib_node->n.vx, ib_node->n.vy, ib_node->n.vz, ib_node->boundary_coef);

		for (gg=0; gg<3; gg++) {
			fprintf(file, "Immersed index(%d,%d,%d) fluid index(%d,%d,%d) coef:%f\n",
					ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index,
					ib_node->fluid_index[gg].x_index, ib_node->fluid_index[gg].y_index,
					ib_node->fluid_index[gg].z_index, ib_node->fluid_coef[gg]);
		}
		fprintf(file, "\n*********************************************************************\n");
	} /* for g */

	fclose(file);
}




/******************************************************************************/
/*
 This function writes the information for the immersed nodes to the files called
 from gvg.c
 */
/******************************************************************************/
void Output_immersed_boundary(MAC_grid *grid, Parameters *params, char which_quantity) {

	FILE *file1, *file2;
	int g;
	int rank;
	char filename[100];
	ImmersedNode *ib_node;
	Immersed *q_immersed;

	/* rank of the current processor */
	rank = params->rank;

	/* filename */
	sprintf(filename, "%c_immersed_boundary_P%d.dat", which_quantity, rank);
	file1 = fopen(filename, "w");
	if (file1 == NULL) {

		printf("Output.c/ Error opening %s file\n", filename);
	} /* if */
	sprintf(filename, "%c_immersed_im_P%d.dat", which_quantity, rank);
	file2 = fopen(filename, "w");
	if (file2 == NULL) {

		printf("Output.c/ Error opening %s file\n", filename);
	} /* if */

	switch (which_quantity) {

		case 'u':
			q_immersed = grid->u_immersed;
			break;
		case 'v':
			q_immersed = grid->v_immersed;
			break;
		case 'w':
			q_immersed = grid->w_immersed;
			break;
		case 'c':
			q_immersed = grid->c_immersed;
			break;
	} /* switch */

	for (g=0; g<q_immersed->N; g++) {

		ib_node = Immersed_get_ib_node(q_immersed, g);
		fprintf(file1, "%f %f %f %f %f %f\n", ib_node->boundary_point.x,
				ib_node->boundary_point.y, ib_node->boundary_point.z,
				ib_node->n.vx, ib_node->n.vy, ib_node->n.vz);
		fprintf(file2, "%f %f %f\n", ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z);
	} /* for g */

	fclose(file1);
	fclose(file2);
}




/******************************************************************************/
/*
 This function writes 2D deposit height to ASCII file for each concentration
 field called from gvg.c
 */
/******************************************************************************/
void Output_write_ascii_deposit_height_dumped(Concentration **c, Parameters *params, MAC_grid *grid, int noutput) {

	FILE *fileout;
	int i, k, NI, NK;
	int iconc, NConc;
	int filenumber;
	char filename[50];

	/* Only write on processor zero */
	if (params->rank == MASTER) {

		NConc      = params->NConc;
		filenumber = 10000 + noutput;
		printf("NConc = %d", NConc);

		NK = grid->NZ-1;
		NI = grid->NX-1;
		printf("NK = %d, NI = %d", NK, NI);

		for (iconc=0; iconc<NConc; iconc++) {
			printf("iconc = %d", iconc);
			/* Only write to file if we have particle, i.e. u_s > 0.0 field */
			if (params->conc_output_dump[iconc]) {

				sprintf(filename, "S_deposit_height_dumped_c%d_%d.dat", iconc, filenumber);
				/* open the new file */
				fileout = fopen(filename, "w");

				/* Number of grid points in z and x directions */
				fprintf(fileout, "%d %d\n", NK, NI);
				fprintf(fileout, "%f %f\n", params->Lz, params->Lx);
				for (k=0; k<NK; k++) {

					for (i=0; i<NI; i++) {

						/* write the deposited height */
						fprintf(fileout, "%2.12f\n", c[iconc]->W_deposit_height_dumped[k][i]);
					} /* for i */
				} /* for k*/
				fclose(fileout);
			} /* if */

		} /* for */
	}/* if */
}




/******************************************************************************/
/*
 This function checks if the file has been opened correctly
 */
/******************************************************************************/
void Output_check_file_open(FILE *file) {

	if (file == NULL) {
		printf("Error opening the file...\n");
		printf("Exiting the program.\n \n");
		Communication_finalize();
		exit(1);
	}
}




/******************************************************************************/
/*
 This Filewriting function is not in final state!
 */
/******************************************************************************/
void Output_timehistory(double time, double frontLocation, double averageHeight,
		Parameters *params) {

/* variables to print out:
 *	time
 *	dkdt
 *	dpEdt
 *	epsKe
 *	epsPe
 *	int epsKe
 *	int epsPe
 *	suspMass
 *	depmass
 *	c[iconc]->W_front_location
 *	c[iconc]->W_total_susp_mass;
 *	c[iconc]->W_ave_height_x
 *
 *
 *
 */
	FILE *fileout;
	char filename[100];

	sprintf(filename, "Output_timehistory.dat");

	if (params->rank == 0) {
		/* append or create new file */
		if (time == 0.0)
			fileout = fopen(filename, "w");
		else
			fileout = fopen(filename, "a");
		Output_check_file_open(fileout);

		fprintf(fileout, "%f %f %f\n", time, frontLocation, averageHeight);

		fclose(fileout);
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
/*
void Output_timehistory_create(Parameters *params) {

	FILE *fileout;
	char filename[100];

	sprintf(filename, "Output_timehistory.dat");
	if (params->rank == 0) {
		fileout = fopen(filename, "w");
		Output_check_file_open(fileout);

		fclose(fileout);
	}
}
*/




/******************************************************************************/
/*
 */
/******************************************************************************/
void Output_copy_from_3dboundary(double ***data3d_b, double ***data3d,
		MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++){
				data3d[k][j][i] = data3d_b[k][j][i];
			}
		}
	}
	return;
}




/******************************************************************************/
/*
 * This function is designed for debugging purposes.
 *
 * It outputs a 3-D variable 'data3d' to an HDF5 file 'filename'.
 *
 * 'node_type' refers to whether the variable is:
 *     'c' - cell-centered
 *     'u' - u-faced
 *     'v' - v-faced
 *     'w' - w-faced
 *
 * 'ghost-type' refers to whether the variable is:
 *     "ghost" - ghost variable
 *     "ng"    - no-ghost variable
 * This is how the variable is stored in memory, not how you want written to
 * file.  Ghost variables will have their ghost nodes written to file, which is
 * useful for debugging boundary conditions.  The grid coordinates written to
 * file will include these ghost nodes.
 */
/******************************************************************************/
void Output_3d_data(double*** data3d, char node_type, char *ghost_type,
		char *filename, Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	// HDF5 handles:
	hid_t fapl_id; // file access property list handle
	hid_t file_id; // file handles
	hsize_t dim_1d[1]; // dimensions of data to be written
	hsize_t count_1d[1];
	herr_t status; // status variable to check return value of HDF5 routines
	int ndim = 3;
	hsize_t dim[ndim];
	hsize_t data_count[ndim], data_offset[ndim];
	hsize_t mem_count[ndim], mem_offset[ndim];
	double *data3d_start;

	char groupname[50];
	char fieldname[50];
	char message[100];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	file_id = Output_h5_open(filename, params, DTRACE("Output_h5_open"));

	double *x = grid -> xc;
	double *y = grid -> yc;
	double *z = grid -> zc;

	int NX = grid -> NX - 1;
	int NY = grid -> NY - 1;
	int NZ = grid -> NZ - 1;

	if (node_type == 'u') {
		x = grid -> xu;
		NX++;
	}
	else if (node_type == 'v') {
		y = grid -> yv;
		NY++;
	}
	else if (node_type == 'w') {
		z = grid -> zw;
		NZ++;
	}
	else if (node_type != 'c') {
		sprintf(message, "Incorrect input node_type = '%c'\n"
				"Must use 'u', 'v', 'w', or 'c'.", node_type);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = min(NX, grid -> G_Ie);
	int Je = min(NY, grid -> G_Je);
	int Ke = min(NZ, grid -> G_Ke);

	if (strcmp(ghost_type, "ghost") == 0) {

		int ghost_nodes = params -> ghost_nodes;

		x -= ghost_nodes;
		y -= ghost_nodes;
		z -= ghost_nodes;

		dim[0] = NZ + 2 * ghost_nodes;
		dim[1] = NY + 2 * ghost_nodes;
		dim[2] = NX + 2 * ghost_nodes;

		if (Is == 0) Is -= ghost_nodes;
		if (Js == 0) Js -= ghost_nodes;
		if (Ks == 0) Ks -= ghost_nodes;

		if (Ie == NX) Ie += ghost_nodes;
		if (Je == NY) Je += ghost_nodes;
		if (Ke == NZ) Ke += ghost_nodes;

		// We have that memspace[0][0][0] == data[L_Ks][L_Js][L_Is]
		data3d_start = &(data3d[grid->L_Ks][grid->L_Js][grid->L_Is]);

		data_count[0] = Ke - Ks;
		data_count[1] = Je - Js;
		data_count[2] = Ie - Is;

		data_offset[0] = Ks + ghost_nodes;
		data_offset[1] = Js + ghost_nodes;
		data_offset[2] = Is + ghost_nodes;

		mem_count[0] = grid->L_Ke - grid->L_Ks;
		mem_count[1] = grid->L_Je - grid->L_Js;
		mem_count[2] = grid->L_Ie - grid->L_Is;

		mem_offset[0] = Ks - grid->L_Ks;
		mem_offset[1] = Js - grid->L_Js;
		mem_offset[2] = Is - grid->L_Is;
	}
	else if (strcmp(ghost_type, "ng") == 0) {

		dim[0] = NZ;
		dim[1] = NY;
		dim[2] = NX;

		// We have that memspace[0][0][0] == data[G_Ks][G_Js][G_Is]
		data3d_start = &(data3d[grid->G_Ks][grid->G_Js][grid->G_Is]);

		data_count[0] = Ke - Ks;
		data_count[1] = Je - Js;
		data_count[2] = Ie - Is;

		data_offset[0] = Ks;
		data_offset[1] = Js;
		data_offset[2] = Is;

		mem_count[0] = grid->G_Ke - grid->G_Ks;
		mem_count[1] = grid->G_Je - grid->G_Js;
		mem_count[2] = grid->G_Ie - grid->G_Is;

		mem_offset[0] = Ks - grid->G_Ks;
		mem_offset[1] = Js - grid->G_Js;
		mem_offset[2] = Is - grid->G_Is;
	}
	else {
		sprintf(message, "Incorrect input ghost_type = %s\n"
		                 "Must use \"ghost\" or \"ng\".", ghost_type);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	//--------------------------------------------------------------------------
	// Write spatial data to *.h5 file
	//--------------------------------------------------------------------------
	sprintf(groupname, "/grid");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	dim_1d[0] = 1;
	sprintf(fieldname, "%s/NX", groupname);
	Output_h5_dataset(&dim[2], H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/NY", groupname);
	Output_h5_dataset(&dim[1], H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/NZ", groupname);
	Output_h5_dataset(&dim[0], H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/x", groupname);
	Output_h5_dataset(x, H5T_NATIVE_DOUBLE, 1, &dim[2], file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/y", groupname);
	Output_h5_dataset(y, H5T_NATIVE_DOUBLE, 1, &dim[1], file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/z", groupname);
	Output_h5_dataset(z, H5T_NATIVE_DOUBLE, 1, &dim[0], file_id, fieldname, params, DTRACE("Output_h5_dataset"));


	//--------------------------------------------------------------------------
	// Write three-dimensional flow field variable
	//--------------------------------------------------------------------------
	sprintf(fieldname, "/data");

	Output_h5_3d_variable(data3d_start, ndim, dim, data_count, data_offset,
	                      mem_count,  mem_offset, file_id, fieldname, params, DTRACE("Output_h5_3d_variable"));

	// Close HDF5 file handles
	assert(H5Fclose(file_id) >= 0);

	sprintf(message, "Output.c: finished writing %s\n", filename);
	Display_progress(params, message);

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_output += T2 - T1;
}




/******************************************************************************/
/*
 * This function is designed for debugging purposes.
 *
 * It outputs a 2-D variable 'data2d' to an HDF5 file 'filename'.
 *
 * 'dim' is an array of length 3 that contains the dimensions of the 2-D slice
 * as a subset of a 3-D slice.  It defines the lengths of 'x', 'y', and 'z'.
 *
 * 'x' - array of x-coordinates of 'data2d' (length dim[0])
 * 'y' - array of y-coordinates of 'data2d' (length dim[1])
 * 'z' - array of z-coordinates of 'data2d' (length dim[2])
 *
 * NOTE: 'data2d' should be initialized as a continuous array and start at
 * 'data2d[0][0]', i.e. use Memory_allocate_2D_double_array()
 *
 * NOTE: Use the function Communication_reduce_2D_arrays() to reduce the 2-D
 * array to a single processor before calling this function
 */
/******************************************************************************/
void Output_2d_data(double** data2d, int *dim, double *x, double *y, double *z,
		char *filename, Cart3d_bag *data_bag, Debug_trace *dtrace) {

	double T1, T2;
	T1 = MPI_Wtime();

	// HDF5 handles:
	hid_t file_id; // file handles
	herr_t status; // status variable to check return value of HDF5 routines
	hsize_t dim_1d[1], dim_2d[2], dim_3d[3];

	char groupname[50];
	char fieldname[50];
	char message[100];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Set dimensions for
	dim_3d[0] = dim[2];
	dim_3d[1] = dim[1];
	dim_3d[2] = dim[0];

	if (dim[0] == 1) {
		dim_2d[0] = dim[2];
		dim_2d[1] = dim[1];
	}
	else if (dim[1] == 1) {
		dim_2d[0] = dim[2];
		dim_2d[1] = dim[0];
	}
	else if (dim[2] == 1) {
		dim_2d[0] = dim[1];
		dim_2d[1] = dim[0];
	}
	else {
		sprintf(message, "No element of 'dim' equals 1\nOne dimension for 'dim' should equal 1\n");
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	if (dim[0] == 1 && dim[1] == 1 || dim[0] == 1 && dim[2] == 1 ||
		dim[1] == 1 && dim[2] == 1) {
		sprintf(message, "Too many elements of 'dim' equal 1\nOnly one dimension for 'dim' should equal 1\n");
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	file_id = Output_h5_open(filename, params, DTRACE("Output_h5_open"));

	//--------------------------------------------------------------------------
	// Write timestamp to *.h5 file
	//--------------------------------------------------------------------------
	sprintf(fieldname, "/time");
	dim_1d[0] = 1;
	Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	//--------------------------------------------------------------------------
	// Write spatial data to *.h5 file
	//--------------------------------------------------------------------------
	sprintf(groupname, "/grid");
	Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

	dim_1d[0] = 1;
	sprintf(fieldname, "%s/NX", groupname);
	Output_h5_dataset(&dim_3d[2], H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/NY", groupname);
	Output_h5_dataset(&dim_3d[1], H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/NZ", groupname);
	Output_h5_dataset(&dim_3d[0], H5T_NATIVE_INT, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/x", groupname);
	Output_h5_dataset(x, H5T_NATIVE_DOUBLE, 1, &dim_3d[2], file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/y", groupname);
	Output_h5_dataset(y, H5T_NATIVE_DOUBLE, 1, &dim_3d[1], file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/z", groupname);
	Output_h5_dataset(z, H5T_NATIVE_DOUBLE, 1, &dim_3d[0], file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	//--------------------------------------------------------------------------
	// Write two-dimensional data
	//--------------------------------------------------------------------------
	sprintf(fieldname, "/data");

	Output_h5_dataset(data2d[0], H5T_NATIVE_DOUBLE, 2, dim_2d,
					 file_id, fieldname, params, DTRACE("Output_h5_dataset"));

	// Close HDF5 file handles
	assert(H5Fclose(file_id) >= 0);

//	sprintf(message, "Output.c: finished writing %s\n", filename);
//	Display_progress(params, message);

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_output += T2 - T1;
}
