
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"

#include "Communication.h"
#include "Conc.h"
#include "Display.h"
#include "Grid.h"
#include "Lagrangian.h"
#include "Memory.h"
#include "Output.h"
#include "Particle.h"
#include "PpEulerian.h"
#include "PpLagrangian.h"
#include "Pressure.h"
#include "Timer.h"
#include "Velocity.h"
#include "Viscosity.h"
#include "Vorticity.h"
#include "Interpolate.h"
#include "Statistics2d.h"

/******************************************************************************/
/*
 * Perform post-processing operations
 */
/******************************************************************************/
void PpEulerian(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i, j, k;
	double temp;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

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

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
#ifdef CONC
	Concentration **c = data_bag -> c;
#endif

	Velocity_cell_center(data_bag);

	Communication_update_ghost_nodes_flow_variable(u->data_bc,  'c', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data_bc,  'c', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data_bc,  'c', 1, data_bag);

	Communication_update_ghost_nodes_flow_variable(u->data,  'u', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data,  'v', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data,  'w', 1, data_bag);

	#ifdef LAG_PARTICLE_RESOLVED
		// Include both local and foreign particles
		Particle_MPI_update(data_bag->lag->p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
		Particle_MPI_update(data_bag->lag->p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));

		// Calculate volume fractions
		Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfc);
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction('c', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));

		#ifdef VOF_PP_VFPRIME

		Memory_reset_flow_variable(grid, params, data_bag->lag->vfu_prime);
		Interpolate_add_to_volume_fraction_vof_prime('u', data_bag->lag->p_mobile_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
		Interpolate_add_to_volume_fraction_vof_prime('u', data_bag->lag->p_fixed_list,
		data_bag, DTRACE("Interpolate_add_to_volume_fraction"));

		#endif

		Particle_list_remove(data_bag->lag->p_mobile_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));
		Particle_list_remove(data_bag->lag->p_fixed_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));

	#endif


	double dA, dy;
#ifdef CONC
	double ***c_data;
	int iconc;
	int NConc = params->NConc;
#endif
	double ***u_data = u->data_bc;
	double ***v_data = v->data_bc;
	double ***w_data = w->data_bc;

#ifdef LAG_PARTICLE_RESOLVED
	double ***vfc = data_bag->lag->ng_vfc;
	#ifdef VOF_PP_VFPRIME
	double ***vfu_p = data_bag->lag->vfu_prime;
	#endif
#endif

	double vf = 0.0;



	double up;
	double vp;
	double wp;
	double cp;

	// Indices start and end on current processor
	int i_start = Is;
	int j_start = Js;
	int k_start = Ks;

	// Exclude the half cell added
	int i_end = min(NX-1, Ie);
	int j_end = min(NY-1, Je);
	int k_end = min(NZ-1, Ke);

	// HDF5 Handles
	hid_t   file_id;  // File handles
	hsize_t dim_1d[1], dim_2d[2]; // Dimensions of data to be written
	herr_t  status;               // Status variable to check return value of HDF5 routines

	int NZslice = 1;
	//--------------------------------------------------------------------------
	// Strings
	//--------------------------------------------------------------------------
	char h5filename[50];    // HDF5 filename to store 3D data
	char groupname[100];    // HDF5 group name to group data
	char subgroupname[100]; // HDF5 subgroup name to group data
	char fieldname[100];    // HDF5 name

	int covar_calc = 0;
	int avg_dir = 1;  // 0=x, 1=y
	int histo_calc = 0;
	int U_mean_calc = 1;
	int vorticity_calc = 0;
	int U_fft_calc = 0;

	double camp = 1.0;
	double cmin = 0.0;
	double eps = 0.0000001;

	#ifdef VOF_PP_VFPRIME

	if (U_mean_calc == 1){


		double U_mean_vof_bl = 0.0;
		double vf_mean = 0.0;



		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					vf = min(vfu_p[k][j][i], 1.0);
					U_mean_vof_bl += (1.0 - vf) * u_data[k][j][i];
					vf_mean += (1.0 - vf);
				}
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, &U_mean_vof_bl, 1, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, &vf_mean, 1, MPI_DOUBLE, MPI_SUM, PCW);
		U_mean_vof_bl = U_mean_vof_bl / vf_mean;

		FILE *fileout;
		char filename[100];
		sprintf(filename, "U_mean_vof_bl.dat");

		if (params->rank == 0) {
			/* append or create new file */
			if (params -> time == 0.0)
				fileout = fopen(filename, "w");
			else
				fileout = fopen(filename, "a");
			Output_check_file_open(fileout);

			fprintf(fileout, "%f %f\n",params -> time, U_mean_vof_bl);

			fclose(fileout);

			}

	}

/*
	sprintf(h5filename, "vfu_prime%d.h5", abs(params->noutput));
	file_id = Output_h5_open(h5filename, params, DTRACE("Output_h5_open"));


	sprintf(fieldname, "/time");
	dim_1d[0] = 1;
	Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));


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


	Output_h5_flow_variable(vfu_p, file_id, "/vfu_prime", grid, params, DTRACE("Output_h5_flow_variable"));
	assert(H5Fclose(file_id) >= 0);
*/

	if (U_fft_calc == 1){





	}
	#endif

	/*

		COVARIANCE CALCULATION

	*/

	if (covar_calc == 1){

		int l;
		int Navg;
		int Ncells;
		if (avg_dir == 0){
			Navg = NX-1;
			Ncells = (NY-1) * (NZ-1);
		}
		else{
			Navg = NY-1;
			Ncells = (NX-1) * (NZ-1);
		}
		// Average v-velocity at each y-coordinate
		double *avg_u = Memory_allocate_1D_array(GVG_DOUBLE, Navg);
		double *avg_v = Memory_allocate_1D_array(GVG_DOUBLE, Navg);
		double *avg_w = Memory_allocate_1D_array(GVG_DOUBLE, Navg);

		// Average conc at each y-coordinate
		double *avg_c = Memory_allocate_1D_array(GVG_DOUBLE, Navg);
		double *avg_mvf = Memory_allocate_1D_array(GVG_DOUBLE, Navg);

		// Average uc
		double *avg_uc = Memory_allocate_1D_array(GVG_DOUBLE, Navg);
		double *avg_vc = Memory_allocate_1D_array(GVG_DOUBLE, Navg);
		double *avg_wc = Memory_allocate_1D_array(GVG_DOUBLE, Navg);

		DSET_ZERO(avg_mvf, Navg);
		DSET_ZERO(avg_c, Navg);
		DSET_ZERO(avg_u, Navg);
		DSET_ZERO(avg_v, Navg);
		DSET_ZERO(avg_w, Navg);
		DSET_ZERO(avg_uc, Navg);
		DSET_ZERO(avg_vc, Navg);
		DSET_ZERO(avg_wc, Navg);

		c_data = c[0]->data;


		// Calculate local fluid flux and volume fraction data

			for (k = k_start; k < k_end; k++) {
				for (j = j_start; j < j_end; j++) {
					for (i = i_start; i < i_end; i++) {

						if(avg_dir==0){
							l=i;
						}
						else{
							l=j;
						}

						#ifdef LAG_PARTICLE_RESOLVED
						vf = min(vfc[k][j][i], 1.0);
						avg_c[l] += (1.0 - vf) * c_data[k][j][i];
						avg_u[l] += (1.0 - vf) * u_data[k][j][i];
						avg_v[l] += (1.0 - vf) * v_data[k][j][i];
						avg_w[l] += (1.0 - vf) * w_data[k][j][i];
						avg_mvf[l] +=  (1.0 - vf);
						#else
						avg_c[l] += c_data[k][j][i];
						avg_u[l] += u_data[k][j][i];
						avg_v[l] += v_data[k][j][i];
						avg_w[l] += w_data[k][j][i];
						avg_mvf[l] += 1.0;
						#endif
					}
				}
			}

		MPI_Allreduce(MPI_IN_PLACE, avg_mvf, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_c, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_u, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_v, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_w, Navg, MPI_DOUBLE, MPI_SUM, PCW);

		for (j = 0; j < Navg; j++) {
			avg_mvf[j] = avg_mvf[j]/Ncells;
			avg_c[j] = avg_c[j]/avg_mvf[j]/Ncells;
			avg_u[j] = avg_u[j]/avg_mvf[j]/Ncells;
			avg_v[j] = avg_v[j]/avg_mvf[j]/Ncells;
			avg_w[j] = avg_w[j]/avg_mvf[j]/Ncells;
		}


		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {

					if(avg_dir==0){
						l=i;
					}
					else{
						l=j;
					}

					#ifdef LAG_PARTICLE_RESOLVED
					vf = min(vfc[k][j][i], 1.0);
					cp = (1.0-vf)*(c_data[k][j][i] -  avg_c[l]);
					up = (1.0-vf)*(u_data[k][j][i] -  avg_u[l]);
					vp = (1.0-vf)*(v_data[k][j][i] -  avg_v[l]);
					wp = (1.0-vf)*(w_data[k][j][i] -  avg_w[l]);
					#else
					cp = c_data[k][j][i] - avg_c[l];
					up = u_data[k][j][i] - avg_u[l];
					vp = v_data[k][j][i] - avg_v[l];
					wp = w_data[k][j][i] - avg_w[l];
					#endif
					avg_uc[l] += up*cp;
					avg_vc[l] += vp*cp;
					avg_wc[l] += wp*cp;
				}
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, avg_uc, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_vc, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_wc, Navg, MPI_DOUBLE, MPI_SUM, PCW);

		for (j = 0; j < Navg; j++) {
			avg_uc[j] = avg_uc[j]/avg_mvf[j]/Ncells;
			avg_vc[j] = avg_vc[j]/avg_mvf[j]/Ncells;
			avg_wc[j] = avg_wc[j]/avg_mvf[j]/Ncells;
		}

		/*------------------------------------------------------------------------*/
		/*
	 	* HDF5 output
	 	*/
		/*------------------------------------------------------------------------*/


		//--------------------------------------------------------------------------
		// Open HDF5 file handle
		//--------------------------------------------------------------------------
		sprintf(h5filename, "Data_uc_%d.h5", abs(params->noutput));
		file_id = Output_h5_open(h5filename, params, DTRACE("Output_h5_open"));

		//--------------------------------------------------------------------------
		// Write timestamp to *.h5 file
		//--------------------------------------------------------------------------
		sprintf(fieldname, "/time");
		dim_1d[0] = 1;
		Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

		dim_1d[0] = Navg;

		if(avg_dir==0){
			sprintf(fieldname, "/x");
			Output_h5_dataset(grid->xc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		}
		else{
			sprintf(fieldname, "/y");
			Output_h5_dataset(grid->yc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		}


		sprintf(groupname, "/covar");
		Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));
		sprintf(fieldname, "%s/uc", groupname);
		Output_h5_dataset(avg_uc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "%s/vc", groupname);
		Output_h5_dataset(avg_vc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "%s/wc", groupname);
		Output_h5_dataset(avg_wc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(groupname, "/avg");
		Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));
		sprintf(fieldname, "%s/c", groupname);
		Output_h5_dataset(avg_c, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "%s/u", groupname);
		Output_h5_dataset(avg_u, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "%s/v", groupname);
		Output_h5_dataset(avg_v, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "%s/w", groupname);
		Output_h5_dataset(avg_w, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "%s/mvf", groupname);
		Output_h5_dataset(avg_mvf, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
		//--------------------------------------------------------------------------
		// Close HDF5 file handle
		//--------------------------------------------------------------------------
		status = H5Fclose(file_id);
		Display_assert_error(status, "Close file_id failed", params, DTRACE("Display_assert_error"));


		free(avg_uc);
		free(avg_vc);
		free(avg_wc);
		free(avg_c);
		free(avg_u);
		free(avg_v);
		free(avg_w);

	}

	if (histo_calc == 1){

		int Nbin = 1000;
		FILE *f;
		int ibin = 0;
		int NTOT = 0;
		double cbd;

		double *histo = (double *)malloc(Nbin * sizeof(double));
		double *range = (double *)malloc(Nbin * sizeof(double));

		for (iconc = 0; iconc < NConc; iconc++) {
			c_data = c[iconc]->data;

			for (ibin=0; ibin<Nbin; ibin++) {
				range[ibin] = cmin + camp*((double) ibin) / Nbin;
				histo[ibin] = 0.0;
			}

			for (k=k_start; k<k_end; k++) {
				for (j=j_start; j<j_end; j++) {
					for (i=i_start; i<i_end; i++) {
						NTOT++;
						cbd = (c_data[k][j][i] - cmin)/camp;

						#ifdef LAG_PARTICLE_RESOLVED
						vf = min(vfc[k][j][i], 1.0);
						//cbd = min(c_data[k][j][i],0.99999);
						//cbd = max(cbd,0.00001);
						cbd = (1-vf)*cbd;
						#endif

						cbd = min(cbd, 1.0-eps);
						cbd = max(cbd, eps);
						ibin = (int) (cbd * Nbin);
						if ( (ibin>=0) && (ibin<Nbin) ) histo[ibin]++;
					}
				}
			}

			MPI_Allreduce(MPI_IN_PLACE, &NTOT, 1, MPI_INT, MPI_SUM, PCW);


			for (ibin=0; ibin<Nbin; ibin++) {
				histo[ibin]*= 1.0/NTOT;
			}

			MPI_Allreduce(MPI_IN_PLACE, &histo[0], Nbin, MPI_DOUBLE, MPI_SUM, PCW);

			sprintf(h5filename, "histo_conc%d_%d.h5", iconc, abs(params->noutput));
			file_id = Output_h5_open(h5filename, params, DTRACE("Output_h5_open"));

			//--------------------------------------------------------------------------
			// Write timestamp to *.h5 file
			//--------------------------------------------------------------------------
			sprintf(fieldname, "/time");
			dim_1d[0] = 1;
			Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

			sprintf(groupname, "/grid");
			Output_h5_create_group(file_id, groupname, params, DTRACE("Output_h5_create_group"));

			dim_1d[0] = NX-1;

			sprintf(fieldname, "%s/x", groupname);
			Output_h5_dataset(grid->xc, H5T_NATIVE_DOUBLE, 1, dim_1d , file_id, fieldname, params, DTRACE("Output_h5_dataset"));

			dim_1d[0] = NY-1;

			sprintf(fieldname, "%s/y", groupname);
			Output_h5_dataset(grid->yc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

			dim_1d[0] = NZ-1;

			sprintf(fieldname, "%s/z", groupname);
			Output_h5_dataset(grid->zc, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

			dim_1d[0] = Nbin;

			sprintf(fieldname, "/histo");
			Output_h5_dataset(histo, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));
			sprintf(fieldname, "/range");
			Output_h5_dataset(range, H5T_NATIVE_DOUBLE, 1, dim_1d, file_id, fieldname, params, DTRACE("Output_h5_dataset"));

			status = H5Fclose(file_id);
			Display_assert_error(status, "Close file_id failed", params, DTRACE("Display_assert_error"));

		}

	free(histo);
	free(range);

	}


}
