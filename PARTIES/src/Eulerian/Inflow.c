#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Display.h"
#include "Inflow.h"
#include "Velocity.h"
#include "Conc.h"
#include "Grid.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/******************************************************************************/
/*
 This function defines the inflow profile
 */
/******************************************************************************/
void Inflow_velocity_profile(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i, j, k;
	double dy, dz, dA;
	double influx, W_influx;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// influx on the current processor
	influx = 0.0;

	// Success in reading inflow data from file
	int statusu = 0;
	int statusv = 0;
	int statusw = 0;
	char message[100];

    // Check if we are running the prescribed advection test.
    if (params->vel_init_type == VEL_INIT_ADVECTION_TEST) {
        // Compute domain lengths. (Assumes these are set; otherwise use xmax-xmin, etc.)
        double Lx = params->xmax - params->xmin;
        double Ly = params->ymax - params->ymin;

        // Decide on the phase: for a time reversal test, reverse after advection_test_time.
        double phase = (params->time < params->advection_test_time) ? 1.0 : -1.0;

        // Loop over the processor’s domain (or the entire domain if appropriate).
        for (int k = grid->G_Ks; k < grid->G_Ke; k++) {
            for (int j = grid->G_Js; j < grid->G_Je; j++) {
                for (int i = grid->G_Is; i < grid->G_Ie; i++) {
                    // For a staggered grid, you might use cell-center coordinates for one velocity
                    // and staggered positions for the other. Here we assume that grid->xc and grid->yc
                    // store the cell-center positions, and grid->xu and grid->yv store the positions
                    // where u and v are defined, respectively.
                    double x_u = grid->xu[i];  // u is defined at xu[]
                    double y_c = grid->yc[j];  // use cell-center for y derivative
                    double x_c = grid->xc[i];  // use cell-center for x derivative
                    double y_v = grid->yv[j];  // v is defined at yv[]

                    // Compute the stream function derivatives for u and v.
                    // For u = -dψ/dy at the u–position:
                    double sin_x_u = sin(PI * x_u / Lx);
                    double sin_y_c = sin(PI * y_c / Ly);
                    double cos_y_c = cos(PI * y_c / Ly);
                    // Derivative with respect to y:
                    double dpsi_dy = 2.0 * sin_y_c * cos_y_c * (PI / Ly) * (sin_x_u * sin_x_u);
                    double u_val = -dpsi_dy;

                    // For v = dψ/dx at the v–position:
                    double sin_x_c = sin(PI * x_c / Lx);
                    double cos_x_c = cos(PI * x_c / Lx);
                    double sin_y_v = sin(PI * y_v / Ly);
                    // Derivative with respect to x:
                    double dpsi_dx = 2.0 * sin_x_c * cos_x_c * (PI / Lx) * (sin_y_v * sin_y_v);
                    double v_val = dpsi_dx;

                    // Set the prescribed velocities (apply phase for reversal)
                    u_data[k][j][i] = phase * u_val;
                    v_data[k][j][i] = phase * v_val;
                    w_data[k][j][i] = 0.0;
                }
            }
        }
        // Update ghost nodes to ensure consistency (NOT NECESSARY I THINK)
        Communication_update_ghost_nodes_flow_variable(u_data, 'u', params->ghost_nodes, data_bag);
        Communication_update_ghost_nodes_flow_variable(v_data, 'v', params->ghost_nodes, data_bag);
        Communication_update_ghost_nodes_flow_variable(w_data, 'w', params->ghost_nodes, data_bag);
        
        return; // Exit the function
    }

	




	//--------------------------------------------------------------------------
	// Inflow at left wall
	//--------------------------------------------------------------------------
#ifdef LEFT_INFLOW
	// Update only on the processors which have the first yz plane
	if (grid->G_Is == 0) {

		double *yc = grid->yc;
		double *zc = grid->zc;
		double Ly  = params->Ly;
		double Lz  = params->Lz;
		double Re  = params->Re;

		// start index on current processor
		int j_start = grid->G_Js;
		int k_start = grid->G_Ks;

		// end index on current processor
		int j_end = min(NY-1, grid->G_Je);
		int k_end = min(NZ-1, grid->G_Ke);

		double  **u_inflow_n, **v_inflow_n, **w_inflow_n;
		double  **u_inflow_o, **v_inflow_o, **w_inflow_o;
		double  **u_inflow,   **v_inflow,   **w_inflow;

		if (params -> vel_init_type == VEL_INIT_PRECURSOR) {
			int read_new_data = params -> read_new_data;

			double time   = params->time;

			u_inflow_n = data_bag -> u -> inflow_n;
			v_inflow_n = data_bag -> v -> inflow_n;
			w_inflow_n = data_bag -> w -> inflow_n;

			u_inflow_o = data_bag -> u -> inflow_o;
			v_inflow_o = data_bag -> v -> inflow_o;
			w_inflow_o = data_bag -> w -> inflow_o;

			u_inflow   = data_bag -> u -> inflow;
			v_inflow   = data_bag -> v -> inflow;
			w_inflow   = data_bag -> w -> inflow;

			// Update inflow data
			read_new_data  = 0;
			if (time > params->t_in_n){
				params->ninflow++;
				params->t_in_o = params->t_in_n;
				read_new_data  = 1;
			}
			statusu = Get_inflow_data(u_inflow, u_inflow_n, u_inflow_o,'u', read_new_data, data_bag);
			statusv = Get_inflow_data(v_inflow, v_inflow_n, v_inflow_o,'v', read_new_data, data_bag);
			statusw = Get_inflow_data(w_inflow, w_inflow_n, w_inflow_o,'w', read_new_data, data_bag);
			params->read_new_data = read_new_data;
		}

		if (statusu >= 0 && statusv >= 0 && statusw >= 0) {
			//------------------------------------------------------------------
			// Set velocity profile at the inlet
			//------------------------------------------------------------------
			i = 0;
			for (k = k_start; k < k_end; k++) {
				for (j = j_start; j < j_end; j++) {
					if (grid->u_status[k][j][i+1] == FLUID) {

						if (params -> vel_init_type == VEL_INIT_UNIFORM) {
							u_data[k][j][i] = params -> ubulk_target;
							v_data[k][j][i-1] = 0.0;
							w_data[k][j][i-1] = 0.0;
						}
						else if (params -> vel_init_type == VEL_INIT_LINEAR) {
							u_data[k][j][i] = 2.0 * params -> ubulk_target * yc[j];
							v_data[k][j][i-1] = 0.0;
							w_data[k][j][i-1] = 0.0;
						}
						else if (params -> vel_init_type == VEL_INIT_POISEUILLE) {
							u_data[k][j][i] = -6.0 * params -> ubulk_target * yc[j] * (yc[j] - 1.0);
							v_data[k][j][i-1] = 0.0;
							w_data[k][j][i-1] = 0.0;
						}
						else if (params -> vel_init_type == VEL_INIT_PRECURSOR) {
							u_data[k][j][i]   =  u_inflow[k][j];
							v_data[k][j][i-1] =  v_inflow[k][j];
							w_data[k][j][i-1] =  w_inflow[k][j];
						}

						else if (params -> vel_init_type == VEL_INIT_JET) {
							double y1 = Ly/2.0 - 0.5;
							double y2 = Ly/2.0 + 0.5;
							double delta = params->delta_jet;
							u_data[k][j][i]   =  0.5*(erf((yc[j]-y1)/delta)-erf((yc[j]-y2)/delta));
							v_data[k][j][i-1] =  0.0;
							w_data[k][j][i-1] = 0.0;
						}

						else if (params -> vel_init_type == VEL_INIT_JET_3D) {
							double r1 = 0.5;
							double ypos = yc[j] - Ly/2.0;
							double zpos = zc[k] - Lz/2.0;

							double rpos = sqrt(ypos*ypos + zpos*zpos);
							double delta = params->delta_jet;
							u_data[k][j][i]   =  0.5*(1-erf((rpos-r1)/delta));
							v_data[k][j][i-1] =  0.0;
							w_data[k][j][i-1] = 0.0;
						}
						else if (params -> vel_init_type == VEL_INIT_LEFT_RIGHT) {
							u_data[k][j][i] = params -> u_in_left;
							v_data[k][j][i-1] = 0.0;
							w_data[k][j][i-1] = 0.0;
						}
	//					// Rotational shear flow
	//					u_data[k][j][0] = 1.0 / params -> G + yc[j] - (0.5 * (yc[grid->NY-2] - yc[0]) + yc[0]);

					} // if

				} // for j
			} // for k



			//------------------------------------------------------------------
			// Calculate the mass flux at the inlet, i.e. u*Area
			//------------------------------------------------------------------
			j_end = min(j_end, grid->NY-1);
			k_end = min(k_end, grid->NZ-1);

			for (k = k_start; k < k_end; k++) {
				for (j = j_start; j < j_end; j++) {
					if (grid->u_status[k][j][1] == FLUID) {

						dy = grid->dy_v[j];
						dz = grid->dz_w[k];
						dA = dy * dz;

						// Calculate the mass flux at the inlet, i.e. u*Area
						influx += u_data[k][j][0] * dA;

					} // if
				} // for j
			} // for k
		} // if successful read
	}  // if

	sprintf(message, "Could not read inflow u data");
	Display_assert_error(statusu, message, params, DTRACE("Display_assert_error"));
	sprintf(message, "Could not read inflow v data");
	Display_assert_error(statusv, message, params, DTRACE("Display_assert_error"));
	sprintf(message, "Could not read inflow w data");
	Display_assert_error(statusw, message, params, DTRACE("Display_assert_error"));

#endif
	//--------------------------------------------------------------------------
	// Inflow at right wall
	//--------------------------------------------------------------------------
#if defined RIGHT_INFLOW

	// Update only on the processors which have the first yz plane
	if (grid->G_Ie == NX) {

		double *yc = grid->yc;
		double Ly  = params->Ly;
		double Re  = params->Re;

		// start index on current processor
		int j_start = grid->G_Js;
		int k_start = grid->G_Ks;

		// end index on current processor
		int j_end = grid->G_Je;
		int k_end = grid->G_Ke;

		//----------------------------------------------------------------------
		// Set velocity profile at the inlet
		//----------------------------------------------------------------------
		i = NX-1;

		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				if (grid->u_status[k][j][i-1] == FLUID) {

					if (params -> vel_init_type == VEL_INIT_UNIFORM) {
						u_data[k][j][i] = params -> ubulk_target;
					}
					else if (params -> vel_init_type == VEL_INIT_LINEAR) {
						u_data[k][j][i] = 2.0 * params -> ubulk_target * yc[j];
					}
					else if (params -> vel_init_type == VEL_INIT_POISEUILLE) {
						u_data[k][j][i] = -6.0 * params -> ubulk_target * yc[j] * (yc[j] - 1.0);
					}
					else if (params -> vel_init_type == VEL_INIT_POISEUILLE) {
						u_data[k][j][i] = -6.0 * params -> ubulk_target * yc[j] * (yc[j] - 1.0);
					}
					else if (params -> vel_init_type == VEL_INIT_LEFT_RIGHT) {
						u_data[k][j][i] = - params -> u_in_right;
					}
/*
					// Rotational shear flow
					u_data[k][j][0] = 1.0 / params -> G + yc[j] - (0.5 * (yc[grid->NY-2] - yc[0]) + yc[0]);
*/
				} // if

				v_data[k][j][i] = 0.0;
				w_data[k][j][i] = 0.0;

			} // for j
		} // for k



		//----------------------------------------------------------------------
		// Calculate the mass flux at the inlet, i.e. u*Area
		//----------------------------------------------------------------------
		j_end = min(j_end, grid->NY-1);
		k_end = min(k_end, grid->NZ-1);

		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				if (grid->u_status[k][j][i-1] == FLUID) {

					dy = grid->dy_v[j];
					dz = grid->dz_w[k];
					dA = dy * dz;

					// Calculate the mass flux at the inlet, i.e. u*Area
					influx += u_data[k][j][i] * dA;

				} // if
			} // for j
		} // for k
	}  // if
#endif


	// Now, send the local influx part and sum them up to get the total influx:
	// W_influx
	MPI_Allreduce ( (void *)&influx, (void *)&W_influx, 1, MPI_DOUBLE, MPI_SUM, PCW);

	// Local part of influx
	data_bag -> u -> G_influx = influx;

	// Total influx (send to all processors)
	data_bag -> u -> W_influx = W_influx;

}




/******************************************************************************/
/*
 This function reads new precursor data for the inflow as time progresses

 Returns 0 on success
 Returns -1 if it could not open the file
 Returns -2 if it had a problem reading data from the file
 */
/******************************************************************************/
int Get_inflow_data(double **inflow, double **inflow_n, double **inflow_o,
		char component, int read_new_data, Cart3d_bag *data_bag) {

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	FILE *fp ;
	double *yc    = grid -> yc;
	double time   = params->time;

	// Start index of bottom-left-back corner on current processor
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid -> G_Ie;
	int Je = grid -> G_Je;
	int Ke = grid -> G_Ke;

	int status;
	int  i, j, k, jj;
	char buff[255];
	char filename[50];
	float  y_in[301], vel_in[301];
	float  t_in_n = params->t_in_n;
	float  t_in_o = params->t_in_o;
	double dvel;

	// Get new data if needed
	if (read_new_data){

		// Shift inflow vector in time
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				inflow_o[k][j] = inflow_n[k][j];
				inflow_o[k][j] = inflow_n[k][j];
				inflow_o[k][j] = inflow_n[k][j];
			}
		}

		// Open precursor file
		if      (component == 'u') {
			sprintf(filename, "./inflow/u_in_%04d.dat", params->ninflow);
		}
		else if (component == 'v') {
			sprintf(filename, "./inflow/v_in_%04d.dat", params->ninflow);
		}
		else if (component == 'w') {
			sprintf(filename, "./inflow/w_in_%04d.dat", params->ninflow);
		}
		else if (component == 'c') {
			sprintf(filename, "./inflow/c_in_%04d.dat", params->ninflow);
		}
//		if(component == 'c') printf("concentration %d %s\n",  params->ninflow, filename);
		// Read data
		fp = fopen(filename, "r");
		if (fp == NULL) {
			return -1;
		}
		status = fscanf(fp, "%s", buff);
		if (status > 0) status = fscanf(fp, "%g", &t_in_n);
		if (status > 0) status = fscanf(fp, "%s", buff);
		if (status > 0) status = fscanf(fp, "%s", buff);

		for (j = 0; j<301; j++){
			if (status > 0) status = fscanf(fp, "%g", &y_in[j]);
			if (status > 0) status = fscanf(fp, "%g", &vel_in[j]);
		}
		fclose(fp);
		params ->t_in_n = t_in_n;

		if (status < 1) {
			return -2;
		}

		// Interpolate to present grid (only y-direction so far!)
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (jj = 0; jj<300; jj++){
					if (yc[j] > y_in[jj] && yc[j]< y_in[jj+1]){
						dvel = (vel_in[jj+1]-vel_in[jj])*(yc[j] - y_in[jj])/(y_in[jj+1] - y_in[jj]);
						inflow_n[k][j] = vel_in[jj] + dvel;
					}
				}
			}
		}
	}


	// Interpolate in time
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			dvel = (inflow_n[k][j]-inflow_o[k][j])*(time - t_in_o)/(t_in_n - t_in_o);
			inflow[k][j] = inflow_o[k][j] + dvel;
		}
	}

	return 0;
}
