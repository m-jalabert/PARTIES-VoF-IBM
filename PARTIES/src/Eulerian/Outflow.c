#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Outflow.h"
#include "Velocity.h"
#include "Conc.h"
#include "Grid.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/******************************************************************************/
/* Outflow boundary condition using the convective condition.

  du_i/dt + U_c ( du_i / dx ) = 0
  dc_i/dt + U_c ( dc_i / dx ) = 0

  where U_c is the convection velocity (set to either the mean or maximum
  velocity at the exit plane). The mass flux is adjusted to keep a constant
  mass flux in the inflow and exit planes.

  ref: I. Orlanski, "A simple boundary condition for unbounded hyperbolic
  flows," J. Comput. Phys. 21, 251 (1976)
*/

/******************************************************************************/
void Outflow_impose_convective_boundary(Cart3d_bag *data_bag) {

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;

	Outflow_vel_impose_convective_boundary(u, v, w, grid, params);

}




/******************************************************************************/
/*
 This function imposes the convective boundary condition at the outflow.
 */
/******************************************************************************/
void Outflow_vel_impose_convective_boundary( Velocity *u, Velocity *v,
		Velocity *w, MAC_grid *grid, Parameters *params) {

	int Js, Ks;
	int Je, Ke;
	int i, j, k;

	double  **du_dn_l_n, **dv_dn_l_n,  **dw_dn_l_n;
	double  **du_dn_l_o, **dv_dn_l_o,  **dw_dn_l_o;

	double  **du_dn_r_n, **dv_dn_r_n,  **dw_dn_r_n;
	double  **du_dn_r_o, **dv_dn_r_o,  **dw_dn_r_o;

	double coef, idx_u, idx_c;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	const double BET[] = {BETA};
	const double GAM[] = {GAMMA};
	const double ZET[] = {ZETA};

	double bet = BET[params -> which_stage];
	double gam = GAM[params -> which_stage];
	double zet = ZET[params -> which_stage];
	double dt  = params -> dt;

	double ***u_data = u -> data;
	double ***v_data = v -> data;
	double ***w_data = w -> data;

	// convective velocity
	double U  = params -> U_conv_outflow;

#ifdef LEFT_OUTFLOW
	// Update only on the processors which have the exit plane (first yz plane)
	if (grid->G_Is == 0) {
		i = 0;

		du_dn_l_n = u -> d_dn_l_n;
		dv_dn_l_n = v -> d_dn_l_n;
		dw_dn_l_n = w -> d_dn_l_n;

		du_dn_l_o = u -> d_dn_l_o;
		dv_dn_l_o = v -> d_dn_l_o;
		dw_dn_l_o = w -> d_dn_l_o;

		idx_u = grid -> idx_u[i+1];
		idx_c = grid -> idx_c[i+1];

		// start index on current processor
		Js = grid -> G_Js;
		Ks = grid -> G_Ks;

		// end index on current processor
		Je = min(grid -> G_Je, NY-1);
		Ke = min(grid -> G_Ke, NZ-1);

		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {

				du_dn_l_n[k][j] = (u_data[k][j][i  ] - u_data[k][j][i+1]) * idx_u;
				dv_dn_l_n[k][j] = (v_data[k][j][i-1] - v_data[k][j][i  ]) * idx_c;
				dw_dn_l_n[k][j] = (w_data[k][j][i-1] - w_data[k][j][i  ]) * idx_c;

				u_data[k][j][i  ] = u_data[k][j][i  ] - U * dt * (gam * du_dn_l_n[k][j] + zet * du_dn_l_o[k][j]);
				v_data[k][j][i-1] = v_data[k][j][i-1] - U * dt * (gam * dv_dn_l_n[k][j] + zet * dv_dn_l_o[k][j]);
				w_data[k][j][i-1] = w_data[k][j][i-1] - U * dt * (gam * dw_dn_l_n[k][j] + zet * dw_dn_l_o[k][j]);

				du_dn_l_o[k][j] = du_dn_l_n[k][j];
				dv_dn_l_o[k][j] = dv_dn_l_n[k][j];
				dw_dn_l_o[k][j] = dw_dn_l_n[k][j];

			} // for j
		} // for k


	} // if
#endif
#ifdef RIGHT_OUTFLOW
	// Update only on the processors which have the exit plane (last yz plane)
	if (grid->G_Ie == grid->NX) {
		i = NX-1;

		du_dn_r_n = u -> d_dn_r_n;
		dv_dn_r_n = v -> d_dn_r_n;
		dw_dn_r_n = w -> d_dn_r_n;

		du_dn_r_o = u -> d_dn_r_o;
		dv_dn_r_o = v -> d_dn_r_o;
		dw_dn_r_o = w -> d_dn_r_o;

		idx_u = grid -> idx_u[NX-2];
		idx_c = grid -> idx_c[NX-2];
		// start index on current processor for u
		Js = grid -> G_Js;
		Ks = grid -> G_Ks;

		// end index on current processor
		Je = min(grid -> G_Je, NY-1);
		Ke = min(grid -> G_Ke, NZ-1);

		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {

				du_dn_r_n[k][j] = (u_data[k][j][i] - u_data[k][j][i-1]) * idx_u;
				dv_dn_r_n[k][j] = (v_data[k][j][i] - v_data[k][j][i-1]) * idx_c;
				dw_dn_r_n[k][j] = (w_data[k][j][i] - w_data[k][j][i-1]) * idx_c;

				u_data[k][j][i] = u_data[k][j][i] - U * dt * (gam * du_dn_r_n[k][j] + zet * du_dn_r_o[k][j]);
				v_data[k][j][i] = v_data[k][j][i] - U * dt * (gam * dv_dn_r_n[k][j] + zet * dv_dn_r_o[k][j]);
				w_data[k][j][i] = w_data[k][j][i] - U * dt * (gam * dw_dn_r_n[k][j] + zet * dw_dn_r_o[k][j]);

				du_dn_r_o[k][j] = du_dn_r_n[k][j];
				dv_dn_r_o[k][j] = dv_dn_r_n[k][j];
				dw_dn_r_o[k][j] = dw_dn_r_n[k][j];

			} // for j
		} // for k


	} // if
#endif
	// Now, update u to conserve global mass flux.
	Outflow_update_u_velocity_to_conserve_mass(u, grid, params) ;

}




/******************************************************************************/
/*
 This function adjusts the mass flux at the exit plane so that it
 matches the mass flux at the inlet
 */
/******************************************************************************/
void Outflow_update_u_velocity_to_conserve_mass(Velocity *u, MAC_grid *grid,
		Parameters *params) {

	int i, j, k;
	int Js, Ks;
	int Je, Ke;
	int j_end, k_end;
	int NX;
	double total_u_outflux;
	double outflow_area;
	double W_outflow_area;
	double W_total_u_outflux;
	double W_total_u_influx;
	double after_u_flux;
	double velocity_correction;
	double dy, dz, dA;
	double send_data[2]={0,0}, W_sum_data[2]={0,0};
	double ***u_data;
	int ***status;

	u_data = u -> data;

	// Update only on the processors which have the last yz plane
	// First, integrate over the exit plane to find the u_outflux
	total_u_outflux = 0.0;
	outflow_area    = 0.0;
#if defined LEFT_OUTFLOW && defined RIGHT_OUTFLOW
	double influx, W_influx;
	influx = 0.0;
#endif


	status = grid -> u_status;

	// start index on current processor
	Js = grid -> G_Js;
	Ks = grid -> G_Ks;

	// end index on current processor
	Je = grid -> G_Je;
	Ke = grid -> G_Ke;
	j_end = min(Je, grid->NY-1);
	k_end = min(Ke, grid->NZ-1);

	NX = grid -> NX;
#ifdef LEFT_OUTFLOW
	if (grid->G_Is == 0) {
		for (k = Ks; k < k_end; k++) {
			for (j = Js; j < j_end; j++) {

				// Check if the node to the left of the outflow node is fluid
//				if (status[k][j][NX-1]  == FLUID) {

					dy = grid -> dy_v[j];
					dz = grid -> dz_w[k];
					dA = dy * dz;
#ifdef RIGHT_OUTFLOW
					// This is done if there are outflow-BC on the left and right
					// Left serves as inflow and right as outflow to conserve mass
					// Calculate the mass flux at the 'inlet', i.e. u*Area
					influx += u_data[k][j][0] * dA;
#else
					// Total outflow area. This is done if there is a geometry
					// (solid) right before the outflow boundary
//Note move this to parameter
					outflow_area += dA;

					// integral over each mesh: u*(dy*dz)
					total_u_outflux += dA * u_data[k][j][0];
#endif
//				} // if
			} // for j
		} // for k
	} // if
#endif
#ifdef RIGHT_OUTFLOW
	if (grid->G_Ie == NX) {

		for (k = Ks; k < k_end; k++) {
			for (j = Js; j < j_end; j++) {

				// Check if the node to the left of the outflow node is fluid
//				if (status[k][j][NX-1]  == FLUID) {

					dy = grid -> dy_v[j];
					dz = grid -> dz_w[k];
					dA = dy * dz;

					// Total outflow area. This is done if there is a geometry
					// (solid) right before the outflow boundary
//Note move this to parameter
					outflow_area += dA;

					// integral over each mesh: u*(dy*dz)
					total_u_outflux += dA * u_data[k][j][NX-1];

//				} // if
			} // for j
		} // for k

	} // if
#endif
#if defined LEFT_OUTFLOW && defined RIGHT_OUTFLOW
	// Now, send the local influx part and sum them up to get the total influx:
	// W_influx
	MPI_Allreduce ( (void *)&influx, (void *)&W_influx, 1, MPI_DOUBLE, MPI_SUM, PCW);

	// Local part of influx
	u -> G_influx = influx;

	// Total influx (send to all processors)
	u -> W_influx = W_influx;
#endif
	// Send two variables to all processors, sum them up and send them back in
	// W_recv_data array
	send_data[0] = total_u_outflux;
	send_data[1] = outflow_area;

	// Now, add all the total_u_outflux and outflow area from all processors
	// and get the sum back on all processors.
	MPI_Allreduce ( (void *)send_data, (void *)W_sum_data, 2, MPI_DOUBLE, MPI_SUM, PCW);

	// Update only on the processors which have the last yz plane
	// Now, add or rescale the velocity at the exit plane to balance the mass-flux
#if defined LEFT_OUTFLOW && !defined RIGHT_OUTFLOW
	// If left and right are defined as OUTFLOW, only right is corrected
	if (grid->G_Is == 0) {

		// (World) Total influx, outflux and outflow area on all processors
		W_total_u_influx  = u -> W_influx;
		W_total_u_outflux = W_sum_data[0];
		W_outflow_area    = W_sum_data[1];

		// Now, calculate the correction amount to be added (substracted) to
		// all the u-nodes at the outflow region
		velocity_correction = (W_total_u_influx - W_total_u_outflux) / W_outflow_area;
//		if (abs(W_total_u_influx) > 1e-3) rescale W_total_u_influx/W_total_u_outflux;
/*					printf("!\n");
					printf("!\n");
					printf("in_out  %g %g \n", W_total_u_influx, W_total_u_outflux);
					printf("!\n");
					printf("!\n");	*/
		total_u_outflux = 0.0;
		for (k = Ks; k < k_end; k++) {
			for (j = Js; j < j_end; j++) {

				// Check if the node to the left of the outflow node is fluid
//				if (status[k][j][NX-1]  == FLUID) {

					dy = grid -> dy_v[j];
					dz = grid -> dz_w[k];
					dA = dy * dz;

					u_data[k][j][0] += velocity_correction;
					total_u_outflux += dA * u_data[k][j][0];

//				} // if
			} // for j
		} // for k

	} // if
#endif
#ifdef RIGHT_OUTFLOW
	if (grid->G_Ie == grid->NX) {

		// (World) Total influx, outflux and outflow area on all processors
		W_total_u_influx  = u -> W_influx;
		W_total_u_outflux = W_sum_data[0];
		W_outflow_area    = W_sum_data[1];

		// Now, calculate the correction amount to be added (substracted) to
		// all the u-nodes at the outflow region
		velocity_correction = (W_total_u_influx - W_total_u_outflux) / W_outflow_area;
//		if (abs(W_total_u_influx) > 1e-3) rescale W_total_u_influx/W_total_u_outflux;

		total_u_outflux = 0.0;
		for (k = Ks; k < k_end; k++) {
			for (j = Js; j < j_end; j++) {

				// Check if the node to the left of the outflow node is fluid
//				if (status[k][j][NX-1]  == FLUID) {

					dy = grid -> dy_v[j];
					dz = grid -> dz_w[k];
					dA = dy * dz;

					u_data[k][j][NX-1] += velocity_correction;
					total_u_outflux += dA * u_data[k][j][NX-1];

//				} // if
			} // for j
		} // for k

	} // if
#endif

//	send_data[0] = total_u_outflux;
//	MPI_Allreduce ( (void *)send_data, (void *)W_sum_data, 2, MPI_DOUBLE, MPI_SUM, PCW);
//	printf("lrank = %d outflux = %e Influx = %e Outflux_bef=%e\n", params->rank, W_sum_data[0],W_total_u_influx, W_total_u_outflux);


}

/******************************************************************************/
