#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Memory.h"
#include "Grid.h"
#include "Pressure.h"
#include "Communication.h"
#include "Cart3d.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

/******************************************************************************/
/*
 This function allocates memory for pressure structure
 */
/******************************************************************************/
Pressure *Pressure_create(MAC_grid *grid, Parameters *params) {

	Pressure *new_p;

	new_p = (Pressure *)malloc(sizeof(Pressure));
	Memory_check_allocation(new_p);

	// Allocate memory for pressure and pressure correction
	new_p->p_data = Memory_allocate_flow_variable(grid, params);
	new_p->deltap = Memory_allocate_flow_variable(grid, params);
	new_p->rhs    = Memory_allocate_flow_variable(grid, params);
	new_p->idxdt = Memory_allocate_1D_array(GVG_DOUBLE, grid->NX);
	new_p->idydt = Memory_allocate_1D_array(GVG_DOUBLE, grid->NY);
	new_p->idzdt = Memory_allocate_1D_array(GVG_DOUBLE, grid->NZ);
#ifdef POST_PROCESS
	new_p->p_data_avg = Memory_allocate_flow_variable(grid, params);
#endif

    // ALLOCATIONS for CG solver 
    new_p->res     = Memory_allocate_flow_variable(grid, params);  // residual array
    new_p->d       = Memory_allocate_flow_variable(grid, params);  // search direction
    new_p->Ad      = Memory_allocate_flow_variable(grid, params);  // operator applied to d
	new_p->M_inv   = Memory_allocate_flow_variable(grid, params);  // Preconditioner

	new_p->project_comm_cpu_time = 0.0;
	new_p->project_update_cpu_time = 0.0;
	new_p->project_u_cpu_time = 0.0;
	new_p->project_v_cpu_time = 0.0;
	new_p->project_w_cpu_time = 0.0;

	new_p->compute_divergence_comm_cpu_time = 0.0;
	new_p->compute_divergence_loop_cpu_time = 0.0;
	new_p->compute_divergence_reduce_cpu_time = 0.0;

	return new_p;
}




/******************************************************************************/
/*
 This function releases the allocated memory for pressure structure
 */
/******************************************************************************/
void Pressure_destroy(Pressure *p, MAC_grid *grid, Parameters *params) {

	Memory_free_flow_variable(grid, params, p->p_data);
	Memory_free_flow_variable(grid, params, p->deltap);
	Memory_free_flow_variable(grid, params, p->rhs);
#ifdef POST_PROCESS
	Memory_free_flow_variable(grid, params, p->p_data_avg);
#endif

    // FREES for CG arrays 
    Memory_free_flow_variable(grid, params, p->res);
    Memory_free_flow_variable(grid, params, p->d);
    Memory_free_flow_variable(grid, params, p->Ad);
	Memory_free_flow_variable(grid, params, p->M_inv);

	free(p->idxdt);
	free(p->idydt);
	free(p->idzdt);
	free(p);

}




/******************************************************************************/
/*
 This function computes the RHS vector {b} for the pressure linear set of
 equations
 */
/******************************************************************************/
void Pressure_set_RHS(Cart3d_bag *data_bag) {

	int i, j, k;
	double dudx, dvdy, dwdz;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Pressure *p = data_bag -> p;

	// Get regular data array for velocities data
	// First, generate the local (including ghost nodes) of velocity data on
	// current processor

	// Now, get the local velocity data on each processor including the ghost
	// nodes
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;

	// Now, got the RHS vector on current processor
	double ***rhs_vec = p -> rhs;

	double *idxdt = p -> idxdt;
	double *idydt = p -> idydt;
	double *idzdt = p -> idzdt;

	const double BET[] = {BETA};
	double a_dt = 1.0 / ( params -> dt * 2.0 * BET[params -> which_stage] );

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	double Lx = params->Lx;
	double Ly = params->Ly;
	double Lz = params->Lz;

	double *xc = grid->xc;
	double *yc = grid->yc;
	double *zc = grid->zc;
	double *yv = grid->yv;

	// Indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	// Exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie);
	int j_end = min(NY-1, grid->G_Je);
	int k_end = min(NZ-1, grid->G_Ke);


#ifdef GRID_UNIFORM

	double ihdt = grid->idx_u[1] * a_dt;
#else


	for (i = 0; i < NX; i++) {
		idxdt[i] = grid->idx_u[i] * a_dt;  // \partial_x/(2 \alpha \Delta t)
	}

	for (j = 0; j < NY; j++) {
		idydt[j] = grid->idy_v[j] * a_dt;
	}

	for (k = 0; k < NZ; k++) {
		idzdt[k] = grid->idz_w[k] * a_dt;
	}
#endif

	/*------------------------------------------------------------------------*/
	/*
	 Go over all the nodes and find divergence.
	 Exclude the last half cell added to the very end.
	 */
	/*------------------------------------------------------------------------*/
//#pragma omp parallel for private(j,i,dudx,dvdy,dwdz)
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				// Div_V_star at cell center p(i,j,k)
				// Regular node
				// Use central differences


#ifdef GRID_UNIFORM
				dudx = ( u_data[k][j][i+1] - u_data[k][j][i] ) * ihdt;
				dvdy = ( v_data[k][j+1][i] - v_data[k][j][i] ) * ihdt;
				dwdz = ( w_data[k+1][j][i] - w_data[k][j][i] ) * ihdt;

#else
				dudx = ( u_data[k][j][i+1] - u_data[k][j][i] ) * idxdt[i];
				dvdy = ( v_data[k][j+1][i] - v_data[k][j][i] ) * idydt[j];
				dwdz = ( w_data[k+1][j][i] - w_data[k][j][i] ) * idzdt[k];
#endif
				rhs_vec[k][j][i] = dudx + dvdy + dwdz; //cos(2*PI*xc[i]/Lx) *  cos(2*PI*yc[j]/Ly) *  cos(2*PI*zc[k]/Lz);//

			} /* for i*/
		} /* for j*/
	} /* for k*/

}




/******************************************************************************/
/*
 This function updates the velocity field using projection method to get a
 divergence free velocity field

 - Single-phase (no #define VOF_PLIC): uses rho_f = 1.
 - VOF multi-phase (#define VOF_PLIC in Boundary.h): uses per-cell rho from data_bag->vof->rho.

 The velocity correction is:
   u^k = u^* - 2 alpha_k dt (1/rho) ∇φ     (VOF case)
   u^k = u^* - 2 alpha_k dt         ∇φ   (single-phase case)

 The pressure update is:
   p^k = p^{k-1} + φ            (VOF)
   p^k = p^{k-1} + φ         (single-phase)
*/
/******************************************************************************/
void Pressure_project_velocity(Cart3d_bag *data_bag) {

	int i, j, k;
	double T1, T2;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Pressure *p = data_bag -> p;

	// Get velocity data. No ghost node is required. So, Global data is
	// retrieved
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***p_data = p -> p_data;

#ifdef VOF_PLIC
    // For VOF (variable density), fetch the local density field
    VolumeFraction *vof = data_bag->vof;
    double ***rho       = vof->rho;
#endif

	// Start index of bottom-left-back corner on current processor
	int Is = grid->G_Is;
	int Js = grid->G_Js;
	int Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid->G_Ie;
	int Je = grid->G_Je;
	int Ke = grid->G_Ke;

	// Same for all quantities
	int NX = grid->NX;
	int NY = grid->NY;
	int NZ = grid->NZ;

	// indices start and end on current processor
	int i_start = Is;
	int j_start = Js;
	int k_start = Ks;

	// exclude the half cell added
	int i_end = min(NX-1, Ie);
	int j_end = min(NY-1, Je);
	int k_end = min(NZ-1, Ke);

	double *idxdt = p->idxdt;
	double *idydt = p->idydt;
	double *idzdt = p->idzdt;

#ifdef POST_PROCESS
	double ***p_data_avg = p -> p_data_avg;
	if (params->which_stage == 0) {
		Memory_reset_flow_variable(grid, params, p_data_avg);
	}
#endif

	// Now, get "delta_p"(pressure correction) (Poisson equation variable)
	//local data since ghost nodes are required
	double ***deltap = p -> deltap;


	// Phi. Get the ghost nodes from the neighboring processors
	T1 = MPI_Wtime();
	Communication_update_ghost_nodes_flow_variable(deltap, CONCENTRATION_PERTURBATION, 3, data_bag);
	T2 = MPI_Wtime();



	p->project_comm_cpu_time += T2-T1;

	const double BET[] = {BETA};
	double mdt = 2.0 * BET[params -> which_stage] * params -> dt;



	for (i=0;i<NX;i++) {
		idxdt[i] = grid->idx_c[i] * mdt;
	}

	for (j=0;j<NY;j++) {
		idydt[j] = grid->idy_c[j] * mdt;
	}

	for (k=0;k<NZ;k++) {
		idzdt[k] = grid->idz_c[k] * mdt;
	}

	//--------------------------------------------------------------------------
	// Update Pressure 
	// p^k = p^{k-1} + phi   (VOF)
	// p^k = p^{k-1} + phi   (single-phase)
	//--------------------------------------------------------------------------
	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

#ifdef VOF_PLIC
                // Multi-phase
                //p_data[k][j][i] += rho[k][j][i] * deltap[k][j][i];
				p_data[k][j][i] += deltap[k][j][i];
#else
                // Single-phase
                p_data[k][j][i] += deltap[k][j][i];
#endif

#ifdef POST_PROCESS
				p_data_avg[k][j][i] += 2.0 * BET[params->which_stage] * p_data[k][j][i];
#endif
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	p->project_update_cpu_time += T2-T1;

	T1 = MPI_Wtime();
	Pressure_apply_BCs(p_data, grid, params);
	Communication_update_ghost_nodes_flow_variable(p_data, CONCENTRATION_PERTURBATION, 3, data_bag);
	T2 = MPI_Wtime();
	p->project_comm_cpu_time += T2-T1;


	//--------------------------------------------------------------------------
	// Update u_star to u_new (divergence free velocity field)
	//    u^k = u^* - ( deltap[i] - deltap[i-1] ) * idxdt[i-1] / rho^k-1  (VOF)
    //    u^k = u^* - ( deltap[i] - deltap[i-1] ) * idxdt[i-1]        (single-phase)
    //
    //    Repeat similarly for v, w in j, k directions.
	//--------------------------------------------------------------------------
	i_start = max(1, Is); // i=0 not included
	j_start = Js;
	k_start = Ks;

	i_end = min(NX-1, Ie);
	j_end = min(NY-1, Je);
	k_end = min(NZ-1, Ke);
#ifdef XPERIODIC
	if (i_end == NX-1)
		i_end = NX;
#endif

	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

#ifdef VOF_PLIC
				double inv_rho_face = 2.0 / (rho[k][j][i] + rho[k][j][i-1]);
                u_data[k][j][i] -= ( deltap[k][j][i] - deltap[k][j][i-1] )
                                   * ( idxdt[i-1] * inv_rho_face );
#else
                u_data[k][j][i] -= ( deltap[k][j][i] - deltap[k][j][i-1] )
                                   * idxdt[i-1];
#endif
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	p->project_u_cpu_time += T2-T1;


	//--------------------------------------------------------------------------
	// Update v_star to v_new (divergence free velocity field)
	//--------------------------------------------------------------------------
	i_start = Is;
	j_start = max(1, Js); // j=0 not included
	k_start = Ks;

	i_end = min(NX-1, Ie);
	j_end = min(NY-1, Je);
	k_end = min(NZ-1, Ke);


#ifdef YPERIODIC
if (j_end == NY-1)
{
	j_end = NY;
}
#endif



	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
				
#ifdef VOF_PLIC
				double inv_rho_face = 2.0 / (rho[k][j][i] + rho[k][j-1][i]);
                v_data[k][j][i] -= ( deltap[k][j][i] - deltap[k][j-1][i] )
                                   * ( idydt[j-1] * inv_rho_face );
#else
                v_data[k][j][i] -= ( deltap[k][j][i] - deltap[k][j-1][i] )
                                   * idydt[j-1];
#endif
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	p->project_v_cpu_time += T2-T1;


	//--------------------------------------------------------------------------
	// Update w_star to w_new (divergence free velocity field)
	//--------------------------------------------------------------------------
	i_start = Is;
	j_start = Js;
	k_start = max(1, Ks); // k=0 not included

	i_end = min(NX-1, Ie);
	j_end = min(NY-1, Je);
	k_end = min(NZ-1, Ke);
#ifdef ZPERIODIC
	if (k_end == NZ-1)
		k_end = NZ;
#endif

	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

#ifdef VOF_PLIC
				double inv_rho_face = 2.0 / (rho[k][j][i] + rho[k-1][j][i]);
                w_data[k][j][i] -= ( deltap[k][j][i] - deltap[k-1][j][i] )
                                   * ( idzdt[k-1] * inv_rho_face );
#else
                w_data[k][j][i] -= ( deltap[k][j][i] - deltap[k-1][j][i] )
                                   * idzdt[k-1];
#endif
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	p->project_w_cpu_time += T2-T1;

}




/******************************************************************************/
/*
 This function finds the maximum divergence of the velocity field to verify the
 solution, i.e. div_max --> 0.0
 */
/******************************************************************************/
double Pressure_compute_velocity_divergence(Cart3d_bag *data_bag) {

	int i, j, k;
	double dudx, dvdy, dwdz;
	double div_cell= 0;
	double div_max = 0.0; // local
	double W_div_max = 0.0; // (World) global
	int i_max, j_max, k_max;
	int i_maxl, j_maxl, k_maxl, iflag;
	int proc_maxl, proc_max;
	double l1_div = 0.0; // local
	double W_l1_div = 0.0; // (World) global
	double l2_div = 0.0; // local
	double W_l2_div = 0.0; // (World) global
	double W_data[2], data[2];
	int W_idata[4], idata[4];
	double T1, T2;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Get velocity data. Including ghost nodes
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;

	/*------------------------------------------------------------------------*/
	/*
	 Generate the local u, v, w data (including the ghost nodes from neighboring
	 processors)
	 */
	/*------------------------------------------------------------------------*/
	T1 = MPI_Wtime();
	// u-vel
//	Communication_update_ghost_nodes_flow_variable(data, 'u', 1, data_bag);
	// v-vel
//	Communication_update_ghost_nodes_flow_variable(data, 'v', 1, data_bag);
	// w-vel
//	Communication_update_ghost_nodes_flow_variable(data, 'w', 1, data_bag);
	T2 = MPI_Wtime();
	data_bag->p->compute_divergence_comm_cpu_time += T2-T1;

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	// exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie);
	int j_end = min(NY-1, grid->G_Je);
	int k_end = min(NZ-1, grid->G_Ke);

	T1 = MPI_Wtime();


	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				// Use central differences
				dudx = (u_data[k][j][i+1] - u_data[k][j][i]) * idx_u[i];
				dvdy = (v_data[k][j+1][i] - v_data[k][j][i]) * idy_v[j];
				dwdz = (w_data[k+1][j][i] - w_data[k][j][i]) * idz_w[k];

				div_cell = fabs(dudx + dvdy + dwdz);
				l1_div = l1_div + div_cell;
				l2_div = l2_div + div_cell*div_cell;

				// Now check the maximum divergence in the whole domain
				if (div_cell >= div_max) {
					div_max = div_cell;
					i_max = i;
					j_max = j;
					k_max = k;
				}

			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	data_bag->p->compute_divergence_loop_cpu_time += T2-T1;


	// Now, find the maximum global divergence on all processors
	T1 = MPI_Wtime();
	MPI_Allreduce (&div_max, &W_div_max, 1, MPI_DOUBLE, MPI_MAX, PCW);
	iflag = -10;

	if (div_max == W_div_max) iflag = 1;
    idata[0] = i_max * iflag;
	idata[1] = j_max * iflag;
	idata[2] = k_max * iflag;
	idata[3] = params->rank * iflag;
	data[0] = l1_div;
	data[1] = l2_div;
	MPI_Reduce(idata, W_idata, 4, MPI_INT, MPI_MAX, 0, PCW);
	MPI_Reduce(data, W_data, 2, MPI_DOUBLE, MPI_SUM, 0, PCW);
	i_max = W_idata[0];
	j_max = W_idata[1];
	k_max = W_idata[2];
	proc_max = W_idata[3];
	W_l1_div = W_data[0];
	W_l2_div = W_data[1];
	T2 = MPI_Wtime();
	data_bag->p->compute_divergence_reduce_cpu_time += T2-T1;

	W_l1_div = W_l1_div / ( (NX-1) * (NY-1) * (NZ-1) );
	W_l2_div = W_l2_div / ( (NX-1) * (NY-1) * (NZ-1) );

	// if (params->rank == 0) {
	// 	printf("Pressure.c/ div:%20.16e at (i,j,k)=(%d,%d,%d) at proc = %d\n",
	// 	       W_div_max, i_max, j_max, k_max,proc_max);
	// 	printf("Pressure.c/ l1_div=%16.12e l2_div= %16.12e\n",
	// 	       W_l1_div, W_l2_div);
	// }
	return (W_div_max);
}

#ifdef VOF_SCALAR

void Pressure_set_RHS_vof(Cart3d_bag *data_bag) {

	int i, j, k;
	double dudx, dvdy, dwdz;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Pressure *p = data_bag -> p;

	// Get regular data array for velocities data
	// First, generate the local (including ghost nodes) of velocity data on
	// current processor

	// Now, get the local velocity data on each processor including the ghost
	// nodes
	double ***u_data = data_bag -> u ->data_vof;
	double ***v_data = data_bag -> v -> data_vof;
	double ***w_data = data_bag -> w -> data_vof;

	// Now, got the RHS vector on current processor
	double ***rhs_vec = p -> rhs;

	double *idxdt = p -> idxdt;
	double *idydt = p -> idydt;
	double *idzdt = p -> idzdt;

	const double BET[] = {BETA};
	double a_dt = 1.0 / ( params -> dt * 2.0 * BET[params -> which_stage] );

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	// Exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie);
	int j_end = min(NY-1, grid->G_Je);
	int k_end = min(NZ-1, grid->G_Ke);

#ifdef GRID_UNIFORM

	double ihdt = grid->idx_u[1] * a_dt;
#else


	for (i = 0; i < NX; i++) {
		idxdt[i] = grid->idx_u[i] * a_dt;  // \partial_x/(2 \alpha \Delta t)
	}

	for (j = 0; j < NY; j++) {
		idydt[j] = grid->idy_v[j] * a_dt;
	}

	for (k = 0; k < NZ; k++) {
		idzdt[k] = grid->idz_w[k] * a_dt;
	}
#endif




	/*------------------------------------------------------------------------*/
	/*
	 Go over all the nodes and find divergence.
	 Exclude the last half cell added to the very end.
	 */
	/*------------------------------------------------------------------------*/
//#pragma omp parallel for private(j,i,dudx,dvdy,dwdz)
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				// Div_V_star at cell center p(i,j,k)
				// Regular node
				// Use central differences

#ifdef GRID_UNIFORM
				dudx = ( u_data[k][j][i+1] - u_data[k][j][i] ) * ihdt;
				dvdy = ( v_data[k][j+1][i] - v_data[k][j][i] ) * ihdt;
				dwdz = ( w_data[k+1][j][i] - w_data[k][j][i] ) * ihdt;

#else
				dudx = ( u_data[k][j][i+1] - u_data[k][j][i] ) * idxdt[i];
				dvdy = ( v_data[k][j+1][i] - v_data[k][j][i] ) * idydt[j];
				dwdz = ( w_data[k+1][j][i] - w_data[k][j][i] ) * idzdt[k];
#endif
				rhs_vec[k][j][i] = dudx + dvdy + dwdz;

			} /* for i*/
		} /* for j*/
	} /* for k*/

}


void Pressure_project_velocity_vof(Cart3d_bag *data_bag) {

	int i, j, k;
	double T1, T2;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Pressure *p = data_bag -> p;

	// Get velocity data. No ghost node is required. So, Global data is
	// retrieved
	double ***u_data = data_bag -> u -> data_vof;
	double ***v_data = data_bag -> v -> data_vof;
	double ***w_data = data_bag -> w -> data_vof;


	// Now, get "delta_p"(pressure correction) (Poisson equation variable)
	//local data since ghost nodes are required
	double ***deltap = p -> deltap;

	// Phi. Get the ghost nodes from the neighboring processors
	T1 = MPI_Wtime();
	Communication_update_ghost_nodes_flow_variable(deltap, CONCENTRATION_PERTURBATION, 1, data_bag);
	T2 = MPI_Wtime();
	p->project_comm_cpu_time += T2-T1;

	const double BET[] = {BETA};
	double mdt = 2.0 * BET[params -> which_stage] * params -> dt;

	// Start index of bottom-left-back corner on current processor
	int Is = grid->G_Is;
	int Js = grid->G_Js;
	int Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid->G_Ie;
	int Je = grid->G_Je;
	int Ke = grid->G_Ke;

	// Same for all quantities
	int NX = grid->NX;
	int NY = grid->NY;
	int NZ = grid->NZ;

	// indices start and end on current processor
	int i_start = Is;
	int j_start = Js;
	int k_start = Ks;

	// exclude the half cell added
	int i_end = min(NX-1, Ie);
	int j_end = min(NY-1, Je);
	int k_end = min(NZ-1, Ke);


#ifdef GRID_UNIFORM

	double ihdt= grid->idx_c[1] * mdt;

#else


	double *idxdt = p->idxdt;
	double *idydt = p->idydt;
	double *idzdt = p->idzdt;

	for (i=0;i<NX;i++) {
		idxdt[i] = grid->idx_c[i] * mdt;
	}

	for (j=0;j<NY;j++) {
		idydt[j] = grid->idy_c[j] * mdt;
	}

	for (k=0;k<NZ;k++) {
		idzdt[k] = grid->idz_c[k] * mdt;
	}
#endif

	//--------------------------------------------------------------------------
	// Update u_star to u_new (divergence free velocity field)
	//--------------------------------------------------------------------------
	i_start = max(1, Is); // i=0 not included
	j_start = Js;
	k_start = Ks;

	i_end = min(NX-1, Ie);
	j_end = min(NY-1, Je);
	k_end = min(NZ-1, Ke);
#ifdef XPERIODIC
	if (i_end == NX-1)
		i_end = NX;
#endif

	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

#ifdef GRID_UNIFORM
				u_data[k][j][i] -= (deltap[k][j][i] - deltap[k][j][i-1]) * ihdt;
#else
				u_data[k][j][i] -= (deltap[k][j][i] - deltap[k][j][i-1]) * idxdt[i-1];
#endif
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	p->project_u_cpu_time += T2-T1;


	//--------------------------------------------------------------------------
	// Update v_star to v_new (divergence free velocity field)
	//--------------------------------------------------------------------------
	i_start = Is;
	j_start = max(1, Js); // j=0 not included
	k_start = Ks;

	i_end = min(NX-1, Ie);
	j_end = min(NY-1, Je);
	k_end = min(NZ-1, Ke);
#ifdef YPERIODIC
	if (j_end == NY-1)
		j_end = NY;
#endif

	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

#ifdef GRID_UNIFORM
				v_data[k][j][i] -= (deltap[k][j][i] - deltap[k][j-1][i]) * ihdt;
#else
				v_data[k][j][i] -= (deltap[k][j][i] - deltap[k][j-1][i]) * idydt[j-1];
#endif
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	p->project_v_cpu_time += T2-T1;


	//--------------------------------------------------------------------------
	// Update w_star to w_new (divergence free velocity field)
	//--------------------------------------------------------------------------
	i_start = Is;
	j_start = Js;
	k_start = max(1, Ks); // k=0 not included

	i_end = min(NX-1, Ie);
	j_end = min(NY-1, Je);
	k_end = min(NZ-1, Ke);
#ifdef ZPERIODIC
	if (k_end == NZ-1)
		k_end = NZ;
#endif

	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
#ifdef GRID_UNIFORM
				w_data[k][j][i] -= (deltap[k][j][i] - deltap[k-1][j][i]) * ihdt;
#else
				w_data[k][j][i] -= (deltap[k][j][i] - deltap[k-1][j][i]) * idzdt[k-1];
#endif
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	p->project_w_cpu_time += T2-T1;

}
#endif

#include "lsolver/psolve_fft.c"
#include "lsolver/psolve_cg.c"
