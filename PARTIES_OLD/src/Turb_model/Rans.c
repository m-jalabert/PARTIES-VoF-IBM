
#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Array.h"
#include "Conc.h"
#include "Communication.h"
#include "Display.h"
#include "Grid.h"
#include "Immersed.h"
#include "Memory.h"
#include "MyMath.h"
#include "Rans.h"
#include "Strain.h"
#include "Two_eqn_rans.h"
#include "Velocity.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/******************************************************************************/
/*
 This function allocates enough memory for the Rans_var structure based on
 parameters defined in "*params"
 */
/******************************************************************************/
Rans *Rans_create(MAC_grid *grid, Parameters *params) {
	Rans *rans;

	rans = (Rans *) malloc(sizeof(Rans));
	Memory_check_allocation(rans);

	rans->nut  = Memory_allocate_flow_variable(grid, params);
	rans->nutc = Memory_allocate_flow_variable(grid, params);

	rans->st_rate = (Strain_rate *) malloc(sizeof(Strain_rate));
	rans->st_rate->strain = Memory_allocate_flow_variable(grid, params);

#ifdef TWO_EQUATION_MODEL
    rans->two_eqn_rans = (Two_equation_rans **)malloc(2 * sizeof(Two_equation_rans *));
	rans->two_eqn_rans[0] = Two_equation_rans_create(grid, params);
	rans->two_eqn_rans[1] = Two_equation_rans_create(grid, params);
#endif

	rans->C_mu = 0.09;
	rans->C_muc = 0.09;

	rans->sigma_1 = 1.0;
	rans->sigma_2 = 1.3;

	rans->C_epsilon1 = 1.44;
	rans->C_epsilon2 = 1.92;
	rans->C_epsilon3 = 0.80;

	// Log law constants
	rans->kappa = 0.41;
	rans->b = 5.0;

	return rans;
}




/******************************************************************************/
/*
 This function releases the allocated memory for Rans structure.
 */
/******************************************************************************/
void Rans_destroy(Rans *rans, Parameters *params, MAC_grid *grid) {

#ifdef TWO_EQUATION_MODEL
	Two_equation_rans_destroy(rans->two_eqn_rans[0], params, grid);
	Two_equation_rans_destroy(rans->two_eqn_rans[1], params, grid);
	free(rans->two_eqn_rans);
#endif

	Memory_free_flow_variable(grid, params, rans->st_rate->strain);
	free(rans->st_rate);

	Memory_free_flow_variable(grid, params, rans->nut);
	Memory_free_flow_variable(grid, params, rans->nutc);

	free(rans);

}




/******************************************************************************/
/*
 This function integrates RANS transport equation for one time step
 */
/******************************************************************************/
void Rans_int_equations(Cart3d_bag *data_bag, double dt, Debug_trace *dtrace) {

	int i, j, k;
	int ierr;
	double T1, T2;
	FILE *fid;
	int pnodes;
	double ***nut;
	Wall_model *log_law;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Timer      *timer  = data_bag -> timer;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
	Pressure *p = data_bag -> p;
#ifdef CONC
	int iconc;
	int NConc = params->NConc;
	Concentration **c = data_bag -> c;
#endif

	Rans *rans = data_bag -> rans;

	pnodes = params->ghost_nodes;

	Display_progress(params," RANS start\n");

	// If which_stage == 0, velocity at cell center is calculated before the
	// to Dtime_cfl

	/*------------------------------------------------------------------------*/
	/*
	 Cell centered velocity: Communication is done for local ghost nodes within
	 the function
	 */
	/*------------------------------------------------------------------------*/
	if (params -> which_stage > 0) {
#ifdef INFLOW
		Inflow_velocity_profile(data_bag, DTRACE("Inflow_velocity_profile"));
#endif
		Velocity_cell_center(data_bag);
	}

#if defined CONC && defined INFLOW
	Inflow_conc_profile(data_bag, DTRACE("Inflow_conc_profile"));
#endif


	/*------------------------------------------------------------------------*/
	/*
	 Communicate cell-centered velocity
	 */
	/*------------------------------------------------------------------------*/
	Communication_update_ghost_nodes_flow_variable(u->data_bc, 'c', pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data_bc, 'c', pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data_bc, 'c', pnodes, data_bag);
#ifdef CONC
	for (iconc=0; iconc<NConc; iconc++) {
		Communication_update_ghost_nodes_flow_variable(c[iconc]->data, 'c', pnodes, data_bag);
	}
#endif


	/*------------------------------------------------------------------------*/
	/*
	 RANS stuff
	 */
	/*------------------------------------------------------------------------*/
	nut = rans->nut;
	log_law = rans->log_law;

#ifdef SCHUMANN
	Schumann_calc_shear(u->data, v->data, w->data, rans->log_law, grid, params);
#endif

#ifdef MIXING_LENGTH
	Rans_eddy_viscosity(rans, u, params, grid);
	Communication_update_ghost_nodes_flow_variable(rans->nut, 'c', pnodes, data_bag);
	Immersed_nut_interpolation(grid, params, rans->nut);
	Communication_update_ghost_nodes_flow_variable(rans->nut, 'c', pnodes, data_bag);
	return;
#endif

	T1 = MPI_Wtime();
	Rans_set_conv_viscous(u, v, w, rans, grid, params, 1);
	Rans_set_conv_viscous(u, v, w, rans, grid, params, 2);
	T2 = MPI_Wtime();
	timer->Wtime_rans_convective += T2 - T1;

	T1 = MPI_Wtime();
	Rans_set_RHS(rans->two_eqn_rans[0], grid, params, 1, dt);
	Rans_set_RHS(rans->two_eqn_rans[1], grid, params, 2, dt);
	Rans_add_source_RHS(data_bag);
	T2 = MPI_Wtime();
	timer->Wtime_rans_rhs += T2 - T1;

	T1 = MPI_Wtime();
#ifdef FULLY_EXPLICIT
	ierr = Two_equation_solve(data_bag, dt);
#endif
	Two_equation_set_boundary_values(rans, grid, params);
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[0]->data, 'c', pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[1]->data, 'c', pnodes, data_bag);
	T2 = MPI_Wtime();
	timer->Wtime_rans_solve += T2 - T1;

	Rans_eddy_viscosity(rans, u, params, grid);

#ifdef IMMERSED_BOUNDARY
	Communication_update_ghost_nodes_flow_variable(rans->nut, 'c', pnodes, data_bag);
	Immersed_nut_interpolation(grid, params, rans->nut);
	Immersed_keps_interpolation(rans, grid, params);
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[0]->data, 'c', pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[1]->data, 'c', pnodes, data_bag);
#endif
	Communication_update_ghost_nodes_flow_variable(rans->nut, 'c', pnodes, data_bag);

	return;
}




/******************************************************************************/
/*
 This function calculates the eddy viscosity and the eddy diffusivity
 */
/******************************************************************************/
void Rans_eddy_viscosity(Rans *rans, Velocity *u, Parameters *params, MAC_grid *grid) {

	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	int i, j, k;
	int NX, NY, NZ;
	double ***tke, ***diss;
	double ***nut, ***nutc;
	double C_mu, C_muc, ksqbyeps;
	double ***u_data;
	double **dudy_wm_b, **dwdy_wm_b, dUdy;
	double **dudy_wm_t, **dwdy_wm_t;
	double *yc;
	double ywall;
	double iRe;


	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	// indices start and end on current processor
	i_start = Is;
	j_start = Js;
	k_start = Ks;

	// exclude the half cell added
	i_end   = min(NX-1, Ie);
	j_end   = min(NY-1, Je);
	k_end   = min(NZ-1, Ke);

#ifdef TWO_EQUATION_MODEL
	tke  = rans->two_eqn_rans[0]->data;
	diss = rans->two_eqn_rans[1]->data;
#endif
	nut = rans->nut;
	u_data = u->data;
	yc = grid->yc;
	iRe = 1.0/params->Re;

	dudy_wm_b = rans->log_law->dudy_wm_bottom;
	dwdy_wm_b = rans->log_law->dwdy_wm_bottom;
	dudy_wm_t = rans->log_law->dudy_wm_top;
	dwdy_wm_t = rans->log_law->dwdy_wm_top;

	C_mu = rans->C_mu;
	C_muc = rans->C_muc;
#ifdef CONC
	nutc = rans->nutc;
#endif
	dUdy = 0.0;

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			ywall = yc[j]-grid->yv[0];
//			ywall = yc[j] -  grid->finer_1d_interface_y[i*grid->fine_step + grid->fine_step/2];
			ywall = min(ywall, grid->yv[NY-1]-grid->yc[j]);
			for (i=i_start; i<i_end; i++){

#ifdef TWO_EQUATION_MODEL

				ksqbyeps = tke[k][j][i]*tke[k][j][i]/diss[k][j][i];
				nut[k][j][i] = max(C_mu*ksqbyeps,iRe);
#ifdef CONC
				nutc[k][j][i] = C_muc*ksqbyeps;
#endif
#endif
			//nut[k][j][i] = 0.001;
		//	nutc[k][j][i] = 0.001;

#ifdef MIXING_LENGTH
// Mixing length model
				if ( (j > 0) && (j<NY-2))
					dUdy = (u_data[k][j+1][i]-u_data[k][j-1][i])/(yc[j+1]-yc[j-1]);


                if (j==0) {
                    dUdy = pow(dudy_wm_b[k][i]*dudy_wm_b[k][i] +
                            dwdy_wm_b[k][i]*dwdy_wm_b[k][i], 0.5);
					dUdy = (u_data[k][j+1][i]-u_data[k][j][i])/(yc[j+1]-yc[j]);
                }


#ifndef TOP_WALL_VELOCITY_FREESLIP
                if (j==NY-2) {
                    dUdy =  pow(dudy_wm_t[k][i]*dudy_wm_t[k][i] +
                            dwdy_wm_t[k][i]*dwdy_wm_t[k][i],0.5);
					dUdy = (u_data[k][j][i]-u_data[k][j-1][i])/(yc[j]-yc[j-1]);
                }
#endif


				nut[k][j][i] = pow(0.41*ywall,2)*fabs(dUdy);
	#ifdef CONC
				nutc[k][j][i] = nut[k][j][i];
	#endif
#endif



			}// for i
		} // for j
	}// for k

	return;

}

/******************************************************************************/
/*
 This function sets value of nut at the boundaries
 */
/******************************************************************************/
void Rans_set_boundary_values(double ***nut, MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ;
	int i, j, k;
	int ierr;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	// Same for all quantities
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

#ifdef TOP_WALL_VELOCITY_NOSLIP
	if (Je==NY) {
		j = NY-1;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				nut[k][j][i] = -nut[k][j-1][i];
			}
		}
	}
#elif defined TOP_WALL_VELOCITY_FREESLIP
	if (Je==NY) {
		j = NY-1;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				nut[k][j][i] = nut[k][j-1][i];
			}
		}
	}
#endif

#ifdef FRONT_WALL_VELOCITY_FREESLIP
	if (Ke==NZ) {
		k = NZ-1;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				nut[k][j][i] = nut[k-1][j][i];
			}
		}
	}
#elif defined FRONT_WALL_VELOCITY_NOSLIP
	if (Ke==NZ) {
		k = NZ-1;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				nut[k][j][i] = -nut[k-1][j][i];
			}
		}
	}
#endif

#if defined OUTFLOW || defined RIGHT_WALL_FREESLIP
	if (Ie==NX) {
		i = NX-1;
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				nut[k][j][i] = nut[k][j][i-1];
			}
		}
	}
#elif defined RIGHT_WALL_NOSLIP
	if (Ie==NX) {
		i = NX-1;
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				nut[k][j][i] = -nut[k][j][i-1];
			}
		}
	}
#endif

}
/******************************************************************************/
