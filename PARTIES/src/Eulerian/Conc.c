#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Conc.h"
#include "Velocity.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Immersed.h"
#include "Display.h"
#include "Array.h"
#include "Cart3d.h"
#include "Lagrangian.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>





static double Conc_innerProd(double ***vec1, double ***vec2,  MAC_grid *grid, Parameters *params) {

	int i, j, k;

	double partialSum, totalSum;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = min( grid->G_Ie, NX-1);  // No ghost nodes
	int Je = min( grid->G_Je, NY-1);
	int Ke = min( grid->G_Ke, NZ-1);


	// Calculate inner product of local nodes
	partialSum = 0;
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {

				partialSum += vec1[k][j][i] * vec2[k][j][i];

			}
		}
	}

	// Communicate and sum inner product among all processors
	MPI_Allreduce(&partialSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, PCW);

	return totalSum;
}


/******************************************************************************/
/*
 Completes matrix-vector multiplication between 'A' matrix and input vector 'x'
 */
/******************************************************************************/
static void Conc_laplacian(double ***Ax, double ***x,  MAC_grid *grid, Parameters *params, int iconc) {

	int i, j, k;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;


	int Ie = min( grid->G_Ie, NX-1);  // No ghost nodes
	int Je = min( grid->G_Je, NY-1);
	int Ke = min( grid->G_Ke, NZ-1);




	// Grid spacing (assuming uniform)
	double iddx = grid -> idx_c[1];
	double iddy = grid -> idy_c[1];
	double iddz = grid -> idz_c[1];
	iddx = iddx * iddx;
	iddy = iddy * iddy;
	iddz = iddz * iddz;

	double iPe = 1.0 / params -> Pe[iconc];
	const double BET[] = {BETA};
	double idtimeb = 1.0 / (BET[params -> which_stage] * params -> dt);

	double ac = idtimeb + 2.0 * iPe * (iddx + iddy + iddz);
	double ax = -iPe * iddx;
	double ay = -iPe * iddy;
	double az = -iPe * iddz;

	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				Ax[k][j][i] = ac * x[k][j][i]
				            + ax * ( x[k][j][i-1] + x[k][j][i+1] )
				            + ay * ( x[k][j-1][i] + x[k][j+1][i] )
				            + az * ( x[k-1][j][i] + x[k+1][j][i] );
			}
		}
	}
}
/******************************************************************************/
/*
 Completes matrix-vector multiplication between 'A' matrix and input vector 'x' with VOF
 */
/******************************************************************************/

static void Conc_laplacian_vof(double ***Ax, double ***x,  double ***vfu,double ***vfv, double ***vfw,double ***ng_vfc , MAC_grid *grid, Parameters *params, int iconc) {

	int i, j, k;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;


	int Ie = min( grid->G_Ie, NX-1);  // No ghost nodes
	int Je = min( grid->G_Je, NY-1);
	int Ke = min( grid->G_Ke, NZ-1);



	// Grid spacing (assuming uniform)
	double iddx = grid -> idx_c[1];
	double iddy = grid -> idy_c[1];
	double iddz = grid -> idz_c[1];
	iddx = iddx * iddx;
	iddy = iddy * iddy;
	iddz = iddz * iddz;

	double iPe = 1.0 / params -> Pe[iconc];
	const double BET[] = {BETA};
	double idtimeb = 1.0 / (BET[params -> which_stage] * params -> dt);



	double conductivity_s = params->conductivity_s[iconc];
	double heat_cap_s = params->vol_heat_cap_s[iconc];
	double ih= grid -> idx_u[1];
	double ihsq=ih*ih;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double lambdaN, lambdaS, lambdaE, lambdaW, lambdaB, lambdaF ;


	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {

					/*------------------------------------------------------------*/
					/*
					 Calculate diffusive terms
					 */
					/*------------------------------------------------------------*/

					//--------------------------------------------------------------
					// Variable conductivity with the arithmetic average
					//--------------------------------------------------------------
					lambdaE =  vfu[k][j][i+1]*(conductivity_s-1)+1;
					lambdaW =  vfu[k][j][i]*(conductivity_s-1)+1;

					lambdaN = vfv[k][j+1][i]*(conductivity_s-1)+1;
					lambdaS = vfv[k][j][i]*(conductivity_s-1)+1;

					lambdaF = vfw[k+1][j][i]*(conductivity_s-1)+1;
					lambdaB = vfw[k][j][i]*(conductivity_s-1)+1;




					dcdxE = ( x[k][j][i+1] - x[k][j][i] );
					dcdxW = ( x[k][j][i] - x[k][j][i-1] ) ;
					dcdyN = ( x[k][j+1][i] - x[k][j][i] ) ;
					dcdyS = ( x[k][j][i] - x[k][j-1][i] ) ;
					dcdzF = ( x[k+1][j][i]-x[k][j][i] ) ;
					dcdzB = ( x[k][j][i]-x[k-1][j][i] ) ;


					d2cdx2 = ( lambdaE * dcdxE - lambdaW * dcdxW ) ;
					d2cdy2 = ( lambdaN * dcdyN - lambdaS * dcdyS ) ;
					d2cdz2 = ( lambdaF * dcdzF - lambdaB * dcdzB ) ;




				Ax[k][j][i] = idtimeb * x[k][j][i]- ihsq*iPe*(d2cdx2 + d2cdy2 + d2cdz2)/( ng_vfc[k][j][i]*(heat_cap_s-1)+1 );

			}
		}
	}
}


/******************************************************************************/
/*
 This function allocates enough memory for the concentration structure based on
 parameters defined in "*params"
 */
/******************************************************************************/
Concentration *Conc_create(int conc_index, MAC_grid *grid, Parameters *params) {

	Concentration *new_conc;
	int NX, NY, NZ;
	int NX_cell, NY_cell, NZ_cell;
	int Ghost_Nodes;
	int ny, nz;

	// Number of real physical cells
	NX_cell = params->NXM;
	NY_cell = params->NYM;
	NZ_cell = params->NZM;

	// Add one cell in each direction to the far most end
	NX = NX_cell + 1; // Number of grid points in x-direction (half cell added in the far most direction)
	NY = NY_cell + 1; // Number of grid points in y-direction (half cell added in the far most direction)
	NZ = NZ_cell + 1; // Number of grid points in z-direction (half cell added in the far most direction)

	// We extend the Field region by one column to the left and right and also
	// one row below and above the regular region. Reason: ASK BRENDON!!
	new_conc = (Concentration *)malloc(sizeof(Concentration));
	Memory_check_allocation(new_conc);

	new_conc->data = Memory_allocate_flow_variable(grid, params);
	new_conc->ng_data_old = Memory_allocate_noghost_variable(grid, params);
#ifdef CONC_BQUICK
	new_conc->blend = Memory_allocate_flow_variable(grid, params);
	new_conc->data_temp = Memory_allocate_flow_variable(grid, params);
	new_conc->ng_temp1 = Memory_allocate_noghost_variable(grid, params);
#endif


#ifdef CONC_QUICK
	Conc_set_quick_coefficients(new_conc, grid, params);
#endif

	new_conc->ng_conv = Memory_allocate_noghost_variable(grid, params);
	new_conc->ng_conv_old = Memory_allocate_noghost_variable(grid, params);
	new_conc->ng_viscous = Memory_allocate_noghost_variable(grid, params);
	new_conc->ng_visc_explicit = Memory_allocate_noghost_variable(grid, params);
	new_conc->ng_rhs = Memory_allocate_noghost_variable(grid, params);
#ifdef IBM_SCALAR
//	new_conc->heat_ibm = Memory_allocate_noghost_variable(grid, params);
#endif


#ifdef CONC_SOLVE_CG
	if(conc_index == 0){
		new_conc -> d = Memory_allocate_flow_variable(grid, params);
		new_conc -> ng_r = Memory_allocate_noghost_variable(grid, params);
		new_conc -> ng_Ad = Memory_allocate_noghost_variable(grid, params);
	}
#endif

#ifdef CONC_BQUICK
	new_conc->ng_conv_quick = Memory_allocate_noghost_variable(grid, params);
#endif


	// Total concentration: To be used in v-momentum equation. Always store the
	// total concentration in c[0]->G_c_total
	if ( (params->NConc > 1) && (conc_index == 0) ){
		new_conc->c_total = Memory_allocate_flow_variable(grid, params);
	} // else, only one concentration field c[0]


	/* concentration averaged height (only a function of x) */
	if (params->ave_height_output) {

		new_conc->G_ave_height_x = Memory_allocate_1D_array(GVG_DOUBLE, NX);
		new_conc->W_ave_height_x = Memory_allocate_1D_array(GVG_DOUBLE, NX);

	}// if


	// Corresponing Peclet Number
	    new_conc->Pe   = params->Pe[conc_index];

	// PARTICLE: Settling Speed
		new_conc->v_settl0       = params->V_s0[conc_index];


	// particle deposit height
	if ( (params->conc_output_deposit_height) ) {

		new_conc->G_deposit_height = Memory_allocate_2D_double_array(NX, NZ);
		new_conc->W_deposit_height = Memory_allocate_2D_double_array(NX, NZ);

		new_conc->G_deposit_height_dumped = Memory_allocate_2D_double_array(NX, NZ);
		new_conc->W_deposit_height_dumped = Memory_allocate_2D_double_array(NX, NZ);
	}
	else {
		new_conc->G_deposit_height = NULL;
		new_conc->W_deposit_height = NULL;

		new_conc->G_deposit_height_dumped = NULL;
		new_conc->W_deposit_height_dumped = NULL;
	} // else

#ifdef IMMERSED_BOUNDARY
	new_conc->G_conc_immersed = Memory_allocate_2D_double_array(NX, NZ);
	new_conc->W_conc_immersed = Memory_allocate_2D_double_array(NX, NZ);
#endif

if(params->conc_init_type[conc_index] == 1  ){
	new_conc-> inflow_n = Memory_allocate_2D_double_array(NY, NZ);
}



#ifdef CONC_SINGLE_BC_NODIFF
params->BC_AN[conc_index]= 1;
params->BC_AS[conc_index]= 1;
params->BC_AE[conc_index]= 1;
params->BC_AW[conc_index]= 1;
params->BC_AF[conc_index]= 1;
params->BC_AB[conc_index]= 1;

params->BC_BN[conc_index]= 0;
params->BC_BS[conc_index]= 0;
params->BC_BE[conc_index]= 0;
params->BC_BW[conc_index]= 0;
params->BC_BF[conc_index]= 0;
params->BC_BB[conc_index]= 0;

params->BC_CN[conc_index]= 0;
params->BC_CS[conc_index]= 0;
params->BC_CE[conc_index]= 0;
params->BC_CW[conc_index]= 0;
params->BC_CF[conc_index]= 0;
params->BC_CB[conc_index]= 0;

#elif defined CONC_SINGLE_BC_RB

params->BC_AN[conc_index]= 0;
params->BC_AS[conc_index]= 0;
params->BC_AE[conc_index]= 1;
params->BC_AW[conc_index]= 1;
params->BC_AF[conc_index]= 1;
params->BC_AB[conc_index]= 1;

params->BC_BN[conc_index]= 1;
params->BC_BS[conc_index]= 1;
params->BC_BE[conc_index]= 0;
params->BC_BW[conc_index]= 0;
params->BC_BF[conc_index]= 0;
params->BC_BB[conc_index]= 0;

params->BC_CN[conc_index]= 0;
params->BC_CS[conc_index]= 1;
params->BC_CE[conc_index]= 0;
params->BC_CW[conc_index]= 0;
params->BC_CF[conc_index]= 0;
params->BC_CB[conc_index]= 0;

#elif defined CONC_SINGLE_BC_ZERO

	params->BC_AN[conc_index]= 0;
	params->BC_AS[conc_index]= 0;
	params->BC_AE[conc_index]= 0;
	params->BC_AW[conc_index]= 0;
	params->BC_AF[conc_index]= 0;
	params->BC_AB[conc_index]= 0;

	params->BC_BN[conc_index]= 1;
	params->BC_BS[conc_index]= 1;
	params->BC_BE[conc_index]= 1;
	params->BC_BW[conc_index]= 1;
	params->BC_BF[conc_index]= 1;
	params->BC_BB[conc_index]= 1;

	params->BC_CN[conc_index]= 0;
	params->BC_CS[conc_index]= 0;
	params->BC_CE[conc_index]= 0;
	params->BC_CW[conc_index]= 0;
	params->BC_CF[conc_index]= 0;
	params->BC_CB[conc_index]= 0;
#endif









	return new_conc;
}




/******************************************************************************/
/*
 This function releases the allocated memory for velocity structure.
 */
/******************************************************************************/
void Conc_destroy(Concentration *c, int conc_index, MAC_grid *grid,
		Parameters *params) {

	int NZ = grid->NZ;

	Memory_free_flow_variable(grid, params, c->data);
	Memory_free_noghost_variable(grid, params, c->ng_data_old);

	Memory_free_noghost_variable(grid, params, c->ng_conv);
	Memory_free_noghost_variable(grid, params, c->ng_conv_old);
	Memory_free_noghost_variable(grid, params, c->ng_viscous);
	Memory_free_noghost_variable(grid, params, c->ng_visc_explicit);


#ifdef CONC_SOLVE_CG
	if(conc_index == 0){

		Memory_free_flow_variable(grid, params, c->d);
		Memory_free_noghost_variable(grid, params, c->ng_r);
		Memory_free_noghost_variable(grid, params, c->ng_Ad);
		}
#endif



#ifdef CONC_BQUICK
	Memory_free_flow_variable(grid, params, c->blend);
	Memory_free_flow_variable(grid, params, c->data_temp);
	Memory_free_noghost_variable(grid, params, c->ng_temp1);
	Memory_free_noghost_variable(grid, params, c->ng_conv_quick);
#endif

	Memory_free_noghost_variable(grid, params, c->ng_rhs);
#ifdef IBM_SCALAR
//	Memory_free_noghost_variable(grid, params, c->heat_ibm);
#endif

	if ( (params->NConc > 1) && (conc_index == 0) ) {
		Memory_free_flow_variable(grid, params, c->c_total);
	} // if

	if (params->ave_height_output) {
		free(c->G_ave_height_x);
		free(c->W_ave_height_x);
	}// if

	if (c->G_deposit_height != NULL) {
		Memory_free_2D_double_array(NZ, c->G_deposit_height);
		Memory_free_2D_double_array(NZ, c->W_deposit_height);
	} // if

	if (c->G_deposit_height_dumped != NULL) {
		Memory_free_2D_double_array(NZ, c->G_deposit_height_dumped);
		Memory_free_2D_double_array(NZ, c->W_deposit_height_dumped);
	} // if

#ifdef IMMERSED_BOUNDARY
	Memory_free_2D_double_array(NZ, c->G_conc_immersed);
	Memory_free_2D_double_array(NZ, c->W_conc_immersed);
#endif

	free(c);

}




/******************************************************************************/
/*
 This function integrates Conc transport equation for one time step
 */
/******************************************************************************/
void Conc_int_equations(Cart3d_bag *data_bag, Debug_trace *dtrace ) {

	int i, j, k;
	double T1, T2, Tstart, Tend;
	int implicit_loops_count;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Timer      *timer  = data_bag -> timer;


	int iconc;
	int NConc = params -> NConc;
	Concentration **c = data_bag -> c;

	int pnodes = params -> ghost_nodes;

	Tstart = MPI_Wtime();
	//--------------------------------------------------------------------------
	// Set convective and viscous terms
	//--------------------------------------------------------------------------
	T1 = MPI_Wtime();
for (iconc=0; iconc<NConc; iconc++) {

	Conc_set_conv_viscous(iconc, data_bag);

	T2 = MPI_Wtime();
	timer->Wtime_c_convective += T2 - T1;

		//----------------------------------------------------------------------
		// Set RHS
		//----------------------------------------------------------------------
	T1 = MPI_Wtime();
	Conc_set_RHS(c[iconc], grid, params);



#ifdef THERMAL_KADER
		Conc_add_source_RHS(iconc, data_bag);
#endif
	T2 = MPI_Wtime();
	timer->Wtime_c_rhs += T2 - T1;

	//----------------------------------------------------------------------
	// Solve for concentration
	//----------------------------------------------------------------------
		T1 = MPI_Wtime();
#ifdef IBM_SCALAR
		Conc_solve_fully_explicit(iconc, data_bag);
	#elif  CONC_SEMI_IMPLICIT
		Conc_solve_semi_implicit(iconc, data_bag);
	#elif  defined CONC_FULLY_EXPLICIT
		Conc_solve_fully_explicit(iconc, data_bag);
	#elif defined CONC_SOLVE_CG
		Conc_solve_cg(iconc, data_bag);
#endif


#ifdef IBM_SCALAR

	data_bag->c[iconc]->heating_flag=0; // == 0 a new heating period begins

    //// calculate LA forcing and add to RHS
	Tstart = MPI_Wtime();

	Conc_set_RHS_IBM_implicit(c[iconc], data_bag->lag->ng_temp ,grid, params); // overwrite ng_rhs!!!

#ifndef FAST
	Communication_update_ghost_nodes_flow_variable(c[iconc]->data, 'c',params->ghost_nodes , data_bag); // sollte kaum notwendig sein
#endif

	Lagrangian_heat(iconc, -1, data_bag, DTRACE("Lagrangian_heat") );


	Communication_update_ghost_nodes_flow_variable(c[iconc]->data, 'c',params->ghost_nodes , data_bag);
	Conc_set_boundary_values(c[iconc]->data, iconc ,CENTRAL_FULL ,grid , params);


	Conc_solve_cg(iconc, data_bag);
    implicit_loops_count = 1;

    while (implicit_loops_count < params->N_impheating_loops)  {

    	Lagrangian_heat(iconc, -1, data_bag, DTRACE("Lagrangian_heat") );

    	Communication_update_ghost_nodes_flow_variable(c[iconc]->data, 'c',params->ghost_nodes , data_bag);
    	Conc_set_boundary_values(c[iconc]->data, iconc ,CENTRAL_FULL ,grid , params);

    	Conc_solve_cg(iconc, data_bag);

    	implicit_loops_count++;
    }


	Lagrangian_heat(iconc,params->N_heating_loops, data_bag, DTRACE("Lagrangian_heat")); // changes directly c[iconc->data]
#ifndef FAST
	Communication_update_ghost_nodes_flow_variable(c[iconc]->data, 'c',params->ghost_nodes , data_bag);
	Conc_set_boundary_values(c[iconc]->data, iconc ,CENTRAL_FULL ,grid , params);
#endif



	Tend = MPI_Wtime();
	timer->Wtime_particle_total += Tend - Tstart;
#endif  // IBM_SCALAR


#ifdef IMMERSED_BOUNDARY
	 	Immersed_conc_interpolation(grid, params, c[iconc]);
		Immersed_conc_solid(c[iconc], data_bag);
#endif
		T2 = MPI_Wtime();
		timer->Wtime_c_solve += T2 - T1;

		//----------------------------------------------------------------------
		// Solve for concentration out of bounds
		//----------------------------------------------------------------------
#ifndef CONC_BQUICK

		int W_outofbounds;
		double tol = 1e-5;

		T1 = MPI_Wtime();
		W_outofbounds = Conc_calc_just_outofbounds(iconc, tol, data_bag);
		T2 = MPI_Wtime();
		timer->Wtime_c_outofbounds += T2 - T1;

#else // CONC_BQUICK

		T1 = MPI_Wtime();
		Conc_set_boundary_values_fortemp(c[iconc], grid, params);

		Communication_update_ghost_nodes_flow_variable(c[iconc]->data_temp, 'c', pnodes, data_bag);

		int W_outofbounds;
		double tol = 1e-5;
		int bquick_iter = 0;
		int tot_interval = 4;
		int old_outofbounds = -1;
		Array_set_noghost(c[iconc]->blend, 0.0, grid, params);

		W_outofbounds = Conc_calc_outofbounds(iconc, old_outofbounds,
											  &tot_interval, tol, data_bag);
		old_outofbounds = W_outofbounds;
		T2 = MPI_Wtime();
		timer->Wtime_c_outofbounds += T2 - T1;

		while (W_outofbounds != 0) {
			T1 = MPI_Wtime();
			Conc_set_conv_outofbounds(iconc, data_bag);
			T2 = MPI_Wtime();
			timer->Wtime_c_convective += T2 - T1;

			//---------
			// Set RHS
			//---------
			T1 = MPI_Wtime();
			Conc_set_RHS(c[iconc], grid, params);
			T2 = MPI_Wtime();
			timer->Wtime_c_rhs += T2 - T1;

			//-------------------------
			// Solve for concentration
			//-------------------------
			T1 = MPI_Wtime();


#ifdef  CONC_SEMI_IMPLICIT
		Conc_solve_semi_implicit(iconc, data_bag);
	#elif  defined CONC_FULLY_EXPLICIT
		Conc_solve_fully_explicit(iconc, data_bag);
	#elif defined CONC_SOLVE_CG
		Conc_solve_cg(iconc, data_bag);
#endif


	#ifdef IMMERSED_BOUNDARY
			Immersed_conc_interpolation(grid, params, c[iconc]);
			Immersed_conc_solid(c[iconc], data_bag);
	#endif
			T2 = MPI_Wtime();
			timer->Wtime_c_solve += T2 - T1;

			T1 = MPI_Wtime();
			Conc_set_boundary_values_fortemp(c[iconc], grid, params);
			Communication_update_ghost_nodes_flow_variable(c[iconc]->data_temp,
														   'c', pnodes, data_bag);
			tol = tol * 2.0;
			tol = min(tol,1e-2);
			W_outofbounds = Conc_calc_outofbounds(iconc, old_outofbounds,
												  &tot_interval, tol, data_bag);
			old_outofbounds = W_outofbounds;
			if ( (bquick_iter == 10) && (W_outofbounds != 0) ) {

				if (W_outofbounds > 10) {
					Conc_show_outofbounds(iconc, old_outofbounds, tol, data_bag);
					Display_throw_error("Temporal_int_all_the_equations()",
							"BQUICK iterations did not converge", params);
				}
				Conc_average_outofbounds(c[iconc], grid, params);
				Communication_update_ghost_nodes_flow_variable(c[iconc]->data_temp,
						'c', pnodes, data_bag);
				W_outofbounds = 0;

			}
			bquick_iter += 1;
			T2 = MPI_Wtime();
			timer->Wtime_c_outofbounds += T2 - T1;
		} // while W_outofbounds

#endif // CONC_BQUICK

		T1 = MPI_Wtime();
#ifdef CONC_BQUICK
		Conc_copy_fromtemp(c[iconc],grid,params);
#endif

#ifdef CONC_CLIP
		//----------------------------------------------------------------------
		// Clip concentration values to be between 0 and 1
		//----------------------------------------------------------------------
		for (k = grid->G_Ks; k < grid->G_Ke; k++) {
			for (j = grid->G_Js; j < grid->G_Je; j++) {
				for (i = grid->G_Is; i < grid->G_Ie; i++) {
					c[iconc]->data[k][j][i] =  max(c[iconc]->data[k][j][i],0.0);
					c[iconc]->data[k][j][i] =  min(c[iconc]->data[k][j][i],1.0);
				}
			}
		}
#endif // CONC_CLIP

		Conc_set_boundary_values(c[iconc]->data, iconc ,CENTRAL_FULL, grid , params); //  this should be dispensable

		T2 = MPI_Wtime();
		timer->Wtime_c_outofbounds += T2 - T1;

	} // for iconc


	Tend = MPI_Wtime();
	timer->Wtime_c_total += Tend - Tstart;
}




/******************************************************************************/
/*
 This function returns the linear system index of the conc. node, or its
 neighbors, depending on the value of parameter 'which'.

 which:
     = 'p'  - return lsys_index of the point itself.
     = 'w'	 - return lsys_index of west neighbour.  (x-direction)
     = 'e'  - return lsys_index of east neighbour.  (x+direction)
     = 's'  - return lsys_index of south neighbour. (y-direction)
     = 'n'  - return lsys_index of north neighbour. (y+direction)
     = 'b'  - return lsys_index of back neighbour.  (z-direction)
     = 'f'  - return lsys_index of front neighbour. (z+direction)
 */
/******************************************************************************/
/*
int Conc_get_lsys_index(Concentration *c, MAC_grid *grid, int x_index,
		int y_index, int z_index, char which) {

	int p_index;

	int NX = grid->NX;
	int NY = grid->NY;
	int NZ = grid->NZ;

	// Index of the point itself
	p_index = x_index + y_index * NX + z_index * (NX * NY);

	switch (which) {

		case 'p':
			return (p_index);
		case 'w':
			return (p_index - 1);
		case 'e':
			return (p_index + 1);
		case 's':
			return (p_index - NX);
		case 'n':
			return (p_index + NX);
		case 'b':
			return (p_index - NX*NY);
		case 'f':
			return (p_index + NX*NY);
		default:
			return -2;
	} // switch
}
*/




/******************************************************************************/
/*
 This Function, calculates the settling speed due to hindered settling effect.

 Various correlations could be used:

     1- Batchelers' model: 1977 (Theoretical equation)
         V_s / V_s0 = (1.0 - 6.55 * phi ) :: Valid for phi <= 0.02

     2- Ham, Homsy: 1988 (experimental relation)
         V_s / V_s0 = (1.0 - 4.0 * phi + 8 * phi^2) :: Valid for phi <= 0.10

     3- Zaki, Richardson type: 1954 (experimental relation)
         V_s / V_s0 = (1.0 - phi)^n :: Valid for phi <= 0.60 (To be verified) ::
         2.39 <= n <= 5.  Typical value n = 4.
 */
/******************************************************************************/
double Conc_settling_speed_function(double conc, double phi_max, double V_s0) {

	double phi = conc * phi_max; /* particle volume fraction */
	double n   = 4.0;

//	return ( V_s0 * ( 1.0 - 6.55 * phi) );
//	return ( V_s0 * ( 1.0 - 4.0 * phi + 8.0 * phi * phi ) );
	return ( V_s0 * pow( 1.0 - phi , n) );
}





/******************************************************************************/
/*
 */
/******************************************************************************/
void Conc_set_conv_viscous_quick(int iconc, Cart3d_bag *data_bag) {

	int i, j, k;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

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

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;

	double *aeW  = c -> aeW;
	double *aeE  = c -> aeE;
	double *aeEE = c -> aeEE;
	double *awWW = c -> awWW;
	double *awW  = c -> awW;
	double *awE  = c -> awE;
	double *anS  = c -> anS;
	double *anN  = c -> anN;
	double *anNN = c -> anNN;
	double *asSS = c -> asSS;
	double *asS  = c -> asS;
	double *asN  = c -> asN;
	double *afB  = c -> afB;
	double *afF  = c -> afF;
	double *afFF = c -> afFF;
	double *abBB = c -> abBB;
	double *abB  = c -> abB;
	double *abF  = c -> abF;

	// Get the local velocities at the location where they are defined
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***c_data = c -> data;

#ifdef VAR_VISC
	#ifdef LES
	double ***nut = data_bag -> smag -> cdev[iconc] -> mSct;
	#elif defined RANS
	double ***nut = data_bag -> rans -> nutc;
	#endif
#endif

	double ***conv = c -> ng_conv;
	double ***visc = c -> ng_viscous;
	double ***conv_old   = c -> ng_conv_old;
	double ***conv_quick = c -> ng_conv_quick;
	double ***visc_explicit = c -> ng_visc_explicit;

	// Constant settling speed of particle. Zero for temperature and Salinity
	double V_s0 = c -> v_settl0;

	if ( (params->which_stage == 1) || (params->which_stage == 2) ) {
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					conv_old[k][j][i] = conv[k][j][i] - visc_explicit[k][j][i];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------
	double iPe = 1.0 / c->Pe;
	double nuE = iPe;
	double nuW = iPe;
	double nuN = iPe;
	double nuS = iPe;
	double nuF = iPe;
	double nuB = iPe;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				/*------------------------------------------------------------*/
				/*
				 Calculate viscosity
				 */
				/*------------------------------------------------------------*/
#ifdef VAR_VISC
				//--------------------------------------------------------------
				// Variable viscosity - eddy viscosity for LES and RANS
				//--------------------------------------------------------------
				nuE = 0.5 * ( nut[k][j][i+1] + nut[k][j][i] );
				nuW = 0.5 * ( nut[k][j][i-1] + nut[k][j][i] );

				nuN = 0.5 * ( nut[k][j+1][i] + nut[k][j][i] );
				nuS = 0.5 * ( nut[k][j-1][i] + nut[k][j][i] );

				nuF = 0.5 * ( nut[k+1][j][i] + nut[k][j][i] );
				nuB = 0.5 * ( nut[k-1][j][i] + nut[k][j][i] );

	#if defined LES || defined RANS
				//--------------------------------------------------------------
				// LES and RANS 'nu' don't account for Peclet number
				//--------------------------------------------------------------
				nuE += iPe;
				nuW += iPe;
				nuN += iPe;
				nuS += iPe;
				nuF += iPe;
				nuB += iPe;
	#endif // LES or RANS
#endif // VAR_VISC


				/*------------------------------------------------------------*/
				/*
				 Calculate viscous terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// dc/dx East/West
				//--------------------------------------------------------------
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * idx_c[i];

				if (i != 0 ) {
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i-1];
				}
				else {
#if defined XPERIODIC || defined INFLOW
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i];
#else
					dcdxW = 0.;
#endif
				}

				//--------------------------------------------------------------
				// dc/dy North/South
				//--------------------------------------------------------------
				dcdyN = ( c_data[k][j+1][i] - c_data[k][j][i] ) * idy_c[j];

				if ( j != 0) {
					dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] ) * idy_c[j-1];
				}
				else {
					// Neumann b/c dcdy = 0
					dcdyS = 0.0;
#ifdef THERMAL_KADER
					dcdyS = 2.*c_data[k][j][i] * idy_v[j];
#endif
				}

				//--------------------------------------------------------------
				// dc/dz Front/Back
				//--------------------------------------------------------------
				dcdzF = ( c_data[k+1][j][i]-c_data[k][j][i] ) * idz_c[k];
				if ( k != 0) {
					dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) * idz_c[k-1];
				}
				else {
#ifdef ZPERIODIC
					dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) * idz_c[k];
#else
					// Neumann b/c dcdz = 0
					dcdzB = 0.;
#endif

				}

				//--------------------------------------------------------------
				// d2c/dx2, d2c/dy2 and d2c/dz2
				//--------------------------------------------------------------
				d2cdx2 = ( nuE * dcdxE - nuW * dcdxW ) * idx_u[i];
				d2cdy2 = ( nuN * dcdyN - nuS * dcdyS ) * idy_v[j];
				d2cdz2 = ( nuF * dcdzF - nuB * dcdzB ) * idz_w[k];


				/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// Convective velocities on the faces
				//--------------------------------------------------------------
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];

				if (grid->c_status[k][j][i] != SOLID) {

						if (j != NY-2) vN = vN + V_s0;
						vS = vS + V_s0;

				}

				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];

				//--------------------------------------------------------------
				// Concentration value on the faces
				// For explanation of various "if conditions", check cN and cS
				// calculation, which is representative of cE and cW, cF and cB
				//--------------------------------------------------------------
				// Calculate  cE
				if (uE > 0) {
					if (i==0) {
#ifdef XPERIODIC
						cE = awWW[i]*c_data[k][j][i-1] + awW[i]*c_data[k][j][i] + awE[i]*c_data[k][j][i+1];
#else
						cE = c_data[k][j][i];
#endif
					}
					else {
						cE = awWW[i]*c_data[k][j][i-1] + awW[i]*c_data[k][j][i] + awE[i]*c_data[k][j][i+1];
					}
				}
				else {
					if (i==NX-2) {
#ifdef XPERIODIC
						cE = aeW[i]*c_data[k][j][i] + aeE[i]*c_data[k][j][i+1] + aeEE[i]*c_data[k][j][i+2];
#else
						cE = c_data[k][j][i+1];
#endif
					}
					else if ( (grid->c_status[k][j][i+2] == FLUID) || (grid->c_status[k][j][i+2] == IMMERSED) ) {
						cE = aeW[i]*c_data[k][j][i] + aeE[i]*c_data[k][j][i+1] + aeEE[i]*c_data[k][j][i+2];
					}
					else {
						cE = c_data[k][j][i+1];
					}
				}

				// Calculate  cW
				if (i==0) {
#ifdef XPERIODIC
					if (uW > 0) {
						cW = awWW[0]*c_data[k][j][i-2] + awW[0]*c_data[k][j][i-1] + awE[0]*c_data[k][j][i];
					}
					else {
						cW = aeW[0]*c_data[k][j][i-1] + aeE[0]*c_data[k][j][i] + aeEE[0]*c_data[k][j][i+1];
					}
#else // not XPERIODIC
	#ifdef INFLOW
					cW = c_data[k][j][-1];
	#else
					//Otherwise impose Neumann b/c i.e dc/dx = 0.0
					cW = c_data[k][j][i];
	#endif
#endif // not XPERIODIC
				}

				else {
					if (uW > 0) {
						if (i==1) {
#ifdef XPERIODIC
							cW = awWW[i-1]*c_data[k][j][i-2] + awW[i-1]*c_data[k][j][i-1] + awE[i-1]*c_data[k][j][i];
#else
							cW = c_data[k][j][i-1];
#endif
						}
						else if ( (grid->c_status[k][j][i-2] == FLUID) || (grid->c_status[k][j][i-2] == IMMERSED) ) {
							cW = awWW[i-1]*c_data[k][j][i-2] + awW[i-1]*c_data[k][j][i-1] + awE[i-1]*c_data[k][j][i];
						}
						else {
							cW = c_data[k][j][i-1];
						}
					}
					else {
						cW = aeW[i-1]*c_data[k][j][i-1] + aeE[i-1]*c_data[k][j][i] + aeEE[i-1]*c_data[k][j][i+1];
					}
				}

				if (vN > 0) {
					if (j==0) {
						// Near the boundary use first order upwind
						cN = c_data[k][j][i];
					}
					else {
						cN = asSS[j]*c_data[k][j-1][i] + asS[j]*c_data[k][j][i] + asN[j]*c_data[k][j+1][i];
					}
				}
				else {
					if (j == NY-2) {
						// Near the boundary use first order upwind
						cN = c_data[k][j+1][i];
					}
					else if ( (grid->c_status[k][j+2][i] == FLUID) || (grid->c_status[k][j+2][i] == IMMERSED) ) {
						// Use QUICK when the stencil  permits (i.e all the
						// points have physical value)
						cN = anS[j]*c_data[k][j][i] + anN[j]*c_data[k][j+1][i] + anNN[j]*c_data[k][j+2][i];
					}
					else {
						// If QUICK stencil has solid nodes, use first order
						// upwind
						cN = c_data[k][j+1][i];
					}
				}

				if (vS > 0) {
					if (j==0) {
						// Impose Neumann b/c i.e dc/dy = 0.0
						cS = c_data[k][j][i];
					}
					else if (j==1) {
						cS = c_data[k][j-1][i];
					}
					else  if ( (grid->c_status[k][j-2][i] == FLUID) || (grid->c_status[k][j-2][i] == IMMERSED) ) {
						// Use QUICK when the stencil  permits (i.e all the
						// points have physical value)
						cS = asSS[j-1]*c_data[k][j-2][i] + asS[j-1]*c_data[k][j-1][i] + asN[j-1]*c_data[k][j][i];
					}
					else {
						// If QUICK stencil has solid nodes, use first order
						// upwind
						cS = c_data[k][j-1][i];
					}
				}
				else {
					if (j==0) {
						// Impose Neumann b/c i.e dc/dy = 0.0
						cS = c_data[k][j][i];
					}
					else {
						cS = anS[j-1]*c_data[k][j-1][i] + anN[j-1]*c_data[k][j][i] + anNN[j-1]*c_data[k][j+1][i];
					}
				}

				// Calculate  cF
				if (wF > 0) {
					if (k==0) {
#ifdef ZPERIODIC
						cF = abBB[k]*c_data[k-1][j][i] + abB[k]*c_data[k][j][i] + abF[k]*c_data[k+1][j][i];
#else
						cF = c_data[k][j][i];
#endif
					}
					else {
						cF = abBB[k]*c_data[k-1][j][i] + abB[k]*c_data[k][j][i] + abF[k]*c_data[k+1][j][i];
					}
				}
				else {
					if (k==NZ-2) {
#ifdef ZPERIODIC
						cF = afB[k]*c_data[k][j][i] + afF[k]*c_data[k+1][j][i] + afFF[k]*c_data[k+2][j][i];
#else
						cF = c_data[k+1][j][i];
#endif
					}
					else if ( (grid->c_status[k+2][j][i] == FLUID) || (grid->c_status[k+2][j][i] == IMMERSED) ) {
						cF = afB[k]*c_data[k][j][i] + afF[k]*c_data[k+1][j][i] + afFF[k]*c_data[k+2][j][i];
					}
					else {
						cF = c_data[k+1][j][i];
					}
				}

				// Calculate  cB
				if (k==0) {
#ifdef ZPERIODIC
					if (wB > 0) {
						cB = abBB[0]*c_data[k-2][j][i] + abB[0]*c_data[k-1][j][i] + abF[0]*c_data[k][j][i];
					}
					else {
						cB = afB[0]*c_data[k-1][j][i] + afF[0]*c_data[k][j][i] + afFF[0]*c_data[k+1][j][i];
					}
#else
					// Otherwise impose Neumann b/c i.e dc/dx = 0.0
					cB = c_data[k][j][i];
#endif
				}

				else {

					if (wB > 0) {
						if (k==1) {
#ifdef ZPERIODIC
							cB = abBB[k-1]*c_data[k-2][j][i] + abB[k-1]*c_data[k-1][j][i] + abF[k-1]*c_data[k][j][i];
#else
							cB = c_data[k-1][j][i];
#endif
						}
						else if ( (grid->c_status[k-2][j][i] == FLUID) || (grid->c_status[k-2][j][i] == IMMERSED) ) {
							cB = abBB[k-1]*c_data[k-2][j][i] + abB[k-1]*c_data[k-1][j][i] + abF[k-1]*c_data[k][j][i];
						}
						else {
							cB = c_data[k-1][j][i];
						}
					}
					else {
						cB = afB[k-1]*c_data[k-1][j][i] + afF[k-1]*c_data[k][j][i] + afFF[k-1]*c_data[k+1][j][i];
					}
				}

				ucE = uE*cE;
				ucW = uW*cW;

				vcN = vN*cN;
				vcS = vS*cS;

				wcF = wF*cF;
				wcB = wB*cB;

				//--------------------------------------------------------------
				// d/dx(uc), d/dy(vc), and d/dz(wc)
				//--------------------------------------------------------------
				dcudx = ( ucE - ucW ) * idx_u[i];
				dcvdy = ( vcN - vcS ) * idy_v[j];
				dcwdz = ( wcF - wcB ) * idz_w[k];


				/*------------------------------------------------------------*/
				/*
				 Store convective and viscous terms
				 */
				/*------------------------------------------------------------*/
				conv[k][j][i] = dcudx + dcvdy + dcwdz;
#ifdef CONC_FULLY_EXPLICIT
				visc[k][j][i] = 0.0;
				visc_explicit[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
#elif defined CONC_FULLY_IMPLICIT
				visc[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
				visc_explicit[k][j][i] = 0.0;
#elif defined CONC_SEMI_IMPLICIT
				visc[k][j][i] = d2cdy2;
				visc_explicit[k][j][i] = d2cdx2 + d2cdz2;
#else
				#error choose a diffusion scheme
#endif

#ifdef BQUICK
				conv_quick[k][j][i] = conv[k][j][i];
#endif

			}
		}
	}

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Conc_set_conv_viscous_ftupwind(int iconc, Cart3d_bag *data_bag) {

	int i, j, k;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;

	double xWW, xW, xE, xEE, xface;
	double ySS, yS, yN, yNN, yface;
	double zBB, zB, zF, zFF, zface;
	FILE *fid;
	char bin_filename[50];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

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

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;
	double *xc = grid -> xc;
	double *yc = grid -> yc;
	double *zc = grid -> zc;
	double *xu = grid -> xu;
	double *yv = grid -> yv;
	double *zw = grid -> zw;

	// Get the local velocities at the location where they are defined
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***c_data = c -> data;

#ifdef VAR_VISC
	#ifdef LES
	double ***nut = data_bag -> smag -> cdev[iconc] -> mSct;
	#elif defined RANS
	double ***nut = data_bag -> rans -> nutc;
	#endif
#endif

	double ***conv = c -> ng_conv;
	double ***visc = c -> ng_viscous;
	double ***conv_old = c -> ng_conv_old;
	double ***visc_explicit = c -> ng_visc_explicit;

	// Constant settling speed of particle. Zero for temperature and Salinity
	double V_s0 = c -> v_settl0;

	if ( (params->which_stage == 1) || (params->which_stage == 2)) {
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					conv_old[k][j][i] = conv[k][j][i] - visc_explicit[k][j][i];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------
	double iPe = 1.0 / c->Pe;
	double nuE = iPe;
	double nuW = iPe;
	double nuN = iPe;
	double nuS = iPe;
	double nuF = iPe;
	double nuB = iPe;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				/*------------------------------------------------------------*/
				/*
				 Calculate viscosity
				 */
				/*------------------------------------------------------------*/
#ifdef VAR_VISC
				//--------------------------------------------------------------
				// Variable viscosity - eddy viscosity for LES and RANS
				//--------------------------------------------------------------
				nuE = 0.5 * ( nut[k][j][i+1] + nut[k][j][i] );
				nuW = 0.5 * ( nut[k][j][i-1] + nut[k][j][i] );

				nuN = 0.5 * ( nut[k][j+1][i] + nut[k][j][i] );
				nuS = 0.5 * ( nut[k][j-1][i] + nut[k][j][i] );

				nuF = 0.5 * ( nut[k+1][j][i] + nut[k][j][i] );
				nuB = 0.5 * ( nut[k-1][j][i] + nut[k][j][i] );

	#if defined LES || defined RANS
				//--------------------------------------------------------------
				// LES and RANS 'nu' don't account for Peclet number
				//--------------------------------------------------------------
				nuE += iPe;
				nuW += iPe;
				nuN += iPe;
				nuS += iPe;
				nuF += iPe;
				nuB += iPe;
	#endif // LES or RANS
#endif // VAR_VISC


				/*------------------------------------------------------------*/
				/*
				 Calculate viscous terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// dc/dx East/West
				//--------------------------------------------------------------
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * idx_c[i];
				if (i != 0 ) {
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i-1];
				}
				else {
#if defined XPERIODIC || defined INFLOW
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i];
#else
					dcdxW = 0.;
#endif
				}

				//--------------------------------------------------------------
				// dc/dy North/South
				//--------------------------------------------------------------
				dcdyN = ( c_data[k][j+1][i] - c_data[k][j][i] ) * idy_c[j];
				if ( j != 0) {
					dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] ) * idy_c[j-1];
				}
				else {
					// Neumann b/c dcdy = 0
					dcdyS = 0.0;
#ifdef THERMAL_KADER
					dcdyS = 2.*c_data[k][j][i] * idy_v[j];
#endif
				}

				//--------------------------------------------------------------
				// dc/dz Front/Back
				//--------------------------------------------------------------
				dcdzF = ( c_data[k+1][j][i]-c_data[k][j][i] ) * idz_c[k];
				if ( k != 0) {
					dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) * idz_c[k-1];
				}
				else {
#ifdef ZPERIODIC
					dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) * idz_c[k];
#else
					// Neumann b/c dcdz = 0
					dcdzB = 0.;
#endif

				}

				//--------------------------------------------------------------
				// d2c/dx2, d2c/dy2 and d2c/dz2
				//--------------------------------------------------------------
				d2cdx2 = ( nuE * dcdxE - nuW * dcdxW ) * idx_u[i];
				d2cdy2 = ( nuN * dcdyN - nuS * dcdyS ) * idy_v[j];
				d2cdz2 = ( nuF * dcdzF - nuB * dcdzB ) * idz_w[k];


				/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// Convective velocities on the faces
				//--------------------------------------------------------------
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];


					if (j != NY-2) vN = vN + V_s0;
					vS = vS + V_s0;


				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];

				//--------------------------------------------------------------
				// Concentration value on the faces
				// For explanation of various "if conditions", check cN and cS
				// calculation, which is representative of cE and cW, cF and cB
				//--------------------------------------------------------------
				// Calculate  cE
				if (uE > 0) {
					cE = c_data[k][j][i];
				}
				else {
					cE = c_data[k][j][i+1];
				}

				// Calculate  cW
				if (i==0) {
#ifdef INFLOW
					cW = c_data[k][j][-1];
#else
					// Otherwise impose Neumann b/c i.e dc/dx = 0.0
					cW = c_data[k][j][i];
#endif
				}

				else {
					if (uW > 0) {
						cW = c_data[k][j][i-1];
					}
					else {
						cW = c_data[k][j][i];
					}
				}

				if (vN > 0) {
					// Near the boundary use first order upwind
					cN = c_data[k][j][i];
				}
				else {
					cN = c_data[k][j+1][i];
				}

				if (vS > 0) {
					if (j==0) {
						// Impose Neumann b/c i.e dc/dy = 0.0
						cS = c_data[k][j][i];
					}
					else {
						// If QUICK stencil has solid nodes, use first order
						// upwind
						cS = c_data[k][j-1][i];
					}
				}
				else {
					cS = c_data[k][j][i];
				}

				// Calculate  cF
				if (wF > 0) {
					cF = c_data[k][j][i];
				}
				else {
					cF = c_data[k+1][j][i];
				}

				// Calculate  cB
				if (k==0) {
					// Otherwise impose Neumann b/c i.e dc/dx = 0.0
					cB = c_data[k][j][i];
				}

				else {

					if (wB > 0) {
						cB = c_data[k-1][j][i];
					}
					else {
						cB = c_data[k][j][i];
					}
				}

				ucE = uE*cE;
				ucW = uW*cW;

				vcN = vN*cN;
				vcS = vS*cS;

				wcF = wF*cF;
				wcB = wB*cB;

				//--------------------------------------------------------------
				// d/dx(uc), d/dy(vc), and d/dz(wc)
				//--------------------------------------------------------------
				dcudx = ( ucE - ucW ) * idx_u[i];
				dcvdy = ( vcN - vcS ) * idy_v[j];
				dcwdz = ( wcF - wcB ) * idz_w[k];

				/*------------------------------------------------------------*/
				/*
				 Store convective and viscous terms
				 */
				/*------------------------------------------------------------*/
				conv[k][j][i] = dcudx + dcvdy + dcwdz;
				#ifdef CONC_FULLY_EXPLICIT
								visc[k][j][i] = 0.0;
								visc_explicit[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
				#elif defined CONC_FULLY_IMPLICIT
								visc[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
								visc_explicit[k][j][i] = 0.0;
				#elif defined CONC_SEMI_IMPLICIT
								visc[k][j][i] = d2cdy2;
								visc_explicit[k][j][i] = d2cdx2 + d2cdz2;
				#else
								#error  choose a diffusion scheme
				#endif
			}
		}
	}

	return;
}




// #ifdef CONC_BQUICK
/******************************************************************************/
/*
 */
/******************************************************************************/
void Conc_set_conv_outofbounds(int iconc, Cart3d_bag *data_bag) {

	int i, j, k;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;

	double xWW, xW, xE, xEE, xface;
	double ySS, yS, yN, yNN, yface;
	double zBB, zB, zF, zFF, zface;
	FILE *fid;
	char bin_filename[50];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

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

	double *idx_u = grid->idx_u;
	double *idy_v = grid->idy_v;
	double *idz_w = grid->idz_w;
	double *idx_c = grid->idx_c;
	double *idy_c = grid->idy_c;
	double *idz_c = grid->idz_c;
	double *xc = grid->xc;
	double *yc = grid->yc;
	double *zc = grid->zc;
	double *xu = grid->xu;
	double *yv = grid->yv;
	double *zw = grid->zw;

	// Get the local velocities at the location where they are defined
	double ***u_data = data_bag->u->data;
	double ***v_data = data_bag->v->data;
	double ***w_data = data_bag->w->data;
	double ***c_data = c->data;
	double ***c_blend = c->blend;
	double ***c_data_temp = c->data_temp;

	double ***conv = c->ng_conv;
	double ***conv_quick = c->ng_conv_quick;

	// Constant settling speed of particle. Zero for temperature and Salinity
	double V_s0 = c->v_settl0;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				// Check if concentration value is out of bounds
//				if ( ((c_data_temp[k][j][i]-1.) > 1e-8) || ((c_data_temp[k][j][i]+1.) < 1e-8) ) {
				if (c_blend[k][j][i] > 1e-8) {

					// Calculate the convective terms
					// Convective velocity on the faces
					uE = u_data[k][j][i+1];
					uW = u_data[k][j][i];

					vN = v_data[k][j+1][i];
					vS = v_data[k][j][i];


						if (j != NY-2) vN = vN + V_s0;
						vS = vS + V_s0;


					wF = w_data[k+1][j][i];
					wB = w_data[k][j][i];

					// Concentration value on the faces
					// For explanation of various "if conditions", check cN and
					// cS calculation, which is representative of cE and cW, cF
					// and cB

					// Calculate  cE
					if (uE > 0)
						cE = c_data[k][j][i];
					else
						cE = c_data[k][j][i+1];

					// Calculate  cW
					if (i==0) {
	#ifdef INFLOW
						cW = c_data[k][j][-1];
	#else
						// Otherwise impose Neumann b/c i.e dc/dx = 0.0
						cW = c_data[k][j][i];
	#endif
					}
					else {
						if (uW > 0)
							cW = c_data[k][j][i-1];
						else
							cW = c_data[k][j][i];
					}

					if (vN > 0)
						// Near the boundary use first order upwind
						cN = c_data[k][j][i];
					else
						cN = c_data[k][j+1][i];

					if (vS > 0) {
						if (j==0)
							// Impose Neumann b/c i.e dc/dy = 0.0
							cS = c_data[k][j][i];
						else
							cS = c_data[k][j-1][i];
					}
					else {
						cS = c_data[k][j][i];
					}

					// Calculate  cF
					if (wF > 0)
						cF = c_data[k][j][i];
					else
						cF = c_data[k+1][j][i];

					// Calculate  cB
					if (k==0)
						// Otherwise impose Neumann b/c i.e dc/dx = 0.0
						cB = c_data[k][j][i];
					else {
						if (wB > 0)
							cB = c_data[k-1][j][i];
						else
							cB = c_data[k][j][i];
					}

					ucE = uE*cE;
					ucW = uW*cW;

					vcN = vN*cN;
					vcS = vS*cS;

					wcF = wF*cF;
					wcB = wB*cB;

					dcudx = ( ucE - ucW ) * idx_u[i];
					dcvdy = ( vcN - vcS ) * idy_v[j];
					dcwdz = ( wcF - wcB ) * idz_w[k];

//					conv[k][j][i] = dcudx + dcvdy + dcwdz;
					conv[k][j][i] = c_blend[k][j][i]*(dcudx + dcvdy + dcwdz) + (1.-c_blend[k][j][i])*conv_quick[k][j][i];
				}

			}
		}
	}

	return;
}
// #endif // CONC_BQUICK




/******************************************************************************/
/*
 This function computes the RHS of the u-momentum linear system. Convective
 terms are treated explicitly
 */
/******************************************************************************/
void Conc_set_RHS(Concentration *c, MAC_grid *grid, Parameters *params) {

	int i, j, k;
	double an, as, ae, aw;
	double rhs_value;
	double value;

	// Pointer to a function
	double ***data = c -> data;

	// Now, got the RHS vector on current processor
	double ***rhs_vec = c -> ng_rhs;

	double ***conv = c -> ng_conv;
	double ***visc = c -> ng_viscous;
	double ***conv_old = c -> ng_conv_old;
	double ***visc_explicit = c -> ng_visc_explicit;

	int which_stage = params -> which_stage;

	const double BET[] = {BETA};
	const double GAMB[] = {GAMBETA};
	const double ZETB[] = {ZETBETA};
	double a_dt = 1.0 / (params -> dt * BET[which_stage]);

	double *idx_u = grid->idx_u;
	double *idx_c = grid->idx_c;
	double *idy_v = grid->idy_v;
	double *idz_w = grid->idz_w;



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

	// indices start and end on current processor
	int i_start = Is;
	int j_start = Js;
	int k_start = Ks;

	// exclude the half cell added
	int i_end   = min(NX-1, Ie);
	int j_end   = min(NY-1, Je);
	int k_end   = min(NZ-1, Ke);

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				rhs_vec[k][j][i] = data[k][j][i]*a_dt;
			}// for i
		} // for
	}// for k


	double *conv_1d          = &conv[Ks][Js][Is];
	double *conv_old_1d      = &conv_old[Ks][Js][Is];
	double *visc_1d          = &visc[Ks][Js][Is];
	double *rhs_1d           = &rhs_vec[Ks][Js][Is];
	double *visc_explicit_1d = &visc_explicit[Ks][Js][Is];

	int index = 0;
	int total_nodes = (Ie - Is) * (Je - Js) * (Ke - Ks);



	if (which_stage == 0) {
		for ( index=0; index < total_nodes; index++ ) {
#ifdef IBM_SCALAR
			rhs_1d[index] += GAMB[0]*(-conv_1d[index])+ 2*visc_1d[index];

#elif defined CONC_FULLY_EXPLICIT
			rhs_1d[index] += GAMB[0]*(visc_explicit_1d[index] - conv_1d[index]);

#elif defined CONC_SEMI_IMPLICIT
			rhs_1d[index] += GAMB[0]*(visc_explicit_1d[index] - conv_1d[index])+ visc_1d[index];

#elif defined CONC_FULLY_IMPLICIT
			rhs_1d[index] += GAMB[0]*(-conv_1d[index])+ visc_1d[index];

#else
			#error
#endif
		}
	}
	else if (which_stage == 1){

		for ( index=0; index < total_nodes; index++ ) {
#ifdef IBM_SCALAR

			rhs_1d[index] += GAMB[1]*(-conv_1d[index]) - ZETB[1]*conv_old_1d[index] + 2*visc_1d[index];

#elif defined CONC_FULLY_EXPLICIT
			rhs_1d[index] += GAMB[1]*(visc_explicit_1d[index] - conv_1d[index]) - ZETB[1]*conv_old_1d[index];
#elif defined CONC_SEMI_IMPLICIT
			rhs_1d[index] += GAMB[1]*(visc_explicit_1d[index] - conv_1d[index]) - ZETB[1]*conv_old_1d[index] + visc_1d[index];
#elif  defined CONC_FULLY_IMPLICIT
			rhs_1d[index] += GAMB[1]*(-conv_1d[index]) - ZETB[1]*conv_old_1d[index] + visc_1d[index];
#else
			#error
#endif
		}
	}
	else if (which_stage == 2){
		for ( index=0; index < total_nodes; index++ ) {
#ifdef IBM_SCALAR

			rhs_1d[index] += GAMB[2]*(-conv_1d[index]) - ZETB[2]*conv_old_1d[index] + 2*visc_1d[index];

#elif defined CONC_FULLY_EXPLICIT
			rhs_1d[index] += GAMB[2]*( visc_explicit_1d[index] - conv_1d[index]) - ZETB[2]*conv_old_1d[index];
#elif defined  CONC_SEMI_IMPLICIT
			rhs_1d[index] += GAMB[2]*(visc_explicit_1d[index] - conv_1d[index]) - ZETB[2]*conv_old_1d[index] + visc_1d[index];
#elif defined CONC_FULLY_IMPLICIT
			rhs_1d[index] += GAMB[2]*(-conv_1d[index]) - ZETB[2]*conv_old_1d[index] + visc_1d[index];
#else
			#error
#endif
		}
	}

#ifdef IMMERSED_BOUNDARY
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

				if (grid->c_status[k][j][i] == SOLID)
					rhs_vec[k][j][i] = 0.0;
			}
		}
	}
#endif


}

void Conc_set_RHS_IBM_implicit(Concentration *c, double ***temp_f ,MAC_grid *grid, Parameters *params) {

	int i, j, k;
	double an, as, ae, aw;
	double rhs_value;
	double value;

	// Pointer to a function
	double ***data = c -> data;

	// Now, got the RHS vector on current processor
	double ***rhs_vec = c -> ng_rhs;
	double ***visc = c -> ng_viscous;



	int which_stage = params -> which_stage;

	const double BET[] = {BETA};

	double a_dt = 1.0 / (params -> dt * BET[which_stage]);

	double *idx_u = grid->idx_u;
	double *idx_c = grid->idx_c;
	double *idy_v = grid->idy_v;
	double *idz_w = grid->idz_w;



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

	// indices start and end on current processor
	int i_start = Is;
	int j_start = Js;
	int k_start = Ks;

	// exclude the half cell added
	int i_end   = min(NX-1, Ie);
	int j_end   = min(NY-1, Je);
	int k_end   = min(NZ-1, Ke);

	double *visc_1d          = &visc[Ks][Js][Is];
	double *rhs_1d           = &rhs_vec[Ks][Js][Is];
    double *temp_f_1d		 = &temp_f[Ks][Js][Is];


	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				rhs_vec[k][j][i] = data[k][j][i]*a_dt;
			}// for i
		} // for j
	}// for k




	int index = 0;
	int total_nodes = (Ie - Is) * (Je - Js) * (Ke - Ks);


	for ( index=0; index < total_nodes; index++ ) {

			rhs_1d[index] +=  -visc_1d[index]; //+  temp_f_1d[index];  // add IBM forcing here

		}




#ifdef IMMERSED_BOUNDARY
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

				if (grid->c_status[k][j][i] == SOLID)
					rhs_vec[k][j][i] = 0.0;
			}
		}
	}
#endif

}


/******************************************************************************/
/*
 This function computes the RHS of the u-momentum linear system.  Convective
 terms are treated explicitly
 */
/******************************************************************************/
void Conc_add_source_RHS(int iconc, Cart3d_bag *data_bag) {

	int i, j, k;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	double ubulk = Velocity_ubulk(data_bag->u, data_bag);
	double c_source = 2.0/ubulk;

	// Now, got the RHS vector on current processor
	double ***rhs_vec   = data_bag -> c[iconc] -> ng_rhs;
	double ***u_data_bc = data_bag -> u -> data_bc;

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

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++){
				rhs_vec[k][j][i] += c_source * u_data_bc[k][j][i];
			}// for i
		} // for j
	}// for k

}




/******************************************************************************/
/*
 This function copies the value of c->data into c->data_old. This is done for
 the time integration purposes.
 */
/******************************************************************************/
void Conc_store_old_data(Concentration *c, MAC_grid *grid, Parameters *params) {

	// copy data into data_old
	Array_copy_noghost(c->data, c->ng_data_old, grid, params);
}




/******************************************************************************/
/*
 This function sets value of conc  at the boundaries
 */
/******************************************************************************/
void Conc_set_boundary_values(double ***data , int iconc ,int type, MAC_grid *grid , Parameters *params) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	// Get regular data array for u vel data
	double ***c_data = data;

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




double AN = params->BC_AN[iconc];  // this set A \partial_x c+ B c= C at the North Boundary
double AS = params->BC_AS[iconc];
double AE = params->BC_AE[iconc];
double AW = params->BC_AW[iconc];
double AF = params->BC_AF[iconc];
double AB = params->BC_AB[iconc];


double BN = params->BC_BN[iconc];  // this set A \partial_x c+ B c= C at the North Boundary
double BS = params->BC_BS[iconc];
double BE = params->BC_BE[iconc];
double BW = params->BC_BW[iconc];
double BF = params->BC_BF[iconc];
double BB = params->BC_BB[iconc];


double CN = params->BC_CN[iconc];  // this set A \partial_x c+ B c= C at the North Boundary
double CS = params->BC_CS[iconc];
double CE = params->BC_CE[iconc];
double CW = params->BC_CW[iconc];
double CF = params->BC_CF[iconc];
double CB = params->BC_CB[iconc];



if (type==CENTRAL_PERTURBATION){
CN=0;
CS=0;
CE=0;
CW=0;
CB=0;
CF=0;
}

double *idx_c = grid -> idx_c;
double *idy_c = grid -> idy_c;
double *idz_c = grid -> idz_c;

double idx_E= idx_c[Is];
double idx_W= idx_c[Ie-1];

//printf("The East boundayr gridsize is %f \n", idx_E);


double idy_S= idy_c[Js];
double idy_N= idy_c[Je-1];

double idz_B= idz_c[Ks];
double idz_F= idz_c[Ke-1];



#ifndef XPERIODIC_CONC
	if (Ie==NX) {  // East Face
		i = NX-1;
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				c_data[k][j][i] = ((AE*idx_E-BE/2.)*c_data[k][j][i-1]+CE)/(AE*idx_E+BE/2.) ;
			}
		}
	}

	if (Is==0) {  // West Face
			i = -1;
			for (k=Ks; k<Ke; k++) {
				for (j=Js; j<Je; j++) {
					c_data[k][j][i] = (-(AW*idx_W+BW/2.)*c_data[k][j][i+1]+CW)/(-AW*idx_W+BW/2.) ;

				}
			}
		}

#endif




	if (Je==NY) {  // North Face

		j = NY-1;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				c_data[k][j][i] = ((AN*idy_N-BN/2.)*c_data[k][j-1][i]+CN)/(AN*idy_N+BN/2.);

			}
		}
	}


	if (Js==0) {  // South Face
			j = -1;
			for (k=Ks; k<Ke; k++) {
				for (i=Is; i<Ie; i++) {
					c_data[k][j][i] = (-(AS/idy_S+BS/2.)*c_data[k][j+1][i]+CS)/(-AS/idy_S+BS/2.);

				}
			}
		}




#ifndef ZPERIODIC_CONC
	if (Ke==NZ) {		// Front
		k = NZ-1;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				c_data[k][j][i] = ((AF*idz_F-BF/2.)*c_data[k-1][j][i]+CF)/(AF*idz_F+BF/2.);

			}
		}
	}

	if (Ks==0) {  // Back
			k = -1;
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					c_data[k][j][i] = (-(AB*idz_B+BB/2.)*c_data[k+1][j][i]+CB)/(-AB*idz_B+BB/2.);

				}
			}
		}

#endif




}




/******************************************************************************/
/*
 This function sets value of   at the boundaries
 */
/******************************************************************************/
void Conc_set_boundary_values_fortemp(Concentration *c, MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ;
	int i, j, k;
	double ***c_data;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	c_data = c->data_temp;

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

	if (Je==NY) {
		j = NY-1;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				c_data[k][j][i] = max(c_data[k][j-1][i],0.0);
				c_data[k][j][i] = min(c_data[k][j][i],1.0);
			}
		}
	}

#ifdef FRONT_WALL_VELOCITY_FREESLIP
	if (Ke==NZ) {
			k = NZ-1;
			for (j=Js; j<Je; j++) {
					for (i=Is; i<Ie; i++) {
							c_data[k][j][i] = max(c_data[k-1][j][i],0.0);
							c_data[k][j][i] = min(c_data[k][j][i],1.0);
					}
			}
	}
#endif
#ifdef FRONT_WALL_VELOCITY_NOSLIP
	if (Ke==NZ) {
		k = NZ-1;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				c_data[k][j][i] = max(c_data[k-1][j][i],0.0);
				c_data[k][j][i] = min(c_data[k][j][i],1.0);
			}
		}
	}
#endif

#ifndef XPERIODIC
	if (Ie==NX) {
		i = NX-1;
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				c_data[k][j][i] = max(c_data[k][j][i-1],0.0);
				c_data[k][j][i] = min(c_data[k][j][i],1.0);
			}
		}
	}
#endif

}




/******************************************************************************/
/*
 This function calculates the number of nodes where scalar value is out of
 bounds and it also sets the blending factor that is used to blend QUICK and
 first order upwind
 */
/******************************************************************************/
int Conc_calc_just_outofbounds(int iconc, double tol, Cart3d_bag *data_bag) {

	int NX, NY, NZ;
	int i, j, k;
	double ***c_data;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int l_outofbounds, W_outofbounds;
	int n_outofbounds, p_outofbounds;
	double nc_value,pc_value;
	double nc_rms, pc_rms;
	double nc_min, pc_max;
	FILE *fid;
	char bin_filename[50];
	double sendc[7], recvc[7];
	int W_n_outofbounds, W_p_outofbounds;
	double W_nc_min, W_nc_rms, W_nc_value;
	double W_pc_max, W_pc_rms, W_pc_value;
	const double BET[] = {BETA};
	double actual_time;
	double ptol, ntol;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	int which_stage = params -> which_stage;

	// Get regular data array for u vel data
	c_data = data_bag -> c[iconc] -> data;

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

	if (Je == NY)
		Je=Je-1;
#ifndef XPERIODIC
	if (Ie == NX)
		Ie=Ie-1;
#endif
#ifndef ZPERIODIC
	if (Ke == NZ)
		Ke=Ke-1;
#endif

	l_outofbounds = 0;
	W_outofbounds = 0;
	n_outofbounds = 0;
	p_outofbounds = 0;
	nc_value = 0;
	pc_value = 0;
	nc_min = 1.0;
	pc_max = -1.0;
	nc_rms = 0.0;
	pc_rms = 0.0;
	ptol = 1.0 + tol;
	ntol = -tol;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				if ( c_data[k][j][i] > ptol) {
					p_outofbounds += 1;
					pc_value = pc_value + c_data[k][j][i];
					pc_rms = pc_rms + c_data[k][j][i]*c_data[k][j][i];
					pc_max = max(pc_max, c_data[k][j][i]);

				}
				if ( c_data[k][j][i] < ntol ) {
					n_outofbounds += 1;
					nc_value = nc_value + c_data[k][j][i];
					nc_rms = nc_rms + c_data[k][j][i]*c_data[k][j][i];
					nc_min = min(nc_min, c_data[k][j][i]);
				}
			}
		}
	}

// 	MPI_Allreduce(&l_outofbounds, &W_outofbounds, 1, MPI_INT, MPI_SUM, PCW);
	l_outofbounds = p_outofbounds + n_outofbounds;

	sendc[0] = l_outofbounds;
	sendc[1] = p_outofbounds;
	sendc[2] = n_outofbounds;
	sendc[3] = pc_value;
	sendc[4] = pc_rms;
	sendc[5] = nc_value;
	sendc[6] = nc_rms;
	MPI_Allreduce(&sendc[0], &recvc[0], 7, MPI_DOUBLE, MPI_SUM, PCW);
	W_outofbounds = recvc[0];
	W_p_outofbounds = recvc[1];
	W_n_outofbounds = recvc[2];
	W_pc_value = recvc[3];
	W_pc_rms   = recvc[4];
	W_nc_value = recvc[5];
	W_nc_rms   = recvc[6];

	W_nc_value = W_nc_value/(W_n_outofbounds + 1e-16);
	W_pc_value = W_pc_value/(W_p_outofbounds + 1e-16);
	W_nc_rms = W_nc_rms/(W_n_outofbounds + 1e-16)-W_nc_value*W_nc_value;
	W_pc_rms = W_pc_rms/(W_p_outofbounds + 1e-16)-W_pc_value*W_pc_value;
	W_nc_rms = sqrt(W_nc_rms);
	W_pc_rms = sqrt(W_pc_rms);
	MPI_Allreduce(&nc_min, &W_nc_min, 1, MPI_DOUBLE, MPI_MIN, PCW);
	MPI_Allreduce(&pc_max, &W_pc_max, 1, MPI_DOUBLE, MPI_MAX, PCW);

	if (params->rank==0) {
		sprintf(bin_filename,"outofbounds_stat%d.dat",iconc);
		fid = fopen(bin_filename,"a");
		if (which_stage == 0)
			actual_time = params->time + params->dt*(BET[0]);
		if (which_stage == 1)
			actual_time = params->time + params->dt*(BET[0] + BET[1]);
		if (which_stage == 2)
			actual_time = params->time + params->dt*(BET[0] + BET[1] + BET[2]);
		fprintf(fid, "%6d %4d %16.8e %16.8e %6d %6d %16.8e %16.8e %16.8e %6d %12.8e %12.8e %12.8e\n"
				,params->ntime,which_stage,params->time,actual_time,W_outofbounds,
				W_p_outofbounds,W_pc_value,W_pc_rms, W_pc_max, W_n_outofbounds,
				W_nc_value, W_nc_rms, W_nc_min);
		fclose(fid);
	}

	return W_outofbounds;

}




/******************************************************************************/
/*
 This function calculates the number of nodes where scalar value is out of
 bounds and it also sets the blending factor that is used to blend QUICK and
 first order upwind
 */
/******************************************************************************/
int Conc_calc_outofbounds(int iconc, int old_outofbounds, int* tot_interval,
		double tol, Cart3d_bag *data_bag) {

	int i, j, k;
	double ***c_data_temp, ***c_blend;
	int l_outofbounds, W_outofbounds;
	int n_outofbounds, p_outofbounds;
	double nc_value,pc_value;
	double nc_rms, pc_rms;
	double nc_min, pc_max;
	double factor;
	FILE *fid;
	char bin_filename[50];
	double sendc[7], recvc[7];
	int W_n_outofbounds, W_p_outofbounds;
	double W_nc_min, W_nc_rms, W_nc_value;
	double W_pc_max, W_pc_rms, W_pc_value;
	const double BET[] = {BETA};
	double actual_time;
	double ptol, ntol;
	int fac_interval;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

	int which_stage = params -> which_stage;

	// Get regular data array for u vel data
	c_data_temp = c->data_temp;
	c_blend = c->blend;

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

	l_outofbounds = 0;
	W_outofbounds = 0;
	n_outofbounds = 0;
	p_outofbounds = 0;
	nc_value = 0;
	pc_value = 0;
	nc_min = 1.0;
	pc_max = -1.0;
	nc_rms = 0.0;
	pc_rms = 0.0;
	ptol = 1.0 + tol;
	ntol = -tol;

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				if ( c_data_temp[k][j][i] > ptol) {
					p_outofbounds += 1;
					pc_value +=  c_data_temp[k][j][i];
					pc_rms +=  c_data_temp[k][j][i]*c_data_temp[k][j][i];
//					pc_max = max(pc_max, c_data_temp[k][j][i]);
					c_blend[k][j][i] = 1.0;

				}
				if ( c_data_temp[k][j][i] < ntol ) {
					n_outofbounds +=1;
					nc_value += c_data_temp[k][j][i];
					nc_rms += c_data_temp[k][j][i]*c_data_temp[k][j][i];
//					nc_min = min(nc_min, c_data_temp[k][j][i]);
					c_blend[k][j][i] = 1.0;
				}
			}
		}
	}

	l_outofbounds = p_outofbounds + n_outofbounds;

	sendc[0] = l_outofbounds;
	sendc[1] = p_outofbounds;
	sendc[2] = n_outofbounds;
	sendc[3] = pc_value;
	sendc[4] = pc_rms;
	sendc[5] = nc_value;
	sendc[6] = nc_rms;
	MPI_Allreduce(&sendc[0], &recvc[0], 7, MPI_DOUBLE, MPI_SUM, PCW);
	W_outofbounds = recvc[0];
	W_p_outofbounds = recvc[1];
	W_n_outofbounds = recvc[2];
	W_pc_value = recvc[3];
	W_pc_rms   = recvc[4];
	W_nc_value = recvc[5];
	W_nc_rms   = recvc[6];

	W_nc_value = W_nc_value/(W_n_outofbounds + 1e-16);
	W_pc_value = W_pc_value/(W_p_outofbounds + 1e-16);
	W_nc_rms = W_nc_rms/(W_n_outofbounds + 1e-16)-W_nc_value*W_nc_value;
	W_pc_rms = W_pc_rms/(W_p_outofbounds + 1e-16)-W_pc_value*W_pc_value;
	W_nc_rms = sqrt(W_nc_rms);
	W_pc_rms = sqrt(W_pc_rms);
	MPI_Allreduce(&nc_min, &W_nc_min, 1, MPI_DOUBLE, MPI_MIN, PCW);
	MPI_Allreduce(&pc_max, &W_pc_max, 1, MPI_DOUBLE, MPI_MAX, PCW);

	if (params->rank==0) {
		sprintf(bin_filename,"outofbounds_stat%d.dat",iconc);
		fid = fopen(bin_filename,"a");
		if (which_stage == 0)
			actual_time = params->time + params->dt*(BET[0]);
		if (which_stage == 1)
			actual_time = params->time + params->dt*(BET[0] + BET[1]);
		if (which_stage == 2)
			actual_time = params->time + params->dt*(BET[0] + BET[1] + BET[2]);
		fprintf(fid, "%6d %4d %16.8e %16.8e %6d %6d %16.8e %16.8e %16.8e %6d %12.8e %12.8e %12.8e\n",params->ntime,which_stage,
				params->time,actual_time,W_outofbounds,W_p_outofbounds,W_pc_value,
				W_pc_rms, W_pc_max, W_n_outofbounds, W_nc_value, W_nc_rms, W_nc_min);
		fclose(fid);
	}

	if (W_outofbounds == 0)
		return W_outofbounds;

	if (W_outofbounds == old_outofbounds )
		*tot_interval = *tot_interval + 1;

	if (params->rank==0)
		printf("tot_interval= %d\n",*tot_interval);

	fac_interval = *tot_interval;

	if (*tot_interval > 10) {
		fac_interval = 10;
		for (i=0;i<*tot_interval-fac_interval;i++) {
			factor = 1.0;
			Communication_update_ghost_nodes_flow_variable(c->blend, 'c', 1, data_bag);
			Conc_set_blendingfactor(c, grid, params, factor);
		}
	}

	for (i=0;i<fac_interval;i++) {
		factor = 1.0 - (i + 1.0) / (*tot_interval + 1);
		Communication_update_ghost_nodes_flow_variable(c->blend, 'c', 1, data_bag);
		Conc_set_blendingfactor(c, grid, params, factor);
	}

	return W_outofbounds;
}




/******************************************************************************/
/*
 This function calculates the number of nodes where scalar value is out of
 bounds and it also sets the blending factor that is used to blend QUICK and
 first order upwind
 */
/******************************************************************************/
void Conc_show_outofbounds(int iconc, int old_outofbounds, double tol,
		Cart3d_bag *data_bag) {

	int i, j, k;
	double ***c_data_temp, ***c_blend;
	int ***c_status;
	int l_outofbounds, W_outofbounds;
	int n_outofbounds, p_outofbounds;
	double nc_value,pc_value;
	double nc_rms, pc_rms;
	double nc_min, pc_max;
	double factor;
	FILE *fid, *fid1;
	char bin_filename[50];
	double sendc[7], recvc[7];
	int W_n_outofbounds, W_p_outofbounds;
	double W_nc_min, W_nc_rms, W_nc_value;
	double W_pc_max, W_pc_rms, W_pc_value;
	const double BET[] = {BETA};
	double actual_time;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

	// Get regular data array for u vel data
	c_blend = c->blend;
	c_data_temp = c->data_temp;
	c_status = grid->c_status;

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

	l_outofbounds = 0;
	W_outofbounds = 0;
	n_outofbounds = 0;
	p_outofbounds = 0;
	nc_value = 0;
	pc_value = 0;
	nc_min = 1.0;
	pc_max = -1.0;
	nc_rms = 0.0;
	pc_rms = 0.0;
	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				if ( (c_data_temp[k][j][i]-1.) > tol) {
					l_outofbounds += 1;
					p_outofbounds += 1;
					pc_value = pc_value + c_data_temp[k][j][i];
					pc_rms = pc_rms + c_data_temp[k][j][i]*c_data_temp[k][j][i];
					pc_max = max(pc_max, c_data_temp[k][j][i]);

					sprintf(bin_filename,"overbound%d.dat",params->rank);
					fid = fopen(bin_filename,"a");
					fprintf(fid,"%e %e %e %e %e \n",grid->xc[i],grid->yc[j],
							grid->zc[k],c_data_temp[k][j][i],c_blend[k][j][i]);
					if  ( (i>0) && (j>0) && (k>0) ) {
						fprintf(fid,"%e %e %e %e %e %e \n",
								c_data_temp[k][j][i-1],c_data_temp[k][j][i+1],
								c_data_temp[k][j-1][i],c_data_temp[k][j+1][i],
								c_data_temp[k-1][j][i],c_data_temp[k+1][j][i]);
						fprintf(fid,"%e %e %e %e %e %e \n",
								c_blend[k][j][i-1],c_blend[k][j][i+1],
								c_blend[k][j-1][i],c_blend[k][j+1][i],
								c_blend[k-1][j][i],c_blend[k+1][j][i]);
						fprintf(fid,"status %d %d %d %d %d %d %d \n",
								c_status[k][j][i],c_status[k][j][i-1],c_status[k][j][i+1],
								c_status[k][j-1][i],c_status[k][j+1][i],
								c_status[k-1][j][i],c_status[k+1][j][i]);
					}
					fprintf(fid,"\n");
					fclose(fid);
				}
				if ( c_data_temp[k][j][i] < -tol ) {
					l_outofbounds += 1;
					n_outofbounds +=1;
					nc_value = nc_value + c_data_temp[k][j][i];
					nc_rms = nc_rms + c_data_temp[k][j][i]*c_data_temp[k][j][i];
					nc_min = min(nc_min, c_data_temp[k][j][i]);

					sprintf(bin_filename,"underbound%d.dat",params->rank);
					fid1 = fopen(bin_filename,"a");
					fprintf(fid1,"%e %e %e %e %e \n",grid->xc[i],grid->yc[j],
							grid->zc[k],c_data_temp[k][j][i],c_blend[k][j][i]);
					if  ( (i>0) && (j>0) && (k>0) ) {
						fprintf(fid1,"%e %e %e %e %e %e \n",
								c_data_temp[k][j][i-1],c_data_temp[k][j][i+1],
								c_data_temp[k][j-1][i],c_data_temp[k][j+1][i],
								c_data_temp[k-1][j][i],c_data_temp[k+1][j][i]);
						fprintf(fid1,"%e %e %e %e %e %e \n",
								c_blend[k][j][i-1],c_blend[k][j][i+1],
								c_blend[k][j-1][i],c_blend[k][j+1][i],
								c_blend[k-1][j][i],c_blend[k+1][j][i]);
					}
					fprintf(fid1,"\n");
					fclose(fid1);
				}
			}
		}
	}

	return;
}




/******************************************************************************/
/*
 This function sets the blending factor that is used to blend QUICK and first
 order upwind
 */
/******************************************************************************/
void Conc_set_blendingfactor(Concentration *c, MAC_grid *grid,
		Parameters *params, double blend_fac) {

	int i, j, k;
	double blend_level, tol;

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

	double ***c_blend  = c -> ng_temp1;
	double ***c_lblend = c -> blend;

	if (Je == NY)
		Je = Je-1;
	if (Js == 0)
		Js = 1;
#ifndef XPERIODIC
	if (Ie == NX)
		Ie=Ie-1;
	if (Is == 0 )
		Is=1;
#endif

#ifndef ZPERIODIC
	if (Ke == NZ)
		Ke=Ke-1;
	if (Ks == 0 )
		Ks=1;
#endif

	for (k=grid->G_Ks; k<grid->G_Ke; k++) {
		for (j=grid->G_Js; j<grid->G_Je; j++) {
			for (i=grid->G_Is; i<grid->G_Ie; i++) {
				c_blend[k][j][i] = c_lblend[k][j][i];
			}
		}
	}

	tol = 1e-8;
	blend_level = blend_fac;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				if (c_lblend[k][j][i] < blend_level) {
					if ( c_lblend[k-1][j][i]  || c_lblend[k+1][j][i]   ||
					     c_lblend[k][j-1][i]  || c_lblend[k][j+1][i]   ||
					     c_lblend[k][j][i-1]  || c_lblend[k][j][i+1]   ) {

						c_blend[k][j][i] = blend_fac;
					}
				}
			}
		}
	}

	if (Js == 1) {
		j = 0;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				if (c_lblend[k][j][i] < blend_level) {
					if ( c_lblend[k-1][j][i]  || c_lblend[k+1][j][i]  ||
					     c_lblend[k][j+1][i]  ||
					     c_lblend[k][j][i-1]  || c_lblend[k][j][i+1]  ) {

						c_blend[k][j][i] = blend_fac;
					}
				}
			}
		}
	}

	if (Ks == 1) {
		k = 0;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				if (c_lblend[k][j][i] < blend_level) {
					if ( c_lblend[k+1][j][i]  ||
					     c_lblend[k][j-1][i]  || c_lblend[k][j+1][i] > tol  ||
					     c_lblend[k][j][i-1]  || c_lblend[k][j][i+1] > tol  ) {

						c_blend[k][j][i] = blend_fac;
					}
				}
			}
		}
	}

	if (Is == 1) {
		i = 0;
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				if (c_lblend[k][j][i] < blend_level) {
					if ( c_lblend[k-1][j][i]   || c_lblend[k+1][j][i]  ||
					     c_lblend[k][j-1][i]   || c_lblend[k][j+1][i]  ||
					     c_lblend[k][j][i+1]  ) {

						c_blend[k][j][i] = blend_fac;
					}
				}
			}
		}
	}

	if ((Is==1) && (Js == 1)) {
		j = 0;
		i = 0;
		for (k=Ks; k<Ke; k++) {
			if (c_lblend[k][j][i] < blend_level) {
					if ( c_lblend[k-1][j][i]  || c_lblend[k+1][j][i]  ||
					     c_lblend[k][j+1][i]  ||
					     c_lblend[k][j][i+1]   ) {

						c_blend[k][j][i] = blend_fac;
					}
			}
		}
	}

	if ((Is==1) && (Ks == 1)) {
		i = 0;
		k = 0;
		for (j=Js; j<Je; j++) {
			if (c_lblend[k][j][i] < blend_level) {
					if ( c_lblend[k+1][j][i]  ||
					     c_lblend[k][j-1][i]  || c_lblend[k][j+1][i]  ||
					     c_lblend[k][j][i+1]  ) {

						c_blend[k][j][i] = blend_fac;
					}
			}
		}
	}

	if ( (Js==1) && (Ks == 1)) {
		j = 0;
		k = 0;
			for (i=Is; i<Ie; i++) {
			if (c_lblend[k][j][i] < blend_level) {
					if ( c_lblend[k+1][j][i]  ||
					     c_lblend[k][j+1][i]  ||
					     c_lblend[k][j][i-1]  || c_lblend[k][j][i+1]  ) {

						c_blend[k][j][i] = blend_fac;
					}
			}
		}
	}

	if ( (Is==1) && (Js==1) && (Ks==1) ){
		i = 0;
		j = 0;
		k = 0;
		if (c_lblend[k][j][i] < blend_level) {
			if ( c_lblend[k+1][j][i]  ||
			     c_lblend[k][j+1][i]  ||
			     c_lblend[k][j][i+1]  ) {

				c_blend[k][j][i] = blend_fac;
			}
		}
	}

	for (k=grid->G_Ks; k<grid->G_Ke; k++) {
		for (j=grid->G_Js; j<grid->G_Je; j++) {
			for (i=grid->G_Is; i<grid->G_Ie; i++) {
				c_lblend[k][j][i] = c_blend[k][j][i];
			}
		}
	}

	return;
}




/******************************************************************************/
/*
 This function sets the blending factor that is used to blend QUICK and first
 order upwind
 */
/******************************************************************************/
void Conc_average_outofbounds(Concentration *c, MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ;
	int i, j, k;
	double ***c_data, ***c_ldata;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double fac;

	// Same for all quantities
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	c_data  = c->ng_temp1;
	c_ldata = c->data_temp;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	if (Je == NY)
		Je = Je-1;
	if (Js == 0)
		Js = 1;
#ifndef XPERIODIC
	if (Ie == NX)
		Ie=Ie-1;
	if (Is == 0 )
		Is=1;
#endif

#ifndef ZPERIODIC
	if (Ke == NZ)
		Ke=Ke-1;
	if (Ks == 0 )
		Ks=1;
#endif

	for (k=grid->G_Ks; k<grid->G_Ke; k++) {
		for (j=grid->G_Js; j<grid->G_Je; j++) {
			for (i=grid->G_Is; i<grid->G_Ie; i++) {
				c_data[k][j][i] = c_ldata[k][j][i];
			}
		}
	}


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
					c_data[k][j][i] = 1./6.*(c_ldata[k-1][j][i] + c_ldata[k+1][j][i] +
											 c_ldata[k][j-1][i] + c_ldata[k][j+1][i] +
											 c_ldata[k][j][i-1] + c_ldata[k][j][i+1] );
				}
			}
		}
	}

	if (Js == 1) {
		j = 0;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
					c_data[k][j][i] = 1./4.*(c_ldata[k-1][j][i] + c_ldata[k+1][j][i] +
											 c_ldata[k][j][i-1] + c_ldata[k][j][i+1] );
				}
			}
		}
	}

	if (Ks == 1) {
		k = 0;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
					c_data[k][j][i] = 1./4.*(c_ldata[k][j-1][i] + c_ldata[k][j+1][i] +
											 c_ldata[k][j][i-1] + c_ldata[k][j][i+1] );
				}
			}
		}
	}

	if (Is == 1) {
		i = 0;
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
					c_data[k][j][i] = 1./4.*(c_ldata[k-1][j][i] + c_ldata[k+1][j][i] +
											 c_ldata[k][j-1][i] + c_ldata[k][j+1][i] );

				}
			}
		}
	}

	if ((Is==1) && (Js == 1)) {
		j = 0;
		i = 0;
		for (k=Ks; k<Ke; k++) {
			if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
				c_data[k][j][i] = 1./2.*(c_ldata[k-1][j][i] + c_ldata[k+1][j][i] );
			}
		}
	}

	if ((Is==1) && (Ks == 1)) {
		i = 0;
		k = 0;
		for (j=Js; j<Je; j++) {
			if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
				c_data[k][j][i] = 1./2.*(c_ldata[k][j-1][i] + c_ldata[k][j+1][i] );
			}
		}
	}

	if ( (Js==1) && (Ks == 1)) {
		j = 0;
		k = 0;
		for (i=Is; i<Ie; i++) {
			if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
				c_data[k][j][i] = 1./2.*(c_ldata[k][j][i-1] + c_ldata[k][j][i+1] );
			}
		}
	}

	if ( (Is==1) && (Js==1) && (Ks==1) ){
		i = 0;
		j = 0;
		k = 0;
		if ((c_ldata[k][j][i] < 0) || (c_ldata[k][j][i]>1.0) ) {
			c_data[k][j][i] = 1./2.*(c_ldata[k][j][i+1] + c_ldata[k+1][j][i] );
		}
	}

	for (k=grid->G_Ks; k<grid->G_Ke; k++) {
		for (j=grid->G_Js; j<grid->G_Je; j++) {
			for (i=grid->G_Is; i<grid->G_Ie; i++) {
				c_ldata[k][j][i] = c_data[k][j][i];
			}
		}
	}

	return;
}






/******************************************************************************/
/*
 This function copies the pressure vector from "smaller domain" to the actual domain
 */
/******************************************************************************/
void Conc_copy_fromtemp(Concentration *c, MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double ***data_temp, ***data;


	// Now, get the pressure data on each processor
	data_temp = c->data_temp;
	data = c->data;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

//	if (Ie == params->NXM+1) Ie = Ie-1;
//	if (Je == params->NYM+1) Je = Je-1;
//	if (Ke == params->NZM+1) Ke = Ke-1;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				data[k][j][i] = data_temp[k][j][i];
//				data[k][j][i] = fmin(data[k][j][i], 0.0);
//				data[k][j][i] = fmax(data[k][j][i], 1.0);
			}
		}
	}

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Conc_set_quick_coefficients(Concentration *c, MAC_grid *grid, Parameters *params) {
// this is needed for non uniform grid
	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	double *aeW, *aeE, *aeEE;
	double *awWW, *awW, *awE;
	double *anS, *anN, *anNN;
	double *asSS, *asS, *asN;
	double *afB, *afF, *afFF;
	double *abBB, *abB, *abF;
	double *xc, *yc, *zc;
	double *xu, *yv, *zw;
	double xWW, xW, xE, xEE, xface;
	double ySS, yS, yN, yNN, yface;
	double zBB, zB, zF, zFF, zface;

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

	// indices start and end on current processor
	i_start = Is;
	j_start = Js;
	k_start = Ks;

	// exclude the half cell added
	i_end   = min(NX-1, Ie);
	j_end   = min(NY-1, Je);
	k_end   = min(NZ-1, Ke);

	xc = grid->xc;
	yc = grid->yc;
	zc = grid->zc;
	xu = grid->xu;
	yv = grid->yv;
	zw = grid->zw;

	aeW  = Memory_allocate_1D_array(GVG_DOUBLE, NX);
	aeE  = Memory_allocate_1D_array(GVG_DOUBLE, NX);
	aeEE = Memory_allocate_1D_array(GVG_DOUBLE, NX);
	awWW = Memory_allocate_1D_array(GVG_DOUBLE, NX);
	awW  = Memory_allocate_1D_array(GVG_DOUBLE, NX);
	awE  = Memory_allocate_1D_array(GVG_DOUBLE, NX);
	anS  = Memory_allocate_1D_array(GVG_DOUBLE, NY);
	anN  = Memory_allocate_1D_array(GVG_DOUBLE, NY);
	anNN = Memory_allocate_1D_array(GVG_DOUBLE, NY);
	asSS = Memory_allocate_1D_array(GVG_DOUBLE, NY);
	asS  = Memory_allocate_1D_array(GVG_DOUBLE, NY);
	asN  = Memory_allocate_1D_array(GVG_DOUBLE, NY);
	afB  = Memory_allocate_1D_array(GVG_DOUBLE, NZ);
	afF  = Memory_allocate_1D_array(GVG_DOUBLE, NZ);
	afFF = Memory_allocate_1D_array(GVG_DOUBLE, NZ);
	abBB = Memory_allocate_1D_array(GVG_DOUBLE, NZ);
	abB  = Memory_allocate_1D_array(GVG_DOUBLE, NZ);
	abF  = Memory_allocate_1D_array(GVG_DOUBLE, NZ);

	for (i=1;i<NX-1;i++){
		xWW = xc[i-1];
		xW  = xc[i];
		xE  = xc[i+1];
		xface  = xu[i+1];
		awWW[i] = ( (xface-xW)*(xface-xE))  / ( (xWW-xW)*(xWW-xE));
		awW[i]  = ( (xface-xWW)*(xface-xE)) / ( (xW-xWW)*(xW-xE));
		awE[i]  = ( (xface-xWW)*(xface-xW)) / ( (xE-xWW)*(xE-xW));
	}

	// coefficient at i=0 & i=NX-1 are needed when using xperiodic
	awWW[0] = awWW[1];
	awW[0]  = awW[1];
	awE[0]  = awE[1];
	awWW[NX-1] = awWW[NX-2];
	awW[NX-1]  = awW[NX-2];
	awE[NX-1]  = awE[NX-2];

	for (i=0;i<NX-2;i++){
		xW  = xc[i];
		xE  = xc[i+1];
		xEE = xc[i+2];
		xface  = xu[i+1];
		aeW[i]  = ( (xface-xEE)*(xface-xE)) / ( (xW-xEE)*(xW-xE));
		aeE[i]  = ( (xface-xEE)*(xface-xW)) / ( (xE-xEE)*(xE-xW));
		aeEE[i] = ( (xface-xE)*(xface-xW))  / ( (xEE-xE)*(xEE-xW));
	}

	// coefficient at i=NX-2 & i=NX-1 are needed when using xperiodic
	aeW[NX-2]  = aeW[NX-3];
	aeE[NX-2]  = aeE[NX-3];
	aeEE[NX-2] = aeEE[NX-3];
	aeW[NX-1]  = aeW[NX-3];
	aeE[NX-1]  = aeE[NX-3];
	aeEE[NX-1] = aeEE[NX-3];

	for (j=1;j<NY-1;j++){
		ySS = yc[j-1];
		yS  = yc[j];
		yN  = yc[j+1];
		yface  = yv[j+1];
		asSS[j] = (  (yface-yS)*(yface-yN)) / ( (ySS-yS)*(ySS-yN));
		asS[j]  = ( (yface-ySS)*(yface-yN)) / ( (yS-ySS)*(yS-yN));
		asN[j]  = ( (yface-ySS)*(yface-yS)) / ( (yN-ySS)*(yN-yS));
	}

	// Following values should never be used, but they are initialized anyways;
	asSS[0] = 0.0;
	asS[0]  = 1.0;
	asN[0]  = 0.0;
	asSS[NY-1] = asSS[NY-2];
	asS[NY-1]  = asS[NY-2];
	asN[NY-1]  = asN[NY-2];

	for (j=0;j<NY-2;j++){
		yS  = yc[j];
		yN  = yc[j+1];
		yNN = yc[j+2];
		yface  = yv[j+1];
		anS[j]  = (  (yface-yN)*(yface-yNN)) / ( (yS-yN)*(yS-yNN));
		anN[j]  = (  (yface-yS)*(yface-yNN)) / ( (yN-yS)*(yN-yNN));
		anNN[j] = (  (yface-yS)*(yface-yN))  / ( (yNN-yS)*(yNN-yN));
	}

	// Following values should never be used, but they are initialized anyways;
	anS[NY-2]  = anS[NY-3];
	anN[NY-2]  = anN[NY-3];
	anNN[NY-2] = anNN[NY-3];
	anS[NY-1]  = anS[NY-3];
	anN[NY-1]  = anN[NY-3];
	anNN[NY-1] = anNN[NY-3];

	for (k=1;k<NZ-1;k++){
		zBB = zc[k-1];
		zB  = zc[k];
		zF  = zc[k+1];
		zface  = zw[k+1];
		abBB[k] = (  (zface-zB)*(zface-zF)) / ( (zBB-zB)*(zBB-zF));
		abB[k]  = ( (zface-zBB)*(zface-zF)) / ( (zB-zBB)*(zB-zF));
		abF[k]  = ( (zface-zBB)*(zface-zB)) / ( (zF-zBB)*(zF-zB));
	}

	// coefficient at k=0 & k=NZ-1 are needed when using zperiodic
	abBB[0] = abBB[1];
	abB[0]  = abB[1];
	abF[0]  = abF[1];
	abBB[NZ-1] = abBB[NZ-2];
	abB[NZ-1]  = abB[NZ-2];
	abF[NZ-1]  = abF[NZ-2];

	for (k=0;k<NZ-2;k++){
		zB  = zc[k];
		zF  = zc[k+1];
		zFF = zc[k+2];
		zface  = zw[k+1];
		afB[k]  = (  (zface-zF)*(zface-zFF)) / ( (zB-zF)*(zB-zFF));
		afF[k]  = (  (zface-zB)*(zface-zFF)) / ( (zF-zB)*(zF-zFF));
		afFF[k] = (  (zface-zB)*(zface-zF))  / ( (zFF-zB)*(zFF-zF));
	}

	// coefficient at k=NZ-2 & k=NZ-1 are needed when using zperiodic
	afB[NZ-2]  = afB[NZ-3];
	afF[NZ-2]  = afF[NZ-3];
	afFF[NZ-2] = afFF[NZ-3];
	afB[NZ-1]  = afB[NZ-3];
	afF[NZ-1]  = afF[NZ-3];
/*
	sprintf(bin_filename, "quickcoeff%d.dat",params->rank);
	fid = fopen(bin_filename, "w");
	for (i=0;i<NX;i++) {
		fprintf(fid,"%f  %f  %f %f  %f  %f\n",aeW[i],aeE[i],aeEE[i],awWW[i],awW[i],awE[i]);
	}
	fprintf(fid,"\n\n\n");
	for (j=0;j<NY;j++) {
		fprintf(fid,"%f  %f  %f %f  %f  %f\n",anS[j],anN[j],anNN[j],asSS[j],asS[j],asN[j]);
	}
	fprintf(fid,"\n\n\n");
	for (k=0;k<NZ;k++) {
		fprintf(fid,"%f  %f  %f %f  %f  %f\n",afB[k],afF[k],afFF[k],abBB[k],abB[k],abF[k]);
	}
	fprintf(fid,"\n\n\n");
*/

	c->aeW  = aeW;
	c->aeE  = aeE;
	c->aeEE = aeEE;
	c->awWW = awWW;
	c->awW  = awW;
	c->awE  = awE;
	c->anS  = anS ;
	c->anN  = anN;
	c->anNN = anNN;
	c->asSS = asSS;
	c->asS  = asS;
	c->asN  = asN;
	c->afB  = afB;
	c->afF  = afF;
	c->afFF = afFF;
	c->abBB = abBB;
	c->abB  = abB;
	c->abF  = abF;

	return;
}


#include "lsolver/csolve_direct.c"



/******************************************************************************/
/*
 Solve the system Ax = b		using the conjugate-gradient method
with diagonal preconditioner see Saad 2003 iterative methods
!!! For the iterative method we use the memory of conc[0] !!!!
 */
/******************************************************************************/


int Conc_solve_cg(int iconc, Cart3d_bag *data_bag) {

	int iters;
	int i, j, k;
	char statement[100];

	double EPS = 1e-16; // Small number to prevent division by zero

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *conc = data_bag -> c[iconc];
	double ***data = data_bag -> c[iconc]->data;
	double ***rhs = data_bag -> c[iconc]-> ng_rhs;
	double ***visc = data_bag -> c[iconc]->ng_viscous;


	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;
    int nghost= params->ghost_nodes;

	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = min( grid->G_Ie, NX-1);  // No ghost nodes
	int Je = min( grid->G_Je, NY-1);
	int Ke = min( grid->G_Ke, NZ-1);



	//--------------------------------------------------------------------------
	// Conjugate-gradient variables
	//--------------------------------------------------------------------------

	double ***r = data_bag -> c[0]->ng_r;  ////c -> ng_r;  // Residual
	double ***d = data_bag -> c[0]->d;  //  c -> d;
	double ***Ad = data_bag -> c[0]->ng_Ad;      //c -> ng_Ad;  // A * d


#ifdef VOF_SCALAR
	double ***vfu = data_bag->lag->vfu;
	double ***vfv = data_bag->lag->vfv;
	double ***vfw = data_bag->lag->vfw;
	double ***ng_vfc = data_bag->lag->ng_vfc;
#endif

	// Inner products
	double rr, rr_old, rr0, dAd, RMS, in_tot;

	// Global  coefficients
	double alpha, beta;

	in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-1));


	double iddx = grid -> idx_c[1];
	double iddy = grid -> idy_c[1];
	double iddz = grid -> idz_c[1];

	iddx = iddx * iddx;
	iddy = iddy * iddy;
	iddz = iddz * iddz;


#ifdef VOF_SCALAR
	Conc_laplacian_vof(Ad, data, vfu, vfv, vfw, ng_vfc,grid, params, iconc);
#else
	Conc_laplacian(Ad, data, grid, params, iconc);
#endif

	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				r[k][j][i] = rhs[k][j][i] - Ad[k][j][i];
				d[k][j][i] = r[k][j][i];


			}
		}
	}


	/*------------------------------------------------------------------------*/
	/*
	 Conjugate-gradient solve
	 */
	/*------------------------------------------------------------------------*/

	rr = Conc_innerProd(r, r, grid, params);
	rr0 = rr + EPS;
	iters  = 0;



	while (iters < params -> CG_MAXIT) {

		iters++;


		// Update boundary cells for 'd'

		Communication_update_ghost_nodes_flow_variable(d, 'c',params->ghost_nodes , data_bag);
		Conc_set_boundary_values(d , iconc ,CENTRAL_PERTURBATION ,grid , params);


		// Ad = A * d
#ifdef VOF_SCALAR
		Conc_laplacian_vof(Ad, d , vfu, vfv, vfw, ng_vfc, grid, params, iconc);
#else
		Conc_laplacian(Ad, d ,grid, params, iconc);
#endif

		dAd = Conc_innerProd(d, Ad, grid, params);

		alpha = rr / (dAd + EPS);

		// x = x + alpha * d;
		// r = r - alpha * Ad;
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
					data[k][j][i] += alpha * d[k][j][i];
					r[k][j][i] -= alpha * Ad[k][j][i];
				}
			}
		}

		rr_old = rr;
		rr = Conc_innerProd(r, r, grid, params);
 		RMS = sqrt(in_tot * rr);
		if (RMS < params -> CG_ETOL)
			break;

		beta = rr / (rr_old + EPS);

		// d = beta * d + r
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

						d[k][j][i] = beta * d[k][j][i] + r[k][j][i];

				}
			}
		}
	}

	if(params->rank == 0) printf("\n Concentration converged to %g after %d iterations \n ", RMS, iters);
	//Display_progress(params, statement);




	return iters;
}


/******************************************************************************/
void Conc_set_conv_viscous_quick_mixed(int iconc, Cart3d_bag *data_bag) {
// has to be tested carefully is no compatible with the ghost cell method
	int i, j, k;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

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

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;

	double *aeW  = c -> aeW;
	double *aeE  = c -> aeE;
	double *aeEE = c -> aeEE;
	double *awWW = c -> awWW;
	double *awW  = c -> awW;
	double *awE  = c -> awE;
	double *anS  = c -> anS;
	double *anN  = c -> anN;
	double *anNN = c -> anNN;
	double *asSS = c -> asSS;
	double *asS  = c -> asS;
	double *asN  = c -> asN;
	double *afB  = c -> afB;
	double *afF  = c -> afF;
	double *afFF = c -> afFF;
	double *abBB = c -> abBB;
	double *abB  = c -> abB;
	double *abF  = c -> abF;

	// Get the local velocities at the location where they are defined
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***c_data = c -> data;

#ifdef VAR_VISC
	#ifdef LES
	double ***nut = data_bag -> smag -> cdev[iconc] -> mSct;
	#elif defined RANS
	double ***nut = data_bag -> rans -> nutc;
	#endif
#endif

	double ***conv = c -> ng_conv;
	double ***visc = c -> ng_viscous;
	double ***conv_old   = c -> ng_conv_old;
	double ***conv_quick = c -> ng_conv_quick;
	double ***visc_explicit = c -> ng_visc_explicit;

	// Constant settling speed of particle. Zero for temperature and Salinity
	double V_s0 = c -> v_settl0;

	if ( (params->which_stage == 1) || (params->which_stage == 2) ) {
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					conv_old[k][j][i] = conv[k][j][i] - visc_explicit[k][j][i];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------
	double iPe = 1.0 / c->Pe;
	double nuE = iPe;
	double nuW = iPe;
	double nuN = iPe;
	double nuS = iPe;
	double nuF = iPe;
	double nuB = iPe;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				/*------------------------------------------------------------*/
				/*
				 Calculate viscosity
				 */
				/*------------------------------------------------------------*/
#ifdef VAR_VISC
				//--------------------------------------------------------------
				// Variable viscosity - eddy viscosity for LES and RANS
				//--------------------------------------------------------------
				nuE = 0.5 * ( nut[k][j][i+1] + nut[k][j][i] );
				nuW = 0.5 * ( nut[k][j][i-1] + nut[k][j][i] );

				nuN = 0.5 * ( nut[k][j+1][i] + nut[k][j][i] );
				nuS = 0.5 * ( nut[k][j-1][i] + nut[k][j][i] );

				nuF = 0.5 * ( nut[k+1][j][i] + nut[k][j][i] );
				nuB = 0.5 * ( nut[k-1][j][i] + nut[k][j][i] );

	#if defined LES || defined RANS
				//--------------------------------------------------------------
				// LES and RANS 'nu' don't account for Peclet number
				//--------------------------------------------------------------
				nuE += iPe;
				nuW += iPe;
				nuN += iPe;
				nuS += iPe;
				nuF += iPe;
				nuB += iPe;
	#endif // LES or RANS
#endif // VAR_VISC


				/*------------------------------------------------------------*/
				/*
				 Calculate viscous terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// dc/dx East/West
				//--------------------------------------------------------------
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * idx_c[i];


				dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i-1];


				//--------------------------------------------------------------
				// dc/dy North/South
				//--------------------------------------------------------------
				dcdyN = ( c_data[k][j+1][i] - c_data[k][j][i] ) * idy_c[j];


				dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] ) * idy_c[j-1];


				//--------------------------------------------------------------
				// dc/dz Front/Back
				//--------------------------------------------------------------
				dcdzF = ( c_data[k+1][j][i]-c_data[k][j][i] ) * idz_c[k];

				dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) * idz_c[k-1];


				//--------------------------------------------------------------
				// d2c/dx2, d2c/dy2 and d2c/dz2
				//--------------------------------------------------------------
				d2cdx2 = ( nuE * dcdxE - nuW * dcdxW ) * idx_u[i];
				d2cdy2 = ( nuN * dcdyN - nuS * dcdyS ) * idy_v[j];
				d2cdz2 = ( nuF * dcdzF - nuB * dcdzB ) * idz_w[k];


				/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// Convective velocities on the faces
				//--------------------------------------------------------------
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];

				if (grid->c_status[k][j][i] != SOLID) {

						if (j != NY-2) vN = vN + V_s0;
						vS = vS + V_s0;

				}

				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];



				//--------------------------------------------------------------
				// Concentration value on the faces
				// For explanation of various "if conditions", check cN and cS
				// calculation, which is representative of cE and cW, cF and cB
				//--------------------------------------------------------------
				// Calculate  cE

				//// There is an error included here!!


				if (uE > 0) {
					if (i==0) {
#ifdef XPERIODIC_CONC
						cE = awWW[i]*c_data[k][j][i-1] + awW[i]*c_data[k][j][i] + awE[i]*c_data[k][j][i+1];
#else
						cE = c_data[k][j][i];
#endif
					}
					else {
						cE = awWW[i]*c_data[k][j][i-1] + awW[i]*c_data[k][j][i] + awE[i]*c_data[k][j][i+1];
					}
				}
				else {
					if (i==NX-2) {
#ifdef XPERIODIC_CONC
						cE = aeW[i]*c_data[k][j][i] + aeE[i]*c_data[k][j][i+1] + aeEE[i]*c_data[k][j][i+2];
#else
						cE = c_data[k][j][i+1];
#endif
					}
					else if ( (grid->c_status[k][j][i+2] == FLUID) || (grid->c_status[k][j][i+2] == IMMERSED) ) {
						cE = aeW[i]*c_data[k][j][i] + aeE[i]*c_data[k][j][i+1] + aeEE[i]*c_data[k][j][i+2];
					}
					else {
						cE = c_data[k][j][i+1];
					}
				}

				// Calculate  cW
				if (i==0) {
#ifdef XPERIODIC_CONC
					if (uW > 0) {
						cW = awWW[0]*c_data[k][j][i-2] + awW[0]*c_data[k][j][i-1] + awE[0]*c_data[k][j][i];
					}
					else {
						cW = aeW[0]*c_data[k][j][i-1] + aeE[0]*c_data[k][j][i] + aeEE[0]*c_data[k][j][i+1];
					}
#else // not XPERIODIC

					cW = c_data[k][j][-1];

#endif // not XPERIODIC
				}

				else {
					if (uW > 0) {
						if (i==1) {
#ifdef XPERIODIC_CONC
							cW = awWW[i-1]*c_data[k][j][i-2] + awW[i-1]*c_data[k][j][i-1] + awE[i-1]*c_data[k][j][i];
#else
							cW = c_data[k][j][i-1];
#endif
						}
						else if ( (grid->c_status[k][j][i-2] == FLUID) || (grid->c_status[k][j][i-2] == IMMERSED) ) {
							cW = awWW[i-1]*c_data[k][j][i-2] + awW[i-1]*c_data[k][j][i-1] + awE[i-1]*c_data[k][j][i];
						}
						else {
							cW = c_data[k][j][i-1];
						}
					}
					else {
						cW = aeW[i-1]*c_data[k][j][i-1] + aeE[i-1]*c_data[k][j][i] + aeEE[i-1]*c_data[k][j][i+1];
					}
				}

				if (vN > 0) {
					if (j==0) {
						// Near the boundary use first order upwind
						cN = c_data[k][j][i];
					}
					else {
						cN = asSS[j]*c_data[k][j-1][i] + asS[j]*c_data[k][j][i] + asN[j]*c_data[k][j+1][i];
					}
				}
				else {
					if (j == NY-2) {
						// Near the boundary use first order upwind
						cN = c_data[k][j+1][i];
					}
					else if ( (grid->c_status[k][j+2][i] == FLUID) || (grid->c_status[k][j+2][i] == IMMERSED) ) {
						// Use QUICK when the stencil  permits (i.e all the
						// points have physical value)
						cN = anS[j]*c_data[k][j][i] + anN[j]*c_data[k][j+1][i] + anNN[j]*c_data[k][j+2][i];
					}
					else {
						// If QUICK stencil has solid nodes, use first order
						// upwind
						cN = c_data[k][j+1][i];
					}
				}

				if (vS > 0) {
					if (j==0) {
						// Impose Neumann b/c i.e dc/dy = 0.0
						cS = c_data[k][j][i];
					}
					else if (j==1) {
						cS = c_data[k][j-1][i];
					}
					else  if ( (grid->c_status[k][j-2][i] == FLUID) || (grid->c_status[k][j-2][i] == IMMERSED) ) {
						// Use QUICK when the stencil  permits (i.e all the
						// points have physical value)
						cS = asSS[j-1]*c_data[k][j-2][i] + asS[j-1]*c_data[k][j-1][i] + asN[j-1]*c_data[k][j][i];
					}
					else {
						// If QUICK stencil has solid nodes, use first order
						// upwind
						cS = c_data[k][j-1][i];
					}
				}
				else {
					if (j==0) {
						// Impose Neumann b/c i.e dc/dy = 0.0
						cS = c_data[k][j][i];
					}
					else {
						cS = anS[j-1]*c_data[k][j-1][i] + anN[j-1]*c_data[k][j][i] + anNN[j-1]*c_data[k][j+1][i];
					}
				}

				// Calculate  cF
				if (wF > 0) {
					if (k==0) {
#ifdef ZPERIODIC_CONC
						cF = abBB[k]*c_data[k-1][j][i] + abB[k]*c_data[k][j][i] + abF[k]*c_data[k+1][j][i];
#else
						cF = c_data[k][j][i];
#endif
					}
					else {
						cF = abBB[k]*c_data[k-1][j][i] + abB[k]*c_data[k][j][i] + abF[k]*c_data[k+1][j][i];
					}
				}
				else {
					if (k==NZ-2) {
#ifdef ZPERIODIC_CONC
						cF = afB[k]*c_data[k][j][i] + afF[k]*c_data[k+1][j][i] + afFF[k]*c_data[k+2][j][i];
#else
						cF = c_data[k+1][j][i];
#endif
					}
					else if ( (grid->c_status[k+2][j][i] == FLUID) || (grid->c_status[k+2][j][i] == IMMERSED) ) {
						cF = afB[k]*c_data[k][j][i] + afF[k]*c_data[k+1][j][i] + afFF[k]*c_data[k+2][j][i];
					}
					else {
						cF = c_data[k+1][j][i];
					}
				}

				// Calculate  cB
				if (k==0) {
#ifdef ZPERIODIC_CONC
					if (wB > 0) {
						cB = abBB[0]*c_data[k-2][j][i] + abB[0]*c_data[k-1][j][i] + abF[0]*c_data[k][j][i];
					}
					else {
						cB = afB[0]*c_data[k-1][j][i] + afF[0]*c_data[k][j][i] + afFF[0]*c_data[k+1][j][i];
					}
#else
					// Otherwise impose Neumann b/c i.e dc/dx = 0.0
					cB = c_data[k][j][i];
#endif
				}

				else {

					if (wB > 0) {
						if (k==1) {
#ifdef ZPERIODIC_CONC
							cB = abBB[k-1]*c_data[k-2][j][i] + abB[k-1]*c_data[k-1][j][i] + abF[k-1]*c_data[k][j][i];
#else
							cB = c_data[k-1][j][i];
#endif
						}
						else if ( (grid->c_status[k-2][j][i] == FLUID) || (grid->c_status[k-2][j][i] == IMMERSED) ) {
							cB = abBB[k-1]*c_data[k-2][j][i] + abB[k-1]*c_data[k-1][j][i] + abF[k-1]*c_data[k][j][i];
						}
						else {
							cB = c_data[k-1][j][i];
						}
					}
					else {
						cB = afB[k-1]*c_data[k-1][j][i] + afF[k-1]*c_data[k][j][i] + afFF[k-1]*c_data[k+1][j][i];
					}
				}










				ucE = uE*cE;
				ucW = uW*cW;

				vcN = vN*cN;
				vcS = vS*cS;

				wcF = wF*cF;
				wcB = wB*cB;

				//--------------------------------------------------------------
				// d/dx(uc), d/dy(vc), and d/dz(wc)
				//--------------------------------------------------------------
				dcudx = ( ucE - ucW ) * idx_u[i];
				dcvdy = ( vcN - vcS ) * idy_v[j];
				dcwdz = ( wcF - wcB ) * idz_w[k];


				/*------------------------------------------------------------*/
				/*
				 Store convective and viscous terms
				 */
				/*------------------------------------------------------------*/
				conv[k][j][i] = dcudx + dcvdy + dcwdz;

#ifdef CONC_FULLY_EXPLICIT
				visc[k][j][i] = 0.0;
				visc_explicit[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
#elif defined CONC_FULLY_IMPLICIT
				visc[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
				visc_explicit[k][j][i] = 0.0;
#elif defined CONC_SEMI_IMPLICIT
				visc[k][j][i] = d2cdy2;
				visc_explicit[k][j][i] = d2cdx2 + d2cdz2;
#else
				#error choose a diffusion scheme
#endif

#ifdef BQUICK
				conv_quick[k][j][i] = conv[k][j][i];
#endif

			}
		}
	}

	return;
}
void Conc_set_conv_viscous_central_mixed(int iconc, Cart3d_bag *data_bag) {
/// The boundary condition are inforced by the boundary nodes!!

	int i, j, k;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	// exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie); // points on the first ghost sell
	int j_end = min(NY-1, grid->G_Je);
	int k_end = min(NZ-1, grid->G_Ke);

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;


	// Get the local velocities at the location where they are defined
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***c_data = c -> data;

#ifdef VAR_VISC
	#ifdef LES
	double ***nut = data_bag -> smag -> cdev[iconc] -> mSct;
	#elif defined RANS
	double ***nut = data_bag -> rans -> nutc;
	#endif
#endif

	double ***conv = c -> ng_conv;
	double ***visc = c -> ng_viscous;
	double ***conv_old = c -> ng_conv_old;
	double ***visc_explicit = c -> ng_visc_explicit;

	// Constant settling speed of particle. Zero for temperature and Salinity
	double V_s0 = c -> v_settl0;


#ifdef GRID_UNIFORM
	double ih= grid -> idx_u[i_start+1];
#else
	double *wc2uW = grid -> wc2uW;
	double *wc2uE = grid -> wc2uE;
	double *wc2vN = grid -> wc2vN;
	double *wc2vS = grid -> wc2vS;
	double *wc2wF = grid -> wc2wF;
	double *wc2wB = grid -> wc2wB;
#endif




	if ( (params->which_stage == 1) || (params->which_stage == 2) ) {
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					conv_old[k][j][i] = conv[k][j][i] - visc_explicit[k][j][i];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------
	double iPe = 1.0 / c->Pe;
	double nuE = iPe;
	double nuW = iPe;
	double nuN = iPe;
	double nuS = iPe;
	double nuF = iPe;
	double nuB = iPe;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				/*------------------------------------------------------------*/
				/*
				 Calculate viscosity
				 */
				/*------------------------------------------------------------*/
#ifdef VAR_VISC
				//--------------------------------------------------------------
				// Variable viscosity - eddy viscosity for LES and RANS
				//--------------------------------------------------------------
				nuE = 0.5 * ( nut[k][j][i+1] + nut[k][j][i] );
				nuW = 0.5 * ( nut[k][j][i-1] + nut[k][j][i] );

				nuN = 0.5 * ( nut[k][j+1][i] + nut[k][j][i] );
				nuS = 0.5 * ( nut[k][j-1][i] + nut[k][j][i] );

				nuF = 0.5 * ( nut[k+1][j][i] + nut[k][j][i] );
				nuB = 0.5 * ( nut[k-1][j][i] + nut[k][j][i] );

	#if defined LES || defined RANS
				//--------------------------------------------------------------
				// LES and RANS 'nu' don't account for Peclet number
				//--------------------------------------------------------------
				nuE += iPe;
				nuW += iPe;
				nuN += iPe;
				nuS += iPe;
				nuF += iPe;
				nuB += iPe;
	#endif // LES or RANS
#endif // VAR_VISC


				/*------------------------------------------------------------*/
				/*
				 Calculate viscous terms
				 */
				/*------------------------------------------------------------*/
#ifdef GRID_UNIFORM
				//--------------------------------------------------------------
				// dc/dx East/West
				//--------------------------------------------------------------
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * ih;


				dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * ih;


				//--------------------------------------------------------------
				// dc/dy North/South
				//--------------------------------------------------------------
				dcdyN = ( c_data[k][j+1][i] - c_data[k][j][i] ) * ih;


				dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] ) * ih;


				//--------------------------------------------------------------
				// dc/dz Front/Back
				//--------------------------------------------------------------
				dcdzF = ( c_data[k+1][j][i]-c_data[k][j][i] ) * ih;


				dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) * ih;



				//--------------------------------------------------------------
				// d2c/dx2, d2c/dy2 and d2c/dz2
				//--------------------------------------------------------------
				d2cdx2 = ( nuE * dcdxE - nuW * dcdxW ) * ih;
				d2cdy2 = ( nuN * dcdyN - nuS * dcdyS ) * ih;
				d2cdz2 = ( nuF * dcdzF - nuB * dcdzB ) * ih;

#else

				//--------------------------------------------------------------
				// dc/dx East/West
				//--------------------------------------------------------------
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * idx_c[i];


				dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i-1];


				//--------------------------------------------------------------
				// dc/dy North/South
				//--------------------------------------------------------------
				dcdyN = ( c_data[k][j+1][i] - c_data[k][j][i] ) * idy_c[j];


				dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] ) * idy_c[j-1];


				//--------------------------------------------------------------
				// dc/dz Front/Back
				//--------------------------------------------------------------
				dcdzF = ( c_data[k+1][j][i]-c_data[k][j][i] ) * idz_c[k];

				dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) * idz_c[k-1];


				//--------------------------------------------------------------
				// d2c/dx2, d2c/dy2 and d2c/dz2
				//--------------------------------------------------------------
				d2cdx2 = ( nuE * dcdxE - nuW * dcdxW ) * idx_u[i];
				d2cdy2 = ( nuN * dcdyN - nuS * dcdyS ) * idy_v[j];
				d2cdz2 = ( nuF * dcdzF - nuB * dcdzB ) * idz_w[k];




#endif
				/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// Convective velocities on the faces
				//--------------------------------------------------------------
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];


				if (j != NY-2) vN = vN + V_s0;
				vS = vS + V_s0;


				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];

				//--------------------------------------------------------------
				// Concentration value on the faces
				//--------------------------------------------------------------
//#ifdef GRID_UNIFORM  // I think the following formulation yields also a second order convergence for a non-uniform grid??
				cE = 0.5 * ( c_data[k][j][i] + c_data[k][j][i+1] );

				cW = 0.5 * ( c_data[k][j][i] + c_data[k][j][i-1] );



				cN = 0.5 * ( c_data[k][j][i] + c_data[k][j+1][i] );

				cS = 0.5 * ( c_data[k][j][i] + c_data[k][j-1][i] );



				cF = 0.5 * ( c_data[k][j][i] + c_data[k+1][j][i] );

				cB = 0.5 * ( c_data[k][j][i] + c_data[k-1][j][i] );
/*#else
				cE = ( c_data[k][j][i]*wc2uW[i+1]  + c_data[k][j][i+1]*wc2uE[i+1] );
				cW = ( c_data[k][j][i]*wc2uE[i] + c_data[k][j][i-1]*wc2uW[i] );

				cN = ( c_data[k][j][i]*wc2vS[j+1] + c_data[k][j+1][i]*wc2vN[j+1] );
				cS = ( c_data[k][j][i]*wc2vN[j] + c_data[k][j-1][i]*wc2vS[j] );

				cF = ( c_data[k][j][i]*wc2wB[k+1] + c_data[k+1][j][i]*wc2wF[k+1] );
				cB = ( c_data[k][j][i]*wc2wF[k] + c_data[k-1][j][i]*wc2wB[k] );

#endif*/

				ucE = uE*cE;
				ucW = uW*cW;

				vcN = vN*cN;
				vcS = vS*cS;

				wcF = wF*cF;
				wcB = wB*cB;

				//--------------------------------------------------------------
				// d/dx(uc), d/dy(vc), and d/dz(wc)
				//--------------------------------------------------------------
				dcudx = ( ucE - ucW ) * idx_u[i];
				dcvdy = ( vcN - vcS ) * idy_v[j];
				dcwdz = ( wcF - wcB ) * idz_w[k];



				/*------------------------------------------------------------*/
				/*
				 Store convective and viscous terms
				 */
				/*------------------------------------------------------------*/
				conv[k][j][i] = dcudx + dcvdy + dcwdz;
#ifdef CONC_FULLY_EXPLICIT
				visc[k][j][i] = 0.0;
				visc_explicit[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
#elif defined CONC_FULLY_IMPLICIT
				visc[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
				visc_explicit[k][j][i] = 0.0;
#elif defined CONC_SEMI_IMPLICIT
				visc[k][j][i] = d2cdy2;
				visc_explicit[k][j][i] = d2cdx2 + d2cdz2;
#else
				#error  choose a diffusion scheme
#endif
			}
		}
	}

	return;
}

#ifdef VOF_SCALAR


void Conc_set_conv_viscous_central_mixed_vof(int iconc, Cart3d_bag *data_bag) {
/// The boundary condition are enforced by the boundary nodes!!

	int i, j, k;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	// exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie); // points on the first ghost cell
	int j_end = min(NY-1, grid->G_Je);
	int k_end = min(NZ-1, grid->G_Ke);

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;


	// Get the local velocities at the location where they are defined

	double ***c_data = c -> data;

	double ***conv = c -> ng_conv;
	double ***visc = c -> ng_viscous;
	double ***conv_old = c -> ng_conv_old;
	double ***visc_explicit = c -> ng_visc_explicit;

	double ***u_vof= data_bag->u->data_vof;
	double ***v_vof= data_bag->v->data_vof;
	double ***w_vof= data_bag->w->data_vof;

	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;


	double ***vfu = data_bag->lag->vfu;
	double ***vfv = data_bag->lag->vfv;
	double ***vfw = data_bag->lag->vfw;
	double ***ng_vfc = data_bag->lag->ng_vfc;

	double iPe = 1.0 / c->Pe;
	double conductivity_s = params->conductivity_s[iconc];
	double heat_cap_s = params->vol_heat_cap_s[iconc];


	// Constant settling speed of particle. Zero for temperature and Salinity
	double ih= grid -> idx_u[i_start];
	double ihsq= ih*ih;

#ifdef CONC_QUICK

	double *aeW  = c -> aeW;
	double *aeE  = c -> aeE;
	double *aeEE = c -> aeEE;
	double *awWW = c -> awWW;
	double *awW  = c -> awW;
	double *awE  = c -> awE;
	double *anS  = c -> anS;
	double *anN  = c -> anN;
	double *anNN = c -> anNN;
	double *asSS = c -> asSS;
	double *asS  = c -> asS;
	double *asN  = c -> asN;
	double *afB  = c -> afB;
	double *afF  = c -> afF;
	double *afFF = c -> afFF;
	double *abBB = c -> abBB;
	double *abB  = c -> abB;
	double *abF  = c -> abF;
#endif




	if ( (params->which_stage == 1) || (params->which_stage == 2) ) {
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					conv_old[k][j][i] = conv[k][j][i] - visc_explicit[k][j][i];
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------

	double lambdaN, lambdaS, lambdaE, lambdaW, lambdaB, lambdaF ;


	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				/*------------------------------------------------------------*/
				/*
				 Calculate diffusive terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// Variable conductivity with the arithmetic average
				//--------------------------------------------------------------
				lambdaE =  vfu[k][j][i+1]*(conductivity_s-1)+1;
				lambdaW =  vfu[k][j][i]*(conductivity_s-1)+1;

				lambdaN = vfv[k][j+1][i]*(conductivity_s-1)+1;
				lambdaS = vfv[k][j][i]*(conductivity_s-1)+1;

				lambdaF = vfw[k+1][j][i]*(conductivity_s-1)+1;
				lambdaB = vfw[k][j][i]*(conductivity_s-1)+1;


				//--------------------------------------------------------------
				// dc/dx East/West
				//--------------------------------------------------------------
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) ; // we divide by the cell size when everything is collected
				dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) ;


				//--------------------------------------------------------------
				// dc/dy North/South
				//--------------------------------------------------------------
				dcdyN = ( c_data[k][j+1][i] - c_data[k][j][i] ) ;
				dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] ) ;


				//--------------------------------------------------------------
				// dc/dz Front/Back
				//--------------------------------------------------------------
				dcdzF = ( c_data[k+1][j][i]-c_data[k][j][i] ) ;
				dcdzB = ( c_data[k][j][i]-c_data[k-1][j][i] ) ;



				//--------------------------------------------------------------
				// d2c/dx2, d2c/dy2 and d2c/dz2
				//--------------------------------------------------------------
				d2cdx2 = ( lambdaE * dcdxE - lambdaW * dcdxW ) ;
				d2cdy2 = ( lambdaN * dcdyN - lambdaS * dcdyS ) ;
				d2cdz2 = ( lambdaF * dcdzF - lambdaB * dcdzB ) ;

	/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// Convective velocities on the faces
				//--------------------------------------------------------------

#ifdef VOF_VELOCITY
				uE = u_vof[k][j][i+1];
				uW = u_vof[k][j][i];

				vN = v_vof[k][j+1][i];
				vS = v_vof[k][j][i];

				wF = w_vof[k+1][j][i];
				wB = w_vof[k][j][i];


#else
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];

				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];
#endif


				//--------------------------------------------------------------
				// Concentration value on the faces
				//--------------------------------------------------------------
#ifdef CONC_QUICK

				if (uE > 0) {
					if (i==0) {
#ifdef XPERIODIC_CONC
						cE = awWW[i]*c_data[k][j][i-1] + awW[i]*c_data[k][j][i] + awE[i]*c_data[k][j][i+1];
#else
						cE = c_data[k][j][i];
#endif
					}
					else {
						cE = awWW[i]*c_data[k][j][i-1] + awW[i]*c_data[k][j][i] + awE[i]*c_data[k][j][i+1];
					}
				}
				else {
					if (i==NX-2) {
#ifdef XPERIODIC_CONC
						cE = aeW[i]*c_data[k][j][i] + aeE[i]*c_data[k][j][i+1] + aeEE[i]*c_data[k][j][i+2];
#else
						cE = c_data[k][j][i+1];
#endif
					}
					else if ( (grid->c_status[k][j][i+2] == FLUID) || (grid->c_status[k][j][i+2] == IMMERSED) ) {
						cE = aeW[i]*c_data[k][j][i] + aeE[i]*c_data[k][j][i+1] + aeEE[i]*c_data[k][j][i+2];
					}
					else {
						cE = c_data[k][j][i+1];
					}
				}

				// Calculate  cW
				if (i==0) {
#ifdef XPERIODIC_CONC
					if (uW > 0) {
						cW = awWW[0]*c_data[k][j][i-2] + awW[0]*c_data[k][j][i-1] + awE[0]*c_data[k][j][i];
					}
					else {
						cW = aeW[0]*c_data[k][j][i-1] + aeE[0]*c_data[k][j][i] + aeEE[0]*c_data[k][j][i+1];
					}
#else // not XPERIODIC
					cW = c_data[k][j][-1];

#endif // not XPERIODIC
				}

				else {
					if (uW > 0) {
						if (i==1) {
#ifdef XPERIODIC_CONC
							cW = awWW[i-1]*c_data[k][j][i-2] + awW[i-1]*c_data[k][j][i-1] + awE[i-1]*c_data[k][j][i];
#else
							cW = c_data[k][j][i-1];
#endif
						}
						else if ( (grid->c_status[k][j][i-2] == FLUID) || (grid->c_status[k][j][i-2] == IMMERSED) ) {
							cW = awWW[i-1]*c_data[k][j][i-2] + awW[i-1]*c_data[k][j][i-1] + awE[i-1]*c_data[k][j][i];
						}
						else {
							cW = c_data[k][j][i-1];
						}
					}
					else {
						cW = aeW[i-1]*c_data[k][j][i-1] + aeE[i-1]*c_data[k][j][i] + aeEE[i-1]*c_data[k][j][i+1];
					}
				}

				if (vN > 0) {
					if (j==0) {
						// Near the boundary use first order upwind
						cN = c_data[k][j][i];
					}
					else {
						cN = asSS[j]*c_data[k][j-1][i] + asS[j]*c_data[k][j][i] + asN[j]*c_data[k][j+1][i];
					}
				}
				else {
					if (j == NY-2) {
						// Near the boundary use first order upwind
						cN = c_data[k][j+1][i];
					}
					else if ( (grid->c_status[k][j+2][i] == FLUID) || (grid->c_status[k][j+2][i] == IMMERSED) ) {
						// Use QUICK when the stencil  permits (i.e all the
						// points have physical value)
						cN = anS[j]*c_data[k][j][i] + anN[j]*c_data[k][j+1][i] + anNN[j]*c_data[k][j+2][i];
					}
					else {
						// If QUICK stencil has solid nodes, use first order
						// upwind
						cN = c_data[k][j+1][i];
					}
				}

				if (vS > 0) {
					if (j==0) {
						// Impose Neumann b/c i.e dc/dy = 0.0
						cS = c_data[k][j][i];
					}
					else if (j==1) {
						cS = c_data[k][j-1][i];
					}
					else  if ( (grid->c_status[k][j-2][i] == FLUID) || (grid->c_status[k][j-2][i] == IMMERSED) ) {
						// Use QUICK when the stencil  permits (i.e all the
						// points have physical value)
						cS = asSS[j-1]*c_data[k][j-2][i] + asS[j-1]*c_data[k][j-1][i] + asN[j-1]*c_data[k][j][i];
					}
					else {
						// If QUICK stencil has solid nodes, use first order
						// upwind
						cS = c_data[k][j-1][i];
					}
				}
				else {
					if (j==0) {
						// Impose Neumann b/c i.e dc/dy = 0.0
						cS = c_data[k][j][i];
					}
					else {
						cS = anS[j-1]*c_data[k][j-1][i] + anN[j-1]*c_data[k][j][i] + anNN[j-1]*c_data[k][j+1][i];
					}
				}

				// Calculate  cF
				if (wF > 0) {
					if (k==0) {
#ifdef ZPERIODIC_CONC
						cF = abBB[k]*c_data[k-1][j][i] + abB[k]*c_data[k][j][i] + abF[k]*c_data[k+1][j][i];
#else
						cF = c_data[k][j][i];
#endif
					}
					else {
						cF = abBB[k]*c_data[k-1][j][i] + abB[k]*c_data[k][j][i] + abF[k]*c_data[k+1][j][i];
					}
				}
				else {
					if (k==NZ-2) {
#ifdef ZPERIODIC_CONC
						cF = afB[k]*c_data[k][j][i] + afF[k]*c_data[k+1][j][i] + afFF[k]*c_data[k+2][j][i];
#else
						cF = c_data[k+1][j][i];
#endif
					}
					else if ( (grid->c_status[k+2][j][i] == FLUID) || (grid->c_status[k+2][j][i] == IMMERSED) ) {
						cF = afB[k]*c_data[k][j][i] + afF[k]*c_data[k+1][j][i] + afFF[k]*c_data[k+2][j][i];
					}
					else {
						cF = c_data[k+1][j][i];
					}
				}

				// Calculate  cB
				if (k==0) {
#ifdef ZPERIODIC_CONC
					if (wB > 0) {
						cB = abBB[0]*c_data[k-2][j][i] + abB[0]*c_data[k-1][j][i] + abF[0]*c_data[k][j][i];
					}
					else {
						cB = afB[0]*c_data[k-1][j][i] + afF[0]*c_data[k][j][i] + afFF[0]*c_data[k+1][j][i];
					}
#else
					// Otherwise impose Neumann b/c i.e dc/dx = 0.0
					cB = c_data[k][j][i];
#endif
				}

				else {

					if (wB > 0) {
						if (k==1) {
#ifdef ZPERIODIC_CONC
							cB = abBB[k-1]*c_data[k-2][j][i] + abB[k-1]*c_data[k-1][j][i] + abF[k-1]*c_data[k][j][i];
#else
							cB = c_data[k-1][j][i];
#endif
						}
						else if ( (grid->c_status[k-2][j][i] == FLUID) || (grid->c_status[k-2][j][i] == IMMERSED) ) {
							cB = abBB[k-1]*c_data[k-2][j][i] + abB[k-1]*c_data[k-1][j][i] + abF[k-1]*c_data[k][j][i];
						}
						else {
							cB = c_data[k-1][j][i];
						}
					}
					else {
						cB = afB[k-1]*c_data[k-1][j][i] + afF[k-1]*c_data[k]
																		  [j][i] + afFF[k-1]*c_data[k+1][j][i];
					}
				}

#else





				cE = 0.5 * ( c_data[k][j][i] + c_data[k][j][i+1] );
				cW = 0.5 * ( c_data[k][j][i] + c_data[k][j][i-1] );

				cN = 0.5 * ( c_data[k][j][i] + c_data[k][j+1][i] );
				cS = 0.5 * ( c_data[k][j][i] + c_data[k][j-1][i] );

				cF = 0.5 * ( c_data[k][j][i] + c_data[k+1][j][i] );
				cB = 0.5 * ( c_data[k][j][i] + c_data[k-1][j][i] );
#endif


#if !defined VOF_SCALAR_NODIV
				ucE = uE*cE;
				ucW = uW*cW;

				vcN = vN*cN;
				vcS = vS*cS;

				wcF = wF*cF;
				wcB = wB*cB;

#else
				ucE = 0.5*(uE+uW)*cE;
				ucW = 0.5*(uE+uW)*cW;

				vcN = 0.5*(vN+vS)*cN;
				vcS = 0.5*(vN+vS)*cS;

				wcF = 0.5*(wB+wF)*cF;
				wcB = 0.5*(wB+wF)*cB;


#endif
				//--------------------------------------------------------------
				// d/dx(uc), d/dy(vc), and d/dz(wc)
				//--------------------------------------------------------------
				dcudx = ( ucE - ucW ) * ih;
				dcvdy = ( vcN - vcS ) * ih;
				dcwdz = ( wcF - wcB ) * ih;




				/*------------------------------------------------------------*/
				/*
				 Store convective and viscous terms
				 */
				/*------------------------------------------------------------*/
				conv[k][j][i] = dcudx + dcvdy + dcwdz;


			//	printf("Index (%d, %d, %d ) my ranks is %d \n", i,j,k, params->rank); fflush(stdout);

#if defined CONC_FULLY_IMPLICIT
				visc[k][j][i] = ihsq*iPe*(d2cdx2 + d2cdy2 + d2cdz2)/( ng_vfc[k][j][i]*(heat_cap_s-1)+1 );
				visc_explicit[k][j][i] = 0.0;
#else
		#error  choose the fully implicity diffusion scheme
#endif
			}
		}
	}




	return;

}
#endif
