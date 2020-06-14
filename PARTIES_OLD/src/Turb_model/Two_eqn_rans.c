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
#include "Strain.h"
#include "Cart3d.h"
#include "Rans.h"
#include "Two_eqn_rans.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/******************************************************************************/
/*
 This function allocates enough memory for the Two_equation_rans structure based on
 parameters defined in "*params"
 */
/******************************************************************************/
Two_equation_rans *Two_equation_rans_create(MAC_grid *grid, Parameters *params) {

	Two_equation_rans *rans;

	rans = (Two_equation_rans *) malloc(sizeof(Two_equation_rans));
	Memory_check_allocation(rans);

	rans->data = Memory_allocate_flow_variable(grid, params);

	rans->ng_explicit = Memory_allocate_noghost_variable(grid, params);
	rans->ng_implicit = Memory_allocate_noghost_variable(grid, params);
	rans->ng_rhs = Memory_allocate_noghost_variable(grid, params);

#ifdef RANS_QUICK
	Rans_set_quick_coefficients(rans, grid, params);
#endif

/*
	new_rans->solve_cpu_time = 0.0;
	new_rans->solve_first_copy_cpu_time = 0.0;
	new_rans->solve_second_copy_cpu_time = 0.0;
	new_rans->solve_first_remap_cpu_time = 0.0;
	new_rans->solve_second_remap_cpu_time = 0.0;
	new_rans->solve_tridiag_cpu_time = 0.0;
*/


	return rans;
}




/******************************************************************************/
/*
 This function releases the allocated memory for Rans structure.
 */
/******************************************************************************/
void Two_equation_rans_destroy(Two_equation_rans *two_eqn_rans,
	Parameters *params, MAC_grid *grid) {

	Memory_free_flow_variable(grid, params, two_eqn_rans->data);
	Memory_free_noghost_variable(grid, params, two_eqn_rans->ng_explicit);
	Memory_free_noghost_variable(grid, params, two_eqn_rans->ng_implicit);
	Memory_free_noghost_variable(grid, params, two_eqn_rans->ng_rhs);

	free(two_eqn_rans);

}



/******************************************************************************/
/*
 This function sets the initial lock profile
 */
/******************************************************************************/
void Rans_initialize (Rans *rans, Velocity *u, Cart3d_bag *data_bag) {

	int i, j, k;
	int NI, NJ, NK;
	int lock_smooth;
	double x_fr, y_fr, z_fr;
	double *xc, *yc, *zc;
	double ***conc;
	double Re;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NY;
	Two_equation_rans **two_eqn_rans;
	double ***u_data;
	double tke, diss, ywall, kappa, C_mu;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Re = params->Re;

	xc   = grid->xc;
	yc   = grid->yc;
	zc   = grid->zc;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	NY = grid->NY;
	kappa = rans->kappa;
	C_mu = rans->C_mu;

	two_eqn_rans = rans->two_eqn_rans;
	u_data = u->data;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			ywall = grid->yc[j]-grid->yv[0];
#ifndef TOP_WALL_VELOCITY_FREESLIP
			ywall = min(ywall, grid->yv[NY-1]-grid->yc[j]);
#endif
			for (i=Is; i<Ie; i++) {

				//10% velocity as u_rms to calculate the TKE
				tke = 0.01*u_data[k][j][i]*u_data[k][j][i];
				if (ywall < 0.1) {
					tke = 1.0;
					diss = 0.25;
				}
				else if (ywall < 0.25) {
					tke = 1.0 - 0.5*(ywall-0.2)/0.3;
					diss = 0.25 - 0.05*(ywall-0.2)/0.3;
				}
				else {
					tke = 0.5;
					diss = 0.2;
				}
				two_eqn_rans[0]->data[k][j][i] = tke/400.;
				two_eqn_rans[1]->data[k][j][i] = pow(tke,1.5)*pow(C_mu,0.75)/(kappa*ywall);
				two_eqn_rans[1]->data[k][j][i] = diss/400.;
//				two_eqn_rans[0]->data[k][j][i] = u_data[k][j][i];
//				two_eqn_rans[1]->data[k][j][i] = u_data[k][j][i];
				if ( (k==2) && (i==2) ) {
					printf("Diss %e %e %e %e %e %e %e\n",grid->yc[j],tke,two_eqn_rans[1]->data[k][j][i],
							pow(tke,1.5),pow(C_mu,0.75),kappa,ywall);
				}

			}
		}
	}
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[0]->data,
										'c', params->ghost_nodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(rans->two_eqn_rans[1]->data,
										'c', params->ghost_nodes, data_bag);


}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Rans_set_conv_viscous_central( Velocity *u, Velocity *v, Velocity *w,
		Rans *rans, MAC_grid *grid, Parameters *params, int eqn_no) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	double iRe;
	double nutE, nutW, nutN, nutS, nutF, nutB;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;
	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double ***nut;
	double ***u_data, ***v_data, ***w_data;
	double ***c_data;
	double ***explicit, ***explicit_old, ***implicit;
	double sigma;
	Two_equation_rans *two_eqn_rans;
	FILE *fid;
	char fname[50];

	int ntime = params -> ntime;

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

	iRe = 1.e0/params->Re;

	idx_u = grid->idx_u;
	idy_v = grid->idy_v;
	idz_w = grid->idz_w;
	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;


	if (eqn_no == 1) two_eqn_rans = rans->two_eqn_rans[0];
	if (eqn_no == 2) two_eqn_rans = rans->two_eqn_rans[1];
	sprintf(fname, "conv_viscos_eq%d_r%d_s%d_ts%d.dat",eqn_no, params->rank,params->which_stage,ntime);
	int ifreq=100000;
//	ifreq = 1;
	if (ntime%ifreq ==0) fid = fopen(fname,"w");

	// Get the local velocities at the location where they are defined
	u_data = u->data;
	v_data = v->data;
	w_data = w->data;
	c_data = two_eqn_rans->data;

	nut = rans->nut;
	if (eqn_no == 1) sigma = rans->sigma_1;
	if (eqn_no == 2) sigma = rans->sigma_2;

	explicit = two_eqn_rans->ng_explicit;
	implicit = two_eqn_rans->ng_implicit;


	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

				// Eddy viscosity
				nutE = 0.5*(nut[k][j][i+1] + nut[k][j][i]);
		#ifdef XPERIODIC
				nutW = 0.5*(nut[k][j][i-1] + nut[k][j][i]);
		#else
				if (i != 0) {
					nutW = 0.5*(nut[k][j][i-1] + nut[k][j][i]);
				}
				else {
					nutW = nut[k][j][i];
				}
		#endif // XPERIODIC
				nutN = 0.5*(nut[k][j+1][i] + nut[k][j][i]);
				if (j == NY-2)  nutN = 0.0;
				if (j !=0 ) {
					nutS = 0.5*(nut[k][j-1][i] + nut[k][j][i]);
				}
				else {
					nutS = 0.0;
				}
				nutF = 0.5*(nut[k+1][j][i] + nut[k][j][i]);
		#ifdef ZPERIODIC
				nutB = 0.5*(nut[k-1][j][i] + nut[k][j][i]);
		#else
				if (k !=0 ) {
					nutB = 0.5*(nut[k-1][j][i] + nut[k][j][i]);
				}
				else {
					nutB = 0.0;
				}
		#endif // ZPERIODIC
				nutE = nutE/sigma + iRe;
				nutW = nutW/sigma + iRe;
				nutN = nutN/sigma + iRe;
				nutS = nutS/sigma + iRe;
				nutF = nutF/sigma + iRe;
				nutB = nutB/sigma + iRe;


				// d2c/dx2, d2c/dy2 and d2c/dz2
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * idx_c[i];

				if (i != 0 ) {
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i-1];
				}
				else {
	#ifdef XPERIODIC
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i];
	#else
		#ifdef INFLOW
					dcdxW = ( c_data[k][j][i] - c_data[k][j][-1] ) * idx_c[i];
		#else
					dcdxW = 0.;
		#endif
	#endif // XPERIODIC
				}

				dcdyN = ( c_data[k][j+1][i] - c_data[k][j][i] ) * idy_c[j];
//if (j==NY-2) dcdyN = rans->log_law->dudy_wm_top[k][i];

				if ( j != 0) {
					dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] ) * idy_c[j-1];
				}
				else {
					// Neumann b/c dcdy = 0
					dcdyS = 0.0;
	#ifdef THERMAL_KADER
					dcdyS = 2.*c_data[k][j][i] * idy_v[j];
	#endif
//dcdyS = rans->log_law->dudy_wm_bottom[k][i];
				}

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

				d2cdx2 = ( nutE * dcdxE - nutW * dcdxW ) * idx_u[i];
				d2cdy2 = ( nutN * dcdyN - nutS * dcdyS ) * idy_v[j];
				d2cdz2 = ( nutF * dcdzF - nutB * dcdzB ) * idz_w[k];


				// Calculate the convective terms
				// Convective velocity on the faces
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];

				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];

				// Concentration value on the faces
				cE = 0.5*(c_data[k][j][i] + c_data[k][j][i+1]);

				if (i != 0 ) {
					cW = 0.5*(c_data[k][j][i] + c_data[k][j][i-1]);
				}
				else {
	#ifdef XPERIODIC
					cW = 0.5*(c_data[k][j][i] + c_data[k][j][i-1]);
	#else
		#ifdef INFLOW
					cW = 0.5*(c_data[k][j][i] + c_data[k][j][i-1]);
		#else
					// Otherwise impose Neumann b/c i.e dc/dx = 0.0
					cW = c_data[k][j][i];
		#endif
	#endif // XPERIODIC
				}

				cN = 0.5*(c_data[k][j][i] + c_data[k][j+1][i]);

				if (j != 0 ) {
					cS = 0.5*(c_data[k][j][i] + c_data[k][j-1][i]);
				}
				else {
					// Impose Neumann b/c i.e dc/dy = 0.0
					cS = c_data[k][j][i];
				}

				cF = 0.5*(c_data[k][j][i] + c_data[k+1][j][i]);

				if (k != 0) {
					cB = 0.5*(c_data[k][j][i] + c_data[k-1][j][i]);
				}
				else {
	#ifdef ZPERIODIC
					cB = 0.5*(c_data[k][j][i] + c_data[k-1][j][i]);
	#else
					// Otherwise impose Neumann b/c i.e dc/dz = 0.0
					cB = c_data[k][j][i];
	#endif
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

	#ifdef csolver_PETSC
				explicit[k][j][i] = dcudx + dcvdy + dcwdz;
				implicit[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
	#else
		#ifdef FULLY_EXPLICIT
				explicit[k][j][i] = dcudx + dcvdy + dcwdz - (d2cdx2 + d2cdy2 + d2cdz2);
				implicit[k][j][i] = 0.0;
		#else
				explicit[k][j][i] = dcudx + dcvdy + dcwdz - (d2cdx2 + d2cdz2);
				implicit[k][j][i] = d2cdy2;
		#endif
	#endif
				if ( (ntime%ifreq==0) && (k == 2))
					fprintf(fid,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",grid->yc[j], dcudx, dcvdy,
           dcwdz, d2cdx2, d2cdy2, d2cdz2,explicit[k][j][i],implicit[k][j][i],uE,uW,cE,cW,ucE,ucW);

			}
		}
	}
	if (ntime%ifreq ==0) fclose(fid);

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Rans_set_conv_viscous_quick( Velocity *u, Velocity *v, Velocity *w,
		Rans *rans, MAC_grid *grid, Parameters *params, int eqn_no) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	double iRe;
	double nutE, nutW, nutN, nutS, nutF, nutB;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;

	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double ***nut;
	double ***u_data, ***v_data, ***w_data;
	double ***c_data;
	double ***explicit, ***explicit_old, ***implicit;
	double *aeW, *aeE, *aeEE;
	double *awWW, *awW, *awE;
	double *anS, *anN, *anNN;
	double *asSS, *asS, *asN;
	double *afB, *afF, *afFF;
	double *abBB, *abB, *abF;
	double sigma;
	Two_equation_rans *two_eqn_rans;

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

	iRe = 1.e0/params->Re;

	idx_u = grid->idx_u;
	idy_v = grid->idy_v;
	idz_w = grid->idz_w;
	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;

	if (eqn_no == 1) two_eqn_rans = rans->two_eqn_rans[0];
	if (eqn_no == 2) two_eqn_rans = rans->two_eqn_rans[1];

	aeW  = two_eqn_rans->aeW;
	aeE  = two_eqn_rans->aeE;
	aeEE = two_eqn_rans->aeEE;
	awWW = two_eqn_rans->awWW;
	awW  = two_eqn_rans->awW;
	awE  = two_eqn_rans->awE;
	anS  = two_eqn_rans->anS;
	anN  = two_eqn_rans->anN;
	anNN = two_eqn_rans->anNN;
	asSS = two_eqn_rans->asSS;
	asS  = two_eqn_rans->asS ;
	asN  = two_eqn_rans->asN ;
	afB  = two_eqn_rans->afB;
	afF  = two_eqn_rans->afF;
	afFF = two_eqn_rans->afFF;
	abBB = two_eqn_rans->abBB;
	abB  = two_eqn_rans->abB ;
	abF  = two_eqn_rans->abF ;

	// Get the local velocities at the location where they are defined
	u_data = u->data;
	v_data = v->data;
	w_data = w->data;
	c_data = two_eqn_rans->data;

	nut = rans->nut;
	if (eqn_no == 1) sigma = rans->sigma_1;
	if (eqn_no == 2) sigma = rans->sigma_2;

	explicit = two_eqn_rans->ng_explicit;
	implicit = two_eqn_rans->ng_implicit;

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

				// Eddy viscosity
				nutE = 0.5*(nut[k][j][i+1] + nut[k][j][i]);
		#ifdef XPERIODIC
				nutW = 0.5*(nut[k][j][i-1] + nut[k][j][i]);
		#else
				if (i != 0) {
					nutW = 0.5*(nut[k][j][i-1] + nut[k][j][i]);
				}
				else {
					nutW = nut[k][j][i];
				}
		#endif // XPERIODIC
				nutN = 0.5*(nut[k][j+1][i] + nut[k][j][i]);

				if (j !=0 ) {
					nutS = 0.5*(nut[k][j-1][i] + nut[k][j][i]);
				}
				else {
					nutS = 0.0;
				}

				nutF = 0.5*(nut[k+1][j][i] + nut[k][j][i]);
		#ifdef ZPERIODIC
				nutB = 0.5*(nut[k-1][j][i] + nut[k][j][i]);
		#else
				if (k !=0 ) {
					nutB = 0.5*(nut[k-1][j][i] + nut[k][j][i]);
				}
				else {
					nutB = 0.0;
				}
		#endif // ZPERIODIC
				nutE = nutE/sigma + iRe;
				nutW = nutW/sigma + iRe;
				nutN = nutN/sigma + iRe;
				nutS = nutS/sigma + iRe;
				nutF = nutF/sigma + iRe;
				nutB = nutB/sigma + iRe;

				// d2c/dx2, d2c/dy2 and d2c/dz2
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * idx_c[i];

				if (i != 0 ) {
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i-1];
				}
				else {
	#ifdef XPERIODIC
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i];
	#else
		#ifdef INFLOW
					dcdxW = ( c_data[k][j][i] - c_data[k][j][-1] ) * idx_c[i];
		#else
					dcdxW = 0.;
		#endif
	#endif // XPERIODIC
				}

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

				d2cdx2 = ( nutE * dcdxE - nutW * dcdxW ) * idx_u[i];
				d2cdy2 = ( nutN * dcdyN - nutS * dcdyS ) * idy_v[j];
				d2cdz2 = ( nutF * dcdzF - nutB * dcdzB ) * idz_w[k];


				// Calculate the convective terms
				// Convective velocity on the faces
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];

				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];

				// Concentration value on the faces
				// For explanation of various "if conditions", check cN and cS
				// calculation, which is representative of cE and cW, cF and cB
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
	#endif // XPERIODIC
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

				dcudx = ( ucE - ucW ) * idx_u[i];
				dcvdy = ( vcN - vcS ) * idy_v[j];
				dcwdz = ( wcF - wcB ) * idz_w[k];
				//printf("rank = %d %d %d %d \n",params->rank,k,j,i);

	#ifdef csolver_PETSC
				explicit[k][j][i] = dcudx + dcvdy + dcwdz;
				implicit[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;





	#else
		#ifdef CONC_FULLY_EXPLICIT
				explicit[k][j][i] = dcudx + dcvdy + dcwdz - (d2cdx2 + d2cdy2 + d2cdz2);
				implicit[k][j][i] = 0.0;
		#elif defined  CONC_SEMI_IMPLICIT
				explicit[k][j][i] = dcudx + dcvdy + dcwdz - (d2cdx2 + d2cdz2);
				implicit[k][j][i] = d2cdy2;
		#elif  defined CONC_FULLY_IMPLICIT
				printf("RANS  NOOOOOOOOOOOOOOOOOOOOOOOOOOOTTTTTTTTTTT Implemented for fully implicit conc!!!!!!!");
		#endif
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
void Rans_set_conv_viscous_ftupwind( Velocity *u, Velocity *v, Velocity *w,
		Rans *rans, MAC_grid *grid, Parameters *params, int eqn_no) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	double iRe;
	double nutE, nutW, nutN, nutS, nutF, nutB;
	double dcdxE, dcdxW, dcdyN, dcdyS, dcdzF, dcdzB;
	double d2cdx2, d2cdy2, d2cdz2;
	double uE, uW, vN, vS, wF, wB;
	double cE, cW, cN, cS, cF, cB;
	double ucE, ucW, vcN, vcS, wcF, wcB;
	double dcudx, dcvdy, dcwdz;
	double rhs;

	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double ***nut;
	double ***u_data, ***v_data, ***w_data;
	double ***c_data;
	double ***explicit, ***explicit_old, ***implicit;
	double *xc, *yc, *zc;
	double *xu, *yv, *zw;
	double xWW, xW, xE, xEE, xface;
	double ySS, yS, yN, yNN, yface;
	double zBB, zB, zF, zFF, zface;
	double sigma;
	Two_equation_rans *two_eqn_rans;

	FILE *fid;
	char bin_filename[50];

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

	iRe = 1.e0/params->Re;

	idx_u = grid->idx_u;
	idy_v = grid->idy_v;
	idz_w = grid->idz_w;
	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;
	xc = grid->xc;
	yc = grid->yc;
	zc = grid->zc;
	xu = grid->xu;
	yv = grid->yv;
	zw = grid->zw;

	if (eqn_no == 1) two_eqn_rans = rans->two_eqn_rans[0];
	if (eqn_no == 2) two_eqn_rans = rans->two_eqn_rans[1];

	// Get the local velocities at the location where they are defined
	u_data = u->data;
	v_data = v->data;
	w_data = w->data;
	c_data = two_eqn_rans->data;

	nut = rans->nut;
	if (eqn_no == 1) sigma = rans->sigma_1;
	if (eqn_no == 2) sigma = rans->sigma_2;

	explicit = two_eqn_rans->ng_explicit;
	implicit = two_eqn_rans->ng_implicit;

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				// Eddy viscosity
				nutE = 0.5*(nut[k][j][i+1] + nut[k][j][i]);
		#ifdef XPERIODIC
				nutW = 0.5*(nut[k][j][i-1] + nut[k][j][i]);
		#else
				if (i != 0) {
					nutW = 0.5*(nut[k][j][i-1] + nut[k][j][i]);
				}
				else {
					nutW = nut[k][j][i];
				}
		#endif
				nutN = 0.5*(nut[k][j+1][i] + nut[k][j][i]);
				if (j !=0 ) {
					nutS = 0.5*(nut[k][j-1][i] + nut[k][j][i]);
				}
				else {
					nutS = 0.0;
				}
				nutF = 0.5*(nut[k+1][j][i] + nut[k][j][i]);
		#ifdef ZPERIODIC
				nutB = 0.5*(nut[k-1][j][i] + nut[k][j][i]);
		#else
				if (k !=0 ) {
					nutB = 0.5*(nut[k-1][j][i] + nut[k][j][i]);
				}
				else {
					nutB = 0.0;
				}
		#endif // ZPERIODIC
				nutE = nutE/sigma + iRe;
				nutW = nutW/sigma + iRe;
				nutN = nutN/sigma + iRe;
				nutS = nutS/sigma + iRe;
				nutF = nutF/sigma + iRe;
				nutB = nutB/sigma + iRe;

				// d2c/dx2, d2c/dy2 and d2c/dz2
				dcdxE = ( c_data[k][j][i+1] - c_data[k][j][i] ) * idx_c[i];
				if (i != 0 ) {
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i-1];
				}
				else {
	#ifdef XPERIODIC
					dcdxW = ( c_data[k][j][i] - c_data[k][j][i-1] ) * idx_c[i];
	#else
		#ifdef INFLOW
					dcdxW = ( c_data[k][j][i] - c_data[k][j][-1] ) * idx_c[i];
		#else
					dcdxW = 0.;
		#endif
	#endif // XPERIODIC
				}

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

				d2cdx2 = ( nutE * dcdxE - nutW * dcdxW ) * idx_u[i];
				d2cdy2 = ( nutN * dcdyN - nutS * dcdyS ) * idy_v[j];
				d2cdz2 = ( nutF * dcdzF - nutB * dcdzB ) * idz_w[k];


				// Calculate the convective terms
				// Convective velocity on the faces
				uE = u_data[k][j][i+1];
				uW = u_data[k][j][i];

				vN = v_data[k][j+1][i];
				vS = v_data[k][j][i];


				wF = w_data[k+1][j][i];
				wB = w_data[k][j][i];

				// Concentration value on the faces
				// For explanation of various "if conditions", check cN and cS
				// calculation, which is representative of cE and cW, cF and cB

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

				dcudx = ( ucE - ucW ) * idx_u[i];
				dcvdy = ( vcN - vcS ) * idy_v[j];
				dcwdz = ( wcF - wcB ) * idz_w[k];

	#ifdef csolver_PETSC
				explicit[k][j][i] = dcudx + dcvdy + dcwdz;
				implicit[k][j][i] = d2cdx2 + d2cdy2 + d2cdz2;
	#else
		#ifdef CONC_FULLY_EXPLICIT
				explicit[k][j][i] = dcudx + dcvdy + dcwdz - (d2cdx2 + d2cdy2 + d2cdz2);
				implicit[k][j][i] = 0.0;
		#elif defined CONC_SEMI_IMPLICIT
				explicit[k][j][i] = dcudx + dcvdy + dcwdz - (d2cdx2 + d2cdz2);
				implicit[k][j][i] = d2cdy2;
		#elif defined CONC_FULLY_IMPLICIT
				printf("RANS  NOOOOOOOOOOOOOOOOOOOOOOOOOOOTTTTTTTTTTT Implemented for fully implicit conc!!!!!!!");
		#endif
	#endif


			}
		}
	}

	return;
}




/******************************************************************************/
/*
 This function computes the RHS of the linear system. Convective and some(or all) viscous
 terms are treated explicitly
 */
/******************************************************************************/
void Rans_set_RHS(Two_equation_rans *two_eqn_rans, MAC_grid *grid,
		Parameters *params, int eqn_no, double dt) {

	int i, j, k;
	double ***explicit, ***explicit_old, ***implicit;
	double ***data;
	double ***rhs_vec;
	double a_dt;
	double rhs_value;
	double value;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	int status;
	const double BET[] = {BETA};
	const double GAMB[] = {GAMBETA};
	const double ZETB[] = {ZETBETA};
	int NX, NY, NZ;
	FILE *fid;
	char fname[50];

	// Pointer to a function
	data = two_eqn_rans->data;

	// Now, got the RHS vector on current processor
	rhs_vec = two_eqn_rans->ng_rhs;

	explicit = two_eqn_rans->ng_explicit;
	implicit = two_eqn_rans->ng_implicit;

	int ntime   = params -> ntime;
	int which_stage = params -> which_stage;
	a_dt = 1.0 / ( dt * BET[which_stage] );


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
	i_end = min(NX-1, Ie);
	j_end = min(NY-1, Je);
	k_end = min(NZ-1, Ke);

	sprintf(fname, "rhs_eq%d_r%d_s%d_ts%d.dat",eqn_no, params->rank,which_stage,ntime);
	int ifreq=100000;
	ifreq=100000;
	if (ntime%ifreq ==0) fid = fopen(fname,"w");

	double sum = 0;
	if (which_stage == 0) {
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					rhs_vec[k][j][i] = data[k][j][i]*a_dt
						- GAMB[0]*explicit[k][j][i] + implicit[k][j][i];
					explicit_old[k][j][i] = explicit[k][j][i];
				if ( (i==2) && (k==2) ) sum += GAMB[0]*explicit[k][j][i];
				if ( (ntime%ifreq==0) && (k == 2))
					fprintf(fid, "%e %e %e %e %e %e %e\n",grid->yc[j],data[k][j][i],
		explicit[k][j][i], implicit[k][j][i], rhs_vec[k][j][i],2*a_dt,sum);
				}// for i
			} // for
		}// for k
	}
	else if (which_stage == 1) {
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					rhs_vec[k][j][i] = data[k][j][i]*a_dt
						- GAMB[1]*explicit[k][j][i] - ZETB[1]*explicit_old[k][j][i] + implicit[k][j][i];
					explicit_old[k][j][i] = explicit[k][j][i];
				}// for i
			} // for
		}// for k
	}
	else if (which_stage == 2) {
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					rhs_vec[k][j][i] = data[k][j][i]*a_dt
						- GAMB[2]*explicit[k][j][i] - ZETB[2]*explicit_old[k][j][i] + implicit[k][j][i];
				}// for i
			} // for
		}// for k
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
	if (ntime%ifreq ==0) fclose(fid);

}

/******************************************************************************/
/*
 This function computes the RHS of the u-momentum linear system.  Convective
 terms are treated explicitly
 */
/******************************************************************************/
void Rans_add_source_RHS(Cart3d_bag *data_bag) {

	int i,j, k;
	double ***tke, ***diss;
	double ***rhs_tke, ***rhs_diss;
	double ***nut, ***nutc;
	double ***strain;
	double ***c_total;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	int status;
	int NX, NY, NZ;
	double C_epsilon1, C_epsilon2, C_epsilon3;
	double prod, b_prod, dissbytke;
	double cN, cS, dcdy;
	double *idy_v;
	double **dudy_wm_b, **dwdy_wm_b;
	double **dudy_wm_t, **dwdy_wm_t;
	double dudy, dwdy, dVdy;
	double iRe, prod_denom_b, prod_denom_t, ywall;
	int NConc;
	FILE *fid;
	char fname[50];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Rans *rans  = data_bag -> rans;
	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
	Concentration **c = data_bag -> c;

	int ntime = params -> ntime;

	Strain_rate_magnitude(u, v, w, grid, params, rans->st_rate);

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


	idy_v = grid->idy_v;

	tke  = rans->two_eqn_rans[0]->data;
	diss = rans->two_eqn_rans[1]->data;

	nut = rans->nut;
	nutc = rans->nutc;
	strain = rans->st_rate->strain;
	dudy_wm_b = rans->log_law->dudy_wm_bottom;
	dwdy_wm_b = rans->log_law->dwdy_wm_bottom;
	dudy_wm_t = rans->log_law->dudy_wm_top;
	dwdy_wm_t = rans->log_law->dwdy_wm_top;

	rhs_tke  = rans->two_eqn_rans[0]->ng_rhs;
	rhs_diss = rans->two_eqn_rans[1]->ng_rhs;

	b_prod  = 0.0;

#ifdef CONC
	NConc = params->NConc;
	if (NConc > 1) {
		Communication_update_ghost_nodes_flow_variable(c[0]->c_total, 'c', 1, data_bag);
		c_total = c[0]->c_total;

	}
	else {

		Communication_update_ghost_nodes_flow_variable(c[0]->data, 'c', 1, data_bag);
		c_total = c[0]->data;
	} // else
#endif

	C_epsilon1 = rans->C_epsilon1;
	C_epsilon2 = rans->C_epsilon2;
	C_epsilon3 = rans->C_epsilon3;

	sprintf(fname, "rhs_source%d_s%d_ts%d.dat",params->rank,params->which_stage,ntime);
	int ifreq=100;
	ifreq=100000;
	if (ntime%ifreq ==0) fid = fopen(fname,"w");

	iRe = 1.0/params->Re;

	ywall = grid->yc[0] - grid->yv[0];
	prod_denom_b = pow(rans->C_mu,0.25)*rans->kappa*ywall;

	ywall = grid->yv[NY-1] - grid->yc[NY-2];
	prod_denom_t = pow(rans->C_mu,0.25)*rans->kappa*ywall;


	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

#ifdef CONC
				cN = 0.5*(c_total[k][j][i] + c_total[k][j+1][i]);
				if (j == NY-2) cN = c_total[k][j][i];


				if (j != 0 ) {
					cS = 0.5*(c_total[k][j][i] + c_total[k][j-1][i]);
				}
				else {
					// Impose Neumann b/c i.e dc/dy = 0.0
					cS = c_total[k][j][i];
				}

				dcdy = (cN - cS) * idy_v[j];

				//Buoyancy Production
				b_prod = nutc[k][j][i]*dcdy;
#endif

				// Turbulence Production
				prod = nut[k][j][i]*strain[k][j][i]*strain[k][j][i];

				dissbytke = diss[k][j][i]/tke[k][j][i];




				if (j==0) {
					dudy = 0.5*(dudy_wm_b[k][i] + dudy_wm_b[k][i]);
					dwdy = 0.5*(dwdy_wm_b[k][i] + dwdy_wm_b[k+1][i]);
					dVdy = pow(dudy*dudy + dwdy*dwdy, 0.5);

					prod = iRe*dVdy/prod_denom_b/pow(tke[k][j][i],0.5);
				}

#ifndef TOP_WALL_VELOCITY_FREESLIP
				if (j==NY-2) {
					dudy = 0.5*(dudy_wm_b[k][i] + dudy_wm_b[k][i]);
					dwdy = 0.5*(dwdy_wm_b[k][i] + dwdy_wm_b[k+1][i]);
					dVdy = pow(dudy*dudy + dwdy*dwdy, 0.5);

					prod = iRe*dVdy/prod_denom_t/pow(tke[k][j][i],0.5);
				}
#endif




#ifdef RANS_SOURCE_IMPLICIT
				if ( (j==0) || (j==NY-2) ) {
					rhs_tke[k][j][i] += 2.*(prod + b_prod - diss[k][j][i]);
					rhs_diss[k][j][i] += 2.*dissbytke*( C_epsilon1*prod
								+ C_epsilon3*b_prod - C_epsilon2*diss[k][j][i]);

				}
				else {
				rhs_tke[k][j][i]  += -2.*(prod + b_prod);
				rhs_diss[k][j][i] +=  2.*(2.*C_epsilon2*diss[k][j][i]*dissbytke);

/*					rhs_tke[k][j][i] += 2.*b_prod + (prod - diss[k][j][i]);
					rhs_diss[k][j][i] += 2.*dissbytke*C_epsilon3*b_prod
								+ dissbytke*(C_epsilon1*prod - C_epsilon2*diss[k][j][i]);
					rhs_diss[k][j][i] += 2.*dissbytke*( C_epsilon1*prod
								+ C_epsilon3*b_prod - C_epsilon2*diss[k][j][i]);

*/
				}
#else

				rhs_tke[k][j][i] += 2.*(prod + b_prod - diss[k][j][i]);
				rhs_diss[k][j][i] += 2.*dissbytke*( C_epsilon1*prod
								+ C_epsilon3*b_prod - C_epsilon2*diss[k][j][i]);

#endif
#ifdef IMMERSED_BOUNDARY
				if (grid->c_status[k][j][i] == SOLID) {
					rhs_tke[k][j][i] = 0.0;
					rhs_diss[k][j][i] = 0.0;
				}
#endif

				if ( (ntime%ifreq==0) && (k == 2))
						fprintf(fid,"%e %e %e %e %e %e \n",grid->yc[j], tke[k][j][i],
						diss[k][j][i], prod,nut[k][j][i],strain[k][j][i]);

			}// for i
		} // for j
	}// for k
	if (ntime%ifreq ==0) fclose(fid);

	return;

}



/******************************************************************************/
/*
 */
/******************************************************************************/

void Rans_set_quick_coefficients(Two_equation_rans *two_eqn_rans,
		MAC_grid *grid, Parameters *params) {

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

	two_eqn_rans->aeW  = aeW;
	two_eqn_rans->aeE  = aeE;
	two_eqn_rans->aeEE = aeEE;
	two_eqn_rans->awWW = awWW;
	two_eqn_rans->awW  = awW;
	two_eqn_rans->awE  = awE;
	two_eqn_rans->anS  = anS ;
	two_eqn_rans->anN  = anN;
	two_eqn_rans->anNN = anNN;
	two_eqn_rans->asSS = asSS;
	two_eqn_rans->asS  = asS;
	two_eqn_rans->asN  = asN;
	two_eqn_rans->afB  = afB;
	two_eqn_rans->afF  = afF;
	two_eqn_rans->afFF = afFF;
	two_eqn_rans->abBB = abBB;
	two_eqn_rans->abB  = abB;
	two_eqn_rans->abF  = abF;

	return;
}




/******************************************************************************/
/*
 This function solves the linear system to find the RANS related variables
 */
/******************************************************************************/
int Two_equation_solve(Cart3d_bag *data_bag, double dt) {

	int iters;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int j_start, j_end;
	int i, j, k;
	double ywall, kappa;
	double turbke;
	int NY;
	double C_epsilon1, C_epsilon2, C_epsilon3;
	double C_mu, C_muc;
	double ***tke, ***diss;
	double ***strain;
	double k_coeff, diss_coeff, dissbytke, tkebydiss, dcdy, strainsq;
	double ak_kcoeff, aeps_kcoeff, bk_epscoeff, beps_epscoeff;
	double ***rhs_tke, ***rhs_diss, ***c_total;
	double ***nut;
	double *idy_v;
	double cN, cS;
	int NConc;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Rans *rans = data_bag -> rans;
	Concentration **c = data_bag -> c;

	const double BET[] = {BETA};
	double dtimeb = BET[params -> which_stage] * dt;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	kappa = rans->kappa;
	C_mu = rans->C_mu;
	C_muc = rans->C_muc;

	C_epsilon1 = rans->C_epsilon1;
	C_epsilon2 = rans->C_epsilon2;
	C_epsilon3 = rans->C_epsilon3;

	NY = grid->NY;
	tke = rans->two_eqn_rans[0]->data;
	diss = rans->two_eqn_rans[1]->data;
	nut = rans->nut;

#ifndef RANS_SOURCE_IMPLICIT
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				rans->two_eqn_rans[0]->data[k][j][i] = rans->two_eqn_rans[0]->ng_rhs[k][j][i]*dtimeb;
				rans->two_eqn_rans[1]->data[k][j][i] = rans->two_eqn_rans[1]->ng_rhs[k][j][i]*dtimeb;
			}
		}
	}

	if (Js == 0) {
		j = Js;
		ywall = grid->yc[j]-grid->yv[0];
		for (k=Ks;k<Ke;k++){
			for (i=Is;i<Ie;i++){
				turbke = rans->two_eqn_rans[0]->data[k][j][i];
				rans->two_eqn_rans[1]->data[k][j][i] = pow(turbke,1.5)*pow(C_mu,0.75)/(kappa*ywall);
			}
		}
	}

	if (Je == NY) {
		j = NY-2;
		ywall = grid->yv[NY-1]-grid->yc[NY-2];
		for (k=Ks;k<Ke;k++){
			for (i=Is;i<Ie;i++){
				turbke = rans->two_eqn_rans[0]->data[k][j][i];
				rans->two_eqn_rans[1]->data[k][j][i] = pow(turbke,1.5)*pow(C_mu,0.75)/(kappa*ywall);
			}
		}
	}

#else  // RANS_SOURCE_IMPLICIT


	j_start = Js;
	j_end = Je;
	if (Js == 0) j_start = 1;
	if (Je == NY) j_end = NY-1;

	strain = rans->st_rate->strain;
	rhs_tke  = rans->two_eqn_rans[0]->ng_rhs;
	rhs_diss = rans->two_eqn_rans[1]->ng_rhs;

	#ifdef CONC
	NConc = params->NConc;
	if (NConc > 1) {
		Communication_update_ghost_nodes_flow_variable(c[0]->c_total, 'c', 1, data_bag);
		c_total = c[0]->c_total;

	}
	else {

		Communication_update_ghost_nodes_flow_variable(c[0]->data, 'c', 1, data_bag);
		c_total = c[0]->data;
	} // else
	idy_v = grid->idy_v;
	#endif


	for (k=Ks;k<Ke;k++){
		for (j=j_start;j<j_end;j++){
			for (i=Is;i<Ie;i++){

/*				dissbytke =  diss[k][j][i]/tke[k][j][i];
				diss_coeff =  1.0/dtimeb + C_epsilon2*dissbytke;
				rans->two_eqn_rans[1]->data[k][j][i] = rans->two_eqn_rans[1]->ng_rhs[k][j][i]/diss_coeff;

				rans->two_eqn_rans[0]->ng_rhs[k][j][i] = rans->two_eqn_rans[0]->ng_rhs[k][j][i] - diss[k][j][i];
				k_coeff =  1.0/dtimeb -  C_epsilon1*C_mu*strain[k][j][i]*strain[k][j][i]/dissbytke;
				rans->two_eqn_rans[0]->data[k][j][i] = rans->two_eqn_rans[0]->ng_rhs[k][j][i]/k_coeff;
*/

	#ifdef CONC
				cN = 0.5*(c_total[k][j][i] + c_total[k][j+1][i]);
				if (j == NY-2) cN = c_total[k][j][i];


				if (j != 0 ) {
					cS = 0.5*(c_total[k][j][i] + c_total[k][j-1][i]);
				}
				else {
					// Impose Neumann b/c i.e dc/dy = 0.0
					cS = c_total[k][j][i];
				}

				dcdy = (cN - cS) * idy_v[j];
	#endif


				dissbytke =  diss[k][j][i]/tke[k][j][i];
				tkebydiss = tke[k][j][i]/diss[k][j][i];
				strainsq = strain[k][j][i]*strain[k][j][i];

				ak_kcoeff = 1.0/dtimeb - 2.*C_mu*(strainsq + dcdy)*2.*tkebydiss;
				aeps_kcoeff = 2.0;

				bk_epscoeff = -2.*( C_epsilon1*C_mu*(strainsq + C_epsilon3*dcdy) + C_epsilon2*dissbytke*diss[k][j][i]);
				beps_epscoeff = 1.0/dtimeb +  2.*C_epsilon2*2.*dissbytke;

				tke[k][j][i] = (rhs_tke[k][j][i]*beps_epscoeff - rhs_diss[k][j][i]*aeps_kcoeff)/
							(ak_kcoeff*beps_epscoeff - bk_epscoeff*aeps_kcoeff);
				diss[k][j][i] = (rhs_tke[k][j][i]*bk_epscoeff - rhs_diss[k][j][i]*ak_kcoeff)/
							(aeps_kcoeff*bk_epscoeff - beps_epscoeff*ak_kcoeff);

/*				ak_kcoeff = 1.0/dtimeb - C_mu*strain[k][j][i]*strain[k][j][i]/dissbytke;
				aeps_kcoeff = 1.0;

				bk_epscoeff = -C_epsilon1*C_mu*strain[k][j][i]*strain[k][j][i];
				beps_epscoeff = 1.0/dtimeb +  C_epsilon2*dissbytke;
*/

			}
		}
	}

	if (Js == 0) {
		j = Js;
		ywall = grid->yc[j]-grid->yv[0];
		for (k=Ks;k<Ke;k++){
			for (i=Is;i<Ie;i++){
				rans->two_eqn_rans[0]->data[k][j][i] = rans->two_eqn_rans[0]->ng_rhs[k][j][i]*dtimeb;
				turbke = rans->two_eqn_rans[0]->data[k][j][i];
				rans->two_eqn_rans[1]->data[k][j][i] = pow(turbke,1.5)*pow(C_mu,0.75)/(kappa*ywall);
			}
		}
	}

	if (Je == NY) {
		j = NY-2;
		ywall = grid->yv[NY-1]-grid->yc[NY-2];
		for (k=Ks;k<Ke;k++){
			for (i=Is;i<Ie;i++){
				rans->two_eqn_rans[0]->data[k][j][i] = rans->two_eqn_rans[0]->ng_rhs[k][j][i]*dtimeb;
				turbke = rans->two_eqn_rans[0]->data[k][j][i];
				rans->two_eqn_rans[1]->data[k][j][i] = pow(turbke,1.5)*pow(C_mu,0.75)/(kappa*ywall);
			}
		}
	}


#endif  // RANS_SOURCE_IMPLICIT

	double tol;
	tol = 1e-6;
	tol = 0.0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				if (diss[k][j][i] < tol) {
//				Set dissipation to a small value
					diss[k][j][i] = 1e-6;
					if (tke[k][j][i] < tol)
//				Adjust tke such that nut_new = 1.1*nut_old
						tke[k][j][i] = 1.1*pow(nut[k][j][i]*diss[k][j][i]/C_mu, 0.5);
					else
//				Adjust diss such that nut_new = 1.1*nut_old
						diss[k][j][i] = C_mu*tke[k][j][i]*tke[k][j][i]/(1.1*nut[k][j][i]);
				}
				if (tke[k][j][i] < 0)
					tke[k][j][i] = 1.1*pow(nut[k][j][i]*diss[k][j][i]/C_mu, 0.5);

			}
		}
	}

#ifdef IMMERSED_BOUNDARY
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				if (grid->c_status[k][j][i] == SOLID) {
					diss[k][j][i] = 1e-4;
					tke[k][j][i]  = 1e-5;

				}
			}
		}
	}

#endif
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				diss[k][j][i] = min(diss[k][j][i], 1e-1);

			}
		}
	}

	iters = 0;
	return iters;
}




/******************************************************************************/
/*
 This function sets value of tke and diss at the boundaries
 */
/******************************************************************************/
void Two_equation_set_boundary_values(Rans *rans, MAC_grid *grid,
				Parameters *params) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double ***tke, ***diss;

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

	tke = rans->two_eqn_rans[0]->data;
	diss = rans->two_eqn_rans[1]->data;

	if (Je==NY) {
		j = NY-1;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {
				tke[k][j][i] = tke[k][j-1][i];
				diss[k][j][i] = diss[k][j-1][i];
			}
		}
	}

#ifndef ZPERIODIC
	if (Ke==NZ) {
		k = NZ-1;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				tke[k][j][i] = tke[k-1][j][i];
				diss[k][j][i] = diss[k-1][j][i];
			}
		}
	}
#endif

#ifndef XPERIODIC
	if (Ie==NX) {
		i = NX-1;
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				tke[k][j][i] = tke[k][j][i-1];
				diss[k][j][i] = diss[k][j][i-1];
			}
		}
	}
#endif


}
/******************************************************************************/
