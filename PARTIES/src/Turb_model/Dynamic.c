#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Array.h"
#include "Cart3d.h"
#include "Communication.h"
#include "Dynamic_smag.h"
#include "Grid.h"
#include "Lesfilter.h"
#include "Memory.h"
#include "MyMath.h"
#include "Subgrid.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/******************************************************************************/
// Dynamic eddy viscosity model," Phys. Fluids A 3, 1760 (1991)
/******************************************************************************/

/******************************************************************************/
/*
 This function allocates enough memory for the dynamic model related
 variables
 */
/******************************************************************************/
void Dynamic_smag_init(Subgrid *smag, MAC_grid *grid, Parameters *params) {

	double one;

	// Allocate Memory for strain rate tensor
	smag->st_rate->s11 = Memory_allocate_flow_variable(grid, params);
	smag->st_rate->s22 = Memory_allocate_flow_variable(grid, params);
	smag->st_rate->s33 = Memory_allocate_flow_variable(grid, params);
	smag->st_rate->s12 = Memory_allocate_flow_variable(grid, params);
	smag->st_rate->s13 = Memory_allocate_flow_variable(grid, params);
	smag->st_rate->s23 = Memory_allocate_flow_variable(grid, params);

	// Allocate Memory for coefficient etc
	smag->ng_Cev = Memory_allocate_noghost_variable(grid, params);
	smag->ng_LM = Memory_allocate_noghost_variable(grid, params);
	smag->ng_MM = Memory_allocate_noghost_variable(grid, params);
#ifdef LES_LAG_AVE
	smag->ILM = Memory_allocate_flow_variable(grid, params);
	smag->IMM = Memory_allocate_flow_variable(grid, params);

	smag->ng_ILMO = Memory_allocate_noghost_variable(grid, params);
	smag->ng_IMMO = Memory_allocate_noghost_variable(grid, params);
#endif

	// Allocate Memory for filtered cell-centered velocity components
	smag->ucf = Memory_allocate_flow_variable(grid, params);
	smag->vcf = Memory_allocate_flow_variable(grid, params);
	smag->wcf = Memory_allocate_flow_variable(grid, params);

	// Allocate Memory for Intermediate variables
	smag->strainf = Memory_allocate_flow_variable(grid, params);
	smag->sijf = Memory_allocate_flow_variable(grid, params);
	smag->strainsijf = Memory_allocate_flow_variable(grid, params);
	smag->uu = Memory_allocate_flow_variable(grid, params);

	// Allocate Memory for "work" variables
	smag->work1 = Memory_allocate_flow_variable(grid, params);
	smag->work2 = Memory_allocate_flow_variable(grid, params);
	smag->work3 = Memory_allocate_flow_variable(grid, params);
	smag->uk = Memory_allocate_variable_jik(grid, params);

	smag->filter_cpu_time = 0.0;
	smag->filter_comm1_cpu_time = 0.0;
	smag->filter_comm2_cpu_time = 0.0;
	smag->filter_comm3_cpu_time = 0.0;
	smag->filter_comm4_cpu_time = 0.0;
	smag->filter_x_cpu_time = 0.0;
	smag->filter_y_cpu_time = 0.0;
	smag->filter_z_cpu_time = 0.0;

	smag->dynamic_cpu_time = 0.0;
	smag->dynamic_1_cpu_time = 0.0;
	smag->dynamic_2_cpu_time = 0.0;
	smag->dynamic_3_cpu_time = 0.0;
	smag->dynamic_4_cpu_time = 0.0;
	smag->dynamic_5_cpu_time = 0.0;
	smag->dynamic_6_cpu_time = 0.0;
	smag->dynamic_7_cpu_time = 0.0;
	smag->dynamic_8_cpu_time = 0.0;
	smag->dynamic_9_cpu_time = 0.0;
	smag->dynamic_10_cpu_time = 0.0;
	smag->dynamic_11_cpu_time = 0.0;
	smag->dynamic_12_cpu_time = 0.0;
	smag->dynamic_13_cpu_time = 0.0;
	smag->dynamic_14_cpu_time = 0.0;
	smag->dynamic_15_cpu_time = 0.0;

#ifdef CONC
	int iconc, NConc;
	ConcDev **cdev;

	// Allocate Memory concentration related dynamic eddy viscosity variables
	NConc = params->NConc;
	cdev = (ConcDev **)malloc(NConc * sizeof(ConcDev *));
	Memory_check_allocation(cdev);

	for (iconc=0; iconc<NConc; iconc++) {
		cdev[iconc] = (ConcDev *)malloc(sizeof(ConcDev));
		Memory_check_allocation(cdev[iconc]);
		cdev[iconc]->mSct = Memory_allocate_flow_variable(grid, params);
		cdev[iconc]->Sct  = Memory_allocate_flow_variable(grid, params);
		cdev[iconc]->ng_TT  = Memory_allocate_noghost_variable(grid, params);
		cdev[iconc]->ng_KT  = Memory_allocate_noghost_variable(grid, params);
		cdev[iconc]->ccf  = Memory_allocate_flow_variable(grid, params);
	#ifdef LES_LAG_AVE
		cdev[iconc]->ITT = Memory_allocate_flow_variable(grid, params);
		cdev[iconc]->IKT = Memory_allocate_flow_variable(grid, params);

		cdev[iconc]->ng_ITTO = Memory_allocate_noghost_variable(grid, params);
		cdev[iconc]->ng_IKTO = Memory_allocate_noghost_variable(grid, params);

	#endif

		one = 1.0;
		Array_set_withghost(cdev[iconc]->mSct, one, grid, params);
		Array_set_withghost(cdev[iconc]->Sct,  one, grid, params);
	} // for iconc
	smag->cdev = cdev;

#endif // CONC
	return;
}
/******************************************************************************/
/*
 * 	Dynamic Smagorinsky model Coefficient calculation
 */
/******************************************************************************/
void Dynamic_smag_coeff(Cart3d_bag *data_bag){

	double pos1, neg1, pos2, filratio_sq, zero;
	double ***LM, ***MM, ***Cev, ***deltasq, ***strainf;
	double ***ILM, ***IMM,***ILMO,***IMMO;
	double ****TT, ****KT, ****ITT, ****IKT, ****ITTO, ****IKTO;
	double ****Sct;
	double numer, denom;
	double base, eps, T;
	double *g_LM, *g_MM;
	double **g_numer, **g_denom, **g_work;
	int NConc, iconc;

	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	int ierr;
	int i, j, k;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
	Concentration **c = data_bag -> c;
	Subgrid *smag = data_bag -> smag;

	double dt = params -> dt;
/*
	// Generate the name of the file based on the current time step
	sprintf(bin_filename, "ilm_P%d.dat", params->rank);
	fid = fopen(bin_filename,"w");
	sprintf(bin_filename, "imm_P%d.dat", params->rank);
	fid1 = fopen(bin_filename,"w");
*/
	int NX, NY, NZ;
	double Tstart, Tend;

	Tstart = MPI_Wtime();

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


	pos1 = 1.0;
	pos2 = 2.0;
	neg1 = -1.0;
	zero = 0.0;
//	filratio_sq = pow(6.0,0.5);
	filratio_sq = 6.0;
#ifdef FILTER_2D
	filratio_sq = pow(6.0, 2.0/3.0);
#endif

	Array_set_noghost(smag->ng_LM, zero, grid, params);
	Array_set_noghost(smag->ng_MM, zero, grid, params);
	Array_set_noghost(smag->strainf, zero, grid, params);


	// Filter the cell centered velocities and strain rate tensor magnitude
	lesfilter(u->data_bc, smag->ucf, smag->work1, data_bag);
	lesfilter(v->data_bc, smag->vcf, smag->work1, data_bag);
	lesfilter(w->data_bc, smag->wcf, smag->work1, data_bag);

	// Calculate the magnitide of the filtered strain rate
	// \hat{Sxx}
	lesfilter(smag->st_rate->s11, smag->sijf, smag->work1, data_bag);
	Array_XXPY_noghost(smag->ng_LM, smag->sijf, grid, params);
	// \hat{Syy}
	lesfilter(smag->st_rate->s22, smag->sijf, smag->work1, data_bag);
	Array_XXPY_noghost(smag->ng_LM, smag->sijf, grid, params);
	// \hat{Szz}
	lesfilter(smag->st_rate->s33, smag->sijf, smag->work1, data_bag);
	Array_XXPY_noghost(smag->ng_LM, smag->sijf, grid, params);
	// \hat{Sxy}
	lesfilter(smag->st_rate->s12, smag->sijf, smag->work1, data_bag);
	Array_AXXPY_noghost(smag->ng_LM, pos2, smag->sijf, grid, params);
	// \hat{Sxz}
	lesfilter(smag->st_rate->s13, smag->sijf, smag->work1, data_bag);
	Array_AXXPY_noghost(smag->ng_LM, pos2, smag->sijf, grid, params);
	// \hat{Syz}
	lesfilter(smag->st_rate->s23, smag->sijf, smag->work1, data_bag);
	Array_AXXPY_noghost(smag->ng_LM, pos2, smag->sijf, grid, params);

	LM = smag->ng_LM;
	strainf = smag->strainf;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++){
				strainf[k][j][i] = sqrt(2.*LM[k][j][i]);
				LM[k][j][i] = 0.0;
			}
		}
	}


	/*------------------------------------------------------------------------*/
	/*
	 Component tau_xx, Lij = \hat{uu} - \hat{u}\hat{u}
	 Mij = \hat{|S|}\hat{S_xx} - \hat{|S|S_xx}
	 uu = \hat{uu}
	 */
	/*------------------------------------------------------------------------*/
	Array_pointwisemult_noghost(smag->work2, u->data_bc, u->data_bc, grid, params);
	lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

	// strainsijf = \hat{|S|Sxx}
	Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->st_rate->s11, grid, params);
	lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

	// Sxx = \hat{Sxx}
	lesfilter(smag->st_rate->s11, smag->sijf, smag->work1, data_bag);

	// Accumulate into LM and MM
	Array_AXXPY_noghost(smag->uu, neg1, smag->ucf, grid, params);					//Lij = G_uu = \hat{u_c u_c} - \hat{u_c}\hat{u_c}
	Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);            //work3 = \hat{|S|}\hat{S_xx}
	Array_AXPBY_noghost(smag->work3, neg1, filratio_sq, smag->strainsijf, grid, params);            //Mij = G_work3 = filratio_sq*work3 - \hat{|S|S_xx} = filratio*work3 - work2
	Array_XZPY_noghost(smag->ng_LM, smag->uu, smag->work3, grid, params);					//LM = G_LM + Lij*Mij = ng_LM + smag->uu*samg->work3
	Array_XXPY_noghost(smag->ng_MM, smag->work3, grid, params);					//MM  = Mij*Mij = ng_MM + smag->work3 * smag->G_work3


	/*------------------------------------------------------------------------*/
	/*
	 Component tau_yy, Lij = \hat{vv} - \hat{v}\hat{v}
	 Mij = \hat{|S|}\hat{S_yy} - \hat{|S|S_yy}
	 uu = \hat{vv}
	 */
	/*------------------------------------------------------------------------*/
	Array_pointwisemult_noghost(smag->work2, v->data_bc, v->data_bc, grid, params);
	lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

	// strainsijf = \hat{|S|Syy}
	Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->st_rate->s22, grid, params);
	lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

	// Syy = \hat{Syy}
	lesfilter(smag->st_rate->s22, smag->sijf, smag->work1, data_bag);

	// Accumulate into LM and MM
	Array_AXXPY_noghost(smag->uu, neg1, smag->vcf, grid, params);					//Lij = G_uu = \hat{v_c v_c} - \hat{v_c}\hat{v_c}
	Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);            //work3 = \hat{|S|}\hat{S_yy}
	Array_AXPBY_noghost(smag->work3, neg1, filratio_sq, smag->strainsijf, grid, params);            //Mij = G_work3 = filratio_sq*work3 - \hat{|S|S_yy} = filratio*work3 - work2
	Array_XZPY_noghost(smag->ng_LM, smag->uu, smag->work3, grid, params);					//LM = G_LM + Lij*Mij = ng_LM + smag->uu*samg->work3
	Array_XXPY_noghost(smag->ng_MM, smag->work3, grid, params);					//MM  = Mij*Mij = ng_MM + smag->work3 * smag->G_work3


	/*------------------------------------------------------------------------*/
	/*
	 Component tau_zz, Lij = \hat{ww} - \hat{w}\hat{w}
	 Mij = \hat{|S|}\hat{S_zz} - \hat{|S|S_zz}
	 uu = \hat{ww}
	 */
	/*------------------------------------------------------------------------*/
	Array_pointwisemult_noghost(smag->work2, w->data_bc, w->data_bc, grid, params);
	lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

	// strainsijf = \hat{|S|Szz}
	Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->st_rate->s33, grid, params);
	lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

	// Szz = \hat{Szz}
	lesfilter(smag->st_rate->s33, smag->sijf, smag->work1, data_bag);

	// Accumulate into LM and MM
	Array_AXXPY_noghost(smag->uu, neg1, smag->wcf, grid, params);					//Lij = G_uu = \hat{w_c w_c} - \hat{w_c}\hat{w_c}
	Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);            //work3 = \hat{|S|}\hat{S_zz}
	Array_AXPBY_noghost(smag->work3, neg1, filratio_sq, smag->strainsijf, grid, params);            //Mij = G_work3 = filratio_sq*work3 - \hat{|S|S_zz} = filratio*work3 - work2
	Array_XZPY_noghost(smag->ng_LM, smag->uu, smag->work3, grid, params);					//LM = G_LM + Lij*Mij = ng_LM + smag->uu*samg->work3
	Array_XXPY_noghost(smag->ng_MM, smag->work3, grid, params);					//MM  = Mij*Mij = ng_MM + smag->work3 * smag->G_work3


	/*------------------------------------------------------------------------*/
	/*
	 Component tau_xy, Lij = \hat{uv} - \hat{u}\hat{v}
	 Mij = \hat{|S|}\hat{S_xy} - \hat{|S|S_xy}
	 uu = \hat{uv}
	 */
	/*------------------------------------------------------------------------*/
	Array_pointwisemult_noghost(smag->work2, u->data_bc, v->data_bc, grid, params);
	lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

	// strainsijf = \hat{|S|Sxy}
	Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->st_rate->s12, grid, params);
	lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

	// Sxy = \hat{Sxy}
	lesfilter(smag->st_rate->s12, smag->sijf, smag->work1, data_bag);


	// Accumulate into LM and MM
	Array_AXZPY_noghost(smag->uu, neg1, smag->ucf, smag->vcf, grid, params);  			//Lij = uu = \hat{u_c v_c} - \hat{u_c}\hat{v_c}
	Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);            //work3 = \hat{|S|}\hat{S_xy}
	Array_AXPBY_noghost(smag->work3, neg1, filratio_sq, smag->strainsijf, grid, params);            //Mij = G_work3 = filratio_sq*work3 - \hat{|S|S_xy} = filratio*work3 - work2
	Array_AXZPY_noghost(smag->ng_LM, pos2, smag->uu, smag->work3, grid, params);			//LM = G_LM + 2.*Lij*Mij = ng_LM + 2.*smag->uu*samg->work3
	Array_AXXPY_noghost(smag->ng_MM, pos2, smag->work3, grid, params);					//MM  = Mij*Mij = ng_MM + 2.*smag->work3 * smag->G_work3


	/*------------------------------------------------------------------------*/
	/*
	 Component tau_xz, Lij = \hat{uw} - \hat{u}\hat{w}
	 Mij = \hat{|S|}\hat{S_xz} - \hat{|S|S_xz}
	 uu = \hat{uw}
	 */
	/*------------------------------------------------------------------------*/
	Array_pointwisemult_noghost(smag->work2, u->data_bc, w->data_bc, grid, params);
	lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

	// strainsijf = \hat{|S|Sxz}
	Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->st_rate->s13, grid, params);
	lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

	// Sxz = \hat{Sxz}
	lesfilter(smag->st_rate->s13, smag->sijf, smag->work1, data_bag);

	// Accumulate into LM and MM
	Array_AXZPY_noghost(smag->uu, neg1, smag->ucf, smag->wcf, grid, params);  			//Lij = uu = \hat{u_c w_c} - \hat{u_c}\hat{w_c}
	Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);            //work3 = \hat{|S|}\hat{S_xz}
	Array_AXPBY_noghost(smag->work3, neg1, filratio_sq, smag->strainsijf, grid, params);            //Mij = G_work3 = filratio_sq*work3 - \hat{|S|S_xz} = filratio*work3 - work2
	Array_AXZPY_noghost(smag->ng_LM, pos2, smag->uu, smag->work3, grid, params);			//LM = G_LM + 2.*Lij*Mij = ng_LM + 2.*smag->uu*samg->work3
	Array_AXXPY_noghost(smag->ng_MM, pos2, smag->work3, grid, params);					//MM  = Mij*Mij = ng_MM + 2.*smag->work3 * smag->G_work3


	/*------------------------------------------------------------------------*/
	/*
	 Component tau_yz, Lij = \hat{vw} - \hat{v}\hat{w}
	 Mij = \hat{|S|}\hat{S_yz} - \hat{|S|S_yz}
	 uu = \hat{vw}
	 */
	/*------------------------------------------------------------------------*/
	Array_pointwisemult_noghost(smag->work2, v->data_bc, w->data_bc, grid, params);
	lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

	// strainsijf = \hat{|S|Syz}
	Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->st_rate->s23, grid, params);
	lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

	// Syz = \hat{Syz}
	lesfilter(smag->st_rate->s23, smag->sijf, smag->work1, data_bag);


	// Accumulate into LM and MM
	Array_AXZPY_noghost(smag->uu, neg1, smag->vcf, smag->wcf, grid, params);  			//Lij = uu = \hat{v_c w_c} - \hat{v_c}\hat{w_c}
	Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);            //work3 = \hat{|S|}\hat{S_yz}
	Array_AXPBY_noghost(smag->work3, neg1, filratio_sq, smag->strainsijf, grid, params);            //Mij = G_work3 = filratio_sq*work3 - \hat{|S|S_xz} = filratio*work3 - work2
	Array_AXZPY_noghost(smag->ng_LM, pos2, smag->uu, smag->work3, grid, params);			//LM = G_LM + 2.*Lij*Mij = ng_LM + 2.*smag->uu*samg->work3
	Array_AXXPY_noghost(smag->ng_MM, pos2, smag->work3, grid, params);					//MM  = Mij*Mij = ng_MM + 2.*smag->work3 * smag->G_work3

	NConc = params->NConc;

#ifdef CONC
	// Filter the cell centered conentration flow field
	for (iconc=0;iconc<NConc;iconc++){

		Array_set_noghost(smag->cdev[iconc]->ng_KT, zero, grid, params);
		Array_set_noghost(smag->cdev[iconc]->ng_TT, zero, grid, params);

		lesfilter(c[iconc]->data, smag->cdev[iconc]->ccf, smag->work1, data_bag);

		/*--------------------------------------------------------------------*/
		/*
		 Component dcdx, Kj = \hat{uc} - \hat{u}\hat{c}
		 Tj = -filratio*\hat{|S|}\hat{dcdx} + \hat{|S|dcdx}
		 uu = \hat{uc}
		/*--------------------------------------------------------------------*/
		Array_pointwisemult_noghost(smag->work2, u->data_bc, c[iconc]->data, grid, params);
		lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

		// strainsijf = \hat{|S|dcdx}
		Subgrid_concentration_derivative(c[iconc], grid, params, smag, 1);		// Caculate dcdx and store it in work3
		Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->work3, grid, params);
//		Array_pointwisemult_noghost(smag->G_work2, smag->G_strain, smag->cdev[iconc]->G_dcdx);
		lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

		// \hat{dcdx}
		lesfilter(smag->work3, smag->sijf, smag->work1, data_bag);

		// Accumulate into TT and KT
		Array_AXZPY_noghost(smag->uu, neg1, smag->ucf, smag->cdev[iconc]->ccf, grid, params); 		//Kj = uu = \hat{u_c c_c} - \hat{u_c}\hat{c_c}
		Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);          //work3 = \hat{|S|}\hat{dcdx}
		Array_AXPBY_noghost(smag->work3, pos1, -filratio_sq, smag->strainsijf, grid, params);         //Tj = G_work3 = -filratio_sq*work3 - \hat{|S|dcdx} = -filratio*work3 + work2
		Array_XZPY_noghost(smag->cdev[iconc]->ng_KT, smag->uu, smag->work3, grid, params);		//KT = ng_KT + Kj*Tj = ng_KT + smag->uu*samg->work3
		Array_XXPY_noghost(smag->cdev[iconc]->ng_TT, smag->work3, grid, params);			//TT  = Tj*Tj = ng_TT + smag->work3 * smag->G_work3
		// Component dcdy, Kj = \hat{vc} - \hat{v}\hat{c}
		// Tj = -filratio*\hat{|S|}\hat{dcdy} + \hat{|S|dcdy}

		// uu = \hat{vc}
		Array_pointwisemult_noghost(smag->work2, v->data_bc, c[iconc]->data, grid, params);
		lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

		// strainsijf = \hat{|S|dcdy}
		Subgrid_concentration_derivative(c[iconc], grid, params, smag, 2);		// Caculate dcdy and store it in work3
		Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->work3, grid, params);
		lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

		// \hat{dcdy}
		lesfilter(smag->work3, smag->sijf, smag->work1, data_bag);

		// Accumulate into TT and KT
		Array_AXZPY_noghost(smag->uu, neg1, smag->vcf, smag->cdev[iconc]->ccf, grid, params); 		//Kj = uu = \hat{v_c c_c} - \hat{v_c}\hat{c_c}
		Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);          //work3 = \hat{|S|}\hat{dcdy}
		Array_AXPBY_noghost(smag->work3, pos1, -filratio_sq, smag->strainsijf, grid, params);         //Tj = G_work3 = -filratio_sq*work3 + \hat{|S|dcdy} = -filratio*work3 + work2
		Array_XZPY_noghost(smag->cdev[iconc]->ng_KT, smag->uu, smag->work3, grid, params);		//KT = ng_KT + Kj*Tj = ng_KT + smag->uu*samg->work3
		Array_XXPY_noghost(smag->cdev[iconc]->ng_TT, smag->work3, grid, params);			//TT  = Tj*Tj = ng_TT + smag->work3 * smag->G_work3


		/*--------------------------------------------------------------------*/
		/*
		 Component dcdz, Kj = \hat{wc} - \hat{w}\hat{c}
		 Tj = -filratio*\hat{|S|}\hat{dcdz} + \hat{|S|dcdz}
		 uu = \hat{wc}
		 */
		/*--------------------------------------------------------------------*/
		Array_pointwisemult_noghost(smag->work2, w->data_bc, c[iconc]->data, grid, params);
		lesfilter(smag->work2, smag->uu, smag->work1, data_bag);

		// strainsijf = \hat{|S|dcdz}
		Subgrid_concentration_derivative(c[iconc], grid, params, smag, 3);		// Caculate dcdz and store it in work3
		Array_pointwisemult_noghost(smag->work2, smag->st_rate->strain, smag->work3, grid, params);
		lesfilter(smag->work2, smag->strainsijf, smag->work1, data_bag);

		// \hat{dcdz}
		lesfilter(smag->work3, smag->sijf, smag->work1, data_bag);

		// Accumulate into TT and KT
		Array_AXZPY_noghost(smag->uu, neg1, smag->wcf, smag->cdev[iconc]->ccf, grid, params); //Kj = uu = \hat{w_c c_c} - \hat{w_c}\hat{c_c}
		Array_pointwisemult_noghost(smag->work3, smag->strainf, smag->sijf, grid, params);    //work3 = \hat{|S|}\hat{dcdz}
		Array_AXPBY_noghost(smag->work3, pos1, -filratio_sq, smag->strainsijf, grid, params); //Tj = G_work3 = -filratio_sq*work3 + \hat{|S|dcdz} = -filratio*work3 + work2
		Array_XZPY_noghost(smag->cdev[iconc]->ng_KT, smag->uu, smag->work3, grid, params);    //KT = ng_KT + Kj*Tj = ng_KT + smag->uu*samg->work3
		Array_XXPY_noghost(smag->cdev[iconc]->ng_TT, smag->work3, grid, params);              //TT  = Tj*Tj = ng_TT + smag->work3 * smag->G_work3

	}  // for iconc
#endif // CONC

//	Array_pointwisemult_noghost(smag->G_LM, smag->G_lengthscale_sq, smag->G_LM, grid, params);
//	Array_pointwisemult_noghost(smag->G_MM, smag->G_lengthscale_sq, smag->G_MM, grid, params);
//	Array_pointwisemult_noghost(smag->G_MM, smag->G_lengthscale_sq, smag->G_MM, grid, params);

	LM = smag->ng_LM;
	MM = smag->ng_MM;
	Cev = smag->ng_Cev;
	deltasq = smag->ng_lengthscale_sq;
#ifdef LES_LAG_AVE
	ILM = smag->ILM;
	IMM = smag->IMM;
	ILMO = smag->ng_ILMO;
	IMMO = smag->ng_IMMO;
#endif


#ifdef CONC
	Sct = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(Sct);
	TT = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(TT);
	KT = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(KT);
	#ifdef LES_LAG_AVE
	ITT = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(ITT);
	IKT = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(IKT);
	ITTO = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(ITTO);
	IKTO = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(IKTO);
	#endif

	for (iconc=0;iconc<NConc;iconc++){
		Sct[iconc]= smag->cdev[iconc]->Sct;
		TT[iconc] = smag->cdev[iconc]->ng_TT;
		KT[iconc] = smag->cdev[iconc]->ng_KT;
	#ifdef LES_LAG_AVE
		ITT[iconc] = smag->cdev[iconc]->ITT;
		IKT[iconc] = smag->cdev[iconc]->IKT;
		ITTO[iconc] = smag->cdev[iconc]->ng_ITTO;
		IKTO[iconc] = smag->cdev[iconc]->ng_IKTO;
	#endif
	}
#endif // CONC

#ifdef LES_PLANE_AVE
	g_LM = (double *) calloc(NY, sizeof(double));
	g_MM = (double *) calloc(NY, sizeof(double));
	Dynamic_smag_global_horizontal_mean(LM, g_LM, grid);
	Dynamic_smag_global_horizontal_mean(MM, g_MM, grid);
	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				Cev[k][j][i] = -0.5*g_LM[j]/(g_MM[j] + 1e-16)/deltasq[k][j][i];
				Cev[k][j][i] = min(Cev[k][j][i],0.05);
				Cev[k][j][i] = max(Cev[k][j][i],0.0);
			}
		}
	}
	#ifdef CONC
	for (iconc=0;iconc<NConc;iconc++){
		Dynamic_smag_global_horizontal_mean(TT, g_LM, grid);
		Dynamic_smag_global_horizontal_mean(KT, g_MM, grid);
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					Sct[iconc][k][j][i] = Cev[k][j][i]*g_LM[j]*deltasq[k][j][i]/(g_MM[j] + 1e-16);
					Sct[iconc][k][j][i] = max(Sct[iconc][k][j][i],0.1);
					Sct[iconc][k][j][i] = min(Sct[iconc][k][j][i],2.0);
				}
			}
		}
	}
	free(g_LM);
	free(g_MM);
	#endif // CONC
#endif // LES_PLANE_AVE


#ifdef LES_SPAN_AVE
	g_numer = Memory_allocate_2D_double_array(NX, NY);
	g_denom = Memory_allocate_2D_double_array(NX, NY);
	g_work  = Memory_allocate_2D_double_array(NX, NY);

	Dynamic_smag_span_ave(LM, g_numer, g_work, grid, params);
	Dynamic_smag_span_ave(MM, g_denom, g_work, grid, params);

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				Cev[k][j][i] = -0.5*g_numer[j][i]/(g_denom[j][i] + 1e-16)/deltasq[k][j][i];
				Cev[k][j][i] = min(Cev[k][j][i],0.05);
				Cev[k][j][i] = max(Cev[k][j][i],0.0);
			}
		}
	}

	#ifdef CONC
	for (iconc=0;iconc<NConc;iconc++){
		Dynamic_smag_span_ave(TT[iconc], g_numer, g_work, grid, params);
		Dynamic_smag_span_ave(KT[iconc], g_denom, g_work, grid, params);
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					Sct[iconc][k][j][i] = Cev[k][j][i]*g_numer[j][i]*deltasq[k][j][i]/(g_denom[j][i] + 1e-16);
					Sct[iconc][k][j][i] = max(Sct[iconc][k][j][i],0.1);
					Sct[iconc][k][j][i] = min(Sct[iconc][k][j][i],2.0);
				}
			}
		}
	}
	#endif // CONC
	Memory_free_2D_double_array(NY, g_numer);
	Memory_free_2D_double_array(NY, g_denom);
	Memory_free_2D_double_array(NY, g_work);

#endif // LES_SPAN_AVE

#ifdef LES_LAG_AVE
	if (params->iter < 10){

		for (i=i_start; i<i_end; i++){
			for (j=j_start; j<j_end; j++) {

				numer = 0.;
				denom = 0.;

				for (k=k_start; k<k_end; k++) {
					numer = numer + LM[k][j][i];
					denom = denom + MM[k][j][i];
				}

				numer = numer/(k_end - k_start + 1. );
				denom = denom/(k_end - k_start + 1. );
				numer = min(numer, 0.);

				for (k=k_start; k<k_end; k++) {
					ILM[k][j][i] = numer;
					IMM[k][j][i] = denom;
					Cev[k][j][i] = -0.5*numer/(denom + 1e-20)/deltasq[k][j][i];
					Cev[k][j][i] = -0.5*g_LM[j]/(g_MM[j] + 1e-16)/deltasq[k][j][i];
					Cev[k][j][i] = -g_LM[j]/(g_MM[j] + 1e-16)/deltasq[k][j][i];
					Cev[k][j][i] = -0.5*g_LM[j]/(g_MM[j] + 1e-16);
					Cev[k][j][i] = min(Cev[k][j][i],0.05);
					Cev[k][j][i] = max(Cev[k][j][i],0.0);
				}
	#ifdef CONC
				for (iconc=0;iconc<NConc;iconc++){
					numer = 0.;
					denom = 0.;
					for (k=k_start; k<k_end; k++) {
						numer = numer + TT[iconc][k][j][i];
						denom = denom + KT[iconc][k][j][i];
					}
					numer = numer/(k_end - k_start + 1. );
					denom = denom/(k_end - k_start + 1. );
					numer = min(numer, 0.);
					for (k=k_start; k<k_end; k++) {
						ITT[iconc][k][j][i] = numer;
						IKT[iconc][k][j][i] = denom;
						Sct[iconc][k][j][i] = Cev[k][j][i]*ITT[iconc][k][j][i]*deltasq[k][j][i]/(IKT[iconc][k][j][i] + 1e-16);
						Sct[iconc][k][j][i] = max(Sct[iconc][k][j][i],0.1);
						Sct[iconc][k][j][i] = min(Sct[iconc][k][j][i],2.0);
/*
						fprintf(fid,"%3d %3d %3d %+-20.16e %+-20.16e %+-20.16e %+-20.16e \n",
								i,j,k,ILM[k][j][i],IMM[k][j][i],deltasq[k][j][i],Cev[k][j][i]);
						fprintf(fid1,"%3d %3d %3d %+-20.16e %+-20.16e %+-20.16e %+-20.16e \n",
								i,j,k,IKT[0][k][j][i],ITT[0][k][j][i],deltasq[k][j][i],Sct[0][k][j][i]);
*/
					}
				}
	#endif // CONC
			}  // for j
		}  // for i
	}  // if params->iter
	else {
		Dynamic_smag_interpolate_xminusdt(u, v, w, smag, grid, params);
		for (k=k_start; k<k_end; k++){
			for (j=j_start; j<j_end; j++) {
				for (i=i_start; i<i_end; i++){
					base = -8.*ILM[k][j][i]*IMM[k][j][i];
					if (base <= 0.) {
						eps=0.;
					}
					else {
						T = 1.5*sqrt(deltasq[k][j][i])/pow(base,0.125);
						eps = dt/(T+dt);
					}
					ILM[k][j][i] = eps*LM[k][j][i] + (1.-eps)*ILMO[k][j][i];
					ILM[k][j][i] = min(ILM[k][j][i], 0.0);
					IMM[k][j][i] = eps*MM[k][j][i] + (1.-eps)*IMMO[k][j][i];
					Cev[k][j][i] = -0.5*ILM[k][j][i]/(IMM[k][j][i] + 1e-16)/deltasq[k][j][i];
					Cev[k][j][i] = -0.5*g_LM[j]/(g_MM[j] + 1e-16)/deltasq[k][j][i];
					Cev[k][j][i] = -0.5*g_LM[j]/(g_MM[j] + 1e-16);
					Cev[k][j][i] = min(Cev[k][j][i],0.05);
					Cev[k][j][i] = max(Cev[k][j][i],0.0);
	#ifdef CONC
					for (iconc=0;iconc<NConc;iconc++){
						ITT[iconc][k][j][i] = eps*TT[iconc][k][j][i] + (1.-eps)*ITTO[iconc][k][j][i];
						IKT[iconc][k][j][i] = eps*KT[iconc][k][j][i] + (1.-eps)*IKTO[iconc][k][j][i];
						Sct[iconc][k][j][i] = Cev[k][j][i]*ITT[iconc][k][j][i]*deltasq[k][j][i]/(IKT[iconc][k][j][i] + 1e-16);
						Sct[iconc][k][j][i] = max(Sct[iconc][k][j][i],0.1);
						Sct[iconc][k][j][i] = min(Sct[iconc][k][j][i],2.0);

					}
	#endif
/*
					fprintf(fid,"%3d %3d %3d %+-20.16e %+-20.16e %+-20.16e %+-20.16e %+-20.16e %+-20.16e %+-20.16e %+-20.16e \n",
							i,j,k,LM[k][j][i],MM[k][j][i],ILM[k][j][i],IMM[k][j][i],eps,Cev[k][j][i],ILMO[k][j][i],IMMO[k][j][i]);
					fprintf(fid1,"%3d %3d %3d %+-20.16e %+-20.16e %+-20.16e %+-20.16e \n",
							i,j,k,ILM[k][j][i],IMM[k][j][i],deltasq[k][j][i],Cev[k][j][i]);
*/
/*
					fprintf(fid,"%3d %3d %3d %+-20.16e %+-20.16e %+-20.16e %+-20.16e \n",
							i,j,k,ILM[k][j][i],IMM[k][j][i],deltasq[k][j][i],Cev[k][j][i]);
					fprintf(fid1,"%3d %3d %3d %+-20.16e %+-20.16e %+-20.16e %+-20.16e \n",
							i,j,k,IKT[0][k][j][i],ITT[0][k][j][i],deltasq[k][j][i],Sct[0][k][j][i]);
*/

				}  // for i
			}  // for j
		}  // for k
/*
		fclose(fid);
		fclose(fid1);
		exit(0);
*/

	}  // else
#endif  // LES_LAG_AVE

#ifdef CONC
	free(Sct);
	free(TT);
	free(KT);
	#ifdef LES_LAG_AVE
	free(ITT);
	free(IKT);
	free(ITTO);
	free(IKTO);
	#endif
#endif // CONC

#ifdef LES_LAG_AVE
	Communication_update_ghost_nodes_flow_variable(smag->ILM, grid, params, 'c', 1);
	Communication_update_ghost_nodes_flow_variable(smag->IMM, grid, params, 'c', 1);
	#ifdef CONC
	for (iconc=0;iconc<NConc;iconc++){
		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->ITT, grid, params, 'c', 1);
		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->IKT, grid, params, 'c', 1);
		Communication_update_ghost_nodes_flow_variable(smag->cdev[iconc]->Sct, grid, params, 'c', 1);
	}
	#endif
#endif // LES_LAG_AVE



/*
	for (i=Is; i<Ie; i++){
		for (j=Js; j<Je; j++) {
			numer = 0.;
			denom = 0.;
			for (k=Ks; k<Ke; k++) {
				for (k=1; k<Ke-1; k++) {
					numer = numer + LM[k][j][i];
					denom = denom + MM[k][j][i];
				}
				numer = numer/(Ke - 1. );
				denom = denom/(Ke - 1. );
				for (k=Ks; k<Ke; k++) {
					Cev[k][j][i] = -0.5*numer/(denom + 1e-20)/deltasq[k][j][i];
				}
			}
		}
	}
*/

	Tend = MPI_Wtime();
	data_bag->timer->Wtime_sgs_dynamic += Tend - Tstart;

	return;
}




/******************************************************************************/
/*
 This function evaluates the quantities needed in coefficient calculation
 (ILM, IMM etc) at the previous particle position. (i.e at x-udt)
*/
/******************************************************************************/
void Dynamic_smag_interpolate_xminusdt(Velocity *u, Velocity *v, Velocity *w, Subgrid *smag,
		MAC_grid *grid, Parameters *params){

	int NX, NY, NZ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;
	int ierr;
	int i, j, k;
	int NConc, iconc;
	double ***ILM, ***IMM, ***uc, ***vc, ***wc;
	double ***ILMO, ***IMMO;
	double ****ITT, ****IKT, ****ITTO, ****IKTO;
	double *xc, *yc, *zc;
	double *xu, *yv, *zw;

	double xp, yp, zp;
	double dx, dz;
	int i0, j0, k0, i1, j1, k1;
	double csi, eta, zeta, csim, etam, zetam;
	double csimetamzetam, csietamzetam, csimetazetam, csietazetam;
	double csimetamzeta, csietamzeta, csimetazeta, csietazeta;

	double dt = params -> dt;

//	FILE *fid, *fid1;
//	char bin_filename[50];

/*
	// Generate the name of the file based on the current time step
	sprintf(bin_filename, "interpolate_P%d.dat", params->rank);
	fid = fopen(bin_filename,"w");
*/

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
//	i_start = max(1,Is); // i=0, i=NX-1 are not included
//	j_start = max(1,Js);
//	k_start = max(1,Ks);
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
	dx = params->Lx/(grid->NX-1);
	dz = params->Lz/(grid->NZ-1);

	ILM = smag->ILM;
	IMM = smag->IMM;
	ILMO = smag->ng_ILMO;
	IMMO = smag->ng_IMMO;
	uc = u->data_bc;
	vc = v->data_bc;
	wc = w->data_bc;

#ifdef CONC
	NConc = params->NConc;
	ITT = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(ITT);
	IKT = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(IKT);
	ITTO = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(ITTO);
	IKTO = (double ****)malloc(NConc * sizeof(double ***));
	Memory_check_allocation(IKTO);

	for (iconc=0;iconc<NConc;iconc++){
		ITT[iconc] = smag->cdev[iconc]->ITT;
		IKT[iconc] = smag->cdev[iconc]->IKT;
		ITTO[iconc] = smag->cdev[iconc]->ng_ITTO;
		IKTO[iconc] = smag->cdev[iconc]->ng_IKTO;
	}
#endif


	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

				// Find the particle location at the previous timestep
				xp = xc[i]- dt*uc[k][j][i];
				yp = yc[j]- dt*vc[k][j][i];
				zp = zc[k]- dt*wc[k][j][i];

#ifndef XPERIODIC
				if (xp < xc[0])
					xp = xc[0];
				if (xp > xc[NX-2])
					xp = xc[NX-2];
				i0 = Grid_get_x_index(xp, grid, 'c');
				i1 = i0 + 1;
				csi  = (xp - xc[i0])/(xc[i0+1]-xc[i0]);
#else
				if (xp<xc[0]) {
					i0 = -1;
					i1 = i0 + 1;
					csim  = (xc[i1]-xp)/dx;
					csi   = 1-csim;
				}
				else if (xp > xc[NX-2]) {
					i0 = NX-2;
					i1 = i0 + 1;
					csi  = (xp - xc[i0])/dx;
				}
				else {
					i0 = Grid_get_x_index(xp, grid, 'c');
					i1 = i0 + 1;
					csi  = (xp - xc[i0])/(xc[i0+1]-xc[i0]);
				}
#endif  // else notXPERIODIC

				if (yp < yc[0])
					yp = yc[0];
				if (yp > yc[NY-2])
					yp = yc[NY-2];
				j0 = Grid_get_y_index(yp, grid, 'c');
				j1 = j0 + 1;
				eta  = (yp - yc[j0])/(yc[j0+1]-yc[j0]);

#ifndef ZPERIODIC
				if (zp < zc[0])
					zp = zc[0];
				if (zp > zc[NZ-2])
					zp = zc[NZ-2];
				k0 = Grid_get_z_index(zp, grid, 'c');
				k1 = k0 + 1;
				zeta = (zp - zc[k0])/(zc[k0+1]-zc[k0]);
#else
				if (zp < zc[0]) {
					k0 = -1;
					k1 = k0 + 1;
					zetam = (zc[k1]-zp)/dz;
					zeta  = 1-zetam;
				}
				if (zp > zc[NZ-2]) {
					k0 = NZ-2;
					k1 = k0+1;
					zeta = (zp - zc[k0])/dz;
				}
				else {
					k0 = Grid_get_z_index(zp, grid, 'c');
					k1 = k0 + 1;
					zeta = (zp - zc[k0])/(zc[k0+1]-zc[k0]);
				}
#endif  // else notZPERIODIC

				csim = 1.-csi;
				etam = 1.-eta;
				zetam = 1.-zeta;

				csimetamzetam = csim*etam*zetam;
				csietamzetam  = csi* etam*zetam;

				csimetazetam = csim*eta*zetam;
				csietazetam  =  csi*eta*zetam;

				csimetamzeta = csim*etam*zeta;
				csietamzeta  = csi*etam*zeta;

				csimetazeta = csim*eta*zeta;
				csietazeta  = csi*eta*zeta;

				ILMO[k][j][i] = ILM[k0][j0][i0]*csimetamzetam + ILM[k0][j0][i1]*csietamzetam
				              + ILM[k0][j1][i0]*csimetazetam  + ILM[k0][j1][i1]*csietazetam
				              + ILM[k1][j0][i0]*csimetamzeta  + ILM[k1][j0][i1]*csietamzeta
				              + ILM[k1][j1][i0]*csimetazeta   + ILM[k1][j1][i1]*csietazeta;


				IMMO[k][j][i] = IMM[k0][j0][i0]*csimetamzetam + IMM[k0][j0][i1]*csietamzetam
				              + IMM[k0][j1][i0]*csimetazetam  + IMM[k0][j1][i1]*csietazetam
				              + IMM[k1][j0][i0]*csimetamzeta  + IMM[k1][j0][i1]*csietamzeta
				              + IMM[k1][j1][i0]*csimetazeta   + IMM[k1][j1][i1]*csietazeta;

//				ILMO[k][j][i] = ILM[k][j][i];
//				IMMO[k][j][i] = IMM[k][j][i];
//				fprintf(fid,"%5d %5d %5d  %+-20.16e %+-20.16e %+-20.16e %+-20.16e \n",
//						i,j,k,ILM[k][j][i],IMM[k][j][i],ILMO[k][j][i],IMMO[k][j][i]);
#ifdef CONC
				for (iconc=0;iconc<NConc;iconc++){
					ITTO[iconc][k][j][i] = ITT[iconc][k0][j0][i0]*csimetamzetam + ITT[iconc][k0][j0][i1]*csietamzetam
					                     + ITT[iconc][k0][j1][i0]*csimetazetam  + ITT[iconc][k0][j1][i1]*csietazetam
					                     + ITT[iconc][k1][j0][i0]*csimetamzeta  + ITT[iconc][k1][j0][i1]*csietamzeta
					                     + ITT[iconc][k1][j1][i0]*csimetazeta   + ITT[iconc][k1][j1][i1]*csietazeta;

					IKTO[iconc][k][j][i] = IKT[iconc][k0][j0][i0]*csimetamzetam + IKT[iconc][k0][j0][i1]*csietamzetam
					                     + IKT[iconc][k0][j1][i0]*csimetazetam  + IKT[iconc][k0][j1][i1]*csietazetam
					                     + IKT[iconc][k1][j0][i0]*csimetamzeta  + IKT[iconc][k1][j0][i1]*csietamzeta
					                     + IKT[iconc][k1][j1][i0]*csimetazeta   + IKT[iconc][k1][j1][i1]*csietazeta;
				}
#endif // CONC
			}  // for i
		}  // for j
	}  // for k
//	fclose(fid);

#ifdef CONC
	free(ITT);
	free(IKT);
	free(ITTO);
	free(IKTO);
#endif
	return;

}



/******************************************************************************/
/*
 Compute horizontal average.
 */
/******************************************************************************/
void Dynamic_smag_global_horizontal_mean(double ***quantity, double *global_mean, MAC_grid *grid) {

	double *mean;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double scale;
	int ierr;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;
	mean = (double *) calloc(NY, sizeof(double));

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
				mean[j] += quantity[k][j][i];
			}
		}
	}

	scale = 1.0 / ((NX - 1) * (NZ - 1));

	for(j = Js; j < Je; j++) {
		mean[j] *= scale;
	}

	for (j=0;j<NY;j++) {
		global_mean[j] = 0.0;
	}

	MPI_Allreduce(mean, global_mean, NY, MPI_DOUBLE, MPI_SUM, PCW);
	free(mean);


	return;
}




/******************************************************************************/
/*
 Compute spanwise average.
 */
/******************************************************************************/
void Dynamic_smag_span_ave(double ***quantity, double **global_mean, double **work,
		MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double scale;
	double *array, *global_array;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);

	for (j=0; j<NY; j++) {
		for (i=Is; i<Ie; i++) {
			work[j][i] = 0.0;
			global_mean[j][i] = 0.0;
		}
	}

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
				work[j][i] += quantity[k][j][i];
			}
		}
	}

	scale = 1.0 / (NZ - 1);

	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			work[j][i] *= scale;
		}
	}

	array = (double *) malloc( (Je-Js)*(Ie-Is)*sizeof(double));
	global_array = (double *) malloc( (Je-Js)*(Ie-Is)*sizeof(double));

	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			array[(j-Js)*(Ie-Is) + i-Is] = work[j][i];
			global_array[(j-Js)*(Ie-Is) + i-Is] = 0.0;
		}
	}


	MPI_Allreduce(array, global_array, (Je-Js)*(Ie-Is),
				  MPI_DOUBLE, MPI_SUM, params->commxy);
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			global_mean[j][i] = global_array[(j-Js)*(Ie-Is) + i-Is];
		}
	}

	free(array);
	free(global_array);
	return;
}
