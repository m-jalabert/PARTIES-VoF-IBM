#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"

#include "Communication.h"
#include "Display.h"
#include "Memory.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Statistics2d.h"
#include "Output.h"
#include "hdf5.h"
#include "Particle.h"
#include "Interpolate.h"

extern MPI_Comm comm3d;

/******************************************************************************/
/*
 Compute statistics of
     - velocities (u, v, w),
     - pressure,
     - concentrations (sediment concentration and salinity),
 at a given output iteration and corresponding time.
 See below what the statistics are.
 */
/******************************************************************************/
Statistics2d *Statistics2d_create(MAC_grid *grid, Parameters *params) {

	Statistics2d *new_statistics2d;

	double ***slice_c, ***mean_c, **mean_u, **mean_v, **mean_w, **mean_p;
	double **turbStress_uu, **turbStress_vv, **turbStress_ww;
	double **turbStress_uv, **turbStress_uw, **turbStress_vw;
	double ***turbStressC;
	double **flux_uu, **flux_vv, **flux_ww;
	double **flux_uv, **flux_uw, **flux_vw;
	double ***fluxC;
	double **turbKinE;
	double **kinE, **Urms;
	double ***potEnergyC ;
	double **potEnergyP;
	double **pwork;
	double **viscDiss, **viscTran, **sgsDiss, **sgsTran;
	double **work, **work1, **work2, **work3, **work4, **vf_avg;

	double iRe;
	int NX, NY, NZ;
	int rank, i, j, ierr;
	char filename[50];
	int iconc, NConc;
	int Nbin;


	Nbin = params->Nbin;

	Display_progress(params,"Compute statistics \n");

	NConc = params->NConc;

	// Initialize variables
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	rank  = params->rank;
	new_statistics2d = (Statistics2d *) malloc(sizeof(Statistics2d));

	mean_u = Memory_allocate_2D_double_array(NX, NY);
	mean_v = Memory_allocate_2D_double_array(NX, NY);
	mean_w = Memory_allocate_2D_double_array(NX, NY);
	mean_p = Memory_allocate_2D_double_array(NX, NY);
	mean_c = Memory_allocate_3D_double_array(NX, NY, NConc);
	slice_c = Memory_allocate_3D_double_array(NX, NY, NConc);


	turbStress_uu = Memory_allocate_2D_double_array(NX, NY);
	turbStress_vv = Memory_allocate_2D_double_array(NX, NY);
	turbStress_ww = Memory_allocate_2D_double_array(NX, NY);
	turbStress_uv = Memory_allocate_2D_double_array(NX, NY);
	turbStress_uw = Memory_allocate_2D_double_array(NX, NY);
	turbStress_vw = Memory_allocate_2D_double_array(NX, NY);
	turbStressC = Memory_allocate_3D_double_array(NX, NY, 4*NConc);

	flux_uu = Memory_allocate_2D_double_array(NX, NY);
	flux_vv = Memory_allocate_2D_double_array(NX, NY);
	flux_ww = Memory_allocate_2D_double_array(NX, NY);
	flux_uv = Memory_allocate_2D_double_array(NX, NY);
	flux_uw = Memory_allocate_2D_double_array(NX, NY);
	flux_vw = Memory_allocate_2D_double_array(NX, NY);
	fluxC = Memory_allocate_3D_double_array(NX, NY, 4*NConc);

	kinE = Memory_allocate_2D_double_array(NX, NY);
	Urms = Memory_allocate_2D_double_array(NX, NY);
	turbKinE = Memory_allocate_2D_double_array(NX, NY);
	pwork = Memory_allocate_2D_double_array(NX, NY);

	viscDiss = Memory_allocate_2D_double_array(NX, NY);
	viscTran = Memory_allocate_2D_double_array(NX, NY);
	sgsDiss = Memory_allocate_2D_double_array(NX, NY);
	sgsTran = Memory_allocate_2D_double_array(NX, NY);
	potEnergyC = Memory_allocate_3D_double_array(NX, NY, NConc);
	potEnergyP = Memory_allocate_2D_double_array(NX, NY);

	work  = Memory_allocate_2D_double_array(NX, NY);
	work1 = Memory_allocate_2D_double_array(NX, NY);
	work2 = Memory_allocate_2D_double_array(NX, NY);
	work3 = Memory_allocate_2D_double_array(NX, NY);
	work4 = Memory_allocate_2D_double_array(NX, NY);

	vf_avg = Memory_allocate_2D_double_array(NX, NY);


	new_statistics2d->mean_c = mean_c;
	new_statistics2d->slice_c = slice_c;
	new_statistics2d->mean_u = mean_u;
	new_statistics2d->mean_v = mean_v;
	new_statistics2d->mean_w = mean_w;
	new_statistics2d->mean_p = mean_p;

	new_statistics2d->mean_nut = Memory_allocate_2D_double_array(NX, NY);
	new_statistics2d->mean_sct = Memory_allocate_2D_double_array(NX, NY);
	new_statistics2d->mean_cev = Memory_allocate_2D_double_array(NX, NY);
if(params->turbulent_terms_switch){
	new_statistics2d->turbStress_uu = turbStress_uu;
	new_statistics2d->turbStress_vv = turbStress_vv;
	new_statistics2d->turbStress_ww = turbStress_ww;
	new_statistics2d->turbStress_uv = turbStress_uv;
	new_statistics2d->turbStress_uw = turbStress_uw;
	new_statistics2d->turbStress_vw = turbStress_vw;
	new_statistics2d->turbStressC   = turbStressC;
}
	new_statistics2d->flux_uu = flux_uu;
	new_statistics2d->flux_vv = flux_vv;
	new_statistics2d->flux_ww = flux_ww;
	new_statistics2d->flux_uv = flux_uv;
	new_statistics2d->flux_uw = flux_uw;
	new_statistics2d->flux_vw = flux_vw;
	new_statistics2d->fluxC = fluxC;

	new_statistics2d->turbKinE = turbKinE;
	new_statistics2d->kinE = kinE;
	new_statistics2d->Urms = Urms;
	new_statistics2d->potEnergyC = potEnergyC;
	new_statistics2d->potEnergyP = potEnergyP;
	new_statistics2d->pwork = pwork;

	new_statistics2d->viscDiss = viscDiss;
	new_statistics2d->viscTran = viscTran;
	new_statistics2d->sgsDiss = sgsDiss;
	new_statistics2d->sgsTran = sgsTran;

	new_statistics2d->work  = work;
	new_statistics2d->work1 = work1;
	new_statistics2d->work2 = work2;
	new_statistics2d->work3 = work3;
	new_statistics2d->work4 = work4;

	new_statistics2d->vf_avg = vf_avg;

	new_statistics2d->histo = Memory_allocate_2D_double_array(Nbin,NConc);
	new_statistics2d->range = Memory_allocate_2D_double_array(Nbin,NConc);

	new_statistics2d->mean_tke = Memory_allocate_2D_double_array(NX, NY);
	new_statistics2d->mean_eps = Memory_allocate_2D_double_array(NX, NY);


	// 1D variables in x direction
	new_statistics2d->current_height  = Memory_allocate_2D_double_array(NX, NConc);
	new_statistics2d->current_height_fluid  = Memory_allocate_2D_double_array(NX, NConc);

	// 0D variables
	new_statistics2d->front_position  = Memory_allocate_1D_array(GVG_DOUBLE, NConc);
	new_statistics2d->integral_potEnergy_fluid  = Memory_allocate_1D_array(GVG_DOUBLE, NConc);
	new_statistics2d->integral_nusselt  = Memory_allocate_1D_array(GVG_DOUBLE, NConc);
	new_statistics2d->integral_conc_fluid  = Memory_allocate_1D_array(GVG_DOUBLE, NConc);
	new_statistics2d->integral_conc_full  = Memory_allocate_1D_array(GVG_DOUBLE, NConc);
	new_statistics2d->integral_buoyantwork = Memory_allocate_1D_array(GVG_DOUBLE, NConc);

	new_statistics2d->grad_y_c = Memory_allocate_3D_double_array(NX, NY, NConc);
	new_statistics2d->nusselt= Memory_allocate_3D_double_array(NX, NY, NConc);

	int Navg;
	if (params->avg_dir == 0){
		Navg = NX-1;
	}
	else{
		Navg = NY-1;
	}



	new_statistics2d->work_avg = Memory_allocate_1D_array(GVG_DOUBLE, Navg);
	new_statistics2d->avg_vf = Memory_allocate_1D_array(GVG_DOUBLE, Navg);
	new_statistics2d->avg_c = Memory_allocate_2D_double_array(Navg,NConc);
	// Average uc
	new_statistics2d->avg_uc = Memory_allocate_2D_double_array(Navg,NConc);
	new_statistics2d->avg_vc = Memory_allocate_2D_double_array(Navg,NConc);
	new_statistics2d->avg_wc = Memory_allocate_2D_double_array(Navg,NConc);

	return new_statistics2d;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Statistics2d_destroy(Statistics2d *st2d, MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ;
	int rank, i, j, ierr;
	int iconc, NConc, Nbin;

	NConc = params->NConc;
	Nbin = params->Nbin;

	// Initialize variables
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;


	// free memory
	Memory_free_2D_double_array(NY, st2d->mean_u);
	Memory_free_2D_double_array(NY, st2d->mean_v);
	Memory_free_2D_double_array(NY, st2d->mean_w);
	Memory_free_2D_double_array(NY, st2d->mean_p);
	Memory_free_3D_double_array(NY, NConc, st2d->mean_c);
	Memory_free_3D_double_array(NY, NConc, st2d->slice_c);
	Memory_free_2D_double_array(NY, st2d->mean_nut);
	Memory_free_2D_double_array(NY, st2d->mean_sct);
	Memory_free_2D_double_array(NY, st2d->mean_cev);
if(params->turbulent_terms_switch){
	Memory_free_2D_double_array(NY, st2d->turbStress_uu);
	Memory_free_2D_double_array(NY, st2d->turbStress_vv);
	Memory_free_2D_double_array(NY, st2d->turbStress_ww);
	Memory_free_2D_double_array(NY, st2d->turbStress_uv);
	Memory_free_2D_double_array(NY, st2d->turbStress_uw);
	Memory_free_2D_double_array(NY, st2d->turbStress_vw);
	Memory_free_3D_double_array(NY, 4*NConc, st2d->turbStressC);
	Memory_free_2D_double_array(NY, st2d->turbKinE);
}

	Memory_free_2D_double_array(NY, st2d->flux_uu);
	Memory_free_2D_double_array(NY, st2d->flux_vv);
	Memory_free_2D_double_array(NY, st2d->flux_ww);
	Memory_free_2D_double_array(NY, st2d->flux_uv);
	Memory_free_2D_double_array(NY, st2d->flux_uw);
	Memory_free_2D_double_array(NY, st2d->flux_vw);
	Memory_free_3D_double_array(NY, 4*NConc, st2d->fluxC);
	Memory_free_2D_double_array(NY, st2d->kinE);
	Memory_free_2D_double_array(NY, st2d->Urms);
	Memory_free_2D_double_array(NY, st2d->pwork);
	Memory_free_2D_double_array(NY, st2d->potEnergyP);

	Memory_free_3D_double_array(NY, NConc, st2d->potEnergyC);
	Memory_free_3D_double_array(NY, NConc, st2d->nusselt);
	Memory_free_3D_double_array(NY, NConc, st2d->grad_y_c);

	Memory_free_2D_double_array(NY, st2d->viscDiss);
	Memory_free_2D_double_array(NY, st2d->viscTran);
	Memory_free_2D_double_array(NY, st2d->sgsDiss);
	Memory_free_2D_double_array(NY, st2d->sgsTran);

	Memory_free_2D_double_array(NY, st2d->work);
	Memory_free_2D_double_array(NY, st2d->work1);
	Memory_free_2D_double_array(NY, st2d->work2);
	Memory_free_2D_double_array(NY, st2d->work3);
	Memory_free_2D_double_array(NY, st2d->work4);

	Memory_free_2D_double_array(NY, st2d->vf_avg);

	Memory_free_2D_double_array(NY, st2d->mean_tke);
	Memory_free_2D_double_array(NY, st2d->mean_eps);


	// 1D variables in x direction
	Memory_free_2D_double_array(NX, st2d->current_height );
	Memory_free_2D_double_array(NX, st2d->current_height_fluid);


	Memory_free_2D_double_array(NConc, st2d->histo);
	Memory_free_2D_double_array(NConc, st2d->range);

	Memory_free_2D_double_array(NConc, st2d->avg_uc);
	Memory_free_2D_double_array(NConc, st2d->avg_vc);
	Memory_free_2D_double_array(NConc, st2d->avg_wc);
	Memory_free_2D_double_array(NConc, st2d->avg_c);

		// 0D variables
	free(st2d->front_position);
	free(st2d->integral_potEnergy_fluid);
	free(st2d->integral_nusselt);
	free(st2d->integral_conc_fluid);
	free(st2d->integral_conc_full);
	free(st2d->integral_buoyantwork);
	free(st2d->avg_vf);
	free(st2d->work_avg);

	free(st2d);

	return;
}

/******************************************************************************/
/*
 */
/******************************************************************************/
void Statistics2d_computeStatistics(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	// FLow variables
	double T1, T2;
	T1 = MPI_Wtime();

	double ***u_cell, ***v_cell, ***w_cell;
	double ***u_face, ***v_face, ***w_face;
	double ***u_forc, ***v_forc, ***w_forc;
	double ****c, ***p;
	double ***nut;
	double tmp;
	// Statistics
	double ***slice_c, ***mean_c, **mean_u, **mean_v, **mean_w, **mean_p;
	double **turbStress_uu, **turbStress_vv, **turbStress_ww;
	double **turbStress_uv, **turbStress_uw, **turbStress_vw;
	double ***turbStressC;
	double **flux_uu, **flux_vv, **flux_ww;
	double **flux_uv, **flux_uw, **flux_vw;
	double ***fluxC;
	double **turbKinE;
	double ***potEnergyC;
	double **kinE, **Urms;
	double **pwork;
	double **potEnergyP ;
	double **viscDiss, **viscTran, **sgsDiss, **sgsTran;
	double **work, **work1, **work2, **work3, **work4, **vf_avg;
	double **mean_tke, **mean_eps;
	double kinEp;
	double potEp;

	Wall_model *log_law;

	double *y;
	double iRe;

	int NX, NY, NZ;
	int rank, i, j;
	char filename[50];
	int iconc, NConc;

	double **mean_nut;
	double tetime;
	FILE *fid;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Display_progress(params,"Compute statistics \n");

	Statistics2d *st2d = data_bag -> st2d;
	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;
	Pressure *pressure = data_bag -> p;
	Concentration **conc = data_bag -> c;
	Subgrid *smag = data_bag -> smag;
	Rans *rans = data_bag -> rans;
	Lagrangian *lag = data_bag->lag;

	double *Ri=params->richardson;
	NConc = params->NConc;

	Communication_update_ghost_nodes_flow_variable(pressure->p_data,  'c', 1, data_bag);
#ifdef CONC
	for (iconc=0; iconc<NConc; iconc++) {
		Communication_update_ghost_nodes_flow_variable(conc[iconc]->data,  'c', 1, data_bag);

	}
#endif
	Communication_update_ghost_nodes_flow_variable(u->data_bc,  'c', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data_bc,  'c', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data_bc,  'c', 1, data_bag);

	Communication_update_ghost_nodes_flow_variable(u->data,  'u', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data,  'v', 1, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data,  'w', 1, data_bag);

	// get data arrays
	u_face = u->data;
	v_face = v->data;
	w_face = w->data;

	u_cell = u->data_bc;
	v_cell = v->data_bc;
	w_cell = w->data_bc;

	#ifdef TURB_FORCING
	u_forc = u->fturb;
	v_forc = v->fturb;
	w_forc = w->fturb;
	#endif

	p = pressure->p_data;

#ifdef CONC
	c = (double ****)calloc(NConc,sizeof(double ***));

	for (iconc = 0; iconc < NConc; iconc++) {
		c[iconc]  = conc[iconc]->data;
	}
#endif

#ifdef LES
	Communication_update_ghost_nodes_flow_variable(smag->nut,'c', 1, data_bag);
	nut = smag->nut;
	#ifdef SCHUMANN
	log_law = smag->log_law;
	#endif
#endif

#ifdef RANS
	Communication_update_ghost_nodes_flow_variable(rans->nut,'c', 1, data_bag);
	nut = rans->nut;
	#ifdef SCHUMANN
	log_law = rans->log_law;
	#endif
#endif

	// Calculate cell centered volume fraction
#ifdef LAG_PARTICLE_RESOLVED
	// Include both local and foreign particles
if (params -> post_processing_switch){	// if this is to be used as post processing tool we need volume fractions
	Particle_MPI_update(data_bag->lag->p_mobile_list, data_bag, DTRACE("Particle_MPI_update"));
	Particle_MPI_update(data_bag->lag->p_fixed_list, data_bag, DTRACE("Particle_MPI_update"));

	// Calculate volume fractions
	Memory_reset_noghost_variable(grid, params, data_bag->lag->ng_vfc);
	Interpolate_add_to_volume_fraction('c', data_bag->lag->p_mobile_list,
	data_bag, DTRACE("Interpolate_add_to_volume_fraction"));
	Interpolate_add_to_volume_fraction('c', data_bag->lag->p_fixed_list,
	data_bag, DTRACE("Interpolate_add_to_volume_fraction"));





	Particle_list_remove(data_bag->lag->p_mobile_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));
	Particle_list_remove(data_bag->lag->p_fixed_list, FOREIGN, grid, params, DTRACE("Particle_list_remove"));
}
#endif

	// Initialize variables
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;
	y  = grid->yc;

	rank  = params->rank;

	mean_u = st2d->mean_u;
	mean_v = st2d->mean_v;
	mean_w = st2d->mean_w;
	mean_p = st2d->mean_p;
	mean_c = st2d->mean_c;
	slice_c = st2d->slice_c;

	turbStress_uu = st2d->turbStress_uu;
	turbStress_vv = st2d->turbStress_vv;
	turbStress_ww = st2d->turbStress_ww;
	turbStress_uv = st2d->turbStress_uv;
	turbStress_uw = st2d->turbStress_uw;
	turbStress_vw = st2d->turbStress_vw;
	turbStressC   = st2d->turbStressC;

	flux_uu = st2d->flux_uu;
	flux_vv = st2d->flux_vv;
	flux_ww = st2d->flux_ww;
	flux_uv = st2d->flux_uv;
	flux_uw = st2d->flux_uw;
	flux_vw = st2d->flux_vw;
	fluxC   = st2d->fluxC;

	kinE     = st2d->kinE;
	Urms     = st2d->Urms;
	turbKinE = st2d->turbKinE;
	pwork    = st2d->pwork;
	potEnergyC = st2d->potEnergyC;
	potEnergyP =  st2d->potEnergyP;

	viscDiss = st2d->viscDiss;
	viscTran = st2d->viscTran;
	sgsDiss  = st2d->sgsDiss;
	sgsTran  = st2d->sgsTran;

	work  = st2d->work;
	work1 = st2d->work1;
	work2 = st2d->work2;
	work3 = st2d->work3;
	work4 = st2d->work4;

	vf_avg = st2d->vf_avg; // fluid volume fraction

	// 1D stuff

	double **current_height = st2d ->current_height;
	double **current_height_fluid = st2d ->current_height_fluid;
	double *front_position = st2d ->front_position;

	double **histo = st2d ->histo;
	double **range = st2d ->range;

	double **avg_uc = st2d ->avg_uc;
	double **avg_vc = st2d ->avg_vc;
	double **avg_wc = st2d ->avg_wc;
	double **avg_c = st2d ->avg_c;

	double *work_avg = st2d ->work_avg;
	double *avg_vf = st2d ->avg_vf;


	double Lx=params->Lx;
	double Ly=params->Ly;
	double Lz=params->Lz;


	/*------------------------------------------------------------------------*/
	/*
	 Start computing statistics
	 */
	/*------------------------------------------------------------------------*/



	//--------------------------------------------------------------------------
	// Horizontally averaged means
	//--------------------------------------------------------------------------
	Statistics2d_horizontalMean2d(u_cell, mean_u, 'c', data_bag);
	Statistics2d_horizontalMean2d(v_cell, mean_v, 'c', data_bag);
	Statistics2d_horizontalMean2d(w_cell, mean_w, 'c', data_bag);
	Statistics2d_horizontalMean2d(p, mean_p, 'c', data_bag);


#ifdef LAG_PARTICLE_RESOLVED
	Statistics2d_vf_avg(vf_avg, data_bag);
#endif

#ifdef CONC
	for (iconc = 0; iconc < NConc; iconc++) {
		Statistics2d_horizontalMean2d(c[iconc], mean_c[iconc], 'c', data_bag);
		st2d->integral_conc_full[iconc]= Statistics2d_volume_integral(c[iconc],  FULL_DOMAIN, data_bag)/Lx/Ly/Lz;
		Statistics2d_histocalc(c[iconc], histo[iconc], range[iconc], data_bag);

		Statistics2d_covarcalc(avg_uc[iconc], c[iconc],u_cell, avg_vf, avg_c[iconc], work_avg, data_bag);
		Statistics2d_covarcalc(avg_vc[iconc], c[iconc],v_cell, avg_vf, avg_c[iconc], work_avg, data_bag);
		Statistics2d_covarcalc(avg_wc[iconc], c[iconc],w_cell, avg_vf, avg_c[iconc], work_avg, data_bag);


#ifdef LAG_PARTICLE_RESOLVED
		st2d->integral_conc_fluid[iconc]= Statistics2d_volume_integral(c[iconc], FLUID_DOMAIN, data_bag)/Lx/Ly/Lz;
#endif

		Statistics2d_horizontalSlice2d(c[iconc], slice_c[iconc], 'c', data_bag);
	}
#endif

	//--------------------------------------------------------------------------
	//Horizontally averaged momentum fluxes
	//--------------------------------------------------------------------------
	Statistics2d_fluxes(flux_uu, u_face, u_face, 'u', data_bag);
	Statistics2d_uface_to_cellcenter(flux_uu, flux_uu, grid);

	Statistics2d_fluxes(flux_vv, v_face, v_face, 'v', data_bag);
	Statistics2d_vface_to_cellcenter(flux_vv, flux_vv, grid);

	Statistics2d_fluxes(flux_ww, w_face, w_face, 'w', data_bag);


	Statistics2d_fluxes(flux_uv, u_cell, v_cell, 'c', data_bag);
	Statistics2d_fluxes(flux_uw, u_cell, w_cell, 'c', data_bag);
	Statistics2d_fluxes(flux_vw, v_cell, w_cell, 'c', data_bag);



	//--------------------------------------------------------------------------
	// Turbulent stresses (Reynolds stresses) <u'u'>, <u'v'>, ...
	//--------------------------------------------------------------------------
if (params-> turbulent_terms_switch == 1) {
	Statistics2d_turbulentStress(turbStress_uu, u_face, u_face, 'u', data_bag);
	Statistics2d_uface_to_cellcenter(turbStress_uu, turbStress_uu, grid);

	Statistics2d_turbulentStress(turbStress_vv, v_face, v_face, 'v', data_bag);
	Statistics2d_vface_to_cellcenter(turbStress_vv, turbStress_vv, grid);

	Statistics2d_turbulentStress(turbStress_ww, w_face, w_face, 'w', data_bag);
	Statistics2d_turbulentStress(turbStress_uv, u_cell, v_cell, 'c', data_bag);
	Statistics2d_turbulentStress(turbStress_uw, u_cell, w_cell, 'c', data_bag);
	Statistics2d_turbulentStress(turbStress_vw, v_cell, w_cell, 'c', data_bag);
}
#ifdef CONC
	int indexC=0;

	if(params->nusselt_switch){
		Statistics2d_compute_nusselt(data_bag);  // calculates nusselt number for all conc fields
	//	printf("Nusselt = %g \n ", st2d->integral_nusselt[0]);
	}

	for (iconc = 0; iconc < NConc; iconc++) {
		st2d->integral_buoyantwork[iconc]=0;
		indexC=iconc*4;
		Statistics2d_fluxes(fluxC[indexC+0], u_cell, c[iconc], 'c', data_bag);
		Statistics2d_fluxes(fluxC[indexC+1], v_cell, c[iconc], 'c', data_bag);
		Statistics2d_fluxes(fluxC[indexC+2], w_cell, c[iconc], 'c', data_bag);
		Statistics2d_fluxes(fluxC[indexC+3], c[iconc], c[iconc], 'c', data_bag);
		if(params->nusselt_switch){
			Statistics2d_vface_to_cellcenter(&st2d->nusselt[iconc][0], &st2d->nusselt[iconc][0], grid);
		}

		if (params-> turbulent_terms_switch == 1) {
			Statistics2d_turbulentStress(turbStressC[indexC+0], u_cell, c[iconc], 'c', data_bag);
			Statistics2d_turbulentStress(turbStressC[indexC+1], v_cell, c[iconc], 'c', data_bag);
			Statistics2d_turbulentStress(turbStressC[indexC+2], w_cell, c[iconc], 'c', data_bag);
			Statistics2d_turbulentStress(turbStressC[indexC+3], c[iconc], c[iconc], 'c', data_bag);
		}


		for (j=0; j < NY-1; j++){
				for (i=0; i < NX-1; i++){
#ifdef LAG_PARTICLE_RESOLVED
					st2d->integral_buoyantwork[iconc] += fluxC[indexC+1][j][i]*Ri[iconc]*vf_avg[j][i];
#else
					st2d->integral_buoyantwork[iconc] += fluxC[indexC+1][j][i]*Ri[iconc];
#endif
				}
		}

		st2d->integral_buoyantwork[iconc] *= 1./(NX-1.)/(NY-1.);

	}




#endif





	if (params-> turbulent_terms_switch == 1) {
		for(j = 0; j < NY; j++) {
			for(i = 0; i < NX; i++) {
				turbKinE[j][i] = 0.5*(turbStress_uu[j][i] + turbStress_vv[j][i] + turbStress_ww[j][i]);
			}
		}
	}

	//--------------------------------------------------------------------------
	// computing kinetic Energy
	//--------------------------------------------------------------------------

	st2d->integral_kinEnergy_fluid=Statistics2d_totKineticEnergy(kinE, u_cell, v_cell, w_cell, data_bag);

	//--------------------------------------------------------------------------
	// computing root mean square of velocity
	//--------------------------------------------------------------------------

	st2d->integral_RMS_fluid=Statistics2d_RMS(Urms, u_cell, v_cell, w_cell, data_bag);

	#ifdef TURB_FORCING
	//--------------------------------------------------------------------------
	// computing the work of turbulent forcing
	//--------------------------------------------------------------------------

	st2d->integral_forcing_fluid=Statistics2d_forcing(u_face, v_face, w_face,	u_forc, v_forc, w_forc, data_bag);
	#endif

#ifdef CONC
	//--------------------------------------------------------------------------
	// computing potential Energy
	//--------------------------------------------------------------------------
	for (iconc = 0; iconc < NConc; iconc++) {
		Statistics2d_potEnergy_y(potEnergyC[iconc], c[iconc], params->richardson[iconc] , &(st2d->integral_potEnergy_fluid[iconc]), data_bag);
		Statistics2d_partial_y(st2d->grad_y_c[iconc], c[iconc], data_bag);
	}
#endif


#ifdef LAG_PARTICLE_RESOLVED
	st2d->integral_potEnergy_part = Statistics2d_potE_part(data_bag, DTRACE("Statistics2d_kinE_part"));
	st2d->integral_kinEnergy_part = Statistics2d_kinE_part(data_bag, DTRACE("Statistics2d_kinE_part"));

	Statistics2d_potEnergyPart(potEnergyP, data_bag);



/*	if (params->rank == 0) {
		if (params->time == 0)
			fid = fopen("kinEp.dat", "w");
		else
			fid = fopen("kinEp.dat", "a");
			fprintf(fid, "%+-20.16e %+-20.16e \n", params->time, st2d->integral_kinEnergy_part);
			fclose(fid);
		if (params->time == 0)
			fid = fopen("potEp.dat", "w");
		else
			fid = fopen("potEp.dat", "a");
			fprintf(fid, "%+-20.16e %+-20.16e \n", params->time, st2d->integral_potEnergy_part);
			fclose(fid);
	}*/ // This stuff is now located inside the hdf5 file

	#endif


	Statistics2d_pressureWork(u_face, v_face, w_face, pressure, pwork, data_bag);

	iRe = 1.e0/params->Re;
	//--------------------------------------------------------------------------
	//Compute kinetic energy budget
	//--------------------------------------------------------------------------
	Statistics2d_ke_budget(u_face, v_face, w_face, u_cell, v_cell, w_cell, nut,
			log_law, viscDiss, viscTran, sgsDiss, sgsTran, data_bag, iRe);


//#ifdef CONC
//	Statistics0d_Frontposition(conc[0]->data, &current_height[0][0], &current_height_fluid[0][0] , front_position  ,  data_bag);
//#endif
	//--------------------------------------------------------------------------
	// write statistics to file
	//--------------------------------------------------------------------------


	Statistics2d_writeStatistics2dToH5(st2d, grid, params, data_bag ,params->noutput_2d, DTRACE("writeStatistics2dToH5"));


#ifdef LES
	mean_nut = st2d->mean_nut;
	for(j = 0; j < NY; j++) {
		for(i = 0; i < NX; i++) {
			mean_nut[j][i] = 0.0;
		}
	}
	Statistics2d_horizontalMean2d(nut, mean_nut, 'c', data_bag);

	#ifdef CONC_DYNAMIC
	double ***Sct, **mean_sct;
	double ***Cev, **mean_cev;
	Sct = smag->cdev[0]->Sct;
	Cev = smag->ng_Cev;
	mean_sct = st2d->mean_sct;
	mean_cev = st2d->mean_cev;
	for(j = 0; j < NY; j++) {
		for(i = 0; i < NX; i++) {
			mean_sct[j][i] = 0.0;
			mean_cev[j][i] = 0.0;
		}
	}
	Statistics2d_horizontalMean2d(Sct, mean_sct,'c', data_bag);
	Statistics2d_horizontalMean2d(Cev, mean_cev,'c', data_bag);

	#endif // CONC_DYNAMIC
#endif // LES

#ifdef RANS
	mean_nut = st2d->mean_nut;
	mean_tke = st2d->mean_tke;
	mean_eps = st2d->mean_eps;

	for(j = 0; j < NY; j++) {
		for(i = 0; i < NX; i++) {
			mean_nut[j][i] = 0.0;
			mean_tke[j][i] = 0.0;
			mean_eps[j][i] = 0.0;
		}
	}
	Statistics2d_horizontalMean2d(nut, mean_nut, 'c', data_bag);

	#ifdef TWO_EQUATION_MODEL
	Statistics2d_horizontalMean2d(rans->two_eqn_rans[0]->data, 'c', data_bag);
	Statistics2d_horizontalMean2d(rans->two_eqn_rans[1]->data, mean_eps, 'c', data_bag);
	#endif
#endif

#ifdef CONC
	free(c);
#endif
	T2 = MPI_Wtime();
	data_bag->timer->Wtime_output_2d += T2 - T1;
	return;
}


/******************************************************************************/
/*
 Compute horizontally average conditional on the fluid phase = < A(x,y,z) phi_f(x,y,z) >_z / < phi_f>_Z
 */
/******************************************************************************/
void Statistics2d_horizontalMean2d(double ***quantity, double **global_mean, char which_component, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double ***vf;
	double **vf_avg;
	double vfrac = 0.0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);

	// Use work buffer for avg volume fraction
	vf_avg = data_bag->st2d->work;

	// Set global mean to 0 before MPI_Allreduce is called
	DSET_ZERO(global_mean[0],NX*NY);
	// Set buffer for vf_avg to 0
	DSET_ZERO(vf_avg[0],NX*NY);

	// Choose the right component of vf if particles are resolved
#ifdef LAG_PARTICLE_RESOLVED
	if (which_component == 'c'){
		vf = data_bag->lag->ng_vfc;
	}
	else if (which_component == 'u'){
		vf = data_bag->lag->ng_vfu;
	}
	else if (which_component == 'v'){
		vf = data_bag->lag->ng_vfv;
	}
	else if (which_component == 'w'){
		vf = data_bag->lag->ng_vfw;
	}
#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				global_mean[j][i] += quantity[k][j][i] * (1.0-vfrac);
				vf_avg[j][i] += (1.0-vfrac);
			}
		}
	}
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			global_mean[j][i] = global_mean[j][i] / vf_avg[j][i] / params->NPZ;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &global_mean[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);

	return;
}

/******************************************************************************/
/*
 Compute horizontally average.
 */
/******************************************************************************/
void Statistics2d_horizontalSlice2d(double ***quantity, double **global_slice, char which_component, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double vfrac = 0.0;


	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

int N_slice = round(NZ/2);

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);
	//Ke = grid->G_Ke;

	//#ifdef XPERIODIC
	//	Is = max(1,Is);
	//	if (Ie == NX-1)
	//		Ie = NX;
	//#endif

	//#ifdef ZPERIODIC
	//	Ks = max(1,Ks);
	//	if (Ke == NZ-1)
	//		Ke = NZ;
	//#endif

	DSET_ZERO(global_slice[0],NX*NY);

	if ( ( N_slice >= Ks) && (N_slice < Ke) ){
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				global_slice[j][i] = quantity[N_slice][j][i];
			}
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &global_slice[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);

	return;
}

void Statistics2d_horizontalSlice2d_near_Wall(double ***quantity, double **global_slice, char which_component, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double vfrac = 0.0;


	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

int N_slice = floor(NZ/params->Lz*0.1); // we output the slice at z=0.1

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);
	//Ke = grid->G_Ke;

	//#ifdef XPERIODIC
	//	Is = max(1,Is);
	//	if (Ie == NX-1)
	//		Ie = NX;
	//#endif

	//#ifdef ZPERIODIC
	//	Ks = max(1,Ks);
	//	if (Ke == NZ-1)
	//		Ke = NZ;
	//#endif

	DSET_ZERO(global_slice[0],NX*NY);

	if ( ( N_slice >= Ks) && (N_slice < Ke) ){
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				global_slice[j][i] = quantity[N_slice][j][i];
			}
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &global_slice[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);

	return;
}

/******************************************************************************/
/*
 Compute turbulent stress <q1' q2'>
 */
/******************************************************************************/
void Statistics2d_turbulentStress(double **turbStress, double ***quantity1,
	 double ***quantity2, char which_component, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double ***vf;
	double **vf_avg, **q1_mean, **q2_mean;
	double vfrac = 0.0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);

	// Use work buffer for avg volume fraction
	vf_avg  = data_bag->st2d->work;
	q1_mean = data_bag->st2d->work1;
	q2_mean = data_bag->st2d->work2;

	// Set global mean to 0 before MPI_Allreduce is called
	DSET_ZERO(turbStress[0],NX*NY);
	// Set buffers to 0
	DSET_ZERO(vf_avg[0],NX*NY);
	DSET_ZERO(q1_mean[0],NX*NY);
	DSET_ZERO(q2_mean[0],NX*NY);

	// Choose the right component of vf if particles are resolved
#ifdef LAG_PARTICLE_RESOLVED
	if (which_component == 'c'){
		vf = data_bag->lag->ng_vfc;
	}
	else if (which_component == 'u'){
		vf = data_bag->lag->ng_vfu;
	}
	else if (which_component == 'v'){
		vf = data_bag->lag->ng_vfv;
	}
	else if (which_component == 'w'){
		vf = data_bag->lag->ng_vfw;
	}
#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				q1_mean[j][i] += quantity1[k][j][i] * (1.0-vfrac);
				q2_mean[j][i] += quantity2[k][j][i] * (1.0-vfrac);
				turbStress[j][i] += quantity1[k][j][i] * quantity2[k][j][i]  * (1.0-vfrac) ;
				vf_avg[j][i] += (1.0-vfrac);
			}
		}
	}

	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			turbStress[j][i] = turbStress[j][i] / vf_avg[j][i];
			q1_mean[j][i] = q1_mean[j][i] / vf_avg[j][i];
			q2_mean[j][i] = q2_mean[j][i] / vf_avg[j][i];
			turbStress[j][i] = ( turbStress[j][i] - (q1_mean[j][i] * q2_mean[j][i])) / params->NPZ;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &turbStress[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	return;
}



/******************************************************************************/
/*
 */
/******************************************************************************/

void Statistics2d_fluxes(double **flux, double ***quantity1, double ***quantity2, char which_component, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double ***vf;
	double **vf_avg;
	double vfrac = 0.0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);

	// Use work buffer for avg volume fraction
	vf_avg = data_bag->st2d->work;

	// Set flux to 0 before MPI_Allreduce is called
	DSET_ZERO(flux[0],NX*NY);
	DSET_ZERO(vf_avg[0],NX*NY);

	// Choose the right component of vf if particles are resolved

#ifdef LAG_PARTICLE_RESOLVED
		if (which_component == 'c'){
			vf = data_bag->lag->ng_vfc;
		}
		else if (which_component == 'u'){
			vf = data_bag->lag->ng_vfu;
		}
		else if (which_component == 'v'){
			vf = data_bag->lag->ng_vfv;
		}
		else if (which_component == 'w'){
			vf = data_bag->lag->ng_vfw;
		}
		else {

		}
#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				flux[j][i] += quantity1[k][j][i] * quantity2[k][j][i] * (1.0-vfrac);
				vf_avg[j][i] += (1.0-vfrac);
			}
		}
	}

	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			flux[j][i] = flux[j][i] / vf_avg[j][i] / params->NPZ;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &flux[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
}



void Statistics2d_partial_y(double **flux, double ***quantity1, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double ***vf;
	double **vf_avg;
	double vfrac = 0.0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

	double ih= grid -> idx_u[Is];
	// Use work buffer for avg volume fraction
	vf_avg = data_bag->st2d->work;

	// Set flux to 0 before MPI_Allreduce is called
	DSET_ZERO(flux[0],NX*NY);
	DSET_ZERO(vf_avg[0],NX*NY);

	// Choose the right component of vf if particles are resolved

#ifdef LAG_PARTICLE_RESOLVED

	vf = data_bag->lag->ng_vfc;

#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				flux[j][i] += 0.5*ih*(quantity1[k][j+1][i]-quantity1[k][j-1][i])*(1.0-vfrac);
				vf_avg[j][i] += (1.0-vfrac);
			}
		}
	}

	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			flux[j][i] = flux[j][i]/(NX-1.)/(NY-1.) ;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &flux[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
}



/******************************************************************************/
/*
 Total kinetic Energy
 Compute kinetic energy:
     \int 0.5 (u^2 + v^2 + w^2) dx dy dz
 using velocities at cell centers.
 */
/******************************************************************************/
double Statistics2d_totKineticEnergy(double **kinE, double ***u_cell,
		double ***v_cell, double ***w_cell, Cart3d_bag *data_bag) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double u,v,w;
	double scale;
	double ***vf;
	double vfrac = 0.0;
	double integral_kinEnergy=0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

	DSET_ZERO(kinE[0],NX*NY);



#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfc;
#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				u = u_cell[k][j][i] * (1.0-vfrac);
				v = v_cell[k][j][i] * (1.0-vfrac);
				w = w_cell[k][j][i] * (1.0-vfrac);
				kinE[j][i] += 0.5 * (u*u + v*v + w*w);
			}
		}
	}

	scale = 1.0/(NZ-1.0);

	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			kinE[j][i] *= scale;
		}
	}


	MPI_Allreduce(MPI_IN_PLACE, &kinE[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);


	for(j = 0; j < NY-1; j++) {
		for(i = 0; i < NX-1; i++) {
			integral_kinEnergy += kinE[j][i];
		}
	}


	return integral_kinEnergy*1.0/(NY-1.0)/(NX-1.0);
}

double Statistics2d_kinE_part(Cart3d_bag *data_bag, Debug_trace *dtrace){
	Particle *p;
	double up,vp,wp, Usq;
	double vol;
	double kinEp = 0.;
	double scale;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	int NZ = grid->NZ;


	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	p = p_mobile_list -> start;
	while (p != NULL) {

		vol = 4/3 * PI * p->R * p->R * p->R;
		up = p->U[0];
		vp = p->U[1];
		wp = p->U[2];
		Usq = up*up + vp*vp + wp*wp;

		kinEp += 0.5 * params->rho_s * vol * Usq;

		p = p -> next;
	}

	scale = 1.0/(NZ-1.0);
	kinEp *= scale;

	MPI_Allreduce(MPI_IN_PLACE, &kinEp, 1, MPI_DOUBLE, MPI_SUM, PCW);

	return kinEp;
}

double Statistics2d_potE_part(Cart3d_bag *data_bag, Debug_trace *dtrace){
	Particle *p;
	double yp;
	double vol;
	double potEp = 0.0;
	double scale;


	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	int NZ = grid->NZ;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	p = p_mobile_list -> start;
	while (p != NULL) {
		vol = 4/3 * PI * p->R * p->R * p->R;
		yp = p->X[1];
		potEp += - params->grav[1] * (params->rho_s - 1.0) * vol * yp;

		p = p -> next;
	}

	scale = 1.0/(NZ-1.0);
	potEp *= scale;

	MPI_Allreduce(MPI_IN_PLACE, &potEp, 1, MPI_DOUBLE, MPI_SUM, PCW);

	return potEp;
}


/******************************************************************************/
/*
 Root Mean Square of Velocity (useful for turbulent flow with zero mean velocity)
 Compute RMS:
    1/V \int sqrt(u^2 + v^2 + w^2) dx dy dz
 using velocities at cell centers.
 */
/******************************************************************************/
double Statistics2d_RMS(double **Urms, double ***u_cell,
		double ***v_cell, double ***w_cell, Cart3d_bag *data_bag) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double u,v,w;
	double scale;
	double ***vf;
	double vfrac = 0.0;
	double integral_RMS=0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

	DSET_ZERO(Urms[0],NX*NY);



#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfc;
#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				u = u_cell[k][j][i] * (1.0-vfrac);
				v = v_cell[k][j][i] * (1.0-vfrac);
				w = w_cell[k][j][i] * (1.0-vfrac);
				Urms[j][i] +=  sqrt((u*u + v*v + w*w)/3);
			}
		}
	}

	scale = 1.0/(NZ-1.0);

	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			Urms[j][i] *= scale;
		}
	}


	MPI_Allreduce(MPI_IN_PLACE, &Urms[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);


	for(j = 0; j < NY-1; j++) {
		for(i = 0; i < NX-1; i++) {
			integral_RMS += Urms[j][i];
		}
	}


	return integral_RMS*1.0/(NY-1.0)/(NX-1.0);
}


/******************************************************************************/
/*
 Work of turbulent forcing (with no particles)
 Compute :
    1/V \int u.f_turb dx dy dz
 */
/******************************************************************************/
double Statistics2d_forcing(double ***u_face, double ***v_face, double ***w_face,
					double ***u_forc, double ***v_forc, double ***w_forc, Cart3d_bag *data_bag) {

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double integral_forcing =0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;


	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

	double dx = grid->dx_u[Is+1];

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
				integral_forcing +=  u_face[k][j][i]*u_forc[k][j][i] + v_face[k][j][i]*v_forc[k][j][i] + w_face[k][j][i]*w_forc[k][j][i];
			}
		}
	}


	MPI_Allreduce(MPI_IN_PLACE, &integral_forcing, 1, MPI_DOUBLE, MPI_SUM, PCW);


	return integral_forcing*1.0/(NX-1.0)/(NY-1.0)/(NZ-1.0);
}


/******************************************************************************/
/*
 Compute potential Energy of concentration field c  E = \int c * Ri* y  dV.
 */
/******************************************************************************/
void Statistics2d_potEnergy_y(double **potE, double ***c, double Ri , double *integral_potEnergy_fluid ,Cart3d_bag *data_bag){

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double scale;
	double u,v,w;
	double ***vf;
	double vfrac = 0.0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
    Statistics2d *st2d = data_bag ->st2d;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;


	double *dx = grid->dx_u;
	double *dy = grid->dy_v;
	double *dz = grid->dz_w;

	double cellvolume= dx[Is] * dy[Js] * dz[Ks];

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

#if(!defined GRID_UNIFORM)
	printf('Function __PRETTY_FUNCTION__ is not ready for nonuniform grid');
#endif
	DSET_ZERO(potE[0],NX*NY);

	double integral_pot=0;

	#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfc;
	#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
	#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
	#endif

				potE[j][i] += grid->yc[j] * c[k][j][i] * Ri * (1.0-vfrac);
				integral_pot += grid->yc[j] * c[k][j][i] * Ri * (1.0-vfrac);

			}
		}
	}


	scale = 1.0/(NZ-1.0);
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			potE[j][i] *= scale;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &potE[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	MPI_Allreduce(MPI_IN_PLACE, &integral_pot, 1, MPI_DOUBLE, MPI_SUM, PCW);

	*integral_potEnergy_fluid =  integral_pot/(NX-1)/(NY-1)/(NZ-1);




	return;
}

/******************************************************************************/
/*
 Compute potential Energy of concentration field c  E = m * g * alpha * y.
 */
/******************************************************************************/
void Statistics2d_potEnergyPart(double **potEpart, Cart3d_bag *data_bag){

	int NX, NY, NZ;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double scale;
	double u,v,w;
	double ***vf;
	double vfrac = 0.0;
	double rho;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	rho = params->rho_s;


	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

	DSET_ZERO(potEpart[0],NX*NY);

	#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfc;
	#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
				vfrac = min(vf[k][j][i],1.0);
				potEpart[j][i] += grid->yc[j] * (rho-1.0) * vfrac;
			}
		}
	}
	scale = 1.0/(NZ-1.0);
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			potEpart[j][i] *= scale;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &potEpart[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);

	return;
}




/******************************************************************************/
/*
 Compute total mass: \int{ c } dx dy dz
 (should be time invariant)
 */
/******************************************************************************/
double Statistics2d_volume_integral(double ***c,   int domain, Cart3d_bag *data_bag) {

	double totalMass = 0;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double *dx, *dy, *dz;
	double cellVolume;
	int NX, NY, NZ;

	MAC_grid   *grid   	= data_bag -> grid;
	Parameters *params 	= data_bag -> params;
	Lagrangian *lag 	= data_bag ->lag;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	dx = grid->dx_u;
	dy = grid->dy_v;
	dz = grid->dz_w;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

#ifdef LAG_PARTICLE_RESOLVED
	double ***ng_vfc= lag->ng_vfc;
#endif


	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
				cellVolume = dx[i] * dy[j] * dz[k];
#ifdef LAG_PARTICLE_RESOLVED

				if ( domain == FLUID_DOMAIN) {
					totalMass += c[k][j][i] * cellVolume*(1-ng_vfc[k][j][i]);
				}
				else if(  domain == PARTICLE_DOMAIN)
				{
					totalMass += c[k][j][i] * cellVolume*ng_vfc[k][j][i];
				}
				else
				{
					totalMass += c[k][j][i] * cellVolume;
				}
#else
				totalMass += c[k][j][i] * cellVolume;
#endif
			}
		}
	}


	MPI_Allreduce(MPI_IN_PLACE, &totalMass,1, MPI_DOUBLE, MPI_SUM, PCW);

	return totalMass;
}



/******************************************************************************/
/*
 Compute pressure work: \int{u_i dp/dx_i}
 */
/******************************************************************************/
void Statistics2d_pressureWork(double ***u, double ***v, double ***w, Pressure *p,
		double **pwork, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double scale;
	double ***dp_dx, ***dp_dy, ***dp_dz;
	double ***p_data;
	double dpdx, dpdy, dpdz;
	double dpdx_ip, dpdy_jp, dpdz_kp;
	double *idx_c, *idy_c, *idz_c;
	double ***vf;
	double vfrac = 0.0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Ie = min(Ie, NX-1);
	Je = min(Je, NY-1);
	Ke = min(Ke, NZ-1);

	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;

	// get p-gradients
	p_data = p->p_data;

	DSET_ZERO(pwork[0],NX*NY);

#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfc;
#endif

	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif

#ifdef XPERIODIC
				if (i==0)
					dpdx = (p_data[k][j][i]-p_data[k][j][i-1]) * idx_c[i-1];
				else
					dpdx = (p_data[k][j][i]-p_data[k][j][i-1]) * idx_c[0];
#else
				if (i==0)
					dpdx = 0.0;
				else
					dpdx = (p_data[k][j][i]-p_data[k][j][i-1]) * idx_c[i-1];
#endif
				dpdx_ip = (p_data[k][j][i+1]-p_data[k][j][i]) * idx_c[i];

				if (j==0)
					dpdy = 0.0;
				else
					dpdy = (p_data[k][j][i]-p_data[k][j-1][i]) * idy_c[j-1];

				dpdy_jp = (p_data[k][j+1][i]-p_data[k][j][i]) * idy_c[j];

#ifdef XPERIODIC
				if (k==0)
					dpdz = (p_data[k][j][i]-p_data[k-1][j][i]) * idz_c[0];
				else
					dpdz = (p_data[k][j][i]-p_data[k-1][j][i]) * idz_c[k-1];
#else
				if (k==0)
					dpdz = 0.0;
				else
					dpdz = (p_data[k][j][i]-p_data[k-1][j][i]) * idz_c[k-1];
#endif
				dpdz_kp = (p_data[k+1][j][i]-p_data[k][j][i]) * idz_c[k];

				pwork[j][i] += (dpdx*u[k][j][i] + dpdx_ip*u[k][j][i+1]
				              + dpdy*v[k][j][i] + dpdy_jp*v[k][j+1][i]
				              + dpdz*w[k][j][i] + dpdz_kp*w[k+1][j][i]);
				pwork[j][i] *= (1.0-vfrac);
			}
		}
	}
	scale = 1.0/(NZ-1.0);
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			pwork[j][i] *= 0.5 * scale;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &pwork[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Statistics2d_uface_to_cellcenter(double **face, double **cell, MAC_grid *grid) {

	int NX, NY, NZ;
	int i, j, k;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	for (j=0;j<NY-1;j++){
		for (i=0;i<NX-1;i++){
			cell[j][i] = 0.5*(face[j][i]+face[j][i+1]);
		}
	}

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Statistics2d_vface_to_cellcenter(double **face, double **cell, MAC_grid *grid) {

	int NX, NY, NZ;
	int i, j, k;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	for (j=0;j<NY-1;j++){
		for (i=0;i<NX-1;i++){
			cell[j][i] = 0.5*(face[j][i]+face[j+1][i]);
		}
	}

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Statistics2d_add_array(double **work, double **work1, MAC_grid *grid) {

	int NX, NY, NZ;
	int i, j, k;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	for (j=0;j<NY-1;j++){
		for (i=0;i<NX-1;i++){
			work1[j][i] += work[j][i];
		}
	}

	return;
}




/******************************************************************************/
/*
 Calculate the all the terms in the KE transport equation at the cell center
 */
/******************************************************************************/
void Statistics2d_ke_budget(double ***u_data, double ***v_data, double ***w_data,
		double ***u_data_bc, double ***v_data_bc, double ***w_data_bc,
		double ***nut, Wall_model *log_law, double **global_viscDiss, double **global_viscTran,
		double **global_sgsDiss, double **global_sgsTran, Cart3d_bag *data_bag, double iRe) {



	int i, j, k;
	int NX, NY, NZ;
	double **viscDiss, **viscTran, **sgsDiss, **sgsTran;


	Statistics2d *st2d = data_bag -> st2d;
	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	// Use work buffer for avg volume fraction
	viscDiss = data_bag->st2d->work1;
	viscTran = data_bag->st2d->work2;
	sgsDiss = data_bag->st2d->work3;
	sgsTran = data_bag->st2d->work4;

	// Set global data to 0 before MPI_Allreduce is called
	DSET_ZERO(global_viscDiss[0],NX*NY);
	DSET_ZERO(global_viscTran[0],NX*NY);
	DSET_ZERO(global_sgsDiss[0],NX*NY);
	DSET_ZERO(global_sgsTran[0],NX*NY);

	//--------------------------------------------------------------------------
	// Calculate the contribution to KE transport from uu/2 transport equation
	//--------------------------------------------------------------------------
	Statistics2d_ke_budget_u(u_data, v_data, w_data, u_data_bc, nut, log_law, viscDiss, viscTran, sgsDiss, sgsTran, data_bag, iRe);

	MPI_Allreduce(MPI_IN_PLACE, &viscDiss[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(viscDiss, viscDiss, grid);
	Statistics2d_add_array(viscDiss, global_viscDiss, grid);

	MPI_Allreduce(MPI_IN_PLACE, &viscTran[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(viscTran, viscTran, grid);
	Statistics2d_add_array(viscTran, global_viscTran, grid);

#ifdef LES
	MPI_Allreduce(MPI_IN_PLACE, &sgsDiss[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(sgsDiss, sgsDiss, grid);
	Statistics2d_add_array(sgsDiss, global_sgsDiss, grid);

	MPI_Allreduce(MPI_IN_PLACE, &sgsTran[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(sgsTran, sgsTran, grid);
	Statistics2d_add_array(sgsTran, global_sgsTran, grid);
#endif
	//--------------------------------------------------------------------------
	// Calculate the contribution to KE transport from vv/2 transport equation
	//--------------------------------------------------------------------------
	Statistics2d_ke_budget_v(u_data, v_data, w_data, v_data_bc, nut, log_law, viscDiss, viscTran, sgsDiss, sgsTran, data_bag, iRe);

	MPI_Allreduce(MPI_IN_PLACE, &viscDiss[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(viscDiss, viscDiss, grid);
	Statistics2d_add_array(viscDiss, global_viscDiss, grid);

	MPI_Allreduce(MPI_IN_PLACE, &viscTran[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(viscTran, viscTran, grid);
	Statistics2d_add_array(viscTran, global_viscTran, grid);
#ifdef LES
	MPI_Allreduce(MPI_IN_PLACE, &sgsDiss[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(sgsDiss, sgsDiss, grid);
	Statistics2d_add_array(sgsDiss, global_sgsDiss, grid);

	MPI_Allreduce(MPI_IN_PLACE, &sgsTran[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(sgsTran, sgsTran, grid);
	Statistics2d_add_array(sgsTran, global_sgsTran, grid);
#endif
	//--------------------------------------------------------------------------
	// Calculate the contribution to KE transport from ww/2 transport equation
	//--------------------------------------------------------------------------
	Statistics2d_ke_budget_w(u_data, v_data, w_data, w_data_bc, nut, log_law, viscDiss, viscTran, sgsDiss, sgsTran, data_bag, iRe);

	MPI_Allreduce(MPI_IN_PLACE, &viscDiss[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(viscDiss, viscDiss, grid);
	Statistics2d_add_array(viscDiss, global_viscDiss, grid);

	MPI_Allreduce(MPI_IN_PLACE, &viscTran[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(viscTran, viscTran, grid);
	Statistics2d_add_array(viscTran, global_viscTran, grid);

#ifdef LES
	MPI_Allreduce(MPI_IN_PLACE, &sgsDiss[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(sgsDiss, sgsDiss, grid);
	Statistics2d_add_array(sgsDiss, global_sgsDiss, grid);

	MPI_Allreduce(MPI_IN_PLACE, &sgsTran[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);
	Statistics2d_uface_to_cellcenter(sgsTran, sgsTran, grid);
	Statistics2d_add_array(sgsTran, global_sgsTran, grid);
#endif
	for (i=0;i<NX-1;i++){
		global_viscDiss[NY-1][i] = global_viscDiss[NY-2][i];
		global_viscTran[NY-1][i] = global_viscTran[NY-2][i];
#ifdef LES
		global_sgsDiss[NY-1][i] = global_sgsDiss[NY-2][i];
		global_sgsTran[NY-1][i] = global_sgsTran[NY-2][i];
#endif
	}

	for (j=0;j<NY;j++){
		global_viscDiss[j][NX-1] = global_viscDiss[j][NX-2];
		global_viscTran[j][NX-1] = global_viscTran[j][NX-2];


#ifdef LES
		global_sgsDiss[j][NX-1] = global_sgsDiss[j][NX-2];
		global_sgsTran[j][NX-1] = global_sgsTran[j][NX-2];
#endif
	}

	st2d->integral_viscDiss=0;

	for (j=0;j<NY-1;j++){
		for (i=0;i<NX-1;i++){
			st2d->integral_viscDiss  += global_viscDiss[j][i];
		}
	}

	st2d->integral_viscDiss *= 1./(NX-1.)/(NY-1.);

	return;

}




/******************************************************************************/
/*
 Compute dissipation & transport (function of x & y):

     vis_eps_ke_due_to_u  = d(\tau_1j u)/dx_j - {u (d2udx2 + d2udy2 + d2udz2) }
     vis_tran_ke_due_to_u = d(\tau_1j u)/dx_j
 \tau_1j in the above two lines is the viscous stress

     sgs_eps_ke_due_to_u  = d(\tau_1j u)/dx_j - {u (d(\tau_1j)/dx_j) }
     sgs_tran_ke_due_to_u = d(\tau_1j u)/dx_j
 \tau_1j in the above two lines is the SGS stress

 The second derivatives are approximated by finite differences.
*/
/******************************************************************************/
void Statistics2d_ke_budget_u(double ***u_data, double ***v_data, double ***w_data,
		double ***u_data_bc, double ***nut, Wall_model *log_law,
		double **viscDiss, double **viscTran, double **sgsDiss,
		double **sgsTran, Cart3d_bag *data_bag, double iRe) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double dudxE, dudxW, dudyN, dudyS, dudzF, dudzB;
	double d2udx2, d2udy2, d2udz2, ddxdudx, ddydvdx,  ddzdwdx;
	double dvdxN, dvdxS, dwdxF, dwdxB;
	double uE, uW, uN, uS, uF, uB;
	double nutE, nutW, nutN, nutS, nutF, nutB;
	double scale;
	double ***vf;
	double vfrac = 0.0;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;
	double vd2udx2, vd2udy2, vd2udz2, vddxdudx, vddydvdx,  vddzdwdx;
	double sd2udx2, sd2udy2, sd2udz2, sddxdudx, sddydvdx,  sddzdwdx;
	double vtd2udx2, vtd2udy2, vtd2udz2, vtddxdudx, vtddydvdx,  vtddzdwdx;
	double std2udx2, std2udy2, std2udz2, stddxdudx, stddydvdx,  stddzdwdx;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

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
	i_start = max(1,Is); // i=0, i=NX-1 are not included
	j_start = Js;
	k_start = Ks;

	// exclude the half cell added
	i_end   = min(NX-1, Ie);
	j_end   = min(NY-1, Je);
	k_end   = min(NZ-1, Ke);

#ifdef XPERIODIC
	if (i_end == NX-1)
		i_end = i_end + 1;
#endif

	idx_u = grid->idx_u;
	idy_v = grid->idy_v;
	idz_w = grid->idz_w;
	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;

	DSET_ZERO(viscDiss[0],NX*NY);
	DSET_ZERO(viscTran[0],NX*NY);
	DSET_ZERO(sgsDiss[0],NX*NY);
	DSET_ZERO(sgsTran[0],NX*NY);

#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfu;
#endif

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif

#ifdef LES
				//--------------------------------------------------------------
				// Eddy viscosity
				//--------------------------------------------------------------
				nutE = nut[k][j][i];
				nutW = nut[k][j][i-1];
				nutN = ( 0.25 * ( nut[k][j][i] + nut[k][j][i-1] + nut[k][j+1][i] + nut[k][j+1][i-1] ) );
				if ( j != 0) {
					nutS = ( 0.25 * ( nut[k][j][i] + nut[k][j][i-1] + nut[k][j-1][i] + nut[k][j-1][i-1] ) );
				}
				else {
	#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
					nutS =  0.0;    //Zero Eddy viscosity on the wall
	#else
					nutS = ( 0.5 * ( nut[k][j][i] + nut[k][j][i-1] ) );
	#endif
				}

				nutF = ( 0.25 * ( nut[k][j][i] + nut[k][j][i-1] + nut[k+1][j][i] + nut[k+1][j][i-1] ) );

	#ifdef ZPERIODIC
				nutB = ( 0.25 * ( nut[k][j][i] + nut[k][j][i-1] + nut[k-1][j][i] + nut[k-1][j][i-1] ) );
	#else
				if ( k != 0) {
					nutB = ( 0.25 * ( nut[k][j][i] + nut[k][j][i-1] + nut[k-1][j][i] + nut[k-1][j][i-1] ) );
				}
				else {
					nutB = ( 0.5 * ( nut[k][j][i] + nut[k][j][i-1]  ) );
				}
	#endif
#else // not LES
				nutE = 0.0;
				nutW = 0.0;
				nutN = 0.0;
				nutS = 0.0;
				nutF = 0.0;
				nutB = 0.0;
#endif // LES

				//--------------------------------------------------------------
				// d2u/dx2, d2u/dy2 and d2u/dz2
				//--------------------------------------------------------------
				dudxE = ( u_data[k][j][i+1] - u_data[k][j][i] ) * idx_u[i];
				dudxW = ( u_data[k][j][i] - u_data[k][j][i-1] ) * idx_u[i-1];

				if (j != NY-2 ) {
					dudyN = ( u_data[k][j+1][i] - u_data[k][j][i] ) * idy_c[j];
				}
				else {
#ifdef TOP_WALL_VELOCITY_NOSLIP
					// Calculate dudyN with ghost cell values (i.e u(NY-1) = -u(NY-2))
					dudyN = -2.*u_data[k][j][i]*idy_v[j];
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
					dudyN = 0.;
#endif
#ifdef TOP_WALL_SCHUMANN
					dudyN = log_law->dudy_wm_top[k][i];
#endif

				}

				if (j !=0 ) {
					dudyS = ( u_data[k][j][i] - u_data[k][j-1][i] ) * idy_c[j-1];
				}
				else {
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
					// Calculate dudyS with ghost cell values (i.e u(-1) = -u(0))
					dudyS =  2.*u_data[k][j][i]*idy_v[j];
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
					dudyS = 0.;
#endif
#ifdef BOTTOM_WALL_SCHUMANN
					dudyS = log_law->dudy_wm_bottom[k][i];
#endif
				}

				if (k != NZ-2) {
					dudzF = ( u_data[k+1][j][i] - u_data[k][j][i] ) * idz_c[k];
				}
				else {
#ifdef FRONT_WALL_VELOCITY_NOSLIP
					dudzF = -2.*u_data[k][j][i]*idz_w[k];
#endif
#ifdef FRONT_WALL_VELOCITY_FREESLIP
					dudzF = 0.;
#endif
#ifdef ZPERIODIC
					dudzF = ( u_data[k+1][j][i] - u_data[k][j][i] ) * idz_c[k];
#endif
				}

				if (k !=0 ) {
					dudzB = ( u_data[k][j][i] - u_data[k-1][j][i] ) * idz_c[k-1];
				}
				else {
#ifdef BACK_WALL_VELOCITY_NOSLIP
					dudzB = 2.*( u_data[k][j][i] ) * idz_w[k];
#endif
#ifdef BACK_WALL_VELOCITY_FREESLIP
					dudzB = 0.;
#endif
#ifdef ZPERIODIC
					dudzB = ( u_data[k][j][i] - u_data[k-1][j][i] ) * idz_c[k];
#endif
				}

				sd2udx2 = ( nutE * dudxE - nutW * dudxW ) * idx_c[i-1];
				sd2udy2 = ( nutN * dudyN - nutS * dudyS ) * idy_v[j];
				sd2udz2 = ( nutF * dudzF - nutB * dudzB ) * idz_w[k];

				vd2udx2 = iRe * (  dudxE -  dudxW ) * idx_c[i-1];
				vd2udy2 = iRe * (  dudyN -  dudyS ) * idy_v[j];
				vd2udz2 = iRe * (  dudzF -  dudzB ) * idz_w[k];

				//--------------------------------------------------------------
				// d/dx(du/dx), d/dy(dv/dx) and d/dz(dw/dx)
				//--------------------------------------------------------------
				dvdxN = ( v_data[k][j+1][i] - v_data[k][j+1][i-1] ) * idx_c[i-1];
				dvdxS = ( v_data[k][j][i] - v_data[k][j][i-1] ) * idx_c[i-1];
				dwdxF = ( w_data[k+1][j][i] - w_data[k+1][j][i-1] ) * idx_c[i-1];
				dwdxB = ( w_data[k][j][i] - w_data[k][j][i-1] ) * idx_c[i-1];

				sddxdudx = sd2udx2;
				sddydvdx = ( nutN * dvdxN - nutS * dvdxS ) * idy_v[j];
				sddzdwdx = ( nutF * dwdxF - nutB * dwdxB ) * idz_w[k];

				vddxdudx = vd2udx2;
				vddydvdx = iRe * (  dvdxN -  dvdxS ) * idy_v[j];
				vddzdwdx = iRe * (  dwdxF -  dwdxB ) * idz_w[k];

				viscDiss[j][i] += u_data[k][j][i] * (vd2udx2 + vd2udy2 + vd2udz2) * (1.0-vfrac);
				viscDiss[j][i] += u_data[k][j][i] * (vddxdudx + vddydvdx + vddzdwdx) * (1.0-vfrac);
#ifdef LES
				sgsDiss[j][i] += u_data[k][j][i] * (sd2udx2 + sd2udy2 + sd2udz2) * (1.0-vfrac);
				sgsDiss[j][i] += u_data[k][j][i] * (sddxdudx + sddydvdx + sddzdwdx) * (1.0-vfrac);
#endif

				//--------------------------------------------------------------
				// Calculate the transport terms
				//--------------------------------------------------------------
				uE = u_data_bc[k][j][i];
				uW = u_data_bc[k][j][i-1];
				uN = 0.5*(u_data[k][j][i] + u_data[k][j+1][i]);
				if (j!=0) {
					uS = 0.5*(u_data[k][j][i] + u_data[k][j-1][i]);
				}
				else {
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
					uS = u_data[k][j][i];
#else
					uS = 0.0;
#endif
				}
				uF = 0.5*(u_data[k][j][i] + u_data[k+1][j][i]);
				if (k!=0) {
					uB = 0.5*(u_data[k][j][i] + u_data[k-1][j][i]);
				}
				else {
#ifdef BACK_WALL_VELOCITY_FREESLIP
					uB = u_data[k][j][i];
#else
					uB = 0.0;
#endif
				}

				std2udx2 = ( uE*nutE * dudxE - uW*nutW * dudxW ) * idx_c[i-1];
				std2udy2 = ( uN*nutN * dudyN - uS*nutS * dudyS ) * idy_v[j];
				std2udz2 = ( uF*nutF * dudzF - uB*nutB * dudzB ) * idz_w[k];

				stddxdudx = std2udx2;
				stddydvdx = ( uN*nutN * dvdxN - uS*nutS * dvdxS ) * idy_v[j];
				stddzdwdx = ( uF*nutF * dwdxF - uB*nutB * dwdxB ) * idz_w[k];

				vtd2udx2 = iRe * ( uE* dudxE - uW* dudxW ) * idx_c[i-1];
				vtd2udy2 = iRe * ( uN* dudyN - uS* dudyS ) * idy_v[j];
				vtd2udz2 = iRe * ( uF* dudzF - uB* dudzB ) * idz_w[k];

				vtddxdudx = vtd2udx2;
				vtddydvdx = iRe * ( uN* dvdxN - uS* dvdxS ) * idy_v[j];
				vtddzdwdx = iRe * ( uF* dwdxF - uB* dwdxB ) * idz_w[k];

				viscTran[j][i] += (vtd2udx2 + vtd2udy2 + vtd2udz2) * (1.0-vfrac);
				viscTran[j][i] += (vtddxdudx + vtddydvdx + vtddzdwdx) * (1.0-vfrac);
#ifdef LES
				sgsTran[j][i] += (std2udx2 + std2udy2 + std2udz2) * (1.0-vfrac) * (1.0-vfrac);
				sgsTran[j][i] += (stddxdudx + stddydvdx + stddzdwdx) * (1.0-vfrac) * (1.0-vfrac);
#endif

			} // for i
		} // for j
	} // for k

	scale = 1.0/(NZ-1.0);
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			viscDiss[j][i] *= scale;
			viscTran[j][i] *= scale;
			sgsDiss[j][i]  *= scale;
			sgsTran[j][i]  *= scale;
			viscDiss[j][i] = viscTran[j][i] - viscDiss[j][i];
			sgsDiss[j][i]  = sgsTran[j][i]  - sgsDiss[j][i];
		}
	}

	if (Is == 1) {
		for(j=Js; j<Je; j++) {
			viscDiss[j][0] = viscDiss[j][1];
			viscTran[j][0] = viscTran[j][1];
			sgsDiss[j][0] = sgsDiss[j][1];
			sgsTran[j][0] = sgsTran[j][1];
		}
	}

	if (Ie == NX-1) {
		for (j=Js; j<Je; j++) {
			viscDiss[j][NX-1] = viscDiss[j][NX-2];
			viscTran[j][NX-1] = viscTran[j][NX-2];
			sgsDiss[j][NX-1] = sgsDiss[j][NX-2];
			sgsTran[j][NX-1] = sgsTran[j][NX-2];
		}
	}

	if (Je == NY-1) {
		for (i=grid->G_Is; i<Ie; i++) {
			viscDiss[NY-1][i] = viscDiss[NY-2][i];
			viscTran[NY-1][i] = viscTran[NY-2][i];
			sgsDiss[NY-1][i] = sgsDiss[NY-2][i];
			sgsTran[NY-1][i] = sgsTran[NY-2][i];
		}
	}

	return;
}




/******************************************************************************/
/*
 Compute dissipation & transport (function of x & y):

     vis_eps_ke_due_to_v   = d(\tau_2j v)/dx_j - {v (d2vdx2 + d2vdy2 + d2vdz2) }
     vis_tran_ke_due_to_v  = d(\tau_2j v)/dx_j
 \tau_2j in the above two lines is the viscous stress

     sgs_eps_ke_due_to_v   = d(\tau_2j v)/dx_j - {v (d(\tau_2j)/dx_j) }
     sgs_tran_ke_due_to_v  = d(\tau_2j v)/dx_j
 \tau_2j in the above two lines is the SGS stress

 The second derivatives are approximated by finite differences.
 */
/******************************************************************************/
void Statistics2d_ke_budget_v(double ***u_data, double ***v_data, double ***w_data,
		double ***v_data_bc, double ***nut, Wall_model *log_law,
		double **viscDiss, double **viscTran, double **sgsDiss,
		double **sgsTran, Cart3d_bag *data_bag, double iRe) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double dvdxE, dvdxW, dvdyN, dvdyS, dvdzF, dvdzB;
	double dudyE, dudyW, dwdyF, dwdyB;
	double vd2vdx2, vd2vdy2, vd2vdz2, vddxdudy, vddydvdy,  vddzdwdy;
	double sd2vdx2, sd2vdy2, sd2vdz2, sddxdudy, sddydvdy,  sddzdwdy;
	double vtd2vdx2, vtd2vdy2, vtd2vdz2, vtddxdudy, vtddydvdy,  vtddzdwdy;
	double std2vdx2, std2vdy2, std2vdz2, stddxdudy, stddydvdy,  stddzdwdy;
	double vE, vW, vN, vS, vF, vB;
	double nutE, nutW, nutN, nutS, nutB, nutF;
	double scale;
	double ***vf;
	double vfrac = 0.0;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

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
	j_start = max(1, Js);// j=0, j=NY-1 are not included
	k_start = Ks;

	// exclude the half cell added
	i_end   = min(NX-1, Ie);
	j_end   = min(NY-1, Je);
	k_end   = min(NZ-1, Ke);

	idx_u = grid->idx_u;
	idy_v = grid->idy_v;
	idz_w = grid->idz_w;
	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;

	DSET_ZERO(viscDiss[0],NX*NY);
	DSET_ZERO(viscTran[0],NX*NY);
	DSET_ZERO(sgsDiss[0],NX*NY);
	DSET_ZERO(sgsTran[0],NX*NY);

#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfv;
#endif


	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				//--------------------------------------------------------------
				// Eddy viscosity at boundary points
				//--------------------------------------------------------------
#ifdef LES
				nutE = 0.25*( nut[k][j][i] + nut[k][j][i+1] + nut[k][j-1][i] + nut[k][j-1][i+1] );
	#ifdef XPERIODIC
				if (i != 0) {
					nutW = 0.25*( nut[k][j][i] + nut[k][j][i-1] + nut[k][j-1][i] + nut[k][j-1][i-1] );
				}
				else {
					nutW = 0.5*( nut[k][j][i] +  nut[k][j-1][i] );
				}
	#else
				nutW = 0.25*( nut[k][j][i] + nut[k][j][i-1] + nut[k][j-1][i] + nut[k][j-1][i-1] );
	#endif
				nutN = nut[k][j][i];
				nutS = nut[k][j-1][i];

				nutF = ( 0.25 * ( nut[k][j][i] + nut[k+1][j][i] + nut[k][j-1][i] + nut[k+1][j-1][i] ) );

				if ( k != 0) {
					nutB = ( 0.25 * ( nut[k][j][i] + nut[k-1][j][i] + nut[k][j-1][i] + nut[k-1][j-1][i] ) );
				}
				else {
					nutB = ( 0.5 * ( nut[k][j][i] +  nut[k][j-1][i] ) );
				}

#else // not LES
				nutE = 0.0;
				nutW = 0.0;
				nutN = 0.0;
				nutS = 0.0;
				nutF = 0.0;
				nutB = 0.0;
#endif // LES

				//--------------------------------------------------------------
				// d2v/dx2, d2v/dy2 and d2v/dz2
				//--------------------------------------------------------------
				dvdxE = ( v_data[k][j][i+1] - v_data[k][j][i] ) * idx_c[i];
				if (i != 0 ) {
					dvdxW = ( v_data[k][j][i] - v_data[k][j][i-1] ) * idx_c[i-1];
				}
				else {
#ifdef XPERIODIC
					dvdxW = ( v_data[k][j][i] - v_data[k][j][i-1] ) * idx_c[i];
#else
					// v = 0 at i=0 face
					dvdxW = 2*( v_data[k][j][i] ) * idx_u[i];
#endif
				}

				// No special treatment needed for dvdyN, dvdyS (unlike u),
				// since v(0)=v(NY-1)=0 as part of the solution
				dvdyN = ( v_data[k][j+1][i] - v_data[k][j][i] ) * idy_v[j];
				dvdyS = ( v_data[k][j][i] - v_data[k][j-1][i] ) * idy_v[j-1];

#ifdef ZPERIODIC
				dvdzF = ( v_data[k+1][j][i] - v_data[k][j][i] ) * idz_c[k];
#else
				if (k != NZ-2) {
					dvdzF = ( v_data[k+1][j][i] - v_data[k][j][i] ) * idz_c[k];
				}
				else {
	#ifdef FRONT_WALL_VELOCITY_NOSLIP
					dvdzF = - 2.*v_data[k][j][i] * idz_w[k];
	#endif
	#ifdef FRONT_WALL_VELOCITY_FREESLIP
					dvdzF = 0.;
	#endif
				}
#endif // ZPERIODIC

				if (k !=0 ) {
					dvdzB = ( v_data[k][j][i] - v_data[k-1][j][i] ) * idz_c[k-1];
				}
				else {
#ifdef BACK_WALL_VELOCITY_NOSLIP
					dvdzB = 2.*( v_data[k][j][i] ) * idz_w[k];
#endif
#ifdef BACK_WALL_VELOCITY_FREESLIP
					dvdzB = 0.;
#endif
#ifdef ZPERIODIC
					dvdzB = ( v_data[k][j][i] - v_data[k-1][j][i] ) * idz_c[k];
#endif
				}

				sd2vdx2 = ( nutE * dvdxE - nutW * dvdxW ) * idx_u[i];
				sd2vdy2 = ( nutN * dvdyN - nutS * dvdyS ) * idy_c[j-1];
				sd2vdz2 = ( nutF * dvdzF - nutB * dvdzB ) * idz_w[k];

				vd2vdx2 = iRe * ( dvdxE -  dvdxW ) * idx_u[i];
				vd2vdy2 = iRe * ( dvdyN -  dvdyS ) * idy_c[j-1];
				vd2vdz2 = iRe * ( dvdzF -  dvdzB ) * idz_w[k];

				//--------------------------------------------------------------
				// d/dx(du/dy), d/dy(dv/dy), d/dz(dw/dy)
				//--------------------------------------------------------------
				dudyE = ( u_data[k][j][i+1] - u_data[k][j-1][i+1] ) * idy_c[j-1];
				dudyW = ( u_data[k][j][i] - u_data[k][j-1][i] ) * idy_c[j-1];
				dwdyF = ( w_data[k+1][j][i] - w_data[k+1][j-1][i] ) * idy_c[j-1];
				dwdyB = ( w_data[k][j][i] - w_data[k][j-1][i] ) * idy_c[j-1];

				sddydvdy = sd2vdy2;
				sddxdudy = ( nutE * dudyE - nutW * dudyW ) * idx_u[i];
				sddzdwdy = ( nutF * dwdyF - nutB * dwdyB ) * idz_w[k];

				vddydvdy = vd2vdy2;
				vddxdudy = iRe * ( dudyE -  dudyW ) * idx_u[i];
				vddzdwdy = iRe * ( dwdyF -  dwdyB ) * idz_w[k];

				viscDiss[j][i] += v_data[k][j][i] * (vd2vdx2 + vd2vdy2 + vd2vdz2) * (1.0-vfrac) ;
				viscDiss[j][i] += v_data[k][j][i] * (vddydvdy + vddxdudy + vddzdwdy) * (1.0-vfrac);
#ifdef LES
				sgsDiss[j][i] += v_data[k][j][i] * (sd2vdx2 + sd2vdy2 + sd2vdz2) * (1.0-vfrac) ;
				sgsDiss[j][i] += v_data[k][j][i] * (sddydvdy + sddxdudy + sddzdwdy) * (1.0-vfrac) ;
#endif

				//--------------------------------------------------------------
				// Calculate the transport terms
				//--------------------------------------------------------------
				vE = 0.5*(v_data[k][j][i] + v_data[k][j][i+1]);
#ifdef XPERIODIC
				vW = 0.5*(v_data[k][j][i] + v_data[k][j][i-1]);
#else
				if (i != 0) {
					vW = 0.5*(v_data[k][j][i] + v_data[k][j][i-1]);
				}
				else {
					vW = 0.0;
				}
#endif

				vN = v_data_bc[k][j][i];
				vS = v_data_bc[k][j-1][i];
				vF = 0.5*(v_data[k][j][i] + v_data[k+1][j][i]);
#ifdef ZPERIODIC
				vB = 0.5*(v_data[k][j][i] + v_data[k-1][j][i]);
#else
				if (k != 0) {
						vB = 0.5*(v_data[k][j][i] + v_data[k-1][j][i]);
				}
				else {
	#ifdef BACK_WALL_VELOCITY_FREESLIP
					vB = 0.0;
	#endif
	#ifdef BACK_WALL_VELOCITY_NOSLIP
					vB = 0.0;
	#endif
				}
#endif // ZPERIODIC

				vtd2vdx2 = iRe * ( vE* dvdxE - vW* dvdxW ) * idx_u[i];
				vtd2vdy2 = iRe * ( vN* dvdyN - vS* dvdyS ) * idy_c[j-1];
				vtd2vdz2 = iRe * ( vF* dvdzF - vB* dvdzB ) * idz_w[k];

				vtddydvdy = vtd2vdy2;
				vtddxdudy = iRe * ( vE* dudyE - vW* dudyW ) * idx_u[i];
				vtddzdwdy = iRe * ( vF* dwdyF - vB* dwdyB ) * idz_w[k];

				viscTran[j][i] += (vtd2vdx2 + vtd2vdy2 + vtd2vdz2) * (1.0-vfrac);
				viscTran[j][i] += (vtddydvdy + vtddxdudy + vtddzdwdy) * (1.0-vfrac);
#ifdef LES
				std2vdx2 = ( vE*nutE * dvdxE - vW*nutW * dvdxW ) * idx_u[i];
				std2vdy2 = ( vN*nutN * dvdyN - vS*nutS * dvdyS ) * idy_c[j-1];
				std2vdz2 = ( vF*nutF * dvdzF - vB*nutB * dvdzB ) * idz_w[k];

				stddydvdy = std2vdy2;
				stddxdudy = ( vE*nutE * dudyE - vW*nutW * dudyW ) * idx_u[i];
				stddzdwdy = ( vF*nutF * dwdyF - vB*nutB * dwdyB ) * idz_w[k];

				sgsTran[j][i]  += (std2vdx2 + std2vdy2 + std2vdz2) * (1.0-vfrac) ;
				sgsTran[j][i]  += (stddydvdy + stddxdudy + stddzdwdy) * (1.0-vfrac) ;
#endif

			}
		}
	}

	scale = 1.0/(NZ-1.0);
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			viscDiss[j][i] *= scale;
			viscTran[j][i] *= scale;
			sgsDiss[j][i]  *= scale;
			sgsTran[j][i]  *= scale;
			viscDiss[j][i] = viscTran[j][i] - viscDiss[j][i];
			sgsDiss[j][i]  = sgsTran[j][i]  - sgsDiss[j][i];
		}
	}

	if (Js == 1) {
		for(i=Is; i<Ie; i++) {
			viscDiss[0][i] = viscDiss[1][i];
			viscTran[0][i] = viscTran[1][i];
			sgsDiss[0][i] = sgsDiss[1][i];
			sgsTran[0][i] = sgsTran[1][i];
		}
	}

	if (Je == NY-1) {
		for (i=Is; i<Ie; i++) {
			viscDiss[NY-1][i] = viscDiss[NY-2][i];
			viscTran[NY-1][i] = viscTran[NY-2][i];
			sgsDiss[NY-1][i] = sgsDiss[NY-2][i];
			sgsTran[NY-1][i] = sgsTran[NY-2][i];
		}
	}

	if (Ie == NX-1) {
		for (j=grid->G_Js; j<Je; j++) {
			viscDiss[j][NX-1] = viscDiss[j][NX-2];
			viscTran[j][NX-1] = viscTran[j][NX-2];
			sgsDiss[j][NX-1] = sgsDiss[j][NX-2];
			sgsTran[j][NX-1] = sgsTran[j][NX-2];
		}
	}

	return;
}




/******************************************************************************/
/*
 Compute dissipation & transport (function of x & y):

     vis_eps_ke_due_to_w   = d(\tau_3j w)/dx_j - {w (d2wdx2 + d2wdy2 + d2wdz2) }
     vis_tran_ke_due_to_w  = d(\tau_3j w)/dx_j
 \tau_3j in the above two lines is the viscous stress

     sgs_eps_ke_due_to_w   = d(\tau_3j w)/dx_j - {w (d(\tau_3j)/dx_j) }
     sgs_tran_ke_due_to_w  = d(\tau_3j w)/dx_j
 \tau_3j in the above two lines is the SGS stress

 The second derivatives are approximated by finite differences.
 */
/******************************************************************************/
void Statistics2d_ke_budget_w(double ***u_data, double ***v_data, double ***w_data,
		double ***w_data_bc, double ***nut, Wall_model *log_law,
		double **viscDiss, double **viscTran, double **sgsDiss,
		double **sgsTran, Cart3d_bag *data_bag, double iRe) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double dwdxE, dwdxW, dwdyN, dwdyS, dwdzF, dwdzB;
	double dudzE, dudzW, dvdzN, dvdzS;
	double vd2wdx2, vd2wdy2, vd2wdz2, vddxdudz, vddydvdz,  vddzdwdz;
	double sd2wdx2, sd2wdy2, sd2wdz2, sddxdudz, sddydvdz,  sddzdwdz;
	double vtd2wdx2, vtd2wdy2, vtd2wdz2, vtddxdudz, vtddydvdz,  vtddzdwdz;
	double std2wdx2, std2wdy2, std2wdz2, stddxdudz, stddydvdz,  stddzdwdz;
	double scale;
	double ***vf;
	double vfrac = 0.0;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	double nutE, nutW, nutN, nutS, nutF, nutB;
	double wE, wW, wN, wS, wF, wB;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

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
	j_start = Js;// j=0, j=NY-1 are not included
	k_start = max(1,Ks);

	// exclude the half cell added
	i_end   = min(NX-1, Ie);
	j_end   = min(NY-1, Je);
	k_end   = min(NZ-1, Ke);
#ifdef ZPERIODIC
	if (k_end == NZ-1)
		k_end = k_end + 1;
#endif

	idx_u = grid->idx_u;
	idy_v = grid->idy_v;
	idz_w = grid->idz_w;
	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;

	DSET_ZERO(viscDiss[0],NX*NY);
	DSET_ZERO(viscTran[0],NX*NY);
	DSET_ZERO(sgsDiss[0],NX*NY);
	DSET_ZERO(sgsTran[0],NX*NY);

#ifdef LAG_PARTICLE_RESOLVED
	vf = data_bag->lag->ng_vfw;
#endif

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

#ifdef LAG_PARTICLE_RESOLVED
				vfrac = min(vf[k][j][i],1.0);
#endif
				//--------------------------------------------------------------
				// Eddy viscosity at boundary points
				//--------------------------------------------------------------
#ifdef LES
				nutE = 0.25*( nut[k][j][i] + nut[k][j][i+1] + nut[k-1][j][i] + nut[k-1][j][i+1] );
	#ifdef XPERIODIC
				nutW = 0.25*( nut[k][j][i] + nut[k][j][i-1] + nut[k-1][j][i] + nut[k-1][j][i-1] );
	#else
				if (i != 0) {
					nutW = 0.25*( nut[k][j][i] + nut[k][j][i-1] + nut[k-1][j][i] + nut[k-1][j][i-1] );
				}
				else {
					nutW = 0.5*( nut[k][j][i] +  nut[k-1][j][i] );
				}
	#endif
				nutN = ( 0.25 * ( nut[k][j][i] + nut[k-1][j][i] + nut[k][j+1][i] + nut[k-1][j+1][i] ) );

				if ( j != 0) {
					nutS = ( 0.25 * ( nut[k][j][i] + nut[k-1][j][i] + nut[k][j-1][i] + nut[k-1][j-1][i] ) );
				}
				else {
	#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
					nutS = 0.0;
	#else
					nutS = ( 0.5 * ( nut[k][j][i] + nut[k-1][j][i] ) );
	#endif
				}

				nutF = nut[k][j][i];
				nutB = nut[k-1][j][i];

#else // not LES
				nutE = 0.0;
				nutW = 0.0;
				nutN = 0.0;
				nutS = 0.0;
				nutF = 0.0;
				nutB = 0.0;
#endif // LES

				//--------------------------------------------------------------
				// d2w/dx2, d2w/dy2 and d2w/dz2
				//--------------------------------------------------------------
				dwdxE = ( w_data[k][j][i+1] - w_data[k][j][i] ) * idx_c[i];

				if (i != 0 ) {
					dwdxW = ( w_data[k][j][i] - w_data[k][j][i-1] ) * idx_c[i-1];
				}
				else {
#ifdef XPERIODIC
					dwdxW = ( w_data[k][j][i] - w_data[k][j][i-1] ) * idx_c[i];
#endif
#ifdef LEFT_WALL_VELOCITY_NOSLIP
					dwdxW = 2*( w_data[k][j][i] ) * idx_u[i];
#endif
#ifdef LEFT_WALL_VELOCITY_FREESLIP
					dwdxW = 0.0;
#endif
				}

				if (j != NY-2) {
					dwdyN = ( w_data[k][j+1][i] - w_data[k][j][i] ) * idy_c[j];
				}
				else {
#ifdef TOP_WALL_VELOCITY_NOSLIP
					dwdyN =  -2.*w_data[k][j][i]*idy_v[j];
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
					dwdyN =  0.0;
#endif
#ifdef TOP_WALL_SCHUMANN
					dwdyN = log_law->dwdy_wm_top[k][i];
#endif
				}

				if ( j != 0) {
					dwdyS = ( w_data[k][j][i] - w_data[k][j-1][i] ) * idy_c[j-1];
				}
				else {
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
					dwdyS =  2.*w_data[k][j][i]*idy_v[j];
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
					dwdyS = 0.;
#endif
#ifdef BOTTOM_WALL_SCHUMANN
					dwdyS = log_law->dwdy_wm_bottom[k][i];
#endif
				}

				dwdzF = ( w_data[k+1][j][i]-w_data[k][j][i] ) * idz_w[k];
				dwdzB = ( w_data[k][j][i]-w_data[k-1][j][i] ) * idz_w[k-1];

				sd2wdx2 = ( nutE * dwdxE - nutW * dwdxW ) * idx_u[i];
				sd2wdy2 = ( nutN * dwdyN - nutS * dwdyS ) * idy_v[j];
				sd2wdz2 = ( nutF * dwdzF - nutB * dwdzB ) * idz_c[k-1];

				vd2wdx2 = iRe * (  dwdxE -  dwdxW ) * idx_u[i];
				vd2wdy2 = iRe * (  dwdyN -  dwdyS ) * idy_v[j];
				vd2wdz2 = iRe * (  dwdzF -  dwdzB ) * idz_c[k-1];

				//--------------------------------------------------------------
				//   d/dx(du/dz), d/dy(dv/dz) and d/dz(dw/dz)
				//--------------------------------------------------------------
				dudzE = ( u_data[k][j][i+1] - u_data[k-1][j][i+1] ) * idz_c[k-1];
				dudzW = ( u_data[k][j][i] - u_data[k-1][j][i] ) * idz_c[k-1];

				dvdzN = ( v_data[k][j+1][i] - v_data[k-1][j+1][i] ) * idz_c[k-1];
				dvdzS = ( v_data[k][j][i]   - v_data[k-1][j][i]   ) * idz_c[k-1];

				sddzdwdz = sd2wdz2;
				sddxdudz = ( nutE * dudzE - nutW * dudzW ) * idx_u[i];
				sddydvdz = ( nutN * dvdzN - nutS * dvdzS ) * idy_v[j];

				vddzdwdz = vd2wdz2;
				vddxdudz = iRe * ( dudzE - dudzW ) * idx_u[i];
				vddydvdz = iRe * ( dvdzN - dvdzS ) * idy_v[j];

				viscDiss[j][i]  += w_data[k][j][i] * (vd2wdx2 + vd2wdy2 + vd2wdz2) * (1.0-vfrac);
				viscDiss[j][i]  += w_data[k][j][i] * (vddzdwdz + vddxdudz + vddydvdz) * (1.0-vfrac);
#ifdef LES
				sgsDiss[j][i]  += w_data[k][j][i] * (sd2wdx2 + sd2wdy2 + sd2wdz2) * (1.0-vfrac) ;
				sgsDiss[j][i]  += w_data[k][j][i] * (sddzdwdz + sddxdudz + sddydvdz) * (1.0-vfrac) ;
#endif

				//--------------------------------------------------------------
				// Calculate the transport terms
				//--------------------------------------------------------------
				wE = 0.5*(w_data[k][j][i] + w_data[k][j][i+1]);
#ifdef XPERIODIC
				wW = 0.5*(w_data[k][j][i] + w_data[k][j][i-1]);
#else
				if (i !=0 ) {
					wW = 0.5*(w_data[k][j][i] + w_data[k][j][i-1]);
				}
				else {
	#ifdef LEFT_WALL_VELOCITY_FREESLIP
					wW = w_data[k][j][i];
	#else
					wW = 0.0;
	#endif
				}
#endif // XPERIODIC

				wN = 0.5*(w_data[k][j][i] + w_data[k][j+1][i]);
				if (j != 0) {
					wS = 0.5*(w_data[k][j][i] + w_data[k][j-1][i]);
				}
				else {
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
					wS = 0.0;
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
					wS = w_data[k][j][i];
#endif
				}
				wF = w_data_bc[k][j][i];
				wB = w_data_bc[k-1][j][i];

				vtd2wdx2 = iRe * ( wE* dwdxE - wW* dwdxW ) * idx_u[i];
				vtd2wdy2 = iRe * ( wN* dwdyN - wS* dwdyS ) * idy_v[j];
				vtd2wdz2 = iRe * ( wF* dwdzF - wB* dwdzB ) * idz_c[k-1];

				vtddzdwdz = vtd2wdz2;
				vtddxdudz = iRe * ( wE* dudzE - wW* dudzW ) * idx_u[i];
				vtddydvdz = iRe * ( wN* dvdzN - wS* dvdzS ) * idy_v[j];

				viscTran[j][i]  += (vtd2wdx2 + vtd2wdy2 + vtd2wdz2) * (1.0-vfrac);
				viscTran[j][i]  += (vtddzdwdz + vtddxdudz + vtddydvdz) * (1.0-vfrac);
#ifdef LES
				std2wdx2 = ( wE*nutE * dwdxE - wW*nutW * dwdxW ) * idx_u[i];
				std2wdy2 = ( wN*nutN * dwdyN - wS*nutS * dwdyS ) * idy_v[j];
				std2wdz2 = ( wF*nutF * dwdzF - wB*nutB * dwdzB ) * idz_c[k-1];

				stddzdwdz = std2wdz2;
				stddxdudz = ( wE*nutE * dudzE - wW*nutW * dudzW ) * idx_u[i];
				stddydvdz = ( wN*nutN * dvdzN - wS*nutS * dvdzS ) * idy_v[j];

				sgsTran[j][i]  += (std2wdx2 + std2wdy2 + std2wdz2) * (1.0-vfrac);
				sgsTran[j][i]  += (stddzdwdz + stddxdudz + stddydvdz) * (1.0-vfrac);
#endif

			}
		}
	}

	scale = 1.0/(NZ-1.0);
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			viscDiss[j][i] *= scale;
			viscTran[j][i] *= scale;
			sgsDiss[j][i]  *= scale;
			sgsTran[j][i]  *= scale;
			viscDiss[j][i] = viscTran[j][i] - viscDiss[j][i];
			sgsDiss[j][i]  = sgsTran[j][i]  - sgsDiss[j][i];
		}
	}

	if (Ie == NX-1) {
		for (j=grid->G_Js; j<Je; j++) {
			viscDiss[j][NX-1] = viscDiss[j][NX-2];
			viscTran[j][NX-1] = viscTran[j][NX-2];
			sgsDiss[j][NX-1] = sgsDiss[j][NX-2];
			sgsTran[j][NX-1] = sgsTran[j][NX-2];
		}
	}

	if (Je == NY-1) {
		for (i=grid->G_Is; i<Ie; i++) {
			viscDiss[NY-1][i] = viscDiss[NY-2][i];
			viscTran[NY-1][i] = viscTran[NY-2][i];
			sgsDiss[NY-1][i] = sgsDiss[NY-2][i];
			sgsTran[NY-1][i] = sgsTran[NY-2][i];
		}
	}
	return;
}

#ifdef LAG_PARTICLE_RESOLVED


void Statistics2d_vf_avg(double **vf_avg, Cart3d_bag *data_bag){
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double ***vf;
	double vfrac = 0.0;

	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);

	// Set global mean to 0 before MPI_Allreduce is called
	DSET_ZERO(vf_avg[0],NX*NY);

	vf = data_bag->lag->ng_vfc;


	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
				vfrac = min(vf[k][j][i],1.0);
				vf_avg[j][i] += (1.0-vfrac);
			}
		}
	}
	for(j = Js; j < Je; j++) {
		for(i = Is; i < Ie; i++) {
			vf_avg[j][i] = vf_avg[j][i]/(NZ-1.);
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &vf_avg[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);

	return;
}


void Statistics2d_vf_slice(double ***quantity, double **global_slice, Cart3d_bag *data_bag) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int NX, NY, NZ;
	double vfrac = 0.0;


	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	int N_slice = round(NZ/2);

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = min(NZ-1, grid->G_Ke);
	//Ke = grid->G_Ke;

	//#ifdef XPERIODIC
	//	Is = max(1,Is);
	//	if (Ie == NX-1)
	//		Ie = NX;
	//#endif

	//#ifdef ZPERIODIC
	//	Ks = max(1,Ks);
	//	if (Ke == NZ-1)
	//		Ke = NZ;
	//#endif

	DSET_ZERO(global_slice[0],NX*NY);

	if ( ( N_slice >= Ks) && (N_slice < Ke) ){
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				global_slice[j][i] = min(1.0,quantity[N_slice][j][i]);
			}
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &global_slice[0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);

	return;
}


#endif


void Statistics2d_globalmean2d(double **mean, double **global_mean, MAC_grid *grid) {

	int i, j, k;
	int NX, NY, NZ;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	for(j = 0; j < NY; j++) {
		MPI_Allreduce (mean[j], global_mean[j], NX, MPI_DOUBLE, MPI_SUM, PCW);
	}
	return;

}





/**********************************************************************/

#ifdef CONC
void Statistics2d_compute_nusselt( Cart3d_bag *data_bag) {


	int i, j, k, iconc;
	double dcdyS;
	double vS;
	double cS;
	double vcS;

	double rhs;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Statistics2d *st2d = data_bag ->st2d;

	int Nconc = params->NConc;
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

	double ***v_data = data_bag -> v -> data;


#ifdef VOF_SCALAR

	double ***v_vof= data_bag->v->data_vof;
	double ***vfv = data_bag->lag->vfv;
	double ***ng_vfc = data_bag->lag->ng_vfc;

#endif



double Pe;
double ih= grid -> idx_u[i_start];
double ihsq= ih*ih;

#ifdef VOF_SCALAR
	double conductivity_s;
	double heat_cap_s;
#endif

double ***c_data;
double lambdaS=1;


for (iconc=0; iconc < Nconc; iconc++ ){
double Pe= params->Pe[iconc];

#ifdef VOF_SCALAR
	 conductivity_s = params->conductivity_s[iconc];
	 heat_cap_s = params->vol_heat_cap_s[iconc];
#endif


	c_data = data_bag->c[iconc]-> data;
	DSET_ZERO(&st2d->nusselt[iconc][0][0], (NX)*(NY));
	st2d->integral_nusselt[iconc]=0;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {


#ifdef VOF_SCALAR
				lambdaS = vfv[k][j][i]*(conductivity_s-1)+1;
#endif
				dcdyS = ( c_data[k][j][i] - c_data[k][j-1][i] )*ih ;


#ifdef VOF_VELOCITY

				vS = v_vof[k][j][i];
				vS = v_data[k][j][i];
#else

				vS = v_data[k][j][i];
#endif

				cS = 0.5 * ( c_data[k][j][i] + c_data[k][j-1][i] );

				vcS = vS*cS;

#ifdef VOF_SCALAR

				st2d->nusselt[iconc][j][i] += vcS*( vfv[k][j][i]*(heat_cap_s-1)+1)*Pe-dcdyS*lambdaS;
#else
				st2d->nusselt[iconc][j][i] += vcS-dcdyS*lambdaS;
#endif
			}
		}
	}



	MPI_Allreduce(MPI_IN_PLACE, &st2d->nusselt[iconc][0][0], NX*NY, MPI_DOUBLE, MPI_SUM, PCW);


	for(j = 0; j < NY-1; j++) {
		for(i = 0; i < NX-1; i++) {
			st2d->integral_nusselt[iconc] += st2d->nusselt[iconc][j][i];
			st2d->nusselt[iconc][j][i] *= 1./(NZ-1);
		}
	}

	st2d->integral_nusselt[iconc] *= 1./((NX-1)*(NY-1)*(NZ-1));

}  // For iconc=0..Nconc-1

	return;

}

#endif

void Statistics2d_fill_ghostnodes_2d( double **quantity, MAC_grid *grid) {


	int i,j;
	int NX = grid->NX;
	int NY = grid->NY;


	int Is = grid->G_Is;
	int Js = grid->G_Js;

	j=NY-1;
		for(i = 0; i < NX-1; i++) {
			quantity[j][i] = quantity[j-1][i];
		}

	i=NX-1;
		for(j = 0; i < NY-1; j++) {
			quantity[j][i] = quantity[j][i-1];
		}


	return;
}

void Statistics2d_covarcalc(double *avg_velc, double ***c_data, double ***vel_data, double *avg_mvf, double *avg_c, double *avg_vel, Cart3d_bag *data_bag){

	int i,j,k;
	int l;
	int Navg = 0;
	int Ncells = 0;
	double cp,up,vf;

	double ***vfc = data_bag->lag->ng_vfc;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

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

	if (params->avg_dir == 0){
		Navg = NX-1;
		Ncells = (NY-1) * (NZ-1);
	}
	else{
		Navg = NY-1;
		Ncells = (NX-1) * (NZ-1);
	}

	DSET_ZERO(avg_velc, Navg);
	DSET_ZERO(avg_mvf, Navg);
	DSET_ZERO(avg_vel, Navg);
	DSET_ZERO(avg_c, Navg);

	// Calculate local fluid flux and volume fraction data

		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {

					if(params->avg_dir==0){
						l=i;
					}
					else{
						l=j;
					}

					#ifdef LAG_PARTICLE_RESOLVED
					vf = min(vfc[k][j][i], 1.0);
					avg_c[l] += (1.0 - vf) * c_data[k][j][i];
					avg_vel[l] += (1.0 - vf) * vel_data[k][j][i];
					avg_mvf[l] +=  (1.0 - vf);
					#else
					avg_c[l] += c_data[k][j][i];
					avg_vel[l] += vel_data[k][j][i];
					avg_mvf[l] += 1.0;
					#endif
				}
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, avg_mvf, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_c, Navg, MPI_DOUBLE, MPI_SUM, PCW);
		MPI_Allreduce(MPI_IN_PLACE, avg_vel, Navg, MPI_DOUBLE, MPI_SUM, PCW);

		for (j = 0; j < Navg; j++) {
			avg_mvf[j] = avg_mvf[j]/Ncells;
			avg_c[j] = avg_c[j]/avg_mvf[j]/Ncells;
			avg_vel[j] = avg_vel[j]/avg_mvf[j]/Ncells;
		}


		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {

					if(params->avg_dir==0){
						l=i;
					}
					else{
						l=j;
					}

					#ifdef LAG_PARTICLE_RESOLVED
					vf = min(vfc[k][j][i], 1.0);
					cp = (1.0-vf)*(c_data[k][j][i] -  avg_c[l]);
					up = (1.0-vf)*(vel_data[k][j][i] -  avg_vel[l]);
					#else
					cp = c_data[k][j][i] - avg_c[l];
					up = vel_data[k][j][i] - avg_vel[l];
					#endif
					avg_velc[l] += up*cp;

				}
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, avg_velc, Navg, MPI_DOUBLE, MPI_SUM, PCW);

		for (j = 0; j < Navg; j++) {
			avg_velc[j] = avg_velc[j]/avg_mvf[j]/Ncells;
		}
}

void Statistics2d_histocalc(double ***c_data, double *histo, double *range, Cart3d_bag *data_bag){

int ibin = 0;
int NTOT = 0;
double eps = 0.000001;
double cbd;
double vf;
int i, j, k;

MAC_grid *grid = data_bag -> grid;
Parameters *params = data_bag -> params;


int Nbin = params->Nbin;
double cbd0 = params->cbd0;
double cbd1 = params->cbd1;
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



#ifdef LAG_PARTICLE_RESOLVED
double ***vfc = data_bag->lag->ng_vfc;
#endif
//double *histo = (double *)malloc(Nbin * sizeof(double));
//double *range = (double *)malloc(Nbin * sizeof(double));

	for (ibin=0; ibin<Nbin; ibin++) {
		range[ibin] = cbd0 + (cbd1-cbd0)*((double) ibin) / Nbin;
		histo[ibin] = 0.0;
	}

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++) {
				NTOT++;
#ifdef LAG_PARTICLE_RESOLVED
				vf = min(vfc[k][j][i], 1.0);
				cbd = ((1-vf)*c_data[k][j][i] - cbd0)/(cbd1-cbd0);
#else
				cbd =(c_data[k][j][i] - cbd0)/(cbd1-cbd0);
#endif

				cbd = min(cbd, 1.0-eps);
				cbd = max(cbd, eps);
				ibin = (int) (cbd * Nbin);
				if ( (ibin>=0) && (ibin<Nbin) ) histo[ibin]++;
			}
		}
	}


	MPI_Allreduce(MPI_IN_PLACE, &NTOT, 1, MPI_INT, MPI_SUM, PCW);
	MPI_Allreduce(MPI_IN_PLACE, &histo[0], Nbin, MPI_DOUBLE, MPI_SUM, PCW);

	for (ibin=0; ibin<Nbin; ibin++) {
		histo[ibin]*= 1.0/(double)NTOT;
	}

}


/******************************************************************************/
/*
 Write Data to H5, assuming that Data2d_XX.h5 is closed only MASTER rank is
 writing
 */
/******************************************************************************/
void Statistics2d_writeStatistics2dToH5(Statistics2d *st2d, MAC_grid *grid,
		Parameters *params, Cart3d_bag *data_bag ,int timestep, Debug_trace *dtrace) {

	int NX = grid->NX;
	int NY = grid->NY;
	int NZ = grid->NZ;
	int NZslice = 1;
	int iconc, NConc;
	char h5filename2d[50];
	char fieldname[50];
	char groupname[50];
	int rank=2; // rank of variables saved into h5
	int ierr;
	int i, j ;
	int Nbin = params->Nbin;
	int Navg;
	if (params->avg_dir == 0){
		Navg = NX-1;
	}
	else{
		Navg = NY-1;
	}
	hid_t file_id_2d; // file handles
	hid_t group_id; // group handle
	hid_t datatype, dataspace, dataset;   /* handles */
	hsize_t dim_2d[2] ={NY, NX}; // dimensions of data to be written
	hsize_t dim_2d_ng[2] ={NY, NX }; // dimensions of data to be written

	hsize_t dim_x[1] = {NX};
	hsize_t dim_y[1] = {NY};
	hsize_t integral[1] = {1};

	hsize_t histosize[1] = {Nbin};
	hsize_t covarsize[1] = {Navg};

	herr_t status; // status variable to check return value of HDF5 routines
	double **work1 = st2d->work1;
	double ***u_cell= data_bag->u->data_bc;
	double ***v_cell = data_bag->v->data_bc;
	double ***w_cell = data_bag->w->data_bc;

	NConc = params->NConc;

	if(params->rank==0){
	printf("Write 2d Statistics\n");
	}
	// Open already existing h5-file in Read-Write mode
	sprintf(h5filename2d, "Data2d_%d.h5", abs(timestep));
	file_id_2d = Output_h5_open(h5filename2d, params, DTRACE("Output_h5_open"));

	// Describe the size of the array and create the data space for fixed size
	// dataset.


	sprintf(fieldname, "/time");

	Output_h5_dataset(&(params->time), H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(groupname, "/grid");
	Output_h5_create_group(file_id_2d, groupname, params, DTRACE("Output_h5_create_group"));

	sprintf(fieldname, "%s/x", groupname);
	Output_h5_dataset(grid->xc, H5T_NATIVE_DOUBLE, 1, dim_x, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "%s/y", groupname);
	Output_h5_dataset(grid->yc, H5T_NATIVE_DOUBLE, 1, dim_y, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


	sprintf(fieldname, "%s/z", groupname);
	Output_h5_dataset(&(grid->zc[NZ/2]), H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


	sprintf(fieldname, "%s/NX", groupname);
	Output_h5_dataset(dim_x, H5T_NATIVE_INT, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/NY", groupname);
	Output_h5_dataset(dim_y, H5T_NATIVE_INT, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "%s/NZ", groupname);
	Output_h5_dataset(&NZslice, H5T_NATIVE_INT, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


	// conc
		sprintf(groupname, "/integral");
		Output_h5_create_group(file_id_2d, groupname, params, DTRACE("Output_h5_create_group"));
		sprintf(groupname, "/histogram");
		Output_h5_create_group(file_id_2d, groupname, params, DTRACE("Output_h5_create_group"));
		sprintf(groupname, "/xdim");  // put quantities that are function of x
		Output_h5_create_group(file_id_2d, groupname, params, DTRACE("Output_h5_create_group"));
		sprintf(groupname, "/covariance");  // put quantities that are function of x
		Output_h5_create_group(file_id_2d, groupname, params, DTRACE("Output_h5_create_group"));

#ifdef LAG_PARTICLE_RESOLVED
	sprintf(fieldname, "/covariance/vf");
	Output_h5_dataset(&st2d->avg_vf[0], H5T_NATIVE_DOUBLE, 1, covarsize, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
#endif

#ifdef CONC
	for (iconc=0; iconc<NConc; iconc++) {
		sprintf(fieldname, "/c%dMean", iconc);
		Output_h5_dataset(st2d->mean_c[iconc][0], H5T_NATIVE_DOUBLE, 2, dim_2d,
		                  file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

  	    sprintf(fieldname, "/c%dSlice", iconc);
  		Output_h5_dataset(st2d->slice_c[iconc][0], H5T_NATIVE_DOUBLE, 2, dim_2d,
  		                  file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	  if(params->nusselt_switch){
  		sprintf(fieldname, "/nusselt_c%d", iconc);
  		Output_h5_dataset(&st2d->nusselt[iconc][0][0], H5T_NATIVE_DOUBLE, 2, dim_2d,file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		}
		sprintf(fieldname, "/grad_y_c%d", iconc);
		Output_h5_dataset(&st2d->grad_y_c[iconc][0][0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

  		//if(params->current_height_switch == 1 && iconc == 0){
  		//	sprintf(fieldname, "/xdim/Current_height_full_c%d", iconc);
  		//	Output_h5_dataset(st2d->current_height[iconc], H5T_NATIVE_DOUBLE, 1, dim_x, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

  		//	sprintf(fieldname, "/xdim/Current_height_fluid_c%d", iconc);
  		//	Output_h5_dataset(st2d->current_height_fluid[iconc], H5T_NATIVE_DOUBLE, 1, dim_x, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

  			//sprintf(fieldname, "/integral/Front_location_c%d", iconc);
  		  	//Output_h5_dataset(&st2d->front_position[iconc], H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


  		//}


		    sprintf(fieldname, "/integral/integral_conc_fluid_c%d", iconc);
		  	Output_h5_dataset(&st2d->integral_conc_fluid[iconc], H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

		  	sprintf(fieldname, "/integral/integral_conc_full_c%d", iconc);
		  	Output_h5_dataset(&st2d->integral_conc_full[iconc], H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


		  	sprintf(fieldname, "/integral/integral_potEnergy_fluid_c%d", iconc);
		  	Output_h5_dataset(&st2d->integral_potEnergy_fluid[iconc], H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

			sprintf(fieldname, "/histogram/histo_%d", iconc);
			Output_h5_dataset(&st2d->histo[iconc][0], H5T_NATIVE_DOUBLE, 1, histosize, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

			sprintf(fieldname, "/histogram/range_%d", iconc);
			Output_h5_dataset(&st2d->range[iconc][0], H5T_NATIVE_DOUBLE, 1, histosize, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

			sprintf(fieldname, "/covariance/uc%d", iconc);
			Output_h5_dataset(&st2d->avg_uc[iconc][0], H5T_NATIVE_DOUBLE, 1, covarsize, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

			sprintf(fieldname, "/covariance/vc%d", iconc);
			Output_h5_dataset(&st2d->avg_vc[iconc][0], H5T_NATIVE_DOUBLE, 1, covarsize, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

			sprintf(fieldname, "/covariance/wc%d", iconc);
			Output_h5_dataset(&st2d->avg_wc[iconc][0], H5T_NATIVE_DOUBLE, 1, covarsize, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

			sprintf(fieldname, "/covariance/c%d", iconc);
			Output_h5_dataset(&st2d->avg_c[iconc][0], H5T_NATIVE_DOUBLE, 1, covarsize, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


  		if(params->nusselt_switch){
  		sprintf(fieldname, "/integral/integral_nusselt_c%d", iconc);
  		Output_h5_dataset(&st2d->integral_nusselt[iconc], H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

  		}

	}
#endif

#ifdef LAG_PARTICLE_RESOLVED
	// vfc
	sprintf(fieldname, "/vf_avg");
	Output_h5_dataset(st2d->vf_avg[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	Statistics2d_horizontalSlice2d(data_bag->lag->ng_vfc, work1, 'c', data_bag);
	sprintf(fieldname, "/vf_slice");
	Output_h5_dataset(work1[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
#endif

	// u
	sprintf(fieldname, "/uMean");
	Output_h5_dataset(st2d->mean_u[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	// v
	sprintf(fieldname, "/vMean");
	Output_h5_dataset(st2d->mean_v[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	// w
	sprintf(fieldname, "/wMean");
	Output_h5_dataset(st2d->mean_w[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	// p
	sprintf(fieldname, "/pMean");
	Output_h5_dataset(st2d->mean_p[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	Statistics2d_horizontalSlice2d(u_cell, work1, 'c', data_bag);
	sprintf(fieldname, "/uslice");
	Output_h5_dataset(work1[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	Statistics2d_horizontalSlice2d(v_cell, work1, 'c', data_bag);
	sprintf(fieldname, "/vslice");
	Output_h5_dataset(work1[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


	if (params-> turbulent_terms_switch == 1) {
		// turbStress_uu and flux
		sprintf(fieldname, "/turbStress_uu");
		Output_h5_dataset(st2d->turbStress_uu[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		// turbStress_vv and flux_vv
		sprintf(fieldname, "/turbStress_vv");
		Output_h5_dataset(st2d->turbStress_vv[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		// turbStress_ww and fluxes
		sprintf(fieldname, "/turbStress_ww");
		Output_h5_dataset(st2d->turbStress_ww[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		// turbStress_uv and fluxes
		sprintf(fieldname, "/turbStress_uv");
		Output_h5_dataset(st2d->turbStress_uv[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		// turbStress_uw and fluxes
		sprintf(fieldname, "/turbStress_uw");
		Output_h5_dataset(st2d->turbStress_uw[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		// turbStress_vw and fluxes
		sprintf(fieldname, "/turbStress_vw");
		Output_h5_dataset(st2d->turbStress_vw[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

		sprintf(fieldname, "/turbKinE");
		Output_h5_dataset(st2d->turbKinE[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	}

	sprintf(fieldname, "/flux_uu");
	Output_h5_dataset(st2d->flux_uu[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "/flux_vv");
	Output_h5_dataset(st2d->flux_vv[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "/flux_ww");
	Output_h5_dataset(st2d->flux_ww[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "/flux_uv");
	Output_h5_dataset(st2d->flux_uv[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "/flux_uw");
	Output_h5_dataset(st2d->flux_uw[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	sprintf(fieldname, "/flux_vw");
	Output_h5_dataset(st2d->flux_vw[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	// turbStress_xcx and fluxes
#ifdef CONC
	int indexC=0;
	for (iconc=0; iconc<NConc; iconc++) {
		indexC=iconc*4;

		if (params-> turbulent_terms_switch == 1) {
			sprintf(fieldname, "/turbStress_uc%d", iconc);
			Output_h5_dataset(st2d->turbStressC[indexC+0][0], H5T_NATIVE_DOUBLE, 2,
						  dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
			sprintf(fieldname, "/turbStress_vc%d", iconc);
			Output_h5_dataset(st2d->turbStressC[indexC+1][0], H5T_NATIVE_DOUBLE, 2,
						  dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
			sprintf(fieldname, "/turbStress_c%dc%d", iconc,iconc);
			Output_h5_dataset(st2d->turbStressC[indexC+3][0], H5T_NATIVE_DOUBLE, 2,
						  dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
			sprintf(fieldname, "/turbStress_wc%d", iconc);
			Output_h5_dataset(st2d->turbStressC[indexC+2][0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

		}
		sprintf(fieldname, "/flux_uc%d", iconc);
		Output_h5_dataset(st2d->fluxC[indexC+0][0], H5T_NATIVE_DOUBLE, 2,dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "/flux_vc%d", iconc);
		Output_h5_dataset(st2d->fluxC[indexC+1][0], H5T_NATIVE_DOUBLE, 2,
								  dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "/flux_wc%d", iconc);
		Output_h5_dataset(st2d->fluxC[indexC+2][0], H5T_NATIVE_DOUBLE, 2,
								  dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
		sprintf(fieldname, "/flux_c%dc%d", iconc,iconc);
		Output_h5_dataset(st2d->fluxC[indexC+3][0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

		sprintf(fieldname, "/integral/integral_buoyantwork_c%d",iconc);
		Output_h5_dataset(&st2d->integral_buoyantwork[iconc], H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	}

#endif

	// kinE totKinE
	sprintf(fieldname, "/totKinE");
	Output_h5_dataset(st2d->kinE[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "/Urms");
	Output_h5_dataset(st2d->Urms[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

#ifdef LAG_PARTICLE_RESOLVED
	sprintf(fieldname, "/totPotEnergy_p");
	Output_h5_dataset(st2d->potEnergyP[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "/integral/integral_potEnergy_part");
	Output_h5_dataset(&st2d->integral_potEnergy_part, H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "/integral/integral_kinEnergy_part");
	Output_h5_dataset(&st2d->integral_kinEnergy_part, H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
#endif

	sprintf(fieldname, "/integral/integral_kinEnergy_fluid");
	Output_h5_dataset(&st2d->integral_kinEnergy_fluid, H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "/integral/integral_RMS_fluid");
	Output_h5_dataset(&st2d->integral_RMS_fluid, H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "/integral/integral_forcing_fluid");
	Output_h5_dataset(&st2d->integral_forcing_fluid, H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	sprintf(fieldname, "/integral/integral_viscDiss");
	Output_h5_dataset(&st2d->integral_viscDiss, H5T_NATIVE_DOUBLE, 1, integral, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));



	double **potEwrite;
	// potEnergy of C

	for (iconc = 0; iconc < NConc; iconc++) {

			potEwrite=st2d->potEnergyC[iconc];
			sprintf(fieldname, "/totPotEnergy_c%d", iconc);
		    //Statistics2d_fill_ghostnodes_2d(st2d->potEnergyC[iconc], grid );
			Output_h5_dataset(&(st2d->potEnergyC[iconc][0][0]), H5T_NATIVE_DOUBLE, 2, dim_2d_ng, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


	}

	// viscous dissipation
	sprintf(fieldname, "/viscDiss");
	Output_h5_dataset(st2d->viscDiss[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


	// viscous transport
	sprintf(fieldname, "/viscTran");
	Output_h5_dataset(st2d->viscTran[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	// viscous transport
	sprintf(fieldname, "/pwork");
	Output_h5_dataset(st2d->pwork[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

#ifdef LES
	// SGS dissipation
	sprintf(fieldname, "/sgsDiss");
	Output_h5_dataset(st2d->sgsDiss[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	// SGS transport
	sprintf(fieldname, "/sgsTran");
	Output_h5_dataset(st2d->sgsTran[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));


	// SGS viscosity
	sprintf(fieldname, "/nut");
	Output_h5_dataset(st2d->mean_nut[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	#ifdef SMAG_DYNAMIC
	// Cev
	sprintf(fieldname, "/cev");
	Output_h5_dataset(st2d->mean_cev[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	#endif

	#ifdef CONC_DYNAMIC
	// Sct
	sprintf(fieldname, "/sct");
	Output_h5_dataset(st2d->mean_sct[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	#endif

#endif

#ifdef RANS
	// RANS viscosity
	sprintf(fieldname, "/nut");
	Output_h5_dataset(st2d->mean_nut[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	#ifdef TWO_EQUATION_MODEL
	// TKE modeled
	sprintf(fieldname, "/tke",);
	Output_h5_dataset(st2d->mean_tke[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));

	// epsilon modeled
	sprintf(fieldname, "/eps");
	Output_h5_dataset(st2d->mean_eps[0], H5T_NATIVE_DOUBLE, 2, dim_2d, file_id_2d, fieldname, params, DTRACE("Output_h5_dataset"));
	#endif

#endif



	//H5Gclose(group_id);
	H5Fclose(file_id_2d);

	return;

}
