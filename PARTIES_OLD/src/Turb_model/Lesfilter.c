#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Conc.h"
#include "Cart3d.h"
#include "Velocity.h"
#include "Subgrid.h"
#include "Lesfilter.h" 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/******************************************************************************/
/*
 This function allocates enough memory for the concentration structure based on
 parameters defined in "*params"
 */
/******************************************************************************/
void Lesfilter_init(Subgrid *smag, MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ;
	int ierr;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end; 
	int j_start, j_end; 
	int k_start, k_end;

	double *wxp,*wxe,*wxw;
	double *wyp,*wys,*wyn;
	double *wzp,*wzb,*wzf;
	double dx, dy, dz; 
	double dxw, dxe, dxt, dxwu, dxwwu; 
	double dys, dyn, dyt, dysv, dynnv; 
	double dzb, dzf, dzt, dzbw, dzffw; 

	FILE *fid;


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
	i_start = Is; // i=0, i=NX-1 are not included
	j_start = Js;
	k_start = Ks; 
		
	// exclude the half cell added
	i_end   = min(NX-1, Ie); 
	j_end   = min(NY-1, Je); 
	k_end   = min(NZ-1, Ke); 

	wxp = (double *) calloc(NX, sizeof(double));
	wxe = (double *) calloc(NX, sizeof(double));
	wxw = (double *) calloc(NX, sizeof(double));
	Memory_check_allocation(wxp);
	Memory_check_allocation(wxe);
	Memory_check_allocation(wxw);

	wyp = (double *) calloc(NY, sizeof(double));
	wyn = (double *) calloc(NY, sizeof(double));
	wys = (double *) calloc(NY, sizeof(double));
	Memory_check_allocation(wyp);
	Memory_check_allocation(wyn);
	Memory_check_allocation(wys);

	wzp = (double *) calloc(NZ, sizeof(double));
	wzb = (double *) calloc(NZ, sizeof(double));
	wzf = (double *) calloc(NZ, sizeof(double));
	Memory_check_allocation(wzp);
	Memory_check_allocation(wzb);
	Memory_check_allocation(wzf);

	for (i=0;i<NX;i++){
		wxp[i] = 1.0;
		wxe[i] = 0.0;
		wxw[i] = 0.0;
	}

	for (i=1;i<NX-1;i++){
		dxw = grid->xc[i]-grid->xc[i-1];
		dxe = grid->xc[i+1]-grid->xc[i];
		dxt = grid->xc[i+1]-grid->xc[i-1];
		wxe[i] = 0.5*dxe/dxt;
		wxw[i] = 0.5*dxw/dxt;
		wxp[i] = 0.5; 
	}

	// Near the inlet
#ifndef XPERIODIC
	i=0;
	dxe = grid->xc[i+1]-grid->xc[i];
	dxwu = grid->xc[i]-grid->xu[i];
	dxwwu = grid->xu[i+2]-grid->xc[i+1];
	dxt = dxe + dxwu + dxwwu;
	wxp[i] =  (dxwu + 0.5*dxe)/dxt;
	wxe[i] =  0.5*(dxe+dxwwu)/dxt;
	wxw[i] = 0.5*dxwwu/dxt;
#else
	i=0;
	wxp[i]=0.5;
	wxe[i]=0.25;
	wxw[i]=0.25;
#endif 

	// Near the Outlet
/*
	i=NX-1;
	dxe = grid->xc[i+1]-grid->xc[i];
	dxwu = grid->xc[i]-grid->xu[i];
	dxwwu = grid->xu[i+2]-grid->xc[i+1];
	dxt = dxe + dxwu + dxeu;
	wxp[i] =  (dxwu + 0.5*dxe)/dxt;
	wxe[i] =  0.5*(dxe+dxwwu)/dxt;
	wxw[i] = 0.5*dxwwu/dxt;
*/
	for (j=0;j<NY;j++){
		wyp[j] = 1.0;
		wyn[j] = 0.0;
		wys[j] = 0.0;
	}

	for (j=1;j<NY-1;j++){
		dys = grid->yc[j]-grid->yc[j-1];
		dyn = grid->yc[j+1]-grid->yc[j];
		dyt = grid->yc[j+1]-grid->yc[j-1];
		wyn[j] = 0.5*dyn/dyt;
		wys[j] = 0.5*dys/dyt;
		wyp[j] = 0.5; 
	}
	// Near the Bottom boundary  
	j=0;
	dyn = grid->yc[j+1]-grid->yc[j];
	dysv = grid->yc[j]-grid->yv[j];
	dynnv = grid->yv[j+2]-grid->yc[j+1];
	dyt = dyn + dysv + dynnv;
	wyp[j] =  (dysv + 0.5*dyn)/dyt;
	wyn[j] =  0.5*(dyn+dynnv)/dyt;
	wys[j] = 0.5*dynnv/dyt;

	dyn = grid->yc[1]-grid->yc[0];
	dynnv = grid->yc[2]-grid->yc[0];
	wyp[j] =  0.5*(dyn + dynnv)/(dynnv);
	wyn[j] =  0.5;
	wys[j] = -0.5*dyn/dynnv;

	j=NY-2;
	dys = grid->yc[j]-grid->yc[j-1];
	dysv = grid->yc[j]-grid->yc[j-2];
	wyp[j] = 0.5*(dys+dysv)/dysv;
	wys[j] = 0.5;
	wyn[j] = -0.5*dys/dysv;

/*
	wyp[0]=0.5;
	wyn[0]=0.25;
	wys[0]=0.25;
*/
	for (k=0;k<NZ;k++){
		wzp[k] = 1.0;
		wzb[k] = 0.0;
		wzf[k] = 0.0;
	}

	for (k=1;k<NZ-1;k++){
		dzb = grid->zc[k]-grid->zc[k-1];
		dzf = grid->zc[k+1]-grid->zc[k];
		dzt = grid->zc[k+1]-grid->zc[k-1];
		wzb[k] = 0.5*dzb/dzt;
		wzf[k] = 0.5*dzf/dzt;
		wzp[k] = 0.5; 
	}
	// Near the Bottom boundary  
#ifndef ZPERIODIC
	k=0;
	dzf = grid->zc[k+1]-grid->zc[k];
	dzbw = grid->zc[k]-grid->zw[k];
	dzffw = grid->zw[k+2]-grid->zc[k+1];
	dzt = dzf + dzbw + dzffw;
	wzp[k] =  (dzbw + 0.5*dzf)/dzt;
	wzf[k] =  0.5*(dzf+dzffw)/dzt;
	wzb[k] = 0.5*dzffw/dzt;
#else
	k=0;
	wzp[k]=0.5;
	wzf[k]=0.25;
	wzb[k]=0.25;
#endif

	smag->wxw = wxw;
	smag->wxp = wxp;
	smag->wxe = wxe;
	smag->wys = wys;
	smag->wyp = wyp;
	smag->wyn = wyn;
	smag->wzb = wzb;
	smag->wzp = wzp;
	smag->wzf = wzf;

	if (params->rank==0) {

		fid = fopen("filter_info.dat","w");
		for (i=0;i<NX;i++){
			fprintf(fid,"wxw = %f wxp = %f wxe = %f \n",wxw[i],wxp[i],wxe[i]);
		}
		fprintf(fid,"\n \n \n");
		for (j=0;j<NY;j++){
			fprintf(fid,"wys = %f wyp = %f wyn = %f \n",wys[j],wyp[j],wyn[j]);
		}
		fprintf(fid,"\n \n \n");
		for (k=0;k<NZ;k++){
			fprintf(fid,"wzb = %f wzp = %f wzf = %f \n",wzb[k],wzp[k],wzf[k]);
		}
		fprintf(fid,"\n \n \n");

		fclose(fid); 
	}

	return;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void lesfilter(double ***u, double ***uf, double ***uf1, Cart3d_bag *data_bag) {

	int NX, NY, NZ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end; 
	int j_start, j_end; 
	int k_start, k_end;
	int ierr; 
	int i, j, k;
	double *wxp,*wxe,*wxw;
	double *wyp,*wys,*wyn;
	double *wzp,*wzb,*wzf;
//	double ***uk;
	double Tstart, Tend;
	
	Tstart = MPI_Wtime();
	
	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Subgrid *smag = data_bag -> smag;
	
	Communication_update_ghost_nodes_x(u, 'c', 1, data_bag);
//	Communication_update_ghost_nodes_flow_variable(u, 'c', 1, data_bag);

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

	i_start = Is; 
	j_start = Js; 
	k_start = Ks; 
	// indices start and end on current processor
#ifndef XPERIODIC
	i_start = max(1,Is); // i=0, i=NX-1 are not included
#endif
	j_start = max(1,Js);
#ifndef ZPERIODIC
	k_start = max(1,Ks); 
#endif
		
	// exclude the half cell added
	i_end   = min(NX-1, Ie); 
	j_end   = min(NY-1, Je); 
	k_end   = min(NZ-1, Ke); 

	wxw = smag->wxw ;
	wxp = smag->wxp ;
	wxe = smag->wxe ;
	wys = smag->wys ;
	wyp = smag->wyp ;
	wyn = smag->wyn ;
	wzb = smag->wzb ;
	wzp = smag->wzp ;
	wzf = smag->wzf ;
	
	// Filter along x-direction
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=i_start; i<i_end; i++){
				uf[k][j][i] = wxw[i]*u[k][j][i-1] + HALF*u[k][j][i] + wxe[i]*u[k][j][i+1];
			}
		}
	}

#ifndef XPERIODIC
	if (Is==0) {
		i=0;	
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				// Use one-sided filter near the boundary 
				uf[k][j][i] = wxw[i]*u[k][j][i+2] + wxp[i]*u[k][j][i] + wxe[i]*u[k][j][i+1];
			}
		}
	}
	if (Ie==NX) {
		i=NX-1;	
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				uf[k][j][i] = u[k][j][i];
			}
		}
	}
#endif // notXPERIODIC

	// Update the values at the ghost cells
	Communication_update_ghost_nodes_y(uf, 'c', 1, data_bag);

#ifdef FILTER_2D
	Array_copy_noghost(uf, uf1, grid, params); 
#else
	// Filter along y-direction 
	for (k=Ks; k<Ke; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=Is; i<Ie; i++){
				uf1[k][j][i] = wys[j]*uf[k][j-1][i] + HALF*uf[k][j][i] + wyn[j]*uf[k][j+1][i];
			}
		}
	}

	if (Js==0) {
		j=0;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++){
				// Use one-sided filter near the boundary 
				uf1[k][j][i] = wyp[j]*uf[k][j][i] + wyn[j]*uf[k][j+1][i] + wys[j]*uf[k][j+2][i];
			}
		}
	}
	if (Je==NY) {
		j=NY-2;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++){
				uf1[k][j][i] = wys[j]*uf[k][j-1][i] + HALF*uf[k][j][i] + wyn[j]*uf[k][j-2][i];
			}
		}
		j=NY-1;
		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++){
				uf1[k][j][i] = uf[k][j][i];
			}
		}
	}

#endif // FILTER_2D
	
	Communication_update_ghost_nodes_z(uf1, 'c', 1, data_bag);
	
/*
	uk = smag->uk; 
	for (k=grid->L_Ks;k<grid->L_Ke;k++) {
		for (j=Js;j<Je;j++) {
			for (i=Is;i<Ie;i++) {
				uk[j][i][k] = uf1[k][j][i];
			}
		}
	}
*/


	// Filter along z-direction
	for (k=k_start; k<k_end; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++){
				uf[k][j][i] = wzb[k]*uf1[k-1][j][i] + HALF*uf1[k][j][i] + wzf[k]*uf1[k+1][j][i];
//				uf[k][j][i] = wzb[k]*uk[j][i][k-1] + HALF*uk[j][i][k] + wzf[k]*uk[j][i][k+1];
			}
		}
	}

#ifndef ZPERIODIC
	if (Ks==0) {
		k=0;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++){
				// Use one-sided filter near the boundary
				uf[k][j][i] = wzp[k]*uf1[k][j][i] + wzf[k]*uf1[k+1][j][i] + wzb[k]*uf1[k+2][j][i];
//				uf[k][j][i] = wzp[k]*uk[j][i][k] + wzf[k]*uk[j][i][k+1] + wzb[k]*uk[j][i][k+2];
			}
		}
	}
	if (Ke==NZ) {
		k=NZ-1;
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++){
				uf[k][j][i] = uf1[k][j][i];
			}
		}
	}
#endif
	
	Tend = MPI_Wtime();
	data_bag->timer->Wtime_sgs_filter += Tend - Tstart;
	return;

}

