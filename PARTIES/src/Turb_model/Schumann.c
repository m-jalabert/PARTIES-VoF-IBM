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
#include "Schumann.h" 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/******************************************************************************/
/*
 This function allocates enough memory for the concentration structure based on
 parameters defined in "*params"
 */
/******************************************************************************/
Wall_model *Schumann_create(MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ;
	int NX_cell, NY_cell, NZ_cell;
	int ierr;
	Wall_model *log_law;

	// Number of real physical cells
	NX_cell = params->NXM;
	NY_cell = params->NYM;
	NZ_cell = params->NZM;

	// Add one cell in each direction to the far most end
	NX = NX_cell + 1; /* Number of grid points in x-direction (half cell added in the far most direction) */
	NY = NY_cell + 1; /* Number of grid points in y-direction (half cell added in the far most direction)*/
	NZ = NZ_cell + 1; /* Number of grid points in z-direction (half cell added in the far most direction)*/

	// Allocate Memory for derivative components that appear in tangential 
	// shear stress at the wall 
	log_law = (Wall_model *) malloc(sizeof(Wall_model));
	log_law->dudy_wm_bottom  = Memory_allocate_2D_double_array(NX, NZ);
	log_law->dwdy_wm_bottom  = Memory_allocate_2D_double_array(NX, NZ);
	log_law->dudy_wm_top  = Memory_allocate_2D_double_array(NX, NZ);
	log_law->dwdy_wm_top  = Memory_allocate_2D_double_array(NX, NZ);

	// Log law constants
	log_law->kappa = 0.41;
	log_law->b = 5.0;

	// Caclulate the intersection point between the liear and the log law
	log_law->ypintr =  Schumann_ypbufl(log_law->kappa, log_law->b);

	return log_law;


}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Schumann_calc_shear(double ***u_data, double ***v_data, double ***w_data, 
		Wall_model *log_law, MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ, NT;
	int NX_cell, NY_cell, NZ_cell;
	int ierr;
	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end; 
	int k_start, k_end;
	int iconc, NConc; 

	double utan, utau, tau, dely, Re, iRe; 
	double uvel, wvel; 

	// Number of real physical cells
	NX_cell = params->NXM;
	NY_cell = params->NYM;
	NZ_cell = params->NZM;

	// Add one cell in each direction to the far most end
	NX = NX_cell + 1; // Number of grid points in x-direction (half cell added in the far most direction)
	NY = NY_cell + 1; // Number of grid points in y-direction (half cell added in the far most direction)
	NZ = NZ_cell + 1; // Number of grid points in z-direction (half cell added in the far most direction)
	NT = NX*NY*NZ;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	if (!( (Js == 0) || (Je == NY-1) )) return;

	Re = params->Re; 
	iRe = 1.e0/Re; 

/*  Calcualte dudy based on tau */
	i_start = Is; 
	k_start = Ks; 

	i_end = Ie;
#ifndef XPERIODIC	
	i_end   = min(NX-1, Ie); 
#endif
	k_end = min(NZ-1, Ke); 

#ifdef BOTTOM_WALL_SCHUMANN
	if (Js == 0){
		dely = grid->yc[0]-grid->yv[0];
		j = 0;
		for (k=Ks; k<Ke; k++) {
			for (i=i_start; i<i_end; i++){

				wvel = 0.25*( w_data[k][j][i-1] + w_data[k][j][i] 
							+ w_data[k+1][j][i-1] + w_data[k+1][j][i]);  
 
				utan = u_data[k][j][i]*u_data[k][j][i] + wvel*wvel;
				utan = pow(utan, 0.5);
				if (fabs(utan) == 0) utan=1e-8;
				utau = Schumann_calc_utau(utan, dely, iRe, log_law->kappa, 
							log_law->b, log_law->ypintr);
				tau = utau*utau;  

				log_law->dudy_wm_bottom[k][i]  = Re*tau*u_data[k][j][i]/utan;

			}
		}
	}
#endif

#ifdef TOP_WALL_SCHUMANN
	if (Je == NY){
		dely = grid->yv[NY-1]-grid->yc[NY-2];
		j = NY-2;
		for (k=Ks; k<Ke; k++) {
			for (i=i_start; i<i_end; i++){

				wvel = 0.25*( w_data[k][j][i-1] + w_data[k][j][i] 
							+ w_data[k+1][j][i-1] + w_data[k+1][j][i]);  

				utan = u_data[k][j][i]*u_data[k][j][i] + wvel*wvel;
				utan = pow(utan, 0.5);
				if (fabs(utan) == 0) utan=1e-8;
				utau = Schumann_calc_utau(utan, dely, iRe, log_law->kappa, 
							log_law->b, log_law->ypintr);
				tau = -utau*utau;  

				log_law->dudy_wm_top[k][i]  = Re*tau*u_data[k][j][i]/utan;


			}
		}
	}
#endif

/*  Calcualte dwdy based on tau */
	i_start = Is; 
	k_start = Ks; 


	i_end   = min(NX-1, Ie); 
	k_end = Ke;
#ifndef ZPERIODIC	
	k_end = min(NZ-1, Ke); 
#endif

#ifdef BOTTOM_WALL_SCHUMANN
	if (Js == 0){
		dely = grid->yc[0]-grid->yv[0];
		j = 0;
		for (k=Ks; k<Ke; k++) {
			for (i=i_start; i<i_end; i++){

				uvel = 0.25*( u_data[k][j][i] + u_data[k][j][i+1] 
							+ u_data[k-1][j][i] + u_data[k-1][j][i+1]);  
 
				utan = uvel*uvel + w_data[k][j][i]*w_data[k][j][i];
				utan = pow(utan, 0.5);
				if (fabs(utan) == 0) utan=1e-8;
				utau = Schumann_calc_utau(utan, dely, iRe, log_law->kappa, 
							log_law->b, log_law->ypintr);
				tau = utau*utau;  

				log_law->dwdy_wm_bottom[k][i]  = Re*tau*w_data[k][j][i]/utan;

			}
		}
	}
#endif

#ifdef TOP_WALL_SCHUMANN
	if (Je == NY){
		dely = grid->yv[NY-1]-grid->yc[NY-2];
		j = NY-2;
		for (k=Ks; k<Ke; k++) {
			for (i=i_start; i<i_end; i++){

				uvel = 0.25*( u_data[k][j][i] + u_data[k][j][i+1] 
							+ u_data[k-1][j][i] + u_data[k-1][j][i+1]);  

				utan = uvel*uvel + w_data[k][j][i]*w_data[k][j][i];
				utan = pow(utan, 0.5);
				if (fabs(utan) == 0) utan=1e-8;
				utau = Schumann_calc_utau(utan, dely, iRe, log_law->kappa, 
							log_law->b, log_law->ypintr);
				tau = -utau*utau;  

				log_law->dwdy_wm_top[k][i]  = Re*tau*w_data[k][j][i]/utan;


			}
		}
	}
#endif

	return;

}


/******************************************************************************/
/*
 Calculates the intersection point between the linear (u+=y+) and the log
 (u+=1/kappa*log(y+)+b) laws
 */
/******************************************************************************/
double Schumann_ypbufl(double kappa, double b)  {

	double ypintr, errmax, f, derf;

	errmax = 1e-8;

	ypintr = 10.5;

	do {
		f = log(ypintr)/kappa - ypintr + b;
		derf = 1.0/(kappa*ypintr) - 1.0;
		ypintr = ypintr - f/derf;
	} while (fabs(f) > errmax);

	return ypintr;

}




/******************************************************************************/
/*
 Calculates the intersection point between the linear (u+=y+) and the log
 (u+=1/kappa*log(y+)+b) laws
 */
/******************************************************************************/
double Schumann_calc_utau(double u, double y, double iRe, double kappa, double bllaw, double ypintr) {
	
	double utau, lutau, f, derf, errmax; 
	double sign;

	errmax = 1e-8;
	
	sign = 1.0;
	if (u < 0) {
		u = fabs(u);
		sign = -1.0;
	} 
	lutau=pow(iRe*u/y,0.5);

	if (u == 0)
		u = 1e-6;
	if ( (lutau*y/iRe) < ypintr) {
		utau = lutau*sign; 
		return utau;
	}
	do {
		f = log(y*lutau/iRe)/kappa - u/lutau + bllaw;
		derf = (1./kappa + u/lutau)/lutau;
		lutau = lutau - f/derf;
	} while (fabs(f) > errmax);

	utau = lutau*sign;
	return utau;

}

/******************************************************************************/


