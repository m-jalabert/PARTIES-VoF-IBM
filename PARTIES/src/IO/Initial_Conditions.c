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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>





/******************************************************************************/
void Conc_init_lock(Concentration *c, Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i, j, k;
	double x_fr, y_fr, z_fr;
	double *xc, *yc, *zc;
	double ***conc;
	double c_test;
	double Re;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
	double xrandom;

	int status;
	char message[500];

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	int Lx = params -> Lx;
	double  **inflow_n=c -> inflow_n;

	conc = c->data;

	Re = params->Re;

	x_fr = params->x_fr;
	y_fr = params->y_fr;
	z_fr = params->z_fr;

	xc   = grid->xc;
	yc   = grid->yc;
	zc   = grid->zc;

	// y-coordinate where concentration profile starts, below concentratin will be one
	double conc_init_y0 = 0.5*params->Ly;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	// Assume the lock is a block. It could be located anywhere within the
	// domain.
	// One has to define the coordinates of the bottom-left-back and
	// top-right-front corner
	xb_s = 0.0;
	xb_e = 0.2 * params->Lx;
	yb_s = 0.0;
	yb_e = 0.5 * params->Ly;
	zb_s = 0.3 * params->Lz;
	zb_e = 0.7 * params->Lz;
	srand(time(NULL)+params->rank);

	if (params -> vel_init_type == VEL_INIT_PRECURSOR) {
		int  jj;
		FILE *fp ;
		char buff[255];
		char filename[50];
		float y_in[301], vel_in[301], t_in;
		double dvel;



		// Open precursor file
		params->ninflow = 500;
		sprintf(filename, "./inflow/c_in_%04d.dat", params->ninflow);


		// Read data
		fp = fopen(filename, "r");
		status = fscanf(fp, "%s", buff);
		if (status > 0) status = fscanf(fp, "%g", &t_in);
		if (status > 0) status = fscanf(fp, "%s", buff);
		if (status > 0) status = fscanf(fp, "%s", buff);

		for (j = 0; j<301; j++){
			if (status > 0) status = fscanf(fp, "%g", &y_in[j]);
			if (status > 0) status = fscanf(fp, "%g", &vel_in[j]);
		}
		fclose(fp);

		params -> t_in_n = t_in;

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

		sprintf(message, "Could not read initial conc data");
		Display_assert_error(status, message, params, DTRACE("Display_assert_error"));
	}


    time_t t;
    /* Intializes random number generator */
     srand((unsigned) time(&t)*params->rank);

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

//				xrandom = rand()/( (double) RAND_MAX);
//				xrandom = (xrandom-0.5)/10;
//				xrandom = 0.0;
//				c_test = 1.0/4.0*
//						(1.0-erf((xc[i]-x_fr+xrandom)*5.)) *
//						(1.0-erf((yc[j]-y_fr)*10.));

//				// Use this if we have a partial lock in z-direction
//				c_test = 1.0/8.0*
//						(1.0-erf((xc[i]-x_fr)*sqrt(Re))) *
//						(1.0-erf((yc[j]-y_fr)*sqrt(Re))) *
//						(1.0-erf((zc[k]-z_fr)*sqrt(Re))) ;

//				if (c_test > 1.0)
//					conc[k][j][i] = 1.0;
//
//				else if (c_test < 0.0)
//					conc[k][j][i] = 0.0;
//				else
//					conc[k][j][i] = c_test;
/*		if (params -> vel_init_type == VEL_INIT_PRECURSOR) {
					conc[k][j][i] = inflow_n[k][j]; //
				}
				else{
					conc[k][j][i] = 1.0/4.0*
							(1.0-erf((xc[i]-x_fr)*0.5 * (double)NY)) *
							(1.0-erf((yc[j]-y_fr)*0.5 * (double)NY));

					conc[k][j][i]+= 1.0/4.0*
							(1.0-erf((xc[i]-2.0 * Lx   )*0.5 * (double)NY)) *
							(1.0-erf((yc[j]-conc_init_y0)*0.5 * (double)NY));


					if (conc[k][j][i] > 1.0)
							conc[k][j][i] = 1.0;  */

				if (xc[i] < x_fr)
					conc[k][j][i] = 1+((double)rand()/(double)RAND_MAX - 0.5)/1e2;
				else
					conc[k][j][i] = 0;



			}
		}
	}

//	i = 10;
//	k = 5;
//	for (j = Js; j<Js+60; j++){
//		printf("%g %g\n", yc[j], conc[k][j][i]);
//	}


#ifdef IMMERSED_BOUNDARY
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				if (grid->c_status[k][j][i] == SOLID) {
					conc[k][j][i] = 0.0;

				}
			}
		}
	}
	Immersed_conc_solid(c, data_bag);
#endif

	return;

}

/******************************************************************************/
/*
 This function sets concentraion to zero
 */
/******************************************************************************/
void Conc_init_zero(Concentration *c, Cart3d_bag *data_bag) {
	int i, j, k;
	double ***conc= c->data;
	int Is, Js, Ks;
	int Ie, Je, Ke;


	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;


	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

					conc[k][j][i] = 0;

			}
		}
	}



	return;
}



/******************************************************************************/
/*
 This function sets a half cosine profile for testing purpose of the diffusion solver
 */
/******************************************************************************/

void Conc_init_cos(Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;

		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;

		printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
		printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
        tmpy= PI/params->Ly;
        tmpz= PI/params->Lz;


		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {


					conc[k][j][i] = (1+cos(tmpx*xc[i])*cos(tmpy*yc[j])*cos(tmpz*zc[k]))/2.;


				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;

	}

void Conc_init_RB_old(Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;

		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;

		printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
		printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
        tmpy= PI/params->Ly;
        tmpz= PI/params->Lz;

        time_t t;
        /* Intializes random number generator */
        srand((unsigned) time(&t)*params->rank);



		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {


					conc[k][j][i] = (1.- yc[j]/params->Ly)+((double)rand()/(double)RAND_MAX)/1e2;  // we add


				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;

	}
// linear varying concentration with zero value at y=0.75 Ly!
void Conc_init_Ardekani(Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;

		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;


        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
        tmpy= PI/params->Ly;
        tmpz= PI/params->Lz;

        time_t t;
        /* Intializes random number generator */
        srand((unsigned) time(&t)*params->rank);



		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {


					conc[k][j][i] = ( yc[j]- params->Ly*0.75)+0*((double)rand()/(double)RAND_MAX);  // we add


				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;

	}
void Conc_init_sin(Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;

		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;

		printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
		printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
        tmpy= PI/params->Ly;
        tmpz= PI/params->Lz;


		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {


					conc[k][j][i] = 1+(sin(tmpx* xc[i])*sin(tmpy*yc[j])*sin(tmpz*zc[k]));


				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;

	}



void Conc_init_fullsin(Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;

		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;

		printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
		printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

        tmpx= 2*PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
        tmpy= PI/params->Ly;
        tmpz= 2*PI/params->Lz;


		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {


					conc[k][j][i] = 1+( sin(tmpx* xc[i])*sin(tmpy*yc[j])*sin(tmpz*zc[k]) );


				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;

	}


void Conc_init_xsin(Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;

		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;

		printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
		printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

        tmpx= 2*2*PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber


		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {


					conc[k][j][i] =(sin(tmpx* xc[i]) );


				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;

}

/******************************************************************************/
/*
 This function sets the initial lock profile
 */
/******************************************************************************/
void Conc_initialize_saltfinger ( Concentration **c, MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int lock_smooth;
	double x_fr, y_fr, z_fr;
	double *xc, *yc, *zc;
	double ***conc, ***salt;
	double c_test;
	double Re;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double erfOfPerturb;
	int NX, NY, NZ;
	double y_int, alpha, l_c;
	// Read the binary 2D perturbation data (array size NX \times  NZ)
	double *delta1d;
	double **delta;  // <----------- dimensions available?
	FILE *file;
	char *filename;
	size_t readSize;
	int index;
	double dmean, drms;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;
	delta = Memory_allocate_2D_double_array(NZ-1, NX-1);
	delta1d = Memory_allocate_1D_array(GVG_DOUBLE, (NZ-1)*(NX-1));

	// Get concentration data
	conc = c[0]->data;
	salt = c[1]->data;

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

	Is = max(Is, 0);
	Js = max(Js, 0);
	Ks = max(Ks, 0);

	Ie = min(Ie, NX-1);
	Je = min(Je, NY-1);
	Ke = min(Ke, NZ-1);

	// Define the parameters related to initial perturbation
	y_int = 0.0;            // salt/sediment interface location
	alpha = 0.25;           // Energy parameter (depth of the perturbation)
	l_c = 0.1*params->ymax; // Error function parameter
	l_c = 2.9102;           // Error function parameter
	l_c = 1.5;              // Error function parameter



	filename = "perturbation.bin";  // <------------- create a parameter?


	if (delta == NULL)
		printf("Error while allocating memory for perturbation field.\n");

	file = fopen(filename, "rb");
	if (file != NULL) {
		readSize = fread(delta1d, sizeof(double), (NX-1) * (NZ-1), file);
		if ( ferror(file) )
			printf("Error while reading file.\n");
		if ( feof(file) )
			printf("Error: reached end of file.\n");
		if (readSize != (NX-1) * (NZ-1))
			printf("Error: read %u chunks from file (expected: %d).\n",
				   (unsigned int) readSize, (NX-1) * (NZ-1));
	}
	else {
			printf("Error: Could not open file of perturbation data.\n");
	}
	fclose(file);

	dmean = 0.0;
	drms = 0.0;
	for (k=0;k<NZ-1;k++) {
		index = k*(NX-1);
		for (i=0;i<NX-1;i++) {
			delta[i][k] = delta1d[index + i];
			dmean = dmean + delta[i][k];
			drms = drms + delta[i][k]*delta[i][k];

		}
	}
	dmean = dmean/((NX-1)*(NZ-1));
	drms = drms/((NX-1)*(NZ-1)) - dmean*dmean;

//	printf("Before the conc field initialization %d\n",params->rank);
	printf("dmean = %16.10e drms = %16.10e rank = %d\n",dmean, drms, params->rank);
	MPI_Barrier(PCW);

	// Calculate erf(f(x,y,z)) where f(x,y,z) = y - y_int + alpha * delta(x,z).
	// Set particle and salinity concentration fields.
	for (k = Ks; k < Ke; k++) {
		for (j = Js;j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				erfOfPerturb = erf( (yc[j] - y_int + params->Ly*0.5*delta[i][k]) / l_c );
				conc[k][j][i] = 0.5 * (1.0 + erfOfPerturb);
				salt[k][j][i] = 0.5 * (1.0 - erfOfPerturb);
			}
		}
	}

	printf("After the conc field initialization %d\n",params->rank);
	MPI_Barrier(PCW);

}
void Conc_init_RB(Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;

		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;

		printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
		printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
        tmpy= PI/params->Ly;
        tmpz= PI/params->Lz;

        time_t t;
        /* Intializes random number generator */
        srand((unsigned) time(&t)*params->rank);



		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {


					conc[k][j][i] = (1.- yc[j]/params->Ly)+sin(2*PI/params->Ly*xc[i] )*1e-4;  // we add


				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;






	}

void Conc_init_layer(int iconc, Concentration *c, Cart3d_bag *data_bag){

		int i, j, k;
		double x_fr, y_fr, z_fr;
		double *xc, *yc, *zc;
		double ***conc;

		double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
		double xrandom, tmpx, tmpy, tmpz;
		int status;
		char message[500];

		MAC_grid *grid = data_bag -> grid;
		Parameters *params = data_bag -> params;

		double layer_height = params->layer_height[iconc];
		double bed_start_x = params->x_fr;


		int NX = grid -> NX;
		int NY = grid -> NY;
		int NZ = grid -> NZ;


		conc = c->data;
		xc	= grid->xc;
		yc	= grid->yc;
		zc	= grid->zc;

		// Start index of bottom-left-back corner on current processor
		int Is = grid->G_Is;
		int Js = grid->G_Js;
		int Ks = grid->G_Ks;

		// End index of top-right-front corner on current processor
		int Ie = grid->G_Ie;  // Also initialize the ghost nodes
		int Je = grid->G_Je;
		int Ke = grid->G_Ke;

		printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
		printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
        tmpy= PI/params->Ly;
        tmpz= PI/params->Lz;


		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {

					if(yc[j]< layer_height && xc[i]>bed_start_x ){
						conc[k][j][i] = 1;
					}else{
						conc[k][j][i] = 0;
					}

				}
			}
		}


	#ifdef IMMERSED_BOUNDARY
		for (k=Ks; k<Ke; k++) {
			for (j=Js; j<Je; j++) {
				for (i=Is; i<Ie; i++) {
					if (grid->c_status[k][j][i] == SOLID) {
						conc[k][j][i] = 0.0;

					}
				}
			}
		}
		Immersed_conc_solid(c, data_bag);
	#endif

		return;

	}




void Conc_init_erf_x(Concentration *c, Cart3d_bag *data_bag){

			int i, j, k;
			double x_fr, y_fr, z_fr;
			double *xc, *yc, *zc;
			double ***conc;
			double xmin, xmax, sig_profile, cbd0, cbd1;

			double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
			double xrandom, tmpx, tmpy, tmpz;

			int status;
			char message[500];

			MAC_grid *grid = data_bag -> grid;
			Parameters *params = data_bag -> params;

			int NX = grid -> NX;
			int NY = grid -> NY;
			int NZ = grid -> NZ;


			xmin = params->xmin;
			xmax = params->xmax;
			sig_profile = params->sig_profile;
			cbd0 = params->cbd0;
			cbd1 = params->cbd1;

			conc = c->data;
			xc	= grid->xc;
			yc	= grid->yc;
			zc	= grid->zc;

			// Start index of bottom-left-back corner on current processor
			int Is = grid->G_Is;
			int Js = grid->G_Js;
			int Ks = grid->G_Ks;

			// End index of top-right-front corner on current processor
			int Ie = grid->G_Ie;  // Also initialize the ghost nodes
			int Je = grid->G_Je;
			int Ke = grid->G_Ke;

			printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
			printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

	        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
	        tmpy= PI/params->Ly;
	        tmpz= PI/params->Lz;


			for (k=Ks; k<Ke; k++) {
				for (j=Js; j<Je; j++) {
					for (i=Is; i<Ie; i++) {


						conc[k][j][i] = cbd0 + (cbd1-cbd0)*0.5* (1.0+erf((xc[i] - (xmax+xmin)/2.0)/sqrt(sig_profile)));


					}
				}
			}


		#ifdef IMMERSED_BOUNDARY
			for (k=Ks; k<Ke; k++) {
				for (j=Js; j<Je; j++) {
					for (i=Is; i<Ie; i++) {
						if (grid->c_status[k][j][i] == SOLID) {
							conc[k][j][i] = 0.0;

						}
					}
				}
			}
			Immersed_conc_solid(c, data_bag);
		#endif

			return;

		}


		void Conc_init_intrusion_DD(int iconc, Concentration *c, Cart3d_bag *data_bag){

					int i, j, k;
					double x_fr, y_fr, z_fr;
					double *xc, *yc, *zc;
					double ***conc;
					double xmin, xmax, sig_profile, cbd0, cbd1, cbd2, cbd3, cbd4, cbd5, Ly, ymin;

					double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
					double xrandom, tmpx, tmpy, tmpz;

					int status;
					char message[500];

					MAC_grid *grid = data_bag -> grid;
					Parameters *params = data_bag -> params;

					int NX = grid -> NX;
					int NY = grid -> NY;
					int NZ = grid -> NZ;

					x_fr = params->x_fr;

					xmin = params->xmin;
					xmax = params->xmax;
					ymin = params->ymin;
					Ly = params->Ly;


					sig_profile = params->sig_profile;
					cbd0 = params->cbd0;
					cbd1 = params->cbd1;
					cbd2 = params->cbd2;
					cbd3 = params->cbd3;
					cbd4 = params->cbd4;
					cbd5 = params->cbd5;


					conc = c->data;
					xc	= grid->xc;
					yc	= grid->yc;
					zc	= grid->zc;

					// Start index of bottom-left-back corner on current processor
					int Is = grid->G_Is;
					int Js = grid->G_Js;
					int Ks = grid->G_Ks;

					// End index of top-right-front corner on current processor
					int Ie = grid->G_Ie;  // Also initialize the ghost nodes
					int Je = grid->G_Je;
					int Ke = grid->G_Ke;

					printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
					printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

			        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
			        tmpy= PI/params->Ly;
			        tmpz= PI/params->Lz;

					if(iconc==0){
						for (k=Ks; k<Ke; k++) {
							for (j=Js; j<Je; j++) {
								for (i=Is; i<Ie; i++) {
									conc[k][j][i] = cbd2*0.5* (1.0-erf((xc[i] - x_fr)/sqrt(sig_profile))) + 0.5*(1.0 + erf((xc[i] - x_fr)/sqrt(sig_profile))) * (cbd1+(yc[j]-ymin)/Ly*(cbd0-cbd1)) ;
								}
							}
						}
					}
					if(iconc==1){
						for (k=Ks; k<Ke; k++) {
							for (j=Js; j<Je; j++) {
								for (i=Is; i<Ie; i++) {
									conc[k][j][i] = cbd5*0.5* (1.0-erf((xc[i] - x_fr)/sqrt(sig_profile))) + 0.5*(1.0 + erf((xc[i] - x_fr)/sqrt(sig_profile))) * (cbd4+(yc[j]-ymin)/Ly*(cbd3-cbd4)) ;
								}
							}
						}
					}


				#ifdef IMMERSED_BOUNDARY
					for (k=Ks; k<Ke; k++) {
						for (j=Js; j<Je; j++) {
							for (i=Is; i<Ie; i++) {
								if (grid->c_status[k][j][i] == SOLID) {
									conc[k][j][i] = 0.0;

								}
							}
						}
					}
					Immersed_conc_solid(c, data_bag);
				#endif

					return;

				}

				void Conc_init_intrusion_DD_partlock(int iconc, Concentration *c, Cart3d_bag *data_bag){

							int i, j, k;
							double x_fr, y_fr, z_fr;
							double *xc, *yc, *zc;
							double ***conc;
							double xmin, xmax, sig_profile, cbd0, cbd1, cbd2, cbd3, cbd4, cbd5, Ly, ymin;

							double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
							double xrandom, tmpx, tmpy, tmpz;

							int status;
							char message[500];

							MAC_grid *grid = data_bag -> grid;
							Parameters *params = data_bag -> params;

							int NX = grid -> NX;
							int NY = grid -> NY;
							int NZ = grid -> NZ;

							x_fr = params->x_fr;

							xmin = params->xmin;
							xmax = params->xmax;
							ymin = params->ymin;
							Ly = params->Ly;


							sig_profile = params->sig_profile;

							cbd0 = params->cbd0;
							cbd1 = params->cbd1;
							cbd2 = params->cbd2;
							cbd3 = params->cbd3;
							cbd4 = params->cbd4;
							cbd5 = params->cbd5;


							conc = c->data;
							xc	= grid->xc;
							yc	= grid->yc;
							zc	= grid->zc;

							// Start index of bottom-left-back corner on current processor
							int Is = grid->G_Is;
							int Js = grid->G_Js;
							int Ks = grid->G_Ks;

							// End index of top-right-front corner on current processor
							int Ie = grid->G_Ie;  // Also initialize the ghost nodes
							int Je = grid->G_Je;
							int Ke = grid->G_Ke;

							double erfprof;
							double hlock = params->hlock;
							double y1 = 0.5 * (Ly-hlock);
							double y2 = 0.5 * (Ly+hlock);

							printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
							printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

							tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
							tmpy= PI/params->Ly;
							tmpz= PI/params->Lz;

							if(iconc==0){
								for (k=Ks; k<Ke; k++) {
									for (j=Js; j<Je; j++) {
										for (i=Is; i<Ie; i++) {
											erfprof = 0.25 * (1.0-erf((xc[i] - x_fr)/sqrt(sig_profile))) * (erf((yc[j] - y1)/sqrt(sig_profile)) + erf(( - yc[j] + y2)/sqrt(sig_profile)));
											conc[k][j][i] = cbd2*erfprof  + (1.0-erfprof) * (cbd1+(yc[j]-ymin)/Ly*(cbd0-cbd1)) ;
										}
									}
								}
							}
							if(iconc==1){
								for (k=Ks; k<Ke; k++) {
									for (j=Js; j<Je; j++) {
										for (i=Is; i<Ie; i++) {
											erfprof = 0.25 * (1.0-erf((xc[i] - x_fr)/sqrt(sig_profile))) * (erf((yc[j] - y1)/sqrt(sig_profile)) + erf(( - yc[j] + y2)/sqrt(sig_profile)));
											conc[k][j][i] = cbd5*erfprof  + (1.0-erfprof) * (cbd4+(yc[j]-ymin)/Ly*(cbd3-cbd4)) ;										}
									}
								}
							}


						#ifdef IMMERSED_BOUNDARY
							for (k=Ks; k<Ke; k++) {
								for (j=Js; j<Je; j++) {
									for (i=Is; i<Ie; i++) {
										if (grid->c_status[k][j][i] == SOLID) {
											conc[k][j][i] = 0.0;

										}
									}
								}
							}
							Immersed_conc_solid(c, data_bag);
						#endif

							return;

						}



void Conc_init_erf_y(Concentration *c, Cart3d_bag *data_bag){

			int i, j, k;
			double x_fr, y_fr, z_fr;
			double *xc, *yc, *zc;
			double ***conc;
			double ymin, ymax, sig_profile, cbd0, cbd1;

			double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
			double xrandom, tmpx, tmpy, tmpz;

			int status;
			char message[500];

			MAC_grid *grid = data_bag -> grid;
			Parameters *params = data_bag -> params;

			int NX = grid -> NX;
			int NY = grid -> NY;
			int NZ = grid -> NZ;


			ymin = params->ymin;
			ymax = params->ymax;
			sig_profile = params->sig_profile;
			cbd0 = params->cbd0;
			cbd1 = params->cbd1;

			conc = c->data;
			xc	= grid->xc;
			yc	= grid->yc;
			zc	= grid->zc;

			// Start index of bottom-left-back corner on current processor
			int Is = grid->G_Is;
			int Js = grid->G_Js;
			int Ks = grid->G_Ks;

			// End index of top-right-front corner on current processor
			int Ie = grid->G_Ie;  // Also initialize the ghost nodes
			int Je = grid->G_Je;
			int Ke = grid->G_Ke;

			printf("Here, NX is %d (I am process %d )\n", grid->NX, params->rank );
			printf("And G_Is is %d and G_IE is %d (I am process %d )\n",grid->G_Is ,grid->G_Ie, params->rank );

	        tmpx= PI/params->Lx;  // one may multiply here with any integer to increase the wavenumber
	        tmpy= PI/params->Ly;
	        tmpz= PI/params->Lz;


			for (k=Ks; k<Ke; k++) {
				for (j=Js; j<Je; j++) {
					for (i=Is; i<Ie; i++) {

						conc[k][j][i] = cbd0 + (cbd1-cbd0)*0.5* (1.0+erf((yc[j] - (ymax+ymin)/2.0)/sqrt(sig_profile)));

					}
				}
			}


		#ifdef IMMERSED_BOUNDARY
			for (k=Ks; k<Ke; k++) {
				for (j=Js; j<Je; j++) {
					for (i=Is; i<Ie; i++) {
						if (grid->c_status[k][j][i] == SOLID) {
							conc[k][j][i] = 0.0;
						}
					}
				}
			}
			Immersed_conc_solid(c, data_bag);
		#endif

			return;

		}

		void Conc_init_unity(Concentration *c, Cart3d_bag *data_bag){

					int i, j, k;
					double x_fr, y_fr, z_fr;
					double *xc, *yc, *zc;
					double ***conc;
					double ymin, ymax, sig_profile, cbd0, cbd1;

					double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;
					double xrandom, tmpx, tmpy, tmpz;

					int status;
					char message[500];

					MAC_grid *grid = data_bag -> grid;
					Parameters *params = data_bag -> params;

					int NX = grid -> NX;
					int NY = grid -> NY;
					int NZ = grid -> NZ;


					ymin = params->ymin;
					ymax = params->ymax;

					conc = c->data;
					xc	= grid->xc;
					yc	= grid->yc;
					zc	= grid->zc;

					// Start index of bottom-left-back corner on current processor
					int Is = grid->G_Is;
					int Js = grid->G_Js;
					int Ks = grid->G_Ks;

					// End index of top-right-front corner on current processor
					int Ie = grid->G_Ie;  // Also initialize the ghost nodes
					int Je = grid->G_Je;
					int Ke = grid->G_Ke;



					for (k=Ks; k<Ke; k++) {
						for (j=Js; j<Je; j++) {
							for (i=Is; i<Ie; i++) {

								conc[k][j][i] = 1.0;

							}
						}
					}


				#ifdef IMMERSED_BOUNDARY
					for (k=Ks; k<Ke; k++) {
						for (j=Js; j<Je; j++) {
							for (i=Is; i<Ie; i++) {
								if (grid->c_status[k][j][i] == SOLID) {
									conc[k][j][i] = 0.0;
								}
							}
						}
					}
					Immersed_conc_solid(c, data_bag);
				#endif

					return;

				}
