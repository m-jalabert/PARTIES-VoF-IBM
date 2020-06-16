#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "EPforcing.h"
#include "DataTypes.h"
#include "Memory.h"
#include "definitions.h"


#ifdef TURB_FORCING   //those fonctions are only useful for turbulent forcing

/*
Ref : A. Chouippe and M. Uhlmann, Phys. Fluids 27, 123301 (2015)
*/


// list of all the forced frequencies
// Kf is the maximum adimensionned force frequency, alloc becomes the final length of the returned list
int ** forced_freq(int *alloc, double Kf){
  int l,m,n;
  int npoints = floor(2.0945*Kf*Kf*Kf+2*Kf+16); //estimates the amount of points in a half sphere
  int **freq_list = (int **)malloc(npoints * sizeof(int *));  // Matrix memory allocation
  for(int i = 0; i <npoints; ++i){
    freq_list[i] = (int *)malloc(3 * sizeof(int));
  }

  *alloc = 0;

  // section 1 : kx>0
  for (l=1; l<=Kf; ++l){
    for (m=-Kf; m<=Kf; ++m){
      for (n=-Kf; n<=Kf; ++n){
        if(l*l+m*m+n*n<=Kf*Kf){
          freq_list[*alloc][0] =l;
          freq_list[*alloc][1] =m;
          freq_list[*alloc][2] =n;
          (*alloc)++;
        }
      }
    }
  }
  // section 2 : kx=0 ky>0
  l=0;
  for (m=1; m<=Kf; ++m){
    for (n=-Kf; n<=Kf; ++n){
      if(l*l+m*m+n*n<=Kf*Kf){
        freq_list[*alloc][0] =l;
        freq_list[*alloc][1] =m;
        freq_list[*alloc][2] =n;
        (*alloc)++;
      }
    }
  }
  // section 3 : kx=0 ky=0 kz>0
  m=0; // and l=0;
  for (n=1; n<=Kf; ++n){
    if(l*l+m*m+n*n<=Kf*Kf){
      freq_list[*alloc][0] =l;
      freq_list[*alloc][1] =m;
      freq_list[*alloc][2] =n;
      (*alloc)++;
    }
  }
  //Memory
  for(l=*alloc; l<npoints; ++l){
    free(freq_list[l]);
  }
  return freq_list;
}
///////////////////////////////////////


// initialize the Fourier struct
void init_turb_forcing(Fourier *fourier, Parameters *params){
  int nb_forced_freq;

  fourier -> k_list = forced_freq(&nb_forced_freq, params->Kf);
  fourier -> nb_forced_freq  = nb_forced_freq;
  fourier -> f_fourier = Memory_allocate_2D_double_complex_array(nb_forced_freq,3);
  fourier -> std_dev = sqrt(2*params->Ds /(params->Re*params->Re*params->Re));

  fourier -> b_random = Memory_allocate_2D_double_complex_array(nb_forced_freq,3);
  // double complex **b = Memory_allocate_2D_double_complex_array(nb_forced_freq,3);
  // for (int axis=0; axis<3; ++axis){
  //   for(int i=0; i<nb_forced_freq; ++i){
  //     printf("%f + i%f\n", creal(b[axis][i]), cimag(b[axis][i]));
  //     b[axis][i] = 0.0;
  //   }
  // }
  // fourier -> b_random = b;
}
///////////////////////////////////////


// function b random process [Ref: eq.(3)]
// and volume turbulent force term in fourier space [Ref: eq.(5)]
void b_and_fourier_update(Fourier *fourier, Parameters *params){
  double normk;
  double complex ratio;
  int axis;
  double complex **b = fourier -> b_random;
  double complex **f_fourier = fourier -> f_fourier;
  int **k_list = fourier -> k_list;
  int nb_freq = fourier->nb_forced_freq;

  double dt = params -> dt;
  double Re = params -> Re;
  double Ds = params -> Ds;

  double cst_b = (1-dt);
  double amplitude = fourier -> std_dev * sqrt(dt);

  // b update
  for(int i=0; i<nb_freq; ++i){
    for (axis=0; axis<3; ++axis){
      b[axis][i] = b[axis][i]*cst_b + (box_muller(0,1) + I * box_muller(0,1))*amplitude;
    }
  }

  //f_fourier update
  for(int i=0; i<nb_freq; ++i){
    ratio = 0.0;
    normk = 0.0;
    for (axis =0; axis<3; ++axis){
      ratio += k_list[i][axis]*b[axis][i];
      normk += k_list[i][axis]*k_list[i][axis];
    }
    ratio = ratio/normk;
    for (axis =0; axis<3; ++axis){
      f_fourier[axis][i] = b[axis][i] - k_list[i][axis]*ratio;
    }
  }
}
////////////////////////////////////////


// volume turbulent force term in physical space [Ref: appendix A]
// Apparently faster than using fft
// f(x) update f for one axis
void ftx_update(Cart3d_bag *data_bag){
  int i, j, k, l, m, intK;

  MAC_grid *grid = data_bag -> grid;
  Parameters *params = data_bag -> params;
  Fourier *fourier = data_bag -> fourier;
  double ***fturb = data_bag -> u -> fturb;

  double complex **f_fourier = fourier -> f_fourier;
  int **k_list = fourier -> k_list;
  int nb_freq = fourier -> nb_forced_freq;

  intK = floor(params->Kf);
  double complex A[intK+1][2*intK+1];
  double complex B[intK+1];
  double wave, k0 = 2*PI/(params->xmax - params->xmin);

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
  // Indices start and end on current processor
  int i_start = Is; // i=0 not included
  int j_start = Js;
  int k_start = Ks;
  // Exclude the half cell added
  int i_end = min(NX-1, Ie);
  int j_end = min(NY-1, Je);
  int k_end = min(NZ-1, Ke);

  double *x = grid -> xu;
  double *y = grid -> yc;
  double *z = grid -> zc;

  for (k = k_start; k < k_end; k++) {
    for(l=0; l<=intK; ++l){
      for(m=0; m<=2*intK; ++m){
        A[l][m] = 0.0;
      }
    }
    wave = z[k]*k0;
    for(i = 0; i<nb_freq; ++i){
      A[k_list[i][0]][k_list[i][1]+intK] += f_fourier[0][i]*cexp(I*k_list[i][2]*wave);  // [Ref: appendix (A3a)] Almz

    }

    for (j = j_start; j < j_end; j++) {
      for (l=0; l<=intK; ++l){
        B[l] = 0.0;
        wave = y[j]*k0;
        for (m=-intK; m<=intK; ++m){  // many 0, amelioration ??
          B[l] += A[l][m+intK]*cexp(I*m*wave);  // [Ref: appendix (A3b)] Blyz
        }
      }

      for (i = i_start; i < i_end; i++) {
        fturb[k][j][i] = 0.0;
        wave = x[i]*k0;
        for (l=0; l<=intK; ++l){
          fturb[k][j][i] += 2*creal(B[l]*cexp(I*l*wave));  // [Ref: appendix (A3c) and f_fourier(n)=conj(f_fourier(-n))]
        }
      }
    }
  }

}
////////////////////////////////////////


// volume turbulent force term in physical space [Ref: appendix A]
// Apparently faster than using fft
// f(x) update f for one axis
void fty_update(Cart3d_bag *data_bag){
  int i, j, k, l, m, intK;

  MAC_grid *grid = data_bag -> grid;
  Parameters *params = data_bag -> params;
  Fourier *fourier = data_bag -> fourier;
  double ***fturb = data_bag -> v -> fturb;

  double complex **f_fourier = fourier -> f_fourier;
  int **k_list = fourier -> k_list;
  int nb_freq = fourier -> nb_forced_freq;

  intK = floor(params->Kf);
  double complex A[intK+1][2*intK+1];
  double complex B[intK+1];
  double wave, k0 = 2*PI/(params->xmax - params->xmin);

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
  // Indices start and end on current processor
  int i_start = Is; // i=0 not included
  int j_start = Js;
  int k_start = Ks;
  // Exclude the half cell added
  int i_end = min(NX-1, Ie);
  int j_end = min(NY-1, Je);
  int k_end = min(NZ-1, Ke);

  double *x = grid -> xc;
  double *y = grid -> yv;
  double *z = grid -> zc;

  for (k = k_start; k < k_end; k++) {
    for(l=0; l<=intK; ++l){
      for(m=0; m<=2*intK; ++m){
        A[l][m] = 0.0;
      }
    }
    wave = z[k]*k0;
    for(i = 0; i<nb_freq; ++i){
      A[k_list[i][0]][k_list[i][1]+intK] += f_fourier[1][i]*cexp(I*k_list[i][2]*wave);  // [Ref: appendix (A3a)] Almz

    }

    for (j = j_start; j < j_end; j++) {
      for (l=0; l<=intK; ++l){
        B[l] = 0.0;
        wave = y[j]*k0;
        for (m=-intK; m<=intK; ++m){  // many 0, amelioration ??
          B[l] += A[l][m+intK]*cexp(I*m*wave);  // [Ref: appendix (A3b)] Blyz
        }
      }

      for (i = i_start; i < i_end; i++) {
        fturb[k][j][i] = 0.0;
        wave = x[i]*k0;
        for (l=0; l<=intK; ++l){
          fturb[k][j][i] += 2*creal(B[l]*cexp(I*l*wave));  // [Ref: appendix (A3c) and f_fourier(n)=conj(f_fourier(-n))]
        }
      }
    }
  }

}
////////////////////////////////////////


// volume turbulent force term in physical space [Ref: appendix A]
// Apparently faster than using fft
// f(x) update f for one axis
void ftz_update(Cart3d_bag *data_bag){
  int i, j, k, l, m, intK;

  MAC_grid *grid = data_bag -> grid;
  Parameters *params = data_bag -> params;
  Fourier *fourier = data_bag -> fourier;
  double ***fturb = data_bag -> w -> fturb;

  double complex **f_fourier = fourier -> f_fourier;
  int **k_list = fourier -> k_list;
  int nb_freq = fourier -> nb_forced_freq;

  intK = floor(params->Kf);
  double complex A[intK+1][2*intK+1];
  double complex B[intK+1];
  double wave, k0 = 2*PI/(params->xmax - params->xmin);

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
  // Indices start and end on current processor
  int i_start = Is; // i=0 not included
  int j_start = Js;
  int k_start = Ks;
  // Exclude the half cell added
  int i_end = min(NX-1, Ie);
  int j_end = min(NY-1, Je);
  int k_end = min(NZ-1, Ke);

  double *x = grid -> xc;
  double *y = grid -> yc;
  double *z = grid -> zw;

  for (k = k_start; k < k_end; k++) {
    for(l=0; l<=intK; ++l){
      for(m=0; m<=2*intK; ++m){
        A[l][m] = 0.0;
      }
    }
    wave = z[k]*k0;
    for(i = 0; i<nb_freq; ++i){
      A[k_list[i][0]][k_list[i][1]+intK] += f_fourier[2][i]*cexp(I*k_list[i][2]*wave);  // [Ref: appendix (A3a)] Almz

    }

    for (j = j_start; j < j_end; j++) {
      for (l=0; l<=intK; ++l){
        B[l] = 0.0;
        wave = y[j]*k0;
        for (m=-intK; m<=intK; ++m){  // many 0, amelioration ??
          B[l] += A[l][m+intK]*cexp(I*m*wave);  // [Ref: appendix (A3b)] Blyz
        }
      }

      for (i = i_start; i < i_end; i++) {
        fturb[k][j][i] = 0.0;
        wave = x[i]*k0;
        for (l=0; l<=intK; ++l){
          fturb[k][j][i] += 2*creal(B[l]*cexp(I*l*wave));  // [Ref: appendix (A3c) and f_fourier(n)=conj(f_fourier(-n))]
        }
      }
    }
  }

}
////////////////////////////////////////



/* boxmuller.c           Implements the Polar form of the Box-Muller Transformation
(c) Copyright 1994, Everett F. Carter Jr.
Permission is granted by the author to use this software for any application provided this  copyright notice is preserved.
*/
double box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;	//uniform on [-1:1]
			x2 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}
////////////////////////////////////////


#endif    // TURB_FORCING
