#include "definitions.h"
#include "DataTypes.h"
#include "Grid.h"
#include "Array.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/******************************************************************************/
/*
 Copy data1 to data2
 */
/******************************************************************************/
void Array_copy_noghost(double ***data1, double ***data2, MAC_grid *grid,
						Parameters *params) {
	
	int i, j, k;    
	int Is, Js, Ks;
	int Ie, Je, Ke;      

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				data2[k][j][i] = data1[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Copy data1 to data2
 */
/******************************************************************************/
void Array_copy_withghost(double ***data1, double ***data2, MAC_grid *grid,
						  Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				data2[k][j][i] = data1[k][j][i];
			}
		}
	}
}





/******************************************************************************/
/*
 Scale array by factor i.e data = data * factor
 */
/******************************************************************************/
void Array_scale_noghost(double ***data, double factor, MAC_grid *grid,
						 Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				data[k][j][i] = data[k][j][i] * factor;
			}
		}
	}
}




/******************************************************************************/
/*
 Scale array by factor i.e data = data * factor
 */
/******************************************************************************/
void Array_scale_withghost(double ***data, double factor, MAC_grid *grid,
						   Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				data[k][j][i] = data[k][j][i] * factor;
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate y = a*x + y
 where x & y are arrays and "a" is a scalar
 */
/******************************************************************************/

void Array_AXPY_noghost(double ***y, double a, double ***x, MAC_grid *grid,
						Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = a * x[k][j][i] + y[k][j][i];
			}
		}
	}
}



/******************************************************************************/
/*
 Calculate y = a*x + y
 where x & y are arrays and "a" is a scalar
 */
/******************************************************************************/
void Array_AXPY_withghost(double ***y, double a, double ***x, MAC_grid *grid,
						  Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = a * x[k][j][i] + y[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate w = x * y
 where w, x & y are arrays
 */
/******************************************************************************/
void Array_pointwisemult_noghost(double ***w, double ***x, double ***y,
								 MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				w[k][j][i] = x[k][j][i] * y[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate w = x * y
 where w, x & y are arrays
 */
/******************************************************************************/
void Array_pointwisemult_withghost(double ***w, double ***x, double ***y,
								   MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				w[k][j][i] = x[k][j][i] * y[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate y = alpha * x + beta * y
 where  x & y are arrays, and alpha and beta are scalars
 */
/******************************************************************************/
void Array_AXPBY_withghost(double ***y, double alpha, double beta, double ***x,
						   MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = alpha * x[k][j][i] + beta * y[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate y = alpha * x + beta * y
 where  x & y are arrays, and alpha and beta are scalars
 */
/******************************************************************************/
void Array_AXPBY_noghost(double ***y, double alpha, double beta, double ***x,
						 MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = alpha * x[k][j][i] + beta * y[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate y = x * x +  y
 where  x & y are arrays
 */
/******************************************************************************/
void Array_XXPY_withghost(double ***y, double ***x, MAC_grid *grid,
						  Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = x[k][j][i] * x[k][j][i] + y[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate y = x * x + y
 where  x & y are arrays
 */
/******************************************************************************/
void Array_XXPY_noghost(double ***y, double ***x, MAC_grid *grid,
						Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = x[k][j][i] * x[k][j][i] + y[k][j][i];
			}
		}
	}
}




/******************************************************************************/
/*
 Calculate y = alpha*x * x + y
 where  x & y are arrays and alpha is a scalar
 */
/******************************************************************************/
void Array_AXXPY_withghost(double ***y, double alpha, double ***x,
						   MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = alpha * x[k][j][i] * x[k][j][i] + y[k][j][i];
			}
		}
	}
	
}




/******************************************************************************/
/*
 Calculate y = alpha*x * x + y
 where  x & y are arrays and alpha is a scalar
 */
/******************************************************************************/
void Array_AXXPY_noghost(double ***y, double alpha, double ***x, MAC_grid *grid,
						 Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = alpha * x[k][j][i] * x[k][j][i] + y[k][j][i];
			}
		}
	}
	
}




/******************************************************************************/
/*
 Calculate y = x * z +  y
 where  x, y & and are arrays
 */
/******************************************************************************/
void Array_XZPY_withghost(double ***y, double ***x, double ***z, MAC_grid *grid,
						  Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = x[k][j][i] * z[k][j][i] + y[k][j][i];
			}
		}
	}
	
}




/******************************************************************************/
/*
 Calculate y = x * z +  y
 where  x, y & and are arrays
 */
/******************************************************************************/
void Array_XZPY_noghost(double ***y, double ***x, double ***z, MAC_grid *grid,
						Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = x[k][j][i] * z[k][j][i] + y[k][j][i];
			}
		}
	}
	
}




/******************************************************************************/
/*
 Calculate y = alpha * x * z +  y
 where  x, y & and are arrays and alpha is a scalar
 */
/******************************************************************************/
void Array_AXZPY_withghost(double ***y, double alpha, double ***x, double ***z,
						   MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = alpha * x[k][j][i] * z[k][j][i] + y[k][j][i];
			}
		}
	}
	
}




/******************************************************************************/
/*
 Calculate y = alpha * x * z +  y
 where  x, y & and are arrays and alpha is a scalar
 */
/******************************************************************************/
void Array_AXZPY_noghost(double ***y, double alpha, double ***x, double ***z,
						 MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;


	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				y[k][j][i] = alpha * x[k][j][i] * z[k][j][i] + y[k][j][i];
			}
		}
	}
	
}

#ifdef PETSC
/******************************************************************************/
/*
 */
/******************************************************************************/
void Array_copy_to_vec(DM *DA_array, Vec *global, double ***data,
					   MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double ***g_data;
	int ierr; 

	DMDAVecGetArray(*DA_array, *global, (void ***)&g_data);

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				g_data[k][j][i] = data[k][j][i];
			}
		}
	}

 	DMDAVecRestoreArray(*DA_array, *global, (void ***)&g_data);
}
#endif  //PETSC




/******************************************************************************/
/*
 Set data to alpha
 */
/******************************************************************************/
void Array_set_noghost(double ***data, double alpha, MAC_grid *grid,
					   Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				data[k][j][i] = alpha;
			}
		}
	}
}




/******************************************************************************/
/*
 Set data to alpha
 */
/******************************************************************************/
void Array_set_withghost(double ***data, double alpha, MAC_grid *grid,
						 Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;

	/* Start index of bottom-left-back corner on current processor */
	Is = grid->L_Is;
	Js = grid->L_Js;
	Ks = grid->L_Ks;
	
	/* End index of top-right-front corner on current processor */
	Ie = grid->L_Ie;
	Je = grid->L_Je;
	Ke = grid->L_Ke;

	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				data[k][j][i] = alpha;
			}
		}
	}
}

