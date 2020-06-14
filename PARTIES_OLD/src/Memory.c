#include "definitions.h"
#include "DataTypes.h"
#include "Memory.h"
#include "Communication.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>


/******************************************************************************/
/*
 This function is to check if the memory allocation was successfull for any
 arbitrary type
 */
/******************************************************************************/
void Memory_check_allocation(void *var) {

	if (var == NULL) {

		printf("**********************************************************\n");
		printf("!!! Error. Not enough Memory. Terminating the program !!!!\n");
		printf("**********************************************************\n");
		Communication_finalize();
		exit(1);
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
double ***Memory_allocate_flow_variable(MAC_grid *grid, Parameters *params) {

	double ***array;
	double **twodarray, *onedarray;
	int istart, iend, jstart, jend, kstart, kend;
	int icount, jcount, kcount;
	int i, j, k;

	istart = grid->L_Is;
	iend   = grid->L_Ie - 1;
	jstart = grid->L_Js;
	jend   = grid->L_Je - 1;
	kstart = grid->L_Ks;
	kend   = grid->L_Ke - 1;

	icount = iend - istart + 1;
	jcount = jend - jstart + 1;
	kcount = kend - kstart + 1;

	array = (double ***) malloc( kcount*(sizeof(double**) + jcount*(sizeof(double*) + icount*sizeof(double)) ) );
	twodarray = (double **) (array + kcount);
	onedarray = (double *)(twodarray + kcount*jcount);

	for (k=0;k<kcount;k++){
		array[k] = &twodarray[k*jcount];
	}

	for (k=0;k<kcount;k++){
		for (j=0;j<jcount;j++){
			array[k][j] = &onedarray[icount*(k*jcount + j)];
		}
	}

	array = array - kstart;
	for (k=kstart;k<=kend;k++){
		array[k] = array[k]- jstart;
	}

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			array[k][j] = array[k][j]- istart;
		}
	}

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			for (i=istart;i<=iend;i++){
				array[k][j][i] = 0.0;
			}
		}
	}

	return array;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Memory_reset_flow_variable(MAC_grid *grid, Parameters *params, double ***array) {

	int i, j, k;
	int istart = grid->L_Is;
	int iend   = grid->L_Ie - 1;
	int jstart = grid->L_Js;
	int jend   = grid->L_Je - 1;
	int kstart = grid->L_Ks;
	int kend   = grid->L_Ke - 1;

	for (k = kstart; k <= kend; k++) {
		for (j = jstart; j <= jend; j++) {
			for (i = istart; i <= iend; i++) {
				array[k][j][i] = 0.0;
			}
		}
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void  Memory_free_flow_variable(MAC_grid *grid, Parameters *params, double ***array) {

//	int istart, iend, jstart, jend, kstart, kend;
//	int i, j, k;
//
//	kstart = grid->L_Ks;
//	free(array+kstart);
//	return;
//
//	istart = grid->L_Is;
//	iend   = grid->L_Ie - 1;
//	jstart = grid->L_Js;
//	jend   = grid->L_Je - 1;
//	kend   = grid->L_Ke - 1;
//
//	for (k=kstart;k<=kend;k++){
//		for (j=jstart;j<=jend;j++){
//			free(array[k][j]+istart);
//		}
//	}
//
//	for (k=kstart;k<=kend;k++) {
//		free(array[k]+jstart);
//	}
//
//	free(array+kstart);

	free(&array[grid->L_Ks]);

	return;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
int ***Memory_allocate_flow_variable_int(MAC_grid *grid, Parameters *params) {

	int ***array;
	int **twodarray, *onedarray;
	int istart, iend, jstart, jend, kstart, kend;
	int icount, jcount, kcount;
	int i, j, k;

	istart = grid->L_Is;
	iend   = grid->L_Ie - 1;
	jstart = grid->L_Js;
	jend   = grid->L_Je - 1;
	kstart = grid->L_Ks;
	kend   = grid->L_Ke - 1;

	icount = iend - istart + 1;
	jcount = jend - jstart + 1;
	kcount = kend - kstart + 1;

	array = (int ***) malloc( kcount*(sizeof(int**) + jcount*(sizeof(int*) + icount*sizeof(int)) ) );
	twodarray = (int **) (array + kcount);
	onedarray = (int *)(twodarray + kcount*jcount);

	for (k=0;k<kcount;k++){
		array[k] = &twodarray[k*jcount];
	}

	for (k=0;k<kcount;k++){
		for (j=0;j<jcount;j++){
			array[k][j] = &onedarray[icount*(k*jcount + j)];
		}
	}

	array = array - kstart;
	for (k=kstart;k<=kend;k++){
		array[k] = array[k]- jstart;
	}

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			array[k][j] = array[k][j]- istart;
		}
	}

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			for (i=istart;i<=iend;i++){
				array[k][j][i] = 0.0;
			}
		}
	}

	return array;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void  Memory_free_flow_variable_int(MAC_grid *grid, Parameters *params, int ***array) {

	int istart, iend, jstart, jend, kstart, kend;
	int i, j, k;

	kstart = grid->L_Ks;
	free(array+kstart);
	return;

	/*------------------------------------------------------------------------*/
	/*
	 End of function ...
	 */
	/*------------------------------------------------------------------------*/

	istart = grid->L_Is;
	iend   = grid->L_Ie - 1;
	jstart = grid->L_Js;
	jend   = grid->L_Je - 1;
	kend   = grid->L_Ke - 1;

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			free(array[k][j]+istart);
		}
	}

	for (k=kstart;k<=kend;k++) {
		free(array[k]+jstart);
	}

	free(array+kstart);

	return;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
double ***Memory_allocate_noghost_variable(MAC_grid *grid, Parameters *params) {

	double ***array;
	double **twodarray, *onedarray;
	int istart, iend, jstart, jend, kstart, kend;
	int icount, jcount, kcount;
	int i, j, k;

	istart = grid->G_Is;
	iend   = grid->G_Ie - 1;
	jstart = grid->G_Js;
	jend   = grid->G_Je - 1;
	kstart = grid->G_Ks;
	kend   = grid->G_Ke - 1;

	icount = iend - istart + 1;
	jcount = jend - jstart + 1;
	kcount = kend - kstart + 1;

	array = (double ***) malloc( kcount*(sizeof(double**) + jcount*(sizeof(double*) + icount*sizeof(double)) ) );
	twodarray = (double **) (array + kcount);
	onedarray = (double *)(twodarray + kcount*jcount);

	for (k=0;k<kcount;k++){
		array[k] = &twodarray[k*jcount];
	}

	for (k=0;k<kcount;k++){
		for (j=0;j<jcount;j++){
			array[k][j] = &onedarray[icount*(k*jcount + j)];
		}
	}

	array = array - kstart;
	for (k=kstart;k<=kend;k++){
		array[k] = array[k]- jstart;
	}

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			array[k][j] = array[k][j]- istart;
		}
	}

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			for (i=istart;i<=iend;i++){
				array[k][j][i] = 0.0;
			}
		}
	}

	return array;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Memory_reset_noghost_variable(MAC_grid *grid, Parameters *params, double ***array) {

	int i, j, k;
	int istart = grid->G_Is;
	int iend   = grid->G_Ie - 1;
	int jstart = grid->G_Js;
	int jend   = grid->G_Je - 1;
	int kstart = grid->G_Ks;
	int kend   = grid->G_Ke - 1;

	for (k = kstart; k <= kend; k++) {
		for (j = jstart; j <= jend; j++) {
			for (i = istart; i <= iend; i++) {
				array[k][j][i] = 0.0;
			}
		}
	}
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void  Memory_free_noghost_variable(MAC_grid *grid, Parameters *params, double ***array) {

	free(&array[grid->G_Ks]);

	return;

}




/******************************************************************************/
/*
 */
/******************************************************************************/
double ***Memory_allocate_variable_jik(MAC_grid *grid, Parameters *params) {

	double ***array;
	int istart, iend, jstart, jend, kstart, kend;
	int i, j, k;

	istart = grid->L_Is;
	iend   = grid->L_Ie - 1;
	jstart = grid->L_Js;
	jend   = grid->L_Je - 1;
	kstart = grid->L_Ks;
	kend   = grid->L_Ke - 1;

	array = malloc((jend-jstart+1)*sizeof(double**));
	array = array - jstart;
	for (j=jstart;j<=jend;j++){
		array[j] = malloc((iend-istart+1)*sizeof(double*));
		array[j] = array[j]- istart;
	}

	for (j=jstart;j<=jend;j++){
		for (i=istart;i<=iend;i++){
			array[j][i] = malloc((kend-kstart+1)*sizeof(double));
			array[j][i] = array[j][i]- kstart;
		}
	}

	for (k=kstart;k<=kend;k++){
		for (j=jstart;j<=jend;j++){
			for (i=istart;i<=iend;i++){
				array[j][i][k] = 0.0;
			}
		}
	}

	return array;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void  Memory_free_variable_jik(MAC_grid *grid, Parameters *params, double ***array) {

	int istart, iend, jstart, jend, kstart, kend;
	int i, j, k;

	istart = grid->L_Is;
	iend   = grid->L_Ie - 1;
	jstart = grid->L_Js;
	jend   = grid->L_Je - 1;
	kstart = grid->L_Ks;
	kend   = grid->L_Ke - 1;

	for (j=jstart;j<=jend;j++){
		for (i=istart;i<=iend;i++){
			free(array[j][i]+kstart);
		}
	}

	for (j=jstart;j<=jend;j++){
		free(array[j]+istart);
	}

	free(array+jstart);

	return;

}




/******************************************************************************/
/*
 This function allocates a 2D array for any arbitrary c-type
 */
/******************************************************************************/
void *Memory_allocate_1D_array(int data_type, int N) {

	switch (data_type) {

		case GVG_DOUBLE: {
			double *array;

			array = (double *)calloc(N, sizeof(double));
			Memory_check_allocation(array);
			return (array);
		}

		case GVG_FLOAT: {
			float *array;

			array = (float *)calloc(N, sizeof(float));
			Memory_check_allocation(array);

			return array;
		}

		case GVG_INT: {
			int *array;

			array = (int *)calloc(N, sizeof(int));
			Memory_check_allocation(array);

			return array;
		}

		case GVG_SHORT_INT: {
			short int *array;

			array = (short int *)calloc(N, sizeof(short int));
			Memory_check_allocation(array);

			return (array);
		}

		case GVG_LONG_INT: {
			long int *array;

			array = (long int *)calloc(N, sizeof(long int));
			Memory_check_allocation(array);

			return (array);
		}

		case GVG_CHAR: {
			char *array;

			array = (char *)calloc(N, sizeof(char));
			Memory_check_allocation(array);

			return (array);
		}
	}
	return NULL;
}




/******************************************************************************/
/*
 This function allocates a 2D double array ::: array[NZ][NY]
 */
/******************************************************************************/
double **Memory_allocate_2D_double_array(int NY, int NZ) {

	int k;

	double **array, *onedarray;

	array = (double **) malloc( NZ * ( sizeof(double*) + NY * sizeof(double) ) );
	Memory_check_allocation(array);

	onedarray = (double *)(array + NZ);

	for (k = 0; k < NZ; k++) {
		array[k] = &onedarray[k * NY];
	}
	memset(onedarray, 0, NZ * NY * sizeof(double));
	return (array);
}



/******************************************************************************/
/*
 This function allocates a 2D double complex array ::: array[NZ][NY]
 */
/******************************************************************************/
double complex **Memory_allocate_2D_double_complex_array(int NY, int NZ) {

	int k;

	double complex **array, *onedarray;

	array = (double complex **) malloc( NZ * ( sizeof(double complex*) + NY * sizeof(double complex) ) );
	Memory_check_allocation(array);

	onedarray = (double complex *)(array + NZ);

	for (k = 0; k < NZ; k++) {
		array[k] = &onedarray[k * NY];
	}
	memset(onedarray, 0, NZ * NY * sizeof(double complex));
	return (array);
}




/******************************************************************************/
/*
 This function allocates a 2D int array ::: array[NZ][NY]
 */
/******************************************************************************/
int **Memory_allocate_2D_int_array(int NY, int NZ) {

	int k;

	int **array, *onedarray;

	array = (int **) malloc( NZ * ( sizeof(int*) + NY * sizeof(int) ) );
	Memory_check_allocation(array);

	onedarray = (int *)(array + NZ);

	for (k = 0; k < NZ; k++) {
		array[k] = &onedarray[k * NY];
	}
	memset(onedarray, 0, NZ * NY * sizeof(int));
	return (array);
}




/******************************************************************************/
/*
 This function allocates a 2D char array ::: array[NZ][NY]
 */
/******************************************************************************/
char **Memory_allocate_2D_char_array(int NY, int NZ) {

	int k;

	char **array, *onedarray;

	array = (char **) malloc( NZ * ( sizeof(char*) + NY * sizeof(char) ) );
	Memory_check_allocation(array);

	onedarray = (char *)(array + NZ);

	for (k = 0; k < NZ; k++) {
		array[k] = &onedarray[k * NY];
	}
	memset(onedarray, 0, NZ * NY * sizeof(char));
	return (array);
}




/******************************************************************************/
/*
 This macro allocates a 3D double array :: array[NZ][NY][NX]
 */
/******************************************************************************/
double ***Memory_allocate_3D_double_array(int NX, int NY, int NZ) {

	int j, k;

	double ***array, **twodarray, *onedarray;

	array = (double ***) malloc( NZ * (sizeof(double**) + NY * (sizeof(double*)+ NX * sizeof(double)) ) );
	twodarray = (double **) (array + NZ);
	onedarray = (double *)(twodarray + NZ * NY);

	for (k = 0; k < NZ; k++){
		array[k] = &twodarray[k * NY];
	}

	for (k = 0; k < NZ; k++){
		for (j = 0; j < NY; j++){
			array[k][j] = &onedarray[NX * (k * NY + j)];
		}
	}
	memset(onedarray, 0, NZ * NY * NX * sizeof(double));
	return (array);
}




/******************************************************************************/
/*
 This fuction frees the allocate memory for a 2D array
 */
/******************************************************************************/
void Memory_free_2D_double_array(int NZ, double **array) {

	free(array);
}




/******************************************************************************/
/*
 This fuction frees the allocate memory for a 2D complex array
 */
/******************************************************************************/
void Memory_free_2D_double_complex_array(int NZ, double **array) {

	free(array);
}




/******************************************************************************/
/*
 This fuction frees the allocate memory for a 2D array
 */
/******************************************************************************/
void Memory_free_2D_int_array(int NZ, int **array) {

	free(array);
}




/******************************************************************************/
/*
 This fuction frees the allocate memory for a 2D array
 */
/******************************************************************************/
void Memory_free_2D_char_array(int NZ, char **array) {

	free(array);
}




/******************************************************************************/
/*
 This fuction frees the allocate memory for a 3D array
 */
/******************************************************************************/
void Memory_free_3D_double_array(int NY, int NZ, double ***array) {

	free(array);
}
