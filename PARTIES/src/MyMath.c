#include "definitions.h"
#include "DataTypes.h"
#include "Grid.h"
#include "MyMath.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

MPI_COMM comm3d;

/******************************************************************************/
/*
 This function interpolates any quantity based on 2 or 4 neighboring nodes's
 given values and coordinates
 */
/******************************************************************************/
double MyMath_interpolate_quantity(Scalar *q, Scalar q_int, int N) {
	
	double q0, q1, q2, q3;
	double a_Dx, a_Dy;
	double Wp0, Wp1;
	double Wp2;
	double m; 
	double q_top, q_bottom;
	double tol = 1.0e-8; 
	
	/*------------------------------------------------------------------------*/
	/*
	 NOTE: If the # of points is equal to two, just a linear interpoleation in
	 X1 direction
	 */
	/*------------------------------------------------------------------------*/
	if (N == 2) {
		
		q0 = q[0].Value;
		q1 = q[1].Value;
		
		m = (q1 - q0) / (q[1].X1 - q[0].X1);
		
		// check if m!= inf
		if (fabs(q[1].X1 - q[0].X1) < tol ) {

			q_int.Value = q[0].Value; // constant value everywhere
		}
		else {
			
			// Linear interpolation: y = m*(x-x1) + y1
			q_int.Value = m*(q_int.X1 - q[0].X1) + q0;   
		} // else
		
		return (q_int.Value);
	} // if N = 2
	
	/*------------------------------------------------------------------------*/
	/* 
	        2---------3
	        |         |
	        |         |
	        |         |
	        0---------1
	 */
	/*------------------------------------------------------------------------*/
	
	if (N == 4) {
		
		q0 = q[0].Value;
		q1 = q[1].Value;
		q2 = q[2].Value;
		q3 = q[3].Value;
		
		a_Dx = 1.0/(q[1].X1 - q[0].X1);
		Wp0  = fabs( (q_int.X1 - q[0].X1) * a_Dx); 
		Wp1  = 1.0 - Wp0;
		
		q_top    = Wp0*q2 + Wp1*q3;
		q_bottom = Wp0*q0 + Wp1*q1;

		a_Dy = 1.0/(q[2].X2 - q[0].X2);
		Wp0  = fabs( (q_int.X2 - q[0].X2) * a_Dy); 
		Wp2  = 1.0 - Wp0;
		
		q_int.Value = Wp0*q_bottom + Wp2*q_top;
		return (q_int.Value);
	}
	
	return (-1.0);
}

// Retruns a random number between [0, max)
double MyMath_random(double max) {

	double rnd;

	rnd = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
	return (max*rnd); 
}




/******************************************************************************/
/*
 This function returns the distance between the two points
 */
/******************************************************************************/
double MyMath_get_point_point_distance(PointType *p1, PointType *p2) {
	
	double d2;
	
	d2 = (p1->x - p2->x)*(p1->x - p2->x) + (p1->y - p2->y)*(p1->y - p2->y) +
	     (p1->z - p2->z)*(p1->z - p2->z);

	return (sqrt(d2)); 
}




/******************************************************************************/
/*
 This function computes the determinent of the matrix
 */
/******************************************************************************/
double MyMath_matrix_determ(MatrixType *mat) {

	int size; 
	double a11, a12, a13; 
	double a21, a22, a23; 
	double a31, a32, a33; 
	double det=0.0; 

	size = mat->size; 

	if (size == 2) {

		a11 = mat->A[0][0]; 
		a12 = mat->A[0][1]; 

		a21 = mat->A[1][0]; 
		a22 = mat->A[1][1]; 

		det = mat->det; 
		mat->det = a11*a22 - a21*a12;
	} // if size == 2

	if (size == 3) {

		a11 = mat->A[0][0]; 
		a12 = mat->A[0][1]; 
		a13 = mat->A[0][2]; 

		a21 = mat->A[1][0]; 
		a22 = mat->A[1][1]; 
		a23 = mat->A[1][2]; 

		a31 = mat->A[2][0]; 
		a32 = mat->A[2][1]; 
		a33 = mat->A[2][2]; 
		
		mat->det = a11*a22*a33 + a21*a32*a13 + a31*a12*a23 - a11*a32*a23 - a31*a22*a13 - a21*a12*a33; 
		det = mat->det; 
	} // if size 3

	if (size > 3) {
	
	} // if size > 3



	return det; 
}




/******************************************************************************/
/*
 This function inverts the matrices of size 2 or 3
 */
/******************************************************************************/
void MyMath_matrix_inv(MatrixType *mat) {

	int size; 
	double tol=1.0e-8; // determinent of the matrix
	double det, a_det; 
	double a11, a12, a13; 
	double a21, a22, a23; 
	double a31, a32, a33; 

	size = mat->size; 

	if (size == 2) {

		a11 = mat->A[0][0]; 
		a12 = mat->A[0][1]; 

		a21 = mat->A[1][0]; 
		a22 = mat->A[1][1]; 

		
		// find the determinent of the current matrix
		det = MyMath_matrix_determ(mat);
		
		a_det = 1.0/det; 

		// inverse of the matrix
		if ( fabs(det) > tol) {

			mat->A_inv[0][0] = a22  * a_det; 
			mat->A_inv[0][1] = -a12 * a_det; 

			mat->A_inv[1][0] = -a21 * a_det; 
			mat->A_inv[1][1] = a11  * a_det; 

		}
		else {
			
			printf("MyMath.c/ Error inverting the matrix. Determinent close to zero:%2.14f\n", det); 
		} // else

	}
	else if (size == 3) {// if size == 2

		a11 = mat->A[0][0]; 
		a12 = mat->A[0][1]; 
		a13 = mat->A[0][2]; 

		a21 = mat->A[1][0]; 
		a22 = mat->A[1][1]; 
		a23 = mat->A[1][2]; 

		a31 = mat->A[2][0]; 
		a32 = mat->A[2][1]; 
		a33 = mat->A[2][2]; 
		
		// find the determinent of the current matrix
		det = MyMath_matrix_determ(mat);
		
		a_det = 1.0/det; 
		// inverse of the matrix
		if (fabs( mat->det) > tol) {

			mat->A_inv[0][0] = (a22*a33 - a23*a32)* a_det; 
			mat->A_inv[0][1] = (a13*a32 - a12*a33)* a_det; 
			mat->A_inv[0][2] = (a12*a23 - a13*a22)* a_det; 

			mat->A_inv[1][0] = (a23*a31 - a21*a33)* a_det; 
			mat->A_inv[1][1] = (a11*a33 - a13*a31)* a_det; 
			mat->A_inv[1][2] = (a13*a21 - a11*a23)* a_det; 

			mat->A_inv[2][0] = (a21*a32 - a22*a31)* a_det; 
			mat->A_inv[2][1] = (a12*a31 - a11*a32)* a_det; 
			mat->A_inv[2][2] = (a11*a22 - a12*a21)* a_det; 
		} else {
			
			printf("MyMath.c/ Error inverting the matrix. Determinent close to zero:%f\n", mat->det); 
		} // else
	}
	else {// if size == 3
		
		MyMath_GaussJordan_matrix_inv(mat); 		
	}

}




/******************************************************************************/
/*
 This function inverses a matrix using Gauss-Jordan method.
 
 Input matrix is in mat->A[][] (size by size)
 
 Inverse matrix will be stored in mat->A_inv[][] (size by size)
 */
/******************************************************************************/
void MyMath_GaussJordan_matrix_inv(MatrixType *mat) {

	int *index_col, *index_row, *index_piv;
	int i, j, k, h, hh;
	double max, temp, a_piv;
	int in_row, in_col; 
	double _temp;
	double tol = 1.0e-8; 
	int n = mat->size; 
	
	index_col = (int *)calloc(n, sizeof(int)); 
	index_row = (int *)calloc(n, sizeof(int)); 
	index_piv = (int *)calloc(n, sizeof(int)); 

	// copy the values
	for (j=0; j<n; j++) {
		for (i=0; i<n; i++) {
			mat->A_inv[j][i] = mat->A[j][i]; 
		} // for i
	} // for j
	
	// pivot index
	for (j=0; j<n; j++) {
		index_piv[j] = -1;
	} // for j
	
	for (i=0; i<n; i++) { 
		max=0.0;

		for (j=0; j<n; j++) {
			if (index_piv[j] != 0) {
				for (k=0; k<n; k++) {
					if (index_piv[k] == -1) {
						if (fabs(mat->A_inv[j][k]) >= max) {
							max = fabs(mat->A_inv[j][k]);
							in_row=j;
							in_col=k;
						} // if
					} // if index_piv
				} // for k
			} // if index_piv
		}// for j
		++(index_piv[in_col]);

		if (in_row != in_col) {

			for (h=0; h<n; h++) { // swap

				_temp = mat->A_inv[in_row][h]; 
				mat->A_inv[in_row][h] = mat->A_inv[in_col][h]; 
				mat->A_inv[in_col][h] = _temp; 
			} // for s
		} // if in_row

		index_row[i] = in_row; 
		index_col[i] = in_col; 
		if ( fabs(mat->A_inv[in_col][in_col]) <= tol ) {

			printf("MyMath.c/ Error inverting the matrix. Singular matrix...\n"); 
		} // if

		a_piv = 1.0/mat->A_inv[in_col][in_col];
		mat->A_inv[in_col][in_col] = 1.0;

		for (h=0; h<n; h++) mat->A_inv[in_col][h] *= a_piv;

		for (hh=0; hh<n; hh++) {

			if (hh != in_col) { 
				temp = mat->A_inv[hh][in_col];
				mat->A_inv[hh][in_col]=0.0;
				for (h=0; h<n; h++) mat->A_inv[hh][h] -= mat->A_inv[in_col][h] * temp;
			}
		} // for hh

	} // for i

	for (h=n-1; h>=0; h--) {
		if (index_row[h] != index_col[h])
		for (k=0; k<n; k++) { // swap

			_temp = mat->A_inv[k][index_row[h]]; 
			mat->A_inv[k][index_row[h]] = mat->A_inv[k][index_col[h]]; 
			mat->A_inv[k][index_col[h]] = _temp; 

		} // for k
	} // for h

	free(index_col); 
	free(index_row); 
	free(index_piv);
}




/******************************************************************************/
/*
 This function does a trilinear interpolation for any given 3D data array and an
 intepolation point at (x_int,y_int,z_int)
 
                 F4------------F6
                /|            /|
               / |           / |
              /  |          /  |
            F5------------F7   |
            |   F0---------|---F2
            |  /     I     |  /
            | /            | /
            |/             |/
            F1-------------F3
 
         z
         |
         |___y
        /
       /
      x
 */
/******************************************************************************/
double MyMath_do_trilinear_inter(double x_int, double y_int, double z_int, double ***data, MAC_grid *grid, char which_quantity) {

	double *xq, *yq, *zq;
	double x0, y0, z0; 
	double x7, y7, z7; 
	double val_int; 
	int i0, j0, k0; 
	int m; 
	double alpha, beta, gamma; 
	double wc[8]; 
	int static dir_x[8] = {0, 1, 0, 1, 0, 1, 0, 1}; 
	int static dir_y[8] = {0, 0, 1, 1, 0, 0, 1, 1}; 
	int static dir_z[8] = {0, 0, 0, 0, 1, 1, 1, 1}; 

	// First the grid coordinates to be used for intepolation
	switch (which_quantity) {
		
		case 'u':
			
			xq = grid->xu;
			yq = grid->yc;
			zq = grid->zc;
			break; 

		case 'v':
			
			xq = grid->xc;
			yq = grid->yv;
			zq = grid->zc;
			break; 

		case 'w':
			
			xq = grid->xc;
			yq = grid->yc;
			zq = grid->zw;
			break; 

		case 'c':
			
			xq = grid->xc; 
			yq = grid->yc; 
			zq = grid->zc; 
			break; 
		
		default: 
			
			printf("MyMath.c/ Error setting coordinates. Unknown quantity %c\n", which_quantity); 
	} // switch

	// Now, get the index of the nodes to the left and bottom of the 
	// interpolation point
	i0  = Grid_get_x_index(x_int, grid, which_quantity); 
	j0  = Grid_get_y_index(y_int, grid, which_quantity); 
	k0  = Grid_get_z_index(z_int, grid, which_quantity); 

	// coordinates of the the two corner nodes, i.e. F0 and F7
	x0 = xq[i0]; 
	y0 = yq[j0];
	z0 = yq[k0];
	
	/* F7 */
	x7 = xq[i0+1]; 
	y7 = yq[j0+1];
	z7 = yq[k0+1];

	// Trilinear interpolation coefficients
	alpha = (x7 - x_int)/(x7 - x0); 
	beta  = (y7 - y_int)/(y7 - y0); 
	gamma = (z7 - z_int)/(z7 - z0); 

	// Get the interpolation coefficients to find the value of the image node 
	// based on the neigboring fluid nodes
	wc[0] = alpha * beta * gamma; 
	wc[1] = (1.0-alpha) * beta * gamma; 
	wc[2] = alpha * (1.0-beta) * gamma;
	wc[3] = (1.0-alpha) * (1.0-beta) * gamma;
	wc[4] = alpha * beta * (1.0-gamma); 
	wc[5] = (1.0-alpha) * beta * (1.0-gamma); 
	wc[6] = alpha * (1.0-beta) * (1.0-gamma); 
	wc[7] = (1.0-alpha)*(1.0-beta)*(1.0-gamma);
	
	val_int = 0.0; 
	// Go through all 8 nodes and add the contribution to the interpolated 
	// value
	for (m=0; m<8; m++) {
		
		val_int += wc[m]*data[ k0+dir_z[m] ][ j0+dir_y[m] ][ i0+dir_x[m] ]; 
	} // for m

	return (val_int); 
} 


