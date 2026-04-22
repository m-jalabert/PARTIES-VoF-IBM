#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Boundary.h"
#include "definitions.h"
#include "DataTypes.h"

#include "Display.h"
#include "Interpolate.h"
#include "Memory.h"
#include "Particle.h"


#define interp_kernel(h,r) (delta(r[0] / h) * delta(r[1] / h) * delta(r[2] / h))  // constant h in all directions

#define lvl_set(x,y,z) ( sqrt( (x - X[0]) * (x - X[0]) + \
                               (y - X[1]) * (y - X[1]) + \
                               (z - X[2]) * (z - X[2]) ) / R - 1.0 )

#define ETA(x,y,z)  R - sqrt( (x - X[0]) * (x - X[0]) + (y - X[1]) * (y - X[1]) + (z - X[2]) * (z - X[2]) )   // distance from the particle surface
//#define SMOOTH(eta, delta_s) erf( max(eta,0 )/delta_s)  // version OI
#define SMOOTH(eta, delta_s) (erf( eta/delta_s)+1)/2  // version OII
#define SMOOTH_BL(eta, delta_s,bl_thick) (erf((bl_thick + eta)/delta_s)+1)/2  // version OII


/******************************************************************************/
/*
 */
/******************************************************************************/
double delta(double r) {
	r = fabs(r);

	if (r <= 0.5)
		return (1 + sqrt(-3 * r * r + 1)) / 3;
	else if (r <= 1.5)
		return (5 - 3 * r - sqrt(-3 * (1 - r) * (1 - r) + 1)) / 6;
	else
		return 0;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
double interp_kernel_nonuniform(double hx, double rx, double hy, double ry, double hz, double rz) {

	return delta(rx/hx) * delta(ry/hy) * delta(rx/hx);
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void inv_3x3(double **M, double **M_inv) {

//  double det = -M[0][0] * M[0][1]*M[0][1]+M[0][0]*M[0][0]*M[1][1] - M[0][2]*M[0][2]*M[1][1] + 2*M[0][1]*M[0][2]*M[1][2] - M[0][0]* M[1][2]*M[1][2];

	double det = -M[0][2]*M[0][2]*M[0][2] + 2*M[0][1]*M[0][2]*M[1][2] - M[0][0]*M[1][2]*M[1][2] - M[0][1]*M[0][1]*M[2][2] + M[0][0]*M[0][2]*M[2][2];

	M_inv[0][0] = (M[1][2] * M[1][2] - M[0][2] * M[2][2]) / det;
	M_inv[0][1] = (-M[0][2] * M[1][2] + M[0][1] * M[2][2]) / det;
	M_inv[0][2] = (M[0][2]*M[0][2] - M[0][1]*M[1][2]) / det;

	M_inv[1][0] = M_inv[0][1];
	M_inv[1][1] = (-M[0][2]*M[0][2] + M[0][0]*M[2][2]) / det;
	M_inv[1][2] = (M[0][1]*M[0][2] - M[0][0]*M[1][2]) / det;

	M_inv[2][0] = M_inv[0][2];
	M_inv[2][1] = M_inv[1][2];
	M_inv[2][2] = (-M[0][1]*M[0][1] + M[0][0]*M[0][2]);

}




/******************************************************************************/
/*
 */
/******************************************************************************/
void mat_vec(int m, int n, double **A, double *B, double *C) { //A[][n], double B[]){

//	double *C = calloc(m, sizeof(double));

	double sum = 0;
	int i,j;


	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {

			sum += A[i][j]*B[j];
		}

		C[i] = sum;
		sum = 0;
	}

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
double delta_mod_3D(double x, double X, double y, double Y, double z, double Z, double dilat_factor_x, double dilat_factor_y, double dilat_factor_z, double *b) {

	return ( b[0] + (y - Y)* b[1] + (y - Y)*(y - Y) * b[2] ) *
	delta((X - x) / dilat_factor_x) * delta((Y - y) / dilat_factor_y) * delta((Z - z) / dilat_factor_z);
}




/******************************************************************************/
/*
 Spread Lagrangian quantities 'A' onto Eulerian quantity 'a' for a uniform grid.

 'a' should be a noghost array
 */
/******************************************************************************/
void Interpolate_Lag_to_Eul_uniform(double *A, double ***a, char which,
        Particle *p, MAC_grid *grid) {
#ifndef GRID_UNIFORM
	printf("Stop __func__ is only working for uniform grids" );
#endif
	int i, j, k;
	int mv;

	double r[3];
	double *X_L, *Y_L, *Z_L;

	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	double *x = grid -> xc;
	double *y = grid -> yc;
	double *z = grid -> zc;

	if (which == 'u' ) {
		x = grid -> xu;
	}
	else if (which == 'v' ) {
		y = grid -> yv;
	}
	else if (which == 'w' ) {
		z = grid -> zw;
	}

	int N_L_local = p -> N_L_local;
	double h = grid -> dx_u[1];
	double dVol = p -> Vol_L / (h * h * h);

	X_L = p -> X_L;
	Y_L = p -> Y_L;
	Z_L = p -> Z_L;

	//--------------------------------------------------------------------------
	// Reset 'a'
	//--------------------------------------------------------------------------
	double *a_1d = &a[grid->G_Ks][grid->G_Js][grid->G_Is];
	for (i = 0; i < grid->ng_total_nodes; i++) {
		a_1d[i] = 0.0;
	}

	//--------------------------------------------------------------------------
	// Let 'A' = 1; used for calculating volume fraction
	//--------------------------------------------------------------------------
	if (A == NULL) {
		for (mv = 0; mv < N_L_local; mv++) {

			// Indices based on uniform grid
			i_start = (int) round((X_L[mv] - x[0]) / h) - 1;
			j_start = (int) round((Y_L[mv] - y[0]) / h) - 1;
			k_start = (int) round((Z_L[mv] - z[0]) / h) - 1;
			i_end = i_start + 3;
			j_end = j_start + 3;
			k_end = k_start + 3;

			// Keep spreading within local processor
			i_start = max(i_start, grid->G_Is);
			i_start = min(i_start, grid->G_Ie);
			i_end = max(i_end, grid->G_Is);
			i_end = min(i_end, grid->G_Ie);

			j_start = max(j_start, grid->G_Js);
			j_start = min(j_start, grid->G_Je);
			j_end = max(j_end, grid->G_Js);
			j_end = min(j_end, grid->G_Je);

			k_start = max(k_start, grid->G_Ks);
			k_start = min(k_start, grid->G_Ke);
			k_end = max(k_end, grid->G_Ks);
			k_end = min(k_end, grid->G_Ke);

			for (k = k_start; k < k_end; k++) {
				r[2] = Z_L[mv] - z[k];

				for (j = j_start; j < j_end; j++) {
					r[1] = Y_L[mv] - y[j];

					for (i = i_start; i < i_end; i++) {
						r[0] = X_L[mv] - x[i];

						a[k][j][i] += interp_kernel(h, r) * dVol;
					}
				}
			}
		}
	}
	//--------------------------------------------------------------------------
	// 'A' != 1; used for interpolating some other value
	//--------------------------------------------------------------------------
	else {
		for (mv = 0; mv < N_L_local; mv++) {

			// Indices based on uniform grid
			i_start = (int) round((X_L[mv] - x[0]) / h) - 1;
			j_start = (int) round((Y_L[mv] - y[0]) / h) - 1;
			k_start = (int) round((Z_L[mv] - z[0]) / h) - 1;
			i_end = i_start + 3;
			j_end = j_start + 3;
			k_end = k_start + 3;

			// Keep spreading within local processor
			i_start = max(i_start, grid->G_Is);
			i_start = min(i_start, grid->G_Ie);
			i_end = max(i_end, grid->G_Is);
			i_end = min(i_end, grid->G_Ie);

			j_start = max(j_start, grid->G_Js);
			j_start = min(j_start, grid->G_Je);
			j_end = max(j_end, grid->G_Js);
			j_end = min(j_end, grid->G_Je);

			k_start = max(k_start, grid->G_Ks);
			k_start = min(k_start, grid->G_Ke);
			k_end = max(k_end, grid->G_Ks);
			k_end = min(k_end, grid->G_Ke);

			for (k = k_start; k < k_end; k++) {
				r[2] = Z_L[mv] - z[k];

				for (j = j_start; j < j_end; j++) {
					r[1] = Y_L[mv] - y[j];

					for (i = i_start; i < i_end; i++) {
						r[0] = X_L[mv] - x[i];

						a[k][j][i] += A[mv] * interp_kernel(h, r) * dVol;
					}
				}
			}
		}
	}
}




#ifndef GRID_UNIFORM
/******************************************************************************/
/*
 Spread Lagrangian quantities 'A' onto Eulerian quantity 'a' for a non-uniform
 grid.

 'a' should be a noghost array
 */
/******************************************************************************/
void Interpolate_Lag_to_Eul_nonuniform(double *A, double ***a, char which,
        Particle *p, MAC_grid *grid) {

	int i, j, k;
	int mv;

	double r[3];
	double *X_L, *Y_L, *Z_L;

	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	double *x = grid -> xc;
	double *y = grid -> yc;
	double *z = grid -> zc;

	if (which == 'u' ) {
		x = grid -> xu;
	}
	else if (which == 'v' ) {
		y = grid -> yv;
	}
	else if (which == 'w' ) {
		z = grid -> zw;
	}

	int N_L_local = p -> N_L_local;

	X_L = p -> X_L;
	Y_L = p -> Y_L;
	Z_L = p -> Z_L;

	int *idx_x_markers = p -> idx_x_markers;
	int *idx_y_markers = p -> idx_y_markers;
	int *idx_z_markers = p -> idx_z_markers;

	double *x_markers = p -> X_L;
	double *y_markers = p -> Y_L;
	double *z_markers = p -> Z_L;

	for (i=0; i< N_L_local; i++) {

		if(fabs(y_markers[i] - y[idx_y_markers[i]+1]) <
		   fabs(y_markers[i] - y[idx_y_markers[i]] ) ) {

			idx_y_markers[i] = idx_y_markers[i] + 1;
		}
	}

	double e1[3] = {1.0, 0.0, 0.0};

	double dx = x[3] - x[2];
	double dz = z[3] - z[2];
	double x_max, y_max, x_min, y_min;
	double **M = Memory_allocate_2D_double_array(3, 3);
	double **M_inv = Memory_allocate_2D_double_array(3, 3);
	double *b;
	double eps_y, dilat_factor_y, delta_y;

	double *Vl = p -> Vl;

	double m_entries[5];
	int n, m, l, idx;
	int idx_y_min, idx_y_max;

	for (m=1; m< N_L_local; m++) {

		y_max = -1;
		y_min = 1024;
		x_max = -1;
		x_min = 1024;

		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				M[i][j] = 0;
			}
		}


		for (i=-1; i<2; i++) {

			if (fabs(y[idx_y_markers[m] + i] - y[idx_y_markers[m] + i - 1]) > y_max) {

				y_max = fabs(y[idx_y_markers[m] + i] - y[idx_y_markers[m] + i - 1]);
			}

			if (fabs(y[idx_y_markers[m] + i] - y[idx_y_markers[m] + i - 1]) < y_min) {

				y_min = fabs(y[idx_y_markers[m] + i] - y[idx_y_markers[m] + i - 1]);
			}
		}


		//-------------------------- calc dilatation factor -------------------
		eps_y = 0.03 * (y[idx_y_markers[m]] - y[idx_y_markers[m]-1]);
		dilat_factor_y = (5 / 6.0 * y_max + 1 / 6.0 * y_min + eps_y);
		//-------------------------------------------------------------------
		//--find idx of nodes inside the influence region of the current marker

		idx = idx_y_markers[m]-1;
		while ( y_markers[m] - y[idx]  <= dilat_factor_y * 1.5 ) idx = idx-1;
		idx_y_min = idx + 1;

		idx = idx_y_markers[m]+1;
		while ( y[idx] - y_markers[m] <= dilat_factor_y * 1.5 ) idx = idx+1;
		idx_y_max = idx - 1;


		for (i=0; i<6; i++) {
			m_entries[i] = 0.0;
		}

		for (l=idx_y_min; l<idx_y_max; l++) {

			delta_y = delta( (y_markers[m] - y[l]) / dilat_factor_y ) * (y[l+1] - y[l-1]) / 2;

			for (i=1; i<6; i++) {
				m_entries[i] = m_entries[i] + pow((y[l] - y_markers[m]), i) * delta_y;
            }

		}

		for (i=0; i<4; i++) {
			for (j=1;j<4;j++) {
				M[i][j] = m_entries[j+i];
			}
		}


		inv_3x3(M, M_inv);
		mat_vec(3, 3, M_inv, e1, b);


		for (k=idx_z_markers[m]-2; k < idx_z_markers[m]+3; k++) {
			for (l=idx_y_min; l<=idx_y_max; l++) {
				//dy = y(l) - y(l-1);
				for (n=idx_x_markers[m]-2; n<idx_x_markers[m]+2; n++) {

					a[k][l][n] = a[k][l][n] + A[m] *
					delta_mod_3D( x[k],x_markers[m], y[l], y_markers[m], z[n], z_markers[m], dx, dilat_factor_y, dz, b ) * Vl[m];

				}
			}
		}


	} // for m
}
#endif




/******************************************************************************/
/*
 Interpolate Eulerian quantity 'a' onto Lagrangian quantities 'A'.
 */
/******************************************************************************/
void Interpolate_Eul_to_Lag(double ***a, double *A, char which, Particle *p,
		MAC_grid *grid) {
#ifndef GRID_UNIFORM
	printf("Stop __func__ is onyl working for uniform grids" );
#endif
	int i, j, k;
	int mv, N_L_local;
	double r[3], h;
	double *X_L, *Y_L, *Z_L;

	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	double *x = grid -> xc;
	double *y = grid -> yc;
	double *z = grid -> zc;

	if (which == 'u') {
		x = grid -> xu;
	}
	else if (which == 'v') {
		y = grid -> yv;
	}
	else if (which == 'w') {
		z = grid -> zw;
	}

	N_L_local = p -> N_L_local;
	h = grid -> dx_u[1];

	X_L = p -> X_L;
	Y_L = p -> Y_L;
	Z_L = p -> Z_L;

#ifdef IBM_SCALAR
	if (which == 'H') {

	X_L = p -> X_H;
	Y_L = p -> Y_H;
	Z_L = p -> Z_H;
	}
#endif

	for (mv = 0; mv < N_L_local; mv++) {

        // TODO: Change bounds for non-uniform grid
		i_start = (int) round((X_L[mv] - x[0]) / h) - 1;
		j_start = (int) round((Y_L[mv] - y[0]) / h) - 1;
		k_start = (int) round((Z_L[mv] - z[0]) / h) - 1;
		i_end = i_start + 3;
		j_end = j_start + 3;
		k_end = k_start + 3;

		// Make sure processor bounds are not exceeded (pathological cases)
		i_start = max(i_start, grid->L_Is);
		j_start = max(j_start, grid->L_Js);
		k_start = max(k_start, grid->L_Ks);
		i_end = min(i_end, grid->L_Ie);
		j_end = min(j_end, grid->L_Je);
		k_end = min(k_end, grid->L_Ke);

#ifdef DEBUG  // TODO: Debugging
		if (i_start < grid -> L_Is) {
			printf("i_start = %d, L_Is = %d, G_Is = %d\n", i_start, grid -> L_Is, grid -> G_Is);
			printf("X_L[mv] = %f, h = %f\n", X_L[mv], h);
			printf("\t%f    %f\n", X_L[mv] + 1.5 * h, grid->xu[grid->G_Is]);
			printf("\t%f   %f   %f\n", x[i_start], X_L[mv], x[i_start+1]);
			printf("\t%f\n", (X_L[mv] - x[0]) / h - 1.5);
			printf("\t%f\n", round((X_L[mv] - x[0]) / h - 1.5));
		}
		if (j_start < grid -> L_Js)
			printf("j_start = %d, L_Js = %d\n", j_start, grid -> L_Js);
		if (k_start < grid -> L_Ks)
			printf("k_start = %d, L_Ks = %d\n", k_start, grid -> L_Ks);
		if (i_end > grid -> L_Ie) {
			printf("i_end = %d, L_Ie = %d\n", i_end, grid -> L_Ie);
			printf("\t%f\n\t%f\n", (X_L[mv] - x[0]) / h + 1.5, round((X_L[mv] - x[0]) / h + 1.5));
		}
		if (j_end > grid -> L_Je)
			printf("j_end = %d, L_Je = %d\n", j_end, grid -> L_Je);
		if (k_end > grid -> L_Ke)
			printf("k_end = %d, L_Ke = %d\n", k_end, grid -> L_Ke);
#endif

		//----------------------------------------------------------------------
		// Reset 'A'
		//----------------------------------------------------------------------
		A[mv] = 0.0;

		//----------------------------------------------------------------------
		// Interpolate
		//----------------------------------------------------------------------
		for (k = k_start; k < k_end; k++) {
			r[2] = Z_L[mv] - z[k];

			for (j = j_start; j < j_end; j++) {
				r[1] = Y_L[mv] - y[j];

				for (i = i_start; i < i_end; i++) {
					r[0] = X_L[mv] - x[i];

#ifdef GRID_UNIFORM
                    A[mv] += a[k][j][i] * interp_kernel(h, r);
#else
					A[mv] += a[k][j][i] * interp_kernel_nonuniform(grid->dx_u[i], r[0], grid->dy_v[j], r[1], grid->dz_w[k], r[2]);  // non-uniform
#endif

				}
			}
		}
	}
}




/******************************************************************************/
/*
 Integrate translational and rotational momentum over each mobile particle
 within the local processor domain and store in Int_U and Int_Omega

 Also calculates volume fraction, which should be reset to zero before calling
 this function for various particle lists.

 Assumes uniform grid
 */
/******************************************************************************/
void Interpolate_integrate_momentum(Velocity *vel, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i, j ,k, ii, jj, kk;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;
	int component, comp_pos, comp_neg;

	double r[3], temp, vf_cell, sum_phi;

	char message[100];

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;

	double *xc = grid -> xc;
	double *yc = grid -> yc;
	double *zc = grid -> zc;

	double *xe = grid -> xu;
	double *ye = grid -> yv;
	double *ze = grid -> zw;

	double ***vf;

	if (vel -> component == 'u') {
		xc = grid -> xu;
		xe = &(grid -> xc[-1]);
		component = 0;
		comp_pos  = 1;
		comp_neg  = 2;
		vf = lag -> ng_vfu;
	}
	else if (vel -> component == 'v') {
		yc = grid -> yv;
		ye = &(grid -> yc[-1]);
		component = 1;
		comp_pos  = 2;
		comp_neg  = 0;
		vf = lag -> ng_vfv;
	}
	else if (vel -> component == 'w') {
		zc = grid -> zw;
		ze = &(grid -> zc[-1]);
		component = 2;
		comp_pos  = 0;
		comp_neg  = 1;
		vf = lag -> ng_vfw;
	}
	else {
		sprintf(message, "Incorrect vel->component = '%c'", vel->component);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	// 'p_list' should have both local and foreign particles
	Display_assert_list_state(p_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	double ***phi  = lag -> temp;
	double ***data = vel -> data;

	// Reset volume fraction  // ASK ED why this line is here
	//double *vf_1d = &vf[grid->G_Ks][grid->G_Js][grid->G_Is];
	//DSET_ZERO(vf_1d, grid->ng_total_nodes);

	// Width and volume of grid cell
	double h = grid -> dx_u[1];
	double dV = h * h * h;

	double R, *X, *Int_U, *Int_Omega;

	Particle *p = p_list -> start;
	while (p != NULL) {

		R = p -> R;
		X = p -> X;
		Int_U     = p -> Int_U;
		Int_Omega = p -> Int_Omega;

		// ---------------------------------------------------------------------
		// Bounds of particle
		// ---------------------------------------------------------------------
		i_start = (int) floor( (X[0] - R - xe[0]) / h);
		j_start = (int) floor( (X[1] - R - ye[0]) / h);
		k_start = (int) floor( (X[2] - R - ze[0]) / h);

		i_end = (int) ceil( (X[0] + R - xe[0]) / h);
		j_end = (int) ceil( (X[1] + R - ye[0]) / h);
		k_end = (int) ceil( (X[2] + R - ze[0]) / h);

		i_start = max(i_start, grid -> G_Is);
		j_start = max(j_start, grid -> G_Js);
		k_start = max(k_start, grid -> G_Ks);

		i_end = min(i_end, grid -> G_Ie);
		j_end = min(j_end, grid -> G_Je);
		k_end = min(k_end, grid -> G_Ke);

		i_end = min(i_end, grid -> NX - 1);
		j_end = min(j_end, grid -> NY - 1);
		k_end = min(k_end, grid -> NZ - 1);

		// ---------------------------------------------------------------------
		// Computer level set function at all nodes
		// ---------------------------------------------------------------------
		for (k = k_start; k <= k_end; k++) {
			for (j = j_start; j <= j_end; j++) {
				for (i = i_start; i <= i_end; i++) {
					phi[k][j][i] = lvl_set(xe[i], ye[j], ze[k]);
				}
			}
		}

		// ---------------------------------------------------------------------
		// Compute integrals of velocity and vorticity
		// ---------------------------------------------------------------------
		for (k = k_start; k < k_end; k++) {
			r[2] = zc[k] - X[2];

			for (j = j_start; j < j_end; j++) {
				r[1] = yc[j] - X[1];

				for (i = i_start; i < i_end; i++) {
					r[0] = xc[i] - X[0];

					vf_cell = 0.0;
					sum_phi = 0.0;

					for (ii = 0; ii < 2; ii++) {
						for (jj = 0; jj < 2; jj++) {
							for (kk = 0; kk < 2; kk++) {

								temp = phi[k+kk][j+jj][i+ii];

								if (temp < 0.0)
									vf_cell -= temp;

								sum_phi += fabs(temp);
							}
						}
					}

					vf_cell = vf_cell / sum_phi;
					vf[k][j][i] += vf_cell;

					Int_U[component] += vf_cell * data[k][j][i] * dV;
					Int_Omega[comp_pos] += vf_cell * r[comp_neg] * data[k][j][i] * dV;
					Int_Omega[comp_neg] -= vf_cell * r[comp_pos] * data[k][j][i] * dV;
				}
			}
		}

		p = p -> next;
	}

}




/******************************************************************************/
/*
 * Add contribution of particles from 'p_list' to the volume fraction specified
 * by 'component' (either 'u', 'v', or 'w').  Volume fraction should be reset to
 * zero before calling this function for various particle lists.
 */
/******************************************************************************/

void Interpolate_add_to_volume_fraction(char component, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i, j ,k, ii, jj, kk;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	double temp, vf_cell, sum_phi;

	char message[100];

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;

	double *xc = grid -> xc;
	double *yc = grid -> yc;
	double *zc = grid -> zc;

	double *xe = grid -> xu;
	double *ye = grid -> yv;
	double *ze = grid -> zw;

	double ***vf;

	if (component == 'u') {
		xc = grid -> xu;
		xe = &(grid -> xc[-1]);
		vf = lag -> ng_vfu;
	}
	else if (component == 'v') {
		yc = grid -> yv;
		ye = &(grid -> yc[-1]);
		vf = lag -> ng_vfv;
	}
	else if (component == 'w') {
		zc = grid -> zw;
		ze = &(grid -> zc[-1]);
		vf = lag -> ng_vfw;
	}
#if (defined LAG_PARTICLE_RESOLVED)
	else if (component == 'c') {
		vf = lag -> ng_vfc;
	}
#endif
else if (component == 'z') {
		xc = grid -> xu;
		xe = &(grid -> xc[-1]);
		yc = grid -> yv;
		ye = &(grid -> yc[-1]);
		vf = lag -> ng_vfz;
	}
	else {
		sprintf(message, "Incorrect component = '%c'", component);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	// 'p_list' should have both local and foreign particles
	Display_assert_list_state(p_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	double ***phi  = lag -> temp;

	// Width of grid cell
	double h = grid -> dx_u[1];

	double R, *X;

	Particle *p = p_list -> start;
	while (p != NULL) {

		R = p -> R;
		X = p -> X;

		// ---------------------------------------------------------------------
		// Bounds of particle
		// ---------------------------------------------------------------------
		i_start = (int) floor( (X[0] - R - xe[0]) / h);
		j_start = (int) floor( (X[1] - R - ye[0]) / h);
		k_start = (int) floor( (X[2] - R - ze[0]) / h);

		i_end = (int) ceil( (X[0] + R - xe[0]) / h);
		j_end = (int) ceil( (X[1] + R - ye[0]) / h);
		k_end = (int) ceil( (X[2] + R - ze[0]) / h);

		i_start = max(i_start, grid -> G_Is);
		j_start = max(j_start, grid -> G_Js);
		k_start = max(k_start, grid -> G_Ks);

		i_end = min(i_end, grid -> G_Ie);
		j_end = min(j_end, grid -> G_Je);
		k_end = min(k_end, grid -> G_Ke);

		i_end = min(i_end, grid -> NX - 1);
		j_end = min(j_end, grid -> NY - 1);
		k_end = min(k_end, grid -> NZ - 1);

		// ---------------------------------------------------------------------
		// Computer level set function at all nodes
		// ---------------------------------------------------------------------
		for (k = k_start; k <= k_end; k++) {
			for (j = j_start; j <= j_end; j++) {
				for (i = i_start; i <= i_end; i++) {
					phi[k][j][i] = lvl_set(xe[i], ye[j], ze[k]);
				}
			}
		}

		// ---------------------------------------------------------------------
		// Compute integrals of velocity and vorticity
		// ---------------------------------------------------------------------
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {

					vf_cell = 0.0;
					sum_phi = 0.0;

					for (ii = 0; ii < 2; ii++) {
						for (jj = 0; jj < 2; jj++) {
							for (kk = 0; kk < 2; kk++) {

								temp = phi[k+kk][j+jj][i+ii];

								if (temp < 0.0)
									vf_cell -= temp;

								sum_phi += fabs(temp);
							}
						}
					}

					vf[k][j][i] += vf_cell / sum_phi;
				}
			}
		}

		p = p -> next;
	}

}


#ifdef VOF_SCALAR
void Interpolate_add_to_volume_fraction_vof(char component, Particle_list *p_list, Cart3d_bag *data_bag, Debug_trace *dtrace) {

	int i, j ,k, ii, jj, kk;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	double temp, vf_cell, sum_phi,eta;

	char message[100];

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;

	double *xc = grid -> xc;
	double *yc = grid -> yc;
	double *zc = grid -> zc;


	double *xe = grid -> xu;
	double *ye = grid -> yv;
	double *ze = grid -> zw;
	double ***vf;
	double ***vf_prime;


	double ***u_vof =  data_bag -> u -> data_vof;
	double ***v_vof =  data_bag -> v -> data_vof;
	double ***w_vof =  data_bag -> w -> data_vof;


	double velo_vof;




/* u- velocity */

	if (component == 'u') {
		xc = grid -> xu;
		xe = &(grid -> xc[-1]);
		vf = lag -> vfu;
#ifdef VOF_SMOOTH_VELO
		vf_prime =lag->vfu_prime;
#endif
	}

	else if (component == 'v') {
		yc = grid -> yv;
		ye = &(grid -> yc[-1]);
		vf = lag -> vfv;
#ifdef VOF_SMOOTH_VELO
		vf_prime =lag->vfv_prime;
#endif
	}
	else if (component == 'w') {
		zc = grid -> zw;
		ze = &(grid -> zc[-1]);
		vf = lag -> vfw;
#ifdef VOF_SMOOTH_VELO
		vf_prime =lag->vfw_prime;
#endif
	}
	else {
		sprintf(message, "Incorrect component = '%c'", component);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	// 'p_list' should have both local and foreign particles
	Display_assert_list_state(p_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	double ***phi  = lag -> temp;

	// Width of grid cell
	double h = grid -> dx_u[1];
	double delta_s=h/2;
	double R, *X, *U, *Omega;

	Particle *p = p_list -> start;
	while (p != NULL) {

		R = p -> R;
		X = p -> X;
		U = p -> U;
		Omega = p->Omega;

		// ---------------------------------------------------------------------
		// Bounds of particle
		// ---------------------------------------------------------------------
		i_start = (int) floor( (X[0] - R - xe[0]) / h);
		j_start = (int) floor( (X[1] - R - ye[0]) / h);
		k_start = (int) floor( (X[2] - R - ze[0]) / h);

		i_end = (int) ceil( (X[0] + R - xe[0]) / h);
		j_end = (int) ceil( (X[1] + R - ye[0]) / h);
		k_end = (int) ceil( (X[2] + R - ze[0]) / h);

		i_start = max(i_start, grid -> G_Is);
		j_start = max(j_start, grid -> G_Js);
		k_start = max(k_start, grid -> G_Ks);

		i_end = min(i_end, grid -> G_Ie);
		j_end = min(j_end, grid -> G_Je);
		k_end = min(k_end, grid -> G_Ke);



		// ---------------------------------------------------------------------
		// Computer level set function at all nodes
		// ---------------------------------------------------------------------
		for (k = k_start; k <= k_end; k++) {
			for (j = j_start; j <= j_end; j++) {
				for (i = i_start; i <= i_end; i++) {
					phi[k][j][i] = lvl_set(xe[i], ye[j], ze[k]);
				}
			}
		}


		if (component == 'u') {

		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {

					vf_cell = 0.0;
					sum_phi = 0.0;
					velo_vof= 0.0;

#ifdef VOF_NO_VOLUME
					eta= ETA(xc[i], yc[j], zc[k]);  // distance to the particle surface positive if inside
					vf[k][j][i] += SMOOTH(eta,delta_s);
					vf_prime[k][j][i] += SMOOTH(eta,delta_s);
					u_vof[k][j][i] += SMOOTH(eta,delta_s)*(U[0]+Omega[1]*(zc[k]-X[2])-Omega[2]*(yc[j]-X[1]));

#else

					for (ii = 0; ii < 2; ii++) {
						for (jj = 0; jj < 2; jj++) {
							for (kk = 0; kk < 2; kk++) {

								temp = phi[k+kk][j+jj][i+ii];

								if (temp < 0.0){  // if inside particle
									vf_cell -= temp;

									velo_vof -=temp*( U[0]+Omega[1]*(ze[k+kk]-X[2])-Omega[2]*(ye[j+jj]-X[1]) );
								}
								sum_phi += fabs(temp);
							}
						}
					}

					vf[k][j][i] += vf_cell / sum_phi;
	#ifdef VOF_SMOOTH_VELO
					eta= ETA(xc[i], yc[j], zc[k]);  // distance to the particle surface positive if inside
					vf_prime[k][j][i] += SMOOTH(eta,delta_s)*vf_cell / sum_phi;
					u_vof[k][j][i] += SMOOTH(eta,delta_s)*velo_vof / sum_phi;
	#else
					u_vof[k][j][i] += velo_vof / sum_phi;
	#endif

#endif
				}
			}
		}





		}// end if component u
		else if (component == 'v') {

			for (k = k_start; k < k_end; k++) {
				for (j = j_start; j < j_end; j++) {
					for (i = i_start; i < i_end; i++) {

						vf_cell = 0.0;
						sum_phi = 0.0;
						velo_vof= 0.0;

#ifdef VOF_NO_VOLUME

						eta= ETA(xc[i], yc[j], zc[k]);  // distance to the particle surface positive if inside
						vf[k][j][i] += SMOOTH(eta,delta_s);
						vf_prime[k][j][i] += SMOOTH(eta,delta_s);
						v_vof[k][j][i] += SMOOTH(eta,delta_s)*( U[1]+Omega[2]*(xc[i]-X[0])-Omega[0]*(zc[k]-X[2]) );

#else

						for (ii = 0; ii < 2; ii++) {
							for (jj = 0; jj < 2; jj++) {
								for (kk = 0; kk < 2; kk++) {

									temp = phi[k+kk][j+jj][i+ii];

									if (temp < 0.0){  // if inside particle
										vf_cell -= temp;
										velo_vof -=temp*( U[1]+Omega[2]*(xe[i+ii]-X[0])-Omega[0]*(ze[k+kk]-X[2]) );
									}
									sum_phi += fabs(temp);
								}
							}
						}

						vf[k][j][i] += vf_cell / sum_phi;
	#ifdef VOF_SMOOTH_VELO

						eta= ETA(xc[i], yc[j], zc[k]);  // distance to the particle surface positive if inside
						vf_prime[k][j][i] += SMOOTH(eta,delta_s)*vf_cell / sum_phi;
						v_vof[k][j][i] += SMOOTH(eta,delta_s)*velo_vof / sum_phi;
	#else
						v_vof[k][j][i] += velo_vof / sum_phi;
	#endif

#endif

					}
				}
			}





		}

		else if (component == 'w') {

			for (k = k_start; k < k_end; k++) {
							for (j = j_start; j < j_end; j++) {
								for (i = i_start; i < i_end; i++) {

									vf_cell = 0.0;
									sum_phi = 0.0;
									velo_vof= 0.0;

#ifdef VOF_NO_VOLUME

									eta= ETA(xc[i], yc[j], zc[k]);  // distance to the particle surface positive if inside
									vf[k][j][i] += SMOOTH(eta,delta_s);
									vf_prime[k][j][i] += SMOOTH(eta,delta_s);
									w_vof[k][j][i] += SMOOTH(eta,delta_s)*( U[2]+Omega[0]*(yc[j]-X[1])-Omega[1]*(xc[i]-X[0]) );


#else
									for (ii = 0; ii < 2; ii++) {
										for (jj = 0; jj < 2; jj++) {
											for (kk = 0; kk < 2; kk++) {

												temp = phi[k+kk][j+jj][i+ii];

												if (temp < 0.0){  // if inside particle
													vf_cell -= temp;
													velo_vof -=temp*( U[2]+Omega[0]*(ye[j+jj]-X[1])-Omega[1]*(xe[i+ii]-X[0]) );
												}
												sum_phi += fabs(temp);
											}
										}
									}
									vf[k][j][i] += vf_cell / sum_phi;
									#ifdef VOF_SMOOTH_VELO
									eta= ETA(xc[i], yc[j], zc[k]);  // distance to the particle surface positive if inside
									vf_prime[k][j][i] += SMOOTH(eta,delta_s)*vf_cell / sum_phi;
									w_vof[k][j][i] += SMOOTH(eta,delta_s)*velo_vof / sum_phi;
#else
									w_vof[k][j][i] += velo_vof / sum_phi;
#endif
#endif
								}
							}
						}

		}







		p = p -> next;
	}






return;
}
#endif

#ifdef VOF_SCALAR
void Interpolate_add_to_volume_fraction_vof_prime(char component, Particle_list *p_list, Cart3d_bag *data_bag, Debug_trace *dtrace) {


    int i, j ,k, ii, jj, kk;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;

	double temp, vf_cell, sum_phi,eta;

	char message[100];

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;
	Lagrangian *lag    = data_bag -> lag;

	double *xc = grid -> xc;
	double *yc = grid -> yc;
	double *zc = grid -> zc;


	double *xe = grid -> xu;
	double *ye = grid -> yv;
	double *ze = grid -> zw;
	double ***vf;
	double ***vf_prime;

    double bl_thick = params->bl_thick;



/* u- velocity */

	if (component == 'u') {
		xc = grid -> xu;
		xe = &(grid -> xc[-1]);
		vf = lag -> vfu;
#ifdef VOF_PP_VFPRIME
		vf_prime =lag->vfu_prime;
#endif
	}

	else if (component == 'v') {
		yc = grid -> yv;
		ye = &(grid -> yc[-1]);
		vf = lag -> vfv;
#ifdef VOF_PP_VFPRIME
		vf_prime =lag->vfv_prime;
#endif
	}
	else if (component == 'w') {
		zc = grid -> zw;
		ze = &(grid -> zc[-1]);
		vf = lag -> vfw;
#ifdef VOF_PP_VFPRIME
		vf_prime =lag->vfw_prime;
#endif
	}
	else {
		sprintf(message, "Incorrect component = '%c'", component);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}


	// 'p_list' should have both local and foreign particles
	Display_assert_list_state(p_list, LIST_STATE_BOTH, params, DTRACE("Display_assert_list_state"));

	double ***phi  = lag -> temp;

	// Width of grid cell
	double h = grid -> dx_u[1];
	double delta_s=h/2;
	double R, *X, *U, *Omega;

	Particle *p = p_list -> start;
	while (p != NULL) {

		R = p -> R;

        R = R + bl_thick;
		X = p -> X;
		U = p -> U;
		Omega = p->Omega;

		// ---------------------------------------------------------------------
		// Bounds of particle
		// ---------------------------------------------------------------------
		i_start = (int) floor( (X[0] - R - xe[0]) / h);
		j_start = (int) floor( (X[1] - R - ye[0]) / h);
		k_start = (int) floor( (X[2] - R - ze[0]) / h);

		i_end = (int) ceil( (X[0] + R - xe[0]) / h);
		j_end = (int) ceil( (X[1] + R - ye[0]) / h);
		k_end = (int) ceil( (X[2] + R - ze[0]) / h);

		i_start = max(i_start, grid -> G_Is);
		j_start = max(j_start, grid -> G_Js);
		k_start = max(k_start, grid -> G_Ks);

		i_end = min(i_end, grid -> G_Ie);
		j_end = min(j_end, grid -> G_Je);
		k_end = min(k_end, grid -> G_Ke);



		// ---------------------------------------------------------------------
		// Computer level set function at all nodes
		// ---------------------------------------------------------------------
		for (k = k_start; k <= k_end; k++) {
			for (j = j_start; j <= j_end; j++) {
				for (i = i_start; i <= i_end; i++) {
					phi[k][j][i] = lvl_set(xe[i], ye[j], ze[k]);
				}
			}
		}

        for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {

					vf_cell = 0.0;
					sum_phi = 0.0;

					for (ii = 0; ii < 2; ii++) {
						for (jj = 0; jj < 2; jj++) {
							for (kk = 0; kk < 2; kk++) {

								temp = phi[k+kk][j+jj][i+ii];

								if (temp < 0.0)
									vf_cell -= temp;

								sum_phi += fabs(temp);
							}
						}
					}

					vf_prime[k][j][i] += vf_cell / sum_phi;
				}
			}
		}


		p = p -> next;
	}

return;
}
#endif










/*
for (k = k_start; k < k_end; k++) {

	if (zc[k] < X[2]) {
		z2_near = ze[k+1];
		z2_far  = ze[k];
	}
	else {
		z2_near = ze[k];
		z2_far  = ze[k+1];
	}
	z2_near = z2_near * z2_near;
	z2_far  = z2_far  * z2_far;

	for (j = j_start; j < j_end; j++) {

		if (yc[j] < X[1]) {
			y2_near = ye[j+1];
			y2_far  = ye[j];
		}
		else {
			y2_near = ye[j];
			y2_far  = ye[j+1];
		}
		y2_near = y2_near * y2_near;
		y2_far  = y2_far  * y2_far;

		for (i = i_start; i < i_end; i++) {

			if (xc[i] < X[0]) {
				x2_near = xe[i+1];
				x2_far  = xe[i];
			}
			else {
				x2_near = xe[i];
				x2_far  = xe[i+1];
			}
			x2_near = x2_near * x2_near;
			x2_far  = x2_far  * x2_far;

			// Cell wholly within sphere
			if (x2_far + y2_far + z2_far <= R2) {
				vf = 1;
			}
			// Cell wholly without sphere
			else if (x2_near + y2_near + z2_near >= R2) {
				vf = 0;
			}
			// Cell partially within sphere
			else {

				vf = 0.0;
				sum_phi = 0.0;

				for (ii = 0; ii < 1; ii++) {
					for (jj = 0; jj < 1; jj++) {
						for (kk = 0; kk < 1; kk++) {

							phi = lvl_set(xe[i+ii], ye[j+jj], ze[k+kk]);

							if (phi > 0)
								vf += phi;

							sum_phi += phi;
						}
					}
				}

				vf = vf / sum_phi;
			}
		}
	}
}
*/


void Interpolate_bound_to_one( double ***vf, Cart3d_bag *data_bag) {
// used for the volume fraction that are used for the calculation of diffusion coefficients

	MAC_grid   *grid   = data_bag -> grid;

	int i,j,k;
	int i_start = grid -> G_Is;
	int j_start =  grid -> G_Js;
	int k_start = grid -> G_Ks;

	int i_end =  grid -> G_Ie;
	int j_end =  grid -> G_Je;
	int k_end =  grid -> G_Ke;




		// ---------------------------------------------------------------------
		// Bound the maximum value of the particle volume fraction to one here
		// ---------------------------------------------------------------------
		for (k = k_start; k <= k_end; k++) {
			for (j = j_start; j <= j_end; j++) {
				for (i = i_start; i <= i_end; i++) {

					vf[k][j][i]= min(1, vf[k][j][i]);

				}
			}
		}





}
