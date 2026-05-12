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

#if defined(LAG_PARTICLE_RESOLVED)

/* Geometry-only view used by the Eulerian solid-volume rasterization. */
typedef struct {
	double X[3];
	double R;
} Interp_solid_geom;

static Interp_solid_geom *
Interpolate_collect_owned_particle_geometry(Particle_list *p_list,
                                            Cart3d_bag *data_bag,
                                            int *n_local)
{
	int n_particles = 0;
	Particle *particles =
	    Particle_collect_owned_overlaps(p_list, data_bag, 0.0, -1.0,
	                                    1, &n_particles);

	if (n_particles == 0) {
		*n_local = 0;
		return NULL;
	}

	Interp_solid_geom *geom =
	    (Interp_solid_geom *)malloc(n_particles * sizeof(Interp_solid_geom));
	Memory_check_allocation(geom);

	for (int n = 0; n < n_particles; n++) {
		geom[n].X[0] = particles[n].X[0];
		geom[n].X[1] = particles[n].X[1];
		geom[n].X[2] = particles[n].X[2];
		geom[n].R = particles[n].R;
	}

	free(particles);
	*n_local = n_particles;
	return geom;
}

static inline double Interpolate_geom_level_set_3d(double x, double y, double z,
                                                   const Interp_solid_geom *g)
{
	const double dx = x - g->X[0];
	const double dy = y - g->X[1];
	const double dz = z - g->X[2];

#if defined(TWOD_CARTESIAN) || defined(AXISYM_RZ)
	return sqrt(dx * dx + dy * dy) / (g->R + 1.0e-30) - 1.0;
#else
	return sqrt(dx * dx + dy * dy + dz * dz) / (g->R + 1.0e-30) - 1.0;
#endif
}

static inline double Interpolate_volume_fraction_from_phi_3d(double phi[2][2][2])
{
	double vf_cell = 0.0;
	double sum_phi = 0.0;

	for (int kk = 0; kk < 2; ++kk)
	for (int jj = 0; jj < 2; ++jj)
	for (int ii = 0; ii < 2; ++ii) {
		const double temp = phi[kk][jj][ii];
		if (temp < 0.0)
			vf_cell -= temp;
		sum_phi += fabs(temp);
	}

	return (sum_phi > 0.0) ? vf_cell / sum_phi : 0.0;
}

static inline double Interpolate_volume_fraction_from_phi_2d(double phi00,
                                                             double phi10,
                                                             double phi01,
                                                             double phi11)
{
	double vf_cell = 0.0;
	double sum_phi = fabs(phi00) + fabs(phi10) + fabs(phi01) + fabs(phi11);

	if (phi00 < 0.0) vf_cell -= phi00;
	if (phi10 < 0.0) vf_cell -= phi10;
	if (phi01 < 0.0) vf_cell -= phi01;
	if (phi11 < 0.0) vf_cell -= phi11;

	return (sum_phi > 0.0) ? vf_cell / sum_phi : 0.0;
}

static void Interpolate_add_to_volume_fraction_global(char component,
                                                      Particle_list *p_list,
                                                      Cart3d_bag *data_bag,
                                                      Debug_trace *dtrace)
{
	(void)dtrace;

	Parameters *params = data_bag->params;
	MAC_grid   *grid   = data_bag->grid;
	Lagrangian *lag    = data_bag->lag;

	double ***vf = NULL;

	double *xe = grid->xu;
	double *ye = grid->yv;
	double *ze = grid->zw;

	if (component == 'u') {
		xe = &(grid->xc[-1]);
		vf = lag->ng_vfu;
	}
	else if (component == 'v') {
		ye = &(grid->yc[-1]);
		vf = lag->ng_vfv;
	}
	else if (component == 'w') {
#if defined(TWOD_MODE)
		return;
#else
		ze = &(grid->zc[-1]);
		vf = lag->ng_vfw;
#endif
	}
	else if (component == 'c') {
		vf = lag->ng_vfc;
	}
	else if (component == 'z') {
		xe = &(grid->xc[-1]);
		ye = &(grid->yc[-1]);
		vf = lag->ng_vfz;
	}
	else {
		char message[100];
		sprintf(message, "Incorrect component = '%c'", component);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	int n_geom = 0;
	Interp_solid_geom *geom =
	    Interpolate_collect_owned_particle_geometry(p_list, data_bag, &n_geom);

	const double h = grid->dx_u[1];

#if defined(TWOD_MODE)

	const int k = grid->G_Ks;

	for (int n = 0; n < n_geom; ++n) {
		const double R = geom[n].R;
		const double *X = geom[n].X;

		int i_start = (int)floor((X[0] - R - xe[0]) / h);
		int j_start = (int)floor((X[1] - R - ye[0]) / h);
		int i_end   = (int)ceil ((X[0] + R - xe[0]) / h);
		int j_end   = (int)ceil ((X[1] + R - ye[0]) / h);

		i_start = max(i_start, grid->G_Is);
		j_start = max(j_start, grid->G_Js);

		i_end = min(min(i_end, grid->G_Ie), grid->NX - 1);
		j_end = min(min(j_end, grid->G_Je), grid->NY - 1);

		for (int j = j_start; j < j_end; ++j) {
			for (int i = i_start; i < i_end; ++i) {
				const double phi00 =
				    Interpolate_geom_level_set_3d(xe[i],     ye[j],     0.0, &geom[n]);
				const double phi10 =
				    Interpolate_geom_level_set_3d(xe[i + 1], ye[j],     0.0, &geom[n]);
				const double phi01 =
				    Interpolate_geom_level_set_3d(xe[i],     ye[j + 1], 0.0, &geom[n]);
				const double phi11 =
				    Interpolate_geom_level_set_3d(xe[i + 1], ye[j + 1], 0.0, &geom[n]);

				vf[k][j][i] += Interpolate_volume_fraction_from_phi_2d(
					phi00, phi10, phi01, phi11);
			}
		}
	}

#else

	for (int n = 0; n < n_geom; ++n) {
		const double R = geom[n].R;
		const double *X = geom[n].X;

		int i_start = (int)floor((X[0] - R - xe[0]) / h);
		int j_start = (int)floor((X[1] - R - ye[0]) / h);
		int k_start = (int)floor((X[2] - R - ze[0]) / h);

		int i_end = (int)ceil((X[0] + R - xe[0]) / h);
		int j_end = (int)ceil((X[1] + R - ye[0]) / h);
		int k_end = (int)ceil((X[2] + R - ze[0]) / h);

		i_start = max(i_start, grid->G_Is);
		j_start = max(j_start, grid->G_Js);
		k_start = max(k_start, grid->G_Ks);

		i_end = min(i_end, grid->G_Ie);
		j_end = min(j_end, grid->G_Je);
		k_end = min(k_end, grid->G_Ke);

		i_end = min(i_end, grid->NX - 1);
		j_end = min(j_end, grid->NY - 1);
		k_end = min(k_end, grid->NZ - 1);

		for (int k = k_start; k < k_end; ++k) {
			for (int j = j_start; j < j_end; ++j) {
				for (int i = i_start; i < i_end; ++i) {
					double phi[2][2][2];

					for (int kk = 0; kk < 2; ++kk)
					for (int jj = 0; jj < 2; ++jj)
					for (int ii = 0; ii < 2; ++ii) {
						phi[kk][jj][ii] =
						    Interpolate_geom_level_set_3d(xe[i + ii],
						                                  ye[j + jj],
						                                  ze[k + kk],
						                                  &geom[n]);
					}

					vf[k][j][i] += Interpolate_volume_fraction_from_phi_3d(phi);
				}
			}
		}
	}

#endif

	free(geom);
}

#endif /* LAG_PARTICLE_RESOLVED */

double delta(double r);

#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
static double interp_kernel_2d(double h, const double *r)
{
	return delta(r[0] / h) * delta(r[1] / h);
}

static double Interpolate_twod_marker_volume(const Particle *p, int mv)
{
	if (p->Vol_L_marker != NULL)
		return p->Vol_L_marker[mv];
	return p->Vol_L;
}

static double Interpolate_twod_cell_volume(const MAC_grid *grid, char which,
                                           const double *x, int i, int j)
{
#ifdef AXISYM_RZ
	/*
	 * A 2D axisymmetric cell represents a full ring.  For u-faces the control
	 * volume is centered on r_u; for all scalar/v/w-style quantities we use the
	 * local interpolation coordinate, which is xc for the physical branches.
	 */
	double radius = (which == 'u') ? fabs(grid->xu[i]) : fabs(x[i]);
	if (radius < TWOD_RADIAL_EPS)
		radius = 0.5 * grid->dx_u[1];
	double ring_weight = (which == 'u' && grid->ring_wt_u != NULL) ?
	                     grid->ring_wt_u[i] :
	                     ((grid->ring_wt_c != NULL) ? grid->ring_wt_c[i] :
	                      TWOD_AXISYM_THETA_SPAN_FULL * radius);
	if (ring_weight < TWOD_RADIAL_EPS)
		ring_weight = TWOD_AXISYM_THETA_SPAN_FULL * radius;
	return grid->dx_u[1] * grid->dy_v[j] * ring_weight;
#else
	double slab = grid->dummy_z_slab_thickness;
	if (slab <= 0.0)
		slab = grid->dx_u[1];
	return grid->dx_u[1] * grid->dy_v[j] * slab;
#endif
}

static double Interpolate_twod_level_set(double x, double y, const Particle *p)
{
	const double *X = p->X;
	const double R = p->R;
	return sqrt((x - X[0]) * (x - X[0]) +
	            (y - X[1]) * (y - X[1])) / R - 1.0;
}

/*
 * Index parity helpers across the left x-boundary.
 *
 * AXISYM_RZ: the left boundary is the symmetry axis r = 0.  u-faces have
 * radial parity (odd), all scalar/v-style quantities have even parity.
 *
 * TWOD_CARTESIAN with LEFT_WALL_VELOCITY_FREESLIP: a half-domain whose left
 * edge is a symmetry plane.  Markers spawned on the negative-x side of the
 * plane fold back into the physical half-domain; u-faces still keep odd
 * parity so the normal velocity at the symmetry face stays zero.  This split
 * matches the "Left-symmetry IBM strategy" rule in the roadmap.
 */
static int Interpolate_axis_index_for_spread(int i, char which,
                                             const MAC_grid *grid, int *sign)
{
	*sign = 1;
#ifdef AXISYM_RZ
	if (grid->G_Is == 0 && i < grid->G_Is) {
		if (which == 'u') {
			*sign = -1;
			return grid->G_Is + (grid->G_Is - i);
		}
		return grid->G_Is + (grid->G_Is - 1 - i);
	}
	if (grid->G_Is == 0 && which == 'u' && i == grid->G_Is)
		return -1;
#elif defined(TWOD_CARTESIAN) && defined(LEFT_WALL_VELOCITY_FREESLIP)
	if (grid->G_Is == 0 && i < grid->G_Is) {
		if (which == 'u') {
			*sign = -1;
			return grid->G_Is + (grid->G_Is - i);
		}
		return grid->G_Is + (grid->G_Is - 1 - i);
	}
	if (grid->G_Is == 0 && which == 'u' && i == grid->G_Is)
		return -1;
#endif
	return i;
}

static int Interpolate_axis_index_for_sample(int i, char which,
                                             const MAC_grid *grid, int *sign)
{
	*sign = 1;
#ifdef AXISYM_RZ
	if (grid->G_Is == 0 && i < grid->G_Is) {
		if (which == 'u') {
			*sign = -1;
			return grid->G_Is + (grid->G_Is - i);
		}
		return grid->G_Is + (grid->G_Is - 1 - i);
	}
	if (grid->G_Is == 0 && which == 'u' && i == grid->G_Is)
		return -1;
#elif defined(TWOD_CARTESIAN) && defined(LEFT_WALL_VELOCITY_FREESLIP)
	if (grid->G_Is == 0 && i < grid->G_Is) {
		if (which == 'u') {
			*sign = -1;
			return grid->G_Is + (grid->G_Is - i);
		}
		return grid->G_Is + (grid->G_Is - 1 - i);
	}
	if (grid->G_Is == 0 && which == 'u' && i == grid->G_Is)
		return -1;
#endif
	return i;
}

static void Interpolate_Lag_to_Eul_uniform_twod(double *A, double ***a,
        char which, Particle *p, MAC_grid *grid)
{
	int i, j, mv;
	double r[3];
	double *x = grid->xc;
	double *y = grid->yc;

	if (which == 'u')
		x = grid->xu;
	else if (which == 'v')
		y = grid->yv;

	double *a_1d = &a[grid->G_Ks][grid->G_Js][grid->G_Is];
	for (i = 0; i < grid->ng_total_nodes; i++)
		a_1d[i] = 0.0;

	if (which == 'w')
		return;

	/*
	 * TWOD_MODE:
	 * The z direction is storage-only.  The Lagrangian marker volume already
	 * contains the full physical slab/ring measure:
	 *
	 *   TWOD_CARTESIAN: dV_L = ds * h * slab_thickness
	 *   AXISYM_RZ:     dV_L = theta_span * r * ds * h
	 *
	 * Therefore the marker force must be spread to one physical k-plane only.
	 * Looping over all stored/periodic k planes duplicates the IBM force.
	 */
	const int k = grid->G_Ks;

	const double h = grid->dx_u[1];

	double *X_L = p->X_L;
	double *Y_L = p->Y_L;

	for (mv = 0; mv < p->N_L_local; mv++) {
		int i_start = (int)round((X_L[mv] - x[0]) / h) - 1;
		int j_start = (int)round((Y_L[mv] - y[0]) / h) - 1;
		int i_end = i_start + 3;
		int j_end = j_start + 3;

#ifdef AXISYM_RZ
		if (grid->G_Is == 0)
			i_start = max(i_start, grid->L_Is);
		else
			i_start = max(i_start, grid->G_Is);
#elif defined(TWOD_CARTESIAN) && defined(LEFT_WALL_VELOCITY_FREESLIP)
		if (grid->G_Is == 0)
			i_start = max(i_start, grid->L_Is);
		else
			i_start = max(i_start, grid->G_Is);
#else
		i_start = max(i_start, grid->G_Is);
#endif

		i_end = min(max(i_end, grid->G_Is), grid->G_Ie);
		j_start = min(max(j_start, grid->G_Js), grid->G_Je);
		j_end = min(max(j_end, grid->G_Js), grid->G_Je);

		for (j = j_start; j < j_end; j++) {
			r[1] = Y_L[mv] - y[j];

			for (i = i_start; i < i_end; i++) {
				int parity;
				int ii = Interpolate_axis_index_for_spread(i, which, grid, &parity);

				if (ii < grid->G_Is || ii >= grid->G_Ie)
					continue;

				r[0] = X_L[mv] - x[i];

				double value = (A == NULL) ? 1.0 : A[mv];

				double dV =
					Interpolate_twod_marker_volume(p, mv) /
					Interpolate_twod_cell_volume(grid, which, x, ii, j);

				a[k][j][ii] +=
					parity * value * interp_kernel_2d(h, r) * dV;
			}
		}
	}
}

static void Interpolate_Eul_to_Lag_twod(double ***a, double *A, char which,
                                        Particle *p, MAC_grid *grid)
{
	int i, j, mv;
	double r[3];
	double *x = grid->xc;
	double *y = grid->yc;
	double *X_L = p->X_L;
	double *Y_L = p->Y_L;

	if (which == 'u')
		x = grid->xu;
	else if (which == 'v')
		y = grid->yv;
	else if (which == 'w') {
		for (mv = 0; mv < p->N_L_local; mv++)
			A[mv] = 0.0;
		return;
	}

#ifdef IBM_SCALAR
	if (which == 'H') {
		X_L = p->X_H;
		Y_L = p->Y_H;
	}
#endif

	const double h = grid->dx_u[1];
	const int k = grid->G_Ks;
	for (mv = 0; mv < p->N_L_local; mv++) {
		int i_start = (int)round((X_L[mv] - x[0]) / h) - 1;
		int j_start = (int)round((Y_L[mv] - y[0]) / h) - 1;
		int i_end = i_start + 3;
		int j_end = j_start + 3;

		i_start = max(i_start, grid->L_Is);
		j_start = max(j_start, grid->L_Js);
		i_end = min(i_end, grid->L_Ie);
		j_end = min(j_end, grid->L_Je);

		A[mv] = 0.0;
		for (j = j_start; j < j_end; j++) {
			r[1] = Y_L[mv] - y[j];
			for (i = i_start; i < i_end; i++) {
				int parity, ii = Interpolate_axis_index_for_sample(i, which, grid, &parity);
				if (ii < grid->L_Is || ii >= grid->L_Ie)
					continue;
				r[0] = X_L[mv] - x[i];
				A[mv] += parity * a[k][j][ii] * interp_kernel_2d(h, r);
			}
		}
	}
}

static double Interpolate_twod_volume_fraction_from_phi(double phi00,
                                                        double phi10,
                                                        double phi01,
                                                        double phi11)
{
	double vf_cell = 0.0;
	double sum_phi = fabs(phi00) + fabs(phi10) + fabs(phi01) + fabs(phi11);
	if (phi00 < 0.0) vf_cell -= phi00;
	if (phi10 < 0.0) vf_cell -= phi10;
	if (phi01 < 0.0) vf_cell -= phi01;
	if (phi11 < 0.0) vf_cell -= phi11;
	return (sum_phi > 0.0) ? vf_cell / sum_phi : 0.0;
}
#endif


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

	return delta(rx/hx) * delta(ry/hy) * delta(rz/hz);
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
	M_inv[2][2] = (-M[0][1]*M[0][1] + M[0][0]*M[0][2]) / det;

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


#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
static void Interpolate_integrate_momentum_twod(Velocity *vel,
                                                Particle_list *p_list,
                                                Cart3d_bag *data_bag,
                                                Debug_trace *dtrace)
{
	int i, j, ii, jj;
	int component = -1;
	double *xc = data_bag->grid->xc;
	double *yc = data_bag->grid->yc;
	double *xe = data_bag->grid->xu;
	double *ye = data_bag->grid->yv;
	double ***vf = NULL;

	Parameters *params = data_bag->params;
	MAC_grid *grid = data_bag->grid;
	Lagrangian *lag = data_bag->lag;
	double ***data = vel->data;

	if (vel->component == 'u') {
		xc = grid->xu;
		xe = &(grid->xc[-1]);
		component = 0;
		vf = lag->ng_vfu;
	} else if (vel->component == 'v') {
		yc = grid->yv;
		ye = &(grid->yc[-1]);
		component = 1;
		vf = lag->ng_vfv;
	} else if (vel->component == 'w') {
		return; // No physical momentum lives in the collapsed direction.
	}

	Display_assert_list_state(p_list, LIST_STATE_BOTH,
	                          params, DTRACE("Display_assert_list_state"));

	const double h = grid->dx_u[1];
	const int k = grid->G_Ks;
	Particle *p = p_list->start;
	while (p != NULL) {
		const double R = p->R;
		const double *X = p->X;

		int i_start = (int)floor((X[0] - R - xe[0]) / h);
		int j_start = (int)floor((X[1] - R - ye[0]) / h);
		int i_end = (int)ceil((X[0] + R - xe[0]) / h);
		int j_end = (int)ceil((X[1] + R - ye[0]) / h);

		i_start = max(i_start, grid->G_Is);
		j_start = max(j_start, grid->G_Js);
		i_end = min(min(i_end, grid->G_Ie), grid->NX - 1);
		j_end = min(min(j_end, grid->G_Je), grid->NY - 1);

		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
				double phi[2][2];
				for (jj = 0; jj < 2; jj++)
					for (ii = 0; ii < 2; ii++)
						phi[jj][ii] = Interpolate_twod_level_set(xe[i + ii],
						                                          ye[j + jj], p);

				double vf_cell = Interpolate_twod_volume_fraction_from_phi(
					phi[0][0], phi[0][1], phi[1][0], phi[1][1]);
				vf[k][j][i] += vf_cell;

#ifdef VOF_IBM
				VolumeFraction *vof = data_bag->vof;
				if (component == 0) {
					double dV_scalar = Interpolate_twod_cell_volume(grid, 'c',
					                                                grid->xc, i, j);
					p->Int_rho_scalar += lag->ng_vfc[k][j][i] * dV_scalar;
				}
#else
				const double rho_loc = 1.0;
#endif
				double r[2] = {xc[i] - X[0], yc[j] - X[1]};
				double dV = Interpolate_twod_cell_volume(grid, vel->component,
				                                         xc, i, j);
				double dP = vf_cell * data[k][j][i] * dV;

#ifdef VOF_IBM
				p->Int_rho[component] += vf_cell * dV;
#endif
				p->Int_U[component] += dP;
#ifdef TWOD_CARTESIAN
				/* Cylinder rotation is about the collapsed z axis only. */
				if (component == 0)
					p->Int_Omega[2] -= dP * r[1];
				else
					p->Int_Omega[2] += dP * r[0];
#endif
			}
		}
		p = p->next;
	}
}

/*
 * Blend rigid-body translation+rotation onto the staggered velocity face
 * inside each particle's level-set support.  Mirrors
 * Interpolate_integrate_momentum_twod's geometry stencil so that a subsequent
 * re-integration sees Int_U == Int_rho * U_rigid for a quiescent exterior,
 * eliminating the cold-start IBM impulse on the first time step.
 */
static void Interpolate_paint_rigid_body_velocity_twod(Velocity *vel,
                                                       Particle_list *p_list,
                                                       Cart3d_bag *data_bag,
                                                       Debug_trace *dtrace)
{
	int i, j, ii, jj;
	int component = -1;
	double *xc = data_bag->grid->xc;
	double *yc = data_bag->grid->yc;
	double *xe = data_bag->grid->xu;
	double *ye = data_bag->grid->yv;

	Parameters *params = data_bag->params;
	MAC_grid *grid = data_bag->grid;
	double ***data = vel->data;

	if (vel->component == 'u') {
		xc = grid->xu;
		xe = &(grid->xc[-1]);
		component = 0;
	} else if (vel->component == 'v') {
		yc = grid->yv;
		ye = &(grid->yc[-1]);
		component = 1;
	} else if (vel->component == 'w') {
		/* Storage-only direction in TWOD_MODE; integrator skips it as well. */
		return;
	}

	Display_assert_list_state(p_list, LIST_STATE_BOTH,
	                          params, DTRACE("Display_assert_list_state"));

	const double h = grid->dx_u[1];
	const int k = grid->G_Ks;
	Particle *p = p_list->start;
	while (p != NULL) {
		const double R = p->R;
		const double *X = p->X;
		const double *Up = p->U;
		const double Omega_z = p->Omega[2]; /* in-plane rotation only */

		int i_start = (int)floor((X[0] - R - xe[0]) / h);
		int j_start = (int)floor((X[1] - R - ye[0]) / h);
		int i_end = (int)ceil((X[0] + R - xe[0]) / h);
		int j_end = (int)ceil((X[1] + R - ye[0]) / h);

		i_start = max(i_start, grid->G_Is);
		j_start = max(j_start, grid->G_Js);
		i_end = min(min(i_end, grid->G_Ie), grid->NX - 1);
		j_end = min(min(j_end, grid->G_Je), grid->NY - 1);

		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
				double phi[2][2];
				for (jj = 0; jj < 2; jj++)
					for (ii = 0; ii < 2; ii++)
						phi[jj][ii] = Interpolate_twod_level_set(xe[i + ii],
						                                          ye[j + jj], p);

				double vf_cell = Interpolate_twod_volume_fraction_from_phi(
					phi[0][0], phi[0][1], phi[1][0], phi[1][1]);
				if (vf_cell <= 0.0)
					continue;

				/* Rigid-body face velocity: u = U + Omega x r.  In TWOD_MODE
				 * only Omega_z (out-of-plane) is meaningful; AXISYM_RZ
				 * additionally forces Omega = 0, so the rotation term is a
				 * no-op there. */
				double u_rigid;
				if (component == 0)
					u_rigid = Up[0] - Omega_z * (yc[j] - X[1]);
				else
					u_rigid = Up[1] + Omega_z * (xc[i] - X[0]);

				data[k][j][i] = vf_cell * u_rigid +
				                (1.0 - vf_cell) * data[k][j][i];
			}
		}
		p = p->next;
	}
}

static void Interpolate_add_to_volume_fraction_twod(char component,
                                                    Particle_list *p_list,
                                                    Cart3d_bag *data_bag,
                                                    Debug_trace *dtrace)
{
	int i, j, ii, jj;
	double *xe = data_bag->grid->xu;
	double *ye = data_bag->grid->yv;
	double ***vf = NULL;

	Parameters *params = data_bag->params;
	MAC_grid *grid = data_bag->grid;
	Lagrangian *lag = data_bag->lag;

	if (component == 'u') {
		xe = &(grid->xc[-1]);
		vf = lag->ng_vfu;
	} else if (component == 'v') {
		ye = &(grid->yc[-1]);
		vf = lag->ng_vfv;
	} else if (component == 'w') {
		return;
	}
#if defined LAG_PARTICLE_RESOLVED
	else if (component == 'c') {
		vf = lag->ng_vfc;
	}
#endif
	else if (component == 'z') {
		xe = &(grid->xc[-1]);
		ye = &(grid->yc[-1]);
		vf = lag->ng_vfz;
	} else {
		char message[100];
		sprintf(message, "Incorrect component = '%c'", component);
		Display_throw_error(message, params, DTRACE("Display_throw_error"));
	}

	Display_assert_list_state(p_list, LIST_STATE_BOTH,
	                          params, DTRACE("Display_assert_list_state"));

	const double h = grid->dx_u[1];
	const int k = grid->G_Ks;
	Particle *p = p_list->start;
	while (p != NULL) {
		const double R = p->R;
		const double *X = p->X;
		int i_start = (int)floor((X[0] - R - xe[0]) / h);
		int j_start = (int)floor((X[1] - R - ye[0]) / h);
		int i_end = (int)ceil((X[0] + R - xe[0]) / h);
		int j_end = (int)ceil((X[1] + R - ye[0]) / h);

		i_start = max(i_start, grid->G_Is);
		j_start = max(j_start, grid->G_Js);
		i_end = min(min(i_end, grid->G_Ie), grid->NX - 1);
		j_end = min(min(j_end, grid->G_Je), grid->NY - 1);

		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {
				double phi[2][2];
				for (jj = 0; jj < 2; jj++)
					for (ii = 0; ii < 2; ii++)
						phi[jj][ii] = Interpolate_twod_level_set(xe[i + ii],
						                                          ye[j + jj], p);
				vf[k][j][i] += Interpolate_twod_volume_fraction_from_phi(
					phi[0][0], phi[0][1], phi[1][0], phi[1][1]);
			}
		}
		p = p->next;
	}
}
#endif



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
#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
    Interpolate_Lag_to_Eul_uniform_twod(A, a, which, p, grid);
    return;
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
#ifdef LEFT_WALL_VELOCITY_FREESLIP
            i_start = max(i_start, grid->L_Is);
            i_end = max(i_end, grid->L_Is); 
#else
            i_start = max(i_start, grid->G_Is);
            i_end = max(i_end, grid->G_Is);
#endif
            i_start = min(i_start, grid->G_Ie);
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

#ifdef LEFT_WALL_VELOCITY_FREESLIP
				double force = A[mv] * interp_kernel(h, r) * dVol;
				if (i < grid->G_Is) { 
					if (which == 'u') { // Odd parity
						int i_mirror = grid->G_Is + (grid->G_Is - i);
						a[k][j][i_mirror] -= force;
					} else { // Even parity
						int i_mirror = grid->G_Is + (grid->G_Is - 1 - i);
						a[k][j][i_mirror] += force;
					}
				} else if (which == 'u' && i == grid->G_Is) {
					// DO NOTHING! Normal velocity at the symmetry plane is strictly zero.
					// Forcing here fights the projection solver and explodes pressure.
				} else {
					a[k][j][i] += force;
				}
#else
                        a[k][j][i] += interp_kernel(h, r) * dVol;
#endif
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
#ifdef LEFT_WALL_VELOCITY_FREESLIP
            i_start = max(i_start, grid->L_Is);
            i_end = max(i_end, grid->L_Is);
#else
            i_start = max(i_start, grid->G_Is);
            i_end = max(i_end, grid->G_Is);
#endif
            i_start = min(i_start, grid->G_Ie);
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

#ifdef LEFT_WALL_VELOCITY_FREESLIP
			double force = A[mv] * interp_kernel(h, r) * dVol;
			if (i < grid->G_Is) { 
				if (which == 'u') { // Odd parity
					int i_mirror = grid->G_Is + (grid->G_Is - i);
					a[k][j][i_mirror] -= force;
				} else { // Even parity
					int i_mirror = grid->G_Is + (grid->G_Is - 1 - i);
					a[k][j][i_mirror] += force;
				}
			} else if (which == 'u' && i == grid->G_Is) {
				// DO NOTHING! Normal velocity at the symmetry plane is strictly zero.
				// Forcing here fights the projection solver and explodes pressure.
			} else {
				a[k][j][i] += force;
			}
#else
                        a[k][j][i] += A[mv] * interp_kernel(h, r) * dVol;
#endif
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

	double m_entries[7];
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
#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
	Interpolate_Eul_to_Lag_twod(a, A, which, p, grid);
	return;
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
	else if (which == 'c') {
		x = grid -> xc;
		y = grid -> yc;
		z = grid -> zc;
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
 * Interpolate_integrate_momentum – VOF-IBM aware
 *
 *  • Integrates translational (Int_U) and rotational (Int_Omega) momentum of
 *    every (mobile) particle in p_list that overlaps the current MPI sub-domain.
 *  • Accumulates the particle volume fraction on the staggered velocity faces
 *    (lag->ng_vfu | ng_vfv | ng_vfw).
 *
 *  ───────────────────────────────────────────────────────────────────────────
 *  NEW (when VOF_IBM is defined)
 *  ───────────────────────────────────────────────────────────────────────────
 *      ‣ Each staggered control-volume contribution to momentum is multiplied
 *        by the local mixture density ρ̃(F)=F+(1-F)(ρ₂/ρ₁) taken from vof->rho.
 *      ‣ Face-centred density for momentum is obtained by arithmetic averaging
 *        of the two adjacent cell-centred values.
 *      ‣ Particle integrated density (p->Int_rho[component]) uses face-centered
 *        density matching the stencil used for Int_U[component]
 *      ‣ When VOF_IBM is *not* defined the routine reverts to the original
 *        single-phase (ρ̃=1) behaviour with zero additional cost.
 *
 *  Assumes a uniform Cartesian grid of spacing h.
 */
/******************************************************************************/
void Interpolate_integrate_momentum(Velocity      *vel,
                                    Particle_list *p_list,
                                    Cart3d_bag    *data_bag,
                                    Debug_trace   *dtrace)
{
#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
    Interpolate_integrate_momentum_twod(vel, p_list, data_bag, dtrace);
    return;
#endif
    int i, j, k, ii, jj, kk;
    int i_start, j_start, k_start;
    int i_end,   j_end,   k_end;
    int component, comp_pos, comp_neg;

    double r[3], temp, vf_cell, sum_phi;

    char message[100];

    /* --------------------------------------------------------------------- */
    /* Pointers to global data                                               */
    /* --------------------------------------------------------------------- */
    Parameters     *params = data_bag->params;
    MAC_grid       *grid   = data_bag->grid;
    Lagrangian     *lag    = data_bag->lag;

#ifdef VOF_IBM
    /* Density field produced by VOF_update_density_viscosity() */
    VolumeFraction *vof       = data_bag->vof;
    double       ***rho_cc    = vof->rho;   /* cell-centred ρ̃ */
#endif

    /* Cell- and face-centred coordinates (will re-map below) */
    double *xc = grid->xc, *yc = grid->yc, *zc = grid->zc;
    double *xe = grid->xu, *ye = grid->yv, *ze = grid->zw;

    /* Face-type-dependent pointers                                         */
    double ***vf;                      /* face-solid volume fraction field */

    if (vel->component == 'u') {       /* ------- u-velocity faces -------- */
        xc        = grid->xu;
        xe        = &(grid->xc[-1]);   /* left cell centres */
        component = 0;                 /* Int_U[x]          */
        comp_pos  = 1;                 /* Int_Omega[y]      */
        comp_neg  = 2;                 /* Int_Omega[z]      */
        vf        = lag->ng_vfu;
    }
    else if (vel->component == 'v') {  /* ------- v-velocity faces -------- */
        yc        = grid->yv;
        ye        = &(grid->yc[-1]);
        component = 1;
        comp_pos  = 2;
        comp_neg  = 0;
        vf        = lag->ng_vfv;
    }
    else if (vel->component == 'w') {  /* ------- w-velocity faces -------- */
        zc        = grid->zw;
        ze        = &(grid->zc[-1]);
        component = 2;
        comp_pos  = 0;
        comp_neg  = 1;
        vf        = lag->ng_vfw;
    }
    else {
        sprintf(message, "Incorrect vel->component = '%c'", vel->component);
        Display_throw_error(message, params, DTRACE("Display_throw_error"));
    }

    /* Ensure list contains local + ghost particles                          */
    Display_assert_list_state(p_list, LIST_STATE_BOTH,
                              params, DTRACE("Display_assert_list_state"));

    /* --------------------------------------------------------------------- */
    /* Short-hand field aliases                                              */
    /* --------------------------------------------------------------------- */
    double ***phi  = lag->temp;       /* temporary level-set field          */
    double ***data = vel->data;       /* staggered velocity component       */

    /* Geometric constants                                                   */
    const double h  = grid->dx_u[1];  /* uniform spacing                    */
    const double dV = h * h * h;

    /* --------------------------------------------------------------------- */
    /* Loop over every particle in the list                                  */
    /* --------------------------------------------------------------------- */
    Particle *p = p_list->start;
    while (p != NULL) {

        const double  R           = p->R;
        const double *X           = p->X;        /* centre */
        double       *Int_U       = p->Int_U;
        double       *Int_Omega   = p->Int_Omega;

        /* ----------------------------------------------------------------- */
        /* 1. Determine the (staggered) index bounds of the particle          */
        /* ----------------------------------------------------------------- */
        i_start = (int)floor((X[0]-R-xe[0])/h);
        j_start = (int)floor((X[1]-R-ye[0])/h);
        k_start = (int)floor((X[2]-R-ze[0])/h);

        i_end   = (int)ceil ((X[0]+R-xe[0])/h);
        j_end   = (int)ceil ((X[1]+R-ye[0])/h);
        k_end   = (int)ceil ((X[2]+R-ze[0])/h);

        /* Clamp to local grid & global domain                               */
        i_start = max(i_start, grid->G_Is);
        j_start = max(j_start, grid->G_Js);
        k_start = max(k_start, grid->G_Ks);

        i_end   = min(i_end  , grid->G_Ie);
        j_end   = min(j_end  , grid->G_Je);
        k_end   = min(k_end  , grid->G_Ke);

        i_end   = min(i_end  , grid->NX-1);
        j_end   = min(j_end  , grid->NY-1);
        k_end   = min(k_end  , grid->NZ-1);

        /* ----------------------------------------------------------------- */
        /* 2. Level-set evaluation at the eight nodes per control volume     */
        /* ----------------------------------------------------------------- */
        for (k = k_start; k <= k_end; k++)
            for (j = j_start; j <= j_end; j++)
                for (i = i_start; i <= i_end; i++)
                    phi[k][j][i] = lvl_set(xe[i], ye[j], ze[k]);

        /* ----------------------------------------------------------------- */
        /* 3. Integrate ρ̃ u and ρ̃ (r×u) over the particle volume            */
        /* ----------------------------------------------------------------- */
        for (k = k_start; k < k_end; k++) {
            r[2] = zc[k] - X[2];

            for (j = j_start; j < j_end; j++) {
                r[1] = yc[j] - X[1];

                for (i = i_start; i < i_end; i++) {
                    r[0] = xc[i] - X[0];

                    /* -- (a) Volume fraction of *this particle* in cell ---- */
                    vf_cell = 0.0;
                    sum_phi = 0.0;
                    for (ii = 0; ii < 2; ii++)
                        for (jj = 0; jj < 2; jj++)
                            for (kk = 0; kk < 2; kk++) {
                                temp = phi[k+kk][j+jj][i+ii];
                                if (temp < 0.0) vf_cell -= temp;
                                sum_phi += fabs(temp);
                            }
                    vf_cell = vf_cell / sum_phi;
                    vf[k][j][i] += vf_cell;

                    /* -- (b) Face-centred density ρ̃^{(c)} ----------------- */
#ifdef VOF_IBM

					// Get CELL-CENTERED density for scalar buoyancy integral
                    double rho_cc_here = rho_cc[k][j][i];
					double vf_center = vof->vfc[k][j][i];
                    
                    // Accumulate scalar density (same for all components)
                    // Only accumulate once per cell, not per velocity component
                    if (component == 0) {  // Only on u-pass to avoid triple counting
                        p->Int_rho_scalar += vf_center * rho_cc_here * dV;
                    }
					
                    double rho_loc;
                    if (component == 0) {     /* u-faces */
                        double rR = rho_cc[k][j][i];
                        double rL = rho_cc[k][j][i-1];
                        rho_loc   = 0.5 * (rR + rL);
                    } else if (component == 1) { /* v-faces */
                        double rT = rho_cc[k][j][i];
                        double rB = rho_cc[k][j-1][i];
                        rho_loc   = 0.5 * (rT + rB);
                    } else {                   /* w-faces */
                        double rF = rho_cc[k][j][i];
                        double rB = rho_cc[k-1][j][i];
                        rho_loc   = 0.5 * (rF + rB);
                    }
#else
                    const double rho_loc = 1.0; /* single-phase */
#endif

                    /* -- (c) Momentum & moment integrals ------------------- */
                    const double u_loc = data[k][j][i];
                    const double dVloc = vf_cell * rho_loc * u_loc * dV;
					
					#ifdef VOF_IBM
						// Integrate density for each component using the SAME face-centered 
						// stencil as used for momentum (ensures consistency: U = Int_U / Int_rho)
						p->Int_rho[component] += vf_cell * rho_loc * dV;
					#endif

                    Int_U[component]       += dVloc;

                    Int_Omega[comp_pos]   += dVloc * r[comp_neg];
                    Int_Omega[comp_neg]   -= dVloc * r[comp_pos];
                }
            }
        }

        p = p->next;
    }
}



/******************************************************************************/
/*
 * Interpolate_paint_rigid_body_velocity
 *
 * Blend rigid-body translation+rotation onto the staggered velocity face
 * (vel->data) inside each particle's level-set support, using the same
 * geometry stencil as Interpolate_integrate_momentum.  After painting and a
 * fresh integration, Int_U[i] == Int_rho[i] * U[i] for a quiescent exterior,
 * which removes the cold-start IBM impulse felt by an imposed-velocity
 * particle on the first time step (see Particle_initialize_velocities).
 *
 * Each rank paints only the cells it owns; callers are responsible for
 * synchronising halos and reapplying physical wall BCs afterwards.
 *
 * NOTE: assumes a uniform Cartesian grid of spacing h (LAG_PARTICLE_RESOLVED
 * already requires GRID_UNIFORM).
 */
/******************************************************************************/
void Interpolate_paint_rigid_body_velocity(Velocity      *vel,
                                           Particle_list *p_list,
                                           Cart3d_bag    *data_bag,
                                           Debug_trace   *dtrace)
{
#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
    Interpolate_paint_rigid_body_velocity_twod(vel, p_list, data_bag, dtrace);
    return;
#endif
    int i, j, k, ii, jj, kk;
    int i_start, j_start, k_start;
    int i_end,   j_end,   k_end;
    int component, comp_pos, comp_neg;

    char message[100];

    Parameters *params = data_bag->params;
    MAC_grid   *grid   = data_bag->grid;
    Lagrangian *lag    = data_bag->lag;

    double *xc = grid->xc, *yc = grid->yc, *zc = grid->zc;
    double *xe = grid->xu, *ye = grid->yv, *ze = grid->zw;

    if (vel->component == 'u') {
        xc        = grid->xu;
        xe        = &(grid->xc[-1]);
        component = 0;
        comp_pos  = 1;
        comp_neg  = 2;
    } else if (vel->component == 'v') {
        yc        = grid->yv;
        ye        = &(grid->yc[-1]);
        component = 1;
        comp_pos  = 2;
        comp_neg  = 0;
    } else if (vel->component == 'w') {
        zc        = grid->zw;
        ze        = &(grid->zc[-1]);
        component = 2;
        comp_pos  = 0;
        comp_neg  = 1;
    } else {
        sprintf(message, "Incorrect vel->component = '%c'", vel->component);
        Display_throw_error(message, params, DTRACE("Display_throw_error"));
        return;
    }

    Display_assert_list_state(p_list, LIST_STATE_BOTH,
                              params, DTRACE("Display_assert_list_state"));

    double ***phi  = lag->temp;
    double ***data = vel->data;
    const double h = grid->dx_u[1];

    Particle *p = p_list->start;
    while (p != NULL) {
        const double  R     = p->R;
        const double *X     = p->X;
        const double *Up    = p->U;
        const double *Omega = p->Omega;

        i_start = (int)floor((X[0]-R-xe[0])/h);
        j_start = (int)floor((X[1]-R-ye[0])/h);
        k_start = (int)floor((X[2]-R-ze[0])/h);

        i_end   = (int)ceil ((X[0]+R-xe[0])/h);
        j_end   = (int)ceil ((X[1]+R-ye[0])/h);
        k_end   = (int)ceil ((X[2]+R-ze[0])/h);

        i_start = max(i_start, grid->G_Is);
        j_start = max(j_start, grid->G_Js);
        k_start = max(k_start, grid->G_Ks);

        i_end   = min(i_end  , grid->G_Ie);
        j_end   = min(j_end  , grid->G_Je);
        k_end   = min(k_end  , grid->G_Ke);

        i_end   = min(i_end  , grid->NX-1);
        j_end   = min(j_end  , grid->NY-1);
        k_end   = min(k_end  , grid->NZ-1);

        for (k = k_start; k <= k_end; k++)
            for (j = j_start; j <= j_end; j++)
                for (i = i_start; i <= i_end; i++)
                    phi[k][j][i] = lvl_set(xe[i], ye[j], ze[k]);

        for (k = k_start; k < k_end; k++) {
            double r[3];
            r[2] = zc[k] - X[2];
            for (j = j_start; j < j_end; j++) {
                r[1] = yc[j] - X[1];
                for (i = i_start; i < i_end; i++) {
                    r[0] = xc[i] - X[0];

                    double vf_cell = 0.0;
                    double sum_phi = 0.0;
                    for (ii = 0; ii < 2; ii++)
                        for (jj = 0; jj < 2; jj++)
                            for (kk = 0; kk < 2; kk++) {
                                double tmp = phi[k+kk][j+jj][i+ii];
                                if (tmp < 0.0) vf_cell -= tmp;
                                sum_phi += fabs(tmp);
                            }
                    if (sum_phi <= 0.0) continue;
                    vf_cell = vf_cell / sum_phi;
                    if (vf_cell <= 0.0) continue;

                    /* u_rigid = U + (Omega x r) on this staggered face.
                     * comp_pos/comp_neg follow the same convention as the
                     * angular-momentum accumulation in the integrator. */
                    const double u_rigid = Up[component]
                        + Omega[comp_pos] * r[comp_neg]
                        - Omega[comp_neg] * r[comp_pos];

                    data[k][j][i] = vf_cell * u_rigid +
                                    (1.0 - vf_cell) * data[k][j][i];
                }
            }
        }

        p = p->next;
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

#if defined(LAG_PARTICLE_RESOLVED)
	/*
	 * Build IBM solid volume fractions from globally gathered particle geometry.
	 * This removes the processor-cut artifact for solids larger than one MPI
	 * neighbor exchange width.
	 */
	Interpolate_add_to_volume_fraction_global(component, p_list, data_bag, dtrace);
	return;
#endif

#if defined(TWOD_MODE) && defined(LAG_PARTICLE_RESOLVED)
	Interpolate_add_to_volume_fraction_twod(component, p_list, data_bag, dtrace);
	return;
#endif
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
