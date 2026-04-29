#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Array.h"
#include "Cart3d.h"
#include "Communication.h"
#include "Conc.h"
#include "Display.h"
#include "Grid.h"
#include "Immersed.h"
#include "Memory.h"
#include "MyMath.h"
#include "TwodOps.h"
#include "Velocity.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/******************************************************************************/
/*
 This function allocates memory for u-velocity structure
 */
/******************************************************************************/
Velocity *Velocity_create(MAC_grid *grid, Parameters *params, char which_velocity,
		Debug_trace *dtrace) {

	Velocity *new_vel;
	int NX, NY, NZ;
	int NX_cell, NY_cell, NZ_cell;
	int Ghost_Nodes;
	int x_s, y_s, z_s;
	int nx, ny, nz;

	// Number of real physical cells
	NX_cell = params->NXM;
	NY_cell = params->NYM;
	NZ_cell = params->NZM;

	// Add one more in each direction. This would help uniforming the for-loops
	// and data distribution layout
	NX = NX_cell + 1;
	NY = NY_cell + 1;
	NZ = NZ_cell + 1;

	new_vel = (Velocity *)malloc(sizeof(Velocity));
	Memory_check_allocation(new_vel);

	new_vel->component  = which_velocity; /* 'u', 'v' or 'w' */
	// Inlfux on current processor
	new_vel->G_influx   = 0.0;

	// Total influx (same on all processors)
	new_vel->W_influx   = 0.0;


	/*------------------------------------------------------------------------*/
	/*
	     - data    : value of the velocity
	     - data_bc : value of the velocity at the cell_center
	 */
	/*------------------------------------------------------------------------*/
	new_vel->data = Memory_allocate_flow_variable(grid, params);

	// Cell centered velocity is mainly used for concentration convective terms
	new_vel->data_bc = Memory_allocate_flow_variable(grid, params);
#ifdef VOF_SCALAR
	new_vel->data_vof = Memory_allocate_flow_variable(grid, params);
#endif


#ifdef LEFT_OUTFLOW
	new_vel-> d_dn_l_n = Memory_allocate_2D_double_array(NY, NZ);
	new_vel-> d_dn_l_o = Memory_allocate_2D_double_array(NY, NZ);
#endif
#ifdef RIGHT_OUTFLOW
	new_vel-> d_dn_r_n = Memory_allocate_2D_double_array(NY, NZ);
	new_vel-> d_dn_r_o = Memory_allocate_2D_double_array(NY, NZ);
#endif

#ifdef LEFT_INFLOW
	if (params -> vel_init_type == VEL_INIT_PRECURSOR) {
		new_vel-> inflow_n = Memory_allocate_2D_double_array(NY, NZ);
		new_vel-> inflow_o = Memory_allocate_2D_double_array(NY, NZ);
		new_vel-> inflow   = Memory_allocate_2D_double_array(NY, NZ);
	}
#endif

	/*------------------------------------------------------------------------*/
	/*
	 Allocate memory for variables that do not use ghost cells:
	     ng_explicit     = explicit terms
	     ng_explicit_old = explicit term at previous sub-step
	     ng_implicit     = implicit terms
	     ng_rhs          = right hand side
		 ng_visc_explicit        = viscous terms
		 ng_visc_explicit_old    = viscous terms at previous sub-step
	 */
	/*------------------------------------------------------------------------*/
	new_vel->ng_explicit = Memory_allocate_noghost_variable(grid, params);
	new_vel->ng_explicit_old = Memory_allocate_noghost_variable(grid, params);
	new_vel->ng_implicit = Memory_allocate_noghost_variable(grid, params);
	new_vel->ng_rhs = Memory_allocate_noghost_variable(grid, params);

	new_vel->ng_visc_explicit = Memory_allocate_noghost_variable(grid, params);
	new_vel->ng_visc_explicit_old = Memory_allocate_noghost_variable(grid, params);

	/*------------------------------------------------------------------------*/
	/*
	 Allocate memory for the shear stress on the bottom:
	     tau = 1/Re sqrt(dudy^2 + dwdy^2)
	 */
	/*------------------------------------------------------------------------*/
	if ( (which_velocity == 'u') && (params->shear_stress_output) ){

		new_vel->G_shear_stress_bottom  = Memory_allocate_2D_double_array(NX, NZ);
		new_vel->W_shear_stress_bottom  = Memory_allocate_2D_double_array(NX, NZ);

		new_vel->G_u_shear = Memory_allocate_2D_double_array(NX, NZ);
		new_vel->G_v_shear = Memory_allocate_2D_double_array(NX, NZ);
		new_vel->G_w_shear = Memory_allocate_2D_double_array(NX, NZ);

		new_vel->W_u_shear = Memory_allocate_2D_double_array(NX, NZ);
		new_vel->W_v_shear = Memory_allocate_2D_double_array(NX, NZ);
		new_vel->W_w_shear = Memory_allocate_2D_double_array(NX, NZ);

	}
	else {

		new_vel->G_shear_stress_bottom = NULL;
		new_vel->W_shear_stress_bottom = NULL;

		new_vel->G_u_shear = NULL;
		new_vel->G_v_shear = NULL;
		new_vel->G_w_shear = NULL;

		new_vel->W_u_shear = NULL;
		new_vel->W_v_shear = NULL;
		new_vel->W_w_shear = NULL;

	} // else

	if  (which_velocity == 'u') {
		new_vel->G_u_streak = Memory_allocate_2D_double_array(NX, NZ);
		new_vel->W_u_streak = Memory_allocate_2D_double_array(NX, NZ);
	}
	else {
		new_vel->G_u_streak = NULL;
		new_vel->W_u_streak = NULL;
	}

#ifdef TURB_FORCING
	new_vel -> fturb = Memory_allocate_flow_variable(grid, params);
#endif

#ifdef CG_SOLVE
	new_vel -> d = Memory_allocate_flow_variable(grid, params);
	new_vel -> ng_r = Memory_allocate_noghost_variable(grid, params);
	new_vel -> ng_Ad = Memory_allocate_noghost_variable(grid, params);
	#ifdef VOF
	new_vel -> M_inv = Memory_allocate_flow_variable(grid, params);
	#endif
#endif

#ifdef BICG_SOLVE
	new_vel -> p = Memory_allocate_flow_variable(grid, params);
	new_vel -> s = Memory_allocate_flow_variable(grid, params);
	new_vel -> ng_r = Memory_allocate_noghost_variable(grid, params);
	new_vel -> ng_r0 = Memory_allocate_noghost_variable(grid, params);
	new_vel -> Ap = Memory_allocate_flow_variable(grid, params);
	new_vel -> ng_As = Memory_allocate_noghost_variable(grid, params);
#endif

	Velocity_nonzero_initialize(new_vel, grid, params, DTRACE("Velocity_nonzero_initialize"));

	new_vel->solve_cpu_time = 0.0;
	new_vel->solve_first_copy_cpu_time = 0.0;
	new_vel->solve_second_copy_cpu_time = 0.0;
	new_vel->solve_first_remap_cpu_time = 0.0;
	new_vel->solve_second_remap_cpu_time = 0.0;
	new_vel->solve_tridiag_cpu_time = 0.0;

	new_vel->cell_center_comm_cpu_time = 0.0;
	new_vel->cell_center_u_cpu_time = 0.0;
	new_vel->cell_center_v_cpu_time = 0.0;
	new_vel->cell_center_w_cpu_time = 0.0;
	new_vel->rhs_p_cpu_time = 0.0;
	new_vel->rhs_loop_cpu_time = 0.0;

	return new_vel;
}




/******************************************************************************/
/*
 This function releases the allocated memory for velocity structure.
 */
/******************************************************************************/
void Velocity_destroy(Velocity *vel, MAC_grid *grid, Parameters *params) {

	int NZ;

	NZ = grid->NZ;

	// Velocity data
	Memory_free_flow_variable(grid, params, vel->data);
	Memory_free_flow_variable(grid, params, vel->data_bc);
#ifdef VOF_SCALAR
	Memory_free_flow_variable(grid, params, vel->data_vof);
#endif

#ifdef TURB_FORCING
	Memory_free_flow_variable(grid, params, vel->fturb);
#endif

	Memory_free_noghost_variable(grid, params, vel->ng_explicit);
	Memory_free_noghost_variable(grid, params, vel->ng_explicit_old);
	Memory_free_noghost_variable(grid, params, vel->ng_implicit);
	Memory_free_noghost_variable(grid, params, vel->ng_rhs);

	Memory_free_noghost_variable(grid, params, vel->ng_visc_explicit);
	Memory_free_noghost_variable(grid, params, vel->ng_visc_explicit_old);

	if (vel->G_u_shear != NULL) {
		Memory_free_2D_double_array(NZ, vel->G_u_shear);
		Memory_free_2D_double_array(NZ, vel->G_v_shear);
		Memory_free_2D_double_array(NZ, vel->G_w_shear);

		Memory_free_2D_double_array(NZ, vel->W_u_shear);
		Memory_free_2D_double_array(NZ, vel->W_v_shear);
		Memory_free_2D_double_array(NZ, vel->W_w_shear);
	}

	if (vel->G_shear_stress_bottom != NULL) {
		Memory_free_2D_double_array(NZ, vel->G_shear_stress_bottom);
		Memory_free_2D_double_array(NZ, vel->W_shear_stress_bottom);
	}

	if (vel->G_u_streak != NULL) {
		Memory_free_2D_double_array(NZ, vel->G_u_streak);
		Memory_free_2D_double_array(NZ, vel->W_u_streak);
	}

#ifdef CG_SOLVE
	Memory_free_flow_variable(grid, params, vel->d);
	Memory_free_noghost_variable(grid, params, vel->ng_r);
	Memory_free_noghost_variable(grid, params, vel->ng_Ad);
	#ifdef VOF
	Memory_free_flow_variable(grid, params, vel->M_inv);
	#endif
#endif

#ifdef BICG_SOLVE
	Memory_free_flow_variable(grid, params, vel->p);
	Memory_free_flow_variable(grid, params, vel->s);
	Memory_free_noghost_variable(grid, params, vel->ng_r);
	Memory_free_noghost_variable(grid, params, vel->ng_r0);
	Memory_free_flow_variable(grid, params, vel->Ap);
	Memory_free_noghost_variable(grid, params, vel->ng_As);
#endif

	free(vel);

}




/******************************************************************************/
/*
 This function calculates velocities at the cell center using linear average of
 values at the two-nodes.
 */
/******************************************************************************/
void Velocity_cell_center(Cart3d_bag *data_bag) {

	int i, j, k;
	double T1, T2;

	T1 = MPI_Wtime();

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;
	int pnodes = params -> ghost_nodes;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;

	// First, update the u, v, w velocity data ghost nodes
	Communication_update_ghost_nodes_flow_variable(u->data, U_VELOCITY, pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data, V_VELOCITY, pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data, W_VELOCITY, pnodes, data_bag);

	// Now, get the array using a local pointer to a 3D array
	double ***u_data = u -> data;
	double ***v_data = v -> data;
	double ***w_data = w -> data;

	// Now get the global array for Velocity at cell center
	double ***u_data_bc = u -> data_bc;
	double ***v_data_bc = v -> data_bc;
	double ***w_data_bc = w -> data_bc;

	// Start index of bottom-left-back corner on current processor
	int Is = grid -> L_Is;
	int Js = grid -> L_Js;
	int Ks = grid -> L_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid -> L_Ie - 1;
	int Je = grid -> L_Je - 1;
	int Ke = grid -> L_Ke - 1;



	// Find u at cell center
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				// Since xc[i] is halfway between xu[i] and xu[i+1]
				u_data_bc[k][j][i] = 0.5 * (u_data[k][j][i] + u_data[k][j][i+1]);
			} // for i
		} // for j
	} // for k

	// Now, find v at cell center
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				// Since yc[j] is halfway between yv[j] and yv[j+1]
				v_data_bc[k][j][i] = 0.5 * (v_data[k][j][i] + v_data[k][j+1][i]);

			} // for i
		} // for j
	} // for k

	// Now, find w at cell center
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				// Since zc[k] is halfway between zw[k] and zw[k+1]
				w_data_bc[k][j][i] = 0.5 * (w_data[k][j][i] + w_data[k+1][j][i]);
			} // for i
		} // for j
	} // for k

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_vel_cell_center += T2 - T1;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Velocity_u_set_implicit_explicit(Cart3d_bag *data_bag) {

	int i, j, k;

	double dudxE, dudxW, dudyN, dudyS, dudzF, dudzB;
	double dvdxN, dvdxS, dwdxF, dwdxB;
	double d2udx2, d2udy2, d2udz2, ddxdudx, ddydvdx,  ddzdwdx;
	double axisym_ur_correction;
	double vCN, vCS, wCF, wCB;
	double uuE, uuW, uvN, uvS, uwF, uwB;
	double duudx, duvdy, duwdz;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	const int collapsed_z = TwodOps_collapsed_component_is_inactive(params);

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Indices start and end on current processor
	int i_start = max(1, grid->G_Is); // i=0 not included
	int j_start = grid -> G_Js;
	int k_start = grid -> G_Ks;

	// Exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie);
	int j_end = min(NY-1, grid->G_Je);
	int k_end = min(NZ-1, grid->G_Ke);

#ifdef XPERIODIC
	if (i_end == NX-1)
		i_end = NX;
#endif

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;

	// Get the local velocities at the location where they are defined
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***u_data_bc = data_bag -> u -> data_bc;

	double ***explicit = data_bag -> u -> ng_explicit;
	double ***implicit = data_bag -> u -> ng_implicit;

	double ***visc_explicit = data_bag -> u -> ng_visc_explicit;

#ifdef VOF
    VolumeFraction *vof = data_bag->vof;
	double ***mu = vof->mu; // cell-centered current/stage viscosity
#endif

#ifdef VAR_VISC
	double ***nu  = data_bag -> viscosity -> nu;
	double ***nuY = data_bag -> viscosity -> nuY;
	double ***nuZ = data_bag -> viscosity -> nuZ;
#endif

#ifdef TOP_WALL_SCHUMANN
	double **dudy_wm_top    = data_bag -> log_law -> dudy_wm_top;
#endif
#ifdef BOTTOM_WALL_SCHUMANN
	double **dudy_wm_bottom = data_bag -> log_law -> dudy_wm_bottom;
#endif

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------
	double iRe = 1.0 / params -> Re;
	double nuE = iRe;
	double nuW = iRe;
	double nuN = iRe;
	double nuS = iRe;
	double nuF = iRe;
	double nuB = iRe;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++){

				/*------------------------------------------------------------*/
				/*
				 Calculate viscosity
				 */
				/*------------------------------------------------------------*/
#ifdef VAR_VISC
				//--------------------------------------------------------------
				// Variable viscosity
				//--------------------------------------------------------------
				nuE = nu[k][j][i];
				nuW = nu[k][j][i-1];

				nuN = nuZ[k][j+1][i];
				nuS = nuZ[k][j][i];

				nuF = nuY[k+1][j][i];
				nuB = nuY[k][j][i];
#endif // VAR_VISC

#ifdef VOF

               
                double muE = mu[k][j][i];
                double muW = mu[k][j][i-1];

                
                int im = i-1;
                int jp = j+1, jm = j-1;
                int kp = k+1, km = k-1;

                
                double mu_yp = 2.0 *
                    ( 0.5*(mu[k][j+1][i-1] + mu[k][j+1][i]) ) *
                    ( 0.5*(mu[k][j  ][i-1] + mu[k][j  ][i]) )
                    / ( (0.5*(mu[k][j+1][i-1] + mu[k][j+1][i])) +
                        (0.5*(mu[k][j  ][i-1] + mu[k][j  ][i])) + 1e-12 );

                double mu_ym = 2.0 *
                    ( 0.5*(mu[k][j  ][i-1] + mu[k][j  ][i]) ) *
                    ( 0.5*(mu[k][j-1][i-1] + mu[k][j-1][i]) )
                    / ( (0.5*(mu[k][j  ][i-1] + mu[k][j  ][i])) +
                        (0.5*(mu[k][j-1][i-1] + mu[k][j-1][i])) + 1e-12 );

                
                double mu_zp = 2.0 *
                    ( 0.5*(mu[k+1][j][i-1] + mu[k+1][j][i]) ) *
                    ( 0.5*(mu[k  ][j][i-1] + mu[k  ][j][i]) )
                    / ( (0.5*(mu[k+1][j][i-1] + mu[k+1][j][i])) +
                        (0.5*(mu[k  ][j][i-1] + mu[k  ][j][i])) + 1e-12 );

                double mu_zm = 2.0 *
                    ( 0.5*(mu[k  ][j][i-1] + mu[k  ][j][i]) ) *
                    ( 0.5*(mu[k-1][j][i-1] + mu[k-1][j][i]) )
                    / ( (0.5*(mu[k  ][j][i-1] + mu[k  ][j][i])) +
                        (0.5*(mu[k-1][j][i-1] + mu[k-1][j][i])) + 1e-12 );

#endif


				/*------------------------------------------------------------*/
				/*
				 Calculate viscous terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// du/dx East/West
				//--------------------------------------------------------------
				dudxE = ( u_data[k][j][i+1] - u_data[k][j][i] ) * idx_u[i];
				dudxW = ( u_data[k][j][i] - u_data[k][j][i-1] ) * idx_u[i-1];

				//--------------------------------------------------------------
				// du/dy North/South
				//--------------------------------------------------------------
				dudyN = ( u_data[k][j+1][i] - u_data[k][j][i] ) * idy_c[j];
#ifdef TOP_WALL_SCHUMANN
				if (j == NY - 2)
					dudyN = dudy_wm_top[k][i];
#endif

				dudyS = ( u_data[k][j][i] - u_data[k][j-1][i] ) * idy_c[j-1];
#ifdef BOTTOM_WALL_SCHUMANN
				if (j == 0)
					dudyS = log_law->dudy_wm_bottom[k][i];
#endif

				//--------------------------------------------------------------
				// du/dz Front/Back
				//--------------------------------------------------------------
				dudzF = 0.0;
				dudzB = 0.0;
				if (!collapsed_z) {
					dudzF = ( u_data[k+1][j][i] - u_data[k][j][i] ) * idz_c[k];
					dudzB = ( u_data[k][j][i] - u_data[k-1][j][i] ) * idz_c[k-1];
				}

#ifdef VOF

                //--------------------------------------------------------------
                // Axial / collapsed-direction diagonal pieces.
                // The x/r piece is rebuilt just below with the shared TwodOps
                // helper so AXISYM_RZ and the CG operator stay in lockstep.
                //--------------------------------------------------------------
                d2udy2 = iRe * ( mu_yp * dudyN - mu_ym * dudyS ) * idy_v[j];

                /*
                 * In 2D mode the stored z slab has no physical diffusion.
                 * Axisymmetric radial corrections are added in the dedicated
                 * mode-specific branch later in Phase 1.
                 */
                d2udz2 = 0.0;
                if (!collapsed_z)
                    d2udz2 = iRe * ( mu_zp * dudzF - mu_zm * dudzB ) * idz_w[k];

#else

				//--------------------------------------------------------------
				// y/z diagonal pieces. The x contribution is rebuilt below with
				// the shared TwodOps helper for the active 2D geometry.
				//--------------------------------------------------------------
				d2udy2 = ( nuN * dudyN - nuS * dudyS ) * idy_v[j];
				d2udz2 = 0.0;
				if (!collapsed_z)
					d2udz2 = ( nuF * dudzF - nuB * dudzB ) * idz_w[k];

#endif

				//--------------------------------------------------------------
				// d/dx(du/dx), d/dy(dv/dx), and d/dz(dw/dx)
				//--------------------------------------------------------------
				dvdxN = ( v_data[k][j+1][i] - v_data[k][j+1][i-1] ) * idx_c[i-1];
				dvdxS = ( v_data[k][j][i]   - v_data[k][j][i-1]   ) * idx_c[i-1];
				dwdxF = 0.0;
				dwdxB = 0.0;
				if (!collapsed_z) {
					dwdxF = ( w_data[k+1][j][i] - w_data[k+1][j][i-1] ) * idx_c[i-1];
					dwdxB = ( w_data[k][j][i]   - w_data[k][j][i-1]   ) * idx_c[i-1];
				}

#ifdef VOF
				/*
				 * The x/r contribution is the one place where AXISYM_RZ differs
				 * from planar Cartesian momentum.  We keep that metric weighting
				 * in TwodOps so the explicit kernel and the implicit CG operator
				 * share the exact same discrete geometry.
				 */
				d2udx2 = iRe * TwodOps_u_cv_x_flux_divergence(
					grid, params, i,
					muE * dudxE,
					muW * dudxW);
				ddxdudx = d2udx2;

				// Cross terms reuse the same face-μ as the diagonal terms:
				ddydvdx = iRe * ( mu_yp * dvdxN - mu_ym * dvdxS ) * idy_v[j];
				ddzdwdx = 0.0;
				if (!collapsed_z)
					ddzdwdx = iRe * ( mu_zp * dwdxF - mu_zm * dwdxB ) * idz_w[k];

				/*
				 * Extra cylindrical linear term in radial momentum:
				 *   -2 μ u_r / (Re r^2)
				 * There is no corresponding cross-term partner, so it belongs
				 * only to the implicit/diagonal side of the split.
				 */
				axisym_ur_correction = TwodOps_axisym_u_radial_linear_term(
					grid, params, i,
					iRe * 0.5 * (muE + muW),
					u_data[k][j][i]);

#else				
				d2udx2 = TwodOps_u_cv_x_flux_divergence(
					grid, params, i,
					nuE * dudxE,
					nuW * dudxW);
				ddxdudx = d2udx2;
				ddydvdx = ( nuN * dvdxN - nuS * dvdxS ) * idy_v[j];
				ddzdwdx = 0.0;
				if (!collapsed_z)
					ddzdwdx = ( nuF * dwdxF - nuB * dwdxB ) * idz_w[k];
				axisym_ur_correction = TwodOps_axisym_u_radial_linear_term(
					grid, params, i,
					iRe,
					u_data[k][j][i]);
#endif


				/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

#if defined(VOF) && !defined(VOF_DIFFUSE)
				double ***mfx = vof->mass_flux_x;
				double ***mfy = vof->mass_flux_y;
				double ***mfz = vof->mass_flux_z;

				/* ---- X-flux: d(ρ u u)/dx ---- */
				double mfx_E = 0.5*(mfx[k][j][i] + mfx[k][j][i+1]);
				double mfx_W = 0.5*(mfx[k][j][i-1] + mfx[k][j][i]);
				double u_E   = 0.5*(u_data[k][j][i] + u_data[k][j][i+1]);
				double u_W   = 0.5*(u_data[k][j][i-1] + u_data[k][j][i]);
				double drhouudx = (mfx_E*u_E - mfx_W*u_W) * idx_c[i-1];

				/* ---- Y-flux: d(ρ u v)/dy ---- */
				/* Mass flux at y-faces of the u-CV: interpolate mfy to x-face i */
				double mfy_N = 0.5*(mfy[k][j+1][i-1] + mfy[k][j+1][i]);
				double mfy_S = 0.5*(mfy[k][j  ][i-1] + mfy[k][j  ][i]);
				double u_N   = 0.5*(u_data[k][j][i]   + u_data[k][j+1][i]);
				double u_S   = 0.5*(u_data[k][j-1][i] + u_data[k][j  ][i]);
				double drhouvdy = (mfy_N*u_N - mfy_S*u_S) * idy_v[j];

				duwdz = 0.0;
				if (!collapsed_z) {
					/* No physical z-flux exists in the storage-only 2D slab. */
					double mfz_F = 0.5*(mfz[k+1][j][i-1] + mfz[k+1][j][i]);
					double mfz_B = 0.5*(mfz[k  ][j][i-1] + mfz[k  ][j][i]);
					double u_F   = 0.5*(u_data[k][j][i]   + u_data[k+1][j][i]);
					double u_B   = 0.5*(u_data[k-1][j][i] + u_data[k  ][j][i]);
					duwdz = (mfz_F*u_F - mfz_B*u_B) * idz_w[k];
				}

				duudx = TwodOps_u_cv_x_flux_divergence(
					grid, params, i,
					mfx_E*u_E,
					mfx_W*u_W);
				duvdy = (mfy_N*u_N - mfy_S*u_S) * idy_v[j];
#else
				/* Velocity-form convection for single-phase and diffuse-VOF runs.
				 * In the VOF_DIFFUSE case, rho-face weighting is applied in RHS. */
				vCN =  0.5 * ( v_data[k][j+1][i] + v_data[k][j+1][i-1] );
				vCS =  0.5 * ( v_data[k][j][i]   + v_data[k][j][i-1]   );
				wCF = 0.0;
				wCB = 0.0;
				if (!collapsed_z) {
					wCF =  0.5 * ( w_data[k+1][j][i] + w_data[k+1][j][i-1] );
					wCB =  0.5 * ( w_data[k][j][i]   + w_data[k][j][i-1]   );
				}

				//--------------------------------------------------------------
				// uu East/West
				//--------------------------------------------------------------
				uuE = u_data_bc[k][j][i]   * u_data_bc[k][j][i];
				uuW = u_data_bc[k][j][i-1] * u_data_bc[k][j][i-1];

				//--------------------------------------------------------------
				// uv North/South
				//--------------------------------------------------------------
				uvN = 0.5 * ( u_data[k][j][i] + u_data[k][j+1][i] ) * vCN;
				uvS = 0.5 * ( u_data[k][j][i] + u_data[k][j-1][i] ) * vCS;

				//--------------------------------------------------------------
				// uw Front/Back
				//--------------------------------------------------------------
				uwF = 0.0;
				uwB = 0.0;
				if (!collapsed_z) {
					uwF = 0.5 * ( u_data[k][j][i] + u_data[k+1][j][i] ) * wCF;
					uwB = 0.5 * ( u_data[k][j][i] + u_data[k-1][j][i] ) * wCB;
				}

				//--------------------------------------------------------------
				// d/dx(uu), d/dy(uv), and d/dz(uw)
				//--------------------------------------------------------------
				duudx = TwodOps_u_cv_x_flux_divergence(
					grid, params, i,
					uuE,
					uuW);
				duvdy = ( uvN - uvS ) * idy_v[j];
				duwdz = 0.0;
				if (!collapsed_z)
					duwdz = ( uwF - uwB ) * idz_w[k];
#endif


				/*------------------------------------------------------------*/
				/*
				 Store explicit and implicit terms
				 */
				/*------------------------------------------------------------*/

				explicit[k][j][i] = -(duudx + duvdy + duwdz);
				//explicit viscous (NOT mixed into 'explicit' for VOF)
				visc_explicit[k][j][i] = ddxdudx + ddydvdx + ddzdwdx; 

#ifdef FULLY_EXPLICIT

				explicit[k][j][i] += d2udx2 + d2udy2 + d2udz2 + axisym_ur_correction;
				implicit[k][j][i] = 0.0;

#elif defined VOF
				// implicit diagonal
				implicit[k][j][i] = d2udx2 + d2udy2 + d2udz2 + axisym_ur_correction;

#elif defined FULLY_IMPLICIT

				implicit[k][j][i] = d2udx2 + d2udy2 + d2udz2 + axisym_ur_correction;
#else

				explicit[k][j][i] += d2udx2 + d2udz2 + axisym_ur_correction;
				implicit[k][j][i] = d2udy2;

#endif

			} // for i
		} // for j
	} // for k

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Velocity_v_set_implicit_explicit(Cart3d_bag *data_bag) {

	int i, j, k;

	double dvdxE, dvdxW, dvdyN, dvdyS, dvdzF, dvdzB;
	double dudyE, dudyW, dwdyF, dwdyB;
	double d2vdx2, d2vdy2, d2vdz2, ddxdudy, ddydvdy,  ddzdwdy;
	double uCE, uCW, wCF, wCB;
	double vuE, vuW, vvN, vvS, vwF, vwB;
	double dvudx, dvvdy, dvwdz;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	const int collapsed_z = TwodOps_collapsed_component_is_inactive(params);

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = max(1, grid->G_Js); // j=0 not included
	int k_start = grid -> G_Ks;

	// Exclude the half cell added
	int i_end = min(NX-1, grid->G_Ie);
	int j_end = min(NY-1, grid->G_Je);
	#ifdef YPERIODIC
	if (j_end == NY-1){
		j_end = NY;
	}
	#endif
	int k_end = min(NZ-1, grid->G_Ke);

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;

	// Get the local velocities at the location where they are defined
	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***v_data_bc = data_bag -> v -> data_bc;

	double ***explicit = data_bag -> v -> ng_explicit;
	double ***implicit = data_bag -> v -> ng_implicit;

	double ***visc_explicit = data_bag -> v -> ng_visc_explicit;

#ifdef VOF
    // Retrieve the current/stage VOF viscosity field
    VolumeFraction *vof = data_bag->vof;
	double ***mu = vof->mu; // cell-centered current/stage viscosity
#endif

#ifdef VAR_VISC
	double ***nu  = data_bag -> viscosity -> nu;
	double ***nuX = data_bag -> viscosity -> nuX;
	double ***nuZ = data_bag -> viscosity -> nuZ;
#endif

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------
	double iRe = 1.0 / params -> Re;
	double nuE = iRe;
	double nuW = iRe;
	double nuN = iRe;
	double nuS = iRe;
	double nuF = iRe;
	double nuB = iRe;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {


				/*------------------------------------------------------------*/
				/*
				 Calculate viscosity
				 */
				/*------------------------------------------------------------*/
#ifdef VAR_VISC
				//--------------------------------------------------------------
				// Variable viscosity - eddy viscosity for LES and RANS
				//--------------------------------------------------------------
				nuE = nuZ[k][j][i+1];
				nuW = nuZ[k][j][i];

				nuN = nu[k][j][i];
				nuS = nu[k][j-1][i];

				nuF = nuX[k+1][j][i];
				nuB = nuX[k][j][i];
#endif // VAR_VISC

#ifdef VOF

                // --- Edge-centered μ for cross terms (keep) ---
                int im = i-1, ip = i+1;
                int jm = j-1;
                int km = k-1, kp = k+1;
				

                

                // x-direction (flux planes at i±1/2), project in y then harmonic in x
                double mu_xp = 2.0 *
                    ( 0.5*(mu[k][j-1][i+1] + mu[k][j][i+1]) ) *
                    ( 0.5*(mu[k][j-1][i  ] + mu[k][j][i  ]) )
                    / ( (0.5*(mu[k][j-1][i+1] + mu[k][j][i+1])) +
                        (0.5*(mu[k][j-1][i  ] + mu[k][j][i  ])) + 1e-12 );

                double mu_xm = 2.0 *
                    ( 0.5*(mu[k][j-1][i  ] + mu[k][j][i  ]) ) *
                    ( 0.5*(mu[k][j-1][i-1] + mu[k][j][i-1]) )
                    / ( (0.5*(mu[k][j-1][i  ] + mu[k][j][i  ])) +
                        (0.5*(mu[k][j-1][i-1] + mu[k][j][i-1])) + 1e-12 );

                // y-direction (normal y): direct harmonic across j
                double mu_yp = mu[k][j][i];
                double mu_ym = mu[k][j-1][i];

                // z-direction (flux planes at k±1/2), project in y then harmonic in z
                double mu_zp = 2.0 *
                    ( 0.5*(mu[k+1][j-1][i] + mu[k+1][j][i]) ) *
                    ( 0.5*(mu[k  ][j-1][i] + mu[k  ][j][i]) )
                    / ( (0.5*(mu[k+1][j-1][i] + mu[k+1][j][i])) +
                        (0.5*(mu[k  ][j-1][i] + mu[k  ][j][i])) + 1e-12 );

                double mu_zm = 2.0 *
                    ( 0.5*(mu[k  ][j-1][i] + mu[k  ][j][i]) ) *
                    ( 0.5*(mu[k-1][j-1][i] + mu[k-1][j][i]) )
                    / ( (0.5*(mu[k  ][j-1][i] + mu[k  ][j][i])) +
                        (0.5*(mu[k-1][j-1][i] + mu[k-1][j][i])) + 1e-12 );

#endif
				/*------------------------------------------------------------*/
				/*
				 Calculate viscous terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// dv/dx East/West
				//--------------------------------------------------------------
				dvdxE = ( v_data[k][j][i+1] - v_data[k][j][i] ) * idx_c[i];
				dvdxW = ( v_data[k][j][i] - v_data[k][j][i-1] ) * idx_c[i-1];

				//--------------------------------------------------------------
				// dv/dy North/South
				//--------------------------------------------------------------
				dvdyN = ( v_data[k][j+1][i] - v_data[k][j][i] ) * idy_v[j];
				dvdyS = ( v_data[k][j][i] - v_data[k][j-1][i] ) * idy_v[j-1];

				//--------------------------------------------------------------
				// dv/dz Front/Back
				//--------------------------------------------------------------
				dvdzF = 0.0;
				dvdzB = 0.0;
				if (!collapsed_z) {
					dvdzF = ( v_data[k+1][j][i] - v_data[k][j][i] ) * idz_c[k];
					dvdzB = ( v_data[k][j][i] - v_data[k-1][j][i] ) * idz_c[k-1];
				}

				//--------------------------------------------------------------
				// d2v/dx2, d2v/dy2, and d2v/dz2
				//--------------------------------------------------------------

#ifdef VOF

                d2vdx2 = iRe * TwodOps_v_cv_x_flux_divergence(
                    grid, params, i,
                    mu_xp * dvdxE,
                    mu_xm * dvdxW);
                d2vdy2 = iRe * ( mu_yp * dvdyN - mu_ym * dvdyS ) * idy_c[j-1];
                d2vdz2 = 0.0;
                if (!collapsed_z)
                    d2vdz2 = iRe * ( mu_zp * dvdzF - mu_zm * dvdzB ) * idz_w[k];

#else
				d2vdx2 = TwodOps_v_cv_x_flux_divergence(
					grid, params, i,
					nuE * dvdxE,
					nuW * dvdxW);
				d2vdy2 = ( nuN * dvdyN - nuS * dvdyS ) * idy_c[j-1];
				d2vdz2 = 0.0;
				if (!collapsed_z)
					d2vdz2 = ( nuF * dvdzF - nuB * dvdzB ) * idz_w[k];
#endif




				//--------------------------------------------------------------
				// d/dx(du/dy), d/dy(dv/dy), d/dz(dw/dy)
				//--------------------------------------------------------------
				dudyE = ( u_data[k][j][i+1] - u_data[k][j-1][i+1] ) * idy_c[j-1];
				dudyW = ( u_data[k][j][i]   - u_data[k][j-1][i]   ) * idy_c[j-1];
				dwdyF = 0.0;
				dwdyB = 0.0;
				if (!collapsed_z) {
					dwdyF = ( w_data[k+1][j][i] - w_data[k+1][j-1][i] ) * idy_c[j-1];
					dwdyB = ( w_data[k][j][i]   - w_data[k][j-1][i]   ) * idy_c[j-1];
				}

#ifdef VOF

				/*
				 * The axial equation uses the same radial metric weighting as
				 * the pressure and scalar operators, but it does not carry the
				 * extra -u_r/r^2 correction that is unique to radial momentum.
				 */
				ddxdudy = iRe * TwodOps_v_cv_x_flux_divergence(
					grid, params, i,
					mu_xp * dudyE,
					mu_xm * dudyW);
                ddydvdy = d2vdy2;
                ddzdwdy = 0.0;
                if (!collapsed_z)
                    ddzdwdy = iRe * ( mu_zp * dwdyF - mu_zm * dwdyB ) * idz_w[k];


#else				

				ddxdudy = TwodOps_v_cv_x_flux_divergence(
					grid, params, i,
					nuE * dudyE,
					nuW * dudyW);
				ddydvdy = d2vdy2;
				ddzdwdy = 0.0;
				if (!collapsed_z)
					ddzdwdy = ( nuF * dwdyF - nuB * dwdyB ) * idz_w[k];

#endif				


				/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

#if defined(VOF) && !defined(VOF_DIFFUSE)
				double ***mfx = vof->mass_flux_x;
				double ***mfy = vof->mass_flux_y;
				double ***mfz = vof->mass_flux_z;

				/* ---- X-flux: d(ρ v u)/dx ----
				* v-CV east/west faces are at x-faces i+1, i.
				* mfx lives at x-faces (cell-centred in y) → interpolate in y
				* across j-1/2 (the v-face): average cells j-1 and j.          */
				double mfx_E = 0.5*(mfx[k][j-1][i+1] + mfx[k][j][i+1]);
				double mfx_W = 0.5*(mfx[k][j-1][i  ] + mfx[k][j][i  ]);
				double v_E   = 0.5*(v_data[k][j][i  ] + v_data[k][j][i+1]);
				double v_W   = 0.5*(v_data[k][j][i-1] + v_data[k][j][i  ]);
				double dvudx = TwodOps_v_cv_x_flux_divergence(
					grid, params, i,
					mfx_E*v_E,
					mfx_W*v_W);

				/* ---- Y-flux: d(ρ v v)/dy ---- (diagonal — native direction)
				* v-CV north/south faces are at cell-centres j, j-1.
				* mfy lives at y-faces → must average pairs to reach cell-centres:
				*   north (cell-centre j)  : avg of y-faces j and j+1
				*   south (cell-centre j-1): avg of y-faces j-1 and j             */
				double mfy_N = 0.5*(mfy[k][j  ][i] + mfy[k][j+1][i]);
				double mfy_S = 0.5*(mfy[k][j-1][i] + mfy[k][j  ][i]);
				double v_N   = 0.5*(v_data[k][j  ][i] + v_data[k][j+1][i]);
				double v_S   = 0.5*(v_data[k][j-1][i] + v_data[k][j  ][i]);
				double dvvdy = (mfy_N*v_N - mfy_S*v_S) * idy_c[j-1];

				dvwdz = 0.0;
				if (!collapsed_z) {
					/*
					 * No physical z transport exists once the third direction is
					 * collapsed to a storage-only slab.
					 */
					double mfz_F = 0.5*(mfz[k+1][j-1][i] + mfz[k+1][j][i]);
					double mfz_B = 0.5*(mfz[k  ][j-1][i] + mfz[k  ][j][i]);
					double v_F   = 0.5*(v_data[k  ][j][i] + v_data[k+1][j][i]);
					double v_B   = 0.5*(v_data[k-1][j][i] + v_data[k  ][j][i]);
					dvwdz = (mfz_F*v_F - mfz_B*v_B) * idz_w[k];
				}
#else
				/* Velocity-form convection for single-phase and diffuse-VOF runs.
				 * In the VOF_DIFFUSE case, rho-face weighting is applied in RHS. */
				uCE = 0.5 * ( u_data[k][j][i+1] + u_data[k][j-1][i+1] );
				uCW = 0.5 * ( u_data[k][j][i]   + u_data[k][j-1][i]   );
				wCF = 0.0;
				wCB = 0.0;
				if (!collapsed_z) {
					wCF = 0.5 * ( w_data[k+1][j][i] + w_data[k+1][j-1][i] );
					wCB = 0.5 * ( w_data[k][j][i]   + w_data[k][j-1][i]   );
				}

				//--------------------------------------------------------------
				// vu East/West
				//--------------------------------------------------------------
				vuE = 0.5 * ( v_data[k][j][i] + v_data[k][j][i+1] ) * uCE;
				vuW = 0.5 * ( v_data[k][j][i] + v_data[k][j][i-1] ) * uCW;

				//--------------------------------------------------------------
				// vv North/South
				//--------------------------------------------------------------
				vvN = v_data_bc[k][j][i]   * v_data_bc[k][j][i];
				vvS = v_data_bc[k][j-1][i] * v_data_bc[k][j-1][i];

				//--------------------------------------------------------------
				// vw Front/Back
				//--------------------------------------------------------------
				vwF = 0.0;
				vwB = 0.0;
				if (!collapsed_z) {
					vwF = 0.5 * ( v_data[k][j][i] + v_data[k+1][j][i] ) * wCF;
					vwB = 0.5 * ( v_data[k][j][i] + v_data[k-1][j][i] ) * wCB;
				}

				//--------------------------------------------------------------
				// d/dx(vu), d/dy(vv), and d/dz(vw)
				//--------------------------------------------------------------
				dvudx = TwodOps_v_cv_x_flux_divergence(
					grid, params, i,
					vuE,
					vuW);
				dvvdy = (vvN - vvS) * idy_c[j-1];
				dvwdz = 0.0;
				if (!collapsed_z)
					dvwdz = (vwF - vwB) * idz_w[k];
#endif

				/*------------------------------------------------------------*/
				/*
				 Store explicit and implicit terms
				 */
				/*------------------------------------------------------------*/

				explicit[k][j][i] = -(dvudx + dvvdy + dvwdz);
				visc_explicit[k][j][i] = ddxdudy + ddydvdy + ddzdwdy;				

#ifdef FULLY_EXPLICIT

				explicit[k][j][i] += d2vdx2 + d2vdy2 + d2vdz2;
				implicit[k][j][i] = 0.0;

#elif defined VOF

				implicit[k][j][i] = d2vdx2 + d2vdy2 + d2vdz2; //implicit terms computed by matVec				

#elif defined FULLY_IMPLICIT
				implicit[k][j][i] = d2vdx2 + d2vdy2 + d2vdz2;

#else


				explicit[k][j][i] += d2vdx2 + d2vdz2;

				implicit[k][j][i] = d2vdy2;

#endif
			}
		}
	}

	return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Velocity_w_set_implicit_explicit(Cart3d_bag *data_bag) {

	int i, j, k;

	double dwdxE, dwdxW, dwdyN, dwdyS, dwdzF, dwdzB;
	double dudzE, dudzW, dvdzN, dvdzS;
	double d2wdx2, d2wdy2, d2wdz2, ddxdudz, ddydvdz,  ddzdwdz;
	double uCE, uCW, vCN, vCS;
	double wuE, wuW, wvN, wvS, wwF, wwB;
	double dwudx, dwvdy, dwwdz;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	if (TwodOps_collapsed_component_is_inactive(params)) {
		Memory_reset_noghost_variable(grid, params, data_bag->w->ng_explicit);
		Memory_reset_noghost_variable(grid, params, data_bag->w->ng_implicit);
		Memory_reset_noghost_variable(grid, params, data_bag->w->ng_visc_explicit);
		return;
	}

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Indices start and end on current processor
	int i_start = grid -> G_Is;
	int j_start = grid -> G_Js;
	int k_start = max(1, grid->G_Ks); // k=0 not included

	// Exclude the half cell added
	int i_end = min(NX-1, grid -> G_Ie);
	int j_end = min(NY-1, grid -> G_Je);
	int k_end = min(NZ-1, grid -> G_Ke);

#ifdef ZPERIODIC
	if (k_end == NZ-1)
		k_end = NZ;
#endif

	double *idx_u = grid -> idx_u;
	double *idy_v = grid -> idy_v;
	double *idz_w = grid -> idz_w;
	double *idx_c = grid -> idx_c;
	double *idy_c = grid -> idy_c;
	double *idz_c = grid -> idz_c;

	double ***u_data = data_bag -> u -> data;
	double ***v_data = data_bag -> v -> data;
	double ***w_data = data_bag -> w -> data;
	double ***w_data_bc = data_bag -> w -> data_bc;

	// Get convective terms global data
	double ***explicit = data_bag -> w -> ng_explicit;
	double ***implicit = data_bag -> w -> ng_implicit;

	double ***visc_explicit = data_bag -> w -> ng_visc_explicit;

#ifdef VOF
    // Retrieve the current/stage VOF viscosity field
    VolumeFraction *vof = data_bag->vof;
	double ***mu = vof->mu; // cell-centered current/stage viscosity
#endif

#ifdef VAR_VISC
	double ***nu  = data_bag -> viscosity -> nu;
	double ***nuX = data_bag -> viscosity -> nuX;
	double ***nuY = data_bag -> viscosity -> nuY;
#endif

#ifdef TOP_WALL_SCHUMANN
	double **dwdy_wm_top    = data_bag -> log_law -> dwdy_wm_top;
#endif
#ifdef BOTTOM_WALL_SCHUMANN
	double **dwdy_wm_bottom = data_bag -> log_law -> dwdy_wm_bottom;
#endif

	//--------------------------------------------------------------------------
	// Constant viscosity
	//--------------------------------------------------------------------------
	double iRe = 1.0 / params -> Re;
	double nuE = iRe;
	double nuW = iRe;
	double nuN = iRe;
	double nuS = iRe;
	double nuF = iRe;
	double nuB = iRe;

	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				/*------------------------------------------------------------*/
				/*
				 Calculate viscosity
				 */
				/*------------------------------------------------------------*/
#ifdef VAR_VISC
				//--------------------------------------------------------------
				// Variable viscosity - eddy viscosity for LES and RANS
				//--------------------------------------------------------------
				nuE = nuY[k][j][i+1];
				nuW = nuY[k][j][i];

				nuN = nuX[k][j+1][i];
				nuS = nuX[k][j][i];

				nuF = nu[k][j][i];
				nuB = nu[k-1][j][i];
#endif // VAR_VISC

#ifdef VOF

                
                int im = i-1, ip = i+1;
                int jm = j-1, jp = j+1;
                int km = k-1, kp = k+1;



                // x-direction (flux planes at i±1/2), project in z then harmonic in x
                double mu_xp = 2.0 *
                    ( 0.5*(mu[k-1][j][i+1] + mu[k][j][i+1]) ) *
                    ( 0.5*(mu[k-1][j][i  ] + mu[k][j][i  ]) )
                    / ( (0.5*(mu[k-1][j][i+1] + mu[k][j][i+1])) +
                        (0.5*(mu[k-1][j][i  ] + mu[k][j][i  ])) + 1e-12 );

                double mu_xm = 2.0 *
                    ( 0.5*(mu[k-1][j][i  ] + mu[k][j][i  ]) ) *
                    ( 0.5*(mu[k-1][j][i-1] + mu[k][j][i-1]) )
                    / ( (0.5*(mu[k-1][j][i  ] + mu[k][j][i  ])) +
                        (0.5*(mu[k-1][j][i-1] + mu[k][j][i-1])) + 1e-12 );

                // y-direction (flux planes at j±1/2), project in z then harmonic in y
                double mu_yp = 2.0 *
                    ( 0.5*(mu[k-1][j+1][i] + mu[k][j+1][i]) ) *
                    ( 0.5*(mu[k-1][j  ][i] + mu[k][j  ][i]) )
                    / ( (0.5*(mu[k-1][j+1][i] + mu[k][j+1][i])) +
                        (0.5*(mu[k-1][j  ][i] + mu[k][j  ][i])) + 1e-12 );

                double mu_ym = 2.0 *
                    ( 0.5*(mu[k-1][j  ][i] + mu[k][j  ][i]) ) *
                    ( 0.5*(mu[k-1][j-1][i] + mu[k][j-1][i]) )
                    / ( (0.5*(mu[k-1][j  ][i] + mu[k][j  ][i])) +
                        (0.5*(mu[k-1][j-1][i] + mu[k][j-1][i])) + 1e-12 );

                // z-direction (normal z): direct harmonic across k
                double mu_zp = mu[k][j][i];
                double mu_zm = mu[k-1][j][i];

#endif


				/*------------------------------------------------------------*/
				/*
				 Calculate viscous terms
				 */
				/*------------------------------------------------------------*/

				//--------------------------------------------------------------
				// dw/dx East/West
				//--------------------------------------------------------------
				dwdxE = ( w_data[k][j][i+1] - w_data[k][j][i] ) * idx_c[i];
				dwdxW = ( w_data[k][j][i] - w_data[k][j][i-1] ) * idx_c[i-1];

				//--------------------------------------------------------------
				// dw/dy North/South
				//--------------------------------------------------------------
				dwdyN = ( w_data[k][j+1][i] - w_data[k][j][i] ) * idy_c[j];
#ifdef TOP_WALL_SCHUMANN
				if (j == NY-2)
					dwdyN = log_law->dwdy_wm_top[k][i];
#endif

				dwdyS = ( w_data[k][j][i] - w_data[k][j-1][i] ) * idy_c[j-1];
#ifdef BOTTOM_WALL_SCHUMANN
				if (j == 0)
					dwdyS = log_law->dwdy_wm_bottom[k][i];
#endif

				//--------------------------------------------------------------
				// dw/dz Front/Back
				//--------------------------------------------------------------
				dwdzF = ( w_data[k+1][j][i] - w_data[k][j][i] ) * idz_w[k];
				dwdzB = ( w_data[k][j][i] - w_data[k-1][j][i] ) * idz_w[k-1];

				//--------------------------------------------------------------
				// d2w/dx2, d2w/dy2 and d2w/dz2
				//--------------------------------------------------------------
#ifdef VOF

                
                d2wdx2 = iRe * ( mu_xp * dwdxE - mu_xm * dwdxW ) * idx_u[i];
                d2wdy2 = iRe * ( mu_yp * dwdyN - mu_ym * dwdyS ) * idy_v[j];
                d2wdz2 = iRe * ( mu_zp * dwdzF - mu_zm * dwdzB ) * idz_c[k-1];

#else
				d2wdx2 = ( nuE * dwdxE - nuW * dwdxW ) * idx_u[i];
				d2wdy2 = ( nuN * dwdyN - nuS * dwdyS ) * idy_v[j];
				d2wdz2 = ( nuF * dwdzF - nuB * dwdzB ) * idz_c[k-1];
#endif

				//--------------------------------------------------------------
				// d/dx(du/dz), d/dy(dv/dz) and d/dz(dw/dz)
				//--------------------------------------------------------------
				dudzE = ( u_data[k][j][i+1] - u_data[k-1][j][i+1] ) * idz_c[k-1];
				dudzW = ( u_data[k][j][i]   - u_data[k-1][j][i]   ) * idz_c[k-1];
				dvdzN = ( v_data[k][j+1][i] - v_data[k-1][j+1][i] ) * idz_c[k-1];
				dvdzS = ( v_data[k][j][i]   - v_data[k-1][j][i]   ) * idz_c[k-1];

#ifdef VOF

				// Cross terms reuse diagonal face-μ:
				ddxdudz = iRe * ( mu_xp * dudzE - mu_xm * dudzW ) * idx_u[i];
				ddydvdz = iRe * ( mu_yp * dvdzN - mu_ym * dvdzS ) * idy_v[j];

#else

				ddxdudz = ( nuE * dudzE - nuW * dudzW ) * idx_u[i];
				ddydvdz = ( nuN * dvdzN - nuS * dvdzS ) * idy_v[j];
#endif
				ddzdwdz = d2wdz2;


				/*------------------------------------------------------------*/
				/*
				 Calculate convective terms
				 */
				/*------------------------------------------------------------*/

#if defined(VOF) && !defined(VOF_DIFFUSE)
				double ***mfx = vof->mass_flux_x;
				double ***mfy = vof->mass_flux_y;
				double ***mfz = vof->mass_flux_z;

				/* ---- X-flux: d(ρ w u)/dx ----
				* w-CV east/west faces are at x-faces i+1, i.
				* mfx lives at x-faces (cell-centred in z) → interpolate in z
				* across k-1/2 (the w-face): average cells k-1 and k.            */
				double mfx_E = 0.5*(mfx[k-1][j][i+1] + mfx[k][j][i+1]);
				double mfx_W = 0.5*(mfx[k-1][j][i  ] + mfx[k][j][i  ]);
				double w_E   = 0.5*(w_data[k][j][i  ] + w_data[k][j][i+1]);
				double w_W   = 0.5*(w_data[k][j][i-1] + w_data[k][j][i  ]);
				double dwudx = (mfx_E*w_E - mfx_W*w_W) * idx_u[i];

				/* ---- Y-flux: d(ρ w v)/dy ----
				* w-CV north/south faces are at y-faces j+1, j.
				* mfy lives at y-faces (cell-centred in z) → interpolate in z
				* across k-1/2: average cells k-1 and k.                          */
				double mfy_N = 0.5*(mfy[k-1][j+1][i] + mfy[k][j+1][i]);
				double mfy_S = 0.5*(mfy[k-1][j  ][i] + mfy[k][j  ][i]);
				double w_N   = 0.5*(w_data[k][j  ][i] + w_data[k][j+1][i]);
				double w_S   = 0.5*(w_data[k][j-1][i] + w_data[k][j  ][i]);
				double dwvdy = (mfy_N*w_N - mfy_S*w_S) * idy_v[j];

				/* ---- Z-flux: d(ρ w w)/dz ---- (diagonal — native direction)
				* w-CV front/back faces are at cell-centres k, k-1.
				* mfz lives at z-faces → must average pairs to reach cell-centres:
				*   front (cell-centre k)  : avg of z-faces k and k+1
				*   back  (cell-centre k-1): avg of z-faces k-1 and k             */
				double mfz_F = 0.5*(mfz[k  ][j][i] + mfz[k+1][j][i]);
				double mfz_B = 0.5*(mfz[k-1][j][i] + mfz[k  ][j][i]);
				double w_F   = 0.5*(w_data[k  ][j][i] + w_data[k+1][j][i]);
				double w_B   = 0.5*(w_data[k-1][j][i] + w_data[k  ][j][i]);
				double dwwdz = (mfz_F*w_F - mfz_B*w_B) * idz_c[k-1];

#else
				/* Velocity-form convection for single-phase and diffuse-VOF runs.
				 * In the VOF_DIFFUSE case, rho-face weighting is applied in RHS. */
				uCE = 0.5 * ( u_data[k][j][i+1] + u_data[k-1][j][i+1] );
				uCW = 0.5 * ( u_data[k][j][i]   + u_data[k-1][j][i]   );
				vCN = 0.5 * ( v_data[k][j+1][i] + v_data[k-1][j+1][i] );
				vCS = 0.5 * ( v_data[k][j][i]   + v_data[k-1][j][i]   );

				//--------------------------------------------------------------
				// wu East/West
				//--------------------------------------------------------------
				wuE = 0.5 * ( w_data[k][j][i] + w_data[k][j][i+1] ) * uCE;
				wuW = 0.5 * ( w_data[k][j][i] + w_data[k][j][i-1] ) * uCW;

				//--------------------------------------------------------------
				// wv North/South
				//--------------------------------------------------------------
				wvN = 0.5 * ( w_data[k][j][i] + w_data[k][j+1][i] ) * vCN;
				wvS = 0.5 * ( w_data[k][j][i] + w_data[k][j-1][i] ) * vCS;

				//--------------------------------------------------------------
				// ww Front/Back
				//--------------------------------------------------------------
				wwF = w_data_bc[k][j][i]   * w_data_bc[k][j][i];
				wwB = w_data_bc[k-1][j][i] * w_data_bc[k-1][j][i];

				//--------------------------------------------------------------
				// d/dx(wu), d/dy(wv), and d/dz(ww)
				//--------------------------------------------------------------
				dwudx = (wuE - wuW) * idx_u[i];
				dwvdy = (wvN - wvS) * idy_v[j];
				dwwdz = (wwF - wwB) * idz_c[k-1];
#endif


				/*------------------------------------------------------------*/
				/*
				 Store explicit and implicit terms
				 */
				/*------------------------------------------------------------*/

				explicit[k][j][i] = -(dwudx + dwvdy + dwwdz);
				visc_explicit[k][j][i] = ddxdudz + ddydvdz + ddzdwdz;				

#ifdef FULLY_EXPLICIT

				explicit[k][j][i] += d2wdx2 + d2wdy2 + d2wdz2;
				implicit[k][j][i] = 0.0;

#elif defined VOF

				implicit[k][j][i] = d2wdx2 + d2wdy2 + d2wdz2; //implicit terms computed by matVec				

#elif defined FULLY_IMPLICIT
				implicit[k][j][i] = d2wdx2 + d2wdy2 + d2wdz2;
#else

				explicit[k][j][i] += d2wdx2 + d2wdz2;
				implicit[k][j][i] = d2wdy2;

#endif

			} // for i
		} // for j
	} // for k

	return;
}




/******************************************************************************/
/*
 This function computes the RHS of the u-momentum linear system. Convective
 terms are treated explicitly
 */
/******************************************************************************/
void Velocity_u_set_RHS(Cart3d_bag *data_bag) {

	int i, j, k;
	double dpdx;
	double T1, T2;

//	int total_nodes, index;
//	double *rhs_1d, *conv_1d, *conv_old_1d, *visc_1d;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Same for all quantities
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	double *xu = grid -> xu;
	double *zc = grid -> zc;
	double *yc = grid -> yc;
	double Ly = params -> Ly;
	double Lz = params -> Lz;

	double xmax = params->xmax;

	// Start index of bottom-left-back corner on current processor
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = grid -> G_Ie;
	int Je = grid -> G_Je;
	int Ke = grid -> G_Ke;

	double time = params -> time;


	// Indices start and end on current processor
	int i_start = max(1,Is); // i=0 not included
	int j_start = Js;
	int k_start = Ks;

	// Exclude the half cell added
	int i_end = min(NX-1, Ie);
	int j_end = min(NY-1, Je);
	int k_end = min(NZ-1, Ke);

#ifdef XPERIODIC
	if (i_end == NX-1)
		i_end = NX;
#endif

	double *idx_c = grid -> idx_c;

	double dp_dx_source = -params -> dp_dx;


#if defined OSCILLATION

	double phase_shift_factor = params -> phase_shift_factor;
	double frequency = params -> frequency;
	double angular_vel = 2 * PI * frequency;	// = omega
	double amplitude;

	int amplitude_mode = params -> amplitude_mode;

	if (amplitude_mode == 1){
		double disp_amplitude = params -> disp_amplitude; 
		amplitude = disp_amplitude * pow(angular_vel,2);		// amplitude = disp_amplitude[L] * omega^2
	} 
	else {
		amplitude = params -> acc_amplitude;		// amplitude = acc_amplitude[L/T^2]
	}

	dp_dx_source = amplitude * sin(time * angular_vel + (phase_shift_factor * PI));

	params -> oscillation = dp_dx_source;

	if (params->oscillation_frame == 1){
		// non-inertial (accelerated) frame
		dp_dx_source = 0;
	} 
	// if inertial (fixed) frame
	// then dp_dx_source remains = oscillation

#endif // OSCILLATION


#ifdef PRESSURE_PULSE
	double ts, te;
	double Tmax  =  0.9;
	double pmax  = 10.0;
	double dt_dx =  0.1;
#endif
	// Velocity and pressure data
	double ***data   = data_bag -> u -> data;
	double ***p_data = data_bag -> p -> p_data;

	// 3-D arrays
	double ***rhs          = data_bag -> u -> ng_rhs;
	double ***explicit     = data_bag -> u -> ng_explicit;
	double ***explicit_old = data_bag -> u -> ng_explicit_old;
	double ***visc_expl    = data_bag -> u -> ng_visc_explicit;      
    double ***visc_expl_old= data_bag -> u -> ng_visc_explicit_old;   
	double ***implicit     = data_bag -> u -> ng_implicit;


#ifdef TURB_FORCING
	double ***fturb   = data_bag -> u -> fturb;
#endif

#ifdef VOF
    /* Stage/current density: rho has already been refreshed after the
     * VOF/CH substep before momentum RHS assembly. */
    double ***rho_data = data_bag->vof->rho;
#endif



	// Runge Kutta coefficients
	int rk = params -> which_stage;
	const double BET[] = {BETA};
	const double GAMB[] = {GAMBETA};
	const double ZETB[] = {ZETBETA};
	double a_dt = 1.0 / (params -> dt * BET[rk]);

	T1 = MPI_Wtime();


	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				dpdx = (p_data[k][j][i] - p_data[k][j][i-1]) * idx_c[i-1];



#ifdef VOF 
				double rho_face = 0.5 * (rho_data[k][j][i] + rho_data[k][j][i-1]);

	#ifdef STATIC_BUBBLE_TESTCASE //inertia terms ignored for static bubble test case
				rhs[k][j][i] = rho_face * a_dt * data[k][j][i] - 2.0 * dpdx
				             + implicit[k][j][i] + (GAMB[rk] * visc_expl[k][j][i] + ZETB[rk] * visc_expl_old[k][j][i]); 


	#else
				{
				double convective_term = GAMB[rk] * explicit[k][j][i]
				                       + ZETB[rk] * explicit_old[k][j][i];
				#ifdef VOF_DIFFUSE
				convective_term *= rho_face;
				#endif

				rhs[k][j][i] = rho_face * a_dt * data[k][j][i]
				             + convective_term
				             + implicit[k][j][i]
				             + (GAMB[rk] * visc_expl[k][j][i] + ZETB[rk] * visc_expl_old[k][j][i])
				             - 2.0 * dpdx;
				}
	#endif

#else

				rhs[k][j][i] = a_dt * data[k][j][i] - 2.0 * dpdx
				             + GAMB[rk] * explicit[k][j][i]
				             + ZETB[rk] * explicit_old[k][j][i]
				             + implicit[k][j][i];

#endif

#ifdef LAG_PARTICLE_RESOLVED
				rhs[k][j][i] += implicit[k][j][i];
#endif

#if defined XPERIODIC && !defined BOUSSINESQ 
	rhs[k][j][i] += 2.0 * dp_dx_source;

#elif defined OSCILLATION
	rhs[k][j][i] -= 2.0 * dp_dx_source;
	
#endif

#ifdef SWIMMERS_JET
				double r1 = 0.5;
				double ypos = yc[j] - Ly/2.0;
				double zpos = zc[k] - Lz/2.0;

				double rpos = sqrt(ypos*ypos + zpos*zpos);

				double delta = params->delta_jet;
				rhs[k][j][i] += (1.0/params->Re)/sqrt(PI)/pow(delta,3.0)*(rpos-r1)*exp(-pow((rpos-r1),2.0)/pow(delta,2.0))*(1.0 - erf((xu[i]-time)/delta));
#endif
#ifdef TURB_FORCING
				rhs[k][j][i] += 2.0  * fturb[k][j][i];
#endif
				explicit_old[k][j][i] = explicit[k][j][i];
				visc_expl_old[k][j][i] = visc_expl[k][j][i];      
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	data_bag->u->rhs_p_cpu_time += T2-T1;


double x_fr = params->x_fr;

double xmin = params->xmin;
double Pe = params->Pe[0];
double sig_profile = params->per_fr*(4*params->time/Pe + params->sig_profile);

if (params->time < params->t_release )
	{
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
									rhs[k][j][i] = rhs[k][j][i] * 0.5* (1.0+erf((xu[i] - x_fr)/sqrt(sig_profile)));
					}
				}
			}
		}

#ifdef IMMERSED_BOUNDARY
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {
				// if point is not fluid, set the rhs to zero
				if (grid->u_status[k][j][i] == SOLID)
					rhs[k][j][i] = 0.0;
			}
		}
	}
#endif

}




/******************************************************************************/
/*
 This function computes the RHS of the v-momentum linear system. Convective
 terms are treated explicitly
 */
/******************************************************************************/
void Velocity_v_set_RHS(Cart3d_bag *data_bag) {

	int i, j, k;
	double dpdy;
	double T1, T2;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Lagrangian *lag = data_bag -> lag;

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

	double *xc = grid -> xc;
	double *zc = grid -> zc;
	double *yv = grid -> yv;

	// Indices start and end on current processor
	int i_start = Is;
	int j_start = max(1, Js); // j=0 not included
	int k_start = Ks;

	// Exclude the half cell added
	int i_end = min(NX-1, Ie);
	int j_end = min(NY-1, Je);
	#ifdef YPERIODIC
	if (j_end==NY-1)
		j_end = NY;
	#endif
	int k_end = min(NZ-1, Ke);

	double *idy_c= grid->idy_c;

	// Velocity and pressure data
	double ***data   = data_bag -> v -> data;
	double ***p_data = data_bag -> p -> p_data;

	// Now, got the RHS vector on current processor
	double ***rhs          = data_bag -> v -> ng_rhs;
	double ***explicit     = data_bag -> v -> ng_explicit;
	double ***explicit_old = data_bag -> v -> ng_explicit_old;
	double ***visc_expl    = data_bag -> v -> ng_visc_explicit;       
    double ***visc_expl_old= data_bag -> v -> ng_visc_explicit_old;   
	double ***implicit     = data_bag -> v -> ng_implicit;
#ifdef TURB_FORCING
	double ***fturb     = data_bag -> v -> fturb;
#endif

#ifdef VOF
    /* Stage/current density: rho has already been refreshed after the
     * VOF/CH substep before momentum RHS assembly. */
    double ***rho_data = data_bag->vof->rho;
#endif

	// Runge Kutta coefficients
	int rk = params -> which_stage;
	const double BET[] = {BETA};
	const double GAMB[] = {GAMBETA};
	const double ZETB[] = {ZETBETA};
	double a_dt = 1.0 / (params -> dt * BET[rk]);

	#ifdef SWIMMERS_SWARM
	double x1 = params->xmin + (params->Lx)/2;
	double z1 = params->zmin + (params->Lz)/2;
	double y1 = params->ymin + params->y_swarm_start + params->V_swarm * params->time ;
	double delta_jet = params->delta_jet;
	double f0 = params->F_jet;
	double rswarm;
	double rpos;
	double eps = 0.1;
	double thet;
	double phi;
	double omega = 0.5;
	double myrand1 = (double)rand()/(double)RAND_MAX;
	double myrand2 = (double)rand()/(double)RAND_MAX;

	#endif

	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				dpdy = (p_data[k][j][i] - p_data[k][j-1][i]) * idy_c[j-1];

		  // -------------------------------------------------------------
          // Multiply the pressure gradient by -2 / rho if VOF is on;
          // else do the old -2.0 * dpdy if single-phase.
          // -------------------------------------------------------------

#ifdef VOF 
				double rho_face = 0.5 * (rho_data[k][j][i] + rho_data[k][j-1][i]);

	#ifdef STATIC_BUBBLE_TESTCASE //inertia terms ignored for static bubble test case
				rhs[k][j][i] = rho_face * a_dt * data[k][j][i] - 2.0 * dpdy
				             + implicit[k][j][i]
							 + (GAMB[rk] * visc_expl[k][j][i] + ZETB[rk] * visc_expl_old[k][j][i]); 


	#else
				{
				double convective_term = GAMB[rk] * explicit[k][j][i]
				                       + ZETB[rk] * explicit_old[k][j][i];
				#ifdef VOF_DIFFUSE
				convective_term *= rho_face;
				#endif

				rhs[k][j][i] = rho_face * a_dt * data[k][j][i]
				             + convective_term
				             + implicit[k][j][i]
				             + (GAMB[rk] * visc_expl[k][j][i] + ZETB[rk] * visc_expl_old[k][j][i])
				             - 2.0 * dpdy;
				}
	#endif

#else

				rhs[k][j][i] = a_dt * data[k][j][i] - 2.0 * dpdy
				             + GAMB[rk] * explicit[k][j][i]
				             + ZETB[rk] * explicit_old[k][j][i]
				             + implicit[k][j][i];

#endif

#ifdef LAG_PARTICLE_RESOLVED
				rhs[k][j][i] += implicit[k][j][i];
#endif
#ifdef TURB_FORCING
				rhs[k][j][i] += 2.0 * fturb[k][j][i];
#endif

#ifdef SWIMMERS_SWARM

				thet = atan( (yv[j] - y1) / (xc[i] - x1));
				phi = atan(sqrt((xc[i] - x1)*(xc[i] - x1) + (yv[j] - y1) *(yv[j] - y1) )/(zc[k] - z1));


				rpos = sqrt((xc[i] - x1)*(xc[i] - x1) + (zc[k] - z1)*(zc[k] - z1) + (yv[j] - y1)*(yv[j] - y1));

				//rswarm  = 1.0 + eps * params->swarm_rand*(-0.5+(double)rand()/(double)RAND_MAX)
				//rswarm  = 1.0 + eps * sin(2*PI*(phi+omega*params->time )) * sin(2*PI*(thet+omega*params->time ));
				rswarm  = 1.0 + eps * sin(2*PI*(phi+myrand1)) * sin(2*PI*(thet+myrand2));
				//rhs[k][j][i] -= (1.0 + params->swarm_rand*(-0.5+(double)rand()/(double)RAND_MAX)) * f0 * 0.5 * (1.0-erf((rpos-1.0)/delta_jet));
				rhs[k][j][i] -= f0 * 0.5 * (1.0-erf((rpos-rswarm)/delta_jet));
				//rhs[k][j][i] -= (1.0 + params->swarm_rand*(-0.5+(double)rand()/(double)RAND_MAX)) * f0 * exp(-rpos*rpos/(2*delta_jet*delta_jet));
#endif


				explicit_old[k][j][i] = explicit[k][j][i];
				visc_expl_old[k][j][i] = visc_expl[k][j][i];
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	data_bag->v->rhs_p_cpu_time += T2-T1;


double xmin = params->xmin;
double xmax = params->xmax;
double x_fr = params->x_fr;

double Pe = params->Pe[0];
double sig_profile = params->per_fr*(4*params->time/Pe + params->sig_profile);

if (params->time < params->t_release )
	{
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
									rhs[k][j][i] = rhs[k][j][i] * 0.5* (1.0+erf((xc[i] - x_fr)/sqrt(sig_profile)));
					}
				}
			}
		}



#ifdef IMMERSED_BOUNDARY
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				// If point is not fluid, set the rhs to zero
				if (grid->v_status[k][j][i] == SOLID)
					rhs[k][j][i] = 0.0;
			}
		}
	}
#endif

}


/******************************************************************************/
/*
 This function computes the RHS of the w-momentum linear system. Convective
 terms are treated explicitly
 */
/******************************************************************************/
void Velocity_w_set_RHS(Cart3d_bag *data_bag) {

	int i, j, k;
	double dpdz;
	double T1, T2;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	if (TwodOps_collapsed_component_is_inactive(params)) {
		Memory_reset_noghost_variable(grid, params, data_bag->w->ng_rhs);
		Memory_reset_noghost_variable(grid, params, data_bag->w->ng_explicit_old);
		Memory_reset_noghost_variable(grid, params, data_bag->w->ng_visc_explicit_old);
		return;
	}

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

	double *xc = grid -> xc;

	// Indices start and end on current processor
	int i_start = Is;
	int j_start = Js;
	int k_start = max(1,Ks); // k=0 not included

	// Exclude the half cell added
	int i_end   = min(NX-1, Ie);
	int j_end   = min(NY-1, Je);
	int k_end   = min(NZ-1, Ke);

#ifdef ZPERIODIC
	if (k_end == NZ-1)
		k_end = NZ;
#endif

	double *idz_c = grid -> idz_c;

	// Velocity and pressure data
	double ***data   = data_bag -> w -> data;
	double ***p_data = data_bag -> p -> p_data;

	double ***rhs          = data_bag -> w -> ng_rhs;
	double ***explicit     = data_bag -> w -> ng_explicit;
	double ***explicit_old = data_bag -> w -> ng_explicit_old;
	double ***visc_expl    = data_bag -> w -> ng_visc_explicit;       
    double ***visc_expl_old= data_bag -> w -> ng_visc_explicit_old;   
	double ***implicit     = data_bag -> w -> ng_implicit;
#ifdef TURB_FORCING
	double ***fturb     = data_bag -> w -> fturb;
#endif

#ifdef VOF
    /* Stage/current density: rho has already been refreshed after the
     * VOF/CH substep before momentum RHS assembly. */
    double ***rho_data = data_bag->vof->rho;  
#endif


	// Runge Kutta coefficients
	int rk = params -> which_stage;
	const double BET[] = {BETA};
	const double GAMB[] = {GAMBETA};
	const double ZETB[] = {ZETBETA};
	double a_dt = 1.0 / (params -> dt * BET[rk]);

	T1 = MPI_Wtime();
	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {

				dpdz = (p_data[k][j][i] - p_data[k-1][j][i]) * idz_c[k-1];

		  // -------------------------------------------------------------
          // Multiply the pressure gradient by -2 / rho if VOF is on;
          // else do the old -2.0 * dpdz if single-phase.
          // -------------------------------------------------------------
	  #ifdef VOF 
			  double rho_face = 0.5 * (rho_data[k][j][i] + rho_data[k-1][j][i]);

		  #ifdef STATIC_BUBBLE_TESTCASE //inertia terms ignored for static bubble test case
				  rhs[k][j][i] = rho_face * a_dt * data[k][j][i] - 2.0 * dpdz
							   + implicit[k][j][i]
							   + (GAMB[rk] * visc_expl[k][j][i] + ZETB[rk] * visc_expl_old[k][j][i]); 
	  
	  
		  #else
				  {
				  double convective_term = GAMB[rk] * explicit[k][j][i]
				                         + ZETB[rk] * explicit_old[k][j][i];
				  #ifdef VOF_DIFFUSE
				  convective_term *= rho_face;
				  #endif

				  rhs[k][j][i] = rho_face * a_dt * data[k][j][i]
				               + convective_term
				               + implicit[k][j][i]
				               + (GAMB[rk] * visc_expl[k][j][i] + ZETB[rk] * visc_expl_old[k][j][i])
				               - 2.0 * dpdz;
				  }
		   #endif

	  #else
	  
					  rhs[k][j][i] = a_dt * data[k][j][i] - 2.0 * dpdz
								   + GAMB[rk] * explicit[k][j][i]
								   + ZETB[rk] * explicit_old[k][j][i]
								   + implicit[k][j][i];
	  
	  #endif


#ifdef LAG_PARTICLE_RESOLVED
				rhs[k][j][i] += implicit[k][j][i];
#endif
#ifdef TURB_FORCING
				rhs[k][j][i] += 2.0 * fturb[k][j][i];
#endif

				explicit_old[k][j][i] = explicit[k][j][i];
				visc_expl_old[k][j][i] = visc_expl[k][j][i];
			} // for i
		} // for j
	} // for k
	T2 = MPI_Wtime();
	data_bag->w->rhs_p_cpu_time += T2-T1;

	double xmin = params->xmin;
	double xmax = params->xmax;
	double Pe = params->Pe[0];
	double sig_profile = params->per_fr*(4*params->time/Pe + params->sig_profile);
	double x_fr = params->x_fr;

	if (params->time < params->t_release )
		{
			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {
										rhs[k][j][i] = rhs[k][j][i] * 0.5* (1.0+erf((xc[i] - x_fr)/sqrt(sig_profile)));
						}
					}
				}
			}


#ifdef IMMERSED_BOUNDARY
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				// If point is not fluid, set the rhs to zero
				if (grid->w_status[k][j][i] == SOLID)
					rhs[k][j][i] = 0.0;
			}
		}
	}
#endif

}

#ifdef BOUSSINESQ
void Velocity_add_buoyancy_2_RHS(Cart3d_bag *data_bag) {

	double T1, T2;
	int i, j, k;
	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Lagrangian *lag = data_bag -> lag;

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
	int i_start = Is;
	int j_start = Js; // j=0 not included
	int k_start = Ks;

	int i_start_u = max(1,Is);
	int j_start_v = max(1,Js); // j=0 not included
	int k_start_w = max(1,Ks);


	// Exclude the half cell added
	int i_end = min(NX-1, Ie);
	int j_end = min(NY-1, Je);
	int k_end = min(NZ-1, Ke);

	int iconc;

#ifndef GRID_UNIFORM
	double *wc2uW = grid->wc2uW;
	double *wc2uE = grid->wc2uE;
	double *wc2vN = grid->wc2vN;
	double *wc2vS = grid->wc2vS;
	double *wc2wF = grid->wc2wF;
	double *wc2wB = grid->wc2wB;
#endif

	// Now, got the RHS vector on current processor

	double ***rhs_u    = data_bag -> u -> ng_rhs;
	double ***rhs_v    = data_bag -> v -> ng_rhs;
	double ***rhs_w    = data_bag -> w -> ng_rhs;


#ifdef   LAG_PARTICLE_RESOLVED
	double ***ng_vfc= lag->ng_vfc;
#ifdef VOF_SCALAR
	double ***ng_vfu= lag->vfu;
	double ***ng_vfv= lag->vfv;
	double ***ng_vfw= lag->vfw;
#else
	double ***ng_vfu= lag->ng_vfu;
	double ***ng_vfv= lag->ng_vfv;
	double ***ng_vfw= lag->ng_vfw;
#endif
#endif

	double c_term_u;
	double c_term_v;
	double c_term_w;

	double ***conc_total;
	int NConc = params -> NConc;
	Concentration **c = data_bag -> c;

	double richardson_0=1;
	double *richardson= params->richardson;
	double *grav= params->grav;
	double grav_magnitude= sqrt(grav[0]*grav[0]+grav[1]*grav[1]+grav[2]*grav[2]);
	double unit_gravity_x, unit_gravity_y, unit_gravity_z;

	if (grav_magnitude > 0) {
		unit_gravity_x= grav[0]/grav_magnitude;
		unit_gravity_y= grav[1]/grav_magnitude;
		unit_gravity_z= grav[2]/grav_magnitude;
	}
	else{
		unit_gravity_x= 0;
		unit_gravity_y= 0;
		unit_gravity_z= 0;
	}

	// Only one concentration field. Use c[0]->conc in v-momentum equation
	if (NConc < 2) {


		conc_total = c[0]->data;
		richardson_0 = richardson[0];
	}
	// If we have more than one concentration field, c_term in y-momentum is
	// stored in c_total (always in c[0])
	else {
		conc_total = c[0]->c_total;

		T1 = MPI_Wtime();
		for (k = k_start; k < k_end; k++) {
			for (j = j_start; j < j_end; j++) {
				for (i = i_start; i < i_end; i++) {
					conc_total[k][j][i]=0;
		}}}

		for (iconc=0; iconc<NConc; iconc++) {



			for (k = k_start; k < k_end; k++) {
				for (j = j_start; j < j_end; j++) {
					for (i = i_start; i < i_end; i++) {
						conc_total[k][j][i] += c[iconc]->data[k][j][i]*richardson[iconc];
			}}}
		}




	}
	Communication_update_ghost_nodes_flow_variable(conc_total, 'c', 1, data_bag);



	T1 = MPI_Wtime();



	for (k = k_start; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start_u; i < i_end; i++) {


				// Just linear interpolation for face centered nodes

#ifdef GRID_UNIFORM
	#ifdef LAG_PARTICLE_RESOLVED
				c_term_u =  (conc_total[k][j][i]+ conc_total[k][j][i-1])*( 1-ng_vfu[k][j][i])*0.5;


	#else
				c_term_u =  (conc_total[k][j][i] + conc_total[k][j][i-1])*0.5;

	#endif
#else

	#ifdef LAG_PARTICLE_RESOLVED
				c_term_u =  ( wc2uE[i]*conc_total[k][j][i] + wc2uW[i]*conc_total[k][j][i-1])*( 1-ng_vfu[k][j][i])  ;



	#else
				c_term_u =  ( wc2uE[i]*conc_total[k][j][i] + wc2uW[i]*conc_total[k][j][i-1] );

	#endif
#endif





				rhs_u[k][j][i] += 2.0 * c_term_u * richardson_0*unit_gravity_x;  //


			} // for i
		} // for j
	} // for k






	for (k = k_start; k < k_end; k++) {
		for (j = j_start_v; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {


				// Just linear interpolation for face centered nodes

#ifdef GRID_UNIFORM
	#ifdef LAG_PARTICLE_RESOLVED
				c_term_v =  (conc_total[k][j][i]+ conc_total[k][j-1][i] )*0.5*( 1-ng_vfv[k][j][i]);
			//	c_term_v =  (conc_total[k][j][i]*( 1-ng_vfc[k][j][i]) + conc_total[k][j-1][i]*( 1-ng_vfc[k][j-1][i]) )*0.5;

	#else
				c_term_v =  (conc_total[k][j][i] + conc_total[k][j-1][i])*0.5;
	#endif
#else

	#ifdef LAG_PARTICLE_RESOLVED
				c_term_v =  ( wc2vN[j]*conc_total[k][j][i]+ wc2vS[j]*conc_total[k][j-1][i])*0.5*( 1-ng_vfv[k][j][i]) ;


	#else
				c_term_v =  ( wc2vN[j]*conc_total[k][j][i] + wc2vS[j]*conc_total[k][j-1][i] );

	#endif
#endif



				rhs_v[k][j][i] += 2.0 * c_term_v * richardson_0*unit_gravity_y;  //

			} // for i
		} // for j
	} // for k




	for (k = k_start_w; k < k_end; k++) {
		for (j = j_start; j < j_end; j++) {
			for (i = i_start; i < i_end; i++) {


				// Just linear interpolation for face centered nodes

#ifdef GRID_UNIFORM
	#ifdef LAG_PARTICLE_RESOLVED
				c_term_w =  (conc_total[k][j][i] + conc_total[k-1][j][i])*0.5*( 1-ng_vfw[k][j][i]) ;

	#else
				c_term_w =  (conc_total[k][j][i] + conc_total[k-1][j][i])*0.5;
	#endif
#else

	#ifdef LAG_PARTICLE_RESOLVED
				c_term_w =  ( wc2wF[k]*conc_total[k][j][i]+ wc2wB[k]*conc_total[k-1][j][i] )*( 1-ng_vfw[k][j][i]) ;


	#else
				c_term_w =  ( wc2wF[k]*conc_total[k][j][i] + wc2wB[k]*conc_total[k-1][j][i] );

	#endif
#endif


				rhs_w[k][j][i] += 2.0 * c_term_w * richardson_0*unit_gravity_z;  //

			} // for i
		} // for j
	} // for k





	T2 = MPI_Wtime();

	data_bag->u->rhs_p_cpu_time += (T2-T1)/3;
	data_bag->v->rhs_p_cpu_time += (T2-T1)/3;
	data_bag->w->rhs_p_cpu_time += (T2-T1)/3;




#ifdef IMMERSED_BOUNDARY
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				// If point is not fluid, set the rhs to zero
				if (grid->v_status[k][j][i] == SOLID)
					rhs_v[k][j][i] = 0.0;
				if (grid->u_status[k][j][i] == SOLID)
					rhs_u[k][j][i] = 0.0;
				if (grid->w_status[k][j][i] == SOLID)
					rhs_w[k][j][i] = 0.0;

			}
		}
	}
#endif

}
#endif




/******************************************************************************/
/*
 Updates boundary values for 'data'

 'component' indicates which velocity component: 'u', 'v', or 'w'
 'type' is currently set to be one of the following:
     - VEL_TYPE_NORMAL: the actual velocity field 'u'
     - VEL_TYPE_CG: the velocity correction field 'd' used by CG and BiCG solvers

 Note: if the velocities are modified at their face boundary (e.g. 'u' at x=0),
 then the ghost cells in the other directions may not reflect the correct value
 (e.g. ghost cells in y-direction at x=0)

 Note: This function is used by Velocity_solve_cg() and Velocity_solve_bicg().
 Be mindful of the boundary conditions for the variables in these methods when
 updating this function.
 */
/******************************************************************************/
void Velocity_update_boundaries(double ***data, char component, int type, Cart3d_bag *data_bag) {

	int i, j, k;
	double T1, T2;

	T1 = MPI_Wtime();

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

#ifdef SHEARED_PERIODIC
	if (component == 'u') {
		if (type == VEL_TYPE_CG) {
			component = U_VELOCITY_PERTURBATION;
		}
		else {
			component = U_VELOCITY;

		}
	}
	else if (component == 'v'){
		component = V_VELOCITY;
	}
	else {
		component = W_VELOCITY;
	}
#endif

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

	// Start index of bottom-left-back corner on current processor
	int Is_g = grid -> L_Is;
	int Js_g = grid -> L_Js;
	int Ks_g = grid -> L_Ks;

	// End index of top-right-front corner on current processor
	int Ie_g = grid -> L_Ie;
	int Je_g = grid -> L_Je;
	int Ke_g = grid -> L_Ke;

	// Top wall velocity, if enabled
#if defined BOTTOM_WALL_VELOCITY && defined TOP_WALL_VELOCITY
	double top_wall_vel = params -> ubulk_target;
#else
	#ifdef STOKES_2ND_PROBLEM
		double u_oscillation;
		double time = params->time;
		double phase_shift_factor = params->phase_shift_factor;
		double frequency = params->frequency;
		double omega = 2 * PI * frequency; 		// angular velocity
		double amplitude = params->disp_amplitude;

		u_oscillation = -amplitude * omega * cos(time * omega + (phase_shift_factor * PI));

		params->u_oscillation = u_oscillation;

		double top_wall_vel = u_oscillation;

	#else
		double top_wall_vel = 2.0 * params -> ubulk_target;
	#endif
#endif

	/*------------------------------------------------------------------------*/
	/*
	 Update non-periodic boundaries
	 */
	/*------------------------------------------------------------------------*/

#ifndef XPERIODIC
	//--------------------------------------------------------------------------
	// Left boundary
	//--------------------------------------------------------------------------
	if (Is == 0) {
		i = 0;
	#ifdef LEFT_INFLOW
		//----------------------------------------------------------------------
		// Left inflow
		//----------------------------------------------------------------------
		// nothing to do

	#elif defined LEFT_OUTFLOW
		//----------------------------------------------------------------------
		// Right outflow
		//----------------------------------------------------------------------
		// nothing to do
	#else
		//----------------------------------------------------------------------
		// Left wall
		//----------------------------------------------------------------------
		if (component == 'u') {
			/*
			 * In AXISYM_RZ the left boundary is the symmetry axis r = 0, so the
			 * radial component must vanish there for both the physical velocity
			 * field and the CG correction field used inside the implicit solve.
			 */
			if (type == VEL_TYPE_NORMAL || params->axisym_rz_enabled) {
				for (k = Ks_g; k < Ke_g; k++) {
					for (j = Js_g; j < Je_g; j++) {
						data[k][j][i] = 0.0;
					}
				}
			}
			}
			else {
				for (k = Ks; k < Ke; k++) {
					for (j = Js; j < Je; j++) {
#ifdef AXISYM_RZ
						data[k][j][i-1] = data[k][j][i];
#endif
			#ifdef LEFT_WALL_VELOCITY_NOSLIP
						data[k][j][i-1] = -data[k][j][i];
			#endif
		#ifdef LEFT_WALL_VELOCITY_FREESLIP
					data[k][j][i-1] = data[k][j][i];
		#endif
				}
			}
		}
	#endif
	}

	//--------------------------------------------------------------------------
	// Right boundary
	//--------------------------------------------------------------------------
	if (Ie == NX) {
		i = NX-1;
	#ifdef RIGHT_OUTFLOW
		//----------------------------------------------------------------------
		// Right outflow
		//----------------------------------------------------------------------
		// nothing to do
	#elif defined RIGHT_INFLOW
		//----------------------------------------------------------------------
		// Right inflow
		//----------------------------------------------------------------------
		// nothing to do
	#else
		//----------------------------------------------------------------------
		// Right wall
		//----------------------------------------------------------------------
		if (component == 'u') {
			// Don't waste time resetting for CG method
			if (type == VEL_TYPE_NORMAL) {
				for (k = Ks_g; k < Ke_g; k++) {
					for (j = Js_g; j < Je_g; j++) {
						data[k][j][i] = 0.0;
					}
				}
			}
		}
		else {
			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
		#ifdef RIGHT_WALL_VELOCITY_NOSLIP
					data[k][j][i] = -data[k][j][i-1];
		#endif
		#ifdef RIGHT_WALL_VELOCITY_FREESLIP
					data[k][j][i] = data[k][j][i-1];
		#endif
				}
			}
		}
	#endif
	}
#endif // not XPERIODIC

#ifndef YPERIODIC
	//--------------------------------------------------------------------------
	// Bottom wall
	//--------------------------------------------------------------------------
	if (Js == 0) {
		j = 0;
		if (component == 'v') {
			// Don't waste time resetting for CG method
			if (type == VEL_TYPE_NORMAL) {
				for (k = Ks_g; k < Ke_g; k++) {
					for (i = Is_g; i < Ie_g; i++) {
						data[k][j][i] = 0.0;
					}
				}
			}
		}
		else {
	#ifdef BOTTOM_WALL_VELOCITY
			if (component == 'u' && type == VEL_TYPE_NORMAL) {
				for (k = Ks; k < Ke; k++) {
					for (i = Is; i < Ie; i++) {
						data[k][j-1][i] = -2.0 * top_wall_vel - data[k][j][i];
					}
				}
			}  // u-component
			// Use this for u-component of CG/BICG methods
			else if (component == 'w' || (component == 'u' && type == VEL_TYPE_CG)) {
				for (k = Ks; k < Ke; k++) {
					for (i = Is; i < Ie; i++) {
						data[k][j-1][i] = -data[k][j][i];
					}
				}
			}  // w-component

	#else  // not lid velocity
			for (k = Ks; k < Ke; k++) {
				for (i = Is; i < Ie; i++) {
		#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
					data[k][j-1][i] = -data[k][j][i];
		#endif
		#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
					data[k][j-1][i] = data[k][j][i];
		#endif
				}
			}
	#endif
		}
	}

	//--------------------------------------------------------------------------
	// Top wall
	//--------------------------------------------------------------------------
	if (Je == NY) {
		j = NY-1;
		if (component == 'v') {

			#ifdef TOP_WALL_OUTFLOW
			if (type == VEL_TYPE_CG) {
				for (k = Ks_g; k < Ke_g; k++) {
					for (i = Is_g; i < Ie_g; i++) {
						data[k][j][i] = 0.0;
					}
				}
			}
			if (type == VEL_TYPE_NORMAL) {
				for (k = Ks_g; k < Ke_g; k++) {
					for (i = Is_g; i < Ie_g; i++) {
						data[k][j][i] = params->voutflow;
					}
				}
			}
			#else
			// Don't waste time resetting for CG method
			if (type == VEL_TYPE_NORMAL) {
				for (k = Ks_g; k < Ke_g; k++) {
					for (i = Is_g; i < Ie_g; i++) {
						data[k][j][i] = 0.0;
					}
				}
			}
			#endif
		}
		else {

	#ifdef TOP_WALL_VELOCITY
			if (component == 'u' && type == VEL_TYPE_NORMAL) {
				for (k = Ks; k < Ke; k++) {
					for (i = Is; i < Ie; i++) {
						data[k][j][i] = 2.0 * top_wall_vel - data[k][j-1][i];
					}
				}
			}  // u-component
			// Use this for u-component of CG/BICG methods
			else if (component == 'w' || (component == 'u' && type == VEL_TYPE_CG)) {
				for (k = Ks; k < Ke; k++) {
					for (i = Is; i < Ie; i++) {
						data[k][j][i] = -data[k][j-1][i];
					}
				}
			}  // w-component

	#else  // not lid velocity
			for (k = Ks; k < Ke; k++) {
				for (i = Is; i < Ie; i++) {
		#if defined TOP_WALL_VELOCITY_NOSLIP || defined TOP_WALL_SCHUMANN
					data[k][j][i] = -data[k][j-1][i];
		#endif
		#ifdef TOP_WALL_VELOCITY_FREESLIP
					data[k][j][i] = data[k][j-1][i];
		#endif
				}
			}
	#endif
		} // not v-component
	}  // Je == NY
#endif // not YPERIODIC

#ifndef ZPERIODIC
	//--------------------------------------------------------------------------
	// Back wall
	//--------------------------------------------------------------------------
	if (Ks == 0) {
		k = 0;
		if (component == 'w') {
			// Don't waste time resetting for CG method
			if (type == VEL_TYPE_NORMAL) {
				for (j = Js_g; j < Je_g; j++) {
					for (i = Is_g; i < Ie_g; i++) {
						data[k][j][i] = 0.0;
					}
				}
			}
		}
		else {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
	#ifdef BACK_WALL_VELOCITY_NOSLIP
					data[k-1][j][i] = -data[k][j][i];
	#endif
	#ifdef BACK_WALL_VELOCITY_FREESLIP
					data[k-1][j][i] = data[k][j][i];
	#endif
				}
			}
		}
	}

	//--------------------------------------------------------------------------
	// Front wall
	//--------------------------------------------------------------------------
	if (Ke == NZ) {
		k = NZ-1;
		if (component == 'w') {
			// Don't waste time resetting for CG method
			if (type == VEL_TYPE_NORMAL) {
				for (j = Js_g; j < Je_g; j++) {
					for (i = Is_g; i < Ie_g; i++) {
						data[k][j][i] = 0.0;
					}
				}
			}
		}
		else {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
	#ifdef FRONT_WALL_VELOCITY_NOSLIP
					data[k][j][i] = -data[k-1][j][i];
	#endif
	#ifdef FRONT_WALL_VELOCITY_FREESLIP
					data[k][j][i] = data[k-1][j][i];
	#endif
				}
			}
		}
	}
#endif // not ZPERIODIC


	/*------------------------------------------------------------------------*/
	/*
	 Update ghost nodes.  This takes care of periodic boundaries
	 */
	/*------------------------------------------------------------------------*/
	Communication_update_ghost_nodes_flow_variable(data, component, params->ghost_nodes, data_bag);

	T2 = MPI_Wtime();
	data_bag->timer->Wtime_vel_boundaries += T2 - T1;
}




/******************************************************************************/
/*
 This function computes the total kinetic energy within the interior of the
 domain:
     K = 0.5 integral (u.u dV)
 */
/******************************************************************************/
void Velocity_compute_total_kinetic_energy( Velocity *u, Velocity *v,
		Velocity *w, MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double ***u_data_bc, ***v_data_bc, ***w_data_bc;
	double dV;
	double G_kinetic_energy;
	double u_, v_, w_;

	// Kinetic energy on the current processor
	G_kinetic_energy = 0.0;

	// Get regular data array for vel data at cell center
	u_data_bc = u->data_bc;
	v_data_bc = v->data_bc;
	w_data_bc = w->data_bc;

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

				// Only include if point is fluid
				if (grid->c_status[k][j][i] == FLUID) {
					/*
					 * Kinetic energy is a physical volume integral, so use the
					 * active 2D/3D control-volume measure here as well.
					 */
					dV = TwodOps_cell_measure_c(grid, params, i, j, k);

					u_ = u_data_bc[k][j][i];
					v_ = v_data_bc[k][j][i];
					w_ = w_data_bc[k][j][i];

					// Add to the kinetic energy on current processor
					G_kinetic_energy += ( u_ * u_ + v_ * v_ + w_ * w_ ) * dV;

				} // if

			} // for i
		} // for j
	} // for k

	G_kinetic_energy *= 0.5;

	u->G_kinetic_energy = G_kinetic_energy;
	v->G_kinetic_energy = G_kinetic_energy;
	w->G_kinetic_energy = G_kinetic_energy;
}




/******************************************************************************/
/*
 This function, summs all the G_.... energies and store them into the W_..
 quantities
 */
/******************************************************************************/
void Velocity_update_world_energies(Velocity *u, Velocity *v, Velocity *w,
		MAC_grid *grid, Parameters *params) {

	double *send_buffer;
	double *recv_buffer;
	int nt;

	// First, compute energies on each processor
//	Velocity_compute_viscous_dissipation_rate(u, v, w, grid, params) ;
	Velocity_compute_total_kinetic_energy(u, v, w, grid, params);

	// Total number of quantities which have to be reduced globally
	nt = 2;
	// First, allocate required memory
	send_buffer = Memory_allocate_1D_array(GVG_DOUBLE, nt);
	recv_buffer = Memory_allocate_1D_array(GVG_DOUBLE, nt);

	// Now, prepare all the data to send at once
	send_buffer[0] = u->G_kinetic_energy;
	send_buffer[1] = u->G_dissipation_rate;

	// Now add up all the local values of the energies and dissipation rate on
	// all processors and send them back to store it in the W_... variables
	MPI_Allreduce( (void *)send_buffer, (void *)recv_buffer, nt, MPI_DOUBLE, MPI_SUM, PCW);

	// Now, store them back into the W_... on each processor
	u->W_kinetic_energy   = recv_buffer[0];
	u->W_dissipation_rate = recv_buffer[1];

	free(send_buffer);
	free(recv_buffer);
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Velocity_compute_bottom_shear_stress(Velocity *u, Velocity *v, Velocity *w,
		Cart3d_bag *data_bag) {

	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i, j_index, k;
	double ***u_data_bc, ***v_data_bc, ***w_data_bc;
	double *yc, *yv;
	double **shear_stress;
	double dudy, dwdy;
	double tau;
	double a_Re;
	double a_Delta_y;
	Indices G_s, G_e, W_e;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	// Get regular data array for vel data
	u_data_bc = u->data_bc;
	v_data_bc = v->data_bc;
	w_data_bc = w->data_bc;

	// Grid coordinates
	yc = grid->yc;
	yv = grid->yv;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	a_Re = 1.0/params->Re;
	shear_stress = u->G_shear_stress_bottom;

	for (k=Ks; k<Ke; k++) {

		for (i=Is; i<Ie; i++) {

			j_index = grid->interface_y_index[k][i];

			// Calculate shear stress only if the current processor contains the
			// bottom boundary
			if ( ( j_index >= Js) && (j_index < Je) ) {

				// Height from the bottom boundary.
				// To avoid unphysical shear stress, we use the second grid (not
				// the first one) above the bottom surface as the velocity
				a_Delta_y = 1.0/(yc[j_index+1] - yv[j_index]);

				dudy = u_data_bc[k][j_index+1][i] * a_Delta_y;

				dwdy = w_data_bc[k][j_index+1][i] * a_Delta_y;

				// Compute total shear stress on the top surface of the bottom
				// boundary
				tau  = -a_Re * sqrt( dudy*dudy + dwdy*dwdy );
			} else {

				tau  = 0.0;
			}

			// Store the shear stress on the current processor
			shear_stress[k][i] = tau;
		} // for i
	} // for k

	// Since the 2D arrays are assumed to have [Y][X] index order, pass always
	// the min, max indices based on this rule.  Start and End indices of the 2D
	// array on the current processor
	G_s.x_index = Is;
	G_s.y_index = Ks;

	G_e.x_index = Ie;
	G_e.y_index = Ke;

	// Total number of grid point in the W_ array
	W_e.x_index = grid->NX;
	W_e.y_index = grid->NZ;

	// Now, sum all the computed shear stresses from all the processors (zero
	// for the ones which don't contain the bottom boundary and send the result
	// to processor zero: Store on W_....
	Communication_reduce_2D_arrays(u->G_shear_stress_bottom, u->W_shear_stress_bottom,
								   &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);

}




/******************************************************************************/
/*
 This function initializes the velocity field to a nonzero value.  Make sure the
 initial velocity field is divergence free
 */
/******************************************************************************/
void Velocity_nonzero_initialize(Velocity *vel, MAC_grid *grid, Parameters *params, Debug_trace *dtrace) {

	int i, j, k;
	double ***data;

	int status;
	char message[500];

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Start index of bottom-left-back corner on current processor
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	// End index of top-right-front corner on current processor
	int Ie = min(NX-1, grid->G_Ie);
	int Je = min(NY-1, grid->G_Je);
	int Ke = min(NZ-1, grid->G_Ke);

#ifdef XPERIODIC
	Ie = grid -> G_Ie;
#endif

	data = vel -> data;

	double *xu = grid -> xu;
	double *yv = grid -> yv;
	double *zw = grid -> zw;

	double *xc = grid -> xc;
	double *yc = grid -> yc;
	double *zc = grid -> zc;

	// Width of y-domain
	double Ly = params->Ly;

	// y-coordinate where velocity profile starts, below velocity will be zero
	double vel_init_y0 = max(params->vel_init_y0, yv[0]);

	// non-dimensional y-coordinate, 0 at vel_init_y0, 1 at top wall
	double yy;

	// Width of velocity profile
	double width_y = fabs(yv[NY-1] - vel_init_y0);

	// Initialize velocity profiles to achieve desired bulk velocity
	double ubulk = params -> ubulk_target;
	if (width_y > 0) ubulk = ubulk * (yv[NY-1] - yv[0]) / width_y;

	// Maximum domain velocity, for setting U_conv_outflow
	double u_max = ubulk;

	//--------------------------------------------------------------------------
	// Uniform flow profile
	//--------------------------------------------------------------------------
	if (params -> vel_init_type == VEL_INIT_UNIFORM) {

		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {
						data[k][j][i] = ubulk;
					}
				}
			}

			u_max = fabs(ubulk);
		}
	}


	if (params -> vel_init_type == VEL_INIT_LEFT_RIGHT) {

		if (vel -> component == 'v' || vel -> component == V_VELOCITY) {

			if (Je==NY-1){
			for (k = Ks; k < Ke; k++) {
					for (i = Is; i < Ie; i++) {
						data[k][Je][i] = params->voutflow;
					}
				}
			}
		}
	}
	//--------------------------------------------------------------------------
	// Linear flow profile
	//--------------------------------------------------------------------------
	else if (params -> vel_init_type == VEL_INIT_LINEAR) {
		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {
#if defined TOP_WALL_VELOCITY && defined BOTTOM_WALL_VELOCITY
						yy = (yc[j] - 0.5*Ly) / Ly;
						data[k][j][i] = 2.0 * ubulk * yy;
#elif defined SHEARED_PERIODIC
						yy = (yc[j] - params->ymin - 0.5*Ly);
						data[k][j][i] = params->shear_sigma_u * yy;// + 0.001*(-0.5+(double)rand()/(double)RAND_MAX);
#else
						yy = (yc[j] - vel_init_y0) / width_y;
						if (yy > 0) {
							data[k][j][i] = 2.0 * ubulk * yy;
						}
#endif
					}
				}
			}

			u_max = 2.0 * fabs(ubulk);
		}
	}
	//--------------------------------------------------------------------------
	// Poiseuille (parabolic) flow profile
	//--------------------------------------------------------------------------
	else if (params -> vel_init_type == VEL_INIT_POISEUILLE) {

		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {

						yy = (yc[j] - vel_init_y0) / width_y;
						if (yy > 0)
							data[k][j][i] = -6.0 * ubulk * yy * (yy - 1.0);
					}
				}
			}

			u_max = 1.5 * fabs(ubulk);
		}
	}

	//
	//

	else if (params -> vel_init_type == VEL_INIT_JET) {

		double y1 = Ly/2.0 - 0.5;
		double y2 = Ly/2.0 + 0.5;
		double delta = params->delta_jet;


		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {

							data[k][j][i] = ((double)rand()/(double)RAND_MAX)/1e2;  //0.5*(erf((yc[j]-y1)/delta)-erf((yc[j]-y2)/delta)) +  ;
					}
				}
			}

		}
	}

	else if (params -> vel_init_type == VEL_INIT_JET_3D) {

		double y1 = Ly/2.0 - 0.5;
		double y2 = Ly/2.0 + 0.5;
		double delta = params->delta_jet;


		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {

							data[k][j][i] = ((double)rand()/(double)RAND_MAX)/1e2;  //0.5*(erf((yc[j]-y1)/delta)-erf((yc[j]-y2)/delta)) +  ;
					}
				}
			}

		}
	}

	else if (params -> vel_init_type == VEL_INIT_NOISE) {

		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {

							data[k][j][i] = ((double)rand()/(double)RAND_MAX)/1e2;  //0.5*(erf((yc[j]-y1)/delta)-erf((yc[j]-y2)/delta)) +  ;
					}
				}
			}

		}
	}


	//--------------------------------------------------------------------------
	// Rotational shear flow profile
	//--------------------------------------------------------------------------
	else if (params -> vel_init_type == VEL_INIT_ROT_SHEAR) {

		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {

						data[k][j][i] = 1.0 / params -> G + yc[j] - (0.5 * (yc[grid->NY-2] - yc[0]) + yc[0]);
					}
				}
			}
		}
	}
	//--------------------------------------------------------------------------
	// Precursor simulation flow profile
	//--------------------------------------------------------------------------
	else if (params -> vel_init_type == VEL_INIT_PRECURSOR) {
		int  jj;
		FILE *fp ;
		char buff[255];
		char filename[50];
		float y_in[301], vel_in[301], t_in;
		double dvel;
		double  **inflow_n;
		inflow_n = vel -> inflow_n;

		// Open precursor file
		params->ninflow = 500;
		if (vel -> component == 'u') {
			sprintf(filename, "./inflow/u_in_%04d.dat", params->ninflow);
		}
		else if (vel -> component == 'v') {
			sprintf(filename, "./inflow/v_in_%04d.dat", params->ninflow);
		}
		else if (vel -> component == 'w') {
			sprintf(filename, "./inflow/w_in_%04d.dat", params->ninflow);
		}

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

		// Initialize field
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
					data[k][j][i] = inflow_n[k][j];
				if (vel -> component == 'v' || vel -> component == 'w') data[k][j][i] = 0.0;
				}
			}
		}

		sprintf(message, "Could not read initial velocity data");
		Display_assert_error(status, message, params, DTRACE("Display_assert_error"));
	}

	else if (params -> vel_init_type == VEL_INIT_TGV) {
		double U0 = 1.0;
		double lam = 1.0;//2*PI;

		if (vel -> component == 'u') {

			if (Ie == NX-1)
				Ie = NX;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {
							data[k][j][i] = U0*sin(lam*xu[i])*cos(lam*yc[j])*cos(lam*zc[k]);
					}
				}
			}
		}
		else if (vel -> component == 'v') {

			if (Je == NY-1)
				Je = NY;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {
							data[k][j][i] = -U0*cos(lam*xc[i])*sin(lam*yv[j])*cos(lam*zc[k]);
					}
				}
			}
			printf("Hello\n");
		}
		else{
			if (Ke == NZ-1)
				Ke = NZ;

			for (k = Ks; k < Ke; k++) {
				for (j = Js; j < Je; j++) {
					for (i = Is; i < Ie; i++) {
							data[k][j][i] = 0.0;
					}
				}
			}
		}
	}


}




/******************************************************************************/
/*
 Calculate ubulk
 */
/******************************************************************************/
double Velocity_ubulk(Velocity *uvel, Cart3d_bag *data_bag) {

	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i, j, k;
	double ***u_data;
	double *yv;
	double ub, W_ub;
	int NX, NY, NZ;

	MAC_grid *grid = data_bag -> grid;

	NX = grid -> NX;
	NY = grid -> NY;
	NZ = grid -> NZ;

	// Start index of bottom-left-back corner on current processor
	Is = max(1,grid -> G_Is); // i=0 not included
	Js = grid -> G_Js;
	Ks = grid -> G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid -> G_Ie;
	Je = min(grid->G_Je, NY-1);
	Ke = min(grid->G_Ke, NZ-1);

	u_data = uvel -> data;

	yv = grid -> yv;
	ub   = 0.0;
	W_ub = 0.0;

#ifdef LAG_PARTICLE_RESOLVED
	double vf_cell;
	double ***vfu = data_bag -> lag -> ng_vfu;
#endif

	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {

#ifdef LAG_PARTICLE_RESOLVED
				vf_cell = min(vfu[k][j][i], 1.0);
				ub += (1 - vf_cell) * u_data[k][j][i] * (yv[j+1] - yv[j]);
#else
				ub += u_data[k][j][i] * (yv[j+1] - yv[j]);
#endif

			} // for i
		} // for j
	} // for k

	ub = ub / (NX * (NZ - 1.0));
	ub = ub / (yv[NY-1] - yv[0]);
	MPI_Allreduce(&ub, &W_ub, 1, MPI_DOUBLE, MPI_SUM, PCW);

	return W_ub;

}




/******************************************************************************/
/*
 Calculate dp_dx
 */
/******************************************************************************/
void Velocity_calculate_dpdx(Velocity *uvel, Cart3d_bag *data_bag) {

	FILE *fid;
	FILE *fid_osc;
	double dt, a_dt;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	params -> ubulk_old = params -> ubulk;
	params -> ubulk     = Velocity_ubulk(uvel, data_bag);

	dt = params -> dt;
	a_dt = 1.0 / dt;
	params -> dp_dx_old = params -> dp_dx;

	#ifdef CONSTANT_MASSFLUX 
		params -> dp_dx     = params -> dp_dx_old
	                      	+ 2.0 * (params->ubulk - params->ubulk_target) / params -> dt
	                      	- (params->ubulk_old - params->ubulk_target) / params -> dt_old;
	#endif

	if (params -> rank == 0) {
		if (params -> time == 0) {
			fid = fopen("dpdx_history.dat","w"); }
		else {
			fid = fopen("dpdx_history.dat","a"); }

		#ifdef STOKES_2ND_PROBLEM
			fprintf(fid, "%8d %e %e %20.12e %20.12e %20.12e %20.12e %20.12e\n", params->ntime,
			params->time,dt, params->dp_dx, params->dp_dx_old, params->ubulk, params->ubulk_old, params->u_oscillation);
			fclose(fid);
		#else
			fprintf(fid, "%8d %e %e %20.12e %20.12e %20.12e %20.12e\n", params->ntime,
				params->time,dt, params->dp_dx, params->dp_dx_old, params->ubulk, params->ubulk_old);
			fclose(fid);
		#endif
	}

	#ifdef OSCILLATION
		if (params -> rank == 0) {
			if (params -> time == 0) {
				fid_osc = fopen("oscillation.dat","w"); }
			else {
				fid_osc = fopen("oscillation.dat","a"); }

		fprintf(fid_osc, "%8f %.20g\n", params->time, params -> oscillation);
		fclose(fid_osc);	
		}
	#endif

return;
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Velocity_wall_shear(Cart3d_bag *data_bag) {

	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i, j, k;
	int in, jn, kn;
	int ib_index;
	double ***u_data_bc, ***v_data_bc, ***w_data_bc;
	double *xc, *yc, *zc;
	double *yv;
	double **shear_stress;
	double tau, d;
	double Re;
	Indices G_s, G_e, W_e;
	double V_dot_n;
	double nx, ny, nz;
	double u_tan, v_tan, w_tan;
	double u_int, v_int, w_int;
	double **u_shear, **v_shear, **w_shear;
	PointType *p_int;
	Immersed *c_immersed;
	ImmersedNode *ib_node;
	Indices index;
	int pnodes;
	int N, g;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	Velocity *u = data_bag -> u;
	Velocity *v = data_bag -> v;
	Velocity *w = data_bag -> w;

	// Allocate memory for the interpolation point
	p_int = (PointType *)calloc(1, sizeof(PointType));

	// Update the values of the ghost nodes based on a BOX stencil
	pnodes = params->ghost_nodes;
	Communication_update_ghost_nodes_flow_variable(u->data_bc, CONCENTRATION_PERTURBATION, pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(v->data_bc, CONCENTRATION_PERTURBATION, pnodes, data_bag);
	Communication_update_ghost_nodes_flow_variable(w->data_bc, CONCENTRATION_PERTURBATION, pnodes, data_bag);

	// Get regular data array for vel data
	u_data_bc = u->data_bc;
	v_data_bc = v->data_bc;
	w_data_bc = w->data_bc;

	// Grid coordinates
	xc = grid->xc;
	yc = grid->yc;
	zc = grid->zc;
	yv = grid->yv;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	Re = params->Re;
	shear_stress = u->G_shear_stress_bottom;

	// Tangential velocity on the bottom
	u_shear = u->G_u_shear;
	v_shear = u->G_v_shear;
	w_shear = u->G_w_shear;

	j = 0;
	for (k=Ks; k<Ke; k++) {
		for (i=Is; i<Ie; i++) {

			u_tan = u_data_bc[k][j][i];
			v_tan = 0.0;
			w_tan = w_data_bc[k][j][i];
			d = yc[j] - yv[j];

			tau = sqrt(u_tan*u_tan + v_tan*v_tan + w_tan*w_tan) / (d * Re);

			// Store the shear stress on the current processor
			shear_stress[k][i] = tau;

			// Tangential velocity
			u_shear[k][i] = u_tan;
			v_shear[k][i] = v_tan;
			w_shear[k][i] = w_tan;

		} // for i
	} // for k


#ifdef IMMERSED_BOUNDARY
	c_immersed = grid->c_immersed;
	N = grid->c_immersed->N;

	for (g=0; g<N; g++) {

		// Get the current immersed node
		ib_node = Immersed_get_ib_node(c_immersed, g);

		p_int->x = ib_node->intersection_point.x;
		p_int->y = ib_node->intersection_point.y;
		p_int->z = ib_node->intersection_point.z;

		i = ib_node->im_index.x_index;
		k = ib_node->im_index.z_index;

		// distance between interpolation point and surface
		d = MyMath_get_point_point_distance(p_int, &ib_node->boundary_point);

		// Use a trilinear interpolation to find the velocity (u,v,w)
		u_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, u_data_bc, grid, 'c');
		v_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, v_data_bc, grid, 'c');
		w_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, w_data_bc, grid, 'c');

		nx = ib_node->n.vx;
		ny = ib_node->n.vy;
		nz = ib_node->n.vz;

		// magnitude of the normal velocity vector
		V_dot_n = u_int*nx + v_int*ny + w_int*nz;

		// Tangential velocity vector = V - V_norm
		u_tan = u_int - V_dot_n*nx;
		v_tan = v_int - V_dot_n*ny;
		w_tan = w_int - V_dot_n*nz;


		tau = sqrt(u_tan*u_tan + v_tan*v_tan + w_tan*w_tan) / (d * Re);

		// Store the shear stress on the current processor
		shear_stress[k][i] = tau;

		// Tangential velocity
		u_shear[k][i] = u_tan;
		v_shear[k][i] = v_tan;
		w_shear[k][i] = w_tan;
	}
#endif // IMMERSED_BOUNDARY

	// Since the 2D arrays are assumed to have [Y][X] index order, pass always
	// the min, max indices based on this rule.  Start and End indices of the 2D
	// array on the current processor
	G_s.x_index = Is;
	G_s.y_index = Ks;

	G_e.x_index = Ie;
	G_e.y_index = Ke;

	// Total number of grid point in the W_ array
	W_e.x_index = grid->NX;
	W_e.y_index = grid->NZ;

	// Now, sum all the computed shear stresses from all the processors (zero
	// for the ones which don't contain the bottom boundary and send the result
	// to processor zero: Store on W_...
	Communication_reduce_2D_arrays(u->G_shear_stress_bottom, u->W_shear_stress_bottom, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);


	// Get the bottom tangential velocity
	Communication_reduce_2D_arrays(u->G_u_shear, u->W_u_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);
	Communication_reduce_2D_arrays(u->G_v_shear, u->W_v_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);
	Communication_reduce_2D_arrays(u->G_w_shear, u->W_w_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);

	free(p_int);
}




/******************************************************************************/
/*
 */
/******************************************************************************/
void Velocity_u_streak(Cart3d_bag *data_bag) {

	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i, j, k;
	int j_plane;
	double ***u_data_bc;
	Indices G_s, G_e, W_e;
	double **u_streak;
	double *u_ave;
	int NX, NY, NZ;
	int pnodes;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;
	pnodes = params->ghost_nodes;

	Velocity *u = data_bag -> u;

	// Update the values of the ghost nodes based on a BOX stencil
	Communication_update_ghost_nodes_flow_variable(u->data_bc, CONCENTRATION_PERTURBATION, pnodes, data_bag);

	// Get regular data array for vel data
	u_data_bc = u->data_bc;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	u_streak = u->G_u_streak;

	j_plane = 5;
	if ((Js <= j_plane) && (Je >=j_plane) ) {

		for (k=Ks; k<Ke; k++) {
			for (i=Is; i<Ie; i++) {

					u_streak[k][i] = u_data_bc[k][j_plane][i];

			} // for i
		} // for k
	}

	// Since the 2D arrays are assumed to have [Y][X] index order, pass always
	// the min, max indices based on this rule.  Start and End indices of the 2D
	// array on the current processor
	G_s.x_index = Is;
	G_s.y_index = Ks;

	G_e.x_index = Ie;
	G_e.y_index = Ke;

	// Total number of grid point in the W_ array
	W_e.x_index = grid->NX;
	W_e.y_index = grid->NZ;

	// Get the bottom velocity
	Communication_reduce_2D_arrays(u->G_u_streak, u->W_u_streak, &G_s, &G_e,
	                               &W_e, REDUCE_TO_MASTER, data_bag);

	if (params->rank == 0) {

		u_ave = (double *) calloc(NX, sizeof(double));
		for (i=0; i<NX; i++) {
			u_ave[i] = 0.0;
		}

		for (k=0; k<NZ-1; k++) {
			for (i=0; i<NX-1; i++) {
				u_ave[i]  += u->W_u_streak[k][i];
			}
		}

		for (i=0; i<NX-1; i++) {
			u_ave[i] = u_ave[i]/(NZ-1.);
		}


		for (k=0; k<NZ-1; k++) {
			for (i=0; i<NX-1; i++) {
				u->W_u_streak[k][i] -= u_ave[i];
			}
		}

		for (i=0; i<NX-1; i++) {
			u->W_u_streak[NZ-1][i] = u->W_u_streak[NZ-2][i];
		}

		for (k=0; k<NZ; k++) {
			u->W_u_streak[k][NX-1] = u->W_u_streak[k][NX-2];
		}

		free(u_ave);
	}

}


#include "lsolver/direct_solver.c"
#include "lsolver/msolve_direct.c"
#ifdef CG_SOLVE
	#include "lsolver/msolve_cg.c"
#endif
#ifdef BICG_SOLVE
	#include "lsolver/msolve_bicg.c"
#endif
