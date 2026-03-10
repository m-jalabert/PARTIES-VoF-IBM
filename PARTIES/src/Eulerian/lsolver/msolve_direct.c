
//#ifdef FULLY_EXPLICIT
//	#define Velocity_solve Velocity_solve_explicit
//#else
//	#define Velocity_solve Velocity_solve_semi_implicit
//#endif

/******************************************************************************/
/*
 Calculate the coefficients of tri-diagonal matrix resulting from dicretization
 of 1-D diffusion equation in y-direction.
 */
/******************************************************************************/
void Velocity_calc_coeff(Velocity *vel, MAC_grid *grid, Parameters *params) {

	int i, j, k;
	int NX, NY, NZ;
	int ierr;
	double *idx_c, *idy_c, *idz_c;
	double *idx_u, *idy_v, *idz_w;
	double *as, *an, *ap;
	double *rhsw, *xw;
	double *modwaveksq, *modwaveisq;
	double idx, idz;
	int NJ;
	double Re, iRe;

	NX  = grid -> NX;
	NY  = grid -> NY;
	NZ  = grid -> NZ;

	/*------------------------------------------------------------------------*/
	/*
	     -ap: diagonal
	     -as: south
	     -an: north
	 */
	/*------------------------------------------------------------------------*/
	idy_v = grid -> idy_v;
	idy_c = grid -> idy_c;
	idx_u = grid -> idx_u;
	idz_w = grid -> idz_w;

	if (vel->component == 'v') {
		NJ = NY;
	}
	else {
		NJ = NY-1;
	}
	as = (double *) malloc( NJ * sizeof(double) );
	an = (double *) malloc( NJ * sizeof(double) );
	ap = (double *) malloc( NJ * sizeof(double) );
	vel -> asw  = (double *) malloc( NJ * sizeof(double) );
	vel -> anw  = (double *) malloc( NJ * sizeof(double) );
	vel -> apw  = (double *) malloc( NJ * sizeof(double) );
	vel -> xw   = (double *) malloc( NJ * sizeof(double) );
	vel -> rhsw = (double *) malloc( NJ * sizeof(double) );

	Re = params -> Re;
	iRe = 1.0 / Re;
#ifdef LES
	// Coefficients are multiplied by (nut+iRe) in Velocity_solve for LES
	iRe = 1.0;
#endif


	// Coefficient from discretization of d2p/dy2 with Neumann boundary
	// condition at both walls/free surfaces

	if (vel -> component == 'v') {
		as[0] = 0.0;
		ap[0] = 1.0;
		an[0] = 0.0;
		for (j = 1; j < NY-1; j++) {
			as[j] = iRe * idy_v[j-1] * idy_c[j-1];
			an[j] = iRe * idy_v[j] * idy_c[j-1];
			ap[j] = -as[j] - an[j];
		}
		as[NY-1] = 0.0;
		ap[NY-1] = 1.0;
		an[NY-1] = 0.0;
	}
	else {
		for (j = 0; j < NY-1; j++) {
			if (j > 0 )   {
				as[j] = iRe * idy_v[j] * idy_c[j-1];
			}
			else {
				as[j] = iRe * idy_v[j] * idy_v[0];
			}
			an[j] = iRe * idy_v[j] * idy_c[j];
			ap[j] = -as[j] - an[j];
		}
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
		ap[0] = -an[0] - 2.0 * as[0];
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
		ap[0] = -an[0];
#endif

#ifdef TOP_WALL_VELOCITY_NOSLIP
		ap[NY-2] = -as[NY-2] - 2.0 * an[NY-2];
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
		ap[NY-2] = -as[NY-2];
#endif
		as[0] = 0.0;
		an[NY-2] = 0.0;
	}

	vel -> as = as;
	vel -> ap = ap;
	vel -> an = an;


	return;

}



/******************************************************************************/
/*
 Create FFTW plans, transpose plans and calculate the coefficient of the matrix
 for the 1-D Helmholtz equation.
 */
/******************************************************************************/
void Velocity_setup_lsys_accounting_geometry(Velocity *vel, MAC_grid *grid,
		Parameters *params) {

	int ierr;
	int maxsize;

	Velocity_calc_coeff(vel, grid, params);

}




/******************************************************************************/
/*
 This function solves the linear system to find velocity
 */
/******************************************************************************/
int Velocity_solve_explicit(Velocity *vel, Cart3d_bag *data_bag) {

	int iters;
	int i, j, k;
	int i_start, j_start ,k_start;
	int i_end, j_end, k_end;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = min(grid->G_Ie, grid->NX-1);
	int Je = min(grid->G_Je, grid->NY-1);
	int Ke = min(grid->G_Ke, grid->NZ-1);

	// We are not solving for the nodes on the boundaries unless periodic
	if (vel -> component == 'u') {
		Is = max(Is, 1);
#ifdef XPERIODIC
		Ie = grid->G_Ie;
#endif
	}
	if (vel -> component == 'v') {
		Js = max(Js, 1);
#ifdef YPERIODIC
		Je = grid->G_Je;
#endif
	}
	if (vel -> component == 'w') {
		Ks = max(Ks, 1);
#ifdef ZPERIODIC
		Ke = grid->G_Ke;
#endif
	}

	const double BET[] = {BETA};
	double dtimeb = BET[params -> which_stage] * params -> dt;

	double ***rhs = vel -> ng_rhs;
	double ***vel_data = vel -> data;

	#ifdef VOF_PLIC
		VolumeFraction *vof = data_bag->vof;
		double ***rho       = vof->rho;   // cell-centered current density (rho^k)
	#endif

	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				#ifndef VOF_PLIC

								/* Original constant-density behavior */
								vel_data[k][j][i] = rhs[k][j][i] * dtimeb;

				#else
								/* 
								* Variable-density correction:
								* rhs is interpreted as (1 / (alpha_k * dt)) * (rho^k u^* - rho^{k-1} u^{k-1}) + ...
								* After multiplying by dtimeb = beta_k * dt, we have something proportional
								* to rho^k u^*. To recover u^*, we must divide by the current density
								* at the velocity location (face).
								*/

								double rho_face = 1.0;

								if (vel->component == 'u') {
									/* u at x-face between cells i-1 and i */
									int icL = i-1;
									int icR = i;
									/* indices should already be valid due to Is >= 1; ghosts hold BCs */
									rho_face = 0.5 * (rho[k][j][icL] + rho[k][j][icR]);
								}
								else if (vel->component == 'v') {
									/* v at y-face between cells j-1 and j */
									int jcB = j-1;
									int jcT = j;
									rho_face = 0.5 * (rho[k][jcB][i] + rho[k][jcT][i]);
								}
								else if (vel->component == 'w') {
									/* w at z-face between cells k-1 and k */
									int kcB = k-1;
									int kcT = k;
									rho_face = 0.5 * (rho[kcB][j][i] + rho[kcT][j][i]);
								}

								/* Safety: avoid division by zero or negative density (should not happen) */
								if (rho_face <= 0.0) {
									rho_face = 1.0;
								}

								vel_data[k][j][i] = (rhs[k][j][i] * dtimeb) / rho_face;
				#endif
			}
		}
	}

	iters = 0;
	return iters;
}




/******************************************************************************/
/*
 This function solves the linear system to find velocity
 */
/******************************************************************************/
int Velocity_solve_semi_implicit(Velocity *vel, Cart3d_bag *data_bag) {

	int iters;
	int ierr;
	double rnorm;
	double *data_t;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int Ijs, Ije;
	int i_start, j_start, k_start;
	int i_end, j_end, k_end;
	int Ijs_start, Ije_end;
	int NX, NY, NZ;
	int count, index, ind1, ind2, ind3, ind4, ind5, ind6;
	double ***rhs, ***vel_data;
	int i, j, k;
	lsolver_transpose *lsolver_remap;
	double nutN, nutS, ***nut;
	double as_top, an_top, as_bot, an_bot;
	double Re, iRe;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;


	lsolver_remap = params->lsolver_remap;
	// Start index of bottom-left-back corner on current processor
	Is = grid -> G_Is;
	Js = grid -> G_Js;
	Ks = grid -> G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid -> G_Ie;
	Je = grid -> G_Je;
	Ke = grid -> G_Ke;

	NX = grid -> NX;
	NY = grid -> NY;
	NZ = grid -> NZ;

	data_t = lsolver_remap->data_transpose;
	Ijs = lsolver_remap -> Ijs;
	Ije = lsolver_remap -> Ije;

	const double BET[] = {BETA};
	double dtimeb = BET[params -> which_stage] * params -> dt;

	Re = params -> Re;
	iRe = 1.0 / Re;
#ifdef LES
	nut = data_bag -> smag -> nut;
#endif
	rhs = vel -> ng_rhs;
	vel_data = vel -> data;


	if (params -> NPY != 1) {
		//----------------------------------------------------------------------
		// Copy the RHS of poisson equation to 1-D array data_t
		//----------------------------------------------------------------------
		count = 0;
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
					data_t[count] = rhs[k][j][i];
					count = count + 1 ;
				}
			}
		}

		//----------------------------------------------------------------------
		// Transpose the data so that tridiagonal equation can be solved along
		// y-direction
		//----------------------------------------------------------------------
		remap_3d(data_t, data_t, NULL, lsolver_remap->remap_plan_origtoyline);
	}

	/*------------------------------------------------------------------------*/
	/*
	 Tridiagonal solve
	 */
	/*------------------------------------------------------------------------*/
	ind1 = (Ije - Ijs) * NY;
	ind2 = (Ije - Ijs);


	//--------------------------------------------------------------------------
	// 'v' velocity
	//--------------------------------------------------------------------------
	if (vel->component=='v') {
		Ije_end = min(Ije, NX-1);
		k_end = min(Ke, NZ-1);

		for (k = Ks; k < k_end; k++) {

			ind3 = ind1 * (k - Ks);

			for (i = Ijs; i < Ije_end; i++) {

				ind4 = ind3 + (i - Ijs);

				for (j = 1; j < NY; j++) {
#ifndef LES
					vel -> asw[j] = -vel -> as[j];
					vel -> anw[j] = -vel -> an[j];
					vel -> apw[j] = 1.0 / dtimeb - vel -> ap[j];
#else
					nutN = nut[k][j][i];
					nutS = nut[k][j-1][i];
					nutN = nutN + iRe;
					nutS = nutS + iRe;
					vel -> anw[j] = -vel -> an[j] * nutN;
					vel -> asw[j] = -vel -> as[j] * nutS;
					vel -> apw[j] = 1.0 / dtimeb - vel -> asw[j] - vel -> anw[j];
#endif
				} // for j

				if (params -> NPY != 1) {
					for (j = 1; j < NY; j++) {
						vel -> rhsw[j] = data_t[ind4 + j * ind2];
					}
				}
				else {
					for (j = 1; j < NY; j++) {
						vel -> rhsw[j] = rhs[k][j][i];
					}
				}

				vel -> rhsw[0]    = 0.0;
				vel -> asw[0]     = 0.0;
				vel -> anw[0]     = 0.0;
				vel -> apw[0]     = 1.0;

				vel -> rhsw[NY-1] = 0.0;
				vel -> asw[NY-1]  = 0.0;
				vel -> anw[NY-1]  = 0.0;
				vel -> apw[NY-1]  = 1.0;
				ltridiag(vel->asw, vel->apw, vel->anw, vel->rhsw, vel->xw, NY);

				if (params -> NPY != 1) {
					for (j = 0; j < NY; j++) {
						data_t[ind4 + j * ind2] = vel -> xw[j];
					}
				}
				else {
					for (j = 0; j < NY; j++) {
						vel_data[k][j][i] = vel -> xw[j];
					}
				}

			} // for i
		} // for k
	} // if v


	//--------------------------------------------------------------------------
	// 'u' velocity
	//--------------------------------------------------------------------------
	else if (vel -> component == 'u') {
		as_bot = grid -> idy_v[0] * grid -> idy_v[0];
		an_bot = grid -> idy_v[0] * grid -> idy_c[0];

		j = NY-2;
		as_top = grid -> idy_v[j] * grid -> idy_c[j-1];
		an_top = grid -> idy_v[j] * grid -> idy_c[j];

		Ijs_start = max(Ijs, 1);
		Ije_end = min(Ije, NX-1);
#ifdef XPERIODIC
		Ije_end = Ije;
#endif
		k_end = min(Ke, NZ-1);

		for (k = Ks; k < k_end; k++) {

			ind3 = ind1 * (k - Ks);

			for (i = Ijs_start; i < Ije_end; i++) {

				ind4 = ind3 + (i - Ijs);

				for (j = 0; j < NY-1; j++) {
#ifndef LES
					vel -> asw[j] = -vel -> as[j];
					vel -> anw[j] = -vel -> an[j];
					vel -> apw[j] = 1.0 / dtimeb - vel -> ap[j];
#else // LES
					nutN = 0.25 * ( nut[k][j][i]   + nut[k][j][i-1]
					              + nut[k][j+1][i] + nut[k][j+1][i-1] );
					if (j != 0) {
						nutS = 0.25 * ( nut[k][j][i]   + nut[k][j][i-1]
						              + nut[k][j-1][i] + nut[k][j-1][i-1] );
					}
					else {
	#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
						// Zero Eddy viscosity on the wall
						nutS =  0.0;
	#else
						nutS = 0.5 * ( nut[k][j][i] + nut[k][j][i-1] );
	#endif
					}

					nutN = nutN + iRe;
					nutS = nutS + iRe;

					vel -> asw[j] = -vel -> as[j] * nutS;
					vel -> anw[j] = -vel -> an[j] * nutN;
					vel -> apw[j] = 1.0 / dtimeb - vel -> asw[j] - vel -> anw[j];
#endif // not LES
					vel -> rhsw[j] = data_t[ind4 + j * ind2];
				} // for j

				if (params -> NPY != 1) {
					for (j = 0; j<NY-1; j++) {
						vel -> rhsw[j] = data_t[ind4 + j * ind2];
					}
				}
				else {
					for (j = 0; j<NY-1; j++) {
						vel -> rhsw[j] = rhs[k][j][i];
					}
				}
#ifdef LES
				j = 0;
				nutN = 0.25 * ( nut[k][j][i]   + nut[k][j][i-1]
				              + nut[k][j+1][i] + nut[k][j+1][i-1] );
				nutN = nutN + iRe;

				vel -> apw[0] = 1.0 / dtimeb + an_bot * nutN;
	#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
				vel -> apw[0] = vel -> apw[0] + 2.0 * as_bot * iRe;
	#endif
				j = NY - 2;
				nutS = 0.25 * ( nut[k][j][i]   + nut[k][j][i-1]
				              + nut[k][j-1][i] + nut[k][j-1][i-1] );
				nutS = nutS + iRe;

				vel -> apw[NY-2] = 1.0 / dtimeb + as_top * nutS;
	#ifdef TOP_WALL_VELOCITY_NOSLIP
				vel -> apw[NY-2] = vel -> apw[NY-2] + 2.0 * an_top * iRe;
	#endif
#endif // LES

				ltridiag(vel->asw, vel->apw, vel->anw, vel->rhsw, vel->xw, NY-1);
				if (params -> NPY != 1) {
					for (j = 0; j < NY-1; j++) {
						data_t[ind4 + j * ind2] = vel -> xw[j];
					}
					data_t[ind4 + (NY-1) * ind2] = vel -> xw[NY-2];
				}
				else {
					for (j = 0; j < NY-1; j++) {
						vel_data[k][j][i] = vel -> xw[j];
					}
					vel_data[k][NY-1][i] = vel -> xw[NY-2];
				}
			} // for i
		} // for k
	} // if u


	//--------------------------------------------------------------------------
	// 'w' velocity
	//--------------------------------------------------------------------------
	else {
		as_bot = grid -> idy_v[0] * grid -> idy_v[0];
		an_bot = grid -> idy_v[0] * grid -> idy_c[0];

		j = NY - 2;
		as_top = grid -> idy_v[j] * grid -> idy_c[j-1];
		an_top = grid -> idy_v[j] * grid -> idy_c[j];

		Ije_end = min(Ije, NX-1);
		k_end = min(Ke, NZ-1);
#ifdef ZPERIODIC
		k_end = Ke;
#endif

		for (k = Ks; k < k_end; k++) {

			ind3 = ind1*(k-Ks);

			for (i = Ijs; i < Ije_end; i++) {

				ind4 = ind3 + (i-Ijs);

				for (j = 0; j < NY-1; j++) {
#ifndef LES
					vel -> asw[j] = -vel -> as[j];
					vel -> anw[j] = -vel -> an[j];
					vel -> apw[j] = 1.0 / dtimeb - vel -> ap[j];
#else // LES

					nutN = 0.25 * ( nut[k][j][i]   + nut[k-1][j][i]
					              + nut[k][j+1][i] + nut[k-1][j+1][i] );

					if (j != 0) {
						nutS = 0.25 * ( nut[k][j][i]   + nut[k-1][j][i]
						              + nut[k][j-1][i] + nut[k-1][j-1][i] );
					}
					else {
						nutS = 0.0;
					}

					nutN = nutN + iRe;
					nutS = nutS + iRe;
					vel -> asw[j] = -vel -> as[j] * nutS;
					vel -> anw[j] = -vel -> an[j] * nutN;
					vel -> apw[j] = 1.0 / dtimeb - vel -> asw[j] - vel -> anw[j];

#endif // not LES
				} // for j

				if (params->NPY != 1) {
					for (j = 0; j < NY-1; j++) {
						vel -> rhsw[j] = data_t[ind4 + j * ind2];
					}
				}
				else {
					for (j = 0; j < NY-1; j++) {
						vel -> rhsw[j] = rhs[k][j][i];
					}
				}

#ifdef LES
				j = 0;
				nutN = 0.25 * ( nut[k][j][i]   + nut[k-1][j][i]
				              + nut[k][j+1][i] + nut[k-1][j+1][i] );
				nutN = nutN + iRe;

				vel -> apw[0] = 1.0 / dtimeb + an_bot * nutN;
	#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
				vel -> apw[0] = vel -> apw[0] + 2.0 * as_bot * iRe;
	#endif

				j = NY-2;
				nutS = 0.25 * ( nut[k][j][i]   + nut[k-1][j][i]
				              + nut[k][j-1][i] + nut[k-1][j-1][i] );
				nutS = nutS + iRe;

				vel -> apw[NY-2] = 1.0 / dtimeb + as_top * nutS;
	#ifdef TOP_WALL_VELOCITY_NOSLIP
				vel -> apw[NY-2] = vel -> apw[NY-2] + 2.0 * an_top * iRe;
	#endif
#endif // LES

				ltridiag(vel->asw, vel->apw, vel->anw, vel->rhsw, vel->xw, NY-1);
				if (params -> NPY != 1) {
					for (j = 0; j < NY-1; j++) {
						data_t[ind4 + j * ind2] = vel -> xw[j];
					}
					data_t[ind4 + (NY-1) * ind2] = vel -> xw[NY-2];
				}
				else {
					for (j = 0; j < NY-1; j++) {
						vel_data[k][j][i] = vel -> xw[j];
					}
					vel_data[k][NY-1][i] = vel -> xw[NY-2];
				}
			} // for i
		} // for k
	} // else (w)


	if (params -> NPY != 1) {
		//----------------------------------------------------------------------
		// Transpose the data to original layout
		//----------------------------------------------------------------------
		remap_3d(data_t, data_t, NULL, lsolver_remap->remap_plan_ylinetoorig);

		i_start = Is;
		i_end = min(Ie, NX-1);
		if (vel -> component == 'u') {
			i_start = max(Is, 1);
#ifdef XPERIODIC
			i_end = Ie;
#endif
		}
		k_start = Ks;
		k_end = min(Ke, NZ-1);
		if (vel -> component == 'w') {
			k_start = max(Ks, 1);
#ifdef ZPERIODIC
			k_end = Ke;
#endif
		}


		//----------------------------------------------------------------------
		// Copy the solution to p_data array
		//----------------------------------------------------------------------
		count = 0;
		if (k_start != Ks)
			count = (Je - Js) * (Ie - Is);

		for (k = k_start; k < k_end; k++) {
			for (j = Js; j < Je; j++) {

				if (Is != i_start)
					count = count + 1 ;

				for (i = i_start; i < i_end; i++) {
					vel_data[k][j][i] = data_t[count];
					count = count + 1 ;
				}

				if (Ie != i_end)
					count = count + 1;
			}
		}
	}

	iters = 0;
	return iters;
}
