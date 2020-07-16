
/******************************************************************************/
/*
 This function calculates the new concentration field in the fully explicit form
 */
/******************************************************************************/
int Conc_solve_fully_explicit(int iconc, Cart3d_bag *data_bag) {

	int iters;
	int i, j, k;
	int i_start, j_start ,k_start;
	int i_end, j_end, k_end;


	double coef, b;

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];

	const double BET[] = {BETA};
	double dtimeb = BET[params -> which_stage] * params -> dt;

	double Pe = c -> Pe;
	double iPe = 1.0 / Pe;



	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = min(grid->G_Ie, grid->NX-1);
	int Je = min(grid->G_Je, grid->NY-1);
	int Ke = min(grid->G_Ke, grid->NZ-1);



	double ***rhs = c-> ng_rhs;
	double ***c_data = c -> data;

	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				c_data[k][j][i] = rhs[k][j][i] * dtimeb;
			}
		}
	}


	return 1;
}


/******************************************************************************/
/*
 Solves the tridiagonal system with coefficients
 
     a, b, c and r for a vector u
 
 where a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1] = r[i]
 
 Assumes a[0] and c[N-1] are zero (i.e does not work for periodic boundary)
 */
/******************************************************************************/
void cltridiag(double *a, double *b, double *c, double *r, double *u, int N){
	
	double *gam, bet;
	int i;
	
	gam = (double *) malloc(N*sizeof(double));

	bet = b[0];
	u[0] = r[0]/bet;

	for (i=1;i<N;i++){
		gam[i] = c[i-1]/bet;
		bet = b[i] - a[i]*gam[i];
		u[i] = (r[i]-a[i]*u[i-1])/bet;
	}
	for (i=N-2;i>=0; i--){
		u[i] = u[i]-gam[i+1]*u[i+1];
	}

	free(gam);
	return;
}




/******************************************************************************/
/*
 Calculate the coefficients of tri-diagonal matrix resulting from dicretization
 of 1-D diffusion equation in y-direction.
 */
/******************************************************************************/
void Conc_calc_coeff(Concentration *c, MAC_grid *grid, Parameters *params) {
	
	int i, j, k;
	int NX, NY, NZ;
	double *idx_c, *idy_c, *idz_c;
	double *idx_u, *idy_v, *idz_w;
	double *as, *an, *ap; 
	double *rhsw, *xw; 
	double *modwaveksq, *modwaveisq; 
	double idx, idz;
	int NJ;
	double Pe, iPe;
	double coef, V_s0, b;
	
	NX  = grid->NX;
	NY  = grid->NY;
	NZ  = grid->NZ;

	/*------------------------------------------------------------------------*/
	/*
	     -ap: diagonal
	     -as: south
	     -an: north
	 */
	/*------------------------------------------------------------------------*/
	idy_v = grid->idy_v; 
	idy_c = grid->idy_c; 
	idx_u = grid->idx_u; 
	idz_w = grid->idz_w; 

	NJ = NY-1;
	as = (double *) malloc( NJ*sizeof(double) );
	an = (double *) malloc( NJ*sizeof(double) );
	ap = (double *) malloc( NJ*sizeof(double) );
	c->asw = (double *) malloc( NJ*sizeof(double) );
	c->anw = (double *) malloc( NJ*sizeof(double) );
	c->apw = (double *) malloc( NJ*sizeof(double) );
	c->xw = (double *) malloc( NJ*sizeof(double) );
	c->rhsw = (double *) malloc( NJ*sizeof(double) );

	Pe = c->Pe;
	iPe = 1.0/Pe;
#ifdef LES
	iPe = 1.0;
#endif
	V_s0 = c->v_settl0;
	
	// Coefficient from discretization of d2p/dy2 with Neumann boundary 
	// condition at both walls/free surfaces
	as[0] = 0.0;
	an[NY-2] = 0.0; 
	for (j=0;j<NY-1;j++) {
		if (j > 0 )   as[j] = iPe*idy_v[j]*idy_c[j-1]; 
		an[j] = iPe*idy_v[j]*idy_c[j];
		ap[j] = -as[j] - an[j]; 
	}
	
#ifndef LES
	if ((int)V_s0) {
//		b    = 2.0 * iPe * idy_v[NY-2];
//		coef = an[NY-2] * (b - fabs(V_s0)) / (b + fabs(V_s0) );
		
		// Note that c[NY-1] >> c[NY-2]
		b = fabs(V_s0 )* 0.5 * Pe /idy_c[NY-1];
		coef = an[NY-2] * (1 + b) / (1 - b );
		coef = an[NY-2];
	}
	else {
		coef = an[NY-2];
	} 
	ap[NY-2] = ap[NY-2] + coef;
	an[NY-2] = 0.0; 
#endif

#ifdef THERMAL_KADER
	j = 0;
	as[j] = iPe*idy_v[j] * idy_v[0];
	ap[0] = -an[0] - 2.*as[0]; 
	as[0] = 0.0; 

	j = NY-2;
	an[NY-2] = iPe*idy_v[j]*idy_c[j];
	ap[NY-2] = -as[NY-2] - 2.*an[NY-2];
	an[NY-2] = 0.0;
#endif

	c->as = as;
	c->ap = ap;
	c->an = an; 

	return; 
	
}




/******************************************************************************/
/*
 Create FFTW plans, transpose plans and calculate the coefficient of the matrix
 for the 1-D Helmholtz equation.
 */
/******************************************************************************/
void Conc_setup_lsys_accounting_geometry(Concentration *c, MAC_grid *grid,
		Parameters *params) {
	
	int maxsize;
	
	Conc_calc_coeff(c, grid, params);
	
}




/******************************************************************************/
/*
 This function solves the linear system to find velocity
 */
/******************************************************************************/
int Conc_solve_semi_implicit(int iconc, Cart3d_bag *data_bag) {
	
	int iters;
	double rnorm;
	int count, index, ind1, ind2, ind3, ind4, ind5, ind6;
	int i, j, k; 
	double as_top, an_top, as_bot, an_bot; 
	lsolver_transpose *lsolver_remap;
	double coef, b;
	
	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	Concentration *c = data_bag -> c[iconc];
	
	const double BET[] = {BETA};
	double dtimeb = BET[params -> which_stage] * params -> dt;
	
	double Pe = c -> Pe;
	double iPe = 1.0 / Pe;
	double V_s0 = c -> v_settl0;
	
	lsolver_remap = params->lsolver_remap;
	
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
	
	int i_end = min(Ie, NX-1);
	int j_end = min(Je, NY-1);
	int k_end = min(Ke, NZ-1);
	
	double ***rhs = c->ng_rhs;
#ifndef CONC_BQUICK
	double ***c_data = c->data;
#else
	double ***c_data = c->data_temp;
#endif

#ifdef LES
	double nutN, nutS;
	double ***nut = data_bag->smag->cdev[iconc]->mSct;
#endif  
	double *data_t = lsolver_remap -> data_transpose;
	int Ijs = lsolver_remap -> Ijs;
	int Ije = lsolver_remap -> Ije;
	int Ije_end = min(Ije, NX-1);




	// Copy the RHS of poisson equation to 1-D array data_t
	if (params->NPY != 1) {
		count = 0;
		for (k=Ks;k<Ke;k++){
			for (j=Js;j<Je;j++){
				for (i=Is;i<Ie;i++){
					data_t[count] = rhs[k][j][i];
					count = count + 1 ;
				}
			}
		}
		
		// Transpose the data so that tridiagonal equation can be solved along 
		// y-direction
		remap_3d(data_t, data_t, NULL, lsolver_remap->remap_plan_origtoyline);
	}

	// Tridiagonal solve
#ifdef LES 
	as_bot = grid->idy_v[0] * grid->idy_v[0];
	an_bot = grid->idy_v[0] * grid->idy_c[0];
	
	j = NY-2;
	as_top = grid->idy_v[j] * grid->idy_c[j-1];
	an_top = grid->idy_v[j] * grid->idy_c[j];
	if (c->Type == PARTICLE) {
		b = fabs(V_s0 )* 0.5 * Pe /grid->idy_c[NY-1];
		coef = an_top * (1 + b) / (1 - b );
		coef = an_top;
	}
	else {
		coef = an_top;
	}
#endif
//	printf("Is= %d Ijs= %d Ie= %d Ije= %d rank= %d\n",Is,Ijs,Ie,Ije,params->rank); 

	ind1 = (Ije-Ijs)*(NY);
	ind2 = (Ije-Ijs);
	for (k=Ks;k<k_end;k++){
		ind3 = ind1*(k-Ks); 
		for (i=Ijs;i<Ije_end;i++){
			ind4 = ind3 + (i-Ijs);
			for (j=0;j<NY-1;j++){
#ifndef LES
				c->asw[j] = -c->as[j]; 
				c->anw[j] = -c->an[j]; 
				c->apw[j] = 1.0/dtimeb-c->ap[j]; 
#else // LES
				nutN = 0.5*(nut[k][j+1][i] + nut[k][j][i]); 
				if (j !=0 ) {
					nutS = 0.5*(nut[k][j-1][i] + nut[k][j][i]);
				}
				else {
	#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
					nutS = 0.0;
	#else
					nutS = nut[k][j][i];
	#endif
				}
				nutN = nutN + iPe;
				nutS = nutS + iPe;
				c->asw[j] = -c->as[j]*nutS; 
				c->anw[j] = -c->an[j]*nutN; 
				c->apw[j] = 1.0/dtimeb-c->asw[j]-c->anw[j]; 
#endif // not LES
			}
			if (params->NPY != 1) {
				for (j=0;j<NY-1;j++){
					c->rhsw[j] = data_t[ind4 + j*ind2]; 
				}
			}
			else {
				for (j=0;j<NY-1;j++){
					c->rhsw[j] = rhs[k][j][i]; 
				}
			}
#ifdef LES
			j = 0;
			nutN = 0.5*(nut[k][j+1][i] + nut[k][j][i]); 
			nutN = nutN + iPe;
			c->anw[0] = -an_bot*nutN;
			c->apw[0] = 1.0/dtimeb + an_bot*nutN;
	#ifdef THERMAL_KADER
			c->apw[0] = c->apw[0] + 2.*as_bot*iPe;
	#endif
			c->asw[0] = 0.0;

			j = NY-2;
			nutS = 0.5*(nut[k][j-1][i] + nut[k][j][i]);
			nutS = nutS + iPe;
			c->asw[NY-2] = -as_top*nutS;
			c->apw[NY-2] = 1.0/dtimeb + as_top*nutS;
	#ifdef THERMAL_KADER
			c->apw[NY-2] = c->apw[NY-2] + 2.*an_top*iPe;
	#else
			c->apw[NY-2] = c->apw[NY-2] + iPe*(an_top-coef); 
	#endif
			c->anw[NY-2] = 0.0;
#endif // LES
			
			cltridiag(c->asw, c->apw, c->anw, c->rhsw, c->xw, NY-1);
			
			if (params->NPY != 1) {
				for (j=0;j<NY-1;j++){
					data_t[ind4 + j*ind2]=c->xw[j]; 
				}
				data_t[ind4 + (NY-1)*ind2]=c->xw[NY-2]; 
			}
			else {
				for (j=0;j<NY-1;j++){
					c_data[k][j][i] = c->xw[j];
				}
				c_data[k][NY-1][i] = c->xw[NY-2];
			}
		}
	}

	if (params->NPY != 1) {
		// Transpose the data to original layout
		remap_3d(data_t, data_t, NULL, lsolver_remap->remap_plan_ylinetoorig);
		
		// Copy the solution to p_data array
		count = 0;
		for (k=Ks;k<k_end;k++){
			for (j=Js;j<Je;j++){
				for (i=Is;i<i_end;i++){
					c_data[k][j][i] = data_t[count];
					count = count + 1 ;
				}
				if (Ie != i_end) count = count + 1;
			}
		}
	}
	
	iters = 0;
	return iters;
}

