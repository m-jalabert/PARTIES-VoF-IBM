

/******************************************************************************/
/*
 Calculates the inner product between two vectors vec1 and vec2
 */
/******************************************************************************/
double innerProd(double ***vec1, double ***vec2, char component, MAC_grid *grid, Parameters *params) {

	int i, j, k;

	double partialSum, totalSum;

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
	if (component == 'u') {
		Is = max(Is, 1);
#ifdef XPERIODIC
		Ie = grid->G_Ie;
#endif
	}
	if (component == 'v') {
		Js = max(Js, 1);
#ifdef YPERIODIC
		Je = grid->G_Je;
#endif
	}
	if (component == 'w') {
		Ks = max(Ks, 1);
#ifdef ZPERIODIC
		Ke = grid->G_Ke;
#endif
	}

	// Calculate inner product of local nodes
	partialSum = 0;
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				partialSum += vec1[k][j][i] * vec2[k][j][i];
			}
		}
	}

	// Communicate and sum inner product among all processors
	MPI_Allreduce(&partialSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, PCW);

	return totalSum;
}




/******************************************************************************/
/*
 Completes matrix-vector multiplication between 'A' matrix and input vector 'x'
 */
/******************************************************************************/
void matVec(double ***Ax, double ***x, char component, MAC_grid *grid, Parameters *params) {

	int i, j, k;

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
	if (component == 'u') {
		Is = max(Is, 1);
#ifdef XPERIODIC
		Ie = grid->G_Ie;
#endif
	}
	if (component == 'v') {
		Js = max(Js, 1);
#ifdef YPERIODIC
		Je = grid->G_Je;
#endif
	}
	if (component == 'w') {
		Ks = max(Ks, 1);
#ifdef ZPERIODIC
		Ke = grid->G_Ke;
#endif
	}

	// Grid spacing (assuming uniform)
	double iddx = grid -> idx_c[1];
	double iddy = grid -> idy_c[1];
	double iddz = grid -> idz_c[1];
	iddx = iddx * iddx;
	iddy = iddy * iddy;
	iddz = iddz * iddz;

	double iRe = 1.0 / params -> Re;
	const double BET[] = {BETA};
	double idtimeb = 1.0 / (BET[params -> which_stage] * params -> dt);

	double ac = idtimeb + 2.0 * iRe * (iddx + iddy + iddz);
	double ax = -iRe * iddx;
	double ay = -iRe * iddy;
	double az = -iRe * iddz;

	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				Ax[k][j][i] = ac * x[k][j][i]
				            + ax * ( x[k][j][i-1] + x[k][j][i+1] )
				            + ay * ( x[k][j-1][i] + x[k][j+1][i] )
				            + az * ( x[k-1][j][i] + x[k+1][j][i] );
			}
		}
	}
}




/******************************************************************************/
/*
 Solve the system
	Ax = b
 using the conjugate-gradient method

 http://en.wikipedia.org/wiki/Conjugate_gradient_method

 http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
 */
/******************************************************************************/
int Velocity_solve_cg(Velocity *vel, Cart3d_bag *data_bag) {

	int iters;
	int i, j, k;
	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;
	char statement[100];
	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;
	double ***data = vel -> data;
	double ***rhs = vel -> ng_rhs;
	char component = vel -> component;
	double EPS = 1e-14; // Small number to prevent division by zero


	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = min(grid->G_Ie, grid->NX-1);
	int Je = min(grid->G_Je, grid->NY-1);
	int Ke = min(grid->G_Ke, grid->NZ-1);

	// We are not solving for the nodes on the boundaries unless periodic
	if (component == 'u') {
		Is = max(Is, 1);
#ifdef XPERIODIC
		Ie = grid->G_Ie;
#endif
	}
	if (component == 'v') {
		Js = max(Js, 1);
#ifdef YPERIODIC
		Je = grid->G_Je;
#endif
	}
	if (component == 'w') {
		Ks = max(Ks, 1);
#ifdef ZPERIODIC
		Ke = grid->G_Ke;
#endif
	}

	//--------------------------------------------------------------------------
	// Conjugate-gradient variables
	//--------------------------------------------------------------------------

	// Residual
	double ***r = vel -> ng_r;

	double ***d = vel -> d;

	// A * d
	double ***Ad = vel -> ng_Ad;

	// Inner products
	double rr, rr_old, rr0, dAd, RMS, in_tot;

	// Global  coefficients
	double alpha, beta;

	if (vel -> component == 'u') {
#ifdef XPERIODIC
		in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-1));
#else
		in_tot = 1.0 / (double)((NX-2)*(NY-1)*(NZ-1));
#endif
	}
	if (vel -> component == 'v') {
#ifdef YPERIODIC
		in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-1));
#else
		in_tot = 1.0 / (double)((NX-1)*(NY-2)*(NZ-1));
#endif
	}
	if (vel -> component == 'w') {
#ifdef ZPERIODIC
		in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-1));
#else
		in_tot = 1.0 / (double)((NX-1)*(NY-1)*(NZ-2));
#endif
	}


	/*------------------------------------------------------------------------*/
	/*
	 Initialize 'r' and 'd' variables
	 */
	/*------------------------------------------------------------------------*/



	matVec(Ad, data, component, grid, params);

	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				r[k][j][i] = rhs[k][j][i] - Ad[k][j][i];
				d[k][j][i] = r[k][j][i];
			}
		}
	}


	/*------------------------------------------------------------------------*/
	/*
	 Conjugate-gradient solve
	 */
	/*------------------------------------------------------------------------*/
/*

	if(component == 'v'){
		printf("My name is r [0][NY-3][0]=%2.3f\n",r[0][NY-3][0]);
		printf("My name is r [0][NY-2][0]=%2.3f\n",r[0][NY-2][0]);
		printf("My name is r [0][NY-1][0]=%2.3f\n",r[0][NY-1][0]);
		printf("My name is r [0][1][0]=%2.3f\n",	r[0][1][0]);
		printf("My name is r [0][2][0]=%2.3f\n",	r[0][2][0]);
		printf("My name is r [0][3][0]=%2.3f\n",	r[0][3][0]);

	}*/
	rr = innerProd(r, r, component, grid, params);

	int yproc = params-> yproccoord;
	int NPY = params->NPY;

int imax = 0;
int jmax =0;
int kmax =0;
double rmax =0.0;


	rr0 = rr + EPS;
	iters = 0;

	while (iters < params -> CG_MAXIT) {

		iters++;

		// Update boundary cells for 'd'
		Velocity_update_boundaries(d, component, VEL_TYPE_CG, data_bag);

/*
		if(component == 'v'){
			printf("My name is d [0][NY-3][0]=%2.3f\n",	d[0][NY-3][0]);
			printf("My name is d [0][NY-2][0]=%2.3f\n",	d[0][NY-2][0]);
			printf("My name is d [0][NY-1][0]=%2.3f\n",	d[0][NY-1][0]);
			printf("My name is d [0][0][0]=%2.3f\n",	d[0][0][0]);
			printf("My name is d [0][1][0]=%2.3f\n",	d[0][1][0]);
			printf("My name is d [0][2][0]=%2.3f\n",	d[0][2][0]);

		}
*/
		// Ad = A * d
		matVec(Ad, d, component, grid, params);

		/* if(component == 'v'){
			if(component == 'v'){
				printf("My name is Ad [0][NY-3][0]=%2.3f\n",	Ad[0][NY-3][0]);
				printf("My name is Ad [0][NY-2][0]=%2.3f\n",	Ad[0][NY-2][0]);
				printf("My name is Ad [0][NY-1][0]=%2.3f\n",	Ad[0][NY-1][0]);
				printf("My name is Ad [0][0][0]=%2.3f\n",	Ad[0][0][0]);
				printf("My name is Ad [0][1][0]=%2.3f\n",	Ad[0][1][0]);
				printf("My name is Ad [0][2][0]=%2.3f\n",	Ad[0][2][0]);

			}

		}*/

		dAd = innerProd(d, Ad, component, grid, params);

		alpha = rr / (dAd + EPS);

		// x = x + alpha * d;
		// r = r - alpha * Ad;


		/* if(component == 'v'){
			printf("d[0][NY-1][0]=%2.3f\n",	d[0][NY-1][0]);
			printf("d[0][NY-2][0]=%2.3f\n",	d[0][NY-2][0]);
			printf("d[0][1][0]=%2.3f\n",	d[0][1][0]);

		}*/

		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
					data[k][j][i] += alpha * d[k][j][i];

					r[k][j][i] -= alpha * Ad[k][j][i];
					if (abs(r[k][j][i])>rmax)
					{rmax = abs(r[k][j][i]);
					imax = i;
					jmax = j;
					kmax = k;}
				}
			}
		}
		//printf("imax = %d, jmax=%d, kmax=%d \n",imax,jmax,kmax);

/*
		if(component == 'v'){
			if (j==NY-1){
				if(k==0){
				if(i==0){
			printf("My name is r [0][0][0]=%e at iteration %d\n",r[k][j][i],iters);}}}}
*/
		rr_old = rr;
		rr = innerProd(r, r, component, grid, params);
 		RMS = sqrt(in_tot * rr);
		if (RMS < params -> CG_ETOL)
			break;

		beta = rr / (rr_old + EPS);

		// d = beta * d + r
		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
					d[k][j][i] = beta * d[k][j][i] + r[k][j][i];
				}
			}
		}
	}

	sprintf(statement, "%c-component of velocity converged to %g after %d iterations\n",
			component, RMS, iters);
	Display_progress(params, statement);

	return iters;
}
