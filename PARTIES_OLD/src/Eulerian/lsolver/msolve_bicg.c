
#include "DataTypes.h"
#include "Velocity.h"
#include "Display.h"
#include "Boundary.h"

#include <stdio.h>

#warning BiCG method not fully tested

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
	// printf("\n\nINNER_PRODUCT_IDX \n Ks= %d,  Ke= %d \n Js= %d,  Je= %d \n Is= %d,  Ie= %d\n\n", Ks,Ke, Js,Je, Is,Ie); //DEBUG
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {
				partialSum += vec1[k][j][i] * vec2[k][j][i];
			}
		}
	}

	/*
	 // Debugging
	 if (component == 'v') {
		fflush(stdout);
		for (i = 0; i < params -> size; i++) {
	 if (params -> rank == i) {
	 printf("\tProcessor %d: rr = %g\n", params -> rank, partialSum);
	 }
	 fflush(stdout);
	 MPI_Barrier(PCW);
		}
	 }
	 */

	// Communicate and sum inner product among all processors
	MPI_Allreduce(&partialSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, PCW);

	return totalSum;
}




/******************************************************************************/
/*
 Completes matrix-vector multiplication between 'A' matrix and input vector 'x'
 */
/******************************************************************************/
void matVec(double ***Ax, double ***x, char which_component, Cart3d_bag *data_bag) {

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

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
	if (which_component == 'u') {
		Is = max(Is, 1);
#ifdef XPERIODIC
		Ie = grid->G_Ie;
#endif
	}
	if (which_component == 'v') {
		Js = max(Js, 1);
#ifdef YPERIODIC
		Je = grid->G_Je;
#endif
	}
	if (which_component == 'w') {
		Ks = max(Ks, 1);
#ifdef ZPERIODIC
		Ke = grid->G_Ke;
#endif
	}

	double iRe = 1.0 / params -> Re;
	const double BET[] = {BETA};
	double idtimeb = 1.0 / (BET[params -> which_stage] * params -> dt);

	// double ***nutN, ***nutS, ***nutE, ***nutW, ***nutF, ***nutB; // TODO: DELETE?

	double *idx2_e, *idx2_w;
	double *idy2_n, *idy2_s;
	double *idz2_f, *idz2_b;

	idx2_e = grid -> idx2_e;
	idx2_w = grid -> idx2_w;
	idy2_n = grid -> idy2_n;
	idy2_s = grid -> idy2_s;
	idz2_f = grid -> idz2_f;
	idz2_b = grid -> idz2_b;


#ifndef VAR_VISC

	// VISCOSITY
	double nut = iRe;

	if (which_component == 'u') {


		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

//					Ax[k][j][i] = (-1) * x[k][j][i] * (nutW[k][j][i] * idx2_e[i-1] +
//													   nutE[k][j][i] * idx2_w[i] +
//													   nutN[k][j][i] * idy2_n[j] +
//													   nutS[k][j][i] * idy2_s[j] +
//													   nutF[k][j][i] * idz2_f[k] +
//													   nutB[k][j][i] * idz2_b[k]) +
//					x[k][j][i-1] * nutW[k][j][i] * idx2_e[i-1] +
//					x[k][j][i+1] * nutE[k][j][i] * idx2_w[i] +
//					x[k][j-1][i] * nutS[k][j][i] * idy2_s[j] +
//					x[k][j+1][i] * nutN[k][j][i] * idy2_n[j] +
//					x[k-1][j][i] * nutB[k][j][i] * idz2_b[k] +
//					x[k+1][j][i] * nutF[k][j][i] * idz2_f[k] ;

					Ax[k][j][i] = x[k][j][i] * ( nut * ( idx2_e[i-1] +
													     idx2_w[i] +
													     idy2_n[j] +
													     idy2_s[j] +
													     idz2_f[k] +
													     idz2_b[k]
													   ) +
												 idtimeb
											   ) -
					nut * ( x[k][j][i-1] * idx2_e[i-1] +
					        x[k][j][i+1] * idx2_w[i] +
					        x[k][j-1][i] * idy2_s[j] +
					        x[k][j+1][i] * idy2_n[j] +
					        x[k-1][j][i] * idz2_b[k] +
					        x[k+1][j][i] * idz2_f[k]
						  );

//					if (idx2_e[i-1] > idx2_w[i], idy2_n[j], idy2_s[j], idz2_f[k], idz2_b[k]) {
//						<#statements#>
//					}
//
//					printf("[%d][%d][%d]\n",k,j,i);
//					printf("X_e(i-1)= %.5g   X_w= %.5g   Y_n= %.5g   Y_s= %.5g   Z_b= %.5g   Z_f= %.5g\n", idx2_e[i-1], idx2_w[i], idy2_n[j], idy2_s[j], idz2_f[k], idz2_b[k]);

//					if (nut != iRe) {
//						printf("\n\n\n******************** NUT different from iRe ***********************");
//						printf("nut= %.4g     iRe= %.4g\n", nut, iRe);
//					}

				}
			}
		}
	}	// IF

	else if (which_component == 'v')
	{


		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

//					Ax[k][j][i] = (-1) * x[k][j][i] * (nutE[k][j][i] * idx2_e[i] +
//													   nutW[k][j][i] * idx2_w[i] +
//													   nutN[k][j][i] * idy2_s[j] +
//													   nutS[k][j][i] * idy2_n[j-1] +
//													   nutF[k][j][i] * idz2_f[k] +
//													   nutB[k][j][i] * idz2_b[k]) +
//					x[k][j][i-1] * nutW[k][j][i] * idx2_w[i] +
//					x[k][j][i+1] * nutE[k][j][i] * idx2_e[i] +
//					x[k][j-1][i] * nutS[k][j][i] * idy2_n[j-1] +
//					x[k][j+1][i] * nutN[k][j][i] * idy2_s[j] +
//					x[k-1][j][i] * nutB[k][j][i] * idz2_b[k] +
//					x[k+1][j][i] * nutF[k][j][i] * idz2_f[k] ;

					Ax[k][j][i] = x[k][j][i] * nut * (idx2_e[i] +
													  idx2_w[i] +
													  idy2_s[j] +
													  idy2_n[j-1] +
													  idz2_f[k] +
													  idz2_b[k]
													 ) +
					x[k][j][i] * idtimeb -
					nut * ( x[k][j][i-1] * idx2_w[i] +
					        x[k][j][i+1] * idx2_e[i] +
					        x[k][j-1][i] * idy2_n[j-1] +
					        x[k][j+1][i] * idy2_s[j] +
					x[k-1][j][i] * idz2_b[k] +
					x[k+1][j][i] * idz2_f[k]);

				}
			}
		}
	} // IF

	else
	{


		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

//					Ax[k][j][i] = (-1) * x[k][j][i] * (nutE[k][j][i] * idx2_e[i] +
//													   nutW[k][j][i] * idx2_w[i] +
//													   nutN[k][j][i] * idy2_n[j] +
//													   nutS[k][j][i] * idy2_s[j] +
//													   nutF[k][j][i] * idz2_b[k] +
//													   nutB[k][j][i] * idz2_f[k-1]) +
//					x[k][j][i-1] * nutW[k][j][i] * idx2_w[i] +
//					x[k][j][i+1] * nutE[k][j][i] * idx2_e[i] +
//					x[k][j-1][i] * nutS[k][j][i] * idy2_s[j] +
//					x[k][j+1][i] * nutN[k][j][i] * idy2_n[j] +
//					x[k-1][j][i] * nutB[k][j][i] * idz2_f[k-1] +
//					x[k+1][j][i] * nutF[k][j][i] * idz2_b[k] ;

					Ax[k][j][i] = x[k][j][i] * nut * (idx2_e[i] +
													  idx2_w[i] +
													  idy2_n[j] +
													  idy2_s[j] +
													  idz2_b[k] +
													  idz2_f[k-1]
													 ) +
					x[k][j][i] * idtimeb -
					nut * (x[k][j][i-1] * idx2_w[i] +
						   x[k][j][i+1] * idx2_e[i] +
						   x[k][j-1][i] * idy2_s[j] +
						   x[k][j+1][i] * idy2_n[j] +
						   x[k-1][j][i] * idz2_f[k-1] +
						   x[k+1][j][i] * idz2_b[k]
						  );

				}
			}
		}
	} // IF

#endif  // not VAR_VISC

#ifdef VAR_VISC

	// VISCOSITY  //TODO: MOVE TO SEPARATE STRUCTURE
	double ***nut, ***nutX, ***nutY, ***nutZ;
	nut = data_bag -> viscosity -> nu;
	nutX = data_bag -> viscosity -> nuX;
	nutY = data_bag -> viscosity -> nuY;
	nutZ = data_bag -> viscosity -> nuZ;

	double AxConstVisc;  //DEBUGGING //TODO: DELETE

	if (which_component == 'u') {


		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

//					Ax[k][j][i] = (-1) * x[k][j][i] * (nutW[k][j][i] * idx2_e[i-1] +
//													   nutE[k][j][i] * idx2_w[i] +
//													   nutN[k][j][i] * idy2_n[j] +
//													   nutS[k][j][i] * idy2_s[j] +
//													   nutF[k][j][i] * idz2_f[k] +
//													   nutB[k][j][i] * idz2_b[k]) +
//					x[k][j][i-1] * nutW[k][j][i] * idx2_e[i-1] +
//					x[k][j][i+1] * nutE[k][j][i] * idx2_w[i] +
//					x[k][j-1][i] * nutS[k][j][i] * idy2_s[j] +
//					x[k][j+1][i] * nutN[k][j][i] * idy2_n[j] +
//					x[k-1][j][i] * nutB[k][j][i] * idz2_b[k] +
//					x[k+1][j][i] * nutF[k][j][i] * idz2_f[k] ;

					Ax[k][j][i] = x[k][j][i] * (nut[k][j][i-1] * idx2_e[i-1] +
												nut[k][j][i] * idx2_w[i] +
												nutZ[k][j+1][i] * idy2_n[j] +
												nutZ[k][j][i] * idy2_s[j] +
												nutY[k+1][j][i] * idz2_f[k] +
												nutY[k][j][i] * idz2_b[k] +
												idtimeb
											   ) -
					x[k][j][i-1] * nut[k][j][i-1] * idx2_e[i-1] -
					x[k][j][i+1] * nut[k][j][i] * idx2_w[i] -
					x[k][j-1][i] * nutZ[k][j][i] * idy2_s[j] -
					x[k][j+1][i] * nutZ[k][j+1][i] * idy2_n[j] -
					x[k-1][j][i] * nutY[k][j][i] * idz2_b[k] -
					x[k+1][j][i] * nutY[k+1][j][i] * idz2_f[k] ;

				}
			}
		}
	}	// IF

	else if (which_component == 'v')
	{


		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

//					Ax[k][j][i] = (-1) * x[k][j][i] * (nutE[k][j][i] * idx2_e[i] +
//													   nutW[k][j][i] * idx2_w[i] +
//													   nutN[k][j][i] * idy2_s[j] +
//													   nutS[k][j][i] * idy2_n[j-1] +
//													   nutF[k][j][i] * idz2_f[k] +
//													   nutB[k][j][i] * idz2_b[k]) +
//					x[k][j][i-1] * nutW[k][j][i] * idx2_w[i] +
//					x[k][j][i+1] * nutE[k][j][i] * idx2_e[i] +
//					x[k][j-1][i] * nutS[k][j][i] * idy2_n[j-1] +
//					x[k][j+1][i] * nutN[k][j][i] * idy2_s[j] +
//					x[k-1][j][i] * nutB[k][j][i] * idz2_b[k] +
//					x[k+1][j][i] * nutF[k][j][i] * idz2_f[k] ;

					Ax[k][j][i] = x[k][j][i] * (nutZ[k][j][i+1] * idx2_e[i] +
												nutZ[k][j][i] * idx2_w[i] +
												nut[k][j][i] * idy2_s[j] +
												nut[k][j-1][i] * idy2_n[j-1] +
												nutX[k+1][j][i] * idz2_f[k] +
												nutX[k][j][i] * idz2_b[k] +
												idtimeb
											   ) -
					x[k][j][i-1] * nutZ[k][j][i] * idx2_w[i] -
					x[k][j][i+1] * nutZ[k][j][i+1] * idx2_e[i] -
					x[k][j-1][i] * nut[k][j-1][i] * idy2_n[j-1] -
					x[k][j+1][i] * nut[k][j][i] * idy2_s[j] -
					x[k-1][j][i] * nutX[k][j][i] * idz2_b[k] -
					x[k+1][j][i] * nutX[k+1][j][i] * idz2_f[k] ;

				}
			}
		}
	} // IF

	else
	{


		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

//					Ax[k][j][i] = (-1) * x[k][j][i] * (nutE[k][j][i] * idx2_e[i] +
//													   nutW[k][j][i] * idx2_w[i] +
//													   nutN[k][j][i] * idy2_n[j] +
//													   nutS[k][j][i] * idy2_s[j] +
//													   nutF[k][j][i] * idz2_b[k] +
//													   nutB[k][j][i] * idz2_f[k-1]) +
//					x[k][j][i-1] * nutW[k][j][i] * idx2_w[i] +
//					x[k][j][i+1] * nutE[k][j][i] * idx2_e[i] +
//					x[k][j-1][i] * nutS[k][j][i] * idy2_s[j] +
//					x[k][j+1][i] * nutN[k][j][i] * idy2_n[j] +
//					x[k-1][j][i] * nutB[k][j][i] * idz2_f[k-1] +
//					x[k+1][j][i] * nutF[k][j][i] * idz2_b[k] ;

					Ax[k][j][i] = x[k][j][i] * (nutY[k][j][i+1] * idx2_e[i] +
												nutY[k][j][i] * idx2_w[i] +
												nutX[k][j+1][i] * idy2_n[j] +
												nutX[k][j][i] * idy2_s[j] +
												nut[k][j][i] * idz2_b[k] +
												nut[k-1][j][i] * idz2_f[k-1] +
												idtimeb) -
					x[k][j][i-1] * nutY[k][j][i] * idx2_w[i] -
					x[k][j][i+1] * nutY[k][j][i+1] * idx2_e[i] -
					x[k][j-1][i] * nutX[k][j][i] * idy2_s[j] -
					x[k][j+1][i] * nutX[k][j+1][i] * idy2_n[j] -
					x[k-1][j][i] * nut[k-1][j][i] * idz2_f[k-1] -
					x[k+1][j][i] * nut[k][j][i] * idz2_b[k] ;

				}
			}
		}
	} // IF

#endif  // VAR_VISC

}




/******************************************************************************/
/*
 Solve the system
	Ax = b
 using the biconjugate-gradient method

 http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
 */
/******************************************************************************/
int Velocity_solve_bicg(Velocity *vel, Cart3d_bag *data_bag) {

	int iters;
	int i, j, k;
	char statement[100];

	double T1,T2;

	double EPS = 1e-14; // Small number to prevent division by zero

	MAC_grid *grid = data_bag -> grid;
	Parameters *params = data_bag -> params;

	double ***data = vel -> data;
	double ***rhs = vel -> ng_rhs;
	char component = vel -> component;

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


	//--------------------------------------------------------------------------
	// Conjugate-gradient variables
	//--------------------------------------------------------------------------

	// Residual
	double ***r = vel -> ng_r;
	double ***p = vel -> p;
	double ***s = vel -> s;
	double ***r0 = vel -> ng_r0;

	// A * p
	// A * s
	double ***Ap = vel -> Ap;
	double ***As = vel -> ng_As;

	// Inner products
	double rr0, rr, r0Ap, AsAs, sAs, RMS, in_tot;

	// Global  coefficients
	double alpha, beta, omega, rho, rho_old;

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


	// see WIKI http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method#Algorithmic_steps

	// >->->->->->->->-> 1., 2., 4. <-<-<-<-<-<-<-<-<-<

	matVec(As, data, component, data_bag);

	T1 = MPI_Wtime();
	for (k = Ks; k < Ke; k++) {			// _g ???
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {

				r[k][j][i] = rhs[k][j][i] - As[k][j][i];
				r0[k][j][i] = r[k][j][i];
				p[k][j][i] = 0;
				Ap[k][j][i] = 0;
			}
		}
	}

	// rr0 used for convergence check
	rr0 = innerProd(r, r, component, grid, params);

	// >->->->->->->->->->-> 3. <-<-<-<-<-<-<-<-<-<-<-<

	rho = alpha = omega = 1;

	// >->->->->->->->->->-> 5. <-<-<-<-<-<-<-<-<-<-<-<

	iters = 0;
	while (iters < params -> CG_MAXIT) {

		iters++;
		//		printf("\n\n ITERATION -->>> %d <<<--\n", iters);  //DEBUG

		// >->->->->->->->->->-> 5.1. <-<-<-<-<-<-<-<-<-<-<-<

		rho_old = rho;
		rho = innerProd(r0,r, component, grid, params);

		// >->->->->->->->->->-> 5.2. <-<-<-<-<-<-<-<-<-<-<-<

		beta = (rho * alpha) / (rho_old * omega + EPS);		// + EPS ???

		// >->->->->->->->->->-> 5.3. <-<-<-<-<-<-<-<-<-<-<-<

		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

					p[k][j][i] = r[k][j][i] + beta * (p[k][j][i] - omega * Ap[k][j][i]);
				}
			}
		}

		// >->->->->->->->->->-> 5.4. <-<-<-<-<-<-<-<-<-<-<-<

		// Update boundary cells for 'p'
		Velocity_update_boundaries(p, component, VEL_TYPE_CG, data_bag);
		// Ap = A * p
		matVec(Ap, p, component, data_bag);

		// >->->->->->->->->->-> 5.5. <-<-<-<-<-<-<-<-<-<-<-<

		r0Ap = innerProd(r0, Ap, component, grid, params);
		alpha = rho / (r0Ap + EPS);

		// >->->->->->->->->->-> 5.6. <-<-<-<-<-<-<-<-<-<-<-<

		for (k = Ks; k < Ke; k++) {
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {

					s[k][j][i] = r[k][j][i] - alpha * Ap[k][j][i];
				}
			}
		}

		// >->->->->->->->->->-> 5.7. <-<-<-<-<-<-<-<-<-<-<-<

		// Update boundary cells for 's'
		Velocity_update_boundaries(s, component, VEL_TYPE_CG, data_bag);
		// As = A * s
		matVec(As, s, component, data_bag);

		// >->->->->->->->->->-> 5.8. <-<-<-<-<-<-<-<-<-<-<-<

		sAs = innerProd(As, s, component, grid, params);
		AsAs = innerProd(As, As, component, grid, params);
		omega = sAs / (AsAs + EPS);

		// >->->->->->->->->->-> 5.9. <-<-<-<-<-<-<-<-<-<-<-<

		// x = x + alpha * Ap;
		for (k = Ks; k < Ke; k++) {				// _g ???
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
					data[k][j][i] = data[k][j][i] + alpha * p[k][j][i] + omega * s[k][j][i];
				}
			}
		}

		// >->->->->->->->->->-> 5.10. <-<-<-<-<-<-<-<-<-<-<-<

		rr = innerProd(r,r, component, grid, params);
		RMS = sqrt(in_tot * rr);
		if ( RMS < params -> CG_ETOL ) // TODO: ACTIVATE Conv check
		{
			// printf("BiCGStab converged on ITERATION >>>> %d <<<<", iters);
			break;
		}

		// >->->->->->->->->->-> 5.11. <-<-<-<-<-<-<-<-<-<-<-<

		// r = s - omega * As;
		for (k = Ks; k < Ke; k++) {				// _g ???
			for (j = Js; j < Je; j++) {
				for (i = Is; i < Ie; i++) {
					r[k][j][i] = s[k][j][i] - omega * As[k][j][i];
				}
			}
		}

		// -------------------------------------------------------------------------------------------------

//				printf("\n\n VELOCITY RESULT AFTER ITERATION >>>%d<<<\n\n", iters);		//TODO: DELETE
//				//printf("U = \n");
//				for (k = Ks; k < Ke; k++) {
//					printf("\n\n    k= %d \n",k);
//					for (j = Js; j < Je; j++) {
//						printf("\n");
//						for (i = Is; i < Ie; i++) {
//							printf("%.6g  ", data[k][j][i]);
//						}
//					}
//				}

	}

	//printf("\n");
	sprintf(statement, "%c-component of velocity converged to %.12g after %d iterations\n",
			component, RMS, iters);
	Display_progress(params, statement);


//	//----------- DEBUG ----------- //TODO: DELETE
//			printf("\n\n Vel = \n\n");
//
//			for (k = Ks; k < Ke; k++) {
//				printf("\n\n k = %d\n", k);
//				for (j = Js; j < Je; j++) {
//					printf("\n");
//					for (i = Is; i < Ie; i++) {
//
//						printf("%.5g  ", data[k][j][i]);
//					}
//				}
//			}
//	//--------------------------------------------



	return iters;
}

//-------------------------------------------------------------
