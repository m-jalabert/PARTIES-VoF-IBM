#include "dsolver.h"

/******************************************************************************/
/*
 Solves the tridiagonal system with coefficients
 
     a, b, c and r for a vector u
 
 where a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1] = r[i]
 
 Assumes a[0] and c[N-1] are zero (i.e does not work for periodic boundary)
 */
/******************************************************************************/
void ltridiag(double *a, double *b, double *c, double *r, double *u, int N){
	
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
 Use Steve Plimpton's (Sandia) functions to tranpose data
 
 Transposing data is needed when we take FFT/FCT or when we solve tri-diagonal
 equation along y-direction. Note that in original MPI decomposition of the
 domain data along x- y- and z-directions are not on the same processor.
 */
/******************************************************************************/
void lsolver_transpose_setup(MAC_grid *grid, Parameters *params) {
	
	struct remap_plan_3d *remap_plan_origtoyline;
	struct remap_plan_3d *remap_plan_ylinetoorig;
	
	int Is, Js, Ks, Ie, Je, Ke;
	int i, j, k;
	int ierr;
	int Ijs, Ije; 
	int NX, NY, NZ; 
	double dnpy; 
	double *data_transpose; 
	int count;
	FILE *fid;
	lsolver_transpose *lsolver_remap; 

	lsolver_remap = (lsolver_transpose *)malloc(sizeof(lsolver_transpose)); 
	params->lsolver_remap = lsolver_remap; 

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;
					
	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	NX  = grid->NX;
	NY  = grid->NY;
	NZ  = grid->NZ;

	dnpy = params->NPY;

	// Starting and ending index along i-direction when all the processors have 
	// all the data (NY) along j-direction
	Ijs = Is + params->yproccoord*(Ie-Is)/dnpy;
	Ije = Is + (params->yproccoord+1)*(Ie-Is)/dnpy;

	// Transpose plan for having all the data in j-direction
	remap_plan_origtoyline = remap_3d_create_plan(
			PCW,
			Is, Ie-1,
			Js, Je-1,
			Ks, Ke-1,
			Ijs, Ije-1,
			0, NY-1,
			Ks, Ke-1,
			1, 0, 1, 2);

	remap_plan_ylinetoorig = remap_3d_create_plan(
			PCW,
			Ijs, Ije-1,
			0, NY-1,
			Ks, Ke-1,
			Is, Ie-1,
			Js, Je-1,
			Ks, Ke-1,
			1, 0, 1, 2);

	count = max((Ie-Is)*(Je-Js)*(Ke-Ks), (Ije-Ijs)*(NY)*(Ke-Ks) );
	data_transpose = (double *) malloc( 2*count*sizeof(double) );
	printf("data_transpose count = %d rank = %d\n",count,params->rank);

	for (i=0;i<params->size;i++){
		if (i==params->rank) {
			if (i==0) {
				fid = fopen("vindex_info.dat","w");
			}
			else {
				fid = fopen("vindex_info.dat","a");
			}
			fprintf(fid,".........Rank=%d........\n",params->rank);
			fprintf(fid,"%4d %4d %4d %4d %4d %4d %9d\n",Is,Ie-1,Js,Je-1,Ks,Ke-1, (Ie-Is)*(Je-Js)*(Ke-Ks));
			fprintf(fid,"%4d %4d %4d %4d %4d %4d %9d\n",Ijs,Ije-1,0,NY-1,Ks,Ke-1, (Ije-Ijs)*(NY)*(Ke-Ks));
			fclose(fid);
		}
		MPI_Barrier(PCW);
	}

	lsolver_remap->remap_plan_origtoyline = remap_plan_origtoyline;
	lsolver_remap->remap_plan_ylinetoorig = remap_plan_ylinetoorig;

	lsolver_remap->Ijs = Ijs; 
	lsolver_remap->Ije = Ije; 
	lsolver_remap->data_transpose = data_transpose;

	return;
}
