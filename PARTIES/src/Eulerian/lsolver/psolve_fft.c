
//#if defined XPERIODIC && defined ZPERIODIC
//	#define Pressure_solve Pressure_solve_xz_period
//#elif defined XPERIODIC && !defined ZPERIODIC
//	#define Pressure_solve Pressure_solve_x_period
//#elif !defined XPERIODIC && defined ZPERIODIC
//	#define Pressure_solve Pressure_solve_z_period
//#else
//	#define Pressure_solve Pressure_solve_no_period
//#endif

/******************************************************************************/
/*
 Creates FFTW plans for transforms and also creates wave-number arrays
 */
/******************************************************************************/
void Pressure_fft_setup(fft *vfft, int N) {

	fftw_plan rplanf, rplanb;
	fftw_plan cplanf, cplanb;
	fftw_plan cosplanf, cosplanb;

	fftw_real *ad;
	fftw_complex *ac;
	int i;
	int *wave, *wave_rehac;

	ad = (fftw_real*)    fftw_malloc( sizeof(fftw_real)    * N );
	ac = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * N );

	/*------------------------------------------------------------------------*/
	/*
	 Create FFTW plans for forward(*f) and backward(*b) real-to-real(*r2r*),
	 complex-to-complex(*dft*) FFT/Fast cosine transforms (FFTW_REDFT10 and
	 FFTW_REDFT01)
	 */
	/*------------------------------------------------------------------------*/
	rplanf   = fftw_plan_r2r_1d(N, ad, ad, FFTW_R2HC, FFTW_PATIENT );
	rplanb   = fftw_plan_r2r_1d(N, ad, ad, FFTW_HC2R, FFTW_PATIENT );
	cplanf   = fftw_plan_dft_1d(N, ac, ac, FFTW_FORWARD, FFTW_PATIENT );
	cplanb   = fftw_plan_dft_1d(N, ac, ac, FFTW_BACKWARD, FFTW_PATIENT );
	cosplanf = fftw_plan_r2r_1d(N, ad, ad, FFTW_REDFT10, FFTW_PATIENT );
	cosplanb = fftw_plan_r2r_1d(N, ad, ad, FFTW_REDFT01, FFTW_PATIENT );

	// Calculate wave number arrays (standard FFTW ordering)
	wave  = (int*) malloc(N*sizeof(int));
	for (i=0;i<N/2;i++) {
		wave[i] = i;
	}
	for (i=N/2;i<N;i++) {
		wave[i] = (N-i);
	}


	// RE-ordered HAlf-Complex(rehac) ordering
	wave_rehac  = (int*) malloc(N*sizeof(int));
	wave_rehac[0] = 0;
	wave_rehac[1] = N/2;
	for (i=2;i<N;i=i+2) {
		wave_rehac[i]   = i/2;
		wave_rehac[i+1] = i/2;
	}

	vfft->N = N;
	vfft->rplanf = rplanf;
	vfft->rplanb = rplanb;
	vfft->cplanf = cplanf;
	vfft->cplanb = cplanb;
	vfft->cosplanf = cosplanf;
	vfft->cosplanb = cosplanb;

	vfft->wave = wave;
	vfft->wave_rehac = wave_rehac;

	return;
}




/******************************************************************************/
/*
 ! fft_module/fftw_halfcomplex_to_rehac
 !
 ! Convert from FFTW's half-complex array to rehac (REordered HAlf Complex)
 ! ordering which is required for parallel processing. The main aim of the
 ! rehac ordering is to keep the real and imaginary parts together in the
 ! array.
 */
/******************************************************************************/
void fftw_halfcomplex_to_rehac(int N, fftw_real *a, fftw_real *b) {

	int i;

	b[0] = a[0];
	b[1] = a[N/2];
	for (i=2;i<N;i=i+2){
		b[i]   = a[i/2];
		b[i+1] = a[N-i/2];
	}
	return;
}




/******************************************************************************/
/*
  ! fft_module/fftw_rehac_to_halfcomplex
  !
  ! Convert from rehac ordering to FFTW's default half-complex ordering.
 */
/******************************************************************************/

void  fftw_rehac_to_halfcomplex(int N, fftw_real *a, fftw_real *b) {

	int i;

	b[0] = a[0];
	b[N/2] = a[1];
	for (i=1; i<N/2; i++) {
		b[i]    = a[2*i];
		b[N-i]  = a[2*i+1];
	}
	return;
}




/******************************************************************************/
/*
 Solves the tridiagonal system with coefficients a, b, c, and r for a vector u
 where

     a[i] * u[i-1] + b[i] * u[i] + c[i] * u[i+1] = r[i]

 Assumes a[0] and c[N-1] are zero (i.e does not work for periodic boundary)
 */
/******************************************************************************/
void tridiag(double *a, double *b, double *c, double *r, double *u, double *gam, int N){

//	double *gam, bet;
	double bet;
	int i;

//	gam = (double *) malloc(N*sizeof(double));

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

//	free(gam);
	return;
}

  void tridiag_yper(double *a, double *b, double *c, double *r, double *u, double *yw, double *qw, double *gam, int N){
//tridiag_yper(p->asw, p->apw, p->anw, p->rhsw, p->u, p->xw, p->qw, p->temp, NY-1);

//	double *gam, bet;
	double bet;
	int i;


//	gam = (double *) malloc(N*sizeof(double));
/*
for (i=1;i<N;i++){
	gam[i] = a[i]/b[i-1];
	b[i] = b[i] - c[i-1]*gam[i];
	r[i] = r[i] - gam[i]*r[i-1];
	u[i] = u[i] - u[i]*u[i-1];
}

yw[N-1] = r[N-1]/b[N-1];
qw[N-1] = u[N-1]/b[N-1];

for (i=N-2;i>=0; i--){
	yw[i] = (r[i]-c[i]*yw[i+1])/b[i];
	qw[i] = (u[i]-c[i]*qw[i+1])/b[i];
}*/

	bet = b[0];
	yw[0] = r[0]/bet;
	qw[0] = u[0]/bet;

	for (i=1;i<N;i++){
		gam[i] = c[i-1]/bet;
		bet = b[i] - a[i]*gam[i];
		yw[i] = (r[i]-a[i]*yw[i-1])/bet;
		qw[i] = (u[i]-a[i]*qw[i-1])/bet;
	}
	for (i=N-2;i>=0; i--){
		yw[i] = yw[i]-gam[i+1]*yw[i+1];
		qw[i] = qw[i]-gam[i+1]*qw[i+1];
	}

//	free(gam);
	return;
}

/******************************************************************************/
/*
 Use Steve Plimpton's (Sandia) functions to tranpose data. Transposing data is
 needed when we take FFT/FCT or when we solve tri-diagonal equation along
 y-direction. Note that in original MPI decomposition of the domain data along
 x- y- and z-directions are not on the same processor.
 */
/******************************************************************************/
void Pressure_transpose_setup(Pressure *p, MAC_grid *grid, Parameters *params) {

	struct remap_plan_3d *remap_plan_origtozline;
	struct remap_plan_3d *remap_plan_zlinetoxline;
	struct remap_plan_3d *remap_plan_xlinetoxyplane;

	struct remap_plan_3d *remap_plan_xyplanetoxline;
	struct remap_plan_3d *remap_plan_xlinetozline;
	struct remap_plan_3d *remap_plan_zlinetoorig;

	int Is, Js, Ks, Ie, Je, Ke;
	int i, j, k;
	int Iks, Ike;
	int Kis, Kie;
	int Kijs, Kije;
	int NX, NY, NZ;
	double dnpx, dnpy, dnpz;
	double *data_transpose;
	int NXM, NZM;
	int PERPROC, REMAINDER, TODISTRI;
	int count;
	FILE *fid;

	NXM = params->NXM;
	if (NXM%2 == 1) {
		printf("NXM should be divisible by 2\n");
		Communication_finalize();
		exit(0);
	}
	NZM = params->NZM;
	if (NZM%2 == 1) {
		printf("NZM should be divisible by 2\n");
		Communication_finalize();
		exit(0);
	}

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

	if (Ie == NX )
		Ie = Ie-1;
	if (Je == NY )
		Je = Je-1;
	if (Ke == NZ )
		Ke = Ke-1;

	dnpz = params->NPZ;
	dnpx = params->NPX;
	dnpy = params->NPY;

	// Starting and ending index along i-direction when all the processors have
	// all the data (NZ-1) along k-direction
	Iks = Is + params->zproccoord*(Ie-Is)/dnpz;
	Ike = Is + (params->zproccoord+1)*(Ie-Is)/dnpz;

//	Kis = Ks + params->xproccoord*(Ke-Ks)/dnpx;
//	Kie = Ks + (params->xproccoord+1)*(Ke-Ks)/dnpx;

	/*------------------------------------------------------------------------*/
	/*
	 Starting and ending index along k-direction when all the processors have
	 all the data (NX-1) along i-direction.

	 Note that data is divided along k-direction so that each processor has even
	 number of data along k-direction.

	 Afer taking fft along k-direction, the data should be complex. However, the
	 data from r2r transform is stored in a real array, where a pair of two
	 numbers (but not necessarily adjacent) in this array form the real and
	 imaginary part of the complex fft transform.
	 */
	/*------------------------------------------------------------------------*/
	i = params->zproccoord*dnpx;
	i = i + params->xproccoord;
	NZM = params->NZM;
	PERPROC = NZM/(dnpx*dnpz);
	REMAINDER = NZM - PERPROC*params->NPX*params->NPZ;
	if (PERPROC%2 != 0) {
		PERPROC -= 1;
		REMAINDER += params->NPX*params->NPZ;
	}
	if (params->rank==0)
		printf("PERPROC = %d REMAINDER = %d\n",PERPROC,REMAINDER);

	TODISTRI = REMAINDER/2;
	i = params->zproccoord*dnpx;
	i = i + params->xproccoord;
	if (TODISTRI == 0) {
		Kis = i*PERPROC;
		Kie = (i+1)*PERPROC;
	}
	else {
		if (i < TODISTRI ) {
			Kis = i*PERPROC + 2*i;
			Kie = (i+1)*PERPROC + 2*(i+1);
		}
		else {
			Kis = i*PERPROC + 2*TODISTRI;
			Kie = (i+1)*PERPROC + 2*TODISTRI;
		}
	}

	/*------------------------------------------------------------------------*/
	/*
	 Starting and ending index along k-direction when all the processors have
	 all the data (NX-1) and (NY-1) along i- and j-direction.
	 */
	/*------------------------------------------------------------------------*/
//	if (dnpy > 1) {
		i = Ks + params->xproccoord*(Ke-Ks)/dnpx;
		j = Ks + (params->xproccoord+1)*(Ke-Ks)/dnpx;
		Kijs = i + params->yproccoord*(j-i)/dnpy;
		Kije = i + (params->yproccoord+1)*(j-i)/dnpy;
/*
		Kis = Ks + params->xproccoord*(Ke-Ks)/dnpx;
		Kie = Ks + (params->xproccoord+1)*(Ke-Ks)/dnpx;
*/
		Kijs = Kis + params->yproccoord*(Kie-Kis)/dnpy;
		Kije = Kis + (params->yproccoord+1)*(Kie-Kis)/dnpy;
//	}

	// Transpose plan for having all the data in k-direction
	remap_plan_origtozline = remap_3d_create_plan(
			PCW,
			Is, Ie-1,
			Js, Je-1,
			Ks, Ke-1,
			Iks, Ike-1,
			Js, Je-1,
			0, NZ-2,
			1, 0, 1, 2);

	// Transpose plan for having all the data in i-direction
	remap_plan_zlinetoxline = remap_3d_create_plan(
			PCW,
			Iks, Ike-1,
			Js, Je-1,
			0, NZ-2,
			0, NX-2,
			Js, Je-1,
			Kis, Kie-1,
			1, 0, 1, 2);

//	if (dnpy > 1) {
		// Transpose plan for having all the data in i- and j-direction
		remap_plan_xlinetoxyplane = remap_3d_create_plan(
				PCW,
				0, NX-2,
				Js, Je-1,
				Kis, Kie-1,
				0, NX-2,
				0, NY-2,
				Kijs, Kije-1,
				1, 0, 1, 2);
		remap_plan_xyplanetoxline = remap_3d_create_plan(
				PCW,
				0, NX-2,
				0, NY-2,
				Kijs, Kije-1,
				0, NX-2,
				Js, Je-1,
				Kis, Kie-1,
				1, 0, 1, 2);
//	}

	remap_plan_xlinetozline = remap_3d_create_plan(
			PCW,
			0, NX-2,
			Js, Je-1,
			Kis, Kie-1,
			Iks, Ike-1,
			Js, Je-1,
			0, NZ-2,
			1, 0, 1, 2);

	remap_plan_zlinetoorig = remap_3d_create_plan(
			PCW,
			Iks, Ike-1,
			Js, Je-1,
			0, NZ-2,
			Is, Ie-1,
			Js, Je-1,
			Ks, Ke-1,
			1, 0, 1, 2);

	count = max((Ie-Is)*(Je-Js)*(Ke-Ks), max((Ike-Iks)*(Je-Js)*(NZ-1),
			max((Kie-Kis)*(Je-Js)*(NX-1), (Kije-Kijs)*(NY-1)*(NX-1)) ) );
	data_transpose = (double *) malloc( count*sizeof(double) );

	for (i=0;i<params->size;i++){
		if (i==params->rank) {
			if (i==0) {
				fid = fopen("index_info.dat","w");
			}
			else {
				fid = fopen("index_info.dat","a");
			}
			fprintf(fid,".........Rank=%d........\n",params->rank);
			fprintf(fid,"%4d %4d %4d %4d %4d %4d %9d\n",Is,Ie-1,Js,Je-1,Ks,Ke-1, (Ie-Is)*(Je-Js)*(Ke-Ks));
			fprintf(fid,"%4d %4d %4d %4d %4d %4d %9d\n",Iks,Ike-1,Js,Je-1,0,NZ-2,(Ike-Iks)*(Je-Js)*(NZ-1));
			fprintf(fid,"%4d %4d %4d %4d %4d %4d %9d\n",0,NX-2,Js,Je-1,Kis,Kie-1,(NX-1)*(Je-Js)*(Kie-Kis));
			fprintf(fid,"%4d %4d %4d %4d %4d %4d %9d\n",0,NX-2,0,NY-2,Kijs,Kije-1,(NX-1)*(NY-1)*(Kije-Kijs));
			fclose(fid);
		}
		MPI_Barrier(PCW);
	}

	p->remap->remap_plan_origtozline     = remap_plan_origtozline;
	p->remap->remap_plan_zlinetoxline    = remap_plan_zlinetoxline;
	p->remap->remap_plan_xlinetoxyplane  = remap_plan_xlinetoxyplane;

	p->remap->remap_plan_xyplanetoxline  = remap_plan_xyplanetoxline;
	p->remap->remap_plan_xlinetozline    = remap_plan_xlinetozline;
	p->remap->remap_plan_zlinetoorig     = remap_plan_zlinetoorig;

	p->remap->Iks = Iks;
	p->remap->Ike = Ike;
	p->remap->Kis = Kis;
	p->remap->Kie = Kie;
	p->remap->Kijs = Kijs;
	p->remap->Kije = Kije;
	p->remap->data_transpose = data_transpose;

	return;
}




/******************************************************************************/
/*
 Calculate the coefficients of tri-diagonal matrix resulting from dicretization
 of 1-D Helomholtz equation in y-direction.

 Note that taking FFT/FCT along x- and z-direction of the original Poisson
 equation gives us series of 1-D Helmholtz equations for each wave-number
 (k_x and k_z).
 */
/******************************************************************************/
void Pressure_calc_coeff(Pressure *p, MAC_grid *grid, Parameters *params) {


	int i, j, k;
	int NX, NY, NZ;
	double *idx_c, *idy_c, *idz_c;
	double *idx_u, *idy_v, *idz_w;
	double *as, *an, *ap;
	double *rhsw, *xw;
	double *modwaveksq, *modwaveisq;
	double idx, idz;


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

	as = (double *) malloc( (NY-1)*sizeof(double) );
	an = (double *) malloc( (NY-1)*sizeof(double) );
	ap = (double *) malloc( (NY-1)*sizeof(double) );
	p->asw = (double *) malloc( (NY-1)*sizeof(double) );
	p->anw = (double *) malloc( (NY-1)*sizeof(double) );
	p->apw = (double *) malloc( (NY-1)*sizeof(double) );
	p->xw = (double *) malloc( (NY-1)*sizeof(double) );


#ifdef YPERIODIC
	p->qw = (double *) malloc( (NY-1)*sizeof(double) );
	p->u = (double *) malloc( (NY-1)*sizeof(double) );
	p->v = (double *) malloc( (NY-1)*sizeof(double) );
#endif

	p->rhsw = (double *) malloc( (NY-1)*sizeof(double) );
	p->temp = (double *) malloc( (NY-1)*sizeof(double) );
	modwaveksq = (double *) malloc( (NZ-1)*sizeof(double) );
	modwaveisq = (double *) malloc( (NX-1)*sizeof(double) );

	// Coefficient from discretization of d2p/dy2 with Neumann boundary
	// condition at both walls/free surfaces

#ifdef YPERIODIC
	as[0] = idy_v[0]*idy_c[NY-2];
	an[NY-2] = idy_v[NY-2]*idy_c[0];
	for (j=0;j<NY-1;j++) {
		if (j > 0 ){
			as[j] = idy_v[j]*idy_c[j-1];
		}
		if (j < NY-2){
			an[j] = idy_v[j]*idy_c[j];
		}
			ap[j] = -as[j] - an[j];
	}
	an[NY-2] = 0.0;
	as[0] = 0.0;
	/*for (j=0;j<NY-1;j++) {
		printf("ap is %2.5f\n", ap[j]);

	}
	for (j=0;j<NY-1;j++) {
		printf("as is %2.5f\n", as[j]);

	}
	for (j=0;j<NY-1;j++) {
		printf("an is %2.5f\n", an[j]);

	}*/
#else
	as[0] = 0.0;
	an[NY-2] = 0.0;
	for (j=0;j<NY-1;j++) {
		if (j > 0 )
			as[j] = idy_v[j]*idy_c[j-1];
			if (j < NY-2)
			an[j] = idy_v[j]*idy_c[j];
			ap[j] = -as[j] - an[j];
		}
		an[NY-2] = 0.0;
#endif

	p->as = as;
	p->ap = ap;
	p->an = an;
	idx = idx_u[0];
	idz = idz_w[0];

	for (k=0;k<NZ-1;k++) {
#ifdef ZPERIODIC
	#ifdef XPERIODIC
		modwaveksq[k] = 2.0*(1-cos(2.*PI*p->zfft->wave_rehac[k]/(NZ-1)))*(idz*idz);
	#else
		modwaveksq[k] = 2.0*(1-cos(2.*PI*p->zfft->wave[k]/(NZ-1)))*(idz*idz);
	#endif
#else
		modwaveksq[k] = 2.0*(1-cos(PI*k/(NZ-1)))*(idz*idz);
#endif

	}

	for (i=0;i<NX-1;i++) {
#ifdef XPERIODIC
		modwaveisq[i] = 2.0*(1-cos(2.*PI*p->xfft->wave[i]/(NX-1)))*(idx*idx);
#else
		modwaveisq[i] = 2.0*(1-cos(PI*i/(NX-1)))*(idx*idx);
#endif
	}

	p->modwaveisq = modwaveisq;
	p->modwaveksq = modwaveksq;


	return;

}




/******************************************************************************/
/*
 Create FFTW plans, transpose plans and calculate the coefficient of the matrix
 for the 1-D Helmholtz equation.
 */
/******************************************************************************/
void Pressure_setup_lsys_accounting_geometry(Pressure *p, MAC_grid *grid,
		Parameters *params) {

	int maxsize;

	p->xfft = (fft *)malloc(sizeof(fft));
	p->zfft = (fft *)malloc(sizeof(fft));
	p->remap = (pres_transpose *)malloc(sizeof(pres_transpose));

	Pressure_fft_setup(p->xfft, params->NXM);
	Pressure_fft_setup(p->zfft, params->NZM);
	if (params->rank==0)
		printf("fft has been set-up successfully...\n");

	Pressure_transpose_setup(p, grid, params);
	if (params->rank==0)
		printf("pressure transpose has been set-up successfully...\n");
	maxsize = max(grid->NX, grid->NZ);

	p->work  = (fftw_real *)malloc(maxsize*sizeof(fftw_real));
	p->work1 = (fftw_real *)malloc(maxsize*sizeof(fftw_real));
	p->cwork = (fftw_complex *)malloc(maxsize*sizeof(fftw_complex));

 	Pressure_calc_coeff(p, grid, params);

}




/******************************************************************************/
/*
 This function solves the linear system to find pressure (phi) at each cell
 center.
 */
/******************************************************************************/
int Pressure_solve_xz_period(Pressure *p, MAC_grid *grid, Parameters *params) {


	int iters;
	double rnorm;
	double *data_t;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int Iks, Ike;
	int Kis, Kie;
	int Kijs, Kije;
	int kstart;
	int NX, NY, NZ;
	int count, index, ind1, ind2, ind3, ind4, ind5, ind6;
	fftw_real *work, *work1;
	fftw_complex *cwork;
	double ***div, ***p_delta;
	double vtq = 0.0;
	int i, j, k;
	double inz, inx;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Ie = min(Ie, NX-1);
	Je = min(Je, NY-1);
	Ke = min(Ke, NZ-1);

	data_t = p->remap->data_transpose;
	work   = p->work;
	work1  = p->work1;
	cwork   = p->cwork;
	Iks = p->remap->Iks;
	Ike = p->remap->Ike;
	Kis = p->remap->Kis;
	Kie = p->remap->Kie;
	Kijs = p->remap->Kijs;
	Kije = p->remap->Kije;
	inx = 1.0/(NX-1.0);
	inz = 1.0/(NZ-1.0);

	// Copy the RHS of poisson equation to 1-D array data_t
	div = p->rhs;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				data_t[count] = div[k][j][i];
				count = count + 1 ;
			}
		}
	}

	/*
	for (j=Js;j<Je;j++){
		printf("%2.5f\n",div[25][j][25]);

	}
*/
	// Transpose the data so that FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_origtozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			fftw_execute_r2r(p->zfft->rplanf, work, work);
			// Rearrange data in REHAC format
			fftw_halfcomplex_to_rehac(NZ-1, work, work1);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work1[k]*inz;
			}
		}
	}

	// Transpose the data so that FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoxline);

	kstart = Kis;
	if (Kis==0) kstart = 2;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);

	/*------------------------------------------------------------------------*/
	/*
	 Take FFT of the complex data

	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)
	 */
	/*------------------------------------------------------------------------*/
	for (k=kstart;k<Kie;k=k+2){
		ind5 = (k-Kis)*ind1;
		ind6 = (k+1-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			ind4 = ind6 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				cwork[i][0] = data_t[ind3+i];
				cwork[i][1] = data_t[ind4+i];

			}
			fftw_execute_dft(p->xfft->cplanf, cwork, cwork);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = cwork[i][0]*inx ;
				data_t[ind4+i] = cwork[i][1]*inx ;
			}
		}
	}

	// Take FFT of the real data
	if (Kis==0) {
		k=0;
		ind5 = k*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->rplanf, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i]*inx ;
			}
		}
		k=1;
		ind5 = k*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->rplanf, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i]*inx ;
			}
		}


	}

	// Transpose the data so that tridiagonal equation can be solved along
	// y-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetoxyplane);

	// Tridiagonal solve
	ind1 = (NX-1)*(NY-1);
	ind2 = (NX-1);
	for (k=Kijs;k<Kije;k++){
		ind3 = ind1*(k-Kijs);
		for (i=0;i<NX-1;i++){
			ind4 = ind3 + i;

#ifdef YPERIODIC

			for (j=0;j<NY-1;j++){
				p->asw[j] = p->as[j];
				p->anw[j] = p->an[j];
				p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
				p->rhsw[j] = data_t[ind4 + j*ind2];
				p->u[j] = 0.0;
				p->v[j] = 0.0;
				/*if (i==0){
					if ( k==0){
				printf("p->rhs[%d] = %e\n",j,p->rhsw[j]);
			}}*/
			}


			p->u[0] = -p->apw[0];
			p->u[NY-2] = p->asw[NY-2];
			p->v[0] = 1.0;
			p->v[NY-2] = -p->anw[0]/p->apw[0];//-p->asw[0]*p->anw[NY-2]/p->apw[0];

			p->apw[NY-2] = p->apw[NY-2] + (p->asw[NY-2]*p->anw[0]/p->apw[NY-2]);
			p->apw[0] = - 2* p->u[0];



			// Since we have periodic bcs. need to fix a point in the solution.
/*
				 if ((i==0) && (k==0)) {
					j =0;
					p->asw[j] = 0.0;
					p->anw[j] = 0.0;
					p->apw[j] = 1.0;
					p->rhsw[j] = 0.0;
					p->u[j] = 0.0;
				}

			printf("u[0] = %f u[NY-2]= %f\n",p->u[0],p->u[NY-2]);
			printf("v[0] = %f v[NY-2]= %f\n",p->v[0],p->v[NY-2]);
*/
			tridiag_yper(p->asw, p->apw, p->anw, p->rhsw, p->u, p->xw, p->qw, p->temp, NY-1);

			/*if ((i==0) && (k==0)) { vtq = 0.0;}
			else{vtq = DOT(p->v,p->xw)/(1+DOT(p->v,p->qw));}
*/
			vtq = (p->v[0]*p->xw[0] + p->v[NY-2]*p->xw[NY-2])/(1.0 + p->v[0]*p->qw[0]+ p->v[NY-2]*p->qw[NY-2]);
			//printf("pxw is %2.2f\n",p->xw[1]);
			//printf("pqw is %2.2f\n",p->qw[1]);

			for (j=0;j<NY-1;j++){
				data_t[ind4 + j*ind2]= p->xw[j] - vtq * p->qw[j] ;
			}
#else
			for (j=0;j<NY-1;j++){
				p->asw[j] = p->as[j];
				p->anw[j] = p->an[j];
				p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
				p->rhsw[j] = data_t[ind4 + j*ind2];
			}
			// Since we have Neumann bcs. need to fix a point in the solution
			if ((i==0) && (k==0)) {
				j = 0;
				p->asw[j]  = 0.0;
				p->anw[j]  = 0.0;
				p->apw[j]  = 1.0;
				p->rhsw[j] = 0.0 ;
			}

			tridiag(p->asw, p->apw, p->anw, p->rhsw, p->xw, p->temp, NY-1);
			for (j = 0; j < NY-1; j++){
				data_t[ind4 + j*ind2]=p->xw[j];
			}
#endif
		}
	}

	// Transpose the data so that Inverse FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xyplanetoxline);

	/*------------------------------------------------------------------------*/
	/*
	 Take Inverse FFT of the complex data

	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)
	 */
	/*------------------------------------------------------------------------*/
	kstart = Kis;
	if (Kis==0) kstart = 2;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);
	for (k=kstart;k<Kie;k=k+2){
		ind5 = (k-Kis)*ind1;
		ind6 = (k+1-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			ind4 = ind6 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				cwork[i][0] = data_t[ind3+i];
				cwork[i][1] = data_t[ind4+i];
			}
			fftw_execute_dft(p->xfft->cplanb, cwork, cwork);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = cwork[i][0] ;
				data_t[ind4+i] = cwork[i][1] ;
			}
		}
	}
	if (Kis==0) {
		k=0;
		ind5 = k*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->rplanb, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i] ;
			}
		}
		k=1;
		ind5 = k*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->rplanb, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i] ;
			}
		}


	}

	// Transpose the data so that Inverse FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			// Rearrange data to FFTW half-complex format
			fftw_rehac_to_halfcomplex(NZ-1, work, work1);
			fftw_execute_r2r(p->zfft->rplanb, work1, work1);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work1[k];
			}
		}
	}

	// Transpose the data to original layout
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoorig);

	// Copy the solution to p_data array
	p_delta = p->deltap;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				p_delta[k][j][i] = data_t[count];
				count = count + 1 ;
			}
		}
	}

	iters = 0;
	return iters;
}




/******************************************************************************/
/*
 This function solves the linear system to find pressure (phi) at each cell
 center.
 */
/******************************************************************************/
int Pressure_solve_z_period(Pressure *p, MAC_grid *grid, Parameters *params) {

	int iters;
	double rnorm;
	double *data_t;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int Iks, Ike;
	int Kis, Kie;
	int Kijs, Kije;
	int kstart;
	int NX, NY, NZ;
	int count, index, ind1, ind2, ind3, ind4, ind5, ind6;
	fftw_real *work, *work1;
	fftw_complex *cwork;
	double ***div, ***p_delta;
	double vtq;
	int i, j, k;
	double inz, inx;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Ie = min(Ie, NX-1);
	Je = min(Je, NY-1);
	Ke = min(Ke, NZ-1);

	data_t = p->remap->data_transpose;
	work   = p->work;
	work1  = p->work1;
	cwork   = p->cwork;
	Iks = p->remap->Iks;
	Ike = p->remap->Ike;
	Kis = p->remap->Kis;
	Kie = p->remap->Kie;
	Kijs = p->remap->Kijs;
	Kije = p->remap->Kije;
	inx = 0.5/(NX-1.0);
	inz = 1.0/(NZ-1.0);

	// Copy the RHS of poisson equation to 1-D array data_t
	div = p->rhs;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				data_t[count] = div[k][j][i];
				count = count + 1 ;
			}
		}
	}

	// Transpose the data so that FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_origtozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			fftw_execute_r2r(p->zfft->rplanf, work, work);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work[k]*inz;
			}
		}
	}

	// Transpose the data so that FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoxline);

	kstart = Kis;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);

	/*------------------------------------------------------------------------*/
	/*
	 Take FCT of the complex data

	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)

	 However, FCT is taken on each k-plane separately, irrespective of whether
	 they are real or imaginary
	 */
	/*------------------------------------------------------------------------*/
	for (k=Kis;k<Kie;k++){
		ind5 = (k-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->cosplanf, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i]*inx ;
			}
		}
	}

	// Transpose the data so that tridiagonal equation can be solved along
	// y-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetoxyplane);

	// Tridiagonal solve
	ind1 = (NX-1)*(NY-1);
	ind2 = (NX-1);
	for (k=Kijs;k<Kije;k++){
		ind3 = ind1*(k-Kijs);
		for (i=0;i<NX-1;i++){
			ind4 = ind3 + i;

			#ifdef YPERIODIC

						for (j=0;j<NY-1;j++){
							p->asw[j] = p->as[j];
							p->anw[j] = p->an[j];
							p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
							p->rhsw[j] = data_t[ind4 + j*ind2];
							p->u[j] = 0.0;
							p->v[j] = 0.0;
						}

						p->u[0] = -p->apw[0];
						p->u[NY-2] = p->asw[0];
						p->v[0] = 1.0;
						p->v[NY-2] = -p->asw[0]*p->anw[NY-2]/p->apw[0];

						p->apw[NY-2] = p->apw[NY-2] + p->asw[0]*p->anw[NY-2]/p->apw[0];
						p->apw[0] = 2.0*p->apw[0];

						// Since we have periodic bcs. need to fix a point in the solution.
						if ((i==0) && (k==0)) {
								j = 0;
								p->asw[j] = 0.0;
								p->anw[j] = 0.0;
								p->apw[j] = 1.0;
								p->rhsw[j] = 0.0;
							}
						//printf("u[0] = %f u[NY-2]= %f\n",p->u[0],p->u[NY-2]);
						//printf("v[0] = %f v[NY-2]= %f\n",p->v[0],p->v[NY-2]);
						//printf("b[0] = %f b[NY-2]= %f b= %f\n",p->apw[0],p->apw[NY-2],p->apw[1]);

						tridiag_yper(p->asw, p->apw, p->anw, p->rhsw, p->u, p->xw, p->qw, p->temp, NY-1);

						if ((i==0) && (k==0)) { vtq = 0.0;}
						else{vtq = DOT(p->v,p->xw)/(1+DOT(p->v,p->qw));}

						//vtq = DOT(p->v,p->xw)/(1+DOT(p->v,p->qw));

						for (j=0;j<NY-1;j++){
							data_t[ind4 + j*ind2]=p->xw[j] - vtq * p->qw[j] ;
						}
			#else
						for (j=0;j<NY-1;j++){
							p->asw[j] = p->as[j];
							p->anw[j] = p->an[j];
							p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
							p->rhsw[j] = data_t[ind4 + j*ind2];
						}
						// Since we have Neumann bcs. need to fix a point in the solution
						if ((i==0) && (k==0)) {
							j = 0;
							p->asw[j]  = 0.0;
							p->anw[j]  = 0.0;
							p->apw[j]  = 1.0;
							p->rhsw[j] = 0.0 ;
						}

						tridiag(p->asw, p->apw, p->anw, p->rhsw, p->xw, p->temp, NY-1);
						for (j=0;j<NY-1;j++){
							data_t[ind4 + j*ind2]=p->xw[j];
						}
			#endif
		}
	}

	// Transpose the data so that Inverse FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xyplanetoxline);

	/*------------------------------------------------------------------------*/
	/*
	 Take Inverse FCT of the complex data
	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)

	 However, Inverse FCT is taken on each k-plane separately, irrespective of
	 whether they are real or imaginary
	 */
	/*------------------------------------------------------------------------*/
	kstart = Kis;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);
	for (k=Kis;k<Kie;k++){
		ind5 = (k-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->cosplanb, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i] ;
			}
		}
	}

	// Transpose the data so that Inverse FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			fftw_execute_r2r(p->zfft->rplanb, work, work);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work[k];
			}
		}
	}

	// Transpose the data to original layout
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoorig);

	// Copy the solution to p_data array
	p_delta = p->deltap;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				p_delta[k][j][i] = data_t[count];
				count = count + 1 ;
			}
		}
	}

	iters = 0;
	return iters;
}




/******************************************************************************/
/*
 This function solves the linear system to find pressure (phi) at each cell
 center.
 */
/******************************************************************************/
int Pressure_solve_no_period(Pressure *p, MAC_grid *grid, Parameters *params) {


	int iters;
	double rnorm;
	double *data_t;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int Iks, Ike;
	int Kis, Kie;
	int Kijs, Kije;
	int kstart;
	int NX, NY, NZ;
	int count, index, ind1, ind2, ind3, ind4, ind5, ind6;
	fftw_real *work, *work1;
	fftw_complex *cwork;
	double ***div, ***p_delta;
	double vtq;
	int i, j, k;
	double inz, inx;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Ie = min(Ie, NX-1);
	Je = min(Je, NY-1);
	Ke = min(Ke, NZ-1);

	data_t = p->remap->data_transpose;
	work   = p->work;
	work1  = p->work1;
	cwork   = p->cwork;
	Iks = p->remap->Iks;
	Ike = p->remap->Ike;
	Kis = p->remap->Kis;
	Kie = p->remap->Kie;
	Kijs = p->remap->Kijs;
	Kije = p->remap->Kije;
	inx = 0.5/(NX-1.0);
	inz = 0.5/(NZ-1.0);

	// Copy the RHS of poisson equation to 1-D array data_t
	div = p->rhs;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				data_t[count] = div[k][j][i];
				count = count + 1 ;
			}
		}
	}

	// Transpose the data so that FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_origtozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			fftw_execute_r2r(p->zfft->cosplanf, work, work);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work[k]*inz;
			}
		}
	}

	// Transpose the data so that FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoxline);

	kstart = Kis;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);

	/*------------------------------------------------------------------------*/
	/*
	 Take FCT of the complex data

	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)

	 However, FCT is taken on each k-plane separately, irrespective of whether
	 they are real or imaginary
	 */
	/*------------------------------------------------------------------------*/
	for (k=Kis;k<Kie;k++){
		ind5 = (k-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->cosplanf, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i]*inx ;
			}
		}
	}

	// Transpose the data so that tridiagonal equation can be solved along
	// y-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetoxyplane);

	// Tridiagonal solve
	ind1 = (NX-1)*(NY-1);
	ind2 = (NX-1);
	for (k=Kijs;k<Kije;k++){
		ind3 = ind1*(k-Kijs);
		for (i=0;i<NX-1;i++){
			ind4 = ind3 + i;

			#ifdef YPERIODIC

						for (j=0;j<NY-1;j++){
							p->asw[j] = p->as[j];
							p->anw[j] = p->an[j];
							p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
							p->rhsw[j] = data_t[ind4 + j*ind2];
							p->u[j] = 0.0;
							p->v[j] = 0.0;
						}

						p->u[0] = -p->apw[0];
						p->u[NY-2] = p->asw[0];
						p->v[0] = 1.0;
						p->v[NY-2] = -p->asw[0]*p->anw[NY-2]/p->apw[0];

						p->apw[NY-2] = p->apw[NY-2] + p->asw[0]*p->anw[NY-2]/p->apw[0];
						p->apw[0] = 2.0*p->apw[0];

						// Since we have periodic bcs. need to fix a point in the solution.
						if ((i==0) && (k==0)) {
								j = 0;
								p->asw[j] = 0.0;
								p->anw[j] = 0.0;
								p->apw[j] = 1.0;
								p->rhsw[j] = 0.0;
							}
						//printf("u[0] = %f u[NY-2]= %f\n",p->u[0],p->u[NY-2]);
						//printf("v[0] = %f v[NY-2]= %f\n",p->v[0],p->v[NY-2]);
						//printf("b[0] = %f b[NY-2]= %f b= %f\n",p->apw[0],p->apw[NY-2],p->apw[1]);

						tridiag_yper(p->asw, p->apw, p->anw, p->rhsw, p->u, p->xw, p->qw, p->temp, NY-1);

						//if ((i==0) && (k==0)) { vtq = 0.0;}
						//else{vtq = DOT(p->v,p->xw)/(1+DOT(p->v,p->qw));}

						vtq = DOT(p->v,p->xw)/(1+DOT(p->v,p->qw));

						for (j=0;j<NY-1;j++){
							data_t[ind4 + j*ind2]=p->xw[j] - vtq * p->qw[j] ;
						}
			#else
						for (j=0;j<NY-1;j++){
							p->asw[j] = p->as[j];
							p->anw[j] = p->an[j];
							p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
							p->rhsw[j] = data_t[ind4 + j*ind2];
						}
						// Since we have Neumann bcs. need to fix a point in the solution
						if ((i==0) && (k==0)) {
							j = 0;
							p->asw[j]  = 0.0;
							p->anw[j]  = 0.0;
							p->apw[j]  = 1.0;
							p->rhsw[j] = 0.0 ;
						}

						tridiag(p->asw, p->apw, p->anw, p->rhsw, p->xw, p->temp, NY-1);
						for (j=0;j<NY-1;j++){
							data_t[ind4 + j*ind2]=p->xw[j];
						}
			#endif
		}
	}

	// Transpose the data so that Inverse FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xyplanetoxline);

	/*------------------------------------------------------------------------*/
	/*
	 Take Inverse FCT of the complex data

	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)

	 However, Inverse FCT is taken on each k-plane separately, irrespective of
	 whether they are real or imaginary
	 */
	/*------------------------------------------------------------------------*/
	kstart = Kis;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);
	for (k=Kis;k<Kie;k++){
		ind5 = (k-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->cosplanb, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i] ;
			}
		}
	}

	// Transpose the data so that Inverse FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			fftw_execute_r2r(p->zfft->cosplanb, work, work);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work[k];
			}
		}
	}

	// Transpose the data to original layout
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoorig);

	// Copy the solution to p_data array
	p_delta = p->deltap;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				p_delta[k][j][i] = data_t[count];
				count = count + 1 ;
			}
		}
	}

	iters = 0;
	return iters;
}




/******************************************************************************/
/*
 This function solves the linear system to find pressure (phi) at each cell
 center.
 */
/******************************************************************************/
int Pressure_solve_x_period(Pressure *p, MAC_grid *grid, Parameters *params) {

	int iters;
	double rnorm;
	double *data_t;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int Iks, Ike;
	int Kis, Kie;
	int Kijs, Kije;
	int kstart;
	int NX, NY, NZ;
	int count, index, ind1, ind2, ind3, ind4, ind5, ind6;
	fftw_real *work, *work1;
	fftw_complex *cwork;
	double ***div, ***p_delta;
	double vtq;
	int i, j, k;
	double inz, inx;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Ie = min(Ie, NX-1);
	Je = min(Je, NY-1);
	Ke = min(Ke, NZ-1);

	data_t = p->remap->data_transpose;
	work   = p->work;
	work1  = p->work1;
	cwork   = p->cwork;
	Iks = p->remap->Iks;
	Ike = p->remap->Ike;
	Kis = p->remap->Kis;
	Kie = p->remap->Kie;
	Kijs = p->remap->Kijs;
	Kije = p->remap->Kije;
	inx = 1.0/(NX-1.0);
	inz = 0.5/(NZ-1.0);



	// Copy the RHS of poisson equation to 1-D array data_t
	div = p->rhs;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				data_t[count] = div[k][j][i];
				count = count + 1 ;
			}
		}
	}

	// Transpose the data so that FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_origtozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			fftw_execute_r2r(p->zfft->cosplanf, work, work);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work[k]*inz;
			}
		}
	}



	// Transpose the data so that FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoxline);

	kstart = Kis;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);

	/*------------------------------------------------------------------------*/
	/*
	 Take FCT of the complex data

	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)

	 However, FCT is taken on each k-plane separately, irrespective of whether
	 they are real or imaginary
	 */
	/*------------------------------------------------------------------------*/
	for (k=Kis;k<Kie;k++){
		ind5 = (k-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->rplanf, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i]*inx ;
			}
		}
	}

	// Transpose the data so that tridiagonal equation can be solved along
	// y-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetoxyplane);

	// Tridiagonal solve
	ind1 = (NX-1)*(NY-1);
	ind2 = (NX-1);
	for (k=Kijs;k<Kije;k++){
		ind3 = ind1*(k-Kijs);
		for (i=0;i<NX-1;i++){
			ind4 = ind3 + i;

			#ifdef YPERIODIC

						for (j=0;j<NY-1;j++){
							p->asw[j] = p->as[j];
							p->anw[j] = p->an[j];
							p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
							p->rhsw[j] = data_t[ind4 + j*ind2];
							p->u[j] = 0.0;
							p->v[j] = 0.0;
						}

						p->u[0] = -p->apw[0];
						p->u[NY-2] = p->asw[0];
						p->v[0] = 1.0;
						p->v[NY-2] = -p->asw[0]*p->anw[NY-2]/p->apw[0];

						p->apw[NY-2] = p->apw[NY-2] + p->asw[0]*p->anw[NY-2]/p->apw[0];
						p->apw[0] = 2.0*p->apw[0];

						// Since we have periodic bcs. need to fix a point in the solution.
						if ((i==0) && (k==0)) {
								j = 0;
								p->asw[j] = 0.0;
								p->anw[j] = 0.0;
								p->apw[j] = 1.0;
								p->rhsw[j] = 0.0;
							}
						//printf("u[0] = %f u[NY-2]= %f\n",p->u[0],p->u[NY-2]);
						//printf("v[0] = %f v[NY-2]= %f\n",p->v[0],p->v[NY-2]);
						//printf("b[0] = %f b[NY-2]= %f b= %f\n",p->apw[0],p->apw[NY-2],p->apw[1]);

						tridiag_yper(p->asw, p->apw, p->anw, p->rhsw, p->u, p->xw, p->qw, p->temp, NY-1);

						//if ((i==0) && (k==0)) { vtq = 0.0;}
						//else{vtq = DOT(p->v,p->xw)/(1+DOT(p->v,p->qw));}

						vtq = DOT(p->v,p->xw)/(1+DOT(p->v,p->qw));

						for (j=0;j<NY-1;j++){
							data_t[ind4 + j*ind2]=p->xw[j] - vtq * p->qw[j] ;
						}
			#else
						for (j=0;j<NY-1;j++){
							p->asw[j] = p->as[j];
							p->anw[j] = p->an[j];
							p->apw[j] = p->ap[j] - p->modwaveksq[k] - p->modwaveisq[i];
							p->rhsw[j] = data_t[ind4 + j*ind2];
						}
						// Since we have Neumann bcs. need to fix a point in the solution
						if ((i==0) && (k==0)) {
							j = 0;
							p->asw[j]  = 0.0;
							p->anw[j]  = 0.0;
							p->apw[j]  = 1.0;
							p->rhsw[j] = 0.0 ;
						}

						tridiag(p->asw, p->apw, p->anw, p->rhsw, p->xw, p->temp, NY-1);
						for (j=0;j<NY-1;j++){
							data_t[ind4 + j*ind2]=p->xw[j];
						}
			#endif
		}
	}

	// Transpose the data so that Inverse FFT/FCT can be taken along i-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xyplanetoxline);

	/*------------------------------------------------------------------------*/
	/*
	 Take Inverse FCT of the complex data

	 Note k=0 and k=1 have real data whereas (k=2 and k=3) and (k=4 and k=5),...
	 have complex data (real and imaginary number)

	 However, Inverse FCT is taken on each k-plane separately, irrespective of
	 whether they are real or imaginary
	 */
	/*------------------------------------------------------------------------*/
	kstart = Kis;
	ind1 = (Je-Js)*(NX-1);
	ind2 = (NX-1);
	for (k=Kis;k<Kie;k++){
		ind5 = (k-Kis)*ind1;
		for (j=Js;j<Je;j++){
			ind3 = ind5 + (j-Js)*ind2;
			for (i=0;i<NX-1;i++){
				work[i] = data_t[ind3+i];
			}
			fftw_execute_r2r(p->xfft->rplanb, work, work);
			for (i=0;i<NX-1;i++){
				data_t[ind3+i] = work[i] ;
			}
		}
	}

	// Transpose the data so that Inverse FFT/FCT can be taken along k-direction
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_xlinetozline);

	ind1 = (Je-Js)*(Ike-Iks);
	ind2 = (Ike-Iks);
	for (j=Js; j<Je; j++){
		ind3 = (j-Js)*ind2;
		for (i=Iks;i<Ike;i++){
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				work[k] = data_t[index];
			}
			fftw_execute_r2r(p->zfft->cosplanb, work, work);
			for (k=0;k<NZ-1;k++){
				index = k*ind1 + ind3 + (i-Iks);
				data_t[index] = work[k];
			}
		}
	}

	// Transpose the data to original layout
	remap_3d(data_t, data_t, NULL, p->remap->remap_plan_zlinetoorig);

	// Copy the solution to p_data array
	p_delta = p->deltap;
	count = 0;
	for (k=Ks;k<Ke;k++){
		for (j=Js;j<Je;j++){
			for (i=Is;i<Ie;i++){
				p_delta[k][j][i] = data_t[count];
				count = count + 1 ;
			}
		}
	}



	iters = 0;
	return iters;
}
