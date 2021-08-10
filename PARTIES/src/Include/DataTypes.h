#ifndef DATATYPES_H
	#define DATATYPES_H

	#include "mpi.h"
	#include "fftw3.h"
	#include "Boundary.h"
	#include "remap_3d.h"
	#include "definitions.h"

	#include <complex.h>

extern MPI_Comm comm3d;
extern MPI_Datatype MPI_PARTICLE;
extern MPI_Datatype MPI_COLLISION;

/******************************************************************************/
/*                             LSOLVER_TRANSPOSE                              */
/******************************************************************************/
struct lsolver_transpose {
	struct remap_plan_3d *remap_plan_origtoyline;
	struct remap_plan_3d *remap_plan_ylinetoorig;
	double *data_transpose;

	int Ijs, Ije;
};
typedef struct lsolver_transpose lsolver_transpose;


/******************************************************************************/
/*                               CONCENTRATION                                */
/******************************************************************************/
struct concentration {

	double Pe;                // Corresponding Peclet number
	double v_settl0;          // Constant settling speed
	double W_total_susp_mass; // Total suspended mass within the fluid region
	double W_front_location;  // x-location of gravity current front

	double ***data;
	double ***ng_data_old;
	double ***ng_rhs;

	// Temp array if BQUICK is used
	double ***data_temp, ***ng_temp1;

	// Blending function if BQUICK is used
	double ***blend;

	// Global and local total concentration fields (it is generated only if NConc > 1 and stored in c[0]
	double ***c_total;

	double ***ng_conv, ***ng_conv_old, ***ng_viscous;
	double ***ng_conv_quick;
	double ***ng_visc_explicit;

	//	Extra array for CG method !!!!! I allocate memory for the first field only !!!!
	double ***d;
	double ***ng_r;
	double ***ng_Ad;


	double *G_ave_height_x; // (concentration) Averaged flow height. Only a function of x
	double *W_ave_height_x; // World ave_height
	double **G_deposit_height; // deposit height for each concentration field on current processor
	double **W_deposit_height; // World, deposit height for each concentration field
	double **G_deposit_height_dumped; /* deposit height for each concentration field on current processor (after dumping whatever is left in the fluid region) */
	double **W_deposit_height_dumped; /* World, deposit height for each concentration field (after dumping whatever is left in the fluid region) */
	double **G_conc_immersed, **W_conc_immersed;

	// Inflow profile from precursor simulation
	double **inflow_n;
	double **inflow_o;
	double **inflow;

	// local and world stokes dissipation rate
	double G_stokes_dissipation_rate;
	double W_stokes_dissipation_rate;

	// Active and passive potential energies
	double G_Ep_active, G_Ep_passive;
	double W_Ep_active, W_Ep_passive;

	// Averaged concentration
	double integral_conc_full; // integral over the full domain
	double integral_conc_fluid; // we subtract the particle
    double dt_integral_conc_full;
    double dt_integral_conc_fluid;

	int heating_flag;   // if == 0 no heating has been applied in the current RK step; used to calculate Q_int
	double Q_int[3]; //  ibm heating for k=1,2,3


	int Type; // PARTICLE, SALINITY or TEMPERATURE

	double *as, *ap, *an;
	double *asw, *apw, *anw;
	double *xw, *rhsw;

	// Following coefficients are needed when QUICK is used
	double *aeW, *aeE, *aeEE;
	double *awWW, *awW, *awE;
	double *anS, *anN, *anNN;
	double *asSS, *asS, *asN;
	double *afB, *afF, *afFF;
	double *abBB, *abB, *abF;


	double solve_cpu_time;
	double solve_first_copy_cpu_time, solve_second_copy_cpu_time;
	double solve_first_remap_cpu_time, solve_second_remap_cpu_time;
	double solve_tridiag_cpu_time;

	double obounds_blend_set_cpu_time, obounds_ifcheck_cpu_time;
	double obounds_mpi1_cpu_time, obounds_mpi2_cpu_time, obounds_mpi3_cpu_time;
	double obounds_comm_cpu_time, obounds_blendingfactor_cpu_time;
	double obounds_total_cpu_time;

	// VOF method
	double vol_heat_cap_particle;


};
typedef struct concentration Concentration;


/******************************************************************************/
/*                                    FFT                                     */
/******************************************************************************/
typedef double fftw_real;
//typedef double fftw_complex[2];
struct fft {

	int N;
	fftw_plan rplanf, rplanb;
	fftw_plan cplanf, cplanb;
	fftw_plan cosplanf, cosplanb;
	int *wave, *wave_rehac;
};
typedef struct fft fft;


/******************************************************************************/
/*                               PRES_TRANSPOSE                               */
/******************************************************************************/
struct pres_transpose {

	struct remap_plan_3d *remap_plan_origtozline;
	struct remap_plan_3d *remap_plan_zlinetoxline;
	struct remap_plan_3d *remap_plan_xlinetoxyplane;

	struct remap_plan_3d *remap_plan_xyplanetoxline;
	struct remap_plan_3d *remap_plan_xlinetozline;
	struct remap_plan_3d *remap_plan_zlinetoorig;
	double *data_transpose;

	int Iks, Ike;
	int Kis, Kie;
	int Kijs, Kije;
};
typedef struct pres_transpose pres_transpose;


/******************************************************************************/
/*                                  PRESSURE                                  */
/******************************************************************************/
struct pressure {

	double ***p_data;
	double ***deltap;
	double ***rhs;
	double *idxdt, *idydt, *idzdt;
#ifdef POST_PROCESS
	double ***p_data_avg;
#endif

	fft *xfft, *zfft;
	pres_transpose *remap;
	fftw_real *work, *work1;
	fftw_complex *cwork;

	double *as, *ap, *an;
	double *modwaveksq, *modwaveisq;
	double *asw, *apw, *anw;
	double *xw, *rhsw;
	double *qw;
	double *u, *v;
	double *temp;

	double project_comm_cpu_time, project_update_cpu_time;
	double project_u_cpu_time, project_v_cpu_time, project_w_cpu_time;

	double compute_divergence_comm_cpu_time;
	double compute_divergence_loop_cpu_time;
	double compute_divergence_reduce_cpu_time;
};
typedef struct pressure Pressure;


/******************************************************************************/
/*                                 PARAMETERS                                 */
/******************************************************************************/
struct parameters {

	/*--------------------------------- MPI ----------------------------------*/
	int rank; // Current processor
	int size; // Total number of processors

	// Number of processors in each direction
	int NPX, NPY, NPZ;

	// Processor coordinates
	int xproccoord, yproccoord, zproccoord;

	// Neighboring processors in each direction, -1 indicates no neighbor
	int npxplus, npxminus, npyplus, npyminus, npzplus, npzminus;

	// Communicators in each direction and plane
	MPI_Comm  new_comm, commx, commy, commz;
	MPI_Comm commxy, commyz, commzx;

	double *send_buffer, *recv_buffer;


	/*------------------------------- GEOMETRY -------------------------------*/
	// Lower and upper bounds of domain
	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;

	// Length of domain in each direction
	double Lx, Ly, Lz;


	/*--------------------------------- GRID ---------------------------------*/
	// Number of grid cells in each direction
	int NXM, NYM, NZM;

	// Boolean to either import a grid or generate a uniform grid
	int ImportGridFromFile;


	/*------------------------------ SIMULATION ------------------------------*/
	// Variables
	double dt, dt_old, dt_substep;
	double time;
	double t_in_n, t_in_o;
	int ntime;
	int noutput;
	int noutput_2d;
	int par_trn;
	double output_time;
	double output_time_2d;
	int which_stage;
	int ninflow;

	// Inputs
	double time_max;
	double output_time_interval;
	double output_time_interval_2d;
	double default_dt;
	double max_dt;
	int constant_dt;
	double cfl;
	int ghost_nodes;
	int read_new_data;

	// Flag which indicates if the simulation starts from t=0 or from a previous
	// run
	int resume;


	/*-------------------------------- LSOLVE --------------------------------*/
	lsolver_transpose *lsolver_remap;

	// Tolerance and maximum number of iterations for CG solver
	double CG_ETOL;
	int CG_MAXIT;


	/*--------------------------------- FLOW ---------------------------------*/
	// Reynolds number
	double Re;

	// Shear parameter
	double G;

	//Convective velocity at the exit plane (needed when with "OUTFLOW" bc)
	double U_conv_outflow;

	// Pressure gradient
	double dp_dx, dp_dx_old;
    double amplitude;
	double tref;
	double startup_time;    // starting time to release the particle (in STARTUP flag)
	double phase_shift_factor;     // phase shift of the oscillating force= phase_shift_factor * Pi

	// Streamwise bulk velocity
	double ubulk, ubulk_old;
	double ubulk_target;
	double u_in_left;
	double u_in_right;
	double shear_sigma_u;
	double shear_sigma_conc;
	double delta_jet;
	double F_jet;
	double V_swarm;
	double y_swarm_start;
	double swarm_rand;
	double voutflow;

	// Parameters for flow velocity initialization
	int vel_init_type;
	double vel_init_y0;


	//Parameters for turbulent forcing
	double Kf;
	double Ds;

	//Lock release
	double t_release;
	double per_fr;


	/*--------------------------------- CONC ---------------------------------*/
	// Number of concentration fields
	int NConc;
	int Nbin;

	// Peclet numbers
	double *Pe;  // For each scalar there is one parameter

	// Contribution of each conc field in the buoyancy term
	double *richardson;

	// Settling velocity of particles: y-direction
	double *V_s0;

	// Type of conc initialization
	int *conc_init_type;

	double *BC_AN;
	double *BC_AS;
	double *BC_AE;
	double *BC_AW;
	double *BC_AF;
	double *BC_AB;

	double *BC_BN;
	double *BC_BS;
	double *BC_BE;
	double *BC_BW;
	double *BC_BF;
	double *BC_BB;

	double *BC_CN;
	double *BC_CS;
	double *BC_CE;
	double *BC_CW;
	double *BC_CF;
	double *BC_CB;


	// Standard deviation of the initial concentration profile
	double sig_profile;

	// min value and amplitude. To be used in as many init profiles as possible
	double cbd0, cbd1, cbd2, cbd3, cbd4, cbd5;

	/*--------------------------------- LOCK ---------------------------------*/
	double x_fr, y_fr, z_fr, hlock;
	int lock_smooth;


	/*------------------------------- IBM-PARTICLE -------------------------------*/
	// Number of additional loops enforcing no-slip after solving viscous terms
	int N_forcing_loops;

	// Number of additional loops enforcing conc BC after solving diffusive terms
	int N_heating_loops;
	int N_impheating_loops;

	// Particle density
	double rho_s;

	// Gravitational acceleration
	double *grav; // used only for the particle momentum balance!!

	// Particle volumetric heat capacity normalize with the fluid's one
	double *vol_heat_cap_s;  // for a solute this is one

	// Particle conductivity normalized with the fluid value
	double *conductivity_s;   // for a solute this is the normalized diffusivity

	// initializes a layer of fluid from ymin to ymin+layer_height
	double *layer_height;

	// Used to smooth velocity field around swimmers to remove viscous BL
	double bl_thick;

#ifdef STARTUP
	int startup_flag;
#endif

	// Collision parameters
	double e_dry_wall, e_dry_particles;  // Restitution coefficients
	double Psi_crit;  // Critical impact angle for ATFM
	double mu_k, mu_s;  // Kinetic and static friction coefficients
	double PoissonsRatio;
	double roughness;
	double lub_range;  // Range of lubrication force divided by grid cell size
	double erm_range;  // Range of electrostatic force divided by grid cell size
	double F_erm;      // Electrostatic repulsion coefficient
	double coh_range;  // Range of cohesive force divided by grid cell size
	double D50;        // Median grain size diameter
	double Co;         // Cohesive Number

	// Number of fluid timesteps per collision
	int Ndt_coll;

	// Squirmer model
	double B1;
	double B2;
	// Swimmers target (light source)
	double *target_coord;
	int target_on;
	double coefA;
	// Vertical cut-off point (above which swimmers get turned off)
	double y_cut;


	/*  						Output-switches								  */
	/* switches to control output  1 == on  anything else  is off */

	int avg_dir;
	int *conc_output_deposit_height;
	int *conc_output_dump;
	int *conc_output_stokes_diss_rate;
	int *conc_output_integral;
	int *conc_output_Q_int;


	int ave_height_output, susp_mass_output;
	int front_location_output, front_speed_output;
	int energies_output;
	int shear_stress_output;
	int sedim_rate_output;
	int ImportBottomInterface;


	int output_vfu;
	int output_vfv;
	int output_vfw;
	int output_vfc;

	// switches for statistics2D

	int nu_bottom_switch;
	int current_height_switch;
	int front_position_switch;
	int turbulent_terms_switch;
	int nusselt_switch;
	int post_processing_switch;
	int near_wall_slice_switch;

};
typedef struct parameters Parameters;


/******************************************************************************/
/*                                  VELOCITY                                  */
/******************************************************************************/
struct velocity {

	/*------------------------------------------------------------------------*/
	/*
	 Which component of the velocity
	    'u', 'v' or 'w'
	 */
	/*------------------------------------------------------------------------*/
	char component;

	/*------------------------------------------------------------------------*/
	/*
	     - data    : value of the velocity
	     - data_bc : value of the velocity at the cell_center
	 */
	/*------------------------------------------------------------------------*/

	// Global values
	double ***data;
	double ***data_bc;
	double ***ng_rhs;
	double ***data_old;

#ifdef VOF_SCALAR
	double ***data_vof;
#endif

	// Convective and viscous terms
	double ***ng_explicit, ***ng_explicit_old, ***ng_implicit;

#ifdef TURB_FORCING
	double ***fturb;
#endif

#ifdef CG_SOLVE
	double ***d;
	double ***ng_r;
	double ***ng_Ad;
#endif

#ifdef BICG_SOLVE
	double ***p;
	double ***s;
	double ***ng_r;
	double ***ng_r0;
	double ***Ap;
	double ***ng_As;
#endif

	// Normal wall velocity gradients (new and old) for left wall outflow
	double **d_dn_l_n;
	double **d_dn_l_o;

	// Normal wall velocity gradients (new and old) for right wall outflow
	double **d_dn_r_n;
	double **d_dn_r_o;

	// Inflow profile from precursor simulation
	double **inflow_n;
	double **inflow_o;
	double **inflow;


	// Total influx to the region
	double G_influx;

	// (World) influx to the region (same on all processors)
	double W_influx;

	// Global and World kinetic energy of the fluid
	double G_kinetic_energy, W_kinetic_energy;

	// Global and World viscous dissipation rate
	double G_dissipation_rate, W_dissipation_rate;

	// Total shear stress on the bottom boundary. Store them only in 'u'
	// velocity
	double **G_shear_stress_bottom, **W_shear_stress_bottom;

	// Tangential velocity on the bottom boundary
	double **G_u_shear, **G_v_shear, **G_w_shear;
	double **W_u_shear, **W_v_shear, **W_w_shear;

	double **G_u_streak, **W_u_streak;

	// total number of points in the linear system (number of fluid points)
	long int lsys_count;

	double *as, *ap, *an;
	double *asw, *apw, *anw;
	double *xw, *rhsw;


	double solve_cpu_time;
	double solve_first_copy_cpu_time, solve_second_copy_cpu_time;
	double solve_first_remap_cpu_time, solve_second_remap_cpu_time;
	double solve_tridiag_cpu_time;

	double cell_center_comm_cpu_time;
	double cell_center_u_cpu_time, cell_center_v_cpu_time, cell_center_w_cpu_time;
	double rhs_p_cpu_time, rhs_loop_cpu_time;
};
typedef struct velocity Velocity;




/******************************************************************************/
/*                                 Fourier Space                              */
/******************************************************************************/
struct fourier {

	int nb_forced_freq;
	int **k_list;
	double std_dev;
	double complex **b_random;
	double complex **f_fourier;

};
typedef struct fourier Fourier;



/******************************************************************************/
/*                                 Vorticity                                  */
/******************************************************************************/
struct vorticity {

	double ***vor_x;
	double ***vor_y;
	double ***vor_z;

	double ***vor_x_bc;
	double ***vor_y_bc;
	double ***vor_z_bc;

};
typedef struct vorticity Vorticity;


/******************************************************************************/
/*                                 VISCOSITY                                  */
/******************************************************************************/
/* Strucutre holding the information for all the immersed node for each quantity */
struct viscosity  {

	double ***nu;
	double ***nuX;
	double ***nuY;
	double ***nuZ;

};
typedef struct viscosity Viscosity;


/******************************************************************************/
/*                                   SCALAR                                   */
/******************************************************************************/
/* S = S(X1, X2, X3) */
struct scalar {

	double X1, X2, X3;
	double Value;
};
typedef struct scalar Scalar;


/******************************************************************************/
/*                                  INDICES                                   */
/******************************************************************************/
/* To keep the 3 index corresponding to variable (u, v, w, c or p) node */
struct indices {

	int x_index;
	int y_index;
	int z_index;
};
typedef struct indices Indices;


/******************************************************************************/
/*                                   VECTOR                                   */
/******************************************************************************/
/* vector (vx,vy,vz); */
struct vector_type {

	double vx, vy, vz;

};
typedef struct vector_type VectorType;


/******************************************************************************/
/*                                 POINT_TYPE                                 */
/******************************************************************************/
/* point (x,y,z); */
struct point_type {

	double x, y, z;

};
typedef struct point_type PointType;


/******************************************************************************/
/*                               IMMERSED_NODE                                */
/******************************************************************************/
/* This sttucture holds the information used for one immersed node */
struct immersed_node {

	// IBM_MAX: 3 Number of fluid nodes used in the interpolation stencil

	// Boundary point (located on the boundary used to enforce the correct
	// boundary conditions
	PointType boundary_point;

	// Coordinates of the immersed node. First fluid node off of the solid wall
	PointType im_point;

	// Data at the intersection_point and boundary_point are used to obtain
	// data at immersed node
	PointType intersection_point;

	// Unit vector normal to the boundary at boundary point
	VectorType n;

	// Coefficients used to relate the value of the IB node to the neigbhoring
	// fluid nodes satisfying the correct boundary conditions
	double fluid_coef[IBM_MAX];

	// (i,j,k) indices of the chosen fluid nodes for interpolation
	Indices fluid_index[IBM_MAX];

	// (i,j,k) index of the immersed node
	Indices im_index;

	// Used to impose any nonzero boundary conditions by taking the nonzero
	// term to the rhs of the lsys
	double boundary_coef;
};
typedef struct immersed_node ImmersedNode;


/******************************************************************************/
/*                                  IMMERSED                                  */
/******************************************************************************/
/* Strucutre holding the information for all the immersed node for each quantity */
struct immersed  {

	// Which flow quantity 'u', 'v', 'w' or 'c'
	char quantity;

	// Total number of the ib_nodes on the current processor
	int N;

	// Array of all the ib_nodes
	ImmersedNode *ib_nodes;

	// Type of the boundary condition to be imposed: NEUMANN, DIRICHLET or ROBIN
	int boundary_condition;
};
typedef struct immersed Immersed;


/******************************************************************************/
/*                                  MAC_GRID                                  */
/******************************************************************************/
/* Structure holding the grid information */
struct mac_grid {

	double *xc, *yc, *zc;
	double *xu, *yv, *zw;
	int NX, NY, NZ, NT;
	int NI, NJ, NK;
	double **interface_position;
	int   **interface_y_index;
	double *finer_1d_interface_x, *finer_1d_interface_y;
	int fine_step;

	// Index of the start and end corner on the current processor
	int G_Is, G_Js, G_Ks;
	int G_Ie, G_Je, G_Ke;

	// Index of the start and end corner on the current processor including
	// ghost nodes
	int L_Is, L_Js, L_Ks;
	int L_Ie, L_Je, L_Ke;

	// Arrays to store the condition of a grid
	int ***u_status;
	int ***v_status;
	int ***w_status;
	int ***c_status;

	// Flow field quantities holding the required information for the
	// adaptation of the immerse boundary method
	Immersed *u_immersed, *v_immersed, *w_immersed, *c_immersed;
	Immersed *nut_immersed;

	// Number of ghost and noghost nodes, for cycling through with a single loop
	int total_nodes, ng_total_nodes;

	double *dx_u, *dy_v, *dz_w;
	double dy_min;			// global minimum
	double *dx_c, *dy_c, *dz_c;

	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double *i2dx_c, *i2dy_c, *i2dz_c;
#ifndef GRID_UNIFORM
	double *wc2vN, *wc2vS, *wc2uW, *wc2uE, *wc2wB, *wc2wF;
#endif
#ifdef BICG_SOLVE
	double *idx2_e, *idx2_w;
	double *idy2_n, *idy2_s;
	double *idz2_f, *idz2_b;
#endif

	// Arrays to store signed distance from the wall when Immersed boundary
	// method is used
	double ***u_sdf, ***v_sdf, ***w_sdf, ***c_sdf;
	// for the sheared periodic config
	double **exchange_slab;
};
typedef struct mac_grid MAC_grid;


/******************************************************************************/
/*                                   TIMER                                    */
/******************************************************************************/
/* Structure holding timing information */
#define N_TIMER_ELEMENTS 34
struct timer {

	// Simulation
	double Wtime_total;
	double Wtime_init;
	double Wtime_intEOM;
	double Wtime_output;
	double Wtime_output_2d;
	double Wtime_cfl;

	// Communication
	double Wtime_comm_3D;
	double Wtime_comm_2D_reduce;

	// Velocity
	double Wtime_vel_convective;
	double Wtime_vel_rhs;
	double Wtime_vel_solve;
	double Wtime_vel_boundaries;
	double Wtime_vel_cell_center;

	// Pressure
	double Wtime_p_rhs;
	double Wtime_p_solve;
	double Wtime_p_project;
	double Wtime_p_divergence;

	// Concentration
	double Wtime_c_total;
	double Wtime_c_convective;
	double Wtime_c_rhs;
	double Wtime_c_solve;
	double Wtime_c_outofbounds;

	// Subgrid
	double Wtime_sgs_total;
	double Wtime_sgs_filter;
	double Wtime_sgs_dynamic;

	// RANS
	double Wtime_rans_total;
	double Wtime_rans_convective;
	double Wtime_rans_rhs;
	double Wtime_rans_solve;

	// Particle
	double Wtime_particle_total;
	double Wtime_particle_comm;
	double Wtime_particle_coll;
	double Wtime_particle_forc;
	double Wtime_particle_int;
};
typedef struct timer Timer;


/******************************************************************************/
/*                                STATISTICS2D                                */
/******************************************************************************/
struct statistics2d {

	double ***slice_c, ***mean_c, **mean_u, **mean_v, **mean_w, **mean_p;
	double **mean_nut, **mean_cev, **mean_sct;
	double **turbStress_uu, **turbStress_vv, **turbStress_ww;
	double **turbStress_uv, **turbStress_uw, **turbStress_vw;
	double ***turbStressC;
	double **flux_uu, **flux_vv, **flux_ww;
	double **flux_uv, **flux_uw, **flux_vw;
	double ***nusselt;
	double ***fluxC;
	double **turbKinE;
	double **kinE, **Urms;

	double **histo, **range;


	double **pwork;

	double ***potEnergyC ;
	double ***grad_y_c;

	double **potEnergyP;


	double **viscDiss, **viscTran, **sgsDiss, **sgsTran;
	double **work, **work1, **work2, **work3, **work4, **vf_avg;
	double **mean_tke, **mean_eps;

	// 1D variables in x
	double **current_height; //  h(x,t)= < c >_{yz}*L_y
	double **current_height_fluid; //  h(x,t)= < c phi_f>_{yz}*L_y


	// 0D variables
	double *front_position;  //  x_F= max{ x: h(x,t)> 0.01  }
	double integral_kinEnergy_fluid;
	double integral_RMS_fluid;
	double integral_forcing_fluid;
	double integral_kinEnergy_part;
	double integral_viscDiss;
	double *integral_buoyantwork;

	double *integral_potEnergy_fluid;
	double *integral_nusselt;
	double integral_potEnergy_part;
	double *integral_conc_fluid;
	double *integral_conc_full;

	double **avg_uc;
	double **avg_vc;
	double **avg_wc;
	double **avg_c;
	double *avg_vf;
	double *work_avg;

};
typedef struct statistics2d Statistics2d;

/******************************************************************************/
/*                          Strain rate tensor                                */
/******************************************************************************/
struct strain_rate {

	double ***strain;

	// Variables to store the strain rate tensor
	double ***s11, ***s22, ***s33, ***s12, ***s13, ***s23;
};
typedef struct strain_rate Strain_rate;

/******************************************************************************/
/*                          Wall-model 		                                  */
/******************************************************************************/
struct wall_model {

	// Schumann Wall model related variables
	double **dudy_wm_bottom, **dudy_wm_top;
	double **dwdy_wm_bottom, **dwdy_wm_top;

	// Log law related variables
	double kappa, b, ypintr;
};
typedef struct wall_model Wall_model;

/******************************************************************************/
/*                                  CONC_DEV                                  */
/******************************************************************************/
/* conc_dev (dynamic eddy viscosity model) holds the ILM, IMM variable etc used in DEV*/
struct conc_dev {

	// SGS Schmidt number
	double ***Sct;

	// Modfied SGS Schmidt number (i.e. mSct = nut/Sct)
	double ***mSct;

	// Varaiables used in Dynamic model calculation
	double ***ng_TT, ***ng_KT;

	// Varaiables used in Lagrangian version of the dynamic model
	double ***ITT, ***IKT;
	double ***ng_ITTO, ***ng_IKTO;

	// Variables to store filtered concentration field
	double ***ccf;
};
typedef struct conc_dev ConcDev;


/******************************************************************************/
/*                                  SUBGRID                                   */
/******************************************************************************/
struct subgrid {

	// Smagorinsky
	double Cs;
	double Sct;

	// Subgrid viscosity
	double ***nut;

	// Variables to store the strain rate tensor
	Strain_rate *st_rate;

	// Length scale
	double ***ng_lengthscale_sq;

	// Dynamic model related variables
	double *wxp,*wxe,*wxw;
	double *wyp,*wys,*wyn;
	double *wzp,*wzb,*wzf;

	// Variables to store coefficient etc
	double ***ILM, ***IMM;
	double ***ng_Cev;
	double ***ng_LM, ***ng_MM;
	double ***ng_ILMO, ***ng_IMMO;

	// Variables to store filtered cell-centered velocity components
	double ***ucf, ***vcf, ***wcf;

	// "Intermediate" variables
	double ***strainf, ***sijf, ***strainsijf, ***uu;

	// "Work" variables
	double ***work1, ***work2, ***work3;
	double ***uk;

	ConcDev **cdev;

	// Schumann Wall model related variables
	Wall_model *log_law;

	double filter_cpu_time;
	double filter_comm1_cpu_time, filter_comm2_cpu_time, filter_comm3_cpu_time;
	double filter_comm4_cpu_time;
	double filter_x_cpu_time, filter_y_cpu_time, filter_z_cpu_time;

	double dynamic_cpu_time;
	double dynamic_1_cpu_time, dynamic_2_cpu_time, dynamic_3_cpu_time;
	double dynamic_4_cpu_time, dynamic_5_cpu_time, dynamic_6_cpu_time;
	double dynamic_7_cpu_time, dynamic_8_cpu_time, dynamic_9_cpu_time;
	double dynamic_10_cpu_time, dynamic_11_cpu_time, dynamic_12_cpu_time;
	double dynamic_13_cpu_time, dynamic_14_cpu_time, dynamic_15_cpu_time;
};
typedef struct subgrid Subgrid;


/******************************************************************************/
/*                             Two Eqn RANS model                             */
/******************************************************************************/
struct two_equation_rans {

	double ***data;
	double ***ng_explicit, ***ng_implicit;
	double ***ng_conv_old;
	double ***ng_rhs;

	// Following coefficients are needed when QUICK is used
	double *aeW, *aeE, *aeEE;
	double *awWW, *awW, *awE;
	double *anS, *anN, *anNN;
	double *asSS, *asS, *asN;
	double *afB, *afF, *afFF;
	double *abBB, *abB, *abF;
};
typedef struct two_equation_rans Two_equation_rans;


/******************************************************************************/
/*                                    RANS                                    */
/******************************************************************************/
struct rans {

	// RANS viscosity
	double ***nut, ***nutc;

	// Two Equation model
	Two_equation_rans **two_eqn_rans;

	// Log law constants
	double kappa, b;

	double C_mu, C_muc;
	double sigma_1, sigma_2;
	double C_epsilon1, C_epsilon2, C_epsilon3;


	Strain_rate *st_rate;

	// Schumann Wall model related variables
	Wall_model *log_law;
};
typedef struct rans Rans;


/******************************************************************************/
/*                                 LAGRANGIAN                                 */
/******************************************************************************/
struct lagrangian {

	struct particle_list *p_mobile_list, *p_fixed_list, *p_release_list;

	// Temporary arrays to be used only within individual functions
	double *Temp_L;
	double *Temp_H;
	double ***ng_temp;
	double ***temp;


#ifdef POST_PROCESS
	// IBM force
	double ***ng_fx_IBM, ***ng_fy_IBM;
#endif

	// Volume fraction of particles
	double ***ng_vfu;
	double ***ng_vfv;
	double ***ng_vfw;
	double ***ng_vfc;
	double ***ng_vfz;

#ifdef VOF_SCALAR
	double ***vfu;
	double ***vfv;
	double ***vfw;

	#if defined(VOF_SMOOTH_VELO) || defined(VOF_PP_VFPRIME)
		double ***vfu_prime;
		double ***vfv_prime;
		double ***vfw_prime;
	#endif
#endif


	// Buffers for MPI communication
	struct particle  *send_part;
	struct particle  *recv_part;
	struct collision *send_coll;
	struct collision *recv_coll;
	int send_part_size, recv_part_size, send_coll_size, recv_coll_size;
	int *send_Nc, *recv_Nc;
	int send_Nc_size, recv_Nc_size;
};
typedef struct lagrangian Lagrangian;


/******************************************************************************/
/*                               PARTICLE_LIST                                */
/******************************************************************************/
struct particle_list {

	// Number of mobile and fixed particles between all processors
//	int Np_mobile;
//	int Np_fixed;
	int Np;

	// Integer used to define type of list.  Current options are (defined in
	// 'definitions.h'):
	//     MOBILE - mobile particles
	//     FIXED  - fixed particles
	int type;

	// Particle ID of first particle in list
	int ID_start;

	/*
	 Character used to identify current state of linked list.  Only meaningful
	 for first particle in linked list:
	     'l' - linked list contains only local particles
	     'f' - linked list contains only foreign particles
	     'b' - linked list contains both local and foreign particles
	     'e' - linked list contains both local and foreign particles, but only
	           near the processor's edges
	 */
//	char p_mobile_state, p_fixed_state;
	char state;

	// Linked list of Particles
	struct particle *start;
};
typedef struct particle_list Particle_list;


/******************************************************************************/
/*                                  PARTICLE                                  */
/******************************************************************************/
struct particle {

	/*------------------------------------------------------------------------*/
	/*
	 !!! Important Note:

	 Elements of this structure must be kept in blocks of the same datatype.
	 When elements are added or removed from this structure, the block counts
	 must be updated in the function Particle_initialize_MPI_datatype()
	 */
	/*------------------------------------------------------------------------*/

	//--------------------------------------------------------------------------
	// Doubles
	//--------------------------------------------------------------------------
	// Position
	double X[3], X_old[3];

	// Rotation matrix
	double Rotn[3][3], Rotn_old[3][3];

	// Translational velocity
	double U[3], U_old[3];

	// Rotational velocity
	double Omega[3], Omega_old[3];

	// Force, fluid
	double F[3];

	// Force, collisional
	double Fc[3], Fc_old[3];

	// Torque, fluid
	double T[3];

	// Torque, collisional
	double Tc[3], Tc_old[3];

	// Volume integral of momentum
	double Int_U[3], Int_U_old[3];

	// Volume integral of angular momentum
	double Int_Omega[3], Int_Omega_old[3];

	// Force and torque from immersed boundary method
	double F_IBM[3], T_IBM[3];

	// Force and torque from rigid body force
	double F_rigid[3], T_rigid[3];

	// Force and torque from collisions
	double F_coll[3], T_coll[3];


#ifdef POST_PROCESS
	// Normal and tangential collision forces, normal and tangential lubrication
	// forces
	double Fc_norm[3], Fc_tan[3], Fl_norm[3], Fl_tan[3];
	double Fc_norm_old[3], Fc_tan_old[3], Fl_norm_old[3], Fl_tan_old[3];
	// Cumulative forces (over entire timestep)
	double Fc_norm_cum[3], Fc_tan_cum[3], Fl_norm_cum[3], Fl_tan_cum[3];
#endif


	// Radius
	double R;

	// Mass divided by particle density
	double M;

	// Moment of inertia divided by particle density
	double I_p;

	double rho_s;

	double Vol_L; // Volume of marker points

	double t_part_release;


	//--------------------------------------------------------------------------
	// Integers
	//--------------------------------------------------------------------------
	// Particle identifier
	//     - unique integer differentiating this from other particles
	int ID;

	int N_L;         // Total number
	int N_L_local;  // Number of all Lagrangian markers inside the subdomain

	//--------------------------------------------------------------------------
	// Characters
	//--------------------------------------------------------------------------
	/*
	 Character used to identify current state of linked list.  Only meaningful
	 for first particle in linked list:
	     'l' - linked list contains only local particles
	     'f' - linked list contains only foreign particles
	     'b' - linked list contains both local and foreign particles
	     'e' - linked list contains both local and foreign particles, but only
	           near the processor's edges
	 */
//	char list_state;

	//--------------------------------------------------------------------------
	// Pointers
	//--------------------------------------------------------------------------
	double *X_L, *Y_L, *Z_L;
#ifdef IBM_SCALAR
	double *X_H, *Y_H, *Z_H;
#endif
#ifdef MASKING_OUTSIDE
double *PVol_L; // transient Volume per Marker point
#endif




#ifndef GRID_UNIFORM
	double *Vl;
	double *N_l_subsec;
	int *idx_subsec_low;
	int *idx_x_markers, *idx_y_markers, *idx_z_markers;
	int N_l_subsec_max;
	double **A, **B;
#endif

	int *flag_L;

	// Linked list identifier
	struct particle *next;

	// Linked list of current particle-particle collisions
	struct collision *particle_collision;

	// Linked list of current particle-wall collisions
	struct collision *wall_collision;
};
typedef struct particle Particle;


/******************************************************************************/
/*                                  COLLISION                                 */
/******************************************************************************/
struct collision {

	/*------------------------------------------------------------------------*/
	/*
	 * Both particles will have their own version of this, but the parameters
	 * should be the same.
	 *
	 * NOTE: If adding or removing elements from this structure, modify the
	 * following functions (instructions in function descriptions):
	 *     - ParticleOutput_h5_collision()
	 *     - ParticleInput_copy_collision_data()
	 */
	/*------------------------------------------------------------------------*/

	// Collision parameters
	double kn, dn;  // Normal stiffness and damping
	double kt, dt;  // Tangential stiffness and damping

	// Stokes number
	double St;

	// Tangential collision force
	double zeta_t[3], zeta_t_old[3];

	// Tangential collision parameters
	int type;

	// Unique collision partner ID: particles (ID >= 0) and walls (ID < 0)
	int other_ID;

	// Collision state for inter-processor communication
	int state;

	// Linked list identifier
	struct collision *next;
};
typedef struct collision Collision;


/******************************************************************************/
/*                                COLLISION_BAG                               */
/******************************************************************************/
struct collision_bag {

	Particle *p, *p2;
	int stage;         // What type of other ('p2') particles are we looking at?

	double n[3];       // Normal unit vector from center of 'p' to that of 'p2'
	double t[3];       // Tangential unit vector

	double g[3];       // Relative translational velocity
	double gn[3];      // Relative translational velocity - normal component
	double gt[3];      // Relative translational velocity - tangential component
	double Om_cross_R[3]; // Relative rotational velocity of two surfaces
	double gt_cp[3];   // Relative tangential velocity of contact point

	double R_cp;       // Distance from particle 'p' to contact point
	double R2_cp;      // Distance from particle 'p2' to contact point

	double g_dot_n;    // Normal component of relative translational velocity
	double gt_cp_norm; // Magnitude of contact point velocity
#ifdef ATFM
	double gt_cp_norm_old; // Magnitude of contact point velocity at older time
#endif

	double surface_distance;
	double normal_force_norm;

};
typedef struct collision_bag Collision_bag;


/******************************************************************************/
/*                                 CART3D_BAG                                 */
/******************************************************************************/
struct cart3d_bag {

	Timer *timer;
	Velocity *u, *v, *w;
	Vorticity *vor;
	Pressure *p;
	Viscosity *viscosity;
	Concentration **c;
	MAC_grid *grid;
	Parameters *params;
	Statistics2d *st2d;
	Lagrangian *lag;
	Subgrid *smag;
	Rans *rans;
	Fourier *fourier;

};
typedef struct cart3d_bag Cart3d_bag;


/******************************************************************************/
/*                               MATRIX_STRUCT                                */
/******************************************************************************/
struct matrix_struct {

	int size;           // size = 2,3,4
	double det;         // determinant
	double A[4][4];     // matrix values
	double A_inv[4][4]; // inverse matrix values

};
typedef struct matrix_struct MatrixType;


/******************************************************************************/
/*                                DEBUG_TRACE                                 */
/******************************************************************************/
struct debug_trace {
    char function[50];
    char file[50];
    int line;
    struct debug_trace *parent;
    struct debug_trace *child;
};
typedef struct debug_trace Debug_trace;

#endif // notDATATYPES_H
