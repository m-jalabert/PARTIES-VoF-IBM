#ifndef DEFINITIONS_H
	#define DEFINITIONS_H

#include "mpi.h"
#include "Boundary.h"

extern MPI_Comm comm3d;

#ifndef PI
	#define PI 3.14159265358979323846
#endif

#define HALF 0.5

#define NO  0
#define YES 1

#define LOCK_EXCHANGE 1

#define CONSTANT_VISCOSITY 0
#define VARIABLE_VISCOSITY 1

#define OPENDX   0
#define TECPLOT  1
#define PARAVIEW 2

/* For outflow velocity */
#define CONVECTIVE_CONSTANT_VELOCITY 0
#define CONVECTIVE_MAXIMUM_VELOCITY 1


/* Memory allocation */
#define GVG_CHAR      0
#define GVG_SHORT_INT 1
#define GVG_INT       2
#define GVG_LONG_INT  3
#define GVG_FLOAT     4
#define GVG_DOUBLE    5


#define maximum(a,b) ((a)<(b)?(b):(a))
#define minimum(a,b) ((a)>(b)?(b):(a))

#define max(a,b) ((a)<(b)?(b):(a))
#define min(a,b) ((a)>(b)?(b):(a))

#ifndef sign
	#define sign(a) ( ((a) > 0.0) ? 1 : ( (a) < 0.0 ? -1 : 0) )
#endif

// Input files
#define PARTIES_INPUT_FILE "parties.inp"
#define P_MOBILE_INPUT_FILE "p_mobile.inp"
#define P_FIXED_INPUT_FILE "p_fixed.inp"
#define P_RELEASE_INPUT_FILE "p_release.inp"

// Radius of delta function used for particle IBM (in grid units)
#define DELTA_FUNC_RADIUS 1.5

// Distance flagged Lagrangian markers need to be to be turned off (in grid units)
#define LAG_FLAG_RANGE 3.0

// Defines some parameter for the scalar VOF method
#define VOF_HARMONIC_MEAN_COND 1 // if == 1 we compute the conductivity(diffusivity resp) by a harmonic mean

/* Stencil of finite difference scheme */
#define STENCIL 7

//#define PCW MPI_COMM_WORLD
#define PCW comm3d

/* Root processor id */
#define MASTER 0

/* output data type */
#define BINARY 1
#define ASCII  2

#define REDUCE_TO_MASTER 0
#define REDUCE_TO_ALL 1

/* states of the nodes */
#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3
#define BACK  4
#define FRONT 5
#define OTHER 6

/* Error analysis */
#define NORM1    1
#define NORM2    2
#define NORM_INF 3

/* Grid status */
#define SOLID    2
#define IMMERSED 3
#define FLUID    4
#define BOUNDARY 5

/*Which Domain*/

#define FULL_DOMAIN 0
#define FLUID_DOMAIN 1
#define MOBILE_PARTICLE_DOMAIN 2
#define FIXED_PARTICLE_DOMAIN	3
#define PARTICLE_DOMAIN 4

/* Number of fluid nodes used in the interpolation stencil */
#define IBM_MAX 3

/* Central-Field-types */
#define CENTRAL_FULL 1
#define CENTRAL_PERTURBATION 0


/* Boundary conditions */
#define DIRICHLET 0
#define NEUMANN   1
#define ROBIN     2

// Velocity initialization options
#define VEL_INIT_ZERO       0
#define VEL_INIT_UNIFORM    1
#define VEL_INIT_LINEAR     2
#define VEL_INIT_POISEUILLE 3
#define VEL_INIT_ROT_SHEAR  4
#define VEL_INIT_PRECURSOR  5
#define VEL_INIT_TGV  6
#define VEL_INIT_JET 7
#define VEL_INIT_JET_3D 8
#define VEL_INIT_LEFT_RIGHT 9
#define VEL_INIT_NOISE 10

// Startup initialization options
#define STUP_INIT_TIME			0
#define STUP_INIT_GOND_ST_27	1
#define STUP_INIT_GOND_ST_152	2

// Concentration initialization options
#define CONC_INIT_TWOLAYERS 0
#define CONC_INIT_LINGRAD 1
#define CONC_INIT_ERF 2
#define GAUSSIAN_BUMP_X 3

// Velocity type, used for setting boundary conditions
#define VEL_TYPE_NORMAL 0
#define VEL_TYPE_CG     1

// Particle information
#define MOBILE  0
#define FIXED   1
#define LOCAL   2
#define FOREIGN 3
#define RELEASE 4

// Particle collision stage information
#define COLL_STAGE_MOBILE  10
#define COLL_STAGE_FIXED   11
#define COLL_STAGE_FOREIGN 12
#define COLL_STAGE_WALL    13

// Options for specifying state of particle linked list
#define LIST_STATE_LOCAL   20
#define LIST_STATE_FOREIGN 21
#define LIST_STATE_BOTH    22
#define LIST_STATE_EDGE    23

// Options for function Lagrangian_collect_forces()
#define LAG_COLLECT_ALL   30
#define LAG_COLLECT_COLL  31
#define LAG_COLLECT_HYDRO 32
#define LAG_COLLECT_NONE  33

// IDs of walls for collisions
#define WALL_ID_XMIN -1
#define WALL_ID_XMAX -2
#define WALL_ID_YMIN -3
#define WALL_ID_YMAX -4
#define WALL_ID_ZMIN -5
#define WALL_ID_ZMAX -6

// Tangential collision types
#define SLIDING 40
#define ROLLING 41

// Collision state flags for passing collisions between processors
#define COLL_STATE_NEUTRAL 50
#define COLL_STATE_OWNER   51
#define COLL_STATE_DESTROY 52

// Coefficients used in RK3 Scheme
#define GAMMA 8.0 / 15.0, 5.0 / 12.0, 3.0 / 4.0  //Explicit part at sub-step level k-1
#define ZETA  0.0, -17.0 / 60.0, -5.0 / 12.0     //Explicit part at sub-step level k-2
#define BETA 4.0 / 15.0, 1.0 / 15.0, 1.0 / 6.0   //Implicit part (at "k-1" and "k")
#define BETA2 8.0 / 15.0, 2.0 / 15.0, 1.0 / 3.0   //Twice the Implicit part (at "k-1" and "k")

#define GAMBETA 2.0, 25.0/4.0, 9.0/2.0 //GAMMA/BETA
#define ZETBETA 0.0, -17.0/4.0, -2.5   //ZETA/BETA

// Useful for short loops
#define FORI3 for (i=0;i<3;i++)


// Set array of doubles to zero
#define DSET_ZERO(array, length) memset(array, 0, (length) * sizeof(double))

// Dot product of two 3-D vectors
#define DOT(a,b) ( a[0] * b[0] + a[1] * b[1] + a[2] * b[2] )

// Debugging trace function call
#define DTRACE(myfun) Display_create_debug_trace_child(dtrace, myfun, __FILE__, __LINE__)

// SCALAR IBM
#define PARTICLE_BC_A 1
#define PARTICLE_BC_B 0
#define PARTICLE_BC_C 0

//
#define CAPACITY(a,b) ( a[0] * b[0] + a[1] * b[1] + a[2] * b[2] )


#define CONCENTRATION 'a'
#define CONCENTRATION_PERTURBATION 'b'

#define U_VELOCITY 'd'
#define U_VELOCITY_PERTURBATION 'e'

#define V_VELOCITY 'f'
#define W_VELOCITY 'g'

#define P_CORRECTION 'h'

// deprecated stuff maybe ready to be deleted
#undef TEST_OUTPUT
#undef  DIRTY_FIX  // if 1 fix is one; to disable put 0
/// some parameter to tune convergence
#undef FAST // disabled some communication terms that may have no influence on the accuracy to be tested



// Momentum solver
#ifdef FULLY_EXPLICIT
	#define Velocity_solve Velocity_solve_explicit
#elif defined CG_SOLVE
	#define Velocity_solve Velocity_solve_cg
#elif defined BICG_SOLVE
	#define Velocity_solve Velocity_solve_bicg
#else
	#define Velocity_solve Velocity_solve_semi_implicit
#endif

// Pressure solver
#if defined XPERIODIC && defined ZPERIODIC
	#define Pressure_solve Pressure_solve_xz_period
#elif defined XPERIODIC && !defined ZPERIODIC
	#define Pressure_solve Pressure_solve_x_period
#elif !defined XPERIODIC && defined ZPERIODIC
	#define Pressure_solve Pressure_solve_z_period
#else
	#define Pressure_solve Pressure_solve_no_period
#endif

// Concentration
#ifdef CONC
		#ifdef VOF_SCALAR
				#define Conc_set_conv_viscous Conc_set_conv_viscous_central_mixed_vof
		#else
			#if defined CONC_QUICK
				#define Conc_set_conv_viscous Conc_set_conv_viscous_quick_mixed
			#elif defined CONC_CENTRAL
				#define	Conc_set_conv_viscous Conc_set_conv_viscous_central_mixed
			#else
				#error 'Choose advection scheme for scalar'
			#endif
		#endif

#else
		#define 	Conc_set_conv_viscous Conc_set_conv_viscous_central_mixed


#endif




// Mergesort for particle linked list
#define MS_STRUCTURE Particle
#define MS_DATA ID

// Particles on non-uniform grid
#ifdef GRID_UNIFORM
	#define Lagrangian_generate_points Lagrangian_generate_points_Leopardi
	#define Interpolate_Lag_to_Eul Interpolate_Lag_to_Eul_uniform
#else
	#define Lagrangian_generate_points Lagrangian_generate_points_Bala
	#define Interpolate_Lag_to_Eul Interpolate_Lag_to_Eul_nonuniform
#endif

// RANS
#ifdef RANS_QUICK
	#define Rans_set_conv_viscous Rans_set_conv_viscous_quick
#elif defined RANS_FTUPWIND
	#define Rans_set_conv_viscous Rans_set_conv_viscous_ftupwind
#else  // RANS_CENTRAL
	#define Rans_set_conv_viscous Rans_set_conv_viscous_central
#endif

#ifdef XPERIODIC
 	#ifndef CONC_INDEPENDENT_BC_X
		#define XPERIODIC_CONC
	#endif
#endif

#ifdef ZPERIODIC
 	#ifndef CONC_INDEPENDENT_BC_Z
		#define ZPERIODIC_CONC
	#endif
#endif

#define HERE printf("!!!!!HERE!!!!!\n"); fflush(stdout);



#endif // DEFINITIONS_H
