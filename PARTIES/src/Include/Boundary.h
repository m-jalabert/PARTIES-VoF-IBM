#ifndef H_BOUNDARY

	#define H_BOUNDARY


/******************************************************************************/
/*                       Predefined Boundary Conditions                       */
/******************************************************************************/

#undef TEST_OUTPUT
#undef  DIRTY_FIX  // if 1 fix is one; to disable put 0


#undef PERIODIC_NOSLIP_BOX
#undef PERIODIC_FREESLIP_BOX
#undef XPERIODIC_FREESLIP_BOX
#undef PERIODIC_SHEAR_FLOW
#undef PERIODIC_DOUBLE_SHEAR_FLOW
#undef DOUBLE_PERIODIC_FREESLIP_BOX
#undef NOSLIP_BOX
#undef FREESLIP_DUCT
#undef NOSLIP_DUCT
#undef GRAVITY_CURRENT
#undef GRAVITY_CURRENT_PERIODIC
#undef TRIPLE_PERIODIC
#undef TRIPLE_PERIODIC_SHEAR
#undef LEFT_RIGHT_INFLOW_TOP_OUTFLOW
#define RISING_BUBBLE

#ifdef LEFT_RIGHT_INFLOW_TOP_OUTFLOW
        #define LEFT_INFLOW
        #define RIGHT_INFLOW
        #define ZPERIODIC
        #define BOTTOM_WALL_VELOCITY_FREESLIP
        #define TOP_WALL_VELOCITY_FREESLIP
		#define TOP_WALL_OUTFLOW

#elif defined RISING_BUBBLE
		#define ZPERIODIC						// The domain is 2D in (x,y), so we set z as periodic
		#define LEFT_WALL_VELOCITY_NOSLIP		// We want no-slip at x=0, x=2
		#define RIGHT_WALL_VELOCITY_NOSLIP		// We want no-slip at x=0, x=2
		#define BOTTOM_WALL_VELOCITY_FREESLIP  // We want slip at y=0, y=0.5
		#define TOP_WALL_VELOCITY_FREESLIP		// We want slip at y=0, y=0.5

#elif defined GRAVITY_CURRENT_PERIODIC
        #define LEFT_WALL_VELOCITY_NOSLIP
        #define RIGHT_WALL_VELOCITY_FREESLIP
        #define ZPERIODIC
        #define BOTTOM_WALL_VELOCITY_NOSLIP
        #define TOP_WALL_VELOCITY_FREESLIP


#elif defined GRAVITY_CURRENT
	#define LEFT_WALL_VELOCITY_NOSLIP
	#define RIGHT_WALL_VELOCITY_FREESLIP
	#define BACK_WALL_VELOCITY_NOSLIP
	#define FRONT_WALL_VELOCITY_NOSLIP
	#define BOTTOM_WALL_VELOCITY_NOSLIP
	#define TOP_WALL_VELOCITY_FREESLIP

#elif defined PERIODIC_NOSLIP_BOX
	#define XPERIODIC
	#define ZPERIODIC
	#define BOTTOM_WALL_VELOCITY_NOSLIP
	#define TOP_WALL_VELOCITY_NOSLIP

#elif defined PERIODIC_FREESLIP_BOX
	#define XPERIODIC
	#define ZPERIODIC
	#define BOTTOM_WALL_VELOCITY_FREESLIP
	#define TOP_WALL_VELOCITY_FREESLIP

#elif defined XPERIODIC_FREESLIP_BOX
	#define XPERIODIC
	#define BACK_WALL_VELOCITY_FREESLIP
	#define FRONT_WALL_VELOCITY_FREESLIP
	#define BOTTOM_WALL_VELOCITY_FREESLIP
	#define TOP_WALL_VELOCITY_FREESLIP

#elif defined PERIODIC_SHEAR_FLOW
	#define XPERIODIC
	#define ZPERIODIC
	#define BOTTOM_WALL_VELOCITY_NOSLIP
	#define TOP_WALL_VELOCITY

#elif defined PERIODIC_DOUBLE_SHEAR_FLOW
		#define XPERIODIC
		#define ZPERIODIC
		#define BOTTOM_WALL_VELOCITY
		#define TOP_WALL_VELOCITY

#elif defined FREESLIP_DUCT
	#define LEFT_INFLOW // any function
	#define RIGHT_OUTFLOW
	#define BACK_WALL_VELOCITY_FREESLIP
	#define FRONT_WALL_VELOCITY_FREESLIP
	#define BOTTOM_WALL_VELOCITY_FREESLIP
	#define TOP_WALL_VELOCITY_FREESLIP


#elif defined NOSLIP_DUCT
	#define LEFT_INFLOW // any function
	#define RIGHT_OUTFLOW
	#define BACK_WALL_VELOCITY_NOSLIP
	#define FRONT_WALL_VELOCITY_NOSLIP
	#define BOTTOM_WALL_VELOCITY_NOSLIP
	#define TOP_WALL_VELOCITY_NOSLIP

#elif defined NOSLIP_BOX
	#define LEFT_WALL_VELOCITY_NOSLIP
	#define RIGHT_WALL_VELOCITY_NOSLIP
	#define BACK_WALL_VELOCITY_NOSLIP
	#define FRONT_WALL_VELOCITY_NOSLIP
	#define BOTTOM_WALL_VELOCITY_NOSLIP
	#define TOP_WALL_VELOCITY_NOSLIP

#elif defined TRIPLE_PERIODIC
		#define XPERIODIC
		#define ZPERIODIC
		#define YPERIODIC

#elif defined TRIPLE_PERIODIC_SHEAR
		#define XPERIODIC
		#define ZPERIODIC
		#define YPERIODIC
		#define SHEARED_PERIODIC
#else
/******************************************************************************/
/*                                X-Boundaries                                */
/******************************************************************************/
	#undef  LEFT_WALL_VELOCITY_NOSLIP
	#undef LEFT_WALL_VELOCITY_FREESLIP
	#define  LEFT_INFLOW
	#undef  LEFT_OUTFLOW

	#undef  RIGHT_WALL_VELOCITY_NOSLIP
	#undef RIGHT_WALL_VELOCITY_FREESLIP
	#undef  RIGHT_INFLOW
	#define  RIGHT_OUTFLOW

	#undef  XPERIODIC


/******************************************************************************/
/*                                Z-Boundaries                                */
/******************************************************************************/
	#undef  BACK_WALL_VELOCITY_NOSLIP
	#undef  BACK_WALL_VELOCITY_FREESLIP

	#undef  FRONT_WALL_VELOCITY_NOSLIP
	#undef  FRONT_WALL_VELOCITY_FREESLIP

	#define ZPERIODIC


/******************************************************************************/
/*                                Y-Boundaries                                */
/******************************************************************************/
	#define BOTTOM_WALL_VELOCITY_NOSLIP
	#undef  BOTTOM_WALL_VELOCITY_FREESLIP
	#undef  BOTTOM_WALL_SCHUMANN
	#undef  BOTTOM_WALL_VELOCITY
	#define TOP_WALL_VELOCITY_NOSLIP
	#undef  TOP_WALL_VELOCITY_FREESLIP
	#undef  TOP_WALL_SCHUMANN
	#undef  TOP_WALL_VELOCITY
	#undef YPERIODIC
#endif

/******************************************************************************/
/*                                Simulation -FlowSolver                      */
/******************************************************************************/

#undef  CONSTANT_MASSFLUX  // Viscous terms solution method (default is semi-implicit FFT)
#undef FLUID_OSCILLATION        // Oscillation force acting on the fluid due to ISS-vibration (corresponds to PARTICLE_OSCILLATION)
#undef  FULLY_EXPLICIT
#define CG_SOLVE
#undef  BICG_SOLVE // nouniform mesh
#undef  TEST


/******************************************************************************/
/*                                   Output                                   */
/******************************************************************************/
#undef  OUTPUT2D


/******************************************************************************/
/* 						Volume of Fluid - PLIC with CSF					      */
/******************************************************************************/
#define VOF_PLIC



/******************************************************************************/
/*                               Concentration                                */
/******************************************************************************/
#undef CONC             // Turn on Concentration
#undef BOUSSINESQ       // Boussinesq approximation
#undef IBM_SCALAR
#undef VOF_SCALAR
#undef VOF_VELOCITY
#undef VOF_SCALAR_DEBUG
#undef VOF_SMOOTH_VELO
#undef VOF_SCALAR_NODIV // if defined we do non-conservative scheme
#undef VOF_NO_VOLUME  // the volume fraction is not used only possible with VOF_SMOOTH_VELO is enabled
#undef VOF_PROJECT
#undef VOF_PP_VFPRIME
//#undef  CONC_CLIP        // not supported Clip concentration to max 1


						/* Concentration solution method*/
#define CONC_CENTRAL     // second order advection scheme
#undef  CONC_QUICK		 // third order advection scheme does only work for a channel flow with zero gradients!!!
//#undef  CONC_BQUICK  // not supported
//#undef  CONC_FTUPWIND  // not supported && not working with periodic conditions  !!!

#define CONC_FULLY_IMPLICIT
#undef  CONC_FULLY_EXPLICIT  // the solutal diffusion is treated fully explicit .. later we may treat temperatures implicit and conc expl !!!! Not tested !!!!!
#undef  CONC_SEMI_IMPLICIT  // the solutal diffusion is treated  explicit for x&z  !!!!!! Not tested !!!!!


// defines some output quantities should be handeled differently
#undef SCALAR_DEBUG	// prints the difference at the Lagrangian points to stdout
#undef POST_PROCESS    //Ed ibm force make field in x direction


						/* Now some predefined BC for a single conc field*/

// We set mixed Boundary condition on the scalar fields
// every wall gets A\partial_x^i c+ c B = C where A,B,C are constant on each face
// constants are defined in parties.inp
#undef CONC_SINGLE_BC_NODIFF  // Zero Gradient at all faces

#undef  CONC_SINGLE_BC_RB      // Rayleigh-Benard convection
#undef  CONC_SINGLE_BC_ZERO   // c=0 at all faces
#undef SCALAR_INIT_MASK  // Scalar is set zero inside particle

#undef CONC_INDEPENDENT_BC_X // if XPERIODIC is defined, but you want to use different BCs in x for the concentration
#undef CONC_INDEPENDENT_BC_Z // if ZPERIODIC is defined,  but you want to use different BCs in z for the concentration


/// some parameter to tune convergence
#undef FAST // disabled some communication terms that may have no influence on the accuracy to be tested


/******************************************************************************/
/*                                 Particles                                  */
/******************************************************************************/
#undef LAG_PARTICLE_RESOLVED  // Turn on Lagrangian particles
#undef  PARTICLE_TRN           // Save Particle_*.h5 files for every timestep in
                               //     subfolder './trn'
#define  SUBSTEP                // Resolve particle collisions with sub-timesteps
#define STARTUP                // Prescribe velocity for particle
#undef  FORCES_DAT             // Print out 'forces.dat' - F acting on particle // error here
#define  DRY_COLLISION          // Turn off fluid forces for large St collisions
#undef ROUGH_COLLISION        // Start collision at surface roughness
#undef  LAG_MARKER_FLAG        // Turn off all competing Lag markers
#undef  LAG_MARKER_PRIORITY    // Turn off only half of competing Lag markers

#undef PARTICLE_OSCILLATION        // Oscillation force acting on the particle due to ISS-vibration (corresponds to FLUID_OSCILLATION)

#undef ONE_WAY					// Turns on one-way coupling, i.e. no feedback from the particles on the fluid
#undef TURB_FORCING			// Turns on the EP turbulent forcing

// Squirmer model
#undef  SQUIRMER_SWIMMER
#undef PARTICLE_RELEASE   // this enables the release feature
#undef SWIMMER_VERTICAL
#undef SWIMMER_HORIZONTAL


// Long-range collision models (surface distance > 0)
#define LUBRICATION_NORMAL     // Turn on lubrication forces (normal direction)
#undef  LUBRICATION_TANGENTIAL // Turn on lubrication forces (tangential direction)
#undef  ELECTROSTATIC_REPULSION // Turn on electrostatic repulsive force
// Normal force collision models
#define ACTM  // Nonlinear spring-dashpot
#undef  DEM   // Linear spring-dashpot
// Tangential force collision models
#undef  ATFM  // Set zero slip between particles (unstable for many particles)
#undef  ENABLE_ATFM_ROLLING
#define LIN_TAN  // Linear spring-dashpot
#undef COHESION // Cohesive force model

// Apply Stokes drag force for a sheared dry suspension (requires a wall
// velocity BC).  Only usable for dry simulations.
#undef  STOKES_DRAG
// Continuously moves top wall from params->ymax to params->vel_init_y0 over the
// time interval {0, time_max}.  Only usable for dry simulations.
#undef  DOWNWARD_MOVING_WALL
#define GRID_UNIFORM
#undef RETRACTION  // give a number here the Lagrangian points are shifted by xx dx inside
#undef MASKING_OUTSIDE // the ibm heating is only applied inside the particle

#undef SWIMMERS_JET
#undef SWIMMERS_SWARM

/****************************conductivity_s**************************************************/
/*                                    LES                                     */
/******************************************************************************/
#undef  LES  // Turn on LES
#undef  SMAG_DYNAMIC
#undef  FILTER_2D
#undef  LES_PLANE_AVE
#undef  LES_SPAN_AVE
#undef  LES_LAG_AVE
#undef  CONC_DYNAMIC
#undef  IMMERSED_BOUNDARY
#undef  THERMAL_KADER


/******************************************************************************/
/*                                    RANS                                    */
/******************************************************************************/
#undef  RANS  // Turn on RANS
#undef  TWO_EQUATION_MODEL
#undef RANS_CENTRAL
#undef  RANS_QUICK
#undef  RANS_BQUICK
#undef  RANS_FTUPWIND


/******************************************************************************/
/*                                Logic checks                                */
/******************************************************************************/

//------------------------------------------------------------------------------
// Oscillation
//------------------------------------------------------------------------------

#if defined FLUID_OSCILLATION || defined PARTICLE_OSCILLATION

	#if defined FLUID_OSCILLATION && defined PARTICLE_OSCILLATION
        #error Both oscillations (acting on fluid and on particle) are specified
	#endif // defined

#endif


//------------------------------------------------------------------------------
// Solution method
//------------------------------------------------------------------------------
#if defined CG_SOLVE || defined BICG_SOLVE

	#if defined CG_SOLVE && defined BICG_SOLVE
		#error Incompatible solvers 'CG_SOLVE' and 'BICG_SOLVE'
	#endif

	#ifdef FULLY_EXPLICIT
		#error Both explicit and implicit solvers specified
	#endif

	#define FULLY_IMPLICIT
#endif


//------------------------------------------------------------------------------
// Variable viscosity
//------------------------------------------------------------------------------
#if defined LES || defined RANS
	#define VAR_VISC
#endif


//------------------------------------------------------------------------------
// Concentration
//------------------------------------------------------------------------------
#ifdef CONC

	#if defined CONC_BQUICK && !defined CONC_QUICK
		#warning 'CONC_BQUICK' defined, but not 'CONC_QUICK'. Defining 'CONC_QUICK'...
		#define CONC_QUICK
	#endif

	#if defined CONC_CENTRAL && defined CONC_QUICK
		#error Incompatible concentration convective methods 'CONC_CENTRAL' and 'CONC_QUICK'
	#elif defined CONC_CENTRAL && defined CONC_FTUPWIND
		#error Incompatible concentration convective methods 'CONC_CENTRAL' and 'CONC_FTUPWIND'
	#elif defined CONC_QUICK && defined CONC_FTUPWIND
		#error Incompatible concentration convective methods 'CONC_QUICK' and 'CONC_FTUPWIND'
	#endif
	#if !defined CONC_CENTRAL && !defined CONC_QUICK && !defined CONC_FTUPWIND
		#error No concentration convective method specified
	#endif

	#ifdef CONC_FULLY_IMPLICIT
    	#define  CONC_SOLVE_CG
	#endif

#else // not CONC
	#undef CONC_DYNAMIC
	#undef CONC_FTUPWIND
	#undef CONC_CENTRAL
	#undef CONC_BQUICK
	#undef CONC_QUICK
	#undef BOUSSINESQ
	#undef CONC_CLIP
	#undef IBM_SCALAR
	#undef VOF_SCALAR
	#undef VOF_SCALAR_DEBUG
#endif // not CONC
#ifndef VOF_SCALAR
	#undef VOF_SCALAR_DEBUG
#endif
//------------------------------------------------------------------------------
// Particles
//------------------------------------------------------------------------------
#if defined ACTM && defined DEM
	#error Incompatible collision models 'ACTM' and 'DEM'
#endif
#if !defined ACTM && !defined DEM
	#error Specify a normal collision model: 'ACTM' or 'DEM'
#endif
#if defined ATFM && defined LIN_TAN
	#error Incompatible collision models 'ATFM' and 'LIN_TAN'
#endif
#if defined ACTM || defined DEM || defined ATFM
	#define STORE_COLLISION
#endif
#ifndef ATFM
	#undef ENABLE_ATFM_ROLLING
#endif
#if defined STOKES_DRAG && !defined TOP_WALL_VELOCITY
	#error 'STOKES_DRAG' requires 'TOP_WALL_VELOCITY' boundary condition
#endif
#if defined STOKES_DRAG && !defined DRY_PARTICLES
	#error 'STOKES_DRAG' only valid for dry simulations
#endif
#if defined DOWNWARD_MOVING_WALL && !defined DRY_PARTICLES
	#error 'DOWNWARD_MOVING_WALL' only valid for dry simulations
#endif
#if defined LAG_MARKER_FLAG && defined LAG_MARKER_PRIORITY
	#error Incompatible options 'LAG_MARKER_FLAG' and 'LAG_MARKER_PRIORITY'
#endif

#if defined LAG_PARTICLE_RESOLVED && !defined GRID_UNIFORM
#error 'The resolved particles only work with a uniform Grid, i.e, 'GRID_UNIFORM' enabled.'
#endif

#if defined SQUIRMER_SWIMMER && defined VOF_SCALAR
	#if !defined VOF_SCALAR_U_FLUID
		#warning 'SQUIRMER_SWIMMER' and 'VOF_SCALAR' defined but not 'VOF_SCALAR_U_FLUID'. Defining 'VOF_SCALAR_U_FLUID'
		#define VOF_SCALAR_U_FLUID
	#endif
#endif
//------------------------------------------------------------------------------
// LES
//------------------------------------------------------------------------------
#ifndef LES
	#undef SMAG_DYNAMIC
	#undef FILTER_2D
	#undef LES_PLANE_AVE
	#undef LES_SPAN_AVE
	#undef LES_LAG_AVE
	#undef CONC_DYNAMIC
#endif


//------------------------------------------------------------------------------
// RANS
//------------------------------------------------------------------------------
#ifdef RANS

	#ifdef LES
		#error Incompatible turbulence models 'LES' and 'RANS'
	#endif

	#if defined RANS_BQUICK && !defined RANS_QUICK
		#warning 'RANS_BQUICK' defined, but not 'RANS_QUICK'. Defining 'RANS_QUICK'...
		#define RANS_QUICK
	#endif

	#if defined RANS_CENTRAL && defined RANS_QUICK
		#error Incompatible RANS convective methods 'RANS_CENTRAL' and 'RANS_QUICK'
	#elif defined RANS_CENTRAL && defined RANS_FTUPWIND
		#error Incompatible RANS convective methods 'RANS_CENTRAL' and 'RANS_FTUPWIND'
	#elif defined RANS_QUICK && defined RANS_FTUPWIND
		#error Incompatible RANS convective methods 'RANS_QUICK' and 'RANS_FTUPWIND'
	#endif
	#if !defined RANS_CENTRAL && !defined RANS_QUICK && !defined RANS_FTUPWIND
		#error No RANS convective method specified
	#endif

#else // not RANS
	#undef TWO_EQUATION_MODEL
	#undef RANS_CENTRAL
	#undef RANS_QUICK
	#undef RANS_BQUICK
	#undef RANS_FTUPWIND
#endif // not RANS


//------------------------------------------------------------------------------
// X-Boundaries
//------------------------------------------------------------------------------
#ifdef XPERIODIC
	#ifdef LEFT_INFLOW
		#error Incompatible boundary conditions 'XPERIODIC' and 'LEFT_INFLOW'
	#elif  RIGHT_INFLOW
		#error Incompatible boundary conditions 'XPERIODIC' and 'RIGHT_INFLOW'
	#elif defined LEFT_OUTFLOW
		#error Incompatible boundary conditions 'XPERIODIC' and 'LEFT_OUTFLOW'
	#elif defined RIGHT_OUTFLOW
		#error Incompatible boundary conditions 'XPERIODIC' and 'RIGHT_OUTFLOW'
	#elif defined LEFT_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'XPERIODIC' and 'LEFT_WALL_VELOCITY_NOSLIP'
	#elif defined RIGHT_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'XPERIODIC' and 'RIGHT_WALL_VELOCITY_NOSLIP'
	#elif defined LEFT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'XPERIODIC' and 'LEFT_WALL_VELOCITY_FREESLIP'
	#elif defined RIGHT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'XPERIODIC' and 'RIGHT_WALL_VELOCITY_FREESLIP'
	#endif
#else
	//-----------
	// Left Wall
	//-----------
	#if defined LEFT_INFLOW && defined LEFT_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'LEFT_INFLOW' and 'LEFT_WALL_VELOCITY_NOSLIP'
	#elif defined LEFT_INFLOW && defined LEFT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'LEFT_INFLOW' and 'LEFT_WALL_VELOCITY_FREESLIP'
	#elif defined LEFT_INFLOW && defined LEFT_OUTFLOW
		#error Incompatible boundary conditions 'LEFT_INFLOW' and 'LEFT_OUTFLOW'
	#elif defined LEFT_OUTFLOW && defined LEFT_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'LEFT_OUTFLOW' and 'LEFT_WALL_VELOCITY_NOSLIP'
	#elif defined LEFT_OUTFLOW && defined LEFT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'LEFT_OUTFLOW' and 'LEFT_WALL_VELOCITY_FREESLIP'
	#elif defined LEFT_WALL_VELOCITY_NOSLIP && defined LEFT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'LEFT_WALL_VELOCITY_NOSLIP' and 'LEFT_WALL_VELOCITY_FREESLIP'
	#endif
	#if !defined LEFT_INFLOW && !defined LEFT_OUTFLOW && !defined LEFT_WALL_VELOCITY_NOSLIP && !defined LEFT_WALL_VELOCITY_FREESLIP
		#error No left wall boundary condition specified
	#endif

	//------------
	// Right Wall
	//------------
	#if defined RIGHT_INFLOW && defined RIGHT_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'RIGHT_INFLOW' and 'RIGHT_WALL_VELOCITY_NOSLIP'
	#elif defined RIGHT_INFLOW && defined RIGHT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'RIGHT_INFLOW' and 'RIGHT_WALL_VELOCITY_FREESLIP'
	#elif defined RIGHT_INFLOW && defined RIGHT_OUTFLOW
		#error Incompatible boundary conditions 'RIGHT_INFLOW' and 'RIGHT_OUTFLOW'
	#elif defined RIGHT_OUTFLOW && defined RIGHT_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'RIGHT_OUTFLOW' and 'RIGHT_WALL_VELOCITY_NOSLIP'
	#elif defined RIGHT_OUTFLOW && defined RIGHT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'RIGHT_OUTFLOW' and 'RIGHT_WALL_VELOCITY_FREESLIP'
	#elif defined RIGHT_WALL_VELOCITY_NOSLIP && defined RIGHT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'RIGHT_WALL_VELOCITY_NOSLIP' and 'RIGHT_WALL_VELOCITY_FREESLIP'
	#endif
	#if !defined RIGHT_INFLOW && !defined RIGHT_OUTFLOW && !defined RIGHT_WALL_VELOCITY_NOSLIP && !defined RIGHT_WALL_VELOCITY_FREESLIP
		#error No right wall boundary condition specified
	#endif

	#if defined LEFT_INFLOW && defined RIGHT_INFLOW
		#if !defined TOP_WALL_OUTFLOW
		#error Incompatible boundary conditions 'LEFT_INFLOW' and 'RIGHT_INFLOW'
		#endif
	#endif
#endif


//------------------------------------------------------------------------------
// Y-Boundaries
//------------------------------------------------------------------------------
#if defined  BOTTOM_WALL_SCHUMANN || defined TOP_WALL_SCHUMANN
	#define SCHUMANN
#endif

//-------------
// Bottom Wall
//-------------
#if defined BOTTOM_WALL_SCHUMANN && defined BOTTOM_WALL_VELOCITY_NOSLIP
	#error Incompatible boundary conditions 'BOTTOM_WALL_SCHUMANN' and 'BOTTOM_WALL_VELOCITY_NOSLIP'
#elif defined BOTTOM_WALL_SCHUMANN && defined BOTTOM_WALL_VELOCITY_FREESLIP
	#error Incompatible boundary conditions 'BOTTOM_WALL_SCHUMANN' and 'BOTTOM_WALL_VELOCITY_FREESLIP'
#elif defined BOTTOM_WALL_VELOCITY_NOSLIP && defined BOTTOM_WALL_VELOCITY_FREESLIP
	#error Incompatible boundary conditions 'BOTTOM_WALL_VELOCITY_NOSLIP' and 'BOTTOM_WALL_VELOCITY_FREESLIP'
#elif defined BOTTOM_WALL_VELOCITY && defined BOTTOM_WALL_VELOCITY_NOSLIP
	#error Incompatible boundary conditions 'BOTTOM_WALL_VELOCITY' and 'BOTTOM_WALL_VELOCITY_NOSLIP'
#elif defined BOTTOM_WALL_VELOCITY && defined BOTTOM_WALL_VELOCITY_FREESLIP
	#error Incompatible boundary conditions 'BOTTOM_WALL_VELOCITY' and 'BOTTOM_WALL_VELOCITY_FREESLIP'
#endif
#if defined BOTTOM_WALL_VELOCITY && !defined TOP_WALL_VELOCITY
	#error Bottom wall velocity boundary condition also requires 'TOP_WALL_VELOCITY'
#endif
#if !defined BOTTOM_WALL_SCHUMANN && !defined BOTTOM_WALL_VELOCITY_NOSLIP && !defined BOTTOM_WALL_VELOCITY_FREESLIP && !defined BOTTOM_WALL_VELOCITY && !defined YPERIODIC
	#error No bottom wall boundary condition specified
#endif

//----------
// Top Wall
//----------
#if defined TOP_WALL_SCHUMANN && defined TOP_WALL_VELOCITY_NOSLIP
	#error Incompatible boundary conditions 'TOP_WALL_SCHUMANN' and 'TOP_WALL_VELOCITY_NOSLIP'
#elif defined TOP_WALL_SCHUMANN && defined TOP_WALL_VELOCITY_FREESLIP
	#error Incompatible boundary conditions 'TOP_WALL_SCHUMANN' and 'TOP_WALL_VELOCITY_FREESLIP'
#elif defined TOP_WALL_VELOCITY_NOSLIP && defined TOP_WALL_VELOCITY_FREESLIP
	#error Incompatible boundary conditions 'TOP_WALL_VELOCITY_NOSLIP' and 'TOP_WALL_VELOCITY_FREESLIP'
#elif defined TOP_WALL_VELOCITY && defined TOP_WALL_SCHUMANN
	#error Incompatible boundary conditions 'TOP_WALL_VELOCITY' and 'TOP_WALL_SCHUMANN'
#elif defined TOP_WALL_VELOCITY && defined TOP_WALL_VELOCITY_NOSLIP
	#error Incompatible boundary conditions 'TOP_WALL_VELOCITY' and 'TOP_WALL_VELOCITY_NOSLIP'
#elif defined TOP_WALL_VELOCITY && defined TOP_WALL_VELOCITY_FREESLIP
	#error Incompatible boundary conditions 'TOP_WALL_VELOCITY' and 'TOP_WALL_VELOCITY_FREESLIP'
#endif
#if !defined TOP_WALL_SCHUMANN && !defined TOP_WALL_VELOCITY_NOSLIP && !defined TOP_WALL_VELOCITY_FREESLIP && !defined TOP_WALL_VELOCITY && !defined YPERIODIC
	#error No top wall boundary condition specified
#endif


//------------------------------------------------------------------------------
// Z-Boundaries
//------------------------------------------------------------------------------
#ifdef ZPERIODIC
	#ifdef BACK_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'ZPERIODIC' and 'BACK_WALL_VELOCITY_NOSLIP'
	#elif defined FRONT_WALL_VELOCITY_NOSLIP
		#error Incompatible boundary conditions 'ZPERIODIC' and 'FRONT_WALL_VELOCITY_NOSLIP'
	#elif defined BACK_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'ZPERIODIC' and 'BACK_WALL_VELOCITY_FREESLIP'
	#elif defined FRONT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'ZPERIODIC' and 'FRONT_WALL_VELOCITY_FREESLIP'
	#endif
#else
	//-----------
	// Back Wall
	//-----------
	#if defined BACK_WALL_VELOCITY_NOSLIP && defined BACK_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'BACK_WALL_VELOCITY_NOSLIP' and 'BACK_WALL_VELOCITY_FREESLIP'
	#endif
	#if !defined BACK_WALL_VELOCITY_NOSLIP && !defined BACK_WALL_VELOCITY_FREESLIP
		#error No back wall boundary condition specified
	#endif

	//------------
	// Front Wall
	//------------
	#if defined FRONT_WALL_VELOCITY_NOSLIP && defined FRONT_WALL_VELOCITY_FREESLIP
		#error Incompatible boundary conditions 'FRONT_WALL_VELOCITY_NOSLIP' and 'FRONT_WALL_VELOCITY_FREESLIP'
	#endif
	#if !defined FRONT_WALL_VELOCITY_NOSLIP && !defined FRONT_WALL_VELOCITY_FREESLIP
		#error No front wall boundary condition specified
	#endif
#endif


	// Check Time stepping scheme for scalar
#ifdef CONC
	#if defined CONC_FULLY_EXPLICIT	 && defined CONC_FULLY_IMPLICIT
		#error Incompatible handling of diffusion 'CONC_FULLY_EXPLICIT' and 'CONC_FULLY_IMPLICIT'
	#endif
	#if defined CONC_FULLY_EXPLICIT && defined CONC_SEMI_IMPLICIT
		#error Incompatible handling of diffusion 'CONC_FULLY_EXPLICIT' and 'CONC_SEMI_IMPLICIT'
	#endif
	#if defined CONC_FULLY_IMPLICIT && defined CONC_SEMI_IMPLICIT
		#error Incompatible handling of diffusion 'CONC_FULLY_EXPLICIT' and 'CONC_SEMI_IMPLICIT'
	#endif
	#if (defined CONC_FULLY_IMPLICIT +defined CONC_SEMI_IMPLICIT + defined CONC_FULLY_EXPLICIT) == 0
		#error One scheme for solutal diffusion has to be choosen
	#endif
	#if defined CONC_FULLY_EXPLICIT	 && defined RANS
		#error RANS and the fully implicit scheme for conc are not implemented
	#endif
	#ifdef CONC_FULLY_IMPLICIT
		#if !defined CONC_SOLVE_CG
		#error  The fully implicit conc scheme is onyl implemented for the CG solver
	#endif
	#endif
	#if  defined CONC_SEMI_IMPLICIT
		#error The semi implicit solver is not yet supported
	#endif
	#if   defined CONC_LEFT_INFLOW
		#error 'CONC_LEFT_INFLOW' has no influence
	#endif
	#if defined IBM_SCALAR  && !defined LAG_PARTICLE_RESOLVED
		#error 'LAG_PARTICLE_RESOLVED' has to be defined for 'IBM_SCALAR'
	#endif
	#if  defined CONC_BQUICK && defined IBM_SCALAR
		#error 'CONC_BQUICK' and 'IBM_SCALAR' do not work togehter
	#endif
	#if defined IBM_SCALAR  && !defined LAG_PARTICLE_RESOLVED
		#error 'LAG_PARTICLE_RESOLVED' has to be defined for 'IBM_SACLAR'
	#endif
	#if defined IBM_SCALAR  &&  !defined CONC_FULLY_IMPLICIT
		#error For 'IBM_SCALAR' 'CONC_FULLY_IMPLICIT' has to be enabled
	#endif
	#if defined IBM_SCALAR && defined VOF_SCALAR
		#error Not both 'IBM_SCALAR' and 'VOF_SCALAR' can be used at the same time
	#endif
	#if (defined VOF_SCALAR && ! defined LAG_PARTICLE_RESOLVED)
		#error If  'VOF_SCALAR' then  'LAG_PARTICLE_RESOLVED' must be enabled
	#endif
	#if (!defined VOF_SCALAR &&  defined VOF_VELOCITY)
		#error If  'VOF_SCALAR' then  'LAG_PARTICLE_RESOLVED' must be enabled
	#endif
	#ifdef CONC_CLIP
		#error  'CONC_CLIP' is not implemented
	#endif
//	#ifdef CONC_QUICK
//		#error  'CONC_QUICK' is not implemented correctly
//	#endif



#endif

#if !(defined CONC && defined LAG_PARTICLE_RESOLVED)
	#undef SCALAR_INIT_MASK
#endif

#ifdef ONE_WAY
	#ifdef BOUSSINESQ
		#error No bouyancy with one-way coupling
	#endif
	#ifdef VOF_SCALAR
		#error No impact on the scalar diffusion by particles
	#endif
#endif

#ifndef LAG_PARTICLE_RESOLVED
	#ifdef ONE_WAY
		#error The particle module has to be switched on for one-way coupling
	#endif
#endif

#ifndef LAG_PARTICLE_RESOLVED
	#undef VOF_SCALAR
#endif

#if  defined VOF_NO_VOLUME  && !defined VOF_SMOOTH_VELO
	#error 'vof_no_velo does not work without vof_smooth_velo'
#endif

#endif // notH_BOUNDARY