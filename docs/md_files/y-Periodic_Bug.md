# y-Periodic Bug

This bug seems to make the code crash, when the YPERIODIC boundary condition is selected for the y-boundary. If the boundary condition TRIPLE_Periodic is selected (which applies the periodic boundary condition to all three boundaries), the code will work just fine. Just, when the periodic boundary condition is selected individually for the y-boundary and other boundary conditions, such as free-slip or no-slip, are selected for the other bodaunries, the simulation will crash. 

Simulation will stop within the first iterations with the following error code:
<pre>
***************************************************************************
!!! Error called by Display_throw_error (Temporal_int.c:1052)           !!!
!!!              at Temporal_int_all_the_equations (Temporal_int.c:473) !!!
!!!              at Temporal_int_rk3 (Cart3d.c:446)                     !!!
!!! Poisson solver did not converge                                     !!!
***************************************************************************
</pre>

This behaviour is likely to be caused by the fact, the initially, the code could not run a simulation with the YPERIODIC boundary condition and the TRIPLE_PERDIODIC boundary condition was added later as a new feature. Within the Display.c-file, which manages among other things, the output during the simulation, no option is implemented for the YPERIODIC boundary condition, but an else-statement that will put "SCHUHMANN" if none of the implemented option is chosen. Because of that, the output will be wrong when the YPERIODIC boundary condition is chosen:

<pre>
	//--------------------------------------------------------------------------
	// Y-Boundaries
	//--------------------------------------------------------------------------
#if defined TOP_WALL_VELOCITY && defined BOTTOM_WALL_VELOCITY
	double wall_vel = params->ubulk_target;
#else
	double wall_vel = 2.0 * params->ubulk_target;
#endif

#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
	printf("Bottom wall: .... no-slip\n");
#elif defined BOTTOM_WALL_VELOCITY_FREESLIP
	printf("Bottom wall: .... free-slip\n");
#elif defined BOTTOM_WALL_VELOCITY
	printf("Bottom wall: .... U = %g\n", -wall_vel);
#else  // SCHUMANN
	printf("Bottom wall: .... Schumann\n");
#endif
#ifdef TOP_WALL_VELOCITY_NOSLIP
	printf("Top wall: ....... no-slip\n");
#elif defined TOP_WALL_VELOCITY_FREESLIP
	printf("Top wall: ....... free-slip\n");
#elif defined TOP_WALL_SCHUMANN
	printf("Top wall: ....... Schumann\n");
#elif defined TOP_WALL_VELOCITY
	printf("Top wall: ....... U = %g\n", wall_vel);
#endif
</pre>

Within the Conc.c file, solely XPERIODIC and ZPERIODIC are specified.

The search for "XPERIODIC" leads to 186 results in 30 files (within PARTIES).<br>
The search for "ZPERIODIC" leads to 159 results in 30 files (within PARTIES).<br>
The search for <code>"YPERIODIC"</code> leads to 57 results in 17 files (within PARTIES).<br>

Indicating that the YPERIODIC boundary condition is not fully implemented but solely so far that it will work for the TRIPLE_PERIODIC boundary condition.
