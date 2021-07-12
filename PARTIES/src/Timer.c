#include <stdio.h>
#include <stdlib.h>

#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Timer.h"

extern MPI_Comm comm3d;

Timer *Timer_create(MAC_grid *grid, Parameters *params) {

	Timer *timer = (Timer *)malloc(sizeof(Timer));

	timer -> Wtime_total = 0.0;
	timer -> Wtime_init = 0.0;
	timer -> Wtime_intEOM = 0.0;
	timer -> Wtime_output = 0.0;
	timer -> Wtime_output_2d = 0.0;
	timer -> Wtime_cfl = 0.0;

	// Communication
	timer -> Wtime_comm_3D = 0.0;
	timer -> Wtime_comm_2D_reduce = 0.0;

	// Velocity
	timer -> Wtime_vel_convective = 0.0;
	timer -> Wtime_vel_rhs = 0.0;
	timer -> Wtime_vel_solve = 0.0;
	timer -> Wtime_vel_boundaries = 0.0;
	timer -> Wtime_vel_cell_center = 0.0;

	// Pressure
	timer -> Wtime_p_rhs = 0.0;
	timer -> Wtime_p_solve = 0.0;
	timer -> Wtime_p_project = 0.0;
	timer -> Wtime_p_divergence = 0.0;

	// Concentration
	timer -> Wtime_c_total = 0.0;
	timer -> Wtime_c_convective = 0.0;
	timer -> Wtime_c_rhs = 0.0;
	timer -> Wtime_c_solve = 0.0;
	timer -> Wtime_c_outofbounds = 0.0;

	// Subgrid
	timer -> Wtime_sgs_total = 0.0;
	timer -> Wtime_sgs_filter = 0.0;
	timer -> Wtime_sgs_dynamic = 0.0;

	// RANS
	timer -> Wtime_rans_total = 0.0;
	timer -> Wtime_rans_convective = 0.0;
	timer -> Wtime_rans_rhs = 0.0;
	timer -> Wtime_rans_solve = 0.0;

	// Particle
	timer -> Wtime_particle_total = 0.0;
	timer -> Wtime_particle_comm = 0.0;
	timer -> Wtime_particle_coll = 0.0;
	timer -> Wtime_particle_forc = 0.0;
	timer -> Wtime_particle_int = 0.0;

	return timer;
}

void Timer_destroy(Timer *timer, MAC_grid *grid, Parameters *params) {

	free(timer);
}
