#ifndef VELOCITY_H
 #define VELOCITY_H

#include "definitions.h"
#include "DataTypes.h"

Velocity *Velocity_create(MAC_grid *grid, Parameters *params, char which_velocity,
		Debug_trace *dtrace);
void Velocity_destroy(Velocity *vel, MAC_grid *grid, Parameters *params);
void Velocity_cell_center(Cart3d_bag *data_bag);
void Velocity_u_set_implicit_explicit(Cart3d_bag *data_bag);
void Velocity_v_set_implicit_explicit(Cart3d_bag *data_bag);
void Velocity_w_set_implicit_explicit(Cart3d_bag *data_bag);
void Velocity_u_set_RHS(Cart3d_bag *data_bag);
void Velocity_v_set_RHS(Cart3d_bag *data_bag);
void Velocity_w_set_RHS(Cart3d_bag *data_bag);
//void Velocity_set_boundary_values(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params);
void Velocity_update_boundaries(double ***data, char component, int type, Cart3d_bag *data_bag);
void Velocity_compute_total_kinetic_energy( Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params);
void Velocity_update_world_energies(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params);
void Velocity_compute_bottom_shear_stress(Velocity *u, Velocity *v, Velocity *w,
		Cart3d_bag *data_bag);
void Velocity_nonzero_initialize(Velocity *vel, MAC_grid *grid, Parameters *params, Debug_trace *dtrace);
double Velocity_ubulk(Velocity *uvel, Cart3d_bag *data_bag);
void Velocity_calculate_dpdx(Velocity *uvel, Cart3d_bag *data_bag);
void Velocity_calc_coeff(Velocity *vel, MAC_grid *grid, Parameters *params);
void Velocity_setup_lsys_accounting_geometry(Velocity *vel, MAC_grid *grid, Parameters *params);
int Velocity_solve(Velocity *vel, Cart3d_bag *data_bag);
void Velocity_wall_shear(Cart3d_bag *data_bag);
void Velocity_u_streak(Cart3d_bag *data_bag);
void Velocity_add_buoyancy_2_RHS(Cart3d_bag *data_bag);

void Velocity_zero_initialize(Cart3d_bag *data_bag);

#ifdef LAG_PARTICLE_RESOLVED
int Velocity_solve_explicit(Velocity *vel, Cart3d_bag *data_bag);
#endif

#endif
