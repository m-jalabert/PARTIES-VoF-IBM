#ifndef OUTFLOW_H
 #define OUTFLOW_H

void Outflow_impose_convective_boundary(Cart3d_bag *data_bag);
void Outflow_vel_impose_convective_boundary(Velocity *u, Velocity *v,
		Velocity *w, MAC_grid *grid, Parameters *params);
void Outflow_conc_impose_convective_boundary(Concentration *c, MAC_grid *grid,
		Parameters *params);
void Outflow_update_u_velocity_to_conserve_mass(Velocity *u, MAC_grid *grid,
		Parameters *params);

#endif
