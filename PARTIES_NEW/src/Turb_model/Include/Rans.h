#ifndef RANS_H
 #define RANS_H

Rans *Rans_create(MAC_grid *grid, Parameters *params);
void Rans_destroy(Rans *rans, Parameters *params, MAC_grid *grid);
void Rans_int_equations(Cart3d_bag *data_bag, double dt, Debug_trace *dtrace);
void Rans_eddy_viscosity(Rans *rans, Velocity *u, Parameters *params, MAC_grid *grid);
void Rans_set_boundary_values(double ***nut, MAC_grid *grid, Parameters *params);

#endif
