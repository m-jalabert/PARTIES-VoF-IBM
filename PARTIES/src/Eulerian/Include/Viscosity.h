#ifndef VISCOSITY_H
 #define VISCOSITY_H

Viscosity *Viscosity_create(MAC_grid *grid, Parameters *params);
void Viscosity_destroy(Viscosity *visc, MAC_grid *grid, Parameters *params);
void Viscosity_set_cell_edges(Cart3d_bag *data_bag);
void Viscosity_update_boundaries(double ***nut, Cart3d_bag *data_bag);

#endif
