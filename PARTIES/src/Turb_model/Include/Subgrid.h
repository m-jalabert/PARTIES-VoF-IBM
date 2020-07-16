#ifndef SUBGRID_H
 #define SUBGRID_H
Subgrid *Subgrid_create(MAC_grid *grid, Parameters *params) ;
void Subgrid_concentration_derivative(Concentration *conc, MAC_grid *grid,
		Parameters *params, Subgrid *smag, int idir);
void Subgrid_smagorinsky_eddy_viscosity(Cart3d_bag *data_bag);
void Subgrid_boundary(double ***var, MAC_grid *grid, Parameters *params);
#endif
