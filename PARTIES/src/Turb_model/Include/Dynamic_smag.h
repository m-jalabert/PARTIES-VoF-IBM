#ifndef DYNAMIC_H
 #define DYNAMIC_H

void Dynamic_smag_init(Subgrid *smag, MAC_grid *grid, Parameters *params);
void Dynamic_smag_coeff(Cart3d_bag *data_bag);
void Dynamic_smag_interpolate_xminusdt(Velocity *u, Velocity *v, Velocity *w, Subgrid *smag, MAC_grid *grid, Parameters *params);
void Dynamic_smag_global_horizontal_mean(double ***quantity, double *global_mean, MAC_grid *grid) ;
void Dynamic_smag_span_ave(double ***quantity, double **global_mean, double **work, MAC_grid *grid, Parameters *params);

#endif
