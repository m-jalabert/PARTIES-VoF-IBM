#ifndef LESFILTER_H
 #define LESFILTER_H
void Lesfilter_init(Subgrid *smag, MAC_grid *grid, Parameters *params);
void lesfilter(double ***u, double ***uf, double ***uf1, Cart3d_bag *data_bag);

#endif
