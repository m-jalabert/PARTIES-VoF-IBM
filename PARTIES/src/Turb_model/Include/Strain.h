#ifndef STRAIN_H
 #define STRAIN_H
void Strain_rate_magnitude(Velocity *uvel, Velocity *vvel, Velocity *wvel, 
		MAC_grid *grid, Parameters *params, Strain_rate *st_rate);
#endif
