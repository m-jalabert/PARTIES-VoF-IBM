#ifndef INTERPOLATE_H
#define INTERPOLATE_H

void Interpolate_Lag_to_Eul(double *A, double ***a, char which, Particle *p,
		MAC_grid *grid);
void Interpolate_Eul_to_Lag(double ***a, double *A, char which, Particle *p,
		MAC_grid *grid);
void Interpolate_integrate_momentum(Velocity *vel, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void Interpolate_add_to_volume_fraction(char component, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void Interpolate_add_to_volume_fraction_vof(char component, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void Interpolate_add_to_volume_fraction_vof_prime(char component, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);

void Interpolate_bound_to_one( double ***vf, Cart3d_bag *data_bag);
#endif
