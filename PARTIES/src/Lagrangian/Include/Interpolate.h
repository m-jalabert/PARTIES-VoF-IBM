#ifndef INTERPOLATE_H
#define INTERPOLATE_H

void Interpolate_Lag_to_Eul(double *A, double ***a, char which, Particle *p,
		MAC_grid *grid);
void Interpolate_Eul_to_Lag(double ***a, double *A, char which, Particle *p,
		MAC_grid *grid);
void Interpolate_integrate_momentum(Velocity *vel, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void Interpolate_paint_rigid_body_velocity(Velocity *vel, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void Interpolate_add_to_volume_fraction(char component, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void Interpolate_add_to_volume_fraction_vof(char component, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void Interpolate_add_to_volume_fraction_vof_prime(char component, Particle_list *p_list,
		Cart3d_bag *data_bag, Debug_trace *dtrace);
void mat_vec(int m, int n, double **A, double *B, double *C);

void Interpolate_bound_to_one( double ***vf, Cart3d_bag *data_bag);

// Add to the existing declarations

#ifdef VOF_IBM
void Interpolate_CCF_to_particle(Particle *p, Cart3d_bag *data_bag);
void Integrate_CCF_to_particle_Eulerian(Particle *p, Cart3d_bag *data_bag);
#endif
#endif
