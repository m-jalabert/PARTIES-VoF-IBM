
#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H
Lagrangian *Lagrangian_create(MAC_grid *grid, Parameters *params);
void Lagrangian_destroy(Lagrangian *lag, MAC_grid *grid, Parameters *params);
void Lagrangian_advect_particles(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Lagrangian_evaluate_fluid_forces(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Lagrangian_integrate_particle_motion(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Lagrangian_collect_forces(Particle_list *p_list, Particle_list *p_list_foreign,
		int type, Parameters *params, Debug_trace *dtrace);

void Lagrangian_force(int force_iter, Cart3d_bag *data_bag, Debug_trace *dtrace);
void Lagrangian_force_individual(Particle *p, int corrector, Cart3d_bag *data_bag);

void Lagrangian_flag_points(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Lagrangian_flag_points_individual(Particle *p, Particle *p2, MAC_grid *grid);
void Lagrangian_flag_points_wall(Particle *p, MAC_grid *grid);
void Lagrangian_generate_points(Particle *p, MAC_grid *grid);
void Lagrangian_heat(int iconc, int heat_iter, Cart3d_bag *data_bag, Debug_trace *dtrace);
void Lagrangian_heat_individual(int iconc, Particle *p, int corrector, Cart3d_bag *data_bag);
void Lagrangian_mask_scalar(int iconc, Cart3d_bag *data_bag);
void Lagrangian_calc_compound_vof_fluid_velo(Cart3d_bag *data_bag);

#endif
