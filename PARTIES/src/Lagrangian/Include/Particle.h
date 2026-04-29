
#ifndef PARTICLE_H
#define PARTICLE_H

#include "hdf5.h"
#include "DataTypes.h"

void Particle_initialize(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Particle_initialize_MPI_datatype();
void Particle_initialize_velocities(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Particle_initialize_nonessential_data(Particle *p);
void Particle_initialize_volume_fraction(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Particle_calc_derived_data(Particle *p, MAC_grid *grid, Parameters *params);
void Particle_destroy(Particle *p);
void Particle_list_destroy(Particle_list *p_list);
void Particle_mobile_turn_off(Particle_list *p_list, Parameters *params);
Particle_list *Particle_list_foreign_create(Particle_list *p_list, Cart3d_bag *data_bag,
		Debug_trace *dtrace);
void Particle_list_copy(Particle_list *p_list_src, Particle_list *p_list_new);
void Particle_release_to_mobile(Particle_list *p_list_release, Particle_list *p_list_mobile, Parameters *params, Debug_trace *dtrace);
void Particle_list_add_array(Particle_list *p_list, Particle *p_array, int p_size,
		Collision *pc_array, int *Nc, MAC_grid *grid);
void Particle_list_remove(Particle_list *p_list, int rm_type, MAC_grid *grid,
		Parameters *params, Debug_trace *dtrace);

void Particle_collision_list_destroy(Particle *p);
void Particle_collision_list_copy(Particle *p_src, Particle *p_new);
int Particle_collision_list_xfer(Particle *p_src, Particle *p_new, char type);
void Particle_collision_list_add_array(Collision **pc_start_ptr,
		Collision *pc_array, int pc_size);
void Particle_collision_list_clean(Particle *p_start);

void Particle_create_internal_arrays(Particle *p);
void Particle_destroy_internal_arrays(Particle *p);

#if defined(LAG_PARTICLE_RESOLVED)
int Particle_center_is_owned_by_rank(const Particle *p, MAC_grid *grid,
		Parameters *params, int rank);
int Particle_center_owner_rank(const Particle *p, MAC_grid *grid,
		Parameters *params);
Particle *Particle_collect_owned_overlaps(Particle_list *p_list,
		Cart3d_bag *data_bag, double extra_range, double min_radius,
		int include_self, int *n_recv);
void Particle_reduce_oversized_forces_to_owner(Particle_list *p_list,
		Cart3d_bag *data_bag);
#endif

void Particle_MPI_update(Particle_list *p_list, Cart3d_bag *data_bag,
		Debug_trace *dtrace);

#endif
