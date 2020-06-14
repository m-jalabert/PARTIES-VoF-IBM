#ifndef H_PARTICLE_INPUT
 #define H_PARTICLE_INPUT
#include "hdf5.h"
#include "DataTypes.h"

void ParticleInput_inp(Particle_list *p_list, MAC_grid *grid, Parameters *params,
		Debug_trace *dtrace);
void ParticleInput_h5(Cart3d_bag *data_bag, Debug_trace *dtrace);
void ParticleInput_h5_data(Particle_list *p_list, hid_t file_id, char *groupname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace);
int ParticleInput_h5_data_element(double ***data, int dataSize, hid_t file_id,
		char *groupname, char *element, Parameters *params, Debug_trace *dtrace);
void ParticleInput_h5_collision(Particle_list *p_list, hid_t file_id,
		char *groupname, Parameters *params, Debug_trace *dtrace);
int ParticleInput_copy_collision_data(Collision **pc_start_ptr, int j,
		int Nc, double *data);

#endif
