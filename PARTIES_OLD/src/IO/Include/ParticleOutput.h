#ifndef PARTICLE_OUTPUT_H
 #define PARTICLE_OUTPUT_H
#include "hdf5.h"
#include "DataTypes.h"

void ParticleOutput_dat(Particle_list *p_list, MAC_grid *grid, Parameters *params,
		Debug_trace *dtrace);
void ParticleOutput_h5(Cart3d_bag *data_bag, int counter, Debug_trace *dtrace);
void ParticleOutput_h5_data(Particle_list *p_list, hid_t file_id,
		char *groupname, Parameters *params, Debug_trace *dtrace);
void ParticleOutput_h5_data_element(Particle_list *p_list, int Np_local,
		hid_t file_id, char *groupname, char *element, int dataSize,
		Parameters *params, Debug_trace *dtrace);
void ParticleOutput_h5_collision(Particle_list *p_list, int Np_local, int Nc_max,
		hid_t file_id, char *groupname, Parameters *params, Debug_trace *dtrace);

#endif
