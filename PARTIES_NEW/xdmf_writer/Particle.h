#ifndef PARTICLE_H
#define PARTICLE_H
#include "DataTypes.h"

void Particle_startFile(DataParticle *p_data_start, char pType, Parameters *params);
void Particle_endFile(DataParticle *p_data_start, char pType, Parameters *params);
void Particle_writeTimestep(DataParticle *p_data_start, char pType, int index,
		double time, int *N_particles, Parameters *params);
void Particle_writeAttribute(FILE *xmf, DataParticle *p_data, char pType, int iter, int N_particles);

#define MOBILE_GROUP "/mobile"
#define FIXED_GROUP "/fixed"

#endif
