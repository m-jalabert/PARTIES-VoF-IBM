#ifndef MARKER_H
#define MARKER_H
#include "DataTypes.h"

void Marker_startFile(DataParticle *p_data_start, int mType, Parameters *params);
void Marker_endFile(DataParticle *p_data_start, int mType, Parameters *params);
void Marker_writeTimestep(DataParticle *p_data_start, int mType, int iter,
		double time, Parameters *params);
void Marker_writeAttribute(FILE *xmf, DataParticle *p_data, int mType, int iter, int N_pts);

#define MARKER_TOPOLOGY_FILE "Marker_topology.h5"
#define MARKER_GROUP_0 "/marker_black"
#define MARKER_GROUP_1 "/marker_red"

#endif
