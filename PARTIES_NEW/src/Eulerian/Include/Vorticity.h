#ifndef VORTICITY_H
 #define VORTICITY_H

#include "definitions.h"
#include "DataTypes.h"

void Vorticity_create(Cart3d_bag *data_bag);
void Vorticity_destroy(Cart3d_bag *data_bag);
void Vorticity_compute(Cart3d_bag *data_bag);

#endif
