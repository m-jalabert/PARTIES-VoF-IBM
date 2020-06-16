#ifndef PARTDRY_H
#define PARTDRY_H

#include "DataTypes.h"

int PartDry_int_rk3(Cart3d_bag *data_bag, Debug_trace *dtrace);
void PartDry_int_all_the_equations(Cart3d_bag *data_bag, Debug_trace *dtrace);
double PartDry_get_dt(Cart3d_bag *data_bag);

#endif
