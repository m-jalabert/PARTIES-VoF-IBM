#ifndef POSTPROCESS_H
#define POSTPROCESS_H
#include "DataTypes.h"

void PostProcess(int argc, char **args, Cart3d_bag *data_bag, Debug_trace *dtrace);
void parse_input(int argc, char **args, int *i_start, int *i_step, int *i_end,
		Parameters *params, Debug_trace *dtrace);

#endif
