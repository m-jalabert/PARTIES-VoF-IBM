#ifndef INFLOW_H
 #define INFLOW_H

void Inflow_velocity_profile(Cart3d_bag *data_bag, Debug_trace *dtrace);

int Get_inflow_data(double **inflow, double **inflow_n, double **inflow_o,
		char component, int read_new_data, Cart3d_bag *data_bag);

#endif
