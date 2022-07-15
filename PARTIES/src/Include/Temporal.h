#ifndef TEMPORAL_H
 #define TEMPORAL_H

int Temporal_int_rk3(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Temporal_int_all_the_equations(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Temporal_int_RANS_equations(Cart3d_bag *data_bag, double dt);
#endif
