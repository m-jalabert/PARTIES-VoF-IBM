#ifndef COMMUNICATION_H
	#define COMMUNICATION_H

#include "DataTypes.h"

void Communication_reduce_2D_arrays(double **G_send_array, double **W_recv_array,
		Indices *G_s, Indices *G_e, Indices *W_e, int status, Cart3d_bag *data_bag);
void Communication_update_ghost_nodes_flow_variable(double ***data,
		char component, int pnodes, Cart3d_bag *data_bag);
void Communication_update_ghost_nodes_mpi(double ***data, char component, int pnodes, Cart3d_bag *data_bag);
void Communication_update_xperiodic_unodes(double ***u_data, int pnodes, Cart3d_bag *data_bag);
void Communication_update_xperiodic_cnodes(double ***c_data, int pnodes, Cart3d_bag *data_bag);
void Communication_update_zperiodic_wnodes(double ***w_data, int pnodes, Cart3d_bag *data_bag);
void Communication_update_zperiodic_cnodes(double ***c_data, int pnodes, Cart3d_bag *data_bag);
void Communication_update_ghost_nodes_x(double ***data, char which_quantity,
		int pnodes, Cart3d_bag *data_bag);
void Communication_update_ghost_nodes_y(double ***data, char which_quantity,
		int pnodes, Cart3d_bag *data_bag);
void Communication_update_ghost_nodes_z(double ***data, char which_quantity,
		int pnodes, Cart3d_bag *data_bag);
void Communication_new_xyz_communicator(MAC_grid *grid, Parameters *params);
int proc_find(int *istartg,int *jstartg,int *kstartg, int size, int ilocal, int jlocal, int klocal);
void bubble_sort(int a[], int size);
void xzperiodic_uvel_ave(Cart3d_bag *data_bag,Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params);
void xzperiodic_saltsediment_ave(Cart3d_bag *data_bag,Concentration **c, MAC_grid *grid, Parameters *params);
void xzperiodic_thermal_ave(Cart3d_bag *data_bag, Concentration **c, MAC_grid *grid, Parameters *params);
void xzperiodic_nut_ave(Cart3d_bag *data_bag, Subgrid *smag, MAC_grid *grid, Parameters *params);
void xzperiodic_rans_nut_ave(Cart3d_bag *data_bag, Rans *rans, MAC_grid *grid, Parameters *params);
void Communication_finalize();
void rans_nut_ave(Rans *rans, MAC_grid *grid, Parameters *params, int iter, int rk);

#endif // notCOMMUNICATION_H
