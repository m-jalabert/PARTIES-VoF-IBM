
#ifndef CONC_H
	#define CONC_H

#include "definitions.h"

Concentration *Conc_create(int conc_index, MAC_grid *grid, Parameters *params);
void Conc_destroy(Concentration *c, int conc_index, MAC_grid *grid, Parameters *params);
void Conc_int_equations(Cart3d_bag *data_bag, Debug_trace *dtrace);
//int Conc_get_lsys_index (Concentration *c, MAC_grid *grid, int x_index, int y_index, int z_index, char which) ;
double Conc_settling_speed_function(double conc, double phi_max, double V_s0);
void Conc_compute_total_concentration(Cart3d_bag *data_bag);

void Conc_set_conv_viscous(int iconc, Cart3d_bag *data_bag);
void Conc_set_conv_outofbounds(int iconc, Cart3d_bag *data_bag);
void Conc_set_RHS(Concentration *c, MAC_grid *grid, Parameters *params);
void Conc_set_RHS_IBM_implicit(Concentration *c, double ***temp_f ,MAC_grid *grid, Parameters *params);
void Conc_add_source_RHS(int iconc, Cart3d_bag *data_bag);
void Conc_store_old_data(Concentration *c, MAC_grid *grid, Parameters *params) ;
void Conc_set_boundary_values(double ***data , int iconc ,int type, MAC_grid *grid , Parameters *params);
void Conc_set_boundary_values_fortemp(Concentration *c, MAC_grid *grid, Parameters *params);
int Conc_calc_just_outofbounds(int iconc, double tol, Cart3d_bag *data_bag);
int Conc_calc_outofbounds(int iconc, int old_outofbounds, int* tot_interval,double tol, Cart3d_bag *data_bag);
void Conc_show_outofbounds(int iconc, int old_outofbounds, double tol,
						   Cart3d_bag *data_bag);
void Conc_set_blendingfactor(Concentration *c, MAC_grid *grid, Parameters *params, double blend_fac) ;
void Conc_average_outofbounds(Concentration *c, MAC_grid *grid, Parameters *params);

void Conc_sdomain_copy(Concentration *c, MAC_grid *grid, Parameters *params);
void Conc_sdomain_copy_totemp(Concentration *c, MAC_grid *grid, Parameters *params);
void Conc_copy_fromtemp(Concentration *c, MAC_grid *grid, Parameters *params);
void Conc_set_quick_coefficients(Concentration *c, MAC_grid *grid, Parameters *params);

void cltridiag(double *a, double *b, double *c, double *r, double *u, int N);
void Conc_calc_coeff(Concentration *c, MAC_grid *grid, Parameters *params);
void Conc_setup_lsys_accounting_geometry(Concentration *c, MAC_grid *grid, Parameters *params);
int Conc_solve_semi_implicit(int iconc, Cart3d_bag *data_bag);
int Conc_solve_fully_explicit(int iconc, Cart3d_bag *data_bag);
int Conc_solve_cg(int iconc, Cart3d_bag *data_bag);
void Conc_set_conv_viscous_central_mixed_vof(int iconc, Cart3d_bag *data_bag);
void Conc_set_conv_viscous_central_mixed(int iconc, Cart3d_bag *data_bag);

#endif // notCONC_H
