#ifndef PRESSURE_H
 #define PRESSURE_H

#include "definitions.h"

Pressure *Pressure_create(MAC_grid *grid, Parameters *params);
void Pressure_destroy(Pressure *p, MAC_grid *grid, Parameters *params);
void Pressure_set_RHS(Cart3d_bag *data_bag);
void Pressure_project_velocity(Cart3d_bag *data_bag);
double Pressure_compute_velocity_divergence(Cart3d_bag *data_bag);
void Pressure_setup_lsys_accounting_geometry(Pressure *p, MAC_grid *grid,
		Parameters *params);
int Pressure_solve(Pressure *p, MAC_grid *grid, Parameters *params);
#endif
