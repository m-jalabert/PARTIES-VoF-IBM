#ifndef TWO_EQN_RANS_H
 #define TWO_EQN_RANS_H

#include "definitions.h"

Two_equation_rans *Two_equation_rans_create(MAC_grid *grid, Parameters *params);
void Two_equation_rans_destroy(Two_equation_rans *two_eqn_rans, 
		Parameters *params, MAC_grid *grid);
void Rans_initialize (Rans *rans, Velocity *u, Cart3d_bag *data_bag);
void Rans_set_conv_viscous( Velocity *u, Velocity *v, Velocity *w,
		Rans *rans, MAC_grid *grid, Parameters *params, int eqn_no);
void Rans_set_RHS(Two_equation_rans *two_eqn_rans, MAC_grid *grid, 
		Parameters *params, int eqn_no, double dt);
void Rans_add_source_RHS(Cart3d_bag *data_bag);
void Rans_set_quick_coefficients(Two_equation_rans *two_eqn_rans, 
		MAC_grid *grid, Parameters *params);
void Two_equation_set_boundary_values(Rans *rans, MAC_grid *grid, 
		Parameters *params);
int Two_equation_solve(Cart3d_bag *data_bag, double dt);
#endif
