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


/******************************************************************************/
// Pressure solver functions for VoF-PLIC
/******************************************************************************/

/******************************************************************************
 * Pressure_compute_preconditioner
 *
 * Builds a **harmonic-average Jacobi** preconditioner for
 *      A φ = ∇·((1/ρ) ∇φ)
 * on a uniform Cartesian grid.  The diagonal entry is
 *
 *   a_ii = (c_E + c_W)/Δx² + (c_N + c_S)/Δy² + (c_T + c_B)/Δz²
 *   with  c_F = 2 / (ρ_i + ρ_F)   (harmonic average of 1/ρ).
 *
 * M_inv = 1/a_ii is stored in p->M_inv.
 ******************************************************************************/
void Pressure_compute_preconditioner(Cart3d_bag *data_bag);




/******************************************************************************/
// Conjugate Gradient solver for the variable-coefficient Poisson system:
//
//     ∇ · ( (1/ρ) ∇φ ) = rhs
//
/******************************************************************************/
int Pressure_solve_cg(Cart3d_bag *data_bag);



/******************************************************************************/
// Pressure_operator_variableCoeff
// Computes: Aphi = ∇·( 1/rho * ∇phi )
/******************************************************************************/
void Pressure_operator_variableCoeff(
    double ***Aphi,         // output: operator(A) * phi
    double ***phi,          // input: phi array
    double ***rho,       // rho at cell centers
    MAC_grid *grid,
    Parameters *params,
    Cart3d_bag *data_bag );



/******************************************************************************/
// Applies physical boundary conditions for pressure.
//
/******************************************************************************/
void Pressure_apply_BCs(double ***phi, MAC_grid *grid, Parameters *params);

/******************************************************************************/
// Initialize the pressure field hydrostatically.
/******************************************************************************/
void Pressure_init_hydrostatic_VOF(Cart3d_bag *data_bag);

/******************************************************************************/
// HYPRE-based pressure solver (PCG + PFMG)
/******************************************************************************/
#ifdef USE_HYPRE
void Pressure_hypre_setup(Cart3d_bag *data_bag);
int  Pressure_solve_hypre(Cart3d_bag *data_bag);
void Pressure_hypre_destroy(void);
void Pressure_hypre_mark_dirty(void);
#endif


#endif
