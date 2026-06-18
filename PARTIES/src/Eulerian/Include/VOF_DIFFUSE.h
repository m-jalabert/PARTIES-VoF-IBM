#ifndef VOF_DIFFUSE_H
#define VOF_DIFFUSE_H

#include "DataTypes.h"

#ifdef VOF_DIFFUSE
void VOF_DIFFUSE_init(Cart3d_bag *db);
void VOF_DIFFUSE_set_boundary_values(double ***f, Cart3d_bag *db);
void VOF_DIFFUSE_update_phase_cache(Cart3d_bag *db);
void VOF_DIFFUSE_remove_disconnected_gas(Cart3d_bag *db);

void VOF_DIFFUSE_compute_C_S(Cart3d_bag *db);
void VOF_DIFFUSE_apply_contact_angle(Cart3d_bag *db);
void VOF_DIFFUSE_extend_psi_LG_contact_angle(Cart3d_bag *db);
void VOF_DIFFUSE_compute_solid_normals_MCL(Cart3d_bag *db);
void VOF_DIFFUSE_compute_laplacian(double ***in, double ***lap, Cart3d_bag *db);
void VOF_DIFFUSE_compute_bulk_S(Cart3d_bag *db);
void VOF_DIFFUSE_compute_psi(Cart3d_bag *db);
void VOF_DIFFUSE_compute_psi_LG(Cart3d_bag *db);

void VOF_DIFFUSE_advect_WENO5(Cart3d_bag *db, double ***rhs_out);
void VOF_DIFFUSE_explicit_diffusion(Cart3d_bag *db, double ***rhs_out);
void VOF_DIFFUSE_solve_implicit_biharmonic(Cart3d_bag *db, double ***rhs_explicit);
void VOF_DIFFUSE_step(Cart3d_bag *db);


void VOF_DIFFUSE_compute_f_sigma(Cart3d_bag *db);
void VOF_DIFFUSE_apply_f_sigma_old(Cart3d_bag *db);

void VOF_DIFFUSE_update_density_viscosity(Cart3d_bag *db);


#endif

#endif
