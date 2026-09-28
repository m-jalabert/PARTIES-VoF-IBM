#ifndef VOF_DIFFUSE_H
#define VOF_DIFFUSE_H

#include "DataTypes.h"

/* Version 2 adds canonical periodic solid geometry to ice transport. */
#define VOF_DIFFUSE_SEDIMENT_TRANSPORT_VERSION 2

/*
 * Threshold above which a cell counts as RESOLVED SEDIMENT rather than fluid.
 * Shared with the ICE_PENALIZATION mask in lsolver/, which must not apply
 * Brinkman damping inside the resolved solid -- the IBM governs there.
 * Was a private #define in VOF_DIFFUSE.c; promoted here so the two cannot drift.
 */
#ifndef DIFFUSE_SOLID_MASS_CUTOFF
#define DIFFUSE_SOLID_MASS_CUTOFF 0.05
#endif

/*
 * Halo around a resolved grain in which the ternary ice indicator
 * phi_s = 1 - C_L - C_S is NOT trustworthy.
 *
 * Brinkman damping suppresses the flow near the grain, so Cahn-Hilliard cannot
 * relax C_L up to 1 - C_S in the wake; the shortfall is then read as ice, which
 * damps harder -- a self-reinforcing feedback.  Measured: phi_s = 0.15-0.30 over
 * R..R+3 cells in an ICE-FREE deck, and LARGER at weaker damping (0.16 at
 * tau=1e-3, 0.30 at tau=3e-3), which is the signature of the feedback rather
 * than of real ice.
 *
 * C_S = 1e-3 sits at R + 4.2 cells for the 3-D shakeout, covering that shell.
 * (1e-12, the geometric support cutoff, would reach R + 26 cells and disable
 * ice rigidity across a large part of the domain -- far too wide.)
 */
#ifndef DIFFUSE_SOLID_HALO_CUTOFF
#define DIFFUSE_SOLID_HALO_CUTOFF 1.0e-3
#endif
/* Above this, a haloed cell is confidently bulk ice and stays rigid. */
#ifndef DIFFUSE_SOLID_HALO_ICE_TRUST
#define DIFFUSE_SOLID_HALO_ICE_TRUST 0.9
#endif

/*
 * Hard ice threshold for ICE_PENALIZATION, matching the reference method:
 * Yang/AFiD use "eta = dt penalty PLUS hard u = 0 for phi_s > 0.9" (roadmap G.8).
 *
 * Why it is needed here: a MOVING resolved grain leaves a wake of cells it has
 * vacated that Cahn-Hilliard has not yet refilled to C_L = 1 - C_S.  The
 * shortfall reads as ice.  Measured in an ICE-FREE settling deck: 97% of cells
 * carry phi_s > 0 at a mean of 0.022 -- ~19% spurious drag EVERYWHERE at
 * tau = 1e-3 -- plus 1,772 wake cells above 0.9.  A graded phi_s therefore
 * brakes the particle no matter what tau is chosen; a hard threshold removes
 * the distributed part (97% -> 0.27% of cells) while leaving bulk ice fully
 * rigid.
 *
 * Enabled by VOF_DIFFUSE_ICE_PENAL_THRESHOLD in Boundary.h; UNDEFINED by
 * default, so Stage A / Yang / melting-RB keep the graded mask they were
 * validated with.
 */
#ifndef DIFFUSE_ICE_PENAL_THRESHOLD_VAL
#define DIFFUSE_ICE_PENAL_THRESHOLD_VAL 0.9
#endif

#ifdef VOF_DIFFUSE
void VOF_DIFFUSE_init(Cart3d_bag *db);
void VOF_DIFFUSE_set_boundary_values(double ***f, Cart3d_bag *db);
void VOF_DIFFUSE_update_phase_cache(Cart3d_bag *db);
void VOF_DIFFUSE_remove_disconnected_gas(Cart3d_bag *db);

void VOF_DIFFUSE_compute_C_S(Cart3d_bag *db);
#ifdef VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT
void VOF_DIFFUSE_move_ice_mask(Cart3d_bag *db);
#endif
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

#if defined(PHASE_CHANGE) && defined(CONC)
void VOF_DIFFUSE_compute_melt_rate(Cart3d_bag *db);
#endif


void VOF_DIFFUSE_compute_f_sigma(Cart3d_bag *db);
void VOF_DIFFUSE_apply_f_sigma_old(Cart3d_bag *db);

void VOF_DIFFUSE_update_density_viscosity(Cart3d_bag *db);


#endif

#endif
