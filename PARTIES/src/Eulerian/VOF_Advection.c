/******************************************************************************
 * VOF_Advection.c
 * Advection of the Volume Fraction field functions for VOF computations.
 ******************************************************************************/

#include "VolumeFraction.h"  
#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Conc.h"
#include "Velocity.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Immersed.h"
#include "Display.h"
#include "Array.h"
#include "Cart3d.h"
#include "Lagrangian.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>


#ifdef VOF_PLIC

/******************************************************************************
 * VoF_set_advection
 *
 * Computes the convective term for the volume-fraction field "vof->F"
 * using a geometry-based approach (adapted from Basilisk's vof.h).
 *
 * PARTIES assumptions:
 *  - "vof->normal_x[k][j][i], normal_y, normal_z, alpha" store the reconstructed
 *    interface plane info (PLIC).
 *  - We do dimension-splitting: computing flux_x, flux_y, flux_z, then summing
 *    up in conv[k][j][i] = dFudx + dFvdy + dFwdz.
 ******************************************************************************/
void VOF_set_advection(Cart3d_bag *data_bag)
{
    /**************************************************************************
     * 1. Basic Declarations
     **************************************************************************/
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;

    // Cell-centered volume fraction and conv array
    double ***F    = vof->F;
    double ***flux_x = vof->flux_x;
    double ***flux_y = vof->flux_y;
    double ***flux_z = vof->flux_z;     
    double ***conv = vof->conv;

    // Face velocities (MAC)
    Velocity *u = data_bag->u; // x-face velocity => [k][j][i]
    Velocity *v = data_bag->v; // y-face velocity => [k][j][i]
    Velocity *w = data_bag->w; // z-face velocity => [k][j][i]

    // Grid extents
    int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    
    double SEPS = 1e-30;
    double dt   = params->dt;

    // Uniform cell spacing
    double *dx_c = grid -> dx_c;
    double *dy_c = grid -> dy_c;
    double *dz_c = grid -> dz_c;

    Memory_reset_flow_variable(grid, params, flux_x);
    Memory_reset_flow_variable(grid, params, flux_y);
    Memory_reset_flow_variable(grid, params, flux_z); 

    /**************************************************************************
     * 2. X-direction face flux
     **************************************************************************/
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) { 

                double face_vel = u->data[k][j][i];  // x-face velocity
                double un = face_vel * dt/(dx_c[i] + SEPS); 
                double s  = (un > 0.) ? 1. : ((un < 0.) ? -1. : 0.);
                // upwind offset
                int ioff = -((int)s + 1)/2; 
                int iup  = i + ioff;

                double c_up = F[k][j][iup];
                double cf   = 0.0;
                if (c_up <= 0.0) {
                    cf = 0.0;
                } else if (c_up >= 1.0) {
                    cf = 1.0;
                } else {
                    // plane normal from upwind cell
                    double nx = vof->normal_x[k][j][iup];
                    double ny = vof->normal_y[k][j][iup];
                    double nz = vof->normal_z[k][j][iup];
                    double alpha = vof->alpha[k][j][iup];

                    // Basilisk logic => plane normal => (-s*nx, ny, nz)
                    PointType plane_n = { -s*nx, ny, nz };

                    // bounding box => from (-0.5, -0.5, -0.5) to (s*un -0.5, 0.5, 0.5)
                    PointType lower = { -0.5, -0.5, -0.5 };
                    PointType upper = {  s*un - 0.5,  0.5,  0.5 };

                    cf = rectangle_fraction(plane_n, alpha, lower, upper);
                }
                flux_x[k][j][i] = cf * face_vel;
            }
        }
    }

    Communication_update_ghost_nodes_flow_variable(flux_x, FLUX_X, 
        params->ghost_nodes, data_bag);
    /**************************************************************************
     * 3. Y-direction face flux
     **************************************************************************/
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {

                double face_vel = v->data[k][j][i]; 
                double un = face_vel * dt/(dy_c[j] + SEPS);
                double s  = (un > 0.) ? 1. : ((un < 0.) ? -1. : 0.);
                int joff  = -((int)s + 1)/2;
                int jup   = j + joff;

                double c_up = F[k][jup][i];
                double cf   = 0.0;
                if (c_up <= 0.0) {
                    cf = 0.0;
                } else if (c_up >= 1.0) {
                    cf = 1.0;
                } else {
                    double nx = vof->normal_x[k][jup][i];
                    double ny = vof->normal_y[k][jup][i];
                    double nz = vof->normal_z[k][jup][i];
                    double alpha = vof->alpha[k][jup][i];

                    // plane normal => (nx, -s*ny, nz)
                    PointType plane_n = { nx, -s*ny, nz };
                    PointType lower   = { -0.5, -0.5, -0.5 };
                    PointType upper   = {  0.5,  s*un - 0.5,  0.5 };

                    cf = rectangle_fraction(plane_n, alpha, lower, upper);
                }
                flux_y[k][j][i] = cf * face_vel;
            }
        }
    }


    Communication_update_ghost_nodes_flow_variable(flux_y, FLUX_Y, 
        params->ghost_nodes, data_bag);

    /**************************************************************************
     * 4. Z-direction face flux
     **************************************************************************/
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {

                double face_vel = w->data[k][j][i]; 
                double un = face_vel * dt/(dz_c[k] + SEPS);
                double s  = (un > 0.) ? 1. : ((un < 0.) ? -1. : 0.);
                int koff  = -((int)s + 1)/2;
                int kup   = k + koff;

                double c_up = F[kup][j][i];
                double cf   = 0.0;
                if (c_up <= 0.0) {
                    cf = 0.0;
                } else if (c_up >= 1.0) {
                    cf = 1.0;
                } else {
                    double nx = vof->normal_x[kup][j][i];
                    double ny = vof->normal_y[kup][j][i];
                    double nz = vof->normal_z[kup][j][i];
                    double alpha = vof->alpha[kup][j][i];

                    // plane normal => (nx, ny, -s*nz)
                    PointType plane_n = { nx, ny, -s*nz };
                    PointType lower   = { -0.5, -0.5, -0.5 };
                    PointType upper   = {  0.5,  0.5,  s*un - 0.5 };

                    cf = rectangle_fraction(plane_n, alpha, lower, upper);
                }
                flux_z[k][j][i] = cf * face_vel;
            }
        }
    }


    Communication_update_ghost_nodes_flow_variable(flux_z, FLUX_Z, 
        params->ghost_nodes, data_bag);

    /**************************************************************************
     * 5. Final Step: conv = dFudx + dFvdy + dFwdz
     **************************************************************************/
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {

                double FxE = flux_x[k][j][i+1];
                double FxW = flux_x[k][j][i];
                double dFudx = (FxE - FxW)/dx_c[i];

                double FyN = flux_y[k][j+1][i];
                double FyS = flux_y[k][j][i];
                double dFvdy = (FyN - FyS)/dy_c[j];

                double FzF = flux_z[k+1][j][i];
                double FzB = flux_z[k][j][i];
                double dFwdz = (FzF - FzB)/dz_c[k];

                conv[k][j][i] = dFudx + dFvdy + dFwdz;
            }
        }
    }


}





/******************************************************************************
 * VOF_update_F
 *
 * Updates the volume-fraction field `F` using a 3-stage Runge–Kutta (RK3) scheme.
 * The update follows:
 *
 *   ng_rhs = GAMMA[stage] * (-conv) + ZETA[stage] * (-conv_old)
 *
 * Then, `F` is updated as:
 *
 *   F = F + dt * ng_rhs
 *
 * - `conv` is the new convective flux divergence from VOF_set_advection().
 * - `conv_old` is from the previous RK stage.
 * - Clamps `F` between [0,1] for physical consistency.
 * - If not at the final RK stage, `conv` is stored in `conv_old` for the next step.
 ******************************************************************************/
void VOF_update_F(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;

    double ***F         = vof->F;
    double ***conv      = vof->conv;
    double ***conv_old  = vof->conv_old;
    double ***rhs_vec   = vof->ng_rhs;

    int which_stage = params->which_stage;
    const double GAM[3] = {GAMMA};  // Explicit term weights
    const double ZET[3] = {ZETA};   // Previous stage correction weights
    double dt = params->dt;

    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    // Initialize RHS to zero
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                rhs_vec[k][j][i] = 0.0;
            }
        }
    }

    // Compute RHS: GAMMA*(-conv) + ZETA*(-conv_old)
    if (which_stage == 0) {
        // Stage 0: Only current convective term
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                for (int i = Is; i < Ie; i++) {
                    rhs_vec[k][j][i] += GAM[0] * (-conv[k][j][i]);
                }
            }
        }
    } else if (which_stage == 1) {
        // Stage 1: Current and previous convective terms
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                for (int i = Is; i < Ie; i++) {
                    rhs_vec[k][j][i] += GAM[1] * (-conv[k][j][i]) + ZET[1] * (-conv_old[k][j][i]);
                }
            }
        }
    } else if (which_stage == 2) {
        // Stage 2: Current and previous convective terms
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                for (int i = Is; i < Ie; i++) {
                    rhs_vec[k][j][i] += GAM[2] * (-conv[k][j][i]) + ZET[2] * (-conv_old[k][j][i]);
                }
            }
        }
    }

    // Update F: F += dt * ng_rhs
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k][j][i] += dt * rhs_vec[k][j][i];
                F[k][j][i] = clampDouble(F[k][j][i], 0.0, 1.0);
            }
        }
    }

    // Update conv_old for next stage (if not final)
    if (which_stage < 2) {
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                for (int i = Is; i < Ie; i++) {
                    conv_old[k][j][i] = conv[k][j][i];
                }
            }
        }
    }
    VOF_set_boundary_values(vof->F, data_bag);
}


/******************************************************************************
 * VOF_set_boundary_values
 *
 * Fill ghost cells of the VOF field F according to the compile-time wall
 * macros, then call the MPI ghost-exchange.  Works for any number of ghost
 * layers (NG = params->ghost_nodes).
 *
 * Conventions
 *   • Zero-gradient (F = interior value)  ⟹  impermeable wall, θ = 90°.
 *   • Dirichlet (F = 1 or user value)     ⟹  inflow / moving wall.
 *   • Outflow (copy one-way)              ⟹  simple convective exit.
 ******************************************************************************/
void VOF_set_boundary_values(double ***F, Cart3d_bag *data_bag)
{
    MAC_grid    *grid   = data_bag->grid;
    Parameters  *params = data_bag->params;
    VolumeFraction *vof = data_bag->vof;

    int NX = grid->NX;
    int NY = grid->NY;
    int NZ = grid->NZ;

    int Is = grid->G_Is; 
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    //-------------------------------------------------------------------------
    // 1) X boundaries — use G_Ks..G_Ke and G_Js..G_Je so corner ghosts are filled
    //-------------------------------------------------------------------------
#ifndef XPERIODIC

    //************************
    // Left boundary (x = 0)
    //************************
    if (Is == 0) {
        int i0 = 0;
      #ifdef LEFT_INFLOW
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i0-1] = 1.0; 
            }
        }
      
      #elif defined LEFT_OUTFLOW
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }

      #elif defined LEFT_WALL_VELOCITY_NOSLIP
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }

      #elif defined LEFT_WALL_VELOCITY_FREESLIP
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }

      #else
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }
      #endif
    }

    //************************
    // Right boundary (x = NX-1)
    //************************
    if (Ie == NX) {
        int i1 = NX - 1;
      #ifdef RIGHT_INFLOW
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i1] = 1.0;
            }
        }

      #elif defined RIGHT_OUTFLOW
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }

      #elif defined RIGHT_WALL_VELOCITY_NOSLIP
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }

      #elif defined RIGHT_WALL_VELOCITY_FREESLIP
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }

      #else
        for (int k = Ks; k < Ke; k++) {
            for (int j = Js; j < Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }
      #endif
    }
#endif // !XPERIODIC


    //-------------------------------------------------------------------------
    // 2) Y boundaries — use G_Ks..G_Ke and G_Is..G_Ie (includes x-ghost overlap)
    //-------------------------------------------------------------------------
#ifndef YPERIODIC

    if (Js == 0) {
#ifdef VOF_WETTING
        if (vof_wetting_manages_array(vof, F)) {
            /* Skip: handled by VOF_WETTING */
        } else
#endif
        {
            int j0 = 0;

          #ifdef BOTTOM_WALL_VELOCITY_NOSLIP
            for (int k = Ks; k < Ke; k++) {
                for (int i = Is; i < Ie; i++) {
                    F[k][j0-1][i] = F[k][j0][i];
                }
            }

          #elif defined BOTTOM_WALL_VELOCITY_FREESLIP
            for (int k = Ks; k < Ke; k++) {
                for (int i = Is; i < Ie; i++) {
                    F[k][j0-1][i] = F[k][j0][i];
                }
            }

          #elif defined BOTTOM_WALL_VELOCITY
            for (int k = Ks; k < Ke; k++) {
                for (int i = Is; i < Ie; i++) {
                    F[k][j0-1][i] = 1.0;
                }
            }

          #elif defined BOTTOM_WALL_SCHUMANN
            for (int k = Ks; k < Ke; k++) {
                for (int i = Is; i < Ie; i++) {
                    F[k][j0-1][i] = F[k][j0][i];
                }
            }

          #else
            for (int k = Ks; k < Ke; k++) {
                for (int i = Is; i < Ie; i++) {
                    F[k][j0-1][i] = F[k][j0][i];
                }
            }
          #endif
        }
    }

    if (Je == NY) {
        int j1 = NY - 1;
      #ifdef TOP_WALL_VELOCITY_NOSLIP
        for (int k = Ks; k < Ke; k++) {
            for (int i = Is; i < Ie; i++) {
                F[k][j1][i] = F[k][j1-1][i];
            }
        }

      #elif defined TOP_WALL_VELOCITY_FREESLIP
        for (int k = Ks; k < Ke; k++) {
            for (int i = Is; i < Ie; i++) {
                F[k][j1][i] = F[k][j1-1][i];
            }
        }

      #elif defined TOP_WALL_VELOCITY
        for (int k = Ks; k < Ke; k++) {
            for (int i = Is; i < Ie; i++) {
                F[k][j1][i] = 1.0; 
            }
        }

      #elif defined TOP_WALL_SCHUMANN
        for (int k = Ks; k < Ke; k++) {
            for (int i = Is; i < Ie; i++) {
                F[k][j1][i] = F[k][j1-1][i];
            }
        }

      #else
        for (int k = Ks; k < Ke; k++) {
            for (int i = Is; i < Ie; i++) {
                F[k][j1][i] = F[k][j1-1][i];
            }
        }
      #endif
    }
#endif // !YPERIODIC


    //-------------------------------------------------------------------------
    // 3) Z boundaries — use G_Js..G_Je and G_Is..G_Ie (includes x,y ghost overlap)
    //-------------------------------------------------------------------------
#ifndef ZPERIODIC

    if (Ks == 0) {
        int k0 = 0;
    #if defined BACK_WALL_VELOCITY_NOSLIP
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k0-1][j][i] = F[k0][j][i];
            }
        }

    #elif defined BACK_WALL_VELOCITY_FREESLIP
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k0-1][j][i] = F[k0][j][i];
            }
        }

    #elif defined BACK_WALL_VELOCITY
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k0-1][j][i] = 1.0;
            }
        }

    #else
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k0-1][j][i] = F[k0][j][i];
            }
        }
    #endif
    }

    if (Ke == NZ) {
        int k1 = NZ - 1;
    #if defined FRONT_WALL_VELOCITY_NOSLIP
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k1][j][i] = F[k1-1][j][i];
            }
        }

    #elif defined FRONT_WALL_VELOCITY_FREESLIP
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k1][j][i] = F[k1-1][j][i];
            }
        }

    #elif defined FRONT_WALL_VELOCITY
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k1][j][i] = 1.0;
            }
        }

    #else
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                F[k1][j][i] = F[k1-1][j][i];
            }
        }
    #endif
    }

#endif // !ZPERIODIC

    //-------------------------------------------------------------------------
    // 4) MPI ghost-cell update
    //-------------------------------------------------------------------------
     Communication_update_ghost_nodes_flow_variable(F,
                                                   VOLUME_FRACTION,
                                                   params->ghost_nodes,
                                                   data_bag);
}


/******************************************************************************
 * VOF_update_density_viscosity
 ******************************************************************************/
void VOF_update_density_viscosity(Cart3d_bag *data_bag)
{
    /**************************************************************************
     * 1. Basic references
     **************************************************************************/
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;

    /* Volume fraction fields:
     *  - F        : raw VOF, used for viscosity mixing
     *  - F_smooth : smoothed VOF, used for density mixing (and CSF/pressure)
     */
    double ***F_raw    = vof->F;

    /* Destination arrays for dimensionless density and viscosity */
    double ***rho_tilde = vof->rho;
    double ***mu_tilde  = vof->mu;

    /* Physical property ratios (from parameters) */
    double rho1 = params->rho1;   /* liquid density   */
    double rho2 = params->rho2;   /* gas density      */
    double mu1  = params->mu1;    /* liquid viscosity */
    double mu2  = params->mu2;    /* gas viscosity    */

#ifdef VOF_IBM
    /* Smoothed solid volume fraction: 0 outside solid, ~1 inside solid */
    double ***vfc = vof->vfc;
#endif

    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {

                /* 1) Fluid-only mixture based on VOF (liquid + gas) */

                /* Use raw F for viscosity */
                double F_ijk  = F_raw[k][j][i];

                double rho_fluid_tilde =
                    F_ijk + (1.0 - F_ijk) * (rho2 / rho1);
                double mu_fluid_tilde =
                    F_ijk + (1.0 - F_ijk) * (mu2 / mu1);

#ifdef VOF_IBM
                double Phi_s = 0.0;
                if (vfc) {
                    Phi_s = vfc[k][j][i];
                    if (Phi_s < 0.0) Phi_s = 0.0;
                    if (Phi_s > 1.0) Phi_s = 1.0;
                }

                /* NEW: Use liquid properties (1.0) for ANY cell with solid presence
                 * This ensures:
                 *   - Consistent IBM forcing (Newton's 3rd law)
                 *   - Correct buoyancy integral (Int_rho_scalar = V_p)
                 *   - Uniform proxy fluid inside solid
                 */
                const double rho_solid_proxy_tilde = 1.0;
                const double mu_solid_proxy_tilde  = 1.0;

                double rhoVal, muVal;
                if (Phi_s > 1e-6) {
                    /* Any solid presence: use proxy density */
                    rhoVal = rho_solid_proxy_tilde;
                    muVal  = mu_solid_proxy_tilde;
                } else {
                    /* Pure fluid cell */
                    rhoVal = rho_fluid_tilde;
                    muVal  = mu_fluid_tilde;
                }
#else
                double rhoVal = rho_fluid_tilde;
                double muVal  = mu_fluid_tilde;
#endif

                /* 3) Store dimensionless NS fields */
                rho_tilde[k][j][i] = rhoVal;
                mu_tilde[k][j][i]  = muVal;
            }
        }
    }

    /* Enforce BCs on the property fields */
    VOF_set_boundary_values(rho_tilde, data_bag);
    VOF_set_boundary_values(mu_tilde, data_bag);
}



#endif // VOF_PLIC