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
#include "VOF_DIFFUSE.h"

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
    double EPS_VOF = 1e-6;  // or even 1e-10; tune to taste
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
                if (c_up <= EPS_VOF) {
                    cf = 0.0;
                } else if (c_up >= (1.0 - EPS_VOF)) {
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
                if (c_up <= EPS_VOF) {
                    cf = 0.0;
                } else if (c_up >= (1.0 - EPS_VOF)) {
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
                if (c_up <= EPS_VOF) {
                    cf = 0.0;
                } else if (c_up >= (1.0 - EPS_VOF)) {
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

    double EPS_FLOT = 1e-6;
    for (int k = Ks; k < Ke; k++)
    for (int j = Js; j < Je; j++)
        for (int i = Is; i < Ie; i++) {
            if (F[k][j][i] < EPS_FLOT)       F[k][j][i] = 0.0;
            else if (F[k][j][i] > 1.0 - EPS_FLOT) F[k][j][i] = 1.0;
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
 * Fill ghost cells of a SCALAR field F according to the compile-time wall
 * macros, then call the MPI ghost-exchange.  Works for any number of ghost
 * layers (NG = params->ghost_nodes).
 *
 * -----------------------------------------
 * - FREE-SLIP walls now use EVEN MIRROR instead of zero-gradient:
 *
 *       F[ghost_g] = F[mirror_interior_g]
 *
 *   This preserves the field profile across the symmetry plane.
 *   Zero-gradient copies only the boundary cell value into ALL ghost
 *   layers, which flattens gradients and corrupts any stencil reading
 *   ghost cells (smoothing kernel, WLS gradient, PLIC reconstruction).
 *
 * - All other wall types remain zero-gradient (unchanged).
 *
 * - Cascading fill order X → Y → Z for edge/corner coverage (unchanged).
 *
 * Mirror index mapping (left-type boundaries, i0 = 0):
 *
 *     ghost layer g    ghost index     mirror interior index
 *     ─────────────    ───────────     ─────────────────────
 *          1             i0 − 1              i0 + 0
 *          2             i0 − 2              i0 + 1
 *          3             i0 − 3              i0 + 2
 *
 * Mirror index mapping (right-type boundaries, i1 = NX−1):
 *
 *     ghost layer g    ghost index     mirror interior index
 *     ─────────────    ───────────     ─────────────────────
 *          0             i1 + 0              i1 − 1
 *          1             i1 + 1              i1 − 2
 *          2             i1 + 2              i1 − 3
 *
 * Conventions (unchanged)
 *   • Zero-gradient (F = interior value)  ⟹  impermeable wall, θ = 90°.
 *   • Even-mirror   (F = mirrored value)  ⟹  symmetry plane (free-slip).
 *   • Dirichlet     (F = 1 or user value) ⟹  inflow / moving wall.
 *   • Outflow       (copy one-way)        ⟹  simple convective exit.
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

    /* Local extents (including ghost cells) */
    int IsL = grid->L_Is;
    int JsL = grid->L_Js;
    int KsL = grid->L_Ks;
    int IeL = grid->L_Ie;
    int JeL = grid->L_Je;
    int KeL = grid->L_Ke;

    const int NG = params->ghost_nodes;

    /*=========================================================================
     * 1) X boundaries
     *=========================================================================*/
#ifndef XPERIODIC

    /*--- Left boundary (x = 0) ---*/
    if (Is == 0) {
        int i0 = 0;

      #ifdef LEFT_INFLOW
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 1; g <= NG; g++)
            F[k][j][i0 - g] = 1.0;

      #elif defined LEFT_OUTFLOW
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 1; g <= NG; g++)
            F[k][j][i0 - g] = F[k][j][i0];

      #elif defined LEFT_WALL_VELOCITY_NOSLIP
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 1; g <= NG; g++)
            F[k][j][i0 - g] = F[k][j][i0];

      #elif defined LEFT_WALL_VELOCITY_FREESLIP
        /* ── EVEN MIRROR: symmetry plane ── */
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 1; g <= NG; g++)
            F[k][j][i0 - g] = F[k][j][i0 + g - 1];

      #else
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 1; g <= NG; g++)
            F[k][j][i0 - g] = F[k][j][i0];
      #endif
    }

    /*--- Right boundary (x = NX-1) ---*/
    if (Ie == NX) {
        int i1 = NX - 1;

      #ifdef RIGHT_INFLOW
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 0; g < NG; g++)
            F[k][j][i1 + g] = 1.0;

      #elif defined RIGHT_OUTFLOW
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 0; g < NG; g++)
            F[k][j][i1 + g] = F[k][j][i1 - 1];

      #elif defined RIGHT_WALL_VELOCITY_NOSLIP
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 0; g < NG; g++)
            F[k][j][i1 + g] = F[k][j][i1 - 1];

      #elif defined RIGHT_WALL_VELOCITY_FREESLIP
        /* ── EVEN MIRROR: symmetry plane ── */
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 0; g < NG; g++)
            F[k][j][i1 + g] = F[k][j][i1 - 1 - g];

      #else
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int g = 0; g < NG; g++)
            F[k][j][i1 + g] = F[k][j][i1 - 1];
      #endif
    }
#endif /* !XPERIODIC */


    /*=========================================================================
     * 2) Y boundaries
     *=========================================================================*/
#ifndef YPERIODIC

    /*--- Bottom boundary (y = 0) ---*/
    if (Js == 0) {
#ifdef VOF_WETTING
        if (vof_wetting_manages_array(vof, F)) {
            /* Skip: handled by VOF_WETTING */
        } else
#endif
        {
            int j0 = 0;

          #ifdef BOTTOM_WALL_VELOCITY_NOSLIP
            for (int k = KsL; k < KeL; k++)
            for (int i = IsL; i < IeL; i++)
            for (int g = 1; g <= NG; g++)
                F[k][j0 - g][i] = F[k][j0][i];

          #elif defined BOTTOM_WALL_VELOCITY_FREESLIP
            /* ── EVEN MIRROR: symmetry plane ── */
            for (int k = KsL; k < KeL; k++)
            for (int i = IsL; i < IeL; i++)
            for (int g = 1; g <= NG; g++)
                F[k][j0 - g][i] = F[k][j0 + g - 1][i];

          #elif defined BOTTOM_WALL_VELOCITY
            for (int k = KsL; k < KeL; k++)
            for (int i = IsL; i < IeL; i++)
            for (int g = 1; g <= NG; g++)
                F[k][j0 - g][i] = 1.0;

          #elif defined BOTTOM_WALL_SCHUMANN
            for (int k = KsL; k < KeL; k++)
            for (int i = IsL; i < IeL; i++)
            for (int g = 1; g <= NG; g++)
                F[k][j0 - g][i] = F[k][j0][i];

          #else
            for (int k = KsL; k < KeL; k++)
            for (int i = IsL; i < IeL; i++)
            for (int g = 1; g <= NG; g++)
                F[k][j0 - g][i] = F[k][j0][i];
          #endif
        }
    }

    /*--- Top boundary (y = NY-1) ---*/
    if (Je == NY) {
        int j1 = NY - 1;

      #ifdef TOP_WALL_VELOCITY_NOSLIP
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k][j1 + g][i] = F[k][j1 - 1][i];

      #elif defined TOP_WALL_VELOCITY_FREESLIP
        /* ── EVEN MIRROR: symmetry plane ── */
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k][j1 + g][i] = F[k][j1 - 1 - g][i];

      #elif defined TOP_WALL_VELOCITY
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k][j1 + g][i] = 1.0;

      #elif defined TOP_WALL_SCHUMANN
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k][j1 + g][i] = F[k][j1 - 1][i];

      #else
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k][j1 + g][i] = F[k][j1 - 1][i];
      #endif
    }
#endif /* !YPERIODIC */


    /*=========================================================================
     * 3) Z boundaries
     *=========================================================================*/
#ifndef ZPERIODIC

    /*--- Back boundary (z = 0) ---*/
    if (Ks == 0) {
        int k0 = 0;

      #if defined BACK_WALL_VELOCITY_NOSLIP
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 1; g <= NG; g++)
            F[k0 - g][j][i] = F[k0][j][i];

      #elif defined BACK_WALL_VELOCITY_FREESLIP
        /* ── EVEN MIRROR: symmetry plane ── */
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 1; g <= NG; g++)
            F[k0 - g][j][i] = F[k0 + g - 1][j][i];

      #elif defined BACK_WALL_VELOCITY
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 1; g <= NG; g++)
            F[k0 - g][j][i] = 1.0;

      #else
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 1; g <= NG; g++)
            F[k0 - g][j][i] = F[k0][j][i];
      #endif
    }

    /*--- Front boundary (z = NZ-1) ---*/
    if (Ke == NZ) {
        int k1 = NZ - 1;

      #if defined FRONT_WALL_VELOCITY_NOSLIP
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k1 + g][j][i] = F[k1 - 1][j][i];

      #elif defined FRONT_WALL_VELOCITY_FREESLIP
        /* ── EVEN MIRROR: symmetry plane ── */
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k1 + g][j][i] = F[k1 - 1 - g][j][i];

      #elif defined FRONT_WALL_VELOCITY
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k1 + g][j][i] = 1.0;

      #else
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int g = 0; g < NG; g++)
            F[k1 + g][j][i] = F[k1 - 1][j][i];
      #endif
    }

#endif /* !ZPERIODIC */

    /*=========================================================================
     * 4) MPI ghost-cell update
     *=========================================================================*/
    Communication_update_ghost_nodes_flow_variable(F,
                                                   VOLUME_FRACTION,
                                                   params->ghost_nodes,
                                                   data_bag);
}



void VOF_compute_conservative_momentum_fluxes(Cart3d_bag *data_bag) {

    VolumeFraction *vof = data_bag->vof;
    MAC_grid *grid      = data_bag->grid;
    Parameters *params  = data_bag->params;

    double rho1 = params->rho1;   // reference (liquid) density
    double rho2 = params->rho2;   // second-phase density  (rho2/rho1 << 1 for gas)
    double r = params->rho2 / params->rho1;
    double dt   = params->dt;

    double ***flux_x = vof->flux_x;   // γ_f * u_f * dt  (from VOF_set_advection)
    double ***flux_y = vof->flux_y;
    double ***flux_z = vof->flux_z;

    double ***u = data_bag->u->data;  // u^{k-1}
    double ***v = data_bag->v->data;
    double ***w = data_bag->w->data;

    double ***mfx = vof->mass_flux_x;
    double ***mfy = vof->mass_flux_y;
    double ***mfz = vof->mass_flux_z;

    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; k++) {
    for (int j = Js; j < Je; j++) {
    for (int i = Is; i < Ie; i++) {

        /* ---- X-face: flux_x[k][j][i] = γ_f * u_f * dt ---- */
        double u_f = u[k][j][i];
        double gamma_x = (fabs(u_f) > 1e-14)
                    ? flux_x[k][j][i] / u_f
                    : 0.5*(vof->F[k][j][i] + vof->F[k][j][i-1]);
        gamma_x = fmax(0.0, fmin(1.0, gamma_x));
        double rho_fx = gamma_x + (1.0 - gamma_x) * r;   // ρ̃_f ∈ [r, 1.0]
        mfx[k][j][i] = rho_fx * u_f;

        /* ---- Y-face ---- */
        double v_f = v[k][j][i];
        double gamma_y = (fabs(v_f) > 1e-14)
                    ? flux_y[k][j][i] / v_f
                    : 0.5*(vof->F[k][j][i] + vof->F[k][j-1][i]);
        gamma_y = fmax(0.0, fmin(1.0, gamma_y));
        double rho_fy = gamma_y + (1.0 - gamma_y) * r;
        mfy[k][j][i] = rho_fy * v_f;

        /* ---- Z-face ---- */
        double w_f = w[k][j][i];
        double gamma_z = (fabs(w_f) > 1e-14)
                    ? flux_z[k][j][i] / w_f
                    : 0.5*(vof->F[k][j][i] + vof->F[k-1][j][i]);
        gamma_z = fmax(0.0, fmin(1.0, gamma_z));
        double rho_fz = gamma_z + (1.0 - gamma_z) * r;
        mfz[k][j][i] = rho_fz * w_f;

    }}} // k,j,i

    /* Ghost node exchange so momentum loops can read halo values */
    Communication_update_ghost_nodes_flow_variable(mfx, FLUX_X, 
        params->ghost_nodes, data_bag); 
    Communication_update_ghost_nodes_flow_variable(mfy, FLUX_Y, 
        params->ghost_nodes, data_bag);
    Communication_update_ghost_nodes_flow_variable(mfz, FLUX_Z, 
        params->ghost_nodes, data_bag);
}

#endif



/******************************************************************************
 * VOF_update_density_viscosity
 ******************************************************************************/
void VOF_update_density_viscosity(Cart3d_bag *data_bag)
{
#ifdef VOF_DIFFUSE
    VOF_DIFFUSE_update_density_viscosity(data_bag);
    return;
#endif

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


                /* 2) Store dimensionless NS fields */
                rho_tilde[k][j][i] = rho_fluid_tilde;
                mu_tilde[k][j][i]  = mu_fluid_tilde;
            }
        }
    }

    /* Enforce BCs on the property fields */
    VOF_set_boundary_values(rho_tilde, data_bag);
    VOF_set_boundary_values(mu_tilde, data_bag);
}




