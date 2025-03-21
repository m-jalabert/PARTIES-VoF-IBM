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
    Communication_update_ghost_nodes_flow_variable(vof->F, VOLUME_FRACTION,
        params->ghost_nodes, data_bag);
}


/******************************************************************************
 * VoF_set_boundary_values
 *
 * Set boundary values for volume fraction, then do MPI ghost updates.
 * This expanded version showcases how you might handle more boundary 
 * condition types (inflow, outflow, contact angle, free-slip, etc.).
 ******************************************************************************/
void VOF_set_boundary_values(double ***F, Cart3d_bag *data_bag)
{
    MAC_grid    *grid   = data_bag->grid;
    Parameters  *params = data_bag->params;

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
    // 1) X boundaries
    //-------------------------------------------------------------------------
#ifndef XPERIODIC

    //************************
    // Left boundary (x = 0)
    //************************
    if (Is == 0) {
        int i0 = 0;
      #ifdef LEFT_INFLOW
        /*
         * Example: if you want a fluid inflow with volume fraction = 1.0 
         * at x=0 (i.e., pure fluid enters),
         * then set the ghost cell to 1.0. 
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i0-1] = 1.0; 
            }
        }
      
      #elif defined LEFT_OUTFLOW
        /*
         * For an "outflow" or "convective" boundary, we might do
         * a zero-gradient approach: F[k][j][i0-1] = F[k][j][i0].
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }

      #elif defined LEFT_WALL_VELOCITY_NOSLIP
        /*
         * If we consider a "wall" for VOF, often we do zero flux 
         * or zero gradient. Possibly also impose contact angle 
         * conditions if relevant. 
         * For now, let's just do zero gradient:
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }

      #elif defined LEFT_WALL_VELOCITY_FREESLIP
        /*
         * Similarly, for a "free-slip" boundary, typically no normal flux 
         * is still correct for VOF. 
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }

      #else
        /*
         * Default fallback: mirror or zero-gradient boundary
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i0-1] = F[k][j][i0];
            }
        }
      #endif // left boundary condition
    } // Is == 0

    //************************
    // Right boundary (x = NX-1)
    //************************
    if (Ie == NX) {
        int i1 = NX - 1;
      #ifdef RIGHT_INFLOW
        /*
         * If fluid inflows from the right side, set volume fraction to 1.
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i1] = 1.0;
            }
        }

      #elif defined RIGHT_OUTFLOW
        /*
         * Zero-gradient outflow, or 'convective' approach
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }

      #elif defined RIGHT_WALL_VELOCITY_NOSLIP
        /*
         * No-slip wall => typically zero flux for VOF
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }

      #elif defined RIGHT_WALL_VELOCITY_FREESLIP
        /*
         * Free-slip => again no normal flux
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }

      #else
        /*
         * Default fallback: mirror boundary
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int j = grid->L_Js; j < grid->L_Je; j++) {
                F[k][j][i1] = F[k][j][i1-1];
            }
        }
      #endif // right boundary condition
    } // Ie == NX
#endif // !XPERIODIC


    //-------------------------------------------------------------------------
    // 2) Y boundaries
    //-------------------------------------------------------------------------
#ifndef YPERIODIC

    //************************
    // Bottom boundary (y=0)
    //************************
    if (Js == 0) {
        int j0 = 0;
      #ifdef BOTTOM_WALL_VELOCITY_NOSLIP
        /*
         * For a no-slip bottom wall, we can do zero flux for VOF
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j0-1][i] = F[k][j0][i];
            }
        }

      #elif defined BOTTOM_WALL_VELOCITY_FREESLIP
        /*
         * Or free-slip also => zero flux approach
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j0-1][i] = F[k][j0][i];
            }
        }

      #elif defined BOTTOM_WALL_VELOCITY
        /*
         * If you had some "moving wall" or "inflow" at bottom, 
         * you might set F=1.0 or a user-specified function.
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j0-1][i] = 1.0; // example
            }
        }

      #elif defined BOTTOM_WALL_SCHUMANN
        /*
         * Some specialized boundary condition. 
         * Typically you'd do a partial flux approach or a special function.
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                // Placeholder for specialized BC 
                F[k][j0-1][i] = F[k][j0][i];
            }
        }

      #else
        /*
         * Default fallback: zero-gradient
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j0-1][i] = F[k][j0][i];
            }
        }
      #endif // bottom boundary condition
    }

    //************************
    // Top boundary (y = NY-1)
    //************************
    if (Je == NY) {
        int j1 = NY - 1;
      #ifdef TOP_WALL_VELOCITY_NOSLIP
        /*
         * Another no-slip boundary => zero flux for VOF
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j1][i] = F[k][j1-1][i];
            }
        }

      #elif defined TOP_WALL_VELOCITY_FREESLIP
        /*
         * Free-slip => zero flux approach
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j1][i] = F[k][j1-1][i];
            }
        }

      #elif defined TOP_WALL_VELOCITY
        /*
         * Possibly set F=1 if fluid is injected from top
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j1][i] = 1.0; 
            }
        }

      #elif defined TOP_WALL_SCHUMANN
        /*
         * Specialized BC again
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                // some specialized approach
                F[k][j1][i] = F[k][j1-1][i];
            }
        }

      #else
        /*
         * Default fallback: mirror
         */
        for (int k = grid->L_Ks; k < grid->L_Ke; k++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k][j1][i] = F[k][j1-1][i];
            }
        }
      #endif // top boundary condition
    }
#endif // !YPERIODIC


//-------------------------------------------------------------------------
// 3) Z boundaries
//-------------------------------------------------------------------------
#ifndef ZPERIODIC

    // Back boundary (z = 0)
    if (Ks == 0) {
        int k0 = 0;
    #if defined BACK_WALL_VELOCITY_NOSLIP
        /*
         * Typical no-slip => zero flux for VOF
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k0-1][j][i] = F[k0][j][i];
            }
        }

    #elif defined BACK_WALL_VELOCITY_FREESLIP
        /*
         * Free-slip => typically zero gradient for volume fraction
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k0-1][j][i] = F[k0][j][i];
            }
        }

    #elif defined BACK_WALL_VELOCITY
        /*
         * Example scenario: a moving or inflow boundary at z=0,
         * e.g. set VOF=1.0 if fluid is injected from the back
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k0-1][j][i] = 1.0;
            }
        }

    #else
        /*
         * Default fallback => mirror or zero gradient
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k0-1][j][i] = F[k0][j][i];
            }
        }
    #endif // back wall macros
    }

    // Front boundary (z = NZ - 1)
    if (Ke == NZ) {
        int k1 = NZ - 1;
    #if defined FRONT_WALL_VELOCITY_NOSLIP
        /*
         * Another no-slip => zero flux for VOF
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k1][j][i] = F[k1-1][j][i];
            }
        }

    #elif defined FRONT_WALL_VELOCITY_FREESLIP
        /*
         * Free-slip => zero gradient approach
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k1][j][i] = F[k1-1][j][i];
            }
        }

    #elif defined FRONT_WALL_VELOCITY
        /*
         * Possibly set F=1.0 if fluid is injected from front
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k1][j][i] = 1.0;
            }
        }

    #else
        /*
         * Default fallback => zero gradient
         */
        for (int j = grid->L_Js; j < grid->L_Je; j++) {
            for (int i = grid->L_Is; i < grid->L_Ie; i++) {
                F[k1][j][i] = F[k1-1][j][i];
            }
        }
    #endif // front wall macros
    }

#endif // !ZPERIODIC



    //-------------------------------------------------------------------------
    // 4) MPI ghost‐cell update
    //-------------------------------------------------------------------------
     Communication_update_ghost_nodes_flow_variable(F,
                                                   VOLUME_FRACTION,
                                                   params->ghost_nodes,
                                                   data_bag);
}

/******************************************************************************
 * VoF_update_density_viscosity
 *
 * Updates the dimensionless density and viscosity fields in the domain,
 * using the piecewise mixing laws:
 *
 *   (1)  rho_tilde(F) = F + (1 - F)* (rho2 / rho1)
 *
 *   (2)  rho_tilde / mu_tilde = F
 *           + (1 - F)* (rho2 * mu1) / (rho1 * mu2)
 *
 * so that:
 *
 *   mu_tilde(F) = rho_tilde(F) / [ F + (1-F)*(rho2*mu1)/(rho1*mu2) ].
 *
 * We store:
 *   vof_density->rho[k][j][i]    = rho_tilde
 *   vof_viscosity->mu[k][j][i]   = mu_tilde
 *
 * where F is the volume fraction in vof->F.
 * This routine is typically called after VoF advection, 
 * ensuring that rho and mu fields remain consistent with the new interface.
 ******************************************************************************/
void VOF_update_density_viscosity(Cart3d_bag *data_bag)
{
    /**************************************************************************
     * 1. Basic references
     **************************************************************************/
    MAC_grid       *grid           = data_bag->grid;
    Parameters     *params         = data_bag->params;
    VolumeFraction *vof            = data_bag->vof;
    //Density_vof    *vof_density   = data_bag->vof_density;
    //Viscosity_vof  *vof_viscosity = data_bag->vof_viscosity;

    // Extract dimensionless volume fraction field
    double ***F = vof->F;

    // Destination arrays for dimensionless density and viscosity
    double ***rho_tilde = vof->rho;
    double ***mu_tilde  = vof->mu;

    // Physical property ratios (from parameters)
    double rho1 = params->rho1;
    double rho2 = params->rho2;
    double mu1  = params->mu1;
    double mu2  = params->mu2;

    // Local domain extents
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    /**************************************************************************
     * 2. Loop over interior cells, compute rho_tilde and mu_tilde
     **************************************************************************/
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {

                // Volume fraction in cell (i,j,k)
                double F_ijk = F[k][j][i];

                // (1) Compute dimensionless density: 
                //     rho_tilde = F + (1-F)*(rho2/rho1)
                double rhoVal = F_ijk + (1.0 - F_ijk)*(rho2 / rho1);

                // (2) Compute ratio [rho_tilde / mu_tilde] 
                //     = F + (1-F)*[ (rho2*mu1)/(rho1*mu2 ) ]
                double rmRatio = F_ijk 
                               + (1.0 - F_ijk)*((rho2*mu1)/(rho1*mu2));

                // (3) => mu_tilde = rho_tilde / rmRatio
                // Be mindful that rmRatio can be near 0 => add tiny eps if needed
                double eps = 1e-30;
                double muVal = (rmRatio > eps) ? (rhoVal / rmRatio) : 0.0;

                // Store dimensionless fields
                rho_tilde[k][j][i] = rhoVal;
                mu_tilde[k][j][i]  = muVal;
            }
        }
    }
}




#endif // VOF_PLIC