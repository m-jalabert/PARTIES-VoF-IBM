/**
 * @file   VOF_InterfaceReconstruction.c
 * @brief  Geometric VOF interface normals and reconstruction (PLIC) 
 *         adapted from Basilisk's methods, integrated into PARTIES.
 *
 * This file provides *only* the routines for:
 *   1) Computing interface normals from the volume fraction field.
 *   2) Reconstructing the plane intercept (alpha) for Piecewise Linear
 *      Interface Construction (PLIC).
 *
 * We base it on Basilisk's vof.h, fractions.h, geometry.h, and myc2d/myc.h
 * (for Mixed-Youngs-Centered approximations). It is presented here in a 
 * single file for PARTIES' style, but in Basilisk these routines are 
 * distributed in multiple files.
 *
 * NOTE:
 *  - We assume a 3D domain. If you are working in 2D, you can adapt 
 *    the code similarly (see Basilisk's myc2d.h).
 *  - We do not implement VOF advection or curvature here. Only 
 *    interface normal calculation and alpha reconstruction are included.
 *  - Many advanced details (adaptive quadtree/octree, face fractions, 
 *    BFS-based stepping, etc.) are omitted or simplified for demonstration.
 *
 * USAGE IN PARTIES:
 *  - The user calls:
 *     VoF_compute_normals(data_bag);
 *     VoF_reconstruct_interface(data_bag);
 *    after having allocated the VolumeFraction arrays (F, normal_x, 
 *    normal_y, normal_z, alpha, etc.).
 *
 *  - #include "VolumeFraction.h" must declare the prototypes:
 *       void VoF_compute_normals(Cart3d_bag *data_bag);
 *       void VoF_reconstruct_interface(Cart3d_bag *data_bag);
 *    and any geometry/utility functions needed.
 *
 * This file is compiled if VOF_PLIC is defined (see Boundary.h).
 */

#include "Boundary.h"         // For VOF_PLIC macro
#include "VolumeFraction.h"   // VolumeFraction struct & function prototypes
#include "definitions.h"
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
#include <assert.h>
#include <time.h>

#ifdef VOF_PLIC

/* ============================================================================
 *  Basic geometric / "Basilisk-inspired" functions
 *  (plane_alpha, plane_volume, etc.)
 *
 *  Basilisk reference: geometry.h/fractions.h.
 *  We only include the minimal subset required for normal computation 
 *  and reconstruction (PLIC).
 * ============================================================================
 */

/**
 * @brief Compute a dimensionless plane intercept alpha (3D).
 *
 * In Basilisk, `plane_alpha()` is used to get the intercept for the 
 * plane n.x = alpha, given the volume fraction c in the unit cell and
 * a normalized normal n. This ensures the plane encloses volume c 
 * inside the cell.
 *
 * @param c     Volume fraction (clamped in [0,1]).
 * @param n     Normalized vector n = (nx, ny, nz).
 * @return      The plane intercept alpha for n.(x - 0.5)= (some shift).
 */
static double plane_alpha_3d(double c, double nx, double ny, double nz)
{
  // In Basilisk, this is a fairly elaborate routine to ensure 
  // consistent geometry. The fully accurate code is ~80 lines. 
  // Here is a simplified version for demonstration:

  // 1) clamp c
  if (c <= 0.) return -0.5*(nx + ny + nz);
  if (c >= 1.) return  0.5*(nx + ny + nz);

  // 2) We assume nx, ny, nz is already normalized
  // 3) For a truly correct approach, Basilisk solves for alpha 
  //    by subdividing or using the plane_volume() approach. 
  //    We do a naive linear guess:
  double alpha = (2.0*c - 1.0) * 0.5 * (nx + ny + nz);

  return alpha;
}

/* ============================================================================
 *  Mixed-Youngs-Centered normal approximation
 *
 *  Basilisk references:
 *    - myc.h (for dimension==3)
 *    - myc2d.h (for dimension==2)
 *
 *  The function "mycs()" or "myc()" typically returns an approximate 
 *  interface normal by evaluating volume fractions in a small stencil.
 *  We replicate a simple version for 3D. 
 *
 *  If you prefer a simpler approach (e.g. central difference of F),
 *  you can skip this. However, Basilisk uses MYC for more accuracy. 
 * ============================================================================
 */

/**
 * @brief Mixed-Youngs-Centered normal in 3D, adapted from Basilisk's myc.h
 * 
 * Basilisk uses a more detailed approach to compute partial face 
 * fractions around the cell. For demonstration, we'll implement a
 * simplified version. 
 *
 * @param F   3D array of volume fractions
 * @param i,j,k  Cell indices
 * @return   Approximate normal as (nx,ny,nz), not necessarily normalized
 */
static inline void myc_3d(double ***F, int i, int j, int k,
                          double *nx, double *ny, double *nz,
                          double dx, double dy, double dz,
                          int Nx, int Ny, int Nz)
{
  // We do a "gradient-based" approach with some weighting 
  // (some lines referencing Basilisk's "Mixed-Youngs" approach).
  // This code is shorter than Basilisk's official myc code, but 
  // captures the concept of a "stencil-based" gradient for the interface.

  // Check bounds 
  // (In Basilisk, we handle adaptivity and edges more carefully.)
  int i_minus = (i>0     ) ? i-1 : i;
  int i_plus  = (i<Nx-1 ) ? i+1 : i;
  int j_minus = (j>0     ) ? j-1 : j;
  int j_plus  = (j<Ny-1 ) ? j+1 : j;
  int k_minus = (k>0     ) ? k-1 : k;
  int k_plus  = (k<Nz-1 ) ? k+1 : k;

  // Weighted difference as a placeholder for MYC logic:
  double dFx = (F[k][j][i_plus] - F[k][j][i_minus]) / 2.0;
  double dFy = (F[k][j_plus][i] - F[k][j_minus][i]) / 2.0;
  double dFz = (F[k_plus][j][i] - F[k_minus][j][i]) / 2.0;

  // Scale by cell size for dimension consistency
  dFx /= dx; 
  dFy /= dy;
  dFz /= dz;

  *nx = dFx;
  *ny = dFy;
  *nz = dFz;
}

/* ============================================================================
 *  1) Compute interface normals
 *
 *  PARTIES style call: VoF_compute_normals(Cart3d_bag *data_bag)
 *
 *  We'll use the Basilisk-inspired MYC routine above, then normalize.
 * ============================================================================
 */
void VoF_compute_normals(Cart3d_bag *data_bag)
{
    VolumeFraction *vof = data_bag->vof;
    MAC_grid       *grid = data_bag->grid;

    int Nx = grid->Nx, Ny = grid->Ny, Nz = grid->Nz;
    double dx = grid->dx, dy = grid->dy, dz = grid->dz;

    double ***F  = vof->F;
    double ***nx = vof->normal_x;
    double ***ny = vof->normal_y;
    double ***nz = vof->normal_z;

    // For interior cells, compute normal using myc_3d
    // For boundary, you can replicate or do a fallback
    for (int k = 1; k < Nz-1; k++) {
      for (int j = 1; j < Ny-1; j++) {
        for (int i = 1; i < Nx-1; i++) {
          double gx, gy, gz;
          myc_3d(F, i, j, k, &gx, &gy, &gz, dx, dy, dz, Nx, Ny, Nz);
          // Now normalize
          double mag = sqrt(gx*gx + gy*gy + gz*gz) + 1e-15;
          nx[k][j][i] = gx / mag;
          ny[k][j][i] = gy / mag;
          nz[k][j][i] = gz / mag;
        }
      }
    }

    // For boundary cells or corners, you can do extension or copy from 
    // nearest interior. Here we do a simple approach:
    for (int k = 0; k < Nz; k++) {
      for (int j = 0; j < Ny; j++) {
        // left/right boundary
        nx[k][j][0]      = nx[k][j][1];
        nx[k][j][Nx-1]   = nx[k][j][Nx-2];
        ny[k][j][0]      = ny[k][j][1];
        ny[k][j][Nx-1]   = ny[k][j][Nx-2];
        nz[k][j][0]      = nz[k][j][1];
        nz[k][j][Nx-1]   = nz[k][j][Nx-2];
      }
    }
    for (int k = 0; k < Nz; k++) {
      for (int i = 0; i < Nx; i++) {
        // bottom/top boundary in y
        nx[k][0][i]      = nx[k][1][i];
        nx[k][Ny-1][i]   = nx[k][Ny-2][i];
        ny[k][0][i]      = ny[k][1][i];
        ny[k][Ny-1][i]   = ny[k][Ny-2][i];
        nz[k][0][i]      = nz[k][1][i];
        nz[k][Ny-1][i]   = nz[k][Ny-2][i];
      }
    }
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        // front/back boundary in z
        nx[0][j][i]      = nx[1][j][i];
        nx[Nz-1][j][i]   = nx[Nz-2][j][i];
        ny[0][j][i]      = ny[1][j][i];
        ny[Nz-1][j][i]   = ny[Nz-2][j][i];
        nz[0][j][i]      = nz[1][j][i];
        nz[Nz-1][j][i]   = nz[Nz-2][j][i];
      }
    }
}

/* ============================================================================
 *  2) Reconstruct interface (PLIC)
 *
 *  We'll fill vof->alpha in each cell with the plane intercept s.t.
 *  plane volume fraction matches F[k][j][i].
 * 
 *  Basilisk approach: 
 *    - For each cell, if 0 < F < 1, we compute alpha = plane_alpha_3d(F, n).
 *    - n is the normal (nx, ny, nz). We assume n is normalized.
 * ============================================================================
 */
void VoF_reconstruct_interface(Cart3d_bag *data_bag)
{
    VolumeFraction *vof = data_bag->vof;
    MAC_grid       *grid = data_bag->grid;

    int Nx = grid->Nx, Ny = grid->Ny, Nz = grid->Nz;

    double ***F    = vof->F;
    double ***nx   = vof->normal_x;
    double ***ny   = vof->normal_y;
    double ***nz   = vof->normal_z;
    double ***alph = vof->alpha;

    for (int k = 0; k < Nz; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
          double fval = F[k][j][i];
          if (fval > 0.0 && fval < 1.0) {
            // normalized normal
            double nnx = nx[k][j][i];
            double nny = ny[k][j][i];
            double nnz = nz[k][j][i];

            // get plane intercept alpha
            double alpha_val = plane_alpha_3d(fval, nnx, nny, nnz);
            alph[k][j][i] = alpha_val;
          } else {
            // Full or empty => alpha not used
            alph[k][j][i] = 0.0;
          }
        }
      }
    }
}

#endif // VOF_PLIC
