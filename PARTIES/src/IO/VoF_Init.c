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
#include "VolumeFraction.h"  

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



#ifdef VOF_PLIC




/******************************************************************************/
/*
 This function initializes the volume fraction field for a bubble, based on 
 Van SintAnnaland (2005) standard advection test case. 
 */
/******************************************************************************/
void VoF_init_bubble(Cart3d_bag *data_bag)
{
    MAC_grid    *grid   = data_bag->grid;
    Parameters  *params = data_bag->params;
    VolumeFraction *vof = data_bag->vof;
    double ***F = vof->F;

    // Bubble parameters (Van SintAnnaland 2005 test case)
    const double x_c = 0.5, y_c = 0.6875, z_c = 0.5;
    const double R = 0.1875;
    
    // Grid parameters
    const double dx = grid->dx_c[0];  // Uniform grid assumed
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];
    const double cell_vol = dx * dy * dz;

    // Subgrid sampling parameters (5x5x5 subcells per grid cell)
    const int samples = 5;
    const double sample_weight = 1.0/(samples*samples*samples);

    // Local domain indices
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; k++) {
        double z_min = grid->zc[k] - 0.5*dz;
        for (int j = Js; j < Je; j++) {
            double y_min = grid->yc[j] - 0.5*dy;
            for (int i = Is; i < Ie; i++) {
                double x_min = grid->xc[i] - 0.5*dx;
                
                // Subgrid sampling to approximate exact volume fraction
                double vol = 0.0;
                for (int sk = 0; sk < samples; sk++) {
                    double z_sub = z_min + (sk + 0.5)*(dz/samples);
                    for (int sj = 0; sj < samples; sj++) {
                        double y_sub = y_min + (sj + 0.5)*(dy/samples);
                        for (int si = 0; si < samples; si++) {
                            double x_sub = x_min + (si + 0.5)*(dx/samples);
                            
                            // Distance from bubble center
                            double dx2 = (x_sub - x_c)*(x_sub - x_c);
                            double dy2 = (y_sub - y_c)*(y_sub - y_c);
                            double dz2 = (z_sub - z_c)*(z_sub - z_c);
                            
                            if (dx2 + dy2 + dz2 < R*R) vol += 1.0;
                        }
                    }
                }
                
                F[k][j][i] = 1.0 - vol * sample_weight;  // Gas phase = 0.0
            }
        }
    }
}




/******************************************************************************/
/*
 This function initializes the volume fraction field for a bubble, for 
 the stationary droplet testcase
 */
/******************************************************************************/
void VoF_init_stationary_droplet(Cart3d_bag *data_bag)
{
    MAC_grid    *grid   = data_bag->grid;
    VolumeFraction *vof = data_bag->vof;
    double ***F = vof->F;

    // Droplet parameters for stationary droplet testcase
    const double x_c = 2.0, y_c = 2.0, z_c = 2.0; 
    const double R = 0.5; 
    
    // Grid parameters
    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    // Subgrid sampling parameters
    const int samples = 5;
    const double sample_weight = 1.0 / (samples * samples * samples);

    // Local domain indices
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k <= Ke; k++) {
        double z_min = grid->zc[k] - 0.5 * dz;
        for (int j = Js; j <= Je; j++) {
            double y_min = grid->yc[j] - 0.5 * dy;
            for (int i = Is; i <= Ie; i++) {
                double x_min = grid->xc[i] - 0.5 * dx;
                
                double vol = 0.0;
                for (int sk = 0; sk < samples; sk++) {
                    double z_sub = z_min + (sk + 0.5) * (dz / samples);
                    for (int sj = 0; sj < samples; sj++) {
                        double y_sub = y_min + (sj + 0.5) * (dy / samples);
                        for (int si = 0; si < samples; si++) {
                            double x_sub = x_min + (si + 0.5) * (dx / samples);
                            
                            // Distance from droplet center
                            double dx2 = (x_sub - x_c) * (x_sub - x_c);
                            double dy2 = (y_sub - y_c) * (y_sub - y_c);
                            double dz2 = (z_sub - z_c) * (z_sub - z_c);
                            
                            if (dx2 + dy2 + dz2 < R * R) vol += 1.0;
                        }
                    }
                }
                
                F[k][j][i] = 1 - vol * sample_weight; 
            }
        }
    }
}

// void VoF_init_stationary_droplet(Cart3d_bag *data_bag)
// {
//     MAC_grid       *grid = data_bag->grid;
//     Parameters     *p    = data_bag->params;
//     VolumeFraction *vof  = data_bag->vof;

//     // 1) your local cell‐center bounds (no ghosts)
//     int Ks = grid->G_Ks, Ke = grid->G_Ke;
//     int Js = grid->G_Js, Je = grid->G_Je;
//     int Is = grid->G_Is, Ie = grid->G_Ie;

//     int ncz = Ke - Ks + 1;
//     int ncy = Je - Js + 1;
//     int ncx = Ie - Is + 1;

//     // 2) droplet geometry
//     const double D  = 0.8, R2 = 0.25*D*D;
//     const double xc = 0.5;
//     const double yc = 0.5;
//     const double zc = 0.5;

//     // 3) build a local Φ[ ncx+1 ][ ncy+1 ][ ncz+1 ] using grid->xc,yc,zc
//     int nvx = ncx + 1, nvy = ncy + 1, nvz = ncz + 1;
//     double ***Phi = malloc(nvx*sizeof(double**));
//     for (int i = 0; i < nvx; i++) {
//       Phi[i] = malloc(nvy*sizeof(double*));
//       for (int j = 0; j < nvy; j++)
//         Phi[i][j] = calloc(nvz, sizeof(double));
//     }

//     for (int kk = 0; kk < nvz; kk++) {
//       int kg = Ks + kk;
//       double z = grid->zc[kg];
//       for (int jj = 0; jj < nvy; jj++) {
//         int jg = Js + jj;
//         double y = grid->yc[jg];
//         for (int ii = 0; ii < nvx; ii++) {
//           int ig = Is + ii;
//           double x = grid->xc[ig];
//           Phi[ii][jj][kk] = R2
//             - ((x - xc)*(x - xc)
//              + (y - yc)*(y - yc)
//              + (z - zc)*(z - zc));
//         }
//       }
//     }

//     // 4) allocate a small contiguous buffer for c_local
//     double *c_data = calloc(ncx*ncy*ncz, sizeof(double));
//     double ***c_local = malloc(ncx*sizeof(double**));
//     for (int i = 0; i < ncx; i++) {
//       c_local[i] = malloc(ncy*sizeof(double*));
//       for (int j = 0; j < ncy; j++)
//         c_local[i][j] = c_data + (i*ncy + j)*ncz;
//     }

//     // 5) run the Basilisk‐style converter on *your* block
//     fractionsLevelSet3D(
//       Phi,       // vertex level‐set
//       c_local,   // will fill c_local[i][j][k]
//       NULL,NULL,NULL,
//       0.0,
//       ncx, ncy, ncz
//     );

//     // 6) copy back into vof->F[k][j][i], looping exactly as your old code did
//     for (int kk = 0; kk < ncz; kk++) {
//       int k = Ks + kk;
//       for (int jj = 0; jj < ncy; jj++) {
//         int j = Js + jj;
//         for (int ii = 0; ii < ncx; ii++) {
//           int i = Is + ii;
//           vof->F[k][j][i] = c_local[ii][jj][kk];
//         }
//       }
//     }

//     // 7) free everything
//     free(c_data);
//     for (int i = 0; i < ncx; i++)
//       free(c_local[i]);
//     free(c_local);

//     for (int i = 0; i < nvx; i++) {
//       for (int j = 0; j < nvy; j++)
//         free(Phi[i][j]);
//       free(Phi[i]);
//     }
//     free(Phi);

//     // 8) one ghost‐exchange so F[k][j][i] is correct everywhere
//     Communication_update_ghost_nodes_flow_variable(
//       vof->F, VOLUME_FRACTION, p->ghost_nodes, data_bag
//     );
// }




/******************************************************************************/  
/*  
  Initialize the VOF field to an ellipsoidal bubble via subgrid sampling.  
  Semi-axes: a (in x), b (in y), c (in z).  
*/  
/******************************************************************************/
void VoF_init_ellipsoid(Cart3d_bag *data_bag)
{
    MAC_grid        *grid   = data_bag->grid;
    VolumeFraction  *vof    = data_bag->vof;
    double         ***F     = vof->F;

    // ------------------------------------------------------------------------
    // Ellipsoid parameters (customize or read from params)
    const double x_c = 0.5, y_c = 0.5, z_c = 0.5;   // centre
    const double a   = 0.3, b   = 0.4, c   = 0.2;   // semi-axes
    // ------------------------------------------------------------------------

    // Grid spacing
    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    // Sub-cell sampling
    const int    samples       = 5;
    const double sample_weight = 1.0/(samples*samples*samples);

    // Local indices (including ghost padding)
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    // Loop over each cell
    for (int k = Ks; k <= Ke; k++) {
        double z0 = grid->zc[k] - 0.5*dz;
        for (int j = Js; j <= Je; j++) {
            double y0 = grid->yc[j] - 0.5*dy;
            for (int i = Is; i <= Ie; i++) {
                double x0 = grid->xc[i] - 0.5*dx;
                double vol = 0.0;

                // Monte-Carlo / uniform subgrid sampling
                for (int sk = 0; sk < samples; sk++) {
                    double z_sub = z0 + (sk + 0.5)*(dz/samples);
                    double dz2   = (z_sub - z_c)*(z_sub - z_c)/(c*c);
                    for (int sj = 0; sj < samples; sj++) {
                        double y_sub = y0 + (sj + 0.5)*(dy/samples);
                        double dy2   = (y_sub - y_c)*(y_sub - y_c)/(b*b);
                        for (int si = 0; si < samples; si++) {
                            double x_sub = x0 + (si + 0.5)*(dx/samples);
                            double dx2   = (x_sub - x_c)*(x_sub - x_c)/(a*a);

                            // inside ellipsoid if dx2 + dy2 + dz2 < 1
                            if (dx2 + dy2 + dz2 < 1.0) 
                                vol += 1.0;
                        }
                    }
                }

                // assign volume fraction
                F[k][j][i] = 1 - vol * sample_weight;
            }
        }
    }
}


/******************************************************************************/
/*
 This function initializes the volume fraction field for a bubble, for 
 the stationary droplet testcase of da Silva et al. MDPI'23
 */
/******************************************************************************/
void VoF_init_rising_bubble(Cart3d_bag *data_bag)
{
    MAC_grid    *grid   = data_bag->grid;
    VolumeFraction *vof = data_bag->vof;
    double ***F = vof->F;

    // Droplet parameters for rising droplet testcase
    const double x_c = 0.5, y_c = 0.5, z_c = 0.5; 
    const double R = 0.25; 
    
    // Grid parameters
    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    // Subgrid sampling parameters
    const int samples = 5;
    const double sample_weight = 1.0 / (samples * samples * samples);

    // Local domain indices
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; k++) {
        double z_min = grid->zc[k] - 0.5 * dz;
        for (int j = Js; j < Je; j++) {
            double y_min = grid->yc[j] - 0.5 * dy;
            for (int i = Is; i < Ie; i++) {
                double x_min = grid->xc[i] - 0.5 * dx;
                
                double vol = 0.0;
                for (int sk = 0; sk < samples; sk++) {
                    double z_sub = z_min + (sk + 0.5) * (dz / samples);
                    for (int sj = 0; sj < samples; sj++) {
                        double y_sub = y_min + (sj + 0.5) * (dy / samples);
                        for (int si = 0; si < samples; si++) {
                            double x_sub = x_min + (si + 0.5) * (dx / samples);
                            
                            // Distance from droplet center
                            double dx2 = (x_sub - x_c) * (x_sub - x_c);
                            double dy2 = (y_sub - y_c) * (y_sub - y_c);
                            double dz2 = (z_sub - z_c) * (z_sub - z_c);
                            
                            if (dx2 + dy2 + dz2 < R * R) vol += 1.0;
                        }
                    }
                }
                
                F[k][j][i] = 1 - vol * sample_weight; 
            }
        }
    }
}

/******************************************************************************
 * VoF_init_two_bubbles_coaxial
 *
 * Initialise two equal, coaxial, spherical bubbles whose centres are
 * separated by three radii (3 R).  Either sphere alone reproduces the
 * single-bubble routine you pasted; the union of both gives the new test.
 *
 *  ───────────── Inputs hard-wired here ─────────────
 *      R        : bubble radius  (non-dimensional)
 *      x_c, y_c : common axial coordinates of both centres
 *      z_c1     : lower-bubble centre
 *      z_c2     : upper-bubble centre  (z_c1 + 3 R)
 *
 *  Feel free to expose these as run-time parameters later; for now they are
 *  constants so that the test case is 100 % reproducible.
 *
 *  ───────────── Conventions ─────────────
 *      • F = 1   → continuous liquid
 *      • F = 0   → gas (bubble interior)
 *      • Union of the two spheres is flagged as gas.
 ******************************************************************************/
void VoF_init_two_bubbles_coaxial(Cart3d_bag *data_bag)
{
    MAC_grid        *grid = data_bag->grid;
    VolumeFraction  *vof  = data_bag->vof;
    double ***F          = vof->F;

    /* ---------- geometry of the pair ---------- */
    const double R   = 0.25;           /* same R you used for a single bubble   */
    const double sep = 3.0 * R;        /* centre-to-centre spacing              */
    const double x_c = 0.5;
    const double y_c = 0.5;
    const double z_c1 = 0.25;          /* lower bubble                          */
    const double z_c2 = z_c1 + sep;    /* upper bubble                          */

    /* ---------- cell metrics ---------- */
    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    /* ---------- sub-cell sampling ---------- */
    const int samples = 5;
    const double w = 1.0 / (samples*samples*samples);

    /* ---------- local domain extents ---------- */
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k <= Ke; k++) {
        double z_min = grid->zc[k] - 0.5*dz;
        for (int j = Js; j <= Je; j++) {
            double y_min = grid->yc[j] - 0.5*dy;
            for (int i = Is; i <= Ie; i++) {
                double x_min = grid->xc[i] - 0.5*dx;

                double vol = 0.0;           /* gas volume in this cell */

                /* triple-nested sub-sampling */
                for (int sk = 0; sk < samples; sk++) {
                    double z_sub = z_min + (sk + 0.5)*(dz/samples);
                    double dz1 = z_sub - z_c1;
                    double dz2 = z_sub - z_c2;

                    for (int sj = 0; sj < samples; sj++) {
                        double y_sub = y_min + (sj + 0.5)*(dy/samples);
                        double dy2 = (y_sub - y_c)*(y_sub - y_c);

                        for (int si = 0; si < samples; si++) {
                            double x_sub = x_min + (si + 0.5)*(dx/samples);
                            double dx2 = (x_sub - x_c)*(x_sub - x_c);

                            /* inside lower or upper sphere? */
                            if (dx2 + dy2 + dz1*dz1 < R*R ||
                                dx2 + dy2 + dz2*dz2 < R*R)
                                vol += 1.0;
                        }
                    }
                }
                /* F = 1 outside gas, 0 in gas; use anti-aliasing fraction */
                F[k][j][i] = 1.0 - vol*w;
            }
        }
    }
}




/******************************************************************************
 * VoF_init_stationary_cylinder
 *
 * Initialize the volume-fraction field for a stationary cylinder of radius R
 * whose axis is parallel to the z-direction and whose centreline passes
 * through (x_c , y_c).  Cells inside the cylinder take F = 0 (liquid-1),
 * outside take F = 1 (liquid-2), matching the convention used for the sphere.
 ******************************************************************************/
void VoF_init_stationary_cylinder (Cart3d_bag *data_bag)
{
    MAC_grid        *grid = data_bag->grid;
    VolumeFraction  *vof  = data_bag->vof;
    double        ***F    = vof->F;

    /* Cylinder parameters -------------------------------------------------- */
    const double x_c = 2.0, y_c = 2.0;   /* axis position (m)            */
    const double R   = 0.5;              /* radius (m)                   */

    /* Grid spacings (assumed uniform per direction) ------------------------ */
    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    /* Sub-cell Monte-Carlo sampling --------------------------------------- */
    const int    samples        = 5;
    const double sample_weight  = 1.0 / (samples * samples * samples);

    /* Local domain indices (ghost-free region) ----------------------------- */
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k <= Ke; k++) {
        double z_min = grid->zc[k] - 0.5 * dz;

        for (int j = Js; j <= Je; j++) {
            double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i <= Ie; i++) {
                double x_min = grid->xc[i] - 0.5 * dx;

                double vol = 0.0;

                /* --- sub-cell samples ------------------------------------ */
                for (int sk = 0; sk < samples; sk++) {
                    double z_sub = z_min + (sk + 0.5) * (dz / samples);

                    for (int sj = 0; sj < samples; sj++) {
                        double y_sub = y_min + (sj + 0.5) * (dy / samples);

                        for (int si = 0; si < samples; si++) {
                            double x_sub = x_min + (si + 0.5) * (dx / samples);

                            /* Radial distance in x–y plane ---------------- */
                            double dx2 = (x_sub - x_c) * (x_sub - x_c);
                            double dy2 = (y_sub - y_c) * (y_sub - y_c);

                            if (dx2 + dy2 < R * R) vol += 1.0;
                        }
                    }
                }

                /* F = 0 inside cylinder, 1 outside ------------------------ */
                F[k][j][i] = 1.0 - vol * sample_weight;
            }
        }
    }
}

/******************************************************************************/
/*
 * VoF_init_vertical_bilayer_Z
 * ---------------------------
 * Fill the bottom half of the box ( z < z_cut ) with heavy fluid (F = 1)
 * and the top half ( z > z_cut ) with light fluid (F = 0).
 *
 *  • Gravity is assumed to act in the −z direction (grav = {0, 0, −1}).
 *  • The interface is placed exactly at the geometric mid-height:
 *        z_cut = 0.5 * (zmin + zmax)
 *  • Anti-alias the interface on coarse grids with a 5×5 sub-sampling
 *    in the x–y plane (same method as your stationary-droplet initialiser).
 *
 *  Ghost layers are **not** touched; call your usual boundary routine
 *  right after this initialisation.
 */
 /******************************************************************************/
void VoF_init_vertical_bilayer_Z (Cart3d_bag *data_bag)
{
    MAC_grid       *grid = data_bag->grid;
    VolumeFraction *vof  = data_bag->vof;
    double       ***F    = vof->F;

    /* -------------------------------------------------------------------- */
    const double z_cut = 0.5 * ( grid->zc[0] + grid->zc[ grid->NZ - 1 ] );

    const int    samples = 5;                       /* sub-samples per cell */
    const double w_samp  = 1.0 / (samples * samples);  /* 2-D weight (x–y) */

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    /* local interior indices (excluding ghosts) -------------------------- */
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    /* -------------------------------------------------------------------- */
    for (int k = Ks; k <= Ke; ++k) {
        double z_min = grid->zc[k] - 0.5 * dz;
        double z_max = z_min + dz;

        for (int j = Js; j <= Je; ++j) {
            double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i <= Ie; ++i) {
                double x_min = grid->xc[i] - 0.5 * dx;

                /* fast paths: fully heavy or fully light ------------------ */
                if (z_max <= z_cut) {           /* cell completely below cut */
                    F[k][j][i] = 1.0;
                    continue;
                }
                if (z_min >= z_cut) {           /* cell completely above cut */
                    F[k][j][i] = 0.0;
                    continue;
                }

                /* interface crosses this cell → sub-sample x–y plane ------ */
                double heavy_frac = 0.0;

                for (int sj = 0; sj < samples; ++sj) {
                    double y_sub = y_min + (sj + 0.5) * (dy / samples);

                    for (int si = 0; si < samples; ++si) {
                        double x_sub = x_min + (si + 0.5) * (dx / samples);

                        /* local column bounds -------------------------------- */
                        double z_sub_min = z_min;
                        double z_sub_max = z_max;

                        if (z_sub_max <= z_cut) {
                            heavy_frac += 1.0;              /* fully heavy */
                        }
                        else if (z_sub_min < z_cut) {       /* partial cell */
                            heavy_frac += (z_cut - z_sub_min) / dz;
                        }
                        /* else (fully light) add 0                                  */
                    }
                }

                F[k][j][i] = heavy_frac * w_samp;          /* 0 ≤ F ≤ 1 */
            }
        }
    }
}

/******************************************************************************/
/*
 * VoF_init_uniform_light_Z
 * ------------------------
 * Initialise the entire computational domain with the light fluid:
 *
 *        F = 0   (light phase)   for every interior cell.
 *
 *  • Designed for non-dimensional VOF runs where the light phase
 *    has density ρ̃ = ρ2/ρ1 and viscosity μ̃ = μ2/μ1.
 *  • Gravity direction and sampling are irrelevant here, so we skip
 *    all sub–cell calculations for speed.
 *
 *  Ghost layers are **not** modified; invoke your usual VOF boundary
 *  routine immediately after this function.
 */
 /******************************************************************************/
void VoF_init_uniform_light_Z (Cart3d_bag *data_bag)
{
    MAC_grid       *grid = data_bag->grid;
    VolumeFraction *vof  = data_bag->vof;
    double       ***F    = vof->F;

    /* local interior indices (excluding ghosts) -------------------------- */
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    /* -------------------------------------------------------------------- */
    for (int k = Ks; k <= Ke; ++k)
        for (int j = Js; j <= Je; ++j)
            for (int i = Is; i <= Ie; ++i)
                F[k][j][i] = 0.0;   /* light phase everywhere */
}



#endif // VOF_PLIC