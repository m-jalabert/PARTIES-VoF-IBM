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
#include "VOF_DIFFUSE.h"
#include "Particle.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



#ifdef VOF

#if defined(VOF_DIFFUSE) && defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)

typedef struct {
    double X[3];
    double R;
} VoFInitSolidGeom;

static double vof_init_solid_support_extra(double cn)
{
    if (cn <= 0.0)
        return 0.0;

    const double shift = sqrt(2.0) * log(19.0) * cn;
    const double denom = 2.0 * sqrt(2.0) * cn;
    const double eps = 1.0e-12;
    const double x = 1.0 - 2.0 * eps;
    const double atanh_x = 0.5 * log((1.0 + x) / (1.0 - x));
    const double extra = -shift + denom * atanh_x;

    return (extra > 0.0) ? extra : 0.0;
}

static VoFInitSolidGeom *
vof_init_collect_diffuse_solid_geometry(Cart3d_bag *data_bag,
                                        double support_extra,
                                        int *n_geom)
{
    Lagrangian *lag = data_bag->lag;
    *n_geom = 0;

    if (lag == NULL)
        return NULL;

    Particle_list *lists[2] = {
        lag->p_mobile_list,
        lag->p_fixed_list
    };

    VoFInitSolidGeom *geom = NULL;
    int n_total = 0;

    for (int list_id = 0; list_id < 2; ++list_id) {
        if (lists[list_id] == NULL)
            continue;

        int n_particles = 0;
        Particle *particles =
            Particle_collect_owned_overlaps(lists[list_id],
                                            data_bag,
                                            support_extra,
                                            -1.0,
                                            1,
                                            &n_particles);

        if (n_particles <= 0) {
            free(particles);
            continue;
        }

        VoFInitSolidGeom *new_geom =
            (VoFInitSolidGeom *)realloc(
                geom,
                (size_t)(n_total + n_particles) * sizeof(VoFInitSolidGeom));
        Memory_check_allocation(new_geom);
        geom = new_geom;

        for (int n = 0; n < n_particles; ++n) {
            geom[n_total + n].X[0] = particles[n].X[0];
            geom[n_total + n].X[1] = particles[n].X[1];
            geom[n_total + n].X[2] = particles[n].X[2];
            geom[n_total + n].R = particles[n].R;
        }

        n_total += n_particles;
        free(particles);
    }

    *n_geom = n_total;
    return geom;
}

static double vof_init_diffuse_solid_fraction(
    const VoFInitSolidGeom *geom,
    int n_geom,
    double x,
    double y,
    double z,
    double shift,
    double denom,
    double support_extra)
{
    double cs_max = 0.0;

    for (int n = 0; n < n_geom; ++n) {
        const double dx = x - geom[n].X[0];
        const double dy = y - geom[n].X[1];
        const double dz = z - geom[n].X[2];

        double dist = sqrt(dx * dx + dy * dy + dz * dz);

#if defined(TWOD_CARTESIAN) || defined(AXISYM_RZ)
        dist = sqrt(dx * dx + dy * dy);
#endif

        if (dist > geom[n].R + support_extra)
            continue;

        const double arg =
            (dist - (geom[n].R - shift)) / (denom + 1.0e-30);
        const double cs = 0.5 - 0.5 * tanh(arg);

        if (cs > cs_max)
            cs_max = cs;
    }

    if (cs_max < 0.0)
        cs_max = 0.0;
    if (cs_max > 1.0)
        cs_max = 1.0;

    return cs_max;
}

static double vof_init_clip_liquid_by_solid(double cl, double cs)
{
    double cap = 1.0 - cs;

    if (cap < 0.0)
        cap = 0.0;
    if (cap > 1.0)
        cap = 1.0;

    if (cl > cap)
        cl = cap;
    if (cl < 0.0)
        cl = 0.0;
    if (cl > 1.0)
        cl = 1.0;

    return cl;
}

#endif




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
    const double x_c = 4.0, y_c = 4.0, z_c = 4.0; 
    const double R = 1.0; 
    
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
    Parameters  *params = data_bag->params;
    VolumeFraction *vof = data_bag->vof;
    double ***F = vof->F;


    const double Lz = params->zmax - params->zmin;
    const double R = 1.0;
    const double x_c = 0.0;
    const double y_c = 1.6;
    const double z_c = 0.0;

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

    #ifdef VOF_DIFFUSE
    {
        const double cn = params->Cn;
        const double denom = 2.0 * sqrt(2.0) * cn;
        const int samples = 6;
        const double inv_samples = 1.0 / (double)samples;
        const double sample_weight = inv_samples * inv_samples * inv_samples;

        for (int k = Ks; k < Ke; ++k) {
            const double z_min = grid->zc[k] - 0.5 * dz;
            for (int j = Js; j < Je; ++j) {
                const double y_min = grid->yc[j] - 0.5 * dy;
                for (int i = Is; i < Ie; ++i) {
                    const double x_min = grid->xc[i] - 0.5 * dx;
                    double cl = 0.0;

                    /* Cell-average the continuous tanh profile so the FV field
                     * starts closer to the discrete diffuse equilibrium. Here
                     * C_L = 0 inside the gas bubble and C_L = 1 in the liquid. */
                    for (int sk = 0; sk < samples; ++sk) {
                        const double z = z_min + (sk + 0.5) * dz * inv_samples;
                        for (int sj = 0; sj < samples; ++sj) {
                            const double y = y_min + (sj + 0.5) * dy * inv_samples;
                            for (int si = 0; si < samples; ++si) {
                                const double x = x_min + (si + 0.5) * dx * inv_samples;
                                const double r = sqrt((x - x_c) * (x - x_c) +
                                                      (y - y_c) * (y - y_c) +
                                                      (z - z_c) * (z - z_c));
                                const double s = r - R;
                                cl += sample_weight * 0.5 * (1.0 + tanh(s / (denom + 1e-30)));
                            }
                        }
                    }

                    if (cl < 0.0) {
                        cl = 0.0;
                    } else if (cl > 1.0) {
                        cl = 1.0;
                    }

                    F[k][j][i] = cl;
                    vof->C_L[k][j][i] = cl;
                }
            }
        }

        VOF_DIFFUSE_set_boundary_values(vof->C_L, data_bag);
        return;
    }
    #endif

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
    const double x_c = 0.5;
    const double y_c = 0.5;
    const double z_c1 = 0.25;          /* lower bubble                          */
    const double z_c2 = 0.625;    /* upper bubble                          */

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

/******************************************************************************/
/*
 * VoF_init_rayleigh_taylor_2d
 * ---------------------------
 * Quasi-2D Rayleigh-Taylor setup extruded in z:
 *
 *   y_int(x) = y_mid + 0.1 Lx cos(2 pi (x - xmin) / Lx)
 *
 * with heavy fluid (F = 1) above the interface and light fluid (F = 0) below.
 * The profile is independent of z so a thin, uniform z-direction can be used
 * to emulate a 2D calculation with the 3D solver.
 */
/******************************************************************************/
void VoF_init_rayleigh_taylor_2d(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double Lx    = params->xmax - params->xmin;
    const double y_mid = 0.5 * (params->ymin + params->ymax);
    const double amp   = 0.1 * Lx;

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];

    const int samples = 6;
    const double inv_samples = 1.0 / (double)samples;
    const double sample_weight = inv_samples * inv_samples;

    const int Is = grid->G_Is;
    const int Ie = grid->G_Ie;
    const int Js = grid->G_Js;
    const int Je = grid->G_Je;
    const int Ks = grid->G_Ks;
    const int Ke = grid->G_Ke;

#ifdef VOF_DIFFUSE
    const double denom = 2.0 * sqrt(2.0) * params->Cn;
#endif

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i < Ie; ++i) {
                const double x_min = grid->xc[i] - 0.5 * dx;
                double f_cell = 0.0;

                for (int sj = 0; sj < samples; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy * inv_samples;

                    for (int si = 0; si < samples; ++si) {
                        const double x = x_min + (si + 0.5) * dx * inv_samples;
                        const double eta = y_mid +
                                           amp * cos(2.0 * PI * (x - params->xmin) / Lx);

#ifdef VOF_DIFFUSE
                        const double s = y - eta;
                        f_cell += sample_weight *
                                  0.5 * (1.0 + tanh(s / (denom + 1e-30)));
#else
                        if (y >= eta) {
                            f_cell += sample_weight;
                        }
#endif
                    }
                }

                F[k][j][i] = f_cell;
            }
        }
    }
}

/******************************************************************************/
/*
 * VoF_init_axisymmetric_rising_bubble_2d
 * --------------------------------------
 * Meridional r-z version of the classic axisymmetric rising-bubble benchmark.
 *
 * The code reuses x as the radial coordinate r and y as the axial coordinate z.
 * The bubble is centred on the symmetry axis at r = xmin and sits 2R above the
 * bottom boundary, matching the standard "bubble on the axis" setup from the
 * roadmap. The third storage direction remains a bookkeeping-only slab.
 *
 * Ding et al. use the same bubble geometry on two domain sizes:
 *   section 4.1 convergence: x in [0, 2R], y in [0, 4R]
 *   section 4.2.2 / Fig. 2: x in [0, 4R], y in [0, 8R]
 */
/******************************************************************************/
void VoF_init_axisymmetric_rising_bubble_2d(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double R   = 1.0;
    const double x_c = 0.0;
    const double y_c = 1.6;

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];

    const int samples = 6;
    const double inv_samples = 1.0 / (double)samples;
    const double sample_weight = inv_samples * inv_samples;

    const int Is = grid->G_Is;
    const int Ie = grid->G_Ie;
    const int Js = grid->G_Js;
    const int Je = grid->G_Je;
    const int Ks = grid->G_Ks;
    const int Ke = grid->G_Ke;

#ifdef VOF_DIFFUSE
    const double denom = 2.0 * sqrt(2.0) * params->Cn;
#endif

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i < Ie; ++i) {
                const double x_min = grid->xc[i] - 0.5 * dx;
                double f_cell = 0.0;

                for (int sj = 0; sj < samples; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy * inv_samples;

                    for (int si = 0; si < samples; ++si) {
                        const double x = x_min + (si + 0.5) * dx * inv_samples;
                        const double r = sqrt((x - x_c) * (x - x_c) +
                                              (y - y_c) * (y - y_c));

#ifdef VOF_DIFFUSE
                        const double s = r - R;
                        f_cell += sample_weight *
                                  0.5 * (1.0 + tanh(s / (denom + 1e-30)));
#else
                        if (r >= R) {
                            f_cell += sample_weight;
                        }
#endif
                    }
                }

                F[k][j][i] = f_cell;
            }
        }
    }
}


/******************************************************************************/
/*
 * VoF_init_planar_rising_bubble_2d
 *
 * Planar 2D rising-bubble benchmark in the x-y plane, extruded through a thin
 * z-direction so the 3D solver behaves quasi-2D. This matches a Cartesian
 * rising-bubble benchmark much better than the axisymmetric surrogate:
 *   domain  [0, 1] x [0, 2]
 *   radius  R = 0.25
 *   center  (0.5, 0.5)
 */
/******************************************************************************/
void VoF_init_planar_rising_bubble_2d(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double Lx  = params->xmax - params->xmin;
    const double R   = 0.25 * Lx;
    const double x_c = 0.5 * (params->xmin + params->xmax);
    const double y_c = params->ymin + 2.0 * R;

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];

    const int samples = 6;
    const double inv_samples = 1.0 / (double)samples;
    const double sample_weight = inv_samples * inv_samples;

    const int Is = grid->G_Is;
    const int Ie = grid->G_Ie;
    const int Js = grid->G_Js;
    const int Je = grid->G_Je;
    const int Ks = grid->G_Ks;
    const int Ke = grid->G_Ke;

#ifdef VOF_DIFFUSE
    const double denom = 2.0 * sqrt(2.0) * params->Cn;
#endif

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i < Ie; ++i) {
                const double x_min = grid->xc[i] - 0.5 * dx;
                double f_cell = 0.0;

                for (int sj = 0; sj < samples; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy * inv_samples;

                    for (int si = 0; si < samples; ++si) {
                        const double x = x_min + (si + 0.5) * dx * inv_samples;
                        const double r = sqrt((x - x_c) * (x - x_c) +
                                              (y - y_c) * (y - y_c));

#ifdef VOF_DIFFUSE
                        const double s = r - R;
                        f_cell += sample_weight *
                                  0.5 * (1.0 + tanh(s / (denom + 1e-30)));
#else
                        if (r >= R) {
                            f_cell += sample_weight;
                        }
#endif
                    }
                }

                F[k][j][i] = f_cell;
            }
        }
    }
}


/*******************************************************************************
 * VoF_init_liu17_sinking_cylinder_2d
 *
 * Liu et al. (2017), section 6.3: horizontal cylinder sinking from a water
 * surface in a planar 2D half-domain.
 ******************************************************************************/
void VoF_init_liu17_sinking_cylinder_2d(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    
    const double y_interface = 3.8;

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];

#ifdef VOF_DIFFUSE
    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double href = (dx < dy) ? dx : dy;
        Cn = 0.75 * href;
    }
    const double bw = 2.0 * sqrt(2.0) * Cn;
#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
    const double shift_S = sqrt(2.0) * log(19.0) * Cn;
    const double denom_S = 2.0 * sqrt(2.0) * Cn;
    const double support_extra = vof_init_solid_support_extra(Cn);
    int n_geom = 0;
    VoFInitSolidGeom *geom =
        vof_init_collect_diffuse_solid_geometry(data_bag,
                                                support_extra,
                                                &n_geom);
#endif
#endif

    const int S = 8;
    const double inv_S = 1.0 / (double)S;
    const double w_samp = inv_S * inv_S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        const double z = grid->zc[k];
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i < Ie; ++i) {
                const double x_min = grid->xc[i] - 0.5 * dx;
                double cl = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy * inv_S;
                    for (int si = 0; si < S; ++si) {
                        const double x = x_min + (si + 0.5) * dx * inv_S;
#ifdef VOF_DIFFUSE
                        const double s = y_interface - y;
                        double cl_samp =
                            0.5 * (1.0 + tanh(s / (bw + 1.0e-30)));
#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
                        const double cs_samp =
                            vof_init_diffuse_solid_fraction(geom,
                                                            n_geom,
                                                            x,
                                                            y,
                                                            z,
                                                            shift_S,
                                                            denom_S,
                                                            support_extra);
                        cl_samp =
                            vof_init_clip_liquid_by_solid(cl_samp, cs_samp);
#endif
                        cl += w_samp * cl_samp;
#else
                        if (y <= y_interface) {
                            cl += w_samp;
                        }
#endif
                    }
                }

                if (cl < 0.0) {
                    cl = 0.0;
                } else if (cl > 1.0) {
                    cl = 1.0;
                }

                F[k][j][i] = cl;
#ifdef VOF_DIFFUSE
                vof->C_L[k][j][i] = cl;
#endif
            }
        }
    }

#if defined(VOF_DIFFUSE) && defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
    free(geom);
#endif

    vof->initial_volume = -1.0;

    {
        char msg[180];
        snprintf(msg, sizeof(msg),
                 "Liu17 section 6.3 IC: water surface y=%.6f, "
                 "solid-aware C_L clipping, Cn=%.6g\n",
                 y_interface, params->Cn);
        Display_progress(params, msg);
    }
}


/*******************************************************************************
 * VoF_init_liu17_axisymmetric_sphere_impact
 *
 * Liu et al. (2017), section 6.4: axisymmetric superhydrophobic sphere impact
 * onto a water pool.  Coordinates are r = x and z = y.  The sphere is supplied
 * by p_mobile.inp; this initializer builds the initially flat water pool and
 * clips it with the same diffuse solid profile used by VOF_DIFFUSE_compute_C_S:
 *
 *   domain        0 <= r <= 4, 0 <= z <= 6
 *   water surface z = 4
 *   sphere center z = 4.5, radius R = 0.5, initially tangent to the surface
 ******************************************************************************/
void VoF_init_liu17_axisymmetric_sphere_impact(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double z_interface = 4.0;

    const double dz = grid->dy_c[0];

#ifdef VOF_DIFFUSE
    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double dr = grid->dx_c[0];
        const double href = (dr < dz) ? dr : dz;
        Cn = 0.75 * href;
    }
    const double bw = 2.0 * sqrt(2.0) * Cn;
#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
    const double shift_S = sqrt(2.0) * log(19.0) * Cn;
    const double denom_S = 2.0 * sqrt(2.0) * Cn;
    const double support_extra = vof_init_solid_support_extra(Cn);
    int n_geom = 0;
    VoFInitSolidGeom *geom =
        vof_init_collect_diffuse_solid_geometry(data_bag,
                                                support_extra,
                                                &n_geom);
#endif
#endif

    const int S = 8;
    const double inv_S = 1.0 / (double)S;
    const double dr = grid->dx_c[0];

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double z_min = grid->yc[j] - 0.5 * dz;

            for (int i = Is; i < Ie; ++i) {
                const double r_min = grid->xc[i] - 0.5 * dr;
                double sum_CL = 0.0;
                double sum_w  = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double z = z_min + (sj + 0.5) * dz * inv_S;
                    for (int si = 0; si < S; ++si) {
                        const double r = r_min + (si + 0.5) * dr * inv_S;
                        const double w_axi = fabs(r);
#ifdef VOF_DIFFUSE
                        const double s = z_interface - z;
                        double cl_samp =
                            0.5 * (1.0 + tanh(s / (bw + 1.0e-30)));
#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
                        const double cs_samp =
                            vof_init_diffuse_solid_fraction(geom,
                                                            n_geom,
                                                            r,
                                                            z,
                                                            grid->zc[k],
                                                            shift_S,
                                                            denom_S,
                                                            support_extra);
                        cl_samp =
                            vof_init_clip_liquid_by_solid(cl_samp, cs_samp);
#endif
                        sum_CL += w_axi * cl_samp;
#else
                        if (z <= z_interface)
                            sum_CL += w_axi;
#endif
                        sum_w += w_axi;
                    }
                }

                double cl = (sum_w > 1.0e-30) ? (sum_CL / sum_w) : 0.0;

                if (cl < 0.0) {
                    cl = 0.0;
                } else if (cl > 1.0) {
                    cl = 1.0;
                }

                F[k][j][i] = cl;
#ifdef VOF_DIFFUSE
                vof->C_L[k][j][i] = cl;
#endif
            }
        }
    }

#if defined(VOF_DIFFUSE) && defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
    free(geom);
#endif

    vof->initial_volume = -1.0;

    {
        char msg[180];
        snprintf(msg, sizeof(msg),
                 "Liu17 section 6.4 IC: axisymmetric water surface z=%.6f, "
                 "solid-aware C_L clipping, sphere center z=4.500000, "
                 "R=0.500000, Cn=%.6g\n",
                 z_interface, params->Cn);
        Display_progress(params, msg);
    }
}


/*******************************************************************************
 * VoF_init_liu17_self_assembly_floating_cylinders_2d
 *
 * Liu et al. (2017), section 6.6: three hydrophilic cylinders initially half
 * immersed at a flat water-air interface.  The cylinders are supplied by
 * p_mobile.inp; this initializer builds the liquid field and clips it with the
 * same diffuse solid profile used by VOF_DIFFUSE_compute_C_S:
 *
 *   domain        0 <= x <= 6, 0 <= y <= 4
 *   water surface y = 2.5
 *   liquid        y < 2.5
 *
 * VOF_DIFFUSE_init() subsequently rebuilds C_S from the mobile cylinders and
 * applies the contact-angle projection, but the t = 0 liquid field is already
 * locally admissible.
 ******************************************************************************/
void VoF_init_liu17_self_assembly_floating_cylinders_2d(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double y_interface = 2.5;

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];

#ifdef VOF_DIFFUSE
    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double href = (dx < dy) ? dx : dy;
        Cn = 0.75 * href;
    }
    const double bw = 2.0 * sqrt(2.0) * Cn;
#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
    const double shift_S = sqrt(2.0) * log(19.0) * Cn;
    const double denom_S = 2.0 * sqrt(2.0) * Cn;
    const double support_extra = vof_init_solid_support_extra(Cn);
    int n_geom = 0;
    VoFInitSolidGeom *geom =
        vof_init_collect_diffuse_solid_geometry(data_bag,
                                                support_extra,
                                                &n_geom);
#endif
#endif

    const int S = 8;
    const double inv_S = 1.0 / (double)S;
    const double w_samp = inv_S * inv_S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        const double z = grid->zc[k];
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i < Ie; ++i) {
                const double x_min = grid->xc[i] - 0.5 * dx;
                double cl = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy * inv_S;
                    for (int si = 0; si < S; ++si) {
                        const double x = x_min + (si + 0.5) * dx * inv_S;
#ifdef VOF_DIFFUSE
                        const double s = y_interface - y;
                        double cl_samp =
                            0.5 * (1.0 + tanh(s / (bw + 1.0e-30)));
#if defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
                        const double cs_samp =
                            vof_init_diffuse_solid_fraction(geom,
                                                            n_geom,
                                                            x,
                                                            y,
                                                            z,
                                                            shift_S,
                                                            denom_S,
                                                            support_extra);
                        cl_samp =
                            vof_init_clip_liquid_by_solid(cl_samp, cs_samp);
#endif
                        cl += w_samp * cl_samp;
#else
                        if (y <= y_interface) {
                            cl += w_samp;
                        }
#endif
                    }
                }

                if (cl < 0.0) {
                    cl = 0.0;
                } else if (cl > 1.0) {
                    cl = 1.0;
                }

                F[k][j][i] = cl;
#ifdef VOF_DIFFUSE
                vof->C_L[k][j][i] = cl;
#endif
            }
        }
    }

#if defined(VOF_DIFFUSE) && defined(VOF_IBM) && defined(LAG_PARTICLE_RESOLVED)
    free(geom);
#endif

    vof->initial_volume = -1.0;

    {
        char msg[180];
        snprintf(msg, sizeof(msg),
                 "Liu17 section 6.6 IC: planar water surface y=%.6f, "
                 "solid-aware C_L clipping, three cylinders initially half "
                 "immersed, Cn=%.6g\n",
                 y_interface, params->Cn);
        Display_progress(params, msg);
    }
}


/*******************************************************************************
 * VoF_droplet_flat_plate_175deg
 *
 * Initialise a spherical-cap droplet (contact angle θ = 175°) resting on the
 * wall (y = 0) whose sphere-of-curvature radius R is chosen so that the cap
 * volume equals that of a full sphere of radius R0:
 *
 *   V_cap(R,θ) = (π/3) R^3 (1 - cosθ)^2 (2 + cosθ) = (4/3) π R0^3
 *     ⇒ R = R0 [ 4 / ((1 - cosθ)^2 (2 + cosθ)) ]^{1/3}
 *
 * Geometry:
 *   • Cap height:   h = R (1 - cosθ)
 *   • Footprint:    a = R sinθ
 *   • Sphere centre: (xc0, yc0, zc0) with yc0 = -R cosθ (below the wall)
 ******************************************************************************/
void VoF_droplet_flat_plate (Cart3d_bag *bag)
{
    MAC_grid       *g   = bag->grid;
    VolumeFraction *vof = bag->vof;
    double       ***F   = vof->F;

    /* user/base parameter: R0 (paper’s initial sphere radius) */
    const double R0 = 1.0;

    /* contact angle and derived sphere-of-curvature R that preserves volume */
    const double theta = 175.0 * PI / 180.0;
    const double c     = cos(theta);
    const double denom = (1.0 - c)*(1.0 - c) * (2.0 + c);
    const double R     = R0 * pow(4.0 / denom, 1.0/3.0);
    const double R2    = R*R;

    /* sphere centre: below the wall by R cosθ; keep it centred in x,z */
    const double xc0 = 3.0;
    const double yc0 = -R * c;   /* centre below y=0 */
    const double zc0 = 3.0;

    const double EPS = 1e-14;

    /* subcell sampler */
    const int    S    = 5;
    const double invS = 1.0 / (double) S;

    /* interior extents (no ghosts) */
    const int Is = g->G_Is, Js = g->G_Js, Ks = g->G_Ks;
    const int Ie = g->G_Ie, Je = g->G_Je, Ke = g->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        const double zc = g->zc[k];
        const double dz = g->dz_c[k];
        const double hz = 0.5 * dz;

        for (int j = Js; j < Je; ++j) {
            const double yc = g->yc[j];
            const double dy = g->dy_c[j];
            const double hy = 0.5 * dy;

            const double y_min = yc - hy;
            const double y_max = yc + hy;

            /* Entire cell strictly below the wall -> gas */
            if (y_max <= 0.0) {
                for (int i = Is; i < Ie; ++i) F[k][j][i] = 0.0;
                continue;
            }

            for (int i = Is; i < Ie; ++i) {
                const double xc = g->xc[i];
                const double dx = g->dx_c[i];
                const double hx = 0.5 * dx;

                /* AABB–sphere quick tests */
                double dxp = fabs(xc0 - xc) - hx; if (dxp < 0.0) dxp = 0.0;
                double dyp = fabs(yc0 - yc) - hy; if (dyp < 0.0) dyp = 0.0;
                double dzp = fabs(zc0 - zc) - hz; if (dzp < 0.0) dzp = 0.0;
                const double dist_min2 = dxp*dxp + dyp*dyp + dzp*dzp;

                const double dxf = fabs(xc0 - xc) + hx;
                const double dyf = fabs(yc0 - yc) + hy;
                const double dzf = fabs(zc0 - zc) + hz;
                const double dist_max2 = dxf*dxf + dyf*dyf + dzf*dzf;

                const int cell_all_above = (y_min >= 0.0);

                /* Fully inside spherical cap (cell entirely above wall AND inside sphere) */
                if (cell_all_above && dist_max2 <= (R + EPS)*(R + EPS)) {
                    F[k][j][i] = 1.0;
                    continue;
                }

                /* Fully outside sphere (no overlap at all) */
                if (dist_min2 >= (R - EPS)*(R - EPS)) {
                    F[k][j][i] = 0.0;
                    continue;
                }

                /* Partial overlap: sample only the slab at y ≥ 0 */
                const double ddx = dx * invS;
                const double ddz = dz * invS;

                const double y_start = (y_min > 0.0) ? y_min : 0.0;
                const double y_span  = fmax(0.0, y_max - y_start);
                if (y_span <= 0.0) { F[k][j][i] = 0.0; continue; }
                const double ddy_eff = y_span / (double) S;

                const double x0 = (xc - hx) + 0.5*ddx;
                const double z0 = (zc - hz) + 0.5*ddz;
                const double y0 = y_start     + 0.5*ddy_eff;

                double cnt = 0.0, tot = 0.0;

                for (int sk = 0; sk < S; ++sk) {
                    const double z   = z0 + sk*ddz;
                    const double dz2 = (z - zc0)*(z - zc0);
                    for (int sj = 0; sj < S; ++sj) {
                        const double y   = y0 + sj*ddy_eff;  /* y ≥ 0 */
                        const double dy2 = (y - yc0)*(y - yc0);
                        for (int si = 0; si < S; ++si) {
                            const double x   = x0 + si*ddx;
                            const double dx2 = (x - xc0)*(x - xc0);
                            ++tot;
                            if (dx2 + dy2 + dz2 <= R2) cnt += 1.0;
                        }
                    }
                }

                double f = (tot > 0.0) ? (cnt / tot) : 0.0;
                if      (f < 0.0) f = 0.0;
                else if (f > 1.0) f = 1.0;
                F[k][j][i] = f;
            }
        }
    }
}





/*******************************************************************************
 * VoF_droplet_flat_plate_90
 *
 * Initialise a hemispherical droplet resting on the wall (y=0) whose radius
 * is Rf = 2^{1/3} * R0, consistent with Patel Eq. (24) for θ = 90°.
 *
 * Why: The post-processing uses R0 in Eq. (24). For θ=90°, the equilibrium
 * spherical cap is a hemisphere with radius-of-curvature Rf = 2^{1/3} R0.
 * Initialising the hemisphere with this Rf makes the volume equal to that of
 * a full sphere of radius R0, so the analytic profile matches and E → 0.
 ******************************************************************************/
void VoF_droplet_flat_plate_90 (Cart3d_bag *bag)
{
    MAC_grid       *g   = bag->grid;
    VolumeFraction *vof = bag->vof;
    double       ***F   = vof->F;

    /* user/base parameter: R0 (paper’s initial sphere radius) */
    const double R0  = 1.0;

    /* hemisphere radius consistent with Eq. (24) at θ=90°: Rf = (2)^{1/3} R0 */
    const double cbrt2 = pow(2.0, 1.0/3.0);
    const double R     = cbrt2 * R0;           /* <-- use this R for the hemisphere */
    const double R2    = R*R;

    /* place the hemisphere on the wall, centered in x,z */
    const double xc0 = 3.0;
    const double yc0 = 0.0;                    /* sphere center lies on the wall plane */
    const double zc0 = 3.0;

    const double EPS = 1e-14;

    /* subcell sampler */
    const int    S    = 5;
    const double invS = 1.0 / (double) S;

    /* interior extents (no ghosts) */
    const int Is = g->G_Is, Js = g->G_Js, Ks = g->G_Ks;
    const int Ie = g->G_Ie, Je = g->G_Je, Ke = g->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        const double zc = g->zc[k];
        const double dz = g->dz_c[k];
        const double hz = 0.5 * dz;

        for (int j = Js; j < Je; ++j) {
            const double yc = g->yc[j];
            const double dy = g->dy_c[j];
            const double hy = 0.5 * dy;

            const double y_min = yc - hy;
            const double y_max = yc + hy;

            /* Entire cell strictly below the wall -> gas */
            if (y_max <= 0.0) {
                for (int i = Is; i < Ie; ++i) F[k][j][i] = 0.0;
                continue;
            }

            for (int i = Is; i < Ie; ++i) {
                const double xc = g->xc[i];
                const double dx = g->dx_c[i];
                const double hx = 0.5 * dx;

                /* AABB–sphere quick tests */
                double dxp = fabs(xc0 - xc) - hx; if (dxp < 0.0) dxp = 0.0;
                double dyp = fabs(yc0 - yc) - hy; if (dyp < 0.0) dyp = 0.0;
                double dzp = fabs(zc0 - zc) - hz; if (dzp < 0.0) dzp = 0.0;
                const double dist_min2 = dxp*dxp + dyp*dyp + dzp*dzp;

                const double dxf = fabs(xc0 - xc) + hx;
                const double dyf = fabs(yc0 - yc) + hy;
                const double dzf = fabs(zc0 - zc) + hz;
                const double dist_max2 = dxf*dxf + dyf*dyf + dzf*dzf;

                const int cell_all_above = (y_min >= 0.0);

                /* Fully inside hemisphere (cell entirely above wall AND inside sphere) */
                if (cell_all_above && dist_max2 <= (R + EPS)*(R + EPS)) {
                    F[k][j][i] = 1.0;
                    continue;
                }

                /* Fully outside sphere (no overlap at all) */
                if (dist_min2 >= (R - EPS)*(R - EPS)) {
                    F[k][j][i] = 0.0;
                    continue;
                }

                /* Partial overlap: sample only the slab at y >= 0 */
                const double ddx = dx * invS;
                const double ddz = dz * invS;

                const double y_start = (y_min > 0.0) ? y_min : 0.0;
                const double y_span  = fmax(0.0, y_max - y_start);
                if (y_span <= 0.0) { F[k][j][i] = 0.0; continue; }
                const double ddy_eff = y_span / (double) S;

                const double x0 = (xc - hx) + 0.5*ddx;
                const double z0 = (zc - hz) + 0.5*ddz;
                const double y0 = y_start     + 0.5*ddy_eff;

                double cnt = 0.0, tot = 0.0;

                for (int sk = 0; sk < S; ++sk) {
                    const double z   = z0 + sk*ddz;
                    const double dz2 = (z - zc0)*(z - zc0);
                    for (int sj = 0; sj < S; ++sj) {
                        const double y   = y0 + sj*ddy_eff;  /* y ≥ 0 */
                        const double dy2 = (y - yc0)*(y - yc0);
                        for (int si = 0; si < S; ++si) {
                            const double x   = x0 + si*ddx;
                            const double dx2 = (x - xc0)*(x - xc0);
                            ++tot;
                            if (dx2 + dy2 + dz2 <= R2) cnt += 1.0;
                        }
                    }
                }

                double f = (tot > 0.0) ? (cnt / tot) : 0.0;
                if      (f < 0.0) f = 0.0;
                else if (f > 1.0) f = 1.0;
                F[k][j][i] = f;
            }
        }
    }
}





/*******************************************************************************
 * VoF_droplet_flat_plate_30deg
 *
 * Initialise a spherical-cap droplet (contact angle θ = 30°) resting on the
 * wall (y = 0) whose sphere-of-curvature radius R is chosen so that the cap
 * volume equals that of a full sphere of radius R0:
 *
 *   V_cap(R,θ) = (π/3) R^3 (1 - cosθ)^2 (2 + cosθ) = (4/3) π R0^3
 *     ⇒ R = R0 [ 4 / ((1 - cosθ)^2 (2 + cosθ)) ]^{1/3}
 *
 * Geometry:
 *   • Cap height:   h = R (1 - cosθ)
 *   • Footprint:    a = R sinθ
 *   • Sphere centre: (xc0, yc0, zc0) with yc0 = -R cosθ.
 *       For θ < 90°, cosθ > 0 ⇒ yc0 < 0 (centre lies below the wall).
 *******************************************************************************/
void VoF_droplet_flat_plate_30 (Cart3d_bag *bag)
{
    MAC_grid       *g   = bag->grid;
    VolumeFraction *vof = bag->vof;
    double       ***F   = vof->F;

    /* user/base parameter: R0 (paper’s initial sphere radius) */
    const double R0 = 1.0;

    /* contact angle and derived sphere-of-curvature R that preserves volume */
    const double theta = 30.0 * PI / 180.0;
    const double c     = cos(theta);
    const double denom = (1.0 - c)*(1.0 - c) * (2.0 + c);
    const double R     = R0 * pow(4.0 / denom, 1.0/3.0);
    const double R2    = R*R;

    /* sphere centre: yc0 = -R cosθ; for θ=30° this is < 0 (below y=0) */
    const double xc0 = 3.0;
    const double yc0 = -R * c;
    const double zc0 = 3.0;

    const double EPS = 1e-14;

    /* subcell sampler */
    const int    S    = 5;
    const double invS = 1.0 / (double) S;

    /* interior extents (no ghosts) */
    const int Is = g->G_Is, Js = g->G_Js, Ks = g->G_Ks;
    const int Ie = g->G_Ie, Je = g->G_Je, Ke = g->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        const double zc = g->zc[k];
        const double dz = g->dz_c[k];
        const double hz = 0.5 * dz;

        for (int j = Js; j < Je; ++j) {
            const double yc = g->yc[j];
            const double dy = g->dy_c[j];
            const double hy = 0.5 * dy;

            const double y_min = yc - hy;
            const double y_max = yc + hy;

            /* Entire cell strictly below the wall -> gas */
            if (y_max <= 0.0) {
                for (int i = Is; i < Ie; ++i) F[k][j][i] = 0.0;
                continue;
            }

            for (int i = Is; i < Ie; ++i) {
                const double xc = g->xc[i];
                const double dx = g->dx_c[i];
                const double hx = 0.5 * dx;

                /* AABB–sphere quick tests */
                double dxp = fabs(xc0 - xc) - hx; if (dxp < 0.0) dxp = 0.0;
                double dyp = fabs(yc0 - yc) - hy; if (dyp < 0.0) dyp = 0.0;
                double dzp = fabs(zc0 - zc) - hz; if (dzp < 0.0) dzp = 0.0;
                const double dist_min2 = dxp*dxp + dyp*dyp + dzp*dzp;

                const double dxf = fabs(xc0 - xc) + hx;
                const double dyf = fabs(yc0 - yc) + hy;
                const double dzf = fabs(zc0 - zc) + hz;
                const double dist_max2 = dxf*dxf + dyf*dyf + dzf*dzf;

                const int cell_all_above = (y_min >= 0.0);

                /* Fully inside spherical cap (cell entirely above wall AND inside sphere) */
                if (cell_all_above && dist_max2 <= (R + EPS)*(R + EPS)) {
                    F[k][j][i] = 1.0;
                    continue;
                }

                /* Fully outside sphere (no overlap at all) */
                if (dist_min2 >= (R - EPS)*(R - EPS)) {
                    F[k][j][i] = 0.0;
                    continue;
                }

                /* Partial overlap: sample only the slab at y ≥ 0 */
                const double ddx = dx * invS;
                const double ddz = dz * invS;

                const double y_start = (y_min > 0.0) ? y_min : 0.0;
                const double y_span  = fmax(0.0, y_max - y_start);
                if (y_span <= 0.0) { F[k][j][i] = 0.0; continue; }
                const double ddy_eff = y_span / (double) S;

                const double x0 = (xc - hx) + 0.5*ddx;
                const double z0 = (zc - hz) + 0.5*ddz;
                const double y0 = y_start     + 0.5*ddy_eff;

                double cnt = 0.0, tot = 0.0;

                for (int sk = 0; sk < S; ++sk) {
                    const double z   = z0 + sk*ddz;
                    const double dz2 = (z - zc0)*(z - zc0);
                    for (int sj = 0; sj < S; ++sj) {
                        const double y   = y0 + sj*ddy_eff;  /* y ≥ 0 */
                        const double dy2 = (y - yc0)*(y - yc0);
                        for (int si = 0; si < S; ++si) {
                            const double x   = x0 + si*ddx;
                            const double dx2 = (x - xc0)*(x - xc0);
                            ++tot;
                            if (dx2 + dy2 + dz2 <= R2) cnt += 1.0;
                        }
                    }
                }

                double f = (tot > 0.0) ? (cnt / tot) : 0.0;
                if      (f < 0.0) f = 0.0;
                else if (f > 1.0) f = 1.0;
                F[k][j][i] = f;
            }
        }
    }
}




/*******************************************************************************
 * VoF_init_droplet_on_sphere
 ******************************************************************************/
void VoF_init_droplet_on_sphere (Cart3d_bag *bag)
{
    MAC_grid       *g   = bag->grid;
    VolumeFraction *vof = bag->vof;
    double       ***F   = vof->F;

    /* ---------------- geometry parameters -------------------------------- */
    const double Rs = 1.0;          /* solid sphere radius                  */
    const double Rd = 1.0;          /* liquid sphere radius                 */

    const double xs = 2.0, ys = 0.5, zs = 2.0;          /* solid centre   */
    const double xd = 2.0,   yd = 2.5 - 0.03,   zd = 2.0;  /* liquid centre  */

    const double Rs2 = Rs*Rs, Rd2 = Rd*Rd;

    /* ---------------- grid & 5×5×5 sampler ------------------------------- */
    const double dx = g->dx_c[0], dy = g->dy_c[0], dz = g->dz_c[0];
    const int    S  = 5;
    const double w  = 1.0 / (S*S*S);

    int Is = g->G_Is, Ie = g->G_Ie;
    int Js = g->G_Js, Je = g->G_Je;
    int Ks = g->G_Ks, Ke = g->G_Ke;
    vof->initial_volume = -1.0; 

    for (int k = Ks; k < Ke; ++k) {
        double zmin = g->zc[k] - 0.5*dz;

        for (int j = Js; j < Je; ++j) {
            double ymin = g->yc[j] - 0.5*dy;

            for (int i = Is; i < Ie; ++i) {
                double xmin = g->xc[i] - 0.5*dx;

                double vol = 0.0;

                for (int sk = 0; sk < S; ++sk) {
                    double z = zmin + (sk+0.5)*dz/S;
                    double dzs2 = (z - zs)*(z - zs);
                    double dzd2 = (z - zd)*(z - zd);

                    for (int sj = 0; sj < S; ++sj) {
                        double y = ymin + (sj+0.5)*dy/S;
                        double dys2 = (y - ys)*(y - ys);
                        double dyd2 = (y - yd)*(y - yd);

                        for (int si = 0; si < S; ++si) {
                            double x = xmin + (si+0.5)*dx/S;
                            double dxs2 = (x - xs)*(x - xs);
                            double dxd2 = (x - xd)*(x - xd);


                            /* count if sample is inside the liquid sphere */
                            if (dxd2 + dyd2 + dzd2 < Rd2) vol += 1.0;
                        }
                    }
                }
                F[k][j][i] = vol * w;      
            }
        }
    }
}


/******************************************************************************/
/*
 * VoF_init_bilayer_at_4D
 * ----------------------
 * Fill the bottom portion of the box ( y < y_cut ) with heavy fluid (F = 1)
 * and the upper portion ( y > y_cut ) with light fluid (F = 0).
 *
 *  • Designed for the particle-laden interface simulation where y_cut = 4D_p.
 *  • With D_p = 1.0 (reference length scale = particle diameter, R = 0.5),
 *    the interface is at y = 4.0.
 *  • Gravity acts in the −y direction (grav = {0, −1, 0}).
 *  • Anti-alias the interface on coarse grids with a 5×5 sub-sampling
 *    in the x–z plane (perpendicular to gravity).
 *
 *  Ghost layers are **not** touched; call your usual boundary routine
 *  right after this initialisation.
 */
/******************************************************************************/
void VoF_init_bilayer_at_4D (Cart3d_bag *data_bag)
{
    MAC_grid       *grid = data_bag->grid;
    VolumeFraction *vof  = data_bag->vof;
    double       ***F    = vof->F;


    const double R     = 0.5;           
                   
    const double y_cut = 3.0; 

    const int    samples = 5;                          /* sub-samples per cell */
    const double w_samp  = 1.0 / (samples * samples);  /* 2-D weight (x–z)     */

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

        for (int j = Js; j <= Je; ++j) {
            double y_min = grid->yc[j] - 0.5 * dy;
            double y_max = y_min + dy;

            for (int i = Is; i <= Ie; ++i) {
                double x_min = grid->xc[i] - 0.5 * dx;

                /* fast paths: fully heavy or fully light ------------------ */
                if (y_max <= y_cut) {           /* cell completely below cut */
                    F[k][j][i] = 1.0;
                    continue;
                }
                if (y_min >= y_cut) {           /* cell completely above cut */
                    F[k][j][i] = 0.0;
                    continue;
                }

                /* interface crosses this cell → sub-sample x–z plane ------ */
                double heavy_frac = 0.0;

                for (int sk = 0; sk < samples; ++sk) {
                    double z_sub = z_min + (sk + 0.5) * (dz / samples);

                    for (int si = 0; si < samples; ++si) {
                        double x_sub = x_min + (si + 0.5) * (dx / samples);

                        /* local column bounds in y direction ----------------- */
                        double y_sub_min = y_min;
                        double y_sub_max = y_max;

                        if (y_sub_max <= y_cut) {
                            heavy_frac += 1.0;              /* fully heavy */
                        }
                        else if (y_sub_min < y_cut) {       /* partial cell */
                            heavy_frac += (y_cut - y_sub_min) / dy;
                        }
                        /* else (fully light) add 0                                  */
                    }
                }

                F[k][j][i] = heavy_frac * w_samp;          /* 0 ≤ F ≤ 1 */
            }
        }
    }
}





#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double clamp_double(double x, double a, double b)
{
    return (x < a) ? a : (x > b) ? b : x;
}

/* Intersection volume of two spheres of radii R and r with center distance d.
   Handles disjoint and containment cases robustly. */
static double sphere_sphere_intersection_volume(double R, double r, double d)
{
    /* disjoint */
    if (d >= R + r) return 0.0;

    /* one contains the other */
    if (d <= fabs(R - r)) {
        const double rmin = (R < r) ? R : r;
        return (4.0/3.0) * M_PI * rmin*rmin*rmin;
    }

    /* partial overlap: standard closed form */
    const double a = R + r - d; /* >0 */
    const double term = (d*d + 2.0*d*r - 3.0*r*r + 2.0*d*R + 6.0*r*R - 3.0*R*R);
    return (M_PI * a*a * term) / (12.0 * d);
}

/* Given Rs, theta, and a candidate Rf, compute l from contact-angle constraint:
      l^2 = Rs^2 + Rf^2 + 2 Rs Rf cos(theta)
   and return volume outside solid after clipping:
      V_out = (4/3)π Rf^3 - V_intersection(Rs, Rf, l)
*/
static double droplet_volume_outside_solid(double Rs, double Rf, double theta_rad)
{
    const double cosT = cos(theta_rad);
    const double l2   = Rs*Rs + Rf*Rf + 2.0*Rs*Rf*cosT;
    const double l    = sqrt(fmax(l2, 0.0));

    const double V_sphere = (4.0/3.0) * M_PI * Rf*Rf*Rf;
    const double V_int    = sphere_sphere_intersection_volume(Rs, Rf, l);

    return V_sphere - V_int;
}

/* Solve for Rf such that V_outside(Rf) = V0 using bisection.
   We expand the upper bracket until sign change (robust). */
static double solve_Rf_for_volume(double Rs, double theta_rad, double V0)
{
    /* Lower bound: small but >0 */
    double a = 1e-6;
    double fa = droplet_volume_outside_solid(Rs, a, theta_rad) - V0;

    /* Choose an initial upper bound near the free sphere radius */
    double b = 1.0;
    double fb = droplet_volume_outside_solid(Rs, b, theta_rad) - V0;

    /* If we already bracket, fine. Otherwise expand b until bracket. */
    int it_expand = 0;
    while (fa * fb > 0.0 && it_expand < 80) {
        b *= 1.5; /* geometric expansion */
        fb = droplet_volume_outside_solid(Rs, b, theta_rad) - V0;
        it_expand++;
    }

    /* If still not bracketed, fall back (should not happen for sane inputs). */
    if (fa * fb > 0.0) {
        /* Return something reasonable rather than NaN */
        return b;
    }

    /* Bisection */
    for (int it = 0; it < 200; ++it) {
        const double m  = 0.5 * (a + b);
        const double fm = droplet_volume_outside_solid(Rs, m, theta_rad) - V0;

        if (fabs(fm) < 1e-12 * V0) return m;

        if (fa * fm <= 0.0) {
            b  = m;
            fb = fm;
        } else {
            a  = m;
            fa = fm;
        }

        if (fabs(b - a) < 1e-12 * fmax(1.0, m)) return 0.5 * (a + b);
    }

    return 0.5 * (a + b);
}

/*******************************************************************************
 * VoF_init_droplet_on_sphere_theta
 *
 * Initializes the droplet as the OUTSIDE-SOLID part of a sphere of radius Rf,
 * where Rf is computed so that the clipped liquid volume equals a free droplet
 * of radius R0=1 (or user-set R0 if you want).
 *
 * Input:
 *   theta_deg : static contact angle through the liquid
 *******************************************************************************/
void VoF_init_droplet_on_sphere_theta(Cart3d_bag *bag)
{
    MAC_grid       *g   = bag->grid;
    VolumeFraction *vof = bag->vof;
    double       ***F   = vof->F;

    /* ---------------- geometry parameters ---------------------------------- */
    const double Rs = 1.0;          /* solid sphere radius */
    const double R0 = 1.0;          /* "starting droplet" radius -> sets volume */

    const double xs = 2.0, ys = 1.30, zs = 2.0;  /* solid center */

    double theta_deg = 180.0;              /* contact angle in degrees */

    double theta = PI - theta_deg * (M_PI / 180.0);
    theta = clamp_double(theta, 1e-12, M_PI - 1e-12); /* avoid extreme degeneracy */

    const double V0 = (4.0/3.0) * M_PI * R0*R0*R0;

    /* Solve curvature radius Rf from volume constraint + contact angle */
    const double Rf = solve_Rf_for_volume(Rs, theta, V0);

    /* Contact-angle geometry -> center distance l */
    const double l2 = Rs*Rs + Rf*Rf + 2.0*Rs*Rf*cos(theta);
    const double l  = sqrt(fmax(l2, 0.0));

    /* Place droplet center along +y from solid center */
    const double xd = xs;
    const double yd = ys + l;
    const double zd = zs;

    const double Rs2 = Rs*Rs;
    const double Rf2 = Rf*Rf;

    /* ---------------- grid & S×S×S sampler -------------------------------- */
    const double dx = g->dx_c[0], dy = g->dy_c[0], dz = g->dz_c[0];

    const int    S  = 5;
    const double w  = 1.0 / (S*S*S);

    const int Is = g->G_Is, Ie = g->G_Ie;
    const int Js = g->G_Js, Je = g->G_Je;
    const int Ks = g->G_Ks, Ke = g->G_Ke;

    vof->initial_volume = -1.0;

    for (int k = Ks; k < Ke; ++k) {
        const double zmin = g->zc[k] - 0.5*dz;

        for (int j = Js; j < Je; ++j) {
            const double ymin = g->yc[j] - 0.5*dy;

            for (int i = Is; i < Ie; ++i) {
                const double xmin = g->xc[i] - 0.5*dx;

                double vol = 0.0;

                for (int sk = 0; sk < S; ++sk) {
                    const double z   = zmin + (sk + 0.5)*dz/S;
                    const double dzd = z - zd;
                    const double dzs = z - zs;

                    for (int sj = 0; sj < S; ++sj) {
                        const double y   = ymin + (sj + 0.5)*dy/S;
                        const double dyd = y - yd;
                        const double dys = y - ys;

                        for (int si = 0; si < S; ++si) {
                            const double x   = xmin + (si + 0.5)*dx/S;
                            const double dxd = x - xd;
                            const double dxs = x - xs;

                            const double r2_d = dxd*dxd + dyd*dyd + dzd*dzd; /* droplet sphere */
                            const double r2_s = dxs*dxs + dys*dys + dzs*dzs; /* solid sphere   */

                            /* liquid only outside solid; inside droplet sphere */
                            if (r2_s >= Rs2 && r2_d < Rf2) vol += 1.0;
                        }
                    }
                }

                F[k][j][i] = vol * w;
            }
        }
    }


}




void VoF_init_all_heavy(Cart3d_bag *data_bag)
{
  MAC_grid *grid = data_bag->grid;
  VolumeFraction *vof = data_bag->vof;
  double ***F = vof->F;

  for (int k=grid->G_Ks; k<=grid->G_Ke; ++k)
  for (int j=grid->G_Js; j<=grid->G_Je; ++j)
  for (int i=grid->G_Is; i<=grid->G_Ie; ++i)
    F[k][j][i] = 1.0;
}

/******************************************************************************/
/*
 This function initializes the volume fraction field for a bubble, for 
 the stationary 3D droplet testcase of Francois et al. (2006) "A balanced-force algorithm
 */
/******************************************************************************/
void VoF_init_stationary_droplet_Francois(Cart3d_bag *data_bag)
{
    MAC_grid    *grid   = data_bag->grid;
    Parameters  *params = data_bag->params;
    VolumeFraction *vof = data_bag->vof;
    double ***F = vof->F;

    // Droplet parameters (Dimensionless)
    // Domain is [0,4], centered at 2, Radius is 1.
    const double x_c = 2.0, y_c = 2.0, z_c = 2.0; 
    const double R = 1.0; 

    // Local domain indices
    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    #ifdef VOF_DIFFUSE
    {
        const double cn = params->Cn;
        const double denom = 2.0 * sqrt(2.0) * cn;
        const double dx = grid->dx_c[0];
        const double dy = grid->dy_c[0];
        const double dz = grid->dz_c[0];
        const int samples = 6;
        const double inv_samples = 1.0 / (double)samples;
        const double sample_weight = inv_samples * inv_samples * inv_samples;

        for (int k = Ks; k <= Ke; k++) {
            const double z_min = grid->zc[k] - 0.5 * dz;

            for (int j = Js; j <= Je; j++) {
                const double y_min = grid->yc[j] - 0.5 * dy;

                for (int i = Is; i <= Ie; i++) {
                    const double x_min = grid->xc[i] - 0.5 * dx;
                    double cl = 0.0;

                    /* Cell-average the continuous tanh profile so the FV field
                     * starts closer to the discrete equilibrium on coarse grids. */
                    for (int sk = 0; sk < samples; ++sk) {
                        const double z = z_min + (sk + 0.5) * dz * inv_samples;
                        for (int sj = 0; sj < samples; ++sj) {
                            const double y = y_min + (sj + 0.5) * dy * inv_samples;
                            for (int si = 0; si < samples; ++si) {
                                const double x = x_min + (si + 0.5) * dx * inv_samples;
                                const double r = sqrt((x - x_c) * (x - x_c) +
                                                      (y - y_c) * (y - y_c) +
                                                      (z - z_c) * (z - z_c));
                                const double s = R - r;
                                cl += sample_weight * 0.5 * (1.0 + tanh(s / (denom + 1e-30)));
                            }
                        }
                    }

                    if (cl < 0.0) {
                        cl = 0.0;
                    } else if (cl > 1.0) {
                        cl = 1.0;
                    }

                    F[k][j][i] = cl;
                    vof->C_L[k][j][i] = cl;
                }
            }
        }

        VOF_DIFFUSE_set_boundary_values(vof->C_L, data_bag);
        return;
    }
    #else
    // Grid parameters
    const double dx = grid->dx_c[0]; // Assuming uniform grid, dx=dy=dz
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    // Subgrid sampling parameters
    // INCREASED samples for better curvature initialization
    const int samples = 10; 
    const double inv_samples = 1.0 / samples;
    const double sample_weight = 1.0 / (samples * samples * samples);

    for (int k = Ks; k <= Ke; k++) {
        // Compute bottom-left-back corner of the cell
        double z_min = grid->zc[k] - 0.5 * dz;
        
        for (int j = Js; j <= Je; j++) {
            double y_min = grid->yc[j] - 0.5 * dy;
            
            for (int i = Is; i <= Ie; i++) {
                double x_min = grid->xc[i] - 0.5 * dx;
                
                double vol_count = 0.0;
                
                // Sub-grid loop
                for (int sk = 0; sk < samples; sk++) {
                    double z_sub = z_min + (sk + 0.5) * (dz * inv_samples);
                    
                    for (int sj = 0; sj < samples; sj++) {
                        double y_sub = y_min + (sj + 0.5) * (dy * inv_samples);
                        
                        for (int si = 0; si < samples; si++) {
                            double x_sub = x_min + (si + 0.5) * (dx * inv_samples);
                            
                            // Distance squared
                            double dist2 = (x_sub - x_c)*(x_sub - x_c) + 
                                           (y_sub - y_c)*(y_sub - y_c) + 
                                           (z_sub - z_c)*(z_sub - z_c);
                            
                            // Check if inside the sphere
                            if (dist2 < R * R) {
                                vol_count += 1.0;
                            }
                        }
                    }
                }
                
                // FIX: F = 1 inside the drop (Fluid 1), F = 0 outside (Fluid 2)
                F[k][j][i] = vol_count * sample_weight; 
            }
        }
    }
    #endif
}


/*
 * VoF_init_meniscus_154deg
 * ------------------------
 * Pre-shapes the gas–liquid interface around a fixed sphere at y_c=4.35
 * to the analytical equilibrium meniscus (θ=154°, Bo≪1).
 *
 * Geometry (all in non-dimensional units, D_p=1, R=0.5):
 *   Sphere center:      (xc, yc, zc) = (4.0, 4.35, 4.0)
 *   Far-field interface: y_inf = 4.0
 *   Contact line:        r_cl = 0.357, y_cl = 4.0
 *   Meniscus parameter:  A    = 0.127
 *   Blend radius:        r_blend = 0.8  (smooth back to flat beyond here)
 */
void VoF_init_meniscus_154deg(Cart3d_bag *data_bag)
{
    MAC_grid       *grid = data_bag->grid;
    VolumeFraction *vof  = data_bag->vof;
    double       ***F    = vof->F;

    /* ── sphere & interface parameters ─────────────────────────────────── */
    const double xc      = 4.0,  yc     = 4.35, zc    = 4.0;
    const double R       = 0.5;
    const double y_inf   = 4.0;   /* far-field interface height             */
    const double r_cl    = 0.357; /* R sin(alpha),  alpha = 45.57 deg       */
    const double y_cl    = 4.0;   /* y_c - R cos(alpha) = contact line y    */
    const double A       = 0.127; /* r_cl * tan(psi),  psi = 19.57 deg      */
    const double r_blend = 0.80;  /* blend meniscus → flat beyond here      */
    const double r_max   = 1.30;  /* domain corner distance (approx)        */

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    int Is = grid->G_Is, Ie = grid->G_Ie;
    int Js = grid->G_Js, Je = grid->G_Je;
    int Ks = grid->G_Ks, Ke = grid->G_Ke;

    const int samples   = 5;
    const double w_samp = 1.0 / ((double)(samples * samples * samples));

    for (int k = Ks; k <= Ke; ++k) {
    for (int j = Js; j <= Je; ++j) {
    for (int i = Is; i <= Ie; ++i) {

        double x0 = grid->xc[i];
        double y0 = grid->yc[j];
        double z0 = grid->zc[k];

        /* ── sub-sample the cell for anti-aliasing ────────────────────── */
        double frac = 0.0;

        for (int sk = 0; sk < samples; ++sk) {
        for (int sj = 0; sj < samples; ++sj) {
        for (int si = 0; si < samples; ++si) {

            double x = (x0 - 0.5*dx) + (si + 0.5) * (dx / samples);
            double y = (y0 - 0.5*dy) + (sj + 0.5) * (dy / samples);
            double z = (z0 - 0.5*dz) + (sk + 0.5) * (dz / samples);

            /* horizontal radius from sphere vertical axis */
            double dx_ = x - xc,  dz_ = z - zc;
            double r   = sqrt(dx_*dx_ + dz_*dz_);
            if (r < 1e-14) r = 1e-14;

            /* ── meniscus height at this (x,z) ───────────────────────── */
            double y_men;

            if (r <= r_cl) {
                /*
                 * Inside the contact-line footprint: the sphere occupies
                 * this column. Use the sphere surface itself as the "bottom"
                 * of the gas region. Below the sphere equator (y < yc) the
                 * submerged cap is water; above is air.
                 * Approximate: treat the sphere bottom as the interface
                 * locally (the IBM will correct this in 1 step).
                 */
                double dy_sph  = sqrt(R*R - r*r);  /* half-height of sphere */
                double y_bot   = yc - dy_sph;       /* bottom of sphere      */
                double y_top   = yc + dy_sph;       /* top of sphere         */

                /* The contact angle condition:  the interface exits the
                 * sphere at y_cl = 4.0. For the sub-cap region (r < r_cl),
                 * the cell is entirely below the sphere contact point, so
                 * it is in the water phase (F = 1) unless it is above the
                 * sphere top (impossible here since y_cl < y_top).
                 * A safe and smooth choice: set the local meniscus to y_cl
                 * everywhere inside the footprint.                          */
                y_men = y_cl;  /* = 4.0 */
            }
            else {
                /* ── log-profile meniscus ─────────────────────────────── */
                double y_log = y_cl + A * log(r_cl / r);

                if (r <= r_blend) {
                    /* pure log profile */
                    y_men = y_log;
                } else {
                    /* smooth blend to flat y_inf beyond r_blend            */
                    double t  = (r - r_blend) / (r_max - r_blend);
                    t = (t > 1.0) ? 1.0 : t;
                    /* smoothstep: 3t²-2t³                                  */
                    double s  = t * t * (3.0 - 2.0 * t);
                    y_men = (1.0 - s) * y_log + s * y_inf;
                }
            }

            /* ── phase assignment: F=1 (water) if sample_y < y_men ───── */
            frac += (y < y_men) ? 1.0 : 0.0;

        }}}  /* end sub-samples */

        F[k][j][i] = frac * w_samp;   /* 0 ≤ F ≤ 1 */

    }}}  /* end cells */
}


/*******************************************************************************
 * VoF_init_planar_droplet_on_static_cylinder_theta
 *
 * Initial condition for the Liu and Ding (2015) section 4.1 half-cylinder
 * wetting test in TWOD_CARTESIAN.  Domain x in [0, 1], y in [0, 2.5];
 * cylinder R_s = 0.5 sitting on the left symmetry plane at (0, 0.8).
 ******************************************************************************/
void VoF_init_planar_droplet_on_static_cylinder_theta(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double Rs = 0.5;
    const double xs = 0.0;
    const double ys = 0.8;
    const double R0 = 0.5;


    double R_f = R0;

    const double l_c = 0.75;
    const double xd  = xs;
    const double yd  = ys + l_c;

    /* Diffuse interface bandwidth follows the Cn used by VOF_DIFFUSE.  If the
     * caller has not yet set Cn (init runs before Input.c finalizes Cn when
     * Cn was given as -1 in parties.inp), fall back to 0.75 * h_ref. */
    const double dx_c = grid->dx_c[0];
    const double dy_c = grid->dy_c[0];
    double Cn = params->Cn;
    if (Cn <= 0.0) {
        double href = (dx_c < dy_c) ? dx_c : dy_c;
        Cn = 0.75 * href;
    }
    const double bw = 2.0 * sqrt(2.0) * Cn;

    const int    S  = 5;
    const double inv_S = 1.0 / (double)S;
    const double w_samp = inv_S * inv_S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy_c;
            for (int i = Is; i < Ie; ++i) {
                const double x_min = grid->xc[i] - 0.5 * dx_c;
                double frac = 0.0;

                /*
                 * Diffuse C_S analytic profile, identical to the one rebuilt
                 * by VOF_DIFFUSE_compute_C_S so the IC limiter clips C_L
                 * inside the cylinder before t = 0.  Without this the H5
                 * snapshot at t = 0 shows liquid bleeding through the
                 * diffuse solid; the next time step would clip it via the
                 * MCL limiter, but the IC frame would still mislead users.
                 */
                const double shift_S = sqrt(2.0) * log(19.0) * Cn;
                const double denom_S = 2.0 * sqrt(2.0) * Cn;

                for (int sj = 0; sj < S; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy_c * inv_S;
                    for (int si = 0; si < S; ++si) {
                        const double x = x_min + (si + 0.5) * dx_c * inv_S;

                        /* Signed distance to the L-G circle, positive inside
                         * the droplet.  This drives the diffuse C_L profile. */
                        double rd = sqrt((x - xd) * (x - xd) +
                                         (y - yd) * (y - yd));
                        double s_lg = R_f - rd;
                        double C_L_pre = 0.5 * (1.0 + tanh(s_lg / bw));

                        /* Diffuse C_S around the cylinder. */
                        double rs = sqrt((x - xs) * (x - xs) +
                                         (y - ys) * (y - ys));
                        double arg_S = (rs - (Rs - shift_S)) / denom_S;
                        double C_S_pre = 0.5 - 0.5 * tanh(arg_S);

                        /* Liu15 mass limiter applied at the sub-sample level
                         * so the diffuse C_L matches the contour the MCL
                         * post-step hook would produce (no liquid leaks into
                         * the diffuse solid). */
                        double cap = 1.0 - C_S_pre;
                        if (cap < 0.0) cap = 0.0;
                        if (C_L_pre > cap) C_L_pre = cap;
                        if (C_L_pre < 0.0) C_L_pre = 0.0;
                        frac += w_samp * C_L_pre;
                    }
                }

                F[k][j][i] = frac;
            }
        }
    }

#ifdef VOF_DIFFUSE
    for (int k = Ks; k < Ke; ++k)
        for (int j = Js; j < Je; ++j)
            for (int i = Is; i < Ie; ++i) {
                vof->C_L[k][j][i] = F[k][j][i];
            }
#endif

    vof->initial_volume = -1.0;
}

/*******************************************************************************
 * Lens volume between two intersecting spheres of radii R and r separated by d.
 * Returns 0 if separated and the full smaller sphere volume if one sphere
 * contains the other.
 ******************************************************************************/
static double mcl_sphere_lens_volume(double R, double r, double d)
{
    if (R <= 0.0 || r <= 0.0)
        return 0.0;

    const double V_R = (4.0 / 3.0) * M_PI * R * R * R;
    const double V_r = (4.0 / 3.0) * M_PI * r * r * r;

    if (d >= R + r)
        return 0.0;

    if (d <= fabs(R - r))
        return (V_R < V_r) ? V_R : V_r;

    if (d < 1.0e-14)
        return (V_R < V_r) ? V_R : V_r;

    const double h = R + r - d;

    double V = M_PI * h * h *
               (d * d + 2.0 * d * (R + r) - 3.0 * (R - r) * (R - r)) /
               (12.0 * d);

    if (V < 0.0) V = 0.0;

    return V;
}

static double axisym_bridge_seed_volume(double Rs, double Rb, double center_gap)
{
    const double V_seed = (4.0 / 3.0) * M_PI * Rb * Rb * Rb;
    return V_seed - 2.0 * mcl_sphere_lens_volume(Rs, Rb, center_gap);
}

static double solve_axisym_bridge_seed_radius(double Rs,
                                              double center_gap,
                                              double V_target)
{
    double lo = 1.0e-8;
    double hi = Rs;

    while (axisym_bridge_seed_volume(Rs, hi, center_gap) < V_target &&
           hi < 10.0 * Rs) {
        hi *= 1.5;
    }

    for (int it = 0; it < 120; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double V_mid = axisym_bridge_seed_volume(Rs, mid, center_gap);

        if (V_mid < V_target)
            lo = mid;
        else
            hi = mid;
    }

    return 0.5 * (lo + hi);
}

static int collect_axisym_bridge_particles(Cart3d_bag *data_bag,
                                           double X_out[2][3],
                                           double R_out[2])
{
    Lagrangian *lag = data_bag->lag;

    if (lag == NULL || lag->p_fixed_list == NULL || lag->p_mobile_list == NULL)
        return 0;

    const int n_total = lag->p_fixed_list->Np + lag->p_mobile_list->Np;
    if (n_total < 2)
        return 0;

    double *send_geom = (double *)calloc((size_t)n_total * 4, sizeof(double));
    double *recv_geom = (double *)calloc((size_t)n_total * 4, sizeof(double));
    int    *send_seen = (int *)calloc((size_t)n_total, sizeof(int));
    int    *recv_seen = (int *)calloc((size_t)n_total, sizeof(int));

    if (send_geom == NULL || recv_geom == NULL ||
        send_seen == NULL || recv_seen == NULL) {
        free(send_geom);
        free(recv_geom);
        free(send_seen);
        free(recv_seen);
        return 0;
    }

    Particle_list *lists[2] = {
        lag->p_fixed_list,
        lag->p_mobile_list
    };

    for (int list_id = 0; list_id < 2; ++list_id) {
        Particle *p = lists[list_id]->start;

        while (p != NULL) {
            if (p->ID >= 0 && p->ID < n_total) {
                send_geom[4 * p->ID + 0] = p->X[0];
                send_geom[4 * p->ID + 1] = p->X[1];
                send_geom[4 * p->ID + 2] = p->X[2];
                send_geom[4 * p->ID + 3] = p->R;
                send_seen[p->ID] = 1;
            }

            p = p->next;
        }
    }

    MPI_Allreduce(send_geom, recv_geom, 4 * n_total, MPI_DOUBLE, MPI_SUM, PCW);
    MPI_Allreduce(send_seen, recv_seen, n_total, MPI_INT, MPI_SUM, PCW);

    int found[2] = {-1, -1};

    for (int id = 0; id < n_total; ++id) {
        if (recv_seen[id] <= 0)
            continue;

        if (found[0] < 0 ||
            recv_geom[4 * id + 1] < recv_geom[4 * found[0] + 1]) {
            found[1] = found[0];
            found[0] = id;
        } else if (found[1] < 0 ||
                   recv_geom[4 * id + 1] < recv_geom[4 * found[1] + 1]) {
            found[1] = id;
        }
    }

    int ok = (found[0] >= 0 && found[1] >= 0);

    if (ok) {
        for (int n = 0; n < 2; ++n) {
            const int id = found[n];
            const double inv_seen = 1.0 / (double)recv_seen[id];
            X_out[n][0] = recv_geom[4 * id + 0] * inv_seen;
            X_out[n][1] = recv_geom[4 * id + 1] * inv_seen;
            X_out[n][2] = recv_geom[4 * id + 2] * inv_seen;
            R_out[n]    = recv_geom[4 * id + 3] * inv_seen;
        }
    }

    free(send_geom);
    free(recv_geom);
    free(send_seen);
    free(recv_seen);

    return ok;
}


/*******************************************************************************
 * VoF_init_axisymmetric_capillary_bridge_two_spheres
 *
 * Rajesh et al. (2026), section 3.1: quasi-static axial force in a liquid
 * bridge between two equal spherical beads.
 *
 * Coordinates:
 *   x -> radial coordinate r
 *   y -> axial coordinate z
 *
 * The bead centers/radii are taken from p_fixed.inp + p_mobile.inp.  The liquid
 * seed is the portion of a small sphere centered halfway between the beads that
 * lies outside both solid spheres.  Its radius is solved from the input
 * bridge_volume_over_r3.  The default V/R^3 = 0.09375 corresponds to
 * V = 0.75 uL and R = 2 mm in the experiments.  The diffuse C_L profile is
 * then clipped by the same diffuse solid profile used by
 * VOF_DIFFUSE_compute_C_S.
 ******************************************************************************/
void VoF_init_axisymmetric_capillary_bridge_two_spheres(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    double Xp[2][3];
    double Rp[2];

    if (!collect_axisym_bridge_particles(data_bag, Xp, Rp)) {
        Xp[0][0] = 0.0; Xp[0][1] = 1.5;  Xp[0][2] = 0.0; Rp[0] = 1.0;
        Xp[1][0] = 0.0; Xp[1][1] = 3.55; Xp[1][2] = 0.0; Rp[1] = 1.0;
    }

    const double R_lower = Rp[0];
    const double R_upper = Rp[1];
    const double R_ref = 0.5 * (R_lower + R_upper);

    const double z_lower = Xp[0][1];
    const double z_upper = Xp[1][1];
    const double z_mid = 0.5 * (z_lower + z_upper);
    const double center_gap = 0.5 * (z_upper - z_lower);
    const double surface_gap = z_upper - z_lower - R_lower - R_upper;

    double V_star = params->bridge_volume_over_r3;
    if (V_star <= 0.0)
        V_star = 0.09375;
    const double V_target = V_star * R_ref * R_ref * R_ref;
    const double R_seed = solve_axisym_bridge_seed_radius(R_ref,
                                                          center_gap,
                                                          V_target);

    const double dr = grid->dx_c[0];
    const double dz = grid->dy_c[0];

    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double href = (dr < dz) ? dr : dz;
        Cn = 0.75 * href;
    }

    const double bw = 2.0 * sqrt(2.0) * Cn;
    const double shift_S = sqrt(2.0) * log(19.0) * Cn;
    const double denom_S = 2.0 * sqrt(2.0) * Cn;

    {
        char msg[320];
        snprintf(msg, sizeof(msg),
                 "Rajesh-Sauret capillary bridge IC: R=%.6f gap/R=%.6f "
                 "V/R^3=%.6f seed_radius/R=%.6f centers z=(%.6f, %.6f) "
                 "theta=%.2fdeg Cn=%.6g\n",
                 R_ref, surface_gap / R_ref, V_star, R_seed / R_ref,
                 z_lower, z_upper, params->contact_angle_deg, Cn);
        Display_progress(params, msg);
    }

    const int S = 6;
    const double inv_S = 1.0 / (double)S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double z_min = grid->yc[j] - 0.5 * dz;

            for (int i = Is; i < Ie; ++i) {
                const double r_min = grid->xc[i] - 0.5 * dr;

                double sum_CL = 0.0;
                double sum_w  = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double z = z_min + (sj + 0.5) * dz * inv_S;

                    for (int si = 0; si < S; ++si) {
                        const double r = r_min + (si + 0.5) * dr * inv_S;
                        const double w_axi = fabs(r);

                        const double dist_seed =
                            sqrt(r * r + (z - z_mid) * (z - z_mid));
                        const double s_lg = R_seed - dist_seed;
                        double C_L_pre =
                            0.5 * (1.0 + tanh(s_lg / (bw + 1.0e-30)));

                        const double dist_lower =
                            sqrt((r - Xp[0][0]) * (r - Xp[0][0]) +
                                 (z - Xp[0][1]) * (z - Xp[0][1]));
                        const double dist_upper =
                            sqrt((r - Xp[1][0]) * (r - Xp[1][0]) +
                                 (z - Xp[1][1]) * (z - Xp[1][1]));

                        const double arg_lower =
                            (dist_lower - (R_lower - shift_S)) /
                            (denom_S + 1.0e-30);
                        const double arg_upper =
                            (dist_upper - (R_upper - shift_S)) /
                            (denom_S + 1.0e-30);

                        const double C_S_lower = 0.5 - 0.5 * tanh(arg_lower);
                        const double C_S_upper = 0.5 - 0.5 * tanh(arg_upper);
                        double C_S_pre = (C_S_lower > C_S_upper) ?
                                         C_S_lower : C_S_upper;

                        double cap = 1.0 - C_S_pre;
                        if (cap < 0.0) cap = 0.0;
                        if (cap > 1.0) cap = 1.0;

                        if (C_L_pre > cap) C_L_pre = cap;
                        if (C_L_pre < 0.0) C_L_pre = 0.0;
                        if (C_L_pre > 1.0) C_L_pre = 1.0;

                        sum_CL += w_axi * C_L_pre;
                        sum_w  += w_axi;
                    }
                }

                if (sum_w > 1.0e-30)
                    F[k][j][i] = sum_CL / sum_w;
                else
                    F[k][j][i] = 0.0;
            }
        }
    }

#ifdef VOF_DIFFUSE
    for (int k = Ks; k < Ke; ++k)
        for (int j = Js; j < Je; ++j)
            for (int i = Is; i < Ie; ++i)
                vof->C_L[k][j][i] = F[k][j][i];
#endif

    vof->initial_volume = -1.0;
}


/*******************************************************************************
 * VoF_init_liu17_two_sinking_cylinders_2d
 *
 * Full-domain planar 2D variant of the Liu17 section 6.3 sinking-cylinder
 * setup.  The two mobile cylinders are supplied by p_mobile.inp; this
 * initializer only builds the initially flat water pool:
 *
 *   domain        0 <= x <= 12.8, 0 <= y <= 9.6
 *   water surface y = 8.6
 *   liquid        y < 8.6
 ******************************************************************************/
void VoF_init_liu17_two_sinking_cylinders_2d(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double y_interface = 8.6;

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];

#ifdef VOF_DIFFUSE
    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double href = (dx < dy) ? dx : dy;
        Cn = 0.75 * href;
    }
    const double bw = 2.0 * sqrt(2.0) * Cn;
#endif

    const int S = 8;
    const double inv_S = 1.0 / (double)S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i < Ie; ++i) {
                double cl = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy * inv_S;
#ifdef VOF_DIFFUSE
                    const double s = y_interface - y;
                    cl += inv_S * 0.5 * (1.0 + tanh(s / (bw + 1.0e-30)));
#else
                    if (y <= y_interface) {
                        cl += inv_S;
                    }
#endif
                }

                if (cl < 0.0) {
                    cl = 0.0;
                } else if (cl > 1.0) {
                    cl = 1.0;
                }

                F[k][j][i] = cl;
#ifdef VOF_DIFFUSE
                vof->C_L[k][j][i] = cl;
#endif
            }
        }
    }

    vof->initial_volume = -1.0;

    {
        char msg[180];
        snprintf(msg, sizeof(msg),
                 "Liu17 two-cylinder IC: planar water surface y=%.6f, "
                 "Cylinders from p_mobile.inp, Cn=%.6g\n",
                 y_interface, params->Cn);
        Display_progress(params, msg);
    }
}


static void VoF_init_liu17_full_domain_sinking_cylinders_2d(
    Cart3d_bag *data_bag,
    const char *case_label)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    const double y_interface = 8.6;

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];

#ifdef VOF_DIFFUSE
    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double href = (dx < dy) ? dx : dy;
        Cn = 0.75 * href;
    }
    const double bw = 2.0 * sqrt(2.0) * Cn;
#endif

    const int S = 8;
    const double inv_S = 1.0 / (double)S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double y_min = grid->yc[j] - 0.5 * dy;

            for (int i = Is; i < Ie; ++i) {
                double cl = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double y = y_min + (sj + 0.5) * dy * inv_S;
#ifdef VOF_DIFFUSE
                    const double s = y_interface - y;
                    cl += inv_S * 0.5 * (1.0 + tanh(s / (bw + 1.0e-30)));
#else
                    if (y <= y_interface) {
                        cl += inv_S;
                    }
#endif
                }

                if (cl < 0.0) {
                    cl = 0.0;
                } else if (cl > 1.0) {
                    cl = 1.0;
                }

                F[k][j][i] = cl;
#ifdef VOF_DIFFUSE
                vof->C_L[k][j][i] = cl;
#endif
            }
        }
    }

    vof->initial_volume = -1.0;

    {
        char msg[180];
        snprintf(msg, sizeof(msg),
                 "%s: planar water surface y=%.6f, "
                 "Cylinders from p_mobile.inp, Cn=%.6g\n",
                 case_label, y_interface, params->Cn);
        Display_progress(params, msg);
    }
}

/*******************************************************************************
 * VoF_init_liu17_three_sinking_cylinders_2d
 *
 * Full-domain planar 2D triplet of sinking cylinders using the same initial
 * water level and particle height as the two-cylinder case. The first mobile
 * particle in p_mobile.inp is the middle cylinder; the second and third are
 * placed symmetrically to the left and right.
 ******************************************************************************/
void VoF_init_liu17_three_sinking_cylinders_2d(Cart3d_bag *data_bag)
{
    VoF_init_liu17_full_domain_sinking_cylinders_2d(
        data_bag, "Liu17 three-cylinder IC");
}





typedef struct {
    double beta;
    double theta;
    double R;
    double a;
    double z_mid;
    double z_half;
    double r_contact;
    double circle_radius;
    double circle_center_r;
    int radial_branch;
    double target_volume;
} Nguyen21_bridge_geom;

static int nguyen21_bridge_geom_from_beta(double R,
                                          double a,
                                          double theta,
                                          double beta,
                                          Nguyen21_bridge_geom *geom)
{
    double delta = theta - beta;
    while (delta <= -0.5 * M_PI)
        delta += M_PI;
    while (delta > 0.5 * M_PI)
        delta -= M_PI;

    const double cos_tb = cos(delta);
    if (R <= 0.0 || a <= R || cos_tb <= 1.0e-12)
        return 0;

    const double z_half = a - R * cos(beta);
    const double r_contact = R * sin(beta);
    const double circle_radius = z_half / cos_tb;
    const double circle_center_r =
        r_contact - circle_radius * sin(delta);
    const int radial_branch =
        (r_contact >= circle_center_r) ? 1 : -1;
    const double r_mid =
        circle_center_r + (double)radial_branch * circle_radius;

    if (z_half <= 0.0 || circle_radius <= 0.0 ||
        z_half > circle_radius * (1.0 + 1.0e-12) ||
        r_mid <= 0.0) {
        return 0;
    }

    geom->beta = beta;
    geom->theta = theta;
    geom->R = R;
    geom->a = a;
    geom->z_half = z_half;
    geom->r_contact = r_contact;
    geom->circle_radius = circle_radius;
    geom->circle_center_r = circle_center_r;
    geom->radial_branch = radial_branch;

    return 1;
}

static double nguyen21_bridge_volume_for_geom(const Nguyen21_bridge_geom *geom)
{
    const int N = 2048;
    const double R = geom->R;
    const double a = geom->a;
    const double dz = 2.0 * geom->z_half / (double)N;
    double integral = 0.0;

    for (int n = 0; n <= N; ++n) {
        const double z = -geom->z_half + (double)n * dz;
        const double circle_arg =
            geom->circle_radius * geom->circle_radius - z * z;
        double r_outer = -1.0;
        double r_inner = 0.0;
        double area = 0.0;

        if (circle_arg >= 0.0) {
            r_outer = geom->circle_center_r +
                      (double)geom->radial_branch * sqrt(circle_arg);
        }

        if (fabs(z + a) < R) {
            const double r_s = sqrt(R * R - (z + a) * (z + a));
            if (r_s > r_inner)
                r_inner = r_s;
        }
        if (fabs(z - a) < R) {
            const double r_s = sqrt(R * R - (z - a) * (z - a));
            if (r_s > r_inner)
                r_inner = r_s;
        }

        if (r_outer > r_inner)
            area = M_PI * (r_outer * r_outer - r_inner * r_inner);

        if (n == 0 || n == N)
            integral += area;
        else if (n % 2 == 1)
            integral += 4.0 * area;
        else
            integral += 2.0 * area;
    }

    return integral * dz / 3.0;
}

static double nguyen21_bridge_volume_for_beta(double R,
                                              double a,
                                              double theta,
                                              double beta)
{
    Nguyen21_bridge_geom geom;

    if (!nguyen21_bridge_geom_from_beta(R, a, theta, beta, &geom))
        return -1.0;

    return nguyen21_bridge_volume_for_geom(&geom);
}

static int solve_nguyen21_bridge_geom(double R,
                                      double a,
                                      double theta,
                                      double target_volume,
                                      Nguyen21_bridge_geom *geom)
{
    const double beta_min = 1.0e-6;
    const double beta_max = M_PI - 1.0e-6;
    const int scan_count = 4096;
    double lo = 0.0;
    double hi = 0.0;
    double prev_beta = 0.0;
    double prev_vol = 0.0;
    int prev_branch = 0;
    int have_prev = 0;
    int bracketed = 0;

    for (int n = 0; n <= scan_count; ++n) {
        const double beta =
            beta_min + (beta_max - beta_min) * (double)n / (double)scan_count;
        Nguyen21_bridge_geom trial_geom;
        double vol = -1.0;

        if (!nguyen21_bridge_geom_from_beta(R, a, theta, beta, &trial_geom))
            continue;

        vol = nguyen21_bridge_volume_for_geom(&trial_geom);
        if (vol < 0.0)
            continue;

        if (have_prev && prev_branch == trial_geom.radial_branch &&
            (prev_vol - target_volume) * (vol - target_volume) <= 0.0) {
            lo = prev_beta;
            hi = beta;
            bracketed = 1;
            break;
        }

        prev_beta = beta;
        prev_vol = vol;
        prev_branch = trial_geom.radial_branch;
        have_prev = 1;
    }

    if (!bracketed)
        return 0;

    for (int it = 0; it < 80; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double vol_lo =
            nguyen21_bridge_volume_for_beta(R, a, theta, lo);
        const double vol =
            nguyen21_bridge_volume_for_beta(R, a, theta, mid);

        if (vol < 0.0 || vol_lo < 0.0)
            return 0;

        if ((vol_lo - target_volume) * (vol - target_volume) <= 0.0)
            hi = mid;
        else
            lo = mid;
    }

    if (!nguyen21_bridge_geom_from_beta(R, a, theta, 0.5 * (lo + hi), geom))
        return 0;

    geom->target_volume = target_volume;
    return 1;
}

/*******************************************************************************
 * VoF_init_nguyen21_axisymmetric_static_bridge
 *
 * Nguyen et al. (2021), APT section 5.3: static liquid bridge force between
 * two stationary equal spheres.  The axisymmetric domain uses x as radius r and
 * y as the axial coordinate.  The bridge seed is a toroidal circular-arc
 * approximation: the contact line is placed so that the sharp bridge volume is
 * params->bridge_volume_over_r3 * R^3 and the lower/upper contact points meet
 * the particle surfaces at params->contact_angle_deg through the liquid.
 ******************************************************************************/
void VoF_init_nguyen21_axisymmetric_static_bridge(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    double Xp[2][3];
    double Rp[2];

    if (!collect_axisym_bridge_particles(data_bag, Xp, Rp)) {
        Xp[0][0] = 0.0; Xp[0][1] = 0.725; Xp[0][2] = 0.0; Rp[0] = 0.5;
        Xp[1][0] = 0.0; Xp[1][1] = 1.975; Xp[1][2] = 0.0; Rp[1] = 0.5;
    }

    const double R_ref = 0.5 * (Rp[0] + Rp[1]);
    const double z_lower = Xp[0][1];
    const double z_upper = Xp[1][1];
    const double z_mid = 0.5 * (z_lower + z_upper);
    const double a = 0.5 * (z_upper - z_lower);
    const double surface_gap = z_upper - z_lower - Rp[0] - Rp[1];
    const double theta = params->contact_angle_deg * M_PI / 180.0;

    double V_star = params->bridge_volume_over_r3;
    if (V_star <= 0.0)
        V_star = 0.3 * (4.0 / 3.0) * M_PI;

    const double V_target = V_star * R_ref * R_ref * R_ref;

    Nguyen21_bridge_geom geom;
    if (!solve_nguyen21_bridge_geom(R_ref, a, theta, V_target, &geom)) {
        Display_throw_warning("Nguyen21 bridge IC could not bracket the toroidal seed; using spherical bridge seed.\n",
                              params);
        VoF_init_axisymmetric_capillary_bridge_two_spheres(data_bag);
        return;
    }
    geom.z_mid = z_mid;

    const double dr = grid->dx_c[0];
    const double dz = grid->dy_c[0];

    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double href = (dr < dz) ? dr : dz;
        Cn = 0.75 * href;
    }

    const double bw = 2.0 * sqrt(2.0) * Cn;
    const double shift_S = sqrt(2.0) * log(19.0) * Cn;
    const double denom_S = 2.0 * sqrt(2.0) * Cn;

    {
        char msg[420];
        snprintf(msg, sizeof(msg),
                 "Nguyen21 static bridge IC: R=%.6f s/dp=%.6f V/R^3=%.9f "
                 "theta=%.2fdeg beta=%.4fdeg rc/R=%.6f rneck/R=%.6f "
                 "branch=%+d Cn=%.6g\n",
                 R_ref, surface_gap / (2.0 * R_ref), V_star,
                 params->contact_angle_deg, geom.beta * 180.0 / M_PI,
                 geom.r_contact / R_ref,
                 (geom.circle_center_r +
                  (double)geom.radial_branch * geom.circle_radius) / R_ref,
                 geom.radial_branch, Cn);
        Display_progress(params, msg);
    }

    const int S = 6;
    const double inv_S = 1.0 / (double)S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double z_min = grid->yc[j] - 0.5 * dz;

            for (int i = Is; i < Ie; ++i) {
                const double r_min = grid->xc[i] - 0.5 * dr;
                double sum_CL = 0.0;
                double sum_w  = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double z =
                        z_min + (sj + 0.5) * dz * inv_S;
                    const double zloc = z - z_mid;

                    for (int si = 0; si < S; ++si) {
                        const double r =
                            r_min + (si + 0.5) * dr * inv_S;
                        const double w_axi = fabs(r);

                        const double circle_arg =
                            geom.circle_radius * geom.circle_radius -
                            zloc * zloc;
                        double s_circle = -1.0e30;
                        if (circle_arg >= 0.0) {
                            const double r_surface =
                                geom.circle_center_r +
                                (double)geom.radial_branch *
                                sqrt(circle_arg);
                            const double dist_circle =
                                sqrt((r - geom.circle_center_r) *
                                     (r - geom.circle_center_r) +
                                     zloc * zloc);
                            double s_signed =
                                (geom.radial_branch > 0) ?
                                geom.circle_radius - dist_circle :
                                dist_circle - geom.circle_radius;

                            if (r <= r_surface)
                                s_circle = fabs(s_signed);
                            else
                                s_circle = -fabs(s_signed);
                        }
                        const double s_axial =
                            geom.z_half - fabs(zloc);
                        const double s_lg =
                            (s_circle < s_axial) ? s_circle : s_axial;

                        double C_L_pre =
                            0.5 * (1.0 + tanh(s_lg / (bw + 1.0e-30)));

                        const double dist_lower =
                            sqrt((r - Xp[0][0]) * (r - Xp[0][0]) +
                                 (z - Xp[0][1]) * (z - Xp[0][1]));
                        const double dist_upper =
                            sqrt((r - Xp[1][0]) * (r - Xp[1][0]) +
                                 (z - Xp[1][1]) * (z - Xp[1][1]));

                        const double arg_lower =
                            (dist_lower - (Rp[0] - shift_S)) /
                            (denom_S + 1.0e-30);
                        const double arg_upper =
                            (dist_upper - (Rp[1] - shift_S)) /
                            (denom_S + 1.0e-30);

                        const double C_S_lower = 0.5 - 0.5 * tanh(arg_lower);
                        const double C_S_upper = 0.5 - 0.5 * tanh(arg_upper);
                        double C_S_pre = (C_S_lower > C_S_upper) ?
                                         C_S_lower : C_S_upper;

                        double cap = 1.0 - C_S_pre;
                        if (cap < 0.0) cap = 0.0;
                        if (cap > 1.0) cap = 1.0;

                        if (C_L_pre > cap) C_L_pre = cap;
                        if (C_L_pre < 0.0) C_L_pre = 0.0;
                        if (C_L_pre > 1.0) C_L_pre = 1.0;

                        sum_CL += w_axi * C_L_pre;
                        sum_w  += w_axi;
                    }
                }

                F[k][j][i] = (sum_w > 1.0e-30) ? sum_CL / sum_w : 0.0;
            }
        }
    }

#ifdef VOF_DIFFUSE
    for (int k = Ks; k < Ke; ++k)
        for (int j = Js; j < Je; ++j)
            for (int i = Is; i < Ie; ++i)
                vof->C_L[k][j][i] = F[k][j][i];
#endif

    vof->initial_volume = -1.0;
}


/*******************************************************************************
 * VoF_init_axisymmetric_droplet_on_static_sphere_theta
 *
 * Axisymmetric static droplet on a fixed sphere.
 *
 * The initial liquid-gas interface is an intersecting sphere chosen so that the
 * contact angle through the liquid is theta_IC = 120 deg.  The liquid field is
 * initialized directly as a diffuse tanh profile, then clipped by the diffuse
 * solid indicator through C_L <= 1 - C_S.
 *
 * Coordinates:
 *   x -> r
 *   y -> z
 *
 * Default geometry:
 *   solid sphere radius Rs = 1
 *   solid sphere center    (r_s, z_s) = (0, 1.5)
 *   target free-droplet diameter D = 1
 *
 * The liquid-gas curvature radius R_f is solved so that the sharp intersecting
 * cap has the same volume as a free sphere of radius D/2.
 ******************************************************************************/
void VoF_init_axisymmetric_droplet_on_static_sphere_theta(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;
    double       ***F      = vof->F;

    /*
     * Keep this consistent with p_fixed.inp.
     * In AXISYM_RZ: x is radius r, y is axial coordinate z.
     */
    const double Rs = 1.0;
    const double r_s = 0.0;
    const double z_s = 1.5;

    /*
     * Target free-droplet volume.  D = 1 follows the Liu-type setup.
     */
    const double D_free = 1.0;
    const double R_free = 0.5 * D_free;
    const double V_target = (4.0 / 3.0) * M_PI * R_free * R_free * R_free;

    /*
     * Desired physical contact angle through the liquid.
     *
     * The two-sphere center-distance formula uses the supplementary geometric
     * angle, exactly like the 2D cylinder initializer:
     *
     *   l^2 = Rs^2 + R_f^2 + 2 Rs R_f cos(pi - theta_liquid)
     *
     * For theta_liquid = 120 deg, the droplet center is above the solid sphere
     * and the liquid-gas sphere intersects the solid sphere.
     */
    const double theta_liquid_deg = params->contact_angle_deg;
    const double theta_liquid     = theta_liquid_deg * M_PI / 180.0;
    const double theta_geom       = M_PI - theta_liquid;
    const double cosTG            = cos(theta_geom);

    /*
     * Solve for the liquid-gas curvature radius R_f such that the sharp
     * cap volume outside the solid sphere equals the target free-droplet volume.
     */
    double R_f = R_free;

    for (int it = 0; it < 200; ++it) {
        const double l = sqrt(Rs * Rs + R_f * R_f + 2.0 * Rs * R_f * cosTG);

        const double V_drop =
            (4.0 / 3.0) * M_PI * R_f * R_f * R_f
            - mcl_sphere_lens_volume(Rs, R_f, l);

        const double res = V_drop - V_target;

        if (fabs(res) < 1.0e-13)
            break;

        const double h  = 1.0e-6 * fmax(R_f, 0.1);
        const double R2 = R_f + h;
        const double l2 = sqrt(Rs * Rs + R2 * R2 + 2.0 * Rs * R2 * cosTG);

        const double V2 =
            (4.0 / 3.0) * M_PI * R2 * R2 * R2
            - mcl_sphere_lens_volume(Rs, R2, l2);

        const double dV = (V2 - V_drop) / h;

        if (fabs(dV) < 1.0e-14)
            break;

        double step = res / dV;

        if (step >  0.5 * R_f) step =  0.5 * R_f;
        if (step < -0.5 * R_f) step = -0.5 * R_f;

        R_f -= step;

        if (R_f < 0.05)
            R_f = 0.05;
    }

    const double l_c = sqrt(Rs * Rs + R_f * R_f + 2.0 * Rs * R_f * cosTG);

    const double r_d = r_s;
    const double z_d = z_s + l_c;

    /*
     * Diffuse interface bandwidth.  This is the same convention as the planar
     * initializer:
     *
     *   C_L = 0.5 * [1 + tanh(s_lg / (2 sqrt(2) Cn))]
     *
     * where s_lg > 0 inside the droplet.
     */
    const double dr = grid->dx_c[0];
    const double dz = grid->dy_c[0];

    double Cn = params->Cn;
    if (Cn <= 0.0) {
        const double href = (dr < dz) ? dr : dz;
        Cn = 0.75 * href;
    }

    const double bw = 2.0 * sqrt(2.0) * Cn;

    /*
     * Diffuse solid profile used only for IC clipping.  This should match the
     * solid C_S construction used later by VOF_DIFFUSE_compute_C_S.
     */
    const double shift_S = sqrt(2.0) * log(19.0) * Cn;
    const double denom_S = 2.0 * sqrt(2.0) * Cn;

    {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "Axisym sphere IC: theta=%.1fdeg R_f=%.8f l_c=%.8f "
                 "center=(%.6f, %.6f) V_target=%.8e Cn=%.4g\n",
                 theta_liquid_deg, R_f, l_c, r_d, z_d, V_target, Cn);
        Display_progress(params, msg);
    }

    const int    S = 5;
    const double inv_S = 1.0 / (double)S;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    for (int k = Ks; k < Ke; ++k) {
        for (int j = Js; j < Je; ++j) {
            const double z_min = grid->yc[j] - 0.5 * dz;

            for (int i = Is; i < Ie; ++i) {
                const double r_min = grid->xc[i] - 0.5 * dr;

                double sum_CL = 0.0;
                double sum_w  = 0.0;

                for (int sj = 0; sj < S; ++sj) {
                    const double z = z_min + (sj + 0.5) * dz * inv_S;

                    for (int si = 0; si < S; ++si) {
                        const double r = r_min + (si + 0.5) * dr * inv_S;

                        /*
                         * Axisymmetric volume weighting.  For normal interior
                         * cells away from the axis this is almost identical to
                         * arithmetic averaging, but it is more consistent near
                         * r = 0.
                         */
                        const double w_axi = fabs(r);

                        /*
                         * Signed distance to the liquid-gas sphere.
                         * Positive inside the droplet.
                         */
                        const double dist_lg =
                            sqrt((r - r_d) * (r - r_d) +
                                 (z - z_d) * (z - z_d));

                        const double s_lg = R_f - dist_lg;

                        double C_L_pre = 0.5 * (1.0 + tanh(s_lg / bw));

                        /*
                         * Diffuse solid sphere.  C_S = 1 inside solid,
                         * C_S = 0 outside solid.
                         */
                        const double dist_s =
                            sqrt((r - r_s) * (r - r_s) +
                                 (z - z_s) * (z - z_s));

                        const double arg_S = (dist_s - (Rs - shift_S)) / denom_S;
                        const double C_S_pre = 0.5 - 0.5 * tanh(arg_S);

                        /*
                         * Enforce no liquid inside the diffuse solid at t = 0.
                         * This avoids a spurious first-step mass correction.
                         */
                        double cap = 1.0 - C_S_pre;
                        if (cap < 0.0) cap = 0.0;
                        if (cap > 1.0) cap = 1.0;

                        if (C_L_pre > cap) C_L_pre = cap;
                        if (C_L_pre < 0.0) C_L_pre = 0.0;
                        if (C_L_pre > 1.0) C_L_pre = 1.0;

                        sum_CL += w_axi * C_L_pre;
                        sum_w  += w_axi;
                    }
                }

                if (sum_w > 1.0e-30)
                    F[k][j][i] = sum_CL / sum_w;
                else
                    F[k][j][i] = 0.0;
            }
        }
    }

#ifdef VOF_DIFFUSE
    for (int k = Ks; k < Ke; ++k)
        for (int j = Js; j < Je; ++j)
            for (int i = Is; i < Ie; ++i)
                vof->C_L[k][j][i] = F[k][j][i];
#endif

    vof->initial_volume = -1.0;
}

#endif // VOF_PLIC
