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
// Single-cell initialization (vof->F) for debugging
/******************************************************************************/
void VoF_init_single_cell(Cart3d_bag *data_bag)
{
    MAC_grid    *grid   = data_bag->grid;
    VolumeFraction *vof = data_bag->vof;

    double ***F = vof->F;  // cell-centered volume fraction array

    // Indices of your local domain portion
    int Is = grid->G_Is; 
    int Ie = grid->G_Ie;
    int Js = grid->G_Js;
    int Je = grid->G_Je;
    int Ks = grid->G_Ks;
    int Ke = grid->G_Ke;

    // We want cell (1,1,1) to have F=0.5, everything else =0
    // for a 3×3×3 domain, global indices go from 0..2
    // so the "center cell" is indeed i=1, j=1, k=1 if it lies in local range

    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                if (i == 3 && j == 3 && k == 3)  {
                    F[k][j][i] = 0.5;  // 
                } else {
                    F[k][j][i] = 0.0;
                }
            }
        }
    }
}
