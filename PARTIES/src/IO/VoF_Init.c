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
 Rider & Kothe (1998) test case. 
 */
/******************************************************************************/
void VoF_init_bubble(Cart3d_bag *data_bag)
{
    
    double x_c, y_c, z_c, R; // bubble center and radius
    
    MAC_grid    *grid   = data_bag->grid;
    Parameters  *params = data_bag->params;
    VolumeFraction *vof = data_bag->vof; // 

    double ***F = vof->F;  // cell-centered volume fraction array


    x_c = 0.5;
    y_c = 0.6875;
    z_c = 0.5;
    R   = 0.1875;


    double *xc = grid->xc;
    double *yc = grid->yc;
    double *zc = grid->zc; 

    // Indices of your local domain portion
    int Is = grid->G_Is; 
    int Ie = grid->G_Ie;
    int Js = grid->G_Js;
    int Je = grid->G_Je;
    int Ks = grid->G_Ks;
    int Ke = grid->G_Ke;

    // Loop over local portion
    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {
                // Coordinates of the cell center (xc[i], yc[j], zc[k])
                double dx = xc[i] - x_c;
                double dy = yc[j] - y_c;
                double dz = zc[k] - z_c;

                double distSq = dx*dx + dy*dy + dz*dz;
                double rSq    = R*R;

                if (distSq < rSq) {
                    // Inside the bubble
                    F[k][j][i] = 0.0;
                } else {
                    // Outside the bubble 
                    F[k][j][i] = 1.0;
                }
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
