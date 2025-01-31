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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



#ifdef VOF_PLIC
#include "VolumeFraction.h"  // So we can see the VoF_create, etc.
#endif


/******************************************************************************/
/*
 This function initializes the volume fraction field for a bubble, based on 
 Basilisk test case. The Basilisk test case is a 2D bubble but built in 3D.
 */
/******************************************************************************/
void VoF_init_bubble(Cart3d_bag *data_bag)
{
    
    double x_c, y_c, z_c, R; // bubble center and radius
    
    MAC_grid    *grid   = data_bag->grid;
    Parameters  *params = data_bag->params;
    VolumeFraction *vof = data_bag->vof; // store convenience pointer

    // If the user didn't enable PLIC or set init_type != 1, do nothing
    if (!params->plic_enabled || params->init_type != 1)
        return;

    double ***F = vof->F;  // cell-centered volume fraction array

    // Retrieve bubble center and radius from params
    x_c = 1.0;
    y_c = 0.0;
    z_c = 0.0;
    R   = 0.25;

    // We assume 2D if (zmax - zmin) = 0 or NZ=1, but code can be 3D if needed
    // We get the cell-centered coordinates from grid->xc, grid->yc, grid->zc
    double *xc = grid->xc;
    double *yc = grid->yc;
    double *zc = grid->zc; // in 2D, you might not even loop over z if NZ=1

    // Indices of your local domain portion
    int Is = grid->G_Is; // or L_Is, depending on your code’s usage
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
                    // Inside the bubble => 100% fluid #1
                    F[k][j][i] = 1.0;
                } else {
                    // Outside the bubble => fluid #2
                    F[k][j][i] = 0.0;
                }
            }
        }
    }

    // Optionally handle boundary or ghost cells, or do synchronization if you do
    // multi-process. For example, you may want to call:
    // VoF_set_boundary_values(vof, grid, params);

    // That’s it for the bubble initialization
    printf("VOF bubble initialized (rank=%d)\n", params->rank);
}
