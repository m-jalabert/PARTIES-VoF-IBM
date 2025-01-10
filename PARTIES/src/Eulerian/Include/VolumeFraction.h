#ifndef VOLUMEFRACTION_H
 #define VOLUMEFRACTION_H

#include "Boundary.h"

#ifdef VOF_PLIC // Only compile if VOF_PLIC is defined in Boundary.h 

#include "definitions.h"
#include "DataTypes.h"

// "Constructor"-style function (similar to Conc_create)
VolumeFraction *VoF_create(MAC_grid *grid, Parameters *params);

// "Destructor"-style function (similar to Conc_destroy)
void VoF_destroy(VolumeFraction *vof, MAC_grid *grid, Parameters *params);

// Main routine to advance the volume fraction (akin to Conc_int_equations)
void VoF_int_equations(Cart3d_bag *data_bag, Debug_trace *dtrace);

// Additional prototypes: boundary updates, advection, etc.
void VoF_set_boundary_values(VolumeFraction *vof, MAC_grid *grid, Parameters *params);
void VoF_set_conv_viscous_vof(Cart3d_bag *data_bag);
// ... or however you want to name them ...

#endif // VOF_PLIC
#endif // VOLUMEFRACTION_H



