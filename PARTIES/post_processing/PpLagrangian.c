
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"

#include "Display.h"
#include "PpLagrangian.h"

/******************************************************************************/
/*
 * Perform post-processing operations
 */
/******************************************************************************/
double PpLagrangian_get_p_flux(Cart3d_bag *data_bag, Debug_trace *dtrace) {

	Particle *p;

	Parameters *params = data_bag -> params;
	MAC_grid   *grid   = data_bag -> grid;

	Particle_list *p_mobile_list = data_bag -> lag -> p_mobile_list;
	Particle_list *p_fixed_list  = data_bag -> lag -> p_fixed_list;

	// Check state of linked lists
	Display_assert_list_state(p_mobile_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));
	Display_assert_list_state(p_fixed_list, LIST_STATE_LOCAL, params, DTRACE("Display_assert_list_state"));

	double p_flux = 0.0;

	p = p_mobile_list -> start;
	while (p != NULL) {

		p_flux += p->M / params->rho_s * p->U[0];

		p = p -> next;
	}

	// Reduce (sum) p_flux to all processors
	MPI_Allreduce(MPI_IN_PLACE, &p_flux, 1, MPI_DOUBLE, MPI_SUM, PCW);

	// Flux per unit width
	p_flux = p_flux / (params->Lx * params->Lz);

	return p_flux;
}
