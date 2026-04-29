/* Add function that perform analysis for output during runtime here */

#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Output.h"
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
#include "post_processing.h"
#include "TwodOps.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>




void Post_processing_conc(Cart3d_bag *data_bag){



	MAC_grid   *grid   = data_bag -> grid;
	Parameters *params = data_bag -> params;
	int int_stop, status, iconc;
	FILE *fid;
	FILE *timeFile;
	char message[500];
	Concentration **c = data_bag -> c;
	int NConc = params->NConc;

	double time, max_time, dt;
	double output_time_interval;
	int ntime;


	for (iconc=0; iconc< NConc; iconc++) {

			if (params->conc_output_deposit_height[iconc]) {

				Conc_integrate_deposited_height(c[iconc], grid, params, params->richardson[iconc]/(params->grav[1]));
			} // if

			Conc_store_old_data(c[iconc], grid, params);




	if ((params->conc_output_dump[iconc]) && (ntime % 20 == 0) && (time > 1.0) ){

		// Write the final deposit profile
		if (params->rank==0)
			printf(" Writing dumped conc to file at time:%f (warning the file index may not be accurate) \n", time);
		Conc_dump_particles(c, data_bag);
		Output_write_ascii_deposit_height_dumped(c, params, grid, params->noutput) ;

	} // if dump_conc
	} // for iconc


	if (params->susp_mass_output)
		Conc_compute_total_suspended_mass(c[0], grid, params);

	if (params->ave_height_output)
		Conc_compute_ave_height_x(c, grid, params);

	// Limit to capture the front location

	double front_limit = params -> Ly; // no-limit
	if (params->front_location_output)
		Conc_find_front_location(c, grid, params, front_limit);




	if (params->front_location_output)
		Output_timehistory(time, c[0]->W_total_susp_mass, c[0]->W_front_location, params);







}





/******************************************************************************/
/*
 This function computes total suspended mass within the fluid region
 */
/******************************************************************************/
void Conc_compute_total_suspended_mass(Concentration *c, MAC_grid *grid, Parameters *params) {

	int NX, NY, NZ;

	int i, j, k;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	double sum, ***c_data;

	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	Ie = min(NX-1, grid->G_Ie);
	Je = min(NY-1, grid->G_Je);
	Ke = min(NZ-1, grid->G_Ke);

	c_data = c->data;
	sum = 0.0;
	for(k = Ks; k < Ke; k++) {
		for(j = Js; j < Je; j++) {
			for(i = Is; i < Ie; i++) {
				double cellVolume = TwodOps_cell_measure_c(grid, params, i, j, k);
				sum += c_data[k][j][i]*cellVolume;
			}
		}
	}

	MPI_Allreduce (&sum, &c->W_total_susp_mass, 1, MPI_DOUBLE, MPI_SUM, PCW);
}






/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*
 Post processing functions
 */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



/******************************************************************************/
/*
This function computes the (concentration based) average height of the flow.
It is only a function of x since it has been averaged in y and z directions.

Note that, it is stored only in processor zero
*/
/******************************************************************************/
void Conc_compute_ave_height_x(Concentration **c, MAC_grid *grid,
	Parameters *params) {


int Is, Js, Ks;
int Ie, Je, Ke;
int i, j, k;
int NConc, iconc;
int NX;
double h_temp, coef;
double dy, dz;
double *yv, *zw;
double *dy_v, *dz_w;
double ***conc;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Js = grid->G_Js;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Je = grid->G_Je;
Ke = grid->G_Ke;

yv = grid->yv;
zw = grid->zw;

// Grid dimensions
dy_v = grid->dy_v;
dz_w = grid->dz_w;

NX = grid->NX;

coef = 1.0/(params->Ly * params->Lz);

NConc = params->NConc;
for (iconc=0; iconc<NConc; iconc++) {

	// Get concentration data
	conc = c[iconc]->data;

	/*--------------------------------------------------------------------*/
	/*
	 At each x location (i-index), find
	 h_ave = 1/(Ly*Lz) sum(over j and k or (y,z) indices) {C(k,j) * (dy*dz)}
	 */
	/*--------------------------------------------------------------------*/

	// This will give us the current height only as a function of
	for (i=Is; i<Ie; i++) {

		h_temp = 0.0;
		for (k=Ks; k<Ke; k++) {

			dz = dz_w[k];
			for (j=Js; j<Je; j++) {

				if (grid->c_status[k][j][i] == FLUID) {
					dy = dy_v[j];
					h_temp += conc[k][j][i]*dy*dz;
				} // if

			} // for
		} // for k

		c[iconc]->G_ave_height_x[i] = h_temp * coef;
	} // for i

	/*--------------------------------------------------------------------*/
	/*
	 Now, reduce all G_ave_height on processor zero (or all processors) and
	 add them all to find the correct ave_height.

	 Note that the trick is that on all processors, we create the full array
	 of ave_height_x[0->NX-1] and set all the parts not located on current
	 processor equal to zero (all the time) and just compute the part
	 located on current processor
	 */
	/*--------------------------------------------------------------------*/

	// Store the final result on processor zero W_ave_height
	MPI_Reduce( (void *)c[iconc]->G_ave_height_x, (void *)c[iconc]->W_ave_height_x,
	           NX, MPI_DOUBLE, MPI_SUM, MASTER, PCW);

	// This will send to all processors after computing ave_height
//		MPI_Allreduce( (void *)c[iconc]->G_ave_height_x, (void *)c[iconc]->W_ave_height_x, NX, MPI_DOUBLE, MPI_SUM, PCW);

} // for iconc
}




/******************************************************************************/
/*
This function finds the front location of each concentration field.  It is done
only on processor zero
*/
/******************************************************************************/
void Conc_find_front_location(Concentration **c, MAC_grid *grid,
	Parameters *params, double front_limit) {

Scalar heights[2];
Scalar height_front;
int iconc, NConc;
int i;
int NX;
double *xc;
double *ave_height;
FILE *fptr;
char front_location_filename[50];

if (params->rank == MASTER) {

	xc    = grid->xc;

	NConc = params->NConc;
	NX    = grid->NX;

	for (iconc=0; iconc<NConc; iconc++) {

		ave_height = c[iconc]->W_ave_height_x;

		for (i=NX-2; i>0; i--) {

			// Found front location
			if ( (ave_height[i] <= front_limit) && (ave_height[i-1] > front_limit) ) {

				/*--------------------------------------------------------*/
				/*
				 For the two point interpolation, Play a little trick since
				 we have front limit, but we want to find the location X.

				 h0 is found at x0. h1 is found at x1. Now, we want to find
				 xf corresponding to the front limit
				 */
				/*--------------------------------------------------------*/

				heights[0].X1  = ave_height[i-1];
				heights[1].X1  = ave_height[i];

				heights[0].Value = xc[i-1];
				heights[1].Value = xc[i];

				height_front.X1 = front_limit;

				// Just linear interpolation for the north and south conc
				// to find c_term at the v_grid node
				c[iconc]->W_front_location = MyMath_interpolate_quantity(heights, height_front, 2);

				break; // Goto next concentration field

			} // if
		} // for i

		// Write front location to file
		sprintf(front_location_filename, "front_location_%d.dat", iconc);
		if (params -> time == params -> default_dt){
			fptr = fopen(front_location_filename, "w");
		}
		else{
			fptr = fopen(front_location_filename, "a");
		}
		fprintf(fptr, "% .10g,% .10g\n",
		params->time, c[iconc]->W_front_location );
		fclose(fptr);
		printf("Write front location:%d %f\n", params->ntime, c[iconc]->W_front_location);
	} // for iconc
} // if rank
}




/******************************************************************************/
/*
This function integrates (in time) the deposited sediment height for the
particle concentration fields.  For now, just use it for cases without hindered
settling.
*/
/******************************************************************************/
void Conc_integrate_deposited_height(Concentration *c, MAC_grid *grid,
	Parameters *params, double weight_factor) {

int Is, Js, Ks;
int Ie, Je, Ke;
int i, k;
int j_index;
int **interface_y_index;
int NX, NZ;
double delta_h;
double dx, dz;
double *xu, *zw;
double *dx_u, *dz_w;
double ***conc, ***conc_old;

double dt = params -> dt;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Js = grid->G_Js;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Je = grid->G_Je;
Ke = grid->G_Ke;

xu = grid->xu;
zw = grid->zw;

// Grid dimensions
dx_u = grid->dx_u;
dz_w = grid->dz_w;

// y-index of first fluid node from bottom, (y=0)
interface_y_index = grid->interface_y_index;

// Get concentration data
conc = c->data;
conc_old = c->ng_data_old;

NX = grid->NX;
NZ = grid->NZ;

// Now go through 2D array and check if the first fluid node lies on
// current processor
for (k=Ks; k<Ke; k++) {

	dz = dz_w[k];
	for (i=Is; i<Ie; i++) {

		// Index of first interior fluid node (on the bottom boundary
		j_index = interface_y_index[k][i];

		dx = dx_u[i];
		// Check if current prcossor includes the bottom geometry
		if ( ( j_index >= Js) && (j_index < Je) ){

			// delta_h added height at each time step: second order in time
			delta_h = dt * (dx*dz) * fabs(c->v_settl0) * weight_factor *
			          0.5 * (conc[k][j_index][i]+conc_old[k][j_index][i]);

			// Update the height of deposited sediment
			c->G_deposit_height[k][i] += delta_h;

		} // if

	} // for
} /* for k*/

}




/******************************************************************************/
/*
This function communicates all the G_deposit_height amongst processors to find
W_deposit height.

Result could be stored either on processor zero or all processors.
*/
/******************************************************************************/
void Conc_update_world_deposited_height(Concentration *c, Cart3d_bag *data_bag) {

int Is, Js, Ks;
int Ie, Je, Ke;
int iconc, NConc;
Indices G_s, G_e, W_e;

MAC_grid *grid = data_bag -> grid;
Parameters *params = data_bag -> params;
NConc = params->NConc;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Ke = grid->G_Ke;

//  Since the 2D arrays are assumed to have [Y][X] index order, pass always
// the min, max indices based on this rule
// Start and End indices of the 2D array on the current processor
G_s.x_index = Is;
G_s.y_index = Ks;

G_e.x_index = Ie;
G_e.y_index = Ke;

// Total number of grid point in the W_ array
W_e.x_index = grid->NX;
W_e.y_index = grid->NZ;

// Go through each particle concentration field and Reduce all the 2D
// arrays on processor zero

		Communication_reduce_2D_arrays(c->G_deposit_height, c->W_deposit_height,
		                               &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);



}




/******************************************************************************/
/*
This function communicates all the G_deposit_height_dumped amongst processors
to find W_deposit height_dumped.

Result could be stored either on processor zero or all processors
*/
/******************************************************************************/
void Conc_update_world_deposited_height_dumped(Concentration **c, Cart3d_bag *data_bag) {

int Is, Js, Ks;
int Ie, Je, Ke;
int iconc, NConc;
Indices G_s, G_e, W_e;

MAC_grid *grid = data_bag -> grid;
Parameters *params = data_bag -> params;
NConc = params->NConc;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Ke = grid->G_Ke;

// Since the 2D arrays are assumed to have [Y][X] index order, pass always
// the min, max indices based on this rule
// Start and End indices of the 2D array on the current processor
G_s.x_index = Is;
G_s.y_index = Ks;

G_e.x_index = Ie;
G_e.y_index = Ke;

// Total number of grid point in the W_ array
W_e.x_index = grid->NX;
W_e.y_index = grid->NZ;

// Go through each particle concentration field and Reduce all the 2D
// arrays on processor zero
for (iconc=0; iconc<NConc; iconc++) {

	if (params->conc_output_deposit_height[iconc]) {

		// Communicate to add all the G_deposit_height[] to processor zero
		// for each concentration field
		Communication_reduce_2D_arrays(c[iconc]->G_deposit_height_dumped, c[iconc]->W_deposit_height_dumped,
		                               &G_s, &G_e, &W_e, REDUCE_TO_MASTER, data_bag);
	} // if
} // for iconc

}




/******************************************************************************/
/*
This function adds all the available particles in the domain to the deposited
height.  This is useful when the velocity field is small and the only important
thing happening in the domain is the settling of particles.  This would somehow
predict how the deposited height would be after a long run calculation.
*/
/******************************************************************************/
void Conc_dump_particles(Concentration **c, Cart3d_bag *data_bag) {

int Is, Js, Ks;
int Ie, Je, Ke;
int i, j, k;
int iconc, NConc;
int NX, NZ;
double ***conc;
double *xu, *yv, *zw;
double *dx_u, *dy_v, *dz_w;
double dx, dy, dz;
double h_temp;
double weight_factor;

MAC_grid *grid = data_bag -> grid;
Parameters *params = data_bag -> params;

xu = grid->xu;
yv = grid->yv;
zw = grid->zw;

// Grid dimensions
dx_u = grid->dx_u;
dy_v = grid->dy_v;
dz_w = grid->dz_w;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Js = grid->G_Js;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Je = grid->G_Je;
Ke = grid->G_Ke;

NX = grid->NX;
NZ = grid->NZ;

NConc = params->NConc;

// Integrate all the particles in y-direction and add them to the
// sed_height for the particle fields
for (iconc=0; iconc<NConc; iconc++) {

	if (params->conc_output_dump[iconc]) {

		conc = c[iconc]->data;

		// At each location (x,z) find c_t = integral (0 to Ly) of
		// c_dx.dy.dz
		weight_factor = params->richardson[iconc]/(params->grav[1] );

		// This will give us the current height only as a function of x
		for (k=Ks; k<Ke; k++) {

			dz = dz_w[k];
			for (i=Is; i<Ie; i++) {

				dx     = dx_u[i];
				h_temp = 0.0;

				c[iconc]->G_deposit_height_dumped[k][i] = c[iconc]->G_deposit_height[k][i];
				for (j=Js; j<Je; j++) {

					if (grid->c_status[k][j][i] == FLUID) {

						dy = dy_v[j];
						h_temp += conc[k][j][i]*dy;
					} // if
				} // for j

				// Now add the integrated value for the shrinked value of
				// concentration to the deposited one and store the result
				// into the two dimensional array locally on each processor
				c[iconc]->G_deposit_height_dumped[k][i] += h_temp *(dx*dz) * weight_factor;

			} // for
		} // for

	} // if particle

} // for iconc


// Now, update the W_dep... on processor zero
Conc_update_world_deposited_height_dumped(c, data_bag);
}




/******************************************************************************/
/*
This function computes the stokes dissipation for any particle concentration
field.  This dissipation is due to the fact that the Stokes flow around
settling particles cause some microscopic dissipation of energy:
 e_s = int( u_s * c dV) domain
*/
/******************************************************************************/
void Conc_compute_stokes_dissipation_rate( int iconc, Concentration **c, MAC_grid *grid,
	Parameters *params) {

int i, j, k;
int Is, Js, Ks;
int Ie, Je, Ke;
double dV;
double ***conc;
double G_stokes_dissipation_rate;

// Get concentration data
conc = c[iconc]->data;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Js = grid->G_Js;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Je = grid->G_Je;
Ke = grid->G_Ke;

G_stokes_dissipation_rate = 0.0;

if ( params->conc_output_stokes_diss_rate[iconc]){
	for (k=Ks; k<Ke; k++) {
		for (j=Js; j<Je; j++) {
			for (i=Is; i<Ie; i++) {

				// only include if point is fluid
				if (grid->c_status[k][j][i] == FLUID) {
					/* Physical cell measure in the active 2D/3D geometry. */
					dV = TwodOps_cell_measure_c(grid, params, i, j, k);

					G_stokes_dissipation_rate += conc[k][j][i] * dV;
				} // if

			} // for i
		} // for j
	}// for k
} // if


// Now, multiply by the constant settling speed
G_stokes_dissipation_rate *= fabs(c[iconc]->v_settl0);

// Store the rate for the current processor
c[iconc]->G_stokes_dissipation_rate = G_stokes_dissipation_rate;

}




/******************************************************************************/
/*
This function computes the stokes dissipation for each particle concentration
field pn each processor and then adds them up to get the W_... values
*/
/******************************************************************************/
void Conc_update_world_stokes_dissipation_rate( Concentration **c,
	MAC_grid *grid, Parameters *params) {

int iconc, NConc;
double G_stokes_dissipation_rate;
double W_stokes_dissipation_rate;

NConc = params->NConc;

G_stokes_dissipation_rate = 0.0;

// Add all the concentration fields dissipation on the current processor
for (iconc=0; iconc<NConc; iconc++) {

	// Compute the stokes dissipation for each conc field
	Conc_compute_stokes_dissipation_rate(iconc, c, grid, params);

	// Now, add them up with the contribution factor
	G_stokes_dissipation_rate += c[iconc]->G_stokes_dissipation_rate * params->richardson[iconc]/(params->grav[1] );
} // for iconc

// Now add them up and store them in the first concentration field
MPI_Allreduce( &G_stokes_dissipation_rate, &W_stokes_dissipation_rate, 1, MPI_DOUBLE, MPI_SUM, PCW);

c[0]->W_stokes_dissipation_rate = W_stokes_dissipation_rate;

}




/******************************************************************************/
/*
This function computes the potential energy for concentration fields within the
interior region. This potential energy could be converted to kinetic energy.
*/
/******************************************************************************/
void Conc_compute_active_potential_energy(Concentration *c, MAC_grid *grid,
	Parameters *params) {

int i, j, k;
int Is, Js, Ks;
int Ie, Je, Ke;
double dV;
double *yc;
double ***conc;
double G_Ep_active;

// Get concentration data
conc = c->data;

yc = grid->yc;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Js = grid->G_Js;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Je = grid->G_Je;
Ke = grid->G_Ke;

G_Ep_active = 0.0;

for (k=Ks; k<Ke; k++) {
	for (j=Js; j<Je; j++) {
		for (i=Is; i<Ie; i++) {

			// only include if point is fluid
			if (grid->c_status[k][j][i] == FLUID) {
				dV = TwodOps_cell_measure_c(grid, params, i, j, k);

				G_Ep_active += conc[k][j][i] * dV * yc[j];
			} // if

		} // for i
	} // for j
}// for k

// Store the rate for the current processor
c->G_Ep_active = G_Ep_active;


}




/******************************************************************************/
/*
This function computes the potential energy for the particles that have settled
out. This potential energy could not be converted to kinetic energy and is
computed just for the cases where we have variable topography.
*/
/******************************************************************************/
void Conc_compute_passive_potential_energy(Concentration *c, MAC_grid *grid,
	Parameters *params) {

double G_Ep_passive;
double dx, dz;
int i, k;
int j_index;
int Is, Js, Ks;
int Ie, Je, Ke;
double *xu, *zw, *yc;
double *dx_u, *dz_w, *dy_v;

// Grid coordinates
xu = grid->xu;
zw = grid->zw;
yc = grid->yc;

// Grid dimensions
dx_u = grid->dx_u;
dy_v = grid->dy_v;
dz_w = grid->dz_w;

// Start index of bottom-left-back corner on current processor
Is = grid->G_Is;
Js = grid->G_Js;
Ks = grid->G_Ks;

// End index of top-right-front corner on current processor
Ie = grid->G_Ie;
Je = grid->G_Je;
Ke = grid->G_Ke;

// Add the potential energy due to the deposited particles * on the bottom
// geometry
G_Ep_passive = 0.0;


	for (k=Ks; k<Ke; k++) {

		dz = dz_w[k];
		for (i=Is; i<Ie; i++) {

			dx = dx_u[i];

			// y-index of first interior fluid node
			j_index = grid->interface_y_index[k][i];

			// Ony add the contribution if the bottom surface lies on the
			// current processor
			if ( (j_index >= Js) && (j_index < Je)) {

				G_Ep_passive += c->G_deposit_height[k][i] * dx*dz * yc[j_index];
			} // if

		} // for i
	}// for k


// Store the rate for the current processor
c->G_Ep_passive = G_Ep_passive;

}




/******************************************************************************/
/*
This function computes the W_... potential energies and sends them back to all
processors
*/
/******************************************************************************/
void Conc_update_world_potential_energies(Concentration **c, MAC_grid *grid,
	Parameters *params) {

int iconc, NConc;
double G_total_Ep_passive;
double G_total_Ep_active;
double *send_buffer;
double *recv_buffer;
int nt;

NConc = params->NConc;
G_total_Ep_passive = 0.0;
G_total_Ep_active  = 0.0;


for (iconc=0; iconc<NConc; iconc++) {

	Conc_compute_active_potential_energy(c[iconc], grid, params) ;
	Conc_compute_passive_potential_energy(c[iconc], grid, params) ;

	G_total_Ep_active  += c[iconc]->G_Ep_active  * params->richardson[iconc]/(params->grav[1] );


	// No need for Conc_alpha since it already was considered in computing
	// the sed_height
	G_total_Ep_passive += c[iconc]->G_Ep_passive;

} // for iconc


// Total number of quantities which have to be reduced globally
nt = 2;
// First, allocate required memory
send_buffer = Memory_allocate_1D_array(GVG_DOUBLE, nt);
recv_buffer = Memory_allocate_1D_array(GVG_DOUBLE, nt);

send_buffer[0] = G_total_Ep_active;
send_buffer[1] = G_total_Ep_passive;

// Now add up all the local values of the energies them back to store it in
// the W_... variables
MPI_Allreduce( (void *)send_buffer, (void *)recv_buffer, nt, MPI_DOUBLE, MPI_SUM, PCW);

// Store the result in the first concentration field
c[0]->W_Ep_active  = recv_buffer[0];
c[0]->W_Ep_passive = recv_buffer[1];

free(send_buffer);
free(recv_buffer);

}




