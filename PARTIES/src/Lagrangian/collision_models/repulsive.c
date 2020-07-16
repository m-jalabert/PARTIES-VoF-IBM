
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"

/******************************************************************************/
/*
 Assumes uniform square grid
 */
/******************************************************************************/
void repulsive(double *normal_force, Collision_bag *bag, Parameters *params) {

	double spring_force;

	double *n = bag -> n;
	double surface_distance = bag -> surface_distance;

	spring_force = params->kn * surface_distance * surface_distance;

	normal_force[0] = -n[0] * spring_force;
	normal_force[1] = -n[1] * spring_force;
	normal_force[2] = -n[2] * spring_force;
}




/******************************************************************************/
/*
 Assumes uniform square grid
 */
/******************************************************************************/
double repulsive_wall(double surface_distance, int dim, int side, Parameters *params) {

	double normal_force;

	normal_force = ((double)(side)) * params->kn * surface_distance * surface_distance;

	return normal_force;
}
