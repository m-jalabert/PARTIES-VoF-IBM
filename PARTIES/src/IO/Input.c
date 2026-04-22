/*
 * config.c
 *
 *  Created on: Oct 17, 2013
 *      Author: daan
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "DataTypes.h"

#include "Display.h"
#include "ini.h"
#include "Input.h"
#include "Memory.h"

static void Input_set_twod_mode_parameters(Parameters *params)
{
	params->twod_mode_enabled = NO;
	params->twod_cartesian_enabled = NO;
	params->axisym_rz_enabled = NO;
	params->axisym_no_swirl = NO;

#ifdef TWOD_CARTESIAN
	params->twod_mode_enabled = YES;
	params->twod_cartesian_enabled = YES;
#endif
#ifdef AXISYM_RZ
	params->twod_mode_enabled = YES;
	params->axisym_rz_enabled = YES;
#endif
#ifdef AXISYM_NO_SWIRL
	params->axisym_no_swirl = YES;
#endif

	if (params->axisym_theta_cells <= 0)
		params->axisym_theta_cells = 1;

	if (params->axisym_theta_span <= 0.0)
		params->axisym_theta_span = TWOD_AXISYM_THETA_SPAN_FULL;

	if (!params->twod_mode_enabled) {
		if (params->twod_slab_thickness <= 0.0)
			params->twod_slab_thickness = params->Lz;
		return;
	}

	{
		double hx = params->Lx / (double)params->NXM;
		double hy = params->Ly / (double)params->NYM;
		double slab_thickness = params->twod_slab_thickness;
		double h_ref_xy = hx;

		if (hy < h_ref_xy)
			h_ref_xy = hy;

		if (slab_thickness <= 0.0)
			slab_thickness = h_ref_xy;
		if (slab_thickness <= 0.0)
			slab_thickness = 1.0;

		params->twod_slab_thickness = slab_thickness;
		params->Lz = slab_thickness;
		params->zmax = params->zmin + slab_thickness;
	}
}



/******************************************************************************/
/*
 Read in default and overriden values to store to Parameters structure
 */
/******************************************************************************/
void Input_set_parameters(Parameters *params, char *infile)
{

	int verbose = 0;

	// Set default values
	#define CFG(s, n, default, reader) params->n = default;
	#define CFGString(s, n, default) sprintf(params->n, default);
	#define CFGArray(s, n, default, data_type) params->n = Input_handle_array(default, data_type, #n);
	if (verbose)
		printf("Attempting to load default values in default.inp\n");
	#include "default.inp"
	#undef CFG
	#undef CFGString
	#undef CFGArray
	if (verbose)
		printf("Loaded default values in default.inp\n");

	if (verbose)
		printf("Attempting to load input file %s\n", infile);
	if (ini_parse(infile, Input_handler, params) < 0) {
		fprintf(stderr, "Can't load input file %s, using defaults\n", infile);
	}
	else if(verbose) {
		printf("Successfully loaded input file %s\n", infile);
	}

	params->Lx = params->xmax - params->xmin;
	params->Ly = params->ymax - params->ymin;
	params->Lz = params->zmax - params->zmin;

	Input_set_twod_mode_parameters(params);

	{
		double hx = params->Lx / (double)params->NXM;
		double hy = params->Ly / (double)params->NYM;
		double hz = params->Lz / (double)params->NZM;
		double h_ref = hx;

		if (hy < h_ref)
			h_ref = hy;
		if (hz < h_ref)
			h_ref = hz;

		if (params->Cn <= 0.0)
			params->Cn = 0.75 * h_ref;

		if (params->Pe_CH <= 0.0 && params->Cn > 0.0)
			params->Pe_CH = 0.9 / params->Cn;
	}

	if (params->front_location_output == 1 && params->ave_height_output == 0) {
		Display_throw_warning("Warning: specified 'front_location_output' but not "
			"'ave_height_output'.\nEnabling 'ave_height_output'...\n", params);
		params->ave_height_output = 1;
	}
}


/******************************************************************************/
/*
 Process a line of the INI file, storing valid values into config struct
 */
/******************************************************************************/
int Input_handler(void *user, const char *section, const char *name,
		const char *value)
{
	int status;
	Parameters *params = (Parameters *)user;

	if (0) {}
	#define CFG(s, n, default, reader) \
		else if (strcmp(section, #s)==0 && strcmp(name, #n)==0) \
			params->n = reader(value);
	#define CFGString(s, n, default) \
		else if (strcmp(section, #s)==0 && strcmp(name, #n)==0) \
			sprintf(params->n, value);
	#define CFGArray(s, n, default, data_type) \
		else if (strcmp(section, #s)==0 && strcmp(name, #n)==0) { \
			free(params->n); \
			params->n = Input_handle_array(value, data_type, name); \
		}
	#include "default.inp"
	#undef CFG
	#undef CFGString
	#undef CFGArray

	return 1;
}

/******************************************************************************/
/*
 Parse a string of numbers 'values' into an array whose pointer is returned.
 The string should be formatted as follows:
     {number1, number2, number3}
 white space does not matter.
 */
/******************************************************************************/
void *Input_handle_array(const char *values, int data_type, const char *name)
{
	int i, N, finished;
	char *start, *end;
	int *temp_int;
	double *temp_double;
	if (data_type == GVG_INT)
		temp_int = (int *)Memory_allocate_1D_array(data_type, MAX_ARRAY_LENGTH);
	else if (data_type == GVG_DOUBLE)
		temp_double = (double *)Memory_allocate_1D_array(data_type, MAX_ARRAY_LENGTH);
	else {
		fprintf(stderr, "Improper dataype for input \"%s\" in input file.  Use GVG_INT or GVG_DOUBLE.\n", name);
		return NULL;
	}


	char values_temp[4 * MAX_ARRAY_LENGTH];
	strcpy(values_temp, values);

	end = values_temp;
	if (*end != '{') {
		fprintf(stderr, "Improper syntax for input \"%s\" in input file.  Use {*,*,*}.\n", name);
		return NULL;
	}
	N = 0;
	finished = 0;

	while (!finished) {
		start = end + 1;
		start = Input_lskip(start);
		end = Input_find_char_next(start);
		if (*end == ',' || *end == '}') {
			if (*end == '}')
				finished = 1;
			*end = '\0';
			if (data_type == GVG_INT)
				temp_int[N] = atoi(start);
			else if (data_type == GVG_DOUBLE)
				temp_double[N] = atof(start);
		}
		else {
			fprintf(stderr, "Improper syntax for input \"%s\" in input file.  Use {*,*,*}.\n", name);
			return NULL;
		}
		N++;
	}

	if (data_type == GVG_INT) {
		int *array = (int *)Memory_allocate_1D_array(data_type, N);
		for (i = 0; i < N; i++)
			array[i] = temp_int[i];
		free(temp_int);
		return array;
	}
	else if (data_type == GVG_DOUBLE) {
		double *array = (double *)Memory_allocate_1D_array(data_type, N);
		for (i = 0; i < N; i++)
			array[i] = temp_double[i];
		free(temp_double);
		return array;
	}
	else {
		fprintf(stderr, "Unrecognized datatype.\n");
		return NULL;
	}
}



/******************************************************************************/
/*
 Return pointer to first non-whitespace char in given string.
 */
/******************************************************************************/
static char* Input_lskip(const char* s)
{
	while (*s && isspace((unsigned char)(*s)))
		s++;
	return (char*)s;
}



/******************************************************************************/
/*
 Return pointer to first ',' or '}' or ';' comment in given string, or pointer
 to null at end of string if neither found. ';' must be prefixed by a whitespace
 character to register as a comment.
 */
/******************************************************************************/
static char* Input_find_char_next(const char* s)
{
	int was_whitespace = 0;
	while (*s && *s != ',' && *s != '}' && !(was_whitespace && *s == ';')) {
		was_whitespace = isspace((unsigned char)(*s));
		s++;
	}
	return (char*)s;
}
