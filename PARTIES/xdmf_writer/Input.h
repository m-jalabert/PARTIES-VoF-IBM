#ifndef INPUT_H
#define INPUT_H
#include "DataTypes.h"
#include "hdf5.h"

int Input_get_values(Parameters *params, Hyperslab *hyperslab);
void Input_set_defaults(Parameters *params, Hyperslab *hyperslab);
void Input_set_values(Parameters *params, Hyperslab *hyperslab);
void Input_set_values_otherData(Parameters *params, Hyperslab *hyperslab);
void Input_set_values_particle(Parameters *params, Hyperslab *hyperslab);
void Input_set_values_marker(Parameters *params, Hyperslab *hyperslab);
int handler_params(void* user, const char* section, const char* name,
		const char* value);
int handler_hs(void* user, const char* section, const char* name,
		const char* value);
herr_t handler_particle(hid_t g_id, const char *name, const H5L_info_t *info,
		void *op_data);
herr_t handler_otherData(hid_t g_id, const char *name, const H5L_info_t *info,
		void *op_data);
#endif
