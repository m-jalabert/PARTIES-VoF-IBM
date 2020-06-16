#ifndef H_RESUME
	#define H_RESUME
#include "hdf5.h"
#include "DataTypes.h"

void Resume_h5_resume(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Resume_h5_data(Cart3d_bag *data_bag, Debug_trace *dtrace);

void Resume_h5_dataset(void *data, hid_t mem_type, int ndim, hsize_t *dim,
		hid_t file_id, char *fieldname, Parameters *params, Debug_trace *dtrace);
void Resume_h5_timer(Timer *timer, hid_t file_id, char *fieldname,
	Parameters *params, Debug_trace *dtrace);
void Resume_h5_flow_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace);
void Resume_h5_noghost_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace);

#endif
