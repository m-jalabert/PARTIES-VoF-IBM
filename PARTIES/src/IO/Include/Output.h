#ifndef OUTPUT_H
 #define OUTPUT_H
#include "hdf5.h"

void Output_h5_resume(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Output_h5_data(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Output_h5_data2d(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Output_h5_timer(Timer *timer, hid_t file_id, char *fieldname,
		Parameters *params, Debug_trace *dtrace);
void Output_h5_flow_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace);
void Output_h5_noghost_variable(double ***data, hid_t file_id, char *fieldname,
		MAC_grid *grid, Parameters *params, Debug_trace *dtrace);
void Output_h5_3d_variable(double *data_start, int ndim, hsize_t *dim,
		hsize_t *data_count, hsize_t *data_offset,
		hsize_t *mem_count,  hsize_t *mem_offset,
		hid_t file_id, char *fieldname, Parameters *params, Debug_trace *dtrace);
void Output_h5_dataset(void *data_start, hid_t datatype, int ndim, hsize_t *dim,
		hid_t file_id, char *fieldname, Parameters *params, Debug_trace *dtrace);
//void assert(herr_t status, Parameters *params, char *message);

hid_t Output_h5_open(char *filename, Parameters *params, Debug_trace *dtrace);
void Output_h5_create_group(hid_t file_id, char *groupname, Parameters *params,
		Debug_trace *dtrace);

void Output_immersed_info(MAC_grid *grid, Parameters *params, char which_quantity);
void Output_immersed_boundary(MAC_grid *grid, Parameters *params, char which_quantity);
void Output_write_ascii_deposit_height_dumped(Concentration **c,
		Parameters *params, MAC_grid *grid, int noutput);
void Output_check_file_open(FILE *file);
void Output_timehistory(double time, double frontLocation,
		double averageHeight, Parameters *params);
void Output_3d_data(double*** data3d, char node_type, char *ghost_type,
		char *filename, Cart3d_bag *data_bag, Debug_trace *dtrace);

void Output_2d_data(double** data2d, int *dim, double *x, double *y, double *z,
		char *filename, Cart3d_bag *data_bag, Debug_trace *dtrace);

void slice_2d_output(double ***data3d, char nme, Cart3d_bag *data_bag, double time, Debug_trace *dtrace);		

#endif
