#ifndef GRID_H
 #define GRID_H

MAC_grid *Grid_create(Parameters *params, Debug_trace *dtrace);
void Grid_destroy(MAC_grid *grid, Parameters *params);
void Grid_generate_mesh_coordinates(MAC_grid *grid, Parameters *params,
		Debug_trace *dtrace);
void Grid_compute_metric_coefficients(MAC_grid *grid, Parameters *params, char which_quantity);
void Grid_describe_interface(MAC_grid *grid, Parameters *params);
int Grid_import_bottom_interface(MAC_grid *grid, Parameters *params, char *filename);
void Grid_tag_q_nodes_using_surface(MAC_grid *grid, Parameters *params, char which_quantity);
void Grid_tag_q_box_boundary_nodes(MAC_grid *grid, char which_quantity, Parameters *params);
void Grid_identify_geometry(Cart3d_bag *data_bag);

int Grid_get_x_index(double x, MAC_grid *grid, char which_quantity);
int Grid_get_y_index(double y, MAC_grid *grid, char which_quantity);
int Grid_get_z_index(double z, MAC_grid *grid, char which_quantity);

int Grid_import_grid_from_file(MAC_grid *grid, Parameters *params);
void Grid_set_interface_y_index(MAC_grid *grid, Parameters *params);
void Grid_find_best_NP(int *NP, double *err_NP, int first_sign, int second_sign, Parameters *params);
#endif
