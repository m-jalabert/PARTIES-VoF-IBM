#ifndef H_IMMERSED
	#define H_IMMERSED 
	
Immersed *Immersed_create(MAC_grid *grid, Parameters *params, char which_quantity);
void Immersed_destroy(Immersed *q_immersed, Parameters *params, MAC_grid *grid) ;
void Immersed_copy_to_nut(MAC_grid *grid, Parameters *params); 
int Immersed_count_immersed_nodes(MAC_grid *grid, Parameters *params, char which_quantity) ;
void Immersed_set_q_immersed_indices(MAC_grid *grid, Parameters *params, char which_quantity);
int Immersed_get_ib_global_index(Immersed *q_immersed, int x_index, int y_index, int z_index);
void Immersed_set_q_immersed_coordinates(MAC_grid *grid, Parameters *params, char which_quantity);
void Immersed_setup_q_immersed_nodes(MAC_grid *grid, Parameters *params, char which_quantity);
void Immersed_set_q_immersed_interpolation(MAC_grid *grid, Parameters *params, char which_quantity);
int Immersed_compute_nearest_plane_interpolation(ImmersedNode *ib_node, MAC_grid *grid, Parameters *params, int boundary_condition, char which_quantity); 
int Immersed_compute_intersecting_xyz_plane_interpolation(ImmersedNode *ib_node, MAC_grid *grid, Parameters *params, int boundary_condition, char which_quantity, int shift); 
ImmersedNode *Immersed_get_ib_node(Immersed *q_immersed, int i);
int Immersed_line_plane_intersection(PointType *line, PointType *plane, double *t, double *u, double *v); 
int Immersed_find_plane_interpolation_coeff(Indices *ind_neighbor,
		PointType *interxn, int idir, char which_quantity, MAC_grid *grid,
		PointType *plane, Indices *ind_plane, double *u, double *v,
		int i_immersed, int j_immersed, int k_immersed, Parameters *params);
void Immersed_line_xplane_intersection_point(PointType *line, double x, PointType *intersection_point);
void Immersed_line_yplane_intersection_point(PointType *line, double y, PointType *intersection_point);
void Immersed_line_zplane_intersection_point(PointType *line, double z, PointType *intersection_point);
int Immersed_point_location_on_xplane(PointType *point, PointType *plane, double *u, double *v);
int Immersed_point_location_on_yplane(PointType *point, PointType *plane, double *u, double *v);
int Immersed_point_location_on_zplane(PointType *point, PointType *plane, double *u, double *v);
int Immersed_point_location_on_xline(PointType *point, PointType *line, double *u);
int Immersed_point_location_on_yline(PointType *point, PointType *line, double *u);
int Immersed_point_location_on_zline(PointType *point, PointType *line, double *u);
void Immersed_find_surface_normal_distance(MAC_grid *grid, Parameters *params, char which_quantity); 
void Immersed_find_boundary_point_and_normal_vector(MAC_grid *grid, Parameters *params, char which_quantity); 
double Immersed_scoordinate_of_boundary(double x, double y, double ax, double bx, double cx, double ay, double by, double cy);  
void Immersed_velocity_interpolation(Velocity *vel, MAC_grid *grid,
		Parameters *params);
void Immersed_conc_interpolation(MAC_grid *grid, Parameters *params,
		Concentration *c);
void Immersed_nut_interpolation(MAC_grid *grid, Parameters *params, Subgrid *smag); 
void Immersed_conc_solid(Concentration *c, Cart3d_bag *data_bag);

#endif
