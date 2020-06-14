#ifndef DISPLAY_H
 #define DISPLAY_H

void Display_DA_3D_info(MAC_grid *grid, Parameters *params);
void Display_2D_outflow(double **outflow, MAC_grid *grid, Parameters *params, char *q_name, char which_quantity);
void Display_2D_inflow(double **inflow, MAC_grid *grid, Parameters *params, char *q_name, char which_quantity);
int Display_parameters(Parameters *params);
void Display_time_results(Cart3d_bag *data_bag);
void Display_flag(int flag, char *name);
void Display_point(PointType *p, char *name);
void Display_immersed_node(ImmersedNode *ib_node);
void Display_grid(double *data, int N, char *name);
void Display_progress(Parameters *params, char *statement);

void Display_assert_list_state(Particle_list *p_list, int state,
		Parameters *params, Debug_trace *dtrace);
void Display_print_list_name(char *name, int state);
void Display_print_list_type(char *name, int type);
void Display_print_dtrace_list(char *message_out, Debug_trace *dtrace);
void Display_assert_error(int status, const char *message_in,
		Parameters *params, Debug_trace *dtrace);
void Display_throw_error(const char *message_in, Parameters *params,
		Debug_trace *dtrace);
void Display_throw_warning(const char *message, Parameters *params);
void Display_error(const char *message_in);

Debug_trace *Display_create_debug_trace_child(Debug_trace *dtrace,
		const char *function, const char *file, const int line);

#endif
