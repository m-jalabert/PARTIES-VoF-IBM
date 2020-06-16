#ifndef MEMORY_H
 #define MEMORY_H

void Memory_check_allocation(void *var);
double*** Memory_allocate_flow_variable(MAC_grid *grid, Parameters *params);
void Memory_reset_flow_variable(MAC_grid *grid, Parameters *params, double ***array);
void Memory_free_flow_variable(MAC_grid *grid, Parameters *params, double ***array);
int ***Memory_allocate_flow_variable_int(MAC_grid *grid, Parameters *params);
void  Memory_free_flow_variable_int(MAC_grid *grid, Parameters *params, int ***array);
double*** Memory_allocate_noghost_variable(MAC_grid *grid, Parameters *params);
void Memory_reset_noghost_variable(MAC_grid *grid, Parameters *params, double ***array);
void Memory_free_noghost_variable(MAC_grid *grid, Parameters *params, double ***array);
double ***Memory_allocate_variable_jik(MAC_grid *grid, Parameters *params);
void  Memory_free_variable_jik(MAC_grid *grid, Parameters *params, double ***array);
void *Memory_allocate_1D_array(int data_type, int N);
double **Memory_allocate_2D_double_array(int NY, int NZ);
double complex **Memory_allocate_2D_double_complex_array(int NY, int NZ);
int **Memory_allocate_2D_int_array(int NY, int NZ);
char **Memory_allocate_2D_char_array(int NY, int NZ);
double ***Memory_allocate_3D_double_array(int NX, int NY, int NZ);
void Memory_free_2D_double_array(int NZ, double **array);
void Memory_free_2D_double_complex_array(int NZ, double **array);
void Memory_free_2D_int_array(int NZ, int **array);
void Memory_free_2D_char_array(int NZ, char **array);
void Memory_free_3D_double_array(int NY, int NZ, double ***array);

#endif
