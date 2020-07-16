#ifndef ARRAY_H
	#define ARRAY_H

void Array_copy_noghost(double ***data1, double ***data2,     MAC_grid *grid, Parameters *params );
void Array_copy_withghost(double ***data1, double ***data2,   MAC_grid *grid, Parameters *params );
void Array_scale_noghost(double ***data, double factor,       MAC_grid *grid, Parameters *params );
void Array_scale_withghost(double ***data, double factor,     MAC_grid *grid, Parameters *params );
void Array_AXPY_noghost(double ***y, double a, double ***x,   MAC_grid *grid, Parameters *params );
void Array_AXPY_withghost(double ***y, double a, double ***x, MAC_grid *grid, Parameters *params );
void Array_pointwisemult_noghost(double ***w, double ***x, double ***y, MAC_grid *grid, Parameters *params);
void Array_pointwisemult_withghost(double ***w, double ***x, double ***y, MAC_grid *grid, Parameters *params);
void Array_AXPBY_withghost(double ***y, double alpha, double beta, double ***x, MAC_grid *grid, Parameters *params);
void Array_AXPBY_noghost(double ***y, double alpha, double beta, double ***x, MAC_grid *grid, Parameters *params);
void Array_XXPY_withghost(double ***y, double ***x, MAC_grid *grid, Parameters *params);
void Array_XXPY_noghost(double ***y, double ***x, MAC_grid *grid, Parameters *params);
void Array_AXXPY_withghost(double ***y, double alpha, double ***x, MAC_grid *grid, Parameters *params);
void Array_AXXPY_noghost(double ***y, double alpha, double ***x, MAC_grid *grid, Parameters *params);
void Array_XZPY_withghost(double ***y, double ***x, double ***z, MAC_grid *grid, Parameters *params);
void Array_XZPY_noghost(double ***y, double ***x, double ***z, MAC_grid *grid, Parameters *params) ;
void Array_AXZPY_withghost(double ***y, double alpha, double ***x, double ***z, MAC_grid *grid, Parameters *params);
void Array_AXZPY_noghost(double ***y, double alpha, double ***x, double ***z, MAC_grid *grid, Parameters *params);


void Array_set_noghost(double ***data, double alpha, MAC_grid *grid, Parameters *params); 
void Array_set_withghost(double ***data, double alpha, MAC_grid *grid, Parameters *params); 

	#ifdef PETSC
void Array_copy_to_vec(DM *DA_array, Vec *global, double ***data, MAC_grid *grid, Parameters *params) ;
	#endif

#endif // notARRAY_H
