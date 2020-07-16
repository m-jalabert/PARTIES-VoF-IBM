#ifndef MYMATH_H
 #define MYMATH_H

double MyMath_interpolate_quantity(Scalar *q, Scalar q_int, int N);
double MyMath_random(double max);
double MyMath_get_point_point_distance(PointType *p1, PointType *p2) ;
double MyMath_matrix_determ(MatrixType *mat) ;
void MyMath_matrix_inv(MatrixType *mat) ;
void MyMath_GaussJordan_matrix_inv(MatrixType *mat) ;
double MyMath_do_trilinear_inter(double x_int, double y_int, double z_int, double ***data, MAC_grid *grid, char which_quantity) ;

#endif
