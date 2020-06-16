#ifndef dsolver_H
 #define dsolver_H

void lsolver_transpose_setup(MAC_grid *grid, Parameters *params); 
void ltridiag(double *a, double *b, double *c, double *r, double *u, int N); 


#endif
