#ifndef SCHUMANN_H
 #define SCHUMANN_H

Wall_model *Schumann_create(MAC_grid *grid, Parameters *params);
void Schumann_calc_shear(double ***u_data, double ***v_data, double ***w_data, 
		Wall_model *log_law, MAC_grid *grid, Parameters *params);
double Schumann_ypbufl(double kappa, double b);
double Schumann_calc_utau(double u, double y, double iRe, double kappa, double bllaw, double ypintr);
#endif
