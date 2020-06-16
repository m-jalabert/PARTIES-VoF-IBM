#ifndef CONC_H
 #define CONC_H



void computeStatistics(Velocity *u, Velocity *v, Velocity *w, Pressure *pressure, Concentration **conc, MAC_grid *grid, Parameters *params);
double kineticEnergy(double ***u_cell, double ***v_cell, double ***w_cell, MAC_grid *grid);
double advectiveFlux(double ***v_cell, double ***conc, MAC_grid *grid);
double* viscousDissipation_u(double ***u, MAC_grid *grid, double iRe);
double* viscousDissipation_v(double ***v, MAC_grid *grid, double iRe);
double* viscousDissipation_w(double ***w, MAC_grid *grid, double iRe);
double pressureWork_u(double ***u, double ***p, MAC_grid *grid);
double pressureWork_v(double ***v, double ***p, MAC_grid *grid);
double pressureWork_w(double ***w, double ***p, MAC_grid *grid);
double totalMass(double ***c, MAC_grid *grid);
double mixingRate(double ***c0, double ***c1, MAC_grid *grid);
double* horizontalMean(double ***quantity, MAC_grid *grid);
double* turbulentStress(double ***quantity1, double ***quantity2, MAC_grid *grid, int rank);
double* maximalVelocity(double ***u_cell, double ***v_cell, double ***w_cell, MAC_grid *grid);
double** min_max_concentration(double ***conc, MAC_grid *grid);




#endif

