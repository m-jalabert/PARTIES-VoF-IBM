#ifndef STATISTIC2D_H
 #define STATISTICS2D_H

//Threshold for front location
#define THRESHOLD_FRONT_LOCATION  0.01

Statistics2d *Statistics2d_create(MAC_grid *grid, Parameters *params);
void Statistics2d_destroy(Statistics2d *st2d, MAC_grid *grid, Parameters *params);
void Statistics2d_computeStatistics(Cart3d_bag *data_bag, Debug_trace *dtrace);
double Statistics2d_totKineticEnergy(double **kinE, double ***u_cell,double ***v_cell, double ***w_cell, Cart3d_bag *data_bag);
double Statistics2d_RMS(double **Urms, double ***u_cell, double ***v_cell, double ***w_cell, Cart3d_bag *data_bag);
double Statistics2d_forcing(double ***u_face, double ***v_face, double ***w_face, double ***u_forc, double ***v_forc, double ***w_forc, Cart3d_bag *data_bag);
double Statistics2d_kinE_part(Cart3d_bag *data_bag, Debug_trace *dtrace);
double Statistics2d_potE_part(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Statistics2d_potEnergy_y(double **potE, double ***c, double Ri ,double *integral_potEnergy_fluid,Cart3d_bag *data_bag);
void Statistics2d_partial_y(double **flux, double ***quantity1, Cart3d_bag *data_bag);
void Statistics2d_potEnergyPart(double **potEpart, Cart3d_bag *data_bag);

double Statistics2d_volume_integral(double ***c,   int domain, Cart3d_bag *data_bag);
void Statistics2d_horizontalMean2d(double ***quantity, double **global_mean, char which_component, Cart3d_bag *data_bag);
void Statistics2d_horizontalSlice2d(double ***quantity, double **global_slice, char which_component, Cart3d_bag *data_bag);
void Statistics2d_turbulentStress(double **turbStress, double ***quantity1,
	           double ***quantity2, char which_component, Cart3d_bag *data_bag);
void Statistics2d_fluxes(double **flux, double ***quantity1, double ***quantity2, char which_component, Cart3d_bag *data_bag);
void Statistics2d_globalmean2d(double **mean, double **global_mean, MAC_grid *grid);
void Statistics2d_pressureWork(double ***u, double ***v, double ***w, Pressure *p, double **pwork, Cart3d_bag *data_bag);
void Statistics2d_uface_to_cellcenter(double **face, double **cell, MAC_grid *grid);
void Statistics2d_vface_to_cellcenter(double **face, double **cell, MAC_grid *grid);
void Statistics2d_add_array(double **work, double **work1, MAC_grid *grid);
void Statistics2d_covarcalc(double *avg_velc, double ***c_data, double ***vel_data, double *avg_mvf, double *avg_c, double *avg_vel, Cart3d_bag *data_bag);
void Statistics2d_histocalc(double ***c_data, double *histo, double *range, Cart3d_bag *data_bag);

void Statistics2d_ke_budget(double ***u_data, double ***v_data, double ***w_data,
		double ***u_data_bc, double ***v_data_bc, double ***w_data_bc,
		double ***nut, Wall_model *log_law, double **global_viscDiss, double **global_viscTran,
		double **global_sgsDiss, double **global_sgsTran, Cart3d_bag *data_bag, double iRe);

void Statistics2d_ke_budget_u(double ***u_data, double ***v_data, double ***w_data,
		double ***w_data_bc, double ***nut, Wall_model *log_law,
		double **viscDiss, double **viscTran, double **sgsDiss,
		double **sgsTran, Cart3d_bag *data_bag, double iRe);

void Statistics2d_ke_budget_v(double ***u_data, double ***v_data, double ***w_data,
		double ***w_data_bc, double ***nut, Wall_model *log_law,
		double **viscDiss, double **viscTran, double **sgsDiss,
		double **sgsTran, Cart3d_bag *data_bag, double iRe);

void Statistics2d_ke_budget_w(double ***u_data, double ***v_data, double ***w_data,
		double ***w_data_bc, double ***nut, Wall_model *log_law,
		double **viscDiss, double **viscTran, double **sgsDiss,
		double **sgsTran, Cart3d_bag *data_bag, double iRe);

void Statistics2d_vf_avg(double **vf_avg, Cart3d_bag *data_bag);

void Statistics2d_writeStatistics2dToH5(Statistics2d *st2d, MAC_grid *grid,
		Parameters *params, Cart3d_bag *data_bag ,int timestep, Debug_trace *dtrace);
//void Statistics1d_MeanYZ(double ***quantity,  double *mean_1d_full, double *mean_1d_fluid , double factor,  Cart3d_bag *data_bag);
void Statistics2d_compute_nusselt( Cart3d_bag *data_bag);
//void Statistics0d_Frontposition (double ***conc, double *current_height, double *current_height_fluid , double *front_position  ,  Cart3d_bag *data_bag);
void Statistics2d_fill_ghostnodes_2d(double **quantity, MAC_grid *grid) ;
//void Statistics2d_horizontalSlice2d_near_Wall(double ***quantity, double **global_slice, char which_component, Cart3d_bag *data_bag);

#endif
