void Post_processing_conc(Cart3d_bag *data_bag);
void Conc_compute_total_suspended_mass(Concentration *c, MAC_grid *grid) ;
void Conc_compute_ave_height_x(Concentration **c, MAC_grid *grid, Parameters *params) ;
void Conc_find_front_location(Concentration **c, MAC_grid *grid, Parameters *params, double front_limit) ;
void Conc_integrate_deposited_height(Concentration *c, MAC_grid *grid, Parameters *params, double weight_factor) ;
void Conc_update_world_deposited_height(Concentration *c, Cart3d_bag *data_bag);

void Conc_update_world_deposited_height_dumped(Concentration **c, Cart3d_bag *data_bag);
void Conc_dump_particles(Concentration **c, Cart3d_bag *data_bag);
void Conc_compute_stokes_dissipation_rate(int iconc, Concentration **c, MAC_grid *grid, Parameters *params);
void Conc_update_world_stokes_dissipation_rate( Concentration **c, MAC_grid *grid, Parameters *params);
void Conc_compute_active_potential_energy(Concentration *c, MAC_grid *grid, Parameters *params);
void Conc_compute_passive_potential_energy(Concentration *c, MAC_grid *grid, Parameters *params) ;
void Conc_update_world_potential_energies(Concentration **c, MAC_grid *grid, Parameters *params) ;
