	#ifndef EPFORCING_H
		#define EPFORCING_H


		#include "DataTypes.h"

	int ** forced_freq(int *alloc, double Kf);
	void init_turb_forcing(Fourier *fourier, Parameters *params);
	void b_and_fourier_update(Fourier *fourier, Parameters *params);
	void ftx_update(Cart3d_bag *data_bag);
	void fty_update(Cart3d_bag *data_bag);
	void ftz_update(Cart3d_bag *data_bag);
	double box_muller(float m, float s);


#endif
