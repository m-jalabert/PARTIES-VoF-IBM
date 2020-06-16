#ifndef DATATYPES_H
#define DATATYPES_H

#define FILE_INDEX(i) (params->iter_start + (i) * params->iter_step)

struct hyperslab{

	// Dimensions of 3-D data (equivalent to those in HDF5 file)
	int NX, NY, NZ;

	// Hyperslab starting indices
	int Is_HS, Js_HS, Ks_HS;

	// Hyperslab stride
	int dI_HS, dJ_HS, dK_HS;

	// Hyperslab dimensions
	int NI_HS, NJ_HS, NK_HS;
};
typedef struct hyperslab Hyperslab;


struct otherData{

	// Filename for other type of data, e.g. "MyFile2D" for MyFile2D_*.h5 files
	char name[50];

	// Dimensions of 'otherData' dataset
	int NX, NY, NZ;


	struct otherDataElement *data_element;

	// Pointer used for linked list functionality
	struct otherData *next;
};
typedef struct otherData OtherData;


struct otherDataElement{

	// Dataname contained within file for 'otherData'
	char dataname[50];

	// Pointer used for linked list functionality
	struct otherDataElement *next;
};
typedef struct otherDataElement OtherDataElement;


struct data3d{
	// 'u', 'v', 'w', or 'c' to denote what grid variable is on
	char nodeType;

	// Name of variable in HDF5 files
	// 3-D data should be in Data_*.h5:/dataname/data
	//     (e.g. Data_0.h5:/u/data)
	char dataname[50];

	// Pointer used for linked list functionality
	struct data3d *next;
};
typedef struct data3d Data3d;


struct dataParticle{
	char dataname[50];
	int length;
	struct dataParticle *next;
};
typedef struct dataParticle DataParticle;


struct parameters{

	// Starting and ending indices of HDF5 output files
	int iter_start, iter_step, iter_end, iter_count;

	// Array of timesteps
	double *timesteps;

	// Array of number of particles for each timestep
	int *N_p_mobile;
	int *N_p_fixed;

	// Array of number of marker points for each timestep
	int *N_pt_marker;
	// Array of number of polygons for each timestep
	int *N_poly_marker;
	// Number of points per marker polygon
	int size_poly_marker;

	int print_u, print_v, print_w;
	int print_uc, print_vc, print_wc;
	int print_p, print_conc, print_vf;
	int print_nut, print_conc_alphat;
	int print_ilm, print_imm, print_conc_itt, print_conc_ikt;
	int NConc;

	char otherData_name[50];
	OtherData *otherData;

	DataParticle *mobile_data, *fixed_data, *marker_data;
};
typedef struct parameters Parameters;

#endif
