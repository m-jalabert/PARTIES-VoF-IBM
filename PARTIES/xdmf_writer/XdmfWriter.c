#include <stdlib.h>
#include <stdio.h>
#include "DataTypes.h"
#include "Input.h"
#include "Marker.h"
#include "Particle.h"

#define CENTER_FILE "Reader_c.xmf"
#define U_FACE_FILE "Reader_u.xmf"
#define V_FACE_FILE "Reader_v.xmf"
#define W_FACE_FILE "Reader_w.xmf"
#define NODE_FILE "Reader_node.xmf"

void initializeData3d(Data3d **data3d_ptr, Data3d **data3d_last_ptr, char *dataname, char nodeType);
void startFile3d(Data3d *data3d_start, char nodeType, Parameters *params, Hyperslab *hyperslab);
void startFileOther(OtherData *otherData, Parameters *params);
void endFile3d(Data3d *data3d_start, char nodeType);
void endFileOther(OtherData *otherData);
void writeGeometry(FILE *xmf, char nodeType, char *indent, int iter, Hyperslab *hyperslab);
void writeGeometryOther(FILE *xmf, OtherData *otherData, char *indent, int iter);
void writeTimestep3d(Data3d *data3d_start, char nodeType, int index,
		double time, Hyperslab *hyperslab, Parameters *params);
void writeTimestepOther(OtherData *otherData, int index, double time, Parameters *params);
void writeAttribute3d(FILE *xmf, Data3d *data3d, int iter, Hyperslab *hyperslab);
int existsNodeType(Data3d *data3d, char nodeType);


/******************************************************************************/
/*
 Writes .xmf files in XDMF format
 (see http://www.xdmf.org/index.php/XDMF_Model_and_Format)

 Default values and explanations for various parameters are found in default.inp
 Structure declarations are in DataTypes.h

 This progam writes up to four files: one each for u-faced, v-faced, w-faced,
 and cell-centered values.

 -------------------------------------------------------------------------------
 To add a new variable to print:
     - add a print flag to the Parameters structure in Datatypes.h
     - add the default value to default.inp
     - add an 'if' statement below to include it in the linked list
 -------------------------------------------------------------------------------
 */
/******************************************************************************/
int main(int argc, char** argv) {

	int i;
	double time;
	char dataname[50];



	//--------------------------------------------------------------------------
	// Read input
	//--------------------------------------------------------------------------
	Parameters *params = (Parameters *)malloc(sizeof(Parameters));
	Hyperslab *hyperslab = (Hyperslab *)malloc(sizeof(Hyperslab));
/*
	// if you pass 2d as argument only the other data is written
	if (argc == 2){
		if ( strncmp( argv[1], "2d", 2 ) == 0){
		params->flag = 1;
	}else{
		params->flag=0;
	}}else{
		params->flag=0;
	}
*/
	Input_get_values(params, hyperslab);


	// Concentration parameters
	int iconc;
	int NConc = params -> NConc;

	// Display parameters
	#define CFG(s, n, default, reader) \
		printf("%s->%s = %d\n", #s, #n, s->n);
	#define CFG_STRING(s, n, default) \
		printf("%s->%s = %s\n", #s, #n, s->n);
	#include "default.inp"
	#undef CFG
	printf("NX = %d\n", hyperslab->NX);
	printf("NY = %d\n", hyperslab->NY);
	printf("NZ = %d\n", hyperslab->NZ);
	for (i = 0; i < params->iter_count; i++) {
		printf("\tt_%d = %f\n", FILE_INDEX(i), params->timesteps[i]);
	}

	// Initialize linked list variables
	Data3d *data3d_start = (Data3d *)malloc(sizeof(Data3d));
	Data3d *data3d, *data3d_last;
	data3d_last = NULL;
	data3d = data3d_start;
	DataParticle *p_data, *p_data_last;
	OtherData *other, *other_last;


	//--------------------------------------------------------------------------
	// Add parameters to linked list
	//--------------------------------------------------------------------------

	// u-faced u velocity
	if (params -> print_u) {
		initializeData3d(&data3d, &data3d_last, "u", 'u');
	}

	// v-faced v velocity
	if (params -> print_v) {
		initializeData3d(&data3d, &data3d_last, "v", 'v');
	}

	// w-faced w velocity
	if (params -> print_w) {
		initializeData3d(&data3d, &data3d_last, "w", 'w');
	}

	// Cell-centered u velocity
	if (params -> print_uc) {
		initializeData3d(&data3d, &data3d_last, "uc", 'c');
	}

	// Cell-centered v velocity
	if (params -> print_vc) {
		initializeData3d(&data3d, &data3d_last, "vc", 'c');
	}

	// Cell-centered w velocity
	if (params -> print_wc) {
		initializeData3d(&data3d, &data3d_last, "wc", 'c');
	}

	// Pressure
	if (params -> print_p) {
		initializeData3d(&data3d, &data3d_last, "p", 'c');
	}

	// Concentration
	if (params -> print_conc) {
		for (iconc = 0; iconc < NConc; iconc++) {
			sprintf(dataname, "Conc/%d", iconc);
			initializeData3d(&data3d, &data3d_last, dataname, 'c');
		}
	}

	// Particle volume fraction (for fully-resolved particles)
	if (params -> print_vf) {
		initializeData3d(&data3d, &data3d_last, "vfu", 'u');
		initializeData3d(&data3d, &data3d_last, "vfv", 'v');
		initializeData3d(&data3d, &data3d_last, "vfw", 'w');
		//initializeData3d(&data3d, &data3d_last, "u_vof", 'u');
		//initializeData3d(&data3d, &data3d_last, "v_vof", 'v');
		//initializeData3d(&data3d, &data3d_last, "w_vof", 'w');
		//initializeData3d(&data3d, &data3d_last, "vfc", 'c');
	}

	// SGS eddy viscosity
	if (params -> print_nut) {
		initializeData3d(&data3d, &data3d_last, "LES/nut", 'c');
	}

	// SGS eddy diffusivity
	if (params -> print_conc_alphat) {
		for (iconc = 0; iconc < NConc; iconc++) {
			sprintf(dataname, "LES/Conc%d_alphat", iconc);
			initializeData3d(&data3d, &data3d_last, dataname, 'c');
		}
	}

	// numerator of dynamic model coefficient
	if (params -> print_ilm) {
		initializeData3d(&data3d, &data3d_last, "LES/ilm", 'c');
	}

	// denominator of dynamic model coefficient
	if (params -> print_imm) {
		initializeData3d(&data3d, &data3d_last, "LES/imm", 'c');
	}

	// numerator of dynamic model coefficient
	if (params -> print_conc_itt) {
		for (iconc = 0; iconc < NConc; iconc++) {
			sprintf(dataname, "LES/Conc%d_itt", iconc);
			initializeData3d(&data3d, &data3d_last, dataname, 'c');
		}
	}

	// denominator of dynamic model coefficient
	if (params -> print_conc_ikt) {
		for (iconc = 0; iconc < NConc; iconc++) {
			sprintf(dataname, "LES/Conc%d_ikt", iconc);
			initializeData3d(&data3d, &data3d_last, dataname, 'c');
		}
	}

	// Check if we are going to create any .xmf files
	if (data3d_last == NULL) {
        if (params->mobile_data == NULL && params->fixed_data == NULL &&
			params->marker_data == NULL && params->otherData == NULL) {
    		printf("Not creating any .xmf files!\n");
    		return EXIT_SUCCESS;
        }
        else {
            free(data3d);
            data3d_start = NULL;
        }
	}
	else {
		free(data3d);
		data3d_last -> next = NULL;
	}


	//--------------------------------------------------------------------------
	// Write XDMF file
	//--------------------------------------------------------------------------

	// Write first part of .xmf file
	startFile3d(data3d_start, 'u', params, hyperslab);
	startFile3d(data3d_start, 'v', params, hyperslab);
	startFile3d(data3d_start, 'w', params, hyperslab);
	startFile3d(data3d_start, 'c', params, hyperslab);
	Particle_startFile(params->mobile_data, 'm', params);
	Particle_startFile(params->fixed_data,  'f', params);
	Marker_startFile(params->marker_data, 0, params);
	Marker_startFile(params->marker_data, 1, params);

	other = params->otherData;
	while (other != NULL) {
		startFileOther(other, params);
		other = other -> next;
	}

	// Write each timestep
	for (i = 0; i < params->iter_count; i++) {

		time = params->timesteps[i];
		//printf("time=%2.3f\n",time);

		writeTimestep3d(data3d_start, 'u', i, time, hyperslab, params);
		writeTimestep3d(data3d_start, 'v', i, time, hyperslab, params);
		writeTimestep3d(data3d_start, 'w', i, time, hyperslab, params);
		writeTimestep3d(data3d_start, 'c', i, time, hyperslab, params);
		Particle_writeTimestep(params->mobile_data, 'm', i, time, params->N_p_mobile, params);
		Particle_writeTimestep(params->fixed_data,  'f', i, time, params->N_p_fixed, params);
		Marker_writeTimestep(params->marker_data, 0, i, time, params);
		Marker_writeTimestep(params->marker_data, 1, i, time, params);

		other = params->otherData;
		while (other != NULL) {
			writeTimestepOther(other, i, time, params);
			other = other -> next;
		}
	}

	// Write end of .xmf file
	endFile3d(data3d_start, 'u');
	endFile3d(data3d_start, 'v');
	endFile3d(data3d_start, 'w');
	endFile3d(data3d_start, 'c');
	Particle_endFile(params->mobile_data, 'm', params);
	Particle_endFile(params->fixed_data,  'f', params);
	Marker_endFile(params->marker_data, 0, params);
	Marker_endFile(params->marker_data, 1, params);

	other = params->otherData;
	while (other != NULL) {
		endFileOther(other);
		other = other -> next;
	}


	//--------------------------------------------------------------------------
	// Free memory
	//--------------------------------------------------------------------------
	data3d = data3d_start;
	while (data3d != NULL) {

		data3d_last = data3d;
		data3d = data3d -> next;
		free(data3d_last);
	}
	p_data = params -> mobile_data;
	while (p_data != NULL) {
		p_data_last = p_data;
		p_data = p_data -> next;
		free(p_data_last);
	}
	p_data = params -> fixed_data;
	while (p_data != NULL) {
		p_data_last = p_data;
		p_data = p_data -> next;
		free(p_data_last);
	}
	p_data = params -> marker_data;
	while (p_data != NULL) {
		p_data_last = p_data;
		p_data = p_data -> next;
		free(p_data_last);
	}
	free(params->timesteps);
	if (params->mobile_data != NULL) {
		free(params->N_p_mobile);
	}
	if (params->fixed_data != NULL) {
		free(params->N_p_fixed);
	}
	if (params->marker_data != NULL) {
		free(params->N_pt_marker);
		free(params->N_poly_marker);
	}
	other = params->otherData;
	while (other != NULL) {
		other_last = other;
		other = other -> next;
		free(other_last);
	}
	free(params);
	free(hyperslab);
	return EXIT_SUCCESS;
}




/******************************************************************************/
/*
 Creates and adds a Data3d structure to the linked list

     - dataname: name of 3-D variable as it appears in the HDF5 file
     - nodeType: specify grid position of variable - 'u', 'v', 'w', 'c', or 'n'
                 for u-faced, v-faced, w-faced, cell-centered, or node-centered
                 (values at cell corners) variable
 */
/******************************************************************************/
void initializeData3d(Data3d **data3d_ptr, Data3d **data3d_last_ptr,
		char *dataname, char nodeType) {

	// Load values into linked list object
	sprintf((*data3d_ptr) -> dataname, "%s", dataname);
	(*data3d_ptr) -> nodeType = nodeType;

	// Advance linked list
	(*data3d_ptr) -> next = (Data3d *)malloc(sizeof(Data3d));
	*data3d_last_ptr = *data3d_ptr;
	*data3d_ptr = (*data3d_ptr) -> next;

}




/******************************************************************************/
/*
 Writes the first part of the file up until the timesteps
 */
/******************************************************************************/
void startFile3d(Data3d *data3d_start, char nodeType, Parameters *params,
		Hyperslab *hyperslab) {

	FILE *xmf;

	// See if we are writing any data for this node type. If not, exit.
	if (!existsNodeType(data3d_start, nodeType))
		return;

	// Write new file
	if (nodeType == 'u')
		xmf = fopen(U_FACE_FILE, "w");
	else if (nodeType == 'v')
		xmf = fopen(V_FACE_FILE, "w");
	else if (nodeType == 'w')
		xmf = fopen(W_FACE_FILE, "w");
	else if (nodeType == 'c')
		xmf = fopen(CENTER_FILE, "w");
	else if (nodeType == 'n')
		xmf = fopen(NODE_FILE, "w");
	else {
		fprintf(stderr, "Error: Invalid input \"nodeType\" = '%c'\n", nodeType);
		fprintf(stderr, "       in function startFile3d(). Use 'u', 'v', 'w', or 'c'.\n\n");
		exit(EXIT_FAILURE);
	}

	// Opening lines of .xmf file
	fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
	fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
	fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
	fprintf(xmf, "  <Domain>\n");
	fprintf(xmf, "\n");

	// Grid geometry
	writeGeometry(xmf, nodeType, "    ", params->iter_start, hyperslab);

	// Opening of temporal grid
	fprintf(xmf, "\n");
	fprintf(xmf, "    <Grid Name=\"TemporalGrid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
	fprintf(xmf, "\n");

	fclose(xmf);
}




/******************************************************************************/
/*
 Writes the first part of the file up until the timesteps
 */
/******************************************************************************/
void startFileOther(OtherData *otherData, Parameters *params) {

	FILE *xmf;
	char filename[50];

	sprintf(filename, "Reader_%s.xmf", otherData->name);
	xmf = fopen(filename, "w");

	// Opening lines of .xmf file
	fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
	fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
	fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
	fprintf(xmf, "  <Domain>\n");
	fprintf(xmf, "\n");

	// Grid geometry
	writeGeometryOther(xmf, otherData, "    ", params->iter_start);

	// Opening of temporal grid
	fprintf(xmf, "\n");
	fprintf(xmf, "    <Grid Name=\"TemporalGrid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
	fprintf(xmf, "\n");

	fclose(xmf);
}




/******************************************************************************/
/*
 Write last few lines to .xmf file and display which files were created
 */
/******************************************************************************/
void endFile3d(Data3d *data3d_start, char nodeType) {

	FILE *xmf;
	char filename[50];

	// See if we are writing any data for this node type. If not, exit.
	if (!existsNodeType(data3d_start, nodeType))
		return;

	// Select which file
	if (nodeType == 'u')
		sprintf(filename, U_FACE_FILE);
	else if (nodeType == 'v')
		sprintf(filename, V_FACE_FILE);
	else if (nodeType == 'w')
		sprintf(filename, W_FACE_FILE);
	else if (nodeType == 'c')
		sprintf(filename, CENTER_FILE);
	else if (nodeType == 'n')
		sprintf(filename, NODE_FILE);
	else {
		fprintf(stderr, "Error: Invalid input \"nodeType\" = '%c'\n", nodeType);
		fprintf(stderr, "       in function endFile3d(). Use 'u', 'v', 'w', or 'c'.\n\n");
		exit(EXIT_FAILURE);
	}

	// Append closing lines
	xmf = fopen(filename, "a");
	fprintf(xmf, "    </Grid>\n");
	fprintf(xmf, "  </Domain>\n");
	fprintf(xmf, "</Xdmf>\n");
	fclose(xmf);

	// Display which file was created
	printf("Wrote %s\n", filename);
}




/******************************************************************************/
/*
 Write last few lines to .xmf file and display which files were created
 */
/******************************************************************************/
void endFileOther(OtherData *otherData) {

	FILE *xmf;
	char filename[50];

	// Append closing lines
	sprintf(filename, "Reader_%s.xmf", otherData->name);
	xmf = fopen(filename, "a");
	fprintf(xmf, "    </Grid>\n");
	fprintf(xmf, "  </Domain>\n");
	fprintf(xmf, "</Xdmf>\n");
	fclose(xmf);

	// Display which file was created
	printf("Wrote %s\n", filename);
}




/******************************************************************************/
/*
 Write specified geometry to file.
     - indent: string of spaces to properly indent geometry writing
     - iter: grid data will be read in from Data_'iter'.h
 */
/******************************************************************************/
void writeGeometry(FILE *xmf, char nodeType, char *indent, int iter, Hyperslab *hyperslab) {

	char xfile[50], yfile[50], zfile[50];

	int NX = hyperslab -> NX;
	int NY = hyperslab -> NY;
	int NZ = hyperslab -> NZ;

	int Is_HS = hyperslab -> Is_HS;
	int Js_HS = hyperslab -> Js_HS;
	int Ks_HS = hyperslab -> Ks_HS;

	int dI_HS = hyperslab -> dI_HS;
	int dJ_HS = hyperslab -> dJ_HS;
	int dK_HS = hyperslab -> dK_HS;

	int NI_HS = hyperslab -> NI_HS;
	int NJ_HS = hyperslab -> NJ_HS;
	int NK_HS = hyperslab -> NK_HS;

	// This takes care of cell-centered data
	sprintf(xfile, "Data_%d.h5:/grid/xc", iter);
	sprintf(yfile, "Data_%d.h5:/grid/yc", iter);
	sprintf(zfile, "Data_%d.h5:/grid/zc", iter);

	if (nodeType == 'u' || nodeType == 'n') {
		NI_HS += 1;
		sprintf(xfile, "Data_%d.h5:/grid/xu", iter);
	}
	if (nodeType == 'v' || nodeType == 'n') {
		NJ_HS += 1;
		sprintf(yfile, "Data_%d.h5:/grid/yv", iter);
	}
	if (nodeType == 'w' || nodeType == 'n') {
		NK_HS += 1;
		sprintf(zfile, "Data_%d.h5:/grid/zw", iter);
	}

	fprintf(xmf, "%s<Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n", indent, NK_HS, NJ_HS, NI_HS);
	fprintf(xmf, "%s<Geometry GeometryType=\"VXVYVZ\">\n", indent);
	fprintf(xmf, "%s  <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d\" Type=\"HyperSlab\">\n", indent, NI_HS);
	fprintf(xmf, "%s    <DataItem Dimensions=\"3\" Format=\"XML\">\n", indent);
	fprintf(xmf, "%s      %-4d %-4d %-4d\n", indent, Is_HS, dI_HS, NI_HS);
	fprintf(xmf, "%s    </DataItem>\n", indent);
	fprintf(xmf, "%s    <DataItem Format=\"HDF\" Dimensions=\"%d\">\n", indent, NX);
	fprintf(xmf, "%s      %s\n", indent, xfile);
	fprintf(xmf, "%s    </DataItem>\n", indent);
	fprintf(xmf, "%s  </DataItem>\n", indent);
	fprintf(xmf, "%s  \n", indent);
	fprintf(xmf, "%s  <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d\" Type=\"HyperSlab\">\n", indent, NJ_HS);
	fprintf(xmf, "%s    <DataItem Dimensions=\"3\" Format=\"XML\">\n", indent);
	fprintf(xmf, "%s      %-4d %-4d %-4d\n", indent, Js_HS, dJ_HS, NJ_HS);
	fprintf(xmf, "%s    </DataItem>\n", indent);
	fprintf(xmf, "%s    <DataItem Format=\"HDF\" Dimensions=\"%d\">\n", indent, NY);
	fprintf(xmf, "%s      %s\n", indent, yfile);
	fprintf(xmf, "%s    </DataItem>\n", indent);
	fprintf(xmf, "%s  </DataItem>\n", indent);
	fprintf(xmf, "%s  \n", indent);
	fprintf(xmf, "%s  <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d\" Type=\"HyperSlab\">\n", indent, NK_HS);
	fprintf(xmf, "%s    <DataItem Dimensions=\"3\" Format=\"XML\">\n", indent);
	fprintf(xmf, "%s      %-4d %-4d %-4d\n", indent, Ks_HS, dK_HS, NK_HS);
	fprintf(xmf, "%s    </DataItem>\n", indent);
	fprintf(xmf, "%s    <DataItem Format=\"HDF\" Dimensions=\"%d\">\n", indent, NZ);
	fprintf(xmf, "%s      %s\n", indent, zfile);
	fprintf(xmf, "%s    </DataItem>\n", indent);
	fprintf(xmf, "%s  </DataItem>\n", indent);
	fprintf(xmf, "%s</Geometry>\n", indent);
}




/******************************************************************************/
/*
 Write specified geometry to file.
     - indent: string of spaces to properly indent geometry writing
     - iter: grid data will be read in from Data_'iter'.h
 */
/******************************************************************************/
void writeGeometryOther(FILE *xmf, OtherData *otherData, char *indent, int iter) {

	char xfile[50], yfile[50], zfile[50];

	// This takes care of cell-centered data
	sprintf(xfile, "%s_%d.h5:/grid/x", otherData->name, iter);
	sprintf(yfile, "%s_%d.h5:/grid/y", otherData->name, iter);
	sprintf(zfile, "%s_%d.h5:/grid/z", otherData->name, iter);

	fprintf(xmf, "%s<Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n", indent, otherData->NZ, otherData->NY, otherData->NX);
	fprintf(xmf, "%s<Geometry GeometryType=\"VXVYVZ\">\n", indent);
	fprintf(xmf, "%s  <DataItem Format=\"HDF\" Dimensions=\"%d\">\n", indent, otherData->NX);
	fprintf(xmf, "%s    %s\n", indent, xfile);
	fprintf(xmf, "%s  </DataItem>\n", indent);
	fprintf(xmf, "%s  \n", indent);
	fprintf(xmf, "%s  <DataItem Format=\"HDF\" Dimensions=\"%d\">\n", indent, otherData->NY);
	fprintf(xmf, "%s    %s\n", indent, yfile);
	fprintf(xmf, "%s  </DataItem>\n", indent);
	fprintf(xmf, "%s  \n", indent);
	fprintf(xmf, "%s  <DataItem Format=\"HDF\" Dimensions=\"%d\">\n", indent, otherData->NZ);
	fprintf(xmf, "%s    %s\n", indent, zfile);
	fprintf(xmf, "%s  </DataItem>\n", indent);
	fprintf(xmf, "%s</Geometry>\n", indent);
}




/******************************************************************************/
/*
 Writes data for entire linked list data3d_start for a single timestep
 */
/******************************************************************************/
void writeTimestep3d(Data3d *data3d_start, char nodeType, int index,
		double time, Hyperslab *hyperslab, Parameters *params) {

	Data3d *data3d;
	FILE *xmf;

	// See if we are writing any data for this node type. If not, exit.
	if (!existsNodeType(data3d_start, nodeType))
		return;

	// Append to file
	if (nodeType == 'u')
		xmf = fopen(U_FACE_FILE, "a");
	else if (nodeType == 'v')
		xmf = fopen(V_FACE_FILE, "a");
	else if (nodeType == 'w')
		xmf = fopen(W_FACE_FILE, "a");
	else if (nodeType == 'c')
		xmf = fopen(CENTER_FILE, "a");
	else if (nodeType == 'n')
		xmf = fopen(NODE_FILE, "a");
	else {
		fprintf(stderr, "Error: Invalid input \"nodeType\" = '%c'\n", nodeType);
		fprintf(stderr, "       in function writeTimestep(). Use 'u', 'v', 'w', or 'c'.\n\n");
		exit(EXIT_FAILURE);
	}

	// Describe time value and grid
	fprintf(xmf, "      <Grid Name=\"SpatialGrid_%d\" GridType=\"Uniform\">\n", index);
	fprintf(xmf, "        <Time Value=\"%f\"/>\n", time);
	fprintf(xmf, "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n");
	fprintf(xmf, "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");

	// Write attributes for any data that matches this node type
	data3d = data3d_start;
	while (data3d != NULL) {

		if (data3d -> nodeType == nodeType)
			writeAttribute3d(xmf, data3d, FILE_INDEX(index), hyperslab);

		data3d = data3d -> next;
	}

	// Close up timestep
	fprintf(xmf, "      </Grid>\n");
	fprintf(xmf, "\n");

	fclose(xmf);
}




/******************************************************************************/
/*
 Writes data for entire linked list data3d_start for a single timestep
 */
/******************************************************************************/
void writeTimestepOther(OtherData *otherData, int index, double time, Parameters *params) {

	FILE *xmf;
	char filename[50];

	OtherDataElement *data_element = otherData->data_element;

	sprintf(filename, "Reader_%s.xmf", otherData->name);
	xmf = fopen(filename, "a");

	// Describe time value and grid
	fprintf(xmf, "      <Grid Name=\"SpatialGrid_%d\" GridType=\"Uniform\">\n", index);
	fprintf(xmf, "        <Time Value=\"%f\"/>\n", time);
	fprintf(xmf, "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n");
	fprintf(xmf, "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");
	while (data_element != NULL) {
		fprintf(xmf, "        <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", data_element->dataname);
		fprintf(xmf, "          <DataItem Format=\"HDF\" NumberType=\"Double\" Precision=\"8\" ");
		fprintf(xmf, "Dimensions=\"%d %d %d\">\n", otherData->NX, otherData->NY, otherData->NZ);
		fprintf(xmf, "            %s_%d.h5:/%s\n", otherData->name, FILE_INDEX(index), data_element->dataname);
		fprintf(xmf, "          </DataItem>\n");
		fprintf(xmf, "        </Attribute>\n");
		data_element = data_element -> next;
	}
	fprintf(xmf, "      </Grid>\n");
	fprintf(xmf, "\n");

	fclose(xmf);
}




/******************************************************************************/
/*
 Writes a single attribute to the file, that is, a single object in the linked
 list.
 */
/******************************************************************************/
void writeAttribute3d(FILE *xmf, Data3d *data3d, int iter, Hyperslab *hyperslab) {

	char *dataname = data3d -> dataname;
	char nodeType  = data3d -> nodeType;

	int NX = hyperslab -> NX;
	int NY = hyperslab -> NY;
	int NZ = hyperslab -> NZ;

	int Is_HS = hyperslab -> Is_HS;
	int Js_HS = hyperslab -> Js_HS;
	int Ks_HS = hyperslab -> Ks_HS;

	int dI_HS = hyperslab -> dI_HS;
	int dJ_HS = hyperslab -> dJ_HS;
	int dK_HS = hyperslab -> dK_HS;

	int NI_HS = hyperslab -> NI_HS;
	int NJ_HS = hyperslab -> NJ_HS;
	int NK_HS = hyperslab -> NK_HS;

	if (nodeType == 'u' || nodeType == 'n')
		NI_HS += 1;
	if (nodeType == 'v' || nodeType == 'n')
		NJ_HS += 1;
	if (nodeType == 'w' || nodeType == 'n')
		NK_HS += 1;

	fprintf(xmf, "        <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", dataname);
	fprintf(xmf, "          <DataItem ItemType=\"HyperSlab\" ");
	fprintf(xmf, "Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", NK_HS, NJ_HS, NI_HS);
	fprintf(xmf, "            <DataItem Dimensions=\"3 3\" Format=\"XML\">\n");
	fprintf(xmf, "              %-4d %-4d %-4d\n", Ks_HS, Js_HS, Is_HS);
	fprintf(xmf, "              %-4d %-4d %-4d\n", dK_HS, dJ_HS, dI_HS);
	fprintf(xmf, "              %-4d %-4d %-4d\n", NK_HS, NJ_HS, NI_HS);
	fprintf(xmf, "            </DataItem>\n");
	fprintf(xmf, "            <DataItem Format=\"HDF\" NumberType=\"Double\" Precision=\"8\" ");
	fprintf(xmf, "Dimensions=\"%d %d %d\">\n", NX, NY, NZ);
	fprintf(xmf, "              Data_%d.h5:/%s\n", iter, dataname);
	fprintf(xmf, "            </DataItem>\n");
	fprintf(xmf, "          </DataItem>\n");
	fprintf(xmf, "        </Attribute>\n");

}




/******************************************************************************/
/*
 Searches Data3d linked list to see if a particular nodeType exists.  This is
 useful to see if we need to create or write a u-faced grid file, for instance.
 */
/******************************************************************************/
int existsNodeType(Data3d *data3d, char nodeType) {

	while (data3d != NULL) {
		if (data3d -> nodeType == nodeType)
			return 1;
		data3d = data3d -> next;
	}

	return 0;
}
