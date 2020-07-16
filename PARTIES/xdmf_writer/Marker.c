#include <stdlib.h>
#include <stdio.h>
#include "DataTypes.h"
#include "Marker.h"

#define MARKER_FILE_0 "Reader_marker_black.xmf"
#define MARKER_FILE_1 "Reader_marker_red.xmf"

/******************************************************************************/
/*
 Writes data for entire linked list data_start for a single timestep
 */
/******************************************************************************/
void Marker_startFile(DataParticle *p_data_start, int mType, Parameters *params) {

	FILE *xmf;

	if (p_data_start == NULL)
		return;

	if (mType == 0)
		xmf = fopen(MARKER_FILE_0, "w");
	else if (mType == 1)
		xmf = fopen(MARKER_FILE_1, "w");
	else {
		fprintf(stderr, "Error: Invalid input \"mType\" = %d\n", mType);
		fprintf(stderr, "       in function Marker_startFile(). Use 0, or 1.\n\n");
		exit(EXIT_FAILURE);
	}

	fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
	fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
	fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
	fprintf(xmf, "  <Domain>\n");
	fprintf(xmf, "\n");
	fprintf(xmf, "    <Grid Name=\"TemporalGrid\" GridType=\"Collection\" CollectionType=\"Temporal\">");
	fprintf(xmf, "\n");
	fclose(xmf);
}




/******************************************************************************/
/*
 Writes data for entire linked list data_start for a single timestep
 */
/******************************************************************************/
void Marker_endFile(DataParticle *p_data_start, int mType, Parameters *params) {

	char filename[50];

	if (p_data_start == NULL)
		return;

	if (mType == 0)
		sprintf(filename, MARKER_FILE_0);
	else if (mType == 1)
		sprintf(filename, MARKER_FILE_1);
	else {
		fprintf(stderr, "Error: Invalid input \"mType\" = %d\n", mType);
		fprintf(stderr, "       in function Marker_endFile(). Use 0, or 1.\n\n");
		exit(EXIT_FAILURE);
	}

	FILE *xmf = fopen(filename, "a");

	fprintf(xmf, "    </Grid>\n");
	fprintf(xmf, "  </Domain>\n");
	fprintf(xmf, "</Xdmf>\n");
	fclose(xmf);

	// Display which file was created
	printf("Wrote %s\n", filename);
}




/******************************************************************************/
/*
 Writes data for entire linked list data_start for a single timestep
 */
/******************************************************************************/
void Marker_writeTimestep(DataParticle *p_data_start, int mType, int index,
		double time, Parameters *params) {

	if (p_data_start == NULL)
		return;

	int iter = FILE_INDEX(index);
	int N_pts = params->N_pt_marker[index];
	int N_poly = params->N_poly_marker[index];
	int size_poly = params->size_poly_marker;

	Data3d *data3d;
	FILE *xmf;
	char gridname[50], groupname[50], typology_type[50];

	// Append to file
	if (mType == 0) {
		xmf = fopen(MARKER_FILE_0, "a");
		sprintf(gridname, "MarkerGrid");
		sprintf(groupname, MARKER_GROUP_0);
	}
	else if (mType == 1) {
		xmf = fopen(MARKER_FILE_1, "a");
		sprintf(gridname, "MarkerGrid");
		sprintf(groupname, MARKER_GROUP_1);
	}
	else {
		fprintf(stderr, "Error: Invalid input \"mType\" = %d\n", mType);
		fprintf(stderr, "       in function Marker_writeTimestep(). Use 0 or 1.\n\n");
		exit(EXIT_FAILURE);
	}

	if (size_poly == 3) {
		sprintf(typology_type, "Triangle");
	}
	else if (size_poly == 4) {
		sprintf(typology_type, "Quadrilateral");
	}
	else {
		fprintf(stderr, "Error: Invalid unsupported polygon size \"size_poly\" = %d\n", size_poly);
		fprintf(stderr, "       in function Marker_writeTimestep(). Use 3 or 4.\n\n");
		exit(EXIT_FAILURE);
	}

	// Describe time value and grid
	fprintf(xmf, "      <Grid Name=\"%s_%d\">\n", gridname, index);
	fprintf(xmf, "        <Time Value=\"%f\"/>\n", time);
	fprintf(xmf, "        <Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n",
		typology_type, N_poly);
	fprintf(xmf, "          <DataItem Format=\"HDF\" NumberType=\"Double\" Dimensions=\"%d %d\">\n",
		N_poly, size_poly);
	fprintf(xmf, "            Marker_%d.h5:%s/topology\n", iter, groupname);
	fprintf(xmf, "          </DataItem>\n");
	fprintf(xmf, "        </Topology>\n");
	fprintf(xmf, "\n");
	fprintf(xmf, "        <Geometry Type=\"XYZ\">\n");
	fprintf(xmf, "          <DataItem Format=\"HDF\" Dimensions=\"%d %d\">\n", N_pts, 3);
	fprintf(xmf, "            Marker_%d.h5:%s/X\n", iter, groupname);
	fprintf(xmf, "          </DataItem>\n");
	fprintf(xmf, "        </Geometry>\n");

	// Write attributes for any data that matches this node type
	DataParticle *p_data = p_data_start;
	while (p_data != NULL) {
		Marker_writeAttribute(xmf, p_data, mType, iter, N_pts);
		p_data = p_data -> next;
	}

	// Close up timestep
	fprintf(xmf, "      </Grid>\n");
	fprintf(xmf, "\n");

	fclose(xmf);
}




/******************************************************************************/
/*
 Writes data for entire linked list data_start for a single timestep
 */
/******************************************************************************/
void Marker_writeAttribute(FILE *xmf, DataParticle *p_data, int mType, int iter, int N_pts) {

	char att_type[50], att_name[50], groupname[50];

	if (mType == 0) {
		sprintf(att_name, "%s", p_data->dataname);
		sprintf(groupname, MARKER_GROUP_0"/%s", p_data->dataname);
	}
	else if (mType == 1) {
		sprintf(att_name, "%s", p_data->dataname);
		sprintf(groupname, MARKER_GROUP_1"/%s", p_data->dataname);
	}
	else {
		fprintf(stderr, "Error: Invalid input \"mType\" = %d\n", mType);
		fprintf(stderr, "       in function startFileParticle(). Use 0 or 1.\n\n");
		exit(EXIT_FAILURE);
	}

	if (p_data->length == 1) {
		sprintf(att_type, "Scalar");
	}
	else if (p_data->length == 3) {
		sprintf(att_type, "Vector");
	}
	else {
		fprintf(stderr, "Error: Invalid data length \"p_data->length = %d\"\n", p_data->length);
		fprintf(stderr, "       for Marker_%d.h5:%s.\n\n", iter, groupname);
		fprintf(stderr, "       in function Particle_writeAttribute().\n\n");
		exit(EXIT_FAILURE);
	}

	fprintf(xmf, "\n");
	fprintf(xmf, "        <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n",
		att_name, att_type);
	fprintf(xmf, "          <DataItem Format=\"HDF\" NumberType=\"Double\" Dimensions=\"%d %d\">\n",
		N_pts, p_data->length);
	fprintf(xmf, "            Marker_%d.h5:%s\n", iter, groupname);
	fprintf(xmf, "          </DataItem>\n");
	fprintf(xmf, "        </Attribute>\n");
}
