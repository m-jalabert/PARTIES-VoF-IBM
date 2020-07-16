#include <stdlib.h>
#include <stdio.h>
#include "DataTypes.h"
#include "Particle.h"

#define MOBILE_FILE "Reader_p_mobile.xmf"
#define FIXED_FILE "Reader_p_fixed.xmf"


/******************************************************************************/
/*
 Writes data for entire linked list data_start for a single timestep
 */
/******************************************************************************/
void Particle_startFile(DataParticle *p_data_start, char pType, Parameters *params) {

	FILE *xmf;

	if (p_data_start == NULL)
		return;

	if (pType == 'm')
		xmf = fopen(MOBILE_FILE, "w");
	else if (pType == 'f')
		xmf = fopen(FIXED_FILE, "w");
	else {
		fprintf(stderr, "Error: Invalid input \"pType\" = '%c'\n", pType);
		fprintf(stderr, "       in function startFileParticle(). Use 'm' or 'f'.\n\n");
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
void Particle_endFile(DataParticle *p_data_start, char pType, Parameters *params) {

	char filename[50];

	if (p_data_start == NULL)
		return;

	if (pType == 'm')
		sprintf(filename, MOBILE_FILE);
	else if (pType == 'f')
		sprintf(filename, FIXED_FILE);
	else {
		fprintf(stderr, "Error: Invalid input \"pType\" = '%c'\n", pType);
		fprintf(stderr, "       in function startFileParticle(). Use 'm' or 'f'.\n\n");
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
void Particle_writeTimestep(DataParticle *p_data_start, char pType, int index,
		double time, int *N_particles, Parameters *params) {

	if (p_data_start == NULL)
		return;

	Data3d *data3d;
	FILE *xmf;
	char gridname[50], groupname[50];

	int iter = FILE_INDEX(index);

	// Append to file
	if (pType == 'm') {
		xmf = fopen(MOBILE_FILE, "a");
		sprintf(gridname, "MobileGrid");
		sprintf(groupname, MOBILE_GROUP);
	}
	else if (pType == 'f') {
		xmf = fopen(FIXED_FILE, "a");
		sprintf(gridname, "FixedGrid");
		sprintf(groupname, FIXED_GROUP);
	}
	else {
		fprintf(stderr, "Error: Invalid input \"pType\" = '%c'\n", pType);
		fprintf(stderr, "       in function startFileParticle(). Use 'm' or 'f'.\n\n");
		exit(EXIT_FAILURE);
	}

	// Describe time value and grid
	fprintf(xmf, "      <Grid Name=\"%s_%d\">\n", gridname, index);
	fprintf(xmf, "        <Time Value=\"%f\"/>\n", time);
	fprintf(xmf, "        <Topology Type=\"Polyvertex\" NumberOfElements=\"%d\" />\n", N_particles[index]);
	fprintf(xmf, "\n");
	fprintf(xmf, "        <Geometry Type=\"XYZ\">\n");
	fprintf(xmf, "          <DataItem Format=\"HDF\" Dimensions=\"%d %d\">\n", N_particles[index], 3);
	fprintf(xmf, "            Particle_%d.h5:%s/X\n", iter, groupname);
	fprintf(xmf, "          </DataItem>\n");
	fprintf(xmf, "        </Geometry>\n");

	// Write attributes for any data that matches this node type
	DataParticle *p_data = p_data_start;
	while (p_data != NULL) {
		Particle_writeAttribute(xmf, p_data, pType, iter, N_particles[index]);
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
void Particle_writeAttribute(FILE *xmf, DataParticle *p_data, char pType, int iter, int N_particles) {

	char att_type[50], att_name[50], groupname[50];

	if (pType == 'm') {
		sprintf(att_name, "%s", p_data->dataname);
		sprintf(groupname, MOBILE_GROUP"/%s", p_data->dataname);
	}
	else if (pType == 'f') {
		sprintf(att_name, "%s", p_data->dataname);
		sprintf(groupname, FIXED_GROUP"/%s", p_data->dataname);
	}
	else {
		fprintf(stderr, "Error: Invalid input \"pType\" = '%c'\n", pType);
		fprintf(stderr, "       in function startFileParticle(). Use 'm' or 'f'.\n\n");
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
		fprintf(stderr, "       in function Particle_writeAttribute().\n\n");
		exit(EXIT_FAILURE);
	}

	fprintf(xmf, "\n");
	fprintf(xmf, "        <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n", att_name, att_type);
	fprintf(xmf, "          <DataItem Format=\"HDF\" NumberType=\"Double\" Dimensions=\"%d %d\">\n", N_particles, p_data->length);
	fprintf(xmf, "            Particle_%d.h5:%s\n", iter, groupname);
	fprintf(xmf, "          </DataItem>\n");
	fprintf(xmf, "        </Attribute>\n");
}
