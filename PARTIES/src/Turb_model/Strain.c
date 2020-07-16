#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Conc.h"
#include "Cart3d.h"
#include "Immersed.h"
#include "Velocity.h"
#include "Subgrid.h"
#include "Strain.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/******************************************************************************/
/*
 * Calculate the magnitude of the strain rate tensor
 * When dynamic Smagorisnky model is used, individual
 * elements in the tensor is also stored
 */
/******************************************************************************/
void Strain_rate_magnitude(Velocity *uvel, Velocity *vvel, Velocity *wvel,
		MAC_grid *grid, Parameters *params, Strain_rate *st_rate) {

	double ***u, ***v, ***w;
	double ***strain; 

	int NX, NY, NZ;
	int Is, Js, Ks;
	int Ie, Je, Ke;
	int i_start, i_end; 
	int j_start, j_end; 
	int k_start, k_end;
	int ierr; 
	int i, j, k;
	double *i2dx_c, *i2dy_c, *i2dz_c;
	double *idx_u, *idy_v, *idz_w;
	double *idx_c, *idy_c, *idz_c;
	double idz;
	double dudyEN, dudyWN, dudyES, dudyWS;
	double dudzET, dudzWT, dudzEB, dudzWB;
	double dvdxNE, dvdxNW, dvdxSE, dvdxSW;
	double dwdxTE, dwdxTW, dwdxBE, dwdxBW;
	double dwdyTN, dwdyTS, dwdyBN, dwdyBS;
	double sxy, sxz, syz; 
	double dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz; 
#ifdef  SMAG_DYNAMIC
	double ***s11, ***s22, ***s33, ***s12, ***s13, ***s23; 
#endif
	char bin_filename[50];
	FILE *fid, *fid1;
/*
	char bin_filename[50];
	FILE *fid, *fid1;

	sprintf(bin_filename,"straininfo%d.dat",params->rank); 
	fid = fopen(bin_filename,"w"); 
	sprintf(bin_filename,"strainmag%d.dat",params->rank); 
	fid1 = fopen(bin_filename,"w"); 
*/

	// Same for all quantities
	NX = grid->NX;
	NY = grid->NY;
	NZ = grid->NZ;

	// Start index of bottom-left-back corner on current processor
	Is = grid->G_Is;
	Js = grid->G_Js;
	Ks = grid->G_Ks;

	// End index of top-right-front corner on current processor
	Ie = grid->G_Ie;
	Je = grid->G_Je;
	Ke = grid->G_Ke;

	// indices start and end on current processor
	i_start = max(1,Is); // i=0, i=NX-1 are not included
	j_start = max(1,Js);
	k_start = max(1,Ks); 

#ifdef XPERIODIC 
	i_start = Is;
#endif
#ifdef ZPERIODIC 
	k_start = Ks;
#endif
		
	// exclude the half cell added
	i_end   = min(NX-1, Ie); 
	j_end   = min(NY-1, Je); 
	k_end   = min(NZ-1, Ke);

	// Get the local velocities at the location where they are defined
	u = uvel->data;
	v = vvel->data;
	w = wvel->data;
	strain = st_rate->strain; 

#ifdef SMAG_DYNAMIC
	s11 = st_rate->s11; 
	s22 = st_rate->s22; 
	s33 = st_rate->s33; 
	s12 = st_rate->s12; 
	s13 = st_rate->s13; 
	s23 = st_rate->s23; 
#endif

	i2dx_c = grid->i2dx_c;
	i2dy_c = grid->i2dy_c;
	i2dz_c = grid->i2dz_c;
	idx_u = grid->idx_u;
	idy_v = grid->idy_v;
	idz_w = grid->idz_w;

	idx_c = grid->idx_c;
	idy_c = grid->idy_c;
	idz_c = grid->idz_c;

	for (k=k_start; k<k_end; k++) {
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){

				// u velocity derivatives, du/dx, du/dy, du/dz
				dudx = ( u[k][j][i+1]-u[k][j][i])*idx_u[i];
				dudyEN = (u[k][j+1][i+1] - u[k][j][i+1]   );
				dudyWN = (u[k][j+1][i]   - u[k][j][i]     );
				dudyES = (u[k][j][i+1]   - u[k][j-1][i+1] );
				dudyWS = (u[k][j][i]     - u[k][j-1][i]   );
				dudy = 0.25*( (dudyEN + dudyWN)*idy_c[j] + (dudyES + dudyWS)*idy_c[j-1] ); 

				dudzET = ( u[k+1][j][i+1] - u[k][j][i+1]) ;  
				dudzWT = ( u[k+1][j][i] - u[k][j][i]) ;  
				if (k != 0) {
					dudzEB = ( u[k][j][i+1] - u[k-1][j][i+1]) ;
					dudzWB = ( u[k][j][i] - u[k-1][j][i]) ;  
					idz = idz_c[k-1]; 
				}
				else {
					dudzEB = ( u[k][j][i+1] - u[k-1][j][i+1]) ;
					dudzWB = ( u[k][j][i] - u[k-1][j][i]) ;  
					idz = idz_c[k]; 
				}
				dudz = 0.25 * ( (dudzET + dudzWT)*idz_c[k] + (dudzEB + dudzWB)*idz );

				// v velocity derivatives, dv/dx, dv/dy, dv/dz
				dvdx = 0.5*( (v[k][j+1][i+1] + v[k][j][i+1] ) - (v[k][j+1][i-1] + v[k][j][i-1]) )*i2dx_c[i];
				dvdxNE =  ( v[k][j+1][i+1] - v[k][j+1][i]   );
				dvdxNW =  ( v[k][j+1][i]   - v[k][j+1][i-1] );
				dvdxSE =  ( v[k][j][i+1] - v[k][j][i]   );
				dvdxSW =  ( v[k][j][i]   - v[k][j][i-1] );
				dvdx   = 0.25*( (dvdxNE + dvdxSE)*idx_c[i] + (dvdxNW + dvdxSW)*idx_c[i-1] ); 

				dvdy = ( v[k][j+1][i] - v[k][j][i] )*idy_v[j];
				dvdz = 0.5*( (v[k+1][j+1][i] + v[k+1][j][i] ) - (v[k-1][j+1][i] + v[k-1][j][i]) )*i2dz_c[k];

				// w velocity derivatives, dw/dx, dw/dy, dw/dz
				dwdx = 0.5*( (w[k+1][j][i+1] + w[k][j][i+1] ) - (w[k+1][j][i-1] + w[k][j][i-1]) )*i2dx_c[i];
				dwdyTN =  (w[k+1][j+1][i] - w[k+1][j][i] ) ;
				dwdyBN =  (w[k][j+1][i] - w[k][j][i] )     ;
				dwdyTS =  (w[k+1][j][i] - w[k+1][j-1][i] ) ;
				dwdyBS =  (w[k][j][i] - w[k][j-1][i]     ) ;
				dwdy = 0.25*( (dwdyTN + dwdyBN)*idy_c[j] + (dwdyTS + dwdyBS)*idy_c[j-1]); 

				dwdz = ( w[k+1][j][i]-w[k][j][i])*idz_w[k];

				// Generate large-scale strain-rate tensor
				sxy = 0.5*(dudy+dvdx);
				sxz = 0.5*(dudz+dwdx);
				syz = 0.5*(dvdz+dwdy);

				strain[k][j][i] = dudx*dudx + dvdy*dvdy + dwdz*dwdz + 2.0*( sxy*sxy + sxz*sxz + syz*syz );
				strain[k][j][i]=sqrt(2.*strain[k][j][i]); 
#ifdef SMAG_DYNAMIC
				s11[k][j][i] = dudx; 
				s22[k][j][i] = dvdy; 
				s33[k][j][i] = dwdz; 
				s12[k][j][i] = sxy; 
				s13[k][j][i] = sxz; 
				s23[k][j][i] = syz; 
#endif

			} // for i
		} // for j
	} // for k

#ifndef XPERIODIC
	if (Is == 0) {
		i=0;
		for (k=k_start; k<k_end; k++) {
			for (j=j_start; j<j_end; j++) {
				
				//--------------------------------------------------------------
				// u velocity derivatives, du/dx, du/dy, du/dz
				//--------------------------------------------------------------
				dudx = ( u[k][j][i+1]-u[k][j][i])*idx_u[i];
				dudy = 0.5*( (u[k][j+1][i+1] + u[k][j+1][i] ) - (u[k][j-1][i+1] + u[k][j-1][i]) )*i2dy_c[j];
				dudz = 0.5*( (u[k+1][j][i+1] + u[k+1][j][i] ) - (u[k-1][j][i+1] + u[k-1][j][i]) )*i2dz_c[k];
				
				//--------------------------------------------------------------
				// v velocity derivatives, dv/dx, dv/dy, dv/dz
				// 	Use one sided finite difference
				//--------------------------------------------------------------
				dvdx = 0.5*( (v[k][j+1][i+1] + v[k][j][i+1] ) - (v[k][j+1][i] + v[k][j][i]) )*i2dx_c[i];
				dvdy = ( v[k][j+1][i] - v[k][j][i] )*idy_v[j];
				dvdz = 0.5*( (v[k+1][j+1][i] + v[k+1][j][i] ) - (v[k-1][j+1][i] + v[k-1][j][i]) )*i2dz_c[k];
				
				//--------------------------------------------------------------
				// w velocity derivatives, dw/dx, dw/dy, dw/dz
				// 	Use one sided finite difference
				//--------------------------------------------------------------
				dwdx = 0.5*( (w[k+1][j][i+1] + w[k][j][i+1] ) - (w[k+1][j][i] + w[k][j][i]) )*i2dx_c[i];
				dwdy = 0.5*( (w[k+1][j+1][i] + w[k][j+1][i] ) - (w[k+1][j-1][i] + w[k][j-1][i]) )*i2dy_c[j]; 
				dwdz = ( w[k+1][j][i]-w[k][j][i])*idz_w[k];
				
				//--------------------------------------------------------------
				// Generate large-scale strain-rate tensor
				//--------------------------------------------------------------
				sxy = 0.5*(dudy+dvdx);
				sxz = 0.5*(dudz+dwdx);
				syz = 0.5*(dvdz+dwdy);

				strain[k][j][i] = dudx*dudx + dvdy*dvdy + dwdz*dwdz + 2.0*( sxy*sxy + sxz*sxz + syz*syz );
				strain[k][j][i]=sqrt(2.*strain[k][j][i]); 
	#ifdef SMAG_DYNAMIC
				s11[k][j][i] = dudx; 
				s22[k][j][i] = dvdy; 
				s33[k][j][i] = dwdz; 
				s12[k][j][i] = sxy; 
				s13[k][j][i] = sxz; 
				s23[k][j][i] = syz; 
	#endif

//				fprintf(fid,"%d %d %d %e %e %e %e %e %e %e %e %e %e\n",i,j,k,
//						dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,strain[k][j][i]);
//				fprintf(fid1,"%e\n",strain[k][j][i]);
			} // for j
		} // for k
	} // if Is==0
#endif // notXPERIODIC

	if (Js == 0) {
		j=0;
		for (k=k_start; k<k_end; k++) {
			for (i=i_start; i<i_end; i++){
				
				//--------------------------------------------------------------
				// u velocity derivatives, du/dx, du/dy, du/dz
				//--------------------------------------------------------------
				dudx = ( u[k][j][i+1]-u[k][j][i])*idx_u[i];
				dudyEN =  (u[k][j+1][i+1] - u[k][j][i+1] ) * idy_c[j];
				dudyWN =  (u[k][j+1][i] - u[k][j][i] )     * idy_c[j];
#if defined BOTTOM_WALL_VELOCITY_NOSLIP || defined BOTTOM_WALL_SCHUMANN
				dudyES =  (2.*u[k][j][i+1]  ) * idy_v[j];
				dudyWS =  (2.*u[k][j][i]    ) * idy_v[j];
#else
				dudyES =  0.0;
				dudyWS =  0.0;
#endif
				dudy = 0.25*(dudyEN + dudyWN + dudyES + dudyWS); 

				dudzET = ( u[k+1][j][i+1] - u[k][j][i+1]) * idz_c[k];  
				dudzWT = ( u[k+1][j][i] - u[k][j][i]) * idz_c[k];  
				if (k != 0) {
					dudzEB = ( u[k][j][i+1] - u[k-1][j][i+1]) * idz_c[k-1];
					dudzWB = ( u[k][j][i] - u[k-1][j][i]) * idz_c[k-1];  
				}
				else {
					dudzEB = ( u[k][j][i+1] - u[k-1][j][i+1]) * idz_c[k];
					dudzWB = ( u[k][j][i] - u[k-1][j][i]) * idz_c[k];  
				}
				dudz = 0.25 * ( dudzET + dudzWT + dudzEB + dudzWB );
				
				//--------------------------------------------------------------
				// v velocity derivatives, dv/dx, dv/dy, dv/dz
				//--------------------------------------------------------------
				dvdx = 0.5*( (v[k][j+1][i+1] + v[k][j][i+1] ) - (v[k][j+1][i-1] + v[k][j][i-1]) )*i2dx_c[i];
				dvdy = ( v[k][j+1][i] - v[k][j][i] )*idy_v[j];
				dvdz = 0.5*( (v[k+1][j+1][i] + v[k+1][j][i] ) - (v[k-1][j+1][i] + v[k-1][j][i]) )*i2dz_c[k];
				
				//--------------------------------------------------------------
				// w velocity derivatives, dw/dx, dw/dy, dw/dz
				//--------------------------------------------------------------
				dwdx = 0.5*( (w[k+1][j][i+1] + w[k][j][i+1] ) - (w[k+1][j][i-1] + w[k][j][i-1]) )*i2dx_c[i];

				dwdyTN =  (w[k+1][j+1][i] - w[k+1][j][i] ) * idy_c[j];
				dwdyBN =  (w[k][j+1][i] - w[k][j][i] )     * idy_c[j];
#if defined BOTTOM_WALL_VELOCITY_NOSLIP || defined BOTTOM_WALL_SCHUMANN
				dwdyTS =  (2.*w[k+1][j][i]    ) * idy_v[j];
				dwdyBS =  (2.*w[k][j][i]      ) * idy_v[j];
#else
				dwdyTS = 0.0;
				dwdyBS = 0.0;
#endif
				dwdy = 0.25*(dwdyTN + dwdyBN + dwdyTS + dwdyBS); 
				dwdz = ( w[k+1][j][i]-w[k][j][i])*idz_w[k];
				
				//--------------------------------------------------------------
				// Generate large-scale strain-rate tensor
				//--------------------------------------------------------------
				sxy = 0.5*(dudy+dvdx);
				sxz = 0.5*(dudz+dwdx);
				syz = 0.5*(dvdz+dwdy);

				strain[k][j][i] = dudx*dudx + dvdy*dvdy + dwdz*dwdz + 2.0*( sxy*sxy + sxz*sxz + syz*syz );
				strain[k][j][i]=sqrt(2.*strain[k][j][i]); 
#ifdef SMAG_DYNAMIC
				s11[k][j][i] = dudx; 
				s22[k][j][i] = dvdy; 
				s33[k][j][i] = dwdz; 
				s12[k][j][i] = sxy; 
				s13[k][j][i] = sxz; 
				s23[k][j][i] = syz; 
#endif 
//				fprintf(fid,"%d %d %d %e %e %e %e %e %e %e %e %e %e\n",i,j,k,
//						dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,strain[k][j][i]);
//				fprintf(fid1,"%e\n",strain[k][j][i]);

			} // for i
		} // for k
	} // if Js==0

#ifndef ZPERIODIC
	if (Ks == 0) {
		k=0;
		for (j=j_start; j<j_end; j++) {
			for (i=i_start; i<i_end; i++){
				
				//--------------------------------------------------------------
				// u velocity derivatives, du/dx, du/dy, du/dz
				//--------------------------------------------------------------
				dudx = ( u[k][j][i+1]-u[k][j][i])*idx_u[i];
				dudy = 0.5*( (u[k][j+1][i+1] + u[k][j+1][i] ) - (u[k][j-1][i+1] + u[k][j-1][i]) )*i2dy_c[j];
				// 	Use one sided finite difference 
				dudz = 0.5*( (u[k+1][j][i+1] + u[k+1][j][i] ) - (u[k][j][i+1] + u[k][j][i]) )*i2dz_c[k];
				
				//--------------------------------------------------------------
				// v velocity derivatives, dv/dx, dv/dy, dv/dz
				//--------------------------------------------------------------
				dvdx = 0.5*( (v[k][j+1][i+1] + v[k][j][i+1] ) - (v[k][j+1][i-1] + v[k][j][i-1]) )*i2dx_c[i];
				dvdy = ( v[k][j+1][i] - v[k][j][i] )*idy_v[j];
				// 	Use one sided finite difference 
				dvdz = 0.5*( (v[k+1][j+1][i] + v[k+1][j][i] ) - (v[k][j+1][i] + v[k][j][i]) )*i2dz_c[k];
				
				//--------------------------------------------------------------
				// w velocity derivatives, dw/dx, dw/dy, dw/dz
				//--------------------------------------------------------------
				dwdx = 0.5*( (w[k+1][j][i+1] + w[k][j][i+1] ) - (w[k+1][j][i-1] + w[k][j][i-1]) )*i2dx_c[i];
				dwdy = 0.5*( (w[k+1][j+1][i] + w[k][j+1][i] ) - (w[k+1][j-1][i] + w[k][j-1][i]) )*i2dy_c[j]; 
				dwdz = ( w[k+1][j][i]-w[k][j][i])*idz_w[k];
				
				//--------------------------------------------------------------
				// Generate large-scale strain-rate tensor
				//--------------------------------------------------------------
				sxy = 0.5*(dudy+dvdx);
				sxz = 0.5*(dudz+dwdx);
				syz = 0.5*(dvdz+dwdy);

				strain[k][j][i] = dudx*dudx + dvdy*dvdy + dwdz*dwdz + 2.0*( sxy*sxy + sxz*sxz + syz*syz );
				strain[k][j][i]=sqrt(2.*strain[k][j][i]); 
	#ifdef SMAG_DYNAMIC
				s11[k][j][i] = dudx; 
				s22[k][j][i] = dvdy; 
				s33[k][j][i] = dwdz; 
				s12[k][j][i] = sxy; 
				s13[k][j][i] = sxz; 
				s23[k][j][i] = syz; 
	#endif

//				fprintf(fid,"%d %d %d %e %e %e %e %e %e %e %e %e %e\n",i,j,k,
//						dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,strain[k][j][i]);
//				fprintf(fid1,"%e\n",strain[k][j][i]);
			}
		}
	}
#endif // notZPERIODIC

	if ( ((Is == 0) && (Js ==0)) && (Is != i_start) )  {
		i=0;
		j=0;
		for (k=k_start; k<k_end; k++) {
			
			//--------------------------------------------------------------
			// u velocity derivatives, du/dx, du/dy, du/dz
			//--------------------------------------------------------------
			dudx = ( u[k][j][i+1]-u[k][j][i])*idx_u[i];
			// Use one sided finite difference 
			dudy = 0.5*( (u[k][j+1][i+1] + u[k][j+1][i] ) - (u[k][j][i+1] + u[k][j][i]) )*i2dy_c[j];
			dudz = 0.5*( (u[k+1][j][i+1] + u[k+1][j][i] ) - (u[k-1][j][i+1] + u[k-1][j][i]) )*i2dz_c[k];
			
			//--------------------------------------------------------------
			// v velocity derivatives, dv/dx, dv/dy, dv/dz
			//--------------------------------------------------------------
			// Use one sided finite difference 
			dvdx = 0.5*( (v[k][j+1][i+1] + v[k][j][i+1] ) - (v[k][j+1][i] + v[k][j][i]) )*i2dx_c[i];
			dvdy = ( v[k][j+1][i] - v[k][j][i] )*idy_v[j];
			dvdz = 0.5*( (v[k+1][j+1][i] + v[k+1][j][i] ) - (v[k-1][j+1][i] + v[k-1][j][i]) )*i2dz_c[k];
			
			//--------------------------------------------------------------
			// w velocity derivatives, dw/dx, dw/dy, dw/dz
			//--------------------------------------------------------------
			// Use one sided finite difference 
			dwdx = 0.5*( (w[k+1][j][i+1] + w[k][j][i+1] ) - (w[k+1][j][i] + w[k][j][i]) )*i2dx_c[i];
			// Use one sided finite difference 
			dwdy = 0.5*( (w[k+1][j+1][i] + w[k][j+1][i] ) - (w[k+1][j][i] + w[k][j][i]) )*i2dy_c[j]; 
			dwdz = ( w[k+1][j][i]-w[k][j][i])*idz_w[k];
			
			//--------------------------------------------------------------
			// Generate large-scale strain-rate tensor
			//--------------------------------------------------------------
			sxy = 0.5*(dudy+dvdx);
			sxz = 0.5*(dudz+dwdx);
			syz = 0.5*(dvdz+dwdy);

			strain[k][j][i] = dudx*dudx + dvdy*dvdy + dwdz*dwdz + 2.0*( sxy*sxy + sxz*sxz + syz*syz );
			strain[k][j][i]=sqrt(2.*strain[k][j][i]); 
#ifdef SMAG_DYNAMIC
			s11[k][j][i] = dudx;
			s22[k][j][i] = dvdy;
			s33[k][j][i] = dwdz; 
			s12[k][j][i] = sxy; 
			s13[k][j][i] = sxz; 
			s23[k][j][i] = syz; 
#endif
//			fprintf(fid,"%d %d %d %e %e %e %e %e %e %e %e %e %e\n",i,j,k,
//					dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,strain[k][j][i]);
//			fprintf(fid1,"%e\n",strain[k][j][i]);
		}
	}

	if ( ((Js == 0) && (Ks ==0)) && (Ks != k_start) )  {
		j=0;
		k=0;
		for (i=i_start; i<i_end; i++){

			//--------------------------------------------------------------
			// u velocity derivatives, du/dx, du/dy, du/dz
			//--------------------------------------------------------------
			dudx = ( u[k][j][i+1]-u[k][j][i])*idx_u[i];
			// Use one sided finite difference 
			dudy = 0.5*( (u[k][j+1][i+1] + u[k][j+1][i] ) - (u[k][j][i+1] + u[k][j][i]) )*i2dy_c[j];
			// Use one sided finite difference 
			dudz = 0.5*( (u[k+1][j][i+1] + u[k+1][j][i] ) - (u[k][j][i+1] + u[k][j][i]) )*i2dz_c[k];
			
			//--------------------------------------------------------------
			// v velocity derivatives, dv/dx, dv/dy, dv/dz
			//--------------------------------------------------------------
			dvdx = 0.5*( (v[k][j+1][i+1] + v[k][j][i+1] ) - (v[k][j+1][i-1] + v[k][j][i-1]) )*i2dx_c[i];
			dvdy = ( v[k][j+1][i] - v[k][j][i] )*idy_v[j];
			// Use one sided finite difference 
			dvdz = 0.5*( (v[k+1][j+1][i] + v[k+1][j][i] ) - (v[k][j+1][i] + v[k][j][i]) )*i2dz_c[k];
			
			//--------------------------------------------------------------
			// w velocity derivatives, dw/dx, dw/dy, dw/dz
			//--------------------------------------------------------------
			dwdx = 0.5*( (w[k+1][j][i+1] + w[k][j][i+1] ) - (w[k+1][j][i-1] + w[k][j][i-1]) )*i2dx_c[i];
			// Use one sided finite difference 
			dwdy = 0.5*( (w[k+1][j+1][i] + w[k][j+1][i] ) - (w[k+1][j][i] + w[k][j][i]) )*i2dy_c[j]; 
			dwdz = ( w[k+1][j][i]-w[k][j][i])*idz_w[k];
			
			//--------------------------------------------------------------
			// Generate large-scale strain-rate tensor
			//--------------------------------------------------------------
			sxy = 0.5*(dudy+dvdx);
			sxz = 0.5*(dudz+dwdx);
			syz = 0.5*(dvdz+dwdy);

			strain[k][j][i] = dudx*dudx + dvdy*dvdy + dwdz*dwdz + 2.0*( sxy*sxy + sxz*sxz + syz*syz );
			strain[k][j][i]=sqrt(2.*strain[k][j][i]); 
#ifdef SMAG_DYNAMIC
			s11[k][j][i] = dudx;
			s22[k][j][i] = dvdy; 
			s33[k][j][i] = dwdz; 
			s12[k][j][i] = sxy; 
			s13[k][j][i] = sxz; 
			s23[k][j][i] = syz; 
#endif 
//			fprintf(fid,"%d %d %d %e %e %e %e %e %e %e %e %e %e\n",i,j,k,
//					dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,strain[k][j][i]);
		}
	}
	if ( ((Ks == 0) && (Is ==0)) &&  ( (Ks != k_start) && (Is != i_start) ) )  {
		k=0;
		i=0;
		for (j=j_start; j<j_end; j++) {

			// u velocity derivatives, du/dx, du/dy, du/dz
			dudx = ( u[k][j][i+1]-u[k][j][i])*idx_u[i];
			dudy = 0.5*( (u[k][j+1][i+1] + u[k][j+1][i] ) - (u[k][j-1][i+1] + u[k][j-1][i]) )*i2dy_c[j];
			// Use one sided finite difference 
			dudz = 0.5*( (u[k+1][j][i+1] + u[k+1][j][i] ) - (u[k][j][i+1] + u[k][j][i]) )*i2dz_c[k];

			// v velocity derivatives, dv/dx, dv/dy, dv/dz
			// Use one sided finite difference 
			dvdx = 0.5*( (v[k][j+1][i+1] + v[k][j][i+1] ) - (v[k][j+1][i] + v[k][j][i]) )*i2dx_c[i];
			dvdy = ( v[k][j+1][i] - v[k][j][i] )*idy_v[j];
			// Use one sided finite difference 
				dvdz = 0.5*( (v[k+1][j+1][i] + v[k+1][j][i] ) - (v[k][j+1][i] + v[k][j][i]) )*i2dz_c[k];

			// w velocity derivatives, dw/dx, dw/dy, dw/dz
			// Use one sided finite difference 
			dwdx = 0.5*( (w[k+1][j][i+1] + w[k][j][i+1] ) - (w[k+1][j][i] + w[k][j][i]) )*i2dx_c[i];
			dwdy = 0.5*( (w[k+1][j+1][i] + w[k][j+1][i] ) - (w[k+1][j-1][i] + w[k][j-1][i]) )*i2dy_c[j]; 
			dwdz = ( w[k+1][j][i]-w[k][j][i])*idz_w[k];

			// Generate large-scale strain-rate tensor
			sxy = 0.5*(dudy+dvdx);
			sxz = 0.5*(dudz+dwdx);
			syz = 0.5*(dvdz+dwdy);

			strain[k][j][i] = dudx*dudx + dvdy*dvdy + dwdz*dwdz + 2.0*( sxy*sxy + sxz*sxz + syz*syz );
			strain[k][j][i]=sqrt(2.*strain[k][j][i]);
#ifdef SMAG_DYNAMIC
			s11[k][j][i] = dudx; 
			s22[k][j][i] = dvdy; 
			s33[k][j][i] = dwdz; 
			s12[k][j][i] = sxy; 
			s13[k][j][i] = sxz; 
			s23[k][j][i] = syz; 
#endif 
//			fprintf(fid,"%d %d %d %e %e %e %e %e %e %e %e %e %e\n",i,j,k,
//					dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,strain[k][j][i]);
//			fprintf(fid1,"%e\n",strain[k][j][i]); 
		}
	}

	if ((Is == 0) && (Js==0) && (Ks==0))
		strain[0][0][0] = 0.5*(strain[1][0][0] + strain[0][1][0] );
//		fclose(fid); 
//		fclose(fid1); 

//	Communication_update_ghost_nodes_flow_variable(smag->strain, grid, params, 'c', 1);

#ifdef SMAG_DYNAMIC
/*
	Communication_update_ghost_nodes_flow_variable(smag->s11, grid, params, 'c', 1);
	Communication_update_ghost_nodes_flow_variable(smag->s22, grid, params, 'c', 1);
	Communication_update_ghost_nodes_flow_variable(smag->s33, grid, params, 'c', 1);
	Communication_update_ghost_nodes_flow_variable(smag->s12, grid, params, 'c', 1);
	Communication_update_ghost_nodes_flow_variable(smag->s13, grid, params, 'c', 1);
	Communication_update_ghost_nodes_flow_variable(smag->s23, grid, params, 'c', 1);
*/
#endif
	return;

}




/******************************************************************************/
