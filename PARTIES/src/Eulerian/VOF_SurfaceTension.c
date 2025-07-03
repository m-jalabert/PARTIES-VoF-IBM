/******************************************************************************
 * VOF_SurfaceTension.c
 * Surface tension force computation for the Volume-Of-Fluid (VOF) method using CSF
 ******************************************************************************/

 #include "VolumeFraction.h"  
 #include "definitions.h"
 #include "Boundary.h"
 #include "DataTypes.h"
 #include "Conc.h"
 #include "Velocity.h"
 #include "MyMath.h"
 #include "Memory.h"
 #include "Grid.h"
 #include "Communication.h"
 #include "Immersed.h"
 #include "Display.h"
 #include "Array.h"
 #include "Cart3d.h"
 #include "Lagrangian.h"
 
 #include <stdlib.h>
 #include <stdbool.h> 
 #include <stdio.h>
 #include <math.h>
 #include <time.h>
 #include <assert.h>
 #include <string.h>
 #include <float.h> 
 
 
 #ifdef VOF_PLIC




#define nodata (1e20)

//------------------------------------------------------------------------------
// Helper: fourth‐order smoothing kernel D(x)
//   D(x) = (15/(16·h))·[ (x/h)^4 – 2·(x/h)^2 + 1 ]   for |x| ≤ h
//        = 0                                       otherwise
//------------------------------------------------------------------------------
static inline double smoothing_kernel(double d, double h) {
  double ad = fabs(d);
  if (ad > h) return 0.0;
  double x = d / h;
  return (15.0 / (16.0 * h)) * (x*x*x*x - 2.0*x*x + 1.0);
}

//------------------------------------------------------------------------------
// VoF_smoothing:
//   Compute a normalised smoothed phase fraction F̃ by convolution with
//   D ⊗ D ⊗ D.  
//
//   NEW: We accumulate both the weighted sum and the total weight so that
//   F̃ = Σ w F  /  Σ w   → guaranteed 0 ≤ F̃ ≤ 1 even near boundaries.
//------------------------------------------------------------------------------
void VoF_smoothing(Cart3d_bag *data_bag)
{
    MAC_grid       *grid = data_bag->grid;
    VolumeFraction *vof  = data_bag->vof;

    double ***F        = vof->F;
    double ***F_smooth = vof->F_smooth;

    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;
    int Ie = min(grid->G_Ie, grid->NX-1);
    int Je = min(grid->G_Je, grid->NY-1);
    int Ke = min(grid->G_Ke, grid->NZ-1);

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    const double h  = 2.0 * dx;            // support radius (Δx = Δy = Δz)
    const int    r  = (int)ceil(h / dx);    // stencil half‑width in cells

    for (int k = Ks; k < Ke; ++k)
        for (int j = Js; j < Je; ++j)
            for (int i = Is; i < Ie; ++i)
            {
                double sum  = 0.0;          // Σ w F
                double wsum = 0.0;          // Σ w

                for (int kk = k - r; kk <= k + r; ++kk)
                {
                    if (kk < 0 || kk >= grid->NK) continue;
                    double Dz = smoothing_kernel((kk - k) * dz, h);
                    if (Dz == 0.0) continue;

                    for (int jj = j - r; jj <= j + r; ++jj)
                    {
                        if (jj < 0 || jj >= grid->NJ) continue;
                        double Dy = smoothing_kernel((jj - j) * dy, h);
                        if (Dy == 0.0) continue;

                        for (int ii = i - r; ii <= i + r; ++ii)
                        {
                            if (ii < 0 || ii >= grid->NI) continue;
                            double Dx = smoothing_kernel((ii - i) * dx, h);
                            if (Dx == 0.0) continue;

                            double w = Dx * Dy * Dz;
                            sum  += w * F[kk][jj][ii];
                            wsum += w;
                        }
                    }
                }

                // Normalise – fall back to unsmoothed F if wsum → 0 (should not happen)
                F_smooth[k][j][i] = (wsum > 0.0) ? (sum / wsum) : F[k][j][i];
            }

    // Apply periodic / physical BCs to smoothed field
    VOF_set_boundary_values(F_smooth, data_bag);
}





//------------------------------------------------------------------------------
// curvature:
//   1) Smooth F → F̃
//   2) Compute interface normal n̂ = ∇F̃/|∇F̃|
//   3) Compute curvature κ = –∇·n̂
//------------------------------------------------------------------------------
void curvature_patel(Cart3d_bag *data_bag) {
  MAC_grid        *grid = data_bag->grid;
  Parameters      *pr   = data_bag->params;
  VolumeFraction  *vof  = data_bag->vof;

  double       ***F_s   = vof->F_smooth;
  double       ***nx_s    = vof->normal_x_smooth;
  double       ***ny_s    = vof->normal_y_smooth;
  double       ***nz_s    = vof->normal_z_smooth;
  double       ***kappa     = vof->kappa;

  double dx = grid->dx_c[0], dy = grid->dy_c[0], dz = grid->dz_c[0];
  int Is = grid->G_Is;
  int Js = grid->G_Js;
  int Ks = grid->G_Ks;
  int Ie = min(grid->G_Ie, grid->NX-1);
  int Je = min(grid->G_Je, grid->NY-1);
  int Ke = min(grid->G_Ke, grid->NZ-1);

  const double eps_if = 1e-8;   /* interface indicator band */





  // 1) Compute unit normals n̂ = ∇F̃ / |∇F̃|
  for (int k = Ks; k < Ke; ++k) {
      for (int j = Js; j < Je; ++j) {
          for (int i = Is; i < Ie; ++i) {
              // central differences for ∇F̃
              double Fx = (F_s[k][j][i+1] - F_s[k][j][i-1]) / (2.0*dx);
              double Fy = (F_s[k][j+1][i] - F_s[k][j-1][i]) / (2.0*dy);
              double Fz = (F_s[k+1][j][i] - F_s[k-1][j][i]) / (2.0*dz);
              double mag = sqrt(Fx*Fx + Fy*Fy + Fz*Fz) + 1e-12;
              nx_s[k][j][i] = Fx / mag;
              ny_s[k][j][i] = Fy / mag;
              nz_s[k][j][i] = Fz / mag;
          }
      }
  }
    VOF_set_boundary_values(nx_s, data_bag);
    VOF_set_boundary_values(ny_s, data_bag);
    VOF_set_boundary_values(nz_s, data_bag);



  // 2) Compute curvature κ = –(∂n_x/∂x + ∂n_y/∂y + ∂n_z/∂z)
  for (int k = Ks; k < Ke; ++k) {
      for (int j = Js; j < Je; ++j) {
          for (int i = Is; i < Ie; ++i) {
             
              const double Fi = F_s[k][j][i];

              /* Outside the interface band → mark as nodata and skip expensive ops */
              if (Fi <= eps_if || Fi >= 1.0 - eps_if) {
                kappa[k][j][i] = nodata;
                continue;
                }

              double dnx_dx = (nx_s[k][j][i+1] - nx_s[k][j][i-1]) / (2.0*dx);
              double dny_dy = (ny_s[k][j+1][i] - ny_s[k][j-1][i]) / (2.0*dy);
              double dnz_dz = (nz_s[k+1][j][i] - nz_s[k-1][j][i]) / (2.0*dz);
              kappa[k][j][i]   = -(dnx_dx + dny_dy + dnz_dz);
          }
      }
  }
  VOF_set_boundary_values(kappa, data_bag);
          
}

/*=============================================================================
   Adds surface‐tension force to the momentum RHS:
     Fσ = (2/We)·κ·∇F
=============================================================================*/
void Velocity_add_surfacetension_2_RHS_patel(Cart3d_bag *data_bag) {
  MAC_grid       *grid   = data_bag->grid;
  VolumeFraction *vof    = data_bag->vof;
  Parameters     *params = data_bag->params;
  Velocity       *u      = data_bag->u;
  Velocity       *v      = data_bag->v;
  Velocity       *w      = data_bag->w;

  double         ***F      = vof->F;
  double         ***F_s     = vof->F_smooth;
  double         ***kappa  = vof->kappa;
  double         ***rho    = vof->rho;
  const double    scale    = 2.0 / params->We;

  int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;
  int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
  int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;


/* ---------- no‑ghost loop bounds ----------------------------------- */
int i_start_u = max(1,Is) , j_start_u = Js,            k_start_u = Ks;
int i_end_u   = min(NX-1, Ie),  j_end_u   = min(NY-1, Je), k_end_u   = min(NZ-1, Ke);
#ifdef XPERIODIC
	if (i_end_u == NX-1)
		i_end_u = NX;
#endif

int i_start_v = Is,            j_start_v = max(1,Js), k_start_v = Ks;
int i_end_v   = min(NX-1, Ie), j_end_v   = min(NY-1, Je),  k_end_v   = min(NZ-1, Ke);
#ifdef YPERIODIC
    if (j_end_v==NY-1)
        j_end_v = NY;
#endif

int i_start_w = Is,            j_start_w = Js,            k_start_w = max(1,Ks);
int i_end_w   = min(NX-1, Ie), j_end_w   = min(NY-1, Je), k_end_w   = min(NZ-1, Ke);
#ifdef ZPERIODIC
	if (k_end_w == NZ-1)
		k_end_w = NZ;
#endif

  double ***rhsu = u->ng_rhs;
  double ***rhsv = v->ng_rhs;
  double ***rhsw = w->ng_rhs;

  // --- u‐momentum (x‐faces) ---
  for (int k = k_start_u; k < k_end_u; ++k) {
    for (int j = j_start_u; j < j_end_u; ++j) {
      for (int i = i_start_u; i < i_end_u; ++i) {
        
        
        double gradF = (F[k][j][i] - F[k][j][i-1]) * grid->idx_u[i-1];
        if (fabs(gradF) < 1e-12) continue; 
        
        /* cell‑centred curvatures ----------------------------------------- */
        double k1 = kappa[k][j][i  ];        /* right  cell  (i  ,j,k) */
        double k2 = kappa[k][j][i-1];        /* left   cell  (i-1,j,k) */

        /* nodata flags ----------------------------------------------------- */
        bool k1_ok = (k1 < 1e19);            /* H defined? */
        bool k2_ok = (k2 < 1e19);

        if (!k1_ok && !k2_ok)                /* both undefined → skip face   */
            continue;

        /* face‑centred curvature κ_{i+1/2,j} ------------------------------ *
        *   – average when both sides valid                                *
        *   – single‑sided when only one valid (rule from the picture)      */
        double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)   /* both defined */
                                    : (k1_ok ? k1 : k2); /* one defined  */
        
        
        rhsu[k][j][i] += scale * kf * gradF;
      }
    }
  }

  // --- v‐momentum (y‐faces) ---
  for (int k = k_start_v; k < k_end_v; ++k) {
    for (int j = j_start_v; j < j_end_v; ++j) {
      for (int i = i_start_v; i < i_end_v; ++i) {
         
        double gradF = (F[k][j][i] - F[k][j-1][i]) * grid->idy_v[j-1];
        if (fabs(gradF) < 1e-12) continue; 
        
        /* cell‑centred curvatures ----------------------------------------- */
        double k1 = kappa[k][j][i  ];       
        double k2 = kappa[k][j-1][i];        

        /* nodata flags ----------------------------------------------------- */
        bool k1_ok = (k1 < 1e19);            /* H defined? */
        bool k2_ok = (k2 < 1e19);

        if (!k1_ok && !k2_ok)                /* both undefined → skip face   */
            continue;

        /* face‑centred curvature κ_{i+1/2,j} ------------------------------ *
        *   – average when both sides valid                                *
        *   – single‑sided when only one valid (rule from the picture)      */
        double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)   /* both defined */
                                    : (k1_ok ? k1 : k2); /* one defined  */
        
        rhsv[k][j][i] += scale * kf * gradF;
      }
    }
  }

  // --- w‐momentum (z‐faces) ---
  for (int k = k_start_w; k < k_end_w; ++k) {
    for (int j = j_start_w; j < j_end_w; ++j) {
      for (int i = i_start_w; i < i_end_w; ++i) {
        
        double gradF = (F[k][j][i] - F[k-1][j][i]) * grid->idz_w[k-1];
        if (fabs(gradF) < 1e-12) continue; 
        
        /* cell‑centred curvatures ----------------------------------------- */
        double k1 = kappa[k][j][i  ];        
        double k2 = kappa[k-1][j][i];        

        /* nodata flags ----------------------------------------------------- */
        bool k1_ok = (k1 < 1e19);            /* H defined? */
        bool k2_ok = (k2 < 1e19);

        if (!k1_ok && !k2_ok)                /* both undefined → skip face   */
            continue;

        /* face‑centred curvature κ_{i+1/2,j} ------------------------------ *
        *   – average when both sides valid                                *
        *   – single‑sided when only one valid (rule from the picture)      */
        double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)   /* both defined */
                                    : (k1_ok ? k1 : k2); /* one defined  */
        
        rhsw[k][j][i] += scale * kf * gradF;
      }
    }
  }
}



void Velocity_add_gravity_2_RHS(Cart3d_bag *data_bag) {
    MAC_grid       *grid   = data_bag->grid;
    VolumeFraction *vof    = data_bag->vof;
    Parameters     *params = data_bag->params;
    Velocity       *u      = data_bag->u;
    Velocity       *v      = data_bag->v;
    Velocity       *w      = data_bag->w;

    double         ***rho  = vof->rho;
    double         *grav   = params->grav;        // gravity vector
    double         *rich   = params->richardson;  // Richardson number(s)
    double          Ri     = rich[0];             // single Richardson
    double          gm     = sqrt(grav[0]*grav[0]
                               + grav[1]*grav[1]
                               + grav[2]*grav[2]);
    double          gx=0, gy=0, gz=0;
    if (gm > 0.0) {
        gx = grav[0]/gm;
        gy = grav[1]/gm;
        gz = grav[2]/gm;
    }

    int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;
    int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
    int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;

    /* no-ghost loop bounds */
    int i_s_u = max(1,Is),       j_s_u = Js,           k_s_u = Ks;
    int i_e_u = min(NX-1,Ie),    j_e_u = min(NY-1,Je), k_e_u = min(NZ-1,Ke);
    #ifdef XPERIODIC
    if (i_e_u==NX-1) i_e_u = NX;
    #endif

    int i_s_v = Is,              j_s_v = max(1,Js),    k_s_v = Ks;
    int i_e_v = min(NX-1,Ie),    j_e_v = min(NY-1,Je), k_e_v = min(NZ-1,Ke);
    #ifdef YPERIODIC
    if (j_e_v==NY-1) j_e_v = NY;
    #endif

    int i_s_w = Is,              j_s_w = Js,           k_s_w = max(1,Ks);
    int i_e_w = min(NX-1,Ie),    j_e_w = min(NY-1,Je), k_e_w = min(NZ-1,Ke);
    #ifdef ZPERIODIC
    if (k_e_w==NZ-1) k_e_w = NZ;
    #endif

    double ***rhsu = u->ng_rhs;
    double ***rhsv = v->ng_rhs;
    double ***rhsw = w->ng_rhs;

    /* u-momentum: add ρ·Ri·g_x at each x-face */
    for (int k = k_s_u; k < k_e_u; ++k)
    for (int j = j_s_u; j < j_e_u; ++j)
    for (int i = i_s_u; i < i_e_u; ++i) {
        double rho_f = 0.5*(rho[k][j][i] + rho[k][j][i-1]);
        rhsu[k][j][i] += 2.0 * rho_f * Ri * gx;
    }

    /* v-momentum: add ρ·Ri·g_y at each y-face */
    for (int k = k_s_v; k < k_e_v; ++k)
    for (int j = j_s_v; j < j_e_v; ++j)
    for (int i = i_s_v; i < i_e_v; ++i) {
        double rho_f = 0.5*(rho[k][j][i] + rho[k][j-1][i]);
        rhsv[k][j][i] += 2.0 * rho_f * Ri * gy;
    }

    /* w-momentum: add ρ·Ri·g_z at each z-face */
    for (int k = k_s_w; k < k_e_w; ++k)
    for (int j = j_s_w; j < j_e_w; ++j)
    for (int i = i_s_w; i < i_e_w; ++i) {
        double rho_f = 0.5*(rho[k][j][i] + rho[k-1][j][i]);
        rhsw[k][j][i] += 2.0 * rho_f * Ri * gz;
    }
}





/*=============================================================================
   Balanced‑force surface‑tension correction to velocity
=============================================================================*/
void SurfaceTension_project_velocity(Cart3d_bag *data_bag) {

    int i, j, k;

    MAC_grid       *grid   = data_bag->grid;
    VolumeFraction *vof    = data_bag->vof;
    Parameters     *params = data_bag->params;
    Velocity       *u      = data_bag->u;
    Velocity       *v      = data_bag->v;
    Velocity       *w      = data_bag->w;
    Pressure        *p      = data_bag->p;

    double         ***F     = vof->F;
    double         ***kappa = vof->kappa;

    const double BET[] = {BETA};

    // combined scale
    const double scale = (2.0 * BET[params->which_stage] * params->dt) / params->We;

    double         *grav   = params->grav;        // gravity vector
    double         *rich   = params->richardson;  // Richardson number(s)
    double          Ri     = rich[0];             // single Richardson
    double          gm     = sqrt(grav[0]*grav[0]
                               + grav[1]*grav[1]
                               + grav[2]*grav[2]);
    double          gx=0, gy=0, gz=0;
    if (gm > 0.0) {
        gx = grav[0]/gm;
        gy = grav[1]/gm;
        gz = grav[2]/gm;
    }

    
    //Communication_update_ghost_nodes_flow_variable(F, VOLUME_FRACTION, params->ghost_nodes, data_bag);

    //VoF_smoothing(data_bag);

    //curvature_patel(data_bag);

      int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;
      int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
      int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;

    double *idxdt = p->idxdt;
	double *idydt = p->idydt;
	double *idzdt = p->idzdt;  

    for (i=0;i<NX;i++) {
		idxdt[i] = grid->idx_c[i] * scale;
	}

	for (j=0;j<NY;j++) {
		idydt[j] = grid->idy_c[j] * scale;
	}

	for (k=0;k<NZ;k++) {
		idzdt[k] = grid->idz_c[k] * scale;
	}



    /* ---------- no‑ghost loop bounds ----------------------------------- */
    int i_start_u = max(1,Is) , j_start_u = Js,            k_start_u = Ks;
    int i_end_u   = min(NX-1, Ie),  j_end_u   = min(NY-1, Je), k_end_u   = min(NZ-1, Ke);
    #ifdef XPERIODIC
        if (i_end_u == NX-1)
            i_end_u = NX;
    #endif

    int i_start_v = Is,            j_start_v = max(1,Js), k_start_v = Ks;
    int i_end_v   = min(NX-1, Ie), j_end_v   = min(NY-1, Je),  k_end_v   = min(NZ-1, Ke);
    #ifdef YPERIODIC
        if (j_end_v==NY-1)
            j_end_v = NY;
    #endif

    int i_start_w = Is,            j_start_w = Js,            k_start_w = max(1,Ks);
    int i_end_w   = min(NX-1, Ie), j_end_w   = min(NY-1, Je), k_end_w   = min(NZ-1, Ke);
    #ifdef ZPERIODIC
        if (k_end_w == NZ-1)
            k_end_w = NZ;
    #endif


    // X‑face correction
    for (int k = k_start_u; k < k_end_u; ++k) {
        for (int j = j_start_u; j < j_end_u; ++j) {
            for (int i = i_start_u; i < i_end_u; ++i) {
                
                double gradF = (F[k][j][i] - F[k][j][i-1]);
                if (fabs(gradF) < 1e-12) continue; 
                        
                /* cell‑centred curvatures ----------------------------------------- */
                double k1 = kappa[k][j][i  ];        
                double k2 = kappa[k][j][i-1];        

                /* nodata flags ----------------------------------------------------- */
                bool k1_ok = (k1 < 1e19);            /* H defined? */
                bool k2_ok = (k2 < 1e19);

                if (!k1_ok && !k2_ok)                /* both undefined → skip face   */
                    continue;

                /* face‑centred curvature κ_{i+1/2,j} ------------------------------ *
                *   – average when both sides valid                                *
                *   – single‑sided when only one valid (rule from the picture)      */
                double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)   /* both defined */
                                            : (k1_ok ? k1 : k2); /* one defined  */

                double inv_rho = 2.0/(vof->rho[k][j][i]+vof->rho[k][j][i-1]);

                u->data[k][j][i] += inv_rho * kf * gradF * idxdt[i-1];
            }
        }
    }

    // Y‑face correction
    for (int k = k_start_v; k < k_end_v; ++k) {
        for (int j = j_start_v; j < j_end_v; ++j) {
            for (int i = i_start_v; i < i_end_v; ++i) {

                double gradF = (F[k][j][i] - F[k][j-1][i]);
                if (fabs(gradF) < 1e-12) continue; 
                        
                /* cell‑centred curvatures ----------------------------------------- */
                double k1 = kappa[k][j][i  ];        
                double k2 = kappa[k][j-1][i];        

                /* nodata flags ----------------------------------------------------- */
                bool k1_ok = (k1 < 1e19);            /* H defined? */
                bool k2_ok = (k2 < 1e19);

                if (!k1_ok && !k2_ok)                /* both undefined → skip face   */
                    continue;

                /* face‑centred curvature κ_{i+1/2,j} ------------------------------ *
                *   – average when both sides valid                                *
                *   – single‑sided when only one valid (rule from the picture)      */
                double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)   /* both defined */
                                            : (k1_ok ? k1 : k2); /* one defined  */

                double inv_rho = 2.0/(vof->rho[k][j][i]+vof->rho[k][j-1][i]);

                v->data[k][j][i] += inv_rho * kf * gradF * idydt[j-1];
            }
        }
    }

    // Z‑face correction
    for (int k = k_start_w; k < k_end_w; ++k) {
        for (int j = j_start_w; j < j_end_w; ++j) {
            for (int i = i_start_w; i < i_end_w; ++i) {

                double gradF = (F[k][j][i] - F[k-1][j][i]);
                if (fabs(gradF) < 1e-12) continue; 
                        
                /* cell‑centred curvatures ----------------------------------------- */
                double k1 = kappa[k][j][i  ];        
                double k2 = kappa[k-1][j][i];        

                /* nodata flags ----------------------------------------------------- */
                bool k1_ok = (k1 < 1e19);            /* H defined? */
                bool k2_ok = (k2 < 1e19);

                if (!k1_ok && !k2_ok)                /* both undefined → skip face   */
                    continue;

                /* face‑centred curvature κ_{i+1/2,j} ------------------------------ *
                *   – average when both sides valid                                *
                *   – single‑sided when only one valid (rule from the picture)      */
                double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)   /* both defined */
                                            : (k1_ok ? k1 : k2); /* one defined  */

                double inv_rho = 2.0/(vof->rho[k][j][i]+vof->rho[k-1][j][i]);

                w->data[k][j][i] += inv_rho * kf * gradF * idzdt[k-1];
            }
        }
    }


}





#endif // VOF_PLIC