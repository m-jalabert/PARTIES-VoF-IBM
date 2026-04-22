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
 #include <mpi.h>
 
 


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
//   D ⊗ D ⊗ D over the *local array including ghosts*.
//
//   Rationale:
//   • Curvature uses central differences of F̃ and explicitly touches j±1,
//     so F̃ must exist and be *physically consistent* in ghost layers.
//   • We therefore smooth F everywhere on the local domain [L_*s, L_*e),
//     not just global interior [G_*s, G_*e), so the first interior row (j=G_Js)
//     sees symmetric support (ghosts carry the θ-extended content).
//
//   Normalisation:
//   • Accumulate Σ w F and Σ w so that F̃ = (ΣwF)/(Σw) → guaranteed 0≤F̃≤1
//     and robust even near walls/partitions.
//
//   IMPORTANT ordering (already in your time loop):
//   1) VOF_wall_contact_angle_extend()   // encode θ in F ghosts
//   2) VoF_smoothing()                   // now uses extended ghosts
//   3) curvature_patel()                 // central differences on F̃
//------------------------------------------------------------------------------ 
void VoF_smoothing(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    VolumeFraction *vof    = data_bag->vof;
    Parameters     *params = data_bag->params;

    #ifdef VOF_IBM
    double ***F        = vof->F;         // input (already extended into ghosts)
    #else
    double ***F        = vof->F;                  // input
    #endif
    double ***F_smooth = vof->F_smooth;  // output

    #if defined(VOF_WETTING) || defined(VOF_IBM)
    // ----- Use LOCAL extents (include ghosts) so we compute F̃ in ghosts, too -----
    const int IsL = grid->L_Is;
    const int JsL = grid->L_Js;
    const int KsL = grid->L_Ks;
    const int IeL = grid->L_Ie;   // one-past-end
    const int JeL = grid->L_Je;
    const int KeL = grid->L_Ke;
    #else
    // ----- Use GLOBAL interior extents (no ghosts) -----
    const int IsL = grid->G_Is;
    const int JsL = grid->G_Js;
    const int KsL = grid->G_Ks;
    const int IeL = grid->G_Ie;   // one-past-end
    const int JeL = grid->G_Je;
    const int KeL = grid->G_Ke;
    #endif

    // (We still clamp donor indices to [L_*s, L_*e) below.)

    const double dx = grid->dx_c[0];
    const double dy = grid->dy_c[0];
    const double dz = grid->dz_c[0];

    // Support radius and integer stencil half-width (uniform grid assumed here)
    const double h  = 2.0 * dx;
    const int    r  = (int)ceil(h / dx);

    // ----- Smooth EVERYWHERE on local domain (including ghosts) -----
    for (int k = KsL; k < KeL; ++k)
    for (int j = JsL; j < JeL; ++j)
    for (int i = IsL; i < IeL; ++i)
    {
        double sum  = 0.0;   // Σ w F
        double wsum = 0.0;   // Σ w

        // Accumulate donor contributions within radius h in index space.
        // Donors are clamped to LOCAL bounds so we only read locally available data.
        for (int kk = k - r; kk <= k + r; ++kk)
        {
            #if defined(VOF_WETTING) || defined(VOF_IBM)
            if (kk < KsL || kk >= KeL) continue; 
            #endif
            const double Dz = smoothing_kernel((kk - k) * dz, h);
            if (Dz == 0.0) continue;

            for (int jj = j - r; jj <= j + r; ++jj)
            {
                #if defined(VOF_WETTING) || defined(VOF_IBM)
                if (jj < JsL || jj >= JeL) continue;
                #endif
                const double Dy = smoothing_kernel((jj - j) * dy, h);
                if (Dy == 0.0) continue;

                for (int ii = i - r; ii <= i + r; ++ii)
                {
                    #if defined(VOF_WETTING) || defined(VOF_IBM)
                    if (ii < IsL || ii >= IeL) continue;
                    #endif
                    const double Dx = smoothing_kernel((ii - i) * dx, h);
                    if (Dx == 0.0) continue;

                    const double w = Dx * Dy * Dz;
                    sum  += w * F[kk][jj][ii];
                    wsum += w;
                }
            }
        }

        // Normalise – fall back to unsmoothed F if wsum → 0 (should not happen with any overlap)
        F_smooth[k][j][i] = (wsum > 0.0) ? (sum / wsum) : F[k][j][i];


    }

    // Exchange F̃ across MPI subdomains.
    VOF_set_boundary_values(vof->F_smooth, data_bag);
    
}



#define MAX_UNKNOWNS 10

//------------------------------------------------------------------------------
// Linear Algebra Helper: Gaussian Elimination for 3 Right-Hand Sides
// Solves A * [x_vec, y_vec, z_vec] = [bx, by, bz]
// Returns 0 on success, -1 if singular.
//------------------------------------------------------------------------------
static int solve_linear_system_3RHS(int N, double A[MAX_UNKNOWNS][MAX_UNKNOWNS], 
                                    double bx[MAX_UNKNOWNS], double by[MAX_UNKNOWNS], double bz[MAX_UNKNOWNS],
                                    double x_out[MAX_UNKNOWNS], double y_out[MAX_UNKNOWNS], double z_out[MAX_UNKNOWNS]) 
{
    int i, j, k, max_row;
    double factor, tmp;

    // Forward Elimination
    for (i = 0; i < N; i++) {
        // Pivot
        max_row = i;
        for (k = i + 1; k < N; k++) {
            if (fabs(A[k][i]) > fabs(A[max_row][i])) max_row = k;
        }

        // Swap rows in A and all B vectors
        if (max_row != i) {
            for (k = i; k < N; k++) {
                tmp = A[max_row][k]; A[max_row][k] = A[i][k]; A[i][k] = tmp;
            }
            tmp = bx[max_row]; bx[max_row] = bx[i]; bx[i] = tmp;
            tmp = by[max_row]; by[max_row] = by[i]; by[i] = tmp;
            tmp = bz[max_row]; bz[max_row] = bz[i]; bz[i] = tmp;
        }

        // Singular check
        if (fabs(A[i][i]) < 1e-25) return -1; 

        // Eliminate
        double inv_diag = 1.0 / A[i][i];
        for (k = i + 1; k < N; k++) {
            factor = A[k][i] * inv_diag;
            for (j = i; j < N; j++) A[k][j] -= factor * A[i][j];
            bx[k] -= factor * bx[i];
            by[k] -= factor * by[i];
            bz[k] -= factor * bz[i];
        }
    }

    // Back Substitution for all 3 vectors
    for (i = N - 1; i >= 0; i--) {
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
        for (j = i + 1; j < N; j++) {
            sum_x += A[i][j] * x_out[j];
            sum_y += A[i][j] * y_out[j];
            sum_z += A[i][j] * z_out[j];
        }
        double inv_diag = 1.0 / A[i][i];
        x_out[i] = (bx[i] - sum_x) * inv_diag;
        y_out[i] = (by[i] - sum_y) * inv_diag;
        z_out[i] = (bz[i] - sum_z) * inv_diag;
    }
    return 0;
}

//------------------------------------------------------------------------------
// Weighted Least Squares (WLS) for Divergence (Curvature)
// Solves partial derivatives for nx, ny, and nz simultaneously.
//
// Inputs:
//   nx_in, ny_in, nz_in: The vector field arrays
//   vfc: Solid volume fraction (used for masking and cut-cell detection)
//   is_curvature_mode: If true, enables O'Brien's boundary injection logic
// Outputs:
//   div_out: The divergence result (-div(n))
//------------------------------------------------------------------------------
static void compute_WLS_divergence_3D(int i, int j, int k, 
                                      double ***nx_in, double ***ny_in, double ***nz_in,
                                      double dx, double dy, double dz,
                                      double *div_out) 
{
    const int N_UNKNOWNS = 10; // Quadratic (1, x, y, z, x2, y2, z2, xy, yz, xz)
    double A[10][10];
    double bx[10], by[10], bz[10];
    double sol_x[10], sol_y[10], sol_z[10];
    
    // Clear systems
    memset(A, 0, sizeof(A));
    memset(bx, 0, sizeof(bx));
    memset(by, 0, sizeof(by));
    memset(bz, 0, sizeof(bz));

    int stencil_count = 0;

    // Iterate 3x3x3 Stencil
    for (int kk = -1; kk <= 1; kk++) {
        for (int jj = -1; jj <= 1; jj++) {
            for (int ii = -1; ii <= 1; ii++) {
                int ni = i + ii;
                int nj = j + jj;
                int nk = k + kk;


                double val_x, val_y, val_z;

                val_x = nx_in[nk][nj][ni];
                val_y = ny_in[nk][nj][ni];
                val_z = nz_in[nk][nj][ni];
                
              

                // --- 2. Gaussian Weighting ---
                // r^2 in index space. Center (0,0,0) has dist=0.
                double dist_sq = (double)(ii*ii + jj*jj + kk*kk);
                double weight = exp(-dist_sq / 2.0); // Gaussian decay sigma=1.41

                // --- 3. Basis Construction ---
                double lx = ii * dx;
                double ly = jj * dy;
                double lz = kk * dz;

                // Basis: [1, x, y, z, x^2, y^2, z^2, xy, yz, xz]
                double basis[10];
                basis[0] = 1.0;
                basis[1] = lx; basis[2] = ly; basis[3] = lz;
                basis[4] = lx*lx; basis[5] = ly*ly; basis[6] = lz*lz;
                basis[7] = lx*ly; basis[8] = ly*lz; basis[9] = lx*lz;

                // Fill Normal Equations: A^T * W * A
                for (int r = 0; r < N_UNKNOWNS; r++) {
                    double w_br = weight * basis[r];
                    for (int c = 0; c < N_UNKNOWNS; c++) {
                        A[r][c] += w_br * basis[c];
                    }
                    // Fill RHS: A^T * W * b
                    bx[r] += w_br * val_x;
                    by[r] += w_br * val_y;
                    bz[r] += w_br * val_z;
                }
                stencil_count++;
            }
        }
    }

    // --- 4. Solve and Extract Derivatives ---
    // We need at least 10 points for Quadratic. Fallback to FD if ill-conditioned.
    if (stencil_count < N_UNKNOWNS || solve_linear_system_3RHS(N_UNKNOWNS, A, bx, by, bz, sol_x, sol_y, sol_z) != 0) {
        // Fallback: Standard Central Difference
        // Note: This relies on neighbors being valid. If solid neighbors were skipped above,
        // this might read garbage, but usually we only fallback in narrow gaps.
        // A safer fallback would be 0.0 or one-sided diff, but let's stick to simple CD.
        double dnx_dx = (nx_in[k][j][i+1] - nx_in[k][j][i-1]) / (2.0*dx);
        double dny_dy = (ny_in[k][j+1][i] - ny_in[k][j-1][i]) / (2.0*dy);
        double dnz_dz = (nz_in[k+1][j][i] - nz_in[k-1][j][i]) / (2.0*dz);
        *div_out = -(dnx_dx + dny_dy + dnz_dz);
    } else {
        // Solution vector x = [c, dx, dy, dz, ...]
        // dnx/dx is index 1 of sol_x
        // dny/dy is index 2 of sol_y
        // dnz/dz is index 3 of sol_z
        *div_out = -(sol_x[1] + sol_y[2] + sol_z[3]);
    }
}

//------------------------------------------------------------------------------
// Helper: Scalar Gradient via WLS (for calculating normal from F_smooth)
// Reuses the logic but for a single scalar field.
//------------------------------------------------------------------------------
static void compute_WLS_gradient_scalar(int i, int j, int k, 
                                        double ***phi,
                                        double dx, double dy, double dz,
                                        double grad[3]) 
{
    const int N_UNKNOWNS = 4; // Linear is sufficient for normals (1, x, y, z)
    double A[10][10] = {0};   // Use slightly larger buffer to reuse solver signature
    double b[10] = {0}, dum[10] = {0};
    double x[10] = {0}, dum_out[10] = {0};
    int stencil_count = 0;

    for (int kk = -1; kk <= 1; kk++) {
        for (int jj = -1; jj <= 1; jj++) {
            for (int ii = -1; ii <= 1; ii++) {
                int ni = i + ii; int nj = j + jj; int nk = k + kk;
                //if (ni < IsL || ni >= IeL || nj < JsL || nj >= JeL || nk < KsL || nk >= KeL) continue;

                

                double weight = exp(-(double)(ii*ii + jj*jj + kk*kk) / 2.0);
                double lx = ii * dx, ly = jj * dy, lz = kk * dz;
                double basis[4] = {1.0, lx, ly, lz};
                double val = phi[nk][nj][ni];

                for (int r = 0; r < N_UNKNOWNS; r++) {
                    double w_br = weight * basis[r];
                    for (int c = 0; c < N_UNKNOWNS; c++) A[r][c] += w_br * basis[c];
                    b[r] += w_br * val;
                }
                stencil_count++;
            }
        }
    }

    if (stencil_count < N_UNKNOWNS || solve_linear_system_3RHS(N_UNKNOWNS, A, b, dum, dum, x, dum_out, dum_out) != 0) {
        grad[0] = (phi[k][j][i+1] - phi[k][j][i-1]) / (2.0*dx);
        grad[1] = (phi[k][j+1][i] - phi[k][j-1][i]) / (2.0*dy);
        grad[2] = (phi[k+1][j][i] - phi[k-1][j][i]) / (2.0*dz);
    } else {
        grad[0] = x[1]; grad[1] = x[2]; grad[2] = x[3];
    }
}


/*=============================================================================
 *  CELESTE 3D Curvature Method
 *  
 *  Extension of O'Brien's 2D CELESTE scheme to 3D:
 *    Step 1: Interface normals from quadratic difference-form WLS on F_smooth
 *            with 1/Δs² weighting
 *    Step 2: Curvature κ = -∇·n̂ from UNWEIGHTED quadratic difference-form WLS
 *            on the unit normal components, solved simultaneously for nx,ny,nz
 *
 *  Key differences from the previous implementation:
 *    - Difference form (no constant term): δF = F_nb - F_c, not F_nb
 *    - Quadratic basis for BOTH steps (9 unknowns in 3D)
 *    - 1/Δs² weighting for normals (not Gaussian)
 *    - UNWEIGHTED for curvature (not Gaussian)
 *    - On uniform grids the curvature matrix A is constant → precomputable
 *
 *  References:
 *    O'Brien, A.R. (2019) PhD Thesis, Section 3.4
 *    Denner & van Wachem, "CELESTE" scheme
 *=============================================================================*/


/* Number of unknowns for the difference-form quadratic basis in 3D:
 *   (Δx, Δy, Δz, Δx²/2, Δy²/2, Δz²/2, ΔxΔy, ΔyΔz, ΔxΔz)        */
#define CELESTE_N  9


/*=============================================================================
 *  Build the 9-component quadratic difference-form basis vector
 *
 *  basis · q  =  F_x Δx + F_y Δy + F_z Δz
 *              + F_xx Δx²/2 + F_yy Δy²/2 + F_zz Δz²/2
 *              + F_xy ΔxΔy  + F_yz ΔyΔz  + F_xz ΔxΔz
 *
 *  So the solution vector q maps as:
 *    q[0] = ∂/∂x    q[1] = ∂/∂y    q[2] = ∂/∂z
 *    q[3] = ∂²/∂x²  q[4] = ∂²/∂y²  q[5] = ∂²/∂z²
 *    q[6] = ∂²/∂x∂y q[7] = ∂²/∂y∂z q[8] = ∂²/∂x∂z
 *=============================================================================*/
static inline void celeste_basis(double lx, double ly, double lz,
                                 double basis[CELESTE_N])
{
    basis[0] = lx;
    basis[1] = ly;
    basis[2] = lz;
    basis[3] = 0.5 * lx * lx;
    basis[4] = 0.5 * ly * ly;
    basis[5] = 0.5 * lz * lz;
    basis[6] = lx * ly;
    basis[7] = ly * lz;
    basis[8] = lx * lz;
}




/*=============================================================================
 *  STEP 1:  Compute interface normal at cell (i,j,k) from F_smooth
 *
 *  Method:
 *    Quadratic difference-form WLS with 1/Δs² weighting on a 3×3×3 stencil
 *    (26 neighbours, excluding center).
 *
 *    The difference form  δF = F_nb − F_c  with no constant term forces the
 *    polynomial to pass exactly through the center point, which improves
 *    gradient accuracy compared to including a constant.
 *
 *
 *  Output:
 *    n_out[3] = unit normal  ∇F̃ / |∇F̃|
 *    Returns |∇F̃| (useful for interface detection / masking)
 *=============================================================================*/
double compute_CELESTE_normal_3D(int i, int j, int k,
                                 double ***F_s,
                                 double dx, double dy, double dz,
                                 double n_out[3])
{
    const int N = CELESTE_N;
    double A[MAX_UNKNOWNS][MAX_UNKNOWNS];
    double b[MAX_UNKNOWNS], dum1[MAX_UNKNOWNS], dum2[MAX_UNKNOWNS];
    double sol[MAX_UNKNOWNS], dum_o1[MAX_UNKNOWNS], dum_o2[MAX_UNKNOWNS];

    memset(A,    0, sizeof(double) * MAX_UNKNOWNS * MAX_UNKNOWNS);
    memset(b,    0, sizeof(double) * MAX_UNKNOWNS);
    memset(dum1, 0, sizeof(double) * MAX_UNKNOWNS);
    memset(dum2, 0, sizeof(double) * MAX_UNKNOWNS);

    const double Fc = F_s[k][j][i];
    int count = 0;

    for (int kk = -1; kk <= 1; kk++)
    for (int jj = -1; jj <= 1; jj++)
    for (int ii = -1; ii <= 1; ii++)
    {
        if (ii == 0 && jj == 0 && kk == 0) continue;   /* skip centre */

        double lx  = ii * dx;
        double ly  = jj * dy;
        double lz  = kk * dz;
        double ds2 = lx*lx + ly*ly + lz*lz;             /* Δs² */
        double w = 1.0 / (ds2 * ds2);                          /* 1/Δs^4 weight */

        double basis[CELESTE_N];
        celeste_basis(lx, ly, lz, basis);

        double dF = F_s[k+kk][j+jj][i+ii] - Fc;         /* difference form */

        for (int r = 0; r < N; r++) {
            double w_br = w * basis[r];
            for (int c = 0; c < N; c++)
                A[r][c] += w_br * basis[c];
            b[r] += w_br * dF;
        }
        count++;
    }

    /* Solve  A q = b  (use 3-RHS solver with dummy vectors) */
    if (count < N ||
        solve_linear_system_3RHS(N, A, b, dum1, dum2,
                                 sol, dum_o1, dum_o2) != 0)
    {
        /* Fallback: central differences */
        sol[0] = (F_s[k][j][i+1] - F_s[k][j][i-1]) / (2.0 * dx);
        sol[1] = (F_s[k][j+1][i] - F_s[k][j-1][i]) / (2.0 * dy);
        sol[2] = (F_s[k+1][j][i] - F_s[k-1][j][i]) / (2.0 * dz);
    }

    /* sol[0..2] = (F̃_x, F̃_y, F̃_z) */
    double mag = sqrt(sol[0]*sol[0] + sol[1]*sol[1] + sol[2]*sol[2]);
    if (mag < 1e-30) {
        n_out[0] = n_out[1] = n_out[2] = 0.0;
        return 0.0;
    }
    n_out[0] = sol[0] / mag;
    n_out[1] = sol[1] / mag;
    n_out[2] = sol[2] / mag;
    return mag;
}


/*=============================================================================
 *  STEP 2:  Compute curvature  κ = −∇·n̂  at cell (i,j,k)
 *
 *  Method:
 *    UNWEIGHTED quadratic difference-form WLS on the unit normal field,
 *    solving all three components (nx, ny, nz) simultaneously.
 *
 *    "More accurate curvatures are achieved when using an unweighted
 *     least-squares procedure for the computation of interface normal
 *     divergence."  — O'Brien (2019), Section 3.4.2
 *
 *    The quadratic terms in the basis absorb second-order truncation error
 *    even though only the linear coefficients enter κ.
 *
 *  NOTE: On a uniform Cartesian grid with a symmetric 3×3×3 stencil and
 *  unit weights, the matrix A = Σ φφᵀ is identical for every cell.
 *  See compute_CELESTE_curvature_3D_precomp() below for the optimised
 *  version that precomputes and LU-factors A once.
 *=============================================================================*/
void compute_CELESTE_curvature_3D(int i, int j, int k,
                                  double ***nx_s, double ***ny_s, double ***nz_s,
                                  double dx, double dy, double dz,
                                  double *kappa_out)
{
    const int N = CELESTE_N;
    double A[MAX_UNKNOWNS][MAX_UNKNOWNS];
    double bx[MAX_UNKNOWNS], by[MAX_UNKNOWNS], bz[MAX_UNKNOWNS];
    double sx[MAX_UNKNOWNS], sy[MAX_UNKNOWNS], sz[MAX_UNKNOWNS];

    memset(A,  0, sizeof(double) * MAX_UNKNOWNS * MAX_UNKNOWNS);
    memset(bx, 0, sizeof(double) * MAX_UNKNOWNS);
    memset(by, 0, sizeof(double) * MAX_UNKNOWNS);
    memset(bz, 0, sizeof(double) * MAX_UNKNOWNS);

    const double nxc = nx_s[k][j][i];
    const double nyc = ny_s[k][j][i];
    const double nzc = nz_s[k][j][i];
    int count = 0;

    for (int kk = -1; kk <= 1; kk++)
    for (int jj = -1; jj <= 1; jj++)
    for (int ii = -1; ii <= 1; ii++)
    {
        if (ii == 0 && jj == 0 && kk == 0) continue;

        double lx = ii * dx;
        double ly = jj * dy;
        double lz = kk * dz;

        double basis[CELESTE_N];
        celeste_basis(lx, ly, lz, basis);

        /* Difference form: δn = n_nb − n_c */
        double dnx = nx_s[k+kk][j+jj][i+ii] - nxc;
        double dny = ny_s[k+kk][j+jj][i+ii] - nyc;
        double dnz = nz_s[k+kk][j+jj][i+ii] - nzc;

        /* Unweighted:  w = 1 */
        for (int r = 0; r < N; r++) {
            for (int c = 0; c < N; c++)
                A[r][c] += basis[r] * basis[c];
            bx[r] += basis[r] * dnx;
            by[r] += basis[r] * dny;
            bz[r] += basis[r] * dnz;
        }
        count++;
    }

    /* Solve  A [qx qy qz] = [bx by bz]  */
    if (count < N ||
        solve_linear_system_3RHS(N, A, bx, by, bz, sx, sy, sz) != 0)
    {
        /* Fallback: central differences */
        double dnx_dx = (nx_s[k][j][i+1] - nx_s[k][j][i-1]) / (2.0 * dx);
        double dny_dy = (ny_s[k][j+1][i] - ny_s[k][j-1][i]) / (2.0 * dy);
        double dnz_dz = (nz_s[k+1][j][i] - nz_s[k-1][j][i]) / (2.0 * dz);
        *kappa_out = -(dnx_dx + dny_dy + dnz_dz);
        return;
    }

    /*  sx[0] = ∂nx/∂x,  sy[1] = ∂ny/∂y,  sz[2] = ∂nz/∂z  */
    *kappa_out = -(sx[0] + sy[1] + sz[2]);
}


/*=============================================================================
 *  OPTIMISED VERSION: Precomputed matrix for the curvature step
 *
 *  On a uniform Cartesian grid, the unweighted matrix  A = Σ φφᵀ  over the
 *  26-neighbour stencil depends only on (dx, dy, dz) and is the SAME for
 *  every cell.  We precompute its LU factorisation once and then do
 *  back-substitution per cell, avoiding 26×9×9 assembly + pivot search.
 *
 *  Usage:
 *    1) At initialisation:  CELESTE_precompute_curvature_matrix(dx,dy,dz,&ctx);
 *    2) Per cell:           CELESTE_curvature_from_precomp(&ctx, i,j,k, ...);
 *=============================================================================*/

typedef struct {
    double A_LU[CELESTE_N][CELESTE_N];   /* LU-factored matrix           */
    int    piv[CELESTE_N];               /* pivot indices                 */
    double dx, dy, dz;                   /* grid spacing (for basis)      */
    int    valid;                        /* 0 if factorisation failed     */
} CELESTE_CurvatureCtx;


/* ---- LU factorisation with partial pivoting (in-place, 9×9) ---- */
static int lu_factorise_9(double A[CELESTE_N][CELESTE_N], int piv[CELESTE_N])
{
    const int N = CELESTE_N;
    for (int i = 0; i < N; i++) piv[i] = i;

    for (int i = 0; i < N; i++) {
        /* partial pivot */
        int max_row = i;
        for (int k = i + 1; k < N; k++)
            if (fabs(A[k][i]) > fabs(A[max_row][i])) max_row = k;

        if (max_row != i) {
            int tmp_p = piv[i]; piv[i] = piv[max_row]; piv[max_row] = tmp_p;
            for (int c = 0; c < N; c++) {
                double tmp = A[i][c]; A[i][c] = A[max_row][c]; A[max_row][c] = tmp;
            }
        }
        if (fabs(A[i][i]) < 1e-25) return -1;   /* singular */

        double inv = 1.0 / A[i][i];
        for (int k = i + 1; k < N; k++) {
            A[k][i] *= inv;
            for (int j = i + 1; j < N; j++)
                A[k][j] -= A[k][i] * A[i][j];
        }
    }
    return 0;
}


/* ---- Solve  LU x = Pb  using pre-factored LU and pivots ---- */
static void lu_solve_9(const double LU[CELESTE_N][CELESTE_N],
                        const int piv[CELESTE_N],
                        const double b_in[CELESTE_N],
                        double x_out[CELESTE_N])
{
    const int N = CELESTE_N;
    double y[CELESTE_N];

    /* Apply permutation and forward substitution (L y = P b) */
    for (int i = 0; i < N; i++) {
        double sum = b_in[piv[i]];
        for (int j = 0; j < i; j++)
            sum -= LU[i][j] * y[j];
        y[i] = sum;
    }
    /* Back substitution (U x = y) */
    for (int i = N - 1; i >= 0; i--) {
        double sum = y[i];
        for (int j = i + 1; j < N; j++)
            sum -= LU[i][j] * x_out[j];
        x_out[i] = sum / LU[i][i];
    }
}


/* ---- Precompute the unweighted curvature matrix ---- */
void CELESTE_precompute_curvature_matrix(double dx, double dy, double dz,
                                         CELESTE_CurvatureCtx *ctx)
{
    const int N = CELESTE_N;
    ctx->dx = dx;  ctx->dy = dy;  ctx->dz = dz;

    memset(ctx->A_LU, 0, sizeof(ctx->A_LU));

    for (int kk = -1; kk <= 1; kk++)
    for (int jj = -1; jj <= 1; jj++)
    for (int ii = -1; ii <= 1; ii++)
    {
        if (ii == 0 && jj == 0 && kk == 0) continue;

        double basis[CELESTE_N];
        celeste_basis(ii * dx, jj * dy, kk * dz, basis);

        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                ctx->A_LU[r][c] += basis[r] * basis[c];
    }

    ctx->valid = (lu_factorise_9(ctx->A_LU, ctx->piv) == 0) ? 1 : 0;
}


/* ---- Compute curvature at one cell using precomputed LU ---- */
void CELESTE_curvature_from_precomp(const CELESTE_CurvatureCtx *ctx,
                                    int i, int j, int k,
                                    double ***nx_s, double ***ny_s, double ***nz_s,
                                    double *kappa_out)
{
    const int N = CELESTE_N;
    const double dx = ctx->dx, dy = ctx->dy, dz = ctx->dz;

    if (!ctx->valid) {
        /* Matrix was singular — fall back to central differences */
        double dnx_dx = (nx_s[k][j][i+1] - nx_s[k][j][i-1]) / (2.0 * dx);
        double dny_dy = (ny_s[k][j+1][i] - ny_s[k][j-1][i]) / (2.0 * dy);
        double dnz_dz = (nz_s[k+1][j][i] - nz_s[k-1][j][i]) / (2.0 * dz);
        *kappa_out = -(dnx_dx + dny_dy + dnz_dz);
        return;
    }

    /* Assemble 3 RHS vectors */
    double bx[CELESTE_N] = {0}, by[CELESTE_N] = {0}, bz[CELESTE_N] = {0};
    const double nxc = nx_s[k][j][i];
    const double nyc = ny_s[k][j][i];
    const double nzc = nz_s[k][j][i];

    for (int kk = -1; kk <= 1; kk++)
    for (int jj = -1; jj <= 1; jj++)
    for (int ii = -1; ii <= 1; ii++)
    {
        if (ii == 0 && jj == 0 && kk == 0) continue;

        double basis[CELESTE_N];
        celeste_basis(ii * dx, jj * dy, kk * dz, basis);

        double dnx = nx_s[k+kk][j+jj][i+ii] - nxc;
        double dny = ny_s[k+kk][j+jj][i+ii] - nyc;
        double dnz = nz_s[k+kk][j+jj][i+ii] - nzc;

        for (int r = 0; r < N; r++) {
            bx[r] += basis[r] * dnx;
            by[r] += basis[r] * dny;
            bz[r] += basis[r] * dnz;
        }
    }

    /* Solve using pre-factored LU (one back-sub per RHS) */
    double sx[CELESTE_N], sy[CELESTE_N], sz[CELESTE_N];
    lu_solve_9(ctx->A_LU, ctx->piv, bx, sx);
    lu_solve_9(ctx->A_LU, ctx->piv, by, sy);
    lu_solve_9(ctx->A_LU, ctx->piv, bz, sz);

    *kappa_out = -(sx[0] + sy[1] + sz[2]);
}


/*=============================================================================
 *  Similarly precompute the NORMAL matrix (1/Δs²-weighted, also constant
 *  on a uniform grid).
 *=============================================================================*/
typedef struct {
    double A_LU[CELESTE_N][CELESTE_N];
    int    piv[CELESTE_N];
    double dx, dy, dz;
    int    valid;
} CELESTE_NormalCtx;


void CELESTE_precompute_normal_matrix(double dx, double dy, double dz,
                                      CELESTE_NormalCtx *ctx)
{
    const int N = CELESTE_N;
    ctx->dx = dx;  ctx->dy = dy;  ctx->dz = dz;

    memset(ctx->A_LU, 0, sizeof(ctx->A_LU));

    for (int kk = -1; kk <= 1; kk++)
    for (int jj = -1; jj <= 1; jj++)
    for (int ii = -1; ii <= 1; ii++)
    {
        if (ii == 0 && jj == 0 && kk == 0) continue;

        double lx  = ii * dx, ly = jj * dy, lz = kk * dz;
        double ds2 = lx*lx + ly*ly + lz*lz;
        double w = 1.0 / (ds2 * ds2);

        double basis[CELESTE_N];
        celeste_basis(lx, ly, lz, basis);

        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                ctx->A_LU[r][c] += w * basis[r] * basis[c];
    }

    ctx->valid = (lu_factorise_9(ctx->A_LU, ctx->piv) == 0) ? 1 : 0;
}


/* ---- Compute normal at one cell using precomputed LU ---- */
double CELESTE_normal_from_precomp(const CELESTE_NormalCtx *ctx,
                                   int i, int j, int k,
                                   double ***F_s,
                                   double n_out[3])
{
    const int N = CELESTE_N;
    const double dx = ctx->dx, dy = ctx->dy, dz = ctx->dz;

    if (!ctx->valid) {
        double gx = (F_s[k][j][i+1] - F_s[k][j][i-1]) / (2.0 * dx);
        double gy = (F_s[k][j+1][i] - F_s[k][j-1][i]) / (2.0 * dy);
        double gz = (F_s[k+1][j][i] - F_s[k-1][j][i]) / (2.0 * dz);
        double mag = sqrt(gx*gx + gy*gy + gz*gz);
        if (mag < 1e-30) { n_out[0]=n_out[1]=n_out[2]=0; return 0; }
        n_out[0]=gx/mag; n_out[1]=gy/mag; n_out[2]=gz/mag;
        return mag;
    }

    /* Assemble weighted RHS */
    double b[CELESTE_N] = {0};
    const double Fc = F_s[k][j][i];

    for (int kk = -1; kk <= 1; kk++)
    for (int jj = -1; jj <= 1; jj++)
    for (int ii = -1; ii <= 1; ii++)
    {
        if (ii == 0 && jj == 0 && kk == 0) continue;

        double lx  = ii * dx, ly = jj * dy, lz = kk * dz;
        double ds2 = lx*lx + ly*ly + lz*lz;
        double w = 1.0 / (ds2 * ds2);

        double basis[CELESTE_N];
        celeste_basis(lx, ly, lz, basis);

        double dF = F_s[k+kk][j+jj][i+ii] - Fc;

        for (int r = 0; r < N; r++)
            b[r] += w * basis[r] * dF;
    }

    double sol[CELESTE_N];
    lu_solve_9(ctx->A_LU, ctx->piv, b, sol);

    double mag = sqrt(sol[0]*sol[0] + sol[1]*sol[1] + sol[2]*sol[2]);
    if (mag < 1e-30) { n_out[0]=n_out[1]=n_out[2]=0; return 0; }
    n_out[0] = sol[0] / mag;
    n_out[1] = sol[1] / mag;
    n_out[2] = sol[2] / mag;
    return mag;
}

/*============================================================================
 *  VOF_set_boundary_values_normal_vector
 *
 *  Fill physical-boundary ghost cells for a three-component VECTOR field
 *  (nx, ny, nz), then call the MPI ghost exchange.
 *
 *  At FREE-SLIP (symmetry) walls the vector is reflected:
 *    - component normal to the wall  → antisymmetric (sign flip)
 *    - components tangential to wall → symmetric     (same sign)
 *
 *  At all other walls: zero-gradient for every component.
 *
 *  Cascading fill order X → Y → Z (same as VOF_set_boundary_values) so
 *  edge and corner ghost cells are correctly populated.
 *============================================================================*/
void VOF_set_boundary_values_normal_vector(double ***nx,
                                           double ***ny,
                                           double ***nz,
                                           Cart3d_bag *data_bag)
{
    MAC_grid   *grid   = data_bag->grid;
    Parameters *params = data_bag->params;

    const int NX = grid->NX;
    const int NY = grid->NY;
    const int NZ = grid->NZ;
    const int NG = params->ghost_nodes;

    const int Is = grid->G_Is, Ie = grid->G_Ie;
    const int Js = grid->G_Js, Je = grid->G_Je;
    const int Ks = grid->G_Ks, Ke = grid->G_Ke;

    const int IsL = grid->L_Is, IeL = grid->L_Ie;
    const int JsL = grid->L_Js, JeL = grid->L_Je;
    const int KsL = grid->L_Ks, KeL = grid->L_Ke;

    /*=====================================================================
     *  1)  X boundaries   (use local j,k ranges)
     *=====================================================================*/
#ifndef XPERIODIC

    /* --- Left boundary (x = 0) --- */
    if (Is == 0) {
        int i0 = 0;

      #ifdef LEFT_WALL_VELOCITY_FREESLIP
        /*  Symmetry plane:  x-component antisymmetric, y/z symmetric
         *
         *    ghost  i0 - gl   ←   mirror interior  i0 + gl - 1
         *      gl=1: ghost −1  ←  cell 0
         *      gl=2: ghost −2  ←  cell 1
         *      gl=3: ghost −3  ←  cell 2                            */
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int gl = 1; gl <= NG; gl++) {
            nx[k][j][i0 - gl] = -nx[k][j][i0 + gl - 1];   /* antisymmetric */
            ny[k][j][i0 - gl] =  ny[k][j][i0 + gl - 1];   /* symmetric     */
            nz[k][j][i0 - gl] =  nz[k][j][i0 + gl - 1];   /* symmetric     */
        }

      #else
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int gl = 1; gl <= NG; gl++) {
            nx[k][j][i0 - gl] = nx[k][j][i0];
            ny[k][j][i0 - gl] = ny[k][j][i0];
            nz[k][j][i0 - gl] = nz[k][j][i0];
        }
      #endif
    }

    /* --- Right boundary (x = NX-1) --- */
    if (Ie == NX) {
        int i1 = NX - 1;

      #ifdef RIGHT_WALL_VELOCITY_FREESLIP
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int gl = 0; gl < NG; gl++) {
            nx[k][j][i1 + gl] = -nx[k][j][i1 - 1 - gl];
            ny[k][j][i1 + gl] =  ny[k][j][i1 - 1 - gl];
            nz[k][j][i1 + gl] =  nz[k][j][i1 - 1 - gl];
        }

      #else
        for (int k = KsL; k < KeL; k++)
        for (int j = JsL; j < JeL; j++)
        for (int gl = 0; gl < NG; gl++) {
            nx[k][j][i1 + gl] = nx[k][j][i1 - 1];
            ny[k][j][i1 + gl] = ny[k][j][i1 - 1];
            nz[k][j][i1 + gl] = nz[k][j][i1 - 1];
        }
      #endif
    }
#endif /* !XPERIODIC */


    /*=====================================================================
     *  2)  Y boundaries   (use local i,k ranges — x-ghosts already filled)
     *=====================================================================*/
#ifndef YPERIODIC

    /* --- Bottom boundary (y = 0) --- */
    if (Js == 0) {
        int j0 = 0;

      #ifdef BOTTOM_WALL_VELOCITY_FREESLIP
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 1; gl <= NG; gl++) {
            nx[k][j0 - gl][i] =  nx[k][j0 + gl - 1][i];
            ny[k][j0 - gl][i] = -ny[k][j0 + gl - 1][i];   /* y is the wall-normal */
            nz[k][j0 - gl][i] =  nz[k][j0 + gl - 1][i];
        }

      #else
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 1; gl <= NG; gl++) {
            nx[k][j0 - gl][i] = nx[k][j0][i];
            ny[k][j0 - gl][i] = ny[k][j0][i];
            nz[k][j0 - gl][i] = nz[k][j0][i];
        }
      #endif
    }

    /* --- Top boundary (y = NY-1) --- */
    if (Je == NY) {
        int j1 = NY - 1;

      #ifdef TOP_WALL_VELOCITY_FREESLIP
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 0; gl < NG; gl++) {
            nx[k][j1 + gl][i] =  nx[k][j1 - 1 - gl][i];
            ny[k][j1 + gl][i] = -ny[k][j1 - 1 - gl][i];
            nz[k][j1 + gl][i] =  nz[k][j1 - 1 - gl][i];
        }

      #else
        for (int k = KsL; k < KeL; k++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 0; gl < NG; gl++) {
            nx[k][j1 + gl][i] = nx[k][j1 - 1][i];
            ny[k][j1 + gl][i] = ny[k][j1 - 1][i];
            nz[k][j1 + gl][i] = nz[k][j1 - 1][i];
        }
      #endif
    }
#endif /* !YPERIODIC */


    /*=====================================================================
     *  3)  Z boundaries   (use local i,j ranges — x/y ghosts already filled)
     *=====================================================================*/
#ifndef ZPERIODIC

    /* --- Back boundary (z = 0) --- */
    if (Ks == 0) {
        int k0 = 0;

      #ifdef BACK_WALL_VELOCITY_FREESLIP
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 1; gl <= NG; gl++) {
            nx[k0 - gl][j][i] =  nx[k0 + gl - 1][j][i];
            ny[k0 - gl][j][i] =  ny[k0 + gl - 1][j][i];
            nz[k0 - gl][j][i] = -nz[k0 + gl - 1][j][i];   /* z is the wall-normal */
        }

      #else
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 1; gl <= NG; gl++) {
            nx[k0 - gl][j][i] = nx[k0][j][i];
            ny[k0 - gl][j][i] = ny[k0][j][i];
            nz[k0 - gl][j][i] = nz[k0][j][i];
        }
      #endif
    }

    /* --- Front boundary (z = NZ-1) --- */
    if (Ke == NZ) {
        int k1 = NZ - 1;

      #ifdef FRONT_WALL_VELOCITY_FREESLIP
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 0; gl < NG; gl++) {
            nx[k1 + gl][j][i] =  nx[k1 - 1 - gl][j][i];
            ny[k1 + gl][j][i] =  ny[k1 - 1 - gl][j][i];
            nz[k1 + gl][j][i] = -nz[k1 - 1 - gl][j][i];
        }

      #else
        for (int j = JsL; j < JeL; j++)
        for (int i = IsL; i < IeL; i++)
        for (int gl = 0; gl < NG; gl++) {
            nx[k1 + gl][j][i] = nx[k1 - 1][j][i];
            ny[k1 + gl][j][i] = ny[k1 - 1][j][i];
            nz[k1 + gl][j][i] = nz[k1 - 1][j][i];
        }
      #endif
    }
#endif /* !ZPERIODIC */

    /*=====================================================================
     *  4)  MPI ghost exchange
     *=====================================================================*/
    Communication_update_ghost_nodes_flow_variable(
        nx, VOLUME_FRACTION, params->ghost_nodes, data_bag);
    Communication_update_ghost_nodes_flow_variable(
        ny, VOLUME_FRACTION, params->ghost_nodes, data_bag);
    Communication_update_ghost_nodes_flow_variable(
        nz, VOLUME_FRACTION, params->ghost_nodes, data_bag);
}


/*============================================================================
 *  VOF_smooth_contact_angle  —  UPDATED
 *
 *  Change: replaced the three separate Communication_update calls at the
 *  end with a single call to VOF_set_boundary_values_normal_vector, which
 *  fills physical-boundary ghosts (with proper symmetry) AND then calls
 *  the MPI communicator.
 *============================================================================*/
#ifdef VOF_IBM
void VOF_smooth_contact_angle(Cart3d_bag *data_bag)
{
    MAC_grid       *g   = data_bag->grid;
    VolumeFraction *v   = data_bag->vof;
    Parameters     *p   = data_bag->params;

    const int IsL = g->L_Is, IeL = g->L_Ie;
    const int JsL = g->L_Js, JeL = g->L_Je;
    const int KsL = g->L_Ks, KeL = g->L_Ke;

    const int Is = g->G_Is, Ie = g->G_Ie;
    const int Js = g->G_Js, Je = g->G_Je;
    const int Ks = g->G_Ks, Ke = g->G_Ke;

    const double dx = g->dx_c[0], dy = g->dy_c[0], dz = g->dz_c[0];
    const double theta = p->contact_angle_deg * PI / 180.0;
    const double sin_t = sin(theta), cos_t = cos(theta);
    const double epsN  = 1e-6;

    double ***Fs   = v->F_smooth;
    double ***nxS  = v->normal_x_smooth;
    double ***nyS  = v->normal_y_smooth;
    double ***nzS  = v->normal_z_smooth;

    CELESTE_NormalCtx n_ctx;
    CELESTE_precompute_normal_matrix(dx, dy, dz, &n_ctx);

    double ***nx_IBM = v->nx_IBM;
    double ***ny_IBM = v->ny_IBM;
    double ***nz_IBM = v->nz_IBM;
    double ***vfc    = v->vfc_smooth;

    for (int k = Ks; k < Ke; ++k)
    for (int j = Js; j < Je; ++j)
    for (int i = Is; i < Ie; ++i)
    {
        const double vcc = vfc[k][j][i];
        if (vcc < 0.05 || vcc > 0.95) continue;

        /* 1. Compute base fluid normal using WLS */
        double gradF[3];
        CELESTE_normal_from_precomp(&n_ctx, i, j, k, Fs, gradF);

        double nf[3] = { gradF[0], gradF[1], gradF[2] };
        double nf_len = sqrt(nf[0]*nf[0] + nf[1]*nf[1] + nf[2]*nf[2]);
        if (nf_len < epsN) continue;
        nf[0] /= nf_len;  nf[1] /= nf_len;  nf[2] /= nf_len;

        /* 2. Wall normal */
        double ns[3] = { nx_IBM[k][j][i], ny_IBM[k][j][i], nz_IBM[k][j][i] };
        double ns_len = sqrt(ns[0]*ns[0] + ns[1]*ns[1] + ns[2]*ns[2]);
        if (ns_len < epsN) continue;
        ns[0] /= ns_len;  ns[1] /= ns_len;  ns[2] /= ns_len;

        /* 3. Robust tangent construction */
        double nf_dot_ns = nf[0]*ns[0] + nf[1]*ns[1] + nf[2]*ns[2];
        double t[3] = { nf[0] - nf_dot_ns*ns[0],
                        nf[1] - nf_dot_ns*ns[1],
                        nf[2] - nf_dot_ns*ns[2] };
        double t_len = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);

        double nstar[3];
        if (t_len < epsN) {
            nstar[0] = nf[0];  nstar[1] = nf[1];  nstar[2] = nf[2];
        } else {
            t[0] /= t_len;  t[1] /= t_len;  t[2] /= t_len;

            /* 4. Rotation to enforce contact angle */
            nstar[0] = ns[0]*cos_t + t[0]*sin_t;
            nstar[1] = ns[1]*cos_t + t[1]*sin_t;
            nstar[2] = ns[2]*cos_t + t[2]*sin_t;
        }

        nxS[k][j][i] = nstar[0];
        nyS[k][j][i] = nstar[1];
        nzS[k][j][i] = nstar[2];
    }

    /* ── Physical boundary fills + MPI exchange ── */
    VOF_set_boundary_values_normal_vector(v->normal_x_smooth,
                                         v->normal_y_smooth,
                                         v->normal_z_smooth,
                                         data_bag);
}
#endif


/*=============================================================================
 *  curvature_patel  —  UPDATED
 *
 *  Change: replaced the three VOF_set_boundary_values() calls on normals
 *  after STEP 1 with a single VOF_set_boundary_values_normal_vector() call.
 *
 *  Curvature (kappa) boundary fill remains VOF_set_boundary_values() —
 *  it is a scalar, symmetric at the mirror plane, and not fed into any
 *  further stencil, so zero-gradient is correct.
 *=============================================================================*/
void curvature_patel(Cart3d_bag *data_bag)
{
    MAC_grid        *grid   = data_bag->grid;
    Parameters      *params = data_bag->params;
    VolumeFraction  *vof    = data_bag->vof;

    double ***F_s     = vof->F_smooth;
    double ***nx_s    = vof->normal_x_smooth;
    double ***ny_s    = vof->normal_y_smooth;
    double ***nz_s    = vof->normal_z_smooth;
    double ***kappa   = vof->kappa;
    double ***grad_mag = vof->grad_mag;

    const double dx = grid->dx_c[0], dy = grid->dy_c[0], dz = grid->dz_c[0];

    const int Is = grid->G_Is, Ie = min(grid->G_Ie, grid->NX - 1);
    const int Js = grid->G_Js, Je = min(grid->G_Je, grid->NY - 1);
    const int Ks = grid->G_Ks, Ke = min(grid->G_Ke, grid->NZ - 1);

    const int IsL = grid->L_Is, IeL = grid->L_Ie;
    const int JsL = grid->L_Js, JeL = grid->L_Je;
    const int KsL = grid->L_Ks, KeL = grid->L_Ke;

    const double eps_grad = 1e-6;

    CELESTE_NormalCtx    n_ctx;
    CELESTE_CurvatureCtx k_ctx;
    CELESTE_precompute_normal_matrix(dx, dy, dz, &n_ctx);
    CELESTE_precompute_curvature_matrix(dx, dy, dz, &k_ctx);

    /*----------------------------------------------------------------------
     *  STEP 1:  Compute unit normals  n̂ = ∇F̃ / |∇F̃|
     *--------------------------------------------------------------------*/
    const int i0 = IsL + 1, i1 = IeL - 1;
    const int j0 = JsL + 1, j1 = JeL - 1;
    const int k0 = KsL + 1, k1 = KeL - 1;

    for (int k = k0; k < k1; ++k)
    for (int j = j0; j < j1; ++j)
    for (int i = i0; i < i1; ++i)
    {
        double n[3];
        double grad_magnit = CELESTE_normal_from_precomp(&n_ctx, i, j, k, F_s, n);
        nx_s[k][j][i] = n[0];
        ny_s[k][j][i] = n[1];
        nz_s[k][j][i] = n[2];
        grad_mag[k][j][i] = grad_magnit;
    }

    /* ── CHANGED: vector-aware boundary fill for normals ── */
    VOF_set_boundary_values_normal_vector(vof->normal_x_smooth,
                                         vof->normal_y_smooth,
                                         vof->normal_z_smooth,
                                         data_bag);

    /*----------------------------------------------------------------------
     *  STEP 2:  Enforce contact angle in cut-cells
     *--------------------------------------------------------------------*/
#ifdef VOF_IBM
#ifdef SURFACE_TENSION
    VOF_smooth_contact_angle(data_bag);
#endif
#endif

    /*----------------------------------------------------------------------
     *  STEP 3:  Compute curvature  κ = −∇·n̂
     *--------------------------------------------------------------------*/
    const double h_min    = fmin(dx, fmin(dy, dz));
    const double kappa_max = 1.0 / h_min;

    for (int k = Ks; k < Ke; ++k)
    for (int j = Js; j < Je; ++j)
    for (int i = Is; i < Ie; ++i)
    {
        if (grad_mag[k][j][i] < eps_grad) {
            kappa[k][j][i] = nodata;
            continue;
        }

        double kval;
        CELESTE_curvature_from_precomp(&k_ctx, i, j, k,
                                       nx_s, ny_s, nz_s, &kval);

        if (kval >  kappa_max) kval =  kappa_max;
        if (kval < -kappa_max) kval = -kappa_max;

        kappa[k][j][i] = kval;
    }

    /* Kappa is a scalar — zero-gradient is correct (symmetric, no stencil reads) */
    VOF_set_boundary_values(vof->kappa, data_bag);
}




/*=============================================================================
   VOF_compute_f_sigma
   Assembles the face-centred surface tension force
       f_sigma = (2/We) * kappa_f * gradF
   from the CURRENT time-level fields (F^k, kappa^k) and stores the result
   in vof->f_sigma_new_{x,y,z}.

   This function does NOT touch the momentum RHS.  The RHS is updated
   separately by VOF_apply_f_sigma_old(), which uses the PREVIOUS substep's
   stored force to maintain the balanced-force property.

   Call order inside one RK substep:
       VoF_smoothing()          -> F_smooth^k
       curvature_patel()        -> kappa^k  (stored in vof->kappa)
       VOF_compute_f_sigma()    -> f_sigma_new  =  (2/We) kappa^k gradF^k
       VOF_apply_f_sigma_old()  -> RHS  +=  f_sigma_old  (= (2/We) kappa^{k-1} gradF^{k-1})
       ... velocity solve, pressure solve ...
       VOF_swap_f_sigma()       -> f_sigma_old  <-  f_sigma_new
=============================================================================*/
void VOF_compute_f_sigma(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    VolumeFraction *vof    = data_bag->vof;
    Parameters     *params = data_bag->params;

    double ***F     = vof->F;        /* raw (unsmoothed) VOF – same field used for gradF */
    double ***kappa = vof->kappa;    /* curvature computed from F_smooth^k by curvature_patel */

    double ***f_new_x = vof->f_sigma_new_x;
    double ***f_new_y = vof->f_sigma_new_y;
    double ***f_new_z = vof->f_sigma_new_z;

    const double scale = 2.0 / params->We;

    int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;
    int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
    int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;

    /* ---------- loop bounds: identical to Velocity_add_surfacetension_2_RHS_patel --- */
    int i_start_u = max(1, Is),      j_start_u = Js,           k_start_u = Ks;
    int i_end_u   = min(NX-1, Ie),   j_end_u   = min(NY-1, Je), k_end_u   = min(NZ-1, Ke);
#ifdef XPERIODIC
    if (i_end_u == NX-1) i_end_u = NX;
#endif

    int i_start_v = Is,              j_start_v = max(1, Js),   k_start_v = Ks;
    int i_end_v   = min(NX-1, Ie),   j_end_v   = min(NY-1, Je), k_end_v   = min(NZ-1, Ke);
#ifdef YPERIODIC
    if (j_end_v == NY-1) j_end_v = NY;
#endif

    int i_start_w = Is,              j_start_w = Js,            k_start_w = max(1, Ks);
    int i_end_w   = min(NX-1, Ie),   j_end_w   = min(NY-1, Je), k_end_w   = min(NZ-1, Ke);
#ifdef ZPERIODIC
    if (k_end_w == NZ-1) k_end_w = NZ;
#endif

    /* -----------------------------------------------------------------------
       x-faces  (u-momentum staggering: face (i-1/2, j, k) sits between
       cell (i-1,j,k) and cell (i,j,k))
    ----------------------------------------------------------------------- */
    for (int k = k_start_u; k < k_end_u; ++k)
    for (int j = j_start_u; j < j_end_u; ++j)
    for (int i = i_start_u; i < i_end_u; ++i)
    {
        /* Gradient of raw F across the face – same stencil as pressure gradient
           to guarantee the balanced-force cancellation.                       */
        double gradF = (F[k][j][i] - F[k][j][i-1]) * grid->idx_u[i-1];

        if (fabs(gradF) < 1e-6) { f_new_x[k][j][i] = 0.0; continue; }

        double k1 = kappa[k][j][i  ];   /* cell (i  , j, k) */
        double k2 = kappa[k][j][i-1];   /* cell (i-1, j, k) */
        bool k1_ok = (k1 < 1e19);
        bool k2_ok = (k2 < 1e19);

        if (!k1_ok && !k2_ok) { f_new_x[k][j][i] = 0.0; continue; }

        double kf = (k1_ok && k2_ok) ? 0.5*(k1 + k2) : (k1_ok ? k1 : k2);

#ifdef VOF_IBM
        /* Do not apply inside the solid.  Store zero so f_sigma_old is also
           zero at solid faces when it is later used as f_sigma_old.         */
        if (data_bag->lag->ng_vfu[k][j][i] > (1.0 - 1e-6) ) { f_new_x[k][j][i] = 0.0; continue; }
#endif
        f_new_x[k][j][i] = scale * kf * gradF;
    }

    /* -----------------------------------------------------------------------
       y-faces  (v-momentum)
    ----------------------------------------------------------------------- */
    for (int k = k_start_v; k < k_end_v; ++k)
    for (int j = j_start_v; j < j_end_v; ++j)
    for (int i = i_start_v; i < i_end_v; ++i)
    {
        double gradF = (F[k][j][i] - F[k][j-1][i]) * grid->idy_v[j-1];

        if (fabs(gradF) < 1e-6) { f_new_y[k][j][i] = 0.0; continue; }

        double k1 = kappa[k][j  ][i];
        double k2 = kappa[k][j-1][i];
        bool k1_ok = (k1 < 1e19);
        bool k2_ok = (k2 < 1e19);

        if (!k1_ok && !k2_ok) { f_new_y[k][j][i] = 0.0; continue; }

        double kf = (k1_ok && k2_ok) ? 0.5*(k1 + k2) : (k1_ok ? k1 : k2);

#ifdef VOF_IBM
        if (data_bag->lag->ng_vfv[k][j][i] > (1.0 - 1e-6) ) { f_new_y[k][j][i] = 0.0; continue; }
#endif
        f_new_y[k][j][i] = scale * kf * gradF;
    }

    /* -----------------------------------------------------------------------
       z-faces  (w-momentum)
    ----------------------------------------------------------------------- */
    for (int k = k_start_w; k < k_end_w; ++k)
    for (int j = j_start_w; j < j_end_w; ++j)
    for (int i = i_start_w; i < i_end_w; ++i)
    {
        double gradF = (F[k][j][i] - F[k-1][j][i]) * grid->idz_w[k-1];

        if (fabs(gradF) < 1e-6) { f_new_z[k][j][i] = 0.0; continue; }

        double k1 = kappa[k  ][j][i];
        double k2 = kappa[k-1][j][i];
        bool k1_ok = (k1 < 1e19);
        bool k2_ok = (k2 < 1e19);

        if (!k1_ok && !k2_ok) { f_new_z[k][j][i] = 0.0; continue; }

        double kf = (k1_ok && k2_ok) ? 0.5*(k1 + k2) : (k1_ok ? k1 : k2);

#ifdef VOF_IBM
        if (data_bag->lag->ng_vfw[k][j][i] > (1.0 - 1e-6) ) { f_new_z[k][j][i] = 0.0; continue; }
#endif
        f_new_z[k][j][i] = scale * kf * gradF;
    }
}







#ifdef VOF_IBM

//------------------------------------------------------------------------------
// VOF_normals_IBM
//   UPDATED: Uses WLS on the raw vfc field.
//   Physics Note: Because 'vfc' is sharp (0 or 1), the normals here will likely
//   be somewhat grid-aligned (staircased). This is useful for collision detection
//   but suboptimal for contact angles.
//------------------------------------------------------------------------------
void VOF_normals_IBM(Cart3d_bag *data_bag)
{
  MAC_grid       *grid = data_bag->grid;
  VolumeFraction *vof  = data_bag->vof;
  Lagrangian     *lag  = data_bag->lag;

  double ***vfc  = lag->ng_vfc;
  double ***vfs  = vof->vfc; // Local copy used for stencil access

  const double eps = 1e-30;
  const double dx = grid->dx_c[0], dy = grid->dy_c[0], dz = grid->dz_c[0];

  // Global Interior bounds (where we compute results)
  const int Is = grid->G_Is, Ie = grid->G_Ie;
  const int Js = grid->G_Js, Je = grid->G_Je;
  const int Ks = grid->G_Ks, Ke = grid->G_Ke;
  
  // Local bounds (for stencil validity checks)
  const int IsL = grid->L_Is, IeL = grid->L_Ie;
  const int JsL = grid->L_Js, JeL = grid->L_Je;
  const int KsL = grid->L_Ks, KeL = grid->L_Ke;

  // 1. Copy volume fraction field to local array
  //    (We do this to ensure we have a writable copy with ghost layers)
  for (int k = Ks; k < Ke; ++k) {
    for (int j = Js; j < Je; ++j) {
      for (int i = Is; i < Ie; ++i) {
        vfs[k][j][i] = vfc[k][j][i];
      }
    }
  }
  
  // 2. Fill ghosts for vfs so WLS has neighbor data
  VOF_set_boundary_values(vfs, data_bag);

}

/*-----------------------------------------------------------------------------
 * Vfc_smoothing  (Approach A – computes vfc_smooth AND IBM normals)
 *
 *  CHANGE LOG (band-restricted normals):
 *  -------------------------------------
 *  - Added eta_band_fluid / eta_band_solid tunables.
 *  - Step 3: normals are only written when eta falls within the band.
 *            vfc_smooth (eta scratch) is still written for ALL cells
 *            in the bounding box so the tanh indicator is correct everywhere.
 *  - Step 4: normalisation uses a magnitude threshold instead of the old
 *            sqrt(...) + eps trick, which silently normalised zero vectors.
 *            Out-of-band cells keep (0,0,0) from initialisation.
 *---------------------------------------------------------------------------*/
void Vfc_smoothing(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    VolumeFraction *vof    = data_bag->vof;
    Parameters     *params = data_bag->params;
    Lagrangian     *lag    = data_bag->lag;

    double ***vfc_smooth = vof->vfc_smooth;        /* output: smooth indicator */
    double ***nxB        = vof->nx_IBM;      /* output: normal x        */
    double ***nyB        = vof->ny_IBM;      /* output: normal y        */
    double ***nzB        = vof->nz_IBM;      /* output: normal z        */

    /* Grid geometry */
    const double *xc = grid->xc;   /* cell centers */
    const double *yc = grid->yc;
    const double *zc = grid->zc;

    const double *xe = grid->xu;   /* cell edges (for bbox index) */
    const double *ye = grid->yv;
    const double *ze = grid->zw;

    /* Global interior extents for this processor */
    const int Is = grid->G_Is;
    const int Js = grid->G_Js;
    const int Ks = grid->G_Ks;
    const int Ie = grid->G_Ie;   /* one-past-end */
    const int Je = grid->G_Je;
    const int Ke = grid->G_Ke;

    /* Uniform grid spacing */
    const double h = grid->dx_u[1];

    /* =====================================================================
     * TUNABLE: delta_s  (transition half-width)
     * =====================================================================*/
    const double delta_s = 1.0 * h;

    const double cutoff = 3.0 * delta_s;

    /* =====================================================================
     * TUNABLE: normal band extents (relative to interface at eta = 0)
     *
     *   eta_band_fluid : how far into the FLUID side (eta < 0) normals
     *                    are stored.  2*delta_s ≈ 2h covers vfc_smooth
     *                    down to ~0.02, well below the rotation threshold
     *                    of 0.05.
     *
     *   eta_band_solid : how far into the SOLID side (eta > 0) normals
     *                    are stored.  5h covers the 4-cell extension
     *                    depth (8 iters × Δτ = 0.5Δx) with margin.
     * =====================================================================*/
    const double eta_band_fluid = 2.0 * delta_s;   /* ~2h  */
    const double eta_band_solid = 5.0 * h;          /* ~5h  */

    /* =================================================================
     * Step 1: Gather ALL particle (X, R) across all MPI ranks.
     * =================================================================*/

    int n_local = 0;
    {
        Particle_list *p_mobile_list = lag->p_mobile_list;
        Particle_list *p_fixed_list  = lag->p_fixed_list;

        if (p_mobile_list != NULL) {
            Particle *p = p_mobile_list->start;
            while (p != NULL) { n_local++; p = p->next; }
        }
        if (p_fixed_list != NULL) {
            Particle *p = p_fixed_list->start;
            while (p != NULL) { n_local++; p = p->next; }
        }
    }

    double *local_buf = (double *)malloc(4 * n_local * sizeof(double));
    {
        int idx = 0;
        Particle_list *p_mobile_list = lag->p_mobile_list;
        Particle_list *p_fixed_list  = lag->p_fixed_list;

        if (p_mobile_list != NULL) {
            Particle *p = p_mobile_list->start;
            while (p != NULL) {
                local_buf[idx++] = p->X[0];
                local_buf[idx++] = p->X[1];
                local_buf[idx++] = p->X[2];
                local_buf[idx++] = p->R;
                p = p->next;
            }
        }
        if (p_fixed_list != NULL) {
            Particle *p = p_fixed_list->start;
            while (p != NULL) {
                local_buf[idx++] = p->X[0];
                local_buf[idx++] = p->X[1];
                local_buf[idx++] = p->X[2];
                local_buf[idx++] = p->R;
                p = p->next;
            }
        }
    }

    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *counts = (int *)malloc(nprocs * sizeof(int));
    int *displs = (int *)malloc(nprocs * sizeof(int));

    int send_count = 4 * n_local;
    MPI_Allgather(&send_count, 1, MPI_INT,
                  counts,      1, MPI_INT, MPI_COMM_WORLD);

    int n_total_doubles = 0;
    for (int r = 0; r < nprocs; ++r) {
        displs[r] = n_total_doubles;
        n_total_doubles += counts[r];
    }

    double *all_buf = (double *)malloc(n_total_doubles * sizeof(double));

    MPI_Allgatherv(local_buf, send_count,        MPI_DOUBLE,
                   all_buf,   counts, displs,     MPI_DOUBLE,
                   MPI_COMM_WORLD);

    int n_particles_global = n_total_doubles / 4;

    free(local_buf);
    free(counts);
    free(displs);

    /* =================================================================
     * Step 2: Initialise scratch fields.
     *
     * vfc_smooth  ← eta_max sentinel (-1e30)
     * nxB/nyB/nzB ← 0  (cells outside the band will keep this value)
     * =================================================================*/
    for (int k = Ks; k < Ke; ++k)
    for (int j = Js; j < Je; ++j)
    for (int i = Is; i < Ie; ++i) {
        vfc_smooth[k][j][i] = -1.0e30;
        nxB[k][j][i] = 0.0;
        nyB[k][j][i] = 0.0;
        nzB[k][j][i] = 0.0;
    }

    /* =================================================================
     * Step 3: Accumulate max(eta) and store direction to closest
     *         particle — BUT only within the normal band.
     *
     *   vfc_smooth ← eta   (written for ALL cells, converted later)
     *   nxB/nyB/nzB         (written ONLY when eta is in-band)
     *
     * Band:  -eta_band_fluid  <  eta  <  +eta_band_solid
     * =================================================================*/
    for (int ip = 0; ip < n_particles_global; ++ip) {

        const double Xp = all_buf[4*ip + 0];
        const double Yp = all_buf[4*ip + 1];
        const double Zp = all_buf[4*ip + 2];
        const double R  = all_buf[4*ip + 3];

        const double extent = R + cutoff;

        int i_start = (int)floor((Xp - extent - xe[0]) / h);
        int j_start = (int)floor((Yp - extent - ye[0]) / h);
        int k_start = (int)floor((Zp - extent - ze[0]) / h);

        int i_end = (int)ceil((Xp + extent - xe[0]) / h);
        int j_end = (int)ceil((Yp + extent - ye[0]) / h);
        int k_end = (int)ceil((Zp + extent - ze[0]) / h);

        i_start = max(i_start, Is);
        j_start = max(j_start, Js);
        k_start = max(k_start, Ks);

        i_end = min(i_end, Ie - 1);
        j_end = min(j_end, Je - 1);
        k_end = min(k_end, Ke - 1);

        for (int k = k_start; k <= k_end; ++k) {
            const double dz  = zc[k] - Zp;
            const double dz2 = dz * dz;
            for (int j = j_start; j <= j_end; ++j) {
                const double dy  = yc[j] - Yp;
                const double dy2 = dy * dy;
                for (int i = i_start; i <= i_end; ++i) {
                    const double dx_  = xc[i] - Xp;
                    const double dist = sqrt(dx_ * dx_ + dy2 + dz2);
                    const double eta  = R - dist;

                    if (eta > vfc_smooth[k][j][i]) {

                        /* Always update eta scratch — the tanh
                         * indicator must be correct everywhere. */
                        vfc_smooth[k][j][i] = eta;

                        /* Only store normal direction if this cell
                         * falls within the band needed by rotation
                         * and extension.  Cells outside keep (0,0,0)
                         * from Step 2 initialisation.                */
                        if (eta > -eta_band_fluid &&
                            eta <  eta_band_solid   ) {
                            nxB[k][j][i] = -dx_;   /* Xp - xc[i] */
                            nyB[k][j][i] = -dy;    /* Yp - yc[j] */
                            nzB[k][j][i] = -dz;    /* Zp - zc[k] */
                        }
                    }
                }
            }
        }
    }

    free(all_buf);

    /* =================================================================
     * Step 4: Convert eta → vfc_smooth (tanh) and normalise normals.
     *
     * Normalisation uses a magnitude threshold (eps) instead of the
     * old  sqrt(...) + 1e-30  trick, which silently normalised zero
     * vectors into arbitrary directions.  Out-of-band cells whose
     * normal was never written remain exactly (0,0,0).
     * =================================================================*/
    const double eps = 1.0e-30;

    for (int k = Ks; k < Ke; ++k)
    for (int j = Js; j < Je; ++j)
    for (int i = Is; i < Ie; ++i) {

        const double eta = vfc_smooth[k][j][i];

        if (eta < -1.0e29) {
            /* No particle in range → pure fluid, zero normal */
            vfc_smooth[k][j][i] = 0.0;
            /* nxB/nyB/nzB already 0.0 from Step 2 */
        }
        else {
            /* --- tanh profile --- */
            vfc_smooth[k][j][i] = 0.5 * (1.0 + tanh(eta / delta_s));

            /* --- Normalise direction to unit normal --- */
            const double nx = nxB[k][j][i];
            const double ny = nyB[k][j][i];
            const double nz = nzB[k][j][i];
            const double mag = sqrt(nx*nx + ny*ny + nz*nz);

            if (mag > eps) {
                nxB[k][j][i] = nx / mag;
                nyB[k][j][i] = ny / mag;
                nzB[k][j][i] = nz / mag;
            }
            else {
                /* Outside the normal band (or exact particle center):
                 * zero normal — safe sentinel for downstream checks
                 * (rotation and extension both guard on |n| > epsN). */
                nxB[k][j][i] = 0.0;
                nyB[k][j][i] = 0.0;
                nzB[k][j][i] = 0.0;
            }
        }
    }

    /* =================================================================
     * Step 5: Fill ghost cells (parallel exchange + physical BCs).
     * =================================================================*/
    VOF_set_boundary_values(vof->vfc_smooth, data_bag);
    VOF_set_boundary_values_normal_vector(vof->nx_IBM,
                                         vof->ny_IBM,
                                         vof->nz_IBM,
                                         data_bag);
}


/******************************************************************************
 * Continuum Capillary Force (CCF) model for contact-line forces
 * Based on Konstantidinis et al. PF'24
 ******************************************************************************/



#define SEPS 1e-30


//------------------------------------------------------------------------------
// VOF_accumulate_solid_capillary_force
//------------------------------------------------------------------------------
void VOF_accumulate_solid_capillary_force(Particle *p, Cart3d_bag *data_bag) 
{
    MAC_grid       *grid   = data_bag->grid;
    Parameters     *params = data_bag->params;
    VolumeFraction *vof    = data_bag->vof;

    // Solid IBM normals pointing into solid
    double ***nx_s = vof->nx_IBM; 
    double ***ny_s = vof->ny_IBM;
    double ***nz_s = vof->nz_IBM;

    // CCF force density fields (cell-centered)
    double ***f_ccf_x = vof->f_ccf_x;
    double ***f_ccf_y = vof->f_ccf_y;
    double ***f_ccf_z = vof->f_ccf_z;
    
    // Fields
    double ***vfc  = vof->vfc;  // Smoothed solid volume fraction
    double ***F    = vof->F_smooth;      // Smoothed fluid color function

    // Geometry
    double *xc = grid->xc;
    double *yc = grid->yc;
    double *zc = grid->zc;
    
    double dx = grid->dx_c[0];
    double dy = grid->dy_c[0];
    double dz = grid->dz_c[0];
    double dV = dx * dy * dz;

    // Physics
    double sigma = 1.0 / params->We;
    
    // Contact angle parameters
    const double theta = params->contact_angle_deg * PI / 180.0;
    const double sin_t = sin(theta), cos_t = cos(theta);

    double X_p[3] = {p->X[0], p->X[1], p->X[2]};
    
    // Use GLOBAL interior extents
    const int Is = grid->G_Is;
    const int Js = grid->G_Js;
    const int Ks = grid->G_Ks;
    const int Ie = grid->G_Ie;
    const int Je = grid->G_Je;
    const int Ke = grid->G_Ke;
    const double epsN = 1e-14;

    Memory_reset_flow_variable(grid, params, f_ccf_x);
    Memory_reset_flow_variable(grid, params, f_ccf_y);
    Memory_reset_flow_variable(grid, params, f_ccf_z);

    for (int k = Ks; k < Ke; k++) {
        for (int j = Js; j < Je; j++) {
            for (int i = Is; i < Ie; i++) {


                // 1. Calculate Gradient of Fluid Color Function (phi)
                double gphi[3];
                gphi[0] = (F[k][j][i+1] - F[k][j][i-1]) / (2.0*dx);
                gphi[1] = (F[k][j+1][i] - F[k][j-1][i]) / (2.0*dy);
                gphi[2] = (F[k+1][j][i] - F[k-1][j][i]) / (2.0*dz);


                double mag_gphi = sqrt(gphi[0]*gphi[0] + gphi[1]*gphi[1] + gphi[2]*gphi[2]);

                // 2. Calculate Gradient of Solid Volume Fraction (phi_s)
                double gphis[3];
                gphis[0] = (vfc[k][j][i+1] - vfc[k][j][i-1]) / (2.0*dx);
                gphis[1] = (vfc[k][j+1][i] - vfc[k][j-1][i]) / (2.0*dy);
                gphis[2] = (vfc[k+1][j][i] - vfc[k-1][j][i]) / (2.0*dz);
                double mag_gphis = sqrt(gphis[0]*gphis[0] + gphis[1]*gphis[1] + gphis[2]*gphis[2]);

                if (mag_gphi < epsN || mag_gphis < epsN) continue;

                // 3. Define Unit Normals
                double ns[3] = {nx_s[k][j][i], ny_s[k][j][i], nz_s[k][j][i]};
                double ns_len = sqrt(ns[0]*ns[0] + ns[1]*ns[1] + ns[2]*ns[2]);
                if (ns_len < epsN) continue;
                ns[0] /= ns_len; ns[1] /= ns_len; ns[2] /= ns_len;

                // Fluid normal
                double nf[3] = {gphi[0]/mag_gphi, gphi[1]/mag_gphi, gphi[2]/mag_gphi};

                // 4. Fluid normal projection onto wall tangent plane
                double t_wall[3];
                double nf_dot_ns = nf[0]*ns[0] + nf[1]*ns[1] + nf[2]*ns[2];
                t_wall[0] = nf[0] - (nf_dot_ns)*ns[0];
                t_wall[1] = nf[1] - (nf_dot_ns)*ns[1];
                t_wall[2] = nf[2] - (nf_dot_ns)*ns[2];
                double t_len = sqrt(t_wall[0]*t_wall[0] + t_wall[1]*t_wall[1] + t_wall[2]*t_wall[2]);

                double tc[3];
                // Handle singularity
                if ( t_len < 1e-10) {
                    // For physics safety
                    tc[0] = ns[0]; tc[1] = ns[1]; tc[2] = ns[2];
                } else {
                    
                    t_wall[0]/=t_len; t_wall[1]/=t_len; t_wall[2]/=t_len;

                    tc[0] = - (ns[0] - (nf_dot_ns)*nf[0]);
                    tc[1] = - (ns[1] - (nf_dot_ns)*nf[1]);
                    tc[2] = - (ns[2] - (nf_dot_ns)*nf[2]);
                    double tc_len = sqrt(tc[0]*tc[0] + tc[1]*tc[1] + tc[2]*tc[2]);
                    if (tc_len < epsN) continue;
                    tc[0] /= tc_len; tc[1] /= tc_len; tc[2] /= tc_len;

                }

                
                // 7. Calculate Magnitude Term
                
                double term_phi  = gphi[0]*t_wall[0] + gphi[1]*t_wall[1] + gphi[2]*t_wall[2];
                double term_phis = gphis[0]*ns[0] + gphis[1]*ns[1] + gphis[2]*ns[2];

                // 8. Force Calculation
                double factor = sigma * term_phi * term_phis * dV;
                double dF[3] = { factor * tc[0], factor * tc[1], factor * tc[2] };

                // For output
                f_ccf_x[k][j][i] = dF[0] / dV;
                f_ccf_y[k][j][i] = dF[1] / dV;
                f_ccf_z[k][j][i] = dF[2] / dV;

                // 9. Accumulate Force on Particle
                p->F_CCF[0] += dF[0];
                p->F_CCF[1] += dF[1];
                p->F_CCF[2] += dF[2];

                // 10. Accumulate Torque
                double r_vec[3];
                r_vec[0] = xc[i] - X_p[0]; 
                r_vec[1] = yc[j] - X_p[1];
                r_vec[2] = zc[k] - X_p[2];

                p->T_CCF[0] += (dF[1] * r_vec[2] - dF[2] * r_vec[1]);
                p->T_CCF[1] += (dF[2] * r_vec[0] - dF[0] * r_vec[2]);
                p->T_CCF[2] += (dF[0] * r_vec[1] - dF[1] * r_vec[0]);
            }
        }
    }


}



// /******************************************************************************/
// /*
//  * VOF_integrate_CSF_over_solid
//  *
//  * Integrates the CSF surface tension force over all solid particles and
//  * adds the resulting force and torque to each particle's F[] and T[] arrays.
//  *
//  * For a particle p with solid volume fraction α on each staggered face:
//  *   F_σ^(p) = Σ α_f · (2/We) · κ_f · (∇F)_f · ΔV
//  *   T_σ^(p) = Σ α_f · r × [ (2/We) · κ_f · (∇F)_f · ΔV ]
//  *
//  * This must be called AFTER:
//  *   - VOF curvature computation (kappa field valid)
//  *   - Lagrangian volume fractions computed (ng_vfu, ng_vfv, ng_vfw)
//  *   - BEFORE particle force integration
//  */
// /******************************************************************************/
// void VOF_integrate_CSF_over_solid(Particle *p, Cart3d_bag *data_bag) {
    
//     MAC_grid       *grid   = data_bag->grid;
//     VolumeFraction *vof    = data_bag->vof;
//     Parameters     *params = data_bag->params;
//     Lagrangian     *lag    = data_bag->lag;
    
//     double ***F      = vof->F;
//     double ***kappa  = vof->kappa;
    
//     /* Solid volume fractions on staggered faces */
//     double ***vfu = lag->ng_vfu;
//     double ***vfv = lag->ng_vfv;
//     double ***vfw = lag->ng_vfw;
    
//     const double scale = 2.0 / params->We;
//     const double h  = grid->dx_u[1];  /* uniform spacing */
//     const double dV = h * h * h;
    
//     int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;
//     int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
//     int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;
    
//     /* Loop bounds (same as Velocity_add_surfacetension_2_RHS_patel) */
//     int i_start_u = max(1, Is),     j_start_u = Js,           k_start_u = Ks;
//     int i_end_u   = min(NX-1, Ie),  j_end_u   = min(NY-1, Je), k_end_u = min(NZ-1, Ke);
//     #ifdef XPERIODIC
//         if (i_end_u == NX-1) i_end_u = NX;
//     #endif
    
//     int i_start_v = Is,             j_start_v = max(1, Js),   k_start_v = Ks;
//     int i_end_v   = min(NX-1, Ie),  j_end_v   = min(NY-1, Je), k_end_v = min(NZ-1, Ke);
//     #ifdef YPERIODIC
//         if (j_end_v == NY-1) j_end_v = NY;
//     #endif
    
//     int i_start_w = Is,             j_start_w = Js,           k_start_w = max(1, Ks);
//     int i_end_w   = min(NX-1, Ie),  j_end_w   = min(NY-1, Je), k_end_w = min(NZ-1, Ke);
//     #ifdef ZPERIODIC
//         if (k_end_w == NZ-1) k_end_w = NZ;
//     #endif
    
//     /* Face-centered coordinates for torque calculation */
//     double *xu = grid->xu;  /* x-coordinate of u-faces */
//     double *yc = grid->yc;  /* y-coordinate at cell centers (for u-faces) */
//     double *zc = grid->zc;  /* z-coordinate at cell centers (for u-faces) */
    
//     double *xc = grid->xc;  /* for v and w faces */
//     double *yv = grid->yv;
//     double *zw = grid->zw;

//         const double R = p->R;
//         const double *X = p->X;  /* Particle center */
//         double *Fp = p->F_CSF_solid;       /* Particle force accumulator */
//         double *Tp = p->T_CSF_solid;       /* Particle torque accumulator */
        
//         /* Determine bounding box of particle */
//         int pi_start = (int)floor((X[0] - R - xu[0]) / h);
//         int pj_start = (int)floor((X[1] - R - yc[0]) / h);
//         int pk_start = (int)floor((X[2] - R - zc[0]) / h);
        
//         int pi_end = (int)ceil((X[0] + R - xu[0]) / h) + 1;
//         int pj_end = (int)ceil((X[1] + R - yc[0]) / h) + 1;
//         int pk_end = (int)ceil((X[2] + R - zc[0]) / h) + 1;
        
//         /* ================================================================== */
//         /* U-component (x-direction surface tension force)                    */
//         /* ================================================================== */
//         for (int k = max(k_start_u, pk_start); k < min(k_end_u, pk_end); ++k) {
//             for (int j = max(j_start_u, pj_start); j < min(j_end_u, pj_end); ++j) {
//                 for (int i = max(i_start_u, pi_start); i < min(i_end_u, pi_end); ++i) {
                    
//                     double alpha = vfu[k][j][i];
//                     if (alpha < 1e-12) continue;  /* No solid here */
                    
//                     double gradF = (F[k][j][i] - F[k][j][i-1]) * grid->idx_u[i-1];
//                     if (fabs(gradF) < 1e-12) continue;
                    
//                     /* Face-centered curvature */
//                     double k1 = kappa[k][j][i];
//                     double k2 = kappa[k][j][i-1];
//                     bool k1_ok = (k1 < 1e19);
//                     bool k2_ok = (k2 < 1e19);
                    
//                     if (!k1_ok && !k2_ok) continue;
                    
//                     double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)
//                                                  : (k1_ok ? k1 : k2);
                    
//                     /* CSF force at this face */
//                     double f_sigma_x = scale * kf * gradF;
                    
//                     /* Weighted by solid fraction and cell volume */
//                     double dF_x = alpha * f_sigma_x * dV;
                    
//                     /* Position vector from particle center to face center */
//                     double rx = xu[i] - X[0];
//                     double ry = yc[j] - X[1];
//                     double rz = zc[k] - X[2];
                    
//                     /* Add force (same sign convention as IBM: force ON particle) */
//                     Fp[0] += dF_x;
                    
//                     /* Add torque: T = r × F, here F = (dF_x, 0, 0) */
//                     Tp[1] += rz * dF_x;   /* T_y = r_z * F_x - r_x * F_z */
//                     Tp[2] -= ry * dF_x;   /* T_z = r_x * F_y - r_y * F_x */
//                 }
//             }
//         }
        
//         /* ================================================================== */
//         /* V-component (y-direction surface tension force)                    */
//         /* ================================================================== */
//         for (int k = max(k_start_v, pk_start); k < min(k_end_v, pk_end); ++k) {
//             for (int j = max(j_start_v, pj_start); j < min(j_end_v, pj_end); ++j) {
//                 for (int i = max(i_start_v, pi_start); i < min(i_end_v, pi_end); ++i) {
                    
//                     double alpha = vfv[k][j][i];
//                     if (alpha < 1e-12) continue;
                    
//                     double gradF = (F[k][j][i] - F[k][j-1][i]) * grid->idy_v[j-1];
//                     if (fabs(gradF) < 1e-12) continue;
                    
//                     double k1 = kappa[k][j][i];
//                     double k2 = kappa[k][j-1][i];
//                     bool k1_ok = (k1 < 1e19);
//                     bool k2_ok = (k2 < 1e19);
                    
//                     if (!k1_ok && !k2_ok) continue;
                    
//                     double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)
//                                                  : (k1_ok ? k1 : k2);
                    
//                     double f_sigma_y = scale * kf * gradF;
//                     double dF_y = alpha * f_sigma_y * dV;
                    
//                     double rx = xc[i] - X[0];
//                     double ry = yv[j] - X[1];
//                     double rz = zc[k] - X[2];
                    
//                     Fp[1] += dF_y;
                    
//                     /* T = r × F, here F = (0, dF_y, 0) */
//                     Tp[0] -= rz * dF_y;   /* T_x = r_y * F_z - r_z * F_y */
//                     Tp[2] += rx * dF_y;   /* T_z = r_x * F_y - r_y * F_x */
//                 }
//             }
//         }
        
//         /* ================================================================== */
//         /* W-component (z-direction surface tension force)                    */
//         /* ================================================================== */
//         for (int k = max(k_start_w, pk_start); k < min(k_end_w, pk_end); ++k) {
//             for (int j = max(j_start_w, pj_start); j < min(j_end_w, pj_end); ++j) {
//                 for (int i = max(i_start_w, pi_start); i < min(i_end_w, pi_end); ++i) {
                    
//                     double alpha = vfw[k][j][i];
//                     if (alpha < 1e-12) continue;
                    
//                     double gradF = (F[k][j][i] - F[k-1][j][i]) * grid->idz_w[k-1];
//                     if (fabs(gradF) < 1e-12) continue;
                    
//                     double k1 = kappa[k][j][i];
//                     double k2 = kappa[k-1][j][i];
//                     bool k1_ok = (k1 < 1e19);
//                     bool k2_ok = (k2 < 1e19);
                    
//                     if (!k1_ok && !k2_ok) continue;
                    
//                     double kf = (k1_ok && k2_ok) ? 0.5 * (k1 + k2)
//                                                  : (k1_ok ? k1 : k2);
                    
//                     double f_sigma_z = scale * kf * gradF;
//                     double dF_z = alpha * f_sigma_z * dV;
                    
//                     double rx = xc[i] - X[0];
//                     double ry = yc[j] - X[1];
//                     double rz = zw[k] - X[2];
                    
//                     Fp[2] += dF_z;
                    
//                     /* T = r × F, here F = (0, 0, dF_z) */
//                     Tp[0] += ry * dF_z;   /* T_x = r_y * F_z - r_z * F_y */
//                     Tp[1] -= rx * dF_z;   /* T_y = r_z * F_x - r_x * F_z */
//                 }
//             }
//         }
        
     
//   }


#endif // VOF_IBM
#endif

#ifdef SURFACE_TENSION
/*=============================================================================
   VOF_apply_f_sigma_old
   Adds the PREVIOUS substep's stored surface tension force to the momentum RHS:
       RHS += f_sigma_old   (= (2/We) kappa^{k-1} gradF^{k-1})

   Because f_sigma_old was assembled when pressure was p^{k-1}, it exactly
   balances the lagged pressure gradient -2*grad(p^{k-1}) already present in
   the RHS, preventing spurious parasitic currents at the interface.

   This function only reads f_sigma_old_{x,y,z} and writes ng_rhs.
   It never reads vof->kappa or vof->F directly.
=============================================================================*/
void VOF_apply_f_sigma_old(Cart3d_bag *data_bag)
{
    MAC_grid       *grid   = data_bag->grid;
    VolumeFraction *vof    = data_bag->vof;
    Parameters     *params = data_bag->params;
    Velocity       *u      = data_bag->u;
    Velocity       *v      = data_bag->v;
    Velocity       *w      = data_bag->w;

    double ***f_old_x = vof->f_sigma_old_x;
    double ***f_old_y = vof->f_sigma_old_y;
    double ***f_old_z = vof->f_sigma_old_z;

    double ***rhsu = u->ng_rhs;
    double ***rhsv = v->ng_rhs;
    double ***rhsw = w->ng_rhs;

    int NX = grid->NX, NY = grid->NY, NZ = grid->NZ;
    int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
    int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;

    /* Loop bounds must be identical to VOF_compute_f_sigma so that only
       faces where a force was stored are touched.                          */
    int i_start_u = max(1, Is),      j_start_u = Js,            k_start_u = Ks;
    int i_end_u   = min(NX-1, Ie),   j_end_u   = min(NY-1, Je), k_end_u   = min(NZ-1, Ke);
#ifdef XPERIODIC
    if (i_end_u == NX-1) i_end_u = NX;
#endif

    int i_start_v = Is,              j_start_v = max(1, Js),    k_start_v = Ks;
    int i_end_v   = min(NX-1, Ie),   j_end_v   = min(NY-1, Je), k_end_v   = min(NZ-1, Ke);
#ifdef YPERIODIC
    if (j_end_v == NY-1) j_end_v = NY;
#endif

    int i_start_w = Is,              j_start_w = Js,             k_start_w = max(1, Ks);
    int i_end_w   = min(NX-1, Ie),   j_end_w   = min(NY-1, Je), k_end_w   = min(NZ-1, Ke);
#ifdef ZPERIODIC
    if (k_end_w == NZ-1) k_end_w = NZ;
#endif

    /* x-faces */
    for (int k = k_start_u; k < k_end_u; ++k)
    for (int j = j_start_u; j < j_end_u; ++j)
    for (int i = i_start_u; i < i_end_u; ++i)
        rhsu[k][j][i] += f_old_x[k][j][i];

    /* y-faces */
    for (int k = k_start_v; k < k_end_v; ++k)
    for (int j = j_start_v; j < j_end_v; ++j)
    for (int i = i_start_v; i < i_end_v; ++i)
        rhsv[k][j][i] += f_old_y[k][j][i];

    /* z-faces */
    for (int k = k_start_w; k < k_end_w; ++k)
    for (int j = j_start_w; j < j_end_w; ++j)
    for (int i = i_start_w; i < i_end_w; ++i)
        rhsw[k][j][i] += f_old_z[k][j][i];
}


/*=============================================================================
   VOF_swap_f_sigma
   Swaps the new and old surface tension force pointer pairs so that
   f_sigma_old always holds the force that was assembled at the just-completed
   substep k, ready to be used as the lagged force in substep k+1.

   This is a pure pointer swap: no data is copied.

   Call this ONCE per RK substep, AFTER the pressure projection and velocity
   correction are complete, and BEFORE advancing to the next substep.
=============================================================================*/
void VOF_swap_f_sigma(Cart3d_bag *data_bag)
{
    VolumeFraction *vof = data_bag->vof;
    double ***tmp;

    tmp = vof->f_sigma_old_x;
    vof->f_sigma_old_x = vof->f_sigma_new_x;
    vof->f_sigma_new_x = tmp;

    tmp = vof->f_sigma_old_y;
    vof->f_sigma_old_y = vof->f_sigma_new_y;
    vof->f_sigma_new_y = tmp;

    tmp = vof->f_sigma_old_z;
    vof->f_sigma_old_z = vof->f_sigma_new_z;
    vof->f_sigma_new_z = tmp;
}
#endif

#ifdef VOF_GRAVITY
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
#endif