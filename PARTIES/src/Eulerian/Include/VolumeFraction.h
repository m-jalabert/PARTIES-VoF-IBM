#ifndef VOLUMEFRACTION_H
 #define VOLUMEFRACTION_H

#include "DataTypes.h"

// "Constructor"-style function (similar to Conc_create)
VolumeFraction *VoF_create(MAC_grid *grid, Parameters *params);

// "Destructor"-style function (similar to Conc_destroy)
void VOF_destroy(VolumeFraction *vof, MAC_grid *grid, Parameters *params);

// Additional prototypes: boundary updates, advection, etc.
void VOF_set_boundary_values(double ***F, Cart3d_bag *data_bag);

/**
 *    Reconstruct interface (PLIC)
 *    We solve for alpha in each interface cell so that the plane
 *    is consistent with the local volume fraction F[i][j][k] and normal.
 */
void VOF_reconstruct_interface(Cart3d_bag *data_bag);

/**
 *    VoF_set_advection
 *
 * Computes convective terms for VOF advection using PLIC-reconstructed interface
 * Geometry-aware advection based on face fluxes calculated using interface normals
 */
void VOF_set_advection(Cart3d_bag *data_bag); 

/**
 *    VOF_set_RHS
 *
 * Updates the volume-fraction field `F` using a 3-stage Runge–Kutta (RK3) scheme.
 */
void VOF_update_F(Cart3d_bag *data_bag);

/**
 *    Compute curvature (kappa)
 *    Possibly using height-function or simpler approach. Fill vof->kappa.
 */
void VOF_compute_curvature (Cart3d_bag *data_bag);

/**
 *    Apply surface tension (CSF)
 *    Add sigma * kappa * grad(F) / rho to momentum or velocity body force
 */
void VOF_apply_surface_tension (Cart3d_bag *data_bag);

/**
 * Updates the dimensionless density and viscosity fields in the domain,
 * using the piecewise mixing laws:
 *
 *   (1)  rho_tilde(F) = F + (1 - F)* (rho2 / rho1)
 *
 *   (2)  rho_tilde / mu_tilde = F
 *           + (1 - F)* (rho2 * mu1) / (rho1 * mu2)
 */
void VOF_update_density_viscosity (Cart3d_bag *data_bag);




//Below are all the low-level functions implemented from Basilisk's geometry.h file used to compute the interface normals and to reconstruct the interface (PLIC)


// Basic functions
double clampDouble(double val, double lower, double upper);
double minDouble(double a, double b);
double maxDouble(double a, double b);

// Function declarations
double line_alpha_1D(double c, PointType n);
double line_alpha_2D(double c, PointType n);
double plane_alpha_3D(double c, PointType n, double *alpha_out);

double line_area_1D(double nx, double alpha);
double line_area_2D(double nx, double ny, double alpha);
double plane_volume_3D(PointType n, double alpha);

double rectangle_fraction(PointType n, double alpha, PointType a, PointType b);

int facets_2D(PointType n, double alpha, PointType p[2]);
int facets_3D(PointType n, double alpha, PointType v[12], double h);

void line_center(PointType m, double alpha, double a, PointType *p);
void plane_center(PointType m, double alpha, double a, PointType *p);

//Below are all the low-level functions from Basilisk's myc2d.h and myc.h files.

/**
 * Computes a 2D interface normal using the Mixed Youngs and Central (MYC) scheme.
 *
 * @param F A 2D array of volume fractions.
 * @param i Index of the cell in the x-direction.
 * @param j Index of the cell in the y-direction.
 * @return A PointType representing the 2D normal (n.x, n.y). The z-component is 0.
 */
PointType mycs2D(double **F, int i, int j);

/**
 * Computes a 3D interface normal using the Mixed Youngs and Central (MYC) scheme.
 *
 * @param F A 3D array of volume fractions.
 * @param i Index of the cell in the x-direction.
 * @param j Index of the cell in the y-direction.
 * @param k Index of the cell in the z-direction.
 * @return A PointType representing the 3D normal (n.x, n.y, n.z).
 */
PointType mycs3D(double ***F, int i, int j, int k);


//Below are all the low-level functions implemented from Basilisk's fractions.h file used to compute the interface normals and to reconstruct the interface (PLIC)

/**
 * Computes the 2D volume and surface fractions based on a level set function.
 *
 * @param phi The level set field (2D array).
 * @param Nx Number of grid cells in x-direction.
 * @param Ny Number of grid cells in y-direction.
 * @param val The level set threshold value.
 * @param c The output volume fraction array (2D).
 * @param sx Optional surface fraction array in x-direction (2D).
 * @param sy Optional surface fraction array in y-direction (2D).
 */
void fractionsLevelSet2D(double **Phi,    // levelset at vertices
                         double **c,      // volume fraction at cells
                         double val,      // interface threshold
                         int Nx, int Ny);

/**
 * Computes the 3D volume and surface fractions based on a level set function.
 *
 * @param phi The level set field (3D array).
 * @param Nx Number of grid cells in x-direction.
 * @param Ny Number of grid cells in y-direction.
 * @param Nz Number of grid cells in z-direction.
 * @param val The level set threshold value.
 * @param c The output volume fraction array (3D).
 * @param sx Optional surface fraction array in x-direction (3D).
 * @param sy Optional surface fraction array in y-direction (3D).
 * @param sz Optional surface fraction array in z-direction (3D).
 */
void fractionsLevelSet3D(double ***Phi,    // vertex-based levelset
                         double ***c,      // cell-based volume fraction
                         double ***sx,     // face fraction in x-direction or NULL
                         double ***sy,     // face fraction in y-direction or NULL
                         double ***sz,     // face fraction in z-direction or NULL
                         double val, 
                         int Nx, int Ny, int Nz);

/**
 * Computes the interface normal vector in 3D using face fractions.
 *
 * @param c       Volume fraction array (3D).
 * @param sx      Surface fraction in the x-direction (Nx+1)xNyxNz array.
 * @param sy      Surface fraction in the y-direction (NxxNy+1xNz array).
 * @param sz      Surface fraction in the z-direction (NxxNyxNz+1 array).
 * @param i, j, k Indices of the cell for which the normal is to be computed.
 * @param data_bag Pointer to Cart3d_bag for grid and parameters.
 * @return A PointType storing the normalized interface normal vector (n.x, n.y, n.z).
 ******************************************************************************/
PointType VoF_facet_normal_3D(
    double ***c,
    double ***sx,
    double ***sy,
    double ***sz,
    int i, 
    int j, 
    int k,
    Cart3d_bag *data_bag);

    
/**
 * VOF_InterfaceArea_3D
 *
 * Computes the surface area of the interface in 3D. This is analogous to
 * Basilisk's "interface_area()" function but specialized to a 3D uniform-grid
 * PARTIES code. 
 *
 * This does NOT do boundary checks for ghost layers; we assume your code 
 * has them or is guaranteed in-bounds for i±1, j±1, k±1.
 ******************************************************************************/
double VOF_InterfaceArea_3D(Cart3d_bag *data_bag);


#endif // VOLUMEFRACTION_H



