#include "DataTypes.h"
#include "Boundary.h"
#include "Communication.h" // For ghost-node exchanges
#include <stdio.h>
#include <math.h>


/******************************************************************************/
// A simple dot-product function for 3D arrays with parallel sum.
// Then we do an MPI_Allreduce to get the global sum.
/******************************************************************************/
static double Pressure_innerProd(double ***vec1, double ***vec2,  MAC_grid *grid, Parameters *params) {

	int i, j, k;

	double partialSum, totalSum;

	int NX = grid -> NX;
	int NY = grid -> NY;
	int NZ = grid -> NZ;

	// Processor start and end indices
	int Is = grid -> G_Is;
	int Js = grid -> G_Js;
	int Ks = grid -> G_Ks;

	int Ie = min( grid->G_Ie, NX-1);  // No ghost nodes
	int Je = min( grid->G_Je, NY-1);
	int Ke = min( grid->G_Ke, NZ-1);


	// Calculate inner product of local nodes
	partialSum = 0;
	for (k = Ks; k < Ke; k++) {
		for (j = Js; j < Je; j++) {
			for (i = Is; i < Ie; i++) {

				partialSum += vec1[k][j][i] * vec2[k][j][i];

			}
		}
	}

	// Communicate and sum inner product among all processors
	MPI_Allreduce(&partialSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, PCW);

	return totalSum;
}

/******************************************************************************/
// Conjugate Gradient solver for the variable-coefficient Poisson system:
//
//     ∇ · ( (1/ρ) ∇φ ) = rhs
//
// The unknown φ is stored in p->deltap. The right-hand side is p->rhs.
// The residual, direction, and A·d are p->res, p->d, and p->Ad, respectively.
//
// We call Pressure_operator_variableCoeff() to apply the operator A to a vector.
//
// At the end, p->deltap holds the solution φ. 
//
// Return value is the number of iterations performed.
//
/******************************************************************************/
int Pressure_solve_cg(Cart3d_bag *data_bag)
{
    Pressure   *p      = data_bag->p;
    MAC_grid   *grid   = data_bag->grid;
    Parameters *params = data_bag->params;
    VolumeFraction *vof = data_bag->vof;

    // VoF varying density
    double ***rho = vof->rho;  


    // Some references to the arrays used for CG:
    double ***phi = p->deltap; // Our unknown "pressure correction"
    double ***rhs = p->rhs;    // The Poisson RHS
    double ***res = p->res;    // CG residual
    double ***d   = p->d;      // CG direction
    double ***Ad  = p->Ad;     // Operator(A)*d

    // Convergence criteria
    double tolerance = params->CG_ETOL;   // e.g. 1e-6
    int maxIters     = params->CG_MAXIT;  // e.g. 2000

    // We'll store dot products here
    double rr, rr_new, alpha, beta;

    // (1) Build initial residual: r = rhs - A*phi
    //     A*phi is computed by calling Pressure_operator_variableCoeff(...)
    Pressure_operator_variableCoeff(Ad, phi, rho, grid, params, data_bag);
    // r = rhs - Ad
    {
        int i, j, k;
        int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
        int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;
        for (k = Ks; k < Ke; k++) {
            for (j = Js; j < Je; j++) {
                for (i = Is; i < Ie; i++) {
                    res[k][j][i] = rhs[k][j][i] - Ad[k][j][i];
                    d[k][j][i]   = res[k][j][i];
                }
            }
        }
    }

    // (2) Initial norm
    rr = Pressure_innerProd(res, res, grid, params);
    double initRes = sqrt(rr);
    if (params->rank == 0) {
        printf(" Pressure CG solver start: initial residual = %g\n", initRes);
    }

    // If already small, we can stop
    if (initRes < tolerance) {
        return 0;
    }

    // (3) Main CG iteration
    int iters;
    for (iters = 1; iters <= maxIters; iters++) {

        // Ad = A*d
        Pressure_operator_variableCoeff(Ad, d, rho, grid, params, data_bag);

        // dot(d, Ad)
        double dAd = Pressure_innerProd(d, Ad, grid, params);

        // alpha = (r,r) / (d,Ad)
        alpha = rr / (dAd + 1e-30); // small safeguard to avoid /0

        // Update: phi = phi + alpha*d   AND   r = r - alpha*Ad
        {
            int i, j, k;
            int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
            int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    for (i = Is; i < Ie; i++) {
                        phi[k][j][i] += alpha * d[k][j][i];
                        res[k][j][i] -= alpha * Ad[k][j][i];
                    }
                }
            }
        }

        // New residual norm
        rr_new = Pressure_innerProd(res, res, grid, params);
        double resNorm = sqrt(rr_new);

        // Check convergence
        if (resNorm < tolerance) {
            if (params->rank == 0) {
                printf(" Pressure CG converged after %d iterations, resNorm=%g\n", iters, resNorm);
            }
            break;
        }

        // beta = rr_new / rr
        beta = rr_new / (rr + 1e-30);

        // d = r + beta*d
        {
            int i, j, k;
            int Is = grid->G_Is, Js = grid->G_Js, Ks = grid->G_Ks;
            int Ie = grid->G_Ie, Je = grid->G_Je, Ke = grid->G_Ke;
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    for (i = Is; i < Ie; i++) {
                        d[k][j][i] = res[k][j][i] + beta * d[k][j][i];
                    }
                }
            }
        }

        // Prepare for next iteration
        rr = rr_new;

        // Optionally print every e.g. 50 iterations
        if (params->rank == 0 && (iters % 50 == 0)) {
            printf(" Pressure CG iteration %d, residual=%g\n", iters, resNorm);
        }
    } // end for

    // iters is either the iteration we ended on, or maxIters+1 if we never broke
    if (iters > maxIters) {
        double finalRes = sqrt(rr);
        if (params->rank == 0) {
            printf(" Pressure CG reached maxIters=%d with final resNorm=%g\n", maxIters, finalRes);
        }
        iters = maxIters;
    }

    return iters;
}



/******************************************************************************/
// Pressure_operator_variableCoeff
// Computes: Aphi = ∇·( 1/rho * ∇phi )
/******************************************************************************/
void Pressure_operator_variableCoeff(
    double ***Aphi,         // output: operator(A) * phi
    double ***phi,          // input: phi array
    double ***rho,       // rho at cell centers
    MAC_grid *grid,
    Parameters *params,
    Cart3d_bag *data_bag )
{
    // 1) Exchange ghost cells and apply BCs to phi.
    Communication_update_ghost_nodes_flow_variable(phi, 'h', params->ghost_nodes, data_bag);
    Pressure_apply_BCs(phi, grid, params); 
    

    // 2) Retrieve grid size info. For uniform grid, just do:
    int i, j, k;
    int NX = grid->NX;
    int NY = grid->NY;
    int NZ = grid->NZ;

    // Start and end indices on this processor (exclude last ghost cell)
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;
    int Ie = min(grid->G_Ie, NX-1);
    int Je = min(grid->G_Je, NY-1);
    int Ke = min(grid->G_Ke, NZ-1);

    // Grid spacing squared (cell-centered)
    double idx2 = grid->idx_c[1] * grid->idx_c[1];
    double idy2 = grid->idy_c[1] * grid->idy_c[1];
    double idz2 = grid->idz_c[1] * grid->idz_c[1];

    // 3) Loop over interior cells and compute the 3D divergence of (invRho * grad(phi)).
    //    Face-based approach. Cell-centered version with face averaging.
 
    for (k = Ks; k < Ke; k++) {
        for (j = Js; j < Je; j++) {
            for (i = Is; i < Ie; i++) {

                //----- X direction -----
                // East face (i+1/2)
                double axp = 2.0 / (rho[k][j][i] + rho[k][j][i+1]) * (phi[k][j][i+1] - phi[k][j][i]);
                // West face (i-1/2)
                double axm = 2.0 / (rho[k][j][i-1] + rho[k][j][i]) * (phi[k][j][i]   - phi[k][j][i-1]);
                double Ax  = (axp - axm) * idx2;

                //----- Y direction -----
                // North face (j+1/2)
                double ayp = 2.0 / (rho[k][j][i] + rho[k][j+1][i]) * (phi[k][j+1][i] - phi[k][j][i]);
                // South face (j-1/2)
                double aym = 2.0 / (rho[k][j-1][i] + rho[k][j][i]) * (phi[k][j][i]   - phi[k][j-1][i]);
                double Ay  = (ayp - aym) * idy2;

                //----- Z direction -----
                // Front face (k+1/2)
                double azp = 2.0 / (rho[k][j][i] + rho[k+1][j][i]) * (phi[k+1][j][i] - phi[k][j][i]);
                // Back face (k-1/2)
                double azm = 2.0 / (rho[k-1][j][i] + rho[k][j][i]) * (phi[k][j][i]   - phi[k-1][j][i]);
                double Az  = (azp - azm) * idz2;

                Aphi[k][j][i] = Ax + Ay + Az;
            }
        }
    }

}



/******************************************************************************/
// Applies physical boundary conditions for the pressure‐correction array `phi`.
//
//
// Typically for no‐slip or free‐slip velocity boundaries, we do Neumann(∂phi/∂n=0).
// For inflow/outflow, we might impose phi=0 or also do Neumann, depending on the
// solver's approach.
//
// NOTE: We assume one layer of ghost cells. If you have more, adapt loops accordingly.
/******************************************************************************/
void Pressure_apply_BCs(double ***phi, MAC_grid *grid, Parameters *params)
{
    int i, j, k;

    // For convenience
    int NX = grid->NX; // interior cell count in x
    int NY = grid->NY; // interior cell count in y
    int NZ = grid->NZ; // interior cell count in z

    // Local domain corners for interior cells (excluding ghost cells):
    // For ghost cells in x, we might look at i == G_Is-1 or i == G_Ie, etc.
    int Is = grid->G_Is;
    int Ie = grid->G_Ie;
    int Js = grid->G_Js;
    int Je = grid->G_Je;
    int Ks = grid->G_Ks;
    int Ke = grid->G_Ke;

    // -------------------------------------------------------------------------
    // X-DIRECTION BOUNDARIES
    // -------------------------------------------------------------------------
#ifndef XPERIODIC
    // Check if this processor is at the left domain boundary (x=0).
    if (Is == 0) {
        // We have cells i = 0..(Ie-1). The ghost is i = -1.
        // E.g. we can impose ∂phi/∂x=0 => phi(k,j,-1) = phi(k,j,0).

        #if defined(LEFT_WALL_VELOCITY_NOSLIP) || defined(LEFT_WALL_VELOCITY_FREESLIP)
            // Typical approach: Neumann => copy interior
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    phi[k][j][Is - 1] = phi[k][j][Is]; 
                }
            }
        #elif defined(LEFT_INFLOW)
            // Some codes prefer Dirichlet = 0
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    phi[k][j][Is - 1] = 0.0;
                }
            }
        #else
            // fallback: Neumann
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    phi[k][j][Is - 1] = phi[k][j][Is];
                }
            }
        #endif
    }

    // Check if this processor is at the right domain boundary (x=NX-1).
    if (Ie == NX) {
        // The ghost is i = NX.
        #if defined(RIGHT_WALL_VELOCITY_NOSLIP) || defined(RIGHT_WALL_VELOCITY_FREESLIP)
            // Neumann => copy interior
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    phi[k][j][Ie] = phi[k][j][Ie - 1];
                }
            }
        #elif defined(RIGHT_INFLOW)
            // Dirichlet => 0
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    phi[k][j][Ie] = 0.0;
                }
            }
        #elif defined(RIGHT_OUTFLOW)
            // Typically outflow => zero derivative
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    phi[k][j][Ie] = phi[k][j][Ie - 1];
                }
            }
        #else
            // fallback: Neumann
            for (k = Ks; k < Ke; k++) {
                for (j = Js; j < Je; j++) {
                    phi[k][j][Ie] = phi[k][j][Ie - 1];
                }
            }
        #endif
    }
#endif // !XPERIODIC

    // -------------------------------------------------------------------------
    // Y-DIRECTION BOUNDARIES
    // -------------------------------------------------------------------------
#ifndef YPERIODIC
    // Bottom boundary (y=0)
    if (Js == 0) {
        // ghost cell at j=Js-1
        #if defined(BOTTOM_WALL_VELOCITY_NOSLIP) || defined(BOTTOM_WALL_VELOCITY_FREESLIP)
            for (k = Ks; k < Ke; k++) {
                for (i = Is; i < Ie; i++) {
                    phi[k][Js - 1][i] = phi[k][Js][i];
                }
            }
        #else
            // fallback: same approach => Neumann
            for (k = Ks; k < Ke; k++) {
                for (i = Is; i < Ie; i++) {
                    phi[k][Js - 1][i] = phi[k][Js][i];
                }
            }
        #endif
    }

    // Top boundary (y=NY-1)
    if (Je == NY) {
        // ghost cell at j=Je
        #if defined(TOP_WALL_VELOCITY_NOSLIP) || defined(TOP_WALL_VELOCITY_FREESLIP)
            for (k = Ks; k < Ke; k++) {
                for (i = Is; i < Ie; i++) {
                    phi[k][Je][i] = phi[k][Je - 1][i];
                }
            }
        #else
            // fallback: Neumann
            for (k = Ks; k < Ke; k++) {
                for (i = Is; i < Ie; i++) {
                    phi[k][Je][i] = phi[k][Je - 1][i];
                }
            }
        #endif
    }
#endif // !YPERIODIC

    // -------------------------------------------------------------------------
    // Z-DIRECTION BOUNDARIES
    // -------------------------------------------------------------------------
#ifndef ZPERIODIC
    // Back boundary (z=0)
    if (Ks == 0) {
        // ghost cell k=Ks-1
        #if defined(BACK_WALL_VELOCITY_NOSLIP) || defined(BACK_WALL_VELOCITY_FREESLIP)
            for (j = Js; j < Je; j++) {
                for (i = Is; i < Ie; i++) {
                    phi[Ks - 1][j][i] = phi[Ks][j][i];
                }
            }
        #else
            // fallback: Neumann
            for (j = Js; j < Je; j++) {
                for (i = Is; i < Ie; i++) {
                    phi[Ks - 1][j][i] = phi[Ks][j][i];
                }
            }
        #endif
    }

    // Front boundary (z=NZ-1)
    if (Ke == NZ) {
        // ghost cell k=Ke
        #if defined(FRONT_WALL_VELOCITY_NOSLIP) || defined(FRONT_WALL_VELOCITY_FREESLIP)
            for (j = Js; j < Je; j++) {
                for (i = Is; i < Ie; i++) {
                    phi[Ke][j][i] = phi[Ke - 1][j][i];
                }
            }
        #else
            // fallback: Neumann
            for (j = Js; j < Je; j++) {
                for (i = Is; i < Ie; i++) {
                    phi[Ke][j][i] = phi[Ke - 1][j][i];
                }
            }
        #endif
    }
#endif // !ZPERIODIC
}



