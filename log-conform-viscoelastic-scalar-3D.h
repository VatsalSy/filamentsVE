/** Title: log-conform-viscoelastic-3D.h
# Version: 1.1
# Main feature 1: A exists in across the domain and relaxes according to \lambda. The stress only acts according to G.
# Main feature 2: This is the 3D implementation of [log-conform-viscoelastic-scalar-2D.h](log-conform-viscoelastic-scalar-2D.h).

# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Updated: Oct 20, 2024

# change log: Oct 19, 2024 (v1.0)
- 3D implementation
- scalar implementation

# change log: Oct 20, 2024 (v1.1)
- Added a check for negative eigenvalues. If any are found, print the location and value of the offending eigenvalue.
- Please report this bug by opening an issue on the GitHub repository. 
- The code works!!! :) 

*/

/** The code is same as http://basilisk.fr/src/log-conform.h but 
- written with G-\lambda formulation. 
- It also fixes the bug where [\sigma_p] = 0 & [\sigma_s] = \gamma\kappa instead of [\sigma_s+\sigma_p] = \gamma\kappa.
*/ 

/**
 * # TODO: (non-critical, non-urgent)
 * axi compatibility is not there. This will not be fixed. To use axi, please use: [log-conform-viscoelastic-scalar-2D.h](log-conform-viscoelastic-scalar-2D.h) for a scalar formulation, or better yet, use [log-conform-viscoelastic.h](log-conform-viscoelastic.h) which is more efficient.
*/

#if AXI
#error "axi compatibility is not there. To keep the code easy to read, we will not implement axi compatibility just yet."
#endif

#include "bcg.h"

/*
conformation tensor */
// diagonal elements
scalar A11[], A22[], A33[];
// off-diagonal elements
scalar A12[], A13[], A23[];
/*
stress tensor */
// diagonal elements
scalar T11[], T22[], T33[];
// off-diagonal elements
scalar T12[], T13[], T23[];

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;

  /*
  initialize A and T
  */
  for (scalar s in {A11, A22, A33}) {
    foreach () {
      s[] = 1.;
    }
  }
  for (scalar s in {T11, T12, T13, T22, T23, T33, A12, A13, A23}) {
    foreach(){
      s[] = 0.;
    }
  }

  for (scalar s in {A11, A22, A33, T11, T22, T33, A12, A13, A23, T12, T13, T23}) {
    if (s.boundary[left] != periodic_bc) {
      s[left] = neumann(0);
	    s[right] = neumann(0);
    }
    if (s.boundary[top] != periodic_bc) {
      s[top] = neumann(0);
	    s[bottom] = neumann(0);
    }
#if dimension == 3
    if (s.boundary[front] != periodic_bc) {
      s[front] = neumann(0);
	    s[back] = neumann(0);
    }
#endif
  }
}

/**
## Useful functions in 2D

The first step is to implement a routine to calculate the eigenvalues
and eigenvectors of the conformation tensor $\mathbf{A}$.

These structs ressemble Basilisk vectors and tensors but are just
arrays not related to the grid. */

#if dimension == 2
typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;

static void diagonalization_2D (pseudo_v * Lambda, pseudo_t * R, pseudo_t * A)
{
  /**
  The eigenvalues are saved in vector $\Lambda$ computed from the
  trace and the determinant of the symmetric conformation tensor
  $\mathbf{A}$. */

  if (sq(A->x.y) < 1e-15) {
    R->x.x = R->y.y = 1.;
    R->y.x = R->x.y = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y;
    return;
  }

  double T = A->x.x + A->y.y; // Trace of the tensor
  double D = A->x.x*A->y.y - sq(A->x.y); // Determinant

  /**
  The eigenvectors, $\mathbf{v}_i$ are saved by columns in tensor
  $\mathbf{R} = (\mathbf{v}_1|\mathbf{v}_2)$. */

  R->x.x = R->x.y = A->x.y;
  R->y.x = R->y.y = -A->x.x;
  double s = 1.;
  for (int i = 0; i < dimension; i++) {
    double * ev = (double *) Lambda;
    ev[i] = T/2 + s*sqrt(sq(T)/4. - D);
    s *= -1;
    double * Rx = (double *) &R->x;
    double * Ry = (double *) &R->y;
    Ry[i] += ev[i];
    double mod = sqrt(sq(Rx[i]) + sq(Ry[i]));
    Rx[i] /= mod;
    Ry[i] /= mod;
  }
}
#endif

/*
Now this is the 3D implementation.
*/
#if dimension == 3

#include "eigen_decomposition.h"

typedef struct { double x, y, z; }   pseudo_v3d;
typedef struct { pseudo_v3d x, y, z; } pseudo_t3d;

static void diagonalization_3D (pseudo_v3d * Lambda, pseudo_t3d * R, pseudo_t3d * A)
{
  // Check if the matrix is already diagonal
  if (sq(A->x.y) + sq(A->x.z) + sq(A->y.z) < 1e-15) {
    R->x.x = R->y.y = R->z.z = 1.;
    R->y.x = R->x.y = R->z.x = R->x.z = R->z.y = R->y.z = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y; Lambda->z = A->z.z;
    return;
  }

  // Compute eigenvalues using the eigen_decomposition function
  double matrix[3][3] = {
    {A->x.x, A->x.y, A->x.z},
    {A->y.x, A->y.y, A->y.z},
    {A->z.x, A->z.y, A->z.z}
  };
  double eigenvectors[3][3];
  double eigenvalues[3];
  
  compute_eigensystem_symmetric_3x3(matrix, eigenvectors, eigenvalues);

  // Store eigenvalues and eigenvectors
  Lambda->x = eigenvalues[0];
  Lambda->y = eigenvalues[1];
  Lambda->z = eigenvalues[2];

  R->x.x = eigenvectors[0][0]; R->x.y = eigenvectors[0][1]; R->x.z = eigenvectors[0][2];
  R->y.x = eigenvectors[1][0]; R->y.y = eigenvectors[1][1]; R->y.z = eigenvectors[1][2];
  R->z.x = eigenvectors[2][0]; R->z.y = eigenvectors[2][1]; R->z.z = eigenvectors[2][2];

}
#endif

/**
The stress tensor depends on previous instants and has to be
integrated in time. In the log-conformation scheme the advection of
the stress tensor is circumvented, instead the conformation tensor,
$\mathbf{A}$ (or more precisely the related variable $\Psi$) is
advanced in time.

In what follows we will adopt a scheme similar to that of [Hao \& Pan
(2007)](#hao2007). We use a split scheme, solving successively

a) the upper convective term:
$$
\partial_t \Psi = 2 \mathbf{B} + (\Omega \cdot \Psi -\Psi \cdot \Omega)
$$
b) the advection term:
$$
\partial_t \Psi + \nabla \cdot (\Psi \mathbf{u}) = 0
$$
c) the model term (but set in terms of the conformation 
tensor $\mathbf{A}$). In an Oldroyd-B viscoelastic fluid, the model is
$$ 
\partial_t \mathbf{A} = -\frac{\mathbf{f}_r (\mathbf{A})}{\lambda}
$$

The implementation below assumes that the values of $\Psi$ and
$\conform_p$ are never needed simultaneously. This means that $\conform_p$ can
be used to store (temporarily) the values of $\Psi$ (i.e. $\Psi$ is
just an alias for $\conform_p$). */

#if dimension == 2
event tracer_advection(i++)
{
  scalar Psi11 = A11;
  scalar Psi12 = A12;
  scalar Psi22 = A22;

  /**
  ### Computation of $\Psi = \log \mathbf{A}$ and upper convective term */

  foreach() {
    /**
      We assume that the stress tensor $\mathbf{\tau}_p$ depends on the
      conformation tensor $\mathbf{A}$ as follows
      $$
      \mathbf{\tau}_p = G_p (\mathbf{A}) =
      G_p (\mathbf{A} - I)
      $$
    */

    pseudo_t A;

    A.x.x = A11[];
    A.x.y = A12[];
    A.y.y = A22[];

    /**
    The conformation tensor is diagonalized through the
    eigenvector tensor $\mathbf{R}$ and the eigenvalues diagonal
    tensor, $\Lambda$. */

    pseudo_v Lambda;
    pseudo_t R;
    diagonalization_2D (&Lambda, &R, &A);
    
    /**
    $\Psi = \log \mathbf{A}$ is easily obtained after diagonalization, 
    $\Psi = R \cdot \log(\Lambda) \cdot R^T$. */
    
    Psi12[] = R.x.x*R.y.x*log(Lambda.x) + R.y.y*R.x.y*log(Lambda.y);
    Psi11[] = sq(R.x.x)*log(Lambda.x) + sq(R.x.y)*log(Lambda.y);
    Psi22[] = sq(R.y.y)*log(Lambda.y) + sq(R.y.x)*log(Lambda.x);

    /**
    We now compute the upper convective term $2 \mathbf{B} +
    (\Omega \cdot \Psi -\Psi \cdot \Omega)$.

    The diagonalization will be applied to the velocity gradient
    $(\nabla u)^T$ to obtain the antisymmetric tensor $\Omega$ and
    the traceless, symmetric tensor, $\mathbf{B}$. If the conformation
    tensor is $\mathbf{I}$, $\Omega = 0$ and $\mathbf{B}= \mathbf{D}$.  */

    pseudo_t B;
    double OM = 0.;
    if (fabs(Lambda.x - Lambda.y) <= 1e-20) {
      B.x.y = (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(4.*Delta); 
      foreach_dimension() 
        B.x.x = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    } else {
      pseudo_t M;
      foreach_dimension() {
        M.x.x = (sq(R.x.x)*(u.x[1] - u.x[-1]) +
        sq(R.y.x)*(u.y[0,1] - u.y[0,-1]) +
        R.x.x*R.y.x*(u.x[0,1] - u.x[0,-1] + 
        u.y[1] - u.y[-1]))/(2.*Delta);
        M.x.y = (R.x.x*R.x.y*(u.x[1] - u.x[-1]) + 
        R.x.y*R.y.x*(u.y[1] - u.y[-1]) +
        R.x.x*R.y.y*(u.x[0,1] - u.x[0,-1]) +
        R.y.x*R.y.y*(u.y[0,1] - u.y[0,-1]))/(2.*Delta);
      }
      double omega = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
      OM = (R.x.x*R.y.y - R.x.y*R.y.x)*omega;

      B.x.y = M.x.x*R.x.x*R.y.x + M.y.y*R.y.y*R.x.y;
      foreach_dimension()
        B.x.x = M.x.x*sq(R.x.x)+M.y.y*sq(R.x.y);	
    }

    /**
    We now advance $\Psi$ in time, adding the upper convective
    contribution. */

    double s = -Psi12[];
    Psi12[] += dt * (2. * B.x.y + OM * (Psi22[] - Psi11[]));
    s *= -1;
    Psi11[] += dt * 2. * (B.x.x + s * OM);
    s *= -1;
    Psi22[] += dt * 2. * (B.y.y + s * OM);
  }

  /**
  ### Advection of $\Psi$
 
  We proceed with step (b), the advection of the log of the
  conformation tensor $\Psi$. */

  advection ({Psi11, Psi12, Psi22}, uf, dt);

  /**
  ### Convert back to \conform_p */

  foreach() {
    /**
    It is time to undo the log-conformation, again by
    diagonalization, to recover the conformation tensor $\mathbf{A}$
    and to perform step (c).*/

    pseudo_t A = {{Psi11[], Psi12[]}, {Psi12[], Psi22[]}}, R;
    pseudo_v Lambda;
    diagonalization_2D (&Lambda, &R, &A);
    Lambda.x = exp(Lambda.x), Lambda.y = exp(Lambda.y);
    
    A.x.y = R.x.x*R.y.x*Lambda.x + R.y.y*R.x.y*Lambda.y;
    foreach_dimension()
      A.x.x = sq(R.x.x)*Lambda.x + sq(R.x.y)*Lambda.y;

    /**
    We perform now step (c) by integrating 
    $\mathbf{A}_t = -\mathbf{f}_r (\mathbf{A})/\lambda$ to obtain
    $\mathbf{A}^{n+1}$. This step is analytic,
    $$
    \int_{t^n}^{t^{n+1}}\frac{d \mathbf{A}}{\mathbf{I}- \mathbf{A}} = 
    \frac{\Delta t}{\lambda}
    $$
    */

    double intFactor = lambda[] != 0. ? exp(-dt/lambda[]): 0.;
    
    A.x.y *= intFactor;
    foreach_dimension()
      A.x.x = (1. - intFactor) + A.x.x*intFactor;

    /**
      Then the Conformation tensor $\mathcal{A}_p^{n+1}$ is restored from
      $\mathbf{A}^{n+1}$.  */
    
    A12[] = A.x.y;
    T12[] = Gp[]*A.x.y;
    A11[] = A.x.x;
    T11[] = Gp[]*(A.x.x - 1.);
    A22[] = A.y.y;
    T22[] = Gp[]*(A.y.y - 1.);
  }
}

#elif dimension == 3

event tracer_advection(i++)
{
  /**
  ### Computation of $\Psi = \log \mathbf{A}$ and upper convective term */

  // start by declaring the scalar variables that will store the components of $\Psi$
  scalar Psi11 = A11, Psi12 = A12, Psi13 = A13,
         Psi22 = A22, Psi23 = A23, Psi33 = A33;

  foreach() {
    pseudo_t3d A, R;
    pseudo_v3d Lambda;

    A.x.x = A11[]; A.x.y = A12[]; A.x.z = A13[];
    A.y.x = A12[]; A.y.y = A22[]; A.y.z = A23[];
    A.z.x = A13[]; A.z.y = A23[]; A.z.z = A33[];

    // Diagonalize the conformation tensor A to obtain the eigenvalues Lambda and eigenvectors R
    diagonalization_3D (&Lambda, &R, &A);

    /*
    Check for negative eigenvalues -- this should never happen. If it does, print the location and value of the offending eigenvalue.
    Please report this bug by opening an issue on the GitHub repository. 
    */
    if (Lambda.x <= 0. || Lambda.y <= 0. || Lambda.z <= 0.) {
      fprintf(ferr, "Negative eigenvalue detected: Lambda.x = %g, Lambda.y = %g, Lambda.z = %g\n", Lambda.x, Lambda.y, Lambda.z);
      fprintf(ferr, "x = %g, y = %g, z = %g, f = %g\n", x, y, z, f[]);
      exit(1);
    }
    
    // Compute Psi = log(A)
    Psi11[] = R.x.x*R.x.x*log(Lambda.x) + R.y.x*R.y.x*log(Lambda.y) + R.z.x*R.z.x*log(Lambda.z);
    Psi22[] = R.x.y*R.x.y*log(Lambda.x) + R.y.y*R.y.y*log(Lambda.y) + R.z.y*R.z.y*log(Lambda.z);
    Psi33[] = R.x.z*R.x.z*log(Lambda.x) + R.y.z*R.y.z*log(Lambda.y) + R.z.z*R.z.z*log(Lambda.z);

    Psi12[] = R.x.x*R.x.y*log(Lambda.x) + R.y.x*R.y.y*log(Lambda.y) + R.z.x*R.z.y*log(Lambda.z);
    Psi13[] = R.x.x*R.x.z*log(Lambda.x) + R.y.x*R.y.z*log(Lambda.y) + R.z.x*R.z.z*log(Lambda.z);
    Psi23[] = R.x.y*R.x.z*log(Lambda.x) + R.y.y*R.y.z*log(Lambda.y) + R.z.y*R.z.z*log(Lambda.z);

    // Compute B and Omega tensors (3D version)
    pseudo_t3d B, M, Omega;

    // Check if any pair of eigenvalues are numerically equal (within a small tolerance)
    if (fabs(Lambda.x - Lambda.y) <= 1e-20 || fabs(Lambda.y - Lambda.z) <= 1e-20 || fabs(Lambda.z - Lambda.x) <= 1e-20) {
      // In case of equal eigenvalues, the calculations for B and Omega simplify significantly

      // Compute off-diagonal elements of B using central differences
      // These represent the symmetric part of the velocity gradient tensor
      B.x.y = (u.y[1,0,0] - u.y[-1,0,0] + u.x[0,1,0] - u.x[0,-1,0])/(4.*Delta);  // (dv/dx + du/dy)/2
      B.x.z = (u.z[1,0,0] - u.z[-1,0,0] + u.x[0,0,1] - u.x[0,0,-1])/(4.*Delta);  // (dw/dx + du/dz)/2
      B.y.z = (u.z[0,1,0] - u.z[0,-1,0] + u.y[0,0,1] - u.y[0,0,-1])/(4.*Delta);  // (dw/dy + dv/dz)/2

      // Compute diagonal elements of B
      // These represent the normal strain rates
      B.x.x = (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta);  // du/dx
      B.y.y = (u.y[0,1,0] - u.y[0,-1,0])/(2.*Delta);  // dv/dy
      B.z.z = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);  // dw/dz

      // Set all components of Omega to zero
      // This is because Omega represents the antisymmetric part of the velocity gradient tensor,
      // which vanishes when eigenvalues are equal
      Omega.x.y = Omega.x.z = Omega.y.z = Omega.y.x = Omega.z.x = Omega.z.y = 0.;

    } else {
      
      /*
      ### Compute the velocity gradient tensor components using central differences
      - These represent the spatial derivatives of each velocity component
      - These gradients form the velocity gradient tensor (nablaU):
      [ dudx  dudy  dudz ]
      [ dvdx  dvdy  dvdz ]
      [ dwdx  dwdy  dwdz ]
      */

      // Derivatives of u (x-component of velocity)
      double dudx = (u.x[1,0,0] - u.x[-1,0,0])/(2.0*Delta);  // du/dx
      double dudy = (u.x[0,1,0] - u.x[0,-1,0])/(2.0*Delta);  // du/dy
      double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.0*Delta);  // du/dz

      // Derivatives of v (y-component of velocity)
      double dvdx = (u.y[1,0,0] - u.y[-1,0,0])/(2.0*Delta);  // dv/dx
      double dvdy = (u.y[0,1,0] - u.y[0,-1,0])/(2.0*Delta);  // dv/dy
      double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.0*Delta);  // dv/dz

      // Derivatives of w (z-component of velocity)
      double dwdx = (u.z[1,0,0] - u.z[-1,0,0])/(2.0*Delta);  // dw/dx
      double dwdy = (u.z[0,1,0] - u.z[0,-1,0])/(2.0*Delta);  // dw/dy
      double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.0*Delta);  // dw/dz

      /*
      Calculate M tensor: M = R^T * (nablaU)^T * R
      - This transforms the velocity gradient tensor to the eigenvector basis
      */

      // First, compute intermediate products of R^T and (nablaU)^T
      double Rx_gradU_x = R.x.x*dudx + R.x.y*dvdx + R.x.z*dwdx;
      double Rx_gradU_y = R.x.x*dudy + R.x.y*dvdy + R.x.z*dwdy;
      double Rx_gradU_z = R.x.x*dudz + R.x.y*dvdz + R.x.z*dwdz;

      double Ry_gradU_x = R.y.x*dudx + R.y.y*dvdx + R.y.z*dwdx;
      double Ry_gradU_y = R.y.x*dudy + R.y.y*dvdy + R.y.z*dwdy;
      double Ry_gradU_z = R.y.x*dudz + R.y.y*dvdz + R.y.z*dwdz;

      double Rz_gradU_x = R.z.x*dudx + R.z.y*dvdx + R.z.z*dwdx;
      double Rz_gradU_y = R.z.x*dudy + R.z.y*dvdy + R.z.z*dwdy;
      double Rz_gradU_z = R.z.x*dudz + R.z.y*dvdz + R.z.z*dwdz;

      // Now compute M components by multiplying the intermediate products with R
      M.x.x = R.x.x*Rx_gradU_x + R.x.y*Rx_gradU_y + R.x.z*Rx_gradU_z;
      M.x.y = R.x.x*Ry_gradU_x + R.x.y*Ry_gradU_y + R.x.z*Ry_gradU_z;
      M.x.z = R.x.x*Rz_gradU_x + R.x.y*Rz_gradU_y + R.x.z*Rz_gradU_z;

      M.y.x = R.y.x*Rx_gradU_x + R.y.y*Rx_gradU_y + R.y.z*Rx_gradU_z;
      M.y.y = R.y.x*Ry_gradU_x + R.y.y*Ry_gradU_y + R.y.z*Ry_gradU_z;
      M.y.z = R.y.x*Rz_gradU_x + R.y.y*Rz_gradU_y + R.y.z*Rz_gradU_z;

      M.z.x = R.z.x*Rx_gradU_x + R.z.y*Rx_gradU_y + R.z.z*Rx_gradU_z;
      M.z.y = R.z.x*Ry_gradU_x + R.z.y*Ry_gradU_y + R.z.z*Ry_gradU_z;
      M.z.z = R.z.x*Rz_gradU_x + R.z.y*Rz_gradU_y + R.z.z*Rz_gradU_z;

      // Compute the off-diagonal elements of the Omega tensor in the eigenvector basis
      double omega_xy = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
      double omega_xz = (Lambda.z*M.x.z + Lambda.x*M.z.x)/(Lambda.z - Lambda.x);
      double omega_yz = (Lambda.z*M.y.z + Lambda.y*M.z.y)/(Lambda.z - Lambda.y);

      // Calculate intermediate terms for the transformation back to the original coordinate system
      double omega_xy_term_x = R.x.x*omega_xy - R.x.z*omega_yz;
      double omega_xy_term_y = R.y.x*omega_xy + R.y.z*omega_xz;
      double omega_xz_term_x = R.x.x*omega_xz + R.y.x*omega_yz;
      double omega_xy_term_xy = R.x.y*omega_xy - R.x.z*omega_yz;
      double omega_xy_term_yy = R.y.y*omega_xy + R.y.z*omega_xz;
      double omega_xz_term_xy = R.x.y*omega_xz + R.y.y*omega_yz;
      double omega_xy_term_xz = R.x.z*omega_xy - R.x.z*omega_yz;
      double omega_xz_term_yz = R.y.z*omega_xy + R.z.z*omega_xz;
      double omega_yz_term_z = R.x.z*omega_xz + R.y.z*omega_yz;

      // Compute Omega components by transforming back to the original coordinate system: Omega = R * omega * R^T
      Omega.x.x = R.x.y*omega_xy_term_x - R.x.x*omega_xy_term_xy + R.x.z*omega_xy_term_xz;
      Omega.x.y = R.y.y*omega_xy_term_x - R.y.x*omega_xy_term_xy + R.y.z*omega_xy_term_xz;
      Omega.x.z = R.z.y*omega_xy_term_x - R.z.x*omega_xy_term_xy + R.z.z*omega_xy_term_xz;

      Omega.y.x = R.x.y*omega_xy_term_y - R.x.x*omega_xy_term_yy + R.x.z*omega_xz_term_yz;
      Omega.y.y = R.y.y*omega_xy_term_y - R.y.x*omega_xy_term_yy + R.y.z*omega_xz_term_yz;
      Omega.y.z = R.z.y*omega_xy_term_y - R.z.x*omega_xy_term_yy + R.z.z*omega_xz_term_yz;

      Omega.z.x = R.x.y*omega_xz_term_x - R.x.x*omega_xz_term_xy + R.x.z*omega_yz_term_z;
      Omega.z.y = R.y.y*omega_xz_term_x - R.y.x*omega_xz_term_xy + R.y.z*omega_yz_term_z;
      Omega.z.z = R.z.y*omega_xz_term_x - R.z.x*omega_xz_term_xy + R.z.z*omega_yz_term_z;
      
      // Extract diagonal components of M (velocity gradient tensor in eigenvector basis)
      double M_diag_x = M.x.x, M_diag_y = M.y.y, M_diag_z = M.z.z;

      /*
      Compute B tensor: B = R * diag(M) * R^T
      - This transforms the diagonal velocity gradient tensor back to the original coordinate system
      - B is symmetric, so we only need to compute the upper triangle
      */

      // Compute diagonal elements of B
      B.x.x = M_diag_x*sq(R.x.x) + M_diag_y*sq(R.y.x) + M_diag_z*sq(R.z.x);
      B.y.y = M_diag_x*sq(R.x.y) + M_diag_y*sq(R.y.y) + M_diag_z*sq(R.z.y);
      B.z.z = M_diag_x*sq(R.x.z) + M_diag_y*sq(R.y.z) + M_diag_z*sq(R.z.z);

      // Compute off-diagonal elements of B (upper triangle)
      B.x.y = M_diag_x*R.x.x*R.x.y + M_diag_y*R.y.x*R.y.y + M_diag_z*R.z.x*R.z.y;
      B.x.z = M_diag_x*R.x.x*R.x.z + M_diag_y*R.y.x*R.y.z + M_diag_z*R.z.x*R.z.z;
      B.y.z = M_diag_x*R.x.y*R.x.z + M_diag_y*R.y.y*R.y.z + M_diag_z*R.z.y*R.z.z;

      // Fill in lower triangle using symmetry of B
      B.y.x = B.x.y;
      B.z.x = B.x.z;
      B.z.y = B.y.z;
    }

    // Update Psi components
    double old_Psi11 = Psi11[];
    double old_Psi22 = Psi22[];
    double old_Psi33 = Psi33[];
    double old_Psi12 = Psi12[];
    double old_Psi13 = Psi13[];
    double old_Psi23 = Psi23[];

    // Psi11
    Psi11[] += dt * (2.0 * B.x.x + 2 * Omega.x.y * old_Psi12 + 2 * Omega.x.z * old_Psi13);
    // Psi22
    Psi22[] += dt * (2.0 * B.y.y - 2 * Omega.x.y * old_Psi12 + 2 * Omega.y.z * old_Psi23);
    // Psi33
    Psi33[] += dt * (2.0 * B.z.z - 2 * Omega.x.z * old_Psi13 - 2 * Omega.y.z * old_Psi23);
    // Psi12
    Psi12[] += dt * (2.0 * B.x.y - Omega.x.y * old_Psi11 + Omega.x.y * old_Psi22 + Omega.x.z * old_Psi23 + Omega.x.z * old_Psi13);
    // Psi13
    Psi13[] += dt * (2.0 * B.x.z - Omega.x.z * old_Psi11 + Omega.x.z * old_Psi33 + Omega.x.y * old_Psi23 + Omega.x.y * old_Psi12);
    // Psi23
    Psi23[] += dt * (2.0 * B.y.z - Omega.y.z * old_Psi22 + Omega.y.z * old_Psi33 - Omega.x.y * old_Psi13 - Omega.x.z * old_Psi12);

  }

  // Advection of Psi, which is the log-conformation tensor
  advection ({Psi11, Psi12, Psi13, Psi22, Psi23, Psi33}, uf, dt);

  /**
  ### Convert back to A and T

  We now convert the log-conformation tensor Psi back to the conformation tensor A
  and compute the stress tensor T. This process involves diagonalization,
  exponentiation of eigenvalues, and application of the relaxation factor.
  */

  foreach() {
    pseudo_t3d A, R;
    pseudo_v3d Lambda;

    // Reconstruct the log-conformation tensor from its components
    A.x.x = Psi11[]; A.x.y = Psi12[]; A.x.z = Psi13[];
    A.y.x = Psi12[]; A.y.y = Psi22[]; A.y.z = Psi23[];
    A.z.x = Psi13[]; A.z.y = Psi23[]; A.z.z = Psi33[];

    // Diagonalize A to obtain eigenvalues and eigenvectors
    diagonalization_3D(&Lambda, &R, &A);
    
    // Exponentiate eigenvalues
    Lambda.x = exp(Lambda.x);
    Lambda.y = exp(Lambda.y);
    Lambda.z = exp(Lambda.z);

    // Reconstruct A using A = R * diag(Lambda) * R^T
    A11[] = Lambda.x * R.x.x * R.x.x + Lambda.y * R.x.y * R.x.y + Lambda.z * R.x.z * R.x.z;
    A12[] = Lambda.x * R.x.x * R.y.x + Lambda.y * R.x.y * R.y.y + Lambda.z * R.x.z * R.y.z;
    A13[] = Lambda.x * R.x.x * R.z.x + Lambda.y * R.x.y * R.z.y + Lambda.z * R.x.z * R.z.z;
    A22[] = Lambda.x * R.y.x * R.y.x + Lambda.y * R.y.y * R.y.y + Lambda.z * R.y.z * R.y.z;
    A23[] = Lambda.x * R.y.x * R.z.x + Lambda.y * R.y.y * R.z.y + Lambda.z * R.y.z * R.z.z;
    A33[] = Lambda.x * R.z.x * R.z.x + Lambda.y * R.z.y * R.z.y + Lambda.z * R.z.z * R.z.z;

    // Apply relaxation using the relaxation time lambda
    double intFactor = lambda[] != 0. ? exp(-dt/lambda[]) : 0.;
    
    A11[] = 1. + (A11[] - 1.)*intFactor;
    A22[] = 1. + (A22[] - 1.)*intFactor;
    A33[] = 1. + (A33[] - 1.)*intFactor;
    A12[] *= intFactor;
    A13[] *= intFactor;
    A23[] *= intFactor;

    // Compute the stress tensor T using the polymer modulus Gp
    T11[] = Gp[]*(A11[] - 1.);
    T22[] = Gp[]*(A22[] - 1.);
    T33[] = Gp[]*(A33[] - 1.);
    T12[] = Gp[]*A12[];
    T13[] = Gp[]*A13[];
    T23[] = Gp[]*A23[];
  }
}
#endif

/**
## Divergence of the viscoelastic stress tensor

The viscoelastic stress tensor $\mathbf{T}$ is defined at cell centers,
while the corresponding force (acceleration) is defined at cell faces. 
For each component of the momentum equation, we need to compute the 
divergence of the corresponding row of the stress tensor.

For example, for the x-component in 3D:

$$
(\nabla \cdot \mathbf{T})_x = \partial_x T_{xx} + \partial_y T_{xy} + \partial_z T_{xz}
$$

The normal stress gradient (e.g. $\partial_x T_{xx}$) is computed directly 
from cell-centered values. The shear stress gradients (e.g. $\partial_y T_{xy}$) 
are computed using vertex-averaged values to avoid checkerboard instabilities.
*/

event acceleration (i++)
{
  face vector av = a;

#if dimension == 2
  // 2D implementation
  foreach_face(x) {
    if (fm.x[] > 1e-20) {
      // y-gradient of T12 (shear stress)
      double shearX = (T12[0,1]*cm[0,1] + T12[-1,1]*cm[-1,1] - 
                       T12[0,-1]*cm[0,-1] - T12[-1,-1]*cm[-1,-1])/4.;
      // x-gradient of T11 (normal stress)
      double gradX_T11 = cm[]*T11[] - cm[-1]*T11[-1];
      
      av.x[] += (shearX + gradX_T11)*alpha.x[]/(sq(fm.x[])*Delta);
    }
  }
  
  foreach_face(y) {
    if (fm.y[] > 1e-20) {
      // x-gradient of T12 (shear stress)
      double shearY = (T12[1,0]*cm[1,0] + T12[1,-1]*cm[1,-1] - 
                       T12[-1,0]*cm[-1,0] - T12[-1,-1]*cm[-1,-1])/4.;
      // y-gradient of T22 (normal stress)
      double gradY_T22 = cm[]*T22[] - cm[0,-1]*T22[0,-1];
      
      av.y[] += (shearY + gradY_T22)*alpha.y[]/(sq(fm.y[])*Delta);
    }
  }

#elif dimension == 3
  // 3D implementation
  foreach_face(x) {
    if (fm.x[] > 1e-20) {
      // y-gradient of T12
      double shearY = (T12[0,1,0]*cm[0,1,0] + T12[-1,1,0]*cm[-1,1,0] - 
                       T12[0,-1,0]*cm[0,-1,0] - T12[-1,-1,0]*cm[-1,-1,0])/4.;
      // z-gradient of T13
      double shearZ = (T13[0,0,1]*cm[0,0,1] + T13[-1,0,1]*cm[-1,0,1] - 
                       T13[0,0,-1]*cm[0,0,-1] - T13[-1,0,-1]*cm[-1,0,-1])/4.;
      // x-gradient of T11
      double gradX_T11 = cm[]*T11[] - cm[-1,0,0]*T11[-1,0,0];
      
      av.x[] += (shearY + shearZ + gradX_T11)*alpha.x[]/(sq(fm.x[])*Delta);
    }
  }

  foreach_face(y) {
    if (fm.y[] > 1e-20) {
      // x-gradient of T12
      double shearX = (T12[1,0,0]*cm[1,0,0] + T12[1,-1,0]*cm[1,-1,0] - 
                       T12[-1,0,0]*cm[-1,0,0] - T12[-1,-1,0]*cm[-1,-1,0])/4.;
      // z-gradient of T23
      double shearZ = (T23[0,0,1]*cm[0,0,1] + T23[0,-1,1]*cm[0,-1,1] - 
                       T23[0,0,-1]*cm[0,0,-1] - T23[0,-1,-1]*cm[0,-1,-1])/4.;
      // y-gradient of T22
      double gradY_T22 = cm[]*T22[] - cm[0,-1,0]*T22[0,-1,0];
      
      av.y[] += (shearX + shearZ + gradY_T22)*alpha.y[]/(sq(fm.y[])*Delta);
    }
  }

  foreach_face(z) {
    if (fm.z[] > 1e-20) {
      // x-gradient of T13
      double shearX = (T13[1,0,0]*cm[1,0,0] + T13[1,0,-1]*cm[1,0,-1] - 
                       T13[-1,0,0]*cm[-1,0,0] - T13[-1,0,-1]*cm[-1,0,-1])/4.;
      // y-gradient of T23
      double shearY = (T23[0,1,0]*cm[0,1,0] + T23[0,1,-1]*cm[0,1,-1] - 
                       T23[0,-1,0]*cm[0,-1,0] - T23[0,-1,-1]*cm[0,-1,-1])/4.;
      // z-gradient of T33
      double gradZ_T33 = cm[]*T33[] - cm[0,0,-1]*T33[0,0,-1];
      
      av.z[] += (shearX + shearY + gradZ_T33)*alpha.z[]/(sq(fm.z[])*Delta);
    }
  }
#endif
}

