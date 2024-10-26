/** Title: log-conform-viscoelastic-2DNaxi.h
# Version: 1.0
# Main feature 1: A exists in across the domain and relaxes according to \lambda. The stress only acts according to G.
# Main feature 2: This is the 2D+axi **scalar** implementation of https://github.com/VatsalSy/BurstingBubble_VE_coated/blob/main/log-conform-viscoelastic.h.

# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Updated: Oct 17, 2024

# change log: Oct 18, 2024 (v1.0)
- 2D+axi implementation
- scalar implementation
*/

/** The code is same as http://basilisk.fr/src/log-conform.h but 
- written with G-\lambda formulation. 
- It also fixes the bug where [\sigma_p] = 0 & [\sigma_s] = \gamma\kappa instead of [\sigma_s+\sigma_p] = \gamma\kappa.
*/ 

#include "bcg.h"

scalar A11[], A12[], A22[]; // conformation tensor
scalar T11[], T12[], T22[]; // stress tensor
#if AXI
scalar conform_qq[], tau_qq[];
#endif

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;

  foreach() {
    A11[] = 1.; A22[] = 1.;
    A12[] = 0.;
    T11[] = 0.; T22[] = 0.; 
    T12[] = 0.;
#if AXI
    tau_qq[] = 0;
    conform_qq[] = 1.;
#endif
  }

  for (scalar s in {T11, T12, T22}) {
    if (s.boundary[left] != periodic_bc) {
        s[left] = neumann(0);
	      s[right] = neumann(0);
      }
  }

  for (scalar s in {A11, A12, A22}) {
    if (s.boundary[left] != periodic_bc) {
        s[left] = neumann(0);
	      s[right] = neumann(0);
    }
  }

#if AXI
  scalar s1 = T12;
  s1[bottom] = dirichlet (0.);  
#endif 

#if AXI
  scalar s2 = A12;
  s2[bottom] = dirichlet (0.);  
#endif
}

/**
## Useful functions in 2D

The first step is to implement a routine to calculate the eigenvalues
and eigenvectors of the conformation tensor $\mathbf{A}$.

These structs ressemble Basilisk vectors and tensors but are just
arrays not related to the grid. */

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

event tracer_advection(i++)
{
  scalar Psi11 = A11;
  scalar Psi12 = A12;
  scalar Psi22 = A22;
#if AXI
  scalar Psiqq = conform_qq;
#endif

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

    A.x.x = A11[]; A.y.y = A22[];
    A.x.y = A12[];

#if AXI
    double Aqq = conform_qq[]; 
    Psiqq[] = log (Aqq); 
#endif

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
	B.x.y = (u.y[1,0] - u.y[-1,0] +
		 u.x[0,1] - u.x[0,-1])/(4.*Delta); 
	foreach_dimension() 
	  B.x.x = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
      }
      else {
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

    /**
    In the axisymmetric case, the governing equation for $\Psi_{\theta
    \theta}$ only involves that component, 
    $$ 
    \Psi_{\theta \theta}|_t - 2 L_{\theta \theta} = 
    \frac{\mathbf{f}_r(e^{-\Psi_{\theta \theta}})}{\lambda} 
    $$
    with $L_{\theta \theta} = u_y/y$. Therefore step (a) for
    $\Psi_{\theta \theta}$ is */

#if AXI
    Psiqq[] += dt*2.*u.y[]/max(y, 1e-20);
#endif
  }

  /**
  ### Advection of $\Psi$
  
  We proceed with step (b), the advection of the log of the
  conformation tensor $\Psi$. */

#if AXI
  advection ({Psi11, Psi12, Psi22, Psiqq}, uf, dt);
#else
  advection ({Psi11, Psi12, Psi22}, uf, dt);
#endif

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
#if AXI
      double Aqq = exp(Psiqq[]);
#endif

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
     
#if AXI
      Aqq = (1. - intFactor) + intFactor*exp(Psiqq[]);
#endif

    A.x.y *= intFactor;
    foreach_dimension()
      A.x.x = (1. - intFactor) + A.x.x*intFactor;

    /**
      Then the Conformation tensor $\mathcal{A}_p^{n+1}$ is restored from
      $\mathbf{A}^{n+1}$.  */
    
    A12[] = A.x.y;
    T12[] = Gp[]*A.x.y;
#if AXI
      conform_qq[] = Aqq;
      tau_qq[] = Gp[]*(Aqq - 1.);
#endif

    A11[] = A.x.x;
    T11[] = Gp[]*(A.x.x - 1.);
    A22[] = A.y.y;
    T22[] = Gp[]*(A.y.y - 1.);
  }
}

/**
### Divergence of the viscoelastic stress tensor

The viscoelastic stress tensor $\mathbf{\tau}_p$ is defined at cell centers
while the corresponding force (acceleration) will be defined at cell
faces. Two terms contribute to each component of the momentum
equation. For example the $x$-component in Cartesian coordinates has
the following terms: $\partial_x \mathbf{\tau}_{p_{xx}} + \partial_y
\mathbf{\tau}_{p_{xy}}$. The first term is easy to compute since it can be
calculated directly from center values of cells sharing the face. The
other one is harder. It will be computed from vertex values. The
vertex values are obtained by averaging centered values.  Note that as
a result of the vertex averaging cells `[]` and `[-1,0]` are not
involved in the computation of shear. */

event acceleration (i++)
{
  face vector av = a;

  foreach_face(x){
    if (fm.x[] > 1e-20) {
      
      double shearX = (T12[0,1]*cm[0,1] + T12[-1,1]*cm[-1,1] - 
      T12[0,-1]*cm[0,-1] - T12[-1,-1]*cm[-1,-1])/4.;
      
      av.x[] += (shearX + cm[]*T11[] - cm[-1]*T11[-1])*
      alpha.x[]/(sq(fm.x[])*Delta);
    
    }
  }

  foreach_face(y){
    if (fm.y[] > 1e-20) {

      double shearY = (T12[1,0]*cm[1,0] + T12[1,-1]*cm[1,-1] - 
      T12[-1,0]*cm[-1,0] - T12[-1,-1]*cm[-1,-1])/4.;
      
      av.y[] += (shearY + cm[]*T22[] - cm[0,-1]*T22[0,-1])*
      alpha.y[]/(sq(fm.y[])*Delta);

    }
  }

#if AXI
  foreach_face(y)
    if (y > 1e-20)
      av.y[] -= (tau_qq[] + tau_qq[0,-1])*alpha.y[]/sq(y)/2.;
#endif
}
