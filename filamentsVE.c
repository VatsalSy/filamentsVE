/**
 * @file filamentsVE.c
 * @brief This code will use the initial condition from filament_initialCondition.c and simulate the filament retraction.
 * @author Vatsal Sanjay
 * @version 2.0
 * @date Oct 26, 2024
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
#include "two-phaseVE.h"

#include "log-conform-viscoelastic-scalar-2D.h"
#define logFile "logAxi-VE.dat"

#include "navier-stokes/conserving.h"
#include "tension.h"

#define tsnap (1e-2)

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-6)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-2)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define trAErr (1e-3)                                // error tolerance in trace of conformation tensor

#define R2(x,y,z) (sqrt(sq(x) + sq(y)))
#define U0 4.0 // this is an adhoc initial condition to make the filament form!

// boundary conditions
u.n[top] = neumann(0.0);
p[top] = dirichlet(0.0);

int MAXlevel;
// Oh -> Solvent Ohnesorge number
// Oha -> air Ohnesorge number
// De -> Deborah number
// Ec -> Elasto-capillary number
// for now there is no viscoelasticity

double Oh, Oha, De, Ec, tmax;
char nameOut[80], dumpFile[80];

int main(int argc, char const *argv[]) {

  L0 = 16.;
  X0 = -L0/2.;
  
  // Values taken from the terminal
  MAXlevel = 12;
  tmax = 1.75;
  Oh = 1.25e-2;
  Oha = 1e-2 * Oh;
  De = 0.10; // 1e-1;
  Ec = 0.25; // 1e-2;

  init_grid (1 << 6);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart_init_stretched_filament");


  rho1 = 1., rho2 = 1e-3;
  mu1 = Oh, mu2 = Oha;
  lambda1 = De, lambda2 = 0.;
  G1 = Ec, G2 = 0.;
  f.sigma = 1.0;

  run();

}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    fprintf(ferr, "No restart file found. Exiting!\n");
    return 1;
  }
  foreach() {
    foreach_dimension()
      u.x[] = 0.0;
  }
}

/**
## Adaptive Mesh Refinement
*/
scalar KAPPA[], trA[];

event adapt(i++){
  curvature(f, KAPPA);
  foreach() {
    trA[] = (A11[] + A22[] + conform_qq[]);
  }
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, trA},
      (double[]){fErr, VelErr, VelErr, KErr, trAErr},
      MAXlevel, 4);
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Ending Simulation
*/
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, Oh %2.1e\n", MAXlevel, Oh);
}

/**
## Log writing
*/
event logWriting (i++) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])+ sq(u.z[])))*sq(Delta);
  }

  static FILE * fp;
  if (pid() == 0) {
    const char* mode = (i == 0) ? "w" : "a";
    fp = fopen(logFile, mode);
    if (fp == NULL) {
      fprintf(ferr, "Error opening log file\n");
      return 1;
    }

    if (i == 0) {
      fprintf(ferr, "Level %d, Oh %2.1e, Oha %2.1e, De %2.1e, Ec %2.1e\n", MAXlevel, Oh, Oha, De, Ec);
      fprintf(ferr, "i dt t ke\n");
      fprintf(fp, "Level %d, Oh %2.1e, Oha %2.1e, De %2.1e, Ec %2.1e\n", MAXlevel, Oh, Oha, De, Ec);
      fprintf(fp, "i dt t ke\n");
    }

    fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
    fprintf(ferr, "%d %g %g %g\n", i, dt, t, ke);

    fflush(fp);
    fclose(fp);
  }

  assert(ke > -1e-10);

  if (i > 1e4 && pid() == 0) {
    if (ke > 1e2*sq(U0) || ke < 1e-8) {
      const char* message = (ke > 1e2) ? 
        "The kinetic energy blew up. Stopping simulation\n" : 
        "kinetic energy too small now! Stopping!\n";
      
      fprintf(ferr, "%s", message);
      
      fp = fopen("log", "a");
      fprintf(fp, "%s", message);
      fflush(fp);
      fclose(fp);
      
      dump(file=dumpFile);
      return 1;
    }
  }

}