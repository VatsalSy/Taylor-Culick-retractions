/** Title: Encapsulation
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

// 1 is Si Pool, 2 is Water Drop and 3 is air
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "tension.h"
// #include "adapt_wavelet_limited.h"
//
// #define locRef 1e0

#define MAXlevel 10                                              // maximum level
#define MINlevel 3                                              // maximum level
#define tmax 100.0                                                 // maximum time
#define tsnap (0.5)
// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in velocity

#define h0 1e0
#define Ohf 0.1  // 1 micros thickness water is order 0.1 Oh \sim h^(-1/2)
#define Mus (1e-4)
#define Rho31 1.0

// density
#define Rho1 (1.0)                                         // density of phase 1
#define Rho3 Rho31                                         // density of phase 3

// viscosity
#define Mu1 (Ohf)                            // Dynamic viscosity of phase 1
#define Mu3 (Mus*Ohf)                       // Dynamic viscosity of phase 3

// domain
#define Ldomain 100                            // Dimension of the domain

// boundary conditions
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

int main(){
  L0=Ldomain;
  X0=0.0; Y0=0.;
  init_grid (1 << (6));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = Rho1; mu1 = Mu1;
  rho2 = Rho3; mu2 = Mu3;

  f.sigma = 1.0;

  fprintf(ferr, "Level %d tmax %g. Ohf %3.2e, h0 %3.2e\n", MAXlevel, tmax, Ohf, h0);
  run();
}



event init(t = 0){
  if(!restore (file = "dump")){
    /**
    We can now initialize the volume fractions in the domain. */
    refine(x<1.02 && y < h0+1.02 && level<MAXlevel);
    fraction(f, y < h0+1 ? 1-(sq(x)+sq(y-1-h0)) : 1-x);
    f.prolongation = refine_bilinear;
    boundary((scalar *){f});

  }
}

// int refRegion(double x, double y, double z){
//   return (x > 10 ? 0:
//     x > 5 ? MAXlevel-4:
//     y > 20 ? MAXlevel-4:
//     y < 2*locRef ? MAXlevel+4:
//     y < 4*locRef ? MAXlevel+2:
//     MAXlevel
//     );
// }
scalar KAPPA[], omega[];
event adapt(i++) {
  vorticity (u, omega);
  curvature(f, KAPPA);
  foreach(){
    omega[] *= f[];
  }
  boundary ((scalar *){omega, KAPPA});
  // adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA, omega},
  //   (double[]){fErr, VelErr, VelErr, KErr, OmegaErr},
  //   refRegion, MINlevel);
  adapt_wavelet((scalar *){f, u.x, u.y, KAPPA, omega},
    (double[]){fErr, VelErr, VelErr, KErr, OmegaErr},
    MAXlevel, MINlevel);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f[]);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
}
