/** Title: Encapsulation
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Updated: Jan 07, 2021
*/

// 1 is Si Pool, 2 is Water Drop and 3 is air
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "tension.h"
// #include "adapt_wavelet_limited.h"
//

#define MAXlevel 10                                              // maximum level

#define MINlevel 3                                              // maximum level
#define tmax 25.0                                                 // maximum time
#define tsnap (0.5)
// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-3)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in velocity

#define hole0 (0.5e0)
#define h0 (1e0)

// domain
#define Ldomain 50                            // Dimension of the domain

// boundary conditions

double Ohf, Ohs, Rho31;

int main(int argc, char const *argv[]) {
  Ohf = atof(argv[1]);
  Ohs = atof(argv[2]);
  Rho31 = atof(argv[3]);

  L0=Ldomain;
  X0=0.0; Y0=0.;
  init_grid (1 << (6));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.000; mu1 = Ohf;
  rho2 = Rho31; mu2 = Ohs;

  f.sigma = 0.50;

  fprintf(ferr, "Level %d tmax %g. Ohf %3.2e, h0 %3.2e, Ohs %g\n", MAXlevel, tmax, Ohf, h0, Ohs);
  run();
}



event init(t = 0){
  if(!restore (file = "dump")){
    /**
    We can now initialize the volume fractions in the domain. */
    refine(x<(h0/2.0)+0.025 && y < hole0+(h0/2.0)+0.025 && level<MAXlevel);
    fraction(f, y < hole0+(h0/2.0) ? sq(h0/2.0)-(sq(x)+sq(y-(h0/2.0)-hole0)) : (h0/2.0)-x);
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
