/** Title: Encapsulation
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Three Phase Taylor Culick
# Last Updated: Dec 28 2020
*/

// 1 is Si Pool, 2 is Water Film and 3 is air
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase.h"
#include "tension.h"
#include "distance.h"
// #include "adapt_wavelet_limited.h"

#define MINlevel 3                                              // maximum level

#define tsnap (5e-1)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in velocity

#define z0 1e-1
#define h0 1e0
#define Ohf 0.10 // 1 micron water film is 0.16

#define Mus 0.01
#define Ohs (Mus*Ohf)

#define Rho21 1.000
#define Rho31 0.001

// Surface tesnion
#define SIGMA13 (0.50)
#define SIGMA12 (1.00)

// density
#define Rho1 (1.0)                                         // density of phase 1
#define Rho2 Rho21                                         // density of phase 2
#define Rho3 Rho31                                         // density of phase 3

// viscosity
#define Mu1 (Ohp)                            // Dynamic viscosity of phase 1
#define Mu2 (Ohf)                            // Dynamic viscosity of phase 2
#define Mu3 (Ohs)                       // Dynamic viscosity of phase 3

// boundary conditions

double Mur, Ohp, tmax, Ldomain;
int MAXlevel;

int main(int argc, char const *argv[]) {
  Mur = atof(argv[1]);
  Ohp = (Mur*Ohf);
  tmax = atof(argv[2]);
  MAXlevel = atoi(argv[3]);
  Ldomain = atof(argv[4]);

  L0=Ldomain;
  X0=-L0/2.; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = Rho1; mu1 = Mu1;
  rho2 = Rho2; mu2 = Mu2;
  rho3 = Rho3; mu3 = Mu3;

  f1.sigma = SIGMA13;
  f2.sigma = SIGMA12;

  fprintf(ferr, "Level %d tmax %g. Ohf %3.2e, Ohp %3.2e\n", MAXlevel, tmax, Ohf, Ohp);
  run();
}

// int refRegion(double x, double y, double z){
//   return (x < -1 ? 2:
//           x > 1.5 ? 2:
//           y > 2.0 ? 2:
//           sq(x)+sq(y)<sq(z0) ? MAXlevel+1:
//           MAXlevel
//         );
// }

event init(t = 0){
  if(!restore (file = "dump")){
    char filename1[60], filename2[60];
    /**
    Initialization for f1 & f2
    */
    sprintf(filename1, "f1_z0%3.2e_h0%3.2e.dat", z0, h0);
    sprintf(filename2, "f2_z0%3.2e_h0%3.2e.dat", z0, h0);

    FILE * fp1 = fopen(filename1,"rb");
    if (fp1 == NULL){
      fprintf(ferr, "There is no file named %s\n", filename1);
      return 1;
    }
    FILE * fp2 = fopen(filename2,"rb");
    if (fp2 == NULL){
      fprintf(ferr, "There is no file named %s\n", filename2);
      return 1;
    }
    coord* InitialShape1;
    coord* InitialShape2;
    InitialShape1 = input_xy(fp1);
    fclose (fp1);
    InitialShape2 = input_xy(fp2);
    fclose (fp2);
    scalar d1[], d2[];
    distance (d1, InitialShape1);
    distance (d2, InitialShape2);
    while (adapt_wavelet ((scalar *){f1, f2, d1, d2}, (double[]){1e-8, 1e-8, 1e-8, 1e-8}, MAXlevel).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi1[], phi2[];
    foreach_vertex(){
      phi1[] = -(d1[] + d1[-1] + d1[0,-1] + d1[-1,-1])/4.;
      phi2[] = -(d2[] + d2[-1] + d2[0,-1] + d2[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fractions in the domain. */
    fractions (phi1, f1);
    fractions (phi2, f2);
  }
}

scalar KAPPA1[], KAPPA2[], omega[];
event adapt(i++) {
  vorticity (u, omega);
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  foreach(){
    omega[] *= f1[]*(1-f2[]);
  }
  boundary ((scalar *){omega, KAPPA1, KAPPA2});
  adapt_wavelet ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2, omega},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr, OmegaErr},
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
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f1[],f2[]);
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
