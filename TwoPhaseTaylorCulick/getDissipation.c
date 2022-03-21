/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
First Draft
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80], outName[80];
int Np;
double Z0, R0, Rmax, Ohf, Ohp, Ohs, l0, eps1, eps2;
scalar f[];

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  sprintf (outName, "%s", arguments[2]);
  Z0 = atof(arguments[3]); R0 = atof(arguments[4]);
  Rmax = atof(arguments[5]); Np = atoi(arguments[6]);
  Ohf = atof(arguments[7]); Ohs = atof(arguments[8]);

  // boundary conditions

  /*
  Actual run and codes!
  */
  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});


  double Deltap = (Rmax-R0)/(Np-1);

  l0 = 0;
  FILE *fp2;
  fp2 = fopen (outName, "a");

  for (int i = 0; i < Np; i++){
    eps1 = 0., eps2 = 0.;
    l0 += Deltap;
    foreach(){
      if (sq(x-Z0)+sq(y-R0) < sq(l0)){
        double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
        double D22 = (u.y[]/max(y, 1e-30));
        double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
        double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
        double D2 = sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13);

        eps1 += 2*y*(clamp(1.-f[], 0., 1.) * Ohs)*D2*sq(Delta);
        eps2 += 2*y*(clamp(f[], 0., 1.) * Ohf)*D2*sq(Delta);
      }
    }
    fprintf(fp2, "%f %f %5.3e %4.3e\n", t, l0, eps1, eps2);
    // fprintf(ferr, "t = %f, l0 = %f, eps1 = %f, eps2 = %f\n", t, l0, eps1, eps2);
  }
  fclose(fp2);
  boundary((scalar *){f, u.x, u.y});
}
