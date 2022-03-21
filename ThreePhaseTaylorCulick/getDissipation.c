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
int N;
double x0, y0, Rmax, Ohf, Ohp, Ohs, l0, eps;
scalar f1[], f2[], D2c[], vel[];
scalar * list = NULL;

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  sprintf (outName, "%s", arguments[2]);
  x0 = atof(arguments[3]); y0 = atof(arguments[4]);
  Rmax = atof(arguments[5]); N = atoi(arguments[5]);
  Ohf = atof(arguments[6]); Ohp = atof(arguments[7]);
  Ohs = atof(arguments[8]);

  // boundary conditions

  /*
  Actual run and codes!
  */
  restore (file = filename);
  f1.prolongation = fraction_refine;
  f2.prolongation = fraction_refine;
  boundary((scalar *){f1, f2, u.x, u.y});

  double Delta = (Rmax-y0)/(N-1);
  l0 = 0, eps = 0;
  for (int i; i < N; i++){
    l0 += Delta;
    foreach(){
      if (sq(x-x0)+sq(y-y0) < sq(l0)){
        double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
        double D22 = (u.y[]/y);
        double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
        double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
        double D2 = sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13);

        eps += 2*y*(clamp(f1[]*(1-f2[]), 0., 1.) * Ohp + clamp(f1[]*f2[], 0., 1.) * Ohf + clamp((1-f1[]), 0., 1.) * Ohs)*D2;
      }
    }
    fprintf(ferr, "t = %f, l0 = %f, eps = %f\n",t, l0, eps);
  }
  boundary((scalar *){f1, f2, u.x, u.y});
}
