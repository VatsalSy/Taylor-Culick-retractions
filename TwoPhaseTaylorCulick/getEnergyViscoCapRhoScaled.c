/* Title: Getting Energy
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Feb 10 2021
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

trace
double interface_energy (scalar c){
  double se = 0.;
  foreach (reduction(+:se)){
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord p, n = interface_normal (point, c);
      double alpha = plane_alpha (c[], n);
      double len = line_length_center(n, alpha, &p);
      se += 2.*pi*( y + p.y*Delta )*(len*Delta); // 2*pi*\int_l (r_c)dl
    }
  }
  return se;
}

scalar f[];
double ke1, ke, se, rho1, rho2, Rhor, Ohf, mu1, mu2, Ohs, eps1, eps;

char filename[80], nameEnergy[80];



int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf(nameEnergy, "%s", arguments[2]);
  Rhor = atof(arguments[3]);
  Ohf = atof(arguments[4]);
  Ohs = atof(arguments[5]);
  // fprintf(ferr, "Rhor %g, Ohf %3.2e, Ohs %3.2e\n", Rhor, Ohf, Ohs);
  // return 1;
  // boundary conditions

  FILE *fp;
  fp = fopen (nameEnergy, "a");
  restore (file = filename);

  rho1 = 1.0/sq(Ohs); mu1 = Ohf/Ohs;
  rho2 = Rhor/sq(Ohs); mu2 = 1.00;

  f.prolongation = refine_bilinear;

  boundary((scalar *){f, u.x, u.y});

  /*
  Do calculations start
  */
  ke1 = 0., ke = 0., se = 0., eps1 = 0., eps = 0.;

  foreach (){
    double rho = clamp(f[], 0., 1.)*(rho1 - rho2) + rho2;
    ke1 += (2*pi*y)*(0.5*clamp(f[], 0., 1.)*rho1*(sq(u.x[]) + sq(u.y[])))*sq(Delta); // 2*pi*\int_A(0.5*rho*v^2)*r*dr*dz
    ke += (2*pi*y)*(0.5*rho*(sq(u.x[]) + sq(u.y[])))*sq(Delta); // 2*pi*\int_A(0.5*rho*v^2)*r*dr*dz

    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));

    double mu = clamp(f[], 0., 1.)*(mu1 - mu2) + mu2;
    eps1 += (2*pi*y)*( 2*mu1*clamp(f[], 0., 1.)*D2 )*sq(Delta);
    eps += (2*pi*y)*( 2*mu*D2 )*sq(Delta);
  }

  boundary((scalar *){f, u.x, u.y});

  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});
  se = interface_energy (f);

  /*
  Do calculations end
  */
  if (t == 0){
    fprintf(ferr, "Rhor %g, Ohf %3.2e, Ohs %3.2e\n", Rhor, Ohf, Ohs);
    fprintf(ferr, "t ke ke1 se eps eps1\n");
    fprintf(fp, "t ke ke1 se eps eps1\n");    
  }

  fprintf(ferr, "%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n", t, ke, ke1, se, eps, eps1);
  fprintf(fp, "%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n", t, ke, ke1, se, eps, eps1);
  fflush (fp);
  fclose (fp);
}
