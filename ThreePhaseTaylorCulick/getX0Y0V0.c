/* Tracking Three Phase Contact line (apparent)
# Last Update: November 12, 2020
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f1[], f2[];
char filename[80], nameTrack[80];
double rOld, DistCutoff;
double wt, xTP, yTP, vTP;

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  sprintf(nameTrack, "%s", arguments[2]);
  rOld = atof(arguments[3]);
  DistCutoff = atof(arguments[4]);
  restore (file = filename);
  f1.prolongation = fraction_refine;
  f2.prolongation = fraction_refine;
  boundary((scalar *){f1, f2});

  // double x0 = -4., y0 = 0.;
  face vector s[];
  s.x.i = -1;
  wt = 0., xTP = 0., yTP = 0., vTP = 0.;
  foreach(){
    if (f2[] > 1e-6 && f2[] < 1. - 1e-6) {
      coord n2 = facet_normal (point, f2, s);
      double alpha2 = plane_alpha (f2[], n2);
      coord segment2[2];
      if (facets (n2, alpha2, segment2) == 2){
        if (f1[] > 1e-6 && f1[] < 1. - 1e-6) {
          coord n1 = facet_normal (point, f1, s);
          double alpha1 = plane_alpha (f1[], n1);
          coord segment1[2];
          if (facets (n1, alpha1, segment1) == 2){
            double x1 = x + (segment1[0].x+segment1[1].x)*Delta/2.;
            double y1 = y + (segment1[0].y+segment1[1].y)*Delta/2.;
            double x2 = x + (segment2[0].x+segment2[1].x)*Delta/2.;
            double y2 = y + (segment2[0].y+segment2[1].y)*Delta/2.;
            double dist = sqrt(sq(x1-x2)+sq(y1-y2));
            double xc = 0.5*(x1+x2);
            double yc = 0.5*(y1+y2);
            if (dist > DistCutoff && yc > rOld){
              wt += 1.0/dist;
              xTP += xc/dist;
              yTP += yc/dist;
              vTP += u.y[]/dist;
            }
          }
        }
      }
    }
    //
    // if(interfacial(point,f1) & interfacial(point,f2) & x > x0){
    //   x0 = x;
    //   y0 = y;
    // }
  }
  if (wt > 0.){
    xTP /= wt;
    yTP /= wt;
    vTP /= wt;
  } else {
    yTP = rOld;
  }
  fprintf(ferr, "%f %f %f %4.3e\n", t, xTP, yTP, vTP);
  FILE *fp2;
  fp2 = fopen (nameTrack, "a");
  fprintf(fp2, "%f %f %f %4.3e\n", t, xTP, yTP, vTP);
  fclose(fp2);
  // fprintf(ferr, "%f %f\n", x0, y0);
}
