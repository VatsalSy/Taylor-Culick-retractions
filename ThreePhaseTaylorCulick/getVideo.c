/* Title: Saving images with bview
# Authors: Vatsal & Youssef
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Dec 24, 2020
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"

scalar f1[], f2[], vel[];
char filename[80], Imagename[80];
double Vmax;

int main(int a, char const *arguments[]) {
  // boundary conditions
  sprintf (filename, "%s", arguments[1]);
  sprintf (Imagename, "%s", arguments[2]);
  Vmax = atof(arguments[3]);
  restore (file = filename);
  foreach(){
    vel[] = sqrt(sq(u.x[]) + sq(u.y[]));
  }
  boundary((scalar *){f1, f2, u.x, u.y, vel});

  view (fov = 20.0, quat = {0,0,-0.707107,0.707107}, tx = 0.0, ty = 0.0, bg = {1,1,1}, width = 1920, height = 1016, samples = 4);
  box(notics=true);
  draw_vof ("f1", lc = {0.0, 0.0, 0.0}, lw=6);
  draw_vof ("f2", lc = {0.85, 0.3725, 0.01}, lw=6);
  squares ("vel", min=0, max=Vmax,  linear=true);
  mirror ({0,1}) {
    draw_vof ("f1", lc = {0.0, 0.0, 0.0}, lw=6);
    draw_vof ("f2", lc = {0.85, 0.3725, 0.01}, lw=6);
    squares ("f2", min=-1, max=1, map = cool_warm);
    cells(lc = {0.5, 0.5, 0.5});
  }
  save (Imagename);
}
