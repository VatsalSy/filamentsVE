// making videos

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"

scalar f[];
char filename[80], Imagename[80];
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  sprintf (Imagename, "%s",arguments[2]);
  restore (file = filename);
  
view (quat = {-0.707, 0.000, 0.000, 0.707},
      fov = 30, near = 0.01, far = 1000,
      tx = 0.008, ty = -0.011, tz = -1.133,
      width = 1756, height = 1141);
  
      
  // box ();
  //cells ();
  draw_vof (c = "f");
  begin_mirror(n={0,0,1},alpha=0.);
  draw_vof (c = "f");
  end_mirror();
  save (Imagename);
}
