#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* Very simple HF implementation for H2. 
 * Parser for files etc. will be added later.
 * Right now, parameters and geometry is hard-coded.
 * This is so I can focus on the maths. 
 * Once I've got it working for one model, I'll start working on testing as well.
 */

#include "basis.h"
#include "system.h"
#include "THO.h"
#include "utils.h"

int main(int argc, char *argv[])
{
  double e[3] = {3.42525091, 0.62391373, 0.16885540};
  double c[3] = {0.15432897, 0.53532814, 0.44463454};
  int s[3] = { 0, 0, 0 };
  double c1[3] = {0.0, 0.0, 0.7};
  double c2[3] = {0.0, 0.0, -0.7};
  sys* hyd = system_create(2, 2);
  bfn* H1 = bfn_create(3, e, c, c1, s, 1);
  bfn* H2 = bfn_create(3, e, c, c2, s, 1);
  hyd->basisfunctions[0] = *H1;
  hyd->basisfunctions[1] = *H2;
  calculate_integral(hyd);
  bfn_destroy(H1);
  bfn_destroy(H2);
  system_destroy(hyd);
  return 0;
}
