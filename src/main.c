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
  /*
  bfn H1 = {
    3,
    {3.42525091, 0.62391373, 0.16885540},
    {0.15432897, 0.53532814, 0.44463454},
    {0, 0, 0.7},
    {0, 0, 0},
    {0, 0, 0},
    1
  };
  
  bfn H2 = {
    3,
    {3.42525091, 0.62391373, 0.16885540},
    {0.15432897, 0.53532814, 0.44463454},
    {0, 0, -0.7},
    {0, 0, 0},
    {0, 0, 0},
    1
  };
  
  sys hyd = {
    2,
    2,
    {H1, H2}     
  };
  normalise(&H1);
  normalise(&H2);
  for(int i = 0; i < hyd.nbfs; i++)
    normalise(&(hyd.basisfunctions[i]));
  calculate_integral(&hyd);
  */
  return 0;
}
