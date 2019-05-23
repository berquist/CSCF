#include "basis.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

void normalise(bfn* basis){
  double l, m, n;
  l = basis->shell[0];
  m = basis->shell[1];
  n = basis->shell[2];
  for(int i = 0; i < 3; i++){
    double num = pow(2.0, 2.0*(l+m+n) + 3.0/2.0) * pow(basis->exps[i], (l+m+n) + 3.0/2.0);
    double denom = double_factorial(2*l-1) * double_factorial(2*m-1) * double_factorial(2*n-1) * pow(M_PI, 3.0/2.0);
    basis->norm[i] = sqrt(num/denom);
  }
}
