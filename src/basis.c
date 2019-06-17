#include "basis.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

bfn* bfn_create(int nprimitives, double* exps, double* coefs, double origin[3], int shell[3], double* norm, double atomno){
  
  bfn *ret_bfn = malloc(sizeof (bfn));
  if (ret_bfn == NULL)
    return NULL;

  ret_bfn->nprimitives = nprimitives;
  ret_bfn->exps = malloc(sizeof (double) * nprimitives);
  if (ret_bfn->exps == NULL){
    free(ret_bfn);
    return NULL;
  }
  ret_bfn->coefs = malloc(sizeof (double) * nprimitives);
  if (ret_bfn->coefs == NULL){
    free(ret_bfn);
    return NULL;
  }
  ret_bfn->norm = malloc(sizeof (double) * nprimitives);
  if (ret_bfn->norm == NULL){
    free(ret_bfn);
    return NULL;
  }
  for(int i = 0; i < nprimitives; i++){
    ret_bfn->exps[i] = exps[i];
    ret_bfn->coefs[i] = coefs[i];
    ret_bfn->norm[i] = norm[i];
  }
  for(int i = 0; i < 3; i++){
    ret_bfn->shell[i] = shell[i];
    ret_bfn->origin[i] = origin[i];
  }
  ret_bfn->atomno = atomno;
  return ret_bfn;
}

void bfn_destroy(bfn* in_bfn){
  if (in_bfn != NULL) {
    free(in_bfn->coefs);
    free(in_bfn->exps);
    free(in_bfn->norm);
    free(in_bfn);
  }
}
