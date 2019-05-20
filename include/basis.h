#pragma once

struct bfn {
  int nprimitives;
  double* exps;
  double* coefs;
  double* origin;
  int* shell;
  double norm;
  double massno;
};

void bfn_set(bfn*, int, double*, double*, double*, int*, double, double);
