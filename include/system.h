#pragma once
#include "basis.h"

typedef struct system {
  int natoms;
  int nbfs;
  bfn basisfunctions[2];
  double S[4];
  double T[4];
  double V[4];
} sys;

void init_system(sys*, FILE*, char*);
