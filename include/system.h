#pragma once
#include <stdio.h>
#include "basis.h"

typedef struct system {
  int natoms;
  int nbfs;
  bfn* basisfunctions;
  double* S;
  double* T;
  double* V;
  double* ERI;
} sys;

void print_system_info(sys*);

void system_from_file(FILE*, char*);

sys* system_create(int, int);

void system_destroy(sys*);
