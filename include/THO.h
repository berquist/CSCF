#pragma once

#include "basis.h"
#include "system.h"

/**
 *   @file THO.h
 *   @brief A Documented file.
 *
 *  Detailed description
 *
 */

typedef enum integrals {
  S = 0,
  T = 1,
  V = 2,
  ERI = 3
} integrals;

void calculate_integral(sys*, integrals type);

double f(double, double, double, double, double);

double overlap_1d(int, int, double, double, double);

double overlap(int*, double*, double, int*, double*, double);

double Sab(bfn*, bfn*);

