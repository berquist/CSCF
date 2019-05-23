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

void calculate_integral(sys*);

double f(double, double, double, double, double);

double overlap_1d(int, int, double, double, double);

double overlap(int*, double*, double, int*, double*, double);

double Sab(bfn, bfn);

