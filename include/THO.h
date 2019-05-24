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

double f(int, int, int, double, double);

double overlap_1d(int, int, double, double, double);

double overlap(int*, double*, double, int*, double*, double);

double Sab(bfn, bfn);

double kinetic(int*, double*, double, int*, double*, double);

double Tab(bfn, bfn);

double A(int, int, int, int, int, double, double, double, double);

void G_array(int, int, double, double, double, double, double*);

double nuclear(int*, double*, double, int*, double*, double, double*);

double Vab(bfn, bfn, double*);
