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
  S = 0
  T = 1
  V = 2
  ERI = 3
} integrals;

void 1eint(sys*, integrals type);
