#pragma once

/**
 *   @file basis.h
 *
 *   @author Abdullah Ahmad
 *
 *   @brief Header file describing the basis function data structure.
 *
 *  Basis function header defining the structure for the bfn struct as well as any member functions.
 *  Currently set only for H_{2} molecule.
 *
 */

typedef struct bfn {
  int nprimitives; /**< The number of primitives in each basisfunction. */
  double exps[3];
  double coefs[3];
  double origin[3];
  int shell[3];
  double norm;
  double massno;
} bfn;

void bfn_set(bfn*, int, double*, double*, int*, double, double);

void normalise(bfn*);
