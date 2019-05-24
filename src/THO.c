#include "THO.h"
#include "utils.h"
#include <math.h>
#include <cblas.h>


void calculate_integral(sys* mol){
  printf("Computing Overlap Integrals...\n");
  for(int i = 0; i < mol->nbfs; i++){
    for(int j = 0; j < mol->nbfs; j++){
      mol->S[i*mol->nbfs + j] = Sab(mol->basisfunctions[i], mol->basisfunctions[j]);
    }
  }
  print_matrix(mol->S, mol->nbfs, mol->nbfs);
  printf("Computing Kinetic Integrals...\n");
  for(int i = 0; i < mol->nbfs; i++){
    for(int j = 0; j < mol->nbfs; j++){
      mol->T[i*mol->nbfs + j] = Tab(mol->basisfunctions[i], mol->basisfunctions[j]);
    }
  }
  print_matrix(mol->T, mol->nbfs, mol->nbfs);
}

double f(int j, int l, int m, double PA, double PB){
  double total = 0.0;
  for (int k = 0; k < j + 1; k++){
    if ((j - l <= k) && (k <= m)){
      total += binomial(l, j - k) * binomial(m, k) * pow(PA, l - j + k) * pow(PB, m - k);
    }
  }
  return total;
}

double overlap_1d(int l1, int l2, double PAx, double PBx, double gamma){
  double total = 0;
  for (int j = 0; j < floor((l1 + l2) * 0.5) + 1; j++){
    total += f(2*j, l1, l2, PAx, PBx) * double_factorial(2*j - 1) / pow(2*gamma, j);
  }
  return total;
}

double overlap(int* lmn1, double* A, double a, int* lmn2, double* B, double b){
  double gamma = a + b;
  double Q[3];
  /* broadcast_vv(A, B, 3, Q); */
  /* cblas_dscal(3, a*b/gamma, Q, 1); */
  Q[0] = (a * A[0] + b * B[0]) / gamma;
  Q[1] = (a * A[1] + b * B[1]) / gamma;
  Q[2] = (a * A[2] + b * B[2]) / gamma;
  double dist = dist2(A, B);
  double S_x = overlap_1d(lmn1[0], lmn2[0], Q[0] - A[0], Q[0] - B[0], gamma);
  double S_y = overlap_1d(lmn1[1], lmn2[1], Q[1] - A[1], Q[1] - B[1], gamma);
  double S_z = overlap_1d(lmn1[2], lmn2[2], Q[2] - A[2], Q[2] - B[2], gamma);
  return pow(M_PI/gamma, 1.5) * exp(-a * b * dist / gamma) * S_x * S_y * S_z;
}

double Sab(bfn bf1, bfn bf2){
  double total = 0;
  for (int i = 0; i < bf1.nprimitives; i++){
    for (int j = 0; j < bf2.nprimitives; j++){
      total +=  bf1.norm[i] * bf2.norm[j] * bf1.coefs[i] * bf2.coefs[j] * overlap(bf1.shell, bf1.origin, bf1.exps[i], bf2.shell, bf2.origin, bf2.exps[j]);
    }
  }
  return total;
}

double kinetic(int* lmn1, double* A, double a, int* lmn2, double* B, double b){
  int l2 = lmn2[0]; int m2 = lmn2[1]; int n2 = lmn2[2];
  int pl2mn[3] = { l2 + 2, m2, n2 }; int plm2n[3] = { l2, m2 + 2, n2 }; int plmn2[3] = { l2, m2, n2 + 2 };
  int ml2mn[3] = { l2 - 2, m2, n2 }; int mlm2n[3] = { l2, m2 - 2, n2 }; int mlmn2[3] = { l2, m2, n2 - 2 };

  double term1 = b * (2 * (l2 + m2 + n2) + 3) * overlap(lmn1, A, a, lmn2, B, b);
  double term2 = -2 * b * b * (
			       overlap(lmn1, A, a, pl2mn, B, b) +
			       overlap(lmn1, A, a, plm2n, B, b) +
			       overlap(lmn1, A, a, plmn2, B, b) );
  double term3 = -0.5 * (
			 (l2 * (l2 - 1) * overlap(lmn1, A, a, ml2mn, B, b)) +
			 (m2 * (m2 - 1) * overlap(lmn1, A, a, mlm2n, B, b)) +
			 (n2 * (n2 - 1) * overlap(lmn1, A, a, mlmn2, B, b)) );

  return term1 + term2 + term3;
}

double Tab(bfn bf1, bfn bf2){
  double total = 0;
  for (int i = 0; i < bf1.nprimitives; i++){
    for (int j = 0; j < bf2.nprimitives; j++){
      total +=  bf1.norm[i] * bf2.norm[j] * bf1.coefs[i] * bf2.coefs[j] * kinetic(bf1.shell, bf1.origin, bf1.exps[i], bf2.shell, bf2.origin, bf2.exps[j]);
    }
  }
  return total;
}
