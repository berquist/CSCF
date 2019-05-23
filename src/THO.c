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
}

double f(double j, double l, double m, double PA, double PB){
  double total = 0;
  for (int k = 0; k < j + 1; k++){
    if ((j - m <= k) && (k <= l)){
      total += binomial(l, k) * binomial(m, j-k) * pow(PA, l-k) * pow(PB, m+k-j);
    }
  }
  return total;
}

double overlap_1d(int l1, int l2, double PAx, double PBx, double gamma){
  double total = 0;
  for (double j = 0; j < floor((l1 + l2) / 2) + 1; j++){
    total += f(2*j, l1, l2, PAx, PBx) * double_factorial(2*j - 1) / pow(2*gamma, j);
  }
}

double overlap(int* lmn1, double* A, double a, int* lmn2, double* B, double b){
  double gamma = a + b;
  double Q[3];
  broadcast_vv(A, B, 3, Q);
  cblas_dscal(3, a*b/gamma, Q, 1);
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
