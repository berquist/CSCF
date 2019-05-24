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
  printf("Computing Nuclear Integrals...\n");
  for(int i = 0; i < mol->nbfs; i++){
    for(int j = 0; j < mol->nbfs; j++){
      for(int k = 0; k < mol->nbfs; k++){
	mol->V[i*mol->nbfs + j] += -mol->basisfunctions[k].massno * Vab(mol->basisfunctions[i], mol->basisfunctions[j], mol->basisfunctions[k].origin);
      }
    }
  }
  print_matrix(mol->V, mol->nbfs, mol->nbfs);
}

double f(int j, int l, int m, double PA, double PB){
  double total = 0.0;
  for (int k = max(0, j-m); k < min(j, l) + 1; k++){
    total += binomial(l, k) * binomial(m, j-k) * pow(PA, l - k) * pow(PB, m + k - j);
  }
  return total;
}

double overlap_1d(int l1, int l2, double PAx, double PBx, double gamma){
  double total = 0.0;
  for (int j = 0; j < floor((l1 + l2) * 0.5) + 1; j++){
    total += f(2*j, l1, l2, PAx, PBx) * double_factorial(2*j - 1) / pow(2*gamma, j);
  }
  return total;
}

double overlap(int* lmn1, double* A, double a, int* lmn2, double* B, double b){
  double gamma = a + b;
  double Q[3];
  gpc(A, a, B, b, Q); 
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
  int  l2 = lmn2[0]; int m2 = lmn2[1]; int n2 = lmn2[2];
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


double A(int i, int r, int u, int l1, int l2, double PAx, double PBx, double CPx, double gamma){
  double total = pow(-1, i) * f(i, l1, l2, PAx, PBx) * pow(-1, u) * factorial(i)*pow(CPx, i-2*r-2*u) * pow(0.25/gamma, r+u) / factorial(r) / factorial(u) / factorial(i-2*r-2*u);
  return total;
}

void G_array(int l1, int l2, double PAx, double PBx, double CPx, double gamma, double* out){
  int imax = l1 + l2 + 1;
  for(int i = 0; i < imax; i++){
    for(int r = 0; r <= i/2; r++){
      for(int u = 0; u <=(i - 2* r) * 0.5; u++){
	int iI = i - 2 * r - u;
	out[iI] = A(i, r, u, l1, l2, PAx, PBx, CPx, gamma);
      }
    }
  }
}

double nuclear(int* lmn1, double* A, double a, int* lmn2, double* B, double b, double* C){
  int  l1 = lmn1[0]; int m1 = lmn1[1]; int n1 = lmn1[2];
  int  l2 = lmn2[0]; int m2 = lmn2[1]; int n2 = lmn2[2];
  double Q[3];
  gpc(A, a, B, b, Q);
  double gamma = a + b;
  double distPC = dist2(C, Q);
  double distAB = dist2(A, B);
  double Gx[l1+l2+1], Gy[m1+m2+1], Gz[n1+n2+1];
  G_array(l1, l2, Q[0] - A[0], Q[0] - B[0], Q[0] - C[0], gamma, Gx);
  G_array(m1, m2, Q[1] - A[1], Q[1] - B[1], Q[1] - C[1], gamma, Gy);
  G_array(n1, n2, Q[2] - A[2], Q[2] - B[2], Q[2] - C[2], gamma, Gz);

  double sum = 0.0;

  for(int i = 0; i <= l1+l2; i++){
    for(int j = 0; j <= m1+m2; j++){
      for(int k = 0; k <= n1+n2; k++){
	sum += Gx[i] * Gy[j] * Gz[k] * boys(i+j+k, distPC*gamma);
      }
    }
  }
  return sum * exp(-a * b * distAB / gamma) * 2 * M_PI / gamma;
}

double Vab(bfn bf1, bfn bf2, double* C){
  double total = 0;
  for (int i = 0; i < bf1.nprimitives; i++){
    for (int j = 0; j < bf2.nprimitives; j++){
      total +=  bf1.norm[i] * bf2.norm[j] * bf1.coefs[i] * bf2.coefs[j] * nuclear(bf1.shell, bf1.origin, bf1.exps[i], bf2.shell, bf2.origin, bf2.exps[j], C);
    }
  }
  return total;
}

