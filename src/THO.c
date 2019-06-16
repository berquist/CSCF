#include "THO.h"
#include "utils.h"
#include <math.h>
#include <cblas.h>
#include <stdlib.h>


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
	mol->V[i*mol->nbfs + j] += mol->basisfunctions[k].atomno * Vab(mol->basisfunctions[i], mol->basisfunctions[j], mol->basisfunctions[k].origin);
      }
    }
  }
  print_matrix(mol->V, mol->nbfs, mol->nbfs);
  printf("Computing Replusion Integrals...\n");
  for(int i = 0; i < mol->nbfs; i++){
    for(int j = 0; j < mol->nbfs; j++){
      for(int k = 0; k < mol->nbfs; k++){
	for(int l = 0; l < mol->nbfs; l++){
	  mol->ERI[i*mol->nbfs * mol->nbfs * mol->nbfs + j * mol->nbfs * mol->nbfs + k * mol->nbfs + l] += ERIabcd(mol->basisfunctions[i], mol->basisfunctions[j], mol->basisfunctions[k], mol->basisfunctions[l]);
	}
      }
    }
  }
  print_matrix(mol->ERI, 8, 2);
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
  double epsilon = 1 / (4*gamma);
  double num = pow(-1, i) * f(i, l1, l2, PAx, PBx) *  pow(-1, u) * factorial(i) * pow(CPx, i - 2*r - 2*u) * pow(epsilon, r+u);
  double denom = factorial(r) * factorial(u) * factorial(i - 2*r - 2*u);
  return num/denom;
}

double* G_array(int l1, int l2, double PAx, double PBx, double CPx, double gamma){
  int l1l2 = l1 + l2;
  double* out = malloc(sizeof (double) * l1l2+1);
  if (l1l2 == 0){
    out[0] = 1.0;
  } else if (l1l2 == 1){
    out[0] = f(0, l1, l2, PAx, PBx);
    out[1] = -CPx;
  } else if (l1l2 == 2){
    out[0] = f(0, l1, l2, PAx, PBx) + (f(1, l1, l2, PAx, PBx) / (2 * gamma));
    out[1] = -f(1, l1, l2, PAx, PBx) * CPx - (f(1, l1, l2, PAx, PBx) / (2 * gamma));
    out[2] = CPx * CPx;
  }
  /* for(int i = 0; i <= l1 + l2; i++){
    for(int r = 0; r <= floor(i/2); r++){
      for(int u = 0; u <=floor((i - 2* r) * 0.5); u++){
	int iI = i - 2 * r - u;
	out[iI] += A(i, r, u, l1, l2, PAx, PBx, CPx, gamma);
      }
    }
  }*/
  return out;
}

double nuclear(int* lmn1, double* A, double a, int* lmn2, double* B, double b, double* C){
  int  l1 = lmn1[0]; int m1 = lmn1[1]; int n1 = lmn1[2];
  int  l2 = lmn2[0]; int m2 = lmn2[1]; int n2 = lmn2[2];
  double Q[3];
  gpc(A, a, B, b, Q);
  double gamma = a + b;
  double distPC = dist2(C, Q);
  double distAB = dist2(A, B);
  double* Gx = G_array(l1, l2, Q[0] - A[0], Q[0] - B[0], Q[0] - C[0], gamma);
  double* Gy = G_array(m1, m2, Q[1] - A[1], Q[1] - B[1], Q[1] - C[1], gamma);
  double* Gz = G_array(n1, n2, Q[2] - A[2], Q[2] - B[2], Q[2] - C[2], gamma);

  double sum = 0.0;
  for(int i = 0; i <= l1+l2; i++){
    for(int j = 0; j <= m1+m2; j++){
      for(int k = 0; k <= n1+n2; k++){
	sum += Gx[i] * Gy[j] * Gz[k] * boys(i+j+k, distPC*gamma);
      }
    }
  }
  free(Gx);
  free(Gy);
  free(Gz);

  return sum * exp(-a * b * distAB / gamma) * -2 * M_PI / gamma;
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

double* H_array(int l1, int l2, double PAx, double PBx, double gamma){
  int imax = l1+l2+1;
  double* out = malloc(sizeof (double) * imax);
  if (l1 == 0 && l2 == 0){
    out[0] = 1.0;
  } else if (l1 == 1 && l2 == 0){
    out[0] = PAx;
    out[1] = 1/(4*gamma);
  }
  return out;
}

double* C_array(int l1, int l2, double PAx, double PBx, double gamma1, int l3, int l4, double QCx, double QDx, double gamma2, double QPx){
  double delta = ((1 / (4 * gamma1)) + (1 / (4 * gamma2)));
  double* H1 = H_array(l1, l2, PAx, PBx, gamma1);
  double* H2 = H_array(l3, l4, QCx, QDx, gamma2);
  double* out = malloc(sizeof (double) * l1+l2+l3+l4+1);
  for (int L = 0; L <= l1 + l2; L++){
    for(int M = 0; M <= l3 + l4; M++){
      for(int u = 0; u <= floor((L+M)/2); u++){
	int I = L+M-u;
	out[I] = H1[L] * pow(-1, M) * H2[M] * factorial(L+M) * pow(-1, u) * pow(QPx, (L+M) - (2*u)) / (factorial(u) * factorial((L + M) - (2 * u)) * pow(delta, L+M-u));
      }
    }
  }
  return out;
}

double electron(int* lmn1, double* A, double a,
		int* lmn2, double* B, double b,
		int* lmn3, double* C, double c,
		int* lmn4, double* D, double d){
  int  l1 = lmn1[0]; int m1 = lmn1[1]; int n1 = lmn1[2];
  int  l2 = lmn2[0]; int m2 = lmn2[1]; int n2 = lmn2[2];
  int  l3 = lmn3[0]; int m3 = lmn3[1]; int n3 = lmn3[2];
  int  l4 = lmn4[0]; int m4 = lmn4[1]; int n4 = lmn4[2];
  double Q[3]; double P[3];
  gpc(A, a, B, b, P);
  gpc(C, c, D, d, Q);
  double gamma1 = a + b;
  double gamma2 = c + d;
  double delta = ((1 / (4 * gamma1)) + (1 / (4 * gamma2)));
  double distPQ = dist2(P, Q);
  double distAB = dist2(A, B);
  double distCD = dist2(C, D);
  double* Cx = C_array(l1, l2, P[0] - A[0], P[0] - B[0], gamma1, l3, l4, Q[0] - C[0], Q[0] - D[0], gamma2, Q[0] - P[0]);
  double* Cy = C_array(m1, m2, P[1] - A[1], P[1] - B[1], gamma1, m3, m4, Q[1] - C[1], Q[1] - D[1], gamma2, Q[1] - P[1]);
  double* Cz = C_array(n1, n2, P[2] - A[2], P[2] - B[2], gamma1, n3, n4, Q[2] - C[2], Q[2] - D[2], gamma2, Q[2] - P[2]);
  double sum = 0.0;
  for(int i = 0; i <= l1 + l2 + l3 + l4; i++){
    for(int j = 0; j <= m1 + m2 + m3 + m4; j++){
      for(int k = 0; k <= n1 + n2 + n3 + n4; k++){
	sum += Cx[i] * Cy[j] * Cz[k] * boys(i+j+k, distPQ / (4 * delta));
      }
    }
  }
  return sum * exp((-a * b * distAB / gamma1)) * exp((-c * d * distCD / gamma2)) * 2 * pow(M_PI, 2) * (1 / (gamma1 * gamma2)) * pow(M_PI * (1 / (gamma1 + gamma2)), 0.5);
}

double ERIabcd(bfn bfn1, bfn bfn2, bfn bfn3, bfn bfn4){
  double total = 0.0;
  for(int i = 0; i < bfn1.nprimitives; i++){
    for(int j = 0; j < bfn2.nprimitives; j++){
      for(int k = 0; k < bfn3.nprimitives; k++){
	for(int l = 0; l < bfn4.nprimitives; l++){
	  total +=  bfn1.norm[i] * bfn2.norm[j] * bfn3.norm[k] * bfn4.norm[l] *
	            bfn1.coefs[i] * bfn2.coefs[j] * bfn3.coefs[k] * bfn4.coefs[l] *
	    electron(bfn1.shell, bfn1.origin, bfn1.exps[i],
		     bfn2.shell, bfn2.origin, bfn2.exps[j],
		     bfn3.shell, bfn3.origin, bfn3.exps[k],
		     bfn4.shell, bfn4.origin, bfn4.exps[l]);
	}
      }
    }
  }
  return total;
}
