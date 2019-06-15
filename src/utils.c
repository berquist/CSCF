#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_hyperg.h>

double min(double a, double b){
  if (a < b){
    return a;
  } else {
    return b;
  }
}

double max(double a, double b){
  if (a > b){
    return a;
  } else {
    return b;
  }
}

int factorial(int n){
  if (n <= 1){
    return 1;
  } else {
    return n*factorial(n-1);
  }
}

int double_factorial(int n){
  if (n <= 1){
    return 1;
  } else {
    return n*double_factorial(n-2);
  }
}

int binomial(int n, int k){
  if(n == k){
    return 1;
  } else {
    return factorial(n) / (factorial(k) * factorial(n-k));
  }
}

void broadcast_vv(double* vec1, double* vec2, int N, double* out_vec){
  for(int i = 0; i < N; i++){
      out_vec[i] = vec1[i] +  vec2[i];
  }
}

double dist2(double* vec1, double* vec2){
  double pre = pow(vec1[0] - vec2[0], 2.0) + pow(vec1[1] - vec2[1], 2.0) + pow(vec1[2] - vec2[2], 2.0);
  return pre;
}

double distd(double* vec1, double* vec2){
  double pre = pow(vec1[0] - vec2[0], 2.0) + pow(vec1[1] - vec2[1], 2.0) + pow(vec1[2] - vec2[2], 2.0);
  return sqrt(pre);
}
  
void print_matrix(double* mat, int M, int K){
  for(int i = 0; i < M; i++){
    for(int j = 0; j < K; j++){
      printf("% .4f ", mat[i * K + j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_matrixi(int* mat, int M, int K){
  for(int i = 0; i < M; i++){
    for(int j = 0; j < K; j++){
      printf("%d ", mat[i * K + j]);
    }
    printf("\n");
  }
  printf("\n");
}

void gpc(double* A, double a, double* B, double b, double* Q){
  Q[0] = (a*A[0] + b*B[0]) / (a+b);
  Q[1] = (a*A[1] + b*B[1]) / (a+b);
  Q[2] = (a*A[2] + b*B[2]) / (a+b);
}

double boys(double n, double T){
  return gsl_sf_hyperg_1F1(n+0.5, n+1.5, -T) / (2.0*n + 1.0);
}
