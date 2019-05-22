#include <math.h>

double factorial(double n){
  if (n <= 1){
    return 1;
  } else {
    return n*factorial(n-1);
  }
}

double double_factorial(double n){
  if (n <= 1){
    return 1;
  } else {
    return n*double_factorial(n-2);
  }
}

double binomial(double n, double k){
  if(n == k){
    return 1;
  } else {
    return factorial(n) / factorial(k) * factorial(n-k);
  }
}

void broadcast_vv(double* vec1, double* vec2, int N, double* out_vec){
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; i++){
      out_vec[i] = vec1[i] +  vec2[i];
    }
  }
}

double dist2(double* vec1, double* vec2){
  double pre = pow(vec2[0] - vec1[0], 2.0) + pow(vec2[1] - vec1[1], 2.0) + pow(vec2[2] - vec1[1], 2.0);
  return pre;
}

