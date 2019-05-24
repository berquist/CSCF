#include <math.h>
#include "THO.h"
#include "utils.h"

int compare(double computed, double ref, double thresh) {
  const double diff = fabs(computed - ref);
  printf("(%lf, %lf, %lf)\n", computed, ref, diff);
  int retval = 0;
  if (diff > thresh)
    retval = 1;
  return retval;
}

int test_double_factorial() {
  printf("double_factorial\n");
  if (double_factorial(0) != 1) {
    printf("failure in double_factorial(0)\n");
    return 1;
  }
  if (double_factorial(1) != 1) {
    printf("failure in double_factorial(1)\n");
    return 1;
  }
  if (double_factorial(2) != 2) {
    printf("failure in double_factorial(2)\n");
    return 1;
  }
  if (double_factorial(3) != (3 * 1)) {
    printf("failure in double_factorial(3)\n");
    return 1;
  }
  if (double_factorial(4) != (4 * 2)) {
    printf("failure in double_factorial(4)\n");
    return 1;
  }
  if (double_factorial(5) != (5 * 3 * 1)) {
    printf("failure in double_factorial(5)\n");
    return 1;
  }
  return 0;
}

int test_f() {
  printf("f\n");
  double thresh = 1.0e-7;

  if (compare(f(1, 1, 1, 0.1, 0.2), 0.3, thresh))
    return 1;
  if (compare(f(1, 1, 1, 0.3, 0.4), 0.7, thresh))
    return 1;
  if (compare(f(1, 3, 1, 0.1, 0.2), 0.007, thresh))
    return 1;
  if (compare(f(2, 3, 1, 0.1, 0.2), 0.09, thresh))
    return 1;

  return 0;
}

int test_overlap1d() {
  printf("overlap1d\n");
  double thresh = 1.0e-5;

  if (compare(overlap_1d(1, 1, 0.1, 0.2, 1.0), 0.52, thresh))
    return 1;
  if (compare(overlap_1d(3, 1, 0.1, 0.2, 1.0), 0.7952, thresh))
    return 1;

  return 0;
}

int test_full_overlap_1() {
  printf("full_overlap_1\n");
  int lmn1[3] = {0, 0, 0};
  double A[3] = {0.0, 0.0, 0.0};
  double a = 1.8;
  int lmn2[3] = {0, 0, 0};
  double B[3] = {0.5, 0.8, -0.2};
  double b = 2.8;
  double computed = overlap(lmn1, A, a, lmn2, B, b);
  double ref = 0.20373275913014607;
  double thresh = 1.0e-12;
  return compare(computed, ref, thresh);
}

int test_full_overlap_2() {
  printf("full_overlap_2\n");
  int lmn1[3] = {1, 0, 0};
  double A[3] = {0.0, 0.0, 0.0};
  double a = 1.8;
  int lmn2[3] = {0, 0, 0};
  double B[3] = {0.5, 0.8, -0.2};
  double b = 2.8;
  double computed = overlap(lmn1, A, a, lmn2, B, b);
  double ref = 0.062005622343957505;
  double thresh = 1.0e-12;
  return compare(computed, ref, thresh);
}

int main() {
  return test_double_factorial() |
    test_f() |
    test_overlap1d() |
    test_full_overlap_1() |
    test_full_overlap_2() |
    0;
}
