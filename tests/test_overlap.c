#include <math.h>
#include "THO.h"

int compare(double computed, double ref, double thresh) {
    const double diff = fabs(computed - ref);
    printf("(%lf, %lf, %lf)\n", computed, ref, diff);
    int retval = 0;
    if (diff > thresh)
        retval = 1;
    return retval;
}

int overlap0() {
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

int overlap1() {
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
    return overlap0()
        | overlap1()
        | 0;
}