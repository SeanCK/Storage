#ifndef TRANSITION_MATRIX
#define TRANSITION_MATRIX

#include "density.h"

double dgamma(double *R, int i, double *c);

double dd(double *dhat, double *R, double *dgam, double ddd, int N_bath, double delta, double *c);

void setwww(double (*www[2][4][4])(double, double, double, double));

#endif
