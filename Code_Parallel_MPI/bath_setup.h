#ifndef BATH_SETUP
#define BATH_SETUP

#include <math.h>
#include <gsl/gsl_rng.h>

void bath_para(double eta, double w_max, int N_bath, double *c, double *m, double *w);

void gauss_init_W(double *R, double *v, gsl_rng * rr, double * ranVector, int N_bath, double *sig);

void randnums(int rand_dim, double *rand_vec, gsl_rng * rr);

double  gauss1(double sigma_x, int i, gsl_rng * rr, double * ranVector, int N_bath);

#endif