#ifndef DENSITY
#define DENSITY

#include "propagation.h"
#include "bath_setup.h"
#include <iostream>
#include <complex>
#include <vector>

int density(double *x,double *p, double *f, double *abszsum1, double *argzsum1, double *habszsum1, double *hargzsum1, double **realsum, double **imagsum,
            double **hrealsum, double **himagsum, gsl_rng * rr, double * ranVector, int N_bath, int Ncut, double *m, int  N_slice, double TSLICE,
            double timestep, double delta, double ddd, double ddd4, double *mww, double *c, double *sig, double (*dens_init[4])(double*, double*, int, double, double, double *),
            double (*obs[4])(double*, double*, int, double, double, double *), double (*obs1[4])(double*, double*, int, double, double*, double *),
            void (*force[4])(double *, double *, int, double, double*, double*), double (*www[2][4][4])(double, double, double, double), double **full_abszsum1, int n);

#endif
