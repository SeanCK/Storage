#ifndef PROPAGATION
#define PROPAGATION

#include "transition_matrix.h"

double U(double *r,double *v, int Sa, double t, double *f, int N_bath, double timestep, double ddd, double ddd4, double *mww, double *c,
         void (*force[4])(double *, double *, int, double, double*, double*));

void Fb(double *R, double *f, int N_bath, double ddd4, double *mww, double *c);

void F1(double *R, double *f, int N_bath, double ddd4, double *mww, double *c);

void F2(double *R, double *f, int N_bath, double ddd4, double *mww, double *c);

void integ_step(double *r, double *v, double dt, int Sa, double *f, int N_bath, double ddd4, double *mww, double *c,
                void (*force[4])(double *, double *, int, double, double*, double*));

double dE(double *R, int N_bath, double ddd, double *c);

double gam(double *R, int N_bath, double *c);

double G(double *R, int N_bath, double delta, double ddd, double *c);

double dens_init_0(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double dens_init_1(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double dens_init_2(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double dens_init_3(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double wigner_harm_osc(double *x, double *p);

double obs_0(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double obs_1(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double obs_2(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double obs_3(double *x,double *p, int N_bath, double delta, double ddd, double *c);

double H_0(double *x,double *p, int N_bath, double ddd, double *mww, double *c);

double H_1(double *x,double *p, int N_bath, double ddd, double *mww, double *c);

double H_2(double *x,double *p, int N_bath, double ddd, double *mww, double *c);

double H_3(double *x,double *p, int N_bath, double ddd, double *mww, double *c);

double Hb(double *R, double *P, int N_bath, double *mww);

#endif