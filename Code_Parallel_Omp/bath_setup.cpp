#include "bath_setup.h"

using namespace std;

#define PI 3.141592653589793

///////////////////////////////////////////////////////////////////////////////
/// FUNCTIONS
///////////////////////////////////////////////////////////////////////////////

/*! Parameters for bath (corresponding to an ohmic spectral density) */
void bath_para(double eta, double w_max, int N_bath, double *c, double *m, double *w){
    double w_0;
    w_0 = (1 - exp(-w_max))/N_bath;
    for (int i = 0; i < N_bath; ++i){
        m[i] = 1.0;
        w[i] = -log( 1-(i+1)*w_0 );
        c[i] = sqrt(eta*w_0*m[i])*w[i];
    }
}

///////////////////////////////////////////////////////////////////////////////
/// RANDOM NUMBER GENERATOR
///////////////////////////////////////////////////////////////////////////////

/*!< Gaussian number generator for (R,P) */
void gauss_init_W(double *R, double *v, gsl_rng * rr, double * ranVector, int N_bath, double *sig){
    double sigma_x, sigma_v;
    randnums(4*N_bath, ranVector, rr);
    for (int i = 0; i < N_bath; ++i){
        sigma_x = sig[i];
        sigma_v = sig[i+N_bath];
        R[i] = gauss1(sigma_x,i, rr, ranVector, N_bath);
        v[i] = gauss1(sigma_v,i + 2*N_bath, rr, ranVector, N_bath);
    }
}

void randnums(int rand_dim, double *rand_vec, gsl_rng * rr){
    for (int i = 0; i < rand_dim; ++i){
        rand_vec[i] = gsl_rng_uniform (rr);
    }
}

double  gauss1(double sigma_x, int i, gsl_rng * rr, double * ranVector, int N_bath){
    double x1,y1,y2,y3;
    y1 = ranVector[i];
    while (fabs(y1) < 1.0e-200){
        y1 = gsl_rng_uniform (rr);
    }
    y2 = ranVector[i+N_bath];
    y3 = sqrt(-2*log(y1));
    x1 = y3*cos(2*PI*y2);
    return (sigma_x*x1);
}