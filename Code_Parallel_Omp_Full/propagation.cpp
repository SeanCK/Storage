#include "propagation.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/// FUNCTIONS
///////////////////////////////////////////////////////////////////////////////

/*! Adiabatic Propagator */
double U(double *r,double *v, int Sa, double t, double *f, int N_bath, double timestep, double ddd, double ddd4, double *mww, double *c,
         void (*force[4])(double *, double *, int, double, double*, double*)){
    double phase, dt;
    int Nsteps;

    force[Sa](r, f, N_bath, ddd4, mww, c);
    dt = timestep;

    if (t <= timestep){
        dt = t;
        Nsteps = 1;
    }
    else{
        Nsteps = t/dt +1;
        dt = t/Nsteps;
    }

    if ((Sa == 0) || (Sa == 3)){
        for (int i = 0; i < Nsteps; ++i){
            integ_step(r , v,  dt, Sa, f, N_bath, ddd4, mww, c, force);
        }
        return 0.0;
    }
    phase = dE(r, N_bath, ddd, c)*0.5;
    for (int i = 0; i < Nsteps; ++i){
        integ_step(r , v,  dt, Sa, f, N_bath, ddd4, mww, c, force);
        phase += dE(r, N_bath, ddd, c);
    }
    phase -=dE(r, N_bath, ddd, c)*0.5;
    phase *= dt;

    if (Sa == 1)
        phase *= -1.0;

    return phase;
}

/*! Pure Bath Force Field */
void Fb(double *R, double *f, int N_bath, double ddd4, double *mww, double *c){
    for (int i= 0; i < N_bath; ++i)
        f[i] = mww[i]*R[i];
}

/*! 00 force field */
void F1(double *R, double *f, int N_bath, double ddd4, double *mww, double *c){
    double g,h;
    g = gam(R, N_bath, c);
    h = g/sqrt(ddd4 + g*g);
    for (int i = 0; i < N_bath; ++i){
        f[i]  = mww[i]*R[i] - h*c[i];
    }
}

/*! 11 force field */
void F2(double *R, double *f, int N_bath, double ddd4, double *mww, double *c){
    double g,h;
    g = gam(R, N_bath, c);
    h = g/sqrt(ddd4 + g*g);
    for (int i = 0; i < N_bath; ++i)
        f[i] = mww[i]*R[i] + h*c[i];
}

/*! Velocity Verlet */
void integ_step(double *r, double *v, double dt, int Sa, double *f, int N_bath, double ddd4, double *mww, double *c,
                void (*force[4])(double *, double *, int, double, double*, double*)){
    double y;
    y = 0.5*dt*dt;
    for (int i = 0; i < N_bath; ++i)
        r[i] += dt*v[i] + y*f[i];
    y = 0.5*dt;
    for (int i = 0; i < N_bath; ++i)
        v[i] += y*f[i];
    force[Sa](r, f, N_bath, ddd4, mww, c);
    for (int i = 0; i < N_bath; ++i)
        v[i] += y*f[i];
}

/*! Energy difference between adiabatic surface (E1 - E0) */
double dE(double *R, int N_bath, double ddd, double *c){
    double g;
    g = gam(R, N_bath, c);
    g *= 4.0*g;
    return (sqrt(ddd + g));
}

/*! Total coupling between quantum and classical systems */
double gam(double *R, int N_bath, double *c){
    double x = 0.0;
    double asyEps = 0.4;

    for (int i = 0; i < N_bath; ++i)
        x += c[i]*R[i];
    x += asyEps;    // asymmetric spin boson
    return -x;
}

/*! Convenience function for calculating eigenvectors of density matrix in adiabatic basis */
double G(double *R, int N_bath, double delta, double ddd, double *c){
    double x,g;
    g = gam(R, N_bath, c);
    if (fabs(g/delta) < 1.0e-7)
        return (g/delta);
    x = (-delta + sqrt(ddd + 4*g*g))/(2*g);
    return x;
}

///////////////////////////////////////////////////////////////////////////////
/// Observables and initial density Matrices
///////////////////////////////////////////////////////////////////////////////

/*! Definition of initial density matrix element */
double dens_init_0(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g,gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z = 0.5*(1.0 + 2.0*g + gg)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_1(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g,gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_2(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g, gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_3(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g,gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z = 0.5*(gg - 2*g  + 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double wigner_harm_osc(double *x, double *p){
    return 1.0;
}

double obs_0(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g,gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z = 2.0*g/(1 + gg);
    return z;
}

double obs_1(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g,gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z = (gg-1)/(1 + gg);
    return z;
}

double obs_2(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g, gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z =  (gg-1)/(1 + gg);
    return z;
}

double obs_3(double *x,double *p, int N_bath, double delta, double ddd, double *c){
    double z;
    double g,gg;
    g = G(x, N_bath, delta, ddd, c); gg = g*g;
    z = -2.0*g/(1 + gg);
    return z;
}

/*! Matrix elements of the Hamiltonian */
double H_0(double *x,double *p, int N_bath, double ddd, double *mww, double *c){
    double z;
    z = Hb(x, p, N_bath, mww) - dE(x, N_bath, ddd, c)*0.5;
    return z;
}

double H_1(double *x,double *p, int N_bath, double ddd, double *mww, double *c){
    return 0.0;
}

double H_2(double *x,double *p, int N_bath, double ddd, double *mww, double *c){
    return 0.0;
}

double H_3(double *x,double *p, int N_bath, double ddd, double *mww, double *c){
    double z;
    z = Hb(x, p, N_bath, mww) + dE(x, N_bath, ddd, c)*0.5;
    return z;
}

/*! Bath Hamiltonian */
double Hb(double *R, double *P, int N_bath, double *mww){
    double x = 0.0;
    for (int i = 0; i < N_bath; ++i)
        x += P[i]*P[i] - mww[i]*R[i]*R[i];
    return x*0.5;
}