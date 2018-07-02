#include "density.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/// PROCESSING TREE
///////////////////////////////////////////////////////////////////////////////

int  density(double *x,double *p, double *f, double *abszsum1, double *argzsum1, double *habszsum1, double *hargzsum1, double **realsum, double **imagsum,
             double **hrealsum, double **himagsum, gsl_rng * rr, double *ranVector, int N_bath, int Ncut, double *m, int  N_slice, double TSLICE,
             double timestep, double delta, double ddd, double ddd4, double *mww, double *c, double *sig, double (*dens_init[4])(double*, double*, int, double, double, double *),
             double (*obs[4])(double*, double*, int, double, double, double *), double (*obs1[4])(double*, double*, int, double, double*, double *),
             void (*force[4])(double *, double *, int, double, double*, double*), double (*www[2][4][4])(double, double, double, double), double **full_abszsum1, int n){

    int SS0,SS1,SS2,SS3,NNjmp = 0,signPdotdhat;
    double phase0 = 0.0,xx;
    double p0,p1,p2,p3,ap0,ap1,ap2,ap3;
    double dn2;
    double *Pperp;
    double *dhat;
    complex<double> z(1,0);
    complex<double> oldz;
    complex<double> initd(1,0);
    complex<double> I(0,1);
    double (*phi)(double*, double*, int, double, double, double *);
    double (*phi1)(double*, double*, int, double, double*, double *);
    double *RR;
    double *PP;
    double abs_d;
    RR = new double[N_bath];
    PP = new double[N_bath];
    dhat = new double[N_bath];
    Pperp = new double[N_bath];

    double Pdotdhat; /*!< Parallel component of momentum*/
    double sina;
    double cosa;
    double de;

    double *dgam; /*!< */
    dgam = new double[N_bath];


    ///////////////////////////////////////////////////////////////////////////////
    /// INITIALIZATION OF INITIAL SURFACE
    ///////////////////////////////////////////////////////////////////////////////

    gauss_init_W(x, p, rr, ranVector, N_bath, sig); /*!< Creates initial random sampling */
    double yy = 4.0*(gsl_rng_uniform (rr));
    /*! Setting initial surface value */
    if (yy < 1.0)
        SS3 = (SS0 = 0);
    else if (yy < 2.0){
        SS0 = 1;
        SS3 = 2;
    }
    else if (yy < 3.0){
        SS0 = 2;
        SS3 = 1;
    }
    else
        SS3 = (SS0 = 3);
    initd = dens_init[SS3](x, p, N_bath, delta, ddd, c);
    z = 4.0;
    /*! Allocating values for position and momentum from Gaussian */
    for (int l = 0; l < N_bath; ++l){
        RR[l] = x[l];
        PP[l] = p[l];
    }
    SS1 = SS0;

    ///////////////////////////////////////////////////////////////////////////////
    /// ITERATION FOR TREE: CALCULATE EACH PATH SEGMENT
    ///////////////////////////////////////////////////////////////////////////////

    for (int l = 0; l < N_slice; ++l) {
      
        SS0 = SS1; /*!< Sets beginning surface value */

        ///////////////////////////////////////////////////////////////////////////////
        /// ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        phase0 = U(RR, PP, SS0, TSLICE*0.5, f, N_bath, timestep, ddd, ddd4, mww, c, force); // exp(iLd/2) (before jump)
        z *= exp(I * phase0);

        ///////////////////////////////////////////////////////////////////////////////
        /// NON-ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        abs_d = dd(dhat, RR, dgam, ddd, N_bath, delta, c); /*!< calculating non-adiabatic coupling matrix */
        de = dE(RR, N_bath, ddd, c); /*!< calculating energy */
        double alpha = 0.0;
        Pdotdhat = 0;
        for (int i = 0; i < N_bath; ++i) {
            Pdotdhat += PP[i] * dhat[i];
            alpha += PP[i] * dhat[i] / m[i];
        }
        alpha *= abs_d * TSLICE;
        signPdotdhat = (Pdotdhat < 0 ? -1 : 1); // -1 if neg, 1 if pos
        Pdotdhat = fabs(Pdotdhat);
        for (int i = 0; i < N_bath; ++i)
            Pperp[i] = PP[i] - signPdotdhat * Pdotdhat * dhat[i];
        alpha *= 2.0;
        sina = sin(alpha);
        cosa = cos(alpha);

        /*!< Importance Sampling - non-adiabatic coupling matrix gives probabilities */
        ap0 = fabs(p0 = ((www[1][SS0][0](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][0](cosa, sina, de, Pdotdhat)));
        ap1 = fabs(p1 = ((www[1][SS0][1](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][1](cosa, sina, de, Pdotdhat)));
        ap2 = fabs(p2 = ((www[1][SS0][2](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][2](cosa, sina, de, Pdotdhat)));
        ap3 = fabs(p3 = ((www[1][SS0][3](cosa, sina, de, Pdotdhat) < -7775.0) ? 0.0 : www[0][SS0][3](cosa, sina, de, Pdotdhat)));
        dn2 = ap0 + ap1 + ap2 + ap3;
        xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements
        //alpha goes to 0, pdotdhat very small, matrix becomes identiy and prob of jumping goes to 0
        oldz = z;
        SS2 = SS1;
        /*!< choose new surface based on probabilities from matrix above
         * calculates new values for the probability weighting
         */
        if (xx < ap0) {
            SS1 = 0;
            z *= p0 * dn2 / ap0;
//            cout << SS1 << endl;
        } else if (xx < ap0 + ap1) {
            SS1 = 1;
            z *= p1 * dn2 / ap1;
//            cout << SS1 << endl;
        } else if (xx < ap0 + ap1 + ap2) {
            SS1 = 2;
            z *= p2 * dn2 / ap2;
//            cout << SS1 << endl;
        } else {
            SS1 = 3;
            z *= p3 * dn2 / ap3;
//            cout << SS1 << endl;
        }

        /*! increases jump counter if a jump was undergone,
         * and exiting if jump counter too high (past truncation value)
         */
        if (SS0 != SS1)
            NNjmp++;
        if (NNjmp > Ncut)
            return 0;

        /*! updating momentum values */
        if (www[1][SS0][SS1](cosa, sina, de, Pdotdhat) != 9999.0){
            for (int i = 0; i < N_bath; ++i) {
                PP[i] = Pperp[i] + signPdotdhat * www[1][SS0][SS1](cosa, sina, de, Pdotdhat) * dhat[i];
            }
        }

        ///////////////////////////////////////////////////////////////////////////////
        /// ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        phase0 = U(RR,PP,SS1,TSLICE*0.5, f, N_bath, timestep, ddd, ddd4, mww, c, force); // exp(iLd/2) (after jump)
        z *= exp(I*phase0);

        ///////////////////////////////////////////////////////////////////////////////
        /// WRITING SUM VALUES OUT (Solving Integral for Observable and Initial Density)
        ///////////////////////////////////////////////////////////////////////////////
        phi = obs[SS1]; /*!< Observable 1 function */

        abszsum1[l] += real(z * phi(RR, PP, N_bath, delta, ddd, c) * initd);
        argzsum1[l] += imag(z * phi(RR, PP, N_bath, delta, ddd, c) * initd);
        realsum[l][NNjmp] += real(z * phi(RR, PP, N_bath, delta, ddd, c) * initd);
        imagsum[l][NNjmp] += imag(z * phi(RR, PP, N_bath, delta, ddd, c) * initd);

        // Unaveraged output to determine variance
        full_abszsum1[n][l] = real(z * phi(RR, PP, N_bath, delta, ddd, c) * initd);

        phi1 = obs1[SS1]; /*!< Observable 2 (Hamiltonian) function */
        habszsum1[l] += real(z * phi1(RR, PP, N_bath, ddd, mww, c) * initd);
        hargzsum1[l] += imag(z * phi1(RR, PP, N_bath, ddd, mww, c) * initd);
        hrealsum[l][NNjmp] += real(z * phi1(RR, PP, N_bath, ddd, mww, c) * initd);
        himagsum[l][NNjmp] += imag(z * phi1(RR, PP, N_bath, ddd, mww, c) * initd);

    }

    delete [] RR; delete [] PP; delete [] Pperp; delete [] dhat;
    delete [] dgam;

    return 0;
}
