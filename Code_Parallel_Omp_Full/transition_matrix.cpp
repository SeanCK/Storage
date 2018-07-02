#include "transition_matrix.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/// FUNCTIONS
///////////////////////////////////////////////////////////////////////////////

/*! Derivative of gam */
double dgamma(double *R, int i, double *c){
    return -c[i];
}

/*! Non-adiabatic coupling matrix */
double  dd(double *dhat, double *R, double *dgam, double ddd, int N_bath, double delta, double *c){
    double abs_d;
    double x1,x2,x3;
    int i;
    x2 = gam(R, N_bath, c);
    for (int i = 0; i < N_bath; ++i)
        dgam[i] = dgamma(R, i, c);
    if (fabs(x2) < 1.0e-4)
        x3 = 1/delta;
    else {
        x1 = G(R, N_bath, delta, ddd, c);
        x3 = -x1/x2 + 2.0/(delta + 2.0*x2*x1);
        x3 = x3/(1.0 + x1*x1);
    }
    for (i = 0,abs_d = 0.0; i < N_bath; ++i){
        dhat[i] = -dgam[i]*x3;
        abs_d += dhat[i]*dhat[i];
    }
    abs_d = sqrt(abs_d);
    for (i = 0; i < N_bath; ++i)
        dhat[i] /= abs_d;
    return abs_d;
}

///////////////////////////////////////////////////////////////////////////////
/// TRANSITION MATRIX
///////////////////////////////////////////////////////////////////////////////

/* Second order Trotter decomposition of the partially Wigner transformed spin boson hamiltonian */

double wwa0_00(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = 1 + cosa;
    return x*0.5;
}

double wwa0_01(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = -sina;
    return x*0.5;
}

double wwa0_02(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = -sina;
    return x*0.5;
}

double wwa0_03(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = 1.0 - cosa;
    return x*0.5;
}

double wwa0_10(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = sina;
    return x*0.5;
}

double wwa0_11(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = 1.0 + cosa;
    return x*0.5;
}

double wwa0_12(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = -1.0 + cosa;
    return x*0.5;
}

double wwa0_13(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = -sina;
    return x*0.5;
}

double wwa0_20(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = sina;
    return x*0.5;
}

double wwa0_21(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = -1.0 + cosa ;
    return x*0.5;
}

double wwa0_22(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = 1.0 + cosa;
    return x*0.5;
}

double wwa0_23(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = -sina ;
    return x*0.5;
}

double wwa0_30(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = 1 - cosa;
    return x*0.5;
}

double wwa0_31(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x =  sina;
    return x*0.5;
}

double wwa0_32(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = sina;
    return x*0.5;
}

double wwa0_33(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = 1 + cosa;
    return x*0.5;
}

//         W_{a 1}

/* _____________________________________________  */

double wwa1_00(double cosa, double sina, double de, double Pdotdhat){
    return 9999.0;
}

double wwa1_01(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_02(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_03(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat - 2.0*de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_10(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_11(double cosa, double sina, double de, double Pdotdhat){
    return 9999.0;
}

double wwa1_12(double cosa, double sina, double de, double Pdotdhat){
    return 9999.0;
}

double wwa1_13(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_20(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_21(double cosa, double sina, double de, double Pdotdhat){
    return 9999.0;

}

double wwa1_22(double cosa, double sina, double de, double Pdotdhat){
    return 9999.0;
}

double wwa1_23(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_30(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat + 2.0*de;
    return sqrt(x);
}

double wwa1_31(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_32(double cosa, double sina, double de, double Pdotdhat){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_33(double cosa, double sina, double de, double Pdotdhat){
    return 9999.0;
}

/*! Non-adiabatic Coupling Matrix */
void setwww(double (*www[2][4][4])(double, double, double, double)){

    // W_{a0}

    www[0][0][0] = wwa0_00;
    www[0][0][1] = wwa0_01;
    www[0][0][2] = wwa0_02;
    www[0][0][3] = wwa0_03;
    www[0][1][0] = wwa0_10;
    www[0][1][1] = wwa0_11;
    www[0][1][2] = wwa0_12;
    www[0][1][3] = wwa0_13;
    www[0][2][0] = wwa0_20;
    www[0][2][1] = wwa0_21;
    www[0][2][2] = wwa0_22;
    www[0][2][3] = wwa0_23;
    www[0][3][0] = wwa0_30;
    www[0][3][1] = wwa0_31;
    www[0][3][2] = wwa0_32;
    www[0][3][3] = wwa0_33;

    // W_{a1}

    www[1][0][0] = wwa1_00;
    www[1][0][1] = wwa1_01;
    www[1][0][2] = wwa1_02;
    www[1][0][3] = wwa1_03;
    www[1][1][0] = wwa1_10;
    www[1][1][1] = wwa1_11;
    www[1][1][2] = wwa1_12;
    www[1][1][3] = wwa1_13;
    www[1][2][0] = wwa1_20;
    www[1][2][1] = wwa1_21;
    www[1][2][2] = wwa1_22;
    www[1][2][3] = wwa1_23;
    www[1][3][0] = wwa1_30;
    www[1][3][1] = wwa1_31;
    www[1][3][2] = wwa1_32;
    www[1][3][3] = wwa1_33;

}


/* ______________________________________   */

