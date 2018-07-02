/*!
 *  \mainpage
 *  \brief Using a Trotter Approximation to calculate a quantum-classical non-adiabatic approximation to
 * separate the quantum and the quasi-classical degrees of freedom to allow a surface-hopping scheme to
 * be implemented that creates a tree-like structure.
 *
 * This code follows one of the possible paths through the tree structure. The non-adiabatic propagator
 * determines the likelihood of a jump based on importance sampling
 *
 * \author Donal MacKernan, Athina Lange, Philip McGrath, Sean Kelly and Shrinath Kumar
 * \date 29.04.18
 */

//#include "opt_parser.h"
#include <fstream>
#include <algorithm>

#include "density.h"
#include "dmatrix.c"
#include <omp.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/// SYSTEM INPUT
///////////////////////////////////////////////////////////////////////////////

// ================================================================
// MAIN
// ================================================================


int main(int argc, char *argv[]){

    char  datafilename[40] = {"Output"};
    FILE   *stream;

    const gsl_rng_type * TT; /*!< Random Number Generator seed based on Gaussian Distribution */

    double ddd; /*!<  \f$ \delta^2 \f$ */
    double ddd4; /*!< \f$ \delta^2/4 \f$ */
    double *m; /*!< Mass of bath particles */
    double *c; /*!< Coupling coefficients between quantum and classical systems */
    double *w; /*!< Frequencies selected to approximate the ohmic spectral density */

    double *mww; /*!< Variable defined for convenience: \f$ mass \times frequency^2 \f$ */
    double *sig; /*!< Standard deviation of gaussian distribution */
    double *mu; /*!< */
    double TSLICE; /*!< Time per time interval*/

    double *abszsum1;
    double *argzsum1;
    double *habszsum1;
    double *hargzsum1;

    double **full_abszsum1;
    double *variance;

    double **realsum;
    double **imagsum;
    double **hrealsum;
    double **himagsum;

    double (*dens_init[4])(double*, double*, int, double, double, double *); /*!< Initial Density Matrix*/
    double (*obs[4])(double*, double*, int, double, double, double *); /*!< Observable Matrix (population density)*/
    double (*obs1[4])(double*, double*, int, double, double*, double *); /*!< Secondary Observable Matrix (hamiltonian)*/
    void (*force[4])(double *, double*, int, double, double*, double*); /*!< Hellman-Feynman Forces*/
    double (*www[2][4][4])(double , double , double , double); /*!< Transition Matrix*/

    ///////////////////////////////////////////////////////////////////////////////
    /// SYSTEM INPUT
    ///////////////////////////////////////////////////////////////////////////////

    // Read input file
//    string input_path(argv[1]);
//    input_params_t input = read_opt_file(input_path);

//    int N_bath = input.N_bath;
//    int N_slice = input.N_tslice;
//    int Ncut = input.N_cut;
//    double timestep = input.prop_tstep;
//    double T = input.total_time;
//    int Nsample = input.N_sample;
//    double w_max = input.w_max;
//    double eta = input.eta;
//    double beta = input.beta;
//    double delta = input.delta;

    int N_bath = 200;
    int N_slice = 60;
    int Ncut = 10;
    double timestep = 0.05;
    double T = 15;
    int Nsample = 10000;
    double w_max = 3;
    double eta = 0.13;
    double beta = 25;
    double delta = 0.8;

    int Regression_test = 1;


    ///////////////////////////////////////////////////////////////////////////////
    /// ALLOCATING MEMORY
    ///////////////////////////////////////////////////////////////////////////////

    mww = new double[N_bath];
    mu = new double[N_bath];
    sig =  new double[2*N_bath];
    c = new double[N_bath];
    m = new double[N_bath];
    w = new double[N_bath];
    abszsum1  = new double[N_slice];
    argzsum1  = new double[N_slice];
    habszsum1  = new double[N_slice];
    hargzsum1  = new double[N_slice];
    variance = new double[N_slice];
    full_abszsum1 = dmatrix(0,Nsample,0,N_slice);
    realsum = dmatrix(0,N_slice,0,N_slice+1);   // Included the real/imag arrays to compare output to OldCode output
    imagsum = dmatrix(0,N_slice,0,N_slice+1);
    hrealsum = dmatrix(0,N_slice,0,N_slice+1);
    himagsum = dmatrix(0,N_slice,0,N_slice+1);

    ///////////////////////////////////////////////////////////////////////////////
    /// INITIALIZATION OF SYSTEM
    ///////////////////////////////////////////////////////////////////////////////

    dens_init[0] = dens_init_0; dens_init[1] = dens_init_1;
    dens_init[2] = dens_init_2; dens_init[3] = dens_init_3;
    obs[0] = obs_0; obs[1] = obs_1; obs[2] = obs_2; obs[3] = obs_3;
    obs1[0] = H_0; obs1[1] = H_1; obs1[2] = H_2; obs1[3] = H_3;

    ddd4 = delta*delta*0.25;
    ddd =  delta*delta;
    TSLICE  = T/N_slice;

    bath_para(eta, w_max, N_bath, c, m, w);       /*!< Defining Bath Parameters */

    for (int i = 0; i < N_bath; ++i){
        mu[i] = beta * w[i] * 0.5;
    }
    for (int i = 0; i < N_bath; ++i){
        sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
        mww[i] = -m[i]*w[i]*w[i];
    }
    for (int i = 0; i < N_bath; ++i) {
        sig[i + N_bath] = 1.0 * sqrt(w[i] / (2.0 * tanh(mu[i])));
    }

    /*! Defining force field */
    force[0] = F1;
    force[1] = Fb;
    force[2] = Fb;
    force[3] = F2;

    // Create output file
    stream = fopen(datafilename,"w");
    fprintf(stream,"%s\n w_max %lf eta %lf beta %lf delta %lf N_bath %d N_slice %d Nsample %d\n", argv[0], w_max, eta, beta, delta, N_bath, N_slice, Nsample);
    fclose(stream);

    /*! Initializing sum1 counters*/
    for (int i = 0; i < N_slice; ++i){
        abszsum1[i] = 0.0;
        argzsum1[i]  = 0.0;
        habszsum1[i] = 0.0;
        hargzsum1[i] = 0.0;

        variance[i] = 0.0;

        for (int j = 0; j <= N_slice; j++){
            realsum[i][j] = 0.0;
            imagsum[i][j] = 0.0;
            hrealsum[i][j] = 0.0;
            himagsum[i][j] = 0.0;
            }
    }

    for (int i1 = 0; i1 < Nsample; ++i1) {
        for (int j = 0; j < N_slice; ++j) {
            full_abszsum1[i1][j] = 0.0;
        }

    }

    /*! Defining non-adiabatic coupling matrix */
    setwww(www);
    ///////////////////////////////////////////////////////////////////////////////
    /// PROCESSING TREE
    ///////////////////////////////////////////////////////////////////////////////

    // Print column headers for output file
    stream = fopen(datafilename,"a");
    fprintf(stream,"Nens\tDt\tj\tabszsum1\targzsum1\trealsum\timagsum\thabszsum1\thargzsum1\threalsum\thimagsum\n");
    fclose(stream);


    int numthreads = omp_get_max_threads(); // Reads the number of threads from the environment variable set by 'export OMP_NUM_THREADS=...'

    gsl_rng ** rrp;
    gsl_rng_env_setup();
    TT = gsl_rng_default;
    rrp = new gsl_rng * [numthreads];
    int * seed_rrp;
    seed_rrp = new int [numthreads];

    /*
    // The time in nanoseconds is used to seed the random number generator
    CurrentTime Time;

    for (int i = 0; i < numthreads; ++i) {
        rrp[i] = gsl_rng_alloc(TT);
        seed_rrp[i] = Time.nanoseconds();
        gsl_rng_set(rrp[i], seed_rrp[i]);
    }
    */

    // The thread id is used to seed the random number generator
    for (int i = 0; i < numthreads; ++i) {
        rrp[i] = gsl_rng_alloc(TT);
        seed_rrp[i] = i;
        gsl_rng_set(rrp[i], seed_rrp[i]);
    }

#pragma omp parallel for num_threads(numthreads) reduction(+:abszsum1[0:N_slice], argzsum1[0:N_slice], habszsum1[0:N_slice], hargzsum1[0:N_slice])
for (int i = 0; i < Nsample; ++i){

	double *f; /*!< Force on particles */
	double  *R1, *v;

    double ranVector[4 * N_bath];

    int id;
    id = omp_get_thread_num();

    /*! Sets up the use of a Gaussian Random Number Generator from GSL */

    gsl_rng * rlocal;
    rlocal = rrp[id];

    R1 = new double[N_bath];
	v = new double[N_bath];
	f = new double[N_bath];

//        if(i % 100 == 0) {
            printf(" calling density %d %d j %d \n", id, numthreads, i);
//        }
	density(R1, v, f, abszsum1, argzsum1, habszsum1, hargzsum1, realsum, imagsum, hrealsum, himagsum, rlocal, ranVector, N_bath, Ncut, m, N_slice, TSLICE, timestep,
            delta, ddd, ddd4, mww, c, sig, dens_init, obs, obs1, force, www, full_abszsum1, i);

	if (i == 1){numthreads = omp_get_num_threads();} // Detects number of threads in parallel region

	delete [] R1; delete [] v; delete [] f;

}

    // Calculating variance for each timestep
    for (int n = 0; n < Nsample; ++n) {
        for (int k = 0; k < N_slice; ++k) {
            variance[k] += pow((full_abszsum1[n][k] - (abszsum1[k]/Nsample)), 2) / Nsample;
        }
    }

    stream = fopen(datafilename, "a");
    fprintf(stream, "dt %lf T %lf Nsample %d\n", timestep, T, Nsample);
    printf("Time        abszsum1   abs_error   argzsum1   realsum    habszsum1   hargzsum1   hrealsum\n");
    fprintf(stream, "Time   abszsum1   abs_error   argzsum1   realsum   habszsum1   hargzsum1   hrealsum\n");
    for (int i = 0; i < N_slice; ++i) {
//        for (int j = 0; j < Ncut; ++j) {
        int j = 0;
        printf("%lf    %lf   %lf   %lf   %lf   %lf   %lf   %lf\n", TSLICE * (i + 1),
               (abszsum1[i] / Nsample), pow(variance[i]/Nsample, 0.5), (argzsum1[i] / Nsample), realsum[i][j] / Nsample,
               (habszsum1[i] / Nsample), (hargzsum1[i] / Nsample), hrealsum[i][j] / Nsample);
        fprintf(stream, "%.2lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf\n", TSLICE * (i + 1),
                (abszsum1[i] / Nsample), pow(variance[i]/Nsample, 0.5), (argzsum1[i] / Nsample), realsum[i][j] / Nsample,
                (habszsum1[i] / Nsample), (hargzsum1[i] / Nsample), hrealsum[i][j] / Nsample);
//        }
    }
    fprintf(stream, "Number of threads was %d", numthreads);
    fclose(stream);

    cout << "Number of threads was "<< numthreads << endl;


    // Regression Testing
//    if(input.Regression_test == 1){
    if(Regression_test == 1){

    // Reading expected and variance values from file
    double *expected_value = new double[62];
    double *std_dev_expected = new double[62];
    ifstream exp_file, var_file;
    exp_file.open("Regression_testing/ExpectedValue_abszsum1");
    var_file.open("Regression_testing/Variance_abszsum1");
    string exp_line_string, var_line_string;
    double exp_line_double, var_line_double;
    int k = 0;
    while (!exp_file.eof()) {
        getline(exp_file, exp_line_string);         // reads line from file as string
        stringstream converter1(exp_line_string);
        converter1 >> exp_line_double;               // converts string to double
        expected_value[k] = exp_line_double;

        getline(var_file, var_line_string);
        stringstream converter2(var_line_string);
        converter2 >> var_line_double;
        std_dev_expected[k] = pow(var_line_double, 0.5);

        ++k;
    }

    // Checks values of abszum1 up to T = 3 are within 5 standard deviations of the expected value
    vector<bool> Regression;
    double val;
    double lower_lim;
    double upper_lim;

    cout << endl << "Regression Test" << endl;
    cout <<"Time  lower_lim < abszum1 < upper_lim" << endl;
    for (int l = 0; l < 60; ++l) {
        val = abszsum1[l] / Nsample;
        lower_lim = expected_value[l + 1] - (5 * std_dev_expected[l + 1]);
        upper_lim = expected_value[l + 1] + (5 * std_dev_expected[l + 1]);

        if (val >= lower_lim && val <= upper_lim) { Regression.push_back(1); }
        else { Regression.push_back(0); } // 1 = passed,  0 = failed

        printf("%.2lf  %lf  < %lf < %lf", (l + 1) * 0.25, lower_lim, val, upper_lim);
        if(Regression[l] == 0){cout << "  X" << endl;}
        else {cout << endl;}
    }

    // If abszsum1 goes out of bounds at any point the test is failed
    if (find(Regression.begin(), Regression.end(), 0)!=Regression.end()){cout << "Regression test Failed" << endl;}
    else {cout << "Regression test Passed" << endl;}
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// DEALLOCATING MEMEORY
    ///////////////////////////////////////////////////////////////////////////////

    delete [] abszsum1; delete [] argzsum1; delete [] habszsum1; delete [] hargzsum1;
    delete [] realsum; delete [] imagsum; delete [] hrealsum; delete [] himagsum;
    delete [] variance; delete [] full_abszsum1;

    delete [] mww; delete [] mu; delete [] sig;
    delete [] c; delete [] m; delete [] w;

    return 0;
}
