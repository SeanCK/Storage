/*!
 *
 *  types.h
 *  Created on: Mar 22, 2018
 *
 * \brief
 * This file contains datatypes declarations for internal
   use in the code
 *
 */
#ifndef TYPES_H_
#define TYPES_H_

#include <chrono>
#include <complex>
#include <limits>
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>

// Primitive types
typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;
typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;
typedef float float32;
typedef double float64;

// Derived types?
typedef float64 real;
typedef std::complex<float64> dcomplex;
typedef std::vector<real> rvec;						// Change vectors and matrix to pointer here
typedef std::vector< std::vector<real> > matrix;    // If performance starts to become an issue

// Chrono typedef's for convenience
typedef std::chrono::nanoseconds  NANOSECONDS;
typedef std::chrono::microseconds MICROSECONDS;
typedef std::chrono::milliseconds MILLISECONDS;
typedef std::chrono::seconds      SECONDS;
typedef std::chrono::high_resolution_clock hires_clock;

/*!
 *  \brief Gets current time since some epoch (probably the unix one?)
 *         using std::chrono
 */
class CurrentTime {

	static hires_clock h_clock;

public:
	static uint64 seconds() {
		return std::chrono::duration_cast<SECONDS>
				(h_clock.now().time_since_epoch()).count();
	}

	static uint64 milliseconds(){
		return std::chrono::duration_cast<MILLISECONDS>
				(h_clock.now().time_since_epoch()).count();
	}
	static uint64 microseconds(){
		return std::chrono::duration_cast<MICROSECONDS>
				(h_clock.now().time_since_epoch()).count();
	}
	static uint64 nanoseconds(){
		return std::chrono::duration_cast<NANOSECONDS>
				(h_clock.now().time_since_epoch()).count();
	}

};

// Math constants
#define PI 3.14159265358979323846
#define E  2.71828182845904523536
//dcomplex I(0,1);


typedef std::unordered_map<std::string, std::string> OptionsMap;
typedef std::unordered_set<std::string> OptionsSet;

/*! Input options read from input file */
typedef struct {
	uint32 N_bath;    /*!< Size of bath */
	uint32 N_tslice;  /*!< Number of time intervals */
	uint32 N_sample;  /*!< Sample Size (No. of trees calculated) */
	uint32 N_cut;     /*!< Truncation parameter */
	real total_time;  /*!< Total time */     //Consider using just dt = total_time/N_tslice?
	real beta;        /*!< Inverse Temperature */
	real delta;       /*!< Related to energy difference between ground and excited state in QM subsystem \f$ (\delta) \f$*/
	real w_max;       /*!< Omega_max a cut-off frequency characterizing the spectral density \f$ (\omega) \f$ */
	real eta;         /*!< The Kondo parameter characterizing the spectral density \f$ (\eta) \f$ */
	real prop_tstep;   /*!< The time step for adiabatic propagator */
	uint8 Regression_test;

} input_params_t;

/*!
 *  \brief Observable from the evolution of the system
 *         stores real and imaginary sums and partial sums
 */
class Observable {
	uint32 DIM;
public:
	rvec rsum, isum;
	matrix part_rsum, part_isum;

	/*!
	 * \brief Sets up and Initializes observable matrices/vectors
	 * 		  of provided dimension.
	 * \param DIM number of observable points (dt)
	 */
	Observable(uint32 _DIM) : DIM(_DIM) {
		rsum.reserve(DIM);
		isum.reserve(DIM);
		for (rvec::size_type i=0; i<DIM; ++i) {
			rsum.push_back(0);
			isum.push_back(0);
		};

		part_rsum.resize(DIM,rvec(DIM));
		part_isum.resize(DIM,rvec(DIM));
		for (rvec::size_type i=0; i<DIM; ++i) {
			for(rvec::size_type j=0; j<DIM+1; ++j) {
				part_rsum[i][j] = 0;
				part_isum[i][j] = 0;
			};
		};
	}
	~Observable(){
		//no need to do anything since memory from std::vector will be automatically freed once
		//the observable object is destroyed.
		//If you change back to raw pointer based arrays then call delete[] here
	}
};

typedef struct {
	rvec mass;		/*!< Mass of bath particles? */
	rvec c;			/*!< A spectral density parameter?*/
	rvec w;			/*!< A spectral density parameter? */
	rvec sigma;		/*!< The std. deviation sqrt(variance)*/
	rvec mu;		/*!< Convenience parameter in dimensionless units */

} bath_params_t;
#endif /* TYPES_H_ */
