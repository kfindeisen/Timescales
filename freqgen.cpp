/** Frequency generators for Scargle functions
 * @file freqgen.cpp
 * @author Krzysztof Findeisen
 * @date Created February 16, 2011
 * @date Last modified April 13, 2011
 */ 

#include <vector>
#include "timescales.h"

using namespace kpftimes;

/** Creates a frequency grid that can be fed to time series analysis 
 * functions. The grid itself is trivial to compute; this function therefore 
 * exists mainly as a convenient wrapper for the most commonly needed grids.
 * 
 * @param[in] times	Times at which data were taken
 * @param[out] freq	The returned frequency grid
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @post freq is a grid of frequencies from 0 to @link pseudoNyquistFreq() 
 *	pseudoNyquistFreq(times) @endlink, spaced in units of 1/2T, where T is 
 *	the time interval covered.
 * @post freq is sorted in ascending order
 * @exception invalid_argument Thrown if preconditions violated.
 */
void kpftimes::freqGen(const DoubleVec &times, DoubleVec &freq) {
	// Delegate input validation to freqGen and pseudoNyquistFreq
	freqGen(times, freq, 0.0, pseudoNyquistFreq(times), 0.5);
}

/** Creates a frequency grid that can be fed to time series analysis 
 * functions. The grid itself is trivial to compute; this function therefore 
 * exists mainly as a convenient wrapper for the most commonly needed grids.
 * 
 * @param[in] times	Times at which data were taken
 * @param[out] freq	The returned frequency grid
 * @param[in] fMin	The requested frequency range
 * @param[in] fMax	The requested frequency range
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre 0 <= fMin < fMax
 * @post freq is a grid of frequencies from fMin to fMax, spaced in units of 
 *	1/2T, where T is the time interval covered.
 * @post freq is sorted in ascending order
 * @exception invalid_argument Thrown if preconditions violated.
 */
void kpftimes::freqGen(const DoubleVec &times, DoubleVec &freq, 
		double fMin, double fMax) {
	// Delegate input validation to freqGen
	freqGen(times, freq, fMin, fMax, 0.5);
}

/** Creates a frequency grid that can be fed to time series analysis 
 * functions. The grid itself is trivial to compute; this function therefore 
 * exists mainly as a convenient wrapper for the most commonly needed grids.
 * 
 * @param[in] times	Times at which data were taken
 * @param[out] freq	The returned frequency grid
 * @param[in] fMin	The requested frequency range
 * @param[in] fMax	The requested frequency range
 * @param[in] fStep	The interval at which freq will sample the frequency 
 *				range in multiples of 1/(Delta t).
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre fMin < fMax
 * @pre 0 < fStep
 * @post freq is a grid of frequencies from fMin to fMax, spaced in units of 
 *	fStep*1/T, where T is the time interval covered.
 * @post freq is sorted in ascending order
 * @exception invalid_argument Thrown if preconditions violated.
 * 
 * @todo How to test this? At all?
 */
void kpftimes::freqGen(const DoubleVec &times, DoubleVec &freq, 
		double fMin, double fMax, double fStep) {
	if (fMin >= fMax) {
		throw std::invalid_argument("fMin should be less than fMax");
	} else if (fStep <= 0) {
		throw std::invalid_argument("negative frequency steps are not allowed");
	}
	// Delegate remaining input validation to deltaT
	double freqUnit = 1.0/deltaT(times);
   
	freq.clear();
	for(double curFreq = fMin; curFreq < fMax; curFreq += fStep*freqUnit) {
		freq.push_back(curFreq);
	}
	// fMax may not be a multiple of freqUnit, but we should have it anyway
	freq.push_back(fMax);
}
