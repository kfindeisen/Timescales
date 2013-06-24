/** Frequency generators for Scargle functions
 * @file timescales/freqgen.cpp
 * @author Krzysztof Findeisen
 * @date Created February 16, 2011
 * @date Last modified June 20, 2013
 */ 

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "timescales.h"

namespace kpftimes {

using boost::lexical_cast;
using std::string;

/** Creates a frequency grid that can be fed to time series analysis 
 * functions. The grid itself is trivial to compute; this function therefore 
 * exists mainly as a convenient wrapper for the most commonly needed grids.
 * 
 * @param[in] times	Times at which data were taken
 * @param[out] freq	The new frequency grid
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * 
 * @post @p freq is a grid of frequencies from 0 to @ref pseudoNyquistFreq() 
 *	"pseudoNyquistFreq(@p times)", spaced in units of 1/2T, where T is 
 *	the time interval covered.
 * @post @p freq is sorted in ascending order
 * 
 * @exception std::invalid_argument Thrown if preconditions violated.
 * 
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 */
void freqGen(const DoubleVec &times, DoubleVec &freq) {
	// Delegate input validation to freqGen and pseudoNyquistFreq
	freqGen(times, freq, 0.0, pseudoNyquistFreq(times), 0.5);
}

/** Creates a frequency grid that can be fed to time series analysis 
 * functions. The grid itself is trivial to compute; this function therefore 
 * exists mainly as a convenient wrapper for the most commonly needed grids.
 * 
 * @param[in] times	Times at which data were taken
 * @param[out] freq	The returned frequency grid
 * @param[in] fMin, fMax The requested frequency range\
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre 0 &le; @p fMin < @p fMax
 * 
 * @post @p freq is a grid of frequencies from @p fMin, inclusive, 
 *	to @p fMax, exclusive, spaced in units of 1/2T, where T is 
 *	the time interval covered.
 * @post @p freq is sorted in ascending order
 * 
 * @exception std::invalid_argument Thrown if preconditions violated.
 * 
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 */
void freqGen(const DoubleVec &times, DoubleVec &freq, 
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
 * @param[in] fMin, fMax The requested frequency range
 * @param[in] fStep	The interval at which @p freq will sample the 
 *			frequency range in multiples of 1/(&Delta;t), where 
 *			&Delta;t = max_element(@p times) - min_element(@p times)
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre 0 &le; @p fMin < @p fMax
 * @pre 0 < @p fStep
 * 
 * @post @p freq is a grid of frequencies from @p fMin, inclusive, 
 *	to @p fMax, exclusive, spaced in units of <tt>fstep</tt>/2T, 
 *	where T is the time interval covered.
 * @post @p freq is sorted in ascending order
 * 
 * @exception std::invalid_argument Thrown if preconditions violated.
 * 
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 * 
 * @todo How to test this? At all?
 */
void freqGen(const DoubleVec &times, DoubleVec &freq, 
		double fMin, double fMax, double fStep) {
	if (fMin >= fMax) {
		try {
			throw std::invalid_argument("Parameter 'fMin' should be less than parameter 'fMax' in freqGen() (gave " 
			+ lexical_cast<string>(fMin) + " for fMin and " 
			+ lexical_cast<string>(fMax) + " for fMax)");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameter 'fMin' should be less than parameter 'fMax' in freqGen()");
		}
	}
	if (fMin < 0) {
		try {
			throw std::invalid_argument("Parameter 'fMin' should be nonnegative in freqGen() (gave " 
			+ lexical_cast<string>(fMin) + ")");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameter 'fMin' should be nonnegative in freqGen()");
		}
	}
	if (fStep <= 0) {
		try {
			throw std::invalid_argument("Parameter 'fStep' should be positive in freqGen() (gave " 
			+ lexical_cast<string>(fStep) + ")");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameter 'fStep' should be positive in freqGen()");
		}
	}
	// Delegate remaining input validation to deltaT
	double freqUnit = 1.0/deltaT(times);

	DoubleVec tempFreq;
	   
	for(double curFreq = fMin; curFreq < fMax; curFreq += fStep*freqUnit) {
		tempFreq.push_back(curFreq);
	}
	
	// IMPORTANT: no exceptions beyond this point
	using std::swap;
	swap(freq, tempFreq);
}

}		// end kpftimes
