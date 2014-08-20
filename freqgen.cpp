/** Frequency generators for Scargle functions
 * @file timescales/freqgen.cpp
 * @author Krzysztof Findeisen
 * @date Created February 16, 2011
 * @date Last modified November 19, 2013
 */ 

/* Copyright 2014, California Institute of Technology.
 *
 * This file is part of the Timescales library.
 * 
 * The Timescales library is free software: you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version, subject to the following 
 * exception added under Section 7 of the License:
 *	* Neither the name of the copyright holder nor the names of its contributors 
 *	  may be used to endorse or promote products derived from this software 
 *	  without specific prior written permission.
 * 
 * The Timescales library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with the Timescales library. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "timeexcept.h"
#include "timescales.h"

namespace kpftimes {

using boost::lexical_cast;
using std::string;

/** @copybrief kpftimes::freqGen(const DoubleVec&,DoubleVec&,double,double,double)
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
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * 
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 */
void freqGen(const DoubleVec &times, DoubleVec &freq) {
	// Delegate input validation to freqGen and pseudoNyquistFreq
	freqGen(times, freq, 0.0, pseudoNyquistFreq(times), 0.5);
}

/** @copybrief kpftimes::freqGen(const DoubleVec&,DoubleVec&,double,double,double)
 * 
 * @param[in] times	Times at which data were taken
 * @param[out] freq	The returned frequency grid
 * @param[in] fMin, fMax The requested frequency range
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
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeRange Thrown if @p fMin &ge; @p fMax
 * @exception std::invalid_argument Thrown if @p fMin or @p fMax is negative
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
 * functions.
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
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeRange Thrown if @p fMin &ge; @p fMax
 * @exception std::invalid_argument Thrown if @p fMin or @p fMax is negative 
 *	or if @p fStep is non-positive
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
			throw except::NegativeRange("Parameter 'fMin' should be less than parameter 'fMax' in freqGen() (gave " 
			+ lexical_cast<string>(fMin) + " for fMin and " 
			+ lexical_cast<string>(fMax) + " for fMax)");
		} catch (const boost::bad_lexical_cast& e) {
			throw except::NegativeRange("Parameter 'fMin' should be less than parameter 'fMax' in freqGen()");
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
