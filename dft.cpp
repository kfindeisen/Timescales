/** Implements the irregularly-sampled discrete Fourier transform
 * @file timescales/dft.cpp
 * @author Krzysztof Findeisen
 * @date Created February 13, 2011
 * @date Last modified November 18, 2013
 */ 

#include <algorithm>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include "dft.h"
#include "../common/stats_except.h"
#include "timeexcept.h"

namespace kpftimes {

using std::string;
using boost::lexical_cast;

/** Calculates the discrete Fourier transform for a list of times and fluxes
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] fluxes	Flux measurements of a source
 * @param[in] freqs	The frequency grid over which the DFT should 
 *			be calculated. See freqGen() for a quick way to 
 *			generate a grid.
 * @param[out] dft	Fourier transform at each frequency.
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre @p fluxes.size() = @p times.size()
 * @pre @p fluxes[i] is the flux of the source at @p times[i], for all i
 * @pre all elements of @p freqs[i] &gt; 0 for all i
 * 
 * @post @p dft.size() = @p freqs.size()
 * @post @p dft[i] is the discrete Fourier transform evaluated at @p freqs[i], for all i
 *
 * @perform O(NF) time, where N = @p times.size() and F = freqs.size()
 *
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some elements of @p freqs are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p times and @p fluxes have 
 *	different lengths.
 * 
 * @exceptsafe The function arguments are unchanged in the event of an exception.
 *
 * @todo Find a faster implementation.
 * @todo Verify that input validation is worth the cost
 */
void dft(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &freqs, ComplexVec &dft) {
	/* This is a brute-force implementation, with no attempt at efficiency
	 * This will later become the reference implementation when I try to 
	 *	replace this with something more subtle
	 */

	using boost::math::double_constants::two_pi;
	const static std::complex<double> I(0.0, 1.0);
	
	size_t nTimes = times.size();
	size_t nFreqs = freqs.size();
	
	// test for non-uniqueness and sorting
	bool diffValues = false, sortedTimes = true;
	for(size_t i = 0; i < nTimes && (!diffValues || sortedTimes); i++) {
		if (!diffValues && times[i] != times[0]) {
			diffValues = true;
		}
		if (sortedTimes && i > 0 && times[i-1] > times[i]) {
			sortedTimes = false;
		}
	}

	// Verify the preconditions
	if (!diffValues) {
		throw except::BadLightCurve("Argument 'times' to dft() contains only one unique date");
	} else if (!sortedTimes) {
		throw kpfutils::except::NotSorted("Argument 'times' to dft() is not sorted in ascending order");
	} else if (fluxes.size() != nTimes) {
		try {
			throw std::invalid_argument("Arguments 'times' and 'fluxes' to dft() are not the same length (gave " 
			+ lexical_cast<string>(nTimes) + " for times and " 
			+ lexical_cast<string>(fluxes.size()) + " for fluxes)");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Arguments 'times' and 'fluxes' to dft() are not the same length");
		}
	}

	// copy-and-swap
	ComplexVec tempDft(nFreqs, 0.0);

	for(size_t i = 0; i < nFreqs; i++) {
		double omega = two_pi * freqs[i];
		for(size_t j = 0; j < nTimes; j++) {
			dft[i] += fluxes[j] * exp(-I * omega * times[j]);
		}
	}
	
	// IMPORTANT: no exceptions beyond this point
	
	using std::swap;
	swap(dft, tempDft);
}

}		// end kpftimes
