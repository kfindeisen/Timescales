/** Implements the irregularly-sampled discrete Fourier transform
 * @file dft.cpp
 * @author Krzysztof Findeisen
 * @date Created February 13, 2011
 * @date Last modified April 14, 2011
 */ 

#include <complex>
#include <stdexcept>
#include <vector>
#include "dft.h"

/* Implementation note: this code extensively uses std::vectors in its 
 *	function arguments to avoid complex memory allocation problems for the 
 *	client and in the function bodies to benefit from the extra 
 *	functionality of array objects compared to C-style arrays
 * The absolute fastest way to fill a vector of known size is to declare it 
 *	without initializing the array, reserve() an array of the right size, 
 *	then use push_back() to fill the array in order. However, this 
 *	approach is vulnerable to counting errors, since the mechanism for 
 *	filling the array never refers to the known array size.
 * A slightly slower but much more maintainable approach is to initialize the 
 *	array to a fixed size, then use the [] operator to set each array 
 *	element correctly. This is the approach taken here.
 */


/** Calculates the discrete Fourier transform for a list of times and fluxes
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] fluxes	Flux measurements of a source
 * @param[in] freqs	The frequency grid over which the DFT should 
 *			be calculated. See freqGen() for a quick way to 
 *			generate a grid.
 * @param[out] dft	Fourier transform at each frequency.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre fluxes is of the same length as times
 * @pre fluxes[i] is the flux of the source at times[i], for all i
 * @pre all elements of freq are >= 0
 * @post dft is of the same length as freq
 * @post dft[i] is the discrete Fourier transform evaluated at freq[i], for all i
 * @exception domain_error Thrown if negative frequencies are provided
 * @exception invalid_argument Thrown if any of the preconditions on the 
 *	format of times or fluxes are violated.
 *
 * @perform O(times.size() × freqs.size()) time
 *
 * @todo Find a faster implementation of dft.
 * @todo Verify that input validation is worth the cost
 */
void kpftimes::dft(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &freqs, ComplexVec &dft) {
	/* This is a brute-force implementation, with no attempt at efficiency
	 * This will later become the reference implementation when I try to 
	 *	replace this with something more subtle
	 */
	const static std::complex<double> I(0, 1);
	const static double PI = 3.1415927;
	
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
		throw std::invalid_argument("times contains only one unique date");
	} else if (!sortedTimes) {
		throw std::invalid_argument("times is not sorted in ascending order");
	} else if (fluxes.size() != nTimes) {
		throw std::invalid_argument("times and fluxes are not the same length");
	}

	// Start working
	dft.clear();
	dft.resize(nFreqs);

	for(size_t i = 0; i < nFreqs; i++) {
		double omega = 2.0 * PI * freqs[i];
		dft[i] = 0.0;
		for(size_t j = 0; j < nTimes; j++) {
			dft[i] += fluxes[j] * exp(-I * omega * times[j]);
		}
	}
}
