/** Autocorrelation function for unevenly sampled data
 * @file autocorr.cpp
 * @author Krzysztof Findeisen
 * @date Created February 16, 2011
 * @date Last modified June 18, 2013
 */ 

#include <complex>
#include <vector>
#include <gsl/gsl_fft_halfcomplex.h>
#include "dft.h"
#include "timescales.h"
#include <cmath>
#include "../common/stats.tmp.h"

/** Calculates the autocorrelation function for a time series. 
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] fluxes	Flux measurements of a source
 * @param[in] offsets	The time grid over which the autocorrelation function 
 *			should be calculated. 
 * @param[out] acf	The value of the autocorrelation function at each 
 *			offset.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre fluxes is of the same length as times
 * @pre fluxes[i] is the flux of the source at times[i], for all i
 * @pre offsets is uniformly sampled from 0 to some maximum offset. This 
 *	requirement will be relaxed in future versions.
 * @post acf is of the same length as offsets
 * @post acf[i] is the Scargle autocorrelation function evaluated at offsets[i], for all i
 * @exception invalid_argument Thrown if the preconditions on times, 
 *	offsets, or length(fluxes) are violated.
 *
 * @perform O(times.size() × offsets.size()) time
 *
 * @todo Allow autoCorr to run with arbitrary offset grids
 * @todo Verify that input validation is worth the cost
 * @bug Many ACFs show a strong fringing effect. The fringes have a wavelength 
 * corresponding to the pseudo-Nyquist frequency, or to the highest peak 
 * frequency that is below the pseudo-Nyquist frequency if the power spectrum 
 * of the data is strongly peaked.
 */
void kpftimes::autoCorr(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &offsets, DoubleVec &acf) {
	autoCorr(times, fluxes, offsets, acf, pseudoNyquistFreq(times));
}

/** Calculates the autocorrelation function for a time series. 
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] fluxes	Flux measurements of a source
 * @param[in] offsets	The time grid over which the autocorrelation function 
 *			should be calculated.
 * @param[in] maxFreq	The maximum frequency to consider when calculating 
 *			the autocorrelation function.
 * @param[out] acf	The value of the autocorrelation function at each 
 *			offset.
 *
 * @note Increasing maxFreq will increase the time resolution of acf at the 
 *	cost of making the entire function noisier.
 * 
 * @warning This version of autoCorr is highly volatile and may be removed 
 *	from the library in future versions. I recommend its use only for 
 *	testing of the autocorrelation algorithm. Wherever possible, use 
 *	autoCorr(), which has a stable (or at least forward-compatible) spec.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre fluxes is of the same length as times
 * @pre fluxes[i] is the flux of the source at times[i], for all i
 * @pre offsets contains at least two elements
 * @pre offsets contains only nonnegative values
 * @pre offsets is uniformly sampled from 0 to some maximum offset. This 
 *	requirement will be relaxed in future versions.
 * @pre maxFreq is positive
 * @post acf is of the same length as offsets
 * @post acf[i] is the Scargle autocorrelation function evaluated at offsets[i], for all i
 * @exception std::invalid_argument Thrown if the preconditions on times, 
 *	offsets, length(fluxes), or maxFreq are violated.
 *
 * @perform O(times.size() × offsets.size()) time
 *
 * @todo Verify that input validation is worth the cost
 */
void kpftimes::autoCorr(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &offsets, DoubleVec &acf, double maxFreq) {
	size_t nTimes  = times.size();
	size_t nOutput = offsets.size();
	double tRange  = deltaT(times);
//	double pseudoNyquist = 0.5*nTimes/tRange;
	
	// Normalize the data
	// Scargle (1989) argues this is too crude, but a fitting method of 
	//	the kind he proposes is too slow and too inflexible
	// While we're at it, construct the flat source we'll use to get the 
	//	window function
	// While we're at it, also test for non-uniqueness and sorting
	DoubleVec zeroFluxes(nTimes), oneFluxes(nTimes);
	double meanFlux = kpfutils::mean(fluxes.begin(), fluxes.end());
	bool diffValues = false, sortedTimes = true;
	for(size_t i = 0; i < nTimes; i++) {
		zeroFluxes[i] = fluxes[i] - meanFlux;
		oneFluxes[i] = 1.0;

		if (!diffValues && times[i] != times.front()) {
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
	} else if (maxFreq <= 0.0) {
		throw std::invalid_argument("maxFreq must be positive");
	}

	// Verify offsets
	if (offsets.size() < 2) {
		throw std::invalid_argument("need at least two elements in offsets for a meaningful ACF");
	} else if (offsets[0] != 0.0) {
		throw std::invalid_argument("first element of offsets must be zero (for now)");
	}
	double offSpace = offsets[1] - offsets[0];
	if (offSpace <= 0.0) {
		throw std::invalid_argument("offsets must be in ascending order");
	}
	for(size_t i = 2; i < nOutput; i++) {
		if (offsets[i] < 0.0) {
			throw std::invalid_argument("offsets must be nonnegative");
		}
		// abs is acting wacky on Cygwin... fix later
		if ((offsets[i] - offsets[i - 1] - offSpace)/offSpace > 1e-3 
				|| (offsets[i] - offsets[i - 1] - offSpace)/offSpace < -1e-3) {
			throw std::invalid_argument("offsets must have uniform spacing (for now)");
		}
	}
	
	// Generate a regular frequency grid
	// A frequency step of 1/2 Delta T is sufficient to prevent wraparound
	// Going to a finer frequency spacing has no effect other than to 
	// 	introduce a string of zeros at offset > Delta T in the final 
	//	answer
	DoubleVec freq;
	double freqStep = 0.5/tRange;
	kpftimes::freqGen(times, freq, 0.0, 0.5/offSpace, 0.5);

	// Forward transform
	ComplexVec xForm, winXForm;
	kpftimes::dft(times, zeroFluxes, freq,    xForm);
	kpftimes::dft(times,  oneFluxes, freq, winXForm);

	// Zero all the high frequencies
	// Note if maxFreq > 0.5/offSpace, this code has no effect
	//	Nor should it, since this means the output grid is too coarse 
	//	to recover the high-frequency information
	// For now, this is a sharp cutoff
	// Introducing some apodizing may be helpful
	for(size_t i= 0; i < freq.size(); i++) {
		if (freq[i] > maxFreq) {
			   xForm[i] = std::complex<double>(0.0, 0.0);
			winXForm[i] = std::complex<double>(0.0, 0.0);
		}
	}
	
	// Get the power spectrum
	for(ComplexVec::iterator it =    xForm.begin(); it !=    xForm.end(); it++) {
		*it = norm(*it);
	}
	for(ComplexVec::iterator it = winXForm.begin(); it != winXForm.end(); it++) {
		*it = norm(*it);
	}
	
	// Reverse transform -- setup
	size_t gslSize = 2*xForm.size() - 1;
	gsl_fft_halfcomplex_wavetable* theTable = 
		gsl_fft_halfcomplex_wavetable_alloc(gslSize);
	gsl_fft_real_workspace * theSpace = 
		gsl_fft_real_workspace_alloc(gslSize);
		
	double *gslBuffer = new double[gslSize];
	gslBuffer[0] = xForm[0].real();
	for (size_t i = 1; i < xForm.size(); i++) {
		gslBuffer[2*i-1] = xForm[i].real();
		gslBuffer[2*i  ] = xForm[i].imag();
	}
	// assert: gsl_size is odd
	// therefore, gslBuffer follows the odd-transform convention, and only 
	//	xForm[0].imag() need be dropped
	
	// Reverse transform -- action
	gsl_fft_halfcomplex_transform(gslBuffer, 1, gslSize, theTable, theSpace);
	
	// Reinterpret gslBuffer as autocorrelation function
	acf.clear();
	for(size_t i = 0; i < gslSize; i++) {
		double time = static_cast<double>(i)/((gslSize-1) * freqStep);
		// Values above tRange are aliases, so don't record them
		if (time <= tRange) {
			//offsets.push_back(time);
			acf    .push_back(gslBuffer[i]/gslBuffer[0]);
		} else {
			acf    .push_back(0.0);
		}
	}
	// Clip or extend to match length of offsets
	acf.resize(offsets.size(), 0.0);
	
	// Now repeat for the sampling ACF
	gslBuffer[0] = winXForm[0].real();
	for (size_t i = 1; i < winXForm.size(); i++) {
		gslBuffer[2*i-1] = winXForm[i].real();
		gslBuffer[2*i  ] = winXForm[i].imag();
	}
	gsl_fft_halfcomplex_transform(gslBuffer, 1, gslSize, theTable, theSpace);
	DoubleVec winAcf;
	for(size_t i = 0; i < gslSize; i++) {
		double time = static_cast<double>(i)/((gslSize-1) * freqStep);
		// Values above tRange are aliases, so don't record them
		if (time <= tRange) {
			winAcf.push_back(gslBuffer[i]/gslBuffer[0]);
		} else {
			winAcf.push_back(1.0);
		}
	}
	// Clip or extend to match length of offsets
	winAcf.resize(offsets.size(), 1.0);
	
	// Clean up
	delete [] gslBuffer;
	gsl_fft_halfcomplex_wavetable_free(theTable);
	gsl_fft_real_workspace_free(theSpace);
	
	// Normalize
	for(size_t i = 0; i < acf.size(); i++) {
		acf[i] = acf[i] / winAcf[i];
	}
}

/** Calculates the autocorrelation window function for a time sampling. 
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] offsets	The time grid over which the autocorrelation function 
 *			should be calculated. 
 * @param[out] wf	The value of the window function at each offset.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre offsets is uniformly sampled from 0 to some maximum offset. This 
 *	requirement will be relaxed in future versions.
 * @post wf is of the same length as offsets
 * @post wf[i] is the Scargle autocorrelation function evaluated at offsets[i], for all i
 * @exception invalid_argument Thrown if the preconditions on times, 
 *	or offsets are violated.
 *
 * @perform O(times.size() × offsets.size()) time
 *
 * @todo Verify that input validation is worth the cost
 */
void kpftimes::acWindow(const DoubleVec &times, const DoubleVec &offsets, DoubleVec &wf) {
	acWindow(times, offsets, wf, pseudoNyquistFreq(times));
}

/** Calculates the autocorrelation window function for a time sampling. 
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] offsets	The time grid over which the autocorrelation function 
 *			should be calculated.
 * @param[in] maxFreq	The maximum frequency to consider when calculating 
 *			the autocorrelation function.
 * @param[out] wf	The value of the window function at each offset.
 *
 * @warning This version of acWindow is highly volatile and may be removed 
 *	from the library in future versions. I recommend its use only for 
 *	testing of the autocorrelation algorithm. Wherever possible, use 
 *	acWindow(), which has a stable (or at least forward-compatible) spec.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre offsets contains at least two elements
 * @pre offsets contains only nonnegative values
 * @pre offsets is uniformly sampled from 0 to some maximum offset. This 
 *	requirement will be relaxed in future versions.
 * @pre maxFreq is positive
 * @post wf is of the same length as offsets
 * @post wf[i] is the Scargle autocorrelation function evaluated at offsets[i], for all i
 * @exception std::invalid_argument Thrown if the preconditions on times, 
 *	offsets, or maxFreq are violated.
 *
 * @perform O(times.size() × offsets.size()) time
 *
 * @todo Verify that input validation is worth the cost
 */
void kpftimes::acWindow(const DoubleVec &times, const DoubleVec &offsets, DoubleVec &wf, 
		double maxFreq) {
	size_t nTimes  = times.size();
	size_t nOutput = offsets.size();
	double tRange  = deltaT(times);
	
	// Construct the flat source we'll use to get the 
	//	window function
	// While we're at it, also test for non-uniqueness and sorting
	DoubleVec oneFluxes(nTimes);
	bool diffValues = false, sortedTimes = true;
	for(size_t i = 0; i < nTimes; i++) {
		oneFluxes[i] = 1.0;

		if (!diffValues && times[i] != times.front()) {
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
	} else if (maxFreq <= 0.0) {
		throw std::invalid_argument("maxFreq must be positive");
	}

	// Verify offsets
	if (offsets.size() < 2) {
		throw std::invalid_argument("need at least two elements in offsets for a meaningful ACF");
	} else if (offsets[0] != 0.0) {
		throw std::invalid_argument("first element of offsets must be zero (for now)");
	}
	double offSpace = offsets[1] - offsets[0];
	if (offSpace <= 0.0) {
		throw std::invalid_argument("offsets must be in ascending order");
	}
	for(size_t i = 2; i < nOutput; i++) {
		if (offsets[i] < 0.0) {
			throw std::invalid_argument("offsets must be nonnegative");
		}
		// abs is acting wacky on Cygwin... fix later
		if ((offsets[i] - offsets[i - 1] - offSpace)/offSpace > 1e-3 
				|| (offsets[i] - offsets[i - 1] - offSpace)/offSpace < -1e-3) {
			throw std::invalid_argument("offsets must have uniform spacing (for now)");
		}
	}
	
	// Generate a regular frequency grid
	// A frequency step of 1/2 Delta T is sufficient to prevent wraparound
	// Going to a finer frequency spacing has no effect other than to 
	// 	introduce a string of zeros at offset > Delta T in the final 
	//	answer
	DoubleVec freq;
	double freqStep = 0.5/tRange;
	kpftimes::freqGen(times, freq, 0.0, 0.5/offSpace, 0.5);

	// Forward transform
	ComplexVec winXForm;
	kpftimes::dft(times,  oneFluxes, freq, winXForm);

	// Zero all the high frequencies
	// Note if maxFreq > 0.5/offSpace, this code has no effect
	//	Nor should it, since this means the output grid is too coarse 
	//	to recover the high-frequency information
	// For now, this is a sharp cutoff
	// Introducing some apodizing may be helpful
	for(size_t i= 0; i < freq.size(); i++) {
		if (freq[i] > maxFreq) {
			winXForm[i] = std::complex<double>(0.0, 0.0);
		}
	}
	
	// Get the power spectrum
	for(ComplexVec::iterator it = winXForm.begin(); it != winXForm.end(); it++) {
		*it = norm(*it);
	}
	
	// Reverse transform -- setup
	size_t gslSize = 2*winXForm.size() - 1;
	gsl_fft_halfcomplex_wavetable* theTable = 
		gsl_fft_halfcomplex_wavetable_alloc(gslSize);
	gsl_fft_real_workspace * theSpace = 
		gsl_fft_real_workspace_alloc(gslSize);
		
	double *gslBuffer = new double[gslSize];
	gslBuffer[0] = winXForm[0].real();
	for (size_t i = 1; i < winXForm.size(); i++) {
		gslBuffer[2*i-1] = winXForm[i].real();
		gslBuffer[2*i  ] = winXForm[i].imag();
	}
	// assert: gsl_size is odd
	// therefore, gslBuffer follows the odd-transform convention, and only 
	//	xForm[0].imag() need be dropped
	gsl_fft_halfcomplex_transform(gslBuffer, 1, gslSize, theTable, theSpace);
	for(size_t i = 0; i < gslSize; i++) {
		double time = static_cast<double>(i)/((gslSize-1) * freqStep);
		// Values above tRange are aliases, so don't record them
		if (time <= tRange) {
			wf.push_back(gslBuffer[i]/gslBuffer[0]);
		} else {
			wf.push_back(1.0);
		}
	}
	// Clip or extend to match length of offsets
	wf.resize(offsets.size(), 1.0);
	
	// Clean up
	delete [] gslBuffer;
	gsl_fft_halfcomplex_wavetable_free(theTable);
	gsl_fft_real_workspace_free(theSpace);
}
