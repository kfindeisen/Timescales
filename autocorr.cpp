/** Autocorrelation function for unevenly sampled data
 * @file timescales/autocorr.cpp
 * @author Krzysztof Findeisen
 * @date Created February 16, 2011
 * @date Last modified November 19, 2013
 */ 

#include <complex>
#include <string>
#include <vector>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr.hpp>
#include <gsl/gsl_fft_halfcomplex.h>
#include "dft.h"
#include "timeexcept.h"
#include "timescales.h"
#include "../common/alloc.tmp.h"
#include "../common/stats.tmp.h"

namespace kpftimes {

using std::string;
using boost::lexical_cast;
using boost::scoped_array;
using boost::shared_ptr;
using kpfutils::checkAlloc;

/** Calculates the autocorrelation function for a time series. 
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] fluxes	Flux measurements of a source
 * @param[in] offsets	The time grid over which the autocorrelation function 
 *			should be calculated. 
 * @param[out] acf	The value of the autocorrelation function at each 
 *			offset.
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre @p fluxes.size() = @p times.size()
 * @pre @p fluxes[i] is the flux of the source at @p times[i], for all i
 * @pre @p offsets contains at least two unique elements
 * @pre @p offsets contains only nonnegative values
 * @pre @p offsets is uniformly sampled from 0 to some maximum value. This 
 *	requirement will be relaxed in future versions.
 * 
 * @post @p acf.size() = @p offsets.size()
 * @post @p acf[i] is the Scargle autocorrelation function evaluated at @p offsets[i], for all i
 * 
 * @perform O(FN) time, where N = @p times.size() and F = @p offsets.size()
 *
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some offsets are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p times and @p fluxes have 
 *	different lengths, if @p offsets has at most one distinct value, or if 
 *	it is not uniformly sampled.
 * @exception std::bad_alloc Thrown if there is not enough memory to perform 
 *	the calculations.
 *
 * @exceptsafe The function arguments are unchanged in the event of an exception.
 *
 * @todo Allow autoCorr() to run with arbitrary offset grids
 * @todo Verify that input validation is worth the cost
 * @todo Prove performance
 * @bug Many ACFs show a strong fringing effect. The fringes have a wavelength 
 * corresponding to the pseudo-Nyquist frequency, or to the highest peak 
 * frequency that is below the pseudo-Nyquist frequency if the power spectrum 
 * of the data is strongly peaked.
 */
void autoCorr(const DoubleVec &times, const DoubleVec &fluxes, 
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
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre @p fluxes.size() = @p times.size()
 * @pre @p fluxes[i] is the flux of the source at @p times[i], for all i
 * @pre @p offsets contains at least two unique elements
 * @pre @p offsets contains only nonnegative values
 * @pre @p offsets is uniformly sampled from 0 to some maximum value. This 
 *	requirement will be relaxed in future versions.
 * @pre @p maxFreq is positive
 * 
 * @post @p acf.size() = @p offsets.size()
 * @post @p acf[i] is the Scargle autocorrelation function evaluated at @p offsets[i], for all i
 * 
 * @perform O(FN) time, where N = @p times.size() and F = @p offsets.size()
 *
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some offsets are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p times and @p fluxes have 
 *	different lengths, if @p offsets has at most one distinct value, if 
 *	it is not uniformly sampled, or if @p maxFreq is non-positive.
 * @exception std::bad_alloc Thrown if there is not enough memory to perform 
 *	the calculations.
 *
 * @exceptsafe The function arguments are unchanged in the event of an exception.
 *
 * @todo Verify that input validation is worth the cost
 * @todo Prove performance
 * @todo Break up this function.
 */
void autoCorr(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &offsets, DoubleVec &acf, double maxFreq) {
	size_t nTimes  = times.size();
	size_t nOutput = offsets.size();
	double tRange  = deltaT(times);
	
	// Normalize the data
	// Scargle (1989) argues this is too crude, but a fitting method of 
	//	the kind he proposes is too slow and too inflexible
	// While we're at it, also test for non-uniqueness and sorting
	DoubleVec zeroFluxes(nTimes);
	double meanFlux = kpfutils::mean(fluxes.begin(), fluxes.end());
	bool diffValues = false, sortedTimes = true;
	for(size_t i = 0; i < nTimes; i++) {
		zeroFluxes[i] = fluxes[i] - meanFlux;

		if (!diffValues && times[i] != times.front()) {
			diffValues = true;
		}
		if (sortedTimes && i > 0 && times[i-1] > times[i]) {
			sortedTimes = false;
		}
	}
	
	// Construct the flat source we'll use to get the window function
	DoubleVec oneFluxes(nTimes, 1.0);
	
	// Verify the preconditions
	if (!diffValues) {
		throw except::BadLightCurve("Argument 'times' to autoCorr() contains only one unique value");
	} else if (!sortedTimes) {
		throw kpfutils::except::NotSorted("Argument 'times' to autoCorr()  is not sorted in ascending order");
	} else if (fluxes.size() != nTimes) {
		try {
			throw std::invalid_argument("Arguments 'times' and 'fluxes' to autoCorr() are not the same length (gave " 
				+ lexical_cast<string>(nTimes)        + " for times, " 
				+ lexical_cast<string>(fluxes.size()) + " for fluxes)");
		} catch(const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Arguments 'times' and 'fluxes' to autoCorr() are not the same length");
		}
	} else if (maxFreq <= 0.0) {
		try {
			throw std::invalid_argument("Argument 'maxFreq' to autoCorr() must be positive (gave " 
				+ lexical_cast<string>(maxFreq) + ")");
		} catch(const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Argument 'maxFreq' to autoCorr() must be positive");
		}
	}

	// Verify offsets
	if (offsets.size() < 2) {
		throw std::invalid_argument("autoCorr(): need at least two elements in offsets for a meaningful ACF");
	} else if (offsets[0] != 0.0) {
		throw std::invalid_argument("autoCorr(): first element of offsets must be zero (for now)");
	}
	if (offsets[0] < 0.0 || offsets[1] < 0.0) {
		throw except::NegativeFreq("autoCorr(): offsets must be nonnegative");
	}
	double offSpace = offsets[1] - offsets[0];
	if (offSpace <= 0.0) {
		throw std::invalid_argument("autoCorr(): offsets must be in ascending order");
	}
	for(size_t i = 2; i < nOutput; i++) {
		if (offsets[i] < 0.0 || offsets[i] <= offsets[i-1]) {
			throw except::NegativeFreq("autoCorr(): offsets must be nonnegative");
		}
		if (offsets[i] <= offsets[i-1]) {
			throw std::invalid_argument("autoCorr(): offsets must be in ascending order");
		}
		if (fabs(offsets[i] - offsets[i-1] - offSpace)/offSpace > 1e-3) {
			throw std::invalid_argument("autoCorr(): offsets must have uniform spacing (for now)");
		}
	}
	
	// Generate a regular frequency grid
	// A frequency step of 1/2 Delta T is sufficient to prevent wraparound
	// Going to a finer frequency spacing has no effect other than to 
	// 	introduce a string of zeros at offset > Delta T in the final 
	//	answer
	DoubleVec freq;
	double freqStep = 0.5/tRange;
	freqGen(times, freq, 0.0, 0.5/offSpace, 0.5);

	// Forward transform
	ComplexVec xForm, winXForm;
	dft(times, zeroFluxes, freq,    xForm);
	dft(times,  oneFluxes, freq, winXForm);

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
	shared_ptr<gsl_fft_halfcomplex_wavetable> theTable(checkAlloc(
		gsl_fft_halfcomplex_wavetable_alloc(gslSize)), 
		&gsl_fft_halfcomplex_wavetable_free);
	shared_ptr<gsl_fft_real_workspace> theSpace(checkAlloc(
		gsl_fft_real_workspace_alloc(gslSize)), 
		&gsl_fft_real_workspace_free);
	
	scoped_array<double> gslBuffer(new double[gslSize]);
	gslBuffer[0] = xForm[0].real();
	for (size_t i = 1; i < xForm.size(); i++) {
		gslBuffer[2*i-1] = xForm[i].real();
		gslBuffer[2*i  ] = xForm[i].imag();
	}
	// assert: gslSize is odd
	// therefore, gslBuffer follows the odd-transform convention, and only 
	//	xForm[0].imag() need be dropped
	
	// Reverse transform -- action
	gsl_fft_halfcomplex_transform(gslBuffer.get(), 1, gslSize, 
		theTable.get(), theSpace.get());
	
	// Reinterpret gslBuffer as autocorrelation function
	DoubleVec tempAcf;
	for(size_t i = 0; i < gslSize; i++) {
		double time = static_cast<double>(i)/((gslSize-1) * freqStep);
		// Values above tRange are aliases, so don't record them
		if (time <= tRange) {
			//offsets.push_back(time);
			tempAcf.push_back(gslBuffer[i]/gslBuffer[0]);
		} else {
			tempAcf.push_back(0.0);
		}
	}
	// Clip or extend to match length of offsets
	tempAcf.resize(offsets.size(), 0.0);
	
	// Now repeat for the sampling ACF
	gslBuffer[0] = winXForm[0].real();
	for (size_t i = 1; i < winXForm.size(); i++) {
		gslBuffer[2*i-1] = winXForm[i].real();
		gslBuffer[2*i  ] = winXForm[i].imag();
	}
	gsl_fft_halfcomplex_transform(gslBuffer.get(), 1, gslSize, 
			theTable.get(), theSpace.get());
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
	
	// Normalize
	for(size_t i = 0; i < tempAcf.size(); i++) {
		tempAcf[i] = tempAcf[i] / winAcf[i];
	}
	
	// IMPORTANT: no exceptions beyond this point
	
	using std::swap;
	swap(acf, tempAcf);
}

/** Calculates the autocorrelation window function for a time sampling. 
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] offsets	The time grid over which the autocorrelation function 
 *			should be calculated. 
 * @param[out] wf	The value of the window function at each offset.
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre @p offsets contains at least two unique elements
 * @pre @p offsets contains only nonnegative values
 * @pre @p offsets is uniformly sampled from 0 to some maximum value. This 
 *	requirement will be relaxed in future versions.
 * 
 * @post @p wf.size = @p offsets.size()
 * @post @p wf[i] is the Scargle autocorrelation function for a constant 
 *	function evaluated at offsets[i], for all i
 * 
 * @perform O(FN) time, where N = @p times.size() and F = @p offsets.size()
 * 
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some offsets are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p offsets has at most one 
 *	distinct value, or if it is not uniformly sampled.
 * @exception std::bad_alloc Thrown if there is not enough memory to perform 
 *	the calculations.
 *
 * @exceptsafe The function arguments are unchanged in the event of an exception.
 *
 * @todo Verify that input validation is worth the cost
 */
void acWindow(const DoubleVec &times, const DoubleVec &offsets, DoubleVec &wf) {
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
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre @p offsets contains at least two unique elements
 * @pre @p offsets contains only nonnegative values
 * @pre @p offsets is uniformly sampled from 0 to some maximum value. This 
 *	requirement will be relaxed in future versions.
 * @pre @p maxFreq is positive
 * 
 * @post wf is of the same length as offsets
 * @post wf[i] is the Scargle autocorrelation function evaluated at offsets[i], for all i
 * 
 * @perform O(FN) time, where N = @p times.size() and F = @p offsets.size()
 * 
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some offsets are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p offsets has at most one 
 *	distinct value, if it is not uniformly sampled, or if @p maxFreq is 
 *	non-negative.
 * @exception std::bad_alloc Thrown if there is not enough memory to perform 
 *	the calculations.
 *
 * @exceptsafe The function arguments are unchanged in the event of an exception.
 *
 * @todo Verify that input validation is worth the cost
 * @todo Prove performance
 * @todo Break up this function.
 */
void acWindow(const DoubleVec &times, const DoubleVec &offsets, DoubleVec &wf, 
		double maxFreq) {
	size_t nTimes  = times.size();
	size_t nOutput = offsets.size();
	double tRange  = deltaT(times);
	
	// Construct the flat source we'll use to get the 
	//	window function
	DoubleVec oneFluxes(nTimes, 1.0);

	// Test for non-uniqueness and sorting
	bool diffValues = false, sortedTimes = true;
	for(size_t i = 0; i < nTimes; i++) {
		if (!diffValues && times[i] != times.front()) {
			diffValues = true;
		}
		if (sortedTimes && i > 0 && times[i-1] > times[i]) {
			sortedTimes = false;
		}
	}
	
	// Verify the preconditions
	if (!diffValues) {
		throw except::BadLightCurve("Argument 'times' to acWindow() contains only one unique value");
	} else if (!sortedTimes) {
		throw kpfutils::except::NotSorted("Argument 'times' to acWindow() is not sorted in ascending order");
	} else if (maxFreq <= 0.0) {
		try {
			throw std::invalid_argument("Argument 'maxFreq' to acWindow() must be positive (gave " 
				+ lexical_cast<string>(maxFreq) + ")");
		} catch(const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Argument 'maxFreq' to acWindow() must be positive");
		}
	}

	// Verify offsets
	if (offsets.size() < 2) {
		throw std::invalid_argument("acWindow(): need at least two elements in offsets for a meaningful ACF");
	} else if (offsets[0] != 0.0) {
		throw std::invalid_argument("acWindow(): first element of offsets must be zero (for now)");
	}
	if (offsets[0] < 0.0 || offsets[1] < 0.0) {
		throw except::NegativeFreq("acWindow(): offsets must be nonnegative");
	}
	double offSpace = offsets[1] - offsets[0];
	if (offSpace <= 0.0) {
		throw std::invalid_argument("acWindow(): offsets must be in ascending order");
	}
	for(size_t i = 2; i < nOutput; i++) {
		if (offsets[i] < 0.0 || offsets[i] <= offsets[i-1]) {
			throw except::NegativeFreq("acWindow(): offsets must be nonnegative");
		}
		if (offsets[i] <= offsets[i-1]) {
			throw std::invalid_argument("acWindow(): offsets must be in ascending order");
		}
		if (fabs(offsets[i] - offsets[i-1] - offSpace)/offSpace > 1e-3) {
			throw std::invalid_argument("acWindow(): offsets must have uniform spacing (for now)");
		}
	}

	// Generate a regular frequency grid
	// A frequency step of 1/2 Delta T is sufficient to prevent wraparound
	// Going to a finer frequency spacing has no effect other than to 
	// 	introduce a string of zeros at offset > Delta T in the final 
	//	answer
	DoubleVec freq;
	double freqStep = 0.5/tRange;
	freqGen(times, freq, 0.0, 0.5/offSpace, 0.5);

	// Forward transform
	ComplexVec winXForm;
	dft(times,  oneFluxes, freq, winXForm);

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
	shared_ptr<gsl_fft_halfcomplex_wavetable> theTable(checkAlloc(
		gsl_fft_halfcomplex_wavetable_alloc(gslSize)), 
		&gsl_fft_halfcomplex_wavetable_free);
	shared_ptr<gsl_fft_real_workspace> theSpace(checkAlloc(
		gsl_fft_real_workspace_alloc(gslSize)), 
		&gsl_fft_real_workspace_free);
	
	scoped_array<double> gslBuffer(new double[gslSize]);
	gslBuffer[0] = winXForm[0].real();
	for (size_t i = 1; i < winXForm.size(); i++) {
		gslBuffer[2*i-1] = winXForm[i].real();
		gslBuffer[2*i  ] = winXForm[i].imag();
	}
	// assert: gsl_size is odd
	// therefore, gslBuffer follows the odd-transform convention, and only 
	//	xForm[0].imag() need be dropped
	gsl_fft_halfcomplex_transform(gslBuffer.get(), 1, gslSize, 
			theTable.get(), theSpace.get());
	
	DoubleVec tempWf;
	for(size_t i = 0; i < gslSize; i++) {
		double time = static_cast<double>(i)/((gslSize-1) * freqStep);
		// Values above tRange are aliases, so don't record them
		if (time <= tRange) {
			tempWf.push_back(gslBuffer[i]/gslBuffer[0]);
		} else {
			tempWf.push_back(1.0);
		}
	}
	// Clip or extend to match length of offsets
	tempWf.resize(offsets.size(), 1.0);
	
	// IMPORTANT: no exceptions beyond this point
	
	using std::swap;
	swap(wf, tempWf);
}

}	// end kpftimes
