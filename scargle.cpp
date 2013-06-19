/** Computes the Lomb-Scargle periodogram of an unevenly sampled lightcurve
 * @file scargle.cpp
 * @author Krzysztof Findeisen
 * @date Derived from scargle.pro (by Joern Wilms et al.) January 25, 2010
 * @date Last modified April 13, 2011
 */ 

/* Although this code is based on scargle.pro, several bugs have been fixed 
 *	and some features have been replaced with the recommendations from 
 *	Horne & Baliunas (1986):
 *	1. The old code made no distinction between the number of independent 
 *		frequencies and the (usually larger) number of frequencies in 
 *		the grid. This caused only the low-frequency end of a grid to 
 *		get calculated.
 *	2. The number of independent frequencies is intended only for 
 *		statistical analysis. For automatically generating grids, 
 *		I go with Delta f = 1/2T
 *	3. Removed the Horne & Baliunas analytic calculation of the detection 
 *		threshold. Detection thresholds are extremely sensitive to how the data 
 *		are sampled, and analytic fits are useful only for evenly 
 *		sampled data. [Changed Jan 18, 2011]
 *	4. The power spectrum calculated from Monte Carlo simulations now uses 
 *		the sample variance of the simulated data, not the sample 
 *		variance of the original data. The original code caused the 
 *		simulated data to give many more spurious peaks than it should 
 *		have. [Fixed Jan 18, 2011]
 */

/* Implementation note: this code extensively uses std::vectors, in its 
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

// If set to 1, Scargle algorithm runs in its slow (low-memory) form
// If set to 0 or not set, algorithm runs in its fast (high-memory) form
// This flag is best set in the call to the compiler, rather than in 
//	the code itself
#ifndef SCARGLE_SLOW 
#define SCARGLE_SLOW 0
#endif

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "timescales.h"
#include "utils.h"
#include "../common/stats.tmp.h"

#if USELFP
#include <lfp/lfp.h>
#endif

#ifndef PI
#define PI 3.1415927
#endif


/** Calculates the Lomb-Scargle periodogram for a time series.
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] fluxes	Flux measurements of a source
 * @param[in] freq	The frequency grid over which the periodogram should 
 *			be calculated. See freqGen() for a quick way to 
 *			generate a grid.
 * @param[out] power	The periodogram power at each frequency.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre fluxes is of the same length as times
 * @pre fluxes has at least two unique values
 * @pre fluxes[i] is the flux of the source at times[i], for all i
 * @pre all elements of freq are >= 0
 * @post power is of the same length as freq
 * @post power[i] is the Lomb-Scargle periodogram evaluated at freq[i], for all i
 * @exception domain_error Thrown if negative frequencies are provided
 * @exception invalid_argument Thrown if any of the preconditions on the 
 *	format of times or fluxes are violated.
 * 
 * @perform O(times.size() × freq.size()) time
 * @perfmore O(times.size() + freq.size()) memory
 * 
 * @todo Find a faster algorithm
 * @todo Optimize for multiple calls with similar [but not identical] values 
 *	of times and freq, as would happen in a large survey where some 
 *	epochs of some stars are removed for technical reasons.
 * @todo Verify that input validation is worth the cost
 */
void kpftimes::lombScargle(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &freq, DoubleVec &power) {
	size_t i, j;
	// Handy initializations
	size_t nTimes = times.size();
	
	// make times of manageable size (Scargle periodogram is time-shift invariant)
	// while we're at it, test for non-uniqueness
	bool diffValues = false, sortedTimes = true;
	DoubleVec times0(nTimes);
	double t0 = times.front();
	// Shift the times so t0 = 0
	for(i = 0; i < nTimes; i++) {
		times0[i] = times[i] - t0;
		if (!diffValues && times[i] != t0) {
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

	// Full sample variance
	double var   = kpfutils::variance(fluxes.begin(), fluxes.end());
	if (var <= 0.0) {
		throw std::invalid_argument("fluxes contains only one unique date");
	}

	size_t nFreq  = freq.size();
	// Equations are best expressed in angular frequency
	DoubleVec om(nFreq);
	for(i = 0; i < nFreq; i++) {
		if(freq[i] < 0) {
			throw std::domain_error("negative frequencies are not allowed");
		} else {
			om[i] = 2.0*PI*freq[i];
		}
	}

	////////////////////////////////
	// Periodogram
	// Ref.: W.H. Press and G.B. Rybicki, 1989, ApJ 338, 277

	DoubleVec cosOmTau(nFreq);
	DoubleVec sinOmTau(nFreq);
	DoubleVec tc2(nFreq);
	DoubleVec ts2(nFreq);
	for (i = 0; i < nFreq; i++) {
		double s2 = 0.0, c2 = 0.0;
		for (j = 0; j < nTimes; j++) {
			s2 += sin(2.0 * om[i] * times0[j]);
			c2 += cos(2.0 * om[i] * times0[j]);
		}
		
		// Eq. (2): Definition -> tan(2omtau)
		// --- tan(2omtau)  =  s2 / c2
		double omTau = 0.5 * atan2(s2, c2);
		cosOmTau[i] = cos(omTau);
		sinOmTau[i] = sin(omTau);
		
		// Eq. (7); total(cos(t-tau)^2) and total(sin(t-tau)^2) 
		double tmp = c2*cos(2.0*omTau) + s2*sin(2.0*omTau);
		tc2[i] = 0.5*(nTimes+tmp);		// total(cos(t-tau)^2)
		ts2[i] = 0.5*(nTimes-tmp);		// total(sin(t-tau)^2)
	}
   
	// computing the periodogram for the original lc
   
	// Subtract mean from data
	double meanF = kpfutils::mean(fluxes.begin(), fluxes.end());
	
	DoubleVec fluxes0(nTimes);
	for(i = 0; i < nTimes; i++) {
		fluxes0[i] = fluxes[i] - meanF;
	}

	////////////////////////////////
	// These 2-4 vectors depend on the data as well as the 
	//	experimental procedure
	// Eq. (5); sh and ch
	DoubleVec sh(nFreq);
	DoubleVec ch(nFreq);

	for (i=0; i < nFreq; i++) {
		sh[i] = 0.0;
		ch[i] = 0.0;
		for(j = 0; j < nTimes; j++) {
			sh[i] += fluxes0[j]*sin(om[i]*times0[j]);
			ch[i] += fluxes0[j]*cos(om[i]*times0[j]);
		}
	}

	// Clean up what we don't need
	fluxes0.clear();

	////////////////////////////////
	// Finally the periodogram itself

	power.clear();
	power.resize(nFreq);
	for(i=0; i < nFreq; i++) {
		// Eq. (3)
		if (om[i] != 0.0) {
			double cc = ch[i]*cosOmTau[i] + sh[i]*sinOmTau[i];
			double sc = sh[i]*cosOmTau[i] - ch[i]*sinOmTau[i];
			// correct normalization 
			power[i] = 0.5*(cc*cc  / tc2[i] + sc*sc  / ts2[i])/var;
		} else {
			// Use the limit as frequency goes to zero
			power[i] = 0.0;
		}
	}
}

/** Calculates the significance threshold for a Lomb-Scargle periodogram.
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] freq	The frequency grid over which the periodogram was 
 *			calculated.
 * @param[in] fap	Desired false alarm probability
 * @param[in] nSims	Number of white noise simulations to find the FAP 
 *			power level.
 *
 * @return The peak power level that will be reached, with probability fap, 
 *	in a periodogram of white noise.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre all elements of freq are >= 0
 * @pre 0 < fap < 1
 * @pre nSims >= 1
 * @pre fap × nSims >> 1
 * @post The function returns the peak power level observed in the 
 *	periodograms of (1-fap) of observations of an uncorrelated Gaussian 
 *	noise source, if the uncorrelated source is sampled at the cadence 
 *	represented by times and the periodogram is measured at frequencies 
 *	freq.
 * @exception domain_error Thrown if negative frequencies, fap, or nSims are 
 *	provided.
 * @exception invalid_argument Thrown if any of the preconditions on the 
 *	format of times or the values of fap and nSims are violated.
 * 
 * @perform O(times.size() × freq.size() × nSims) time
 * @perfmore O(times.size() × freq.size()) memory by default, or 
 *	O(times.size() + freq.size()) memory if SCARGLE_SLOW is set
 *
 * @remark The intended use is that %lsThreshold() will accompany a call, or multiple 
 *	calls, to lombScargle() with the same values of times and freq. Although 
 *	%lsThreshold() is optimized for multiple calculations to a degree that 
 *	%lombScargle() cannot, it will still run roughly nSims/10 times longer than 
 *	%lombScargle() itself. Don't run it after every real periodogram!
 *
 * @remark The significance threshold is a strong function of the observing 
 *	cadence (times), and a weaker function of the frequency grid 
 *	(freq). The latter dependence arises because in finer grids a 
 *	(random or observed) periodogram peak is more likely to be sampled 
 *	close to its maximum value.
 * 
 * @todo Test the performance advantage AFTER refactoring
 * @todo Verify that input validation is worth the cost
 */
double kpftimes::lsThreshold(const DoubleVec &times, const DoubleVec &freq, 
		double fap, long nSims) {
	size_t i, j;
	// Handy initializations
	size_t nTimes = times.size();
	
	// make times of manageable size (Scargle periodogram is time-shift invariant)
	// while we're at it, test for non-uniqueness
	bool diffValues = false, sortedTimes = true;
	DoubleVec times0(nTimes);
	double t0 = times.front();
	// Shift the times so t0 = 0
	for(i = 0; i < nTimes; i++) {
		times0[i] = times[i] - t0;
		if (!diffValues && times[i] != t0) {
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
	} else if (nSims < 1) {
		throw std::invalid_argument("Must run at least one simulation");
	} else if (fap <= 0.0 || fap >= 1.0) {
		throw std::invalid_argument("False alarm probability must be in the interval (0, 1)");
	} else if (nSims*fap < 10) {
		throw std::invalid_argument("Not enough simulations to get a significant peak at the desired false alarm probability");
	}

	size_t nFreq  = freq.size();
	// Equations are best expressed in angular frequency
	DoubleVec om(nFreq);
	for(i = 0; i < nFreq; i++) {
		if(freq[i] < 0) {
			throw std::domain_error("negative frequencies are not allowed");
		} else {
			om[i] = 2.0*PI*freq[i];
		}
	}

	////////////////////////////////
	// These four vectors are functions only of om, t, and t.size(), 
	//	i.e. they depend on the experimental procedure but 
	//	not the data. We can reuse them for each simulation
	// Eq. (6); s2, c2
	DoubleVec cosOmTau(nFreq);
	DoubleVec sinOmTau(nFreq);
	DoubleVec tc2(nFreq);
	DoubleVec ts2(nFreq);
	for (i = 0; i < nFreq; i++) {
		double s2 = 0.0, c2 = 0.0;
		for (j = 0; j < nTimes; j++) {
			s2 += sin(2.0 * om[i] * times0[j]);
			c2 += cos(2.0 * om[i] * times0[j]);
		}
		
		// Eq. (2): Definition -> tan(2omtau)
		// --- tan(2omtau)  =  s2 / c2
		double omTau = 0.5 * atan2(s2, c2);
		cosOmTau[i] = cos(omTau);
		sinOmTau[i] = sin(omTau);
		
		// Eq. (7); total(cos(t-tau)^2) and total(sin(t-tau)^2) 
		double tmp = c2*cos(2.0*omTau) + s2*sin(2.0*omTau);
		tc2[i] = 0.5*(nTimes+tmp);		// total(cos(t-tau)^2)
		ts2[i] = 0.5*(nTimes-tmp);		// total(sin(t-tau)^2)
	}
   
	////////////////////////////////
	// These 2 tables depend only on om and t as well
	// Eq. (5); sh and ch
	#if SCARGLE_SLOW==0
	FastTable sisi(nFreq, nTimes);
	FastTable coco(nFreq, nTimes);

	for (i=0; i < nFreq; i++) {
		for(j = 0; j < nTimes; j++) {
			sisi.at(i,j) = sin(om[i]*times0[j]);
			coco.at(i,j) = cos(om[i]*times0[j]);
		}
	}
	#endif

	////////////////////////////////
	// These 2 vectors depend on the data as well as the 
	//	experimental procedure
	// Eq. (5); sh and ch
	DoubleVec sh(nFreq);
	DoubleVec ch(nFreq);

	//////////////////////////////////////////////////////////////
	// Simulations
	
	DoubleVec psdPeak(nSims); 
	DoubleVec fluxes(nTimes);
	double meanF;
	gsl_rng* noiseGen = gsl_rng_alloc(gsl_rng_ranlxd2);
	// Not the cleverest seed, but we only choose one per run of lsThreshold()
	// I'm hoping that for any reasonable value of nSims two consecutive 
	//	runs of lsThreshold() should be far enough apart in time that 
	//	the seeds will be different
	time_t foo;
	gsl_rng_set(noiseGen, static_cast<unsigned long>(time(&foo)));

	for(long m = 0; m < nSims; m++) {
		#if USELFP
		PStart(1);
		#endif
		// white noise simulation
		meanF = 0.0;
		for (i = 0; i < nTimes; i++) {
			fluxes[i] = gsl_ran_ugaussian(noiseGen);
			meanF += fluxes[i];
		}
		meanF /= nTimes;
		for (i = 0; i < nTimes; i++) {
			fluxes[i] = fluxes[i]-meanF; 	// .. force OBSERVED count rate to zero
		}

		// Full sample variance
		double var   = kpfutils::variance(fluxes.begin(), fluxes.end());
		//double cnnoise = sqrt(cnVar);

		#if USELFP
		PStart(4);
		#endif
		for (i=0; i < nFreq; i++) {
			sh[i] = 0.0;
			ch[i] = 0.0;
			
			for(j = 0; j < nTimes; j++) {
				#if SCARGLE_SLOW==0
				sh[i] += fluxes[j]*sisi.at(i,j);
				ch[i] += fluxes[j]*coco.at(i,j);
				#else
				sh[i] += fluxes[j]*sin(om[i]*times0[j]);
				ch[i] += fluxes[j]*cos(om[i]*times0[j]);
				#endif
			}
		}
		#if USELFP
		PEnd(4);
		#endif

		// Eq. (3) ; computing the periodogram for each simulation
		//	Since we're only interested in the peak of the 
		//	periodogram, don't calculate the whole thing 
		//	-- just keep a running max in psdPeak[]
		double cc = ch[0]*cosOmTau[0] + sh[0]*sinOmTau[0];
		double sc = sh[0]*cosOmTau[0] - ch[0]*sinOmTau[0];
		psdPeak[m] = cc*cc / tc2[0] + sc*sc / ts2[0];
		for (i = 1; i < nFreq; i++) {
			double pp;
			if (om[i] != 0.0) {
				cc = ch[i]*cosOmTau[i] + sh[i]*sinOmTau[i];
				sc = sh[i]*cosOmTau[i] - ch[i]*sinOmTau[i];
				pp = cc*cc / tc2[i] + sc*sc / ts2[i];
			} else {
				// Use the limit as frequency goes to zero
				pp = 0.0;
			}

			if (pp > psdPeak[m]) {
				psdPeak[m] = pp;
			}
		}
		// correct normalization
		psdPeak[m] = 0.5 * psdPeak[m]/var;
		#if USELFP
		PEnd(1);
		#endif
	}

	// False Alarm Probability according to simulations
	// We have all the peaks, now find the fapth percentile
	// That's the level we expect to be exceeded one fapth of the time
	sort(psdPeak.begin(), psdPeak.end());
	return psdPeak[static_cast<size_t>((1.0 - fap) * (nSims - 1))];
}

/** Calculates the empirical distribution function of false peaks for a 
 *	Lomb-Scargle periodogram. This function is a generalization 
 *	of lsThreshold().
 *
 * @param[in] times	Times at which data were taken
 * @param[in] freqs	The frequency grid over which the periodogram was 
 *			calculated.
 * @param[out] powers	The power levels at which the EDF is measured
 * @param[out] probs	probs[i] contains the estimated probability that a 
 *			periodogram calculated from white noise has a peak 
 *			less than or equal to powers[i]
 * @param[in] nSims	Number of white noise simulations to find the EDF.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre all elements of freq are >= 0
 * @pre nSims >> 1
 * @post powers.size() == probs.size() == nSims
 * @post powers is sorted in ascending order
 * @post probs is sorted in ascending order
 * @post powers and probs represent an empirical distribution function, such 
 *	that probs[i] = EDF(powers[i]). The EDF is that of the peak power 
 *	level observed in the periodograms of observations of an uncorrelated 
 *	Gaussian noise source, if the uncorrelated source is sampled at the 
 *	cadence represented by times and the periodogram is measured at 
 *	frequencies freq.
 * 
 * @exception domain_error Thrown if negative frequencies or nonpositive nSims 
 *	are provided.
 * @exception invalid_argument Thrown if the preconditions on the 
 *	format of times are violated.
 *
 * @perform O(times.size() × freq.size() × nSims) time
 * @perfmore O(times.size() × freq.size()) memory by default, or 
 *	O(times.size() + freq.size()) memory if SCARGLE_SLOW is set
 *
 * @remark The intended use is that %lsRandomEdf() will accompany a call, or multiple 
 *	calls, to lombScargle() with the same values of times and freq. Although 
 *	%lsRandomEdf() is optimized for multiple calculations to a degree that 
 *	%lombScargle() cannot, it will still run roughly nSims/10 times longer than 
 *	%lombScargle() itself. Don't run it after every real periodogram!
 *
 * @test A 1-element time series and nSims = 100. Expected behavior = throw 
 *	invalid_argument
 * @test A 2-element time series, sorted with no duplicates, and nSims = 100. 
 *	Expected behavior = matches result of running lombScargle 100 times.
 * @test A 2-element time series, sorted with duplicates, and nSims = 100. 
 *	Expected behavior = throw invalid_argument
 * @test A 2-element time series, unsorted, and nSims = 100. Expected behavior 
 *	= throw invalid_argument.
 * @test A 100-element nonuniformly sampled time series, sorted with no 
 *	duplicates, and negative frequencies. Expected behavior = throw domain_error.
 * @test A 100-element nonuniformly sampled time series, sorted with no 
 *	duplicates, and multiple zero frequencies. Expected behavior = matches 
 *	result of running lombScargle 100 times.
 * @test A 100-element uniformly sampled time series, sorted with no 
 *	duplicates, and nSims = 100. Expected behavior = matches result of 
 *	running lombScargle 100 times.
 * @test A 100-element nonuniformly sampled time series, sorted with no 
 *	duplicates, and nSims = 100. Expected behavior = matches result of 
 *	running lombScargle 100 times.
 * @test A 100-element nonuniformly sampled time series, sorted with no 
 *	duplicates, and nSims = 0. Expected behavior = throw domain_error.
 * @test A 100-element nonuniformly sampled time series, sorted with no 
 *	duplicates, and nSims = 1. Expected behavior = (probs = 1.0, powers = 
 *	undefined)
 * @test A 100-element nonuniformly sampled time series, sorted with 
 *	duplicates, and nSims = 100. Expected behavior = matches result of 
 *	running lombScargle 100 times.
 * @test A 100-element nonuniformly sampled time series, unsorted. Expected 
 *	behavior = throw invalid_argument.
 */
void kpftimes::lsNormalEdf(const DoubleVec &times, const DoubleVec &freqs, 
		DoubleVec &powers, DoubleVec &probs, long nSims) {
	size_t i, j;
	// Handy initializations
	size_t nTimes = times.size();
	
	// make times of manageable size (Scargle periodogram is time-shift invariant)
	// while we're at it, test for non-uniqueness
	bool diffValues = false, sortedTimes = true;
	DoubleVec times0(nTimes);
	double t0 = times.front();
	// Shift the times so t0 = 0
	for(i = 0; i < nTimes; i++) {
		times0[i] = times[i] - t0;
		if (!diffValues && times[i] != t0) {
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
	} else if (nSims < 1) {
		throw std::domain_error("Must run at least one simulation");
	}

	size_t nFreq  = freqs.size();
	// Equations are best expressed in angular frequency
	DoubleVec om(nFreq);
	for(i = 0; i < nFreq; i++) {
		if(freqs[i] < 0) {
			throw std::domain_error("negative frequencies are not allowed");
		} else {
			om[i] = 2.0*PI*freqs[i];
		}
	}

	////////////////////////////////
	// These four vectors are functions only of om, t, and t.size(), 
	//	i.e. they depend on the experimental procedure but 
	//	not the data. We can reuse them for each simulation
	// Eq. (6); s2, c2
	DoubleVec cosOmTau(nFreq);
	DoubleVec sinOmTau(nFreq);
	DoubleVec tc2(nFreq);
	DoubleVec ts2(nFreq);
	for (i = 0; i < nFreq; i++) {
		double s2 = 0.0, c2 = 0.0;
		for (j = 0; j < nTimes; j++) {
			s2 += sin(2.0 * om[i] * times0[j]);
			c2 += cos(2.0 * om[i] * times0[j]);
		}
		
		// Eq. (2): Definition -> tan(2omtau)
		// --- tan(2omtau)  =  s2 / c2
		double omTau = 0.5 * atan2(s2, c2);
		cosOmTau[i] = cos(omTau);
		sinOmTau[i] = sin(omTau);
		
		// Eq. (7); total(cos(t-tau)^2) and total(sin(t-tau)^2) 
		double tmp = c2*cos(2.0*omTau) + s2*sin(2.0*omTau);
		tc2[i] = 0.5*(nTimes+tmp);		// total(cos(t-tau)^2)
		ts2[i] = 0.5*(nTimes-tmp);		// total(sin(t-tau)^2)
	}

	////////////////////////////////
	// These 2 tables depend only on om and t as well
	// Eq. (5); sh and ch
	#if SCARGLE_SLOW==0
	FastTable sisi(nFreq, nTimes);
	FastTable coco(nFreq, nTimes);

	for (i=0; i < nFreq; i++) {
		for(j = 0; j < nTimes; j++) {
			sisi.at(i,j) = sin(om[i]*times0[j]);
			coco.at(i,j) = cos(om[i]*times0[j]);
		}
	}
	#endif

	////////////////////////////////
	// These 2 vectors depend on the data as well as the 
	//	experimental procedure
	// Eq. (5); sh and ch
	DoubleVec sh(nFreq);
	DoubleVec ch(nFreq);

	//////////////////////////////////////////////////////////////
	// Simulations
	
	DoubleVec fluxes(nTimes);
	powers.clear();
	powers.resize(nSims);
	double meanF;
	gsl_rng* noiseGen = gsl_rng_alloc(gsl_rng_ranlxd2);
	// Not the cleverest seed, but we only choose one per run of lsThreshold()
	// I'm hoping that for any reasonable value of nSims two consecutive 
	//	runs of lsThreshold() should be far enough apart in time that 
	//	the seeds will be different
	time_t foo;
	gsl_rng_set(noiseGen, static_cast<unsigned long>(time(&foo)));

	for(long m = 0; m < nSims; m++) {
		#if USELFP
		PStart(1);
		#endif
		// white noise simulation
		meanF = 0.0;
		for (i = 0; i < nTimes; i++) {
			fluxes[i] = gsl_ran_ugaussian(noiseGen);
			meanF += fluxes[i];
		}
		meanF /= nTimes;
		for (i = 0; i < nTimes; i++) {
			fluxes[i] = fluxes[i]-meanF; 	// .. force OBSERVED count rate to zero
		}

		// Full sample variance
		double var   = kpfutils::variance(fluxes.begin(), fluxes.end());

		#if USELFP
		PStart(4);
		#endif
		for (i=0; i < nFreq; i++) {
			sh[i] = 0.0;
			ch[i] = 0.0;
			
			for(j = 0; j < nTimes; j++) {
				#if SCARGLE_SLOW==0
				sh[i] += fluxes[j]*sisi.at(i,j);
				ch[i] += fluxes[j]*coco.at(i,j);
				#else
				sh[i] += fluxes[j]*sin(om[i]*times0[j]);
				ch[i] += fluxes[j]*cos(om[i]*times0[j]);
				#endif
			}
		}
		#if USELFP
		PEnd(4);
		#endif

		// Eq. (3) ; computing the periodogram for each simulation
		//	Since we're only interested in the peak of the 
		//	periodogram, don't calculate the whole thing 
		//	-- just keep a running max in powers[]
		double cc = ch[0]*cosOmTau[0] + sh[0]*sinOmTau[0];
		double sc = sh[0]*cosOmTau[0] - ch[0]*sinOmTau[0];
		powers[m] = cc*cc / tc2[0] + sc*sc / ts2[0];
		for (i = 1; i < nFreq; i++) {
			double pp;
			if (om[i] != 0.0) {
				cc = ch[i]*cosOmTau[i] + sh[i]*sinOmTau[i];
				sc = sh[i]*cosOmTau[i] - ch[i]*sinOmTau[i];
				pp = cc*cc / tc2[i] + sc*sc / ts2[i];
			} else {
				// Use the limit as frequency goes to zero
				pp = 0.0;
			}

			if (pp > powers[m]) {
				powers[m] = pp;
			}
		}
		// correct normalization
		powers[m] = 0.5 * powers[m]/var;
		#if USELFP
		PEnd(1);
		#endif
	}

	sort(powers.begin(), powers.end());
	probs.clear();
	for(size_t i = 1; i <= powers.size(); i++) {
		probs.push_back(static_cast<double>(i)/powers.size());
	}
}
