/** Computes the Lomb-Scargle periodogram of an unevenly sampled lightcurve
 * @file timescales/scargle.cpp
 * @author Krzysztof Findeisen
 * @date Derived from scargle.pro (by Joern Wilms et al.) January 25, 2010
 * @date Last modified November 19, 2013
 */ 

#include <algorithm>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/smart_ptr.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "timescales.h"
#include "utils.h"
#include "../common/alloc.tmp.h"
#include "../common/stats.tmp.h"
#include "timeexcept.h"

namespace kpftimes {

using std::string;
using boost::lexical_cast;
using boost::math::double_constants::two_pi;
using boost::shared_ptr;
using kpfutils::checkAlloc;

/** Calculates the Lomb-Scargle periodogram for a time series.
 * 
 * @param[in] times	Times at which @p data were taken
 * @param[in] data	Measurements of a time series
 * @param[in] freqs	The frequency grid over which the periodogram should 
 *			be calculated. See freqGen() for a quick way to 
 *			generate a grid.
 * @param[out] power	The periodogram power at each frequency.
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre @p data.size() = @p times.size()
 * @pre @p data has at least two unique values
 * @pre @p data[i] is the measurement of the source at @p times[i], for all i
 * @pre all elements of @p freqs are &ge; 0
 * 
 * @post @p power.size() = @p freqs.size()
 * @post @p power[i] is the Lomb-Scargle periodogram evaluated at @p freqs[i], for all i
 *
 * @perform O(NF) time, where N = @p times.size() and F = @p freqs.size()
 * @perfmore O(N + F) memory
 * 
 * @exception kpftimes::except::BadLightCurve Thrown if @p times or @p data has 
 *	at most one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some elements of @p freqs are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p times and @p data have 
 *	different lengths.
 * @exception std::bad_alloc Thrown if there is not enough memory to do the 
 *	calculations.
 *
 * @exceptsafe The function arguments are unchanged in the event of an exception
 * 
 * @todo Factor this function
 * @todo Find a faster algorithm
 * @todo Optimize for multiple calls with similar (but not identical) values 
 *	of @p times and @p freqs, as would happen in a large survey where some 
 *	epochs of some stars are removed for technical reasons.
 * @todo Verify that input validation is worth the cost
 */
void lombScargle(const DoubleVec &times, const DoubleVec &data, 
		const DoubleVec &freqs, DoubleVec &power) {
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
		throw except::BadLightCurve("Parameter 'times' in lombScargle() contains only one unique date");
	} else if (!sortedTimes) {
		throw kpfutils::except::NotSorted("Parameter 'times' in lombScargle() is not sorted in ascending order");
	} else if (data.size() != nTimes) {
		try {
			throw std::invalid_argument("Parameters 'times' and 'data' in lombScargle() are not the same length (gave " 
			+ lexical_cast<string>(nTimes) + " for times and " 
			+ lexical_cast<string>(data.size()) + " for data)");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameters 'times' and 'data' in lombScargle() are not the same length");
		}
	}

	// Full sample variance
	double var   = kpfutils::variance(data.begin(), data.end());
	if (var <= 0.0) {
		throw except::BadLightCurve("Parameter 'data' in lombScargle() has no variability");
	}

	size_t nFreq  = freqs.size();
	// Equations are best expressed in angular frequency
	DoubleVec om(nFreq);
	for(i = 0; i < nFreq; i++) {
		if(freqs[i] < 0) {
			throw except::NegativeFreq("Parameter 'freqs' in lombScargle() contains negative frequencies");
		} else {
			om[i] = two_pi*freqs[i];
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
	double meanF = kpfutils::mean(data.begin(), data.end());
	
	DoubleVec data0(nTimes);
	for(i = 0; i < nTimes; i++) {
		data0[i] = data[i] - meanF;
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
			sh[i] += data0[j]*sin(om[i]*times0[j]);
			ch[i] += data0[j]*cos(om[i]*times0[j]);
		}
	}

	// Clean up what we don't need
	data0.clear();

	////////////////////////////////
	// Finally the periodogram itself

	// copy-and-swap
	DoubleVec tempPower(nFreq);
	
	for(i=0; i < nFreq; i++) {
		// Eq. (3)
		if (om[i] != 0.0) {
			double cc = ch[i]*cosOmTau[i] + sh[i]*sinOmTau[i];
			double sc = sh[i]*cosOmTau[i] - ch[i]*sinOmTau[i];
			// correct normalization 
			tempPower[i] = 0.5*(cc*cc  / tc2[i] + sc*sc  / ts2[i])/var;
		} else {
			// Use the limit as frequency goes to zero
			tempPower[i] = 0.0;
		}
	}
	
	using std::swap;
	swap(power, tempPower);
}

/** Calculates the significance threshold for a Lomb-Scargle periodogram.
 * 
 * The intended use is that lsThreshold() will accompany a call, or multiple 
 *	calls, to lombScargle() with the same values of @p times and @p freqs. Although 
 *	lsThreshold() is optimized for multiple calculations to a degree that 
 *	lombScargle() cannot, it will still run roughly nSims/10 times longer than 
 *	lombScargle() itself. Don't run it after every real periodogram!
 *
 * The significance threshold is a strong function of the observing 
 *	cadence (@p times), and a weaker function of the frequency grid 
 *	(@p freqs). The latter dependence arises because in finer grids a 
 *	(random or observed) periodogram peak is more likely to be sampled 
 *	close to its maximum value.
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] freqs	The frequency grid over which the periodogram was 
 *			calculated.
 * @param[in] fap	Desired false alarm probability
 * @param[in] nSims	Number of Gaussian white noise simulations to find 
 *			the FAP power level.
 *
 * @return The peak power level that will be reached, with probability @p fap, 
 *	in a periodogram of Gaussian white noise.
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre all elements of @p freqs are &ge; 0
 * @pre 0 < @p fap < 1
 * @pre @p nSims &ge; 1
 * @pre @p fap × @p nSims >> 1
 * 
 * @post The function returns the peak power level observed in the 
 *	periodograms of (1-@p fap) of observations of an uncorrelated Gaussian 
 *	noise source, if the uncorrelated source is sampled at the cadence 
 *	represented by @p times and the periodogram is measured at frequencies 
 *	@p freqs.
 *
 * @perform O(NF × @p nSims) time, where N = @p times.size() and F = @p freqs.size()
 * @perfmore O(NF) memory
 * 
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has 
 *	at most one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some elements of @p freqs are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p fap is outside (0, 1) or 
 *	if @p nSims is nonpositive
 * @exception std::bad_alloc Thrown if there is not enough memory to do the 
 *	calculations.
 * 
 * @exceptsafe The function arguments are unchanged in the event of an exception
 *
 * @todo Factor this function
 * @todo Test the performance advantage AFTER refactoring
 * @todo Verify that input validation is worth the cost
 */
double lsThreshold(const DoubleVec &times, const DoubleVec &freqs, 
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
		throw except::BadLightCurve("Parameter 'times' in lsThreshold() contains only one unique date");
	} else if (!sortedTimes) {
		throw kpfutils::except::NotSorted("Parameter 'times' in lsThreshold() is not sorted in ascending order");
	} else if (nSims < 1) {
		try {
			throw std::invalid_argument("Must run at least one simulation in lsThreshold() (gave " + lexical_cast<string>(nSims) + ")");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Must run at least one simulation in lsThreshold().");
		}
	} else if (fap >= 1.0 || fap <= 0.0) {
		try {
			throw std::invalid_argument("False alarm probability in lsThreshold() must be in the interval (0, 1) (gave " + lexical_cast<string>(fap) + ")");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("False alarm probability in lsThreshold() must be in the interval (0, 1)");
		}
	} else if (nSims*fap < 10) {
		throw std::invalid_argument("Not enough simulations in lsThreshold() to get a significant peak at the desired false alarm probability");
	}

	size_t nFreq  = freqs.size();
	// Equations are best expressed in angular frequency
	DoubleVec om(nFreq);
	for(i = 0; i < nFreq; i++) {
		if(freqs[i] < 0) {
			throw except::NegativeFreq("Parameter 'freqs' in lsThreshold() contains negative frequencies");
		} else {
			om[i] = two_pi*freqs[i];
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
	FastTable sisi(nFreq, nTimes);
	FastTable coco(nFreq, nTimes);

	for (i=0; i < nFreq; i++) {
		for(j = 0; j < nTimes; j++) {
			sisi.at(i,j) = sin(om[i]*times0[j]);
			coco.at(i,j) = cos(om[i]*times0[j]);
		}
	}

	////////////////////////////////
	// These 2 vectors depend on the data as well as the 
	//	experimental procedure
	// Eq. (5); sh and ch
	DoubleVec sh(nFreq);
	DoubleVec ch(nFreq);

	//////////////////////////////////////////////////////////////
	// Simulations
	
	DoubleVec psdPeak(nSims); 
	DoubleVec data(nTimes);
	double meanF;
	shared_ptr<gsl_rng> noiseGen(checkAlloc(
		gsl_rng_alloc(gsl_rng_ranlxd2)), &gsl_rng_free);
	// Not the cleverest seed, but we only choose one per run of lsThreshold()
	// I'm hoping that for any reasonable value of nSims two consecutive 
	//	runs of lsThreshold() should be far enough apart in time that 
	//	the seeds will be different
	time_t foo;
	gsl_rng_set(noiseGen.get(), static_cast<unsigned long>(time(&foo)));

	for(long m = 0; m < nSims; m++) {
		// white noise simulation
		meanF = 0.0;
		for (i = 0; i < nTimes; i++) {
			data[i] = gsl_ran_ugaussian(noiseGen.get());
			meanF += data[i];
		}
		meanF /= nTimes;
		for (i = 0; i < nTimes; i++) {
			data[i] = data[i]-meanF; 	// .. force OBSERVED count rate to zero
		}

		// Full sample variance
		double var   = kpfutils::variance(data.begin(), data.end());

		for (i=0; i < nFreq; i++) {
			sh[i] = 0.0;
			ch[i] = 0.0;
			
			for(j = 0; j < nTimes; j++) {
				sh[i] += data[j]*sisi.at(i,j);
				ch[i] += data[j]*coco.at(i,j);
			}
		}

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
	}

	// False Alarm Probability according to simulations
	// We have all the peaks, now find the fapth percentile
	// That's the level we expect to be exceeded one fapth of the time
	return kpfutils::quantile(psdPeak.begin(), psdPeak.end(), 1.0 - fap);
}

/** Calculates the empirical distribution function of false peaks for a 
 *	Lomb-Scargle periodogram. This function is a generalization 
 *	of lsThreshold().
 *
 * The intended use is that lsNormalEdf() will accompany a call, or multiple 
 *	calls, to lombScargle() with the same values of @p times and @p freqs. Although 
 *	lsNormalEdf() is optimized for multiple calculations to a degree that 
 *	lombScargle() cannot, it will still run roughly @p nSims/10 times longer than 
 *	lombScargle() itself. Don't run it after every real periodogram!
 *
 * @param[in] times	Times at which data were taken
 * @param[in] freqs	The frequency grid over which the periodogram was 
 *			calculated.
 * @param[out] powers	The power levels at which the EDF is measured
 * @param[out] probs	The probability that a periodogram calculated from 
 *			white noise has a peak less than or equal to 
 *			a power level
 * @param[in] nSims	Number of Gaussian white noise simulations to find 
 *			the EDF.
 * 
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre all elements of @p freqs are &ge; 0
 * @pre @p nSims &ge; 1
 * 
 * @post @p powers.size() == @p probs.size() == @p nSims
 * @post @p powers is sorted in ascending order
 * @post @p probs is sorted in ascending order
 * 
 * @post @p powers and @p probs represent an empirical distribution function, 
 *	such that probs[i] = EDF(powers[i]). The EDF is that of the peak power 
 *	level observed in the periodograms of observations of an uncorrelated 
 *	Gaussian noise source, if the uncorrelated source is sampled at the 
 *	cadence represented by @p times and the periodogram is measured at 
 *	frequencies @p freqs.
 * 
 * @perform O(NF × @p nSims) time, where N = @p times.size() and F = @p freqs.size()
 * @perfmore O(NF) memory
 *
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has 
 *	at most one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception kpftimes::except::NegativeFreq Thrown if some elements of @p freqs are 
 *	negative.
 * @exception std::invalid_argument Thrown if @p nSims is nonpositive
 * @exception std::bad_alloc Thrown if there is not enough memory to do the 
 *	calculations.
 *
 * @exceptsafe The function arguments are unchanged in the event of an exception
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
 *
 * @todo Factor this function
 */
void lsNormalEdf(const DoubleVec &times, const DoubleVec &freqs, 
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
		throw except::BadLightCurve("Parameter 'times' in lsNormalEdf() contains only one unique date");
	} else if (!sortedTimes) {
		throw kpfutils::except::NotSorted("Parameter 'times' in lsNormalEdf() is not sorted in ascending order");
	} else if (nSims < 1) {
		try {
			throw std::invalid_argument("Must run at least one simulation in lsNormalEdf() (gave " + lexical_cast<string>(nSims) + ")");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Must run at least one simulation in lsNormalEdf().");
		}
	}

	size_t nFreq  = freqs.size();
	// Equations are best expressed in angular frequency
	DoubleVec om(nFreq);
	for(i = 0; i < nFreq; i++) {
		if(freqs[i] < 0) {
			throw except::NegativeFreq("Parameter 'freq' in lsThreshold() contains negative frequencies");
		} else {
			om[i] = two_pi*freqs[i];
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
	FastTable sisi(nFreq, nTimes);
	FastTable coco(nFreq, nTimes);

	for (i=0; i < nFreq; i++) {
		for(j = 0; j < nTimes; j++) {
			sisi.at(i,j) = sin(om[i]*times0[j]);
			coco.at(i,j) = cos(om[i]*times0[j]);
		}
	}

	////////////////////////////////
	// These 2 vectors depend on the data as well as the 
	//	experimental procedure
	// Eq. (5); sh and ch
	DoubleVec sh(nFreq);
	DoubleVec ch(nFreq);

	//////////////////////////////////////////////////////////////
	// Simulations
	
	DoubleVec data(nTimes);
	DoubleVec tempPowers(nSims);
	double meanF;
	gsl_rng* noiseGen = gsl_rng_alloc(gsl_rng_ranlxd2);
	// Not the cleverest seed, but we only choose one per run of lsThreshold()
	// I'm hoping that for any reasonable value of nSims two consecutive 
	//	runs of lsThreshold() should be far enough apart in time that 
	//	the seeds will be different
	time_t foo;
	gsl_rng_set(noiseGen, static_cast<unsigned long>(time(&foo)));

	for(long m = 0; m < nSims; m++) {
		// white noise simulation
		meanF = 0.0;
		for (i = 0; i < nTimes; i++) {
			data[i] = gsl_ran_ugaussian(noiseGen);
			meanF += data[i];
		}
		meanF /= nTimes;
		for (i = 0; i < nTimes; i++) {
			data[i] = data[i]-meanF; 	// .. force OBSERVED count rate to zero
		}

		// Full sample variance
		double var   = kpfutils::variance(data.begin(), data.end());

		for (i=0; i < nFreq; i++) {
			sh[i] = 0.0;
			ch[i] = 0.0;
			
			for(j = 0; j < nTimes; j++) {
				sh[i] += data[j]*sisi.at(i,j);
				ch[i] += data[j]*coco.at(i,j);
			}
		}

		// Eq. (3) ; computing the periodogram for each simulation
		//	Since we're only interested in the peak of the 
		//	periodogram, don't calculate the whole thing 
		//	-- just keep a running max in tempPowers[]
		double cc = ch[0]*cosOmTau[0] + sh[0]*sinOmTau[0];
		double sc = sh[0]*cosOmTau[0] - ch[0]*sinOmTau[0];
		tempPowers[m] = cc*cc / tc2[0] + sc*sc / ts2[0];
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

			if (pp > tempPowers[m]) {
				tempPowers[m] = pp;
			}
		}
		// correct normalization
		tempPowers[m] = 0.5 * tempPowers[m]/var;
	}

	sort(tempPowers.begin(), tempPowers.end());
	size_t nPowers = tempPowers.size();
	DoubleVec tempProbs;
	for(size_t i = 1; i <= nPowers; i++) {
		tempProbs.push_back(static_cast<double>(i)/nPowers);
	}
	
	// IMPORTANT: no exceptions beyond this point
	
	using std::swap;
	swap(powers, tempPowers);
	swap(probs , tempProbs );
}

}		// end kpftimes
