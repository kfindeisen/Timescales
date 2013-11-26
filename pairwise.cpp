/** Computes delta-T pair plots, including self-correlation functions and gap plots
 * @file timescales/pairwise.cpp
 * @author Krzysztof Findeisen
 * @date Created July 24, 2011
 * @date Last modified October 29, 2013
 */ 

#include <algorithm>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include "../common/stats_except.h"
#include "../common/stats.tmp.h"
#include "timeexcept.h"
#include "timescales.h"

namespace kpftimes {

using std::string;
using boost::lexical_cast;

/** Calculates a &Delta;m&Delta;t plot.
 * 
 * @param[in] times	Times at which @p mags were taken
 * @param[in] mags	Magnitude measurements of a source
 * @param[out] deltaT	A list of the time intervals between all pairs of sources.
 * @param[out] deltaM	A list of the magnitude difference between each pair in @p deltaT.
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * @pre @p mags.size() = @p times.size()
 * @pre @p mags[i] is the magnitude of the source at @p times[i], for all i
 *
 * @post @p deltaT.size() = @p deltaM.size() = N(N-1)/2, where N = @p times.size()
 * @post @p deltaT is sorted in ascending order
 * @post Each element of @p deltaT represents the absolute time separation of 
 * 	a unique pair of values in @p times, and each element of @p deltaM 
 *	represents the absolute magnitude separation of a corresponding pair 
 *	of values in @p mags
 * 
 * @perform O(N<sup>2</sup>) time, where N = times.size()
 *
 * @exception kpftimes::except::BadLightCurve Thrown if @p times has at most 
 *	one distinct value.
 * @exception kpfutils::except::NotSorted Thrown if @p times is not in 
 *	ascending order.
 * @exception std::invalid_argument Thrown if @p times and @p mags have 
 *	different lengths.
 * @exception std::bad_alloc Thrown if there is not enough memory to compute 
 *	the &Delta;m&Delta;t plot
 * 
 * @exceptsafe The function arguments are unchanged in the event of an exception.
 */
void dmdt(const DoubleVec &times, const DoubleVec &mags, 
		DoubleVec &deltaT, DoubleVec &deltaM) {
	size_t nTimes  = times.size();
	
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
		throw except::BadLightCurve("Parameter 'times' in dmdt() contains only one unique date");
	} else if (!sortedTimes) {
		throw kpfutils::except::NotSorted("Parameter 'times' in dmdt() is not sorted in ascending order");
	} else if (mags.size() != nTimes) {
		try {
			throw std::invalid_argument("Parameters 'times' and 'mags' in dmdt() are not the same length (gave " 
			+ lexical_cast<string>(nTimes) + " for times and " 
			+ lexical_cast<string>(mags.size()) + " for mags)");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameters 'times' and 'mags' in dmdt() are not the same length");
		}
	}

	// Needed for sorting by deltaT
	typedef std::vector<std::pair<double, double> > pairVec;
	pairVec sortableVec;
	
	for(size_t i = 0; i < nTimes; i++) {
		for(size_t j = i+1; j < nTimes; j++) {
			// Time must be first element so that the pairs get sorted properly
			sortableVec.push_back(std::make_pair(fabs(times[i]-times[j]), 
					fabs(mags[i]-mags[j]) ));
		}
	}
	
	// Now we just need to sort it
	// I don't have a parallel sort function yet
	// For now, do it the clumsy way, using pair<>
	std::sort(sortableVec.begin(), sortableVec.end());
	
	// copy-and-swap
	DoubleVec tempTimes, tempMags;
	
	for(pairVec::const_iterator it = sortableVec.begin(); it != sortableVec.end(); it++) {
		tempTimes.push_back(it->first );
		tempMags .push_back(it->second);
	}
	
	using std::swap;
	swap(deltaT, tempTimes);
	swap(deltaM, tempMags );
}

/** Computes the fraction of pairs of magnitudes above some threshold found 
 *	in each &Delta;t bin of a &Delta;m&Delta;t plot.
 *
 * @param[in] deltaT The &Delta;t values of the &Delta;m&Delta;t plot.
 * @param[in] deltaM The corresponding &Delta;m values of the &Delta;m&Delta;t plot.
 * @param[in] binEdges A vector containing the (N+1) boundaries of the N &Delta;t bins 
 *	in which to count high-&Delta;m pairs.
 * @param[out] fracs A vector containing the fraction of &Delta;m values in 
 *	each bin that exceed threshold.
 * @param[in] threshold The characteristic magnitude difference, in magnitudes, 
 *	above which &Delta;m values are to be counted.
 *
 * @pre @p deltaT.size() = @p deltaM.size()
 * @pre @p deltaT is sorted in ascending order
 * @pre @p binEdges is sorted in ascending order
 *
 * @pre @p deltaT does not contain any NaNs
 * @pre @p deltaM does not contain any NaNs
 * @pre @p binEdges does not contain any NaNs
 * 
 * @post @p fracs.size() = @p binEdges.size() - 1
 * @post For all i &isin; [0, @p binEdges.size()-1], @p fracs[i] contains the 
 *	fraction of @p deltaM > @p threshold, given deltaT &isin; 
 *	[@p binEdges[i], @p binEdges[i+1]).
 *
 * @perform O(M + N) time, where N = @p deltaM.size() and M = @p binEdges.size()
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to store 
 *	the bin fractions.
 * @exception kpfutils::except::NotSorted Thrown if either @p deltaT or 
 *	@p binEdges is unsorted.
 * @exception std::invalid_argument Thrown if @p deltaT.size() &ne; @p deltaM.size()
 *
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 */
void hiAmpBinFrac(const DoubleVec &deltaT, const DoubleVec &deltaM, 
		const DoubleVec &binEdges, DoubleVec &fracs, double threshold) {
	using std::swap;

	if (deltaT.size() != deltaM.size()) {
		try {
			throw std::invalid_argument("Parameters 'deltaT' and 'deltaM' in hiAmpBinFrac() are not the same length (gave " 
			+ lexical_cast<string>(deltaT.size()) + " for deltaT and " 
			+ lexical_cast<string>(deltaM.size()) + " for deltaM)");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameters 'deltaT' and 'deltaM' in hiAmpBinFrac() are not the same length");
		}
	}
	if (!kpfutils::isSorted(deltaT.begin(), deltaT.end())) {
		throw kpfutils::except::NotSorted("deltaT is not sorted in hiAmpBinFrac()");
	}
	if (!kpfutils::isSorted(binEdges.begin(), binEdges.end())) {
		throw kpfutils::except::NotSorted("binEdges is not sorted in hiAmpBinFrac()");
	}
	
	DoubleVec tempFracs;
	tempFracs.reserve(binEdges.size()-1);

	DoubleVec::const_iterator curDelta = std::lower_bound(deltaT.begin(), deltaT.end(), 
			binEdges.front());
	// assert: if any deltaT values fall within the bins specified by binEdges, 
	//	then curDelta is the lowest such value
	for(DoubleVec::const_iterator curBin = binEdges.begin(); curBin+1 != binEdges.end(); 
			curBin++) {
		// loop invariant: the pair pointed to by curDelta has not 
		//	been counted for this or any previous bin

		// Restart the count for each bin
		long numPairs = 0, numHighPairs = 0;
		while (curDelta+1 != deltaT.end() && 
				*curDelta >= *curBin && *curDelta < *(curBin+1)) {
			numPairs++;
			
			// what is the magnitude difference for the pair pointed to by curDelta?
			DoubleVec::const_iterator curDeltaM = deltaM.begin() 
					+ (curDelta - deltaT.begin());
			if (*curDeltaM > threshold) {
				numHighPairs++;
			}
			
			curDelta++;
		} // else no new pairs to add
		// assert: *curDelta > *(curBin+1) xor all pairs in bin already added
		
		tempFracs.push_back(numPairs > 0 ? 
				static_cast<double>(numHighPairs)/
				static_cast<double>(numPairs) : 
				std::numeric_limits<double>::signaling_NaN());
	}
	
	// IMPORTANT: no exceptions beyond this point
	
	swap(fracs, tempFracs);
}

/** Computes the quantile of pairs of magnitudes found in each &Delta;t bin 
 *	of a &Delta;m&Delta;t plot.
 *
 * @param[in] deltaT The &Delta;t values of the &Delta;m&Delta;t plot.
 * @param[in] deltaM The corresponding &Delta;m values of the &Delta;m&Delta;t plot.
 * @param[in] binEdges A vector containing the (N+1) boundaries of the N &Delta;t bins 
 *	in which to calculate quantiles.
 * @param[out] quants A vector containing the quantiles within each bin.
 * @param[in] q The quantile to calculate.
 *
 * @pre @p deltaT.size() = @p deltaM.size()
 * @pre @p deltaT is sorted in ascending order
 * @pre @p binEdges is sorted in ascending order
 * @pre 0 < @p q < 1
 *
 * @pre @p deltaT does not contain any NaNs
 * @pre @p deltaM does not contain any NaNs
 * @pre @p binEdges does not contain any NaNs
 * 
 * @post @p quants.size() = @p binEdges.size() - 1
 * @post For all i &isin; [0, @p binEdges.size()-1], @p quants[i] contains the 
 *	<tt>q</tt>th quantile of @p deltaM, given @p deltaT &isin; 
 *	[@p binEdges[i], @p binEdges[i+1])
 *
 * @perform O(N log N) time, where N = @p deltaT.size()
 * @perfmore o(N log N - N log M) time, where N = @p deltaT.size() and M = @p binEdges.size()
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to store 
 *	the bin fractions.
 * @exception std::invalid_argument Thrown if @p q is not in (0, 1) or if 
 *	@p deltaT.size() &ne; @p deltaM.size()
 * @exception kpfutils::except::NotSorted Thrown if either @p deltaT or @p binEdges 
 *	is unsorted.
 *
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 */
void deltaMBinQuantile(const DoubleVec &deltaT, const DoubleVec &deltaM, 
		const DoubleVec &binEdges, DoubleVec &quants, double q) {
	using std::swap;
	
	if (deltaT.size() != deltaM.size()) {
		try {
			throw std::invalid_argument("Parameters 'deltaT' and 'deltaM' in deltaMBinQuantile() are not the same length (gave " 
			+ lexical_cast<string>(deltaT.size()) + " for deltaT and " 
			+ lexical_cast<string>(deltaM.size()) + " for deltaM)");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameters 'deltaT' and 'deltaM' in deltaMBinQuantile() are not the same length");
		}
	}
	if (q <= 0 || q >= 1) {
		throw std::invalid_argument("Quantile must be in (0, 1) (gave " 
			+ lexical_cast<std::string>(q) + ")");
	}
	if (!kpfutils::isSorted(deltaT.begin(), deltaT.end())) {
		throw kpfutils::except::NotSorted("deltaT is not sorted in deltaMBinQuantile()");
	}
	if (!kpfutils::isSorted(binEdges.begin(), binEdges.end())) {
		throw kpfutils::except::NotSorted("binEdges is not sorted in deltaMBinQuantile()");
	}
	
	// copy-and-swap
	DoubleVec tempQuants;
	tempQuants.reserve(binEdges.size()-1);
	
	for(DoubleVec::const_iterator curBin = binEdges.begin(); 
			curBin+1 != binEdges.end(); curBin++) {
		// Want deltam values for all pairs in the current bin
		// Start by finding the deltaT range in the interval [curBin, curBin+1)
		// Note: *_bound works only if deltaT is sorted
		DoubleVec::const_iterator deltaTStart = std::lower_bound(deltaT.begin(), 
				deltaT.end(), *curBin);
		DoubleVec::const_iterator deltaTEnd = std::lower_bound(deltaT.begin(), 
				deltaT.end(), *(curBin+1));
		DoubleVec::const_iterator deltaMStart = deltaM.begin() 
				+ (deltaTStart - deltaT.begin());
		DoubleVec::const_iterator deltaMEnd = deltaM.begin() 
				+ (deltaTEnd - deltaT.begin());
		
		if (deltaMStart != deltaMEnd) {
			double result = kpfutils::quantile(static_cast<DoubleVec::const_iterator>(deltaMStart), 
					deltaMEnd, q);
			tempQuants.push_back(result);
		} else {
			// The bin is empty
			tempQuants.push_back(std::numeric_limits<double>::quiet_NaN());
		}
	}

	// IMPORTANT: no exceptions beyond this point
	
	swap(quants, tempQuants);
}

}		// end kpftimes
