/** Computes delta-T pair plots, including self-correlation functions and gap plots
 * @file timescales/pairwise.cpp
 * @author Krzysztof Findeisen
 * @date Created July 24, 2011
 * @date Last modified June 20, 2013
 */ 

#include <algorithm>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include "timescales.h"

namespace kpftimes {

using std::string;
using boost::lexical_cast;
using boost::math::double_constants::two_pi;

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
 * @exception std::invalid_argument Thrown if the preconditions on 
 *	@p times or @p mags.size() are violated.
 * @exception std::bad_alloc Thrown if there is not enough memory to compute 
 *	the &Delta;m&Delta;t plot
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
		throw std::invalid_argument("Parameter 'times' in lombScargle() contains only one unique date");
	} else if (!sortedTimes) {
		throw std::invalid_argument("Parameter 'times' in lombScargle() is not sorted in ascending order");
	} else if (mags.size() != nTimes) {
		try {
			throw std::invalid_argument("Parameters 'times' and 'mags' in lombScargle() are not the same length (gave " 
			+ lexical_cast<string>(nTimes) + " for times and " 
			+ lexical_cast<string>(mags.size()) + " for mags)");
		} catch (const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Parameters 'times' and 'mags' in lombScargle() are not the same length");
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

}		// end kpftimes
