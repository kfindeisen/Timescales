/** Computes delta-T pair plots, including self-correlation functions and gap plots
 * @file pairwise.cpp
 * @author Krzysztof Findeisen
 * @date Created July 24, 2011
 * @date Last modified July 24, 2011
 */ 

#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>
#include "timescales.h"

#if USELFP
#include <lfp/lfp.h>
#endif

#ifndef PI
#define PI 3.1415927
#endif

using namespace std;

/** Calculates all Delta-m Delta-T pairs. deltaT is returned sorted in 
 *	ascending order. Pairs of observations with the same time separation may be returned 
 *	in any order.
 * 
 * @param[in] times	Times at which data were taken
 * @param[in] fluxes	Flux measurements of a source
 * @param[out] deltaT	A list of the time intervals between all pairs of sources.
 * @param[out] deltaM	A list of the magnitude difference between each pair in deltaT.
 *
 * @pre times contains at least two unique values
 * @pre times is sorted in ascending order
 * @pre fluxes is of the same length as times
 * @pre fluxes[i] is the flux of the source at times[i], for all i
 * @post if N == length(times), length(deltaT) == length(deltaM) == N(N-1)/2
 * @exception invalid_argument Thrown if the preconditions on times or 
 *	length(fluxes) are violated.
 *
 * @perform O(times.size()^2) time
 *
 */
void kpftimes::dmdt(const DoubleVec &times, const DoubleVec &fluxes, 
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
		throw invalid_argument("times contains only one unique date");
	} else if (!sortedTimes) {
		throw invalid_argument("times is not sorted in ascending order");
	} else if (fluxes.size() != nTimes) {
		throw invalid_argument("times and fluxes are not the same length");
	}

	// Start making the output
	deltaT.clear();
	deltaT.reserve(nTimes*(nTimes-1)/2);
	deltaM.clear();
	deltaM.reserve(nTimes*(nTimes-1)/2);
	
	// Needed for sorting by deltaT
	vector<pair<double, double> > sortableVec;
	
	for(size_t i = 0; i < nTimes; i++) {
		for(size_t j = i+1; j < nTimes; j++) {
//			deltaT.push_back(abs(times[i]-times[j]));
//			deltaM.push_back(abs(fluxes[i]-fluxes[j]));

			// Time must be first element so that the pairs get sorted properly
			sortableVec.push_back(pair<double, double>(abs(times[i]-times[j]), 
					abs(fluxes[i]-fluxes[j]) ));
		}
	}
	
	// Now we just need to sort it
	// I don't have a parallel sort function yet
	// For now, do it the clumsy way, using pair<>
	sort(sortableVec.begin(), sortableVec.end());
	
	for(vector<pair<double, double> >::const_iterator it = sortableVec.begin(); 
			it != sortableVec.end(); it++) {
		deltaT.push_back(it->first);
		deltaM.push_back(it->second);
	}
}
