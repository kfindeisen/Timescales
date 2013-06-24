/** Computes characteristic frequencies of a sampling cadence
 * @file timescales/specialfreqs.cpp
 * @author Krzysztof Findeisen
 * @date Created April 13, 2011
 * @date Last modified June 20, 2013
 */ 

#include <vector>
#include "../common/stats.tmp.h"
#include "timescales.h"

namespace kpftimes {

/** Returns the pseudo-Nyquist frequency for a grid of observations.
 * 
 * The pseudo-Nyquist frequency is defined as N/2T, where N is the number of 
 * observations and T is the length of the time interval covered by the data.
 * @ingroup util
 * 
 * @param[in] times	Times at which data were taken
 * 
 * @return The pseudo-Nyquist frequency, in the inverse of whatever units 
 *	times is in.
 *
 * @pre @p times contains at least two unique values
 * 
 * @perform O(N) time, where N = @p times.size()
 * 
 * @exception std::invalid_argument Thrown if preconditions violated.
 * 
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 *
 * @test Regular grid, length 1. Expected behavior: throws invalid_argument
 * @test Regular grid, length 2. Expected behavior: returns PNF = 1/(2*step)
 * @test Regular grid, length 100. Expected behavior: returns PNF = 1/(2*step)
 * 
 * @todo Come up with test cases for an irregular grid.
 */
double pseudoNyquistFreq(const DoubleVec &times) {
	// Delegate input validation to deltaT
	return 0.5 * times.size() / deltaT(times);
}

/** Returns the time interval covered by the data.
 * 
 * @param[in] times	Times at which data were taken
 *
 * @return The length of time between the earliest observation in times and 
 *	the latest observation in times, in whatever units @p times is in.
 *
 * @pre @p times contains at least two unique values
 * 
 * @perform O(N) time, where N = @p times.size()
 * 
 * @exception std::invalid_argument Thrown if preconditions violated.
 * 
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 * 
 * @test Regular grid, length 1. Expected behavior: throws invalid_argument
 * @test Regular grid, length 2. Expected behavior: returns step
 * @test Regular grid, length 100. Expected behavior: returns 99*step
 * @test Irregular grid, 2 values randomly chosen from [0, 1). Expected behavior: returns max_element()-min_element()
 * @test Irregular grid, 100 values randomly chosen from [0, 1). Expected behavior: returns max_element()-min_element()
 */
double deltaT(const DoubleVec &times) {
	if (times.size() < 2) {
		throw std::invalid_argument("Parameter 'times' in deltaT() contains fewer than 2 observations");
	}
	
	// In C++11, this entire body can be replaced with a call to std::minmax_element()

	// Scanning the array to verify that it's sorted would take just as 
	//	long as scanning it for min and max
	double tMin = times.front();
	double tMax = tMin;
	
	for (DoubleVec::const_iterator i = times.begin(); i != times.end(); i++) {
		if (*i < tMin) {
			tMin = *i;
		}
		if (*i > tMax) {
			tMax = *i;
		}
	}
	
	if (tMax <= tMin) {
		throw std::invalid_argument("Parameter 'times' in deltaT() contains only one unique value");
	}
	else {
		return (tMax - tMin);
	}
}

/** Returns the highest frequency that can be probed by the data. This is 
 * defined as 1/2dt, where dt > 0 is the @b smallest time interval between 
 * any two observations.
 * 
 * @param[in] times	Times at which data were taken
 *
 * @return The highest meaningful frequency, in the inverse of whatever units 
 *	@p times is in.
 *
 * @pre @p times contains at least two unique values
 * @pre @p times is sorted in ascending order
 * 
 * @perform O(N) time, where N = @p times.size()
 * 
 * @exception std::invalid_argument Thrown if preconditions violated.
 * 
 * @exceptsafe The function arguments are unchanged in the event 
 *	of an exception.
 * 
 * @test Regular grid, length 1. Expected behavior: throws invalid_argument
 * @test Regular grid, length 2. Expected behavior: returns PNF = 1/(2*step)
 * @test Regular grid, length 100. Expected behavior: returns PNF = 1/(2*step)
 * @todo How to test this for irregular grids?
 */
double maxFreq(const DoubleVec &times) {
	if (times.size() < 2) {
		throw std::invalid_argument("Parameter 'times' in maxFreq() contains fewer than 2 observations");
	}

	// Test for sort in O(N)
	// Faster than sorting, O(N log N), or unsorted test, O(N^2)
	if(!kpfutils::isSorted(times.begin(), times.end())) {
		throw std::invalid_argument("Parameter 'times' in maxFreq() is unsorted");
	}

	// Look for the smallest interval
	DoubleVec::const_iterator t1 = times.begin();
	DoubleVec::const_iterator t2 = t1 + 1;
	double minDeltaT = 0.0;
	
	while(t2 != times.end()) {
		if (*t2 > *t1 && (*t2 - *t1 < minDeltaT || minDeltaT == 0.0)) {
			minDeltaT = *t2 - *t1;
		}
		t1++;
		t2++;
	}
	
	// Report the results
	if (minDeltaT == 0.0) {
		throw std::invalid_argument("Parameter 'times' in maxFreq() contains only one unique value");
	} else {
		return 0.5 / minDeltaT;
	}
}

}		// end kpftimes
