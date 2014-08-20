/** Performs unit testing of the function kpftimes::lsNormalEdf()
 * @file timescales/tests/unit_lsNormalEdf.cpp
 * @author Krzysztof Findeisen
 * @date Created May 19, 2011
 * @date Last modified November 19, 2013
 */

/* Copyright 2014, California Institute of Technology.
 *
 * This file is part of the Timescales library.
 * 
 * The Timescales library is free software: you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version, subject to the following 
 * exception added under Section 7 of the License:
 *	* Neither the name of the copyright holder nor the names of its contributors 
 *	  may be used to endorse or promote products derived from this software 
 *	  without specific prior written permission.
 * 
 * The Timescales library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with the Timescales library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../common/warnflags.h"

// Boost.Test uses C-style casts and non-virtual destructors
#ifdef GNUC_COARSEWARN
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Weffc++"
#endif

// Boost.Test uses C-style casts and non-virtual destructors
#ifdef GNUC_FINEWARN
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/test/unit_test.hpp>

// Re-enable all compiler warnings
#ifdef GNUC_FINEWARN
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <boost/smart_ptr.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "../../common/alloc.tmp.h"
#include "../../common/cerror.h"
#include "../../common/stats.tmp.h"
#include "../timescales.h"

namespace kpftimes { namespace test {

using boost::shared_ptr;
using kpfutils::checkAlloc;

/** Data common to the test cases.
 *
 * Contains generic time and frequency grids
 */
class EdfData {
public: 
	/** Defines the data for each test case.
	 *
	 * @exception std::bad_alloc Thrown if there is not enough memory to 
	 *	store the testing data.
	 *
	 * @exceptsafe Object construction is atomic.
	 */
	EdfData(): times1(), times2sort(), times2unsort(), times2dupe(), 
			times100randsort(), times100randunsort(), times100rand2sort(), times100unifsort(), 
			posFreq(), zero2Freq(), negFreq() {
	
		// Define time grids
		times1      .push_back(30.0);
		times2sort  .push_back(30.0);
		times2sort  .push_back(40.0);
		times2unsort.push_back(40.0);
		times2unsort.push_back(30.0);
		times2dupe  .push_back(30.0);
		times2dupe  .push_back(30.0);
		
		shared_ptr<gsl_rng> timeGen(checkAlloc(gsl_rng_alloc(gsl_rng_taus2)), 
			&gsl_rng_free);
		gsl_rng_set(timeGen.get(), 42);
		
		for(size_t i = 0; i < 100; i++) {
			times100unifsort.push_back(0.452*(i+42));
			
			double randTime = 0.452*(100*gsl_rng_uniform(timeGen.get())+42);
			times100randsort.push_back(randTime);
			
			times100rand2sort.push_back(randTime);	// Make everything duplicated
			times100rand2sort.push_back(randTime);
		}
		
		times100randunsort = times100randsort;
		std::sort(times100randsort .begin(), times100randsort .end());
		std::sort(times100rand2sort.begin(), times100rand2sort.end());
		
		// Define frequency grids
		for(double i = 0.0; i < 2.5; i+=0.01) {
			if (i != 0.0) {
				posFreq.push_back(i);
				negFreq.push_back(i);
			}
			zero2Freq.push_back(i);
			zero2Freq.push_back(i);			// Make all frequencies duplicated
		}
		for(double i = 0.01; i < 0.1; i+=0.01) {
			negFreq.push_back(-i);			// Who says negFreq needs to be sorted?
		}
	}
	
	virtual ~EdfData() {
	}
	
	/** Grid with a single time
	 */
	DoubleVec times1;
	/** Grid with two times in ascending order
	 */
	DoubleVec times2sort;
	/** Grid with two times in descending order
	 */
	DoubleVec times2unsort;
	/** Grid with two identical times
	 */
	DoubleVec times2dupe;
	
	/** Grid with 100 random times in ascending order
	 */
	DoubleVec times100randsort;
	/** Grid with 100 random times, unsorted
	 */
	DoubleVec times100randunsort;
	/** Grid with 100 random times in ascending order, with duplicates
	 */
	DoubleVec times100rand2sort; 
	/** Grid with 100 evenly spaced times in ascending order
	 */
	DoubleVec times100unifsort;
	
	/** Grid with only positive frequencies, in ascending order
	 */
	DoubleVec posFreq;
	/** Grid with zero and positive frequencies, in ascending order, 
	 *	with duplicates
	 */
	DoubleVec zero2Freq;
	/** Grid with negative frequencies, unsorted
	 */
	DoubleVec negFreq;
};

/** Tests the specification of lsNormalEdf() for specific inputs
 *
 * @param[in] times The times array to be passed to lsNormalEdf().
 * @param[in] freqs The frequency array to be passed to lsNormalEdf().
 * @param[in] nSims The number of simulations to request of lsNormalEdf().
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to run the test.
 * @exception std::exception Thrown if incorrect arguments are passed to lsNormalEdf()
 *
 * @exceptsafe Program is in a valid state in the event of an exception.
 *
 * @internal @note Program opens and writes to files, preventing a stronger guarantee
 */
void testEdf(const DoubleVec& times, const DoubleVec& freqs, size_t nSims) {
	const static size_t numRealRuns = 1000;

	DoubleVec powers, probs;
	// Do not wrap exceptions here -- let them be tested by the caller
	lsNormalEdf(times, freqs, powers, probs, nSims);
	
	// Make sure that, if lsNormalEdf did not throw, then the data is valid
	BOOST_REQUIRE(nSims > 0);

	BOOST_REQUIRE_EQUAL(powers.size(), nSims);
	BOOST_REQUIRE_EQUAL(probs .size(), nSims);
	
	// powers and probs should both be sorted
	BOOST_REQUIRE(kpfutils::isSorted(powers.begin(), powers.end()));
	BOOST_REQUIRE(kpfutils::isSorted(probs .begin(), probs .end()));
	
	// powers matches output of lombScargle... at least, in a 
	//	statistical sense
	shared_ptr<gsl_rng> calGen(checkAlloc(gsl_rng_alloc(gsl_rng_mt19937)), &gsl_rng_free);
	gsl_rng_set(calGen.get(), 101);
	
	DoubleVec highestPeak, randomFluxes, truePower;
	for (size_t i = 0; i < numRealRuns; i++) {
		randomFluxes.clear();
		truePower.clear();
		for(DoubleVec::const_iterator it = times.begin(); it != times.end(); 
				it++) {
			randomFluxes.push_back(gsl_ran_gaussian(calGen.get(), 1.0));
		}

		BOOST_REQUIRE_NO_THROW(lombScargle(times, randomFluxes, freqs, truePower) );

		highestPeak.push_back(*std::max_element(truePower.begin(), 
				truePower.end()));
	}

	/** Compare highestPeak to powers in a statistical sense. For 
	 *	now, do this by visual inspection of plots. Eventually 
	 *	we want to implement a KS or (preferably) AD test to 
	 *	automate the comparison. GSL doesn't have either.
	 *
	 * @todo: add automated oracle for EDFs
	 */
	 if (nSims >= 100) {
	 	shared_ptr<FILE> trueEdf = kpfutils::fileCheckOpen("trueedf.txt", "w");
	 	shared_ptr<FILE>  ourEdf = kpfutils::fileCheckOpen( "ouredf.txt", "w");
		
		std::sort(highestPeak.begin(), highestPeak.end());
		for(size_t i = 0; i < highestPeak.size(); i++) {
			if (fprintf(trueEdf.get(), "%.3f %.3f\n", highestPeak[i], 
					static_cast<double>(i+1) / highestPeak.size()) < 0) {
				kpfutils::fileError(trueEdf.get(), 
					"Could not write to trueedf.txt: ");
			}
		}
		for(size_t i = 0; i < powers.size(); i++) {
			if (fprintf(ourEdf.get(), "%.3f %.3f\n", powers[i], probs[i]) < 0) {
				kpfutils::fileError(ourEdf.get(), 
					"Could not write to ouredf.txt: ");
			}
		}
	}
}


/** Test cases for lsNormalEdf()
 * @class BoostTest::test_edf
 */
BOOST_FIXTURE_TEST_SUITE(test_edf, EdfData)

/** Tests whether lsNormalEdf() matches behavior of lombScargle() on white 
 *	Gaussian noise
 *
 * @exceptsafe Does not throw exceptions.
 */
BOOST_AUTO_TEST_CASE(normal) {
	////////////////////////////////////////////////////
	// Minimal-length time series

	/* @test A 1-element time series and nSims = 100. Expected behavior = throw 
	 *	invalid_argument
	 */
	BOOST_CHECK_THROW(testEdf(times1, posFreq, 1000), std::invalid_argument);
	/* @test A 2-element time series, sorted with no duplicates, and nSims = 100. 
	 *	Expected behavior = matches result of running lombScargle 1000 times.
	 */
	BOOST_CHECK_NO_THROW(testEdf(times2sort, posFreq, 1000));
	/* @test A 2-element time series, sorted with duplicates, and nSims = 100. 
	 *	Expected behavior = throw invalid_argument
	 */
	BOOST_CHECK_THROW(testEdf(times2dupe, posFreq, 1000), std::invalid_argument);
	/* @test A 2-element time series, unsorted, and nSims = 100. Expected behavior 
	 *	= throw invalid_argument.
	 */
	BOOST_CHECK_THROW(testEdf(times2unsort, posFreq, 1000), std::invalid_argument);

	////////////////////////////////////////////////////
	// Invalid input

	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 0. Expected behavior = throw invalid_argument.
	 */
	BOOST_CHECK_THROW(testEdf(times100randsort, posFreq, 0), std::invalid_argument);
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 100. Expected behavior = matches result of 
	 *	running lombScargle 1000 times.
	 */
	BOOST_CHECK_NO_THROW(testEdf(times100randsort, posFreq, 1000));
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 1. Expected behavior = (probs = 1.0, powers = 
	 *	undefined)
	 */
	BOOST_CHECK_NO_THROW(testEdf(times100randsort, posFreq, 1));

	////////////////////////////////////////////////////
	// Frequency grid testing

	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and negative frequencies. Expected behavior = throw invalid_argument.
	 */
	BOOST_CHECK_THROW(testEdf(times100randsort, negFreq, 1000), std::invalid_argument);
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and multiple zero frequencies. Expected behavior = matches 
	 *	result of running lombScargle 1000 times.
	 */
	BOOST_CHECK_NO_THROW(testEdf(times100randsort, zero2Freq, 1000));

	////////////////////////////////////////////////////
	// Time grid testing

	/* @test A 100-element uniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 100. Expected behavior = matches result of 
	 *	running lombScargle 1000 times.
	 */
	BOOST_CHECK_NO_THROW(testEdf(times100unifsort, posFreq, 1000));
	/* @test A 100-element nonuniformly sampled time series, sorted with 
	 *	duplicates, and nSims = 100. Expected behavior = matches result of 
	 *	running lombScargle 1000 times.
	 */
	BOOST_CHECK_NO_THROW(testEdf(times100rand2sort, posFreq, 1000));
	/* @test A 100-element nonuniformly sampled time series, unsorted. Expected 
	 *	behavior = throw invalid_argument.
	 */
	BOOST_CHECK_THROW(testEdf(times100randunsort, posFreq, 1000), std::invalid_argument);
	
	BOOST_MESSAGE("WARNING: some tests require manual verification. Please examine trueedf.txt and ouredf.txt.");
}

BOOST_AUTO_TEST_SUITE_END()

}}		// end kpftimes::test
