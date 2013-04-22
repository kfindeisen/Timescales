/** Performs unit testing of the function kpftimes::lsNormalEdf()
 * @file unit_lsNormalEdf.cpp
 * @author Krzysztof Findeisen
 * @date Created May 19, 2011
 * @date Last modified May 24, 2011
 */

#include <cstdio>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "../timescales.h"
#include "../utils.h"

using namespace kpftimes;

/** Carries out a generic test of the lsNormalEdf() function
 *
 * @param[in] times The times array to be passed to lsNormalEdf().
 * @param[in] freqs The frequency array to be passed to lsNormalEdf().
 * @param[in] nSims The number of simulations to request of lsNormalEdf().
 * @param[in] expectInvalid Whether the correct behavior is throwing a \
 *	invalid_argument exception. May be set simultaneously with 
 *	expectDomain, in which case catching either invalid_argument or 
 *	domain_behavior is considered correct behavior.
 * @param[in] expectDomain Whether the correct behavior is throwing a \
 *	domain_error exception. May be set simultaneously with 
 *	expectInvalid, in which case catching either invalid_argument or 
 *	domain_behavior is considered correct behavior.
 *
 * @return true if test passed, false if test failed
 * 
 */
bool testEdf(const DoubleVec& times, const DoubleVec& freqs, int nSims, 
		bool expectInvalid, bool expectDomain) throw () {
	const static size_t numRealRuns = 1000;

	DoubleVec powers, probs;
	try {
		lsNormalEdf(times, freqs, powers, probs, nSims);
		
		// If lsNormalEdf was supposed to throw an exception, but 
		//	instead ran to completion, that's an error
		if (expectInvalid || expectDomain) {
			return false;
		}
		// powers.size == probs.size == nSims
		if(nSims <= 0 || powers.size() != static_cast<size_t>(nSims) 
			|| probs.size() != static_cast<size_t>(nSims)) {
			return false;
		}
		// powers and probs should both be sorted
		if(!isSortedAsc(powers) || !isSortedAsc(probs)) {
			return false;
		}

		// powers matches output of lombScargle... at least, in a 
		//	statistical sense
		gsl_rng * calGen = gsl_rng_alloc(gsl_rng_taus2);
		/* Use a fixed seed for the reference distribution to avoid 
		 * problems from lsNormalEdf and lombScargle using identical seeds
		 * lsNormalEdf uses a variable seed in order to produce non-
		 * reproducible results
		 */
		gsl_rng_set(calGen, 101);
		DoubleVec highestPeak, randomFluxes, truePower;
		for (size_t i = 0; i < numRealRuns; i++) {
			randomFluxes.clear();
			truePower.clear();
			for(DoubleVec::const_iterator it = times.begin(); it != times.end(); 
					it++) {
				randomFluxes.push_back(gsl_ran_gaussian(calGen, 1.0));
			}
			try {
				lombScargle(times, randomFluxes, freqs, truePower);
			} catch (std::exception e) {
				return false;
			}
			highestPeak.push_back(*std::max_element(truePower.begin(), 
					truePower.end()));
		}
		gsl_rng_free(calGen);

		/** Compare highestPeak to powers in a statistical sense. For 
		 *	now, do this by visual inspection of plots. Eventually 
		 *	we want to implement a KS or (preferably) AD test to 
		 *	automate the comparison. GSL doesn't have either, 
		 *	believe it or not.
		 *
		 * @warning: missing automated oracle for EDFs
		 */
		 if (nSims >= 100) {
			 FILE *trueEdf = NULL, *ourEdf = NULL;
			 trueEdf = fopen("trueedf.txt", "w");
			 ourEdf = fopen( "ouredf.txt", "w");
			 if(trueEdf != NULL && ourEdf != NULL) {
			 	std::sort(highestPeak.begin(), highestPeak.end());
				for(size_t i = 0; i < highestPeak.size(); i++) {
					fprintf(trueEdf, "%.3f %.3f\n", highestPeak[i], 
							static_cast<double>(i+1) / 
							highestPeak.size());
				}
				for(size_t i = 0; i < powers.size(); i++) {
					fprintf(ourEdf, "%.3f %.3f\n", powers[i], probs[i]);
				}
			 }
			 // Clean up... but we don't know if file open succeeded
			 if ( ourEdf != NULL) {
			 	fclose( ourEdf);
			 	 ourEdf = NULL;
			 }
			 if (trueEdf != NULL) {
			 	fclose(trueEdf);
			 	trueEdf = NULL;
			 }
		 }
	} catch (std::invalid_argument e) {
		return expectInvalid;
	} catch (std::domain_error e) {
		return expectDomain;
	} catch (std::exception e) {
		return false;
	}
	
	return true;
}

/** Runs all test cases defined for lsNormalEdf(). Results are printed to stdout.
 */
void runLsNormalEdfTests() {
	// Set up common data structures
	DoubleVec times1, times2sort, times2unsort, times2dupe, 
			times100randsort, times100rand2sort, times100unifsort, 
			posFreq, zero2Freq, negFreq;

	// Define time grids
	times1      .push_back(30.0);
	times2sort  .push_back(30.0);
	times2sort  .push_back(40.0);
	times2unsort.push_back(40.0);
	times2unsort.push_back(30.0);
	times2dupe  .push_back(30.0);
	times2dupe  .push_back(30.0);
	gsl_rng * timeGen = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(timeGen, 42);
	for(size_t i = 0; i < 100; i++) {
		times100unifsort.push_back(0.452*(i+42));
		times100randsort.push_back(0.452*(100*gsl_rng_uniform(timeGen)+42));
		times100rand2sort.push_back(0.452*(100*gsl_rng_uniform(timeGen)+42));	// Make everything duplicated
		times100rand2sort.push_back(0.452*(100*gsl_rng_uniform(timeGen)+42));
	}
	gsl_rng_free(timeGen);
	DoubleVec times100randunsort = times100randsort;
	std::sort(times100randsort.begin(), times100randsort.end());
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
	
	std::vector<bool> testResults;
	
	/* @test A 1-element time series and nSims = 100. Expected behavior = throw 
	 *	invalid_argument
	 */
	testResults.push_back(testEdf(times1, posFreq, 1000, true, false));
	/* @test A 2-element time series, sorted with no duplicates, and nSims = 100. 
	 *	Expected behavior = matches result of running lombScargle 1000 times.
	 */
	testResults.push_back(testEdf(times2sort, posFreq, 1000, false, false));
	/* @test A 2-element time series, sorted with duplicates, and nSims = 100. 
	 *	Expected behavior = throw invalid_argument
	 */
	testResults.push_back(testEdf(times2dupe, posFreq, 1000, true, false));
	/* @test A 2-element time series, unsorted, and nSims = 100. Expected behavior 
	 *	= throw invalid_argument.
	 */
	testResults.push_back(testEdf(times2unsort, posFreq, 1000, true, false));
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and negative frequencies. Expected behavior = throw domain_error.
	 */
	testResults.push_back(testEdf(times100randsort, negFreq, 1000, false, true));
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and multiple zero frequencies. Expected behavior = matches 
	 *	result of running lombScargle 1000 times.
	 */
/* FAIL */
	testResults.push_back(testEdf(times100randsort, zero2Freq, 1000, false, false));
	/* @test A 100-element uniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 100. Expected behavior = matches result of 
	 *	running lombScargle 1000 times.
	 */
	testResults.push_back(testEdf(times100unifsort, posFreq, 1000, false, false));
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 100. Expected behavior = matches result of 
	 *	running lombScargle 1000 times.
	 */
	testResults.push_back(testEdf(times100randsort, posFreq, 1000, false, false));
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 0. Expected behavior = throw domain_error.
	 */
/* FAIL */
	testResults.push_back(testEdf(times100randsort, posFreq, 0, false, true));
	/* @test A 100-element nonuniformly sampled time series, sorted with no 
	 *	duplicates, and nSims = 1. Expected behavior = (probs = 1.0, powers = 
	 *	undefined)
	 */
	testResults.push_back(testEdf(times100randsort, posFreq, 1, false, false));
	/* @test A 100-element nonuniformly sampled time series, sorted with 
	 *	duplicates, and nSims = 100. Expected behavior = matches result of 
	 *	running lombScargle 1000 times.
	 */
	testResults.push_back(testEdf(times100rand2sort, posFreq, 1000, false, false));
	/* @test A 100-element nonuniformly sampled time series, unsorted. Expected 
	 *	behavior = throw invalid_argument.
	 */
	testResults.push_back(testEdf(times100randunsort, posFreq, 1000, true, false));
	
	// Print out the results
	
	size_t numPassed = 0;
	for(std::vector<bool>::const_iterator it = testResults.begin(); 
			it != testResults.end(); it++) {
		if (*it) {
			numPassed++;
		}
	}
	
	if (numPassed == testResults.size()) {
		printf("lsNormalEdf: All tests passed.\n");
		printf("\tWARNING: some tests require manual verification. Please examine trueedf.txt and ouredf.txt.\n");
	} else {
		printf("lsNormalEdf: %i/%i tests passed.\n", numPassed, testResults.size());
		for(std::vector<bool>::const_iterator it = testResults.begin(); 
				it != testResults.end(); it++) {
			if (*it) {
				printf("\tTest passed\n");
			} else {
				printf("\tTest failed\n");
			}
		}
	}
}
