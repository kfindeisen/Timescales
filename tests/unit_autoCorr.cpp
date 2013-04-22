/** Performs unit testing of the function kpftimes::autoCorr()
 * @file unit_autoCorr.cpp
 * @author Krzysztof Findeisen
 * @date Created July 15, 2011
 * @date Last modified July 15, 2011
 */

#include <stdexcept>
#include <string>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "../timescales.h"

using namespace kpftimes;

/** Carries out a generic test of the autoCorr() function
 *
 * @param[in] times The times array to be passed to autoCorr().
 * @param[in] fluxes The flux array to be passed to autoCorr().
 * @param[in] offsets The offset array to be passed to autoCorr().
 *
 * @param[in] expectInvalid Whether the correct behavior is throwing a \
 *	invalid_argument exception.
 *
 * @param[in] outFile The name of the file to which to print the ACF. No printing if empty string.
 *
 * @return true if test passed, false if test failed
 * 
 */
bool testAcf(const DoubleVec& times, const DoubleVec& fluxes, const DoubleVec& offsets, 
		bool expectInvalid, const std::string &outFile) throw () {
	DoubleVec acf;
	try {
		autoCorr(times, fluxes, offsets, acf);
		
		// If lsNormalEdf was supposed to throw an exception, but 
		//	instead ran to completion, that's an error
		if (expectInvalid) {
			return false;
		}
		
		/** @warning: missing automated oracle for ACFs
		 */
		 if (outFile.size() > 0) {
			 FILE *hAcfOutput = NULL;
			 hAcfOutput = fopen(outFile.c_str(), "w");
			 if(hAcfOutput != NULL) {
			 	for(size_t i = 0; i < acf.size(); i++) {
					fprintf(hAcfOutput, "%.3f, %.3f\n", offsets[i], acf[i]);
				}
				 // Clean up... 
			 	fclose(hAcfOutput);
			 	hAcfOutput = NULL;
			 }
		 }
	} catch (std::invalid_argument e) {
		return expectInvalid;
	} catch (std::exception e) {
		return false;
	}
	
	return true;
}

/** Runs all test cases defined for autoCorr(). Results are printed to stdout.
 */
void runAutoCorrTests() {
	const static double WAVEFREQ  = 4.52;
	const static double WAVEPHASE = 0.45;
	const static double OBSFREQ   = 0.01;
	const static double OBSPHASE  = 22;
	const static double PI        = 3.1415927;

	// Set up common data structures
	DoubleVec times100randsort, times100rand2sort, times100unifsort;
	DoubleVec flux100noise0, flux100noise001, flux100noise01, flux100noise05, 
		flux100noise1, flux100noise2;
	DoubleVec offsets100;

	// Define time grids
	gsl_rng * timeGen = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(timeGen, 42);
	for(size_t i = 0; i < 100; i++) {
		times100unifsort.push_back( OBSFREQ*(i+OBSPHASE));
		times100randsort.push_back( OBSFREQ*(100*gsl_rng_uniform(timeGen)+OBSPHASE));
		times100rand2sort.push_back(OBSFREQ*(100*gsl_rng_uniform(timeGen)+OBSPHASE));	// Make everything duplicated
		times100rand2sort.push_back(OBSFREQ*(100*gsl_rng_uniform(timeGen)+OBSPHASE));
	}
	DoubleVec times100randunsort = times100randsort;
	std::sort(times100randsort.begin(), times100randsort.end());
	std::sort(times100rand2sort.begin(), times100rand2sort.end());
	
	// Define flux grids
	for(DoubleVec::const_iterator it = times100randsort.begin(); 
			it != times100randsort.end(); it++) {
		double noise = gsl_ran_gaussian(timeGen, 1.0);
		
		flux100noise0  .push_back(sin(2.0*PI*(*it)*WAVEFREQ+WAVEPHASE));
		flux100noise001.push_back(sin(2.0*PI*(*it)*WAVEFREQ+WAVEPHASE) + 0.01*noise);
		flux100noise01 .push_back(sin(2.0*PI*(*it)*WAVEFREQ+WAVEPHASE) + 0.1 *noise);
		flux100noise05 .push_back(sin(2.0*PI*(*it)*WAVEFREQ+WAVEPHASE) + 0.5 *noise);
		flux100noise1  .push_back(sin(2.0*PI*(*it)*WAVEFREQ+WAVEPHASE) + 1.0 *noise);
		flux100noise2  .push_back(sin(2.0*PI*(*it)*WAVEFREQ+WAVEPHASE) + 2.0 *noise);
	}

	{
		FILE *hSine0Output   = fopen("flux.sine0.csv",   "w");
		FILE *hSine001Output = fopen("flux.sine001.csv", "w");
		FILE *hSine01Output  = fopen("flux.sine01.csv",  "w");
		FILE *hSine05Output  = fopen("flux.sine05.csv",  "w");
		FILE *hSine1Output   = fopen("flux.sine1.csv",   "w");
		FILE *hSine2Output   = fopen("flux.sine2.csv",   "w");
		for(size_t i = 0; i < times100randsort.size(); i++) {
			fprintf(hSine0Output, "%.3f, %.3f\n",   times100randsort[i], 
					flux100noise0[i]);
			fprintf(hSine001Output, "%.3f, %.3f\n", times100randsort[i], 
					flux100noise001[i]);
			fprintf(hSine01Output, "%.3f, %.3f\n",  times100randsort[i], 
					flux100noise01[i]);
			fprintf(hSine05Output, "%.3f, %.3f\n",  times100randsort[i], 
					flux100noise05[i]);
			fprintf(hSine1Output, "%.3f, %.3f\n",   times100randsort[i], 
					flux100noise1[i]);
			fprintf(hSine2Output, "%.3f, %.3f\n",   times100randsort[i], 
					flux100noise2[i]);
		}
		 // Clean up... 
		fclose(hSine0Output);
		fclose(hSine001Output);
		fclose(hSine01Output);
		fclose(hSine05Output);
		fclose(hSine1Output);
		fclose(hSine2Output);
	}
	
	// Define offset grids
	for(double i = 0.0; i < 100*OBSFREQ; i+=0.001) {
		offsets100.push_back(i);
	}

	gsl_rng_free(timeGen);
	
	std::vector<bool> testResults;
	
	/* @test A 100 element time series with no noise
	 */
	testResults.push_back(testAcf(times100randsort, flux100noise0, offsets100, 
		false, "acf.noise0.csv"));
	/* @test A 100 element time series with SNR 100
	 */
	testResults.push_back(testAcf(times100randsort, flux100noise001, offsets100, 
		false, "acf.noise001.csv"));
	/* @test A 100 element time series with SNR 10
	 */
	testResults.push_back(testAcf(times100randsort, flux100noise01, offsets100, 
		false, "acf.noise01.csv"));
	/* @test A 100 element time series with SNR 2
	 */
	testResults.push_back(testAcf(times100randsort, flux100noise05, offsets100, 
		false, "acf.noise05.csv"));
	/* @test A 100 element time series with SNR 1
	 */
	testResults.push_back(testAcf(times100randsort, flux100noise1, offsets100, 
		false, "acf.noise1.csv"));
	/* @test A 100 element time series with SNR 0.5
	 */
	testResults.push_back(testAcf(times100randsort, flux100noise2, offsets100, 
		false, "acf.noise2.csv"));
	
	
	// Print out the results
	size_t numPassed = 0;
	for(std::vector<bool>::const_iterator it = testResults.begin(); 
			it != testResults.end(); it++) {
		if (*it) {
			numPassed++;
		}
	}
	
	if (numPassed == testResults.size()) {
		printf("autoCorr: All tests passed.\n");
		printf("\tWARNING: some tests require manual verification. Please examine acf.*.csv.\n");
	} else {
		printf("autoCorr: %i/%i tests passed.\n", numPassed, testResults.size());
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
