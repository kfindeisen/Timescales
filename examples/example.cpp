/** Brief demo of the Timescales capabilities.
 * @file examples/example.cpp
 * @author Krzysztof Findeisen
 * @date Created February 14, 2011
 * @date Last modified November 27, 2013
 *
 * This program computes a single periodogram and single autocorrelation 
 * function of an irregularly sampled sinusoid.
 */

#include <algorithm>
#include <stdexcept>
#include <string>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <boost/math/constants/constants.hpp>
#include <boost/smart_ptr.hpp>
#include "../timescales.h"

using namespace kpftimes;
using boost::shared_ptr;

/** Generates a randomly sampled sine curve
 * 
 * @param[out] time A grid of random time samplings
 * @param[out] signal A sinusoidal function evaluated at @p time
 */
void makeSine(DoubleVec& time, DoubleVec& signal);

/** Prints a two-column file
 *
 * @param[in] fileName The name of a text file to open
 * @param[in] col1 The first column of data to print
 * @param[in] col2 The second column of data to print
 *
 * @pre @p col1.size() = @p col2.size()
 *
 * @post creates or overwrites @p fileName, printing corresponding elements 
 *	of @p col1 and @p col2 in CSV format
 */
void writeTwoColumns(const std::string& fileName, 
		const DoubleVec& col1, const DoubleVec& col2);

/** Test driver for a simple user example.
 *
 * @return Status code 0
 * 
 * @post creates or overwrites the file "test.obs.txt" containing a 
 *	sinusoidal light curve in CSV format
 * @post creates or overwrites the file "test.pgram.txt" containing the 
 *	periodogram of test.obs.txt in CSV format
 * @post creates or overwrites the files "test.acf.txt" and "test.acf2.txt" 
 *	containing the autocorrelation function of test.obs.txt in CSV format
 * @post creates or overwrites the file "test.unwinacf.txt" containing the 
 *	ACF, not corrected for the window function, of test.obs.txt in CSV format
 * @post creates or overwrites the file "test.win.txt" containing the 
 *	ACF window function of test.obs.txt in CSV format
 */
int main() {
	//--------------------------------------------------
	// Generate randomly sampled data
	DoubleVec t, x;
	makeSine(t, x);
	double tMax = *std::max_element(t.begin(), t.end());
	
	writeTwoColumns("test.obs.txt", t, x);

	//--------------------------------------------------
	// Generate periodogram
	try {
		DoubleVec freq, power;
		freqGen(t, freq);			// Fill freq[] with a basic grid
		lombScargle(t, x, freq, power);		// Fill power[] with the periodogram

		writeTwoColumns("test.pgram.txt", freq, power);
	} catch (const std::exception& e) {
		fprintf(stderr, "Periodogram generation failed: %s\n", e.what());
	}
	
	//--------------------------------------------------
	// Generate autocorrelation function
	
	try {
		DoubleVec off, offHalf, acf;
		// Make a regular grid
		// Support for irregular grids will be added later
		for(double i = 0.0; i <= tMax; i+= 0.001) {
			off.push_back(i);
			if(i < tMax/2) {
				offHalf.push_back(i);
			}
		}
		

		autoCorr(t, x, offHalf, acf);
		writeTwoColumns("test.acf2.txt", offHalf, acf);

		autoCorr(t, x, off, acf);
		writeTwoColumns("test.acf.txt", off, acf);
	} catch (const std::exception& e) {
		fprintf(stderr, "Autocorrelation generation failed: %s\n", e.what());
	}
	
	//--------------------------------------------------
	// Generate unwindowed window function
	
	try {
		DoubleVec off, acf, wf;
		// Make a regular grid
		// Support for irregular grids will be added later
		for(double i = 0.0; i <= tMax; i+= 0.001) {
			off.push_back(i);
		}


		autoCorr(t, x, off, acf);
		acWindow(t, off, wf);
		writeTwoColumns("test.win.txt", off, wf);

		DoubleVec acfWf;
		for(size_t i = 0; i < acf.size(); i++) {
			acfWf.push_back(acf[i]*wf[i]);
		}		
		writeTwoColumns("test.unwinacf.txt", off, acfWf);
	} catch (const std::exception& e) {
		fprintf(stderr, "ACF window function generation failed: %s\n", e.what());
	}
	
	return 0;
}

void makeSine(DoubleVec& time, DoubleVec& signal) {
	using boost::math::double_constants::two_pi;
	const double FREQ = 5.432;
	const double DELT = 1.0;
	const size_t N = 100;

	shared_ptr<gsl_rng> generator(gsl_rng_alloc(gsl_rng_ranlxd2), &gsl_rng_free);
	// fixed seed, so program output will be the same from run to run
	gsl_rng_set(generator.get(), 42);
	
	time.clear();
	signal.clear();
	
	for(size_t i = 0; i < N; i++) {
		time.push_back(DELT*gsl_rng_uniform(generator.get()));
	}
	std::sort(time.begin(), time.end());
	
	for(DoubleVec::const_iterator it = time.begin(); it != time.end(); it++) {
		signal.push_back(sin(two_pi*(*it)*FREQ));
	}
}

void writeTwoColumns(const std::string& fileName, 
		const DoubleVec& col1, const DoubleVec& col2) {
	shared_ptr<FILE> hOutput(fopen(fileName.c_str(), "w"), &fclose);
	
	for(size_t i = 0; i < col1.size(); i++) {
		fprintf(hOutput.get(), "%.4f, %.4f\n", col1[i], col2[i]);
	}
}
