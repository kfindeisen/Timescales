/** Brief demo of the Timescales capabilities. This program computes a single 
 *	periodogram and single autocorrelation function of an irregularly 
 *	sampled sinusoid.
 * @file example.cpp
 * @author Krzysztof Findeisen
 * @date Created February 14, 2011
 * @date Last modified April 15, 2011
 */

#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include "../timescales.h"

using namespace kpftimes;

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
 * 
 */
int main()
{
	const double PI   = 3.1415927;
	const double FREQ = 5.432;
	const double DELT = 1.0;
	
	FILE* hOutput = NULL;
	gsl_rng* generator = gsl_rng_alloc(gsl_rng_ranlxd2);
	// fixed seed, so program output will be the same from run to run
	gsl_rng_set(generator, 42);
	
	//--------------------------------------------------
	// Generate randomly sampled data
	size_t N = 100;
	DoubleVec t, x;
	for(size_t i = 0; i < N; i++) {
		t.push_back(DELT*gsl_rng_uniform(generator));
	}
	std::sort(t.begin(), t.end());
	for(DoubleVec::const_iterator it = t.begin(); it != t.end(); it++) {
		x.push_back(sin(2.0*PI*(*it)*FREQ));
	}
	
	hOutput = fopen("test.obs.txt", "w");
	for(size_t i = 0; i < t.size(); i++) {
		fprintf(hOutput, "%.4f, %.4f\n", t[i], x[i]);
	}
	fclose(hOutput);

	//--------------------------------------------------
	// Generate periodogram
	try {
		DoubleVec freq, power;
		freqGen(t, freq);			// Fill freq[] with a basic grid
		lombScargle(t, x, freq, power);		// Fill power[] with the periodogram

		hOutput = fopen("test.pgram.txt", "w");
		for(size_t i = 0; i < freq.size(); i++) {
			fprintf(hOutput, "%.4f, %.4f\n", freq[i], power[i]);
		}
		fclose(hOutput);
	} catch (std::logic_error e) {
		fprintf(stderr, "Periodogram generation failed: %60s\n", e.what());
	}
	
	//--------------------------------------------------
	// Generate autocorrelation function
	
	try {
		DoubleVec off, offHalf, acf;
		// Make a regular grid
		// Support for irregular grids will be added later
		for(double i = 0.0; i <= DELT; i+= 0.001) {
			off.push_back(i);
			if(i < DELT/2) {
				offHalf.push_back(i);
			}
		}
		
		autoCorr(t, x, offHalf, acf);

		hOutput = fopen("test.acf2.txt", "w");
		for(size_t i = 0; i < acf.size(); i++) {
			fprintf(hOutput, "%.4f, %.4f\n", offHalf[i], acf[i]);
		}
		fclose(hOutput);

		autoCorr(t, x, off, acf);

		hOutput = fopen("test.acf.txt", "w");
		for(size_t i = 0; i < acf.size(); i++) {
			fprintf(hOutput, "%.4f, %.4f\n", off[i], acf[i]);
		}
		fclose(hOutput);
	} catch (std::logic_error e) {
		fprintf(stderr, "Autocorrelation generation failed: %60s\n", e.what());
	}
	
	//--------------------------------------------------
	// Generate unwindowed window function
	
	try {
		DoubleVec off, acf, wf;
		// Make a regular grid
		// Support for irregular grids will be added later
		for(double i = 0.0; i <= DELT; i+= 0.001) {
			off.push_back(i);
		}

		autoCorr(t, x, off, acf);
		acWindow(t, off, wf);

		hOutput = fopen("test.win.txt", "w");
		for(size_t i = 0; i < wf.size(); i++) {
			fprintf(hOutput, "%.4f, %.4f\n", off[i], wf[i]);
		}
		fclose(hOutput);
		
		hOutput = fopen("test.unwinacf.txt", "w");
		for(size_t i = 0; i < off.size(); i++) {
			fprintf(hOutput, "%.4f, %.4f\n", off[i], acf[i]*wf[i]);
		}
		fclose(hOutput);
	} catch (std::logic_error e) {
		fprintf(stderr, "Autocorrelation generation failed: %60s\n", e.what());
	}
	
	return 0;
}
