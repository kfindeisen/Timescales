/** Test unit for statistics template routines in utils.h
 * @file unit_stats.cpp
 * @author Krzysztof Findeisen
 * @date Created July 20, 2011
 * @date Last modified July 20, 2011
 */

#include <boost/test/unit_test.hpp>
#include "../utils.h"

#include <list>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>

// Data common to the test cases
// Format: Three trivial lists
//	TEST_COUNT realizations of random collections (array, list, vector) 
//	of TEST_LEN numbers each
class StatsData {
public: 
	const static size_t TEST_LEN       = 100;
	const static size_t TEST_COUNT     =  10;
	const static double TEST_TOLERANCE = 1e-10;

	StatsData() : emptyList(), oneList(1, 42), twoList() {
		twoList.push_back(-10);
		twoList.push_back(27);
	
		// Data come from a random number generator
		gsl_rng * testRng = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(testRng, 42);
		
		// Define the arrays, then the STL versions
		for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
			intArray[nTest] = new int   [TEST_LEN];
			dblArray[nTest] = new double[TEST_LEN];
			dblVec[nTest].reserve(TEST_LEN);
			for (size_t i = 0; i < TEST_LEN; i++) {
				// Large sigma to allow interesting truncation to int
				intArray[nTest][i] = static_cast<int>(gsl_ran_gaussian(testRng, 1000.0));
				intList[nTest].push_back(intArray[nTest][i]);
				
				dblArray[nTest][i] = gsl_ran_ugaussian(testRng);
				dblList[nTest].push_back(dblArray[nTest][i]);
				dblVec[nTest].push_back(dblArray[nTest][i]);
			}
		}
		
		// Clean up
		gsl_rng_free(testRng);
	}
	
	~StatsData() {
		for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
			if (intArray[nTest] != NULL) {
				delete [] intArray[nTest];
			}
			if (dblArray[nTest] != NULL) {
				delete [] dblArray[nTest];
			}
		}
	}

	// Empty list of ints
	std::list<int> emptyList;
	// Length-1 list of ints
	std::list<int> oneList;
	// Length-2 list of ints
	std::list<int> twoList;
	// Length-100 collections of integers
	// Need TEST_COUNT of them to avoid statistical flukes
	int* intArray[TEST_COUNT];
	std::list<int> intList[TEST_COUNT];
	// Length-100 collections of floating point numbers
	double *dblArray[TEST_COUNT];
	std::list<double> dblList[TEST_COUNT];
	std::vector<double> dblVec[TEST_COUNT];
};

BOOST_FIXTURE_TEST_SUITE(test_stats, StatsData)

BOOST_AUTO_TEST_CASE(mean)
{
	//@test List of ints, length 0. Expected behavior: throw invalid_argument.
	BOOST_CHECK_THROW(kpftimes::mean(emptyList.begin(), emptyList.end()), 
			std::invalid_argument);
	
	//@test List of ints, length 1. Expected behavior: return list[0]
	BOOST_CHECK_EQUAL(kpftimes::mean(oneList.begin(), oneList.end()), oneList.front());
	
	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		int trueMean = static_cast<int>(gsl_stats_int_mean(intArray[nTest], 1, 100));
		//@test List of ints, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_mean in 10 out of 10 trials.
		BOOST_CHECK_EQUAL(kpftimes::mean(intList[nTest].begin(), intList[nTest].end()), 
				trueMean);
	}

	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		double trueMean = gsl_stats_mean(dblArray[nTest], 1, 100);
		//@test List of doubles, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_mean to within 1e-10 in 10 out of 10 trials.
		BOOST_CHECK_CLOSE(kpftimes::mean(dblList[nTest].begin(), dblList[nTest].end()), 
				trueMean, TEST_TOLERANCE*100.0);
		
		//@test Vector of doubles, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_mean to within 1e-10 in 10 out of 10 trials.
		BOOST_CHECK_CLOSE(kpftimes::mean(dblVec[nTest].begin(), dblVec[nTest].end()), 
				trueMean, TEST_TOLERANCE*100.0);
		
		//@test Array of doubles, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_mean to within 1e-10 in 10 out of 10 trials.
		BOOST_CHECK_CLOSE(kpftimes::mean(dblArray[nTest], dblArray[nTest]+TEST_LEN), 
				trueMean, TEST_TOLERANCE*100.0);
	}
}

BOOST_AUTO_TEST_CASE(variance)
{
	//@test List of ints, length 0. Expected behavior: throw invalid_argument.
	BOOST_CHECK_THROW(kpftimes::variance(emptyList.begin(), emptyList.end()), 
			std::invalid_argument);
	
	//@test List of ints, length 1. Expected behavior: throw invalid_argument.
	BOOST_CHECK_THROW(kpftimes::variance(oneList.begin(), oneList.end()), 
			std::invalid_argument);
	
	//@test List of ints, length 2. Expected behavior: return 
	//		(list::back()-list::front())^2/2
	{
		int trueVar = twoList.back()-twoList.front();
		trueVar *= trueVar;
		trueVar /= 2;
		BOOST_CHECK_EQUAL(kpftimes::variance(twoList.begin(), twoList.end()), trueVar);
	}
	
	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		int trueVar = static_cast<int>(gsl_stats_int_variance(intArray[nTest], 1, 100));
		//@test List of ints, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_variance in 10 out of 10 trials.
		BOOST_CHECK_EQUAL(kpftimes::variance(intList[nTest].begin(), 
				intList[nTest].end()), trueVar);
	}

	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		double trueVar = gsl_stats_variance(dblArray[nTest], 1, 100);
		//@test List of doubles, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_variance to within 1e-10 in 10 out of 10 trials.
		BOOST_CHECK_CLOSE(kpftimes::variance(dblList[nTest].begin(), 
				dblList[nTest].end()),     trueVar, TEST_TOLERANCE*100.0);
		
		//@test Vector of doubles, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_variance to within 1e-10 in 10 out of 10 trials.
		BOOST_CHECK_CLOSE(kpftimes::variance(dblVec[nTest].begin(), 
				dblVec[nTest].end()),      trueVar, TEST_TOLERANCE*100.0);
		
		//@test Array of doubles, length 100, randomly generated. Expected behavior: 
		//	agrees with gsl_stats_variance to within 1e-10 in 10 out of 10 trials.
		BOOST_CHECK_CLOSE(kpftimes::variance(dblArray[nTest], 
				dblArray[nTest]+TEST_LEN), trueVar, TEST_TOLERANCE*100.0);
	}
}

BOOST_AUTO_TEST_SUITE_END()
