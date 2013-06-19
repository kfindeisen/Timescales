/** Test unit for statistics template routines in utils.h
 * @file timescales/tests/unit_stats.cpp
 * @author Krzysztof Findeisen
 * @date Created July 20, 2011
 * @date Last modified June 17, 2013
 */

#include "../warnflags.h"

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

#include <list>
#include <stdexcept>
#include <vector>
#include <boost/smart_ptr.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>
#include "../utils.h"

namespace kpftimes { namespace test {

using boost::shared_array;
using boost::shared_ptr;
using namespace std;

/** Data common to the test cases.
 *
 * Contains arrays, lists, and vectors with mock data.
 */
class StatsData {
public: 
	const static size_t TEST_LEN       = 100;
	const static size_t TEST_COUNT     =  10;
	const static double TEST_TOLERANCE = 1e-8;

	/** Defines the data for each test case.
	 *
	 * @exception std::bad_alloc Thrown if there is not enough memory to 
	 *	store the testing data.
	 *
	 * @exceptsafe Object construction is atomic.
	 */
	StatsData() : emptyList(), oneList(1, 42), twoList(), 
			intArray(), intList(), dblArray(), dblList(), dblVec() {
		twoList.push_back(-10);
		twoList.push_back(27);
	
		// Data come from a random number generator
		shared_ptr<gsl_rng> testRng(gsl_rng_alloc(gsl_rng_mt19937), &gsl_rng_free);
		if (testRng.get() == NULL) {
			throw bad_alloc();
		}
		gsl_rng_set(testRng.get(), 42);
		
		// Define the arrays, then the STL versions
		for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
			intArray.push_back(shared_array<int   >(new int   [TEST_LEN]));
			intList.push_back(  list<int   >());

			dblArray.push_back(shared_array<double>(new double[TEST_LEN]));
			dblList.push_back(  list<double>());
			dblVec .push_back(vector<double>());
			dblVec.back().reserve(TEST_LEN);
			for (size_t i = 0; i < TEST_LEN; i++) {
				{
					// Large sigma to allow interesting truncation to int
					int intVal = static_cast<int>(
						gsl_ran_gaussian(testRng.get(), 1000.0));
					intArray.back()[i] = intVal;
					intList .back().push_back(intVal);
				}

				{
					double dblVal = gsl_ran_ugaussian(testRng.get());
					dblArray.back()[i] = dblVal;
					dblList .back().push_back(dblVal);
					dblVec  .back().push_back(dblVal);
				}
			}
		}
	}
	
	/** Empty list of ints
	 */
	list<int> emptyList;
	/** Length-1 list of ints
	 */
	list<int> oneList;
	/** Length-2 list of ints
	 */
	list<int> twoList;
	
	/** Length-100 collections of integers
	 */
	vector<shared_array<int   > > intArray;
	vector<        list<int   > > intList;
	/** Length-100 collections of floating point numbers
	 */
	vector<shared_array<double> > dblArray;
	vector<        list<double> > dblList;
	vector<      vector<double> > dblVec;
};

/** Test cases for generic mean and variance functions
 * @class BoostTest::test_stats
 */
BOOST_FIXTURE_TEST_SUITE(test_stats, StatsData)

/** Tests whether mean() works as advertised
 *
 * @exceptsafe Does not throw exceptions.
 */
//BOOST_AUTO_TEST_CASE(mean)
//{
//	/** @test List of ints, length 0. Expected behavior: throw invalid_argument.
//	 */
//	BOOST_CHECK_THROW(kpftimes::mean(emptyList.begin(), emptyList.end()), 
//		std::invalid_argument);
//	
//	/** @test List of ints, length 1. Expected behavior: return list[0]
//	 */
//	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(
//		kpftimes::mean(oneList.begin(), oneList.end()), 
//		oneList.front()));
//	
//	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
//		int trueMean = static_cast<int>(
//			gsl_stats_int_mean(intArray[nTest].get(), 1, TEST_LEN));
//		/** @test List of ints, length 100, randomly generated. Expected behavior: 
//		 * 	agrees with gsl_stats_mean in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(
//			kpftimes::mean(intList[nTest].begin(), intList[nTest].end()), 
//			trueMean));
//		
//		/** @test Array of ints, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_mean in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(
//			kpftimes::mean(intArray[nTest].get(), intArray[nTest].get()+TEST_LEN), 
//			trueMean));
//	}
//
//	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
//		double trueMean = gsl_stats_mean(dblArray[nTest].get(), 1, TEST_LEN);
//		/** @test List of doubles, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_mean to within 1e-8 in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_CLOSE(
//			kpftimes::mean(dblList[nTest].begin(), dblList[nTest].end()), 
//			trueMean, TEST_TOLERANCE));
//		
//		/** @test Vector of doubles, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_mean to within 1e-8 in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_CLOSE(
//			kpftimes::mean(dblVec[nTest].begin(), dblVec[nTest].end()), 
//			trueMean, TEST_TOLERANCE));
//		
//		/** @test Array of doubles, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_mean to within 1e-8 in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_CLOSE(
//			kpftimes::mean(dblArray[nTest].get(), dblArray[nTest].get()+TEST_LEN), 
//			trueMean, TEST_TOLERANCE));
//	}
//}
//
///** Tests whether variance() works as advertised
// *
// * @exceptsafe Does not throw exceptions.
// */
//BOOST_AUTO_TEST_CASE(variance)
//{
//	/** @test List of ints, length 0. Expected behavior: throw invalid_argument.
//	 */
//	BOOST_CHECK_THROW(kpftimes::variance(emptyList.begin(), emptyList.end()), 
//		std::invalid_argument);
//	
//	/** @test List of ints, length 1. Expected behavior: throw invalid_argument.
//	 */
//	BOOST_CHECK_THROW(kpftimes::variance(oneList.begin(), oneList.end()), 
//		std::invalid_argument);
//	
//	/** @test List of ints, length 2. Expected behavior: return 
//	 *		(list::back()-list::front())^2/2
//	 */
//	{
//		int trueVar = twoList.back()-twoList.front();
//		trueVar *= trueVar;
//		trueVar /= 2;
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(
//			kpftimes::variance(twoList.begin(), twoList.end()), 
//			trueVar));
//	}
//	
//	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
//		int trueVar = static_cast<int>(
//			gsl_stats_int_variance(intArray[nTest].get(), 1, TEST_LEN));
//		/** @test List of ints, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_variance in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(
//			kpftimes::variance(intList[nTest].begin(), intList[nTest].end()), 
//			trueVar));
//	}
//
//	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
//		double trueVar = gsl_stats_variance(dblArray[nTest].get(), 1, TEST_LEN);
//		/** @test List of doubles, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_variance to within 1e-8 in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_CLOSE(
//			kpftimes::variance(dblList[nTest].begin(), dblList[nTest].end()),
//			trueVar, TEST_TOLERANCE));
//		
//		/** @test Vector of doubles, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_variance to within 1e-8 in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_CLOSE(
//			kpftimes::variance(dblVec[nTest].begin(), dblVec[nTest].end()), 
//			trueVar, TEST_TOLERANCE));
//		
//		/** @test Array of doubles, length 100, randomly generated. Expected behavior: 
//		 *	agrees with gsl_stats_variance to within 1e-8 in 10 out of 10 trials.
//		 */
//		BOOST_CHECK_NO_THROW(BOOST_CHECK_CLOSE(
//			kpftimes::variance(dblArray[nTest].get(), 
//			dblArray[nTest].get()+TEST_LEN), 
//			trueVar, TEST_TOLERANCE));
//	}
//}

BOOST_AUTO_TEST_SUITE_END()

}}	// end kpftimes::test
