/** Test unit for FastTable class
 * @file timescales/tests/unit_FastTable.cpp
 * @author Krzysztof Findeisen
 * @date Created July 20, 2011
 * @date Last modified August 21, 2013
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

#include <ctime>
#include "../utils.h"

namespace kpftimes { namespace test {

/** Data common to the test cases.
 *
 * Contains @ref FastTable "FastTables" of varying dimensions
 */
class FtData {
public: 
	const static size_t TEST_COUNT   =   10;
	const static size_t TEST_BIGSIZE = 1000;

	/** Defines the data for each test case.
	 *
	 * @exception std::bad_alloc Thrown if there is not enough memory to 
	 *	store the testing data.
	 *
	 * @exceptsafe Object construction is atomic.
	 */
	FtData() : FastTable11(1,1), FastTable45(4,5), FastTable54(5,4) {
		BOOST_CHECK_NO_THROW(FastTable11.at(1, 1) = 42.0);
		
		for(size_t x = 0; x < 4; x++) {
			for(size_t y = 0; y < 5; y++) {
				BOOST_CHECK_NO_THROW(FastTable45.at(x, y) = x + (x+1)*(y+1) - y);
			}
		}
		for(size_t x = 0; x < 5; x++) {
			for(size_t y = 0; y < 4; y++) {
				BOOST_CHECK_NO_THROW(FastTable54.at(x, y) = y + (x+1)*(y+1) - x);
			}
		}
	}
	
	virtual ~FtData() {
	}
	
	/** A one-element table
	 */
	FastTable FastTable11;

	/** A fast twenty-element table
	 */
	FastTable FastTable45;

	/** A slow twenty-element table
	 */
	FastTable FastTable54;
};

/** Tests whether two @ref FastTable "FastTables" are identical but independent copies
 *
 * @param[in]  oldTable The original copy of the table
 * @param[in] copyTable The new copy of the table
 * @param[in] sizeX, sizeY The dimensions of both @p oldTable and @p copyTable
 *
 * @post @p oldTable and @p copyTable may be modified as part of the testing
 *
 * @exceptsafe Does not throw exceptions.
 */
void testCopy(FastTable& oldTable, FastTable& copyTable, size_t sizeX, size_t sizeY) {
	// Overall test strategy:
	//	1. Make a copy
	//	2. Verify that it's equal to the original
	//	3. Edit the copy
	//	4. Verify that it's no longer equal to the original
	
	for(size_t x = 0; x < sizeX; x++) {
		for(size_t y = 0; y < sizeY; y++) {
			BOOST_CHECK_NO_THROW(
				BOOST_CHECK_EQUAL(copyTable.at(x, y), oldTable.at(x, y)) );
		}
	}
	for(size_t x = 0; x < sizeX; x++) {
		for(size_t y = 0; y < sizeY; y++) {
			BOOST_CHECK_NO_THROW(copyTable.at(x, y) += 42.0);
		}
	}
	for(size_t x = 0; x < sizeX; x++) {
		for(size_t y = 0; y < sizeY; y++) {
			BOOST_CHECK_NO_THROW(
				BOOST_CHECK_NE(copyTable.at(x, y), oldTable.at(x, y)) );
		}
	}
}

// Boost.Test uses C-style casts and non-virtual destructors
#ifdef GNUC_FINEWARN
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif

/** Test cases for FastTable construction and other fundamental operations 
 *	that must not rely on FtData
 * 
 * @class BoostTest::test_FastTableBuild
 */
BOOST_AUTO_TEST_SUITE(test_FastTableBuild)

/** Tests whether FastTable constructor accepts appropriate arguments
 *
 * @exceptsafe Does not throw exceptions.
 */
BOOST_AUTO_TEST_CASE(constructors) {
	//@test FastTable(0, 1). Expected behavior: throw invalid_argument
	BOOST_CHECK_THROW(FastTable(0, 1), std::invalid_argument);

	//@test FastTable(1, 0). Expected behavior: throw invalid_argument
	BOOST_CHECK_THROW(FastTable(1, 0), std::invalid_argument);
	
	// Valid constructors mostly tested in FtData
	// Just make sure they don't throw
	BOOST_CHECK_NO_THROW(FastTable(1,1));
	BOOST_CHECK_NO_THROW(FastTable(4,5));
	BOOST_CHECK_NO_THROW(FastTable(5,4));
}

BOOST_AUTO_TEST_SUITE_END()

// Re-enable all compiler warnings
#ifdef GNUC_FINEWARN
#pragma GCC diagnostic pop
#endif

/** Test cases for general FastTable operations
 * 
 * @class BoostTest::test_FastTable
 */
BOOST_FIXTURE_TEST_SUITE(test_FastTable, FtData)

/** Tests whether FastTable element access behaves consistently
 *
 * @exceptsafe Does not throw exceptions.
 */
BOOST_AUTO_TEST_CASE(access) {
	/** @test FastTable(4, 5).at(3, 2). Expected behavior: update reflected by separate .at().
	 */
	BOOST_CHECK_NO_THROW(FastTable45.at(3, 2) = 42.2);
	// Exact equality check, because no math operations should have occurred
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(3, 2), 42.2));
	
	/** @test FastTable(4, 5).at(0, 2). Expected behavior: update reflected by separate 
	 .at().
	 */
	BOOST_CHECK_NO_THROW(FastTable45.at(0, 2) = -10.56);
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(0, 2), -10.56));
	
	/** @test FastTable(4, 5).at(2, 0). Expected behavior: update reflected by separate .at().
	 */
	BOOST_CHECK_NO_THROW(FastTable45.at(2, 0) = 1234.56789);
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(2, 0), 1234.56789));

	/** @test FastTable(4, 5).at(2, 4). Expected behavior: update reflected by separate .at().
	 */
	BOOST_CHECK_NO_THROW(FastTable45.at(2, 4) = 1e-5);
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(2, 4), 1e-5));
	
	// Now check that we haven't overwritten anything
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(3, 2), 42.2));
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(2, 0), 1234.56789));
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(2, 4), 1e-5));
	BOOST_CHECK_NO_THROW(BOOST_CHECK_EQUAL(FastTable45.at(0, 2), -10.56));
}

/** Tests whether the FastTable copy-constructor creates correct copies
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to copy the objects.
 *
 * @exceptsafe The function arguments and the test case are unchanged in 
 *	the event of an exception.
 */
BOOST_AUTO_TEST_CASE(copyconstruction) {
	// This is a nontrivial test, so it's been placed in testCopy()
		
	//@test FastTable(1, 1). Expected behavior: success. Can iterate over 
	//	elements with external loop and edit independently.
	FastTable newTable1(FastTable11);
	testCopy(FastTable11, newTable1, 1, 1);
	
	//@test FastTable(4, 5). Expected behavior: success. Can iterate over 
	//	elements with external loop and edit independently.
	FastTable newTable2(FastTable45);
	testCopy(FastTable45, newTable2, 4, 5);

	//@test FastTable(5, 4). Expected behavior: success. Can iterate over 
	//	elements with external loop and edit independently.
	FastTable newTable3(FastTable54);
	testCopy(FastTable54, newTable3, 5, 4);
}

/** Tests whether the FastTable assignment operator creates correct copies
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to copy the objects.
 *
 * @exceptsafe The function arguments and the test case are unchanged in 
 *	the event of an exception.
 */
BOOST_AUTO_TEST_CASE(assignment) {
	// This is a nontrivial test, so it's been placed in testCopy()
	
	//@test FastTable(1, 1) = FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop and edit independently.
	FastTable newTable1(1, 1);
	BOOST_CHECK_NO_THROW(newTable1 = FastTable11);
	testCopy(FastTable11, newTable1, 1, 1);
	
	//@test FastTable(4, 5) = FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop and edit independently.
	FastTable newTable2(4, 5);
	BOOST_CHECK_NO_THROW(newTable2 = FastTable45);
	testCopy(FastTable45, newTable2, 4, 5);
	
	//@test FastTable(4, 5) = FastTable(5, 4). Expected behavior: throw invalid_argument.
	//@test FastTable(5, 4) = FastTable(4, 5). Expected behavior: throw invalid_argument.
	FastTable newTable3(4, 5), newTable4(5, 4);
	BOOST_CHECK_THROW(newTable3 = FastTable54, std::invalid_argument);
	BOOST_CHECK_THROW(newTable4 = FastTable45, std::invalid_argument);
}

/** Measures the time needed to traverse a table in y, then x
 *
 * @param[in] table The FastTable to update
 * @param[in] dimX, dimY The dimensions of @p table
 *
 * @return The number of clock ticks needed to iterate once over @p table.
 *
 * @exceptsafe Does not throw exceptions.
 */
long clockFast(FastTable& table, size_t dimX, size_t dimY) {
	clock_t start, end;
	start = clock();
	for(size_t x = 0; x < dimX; x++) {
		for(size_t y = 0; y < dimY; y++) {
			table.at(x, y) = x + y + x*y;
		}
	}
	end = clock();
	
	return end - start;
}

/** Measures the time needed to traverse a table in x, then y
 *
 * @param[in] table The FastTable to update
 * @param[in] dimX, dimY The dimensions of @p table
 *
 * @return The number of clock ticks needed to iterate once over @p table.
 *
 * @exceptsafe Does not throw exceptions.
 */
long clockSlow(FastTable& table, size_t dimX, size_t dimY) {
	clock_t start, end;
	start = clock();
	for(size_t y = 0; y < dimY; y++) {
		for(size_t x = 0; x < dimX; x++) {
			table.at(x, y) = x + y + x*y;
		}
	}
	end = clock();
	
	return end - start;
}

/** Tests whether a FastTable allocated on the stack has faster element 
 *	access in the y direction than the x direction
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to run the test.
 *
 * @exceptsafe The function arguments and the test case are unchanged in 
 *	the event of an exception.
 */
BOOST_AUTO_TEST_CASE(stack_speed) {
	/** @test FastTable(1000, 1000).at(x,y) = x+y±x*y iteration in both directions. 
	 *	Expected behavior: y-inner loop faster than x-inner loop on average over 10 
	 *	objects.
	 */
	long fastTotal = 0;
	long slowTotal = 0;
	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		FastTable speedTable(TEST_BIGSIZE, TEST_BIGSIZE);
		
		fastTotal += clockFast(speedTable, TEST_BIGSIZE, TEST_BIGSIZE);

		slowTotal += clockSlow(speedTable, TEST_BIGSIZE, TEST_BIGSIZE);
	}
	// Are we faster on average?
	BOOST_CHECK_LT(fastTotal, slowTotal);
}

/** Tests whether a FastTable allocated on the heap has faster element 
 *	access in the y direction than the x direction
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to run the test.
 *
 * @exceptsafe The function arguments and the test case are unchanged in 
 *	the event of an exception.
 */
BOOST_AUTO_TEST_CASE(heap_speed) {
	/** @test FastTable(1000, 1000).at(x,y) = x+y±x*y iteration in both directions. 
	 *	Expected behavior: y-inner loop faster than x-inner loop on average over 10 
	 *	objects.
	 */
	long fastTotal = 0;
	long slowTotal = 0;
	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		FastTable *speedTable = new FastTable(TEST_BIGSIZE, TEST_BIGSIZE);
		fastTotal += clockFast(*speedTable, TEST_BIGSIZE, TEST_BIGSIZE);
		delete speedTable;		

		speedTable = new FastTable(TEST_BIGSIZE, TEST_BIGSIZE);
		slowTotal += clockSlow(*speedTable, TEST_BIGSIZE, TEST_BIGSIZE);

		delete speedTable;
	}
	// Are we faster on average?
	BOOST_CHECK_LT(fastTotal, slowTotal);
}

BOOST_AUTO_TEST_SUITE_END()

}}		// end kpftimes::test
