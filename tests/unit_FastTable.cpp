/** Test unit for FastTable class
 * @file unit_FastTable.cpp
 * @author Krzysztof Findeisen
 * @date Created July 20, 2011
 * @date Last modified July 20, 2011
 */

#include <ctime>

#include <boost/test/unit_test.hpp>
#include "../utils.h"

using kpftimes::FastTable;

// Data common to the test cases
class FtData {
public: 
	const static size_t TEST_COUNT   =   10;
	const static size_t TEST_BIGSIZE = 1000;

	//@test FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop.
	//@test FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop.
	//@test FastTable(5, 4). Expected behavior: success. Can iterate over elements with external loop.
	FtData() : FastTable11(1,1), FastTable45(4,5), FastTable54(5,4) {
printf("DEBUG: constructing FtData\n");
		FastTable11.at(1, 1) = 42.0;
		
		for(size_t x = 0; x < 4; x++) {
			for(size_t y = 0; y < 5; y++) {
				FastTable45.at(x, y) = x + (x+1)*(y+1) - y;
			}
		}
		for(size_t x = 0; x < 5; x++) {
			for(size_t y = 0; y < 4; y++) {
				FastTable54.at(x, y) = y + (x+1)*(y+1) - x;
			}
		}
printf("DEBUG: FtData complete\n");
	}
	
	~FtData() {
printf("DEBUG: destroying FtData\n");
	}
	
	void testCopy(FastTable& oldTable, FastTable& copyTable, size_t sizeX, size_t sizeY) {
		// Overall test strategy:
		//	1. Make a copy
		//	2. Verify that it's equal to the original
		//	3. Edit the copy
		//	4. Verify that it's no longer equal to the original
		
		for(size_t x = 0; x < sizeX; x++) {
			for(size_t y = 0; y < sizeY; y++) {
				BOOST_CHECK_EQUAL(copyTable.at(x, y), oldTable.at(x, y));
			}
		}
		for(size_t x = 0; x < sizeX; x++) {
			for(size_t y = 0; y < sizeY; y++) {
				copyTable.at(x, y) = copyTable.at(x, y) + 42.0;
			}
		}
		for(size_t x = 0; x < sizeX; x++) {
			for(size_t y = 0; y < sizeY; y++) {
				BOOST_CHECK_NE(copyTable.at(x, y), oldTable.at(x, y));
			}
		}
	}
	
	FastTable FastTable11, FastTable45, FastTable54;
};

BOOST_FIXTURE_TEST_SUITE(test_FastTable, FtData)

BOOST_AUTO_TEST_CASE(constructors) {
printf("DEBUG: beginning constructors check\n");
	//@test FastTable(-1, 1). Expected behavior: throw invalid_argument
	//BOOST_CHECK_THROW(FastTable(-1, 1), std::invalid_argument);
	
	//@test FastTable(0, 1). Expected behavior: throw invalid_argument
	BOOST_CHECK_THROW(FastTable(0, 1), std::invalid_argument);

	//@test FastTable(1, -1). Expected behavior: throw invalid_argument
	//BOOST_CHECK_THROW(FastTable(1, -1), std::invalid_argument);

	//@test FastTable(1, 0). Expected behavior: throw invalid_argument
	BOOST_CHECK_THROW(FastTable(1, 0), std::invalid_argument);
	
	// Valid constructors already used in FtData
printf("DEBUG: constructors check complete\n");
}

BOOST_AUTO_TEST_CASE(access) {
printf("DEBUG: beginning access check\n");
	//@test FastTable(4, 5).at(3, 2). Expected behavior: update reflected by separate .at().
	FastTable45.at(3, 2) = 42.2;
	// Exact equality check, because no math operations should have occurred
	BOOST_CHECK_EQUAL(FastTable45.at(3, 2), 42.2);
	
	//@test FastTable(4, 5).at(0, 2). Expected behavior: update reflected by separate .at().
	FastTable45.at(0, 2) = -10.56;
	BOOST_CHECK_EQUAL(FastTable45.at(0, 2), -10.56);
	
	//@test FastTable(4, 5).at(2, 0). Expected behavior: update reflected by separate .at().
	FastTable45.at(2, 0) = 1234.56789;
	BOOST_CHECK_EQUAL(FastTable45.at(2, 0), 1234.56789);

	//@test FastTable(4, 5).at(2, 4). Expected behavior: update reflected by separate .at().
	FastTable45.at(2, 4) = 1e-5;
	BOOST_CHECK_EQUAL(FastTable45.at(2, 4), 1e-5);
	
	// Now check that we haven't overwritten anything
	BOOST_CHECK_EQUAL(FastTable45.at(3, 2), 42.2);
	BOOST_CHECK_EQUAL(FastTable45.at(2, 0), 1234.56789);
	BOOST_CHECK_EQUAL(FastTable45.at(2, 4), 1e-5);
	BOOST_CHECK_EQUAL(FastTable45.at(0, 2), -10.56);
printf("DEBUG: access check complete\n");
}

/*BOOST_AUTO_TEST_CASE(stack_speed) {
	//@test FastTable(1000, 1000).at(x,y) = x+y±x*y iteration in both directions. 
	//	Expected behavior: y-inner loop faster than x-inner loop for 10 out of 10 
	//	objects.
	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		FastTable speedTable(TEST_BIGSIZE, TEST_BIGSIZE);
		
		time_t startFast, endFast;
		startFast = time(NULL);
		for(size_t x = 0; x < TEST_BIGSIZE; x++) {
			for(size_t y = 0; y < TEST_BIGSIZE; y++) {
				speedTable.at(x, y) = x + y + x*y;
			}
		}
		endFast = time(NULL);

		time_t startSlow, endSlow;
		startSlow = time(NULL);
		for(size_t y = 0; y < TEST_BIGSIZE; y++) {
			for(size_t x = 0; x < TEST_BIGSIZE; x++) {
				speedTable.at(x, y) = x + y + x*y;
			}
		}
		endSlow = time(NULL);

		// Are we faster?
		BOOST_CHECK_LT(difftime(endFast, startFast), difftime(endSlow, startSlow));
	}
}

BOOST_AUTO_TEST_CASE(heap_speed) {
	//@test FastTable(1000, 1000).at(x,y) = x+y±x*y iteration in both directions. 
	//	Expected behavior: y-inner loop faster than x-inner loop for 10 out of 10 
	//	objects.
	for (size_t nTest = 0; nTest < TEST_COUNT; nTest++) {
		FastTable *speedTable = new FastTable(TEST_BIGSIZE, TEST_BIGSIZE);
		
		time_t startFast, endFast;
		startFast = time(NULL);
		for(size_t x = 0; x < TEST_BIGSIZE; x++) {
			for(size_t y = 0; y < TEST_BIGSIZE; y++) {
				speedTable->at(x, y) = x + y + x*y;
			}
		}
		endFast = time(NULL);

		delete speedTable;		

		speedTable = new FastTable(TEST_BIGSIZE, TEST_BIGSIZE);
		
		time_t startSlow, endSlow;
		startSlow = time(NULL);
		for(size_t y = 0; y < TEST_BIGSIZE; y++) {
			for(size_t x = 0; x < TEST_BIGSIZE; x++) {
				speedTable->at(x, y) = x + y + x*y;
			}
		}
		endSlow = time(NULL);

		delete speedTable;
		
		// Are we faster?
		BOOST_CHECK_LT(difftime(endFast, startFast), difftime(endSlow, startSlow));
	}
}

BOOST_AUTO_TEST_CASE(copyconstruction) {
	// This is a nontrivial test, so it's been placed in FtData::testCopy()
		
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

BOOST_AUTO_TEST_CASE(assignment) {
	// This is a nontrivial test, so it's been placed in FtData::testCopy()
		
	//@test FastTable(1, 1) = FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop and edit independently.
	FastTable newTable1(1, 1);
	newTable1 = FastTable11;
	testCopy(FastTable11, newTable1, 1, 1);
	
	//@test FastTable(4, 5) = FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop and edit independently.
	FastTable newTable2(4, 5);
	newTable2 = FastTable45;
	testCopy(FastTable45, newTable2, 4, 5);

	//@test FastTable(4, 5) = FastTable(5, 4). Expected behavior: throw domain_error.
	//@test FastTable(5, 4) = FastTable(4, 5). Expected behavior: throw domain_error.
	FastTable newTable3(4, 5), newTable4(5, 4);
	BOOST_CHECK_THROW(newTable3 = FastTable54, std::domain_error);
	BOOST_CHECK_THROW(newTable4 = FastTable45, std::domain_error);
}*/
	
BOOST_AUTO_TEST_SUITE_END()
