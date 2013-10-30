/** Test code for Timescales
 * @file driver.cpp
 * @author Krzysztof Findeisen
 * @date Created May 23, 2011
 * @date Last modified June 17, 2013
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

/** Name of the unit test module
 */
#define BOOST_TEST_MODULE TimescaleTesting
#include <boost/test/unit_test.hpp>

// Re-enable all compiler warnings
#ifdef GNUC_FINEWARN
#pragma GCC diagnostic pop
#endif

namespace kpftimes { namespace test {

/** This function is a wrapper for a trusted approximate comparison method.
 * 
 * @param[in] val1, val2 The values to compare
 * @param[in] frac The fractional difference allowed between them
 * 
 * @return true iff |val1 - val2|/|val1| and |val1 - val2|/|val2| < @p frac
 *
 * @exceptsafe Does not throw exceptions.
 */
bool isClose(double val1, double val2, double frac) {
	using namespace ::boost::test_tools;
	
	return static_cast<bool>(check_is_close(val1, val2, 
			fraction_tolerance_t<double>(frac)));
}

}}		// end kpftimes::test
