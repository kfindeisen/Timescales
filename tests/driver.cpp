/** Test code for Timescales
 * @file driver.cpp
 * @author Krzysztof Findeisen
 * @date Created May 23, 2011
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

/** Name of the unit test module
 */
#define BOOST_TEST_MODULE TimescaleTesting
#include <boost/test/unit_test.hpp>

// Re-enable all compiler warnings
#ifdef GNUC_FINEWARN
#pragma GCC diagnostic pop
#endif

// All other code in the individual test libraries
