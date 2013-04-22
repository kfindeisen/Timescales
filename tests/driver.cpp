/** Test code for Timescales
 * @file driver.cpp
 * @author Krzysztof Findeisen
 * @date Created May 23, 2011
 * @date Last modified May 23, 2011
 */

// I have absolutely no idea why I need this when I want static linking
// Nor do I know if it's a Cygwin-only weirdness, but I'll assume it until proven otherwise
#ifdef __CYGWIN__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE TimescaleTesting
#include <boost/test/unit_test.hpp>

// All other code in the individual test libraries
