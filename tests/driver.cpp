/** Test code for Timescales
 * @file driver.cpp
 * @author Krzysztof Findeisen
 * @date Created May 23, 2011
 * @date Last modified June 17, 2013
 */

/* Copyright 2014, California Institute of Technology.
 *
 * This file is part of the Timescales library.
 * 
 * The Timescales library is free software: you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version, subject to the following 
 * exception added under Section 7 of the License:
 *	* Neither the name of the copyright holder nor the names of its contributors 
 *	  may be used to endorse or promote products derived from this software 
 *	  without specific prior written permission.
 * 
 * The Timescales library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with the Timescales library. If not, see <http://www.gnu.org/licenses/>.
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
