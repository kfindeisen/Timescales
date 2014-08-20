/** Exception definitions for Timescales library
 * @file timescales/timeexcept.h
 * @author Krzysztof Findeisen
 * @date Created November 18, 2013
 * @date Last modified November 18, 2013
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

#ifndef TIMEEXCEPTH
#define TIMEEXCEPTH

#include <stdexcept>
#include <string>

namespace kpftimes { namespace except {

using std::string;

/** This exception is thrown if a data set does not represent a variable 
 *	time series.
 */
class BadLightCurve : public std::invalid_argument {
public:
	/** Constructs a BadLightCurve object.
	 */
	explicit BadLightCurve(const string& what_arg);
};

/** This exception is thrown if an analysis is run on negative frequencies.
 */
class NegativeFreq : public std::invalid_argument {
public:
	/** Constructs a NegativeFreq object.
	 */
	explicit NegativeFreq(const string& what_arg);
};

/** This exception is thrown if a range is given with minimum greater 
 *	than maximum.
 */
class NegativeRange : public std::invalid_argument {
public:
	/** Constructs a NegativeRange object.
	 */
	explicit NegativeRange(const string& what_arg);
};

}}		// end kpftimes::except

#endif		// end ifndef TIMEEXCEPTH
