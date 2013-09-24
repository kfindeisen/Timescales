/** Exception definitions for Timescales library
 * @file timescales/timeexcept.h
 * @author Krzysztof Findeisen
 * @date Created November 18, 2013
 * @date Last modified November 18, 2013
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
