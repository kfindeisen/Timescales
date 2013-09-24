/** Exceptions for Timescales library: unworkable data
 * @file timescales/baddata.cpp
 * @author Krzysztof Findeisen
 * @date Created November 18, 2013
 * @date Last modified November 18, 2013
 */

#include <stdexcept>
#include <string>
#include "timeexcept.h"

namespace kpftimes { namespace except {

using std::string;

/** Constructs a BadLightCurve object.
 *
 * @param[in] what_arg A string with the same content as the value 
 *	returned by what().
 *
 * @post this->what() = @p what_arg.c_str()
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to 
 *	construct the exception.
 * 
 * @exceptsafe Object construction is atomic.
 */
BadLightCurve::BadLightCurve(const string& what_arg) : invalid_argument(what_arg) {
}

}}		// end kpftimes::except
