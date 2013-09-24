/** Exceptions for Timescales library: invalid function options
 * @file timescales/badoption.cpp
 * @author Krzysztof Findeisen
 * @date Created November 18, 2013
 * @date Last modified November 18, 2013
 */

#include <stdexcept>
#include <string>
#include "timeexcept.h"

namespace kpftimes { namespace except {

using std::string;

/** Constructs a NegativeFreq object.
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
NegativeFreq::NegativeFreq(const string& what_arg) : invalid_argument(what_arg) {
}

/** Constructs a NegativeRange object.
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
NegativeRange::NegativeRange(const string& what_arg) : invalid_argument(what_arg) {
}

}}		// end kpftimes::except
