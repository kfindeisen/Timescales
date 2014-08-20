/** Implements support code for the library. None of these routines are 
 * intended as part of the public API.
 * @file timescales/utils.cpp
 * @author Krzysztof Findeisen
 * @date Created April 13, 2011
 * @date Last modified November 19, 2013
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

#include <algorithm>
#include <stdexcept>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr.hpp>
#include "utils.h"

namespace kpftimes {

using boost::lexical_cast;
using std::string;

/** Constructs but does not initialize a table of fixed size.
 * 
 * @param[in] initX The first dimension of the table.
 * @param[in] initY The first dimension of the table.
 *
 * @pre @p initX > 0 and @p initY > 0
 * @post The object represents an @p initX by @p initY table of values. The table 
 *	elements are not constrained.
 *
 * @exception std::invalid_argument Thrown if the table dimensions aren't valid.
 * @exception std::bad_alloc Thrown if there is not enough memory to construct 
 *	the table.
 * 
 * @exceptsafe Object construction is atomic.
 *
 * @test FastTable(-1, 1). Expected behavior: throw invalid_argument
 * @test FastTable(0, 1). Expected behavior: throw invalid_argument
 * @test FastTable(1, -1). Expected behavior: throw invalid_argument
 * @test FastTable(1, 0). Expected behavior: throw invalid_argument
 * @test FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop.
 * @test FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop.
 * @test FastTable(5, 4). Expected behavior: success. Can iterate over elements with external loop.
 */
FastTable::FastTable(size_t initX, size_t initY)
		: table(NULL), dimX(initX), dimY(initY) {
	if (initX <= 0 || initY <= 0) {
		try {
			throw std::invalid_argument("Cannot create FastTable with zero dimensions (gave " + lexical_cast<string>(initX) + "×" + lexical_cast<string>(initY) + ").");
		} catch(const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Cannot create FastTable with zero dimensions.");
		}
	}
	
	// Allocating table is the only operation that could potentially 
	//	violate strong exception guarantee
	// Do it last to ensure that it's allocated only if nothing else 
	//	throws an exception
	table = new double[initX*initY];
}

/** Constructs an independent copy of the given object, containing identical data.
 * 
 * @param[in] other The object to copy.
 *
 * @post The object is identical to @p other, but may be modified or 
 *	deleted independently.
 *
 * @exception std::bad_alloc Thrown if there is not enough memory to construct 
 *	the table.
 * 
 * @exceptsafe Object construction is atomic.
 *
 * @test FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(5, 4). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 */
FastTable::FastTable(const FastTable &other) : 
		// This allocation is safe because it is the only operation 
		// that can throw an exception
		table(new double[other.dimX*other.dimY]), 
		dimX(other.dimX), dimY(other.dimY) {
	size_t n = dimX*dimY;
	for(size_t i = 0; i < n; i++) {
		this->table[i] = other.table[i];
	}
}

/** Alters the existing object to be identical to an existing object.
 *
 * @param[in] other The object to copy.
 *
 * @return An assignable lvalue to this object.
 *
 * @pre This object has the same dimensions as @p other.
 *
 * @post Any data previously in this object is deleted.
 * @post The object is identical to @p other, but may be modified or 
 *	deleted independently.
 * 
 * @exception std::invalid_argument Thrown if the two tables have mismatched dimensions.
 * @exception std::bad_alloc Thrown if there is not enough memory to construct 
 *	the table.
 *
 * @exceptsafe Neither object is changed in the event of an exception.
 *
 * @test FastTable(1, 1) = FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(4, 5) = FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(4, 5) = FastTable(5, 5). Expected behavior: throw invalid_argument.
 * @test FastTable(5, 4) = FastTable(5, 5). Expected behavior: throw invalid_argument.
 */
FastTable& FastTable::operator=(const FastTable &other) {
	if(this->dimX != other.dimX || this->dimY != other.dimY) {
		try {
			// For some reason lexical_cast crashes if you use it 
			//	more than twice here. No time to isolate the bug now.
			throw std::invalid_argument("Cannot assign to a table of mismatched dimensions (lvalue is " 
			+ lexical_cast<string>(this->dimX) + "×" + lexical_cast<string>(this->dimY) 
			+ ", rvalue is " 
			+ lexical_cast<string>(other.dimX) + "×" + lexical_cast<string>(other.dimY) 
			+ ").");
		} catch(const boost::bad_lexical_cast& e) {
			throw std::invalid_argument("Cannot assign to a table of mismatched dimensions.");
		}
	}

	// copy-and-swap	
	size_t n = dimX*dimY;
	double* newTable = new double[n];
	
	// IMPORTANT: no exceptions beyond this point

	for(size_t i = 0; i < n; i++) {
		newTable[i] = other.table[i];
	}
	
	using std::swap;
	swap(this->table, newTable);
	
	// newTable now contains the previous table pointed to by this->table
	delete [] newTable;
	
	return *this;
}

FastTable::~FastTable() {
	delete [] this->table;
}

/** Access function to allow reads and writes of a table element. Elements 
 *	with the same x coordinate and adjacent y coordinates are adjacent in 
 *	memory.
 * 
 * @param[in] x The first-dimension coordinate to be read.
 * @param[in] y The second-dimension coordinate to be read.
 * 
 * @return An lvalued reference to the element Table[x, y]
 * 
 * @pre 0 <= x < dimX
 * @pre 0 <= y < dimY
 *
 * @warning In the interest of speed at all costs, preconditions are @em not 
 *	checked.
 *
 * @exceptsafe Does not throw exceptions.
 *
 * @test FastTable(4, 5).at(0, 2). Expected behavior: update reflected by separate .at().
 * @test FastTable(4, 5).at(3, 2). Expected behavior: update reflected by separate .at().
 * @test FastTable(4, 5).at(2, 0). Expected behavior: update reflected by separate .at().
 * @test FastTable(4, 5).at(2, 4). Expected behavior: update reflected by separate .at().
 * @test FastTable(1000, 1000).at(x,y) = x+y±x*y iteration in both directions. 
 *	Expected behavior: y-inner loop faster than x-inner loop for 10 out of 10 objects.
 */
double& FastTable::at(size_t x, size_t y) {
	// need adjacent y-indices to be adjacent in memory
	return table[dimY*x + y];
}

}
