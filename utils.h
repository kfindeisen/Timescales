/** Support code for the library. None of these routines are intended as part of the public API.
 * @file timescales/utils.h
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

#include <stdexcept>
#include <vector>
 
/** A convenient shorthand for vectors of doubles.
 */
typedef std::vector<double> DoubleVec;

namespace kpftimes {

//----------------------------------------------------------
/** @defgroup util Utility functions
 *  @{
 */
 
/** Two-dimensional dynamically allocated table. The table is designed to 
 *	allow fast (i.e. localized) scans along the Y axis. It has 
 *	barely any other functionality.
 */
class FastTable {
public:
	/** Constructs but does not initialize a table of fixed size.
	 */
	FastTable(size_t dimX, size_t dimY);
	FastTable(const FastTable& otherTable);
	FastTable& operator=(const FastTable& otherTable);
	~FastTable();

	/** Access function to allow reads and writes of a table element.
	 */
	double& at(size_t x, size_t y);
	
	/** Returns the X (outer) dimension of the FastTable
	 */
	size_t getX() const;

	/** Returns the Y (inner) dimension of the FastTable
	 */
	size_t getY() const;
private:
	double* table;
	const size_t dimX, dimY;
};


/** @} */

}
