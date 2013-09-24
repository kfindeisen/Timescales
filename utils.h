/** Support code for the library. None of these routines are 
 * intended as part of the public API.
 * @file timescales/utils.h
 * @author Krzysztof Findeisen
 * @date Created April 13, 2011
 * @date Last modified November 19, 2013
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
