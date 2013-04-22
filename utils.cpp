/** Implements support code for the library. None of these routines are 
 * intended as part of the public API.
 * @file utils.cpp
 * @author Krzysztof Findeisen
 * @date Created April 13, 2011
 * @date Last modified April 13, 2011
 */ 

#include "utils.h"

using namespace kpftimes;

/** Constructs but does not initialize a table of fixed size.
 * @param[in] initX The first dimension of the table.
 * @param[in] initY The first dimension of the table.
 *
 * @pre initX > 0 and initY > 0
 * @post The object represents an initX by initY table of values. The table 
 *	elements are not constrained.
 *
 * @exception std::invalid_argument Thrown if the table dimensions aren't valid.
 *
 * @note Once the FastTable object has been constructed, its dimensions cannot 
 *	be changed.
 *
	 * @test FastTable(-1, 1). Expected behavior: throw invalid_argument
 * @test FastTable(0, 1). Expected behavior: throw invalid_argument
	 * @test FastTable(1, -1). Expected behavior: throw invalid_argument
 * @test FastTable(1, 0). Expected behavior: throw invalid_argument
 * @test FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop.
 * @test FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop.
 * @test FastTable(5, 4). Expected behavior: success. Can iterate over elements with external loop.
 */
kpftimes::FastTable::FastTable(size_t initX, size_t initY)
		: table(NULL), dimX(initX), dimY(initY) {
	if (initX <= 0 || initY <= 0) {
		throw std::invalid_argument("Cannot create FastTable with negative dimensions");
	}
	table = new double[initX*initY];
}

/** Constructs an independent copy of the given object, containing identical data.
 * @param[in] other The object to copy.
 *
 * @post The object is indistinguishable from other, but may be modified or 
 *	deleted independently.
 *
 * @test FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(5, 4). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 */
kpftimes::FastTable::FastTable(const FastTable &other) : 
		table(new double[other.dimX*other.dimY]), 
		dimX(other.dimX), dimY(other.dimY) {
	size_t n = dimX*dimY;
	for(size_t i = 0; i < n; i++) {
		this->table[i] = other.table[i];
	}
}

/** Alters the existing object to be identical to otherTable.
 *
 * @param[in] other The object to copy.
 *
 * @pre This object has the same dimensions as other.
 * @post The object is indistinguishable from other, but may be modified or 
 *	deleted independently.
 * 
 * @exception domain_error Thrown if the two tables have mismatched dimensions.
 *
 * @test FastTable(1, 1) = FastTable(1, 1). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(4, 5) = FastTable(4, 5). Expected behavior: success. Can iterate over elements with external loop and edit independently.
 * @test FastTable(4, 5) = FastTable(5, 5). Expected behavior: throw domain_error.
 * @test FastTable(5, 4) = FastTAble(5, 5). Expected behavior: throw domain_error.
 */
kpftimes::FastTable& FastTable::operator=(const FastTable &other) {
	if(this->dimX != other.dimX || this->dimY != other.dimY) {
		throw std::domain_error("Mismatched tables; cannot reassign");
	}

	delete [] table;
	size_t n = dimX*dimY;
	this->table = new double[n];

	for(size_t i = 0; i < n; i++) {
		this->table[i] = other.table[i];
	}
	
	return *this;
}

kpftimes::FastTable::~FastTable() {
	delete [] table;
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
 * @warning In the interests of speed-at-all-costs, preconditions are @em not 
 *	checked.
 *
 * @test FastTable(4, 5).at(0, 2). Expected behavior: update reflected by separate .at().
 * @test FastTable(4, 5).at(3, 2). Expected behavior: update reflected by separate .at().
 * @test FastTable(4, 5).at(2, 0). Expected behavior: update reflected by separate .at().
 * @test FastTable(4, 5).at(2, 4). Expected behavior: update reflected by separate .at().
 * @test FastTable(1000, 1000).at(x,y) = x+y±x*y iteration in both directions. 
 *	Expected behavior: y-inner loop faster than x-inner loop for 10 out of 10 objects.
 */
double& kpftimes::FastTable::at(size_t x, size_t y) {
	// need adjacent y-indices to be adjacent in memory
	return table[dimY*x + y];
}

/** Tests whether the argument is sorted
 * 
 * @param[in] list	Times at which data were taken
 * @return TRUE if and only if the argument is sorted in ascending order.
 *
 * @perform O(list.size()) time
 * 
 * @test Empty list. Expected behavior: true.
 * @test List of size 1. Expected behavior: true.
 * @test List of size 2, sorted asc. Expected behavior: true.
 * @test List of size 2, sorted desc. Expected behavior: false.
 * @test List of size 10, list[i] = i^2. Expected behavior: true.
 * @test List of size 10, list[i] = (i-5)^2. Expected behavior: false.
 * @test List of size 10, random elements. Expected behavior: false.
 */
bool kpftimes::isSortedAsc(const DoubleVec &list) {
	if (list.size() <= 1) {
		return true;
	}
	
	bool foundError = false;
	DoubleVec::const_iterator v1 = list.begin();
	DoubleVec::const_iterator v2 = v1 + 1;
	while(!foundError && v2 != list.end()) {
		if (*v1 > *v2) {
			foundError = true;
		}
		else {
			v1++;
			v2++;
		}
	}
	
	return !foundError;
}
