/** Support code for the library. None of these routines are intended as part of the public API.
 * @file utils.h
 * @author Krzysztof Findeisen
 * @date Created April 13, 2011
 * @date Last modified July 24, 2011
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
private:
	double *table;
	const size_t dimX, dimY;
};

/** Tests whether the argument is sorted.
 */
bool isSortedAsc(const DoubleVec &list);

/*----------------------------------------------------------
 * Speed-optimized statistics routines
 * I already had these written, so it was simpler than writing a wrapper to 
 *	convert vectors to C arrays just so GSL can understand them
 */

/** Finds the mean of the values in a generic container object. The container 
 *	class is accessed using first and last iterators, as in the C++ STL 
 *	convention, and the variance is computed over the interval [first, last).
 * 
 * @tparam ForwardIterator The iterator type for the container over which the 
 *	mean is to be calculated
 * @param[in] first Forward iterator marking the first element in the container.
 * @param[out] last Forward iterator marking the position after the last 
 *	element in the container.
 *
 * @return The arithmetic mean of the elements between first, inclusive, and 
 *	last, exclusive. The return type is that of the elements pointed to by the first and 
 *	last iterators.
 *
 * @pre first is "before" last in the sense that incrementing first repeatedly 
 *	would reach last.
 * @pre There is at least one element in the interval [first, last)
 * @exception invalid_argument Thrown if there are not enough elements.
 *
 * @test List of ints, length 0. Expected behavior: throw invalid_argument.
 * @test List of ints, length 1. Expected behavior: return list[0]
 * @test List of ints, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_mean in 10 out of 10 trials.
 * @test List of doubles, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_mean to within 1e-10 in 10 out of 10 trials.
 * @test Vector of doubles, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_mean to within 1e-10 in 10 out of 10 trials.
 * @test Array of doubles, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_mean to within 1e-10 in 10 out of 10 trials.
 */
template <class ForwardIterator> 				// Iterator to use
		typename std::iterator_traits<ForwardIterator>::value_type 	// Type pointed to by iterator
		mean(ForwardIterator first, ForwardIterator last) {
	typedef typename std::iterator_traits<ForwardIterator>::value_type Value;

	Value  sum = 0;
	long count = 0;

	for(; first != last; first++) {
		sum += (*first);
		count++;
	}
	if (count <= 0) {
		throw std::invalid_argument("Not enough data to compute mean");
	}
	return static_cast<Value>(sum / count);
}

/** Finds the variance of the values in a generic container object. The 
 *	container class is accessed using first and last iterators, as in the 
 *	C++ STL convention, and the mean is computed over the interval 
 *	[first, last).
 * 
 * @tparam ForwardIterator The iterator type for the container over which the 
 *	variance is to be calculated
 * @param[in] first Forward iterator marking the first element in the container.
 * @param[out] last Forward iterator marking the position after the last element in the container.
 *
 * @return The (unbiased) sample variance of the elements between first, 
 *	inclusive, and last, exclusive. The return type is that of the 
 *	elements pointed to by the first and last iterators.
 *
 * @pre first is "before" last in the sense that incrementing first repeatedly 
 *	would reach last.
 * @pre There is at least two elements in the interval [first, last)
 * @exception invalid_argument Thrown if there are not enough elements.
 *
 * @test List of ints, length 0. Expected behavior: throw invalid_argument.
 * @test List of ints, length 1. Expected behavior: throw invalid_argument.
 * @test List of ints, length 2. Expected behavior: return (list::back()-list::front())^2/2
 * @test List of ints, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_variance in 10 out of 10 trials.
 * @test List of doubles, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_variance to within 1e-10 in 10 out of 10 trials.
 * @test Vector of doubles, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_variance to within 1e-10 in 10 out of 10 trials.
 * @test Array of doubles, length 100, randomly generated. Expected behavior: 
 *	agrees with gsl_stats_variance to within 1e-10 in 10 out of 10 trials.
 */
template <class ForwardIterator> 				// Iterator to use
		typename std::iterator_traits<ForwardIterator>::value_type 	// Type pointed to by iterator 
		variance(ForwardIterator first, ForwardIterator last) {
	typedef typename std::iterator_traits<ForwardIterator>::value_type Value;

	Value  sum = 0, sumsq = 0;
	long count = 0;

	for(; first != last; first++) {
		sum   += (*first);
		sumsq += (*first)*(*first);
		count++;
	}
	if (count <= 1) {
		throw std::invalid_argument("Not enough data to compute variance");
	}

	// Minimize number of divisions and maximize dividend in case value_type is integral
	return static_cast<Value>((sumsq - sum*sum/count)/(count-1));
}

/** @} */

}
