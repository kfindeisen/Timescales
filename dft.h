/** Computes the irregularly-sampled discrete Fourier transform
 * @file dft.h
 * @author Krzysztof Findeisen
 * @date Created February 13, 2011
 * @date Last modified April 14, 2011
 */ 

#ifndef SCARGDFTH
#define SCARGDFTH

#include <complex>
#include <vector>

/** A convenient shorthand for vectors of doubles.
 */
typedef std::vector<double               >  DoubleVec;
/** A convenient shorthand for vectors of complex numbers.
 */
typedef std::vector<std::complex<double> > ComplexVec;

namespace kpftimes {

/** Calculates the discrete Fourier transform for a list of times and fluxes
 * @ingroup util
 */
void dft(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &freqs, ComplexVec &dft);

}	// end kpftimes::

#endif
