/** Computes the irregularly-sampled discrete Fourier transform
 * @file dft.h
 * @author Krzysztof Findeisen
 * @date Created February 13, 2011
 * @date Last modified April 14, 2011
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
