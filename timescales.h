/** Primary header for Krzysztof's timescales library.
 *  @file timescales.h
 *  @author Krzysztof Findeisen
 *  @date Created January 25, 2010
 *  @date Last modified July 24, 2011
 */
 
/** @mainpage
 *
 * @section intro Introduction
 * 
 * The Timescales library provides basic functions for lightcurve analysis. 
 * The target application is automated reduction of large time-series data sets. 
 *
 * The library is organized as a series of global functions under the 
 * kpftimes namespace, rather than as an object heirarchy. This architecture 
 * was chosen for its simplicity: each function performs a single, narrowly 
 * defined task, making it (hopefully!) easy to chain functions together 
 * into pipelines.
 * 
 * @section metahelp About this Documentation
 * 
 * @htmlonly
 * New users will find the Modules tab at the top of this page the best 
 * starting point for learning about the Timescales API. There they will find 
 * a list of the main functions in the library, organized by category. The 
 * other tabs are more useful for people seeking to understand the code 
 * itself.
 * @endhtmlonly
 * 
 * @latexonly
 * New users will find the Module Documentation chapter the best 
 * starting point for learning about the Timescales API. There they will find 
 * a list of the main functions in the library, organized by category. The 
 * other chapters are more useful for people seeking to understand the code 
 * itself.
 * @endlatexonly
 * 
 * @section install Installation and Use
 * 
 * Timescales itself should compile on any UNIX-like system. In many cases 
 * you need simply unpack the .tar contents into the appropriate directory, 
 * run @c make, and move libtimescales.a into a directory of your choice. If 
 * you do not use GCC, you may need to edit the makefile before you can build 
 * the library.
 *
 * To use Timescales, include <tt>\<timescales.h\></tt> in your source code 
 * (see examples/example.cpp). Timescales relies on the GNU Scientific Library 
 * (<a href="http://www.gnu.org/software/gsl/">GSL</a>) for some of its more 
 * complex mathematics, so you must link your program with @e both Timescales 
 * and GSL for it to run correctly. Check with your system administrator if 
 * you're not sure whether GSL is available on your machine.
 * 
 * @section changelog Recent changes
 *
 * @subsection v022 0.2.2
 *
 * @li Removed exception specifications from all functions
 * @li Added delMDelT()
 *
 * @subsection v021 0.2.1
 *
 * @li Added acWindow()
 * @li Added lsNormalEdf()
 * 
 */

#ifndef TIMESCALEH
#define TIMESCALEH

#include <stdexcept>
#include <vector>

/** A convenient shorthand for vectors of doubles.
 */
typedef std::vector<double> DoubleVec;

/** The kpftimes namespace uniquely identifies member functions of the Timescales library.
 */
namespace kpftimes {

//----------------------------------------------------------
/** @defgroup period Periodogram generation
 *  @{
 */
 
/** @def SCARGLE_SLOW 
 * A compiler flag controlling the strategy used for lsThreshold().
 * If set to 0 or not set [the default], lsThreshold() is optimized for speed 
 *	at the expense of memory. If set to 1, it is instead optimized for 
 *	memory use at the expense of speed. The default setting should 
 *	perform better on all but the oldest systems.
 * 
 * This flag is best set in the call to the compiler, rather than in 
 *	the code itself.
 */

/** Calculates the Lomb-Scargle periodogram for a time series.
 */
void lombScargle(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &freq, DoubleVec &power);

/** Calculates the significance threshold for a Lomb-Scargle periodogram.
 */
double lsThreshold(const DoubleVec &times, const DoubleVec &freq, double fap, long nSims);

/** Calculates the empirical distribution function of false peaks for a 
 *	Lomb-Scargle periodogram.
 */
void lsNormalEdf(const DoubleVec &times, const DoubleVec &freqs, 
		DoubleVec &powers, DoubleVec &probs, long nSims);

/** @} */

//----------------------------------------------------------
/** @defgroup acf Autocorrelation function generation
 *  @{
 */
 
/** Calculates the autocorrelation function for a time series. 
 */
void autoCorr(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &offsets, DoubleVec &acf);

/** Calculates the autocorrelation function for a time series. 
 */
void autoCorr(const DoubleVec &times, const DoubleVec &fluxes, 
		const DoubleVec &offsets, DoubleVec &acf, double maxFreq);

/** Calculates the autocorrelation window function for a time sampling. 
 */
void acWindow(const DoubleVec &times, const DoubleVec &offsets, DoubleVec &wf);

/** Calculates the autocorrelation window function for a time sampling. 
 */
void acWindow(const DoubleVec &times, const DoubleVec &offsets, DoubleVec &wf, 
		double maxFreq);

/** @} */

//----------------------------------------------------------
/** @defgroup dmdt Delta-T pair diagram generation
 *  @{
 */

/** Calculates all Delta-m Delta-T pairs
 */
void dmdt(const DoubleVec &times, const DoubleVec &fluxes, 
		DoubleVec &deltaT, DoubleVec &deltaM);
 
/** @} */

//----------------------------------------------------------
/** @defgroup grid Frequency/offset grid generation
 *  @{
 */
 
/** Returns the time interval covered by the data.
 */
double deltaT(const DoubleVec &times);

/** Returns the pseudo-Nyquist frequency for a grid of observations.
 */
double pseudoNyquistFreq(const DoubleVec &times);

/** Returns the highest frequency that can be probed by the data.
 */
double maxFreq(const DoubleVec &times);

/** Creates a frequency grid that can be fed to time series analysis 
 * functions. 
 */
void freqGen(const DoubleVec &times, DoubleVec &freq);

/** @overload 
 */
void freqGen(const DoubleVec &times, DoubleVec &freq, 
		double fMin, double fMax);

/** @overload 
 */
void freqGen(const DoubleVec &times, DoubleVec &freq, 
		double fMin, double fMax, double fStep);
}

/** @} */

#endif		// defined(TIMESCALEH)
