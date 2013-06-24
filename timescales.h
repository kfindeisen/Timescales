/** Primary header for Krzysztof's timescales library.
 *  @file timescales.h
 *  @author Krzysztof Findeisen
 *  @date Created January 25, 2010
 *  @date Last modified June 21, 2013
 *
 * @todo Fix any bugs related to the \@overload documentation tag
 */
 
/** @mainpage
 *
 * @section intro Overview
 * 
 * The Timescales library provides basic functions for lightcurve analysis. 
 * The target application is automated reduction of large time-series data sets. 
 *
 * The library is organized as a series of global functions under the 
 * @c kpftimes namespace, rather than as an object heirarchy. This architecture 
 * was chosen for its simplicity: each function performs a single, narrowly 
 * defined task, making it (hopefully!) easy to chain functions together 
 * into pipelines.
 * 
 * @section metahelp About this Documentation
 * 
 * New users will find the @htmlonly Modules tab at the top of this page @endhtmlonly 
 * @latexonly Module Documentation chapter @endlatexonly the best 
 * starting point for learning about the Timescales API. There they will find 
 * a list of the main functions in the library, organized by category. The 
 * other @htmlonly tabs @endhtmlonly @latexonly chapters @endlatexonly 
 * are more useful for people seeking to understand the code itself.
 *
 * @section credits Credits
 *
 * Lightcurve MC was primarily written by Krzysztof Findeisen. Please contact 
 * him at krzys, astro caltech edu for questions, feedback, or bug reports.
 * 
 * @page install Installation
 *
 * @page install Installation
 *
 * @section install_reqs Requirements
 *
 * Timescales is written in standard C++, and should build on any 
 * C++98/C++03-compliant compiler. For portability between build 
 * environments, no TR1 or C++11 features are included. However, the 
 * makefile itself is specific to GCC, so users wishing to install 
 * Timescales on a non-Unix platform may need to revise the contents of 
 * the following files: 
 * - @c makefile
 * - @c makefile.inc
 * - @c makefile.common
 * - @c makefile.subdirs
 * 
 * The library may switch to CMake in the future for improved portability.
 *
 * Timescales depends on the following external libraries:
 * - <a href="http://www.boost.org/">Boost</a> 1.31 or later
 * - <a href="http://www.gnu.org/software/gsl/">Gnu Scientific Library (GSL)</a> 
 *	1.3 or later
 * 
 * In addition, (re-)generating this documentation requires 
 * <a href="http://www.doxygen.org/">Doxygen</a> 1.8.1 or later.
 *
 * @section install_make Build Commands
 *
 * Running <tt>make</tt>, with no arguments, will build @c libtimescales.a in 
 * the @c timescales directory. Running <tt>make unittest</tt> will build the 
 * unit test suite at @c tests/test. Running <tt>make autotest</tt> will 
 * build the library and the test suite, if neccessary, before running all 
 * tests. Finally, running <tt>make doc</tt> will 
 * generate this documentation at @c doc/html/index.html and 
 * (if you have LaTeX installed) @c doc/latex/refman.pdf, and running 
 * <tt>make all</tt> will generate the library, the test suite, and the 
 * documentation.
 *
 * Place @c timescales.h in a directory where your C++ compiler can find 
 * header files, such as <tt>/usr/local/include/</tt>. 
 * Place @c libtimescales.a in a directory where your C++ compiler can find 
 * static libraries, such as <tt>/usr/local/lib</tt>. A shared version of the 
 * Timescales library is not available at present.
 * 
 * @section install_use Using Timescales in Your Software
 * 
 * To use Timescales, include <tt>\<timescales.h\></tt> in your source code 
 * (see @c examples/example.cpp). 
 * You must link the final program with @c -lgsl as well as @c -ltimescales. 
 * The library depends on Boost only through template headers, so you do not 
 * need to link against any Boost libraries.
 * 
 * @page use User's Guide
 * 
 * @page changelog Version History
 * 
 * Timescales conforms to 
 * <a href="http://semver.org/spec/v2.0.0-rc.1.html">version 2.0.0-rc.1 of the Semantic Versioning specification</a>. 
 * All version numbers are to be interpreted as described therein. 
 * This documentation constitutes the public API for the library.
 *
 * @section v04 0.4.0-devel
 *
 * @subsection v0_4_0_diff Changes 
 * 
 * - From now on, version numbers will conform to 
 *	<a href="http://semver.org/spec/v2.0.0-rc.1.html">version 2.0.0-rc.1 of the Semantic Versioning specification</a>.
 * - Removed statistics routines from utils.h to common/stats.tmp.h
 * - Improved quality of documentation
 * - Improved quality of source code. All functions are now exception-safe.
 * - Removed the SCARGLE_SLOW compile-time option
 * - Wrote more informative exception strings
 * 
 * @subsection v0_4_0_new New Features 
 * 
 * 
 * @subsection v0_4_0_fix Bug Fixes 
 * 
 * - autoCorr() and acWindow() now correctly test whether the input time 
 *	lag grid is uniform
 * 
 * @section v03 0.3
 *
 * No version information recorded.
 *
 * @section v022 0.2.2
 *
 * @subsection v0_2_2_diff Changes
 * 
 * - Removed C++98-style exception specifications from all functions
 * 
 * @subsection v0_2_2_new New Features 
 * 
 * - Added delMDelT()
 *
 * @section v021 0.2.1
 *
 * @subsection v0_2_1_new New Features 
 * 
 * - Added acWindow()
 * - Added lsNormalEdf()
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

/** @} */	// end Periodogram generation

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

/** @} */	// end Autocorrelation function generation

//----------------------------------------------------------
/** @defgroup dmdt &Delta;m&Delta;t pair diagram generation
 *  @{
 */

/** Calculates a &Delta;m&Delta;t plot.
 */
void dmdt(const DoubleVec &times, const DoubleVec &fluxes, 
		DoubleVec &deltaT, DoubleVec &deltaM);
 
/** @} */	// end &Delta;m&Delta;t generation

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

/** @} */	// end Frequency/offset grid generation

}		// end kpftimes

#endif		// TIMESCALEH
