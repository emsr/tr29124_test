#ifndef ASIN_CPP
#define ASIN_CPP

#include "asin.h"
#include <cmath>

using namespace std ;

/** Complex valued arc sine
* The need for an implementation of this mathematical function arises from
* the fact that g++ (as with 4.1.2) does not yet provide this functionality.
* The alternatives are to use the boost library or the AIX compiler, for
* example.
* @param[in] z complex parameter, the argument
* @return the value of the arc sine on the principal branch
* @author Richard J. Mathar
*/
template <class _FLT> complex<_FLT> asin (const complex<_FLT> & z)
{
#if 0
	_FLT x = real(z) ;
	_FLT y = imag(z) ;
	// Abramowitz-Stegun, 4.4.39
	_FLT alpha = 0.5*(hypot(x+1.,y)+hypot(x-1.,y)) ;
	_FLT beta = 0.5*(hypot(x+1.,y)-hypot(x-1.,y)) ;
	// Abramowitz-Stegun, 4.4.37
	return complex<_FLT> (asin (beta), log (alpha+sqrt(alpha*alpha-1.) ));
#else
	return complex<_FLT>(0.,-1.)*log(complex<_FLT>(0.,1.)*z+sqrt(complex<_FLT>(1.,0.)-z*z)) ;
#endif
}
#if 0
template complex<float> asin (const complex<float> z) ;
template complex<double> asin (const complex<double> z) ;
template complex<long double> asin (const complex<long double> z) ;
#endif

#endif /* ASIN_CPP */
