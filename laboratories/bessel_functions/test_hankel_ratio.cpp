/**
 *
 */

#include <complex>
#include <iostream>
#include <iomanip>

#include <ext/continued_fractions.h>

/**
 * Return the ratio of the Hankel functions
 * @f[
 *   \frac{H_{\nu}^{(1)}(z)}{H_{\nu+1}^{(1)}(z)}
 *   = \frac{\nu}{z} - \frac{H'_{\nu}^{(1)}(z)}{H_{\nu}^{(1)}(z)}
 *   = 
 * @f]
 */
template<typename _Tnu, typename _Tp>
  _Tp
  __cyl_hankel_ratio(_Tnu __nu, _Tp __x)
  {
    using _Cmplx = std::complex<_Tp>;

    auto a_trigint
      = [](std::size_t i, _Tp)
	-> _Cmplx
	{
	  if (i == 1)
	    return _Tp(1);
	  else
	    return -_Tp(i - 1) * _Tp(i - 1);
	};
    using _AFun = decltype(a_trigint);

    auto b_trigint
      = [](std::size_t i, _Tp __x)
	{
	  if (i == 0)
	    return _Cmplx{0};
	  else
	    return _Cmplx{_Tp(2 * i - 1), __x};
	};
    using _BFun = decltype(b_trigint);

    auto w_trigint
      = [](std::size_t, _Tp)
	{ return _Cmplx{}; };
    using _TailFun = decltype(w_trigint);

    _SteedContinuedFraction<_Tp, _AFun, _BFun, _TailFun>
      SiCi(a_trigint, b_trigint, w_trigint);
  }

int
main()
{
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

}
