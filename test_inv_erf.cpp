// $HOME/bin_specfun/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_inv_erf test_inv_erf.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_inv_erf > test_inv_erf.txt

// g++ -std=c++1z -o test_inv_erf test_inv_erf.cpp

// ./test_inv_erf.exe > test_inv_erf.txt

// AAOF pp. 408-409.

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>
#include "float128.h"

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __erf_series(_Tp __x)
    {
      
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __experf_series(_Tp __x)
    {
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};
      const auto _S_max_iter = 200;
      const auto __x2 = __x * __x;
      auto __term = __x;
      auto __sum = __term;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  __term *= _Tp{2} * __x2 / _Tp(2 * __k + 1);
	  __sum += __term;
	}
      return  _Tp{2} * __sum / _S_sqrt_pi;
    }


  /**
   * Return a quick approximation to the experfc function.
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x)erfc(\sqrt{x})
   * @f]
   * The approximation is:
   * @f[
   *   experfc(x) = \frac{2}{\sqrt{\pi x}
   *       + \sqrt{\pi(x+2)-2(\pi-2)exp(-\sqrt(5x/7))}}
   * @f]
   * Surprisingly, this is accurate to within 0.1% over the whole range
   * [0, infty].  It is used to start agm algorithms of the experfc function.
   */
  template<typename _Tp>
    _Tp
    __experfc_approx(_Tp __x)
    {
      const auto _S_pi =  _Tp{3.1415926535897932384626433832795029L};
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};
      auto __experfc = _S_sqrt_pi * std::sqrt(__x)
		     + std::sqrt(_S_pi * (__x + _Tp{2})
			       - _Tp{2} * (_S_pi - _Tp{2})
			       * std::exp(-std::sqrt(_Tp{5} * __x / _Tp{7})));
      return _Tp{2} / __experfc;
    }

  /**
   * Return the experfc function by asymptotic series.
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x^2)erfc(x)
   * @f]
   * The asymptotic series is:
   * @f[
   *   experfc(x) = \frac{1}{\sqrt{\pi} x}
   *       \sum_{k=0}^{\infty}\frac{(2k-1)!!}{(-2x^2)^k}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __experfc_asymp(_Tp __x)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};
      const auto _S_max_iter = 200;
      const auto __xm2 = -_Tp{1} / (_Tp{2} * __x * __x);
      auto __term = _Tp{1} / __x;
      auto __sum = __term;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  __term *= __xm2 * _Tp(2 * __k - 1);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return __sum / _S_sqrt_pi;
    }

  /**
   * Return the experfc function by continued fraction.
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x^2)erfc(x)
   * @f]
   * The continued fraction is
   * @f[
   *   experfc(x) = \cfrac{\sqrt{2}\pi}{\sqrt{2}x
   *              + \cfrac{1}{\sqrt{2}x
   *              + \cfrac{2}{\sqrt{2}x
   *              + \cfrac{3}{\sqrt{2}x + ...} } } }
   * @f]
   */
  template<typename _Tp>
    _Tp
    __experfc_cont_frac(_Tp __x)
    {
      const auto _S_sqrt_2 = _Tp{1.414213562373095048801688724209698078569L};
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};
      const auto __b = std::sqrt(_Tp{2} * __x * __x);
      auto __s = __b;
      for (int __n = 100; __n >= 1; --__n)
	{
	  auto __a = _Tp(__n);
	  __s = __b + __a / __s;
	}
	__s = (_S_sqrt_2 / _S_sqrt_pi) / __s;
      return __s;
    }

  /**
   * exp(x^2) erfc(x).
   */
  template<typename _Tp>
    std::complex<_Tp>
    __experfc_series_wanker(std::complex<_Tp> __z,
		     int __trunc_lo = 32, int __trunc_hi = 193,
		     _Tp __sep = _Tp{8})
    {
      const auto _S_pi =  _Tp{3.1415926535897932384626433832795029L};
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};
      const auto _S_i = std::complex<_Tp>{0, 1};

      if (__sep < _Tp{0} || std::abs(__z) <= __sep)
        {
          // Infinite series approximation, Abramowitz p. 313 7.1.29
          auto __x = std::real(__z);
          auto __y = std::imag(__z);

          auto __ez2 = std::exp(__z * __z);

          auto __s1 = erf(__x);

          auto __k1 = std::exp(-__y * __y);
          auto __k2 = std::exp(_Tp{2} * _S_i * __x * __y);

          std::complex<_Tp> __s2;
          if (__x == _Tp{0})
            __s2 = _S_i * __y * __k1 / _S_pi;
          else
            __s2 = __k1 * (__k2 - _Tp{1} / (_Tp{2} * _S_pi * __x));

          auto __retval = __ez2 * __s1 + __s2;

          if (__y != _Tp{0})
	    {
              auto __sum = std::complex<_Tp>{};
              for (unsigned __n = 1; __n <= __trunc_lo; ++__n)
		{
		  auto __enn = _Tp(__n);
		  auto __s3 = std::exp(-__enn * __enn / _Tp{4}) / (__enn * __enn + _Tp{4} * __x * __x);
                  auto __s4 = _Tp{2} * __x * __k1 * __k2
                      - (__x + _S_i * __enn / _Tp{2}) * std::exp(-__y * (__enn + __y))
                      - (__x - _S_i * __enn / _Tp{2}) * std::exp(+__y * (__enn - __y));
                  __sum += __s3 * __s4;
		}
              __retval += _Tp{2} * __sum / _S_pi;
	    }
          return __ez2 - __retval;
        }
      else
	{
          // Asymptotic expansion, Abramowitz p. 312 7.1.23
          bool __isneg = (std::real(__z) < _Tp{0});
          if (__isneg)
	    __z = -__z;

          std::complex<_Tp> __s = _Tp{1};
          std::complex<_Tp> __y = _Tp{2} * __z * __z;
          for (auto __n = __trunc_hi; __n >= 1; __n -= 2)
	    __s = _Tp{1} - _Tp(__n) * (__s / __y);

          auto __retval = __s / (_S_sqrt_pi * __z);

          if (__isneg)
	    {
              __z = -__z;
              __retval = _Tp{2} - __retval;
            }

          return __retval;
        }
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    __experfc(_Tp __x)
    {
      if (__x < _Tp{3})
	return std::exp(__x * __x) * erfc(__x);
      else
	return __experfc_cont_frac(__x);
    }

  /**
   * Return the inverse error function.
   */
  template<typename _Tp>
    _Tp
    __erfi(_Tp __p)
    {
std::cout.precision(std::numeric_limits<_Tp>::digits10);
auto width = 8 + std::cout.precision();
      constexpr auto _S_eps = _Tp{10} * std::numeric_limits<_Tp>::epsilon();
      // Iterate experfc(x^2).
      if (__p < _Tp{0})
	return -__erfi(-__p);
      else
	{
	  auto __x2 = _Tp{25};
	  auto __x2prev2 = _Tp{0}, __x2prev = _Tp{0};
	  const auto _S_max_iter = 500;
	  auto __iter = 0;
	  while (++__iter < _S_max_iter)
	    {
		    
	      __x2prev2 = __x2prev;
	      __x2prev = __x2;
	      __x2 = std::log(__experfc(__x2) / (_Tp{1} - __p));
//std::cout
// << ' ' << std::setw(width) << __x2
// << ' ' << std::setw(width) << __x2 - __x2prev
// << ' ' << std::setw(width) << __x2 - __x2prev2
// << '\n';
	      // If the fraction jumps < 0 just bop it back.
	      if (__x2 < _Tp{0})
		__x2 = -__x2;
	      if (_S_eps > std::abs(__x2 - __x2prev) / std::abs(__x2))
		break;
	      if (_S_eps > std::abs(__x2 - __x2prev2) / std::abs(__x2))
		break;
	    }
	  return std::sqrt(__x2);
	}
    }

/**
 * Test the inverse error function.
 */
template<typename _Tp>
  void
  test_inv_erf()
  {
    //  Build the series coefficients.
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_pi =  _Tp{3.1415926535897932384626433832795029L};
    const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = 8 + std::cout.precision();

    const int n_max = 250;
    std::vector<_Tp> a;
    a.push_back(1);
    for (int __n = 1; __n < n_max; ++__n)
      {
	auto __atemp = _Tp{0};
	for (int __k = 1; __k <= __n; ++__k)
	  __atemp += _Tp(2 * (__k - 1) + 1) * a[__k - 1]
		 * _Tp(2 * (__n - __k) + 1) * a[__n - __k]
		 / _Tp(__k * (2 * __k - 1));
	__atemp /= _Tp(2 * __n + 1);
	a.push_back(__atemp);
      }

    std::cout << "\n\n " << std::setw(width) << "a_k" << '\n';
    for (auto __aa : a)
      std::cout << ' ' << std::setw(width) << __aa << '\n';

    std::vector<_Tp> c;
    c.push_back(1);
    for (int __n = 1; __n < n_max; ++__n)
      {
	auto __ctemp = _Tp{0};
	for (int __k = 1; __k <= __n; ++__k)
	  __ctemp += c[__k - 1] * c[__n - __k] / _Tp(__k * (2 * __k - 1));
	c.push_back(__ctemp);
      }
    for (int __n = 1; __n < n_max; ++__n)
      c[__n] /= _Tp(2 * __n + 1);

    std::cout << "\n\n " << std::setw(width) << "c_k" << '\n';
    for (auto __cc : c)
      std::cout << ' ' << std::setw(width) << __cc << '\n';

    std::cout << "\n\n"
	      << std::setw(width) << "p"
	      << std::setw(width) << "inv_erf(p)"
	      << std::setw(width) << "erf(inv_erf(p))"
	      << std::setw(width) << "erf(inv_erf(p)) - p"
	      << '\n';
    for (int __i = -100; __i <= 100; ++__i)
      {
	auto __x = __i * _Tp{0.01L};
	auto __chi = _S_sqrt_pi * __x / _Tp{2};
	auto __chi2 = __chi * __chi;
	auto __chip = __chi;
	auto __inverf = _Tp{0};
	for (int __k = 0; __k < n_max; ++__k)
	  {
	    auto __term = a[__k] * __chip;
	    __inverf += __term;
	    if (std::abs(__term) < _S_eps * std::abs(__inverf))
	      break;
	    __chip *= __chi2;
	  }
	std::cout << ' ' << std::setw(width) << __x
		  << ' ' << std::setw(width) << __inverf
		  << ' ' << std::setw(width) << erf(__inverf)
		  << ' ' << std::setw(width) << erf(__inverf) - __x
		  << '\n';
      }

    std::cout << "\n\n"
	      << std::setw(width) << "x"
	      << std::setw(width) << "erf(x)"
	      << std::setw(width) << "inv_erf(erf(x))"
	      << std::setw(width) << "inv_erf(erf(x)) - x"
	      << std::setw(width) << "inv_erf(erf(x))"
	      << std::setw(width) << "inv_erf(erf(x)) - x"
	      << '\n';
    for (int __i = -200; __i <= 200; ++__i)
      {
	auto __x = __i * _Tp{0.01L};
	auto __erfx = erf(__x);
	auto __chi = _S_sqrt_pi * __erfx / _Tp{2};
	auto __chi2 = __chi * __chi;
	auto __chip = __chi;
	auto __inverf = _Tp{0};
	for (int __k = 0; __k < n_max; ++__k)
	  {
	    auto __term = a[__k] * __chip;
	    __inverf += __term;
	    if (std::abs(__term) < _S_eps * std::abs(__inverf))
	      break;
	    __chip *= __chi2;
	  }
	std::cout << ' ' << std::setw(width) << __x
		  << ' ' << std::setw(width) << __erfx
		  << ' ' << std::setw(width) << __inverf
		  << ' ' << std::setw(width) << __inverf - __x
		  << ' ' << std::setw(width) << __erfi(__erfx)
		  << ' ' << std::setw(width) << __erfi(__erfx) - __x
		  << '\n';
      }

    std::cout << "\n\n"
	      << std::setw(width) << "x"
	      << std::setw(width) << "approx experfc(x)"
	      << std::setw(width) << "exp(x)erfc(rt(x))"
	      << std::setw(width) << "cfrac experfc(x)"
	      << std::setw(width) << "asymp experfc(x)"
	      << '\n';
    for (int __i = 0; __i <= 5500; ++__i)
      {
	auto __x = __i * _Tp{0.01L};
	std::cout << ' ' << std::setw(width) << __x
		  << ' ' << std::setw(width) << __experfc_approx(__x * __x)
		  << ' ' << std::setw(width) << std::exp(__x * __x) * erfc(__x)
		  << ' ' << std::setw(width) << __experfc_cont_frac(__x)
		  << ' ' << std::setw(width) << __experfc_asymp(__x)
		  << '\n';
      }
  }

int
main()
{
  std::cout << "\n\n  float\n";
  std::cout << "  =====\n";
  test_inv_erf<float>();

  std::cout << "\n\n  double\n";
  std::cout << "  ======\n";
  test_inv_erf<double>();

  std::cout << "\n\n  long double\n";
  std::cout << "  ===========\n";
  test_inv_erf<long double>();

  //std::cout << "\n\n  __float128\n";
  //std::cout << "  ==========\n";
  //test_inv_erf<__float128>();
}
