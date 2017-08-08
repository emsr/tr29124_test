/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_experfc test_experfc.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_experfc > test_experfc.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_experfc test_experfc.cpp -lquadmath -lmpfr
./test_experfc > test_experfc.txt
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>
#include <bits/float128_io.h>
#include <mpreal.h>

  template<typename _Tp>
    _Tp
    __experfc_func(_Tp __x)
    {
      mpfr::mpreal __X((long double)__x, 256);
      return (long double)(mpfr::exp(__X * __X) * mpfr::erfc(__X));
    }

  /**
   * The experf function is defined by
   * @f[
   *   experf(x) = exp(x^2)erf(x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __experf_series(_Tp __x)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};
      const auto _S_max_iter = 200;
      const auto __x2 = _Tp{2} * __x * __x;
      auto __term = __x;
      auto __sum = __term;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  __term *= __x2 / _Tp(2 * __k + 1);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return  _Tp{2} * __sum / _S_sqrt_pi;
    }

  /**
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x^2)(1 - erf(x))
   * @f]
   */
  template<typename _Tp>
    _Tp
    __experfc_func_bad(_Tp __x)
    { return std::exp(__x * __x) - __experf_series(__x); }

  /**
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x^2)(1 - erf(x))
   * @f]
   */
  template<typename _Tp>
    _Tp
    __experfc_series(_Tp __x)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};
      const auto _S_max_iter = 200;
      auto __gamr_e = _Tp{1};
      auto __gamr_o = _Tp{2} / _S_sqrt_pi;
      auto __term = _Tp{0};
      auto __sum = _Tp{0};
      auto __xk = _Tp{1};

      for (int __k = 0; __k < _S_max_iter; ++__k)
	{
	  if (__k & 1)
	    {
	      __term = __gamr_o * __xk;
	      __sum += __term;
	      __gamr_o /= _Tp(2 + __k) / _Tp{2};
	    }
	  else
	    {
	      __term = __gamr_e * __xk;
	      __sum += __term;
	      __gamr_e /= _Tp(1 + __k / 2);
	    }
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	  __xk *= -__x;
	}
      return __sum;
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
      const auto _S_pi =  _Tp{3.1415926535897932384626433832795029Q};
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};
      auto __experfc = _S_sqrt_pi * std::sqrt(__x)
		     + std::sqrt(_S_pi * (__x + _Tp{2})
			       - _Tp{2} * (_S_pi - _Tp{2})
			       * std::exp(-std::sqrt(_Tp{5} * __x / _Tp{7})));
      return _Tp{2} / __experfc;
    }

  /**
   * Return scaled repeated integrals of the erfc function by asymptotic series.
   * The experfc function is defined by
   * @f[
   *   experfc(k, x) = exp(x^2)I^k erfc(x)
   * @f]
   * where the integral of the comlementary error function is:
   * @f[
   *   I^k erfc(x) = \frac{2}{k!\sqrt{\pi}}
   *       \int_{x}^{\infty}(t-x)^ke^{-t^2}dt
   * @f]
   * @see Cuyt, et.al. 13.3.2
   */
  template<typename _Tp>
    _Tp
    __experfc_asymp(int __k, _Tp __x)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(__x);
      const auto _S_max_iter = 200;
      const auto __2x = _Tp{2} * __x;
      const auto __2xm2 = -_Tp{1} / (__2x * __2x);
      auto __term = _Tp{1};
      auto __sum = __term;
      auto __kfact = std::__detail::__factorial<_Tp>(__k);
      auto __prev_term = std::abs(__term);
      for (int __m = 1; __m < _S_max_iter; ++__m)
	{
	  __term *= __2xm2 * _Tp(__k + 2 * __m) * _Tp(__k + 2 * __m - 1)
		  / _Tp(__m) / __kfact;
	  if (std::abs(__term) > __prev_term)
	    break;
	  __prev_term = std::abs(__term);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      const auto __fact = _Tp{2} / std::pow(__2x, _Tp(__k + 1)) / _S_sqrt_pi;
      return __fact * __sum;
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
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};
      const auto _S_max_iter = 200;
      const auto __xm2 = -_Tp{1} / (_Tp{2} * __x * __x);
      auto __term = _Tp{1} / __x;
      auto __sum = __term;
      auto __prev_term = std::abs(__term);
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  __term *= __xm2 * _Tp(2 * __k - 1);
	  __sum += __term;
	  if (std::abs(__term) > __prev_term)
	    break;
	  __prev_term = std::abs(__term);
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
      const auto _S_sqrt_2 = _Tp{1.414213562373095048801688724209698078569Q};
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};
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
    __experfc_series_aw(std::complex<_Tp> __z,
			int __trunc_lo = 32, int __trunc_hi = 193,
			_Tp __sep = _Tp{8})
    {
      const auto _S_pi =  _Tp{3.1415926535897932384626433832795029Q};
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};
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
		  auto __s3 = std::exp(-__enn * __enn / _Tp{4})
			    / (__enn * __enn + _Tp{4} * __x * __x);
                  auto __s4 = _Tp{2} * __x * __k1 * __k2
			    - (__x + _S_i * __enn / _Tp{2})
			     * std::exp(-__y * (__enn + __y))
			    - (__x - _S_i * __enn / _Tp{2})
			     * std::exp(+__y * (__enn - __y));
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
    __erfc_scaled(_Tp __x)
    {
      const auto _S_cfrac = _Tp{0.025} * std::numeric_limits<_Tp>::digits;
      // The asymptotic series gets good by here but never really beats C.F.
      //const auto _S_asymp = _Tp{0.18} * std::numeric_limits<_Tp>::digits;
      if (__x < _S_cfrac)
	return __experfc_series(__x);
      else
	return __experfc_cont_frac(__x);
    }

/**
 * Test the scaled complementary error function - experfc(x) = exp(x^2)erfc(x).
 */
template<typename _Tp>
  void
  test_experfc()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    decltype(std::cout.precision()) xw = 22;
    auto w = std::max(xw, 8 + std::cout.precision());

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"x\""
	      << ' ' << std::setw(w) << "\"approx experfc(x)\""
	      << ' ' << std::setw(w) << "\"exp(x^2)erfc(x)\""
	      << ' ' << std::setw(w) << "\"cfrac experfc(x)\""
	      << ' ' << std::setw(w) << "\"asymp experfc(x)\""
	      << ' ' << std::setw(w) << "\"series experfc(x)\""
	      << ' ' << std::setw(w) << "\"delta cfrac\""
	      << ' ' << std::setw(w) << "\"delta asymp\""
	      << ' ' << std::setw(w) << "\"delta series\""
	      << '\n';
    for (int __i = 0; __i <= 5500; ++__i)
      {
	auto __x = __i * _Tp{0.01Q};
	auto __experfc_apr = __experfc_approx(__x * __x);
	auto __experfc_fun = __experfc_func(__x);
	auto __experfc_cfr = __experfc_cont_frac(__x);
	auto __experfc_asy = __experfc_asymp(__x);
	auto __experfc_ser = __experfc_series(__x);
	std::cout << ' ' << std::setw(w) << __x
		  << ' ' << std::setw(w) << __experfc_apr
		  << ' ' << std::setw(w) << __experfc_fun
		  << ' ' << std::setw(w) << __experfc_cfr
		  << ' ' << std::setw(w) << __experfc_asy
		  << ' ' << std::setw(w) << __experfc_ser
		  << ' ' << std::setw(w) << (__experfc_cfr - __experfc_fun) / __experfc_fun
		  << ' ' << std::setw(w) << (__experfc_asy - __experfc_fun) / __experfc_fun
		  << ' ' << std::setw(w) << (__experfc_ser - __experfc_fun) / __experfc_fun
		  << '\n';
      }
  }

int
main()
{
  std::cout << "\n\n  float\n";
  std::cout << "  =====\n";
  test_experfc<float>();

  std::cout << "\n\n  double\n";
  std::cout << "  ======\n";
  test_experfc<double>();

  std::cout << "\n\n  long double\n";
  std::cout << "  ===========\n";
  test_experfc<long double>();

  std::cout << "\n\n  __float128\n";
  std::cout << "  ==========\n";
  test_experfc<__float128>();
}
