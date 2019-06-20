/**
 *
 */

#include <ext/math_constants.h>
#include <ext/complex_util.h>
#include <bits/float128_io.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <cmath>

#include <bits/notsospecfun.h> // For complex log1p.
#include <bits/numeric_limits.h>

#include <wrap_gsl.h>

  /**
   *   @brief Return the hypergeometric function @f$ _2F_1(a,b;c;x) @f$
   *   by series expansion.
   *
   *   The hypergeometric function is defined by
   *   @f[
   *     _2F_1(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
   *                      \sum_{n=0}^{\infty}
   *                      \frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)}
   *                      \frac{x^n}{n!}
   *   @f]
   *
   *   This works and it's pretty fast.
   *
   *   @param  __a  The first @a numerator parameter.
   *   @param  __b  The second @a numerator parameter.
   *   @param  __c  The @a denominator parameter.
   *   @param  __x  The argument of the confluent hypergeometric function.
   *   @return  The confluent hypergeometric function.
   */
  template<typename _Tp>
    _Tp
    __hyperg_series(_Tp __a, _Tp __b, _Tp __c, _Tp __x)
    {
      using _Val = std::__detail::__num_traits_t<_Tp>;
      const auto __eps = __gnu_cxx::__epsilon<_Val>();
      const unsigned int _S_max_iter = 100000u;
      auto __aint = __gnu_cxx::__fp_is_integer(__a);
      auto __bint = __gnu_cxx::__fp_is_integer(__b);
      auto __max_iter = _S_max_iter;
      if (__aint && __aint() < 0 && int(_S_max_iter) > -__aint())
	__max_iter = -__aint();
      if (__bint && __bint() < 0 && int(_S_max_iter) > -__bint())
	__max_iter = -__bint();

      auto __term = _Tp{1};
      auto __Fabc = _Tp{1};
      unsigned int __i;
      for (__i = 0; __i < __max_iter; ++__i)
	{
	  __term *= (__a + _Tp(__i)) * (__b + _Tp(__i)) * __x
		  / ((__c + _Tp(__i)) * _Tp(1 + __i));
	  if (std::abs(__term) < __eps)
	    break;
	  __Fabc += __term;
	}
      if (__i == __max_iter)
	std::__throw_runtime_error(__N("Series failed to converge "
				       "in __hyperg_series."));

      return __Fabc;
    }

  /**
   * Do Buhring's analytic continuation.
   * @param __s  The numerator parameter, (a or b)
   * @todo This could be written to grow an existing array.
   */
  template<typename _Tp>
    std::vector<_Tp>
    __hyperg_buhring_d(int __ss, _Tp __a, _Tp __b, _Tp __c, _Tp __z0)
    {
      const unsigned int _N = 20; // ?
      std::vector<_Tp> __d(_N);
      const auto __s = __ss == 0 ? __a : __b;
      unsigned int __n = 1;
      auto __danm2 = _Tp{0};
      auto __danm1 = _Tp{1};
      __d[0] = __danm1;
      auto __dan = (_Tp(__n - 1) + __s)
		 * (__z0 * (_Tp{1} - __z0) * (_Tp(__n - 2) + __s) * __danm2
		  + ((_Tp(__n) + __s) * (_Tp{1} - _Tp{2} * __z0)
		   + (__a + __b + _Tp{1}) * __z0 - __c) * __danm1)
		 / _Tp(__n) / (_Tp(__n) + _Tp{2} * __s - __a - __b);
      __d[1] = __dan;
      for (__n = 2; __n < _N; ++__n)
	{
	  auto __danm2 = __danm1;
	  auto __danm1 = __dan;
	  auto __dan = (_Tp(__n - 1) + __s)
		     * (__z0 * (_Tp{1} - __z0) * (_Tp(__n - 2) + __s) * __danm2
		      + ((_Tp(__n) + __s) * (_Tp{1} - _Tp{2} * __z0)
		       + (__a + __b + _Tp{1}) * __z0 - __c) * __danm1)
		     / _Tp(__n) / (_Tp(__n) + _Tp{2} * __s - __a - __b);
	  __d[__n] = __dan;
	}
      return __d;
    }

  /**
   * Do Buhring's analytic continuation.
   * @param __s  The numerator parameter, (a or b)
   * @todo This could be written to grow an existing array.
   */
  template<typename _Tp>
    _Tp
    __hyperg_buhring(_Tp __a, _Tp __b, _Tp __c, _Tp __z)
    {
      /// Find nearest z0 @f$ z_0 = e^{\plusminus i\pi/3} @f$
      using _Val = std::__detail::__num_traits_t<_Tp>;
      constexpr auto _S_pi_3 = __gnu_cxx::math::__pi_third_v<_Val>;
      const auto __z0p = __z - std::polar(_Val{1}, +_S_pi_3);
      const auto __z0m = __z - std::polar(_Val{1}, -_S_pi_3);
      const auto __z0 = std::abs(__z0m) < std::abs(__z0p) ? __z0m : __z0p;

      auto __dz = __z - __z0;
      auto __rdz = _Tp{1} / __dz;

      auto __da = __hyperg_buhring_d(0, __a, __b, __c, __z0);
      auto __suma = __da[0];
      decltype(__rdz) __terma(1);
      for (unsigned int __n = 1; __n < __da.size(); ++__n)
	{
	  __terma *= __rdz;
	  __suma += __da[__n] * __terma;
	}

      auto __db = __hyperg_buhring_d(1, __a, __b, __c, __z0);
      auto __sumb = __db[0];
      decltype(__rdz) __termb(1);
      for (unsigned int __n = 1; __n < __db.size(); ++__n)
	{
	  __termb *= __rdz;
	  __sumb += __db[__n] * __termb;
	}

      // This is where Buhring's gamma ratio might come in handy.
      const auto _Gama = std::__detail::__gamma(__a);
      const auto _Gamb = std::__detail::__gamma(__b);
      const auto _Gamc = std::__detail::__gamma(__c);
      const auto _Gambma = std::__detail::__gamma(__b - __a);
      const auto _Gamamb = std::__detail::__gamma(__a - __b);
      const auto _Gamcma = std::__detail::__gamma(__c - __a);
      const auto _Gamcmb = std::__detail::__gamma(__c - __b);

      return (_Gamc * _Gambma / _Gamb / _Gamcma) * std::pow(__dz, -__a) * __suma
	  + (_Gamc * _Gamamb / _Gama / _Gamcmb) * std::pow(__dz, -__b) * __sumb;
    }

  /**
   * @brief Return the hypergeometric function @f$ _2F_1(a,b;c;x) @f$.
   */
  template<typename _Tp>
    _Tp
    __hyperg(_Tp __a, _Tp __b, _Tp __c, _Tp __x, _Tp __rho = _Tp{0.5Q})
    {
      auto __aint = __gnu_cxx::__fp_is_integer(__a);
      auto __bint = __gnu_cxx::__fp_is_integer(__b);
      auto __cint = __gnu_cxx::__fp_is_integer(__c);
      auto __d = __c - __a - __b;
      auto __dint = __gnu_cxx::__fp_is_integer(__d);
      const auto __toler = _Tp{1000} * __gnu_cxx::__epsilon(__x);

      if (std::isnan(__a) || std::isnan(__b)
	 || std::isnan(__c) || std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__cint && __cint() <= 0)
	return __gnu_cxx::__infinity(__x);
      else if (std::abs(__c - __b) < __toler
	    || std::abs(__c - __a) < __toler)
	return std::pow(_Tp{1} - __x, __d);
      else if (std::abs(__x) <= __rho)
        return __hyperg_series(__a, __b, __c, __x);
      else if (std::abs(_Tp{1} - __x) <= __rho)
        {
	  if (__dint)
            {
	      auto __m = __dint();
	      auto _Gamm = std::__detail::__gamma(_Tp(__m));
	      auto _Gamabm = std::__detail::__gamma(__a + __b + _Tp(__m));
	      auto _Gamam = std::__detail::__gamma(__a + _Tp(__m));
	      auto _Gambm = std::__detail::__gamma(__b + _Tp(__m));
	      auto __sum = _Tp{1};
	      auto __term = _Tp{1};
	      for (int __k = 0; __k < __m; ++__k)
	        {
		  __term *= (__a + _Tp(__m + __k)) * (__b + _Tp(__m + __k))
			  / _Tp(1 - __m + __k) / _Tp(__k) * (_Tp{1} - __x);
		  __sum += __term;
		}
	      return __sum * _Gamm * _Gamabm / _Gamam / _Gambm;
	    }
	  else
            {
	      // This is where Buhring's gamma ratio might come in handy.
	      auto _Gama = std::__detail::__gamma(__a);
	      auto _Gamb = std::__detail::__gamma(__b);
	      auto _Gamc = std::__detail::__gamma(__c);
	      auto _Gamd = std::__detail::__gamma(__d);
	      auto _Gammd = std::__detail::__gamma(-__d);
	      auto _Gamcma = std::__detail::__gamma(__c - __a);
	      auto _Gamcmb = std::__detail::__gamma(__c - __b);
	      return _Gamc * _Gamd
		   * __hyperg_series(__a, __b, _Tp{1} - __d, _Tp{1} - __x)
		   / _Gamcma / _Gamcmb
		  + _Gamc * _Gammd
		   * __hyperg_series(__c - __a, __c - __b, _Tp{1} + __d, _Tp{1} - __x)
		   / _Gama / _Gamb;
	    }
	}
      else
	return __gnu_cxx::__quiet_NaN(__x);
    }

/**
 * Test harness.
 */
template<typename _Tp>
  void
  test_hyperg(_Tp proto = _Tp{})
  {
    //using _Val = _Tp;
    //using _Real = std::__detail::__num_traits_t<_Val>;
    std::vector<_Tp> parm{_Tp{0.25}, _Tp{0.5}, _Tp{1}, _Tp{2}, _Tp{5}};

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "b"
	      << ' ' << std::setw(width) << "c"
	      << ' ' << std::setw(width) << "z"
	      << ' ' << std::setw(width) << "hyperg0"
	      << ' ' << std::setw(width) << "hyperg"
	      << '\n';
    int i_min = -99;
    for (auto a : parm)
      for (auto b : parm)
	for (auto c : parm)
	  for (int i = i_min; i <= +99; ++i)
	    {
	      auto z = _Tp{0.01Q} * i;
	      try
		{
		  auto hyperg0 = __gnu_cxx::hyperg(a, b, c, z);
		  auto hyperg = __hyperg(a, b, c, z);
		  std::cout << ' ' << std::setw(width) << a
			    << ' ' << std::setw(width) << b
			    << ' ' << std::setw(width) << c
			    << ' ' << std::setw(width) << z
			    << ' ' << std::setw(width) << hyperg0
			    << ' ' << std::setw(width) << hyperg
			    << ' ' << std::setw(width) << hyperg - hyperg0
			    << '\n';
		}
	      catch(std::exception& err)
		{
		  std::cerr << "Error in test_hyperg"
			    << " at z = " << z << "; a = " << a << "; b = " << b << "; c = " << c
			    << ": " << err.what() << '\n';
		}
	    }
  }

template<typename _Tp>
  void
  test_gsl_issue(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    for (int n = 1; n <= 20; ++n)
      {
	for (int i = 1; i <= n; ++i)
	  {
	    for (int j = 1; j <= n; ++j)
	      {
		for (int k = 0; k <= 20; ++k)
		  {
		    auto x = _Tp(k * 0.05L);
		    auto gnu = __gnu_cxx::hyperg(_Tp(-i), _Tp(-n + j), _Tp(1 - i + j), x);
		    auto gsl = gsl::hyperg(_Tp(-i), _Tp(-n + j), _Tp(1 - i + j), x);
		    auto del = gnu - gsl;
		    std::cout << ' ' << std::setw(2) << n
			      << ' ' << std::setw(2) << i
			      << ' ' << std::setw(2) << j
			      << ' ' << std::setw(w) << gnu
			      << ' ' << std::setw(w) << gsl
			      << ' ' << std::setw(w) << del
			      << '\n';
		  }
	      }
	  }
      }
  }

template<typename _Tp>
  void
  test_complex(_Tp proto = _Tp{})
  {
    using namespace std::literals::complex_literals;

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    using cmplx = std::complex<_Tp>;
    constexpr auto _S_pi_3 = __gnu_cxx::math::__pi_third_v<_Tp>;
    const auto z0p = std::polar(_Tp{1}, +_S_pi_3);
    const auto z0m = std::polar(_Tp{1}, -_S_pi_3);

    cmplx a = 1, b = 2, c = 3;
    auto f0p = __hyperg_buhring(a, b, c, z0p);
    std::cout << ' ' << z0p << ' ' << f0p << '\n';
    auto f0m = __hyperg_buhring(a, b, c, z0m);
    std::cout << ' ' << z0m << ' ' << f0m << '\n';

    for (auto aa : {0.25l + 0.25il, 1.0l + 0.0il, 4.0l - 1.0il})
      {
	for (auto bb : {0.25l + 0.25il, 1.0l + 0.0il, 4.0l - 1.0il})
	  {
	    for (auto cc : {0.25l + 0.25il, 1.0l + 0.0il, 4.0l - 1.0il})
	      {
		for (int n = 1; n <= 20; ++n)
		  {
		    auto a = cmplx(aa);
		    auto b = cmplx(bb);
		    auto c = cmplx(cc);
		    auto z = cmplx(n * 0.05l + 0.0il);
		    auto gnu = __gnu_cxx::hyperg(a, b, c, z);
		    std::cout << ' ' << std::setw(w) << a
			      << ' ' << std::setw(w) << b
			      << ' ' << std::setw(w) << c
			      << ' ' << std::setw(w) << z
			      << ' ' << std::setw(w) << gnu
			      << '\n';
		  }
	      }
	  }
      }
    __hyperg_series(a, b, c, z0p);
  }

int
main()
{
  //test_gsl_issue(1.0);

  test_hyperg(1.0F);

  test_hyperg(1.0);

  test_hyperg(1.0L);

  test_complex(1.0L);
}
