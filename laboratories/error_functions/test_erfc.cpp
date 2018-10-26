/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_erfc test_erfc.cpp -lquadmath
./test_erfc > test_erfc.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_erfc test_erfc.cpp -lquadmath
./test_erfc > test_erfc.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <LentzContinuedFraction.tcc>

namespace std
{
namespace __detail
{

  /**
   * Return the complementary error function by series:
   * @f[
   *    erfc(z) = 
   * @f]
   * @see On the Evaluation of Integrals Related to the Error Function
   * of Real and Complex Variable with High Relative Precision, Masatake Mori,
   * Publ. RIMS, Kyoto Univ. 19 (1983), 1081-1094
   */
  template<typename _Tp>
    _Tp
    __erfc_series(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      const auto _S_max_iter = 1000;
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__z));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));

      const auto __h = _Real{1} / _Real{2};

      const auto __z2 = __z * __z;
      auto __sum = _Tp{0};
      auto __a = std::exp(__h * __h);
      auto __b = _Real{1};
      const auto __c = _Real{1} / __a / __a;
      for (auto __k = 1; __k < _S_max_iter; ++__k)
      {
	__a *= __c;
	__b *= __a;
	const auto __kh = __k * __h;
	const auto __kh2 = __kh * __kh;
	const auto __term = __b / (__kh2 + __z2);
	__sum += __term;
	if (std::abs(__term) < _S_eps * std::abs(__sum))
	  break;
      }
      auto __erfc = _Tp{2} * __z * __h * std::exp(-__z2)
		  * (_Tp{1} / __z2 / _Tp{2} + __sum) / _S_pi;

      auto __rez = std::real(__z);
      auto __imz = std::imag(__z);
      if (__rez < _Real{0})
	__erfc += _Real{2};
      if (std::abs(__rez) < _S_pi / __h)
	{
	  auto __sz = _Tp(__rez < _Tp{0} ? -1 : +1) * __z;
	  auto __Rcorr = _Tp{2} / std::expm1(_Tp{2} * _S_pi * __sz / __h);
	  if (__rez >= _Real{0})
	    {
	      if (std::abs(__imz) < _S_pi / _Tp{2} - __rez)
		__erfc -= __Rcorr;
	    }
	  else
	    {
	      if (std::abs(__imz) < _S_pi / _Tp{2} + __rez)
		__erfc += __Rcorr;
	    }
	}

      return __erfc;
    }

  /**
   * Use the continued fraction to evaluate the complementary
   * error function.
   */
  template<typename _Tp>
    _Tp
    __erfc_cont_frac(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__z));
      const auto __a = [](size_t __k, _Tp /*__z*/)
		 ->_Tp
		 { return __k == 1 ? _Tp{1} : _Tp(__k - 1) / _Tp{2}; };
      using _AFun = decltype(__a);
      const auto __b = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 0 ? _Tp{0} : __k & 1 ? __z * __z : _Tp{1}; };
      using _BFun = decltype(__b);
      const auto __w = [__b](size_t __k, _Tp __z)
		 ->_Tp
		 { return __b(__k, __z)
		  / (std::sqrt(_Tp{1} + _Tp(2 * __k) / __z / __z) - _Tp{1})
		  / _Tp{2}; };
      using _WFun = decltype(__w);
      using _CFrac = _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac __cf(__a, __b, __w);
      auto __erfc = std::exp(-__z * __z) * __cf(__z) * __z / _S_sqrt_pi;
      if (std::real(__z) < _Real{0})
	__erfc += _Real{2};
      return __erfc;
    }

  /**
   * Use the continued fraction to evaluate the complementary
   * error function.
   * @f[
   *    w(z) = 
   * @f]
   */
  template<typename _Tp>
    _Tp
    __fadeeva_cont_frac(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_i = _Cmplx{0, 1};
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__z));
      const auto __a = [](size_t __k, _Tp /*__z*/)
		 ->_Tp
		 { return __k == 1 ? _Tp{1} : _Tp(__k - 1) / _Tp{2}; };
      using _AFun = decltype(__a);
      const auto __b = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 0 ? _Tp{0} : __k & 1 ? __z * __z : _Tp{1}; };
      using _BFun = decltype(__b);
      const auto __w = [__b](size_t __k, _Tp __z)
		 ->_Tp
		 { return __b(__k, __z)
		  / (std::sqrt(_Tp{1} + _Tp(2 * __k) / __z / __z) - _Tp{1})
		  / _Tp{2}; };
      using _WFun = decltype(__w);
      using _CFrac = _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac __cf(__a, __b, __w);
      const auto __cfrac = _S_i * __cf(__z) / _S_sqrt_pi;
      if (std::real(__z) < _Real{0})
	__erfc += _Real{2};
      return __cfrac;
    }

  /**
   * Use the even continued fraction to evaluate the complementary
   * error function.
   * @f[
   *    erfc(z) = 
   * @f]
   */
  template<typename _Tp>
    _Tp
    __erfc_cont_frac_even(_Tp __z)
    {
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__z));
      const auto __a = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 1
			 ? _Tp{2} * __z
			 : __k & 1 ? __z * __z : _Tp{1}; };
      using _AFun = decltype(__a);
      const auto __b = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 0
			? _Tp{0}
			: __k & 1 ? _Tp{2} * __z * __z
				  : _Tp{1}; };
      using _BFun = decltype(__b);
      const auto __w = [__b](size_t, _Tp)->_Tp{ return _Tp{0}; };
      using _WFun = decltype(__w);
      using _CFrac = _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac __cf(__a, __b, __w);
      return std::exp(-__z * __z) * __cf(__z) / _S_sqrt_pi;
    }

  /**
   * Return the Fadeeva function.
   * @f[
   *    w(z) = e^{-z^2}erfc(-iz)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __fadeeva(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::real(__z));
      const auto _S_i = _Cmplx{0, 1};
      if (std::isnan(__z))
	return _Cmplx{_S_NaN, _S_NaN};
      else if (std::real(__z) < _Real{0})
	return _Real{2} * std::exp(-__z * __z) - __fadeeva(-__z);
      else if (std::abs(__z) < _Real{15})
	return __erfc_series(__z);
      else
	return __fadeeva_cont_frac(__z);
    }

  /**
   * Return the complementary error function.
   */
  template<typename _Tp>
    _Tp
    __erfc(_Tp __x)
    {
      const auto _S_inf = __gnu_cxx::__infinity(__x);
      const auto _S_cfrac = _Tp{0.025} * std::numeric_limits<_Tp>::digits;

      if (std::isnan(__x))
	return __x;
      else if (__x == -_S_inf)
	return _Tp{2};
      else if (__x == +_S_inf)
	return _Tp{0};
      else if (__x < _S_cfrac)
	return __erfc_series(__x);
      else
	return __erfc_cont_frac(__x);
    }

  /**
   * Return the error function.
   */
  template<typename _Tp>
    _Tp
    __erf(_Tp __x)
    {
      const auto _S_inf = std::numeric_limits<_Tp>::infinity();
      const auto _S_cfrac = _Tp{0.025} * std::numeric_limits<_Tp>::digits;

      if (std::isnan(__x))
	return __x;
      else if (__x == -_S_inf)
	return _Tp{0};
      else if (__x == +_S_inf)
	return _Tp{1};
      else if (__x < _S_cfrac)
	return _Tp{1} - __erfc_series(__x);
      else
	return _Tp{1} - __erfc_cont_frac(__x);
    }

  /**
   * Return scaled repeated integrals of the erfc function by asymptotic series.
   * The experfc function is defined by
   * @f[
   *   erfc(k, x) = i^kerfc(x)
   * @f]
   * where the integral of the comlementary error function is:
   * @f[
   *   i^kerfc(x) = \frac{2}{k!\sqrt{\pi}}
   *       \int_{x}^{\infty}(t-x)^ke^{-t^2}dt
   * @f]
   * @see Cuyt, et.al. 13.3.2
   */
  template<typename _Tp>
    _Tp
    __erfc_asymp(int __k, _Tp __x)
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
      const auto __fact = _Tp{2} * std::exp(__x * __x)
			/ std::pow(__2x, _Tp(__k + 1)) / _S_sqrt_pi;
      return __fact * __sum;
    }

  /**
   * Return scaled repeated integrals of the erfc function by asymptotic series.
   * The experfc function is defined by
   * @f[
   *   erfc(k, x) = i^kerfc(x)
   * @f]
   * where the integral of the comlementary error function is:
   * @f[
   *   i^kerfc(x) = \frac{2}{k!\sqrt{\pi}}
   *       \int_{x}^{\infty}(t-x)^ke^{-t^2}dt
   * @f]
   *
   * The recursion is:
   * @f[
   *    i^kerfc(x) = -\frac{x}{k}i^{k-1}erfc(x) + \frac{1}{2k}i^{k-2}erfc(x)
   * @f]
   * starting with
   * @f[
   *    i^0erfc(x) = erfc(x) \hbox{  and  }
   *        i^{-1}erfc(x) = \frac{2}{\sqrt{\pi}}e^{-x^2}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __erfc_recur(int __k, _Tp __x)
    {
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(__x);

      auto __erfcm2 = _Tp{2} * std::exp(-__x * __x) / _S_sqrt_pi;
      if (__k == -1)
	return __erfcm2;

      auto __erfcm1 = __erfc(__x);
      if (__k == 0)
	return __erfcm1;

      auto __erfcm0 = -__x * __erfcm1 + __erfcm2 / _Tp{2};
      for (int __i = 2; __i <= __k; ++__i)
	{
	  __erfcm2 = __erfcm1;
	  __erfcm1 = __erfcm0;
	  __erfcm0 = (-__x * __erfcm1 + __erfcm2 / _Tp{2}) / __i;
	}
      return __erfcm0;
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    __erfc(int __k, _Tp __x)
    {
      const auto _S_inf = std::numeric_limits<_Tp>::infinity();

      if (std::isnan(__x))
	return __x;
      else if (__x == -_S_inf)
	return +_S_inf;
      else if (__x == +_S_inf)
	return _Tp{0};
      else
	return __erfc_recur(__k, __x);
    }

} // namespace std
} // namespace __detail

template<typename _Tp>
  void
  test_erfc(_Tp proto = _Tp{})
  {
    using namespace std::literals::complex_literals;
    const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto del = _Tp{1} / _Tp{100};
    for (int i = -200; i <= 1000; ++i)
      {
	auto x = del * i;
	auto __erfc = std::erfc(x);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << __erfc;

	try
	  {
	    auto __erfcs = std::__detail::__erfc_series(x);
	    std::cout << ' ' << std::setw(width) << __erfcs;
	  }
	catch (std::runtime_error& err)
	  {
	    std::cout << ' ' << std::setw(width) << _S_NaN;
	    std::cerr << err.what() << '\n';
	  }

	try
	  {
	    auto __erfccf = std::__detail::__erfc_cont_frac(x);
	    std::cout << ' ' << std::setw(width) << __erfccf;
	  }
	catch (std::runtime_error& err)
	  {
	    std::cout << ' ' << std::setw(width) << _S_NaN;
	    std::cerr << err.what() << '\n';
	  }

	std::cout << '\n';
      }
  }

/**
 * Test the scaled complementary error function - experfc(x) = exp(x^2)erfc(x).
 */
template<typename _Tp>
  void
  plot_erfc()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 8 + std::cout.precision();

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "erf(x)"
	      << '\n';
    for (int __k = -200; __k <= 200; ++__k)
      {
	auto __x = __k * _Tp{0.01Q};
	auto __erfx = std::__detail::__erf(__x);
	std::cout << ' ' << std::setw(w) << __x
		  << ' ' << std::setw(w) << __erfx
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "erfc(x)"
	      << '\n';
    for (int __k = -200; __k <= 200; ++__k)
      {
	auto __x = __k * _Tp{0.01Q};
	auto __erfcx = std::__detail::__erfc(__x);
	std::cout << ' ' << std::setw(w) << __x
		  << ' ' << std::setw(w) << __erfcx
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "i^{-1}erfc(x)"
	      << ' ' << std::setw(w) << "i^0erfc(x)"
	      << ' ' << std::setw(w) << "i^1erfc(x)"
	      << ' ' << std::setw(w) << "i^2erfc(x)"
	      << ' ' << std::setw(w) << "i^3erfc(x)"
	      << ' ' << std::setw(w) << "i^4erfc(x)"
	      << ' ' << std::setw(w) << "i^5erfc(x)"
	      << '\n';
    for (int __k = -200; __k <= 500; ++__k)
      {
	auto __x = __k * _Tp{0.01Q};
	std::cout << ' ' << std::setw(w) << __x;
	for (int __n = -1; __n <= 5; ++__n)
	  std::cout << ' ' << std::setw(w) << std::__detail::__erfc(__n, __x);
	std::cout << '\n';
      }
  }


int
main()
{
  plot_erfc<double>();

  std::cout << "\n\n  float\n";
  std::cout << "  =====\n";
  test_erfc(1.0F);

  std::cout << "\n\n  double\n";
  std::cout << "  ======\n";
  test_erfc(1.0);

  std::cout << "\n\n  long double\n";
  std::cout << "  ===========\n";
  test_erfc(1.0L);

  std::cout << "\n\n  __float128\n";
  std::cout << "  ==========\n";
  test_erfc<__float128>();
}

