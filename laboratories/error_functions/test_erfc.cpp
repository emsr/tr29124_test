/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

#include <emsr/math_constants.h>
#include <emsr/continued_fractions.h>
#include <emsr/special_functions.h>

namespace emsr
{
namespace detail
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
    erfc_series(_Tp z)
    {
      using _Real = emsr::num_traits_t<_Tp>;
      const auto s_max_iter = 1000;
      const auto s_eps = emsr::epsilon(std::real(z));
      const auto s_pi = emsr::pi_v<_Real>;

      const auto h = _Real{1} / _Real{2};

      const auto z2 = z * z;
      auto sum = _Tp{0};
      auto a = std::exp(h * h);
      auto b = _Real{1};
      const auto c = _Real{1} / a / a;
      for (auto k = 1; k < s_max_iter; ++k)
      {
	a *= c;
	b *= a;
	const auto kh = k * h;
	const auto kh2 = kh * kh;
	const auto term = b / (kh2 + z2);
	sum += term;
	if (std::abs(term) < s_eps * std::abs(sum))
	  break;
      }
      auto erfc = _Tp{2} * z * h * std::exp(-z2)
		  * (_Tp{1} / z2 / _Tp{2} + sum) / s_pi;

      auto rez = std::real(z);
      auto imz = std::imag(z);
      if (rez < _Real{0})
	erfc += _Real{2};
      if (std::abs(rez) < s_pi / h)
	{
	  auto sz = _Tp(rez < _Tp{0} ? -1 : +1) * z;
	  auto Rcorr = _Tp{2} / std::expm1(_Tp{2} * s_pi * sz / h);
	  if (rez >= _Real{0})
	    {
	      if (std::abs(imz) < s_pi / _Tp{2} - rez)
		erfc -= Rcorr;
	    }
	  else
	    {
	      if (std::abs(imz) < s_pi / _Tp{2} + rez)
		erfc += Rcorr;
	    }
	}

      return erfc;
    }

  /**
   * Use the continued fraction to evaluate the complementary
   * error function.
   */
  template<typename _Tp>
    _Tp
    erfc_cont_frac(_Tp z)
    {
      using _Real = emsr::num_traits_t<_Tp>;
      const auto s_sqrt_pi = emsr::sqrtpi_v<_Real>;
      const auto a = [](size_t k, _Tp /*z*/)
		 ->_Tp
		 { return k == 1 ? _Tp{1} : _Tp(k - 1) / _Tp{2}; };
      using _AFun = decltype(a);
      const auto b = [](size_t k, _Tp z)
		 ->_Tp
		 { return k == 0 ? _Tp{0} : (k & 1) ? z * z : _Tp{1}; };
      using _BFun = decltype(b);
      const auto w = [b](size_t k, _Tp z)
		 ->_Tp
		 { return b(k, z)
		  / (std::sqrt(_Tp{1} + _Tp(2 * k) / z / z) - _Tp{1})
		  / _Tp{2}; };
      using _WFun = decltype(w);
      using _CFrac = emsr::LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac cf(a, b, w);
      auto erfc = std::exp(-z * z) * cf(z) * z / s_sqrt_pi;
      if (std::real(z) < _Real{0})
	erfc += _Real{2};
      return erfc;
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
    fadeeva_cont_frac(_Tp z)
    {
      using _Real = emsr::num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto s_i = _Cmplx{0, 1};
      const auto s_sqrt_pi = emsr::sqrtpi_v<_Real>;
      const auto a = [](size_t k, _Tp /*z*/)
		 ->_Tp
		 { return k == 1 ? _Tp{1} : _Tp(k - 1) / _Tp{2}; };
      using _AFun = decltype(a);
      const auto b = [](size_t k, _Tp z)
		 ->_Tp
		 { return k == 0 ? _Tp{0} : (k & 1) ? z * z : _Tp{1}; };
      using _BFun = decltype(b);
      const auto w = [b](size_t k, _Tp z)
		 ->_Tp
		 { return b(k, z)
		  / (std::sqrt(_Tp{1} + _Tp(2 * k) / z / z) - _Tp{1})
		  / _Tp{2}; };
      using _WFun = decltype(w);
      using _CFrac = emsr::LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac cf(a, b, w);
      const auto cfrac = s_i * cf(z) / s_sqrt_pi;
      if (std::real(z) < _Real{0})
	erfc += _Real{2};
      return cfrac;
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
    erfc_cont_frac_even(_Tp z)
    {
      const auto s_sqrt_pi = emsr::sqrtpi_v<_Tp>;
      const auto a = [](size_t k, _Tp z)
		 ->_Tp
		 { return k == 1
			 ? _Tp{2} * z
			 : (k & 1) ? z * z : _Tp{1}; };
      using _AFun = decltype(a);
      const auto b = [](size_t k, _Tp z)
		 ->_Tp
		 { return k == 0
			? _Tp{0}
			: (k & 1) ? _Tp{2} * z * z
				  : _Tp{1}; };
      using _BFun = decltype(b);
      const auto w = [b](size_t, _Tp)->_Tp{ return _Tp{0}; };
      using _WFun = decltype(w);
      using _CFrac = emsr::LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac cf(a, b, w);
      return std::exp(-z * z) * cf(z) / s_sqrt_pi;
    }

  /**
   * Return the Fadeeva function.
   * @f[
   *    w(z) = e^{-z^2}erfc(-iz)
   * @f]
   */
  template<typename _Tp>
    _Tp
    fadeeva(_Tp z)
    {
      using _Real = emsr::num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto s_NaN = emsr::quiet_NaN(std::real(z));
      const auto s_i = _Cmplx{0, 1};
      if (std::isnan(z))
	return _Cmplx{s_NaN, s_NaN};
      else if (std::real(z) < _Real{0})
	return _Real{2} * std::exp(-z * z) - fadeeva(-z);
      else if (std::abs(z) < _Real{15})
	return erfc_series(z);
      else
	return fadeeva_cont_frac(z);
    }

  /**
   * Return the complementary error function.
   */
  template<typename _Tp>
    _Tp
    erfc(_Tp x)
    {
      const auto s_inf = emsr::infinity(x);
      const auto s_cfrac = _Tp{0.025} * std::numeric_limits<_Tp>::digits;

      if (std::isnan(x))
	return x;
      else if (x == -s_inf)
	return _Tp{2};
      else if (x == +s_inf)
	return _Tp{0};
      else if (x < s_cfrac)
	return erfc_series(x);
      else
	return erfc_cont_frac(x);
    }

  /**
   * Return the error function.
   */
  template<typename _Tp>
    _Tp
    erf(_Tp x)
    {
      const auto s_inf = std::numeric_limits<_Tp>::infinity();
      const auto s_cfrac = _Tp{0.025} * std::numeric_limits<_Tp>::digits;

      if (std::isnan(x))
	return x;
      else if (x == -s_inf)
	return _Tp{0};
      else if (x == +s_inf)
	return _Tp{1};
      else if (x < s_cfrac)
	return _Tp{1} - erfc_series(x);
      else
	return _Tp{1} - erfc_cont_frac(x);
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
    erfc_asymp(int k, _Tp x)
    {
      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();
      const auto s_sqrt_pi = emsr::sqrtpi_v<_Tp>;
      const auto s_max_iter = 200;
      const auto __2x = _Tp{2} * x;
      const auto __2xm2 = -_Tp{1} / (__2x * __2x);
      auto term = _Tp{1};
      auto sum = term;
      auto kfact = emsr::detail::factorial<_Tp>(k);
      auto prev_term = std::abs(term);
      for (int m = 1; m < s_max_iter; ++m)
	{
	  term *= __2xm2 * _Tp(k + 2 * m) * _Tp(k + 2 * m - 1)
		  / _Tp(m) / kfact;
	  if (std::abs(term) > prev_term)
	    break;
	  prev_term = std::abs(term);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      const auto fact = _Tp{2} * std::exp(x * x)
			/ std::pow(__2x, _Tp(k + 1)) / s_sqrt_pi;
      return fact * sum;
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
    erfc_recur(int k, _Tp x)
    {
      const auto s_sqrt_pi = emsr::sqrtpi_v<_Tp>;

      auto erfcm2 = _Tp{2} * std::exp(-x * x) / s_sqrt_pi;
      if (k == -1)
	return erfcm2;

      auto erfcm1 = erfc(x);
      if (k == 0)
	return erfcm1;

      auto erfcm0 = -x * erfcm1 + erfcm2 / _Tp{2};
      for (int i = 2; i <= k; ++i)
	{
	  erfcm2 = erfcm1;
	  erfcm1 = erfcm0;
	  erfcm0 = (-x * erfcm1 + erfcm2 / _Tp{2}) / i;
	}
      return erfcm0;
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    erfc(int k, _Tp x)
    {
      const auto s_inf = std::numeric_limits<_Tp>::infinity();

      if (std::isnan(x))
	return x;
      else if (x == -s_inf)
	return +s_inf;
      else if (x == +s_inf)
	return _Tp{0};
      else
	return erfc_recur(k, x);
    }

} // namespace detail
} // namespace emsr

template<typename _Tp>
  void
  test_erfc(_Tp proto = _Tp{})
  {
    using namespace std::complex_literals;
    const auto s_NaN = std::numeric_limits<_Tp>::quiet_NaN();
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto del = _Tp{1} / _Tp{100};
    for (int i = -200; i <= 1000; ++i)
      {
	auto x = del * i;
	auto erfc = std::erfc(x);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << erfc;

	try
	  {
	    auto erfcs = emsr::detail::erfc_series(x);
	    std::cout << ' ' << std::setw(width) << erfcs;
	  }
	catch (std::runtime_error& err)
	  {
	    std::cout << ' ' << std::setw(width) << s_NaN;
	    std::cerr << err.what() << '\n';
	  }

	try
	  {
	    auto erfccf = emsr::detail::erfc_cont_frac(x);
	    std::cout << ' ' << std::setw(width) << erfccf;
	  }
	catch (std::runtime_error& err)
	  {
	    std::cout << ' ' << std::setw(width) << s_NaN;
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
    for (int k = -200; k <= 200; ++k)
      {
	auto x = k * _Tp{0.01Q};
	auto erfx = emsr::detail::erf(x);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << erfx
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "erfc(x)"
	      << '\n';
    for (int k = -200; k <= 200; ++k)
      {
	auto x = k * _Tp{0.01Q};
	auto erfcx = emsr::detail::erfc(x);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << erfcx
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
    for (int k = -200; k <= 500; ++k)
      {
	auto x = k * _Tp{0.01Q};
	std::cout << ' ' << std::setw(w) << x;
	for (int n = -1; n <= 5; ++n)
	  std::cout << ' ' << std::setw(w) << emsr::detail::erfc(n, x);
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

  //std::cout << "\n\n  __float128\n";
  //std::cout << "  ==========\n";
  //test_erfc<__float128>();
}

