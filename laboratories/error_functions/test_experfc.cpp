/**
 *
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

#include <emsr/float128_io.h>
#include <mpreal.h>

  template<typename Tp>
    Tp
    experfc_func(Tp x)
    {
      mpfr::mpreal X(static_cast<long double>(x), 1024);
      if (x < Tp{20})
	return static_cast<long double>(mpfr::exp(X * X) * mpfr::erfc(X));
      else
	{
return 0;
/*
	  constexpr auto log2_e = 1.442695040888963407359924681001892137427L;
	  auto lg = X * X + mpfr::log(mpfr::erfc(X));
	  if (lg * log2_e < 1022)
	    return (long double)(mpfr::exp(lg));
	  else
	    return Tp{0};
*/
	}
    }

  /**
   * The experf function is defined by
   * @f[
   *   experf(x) = exp(x^2)erf(x)
   *             = \frac{2}{\sqrt{\pi}}
   *               \sum_{k=0}^{\infty}\frac{x^{2k+1}}{(2k+1)!!}
   * @f]
   */
  template<typename Tp>
    Tp
    experf_series(Tp x)
    {
      constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
      const auto s_sqrt_pi = Tp{1.772453850905516027298167483341145182797Q};
      const auto s_max_iter = 200;
      const auto x2 = Tp{2} * x * x;
      auto term = x;
      auto sum = term;
      for (int k = 1; k < s_max_iter; ++k)
	{
	  term *= x2 / Tp(2 * k + 1);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return  Tp{2} * sum / s_sqrt_pi;
    }

  /**
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x^2)erfc(x) = exp(x^2)(1 - erf(x))
   * @f]
   */
  template<typename Tp>
    Tp
    experfc_func_bad(Tp x)
    { return std::exp(x * x) - experf_series(x); }

  /**
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x^2)erfc(x) = exp(x^2)(1 - erf(x))
   * @f]
   * The series is
   * @f[
   *    experfc(x) = \sum_{k=0}^{\infty} \frac{(-x)^k}{\Gamma(1+k/2)}
   * @f]
   */
  template<typename Tp>
    Tp
    experfc_series(Tp x)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      const auto s_sqrt_pi = Tp{1.772453850905516027298167483341145182797Q};
      const auto s_max_iter = 200;
      auto gamr_e = Tp{1};
      auto gamr_o = Tp{2} / s_sqrt_pi;
      auto term = Tp{0};
      auto sum = Tp{0};
      auto xk = Tp{1};

      for (int k = 0; k < s_max_iter; ++k)
	{
	  if (k & 1)
	    {
	      term = gamr_o * xk;
	      sum += term;
	      gamr_o /= Tp(2 + k) / Tp{2};
	    }
	  else
	    {
	      term = gamr_e * xk;
	      sum += term;
	      gamr_e /= Tp(1 + k / 2);
	    }
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	  xk *= -x;
	}
      return sum;
    }

  /**
   * Return a quick approximation to the experfc function.
   * The experfc function is defined by
   * @f[
   *   experfc(x) = exp(x)erfc(\sqrt{x})
   * @f]
   * Note the change @f$ x -> \sqrt{x} @f$ relative to the other functions here.
   * The approximation is:
   * @f[
   *   experfc(x) = \frac{2}{\sqrt{\pi x}
   *       + \sqrt{\pi(x+2)-2(\pi-2)exp(-\sqrt(5x/7))}}
   * @f]
   * Surprisingly, this is accurate to within 0.1% over the whole range
   * [0, infty].  It is used to start AGM algorithms of the experfc function.
   */
  template<typename Tp>
    Tp
    experfc_approx(Tp x)
    {
      const auto s_pi =  Tp{3.1415926535897932384626433832795029Q};
      const auto s_sqrt_pi = Tp{1.772453850905516027298167483341145182797Q};
      auto experfc = s_sqrt_pi * std::sqrt(x)
		     + std::sqrt(s_pi * (x + Tp{2})
			       - Tp{2} * (s_pi - Tp{2})
			       * std::exp(-std::sqrt(Tp{5} * x / Tp{7})));
      return Tp{2} / experfc;
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
  template<typename Tp>
    Tp
    experfc_asymp(Tp x)
    {
      constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
      const auto s_sqrt_pi = Tp{1.772453850905516027298167483341145182797Q};
      const auto s_max_iter = 200;
      const auto xm2 = -Tp{1} / (Tp{2} * x * x);
      auto term = Tp{1} / x;
      auto sum = term;
      auto prev_term = std::abs(term);
      for (int k = 1; k < s_max_iter; ++k)
	{
	  term *= xm2 * Tp(2 * k - 1);
	  sum += term;
	  if (std::abs(term) > prev_term)
	    break;
	  prev_term = std::abs(term);
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return sum / s_sqrt_pi;
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
  template<typename Tp>
    Tp
    experfc_cont_frac(Tp x)
    {
      const auto s_sqrt_2 = Tp{1.414213562373095048801688724209698078569Q};
      const auto s_sqrt_pi = Tp{1.772453850905516027298167483341145182797Q};
      const auto b = std::sqrt(Tp{2} * x * x);
      auto s = b;
      for (int n = 100; n >= 1; --n)
	{
	  auto a = Tp(n);
	  s = b + a / s;
	}
      s = (s_sqrt_2 / s_sqrt_pi) / s;
      return s;
    }

  /**
   * exp(x^2) erfc(x).
   */
  template<typename Tp>
    std::complex<Tp>
    experfc_series_aw(std::complex<Tp> z,
			int trunc_lo = 32, int trunc_hi = 193,
			Tp sep = Tp{8})
    {
      const auto s_pi =  Tp{3.1415926535897932384626433832795029Q};
      const auto s_sqrt_pi = Tp{1.772453850905516027298167483341145182797Q};
      const auto s_i = std::complex<Tp>{0, 1};

      if (sep < Tp{0} || std::abs(z) <= sep)
        {
          // Infinite series approximation, Abramowitz p. 313 7.1.29
          auto x = std::real(z);
          auto y = std::imag(z);

          auto ez2 = std::exp(z * z);

          auto s1 = erf(x);

          auto k1 = std::exp(-y * y);
          auto k2 = std::exp(Tp{2} * s_i * x * y);

          std::complex<Tp> s2;
          if (x == Tp{0})
            s2 = s_i * y * k1 / s_pi;
          else
            s2 = k1 * (k2 - Tp{1} / (Tp{2} * s_pi * x));

          auto retval = ez2 * s1 + s2;

          if (y != Tp{0})
	    {
              auto sum = std::complex<Tp>{};
              for (unsigned n = 1; n <= trunc_lo; ++n)
		{
		  auto enn = Tp(n);
		  auto s3 = std::exp(-enn * enn / Tp{4})
			    / (enn * enn + Tp{4} * x * x);
                  auto s4 = Tp{2} * x * k1 * k2
			    - (x + s_i * enn / Tp{2})
			     * std::exp(-y * (enn + y))
			    - (x - s_i * enn / Tp{2})
			     * std::exp(+y * (enn - y));
                  sum += s3 * s4;
		}
              retval += Tp{2} * sum / s_pi;
	    }
          return ez2 - retval;
        }
      else
	{
          // Asymptotic expansion, Abramowitz p. 312 7.1.23
          bool isneg = (std::real(z) < Tp{0});
          if (isneg)
	    z = -z;

          std::complex<Tp> s = Tp{1};
          std::complex<Tp> y = Tp{2} * z * z;
          for (auto n = trunc_hi; n >= 1; n -= 2)
	    s = Tp{1} - Tp(n) * (s / y);

          auto retval = s / (s_sqrt_pi * z);

          if (isneg)
	    {
              z = -z;
              retval = Tp{2} - retval;
            }

          return retval;
        }
    }

  /**
   *
   */
  template<typename Tp>
    Tp
    erfc_scaled(Tp x)
    {
      const auto s_inf = std::numeric_limits<Tp>::infinity();
      const auto s_cfrac = Tp{0.025} * std::numeric_limits<Tp>::digits;
      // The asymptotic series gets good by here but never really beats C.F.
      //const auto s_asymp = Tp{0.18} * std::numeric_limits<Tp>::digits;

      if (std::isnan(x))
	return x;
      else if (x == -s_inf)
	return +s_inf;
      else if (x == +s_inf)
	return Tp{0};
      else if (x < s_cfrac)
	return experfc_series(x);
      else
	return experfc_cont_frac(x);
    }

/**
 * Test the scaled complementary error function - experfc(x) = exp(x^2)erfc(x).
 */
template<typename Tp>
  void
  plot_experfc()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "experfc(x)"
	      << '\n';
    for (int k = -200; k <= 5500; ++k)
      {
	auto x = k * Tp{0.01Q};
	auto experfc = erfc_scaled(x);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << experfc
		  << '\n';
      }
  }

/**
 * Test the scaled complementary error function - experfc(x) = exp(x^2)erfc(x).
 */
template<typename Tp>
  void
  test_experfc()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
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
    for (int i = 0; i <= 5500; ++i)
      {
	auto x = i * Tp{0.01Q};
	auto experfc_apr = experfc_approx(x * x);
	auto experfc_fun = experfc_func(x);
	auto experfc_cfr = experfc_cont_frac(x);
	auto experfc_asy = experfc_asymp(x);
	auto experfc_ser = experfc_series(x);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << experfc_apr
		  << ' ' << std::setw(w) << experfc_fun
		  << ' ' << std::setw(w) << experfc_cfr
		  << ' ' << std::setw(w) << experfc_asy
		  << ' ' << std::setw(w) << experfc_ser
		  << ' ' << std::setw(w) << (experfc_cfr - experfc_fun) / experfc_fun
		  << ' ' << std::setw(w) << (experfc_asy - experfc_fun) / experfc_fun
		  << ' ' << std::setw(w) << (experfc_ser - experfc_fun) / experfc_fun
		  << '\n';
      }
  }

int
main()
{
  plot_experfc<double>();

  std::cout << "\n\n  float\n";
  std::cout << "  =====\n";
  test_experfc<float>();

  std::cout << "\n\n  double\n";
  std::cout << "  ======\n";
  test_experfc<double>();

  std::cout << "\n\n  long double\n";
  std::cout << "  ===========\n";
  test_experfc<long double>();

  //std::cout << "\n\n  __float128\n";
  //std::cout << "  ==========\n";
  //test_experfc<__float128>();
}
