/**
 *
 */

#include <emsr/math_constants.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/float128_math.h>
#include <emsr/float128_limits.h>

#include <emsr/notsospecfun.h> // For complex log1p.
#include <emsr/numeric_limits.h>
#include <emsr/sf_hyperg.h>

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
   *   @param  a  The first @a numerator parameter.
   *   @param  b  The second @a numerator parameter.
   *   @param  c  The @a denominator parameter.
   *   @param  x  The argument of the confluent hypergeometric function.
   *   @return  The confluent hypergeometric function.
   */
  template<typename _Tp>
    _Tp
    hyperg_series(_Tp a, _Tp b, _Tp c, _Tp x)
    {
      using _Val = emsr::num_traits_t<_Tp>;
      const auto eps = emsr::epsilon<_Val>();
      const unsigned int s_max_iter = 100000u;
      auto aint = emsr::fp_is_integer(a);
      auto bint = emsr::fp_is_integer(b);
      auto max_iter = s_max_iter;
      if (aint && aint() < 0 && int(s_max_iter) > -aint())
	max_iter = -aint();
      if (bint && bint() < 0 && int(s_max_iter) > -bint())
	max_iter = -bint();

      auto term = _Tp{1};
      auto Fabc = _Tp{1};
      unsigned int i;
      for (i = 0; i < max_iter; ++i)
	{
	  term *= (a + _Tp(i)) * (b + _Tp(i)) * x
		  / ((c + _Tp(i)) * _Tp(1 + i));
	  if (std::abs(term) < eps)
	    break;
	  Fabc += term;
	}
      if (i == max_iter)
	throw std::runtime_error("Series failed to converge in hyperg_series.");

      return Fabc;
    }

  /**
   * Do Buhring's analytic continuation.
   * @param s  The numerator parameter, (a or b)
   * @todo This could be written to grow an existing array.
   */
  template<typename _Tp>
    std::vector<_Tp>
    hyperg_buhring_d(int ss, _Tp a, _Tp b, _Tp c, _Tp z0)
    {
      const unsigned int _N = 20; // ?
      std::vector<_Tp> d(_N);
      const auto s = ss == 0 ? a : b;
      unsigned int n = 1;
      auto danm2 = _Tp{0};
      auto danm1 = _Tp{1};
      d[0] = danm1;
      auto dan = (_Tp(n - 1) + s)
		 * (z0 * (_Tp{1} - z0) * (_Tp(n - 2) + s) * danm2
		  + ((_Tp(n) + s) * (_Tp{1} - _Tp{2} * z0)
		   + (a + b + _Tp{1}) * z0 - c) * danm1)
		 / _Tp(n) / (_Tp(n) + _Tp{2} * s - a - b);
      d[1] = dan;
      for (n = 2; n < _N; ++n)
	{
	  auto danm2 = danm1;
	  auto danm1 = dan;
	  auto dan = (_Tp(n - 1) + s)
		     * (z0 * (_Tp{1} - z0) * (_Tp(n - 2) + s) * danm2
		      + ((_Tp(n) + s) * (_Tp{1} - _Tp{2} * z0)
		       + (a + b + _Tp{1}) * z0 - c) * danm1)
		     / _Tp(n) / (_Tp(n) + _Tp{2} * s - a - b);
	  d[n] = dan;
	}
      return d;
    }

  /**
   * Do Buhring's analytic continuation.
   * @param s  The numerator parameter, (a or b)
   * @todo This could be written to grow an existing array.
   */
  template<typename _Tp>
    _Tp
    hyperg_buhring(_Tp a, _Tp b, _Tp c, _Tp z)
    {
      /// Find nearest z0 @f$ z_0 = e^{\plusminus i\pi/3} @f$
      using _Val = emsr::num_traits_t<_Tp>;
      constexpr auto s_pi_3 = emsr::pi_v<_Val> / _Val{3};
      const auto z0p = z - std::polar(_Val{1}, +s_pi_3);
      const auto z0m = z - std::polar(_Val{1}, -s_pi_3);
      const auto z0 = std::abs(z0m) < std::abs(z0p) ? z0m : z0p;

      auto dz = z - z0;
      auto rdz = _Tp{1} / dz;

      auto da = hyperg_buhring_d(0, a, b, c, z0);
      auto suma = da[0];
      decltype(rdz) terma(1);
      for (unsigned int n = 1; n < da.size(); ++n)
	{
	  terma *= rdz;
	  suma += da[n] * terma;
	}

      auto db = hyperg_buhring_d(1, a, b, c, z0);
      auto sumb = db[0];
      decltype(rdz) termb(1);
      for (unsigned int n = 1; n < db.size(); ++n)
	{
	  termb *= rdz;
	  sumb += db[n] * termb;
	}

      // This is where Buhring's gamma ratio might come in handy.
      const auto _Gama = emsr::detail::gamma(a);
      const auto _Gamb = emsr::detail::gamma(b);
      const auto _Gamc = emsr::detail::gamma(c);
      const auto _Gambma = emsr::detail::gamma(b - a);
      const auto _Gamamb = emsr::detail::gamma(a - b);
      const auto _Gamcma = emsr::detail::gamma(c - a);
      const auto _Gamcmb = emsr::detail::gamma(c - b);

      return (_Gamc * _Gambma / _Gamb / _Gamcma) * std::pow(dz, -a) * suma
	  + (_Gamc * _Gamamb / _Gama / _Gamcmb) * std::pow(dz, -b) * sumb;
    }

  /**
   * @brief Return the hypergeometric function @f$ _2F_1(a,b;c;x) @f$.
   */
  template<typename _Tp>
    _Tp
    hyperg(_Tp a, _Tp b, _Tp c, _Tp x, _Tp rho = _Tp{0.5Q})
    {
      auto aint = emsr::fp_is_integer(a);
      auto bint = emsr::fp_is_integer(b);
      auto cint = emsr::fp_is_integer(c);
      auto d = c - a - b;
      auto dint = emsr::fp_is_integer(d);
      const auto toler = _Tp{1000} * emsr::epsilon(x);

      if (std::isnan(a) || std::isnan(b)
	 || std::isnan(c) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (cint && cint() <= 0)
	return emsr::infinity(x);
      else if (std::abs(c - b) < toler
	    || std::abs(c - a) < toler)
	return std::pow(_Tp{1} - x, d);
      else if (std::abs(x) <= rho)
        return hyperg_series(a, b, c, x);
      else if (std::abs(_Tp{1} - x) <= rho)
        {
	  if (dint)
            {
	      auto m = dint();
	      auto _Gamm = emsr::detail::gamma(_Tp(m));
	      auto _Gamabm = emsr::detail::gamma(a + b + _Tp(m));
	      auto _Gamam = emsr::detail::gamma(a + _Tp(m));
	      auto _Gambm = emsr::detail::gamma(b + _Tp(m));
	      auto sum = _Tp{1};
	      auto term = _Tp{1};
	      for (int k = 0; k < m; ++k)
	        {
		  term *= (a + _Tp(m + k)) * (b + _Tp(m + k))
			  / _Tp(1 - m + k) / _Tp(k) * (_Tp{1} - x);
		  sum += term;
		}
	      return sum * _Gamm * _Gamabm / _Gamam / _Gambm;
	    }
	  else
            {
	      // This is where Buhring's gamma ratio might come in handy.
	      auto _Gama = emsr::detail::gamma(a);
	      auto _Gamb = emsr::detail::gamma(b);
	      auto _Gamc = emsr::detail::gamma(c);
	      auto _Gamd = emsr::detail::gamma(d);
	      auto _Gammd = emsr::detail::gamma(-d);
	      auto _Gamcma = emsr::detail::gamma(c - a);
	      auto _Gamcmb = emsr::detail::gamma(c - b);
	      return _Gamc * _Gamd
		   * hyperg_series(a, b, _Tp{1} - d, _Tp{1} - x)
		   / _Gamcma / _Gamcmb
		  + _Gamc * _Gammd
		   * hyperg_series(c - a, c - b, _Tp{1} + d, _Tp{1} - x)
		   / _Gama / _Gamb;
	    }
	}
      else
	return emsr::quiet_NaN(x);
    }

/**
 * Test harness.
 */
template<typename _Tp>
  void
  test_hyperg(_Tp proto = _Tp{})
  {
    //using _Val = _Tp;
    //using _Real = emsr::num_traits_t<_Val>;
    std::vector<_Tp> parm{_Tp{0.25}, _Tp{0.5}, _Tp{1}, _Tp{2}, _Tp{5}};

    std::cout.precision(emsr::digits10(proto));
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
		  auto hyperg0 = emsr::hyperg(a, b, c, z);
		  auto hyperg = ::hyperg(a, b, c, z);
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
    std::cout.precision(emsr::digits10(proto));
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
		    auto gnu = emsr::hyperg(_Tp(-i), _Tp(-n + j), _Tp(1 - i + j), x);
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
    using namespace std::complex_literals;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    using cmplx = std::complex<_Tp>;
    constexpr auto s_pi_3 = emsr::pi_v<_Tp> / _Tp{3};
    const auto z0p = std::polar(_Tp{1}, +s_pi_3);
    const auto z0m = std::polar(_Tp{1}, -s_pi_3);

    cmplx a = 1, b = 2, c = 3;
    auto f0p = hyperg_buhring(a, b, c, z0p);
    std::cout << ' ' << z0p << ' ' << f0p << '\n';
    auto f0m = hyperg_buhring(a, b, c, z0m);
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
		    auto gnu = emsr::hyperg(a, b, c, z);
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
    hyperg_series(a, b, c, z0p);
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
