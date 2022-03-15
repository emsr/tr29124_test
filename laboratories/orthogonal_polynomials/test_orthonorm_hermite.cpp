
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>

#include <emsr/integration.h>
#include <emsr/math_constants.h>
#include <emsr/special_functions.h>
#include <emsr/math_constants.h>

// Normalized Hermite polynomial.
template<typename Tp>
  Tp
  normalized_hermite(int n, Tp x)
  {
    constexpr auto s_lnpi = emsr::lnpi_v<Tp>;
    constexpr auto s_ln2 = emsr::ln2_v<Tp>;
    const auto lnorm = Tp{0.5} * (Tp{0.5} * s_lnpi + Tp(n) * s_ln2 + emsr::lfactorial<Tp>(n));
    return emsr::hermite(n, x) * std::exp(-lnorm);
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp x)
  {
    return std::exp(-x * x)
	 * normalized_hermite(n1, x)
	 * normalized_hermite(n2, x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_hermite()
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    //const auto infty = std::numeric_limits<Tp>::infinity();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2](Tp x)
			-> Tp
			{ return integrand(n1, n2, x); };

	    auto [result, error]
		//= emsr::integrate_singular(func, -infty, infty, abs_precision, rel_precision);
		//= emsr::integrate(func, -infty, infty, abs_precision, rel_precision);
		//= emsr::integrate_infinite(func, abs_precision, rel_precision);
		//= emsr::integrate_singular_infinite(func, abs_precision, rel_precision);
		//= emsr::integrate_oscillatory(map_minf_pinf<Tp, func_t>(func),
		//			  Tp{0}, Tp{1}, abs_precision, rel_precision);
		//= emsr::integrate_clenshaw_curtis(map_minf_pinf<Tp, func_t>(func),
		//			      Tp{0}, Tp{1}, abs_precision, rel_precision);
		= emsr::integrate_sinh_sinh(func, abs_precision, rel_precision);

	    if (std::abs(delta<Tp>(n1, n2) - result) > cmp_precision)
	      {
		std::stringstream ss;
		ss.precision(std::numeric_limits<Tp>::digits10);
		ss << std::showpoint << std::scientific;
		ss << "Integration failed at n1=" << n1 << ", n2=" << n2
		   << ", returning result " << result
		   << ", with error " << error
		   << " instead of the expected " << delta<Tp>(n1, n2)
		   << " with absolute precision " << cmp_precision << '\n';
		throw std::logic_error(ss.str());
	      }
	  }
	std::cout << "Integration successful for hermite polynomials up to n = " << n1
		  << '\n' << std::flush;
      }

    int n1_lower = n1 - 1;
    int n1_upper = 2 * n1_lower;
    int del = 2;
    bool breakout = false;
    while (n1_upper != n1_lower)
      {
	RESTART:
	for (int n2 = 0; n2 <= n1_upper; n2 += del)
	  {
	    auto func = [n1 = n1_upper, n2](Tp x)
			-> Tp
			{ return integrand(n1, n2, x); };

	    auto [result, error]
		//= emsr::integrate_singular(func, -infty, infty,
		//				abs_precision, rel_precision);
		= emsr::integrate_sinh_sinh(func, abs_precision, rel_precision);

	    if (std::abs(delta<Tp>(n1_upper, n2) - result) > cmp_precision)
	      {
		if ((n1_lower + n1_upper) / 2 < n1_upper)
		  {
		    n1_upper = (n1_lower + n1_upper) / 2;
		    goto RESTART;
		  }
		else
		  {
		    breakout = true;
		    break;
		  }
	      }
	  }

	std::cout << "Integration successful for hermite polynomials up to n = " << n1_upper
		  << '\n' << std::flush;

	if (breakout)
	  break;

	n1_lower = n1_upper;
	if (n1_upper > 1000)
	  {
	    std::cout << "\nGood enough!\n" << std::flush;
	    break;
	  }
	else if (n1_upper <= std::numeric_limits<int>::max() / 2)
	  n1_upper *= 2;
	else
	  break;
        del *= 2;
      }
  }

int
main()
{
  std::cout << "\n\nOrthonormality tests for float\n";
  try
    {
      test_hermite<float>();
    }
  catch (emsr::integration_error<float, float>& ierr)
    {
      std::cerr << ierr.what() << '\n';
      std::cerr << " result = " << ierr.result() << " abserr = " << ierr.abserr() << '\n';
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for double\n";
  try
    {
      test_hermite<double>();
    }
  catch (emsr::integration_error<double, double>& ierr)
    {
      std::cerr << ierr.what() << '\n';
      std::cerr << " result = " << ierr.result() << " abserr = " << ierr.abserr() << '\n';
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for long double\n";
  try
    {
      test_hermite<long double>();
    }
  catch (emsr::integration_error<long double, long double>& ierr)
    {
      std::cerr << ierr.what() << '\n';
      std::cerr << " result = " << ierr.result() << " abserr = " << ierr.abserr() << '\n';
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
