
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>

#include <emsr/integration.h>
#include <emsr/special_functions.h>

// Chebyshev polynomials of the first kind are already normalized.

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp x)
  {
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    const auto s_inf = std::numeric_limits<Tp>::infinity();
    if (std::abs(x - Tp{1}) < s_eps)
      return s_inf;
    else if (std::abs(x + Tp{1}) < s_eps)
      return ((n1 + n2) & 1) ? -s_inf : s_inf;
    else
      return emsr::chebyshev_t(n2, x)
	   * emsr::chebyshev_t(n1, x)
	   / std::sqrt(Tp{1} - x * x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_chebyshev_t()
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2](Tp x)->Tp{return integrand(n1, n2, x);};

	    auto [result, error]
		//= emsr::integrate(func, Tp{-1}, Tp{1},
		//			abs_precision, rel_precision);
		//= emsr::integrate_singular(func, Tp{-1}, Tp{1},
		//				abs_precision, rel_precision);
		= emsr::integrate_singular_endpoints(func,
				 Tp{-1}, Tp{1},
				 Tp{-0.5}, Tp{-0.5}, 0, 0,
				 abs_precision, rel_precision);

	    if (std::abs(delta<Tp>(n1, n2) - result) > cmp_precision)
	      {
		std::stringstream ss;
		ss.precision(std::numeric_limits<Tp>::digits10);
		ss << std::showpoint << std::scientific;
		ss << "Integration failed at n1=" << n1 << ", n2=" << n2
		   << ", returning result " << result
		   << ", with error " << error
		   << " instead of the expected " << delta<Tp>(n1, n2) << '\n';
		throw std::logic_error(ss.str());
	      }
	  }
	std::cout << "Integration successful for chebyshev_t polynomials up to n = " << n1
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
	    auto func = [n1 = n1_upper, n2](Tp x)->Tp{return integrand(n1, n2, x);};

	    auto [result, error]
		//= emsr::integrate(func, Tp{-1}, Tp{1},
		//			abs_precision, rel_precision);
		//= emsr::integrate_singular(func, Tp{-1}, Tp{1},
		//				abs_precision, rel_precision);
		= emsr::integrate_singular_endpoints(func,
				 Tp{-1}, Tp{1},
				 Tp{-0.5}, Tp{-0.5}, 0, 0,
				 abs_precision, rel_precision);

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

	std::cout << "Integration successful for chebyshev_t polynomials up to n = " << n1_upper
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
      test_chebyshev_t<float>();
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
      test_chebyshev_t<double>();
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
      test_chebyshev_t<long double>();
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
