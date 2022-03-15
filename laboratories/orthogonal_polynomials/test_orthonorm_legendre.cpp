
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>

#include <emsr/integration.h>
#include <emsr/special_functions.h>

// Normalized Legendre polynomial.
template<typename Tp>
  Tp
  normalized_legendre(int l, Tp x)
  {
    return emsr::legendre(l, x) * std::sqrt(Tp(2 * l + 1) / Tp{2});
  }

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int l1, int l2, Tp x)
  {
    return normalized_legendre(l1,x)
	 * normalized_legendre(l2,x);
  }

template<typename Tp>
  Tp
  delta(int l1, int l2)
  { return l1 == l2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_legendre()
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    int l1 = 0;
    for (; l1 <= 128; ++l1)
      {
	for (int l2 = 0; l2 <= l1; ++l2)
	  {
	    auto func = [l1, l2](Tp x)
			-> Tp
			{ return integrand(l1, l2, x); };

	    auto [result, error]
		//= emsr::integrate(func, Tp{-1}, Tp{1},
		//			abs_precision, rel_precision);
		= emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1},
						 abs_precision, rel_precision, 6);

	    if (std::abs(delta<Tp>(l1, l2) - result) > cmp_precision)
	      {
		std::stringstream ss;
		ss.precision(std::numeric_limits<Tp>::digits10);
		ss << std::showpoint << std::scientific;
		ss << "Integration failed at l1=" << l1 << ", l2=" << l2
		   << ", returning result " << result
		   << ", with error " << error
		   << " instead of the expected " << delta<Tp>(l1, l2) << '\n';
		throw std::logic_error(ss.str());
	      }
	  }
	std::cout << "Integration successful for legendre polynomials up to l = " << l1
		  << '\n' << std::flush;
      }

    int n1_lower = l1 - 1;
    int n1_upper = 2 * n1_lower;
    int del = 2;
    bool breakout = false;
    while (n1_upper != n1_lower)
      {
	RESTART:
	for (int l2 = 0; l2 <= n1_upper; l2 += del)
	  {
	    auto func = [l1 = n1_upper, l2](Tp x)
			-> Tp
			{ return integrand(l1, l2, x); };

	    auto [result, error]
		//= emsr::integrate(func, Tp{-1}, Tp{1},
		//			abs_precision, rel_precision);
		= emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1},
						 abs_precision, rel_precision, 6);

	    if (std::abs(delta<Tp>(n1_upper, l2) - result) > cmp_precision)
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

	std::cout << "Integration successful for legendre polynomials up to l = " << n1_upper
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
      test_legendre<float>();
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
      test_legendre<double>();
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
      test_legendre<long double>();
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
