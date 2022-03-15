
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include <emsr/integration.h>
#include <emsr/special_functions.h>

// Normalized radial polynomial.
template<typename Tp>
  Tp
  normalized_radpoly(int n, int m, Tp rho)
  {
    return std::sqrt(Tp(2 * n + 2)) * emsr::radpoly(n, m, rho);
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int m1, int n2, int m2, Tp rho)
  {
    return rho
	 * normalized_radpoly(n1, m1, rho)
	 * normalized_radpoly(n2, m2, rho);
  }

template<typename Tp>
  Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_radpoly()
  {
    bool full_range = false;

    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    std::vector<int> degree{0,1,2,3,4,5,8,9,16,17,32,33,64,65,128,129};

    int n1_save = 0;
    for (int n1 : degree)
      {
	for (int m1 : degree)
	  {
	    if (m1 > n1 || (n1 - m1) & 1)
	      continue;
	    for (int n2 : degree)
	      {
		if (n2 > n1)
		  continue;
		// The orthonormality only works for m2 == m1.
		int m2 = m1;
		if (m2 > n2 || (n2 - m2) & 1)
		  continue;
		auto func = [n1, m1, n2, m2](Tp x)
			    -> Tp
			    { return integrand(n1, m1, n2, m2, x); };

		auto [result, error]
			= emsr::integrate_tanh_sinh(func, Tp{0}, Tp{1},
					      abs_precision, rel_precision, 8);

		if (std::abs(delta<Tp>(n1, m1, n2, m2) - result) > cmp_precision)
		  {
		    std::stringstream ss;
		    ss.precision(std::numeric_limits<Tp>::digits10);
		    auto w = 8 + ss.precision();
		    ss << std::showpoint << std::scientific;
		    ss << "Integration failed at n1=" << n1 << ", m1=" << m1
		       << ", n2=" << n2 << ", m2=" << m2
		       << ", returning result = " << std::setw(w) << result
		       << ", with error = " << std::setw(w) << error
		       << " instead of the expected " << delta<Tp>(n1, m1, n2, m2) << '\n';
		    std::cerr << ss.str();
		  }
	      }
	  }
	n1_save = n1;
	std::cout << "Integration successful for radpoly polynomials up to n = " << n1_save
		  << '\n' << std::flush;
      }

    if (!full_range)
      return;

    int n1_lower = n1_save;
    int n1_upper = 2 * n1_lower;
    int del = 2;
    bool breakout = false;
    while (n1_upper != n1_lower)
      {
	RESTART:
	for (int m1 = 0; m1 <= n1_upper; ++m1)
	  {
	    if ((n1_upper - m1) & 1)
	      continue;
	    for (int n2 = 0; n2 <= n1_upper; n2 += del)
	      {
		// The orthonormality only works for m2 == m1.
		int m2 = m1;
		if (m2 > n2 || (n2 - m2) & 1)
		  continue;
		auto func = [n1 = n1_upper, m1, n2, m2](Tp x)
			    -> Tp
			    { return integrand(n1, m1, n2, m2, x); };

		auto [result, error]
			= emsr::integrate_tanh_sinh(func, Tp{0}, Tp{1},
					      abs_precision, rel_precision, 8);

		if (std::abs(delta<Tp>(n1_upper, m1, n2, m2) - result) > cmp_precision)
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
	    if (breakout)
	      break;
	  }

	std::cout << "Integration successful for radpoly polynomials up to n = " << n1_upper
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
      test_radpoly<float>();
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
      test_radpoly<double>();
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
      test_radpoly<long double>();
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
