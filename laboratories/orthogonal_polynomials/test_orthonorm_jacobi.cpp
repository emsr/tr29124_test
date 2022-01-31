
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>

#include <emsr/integration.h>
#include <emsr/special_functions.h>

// Try to manage the four-gamma ratio.
// alpha > -1, beta > -1.
template<typename Tp>
  Tp
  gamma_ratio(int n, Tp alpha, Tp beta)
  {
    const auto _S_eps = std::numeric_limits<Tp>::epsilon();
    if (std::abs(Tp(1) + alpha + beta) < _S_eps)
      return Tp(0);
    else
      {
	auto gaman1 = std::tgamma(Tp(1) + alpha);
	auto gambn1 = std::tgamma(Tp(1) + beta);
	auto gamabn1 = std::tgamma(Tp(1) + alpha + beta);
	auto fact = gaman1 * gambn1 / gamabn1;
	for (int k = 1; k <= n; ++k)
	  fact *= (Tp(k) + alpha) * (Tp(k) + beta)
		/ (Tp(k) + alpha + beta) / Tp(k);
	return fact;
      }
  }

// Normalized Jacobi polynomial.
template<typename Tp>
  Tp
  normalized_jacobi(int n, Tp alpha, Tp beta, Tp x)
  {
    auto gam = gamma_ratio(n, alpha, beta);
    auto norm = std::sqrt(std::pow(Tp{2}, Tp{1} + alpha + beta)
	                * gam / (Tp(2 * n + 1) + alpha + beta));
    return emsr::jacobi(n, alpha, beta, x) / norm;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp alpha, Tp beta, Tp x)
  {
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    if (std::abs(x - Tp{1}) < s_eps)
      return Tp{0};
    else if (std::abs(x + Tp{1}) < s_eps)
      return Tp{0};
    else
      {
	return std::pow(Tp{1} - x, alpha) * std::pow(Tp{1} + x, beta)
	     * normalized_jacobi(n1, alpha, beta, x)
	     * normalized_jacobi(n2, alpha, beta, x);
      }
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_jacobi(Tp alpha, Tp beta)
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    const bool singular = (alpha < Tp{0} || beta < Tp{0});

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2, alpha, beta](Tp x)
			-> Tp
			{ return integrand<Tp>(n1, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					       Tp{-1}, Tp{1},
					       alpha, beta, 0, 0,
					       abs_precision, rel_precision)
		//: emsr::integrate(func, Tp{-1}, Tp{1},
		//			abs_precision, rel_precision);
		: emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1},
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
	std::cout << "Integration successful for jacobi polynomials up to n = " << n1
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
	    auto func = [n1_upper, n2, alpha, beta](Tp x)
			-> Tp
			{ return integrand<Tp>(n1_upper, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					       Tp{-1}, Tp{1},
					       alpha, beta, 0, 0,
					       abs_precision, rel_precision)
		//: emsr::integrate(func, Tp{-1}, Tp{1},
		//			abs_precision, rel_precision);
		: emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1},
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

	std::cout << "Integration successful for jacobi polynomials up to n = " << n1_upper
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
      test_jacobi<float>(0.5F, 1.5F);
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
      test_jacobi<double>(0.5, 1.5);
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
      test_jacobi<long double>(0.5L, 1.5L);
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
