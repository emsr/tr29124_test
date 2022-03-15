/**
 *
 */

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath> // FIXME: For isnan for math_util.h

#include <emsr/math_util.h>
#include <emsr/solver_jenkins_traub.h>
#include <emsr/polynomial.h>

#include <emsr/float128_io.h>
#include <emsr/float128_math.h>

#include <sf_gegenbauer_neg_params.tcc>

namespace emsr
{
namespace detail
{

  /**
   * Return the Gegenbauer polynomial as a polynomial:
   * @f[
   *    C_n^{(\lambda)}(x) = \frac{(2\lambda)_n}{n!}
   *      {}_2F_1(-n, n + 2\lambda; \lambda + 1/2; \frac{1-x}{2})
   * @f]
   *
   * @tparam Tp The real type of the argument and degree parameters.
   * @param[in]  n  The degree of the Gegenbauer polynomial
   * @param[in]  lambda  The order of the Gegenbauer polynomial
   */
  template<typename Tp>
    emsr::Polynomial<Tp>
    gegenbauer_poly(unsigned int n, Tp lambda)
    {
      emsr::Polynomial<Tp> poly;

      if (std::isnan(lambda))
	return poly;

      auto term = emsr::Polynomial<Tp>{1};
      poly += term;
      if (n == 0)
	return poly;

      const auto __2lambda = Tp{2} * lambda;
      const auto lambdaph = lambda + Tp{1} / Tp{2};

      auto m = int(n);
      if (const auto pint = emsr::fp_is_integer(n + __2lambda);
	  pint && pint() <= 0 && -pint() < m)
	m = -pint();

      const emsr::Polynomial<Tp> arg({Tp{1}/Tp{2}, Tp{-1}/Tp{2}});

      auto fact = Tp{1};
      for (unsigned int k = 1; k <= n; ++k)
	fact *= Tp(__2lambda + Tp(k - 1)) / Tp(k);

      for (int k = 1; k <= m; ++k)
	{
	  const auto km1 = Tp(k - 1);

	  term *= (Tp(-int(n) + km1) / Tp(k))
		  * (Tp(n + __2lambda + km1) / (lambdaph + km1))
		  * arg;

	  poly += term;
	}

      return fact * poly;
    }

  /**
   * Highest degree term coefficient.
   */
  template<typename Tp>
    Tp
    gegenbauer_norm(unsigned int n, Tp lambda)
    {
      int sgam1, sgam2;
      const auto lgam1 = lgamma_r(Tp(2 * n + 2 * lambda), &sgam1);
      const auto lgam2 = lgamma_r(Tp(n + 2 * lambda), &sgam2);
      return sgam1 * sgam2 * std::exp(lgam1 - std::lgamma(Tp(n + 1))
   				    - lgam2 - Tp(n) * std::log(Tp{2}));
    }

} // namespace detail
} // namespace emsr

template<typename Tp>
  void
  test_neg_parm_gegenbauer_roots(unsigned n, Tp lambda, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; lambda = " << lambda
	      << '\n';

    const auto poly = emsr::detail::gegenbauer_poly(n, lambda);
    auto coef = poly.coefficients();
    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << emsr::detail::gegenbauer_norm(n, lambda) << '\n';
    std::cout << std::flush;

    std::reverse(coef.begin(), coef.end());

    auto jt = emsr::JenkinsTraubSolver(coef);
    auto roots = jt.solve();
    std::cout << "\nThe roots are:\n";
    for (const auto& z : roots)
      {
	if (z.index() == 0)
	  continue;
	else if (z.index() == 1)
	  std::cout << ' ' << std::setw(w) << std::get<1>(z)
		    << ' ' << std::setw(w) << 0.0
		    << ' ' << std::setw(w) << poly(std::get<1>(z))
		    << '\n';
	else
	  std::cout << ' ' << std::setw(w) << std::real(std::get<2>(z))
		    << ' ' << std::setw(w) << std::imag(std::get<2>(z))
		    << ' ' << std::setw(w) << poly(std::get<2>(z))
		    << '\n';
      }

    gp << std::setprecision(prec);
    gp << "\n\n";
    gp << "# n = " << n
       << "; lambda = " << lambda << '\n';
    for (const auto& z : roots)
      {
	if (z.index() == 0)
	  continue;
	else if (z.index() == 1)
	  gp << ' ' << std::setw(w) << std::get<1>(z)
	     << ' ' << std::setw(w) << 0.0
	     << '\n';
	else
	  gp << ' ' << std::setw(w) << std::real(std::get<2>(z))
	     << ' ' << std::setw(w) << std::imag(std::get<2>(z))
	     << '\n';
      }
    gp << std::flush;
  }

/**
 * Numerical Methods for Special Functions, Gil, Segura, Temme, pp. 192.
 */
template<typename Tp>
  void
  run()
  {
    std::ofstream gp("gegenbauer_roots.gp");

    unsigned n = 50;
    Tp lambda;

    lambda = Tp(1 - int(n)) / Tp{2};
    test_neg_parm_gegenbauer_roots(n, lambda, gp);
  }

template<typename Tp>
  void
  test_poly()
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";

    for (int n : {10, 15, 20})
      for (int m : {0, 1, 2, 3})
	{
	  const auto lambda = (1 - n - m) / Tp{2};
	  auto P = emsr::detail::gegenbauer_poly(n, lambda);
	  std::cout << " n = " << n
		    << "; lambda = " << lambda
		    << '\n';
	  for (int i = -10; i <= +10; ++i)
	    {
	      const auto x = Tp(i * 0.1L);
	      std::cout << ' ' << x
			<< ' ' << std::setw(w) << P(x)
			//<< ' ' << std::setw(w) << emsr::gegenbauer(n, lambda, x).C_n
			<< ' ' << std::setw(w) << lab::gegenbauer_recur(n, lambda, x).C_n
			<< '\n';
	    }
	  std::cout << '\n';
	}
  }

int
main()
{
  run<long double>();

  //run<__float128>();

  test_poly<double>();
}
