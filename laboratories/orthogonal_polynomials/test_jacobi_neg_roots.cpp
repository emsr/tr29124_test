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

#include <sf_jacobi_neg_params.tcc>

namespace emsr
{
namespace detail
{

  /**
   * Return the Jacobi polynomial as a polynomial:
   * @f[
   *    P_n^{(\alpha,\beta)}(x) = \frac{(\alpha + 1)_n}{n!}
   *      {}_2F_1(-n, n + 1 + \alpha + \beta; \alpha + 1; \frac{1-x}{2})
   * @f]
   *
   * @tparam Tp The real type of the argument and degree parameters.
   * @param[in]  n  The degree of the Jacobi polynomial
   * @param[in]  alpha1  The first parameter of the Jacobi polynomial
   * @param[in]  beta1  The second parameter of the Jacobi polynomial
   */
  template<typename Tp>
    emsr::Polynomial<Tp>
    jacobi_poly(unsigned int n, Tp alpha1, Tp beta1)
    {
      emsr::Polynomial<Tp> poly;

      if (std::isnan(alpha1) || std::isnan(beta1))
	return poly;

      auto term = emsr::Polynomial<Tp>{1};
      poly += term;
      if (n == 0)
	return poly;

      const auto apb = alpha1 + beta1;

      auto m = int(n);
      if (const auto pint = emsr::fp_is_integer(n + 1 + apb);
	  pint && pint() <= 0 && -pint() < m)
	m = -pint();

      const emsr::Polynomial<Tp> arg({Tp{1}/Tp{2}, Tp{-1}/Tp{2}});

      auto fact = Tp{1};
      for (unsigned int k = 1; k <= n; ++k)
	fact *= Tp(alpha1 + k) / Tp(k);

      for (int k = 1; k <= m; ++k)
	{

	  term *= (Tp(-int(n) + k - 1) / Tp(k))
		  * (Tp(n + k + apb) / Tp(alpha1 + k))
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
    jacobi_norm(unsigned int n, Tp alpha1, Tp beta1)
    {
      int sgam1, sgam2;
      const auto lgam1 = lgamma_r(Tp(2 * n + alpha1 + beta1 + 1), &sgam1);
      const auto lgam2 = lgamma_r(Tp(n + alpha1 + beta1 + 1), &sgam2);
      return sgam1 * sgam2 * std::exp(lgam1 - std::lgamma(Tp(n + 1))
   				    - lgam2 - Tp(n) * std::log(Tp{2}));
    }

} // namespace detail
} // namespace emsr

template<typename Tp>
  void
  test_neg_parm_jacobi_roots(unsigned n, Tp alpha1, Tp beta1, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; alpha = " << alpha1
	      << "; beta = " << beta1
	      << '\n';

    const auto poly = emsr::detail::jacobi_poly(n, alpha1, beta1);
    auto coef = poly.coefficients();
    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << emsr::detail::jacobi_norm(n, alpha1, beta1) << '\n';
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
       << "; alpha = " << alpha1
       << "; beta = " << beta1 << '\n';
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
    std::ofstream gp("jacobi_roots.gp");

    unsigned n = 50;
    Tp alpha1, beta1;

    alpha1 = Tp{2};
    beta1 = Tp{-83} / Tp{2};
    test_neg_parm_jacobi_roots(n, alpha1, beta1, gp);

    alpha1 = Tp{2};
    beta1 = Tp{-52};
    test_neg_parm_jacobi_roots(n, alpha1, beta1, gp);

    alpha1 = Tp{2};
    beta1 = Tp{-127} / Tp{2};
    test_neg_parm_jacobi_roots(n, alpha1, beta1, gp);

    // Flip alpha and beta.

    alpha1 = Tp{2};
    beta1 = Tp{-83} / Tp{2};
    test_neg_parm_jacobi_roots(n, beta1, alpha1, gp);

    alpha1 = Tp{2};
    beta1 = Tp{-52};
    test_neg_parm_jacobi_roots(n, beta1, alpha1, gp);

    alpha1 = Tp{2};
    beta1 = Tp{-127} / Tp{2};
    test_neg_parm_jacobi_roots(n, beta1, alpha1, gp);
  }

/*
 * Test polynomial evaluations against good ol' recursion.
 */
template<typename Tp>
  void
  test_poly()
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
emsr::jacobi(10, 2.0, -12.0, -1.0);
    std::cout << "\n\n";
    for (int n : {10, 15, 20})
      for (Tp alpha : {1, 2})
	for (Tp beta : {1, 2})
	  {
	    auto P = emsr::detail::jacobi_poly(n, alpha, beta);
	    std::cout << " n = " << n
		      << "; alpha = " << alpha
		      << "; beta = " << beta
		      << "; P = " << P
		      << '\n';
	    for (int i = -10; i <= +10; ++i)
	      {
		const auto x = Tp(i * 0.1L);
		std::cout << ' ' << x
			  << ' ' << std::setw(w) << P(x)
			  << ' ' << std::setw(w) << emsr::jacobi(n, alpha, beta, x)
			  << '\n';
	      }
	    std::cout << '\n';
	  }

    std::cout << "\n\n";
    const auto alpha = Tp{2};
    for (int n : {10, 15, 20})
      for (int m : {0, 1, 2, 3})
	{
	  const auto beta = -n - m - alpha;
	  auto P = emsr::detail::jacobi_poly(n, alpha, beta);
	  std::cout << " n = " << n
		    << "; alpha = " << alpha
		    << "; beta = " << beta
		    << "; P = " << P
		    << '\n';
	  for (int i = -10; i <= +10; ++i)
	    {
	      const auto x = Tp(i * 0.1L);
	      std::cout << ' ' << x
			<< ' ' << std::setw(w) << P(x)
			//<< ' ' << std::setw(w) << emsr::jacobi(n, alpha, beta, x)
			<< ' ' << std::setw(w) << lab::jacobi_recur(n, alpha, beta, x).P_n
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
