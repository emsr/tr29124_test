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
#include <emsr/solver_madsen_reid.h>
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
   * @tparam _Tp The real type of the argument and degree parameters.
   * @param[in]  n  The degree of the Gegenbauer polynomial
   * @param[in]  lambda  The order of the Gegenbauer polynomial
   */
  template<typename _Tp>
    emsr::Polynomial<_Tp>
    gegenbauer_poly(unsigned int n, _Tp lambda)
    {
      emsr::Polynomial<_Tp> poly;

      if (std::isnan(lambda))
	return poly;

      auto term = emsr::Polynomial<_Tp>{1};
      poly += term;
      if (n == 0)
	return poly;

      const auto _2lambda = _Tp{2} * lambda;
      const auto lambdaph = lambda + _Tp{1} / _Tp{2};

      auto m = int(n);
      if (const auto pint = emsr::fp_is_integer(n + _2lambda);
	  pint && pint() <= 0 && -pint() < m)
	m = -pint();

      const emsr::Polynomial<_Tp> arg({_Tp{1}/_Tp{2}, _Tp{-1}/_Tp{2}});

      auto fact = _Tp{1};
      for (unsigned int k = 1; k <= n; ++k)
	fact *= _Tp(_2lambda + _Tp(k - 1)) / _Tp(k);

      for (int k = 1; k <= m; ++k)
	{
	  const auto km1 = _Tp(k - 1);

	  term *= (_Tp(-int(n) + km1) / _Tp(k))
		  * (_Tp(n + _2lambda + km1) / (lambdaph + km1))
		  * arg;

	  poly += term;
	}

      return fact * poly;
    }

  /**
   * Highest degree term coefficient.
   */
  template<typename _Tp>
    _Tp
    gegenbauer_norm(unsigned int n, _Tp lambda)
    {
      int sgam1, sgam2;
      const auto lgam1 = lgamma_r(_Tp(2 * n + 2 * lambda), &sgam1);
      const auto lgam2 = lgamma_r(_Tp(n + 2 * lambda), &sgam2);
      return sgam1 * sgam2 * std::exp(lgam1 - std::lgamma(_Tp(n + 1))
   				    - lgam2 - _Tp(n) * std::log(_Tp{2}));
    }

} // namespace detail
} // namespace emsr

/**
 * Write roots.
 */
template<typename Tp>
  void
  write_roots(unsigned n, Tp lambda, const std::vector<Tp>& coef,
              const emsr::Polynomial<Tp>& poly,
              const std::vector<emsr::solution_t<Tp>>& roots, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; lambda = " << lambda
	      << '\n';

    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << emsr::detail::gegenbauer_norm(n, lambda) << '\n';
    std::cout << std::flush;

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
 * Write roots.
 */
template<typename Tp>
  void
  write_roots(unsigned n, Tp lambda, const std::vector<std::complex<Tp>>& coef,
              const emsr::Polynomial<Tp>& poly,
              const std::vector<std::complex<Tp>>& roots, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; lambda = " << lambda
	      << '\n';

    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << emsr::detail::gegenbauer_norm(n, lambda) << '\n';
    std::cout << std::flush;

    std::cout << "\nThe roots are:\n";
    for (const auto& z : roots)
      {
	std::cout << ' ' << std::setw(w) << std::real(z)
		  << ' ' << std::setw(w) << std::imag(z)
		  << ' ' << std::setw(w) << poly(z)
		  << '\n';
      }

    gp << std::setprecision(prec);
    gp << "\n\n";
    gp << "# n = " << n
       << "; lambda = " << lambda << '\n';
    for (const auto& z : roots)
      {
	  gp << ' ' << std::setw(w) << std::real(z)
	     << ' ' << std::setw(w) << std::imag(z)
	     << '\n';
      }
    gp << std::flush;
  }

/**
 * Test negative parm roots.
 */
template<typename Tp>
  void
  test_neg_parm_gegenbauer_roots(unsigned n, Tp lambda, std::ofstream& gp)
  {
    const auto poly = emsr::detail::gegenbauer_poly(n, lambda);
    auto coef = poly.coefficients();

    std::reverse(coef.begin(), coef.end());

    std::cout << "\n  Jenkins-Traub solver...\n";
    gp << "\n#  Jenkins-Traub solver...";
    auto jt = emsr::JenkinsTraubSolver(coef);
    auto roots_jt = jt.solve();
    write_roots(n, lambda, coef, poly, roots_jt, gp);

    // Try complex solvers...
    std::vector<std::complex<Tp>> ccoef(coef.size());
    for (size_t c = 0; c < coef.size(); ++c)
      ccoef[c] = coef[c]; 

    std::cout << "\n  Complex Jenkins-Traub solver...\n";
    gp << "\n#  Complex Jenkins-Traub solver...";
    auto jtc = emsr::JenkinsTraubSolver(ccoef);
    auto roots_jtc = jtc.solve();
    write_roots(n, lambda, ccoef, poly, roots_jtc, gp);

    std::cout << "\n  Madsen-Reid solver...\n";
    gp << "\n#  Madsen-Reid solver...";
    auto mr = SolverMadsenReid(ccoef);
    auto roots_mr = mr.solve();
    write_roots(n, lambda, ccoef, poly, roots_mr, gp);
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
