/**
 *
 */

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <cmath> // FIXME: For isnan for math_util.h
#include <ext/math_util.h>
#include <ext/solver_jenkins_traub.h>
#include <ext/solver_madsen_reid.h>
#include <ext/polynomial.h>

#include <ext/float128_io.h>
#include <ext/float128_math.h>

#include <sf_jacobi_neg_params.tcc>

namespace std
{
namespace __detail
{

  /**
   * Return the Jacobi polynomial as a polynomial:
   * @f[
   *    P_n^{(\alpha,\beta)}(x) = \frac{(\alpha + 1)_n}{n!}
   *      {}_2F_1(-n, n + 1 + \alpha + \beta; \alpha + 1; \frac{1-x}{2})
   * @f]
   *
   * @tparam _Tp The real type of the argument and degree parameters.
   * @param[in]  n  The degree of the Jacobi polynomial
   * @param[in]  alpha1  The first parameter of the Jacobi polynomial
   * @param[in]  beta1  The second parameter of the Jacobi polynomial
   */
  template<typename _Tp>
    __gnu_cxx::_Polynomial<_Tp>
    __jacobi_poly(unsigned int __n, _Tp __alpha1, _Tp __beta1)
    {
      __gnu_cxx::_Polynomial<_Tp> __poly;

      if (std::isnan(__alpha1) || std::isnan(__beta1))
	return __poly;

      auto __term = __gnu_cxx::_Polynomial<_Tp>{1};
      __poly += __term;
      if (__n == 0)
	return __poly;

      const auto __apb = __alpha1 + __beta1;

      auto __m = int(__n);
      if (const auto __pint = __gnu_cxx::__fp_is_integer(__n + 1 + __apb);
	  __pint && __pint() <= 0 && -__pint() < __m)
	__m = -__pint();

      const __gnu_cxx::_Polynomial<_Tp> __arg({_Tp{1}/_Tp{2}, _Tp{-1}/_Tp{2}});

      auto __fact = _Tp{1};
      for (unsigned int __k = 1; __k <= __n; ++__k)
	__fact *= _Tp(__alpha1 + __k) / _Tp(__k);

      for (int __k = 1; __k <= __m; ++__k)
	{

	  __term *= (_Tp(-int(__n) + __k - 1) / _Tp(__k))
		  * (_Tp(__n + __k + __apb) / _Tp(__alpha1 + __k))
		  * __arg;

	  __poly += __term;
	}

      return __fact * __poly;
    }

  /**
   * Highest degree term coefficient.
   */
  template<typename _Tp>
    _Tp
    __jacobi_norm(unsigned int __n, _Tp __alpha1, _Tp __beta1)
    {
      int sgam1, sgam2;
      const auto lgam1 = lgamma_r(_Tp(2 * __n + __alpha1 + __beta1 + 1), &sgam1);
      const auto lgam2 = lgamma_r(_Tp(__n + __alpha1 + __beta1 + 1), &sgam2);
      return sgam1 * sgam2 * std::exp(lgam1 - std::lgamma(_Tp(__n + 1))
   				    - lgam2 - _Tp(__n) * std::log(_Tp{2}));
    }

} // namespace std
} // namespace __detail

/**
 * Write roots.
 */
template<typename Tp>
  void
  write_roots(unsigned n, Tp alpha, Tp beta, const std::vector<Tp>& coef,
              const __gnu_cxx::_Polynomial<Tp>& poly,
              const std::vector<__gnu_cxx::solution_t<Tp>>& roots, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; alpha = " << alpha << "; beta = " << beta
	      << '\n';

    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << std::__detail::__jacobi_norm(n, alpha, beta) << '\n';
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
       << "; alpha = " << alpha << "; beta = " << beta << '\n';
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
  write_roots(unsigned n, Tp alpha, Tp beta, const std::vector<std::complex<Tp>>& coef,
              const __gnu_cxx::_Polynomial<Tp>& poly,
              const std::vector<std::complex<Tp>>& roots, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; alpha = " << alpha << "; beta = " << beta
	      << '\n';

    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << std::__detail::__jacobi_norm(n, alpha, beta) << '\n';
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
       << "; alpha = " << alpha << "; beta = " << beta << '\n';
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
  test_neg_parm_jacobi_roots(unsigned n, Tp alpha1, Tp beta1, std::ofstream& gp)
  {
    const auto poly = std::__detail::__jacobi_poly(n, alpha1, beta1);
    auto coef = poly.coefficients();

    std::reverse(coef.begin(), coef.end());

    std::cout << "\n  Jenkins-Traub solver...\n";
    gp << "\n#  Jenkins-Traub solver...";
    auto jt = __gnu_cxx::_JenkinsTraubSolver(coef);
    auto roots_jt = jt.solve();
    write_roots(n, alpha1, beta1, coef, poly, roots_jt, gp);

    // Try complex solvers...
    std::vector<std::complex<Tp>> ccoef(coef.size());
    for (size_t c = 0; c < coef.size(); ++c)
      ccoef[c] = coef[c];

    std::cout << "\n  Complex Jenkins-Traub solver...\n";
    gp << "\n#  Complex Jenkins-Traub solver...";
    auto jtc = __gnu_cxx::_JenkinsTraubSolver(ccoef);
    auto roots_jtc = jtc.solve();
    write_roots(n, alpha1, beta1, ccoef, poly, roots_jtc, gp);

    std::cout << "\n  Madsen-Reid solver...\n";
    gp << "\n#  Madsen-Reid solver...";
    auto mr = SolverMadsenReid(ccoef);
    auto roots_mr = mr.solve();
    write_roots(n, alpha1, beta1, ccoef, poly, roots_mr, gp);
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
__gnu_cxx::jacobi(10, 2.0, -12.0, -1.0);
    std::cout << "\n\n";
    for (int n : {10, 15, 20})
      for (Tp alpha : {1, 2})
	for (Tp beta : {1, 2})
	  {
	    auto P = std::__detail::__jacobi_poly(n, alpha, beta);
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
			  << ' ' << std::setw(w) << __gnu_cxx::jacobi(n, alpha, beta, x)
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
	  auto P = std::__detail::__jacobi_poly(n, alpha, beta);
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
			//<< ' ' << std::setw(w) << __gnu_cxx::jacobi(n, alpha, beta, x)
			<< ' ' << std::setw(w) << lab::__jacobi_recur(n, alpha, beta, x).__P_n
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
