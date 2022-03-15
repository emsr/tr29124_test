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
#include <ext/polynomial.h>

#include <ext/float128_io.h>
#include <ext/float128_math.h>

#include <sf_gegenbauer_neg_params.tcc>

namespace std
{
namespace __detail
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
    __gnu_cxx::_Polynomial<_Tp>
    __gegenbauer_poly(unsigned int __n, _Tp __lambda)
    {
      __gnu_cxx::_Polynomial<_Tp> __poly;

      if (std::isnan(__lambda))
	return __poly;

      auto __term = __gnu_cxx::_Polynomial<_Tp>{1};
      __poly += __term;
      if (__n == 0)
	return __poly;

      const auto __2lambda = _Tp{2} * __lambda;
      const auto __lambdaph = __lambda + _Tp{1} / _Tp{2};

      auto __m = int(__n);
      if (const auto __pint = __gnu_cxx::__fp_is_integer(__n + __2lambda);
	  __pint && __pint() <= 0 && -__pint() < __m)
	__m = -__pint();

      const __gnu_cxx::_Polynomial<_Tp> __arg({_Tp{1}/_Tp{2}, _Tp{-1}/_Tp{2}});

      auto __fact = _Tp{1};
      for (unsigned int __k = 1; __k <= __n; ++__k)
	__fact *= _Tp(__2lambda + _Tp(__k - 1)) / _Tp(__k);

      for (int __k = 1; __k <= __m; ++__k)
	{
	  const auto __km1 = _Tp(__k - 1);

	  __term *= (_Tp(-int(__n) + __km1) / _Tp(__k))
		  * (_Tp(__n + __2lambda + __km1) / (__lambdaph + __km1))
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
    __gegenbauer_norm(unsigned int __n, _Tp __lambda)
    {
      int sgam1, sgam2;
      const auto lgam1 = lgamma_r(_Tp(2 * __n + 2 * __lambda), &sgam1);
      const auto lgam2 = lgamma_r(_Tp(__n + 2 * __lambda), &sgam2);
      return sgam1 * sgam2 * std::exp(lgam1 - std::lgamma(_Tp(__n + 1))
   				    - lgam2 - _Tp(__n) * std::log(_Tp{2}));
    }

} // namespace std
} // namespace __detail

template<typename _Tp>
  void
  test_neg_parm_gegenbauer_roots(unsigned n, _Tp lambda, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<_Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; lambda = " << lambda
	      << '\n';

    const auto poly = std::__detail::__gegenbauer_poly(n, lambda);
    auto coef = poly.coefficients();
    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << std::__detail::__gegenbauer_norm(n, lambda) << '\n';
    std::cout << std::flush;

    std::reverse(coef.begin(), coef.end());

    auto jt = __gnu_cxx::_JenkinsTraubSolver(coef);
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
template<typename _Tp>
  void
  run()
  {
    std::ofstream gp("gegenbauer_roots.gp");

    unsigned n = 50;
    _Tp lambda;

    lambda = _Tp(1 - int(n)) / _Tp{2};
    test_neg_parm_gegenbauer_roots(n, lambda, gp);
  }

template<typename _Tp>
  void
  test_poly()
  {
    const auto prec = std::numeric_limits<_Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";

    for (int n : {10, 15, 20})
      for (int m : {0, 1, 2, 3})
	{
	  const auto lambda = (1 - n - m) / _Tp{2};
	  auto P = std::__detail::__gegenbauer_poly(n, lambda);
	  std::cout << " n = " << n
		    << "; lambda = " << lambda
		    << '\n';
	  for (int i = -10; i <= +10; ++i)
	    {
	      const auto x = _Tp(i * 0.1L);
	      std::cout << ' ' << x
			<< ' ' << std::setw(w) << P(x)
			//<< ' ' << std::setw(w) << __gnu_cxx::gegenbauer(n, lambda, x).__C_n
			<< ' ' << std::setw(w) << lab::__gegenbauer_recur(n, lambda, x).__C_n
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
