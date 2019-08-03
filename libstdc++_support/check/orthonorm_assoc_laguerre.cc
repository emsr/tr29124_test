
#include <cmath>

#include <ext/integration.h>

/* Orthonormality only works for integer alpha. */

// Try to manage the gamma ratio.
template<typename _Tp>
  _Tp
  gamma_ratio(int n, _Tp alpha)
  {
    auto gaman1 = std::tgamma(_Tp(1) + alpha);
    auto fact = gaman1;
    for (int k = 1; k <= n; ++k)
      fact *= (_Tp(k) + alpha) / _Tp(k);
    return fact;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_assoc_laguerre(int n1, int n2, _Tp alpha, _Tp x)
  {
    auto norm = gamma_ratio(n1, alpha);
    return std::pow(x, alpha) * std::exp(-x)
	 * std::assoc_laguerre(n1, alpha, x)
	 * std::assoc_laguerre(n2, alpha, x) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_assoc_laguerre(_Tp alpha)
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = _Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int fail = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2, alpha](_Tp x)
			-> _Tp
			{ return norm_assoc_laguerre<_Tp>(n1, n2, alpha, x); };

	    auto [result, error]
		= __gnu_cxx::integrate_exp_sinh(func, _Tp{0},
						abs_prec, rel_prec);
	    if (std::abs(delta<_Tp>(n1, n2) - result) > cmp_prec)
	      ++fail;
	  }
      }
    return fail;
  }

int
main()
{
  test_assoc_laguerre<float>(0.0F);
  test_assoc_laguerre<double>(0.0);
  test_assoc_laguerre<long double>(0.0L);
}
