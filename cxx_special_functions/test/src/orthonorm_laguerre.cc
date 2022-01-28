
#include <emsr/integration.h>
#include <emsr/sf_laguerre.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_laguerre(int n1, int n2, _Tp x)
  {
    return std::exp(-x)
	 * emsr::laguerre(n1, x)
	 * emsr::laguerre(n2, x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_laguerre()
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = _Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2](_Tp x)
			-> _Tp
			{ return norm_laguerre<_Tp>(n1, n2, x); };

	    auto [result, error]
		= emsr::integrate_exp_sinh(func, _Tp{0}, abs_prec, rel_prec);

	    if (std::abs(delta<_Tp>(n1, n2) - result) > cmp_prec)
	      ++num_errors;
	  }
      }
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  num_errors += test_laguerre<float>();
  num_errors += test_laguerre<double>();
  num_errors += test_laguerre<long double>();
  return num_errors;
}
