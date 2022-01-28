
#include <emsr/integration.h>
#include <emsr/sf_chebyshev.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_chebyshev_v(int n1, int n2, _Tp x)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_inf = std::numeric_limits<_Tp>::infinity();
    if (std::abs(x - _Tp{1}) < _S_eps)
      return (n1 + n2) & 1 ? -_S_inf : _S_inf;
    else
      return emsr::chebyshev_v(n1, x)
	   * emsr::chebyshev_v(n2, x)
	   * std::sqrt((_Tp{1} + x) / (_Tp{1} - x))
	   / _Tp{2};
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_chebyshev_v()
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
			{return norm_chebyshev_v(n1, n2, x);};

	    auto [result, error]
		= emsr::integrate_singular_endpoints(func,
				 _Tp{-1}, _Tp{1},
				 _Tp{0.5}, _Tp{-0.5}, 0, 0,
				 abs_prec, rel_prec);

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
  num_errors += test_chebyshev_v<float>();
  num_errors += test_chebyshev_v<double>();
  num_errors += test_chebyshev_v<long double>();
  return num_errors;
}
