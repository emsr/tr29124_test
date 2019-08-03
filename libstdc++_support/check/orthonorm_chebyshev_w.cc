#include <cmath>

#include <ext/integration.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_chebyshev_w(int n1, int n2, _Tp x)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_inf = std::numeric_limits<_Tp>::infinity();
    const auto
    _S_pi = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    if (std::abs(x + _Tp{1}) < _S_eps)
      return (n1 + n2) & 1 ? -_S_inf : _S_inf;
    else
      return __gnu_cxx::chebyshev_w(n1, x)
	   * __gnu_cxx::chebyshev_w(n2, x)
	   * std::sqrt((_Tp{1} - x) / (_Tp{1} + x))
	   / _S_pi;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_chebyshev_w()
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
	    auto func = [n1, n2](_Tp x)
			-> _Tp
			{return norm_chebyshev_w(n1, n2, x);};

	    auto [result, error]
		= __gnu_cxx::integrate_singular_endpoints(func,
				 _Tp{-1}, _Tp{1},
				 _Tp{-0.5}, _Tp{0.5}, 0, 0,
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
  test_chebyshev_w<float>();
  test_chebyshev_w<double>();
  test_chebyshev_w<long double>();
}
