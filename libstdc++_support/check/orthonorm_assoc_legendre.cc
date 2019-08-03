#include <cmath>

#include <ext/integration.h>

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_assoc_legendre(int l1, int m1, int l2, int m2, _Tp x)
  {
    return (_Tp(l1 + l2 + 1) / _Tp{2})
	 * std::assoc_legendre(l1, m1, x)
	 * std::assoc_legendre(l2, m2, x);
  }

template<typename _Tp>
  _Tp
  delta(int l1, int l2)
  { return l1 == l2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_assoc_legendre(int m1, int m2)
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = _Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int fail = 0;

    for (const auto l1 : degree)
      {
	for (const auto l2 : degree)
	  {
	    auto func = [l1, m1, l2, m2](_Tp x)
			-> _Tp
			{ return normalized_assoc_legendre(l1, m1, l2, m2, x); };

	    auto [result, error]
		= __gnu_cxx::integrate_tanh_sinh(func, _Tp{-1}, _Tp{1},
						 abs_prec, rel_prec, 6);

	    if (std::abs(delta<_Tp>(l1, l2) - result) > cmp_prec)
	      ++fail;
	  }
      }
    return fail;
  }

int
main()
{
  test_assoc_legendre<float>(0, 0);
  test_assoc_legendre<double>(0, 0);
  test_assoc_legendre<long double>(0, 0);
}
