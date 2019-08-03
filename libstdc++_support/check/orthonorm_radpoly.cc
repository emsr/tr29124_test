#include <cmath>

#include <ext/integration.h>

// Neumann's number
template<typename _Tp>
  _Tp
  epsilon(int m)
  { return m == 0 ? _Tp{2} : _Tp{1}; }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_radpoly(int n1, int m1, int n2, int m2, _Tp rho)
  {
    auto norm = _Tp{1} / std::sqrt(_Tp(2 * n1 + 2) * _Tp(2 * n2 + 2));
    return rho
	 * __gnu_cxx::radpoly(n1, m1, rho)
	 * __gnu_cxx::radpoly(n2, m2, rho) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_radpoly(int m1, int m2)
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
	if ((n1 - m1) & 1)
	  continue;
	for (const auto n2 : degree)
	  {
	    //int m2 = m1;
	    if (m2 > n2)
	      continue;
	    if ((n2 - m2) & 1)
	      continue;

	    auto func = [n1, m1, n2, m2](_Tp x)
			-> _Tp
			{ return norm_radpoly(n1, m1, n2, m2, x); };

	    auto [result, error]
		= __gnu_cxx::integrate_singular(func, _Tp{0}, _Tp{1}, abs_prec, rel_prec);

	    if (std::abs(delta<_Tp>(n1, m1, n2, m2) - result) > cmp_prec)
	      ++fail;
	  }
      }
    return fail;
  }

int
main()
{
  test_radpoly<float>(1, 2);
  test_radpoly<double>(1, 2);
  test_radpoly<long double>(1, 2);
}
