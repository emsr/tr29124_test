/*
$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -o test_charlier test_charlier.cpp -Lwrappers/debug -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=$HOME/bin/lib64:wrappers/debug:$LD_LIBRARY_PATH ./test_charlier > test_charlier.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <wrap_burkhardt.h>

/**
 * Compute the Charlier polynomial by recursion:
 * @f[
 *    -x C_n(x) = a C_{n+1}(x) - (n + a) C_n(x) + n C_{n-1}(x)
 * @f]
 * where @f$ C_n(x) = C_n(x; a) @f$.
 */
template<typename _Tp, typename _TpX>
  _Tp
  __charlier_recur(int n, _Tp a, _TpX x)
  {
    auto Cnm1 = _Tp{1};
    if (n == 0)
      return Cnm1;

    auto Cn = _Tp{1} - _Tp(x) / a;
    if (n == 1)
      return Cn;

    auto Cnp1 = ((_Tp{1} + a - _Tp(x)) * Cn - Cnm1) / a;
    for (int k = 2; k < n; ++k)
      {
	Cnm1 = Cn;
	Cn = Cnp1;
	Cnp1 = ((_Tp(k) + a - _Tp(x)) * Cn - _Tp(k) * Cnm1) / a;
      }

    return Cnp1;
  }

/**
 * Return the Charlier polynomial defined by
 * @f[
 *    C_n(x; a) = {}_2F_0(-n, -x; ; -\frac{1}{a})
 * @f]
 */
template<typename _Tp, typename _TpX>
  _Tp
  __charlier(int n, _Tp a, _TpX x)
  {
    if (std::isnan(a))
      return a;
    else if (std::isnan(x))
      return x;
    else
      return __charlier_recur(n, a, x);
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_charlier(int n_max, _Tp a)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto C = __charlier(n, a, x);
	    auto C_test = burkhardt::charlier(n, a, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << C
		      << ' ' << std::setw(w) << C_test
		      << '\n';
	  }
      }
  }

int
main()
{
  test_charlier<float>(10, 2.0f);
}
