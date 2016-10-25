/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_hypot test_hypot.cpp -L$HOME/bin/lib64 -lquadmath
./test_hypot > test_hypot.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_hypot test_hypot.cpp -L$HOME/bin/lib64 -lquadmath
./test_hypot > test_hypot.txt

g++ -std=gnu++14 -Wall -Wextra -DNO_LOGBQ -I. -o test_hypot test_hypot.cpp -lquadmath
./test_hypot > test_hypot.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <bits/specfun_util.h>

#define __cpp_lib_hypot 201603L

namespace std
{
namespace __detail
{
#if LAMBDA
  /**
   * Return the three-dimensional hypoteneuse of @c x, @c y, @c z
   * @f[
   *    hypot(x,y,z) = \sqrt{x^2 + y^2 + z^2}
   * @f]
   * avoiding underflow/overflow with small/large arguments.
   */
  template<typename _Tp>
    constexpr _Tp
    __hypot(_Tp __x, _Tp __y, _Tp __z)
    {
      if (__isnan(__x) || __isnan(__y) || __isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto __abs_max = [](_Tp __a, _Tp __b)
			    -> bool
			    { return std::abs(__a) < std::abs(__b); };

	  auto __amax = std::max({__x, __y, __z}, __abs_max);
	  if (__amax == _Tp{0})
	    return _Tp{0};
	  else if (__isinf(__amax))
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    {
	      __x /= __amax;
	      __y /= __amax;
	      __z /= __amax;
	      return __amax * std::sqrt(__x * __x + __y * __y + __z * __z);
            }
	}
    }
#else
  /**
   * Return the three-dimensional hypoteneuse of @c x, @c y, @c z
   * @f[
   *    hypot(x,y,z) = \sqrt{x^2 + y^2 + z^2}
   * @f]
   * avoiding underflow/overflow with small/large arguments.
   */
  template<typename _Tp>
    constexpr _Tp
    __hypot(_Tp __x, _Tp __y, _Tp __z)
    {
      if (__isnan(__x) || __isnan(__y) || __isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  __x = std::abs(__x);
	  __y = std::abs(__y);
	  __z = std::abs(__z);
	  auto __amax = std::max({__x, __y, __z});
	  if (__amax == _Tp{0})
	    return _Tp{0};
	  else if (__isinf(__amax))
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    {
	      __x /= __amax;
	      __y /= __amax;
	      __z /= __amax;
	      return __amax * std::sqrt(__x * __x + __y * __y + __z * __z);
            }
	}
    }
#endif
} // namespace __detail

  float
  hypot(float __x, float __y, float __z)
  { return std::__detail::__hypot<float>(__x, __y, __z); }

  double
  hypot(double __x, double __y, double __z)
  { return std::__detail::__hypot<double>(__x, __y, __z); }

  long double
  hypot(long double __x, long double __y, long double __z)
  { return std::__detail::__hypot<long double>(__x, __y, __z); }

  template<typename _Tpx, typename _Tpy, typename _Tpz>
    inline __gnu_cxx::__promote_fp_t<_Tpx, _Tpy, _Tpz>
    hypot(_Tpx __x, _Tpy __y, _Tpz __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpx, _Tpy, _Tpz>;
      return std::__detail::__hypot<__type>(__x, __y, __z);
    }

} // namespace std

template<typename _Tp>
  void
  test()
  {
    constexpr auto tiny = _Tp{8} * std::numeric_limits<_Tp>::lowest();
    constexpr auto huge = _Tp{0.125L} * std::numeric_limits<_Tp>::max();

    constexpr auto m123 = std::hypot(_Tp{1}, _Tp{2}, _Tp{3});
    constexpr auto m1huge = std::hypot(huge, _Tp{2}, _Tp{3});
    constexpr auto m2huge = std::hypot(huge, _Tp{2} * huge, _Tp{3});
    constexpr auto m3huge = std::hypot(huge, _Tp{2} * huge, _Tp{3} * huge);
    constexpr auto m1tiny = std::hypot(tiny, _Tp{2}, _Tp{3});
    constexpr auto m2tiny = std::hypot(tiny, _Tp{2} * tiny, _Tp{3});
    constexpr auto m3tiny = std::hypot(tiny, _Tp{2} * tiny, _Tp{3} * tiny);
  }

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  auto w = 8 + std::cout.precision();
  auto inf = std::numeric_limits<double>::infinity();

  auto m123 = std::hypot(1.0, 2.0, 3.0);
  std::cout << "m123    = " << std::setw(w) << m123 << '\n';

  auto m1inf = std::hypot(1.0, 2.0, inf);
  std::cout << "m1inf   = " << std::setw(w) << m1inf << '\n';

  auto m2inf = std::hypot(inf, 2.0, inf);
  std::cout << "m2inf   = " << std::setw(w) << m2inf << '\n';

  auto m3inf = std::hypot(inf, -inf, inf);
  std::cout << "m3inf   = " << std::setw(w) << m3inf << '\n';

  auto m1big = std::hypot(1.0e300, 2.0, 3.0);
  std::cout << "m1big   = " << std::setw(w) << m1big << '\n';

  auto m2big = std::hypot(1.0e300, 2.0e300, 3.0);
  std::cout << "m2big   = " << std::setw(w) << m2big << '\n';

  auto m3big = std::hypot(1.0e300, 2.0e300, 3.0e300);
  std::cout << "m3big   = " << std::setw(w) << m3big << '\n';

  auto m1small = std::hypot(1.0e-300, 2.0, 3.0);
  std::cout << "m1small = " << std::setw(w) << m1small << '\n';

  auto m2small = std::hypot(1.0e-300, 2.0e-300, 3.0);
  std::cout << "m2small = " << std::setw(w) << m2small << '\n';

  auto m3small = std::hypot(1.0e-300, 2.0e-300, 3.0e-300);
  std::cout << "m3small = " << std::setw(w) << m3small << '\n';

  auto m3zero = std::hypot(0.0, 0.0, 0.0);
  std::cout << "m3zero  = " << std::setw(w) << m3zero << '\n';

  // x^2 + (x+1)^2 + [x(x+1)]^2 = [x(x+1) + 1]^2
  auto m236 = std::hypot(2.0, 3.0, 6.0);
  std::cout << "m236    = " << std::setw(w) << m236 << '\n';
}

/*
m123    =        3.74165738677394
m1inf   =                     inf
m2inf   =                     inf
m3inf   =                     inf
m1big   =                  1e+300
m2big   =   2.23606797749979e+300
m3big   =   3.74165738677394e+300
m1small =        3.60555127546399
m2small =                       3
m3small =   3.74165738677394e-300
m3zero  =                       0
m236    =                       7
*/
