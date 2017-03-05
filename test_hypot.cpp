/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hypot test_hypot.cpp -L$HOME/bin/lib64 -lquadmath
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

//#define __cpp_lib_hypot 201603L

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
    __hypot3(_Tp __x, _Tp __y, _Tp __z)
    {
      if (__isnan(__x) || __isnan(__y) || __isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto __abs_max = [](_Tp __a, _Tp __b)
			    -> bool
			    { return std::abs(__a) < std::abs(__b); };

	  const auto __amax = std::max({__x, __y, __z}, __abs_max);
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

  // Avoid including all of <algorithm>
  template<typename _Tp>
    constexpr _Tp
    __max3(_Tp __x, _Tp __y, _Tp __z)
    {
      return std::max(std::max(__x, __y), std::max(__y, __z));
    }

  /**
   * Return the three-dimensional hypoteneuse of @c x, @c y, @c z
   * @f[
   *    hypot(x,y,z) = \sqrt{x^2 + y^2 + z^2}
   * @f]
   * avoiding underflow/overflow with small/large arguments.
   */
  template<typename _Tp>
    constexpr _Tp
    __hypot3(_Tp __x, _Tp __y, _Tp __z)
    {
      if (__isnan(__x) || __isnan(__y) || __isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  __x = std::abs(__x);
	  __y = std::abs(__y);
	  __z = std::abs(__z);
	  const auto __amax = __max3(__x, __y, __z);
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

  /**
   * Return the two-dimensional hypoteneuse of complex arguments @c x, @c y:
   * @f[
   *    hypot(x,y) = \sqrt{|x|^2 + |y|^2}
   * @f]
   * avoiding underflow/overflow with small/large arguments.
   */
  template<typename _Tp>
    constexpr __gnu_cxx::__promote_fp_t<_Tp>
    __hypot3(std::complex<_Tp>& __x, std::complex<_Tp>& __y)
    {
      if (__isnan(__x) || __isnan(__y))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto __ax = std::abs(__x);
	  auto __ay = std::abs(__y);
	  const auto __amax = std::max<_Tp>(__ax, __ay);
	  if (__amax == _Tp{0})
	    return _Tp{0};
	  else if (__isinf(__amax))
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    {
	      __ax /= __amax;
	      __ay /= __amax;
	      return __amax * std::sqrt(__ax * __ax + __ay * __ay);
            }
	}
    }

  /**
   * Return the three-dimensional hypoteneuse of complex arguments
   * @c x, @c y, @c z:
   * @f[
   *    hypot(x,y) = \sqrt{|x|^2 + |y|^2 + |z|^2}
   * @f]
   * avoiding underflow/overflow with small/large arguments.
   */
  template<typename _Tp>
    constexpr __gnu_cxx::__promote_fp_t<_Tp>
    __hypot3(std::complex<_Tp>& __x, std::complex<_Tp>& __y,
	     std::complex<_Tp>& __z)
    {
      if (__isnan(__x) || __isnan(__y))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto __ax = std::abs(__x);
	  auto __ay = std::abs(__y);
	  auto __az = std::abs(__z);
	  const auto __amax = std::max<_Tp>({__ax, __ay, __az});
	  if (__amax == _Tp{0})
	    return _Tp{0};
	  else if (__isinf(__amax))
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    {
	      __ax /= __amax;
	      __ay /= __amax;
	      __az /= __amax;
	      return __amax
		   * std::sqrt(__ax * __ax + __ay * __ay + __az * __az);
            }
	}
    }

} // namespace __detail
} // namespace std

  constexpr inline float
  hypot(float __x, float __y, float __z)
  { return std::__detail::__hypot3<float>(__x, __y, __z); }

  constexpr inline double
  hypot(double __x, double __y, double __z)
  { return std::__detail::__hypot3<double>(__x, __y, __z); }

  constexpr inline long double
  hypot(long double __x, long double __y, long double __z)
  { return std::__detail::__hypot3<long double>(__x, __y, __z); }

  template<typename _Tpx, typename _Tpy, typename _Tpz>
    constexpr inline __gnu_cxx::__promote_fp_t<_Tpx, _Tpy, _Tpz>
    hypot(_Tpx __x, _Tpy __y, _Tpz __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpx, _Tpy, _Tpz>;
      return std::__detail::__hypot3<__type>(__x, __y, __z);
    }


  constexpr inline float
  hypot(std::complex<float> __x, std::complex<float> __y)
  { return std::__detail::__hypot3<float>(__x, __y); }

  constexpr inline double
  hypot(std::complex<double> __x, std::complex<double> __y)
  { return std::__detail::__hypot3<double>(__x, __y); }

  constexpr inline long double
  hypot(std::complex<long double> __x, std::complex<long double> __y)
  { return std::__detail::__hypot3<long double>(__x, __y); }

  template<typename _Tpx, typename _Tpy>
    constexpr inline __gnu_cxx::__promote_fp_t<_Tpx, _Tpy>
    hypot(std::complex<_Tpx> __x, std::complex<_Tpy> __y)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpx, _Tpy>;
      return std::__detail::__hypot3<__type>(__x, __y);
    }

  constexpr inline float
  hypot(std::complex<float> __x, std::complex<float> __y,
	std::complex<float> __z)
  { return std::__detail::__hypot3<float>(__x, __y, __z); }

  constexpr inline double
  hypot(std::complex<double> __x, std::complex<double> __y,
	std::complex<double> __z)
  { return std::__detail::__hypot3<double>(__x, __y, __z); }

  constexpr inline long double
  hypot(std::complex<long double> __x, std::complex<long double> __y,
	std::complex<long double> __z)
  { return std::__detail::__hypot3<long double>(__x, __y, __z); }

  template<typename _Tpx, typename _Tpy, typename _Tpz>
    constexpr inline __gnu_cxx::__promote_fp_t<_Tpx, _Tpy, _Tpz>
    hypot(std::complex<_Tpx> __x, std::complex<_Tpy> __y,
	  std::complex<_Tpz> __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpx, _Tpy, _Tpz>;
      return std::__detail::__hypot3<__type>(__x, __y, __z);
    }


//} // namespace std

template<typename _Tp>
  void
  test()
  {
    constexpr auto tiny = _Tp{8} * std::numeric_limits<_Tp>::lowest();
    constexpr auto huge = _Tp{0.125L} * std::numeric_limits<_Tp>::max();
    constexpr auto inf = std::numeric_limits<_Tp>::infinity();
    constexpr auto bigexp = std::numeric_limits<_Tp>::max_exponent - 2;
    constexpr auto big1 = std::ldexp(_Tp{1}, bigexp);
    constexpr auto big2 = std::ldexp(_Tp{2}, bigexp);
    constexpr auto big3 = std::ldexp(_Tp{3}, bigexp);
    constexpr auto smallexp = std::numeric_limits<_Tp>::min_exponent + 2;
    constexpr auto small1 = std::ldexp(_Tp{1}, smallexp);
    constexpr auto small2 = std::ldexp(_Tp{2}, smallexp);
    constexpr auto small3 = std::ldexp(_Tp{3}, smallexp);

    constexpr auto m123 = hypot(_Tp{1}, _Tp{2}, _Tp{3});
    constexpr auto m1inf = hypot(_Tp{1}, _Tp{2}, inf);
    constexpr auto m2inf = hypot(inf, _Tp{2}, inf);
    constexpr auto m3inf = hypot(inf, -inf, inf);
    constexpr auto m1huge = hypot(huge, _Tp{2}, _Tp{3});
    constexpr auto m2huge = hypot(huge, _Tp{2} * huge, _Tp{3});
    constexpr auto m3huge = hypot(huge, _Tp{2} * huge, _Tp{3} * huge);
    constexpr auto m1tiny = hypot(tiny, _Tp{2}, _Tp{3});
    constexpr auto m2tiny = hypot(tiny, _Tp{2} * tiny, _Tp{3});
    constexpr auto m3tiny = hypot(tiny, _Tp{2} * tiny, _Tp{3} * tiny);
    if constexpr (std::numeric_limits<_Tp>::has_denorm == std::denorm_present)
    {
      constexpr auto denorm = std::numeric_limits<_Tp>::denorm_min();
      constexpr auto denorm1 = _Tp{1} * denorm;
      constexpr auto denorm2 = _Tp{2} * denorm;
      constexpr auto denorm3 = _Tp{3} * denorm;
    }
  }

template<typename _Tp>
  int
  test_hypot()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 8 + std::cout.precision();

    constexpr auto inf = std::numeric_limits<_Tp>::infinity();
    constexpr auto bigexp = std::numeric_limits<_Tp>::max_exponent - 2;
    const auto big1 = std::ldexp(_Tp{1}, bigexp);
    const auto big2 = std::ldexp(_Tp{2}, bigexp);
    const auto big3 = std::ldexp(_Tp{3}, bigexp);
    constexpr auto smallexp = std::numeric_limits<_Tp>::min_exponent + 2;
    const auto small1 = std::ldexp(_Tp{1}, smallexp);
    const auto small2 = std::ldexp(_Tp{2}, smallexp);
    const auto small3 = std::ldexp(_Tp{3}, smallexp);

    std::cout << '\n';

    auto m123 = hypot(_Tp{1}, _Tp{2}, _Tp{3});
    std::cout << "m123      = " << std::setw(w) << m123 << '\n';

    auto m1inf = hypot(_Tp{1}, _Tp{2}, inf);
    std::cout << "m1inf     = " << std::setw(w) << m1inf << '\n';

    auto m2inf = hypot(inf, _Tp{2}, inf);
    std::cout << "m2inf     = " << std::setw(w) << m2inf << '\n';

    auto m3inf = hypot(inf, -inf, inf);
    std::cout << "m3inf     = " << std::setw(w) << m3inf << '\n';

    auto m1big = hypot(big1, _Tp{2}, _Tp{3});
    std::cout << "m1big     = " << std::setw(w) << m1big << '\n';

    auto m2big = hypot(big1, big2, _Tp{3});
    std::cout << "m2big     = " << std::setw(w) << m2big << '\n';

    auto m3big = hypot(big1, big2, big3);
    std::cout << "m3big     = " << std::setw(w) << m3big << '\n';
    std::cout << "m3big*    = " << std::setw(w) << m3big / big1 / m123 << '\n';

    auto m1small = hypot(small1, _Tp{2}, _Tp{3});
    std::cout << "m1small   = " << std::setw(w) << m1small << '\n';

    auto m2small = hypot(small1, small2, _Tp{3});
    std::cout << "m2small   = " << std::setw(w) << m2small << '\n';

    auto m3small = hypot(small1, small2, small3);
    std::cout << "m3small   = " << std::setw(w) << m3small << '\n';
    std::cout << "m3small*  = " << std::setw(w) << m3small / small1 / m123 << '\n';

    auto m3zero = hypot(_Tp{0}, _Tp{0}, _Tp{0});
    std::cout << "m3zero    = " << std::setw(w) << m3zero << '\n';

    if constexpr (std::numeric_limits<_Tp>::has_denorm == std::denorm_present)
    {
      constexpr auto denorm = std::numeric_limits<_Tp>::denorm_min();
      constexpr auto denorm1 = _Tp{1} * denorm;
      constexpr auto denorm2 = _Tp{2} * denorm;
      constexpr auto denorm3 = _Tp{3} * denorm;
      auto m3denorm = hypot(denorm1, denorm2, denorm3);
      std::cout << "m3denorm  = " << std::setw(w) << m3denorm << '\n';
      std::cout << "m3denorm* = " << std::setw(w) << m3denorm / denorm / m123 << '\n';
    }

    // x^2 + (x+1)^2 + [x(x+1)]^2 = [x(x+1) + 1]^2
    auto m236 = hypot(_Tp{2}, _Tp{3}, _Tp{6});
    std::cout << "m236      = " << std::setw(w) << m236 << '\n';

    return 0;
  }

template<typename _Tp>
  int
  test_hypot_complex()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 8 + std::cout.precision();

    constexpr auto inf = std::numeric_limits<_Tp>::infinity();
    constexpr auto bigexp = std::numeric_limits<_Tp>::max_exponent - 2;
    const auto big1 = std::ldexp(_Tp{1}, bigexp);
    const auto big2 = std::ldexp(_Tp{2}, bigexp);
    const auto big3 = std::ldexp(_Tp{3}, bigexp);
    constexpr auto smallexp = std::numeric_limits<_Tp>::min_exponent + 2;
    const auto small1 = std::ldexp(_Tp{1}, smallexp);
    const auto small2 = std::ldexp(_Tp{2}, smallexp);
    const auto small3 = std::ldexp(_Tp{3}, smallexp);

    std::complex<_Tp> q(1, -2);
    std::complex<_Tp> x(-3, 4);
    std::complex<_Tp> y(-5, 6);
    std::complex<_Tp> z(-8, 8);

    auto hxy = hypot(x, y);
    std::cout << "hxy   = " << std::setw(w) << hxy << '\n';
    auto hyz = hypot(y, z);
    std::cout << "hyz   = " << std::setw(w) << hyz << '\n';
    auto hzx = hypot(z, x);
    std::cout << "hzx   = " << std::setw(w) << hzx << '\n';

    auto hxyz = hypot(x, y, z);
    std::cout << "hxyz   = " << std::setw(w) << hxyz << '\n';
    auto hqxy = hypot(q, x, y);
    std::cout << "hwxy   = " << std::setw(w) << hqxy << '\n';
    auto hqyz = hypot(q, y, z);
    std::cout << "hwyz   = " << std::setw(w) << hqyz << '\n';
    auto hqzx = hypot(q, z, x);
    std::cout << "hqzx   = " << std::setw(w) << hqzx << '\n';

    using _Cmplx = std::complex<_Tp>;
    auto m1rinf = hypot(_Cmplx(_Tp{1}), _Cmplx(_Tp{2}), _Cmplx(inf));
    std::cout << "m1rinf    = " << std::setw(w) << m1rinf << '\n';
    auto m1iinf = hypot(_Cmplx(_Tp{0}, _Tp{1}), _Cmplx(_Tp{0}, _Tp{2}), _Cmplx(_Tp{0}, inf));
    std::cout << "m1iinf    = " << std::setw(w) << m1iinf << '\n';

    auto m2rinf = hypot(_Cmplx(inf), _Cmplx(_Tp{2}), _Cmplx(inf));
    std::cout << "m2rinf    = " << std::setw(w) << m2rinf << '\n';
    auto m2iinf = hypot(_Cmplx(_Tp{0}, inf), _Cmplx(_Tp{0}, _Tp{2}), _Cmplx(_Tp{0}, inf));
    std::cout << "m2iinf    = " << std::setw(w) << m2iinf << '\n';

    auto m3rinf = hypot(_Cmplx(inf), _Cmplx(-inf), _Cmplx(inf));
    std::cout << "m3rinf    = " << std::setw(w) << m3rinf << '\n';
    auto m3iinf = hypot(_Cmplx(_Tp{0}, inf), _Cmplx(_Tp{0}, -inf), _Cmplx(_Tp{0}, inf));
    std::cout << "m3iinf    = " << std::setw(w) << m3iinf << '\n';

    auto m1rbig = hypot(_Cmplx(big1), _Cmplx(_Tp{2}), _Cmplx(_Tp{3}));
    std::cout << "m1rbig    = " << std::setw(w) << m1rbig << '\n';
    auto m1ibig = hypot(_Cmplx(_Tp{0}, big1), _Cmplx(_Tp{0}, _Tp{2}), _Cmplx(_Tp{0}, _Tp{3}));
    std::cout << "m1ibig    = " << std::setw(w) << m1ibig << '\n';
    auto m1abig = hypot(_Cmplx(big1, big1), _Cmplx(_Tp{2}, _Tp{2}), _Cmplx(_Tp{3}, _Tp{3}));
    std::cout << "m1abig    = " << std::setw(w) << m1abig << '\n';

    auto m2rbig = hypot(_Cmplx(big1), _Cmplx(big2), _Cmplx(_Tp{3}));
    std::cout << "m2rbig    = " << std::setw(w) << m2rbig << '\n';
    auto m2ibig = hypot(_Cmplx(_Tp{0}, big1), _Cmplx(_Tp{0}, big2), _Cmplx(_Tp{0}, _Tp{3}));
    std::cout << "m2ibig    = " << std::setw(w) << m2ibig << '\n';
    auto m2abig = hypot(_Cmplx(big1, big1), _Cmplx(big2, big2), _Cmplx(_Tp{3}, _Tp{3}));
    std::cout << "m2abig    = " << std::setw(w) << m2abig << '\n';

    auto m3rbig = hypot(_Cmplx(big1), _Cmplx(big2), _Cmplx(big3));
    std::cout << "m3rbig    = " << std::setw(w) << m3rbig << '\n';
    auto m3ibig = hypot(_Cmplx(_Tp{0}, big1), _Cmplx(_Tp{0}, big2), _Cmplx(_Tp{0}, big3));
    std::cout << "m3ibig    = " << std::setw(w) << m3ibig << '\n';
    auto m3abig = hypot(_Cmplx(big1, big1), _Cmplx(big2, big2), _Cmplx(big3, big3));
    std::cout << "m3abig    = " << std::setw(w) << m3abig << '\n';

    auto m1rsmall = hypot(_Cmplx(small1), _Cmplx(_Tp{2}), _Cmplx(_Tp{3}));
    std::cout << "m1rsmall  = " << std::setw(w) << m1rsmall << '\n';
    auto m1ismall = hypot(_Cmplx(_Tp{0}, small1), _Cmplx(_Tp{0}, _Tp{2}), _Cmplx(_Tp{0}, _Tp{3}));
    std::cout << "m1ismall  = " << std::setw(w) << m1ismall << '\n';
    auto m1asmall = hypot(_Cmplx(small1, small1), _Cmplx(_Tp{2}, _Tp{2}), _Cmplx(_Tp{3}, _Tp{3}));
    std::cout << "m1asmall  = " << std::setw(w) << m1asmall << '\n';

    auto m2rsmall = hypot(_Cmplx(small1), _Cmplx(small2), _Cmplx(_Tp{3}));
    std::cout << "m2rsmall  = " << std::setw(w) << m2rsmall << '\n';
    auto m2ismall = hypot(_Cmplx(_Tp{0}, small1), _Cmplx(_Tp{0}, small2), _Cmplx(_Tp{0}, _Tp{3}));
    std::cout << "m2ismall  = " << std::setw(w) << m2ismall << '\n';
    auto m2asmall = hypot(_Cmplx(small1, small1), _Cmplx(small2, small2), _Cmplx(_Tp{3}, _Tp{3}));
    std::cout << "m2asmall  = " << std::setw(w) << m2asmall << '\n';

    auto m3rsmall = hypot(_Cmplx(small1), _Cmplx(small2), _Cmplx(small3));
    std::cout << "m3rsmall  = " << std::setw(w) << m3rsmall << '\n';
    auto m3ismall = hypot(_Cmplx(_Tp{0}, small1), _Cmplx(_Tp{0}, small2), _Cmplx(_Tp{0}, small3));
    std::cout << "m3ismall  = " << std::setw(w) << m3ismall << '\n';
    auto m3asmall = hypot(_Cmplx(small1, small1), _Cmplx(small2, small2), _Cmplx(small3, small3));
    std::cout << "m3asmall  = " << std::setw(w) << m3asmall << '\n';

    return 0;
  }

int
main()
{
  test_hypot<float>();
  test_hypot<double>();
  test_hypot<long double>();
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  test_hypot<__float128>();
#endif

  //test<float>();
  //test<double>();
  //test<long double>();

  test_hypot_complex<double>();
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
