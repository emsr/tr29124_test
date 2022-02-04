/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/complex_util.h> // isnan

//#define cpp_lib_hypot 201603L

namespace emsr
{
namespace detail
{
#if LAMBDA

  /**
   * Return the three-dimensional hypoteneuse of @c x, @c y, @c z
   * @f[
   *    hypot(x,y,z) = \sqrt{x^2 + y^2 + z^2}
   * @f]
   * avoiding underflow/overflow with small/large arguments.
   */
  template<typename Tp>
    constexpr Tp
    hypot3(Tp x, Tp y, Tp z)
    {
      if (std::isnan(x) || std::isnan(y) || std::isnan(z))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  auto abs_max = [](Tp a, Tp b)
			    -> bool
			    { return std::abs(a) < std::abs(b); };

	  const auto amax = std::max({x, y, z}, abs_max);
	  if (amax == Tp{0})
	    return Tp{0};
	  else if (std::isinf(amax))
	    return std::numeric_limits<Tp>::infinity();
	  else
	    {
	      x /= amax;
	      y /= amax;
	      z /= amax;
	      return amax * std::sqrt(x * x + y * y + z * z);
            }
	}
    }

#else

  // Avoid including all of <algorithm>
  template<typename Tp>
    constexpr Tp
    max3(Tp x, Tp y, Tp z)
    {
      return std::max(std::max(x, y), std::max(y, z));
    }

  /**
   * Return the three-dimensional hypoteneuse of @c x, @c y, @c z
   * @f[
   *    hypot(x,y,z) = \sqrt{x^2 + y^2 + z^2}
   * @f]
   * avoiding underflow/overflow with small/large arguments.
   */
  template<typename Tp>
    constexpr Tp
    hypot3(Tp x, Tp y, Tp z)
    {
      if (std::isnan(x) || std::isnan(y) || std::isnan(z))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  x = std::abs(x);
	  y = std::abs(y);
	  z = std::abs(z);
	  const auto amax = max3(x, y, z);
	  if (amax == Tp{0})
	    return Tp{0};
	  else if (std::isinf(amax))
	    return std::numeric_limits<Tp>::infinity();
	  else
	    {
	      x /= amax;
	      y /= amax;
	      z /= amax;
	      return amax * std::sqrt(x * x + y * y + z * z);
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
  template<typename Tp>
    constexpr emsr::fp_promote_t<Tp>
    hypot3(const std::complex<Tp>& x, const std::complex<Tp>& y)
    {
      if (std::isnan(x) || std::isnan(y))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  auto ax = std::abs(x);
	  auto ay = std::abs(y);
	  const auto amax = std::max<Tp>(ax, ay);
	  if (amax == Tp{0})
	    return Tp{0};
	  else if (std::isinf(amax))
	    return std::numeric_limits<Tp>::infinity();
	  else
	    {
	      ax /= amax;
	      ay /= amax;
	      return amax * std::sqrt(ax * ax + ay * ay);
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
  template<typename Tp>
    constexpr emsr::fp_promote_t<Tp>
    hypot3(const std::complex<Tp>& x, const std::complex<Tp>& y,
	     const std::complex<Tp>& z)
    {
      if (std::isnan(x) || std::isnan(y))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  auto ax = std::abs(x);
	  auto ay = std::abs(y);
	  auto az = std::abs(z);
	  const auto amax = std::max<Tp>({ax, ay, az});
	  if (amax == Tp{0})
	    return Tp{0};
	  else if (std::isinf(amax))
	    return std::numeric_limits<Tp>::infinity();
	  else
	    {
	      ax /= amax;
	      ay /= amax;
	      az /= amax;
	      return amax
		   * std::sqrt(ax * ax + ay * ay + az * az);
            }
	}
    }

} // namespace detail
} // namespace emsr

  constexpr inline float
  hypot(float x, float y, float z)
  { return emsr::detail::hypot3<float>(x, y, z); }

  constexpr inline double
  hypot(double x, double y, double z)
  { return emsr::detail::hypot3<double>(x, y, z); }

  constexpr inline long double
  hypot(long double x, long double y, long double z)
  { return emsr::detail::hypot3<long double>(x, y, z); }

  template<typename Tpx, typename Tpy, typename Tpz>
    constexpr inline emsr::fp_promote_t<Tpx, Tpy, Tpz>
    hypot(Tpx x, Tpy y, Tpz z)
    {
      using type = emsr::fp_promote_t<Tpx, Tpy, Tpz>;
      return emsr::detail::hypot3<type>(x, y, z);
    }


  constexpr inline float
  hypot(const std::complex<float>& x, const std::complex<float>& y)
  { return emsr::detail::hypot3<float>(x, y); }

  constexpr inline double
  hypot(const std::complex<double>& x, const std::complex<double>& y)
  { return emsr::detail::hypot3<double>(x, y); }

  constexpr inline long double
  hypot(const std::complex<long double>& x,
	const std::complex<long double>& y)
  { return emsr::detail::hypot3<long double>(x, y); }

  template<typename Tpx, typename Tpy>
    constexpr inline emsr::fp_promote_t<Tpx, Tpy>
    hypot(const std::complex<Tpx>& x, const std::complex<Tpy>& y)
    {
      using type = emsr::fp_promote_t<Tpx, Tpy>;
      return emsr::detail::hypot3<type>(x, y);
    }

  constexpr inline float
  hypot(const std::complex<float>& x, const std::complex<float>& y,
	const std::complex<float>& z)
  { return emsr::detail::hypot3<float>(x, y, z); }

  constexpr inline double
  hypot(const std::complex<double>& x, const std::complex<double>& y,
	const std::complex<double>& z)
  { return emsr::detail::hypot3<double>(x, y, z); }

  constexpr inline long double
  hypot(const std::complex<long double>& x,
	const std::complex<long double>& y,
	const std::complex<long double>& z)
  { return emsr::detail::hypot3<long double>(x, y, z); }

  template<typename Tpx, typename Tpy, typename Tpz>
    constexpr inline emsr::fp_promote_t<Tpx, Tpy, Tpz>
    hypot(const std::complex<Tpx>& x, const std::complex<Tpy>& y,
	  const std::complex<Tpz>& z)
    {
      using type = emsr::fp_promote_t<Tpx, Tpy, Tpz>;
      return emsr::detail::hypot3<type>(x, y, z);
    }


//} // namespace emsr

template<typename Tp>
  void
  test()
  {
    constexpr auto tiny = Tp{8} * std::numeric_limits<Tp>::lowest();
    constexpr auto huge = Tp{0.125L} * std::numeric_limits<Tp>::max();
    constexpr auto inf = std::numeric_limits<Tp>::infinity();
    constexpr auto bigexp = std::numeric_limits<Tp>::max_exponent - 2;
    constexpr auto big1 = std::ldexp(Tp{1}, bigexp);
    constexpr auto big2 = std::ldexp(Tp{2}, bigexp);
    constexpr auto big3 = std::ldexp(Tp{3}, bigexp);
    constexpr auto smallexp = std::numeric_limits<Tp>::min_exponent + 2;
    constexpr auto small1 = std::ldexp(Tp{1}, smallexp);
    constexpr auto small2 = std::ldexp(Tp{2}, smallexp);
    constexpr auto small3 = std::ldexp(Tp{3}, smallexp);

    constexpr auto m123 = hypot(Tp{1}, Tp{2}, Tp{3});
    constexpr auto m1inf = hypot(Tp{1}, Tp{2}, inf);
    constexpr auto m2inf = hypot(inf, Tp{2}, inf);
    constexpr auto m3inf = hypot(inf, -inf, inf);
    constexpr auto m1huge = hypot(huge, Tp{2}, Tp{3});
    constexpr auto m2huge = hypot(huge, Tp{2} * huge, Tp{3});
    constexpr auto m3huge = hypot(huge, Tp{2} * huge, Tp{3} * huge);
    constexpr auto m1tiny = hypot(tiny, Tp{2}, Tp{3});
    constexpr auto m2tiny = hypot(tiny, Tp{2} * tiny, Tp{3});
    constexpr auto m3tiny = hypot(tiny, Tp{2} * tiny, Tp{3} * tiny);
    if constexpr (std::numeric_limits<Tp>::has_denorm == std::denorm_present)
    {
      constexpr auto denorm = std::numeric_limits<Tp>::denorm_min();
      constexpr auto denorm1 = Tp{1} * denorm;
      constexpr auto denorm2 = Tp{2} * denorm;
      constexpr auto denorm3 = Tp{3} * denorm;
    }
  }

template<typename Tp>
  int
  test_hypot()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();

    constexpr auto inf = std::numeric_limits<Tp>::infinity();
    constexpr auto bigexp = std::numeric_limits<Tp>::max_exponent - 2;
    const auto big1 = std::ldexp(Tp{1}, bigexp);
    const auto big2 = std::ldexp(Tp{2}, bigexp);
    const auto big3 = std::ldexp(Tp{3}, bigexp);
    constexpr auto smallexp = std::numeric_limits<Tp>::min_exponent + 2;
    const auto small1 = std::ldexp(Tp{1}, smallexp);
    const auto small2 = std::ldexp(Tp{2}, smallexp);
    const auto small3 = std::ldexp(Tp{3}, smallexp);

    std::cout << '\n';

    auto m123 = hypot(Tp{1}, Tp{2}, Tp{3});
    std::cout << "m123      = " << std::setw(w) << m123 << '\n';

    auto m1inf = hypot(Tp{1}, Tp{2}, inf);
    std::cout << "m1inf     = " << std::setw(w) << m1inf << '\n';

    auto m2inf = hypot(inf, Tp{2}, inf);
    std::cout << "m2inf     = " << std::setw(w) << m2inf << '\n';

    auto m3inf = hypot(inf, -inf, inf);
    std::cout << "m3inf     = " << std::setw(w) << m3inf << '\n';

    auto m1big = hypot(big1, Tp{2}, Tp{3});
    std::cout << "m1big     = " << std::setw(w) << m1big << '\n';

    auto m2big = hypot(big1, big2, Tp{3});
    std::cout << "m2big     = " << std::setw(w) << m2big << '\n';

    auto m3big = hypot(big1, big2, big3);
    std::cout << "m3big     = " << std::setw(w) << m3big << '\n';
    std::cout << "m3big*    = " << std::setw(w) << m3big / big1 / m123 << '\n';

    auto m1small = hypot(small1, Tp{2}, Tp{3});
    std::cout << "m1small   = " << std::setw(w) << m1small << '\n';

    auto m2small = hypot(small1, small2, Tp{3});
    std::cout << "m2small   = " << std::setw(w) << m2small << '\n';

    auto m3small = hypot(small1, small2, small3);
    std::cout << "m3small   = " << std::setw(w) << m3small << '\n';
    std::cout << "m3small*  = " << std::setw(w) << m3small / small1 / m123 << '\n';

    auto m3zero = hypot(Tp{0}, Tp{0}, Tp{0});
    std::cout << "m3zero    = " << std::setw(w) << m3zero << '\n';

    if constexpr (std::numeric_limits<Tp>::has_denorm == std::denorm_present)
    {
      constexpr auto denorm = std::numeric_limits<Tp>::denorm_min();
      constexpr auto denorm1 = Tp{1} * denorm;
      constexpr auto denorm2 = Tp{2} * denorm;
      constexpr auto denorm3 = Tp{3} * denorm;
      auto m3denorm = hypot(denorm1, denorm2, denorm3);
      std::cout << "m3denorm  = " << std::setw(w) << m3denorm << '\n';
      std::cout << "m3denorm* = " << std::setw(w) << m3denorm / denorm / m123 << '\n';
    }

    // x^2 + (x+1)^2 + [x(x+1)]^2 = [x(x+1) + 1]^2
    auto m236 = hypot(Tp{2}, Tp{3}, Tp{6});
    std::cout << "m236      = " << std::setw(w) << m236 << '\n';

    return 0;
  }

template<typename Tp>
  int
  test_hypot_complex()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();

    constexpr auto inf = std::numeric_limits<Tp>::infinity();
    constexpr auto bigexp = std::numeric_limits<Tp>::max_exponent - 2;
    const auto big1 = std::ldexp(Tp{1}, bigexp);
    const auto big2 = std::ldexp(Tp{2}, bigexp);
    const auto big3 = std::ldexp(Tp{3}, bigexp);
    constexpr auto smallexp = std::numeric_limits<Tp>::min_exponent + 2;
    const auto small1 = std::ldexp(Tp{1}, smallexp);
    const auto small2 = std::ldexp(Tp{2}, smallexp);
    const auto small3 = std::ldexp(Tp{3}, smallexp);

    std::complex<Tp> q(1, -2);
    std::complex<Tp> x(-3, 4);
    std::complex<Tp> y(-5, 6);
    std::complex<Tp> z(-8, 8);

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

    using Cmplx = std::complex<Tp>;
    auto m1rinf = hypot(Cmplx(Tp{1}), Cmplx(Tp{2}), Cmplx(inf));
    std::cout << "m1rinf    = " << std::setw(w) << m1rinf << '\n';
    auto m1iinf = hypot(Cmplx(Tp{0}, Tp{1}), Cmplx(Tp{0}, Tp{2}), Cmplx(Tp{0}, inf));
    std::cout << "m1iinf    = " << std::setw(w) << m1iinf << '\n';

    auto m2rinf = hypot(Cmplx(inf), Cmplx(Tp{2}), Cmplx(inf));
    std::cout << "m2rinf    = " << std::setw(w) << m2rinf << '\n';
    auto m2iinf = hypot(Cmplx(Tp{0}, inf), Cmplx(Tp{0}, Tp{2}), Cmplx(Tp{0}, inf));
    std::cout << "m2iinf    = " << std::setw(w) << m2iinf << '\n';

    auto m3rinf = hypot(Cmplx(inf), Cmplx(-inf), Cmplx(inf));
    std::cout << "m3rinf    = " << std::setw(w) << m3rinf << '\n';
    auto m3iinf = hypot(Cmplx(Tp{0}, inf), Cmplx(Tp{0}, -inf), Cmplx(Tp{0}, inf));
    std::cout << "m3iinf    = " << std::setw(w) << m3iinf << '\n';

    auto m1rbig = hypot(Cmplx(big1), Cmplx(Tp{2}), Cmplx(Tp{3}));
    std::cout << "m1rbig    = " << std::setw(w) << m1rbig << '\n';
    auto m1ibig = hypot(Cmplx(Tp{0}, big1), Cmplx(Tp{0}, Tp{2}), Cmplx(Tp{0}, Tp{3}));
    std::cout << "m1ibig    = " << std::setw(w) << m1ibig << '\n';
    auto m1abig = hypot(Cmplx(big1, big1), Cmplx(Tp{2}, Tp{2}), Cmplx(Tp{3}, Tp{3}));
    std::cout << "m1abig    = " << std::setw(w) << m1abig << '\n';

    auto m2rbig = hypot(Cmplx(big1), Cmplx(big2), Cmplx(Tp{3}));
    std::cout << "m2rbig    = " << std::setw(w) << m2rbig << '\n';
    auto m2ibig = hypot(Cmplx(Tp{0}, big1), Cmplx(Tp{0}, big2), Cmplx(Tp{0}, Tp{3}));
    std::cout << "m2ibig    = " << std::setw(w) << m2ibig << '\n';
    auto m2abig = hypot(Cmplx(big1, big1), Cmplx(big2, big2), Cmplx(Tp{3}, Tp{3}));
    std::cout << "m2abig    = " << std::setw(w) << m2abig << '\n';

    auto m3rbig = hypot(Cmplx(big1), Cmplx(big2), Cmplx(big3));
    std::cout << "m3rbig    = " << std::setw(w) << m3rbig << '\n';
    auto m3ibig = hypot(Cmplx(Tp{0}, big1), Cmplx(Tp{0}, big2), Cmplx(Tp{0}, big3));
    std::cout << "m3ibig    = " << std::setw(w) << m3ibig << '\n';
    auto m3abig = hypot(Cmplx(big1, big1), Cmplx(big2, big2), Cmplx(big3, big3));
    std::cout << "m3abig    = " << std::setw(w) << m3abig << '\n';

    auto m1rsmall = hypot(Cmplx(small1), Cmplx(Tp{2}), Cmplx(Tp{3}));
    std::cout << "m1rsmall  = " << std::setw(w) << m1rsmall << '\n';
    auto m1ismall = hypot(Cmplx(Tp{0}, small1), Cmplx(Tp{0}, Tp{2}), Cmplx(Tp{0}, Tp{3}));
    std::cout << "m1ismall  = " << std::setw(w) << m1ismall << '\n';
    auto m1asmall = hypot(Cmplx(small1, small1), Cmplx(Tp{2}, Tp{2}), Cmplx(Tp{3}, Tp{3}));
    std::cout << "m1asmall  = " << std::setw(w) << m1asmall << '\n';

    auto m2rsmall = hypot(Cmplx(small1), Cmplx(small2), Cmplx(Tp{3}));
    std::cout << "m2rsmall  = " << std::setw(w) << m2rsmall << '\n';
    auto m2ismall = hypot(Cmplx(Tp{0}, small1), Cmplx(Tp{0}, small2), Cmplx(Tp{0}, Tp{3}));
    std::cout << "m2ismall  = " << std::setw(w) << m2ismall << '\n';
    auto m2asmall = hypot(Cmplx(small1, small1), Cmplx(small2, small2), Cmplx(Tp{3}, Tp{3}));
    std::cout << "m2asmall  = " << std::setw(w) << m2asmall << '\n';

    auto m3rsmall = hypot(Cmplx(small1), Cmplx(small2), Cmplx(small3));
    std::cout << "m3rsmall  = " << std::setw(w) << m3rsmall << '\n';
    auto m3ismall = hypot(Cmplx(Tp{0}, small1), Cmplx(Tp{0}, small2), Cmplx(Tp{0}, small3));
    std::cout << "m3ismall  = " << std::setw(w) << m3ismall << '\n';
    auto m3asmall = hypot(Cmplx(small1, small1), Cmplx(small2, small2), Cmplx(small3, small3));
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
  //test_hypot<__float128>();
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
