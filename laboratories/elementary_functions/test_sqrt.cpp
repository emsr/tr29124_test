/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <emsr/math_constants.h> // For math constants.

/**
 * 
 */
template<typename Tp>
  constexpr Tp
  sqrt_newton(Tp x)
  {
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    const auto s_digs = std::numeric_limits<Tp>::digits;
    constexpr auto rsqrt2 = Tp{1} / emsr::sqrt2_v<Tp>;

    if (x == Tp{0})
      return Tp{0};

    int e;
    auto xx = std::frexp(x, &e);

    auto f = 0.4 + 0.6 * xx;

    for (int i = 0; i < s_digs; ++i)
      {
	f -= Tp{0.5L} * (f - xx / f);
	if (std::abs(f * f - xx) < s_eps * f)
	  break;
      }
    if (e & 1)
      return rsqrt2 * std::ldexp(f, (e + 1) >> 1);
    else
      return std::ldexp(f, e >> 1);
  }

/**
 * 
 */
template<typename Tp>
  constexpr Tp
  sqrt_recip_newton(Tp x)
  {
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    const auto s_digs = std::numeric_limits<Tp>::digits;

    if (x == Tp{0})
      return Tp{0};

    // TODO: Argument reduction...

    auto f = 0.01;
    auto ff = f * f;
    for (int i = 0; i < s_digs; ++i)
      {
	f *= (Tp{3} - x * ff) / Tp{2};
	ff = f * f;
	if (std::abs(ff - x) < s_eps * f)
	  break;
      }

    return f;
  }

/**
 * 
 */
template<typename Tp>
  void
  test_sqrt()
  {
    std::cout.precision(std::numeric_limits<Tp>::max_digits10);
    const auto w = 6 + std::cout.precision();

    for (int i = 0; i <= 1000; ++i)
      {
	auto x = Tp{1.0L} * i;
	auto sqrtr = x * sqrt_recip_newton(x);
	auto sqrtcw = sqrt_newton(x);
	auto sqrt = std::sqrt(x);
	auto deltar = sqrtr - sqrt;
	auto deltacw = sqrtcw - sqrt;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << sqrtr
		  << ' ' << std::setw(w) << sqrtcw
		  << ' ' << std::setw(w) << sqrt
		  << ' ' << std::setw(w) << deltar
		  << ' ' << std::setw(w) << deltacw
		  << '\n';
      }
  }

int
main()
{
  test_sqrt<double>();
}
