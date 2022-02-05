/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/special_functions.h>

  /**
   * A struct to store a cosine and a sine value.
   */
  template<typename Tp>
    struct sincos_t
    {
      Tp sin_v;
      Tp cos_v;
    };

  /**
   * Default implementation of sincos.
   */
  template<typename Tp>
    inline sincos_t<Tp>
    sincos(Tp x)
    { return sincos_t<Tp>{std::sin(x), std::cos(x)}; }

  template<>
    inline sincos_t<float>
    sincos(float x)
    {
      float sin, cos;
      __builtin_sincosf(x, &sin, &cos);
      return sincos_t<float>{sin, cos};
    }

  template<>
    inline sincos_t<double>
    sincos(double x)
    {
      double sin, cos;
      __builtin_sincos(x, &sin, &cos);
      return sincos_t<double>{sin, cos};
    }

  template<>
    inline sincos_t<long double>
    sincos(long double x)
    {
      long double sin, cos;
      __builtin_sincosl(x, &sin, &cos);
      return sincos_t<long double>{sin, cos};
    }

#ifdef EMSR_HAVE_FLOAT128
/*
  template<>
    inline sincos_t<__float128>
    sincos(__float128 x)
    {
      __float128 sin, cos;
      ::sincosq(x, &sin, &cos);
      return sincos_t<__float128>{sin, cos};
    }
*/
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

  /**
   * Reperiodized sincos.
   */
  template<typename Tp>
    sincos_t<Tp>
    sincos_pi(Tp x)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_NaN = emsr::quiet_NaN<Tp>(x);
      if (std::isnan(x))
	return sincos_t<Tp>{s_NaN, s_NaN};
      else if (x < Tp{0})
	{
	  sincos_t<Tp> tempsc = sincos_pi(-x);
	  return sincos_t<Tp>{-tempsc.sin_v,
					     tempsc.cos_v};
	}
      else if (x < Tp{0.5L})
	return sincos(s_pi * x);
      else if (x < Tp{1})
	{
	  sincos_t<Tp>
	    tempsc = sincos(s_pi * (Tp{1} - x));
	  return sincos_t<Tp>{tempsc.sin_v,
					   -tempsc.cos_v};
	}
      else
	{
	  auto nu = std::floor(x);
	  auto arg = x - nu;
	  auto sign = (int(nu) & 1) == 1 ? Tp{-1} : Tp{+1};

	  auto sinval = (arg < Tp{0.5L})
			? std::sin(s_pi * arg)
			: std::sin(s_pi * (Tp{1} - arg));
	  auto cosval = std::cos(s_pi * arg);
	  return sincos_t<Tp>{sign * sinval,
					    sign * cosval};
	}
    }

  /**
   * Reperiodized complex constructor.
   */
  template<typename Tp>
    inline std::complex<Tp>
    polar_pi(Tp rho, Tp phi_pi)
    {
      sincos_t<Tp> sc = sincos_pi(phi_pi);
      return std::complex<Tp>(rho * sc.cos_v, rho * sc.sin_v);
    }


template<typename Tp>
  void
  test_sincos(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto pi = emsr::pi_v<Tp>;

    std::cout << '\n';
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sincos.sin"
	      << std::setw(width) << "sincos_pi.sin"
	      << std::setw(width) << "delta sin"
	      << std::setw(width) << "delta sin"
	      << std::setw(width) << "sincos.cos"
	      << std::setw(width) << "sincos_pi.cos"
	      << std::setw(width) << "delta cos"
	      << std::setw(width) << "delta cos"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto del = Tp{1} / Tp{10};
    for (int i = -40; i <= +40; ++i)
      {
	auto x = del * i;
	auto sc = emsr::detail::sincos(pi * x);
	auto scpi = emsr::detail::sincos_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sc.sin_v
		  << std::setw(width) << scpi.sin_v
		  << std::setw(width) << sc.sin_v - scpi.sin_v
		  << std::setw(width) << scpi.sin_v - std::sin(pi * x)
		  << std::setw(width) << sc.cos_v
		  << std::setw(width) << scpi.cos_v
		  << std::setw(width) << sc.cos_v - scpi.cos_v
		  << std::setw(width) << scpi.cos_v - std::cos(pi * x)
		  << '\n';
      }
  }

int
main()
{
  constexpr auto pif = emsr::pi_v<float>;
  constexpr auto pi = emsr::pi_v<double>;
  constexpr auto pil = emsr::pi_v<long double>;
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //constexpr auto piq = emsr::pi_v<__float128>;
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

  auto a1 [[maybe_unused]] = /*emsr::detail::*/sincos(pif * 1.5f);
  auto a2 [[maybe_unused]] = /*emsr::detail::*/sincos_pi(1.5f);

  auto b1 [[maybe_unused]] = /*emsr::detail::*/sincos(pi * 1.5);
  auto b2 [[maybe_unused]] = /*emsr::detail::*/sincos_pi(1.5);

  auto c1 [[maybe_unused]] = /*emsr::detail::*/sincos(pil * 1.5l);
  auto c2 [[maybe_unused]] = /*emsr::detail::*/sincos_pi(1.5l);

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //auto d1 [[maybe_unused]] = /*emsr::detail::*/sincos(piq * 1.5q);
  //auto d2 [[maybe_unused]] = /*emsr::detail::*/sincos_pi(1.5q);
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

  test_sincos<double>();
}
