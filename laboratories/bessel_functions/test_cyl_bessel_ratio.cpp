/**
 *
 */

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#include <emsr/continued_fractions.h>
#include <emsr/fp_type_util.h>
#include <emsr/complex_util.h> // is_complex
#include <emsr/special_functions.h>

  /**
   * Compute ratios of Bessel functions using the S-fraction.
   *
   * @param  nu  The order of the Hankel ratio.
   * @param  x   The argument of the Hankel ratio.
   * @param  zeta A variable encapsulating the regular and irregular
   *           Bessel functions; Use @f$ \zeta = (iz)^2 @f$ for @f$ J_\nu(z) @f$
   *           and @f$ \zeta = z^2 @f$ for @f$ I_\nu(z) @f$.
   */
  template<typename _Tnu, typename Tp, typename _Tzeta>
    std::complex<emsr::num_traits_t<
		 emsr::fp_promote_t<_Tnu, Tp, _Tzeta>>>
    cyl_bessel_ratio_s_frac(_Tnu nu, Tp z, _Tzeta zeta)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;

      auto a_J
	= [nu, z, zeta](std::size_t k, Tp)
	  {
	    using type = decltype(_Tnu{} * _Tzeta{});
	    if (k == 1)
	      return type(z / (_Tnu{2} * nu + _Tnu{2}));
	    else
	      return zeta
		   / (_Tnu{4} * (nu + _Tnu(k - 1)) * (nu + _Tnu(k)));
	  };
      using _AFun = decltype(a_J);

      auto b_J = [](std::size_t, Tp) -> _Real { return _Real{1}; };
      using _BFun = decltype(b_J);

      auto w_J = [](std::size_t, Tp) -> _Real { return _Real{0}; };
      using _WFun = decltype(w_J);

      //emsr::SteedContinuedFraction<Tp, _AFun, _BFun, _WFun>
      emsr::LentzContinuedFraction<Tp, _AFun, _BFun, _WFun>
      _J(a_J, b_J, w_J);

      // b_0 is 0 not 1 so subtract 1.
      return _J(z) - _Real{1};
    }

  /**
   * 
   */
  template<typename _Tnu, typename Tp,
	   typename _Val = emsr::fp_promote_t<_Tnu, Tp>>
    std::conditional_t<emsr::is_complex_v<_Val>,
			std::complex<emsr::num_traits_t<_Val>>,
			_Val>
    cyl_bessel_j_ratio_s_frac(_Tnu nu, Tp z)
    {
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto iz = _Cmplx{0, 1} * z;
      const auto zeta = iz * iz;
      const auto _Jrat = cyl_bessel_ratio_s_frac(nu, z, zeta);

      if constexpr (!emsr::is_complex_v<_Val>)
	return std::real(_Jrat);
      else
	return _Jrat;
    }

  /**
   * 
   */
  template<typename _Tnu, typename Tp,
	   typename _Val = emsr::fp_promote_t<_Tnu, Tp>>
    std::conditional_t<emsr::is_complex_v<_Val>,
			std::complex<emsr::num_traits_t<_Val>>,
			_Val>
    cyl_bessel_i_ratio_s_frac(_Tnu nu, Tp z)
    {
      const auto zeta = z * z;
      const auto _Irat = cyl_bessel_ratio_s_frac(nu, z, zeta);

      if constexpr (!emsr::is_complex_v<_Val>)
	return std::real(_Irat);
      else
	return _Irat;
    }

template<typename Tp>
  Tp
  cyl_bessel_j_cf1(Tp nu, Tp x)
  {
    const auto s_fp_min = emsr::sqrt_min(nu);
    const auto s_eps = emsr::epsilon(x);
    constexpr int s_max_iter = 15000;

    int isign = 1;
    const auto xi = Tp{1} / x;
    const auto xi2 = Tp{2} * xi;
    auto h = std::max(s_fp_min, nu * xi);
    auto b = xi2 * nu;
    auto d = Tp{0};
    auto c = h;
    int i;
    for (i = 1; i <= s_max_iter; ++i)
      {
	b += xi2;
	d = b - d;
	if (std::abs(d) < s_fp_min)
	  d = s_fp_min;
	d = Tp{1} / d;
	c = b - Tp{1} / c;
	if (std::abs(c) < s_fp_min)
	  c = s_fp_min;
	const auto del = c * d;
	h *= del;
	if (d < Tp{0})
	  isign = -isign;
	if (std::abs(del - Tp{1}) < s_eps)
	  break;
      }
    return h;
  }

template<typename Tp>
  Tp
  cyl_bessel_i_cf1(Tp nu, Tp x)
  {
    const auto s_fp_min = emsr::sqrt_min(nu);
    const auto s_eps = emsr::epsilon(x);
    constexpr int s_max_iter = 15000;

    const auto xi = Tp{1} / x;
    const auto xi2 = Tp{2} * xi;
    auto h = std::max(s_fp_min, nu * xi);
    auto b = xi2 * nu;
    auto d = Tp{0};
    auto c = h;
    int i;
    for (i = 1; i <= s_max_iter; ++i)
      {
	b += xi2;
	d = Tp{1} / (b + d);
	c = b + Tp{1} / c;
	const auto del = c * d;
	h *= del;
	if (std::abs(del - Tp{1}) < s_eps)
	  break;
      }
    return h;
  }

template<typename Tp>
  void
  test_cyl_bessel_ratio()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 6 + std::cout.precision();
    using Ret = decltype(cyl_bessel_j_ratio_s_frac(Tp{1}, Tp{1}));

    // What's with the NaN at 2?
    cyl_bessel_j_ratio_s_frac(Tp{1}, Tp{2});

    std::cout << "\n\nRatio J_{\\nu+1}(x) / J_{\\nu}(x)\n";
    std::cout << ' ' << std::setw(w) << "z"
	      << ' ' << std::setw(w) << "new ratio"
	      << ' ' << std::setw(w) << "lab w/deriv ratio"
	      << ' ' << std::setw(w) << "nu / z - cf1"
	      << '\n';
    for (auto nu : {Tp{0}, Tp{1}/Tp{3}, Tp{1}/Tp{2}, Tp{2}/Tp{3}, Tp{1}, Tp{2}, Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = Tp(i) / Tp{10};
	    Ret r{};
	    try
	      {
		r = cyl_bessel_j_ratio_s_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    const auto cf1 = cyl_bessel_j_cf1(nu, z);
	    const auto jn = emsr::detail::cyl_bessel_jn(nu, z);
	    const auto s = nu / z - jn.J_deriv / jn.J_value;
	    std::cout << ' ' << std::setw(w) << z
		      << ' ' << std::setw(w) << r
		      << ' ' << std::setw(w) << s
		      << ' ' << std::setw(w) << nu / z - cf1
		      << '\n';
	  }
      }

    std::cout << "\n\nRatio I_{\\nu+1}(x) / I_{\\nu}(x)\n";
    std::cout << ' ' << std::setw(w) << "z"
	      << ' ' << std::setw(w) << "new ratio"
	      << ' ' << std::setw(w) << "lab w/deriv ratio"
	      << ' ' << std::setw(w) << "cf1 - nu / z"
	      << '\n';
    for (auto nu : {Tp{0}, Tp{1}/Tp{3}, Tp{1}/Tp{2}, Tp{2}/Tp{3}, Tp{1}, Tp{2}, Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = Tp(i) / Tp{10};
	    Ret r{};
	    try
	      {
		r = cyl_bessel_i_ratio_s_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    const auto cf1 = cyl_bessel_i_cf1(nu, z);
	    const auto ik = emsr::detail::cyl_bessel_ik(nu, z);
	    const auto s = -nu / z + ik.I_deriv / ik.I_value;
	    std::cout << ' ' << std::setw(w) << z
		      << ' ' << std::setw(w) << r
		      << ' ' << std::setw(w) << s
		      << ' ' << std::setw(w) << cf1 - nu / z
		      << '\n';
	  }
      }
  }

int
main()
{
  test_cyl_bessel_ratio<double>();

  return 0;
}
