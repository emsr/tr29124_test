
/**
 *
 */

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#include <emsr/continued_fractions.h>
#include <emsr/fp_type_util.h>
#include <emsr/math_constants.h>
#include <emsr/complex_util.h> // is_complex
#include <emsr/special_functions.h>

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename Tp, typename _Tzeta>
    std::complex<emsr::num_traits_t<
		emsr::fp_promote_t<_Tnu, Tp, _Tzeta>>>
    cyl_hankel_ratio_c_frac(_Tnu nu, Tp z, _Tzeta zeta)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp, _Tzeta>;
      using _Real = emsr::num_traits_t<_Val>;

      auto c
	= [nu](std::size_t k)
	  -> _Tnu
	  {
	    if (k % 2 == 0)
	      return _Real(2 * k - 3) - _Real{2} * nu;
	    else
	      return _Real(2 * k + 1) + _Real{2} * nu;
	  };

      auto a_H1
	= [nu, zeta, c](std::size_t k, Tp)
	  {
	    using type = decltype(_Tnu{} * _Tzeta{});
	    if (k == 1)
	      return type{1};
	    else
	      return c(k) * zeta;
	  };
      using _NumFun = decltype(a_H1);

      auto b_H1 = [](std::size_t, Tp) -> _Real { return _Real{1}; };
      using _DenFun = decltype(b_H1);

      auto w_H1
	= [nu, zeta, c](std::size_t k, Tp)
	  {
	    return (_Real{-1}
		    + std::sqrt(_Real{1} + _Real{8} * c(k) * zeta))
		 / _Real{2};
	  };
      using _TailFun = decltype(w_H1);

      emsr::SteedContinuedFraction<Tp, _NumFun, _DenFun, _TailFun>
      _H1(a_H1, b_H1, w_H1);

      // b_0 is 0 not 1 so subtract 1.
      return _H1(z) - _Real{1};
    }

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename Tp>
    std::complex<emsr::num_traits_t<
		emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_1_ratio_c_frac(_Tnu nu, Tp z)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto zeta = _Real{-1} / (_Cmplx{0, 2} * z);
      return -cyl_hankel_ratio_c_frac(nu, z, zeta);
    }

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename Tp>
    std::complex<emsr::num_traits_t<
		emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_2_ratio_c_frac(_Tnu nu, Tp z)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto zeta = _Real{1} / (_Cmplx{0, 2} * z);
      return cyl_hankel_ratio_c_frac(nu, z, zeta);
    }

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename Tp>
    std::complex<emsr::num_traits_t<
		emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_bessel_k_ratio_c_frac(_Tnu nu, Tp z)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;

      const auto zeta = _Real{1} / (_Real{2} * z);
      if constexpr (!emsr::is_complex_v<Tp>)
	return std::real(cyl_hankel_ratio_c_frac(nu, z, zeta));
      else
	return cyl_hankel_ratio_c_frac(nu, z, zeta);
    }


  /**
   * Compute ratios of Hankel functions using the J-fraction.
   */
  template<typename _Tnu, typename Tp>
    std::complex<emsr::num_traits_t<
		 emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_ratio_j_frac(_Tnu nu, Tp z, Tp sgn)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto zeta = _Cmplx{0, 2} * z;
      using _Tzeta = decltype(zeta);

      auto a_H
	= [nu](std::size_t k, Tp)
	  {
	    const auto kk = _Tnu(2 * k - 1) / _Tnu{2};
	    return (nu - kk) * (nu + kk);
	  };
      using _NumFun = decltype(a_H);

      auto b_H
	= [zeta, sgn](std::size_t k, Tp)
	  { return sgn * _Tzeta(2 * k) + zeta; };
      using _DenFun = decltype(b_H);

      auto w_H
	= [zeta](std::size_t k, Tp)
	  { return _Tzeta(k) + zeta / _Tzeta{2}; };
      using _TailFun = decltype(w_H);

      emsr::SteedContinuedFraction<Tp, _NumFun, _DenFun, _TailFun>
      _H(a_H, b_H, w_H);

      return (_Tzeta(2 * nu + 1) + sgn * zeta) / (Tp{2} * z)
	   + sgn * (_H(z) - b_H(0, Tp{})) / z;
    }

  /**
   * Return the Hankel function ratio of the first kind from the J-fraction.
   */
  template<typename _Tnu, typename Tp>
    inline std::complex<emsr::num_traits_t<
			emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_1_ratio_j_frac(_Tnu nu, Tp z)
    { return cyl_hankel_ratio_j_frac(nu, z, Tp{-1}); }

  /**
   * Return the Hankel function ratio of the second kind from the J-fraction.
   */
  template<typename _Tnu, typename Tp>
    inline std::complex<emsr::num_traits_t<
			emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_2_ratio_j_frac(_Tnu nu, Tp z)
    { return cyl_hankel_ratio_j_frac(nu, z, Tp{+1}); }

  /**
   * Return the modified Bessel function ratio of the second kind
   * from the J-fraction ratios of Hankel functions.
   */
  template<typename _Tnu, typename Tp,
	   typename _Val = emsr::fp_promote_t<_Tnu, Tp>>
    std::conditional_t<emsr::is_complex_v<_Val>,
			std::complex<emsr::num_traits_t<_Val>>,
			_Val>
    cyl_bessel_k_ratio_j_frac(_Tnu nu, Tp z)
    {
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto s_i = _Cmplx{0, 1};
      const auto s_pi = emsr::pi_v<_Real>;
      const auto ph = std::arg(z);

      _Cmplx _Krat;
      if (ph > -s_pi && ph <= s_pi / _Real{2})
	_Krat = s_i * cyl_hankel_1_ratio_j_frac(nu, s_i * z);
      else
	_Krat = -s_i * cyl_hankel_2_ratio_j_frac(nu, -s_i * z);

      if constexpr (!emsr::is_complex_v<_Val>)
	return std::real(_Krat);
      else
	return _Krat;
    }


template<typename Tp>
  std::complex<Tp>
  cyl_bessel_j_cf2(Tp nu, Tp x)
  {
    const auto s_i = std::complex<Tp>{0, 1};
    const auto s_fp_min = emsr::sqrt_min(nu);
    const auto s_eps = emsr::epsilon(x);
    constexpr int s_max_iter = 15000;

    //const int n = std::max(0, static_cast<int>(nu - x + Tp{1.5L}));
    //const auto mu = nu - Tp(n);
    const auto mu = nu;
    const auto mu2 = mu * mu;
    const auto xi = Tp{1} / x;
    auto a = Tp{0.25L} - mu2;
    auto pq = std::complex<Tp>(-xi / Tp{2}, Tp{1});
    auto b = std::complex<Tp>(Tp{2} * x, Tp{2});
    auto fact = a * xi / std::norm(pq);
    auto c = b + s_i * fact * std::conj(pq);
    auto d = std::conj(b) / std::norm(b);
    auto dl = c * d;
    pq *= dl;
    int i;
    for (i = 2; i <= s_max_iter; ++i)
      {
	a += Tp(2 * (i - 1));
	b += s_i * Tp{2};
	d = a * d + b;
	if (std::abs(d) < s_fp_min)
	  d = s_fp_min;
	fact = a / std::norm(c);
	c = b + fact * std::conj(c);
	if (std::abs(c) < s_fp_min)
	  c = s_fp_min;
	d = std::conj(d) / std::norm(d);
	dl = c * d;
	pq *= dl;
	if (std::abs(dl - Tp{1}) < s_eps)
	  break;
      }
    return pq;
  }


template<typename Tp>
  Tp
  cyl_bessel_k_cf2(Tp nu, Tp x)
  {
    const auto s_eps = emsr::epsilon(x);
    constexpr int s_max_iter = 15000;

    const auto mu = nu;
    const auto mu2 = mu * mu;
    auto b = Tp{2} * (Tp{1} + x);
    auto d = Tp{1} / b;
    auto delh = d;
    auto h = delh;
    auto q1 = Tp{0};
    auto q2 = Tp{1};
    const auto a1 = Tp{0.25L} - mu2;
    auto c = a1;
    auto q = c;
    auto a = -a1;
    auto s = Tp{1} + q * delh;
    int i;
    for (i = 2; i <= s_max_iter; ++i)
      {
	a -= Tp{2 * (i - 1)};
	c = -a * c / i;
	const auto qnew = (q1 - b * q2) / a;
	q1 = q2;
	q2 = qnew;
	q += c * qnew;
	b += Tp{2};
	d = Tp{1} / (b + a * d);
	delh = (b * d - Tp{1}) * delh;
	h += delh;
	const auto dels = q * delh;
	s += dels;
	if (std::abs(dels / s) < s_eps)
	  break;
      }
    if (i > s_max_iter)
      throw std::runtime_error("cyl_bessel_k_cf2: Steed's method failed");
    h *= a1;

    return h;
  }

template<typename Tp>
  void
  test_cyl_hankel_ratio()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto wr = 6 + std::cout.precision();
    auto wc = 4 + 2 * wr;
    using Ret = decltype(cyl_hankel_1_ratio_j_frac(Tp{1}, Tp{1}));

    std::vector<Tp> nu_vec{Tp{0}, Tp{1}/Tp{3}, Tp{1}/Tp{2}, Tp{2}/Tp{3},
			    Tp{1}, Tp{2}, Tp{5}, Tp{10}, Tp{20}, Tp{50}, Tp{100},
			    Tp{128},
			    Tp{200}, Tp{500}, Tp{1000}};

    std::cout << "\n\nRatio H^{(1)}_{\\nu+1}(z) / H^{(1)}_{\\nu}(z)\n";
    std::cout << ' ' << std::setw(wr) << "z"
	      << ' ' << std::setw(wc) << "ratio"
	      << ' ' << std::setw(wc) << "from hankel&deriv"
	      << ' ' << std::setw(wc) << "from hankels"
	      << ' ' << std::setw(wc) << "old_cf2"
	      << ' ' << std::setw(wc) << "nu / z - old_cf2"
	      << ' ' << std::setw(wr) << "delta_r / r"
	      << '\n';
    for (auto nu : nu_vec)
      {
	std::cout << "\n nu = " << std::setw(wr) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = Tp(i) / Tp{10};
	    Ret r{};
	    try
	      {
		r = cyl_hankel_1_ratio_j_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(wr) << z
			  << ' ' << std::setw(wr) << "FAIL\n";
		continue;
	      }
	    const auto h1h2 = emsr::detail::cyl_hankel_h1h2(nu, z);
	    const auto s1 = nu / z - h1h2.H1_deriv / h1h2.H1_value;
	    const auto h1nup1 = emsr::cyl_hankel_1(nu + 1, z);
	    const auto h1nu = emsr::cyl_hankel_1(nu, z);
	    const auto s2 = h1nup1 / h1nu;
	    const auto cf2 = cyl_bessel_j_cf2(nu, z);
	    const auto cf2x = nu / z - cf2;
	    const auto test = std::abs((r - cf2x) / cf2x);
	    std::cout << ' ' << std::setw(wr) << z
		      << ' ' << std::setw(wc) << r
		      << ' ' << std::setw(wc) << s1
		      << ' ' << std::setw(wc) << s2
		      << ' ' << std::setw(wc) << cf2
		      << ' ' << std::setw(wc) << cf2x
		      << ' ' << std::setw(wr) << test
		      << '\n';
	  }
      }

    std::cout << "\n\nRatio H^{(2)}_{\\nu+1}(z) / H^{(2)}_{\\nu}(z)\n";
    std::cout << ' ' << std::setw(wr) << "z"
	      << ' ' << std::setw(wc) << "ratio"
	      << ' ' << std::setw(wc) << "from hankel&deriv"
	      << ' ' << std::setw(wc) << "from hankels"
	      << ' ' << std::setw(wc) << "old_cf2"
	      << ' ' << std::setw(wc) << "nu / z - old_cf2*"
	      << ' ' << std::setw(wr) << "delta_r / r"
	      << '\n';
    for (auto nu : nu_vec)
      {
	std::cout << "\n nu = " << std::setw(wr) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = Tp(i) / Tp{10};
	    Ret r{};
	    try
	      {
		r = cyl_hankel_2_ratio_j_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(wr) << z
			  << ' ' << std::setw(wr) << "FAIL\n";
		continue;
	      }
	    const auto h1h2 = emsr::detail::cyl_hankel_h1h2(nu, z);
	    const auto s1 = nu / z - h1h2.H2_deriv / h1h2.H2_value;
	    const auto h2nup1 = emsr::cyl_hankel_2(nu + 1, z);
	    const auto h2nu = emsr::cyl_hankel_2(nu, z);
	    const auto s2 = h2nup1 / h2nu;
	    const auto cf2 = cyl_bessel_j_cf2(nu, z);
	    const auto cf2x = nu / z - std::conj(cf2);
	    const auto test = std::abs((r - cf2x) / cf2x);
	    std::cout << ' ' << std::setw(wr) << z
		      << ' ' << std::setw(wc) << r
		      << ' ' << std::setw(wc) << s1
		      << ' ' << std::setw(wc) << s2
		      << ' ' << std::setw(wc) << cf2
		      << ' ' << std::setw(wc) << cf2x
		      << ' ' << std::setw(wr) << test
		      << '\n';
	  }
      }

    std::cout << "\n\nRatio K_{\\nu+1}(x) / K_{\\nu}(x)\n";
    std::cout << ' ' << std::setw(wr) << "z"
	      << ' ' << std::setw(wc) << "calculated ratio"
	      << ' ' << std::setw(wc) << "bessel function"
//	      << ' ' << std::setw(wr) << "old_cf2"
	      << ' ' << std::setw(wr) << "delta_r / r"
	      << '\n';
    for (auto nu : nu_vec)
      {
	std::cout << "\n nu = " << std::setw(wr) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = Tp(i) / Tp{10};
	    Ret r{};
	    try
	      {
		r = cyl_bessel_k_ratio_j_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(wr) << z
			  << ' ' << std::setw(wr) << "FAIL\n";
		continue;
	      }
	    const auto ik = emsr::detail::cyl_bessel_ik(nu, z);
	    const auto s = nu / z - ik.K_deriv / ik.K_value;
//	    const auto h = cyl_bessel_k_cf2(nu, z);
	    std::cout << ' ' << std::setw(wr) << z
		      << ' ' << std::setw(wc) << r
		      << ' ' << std::setw(wc) << s
//		      << ' ' << std::setw(wr) << h
		      << ' ' << std::setw(wr) << std::abs((r - s) / s)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_cyl_hankel_ratio<double>();

  return 0;
}
