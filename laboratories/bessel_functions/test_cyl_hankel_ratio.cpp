
/**
 *
 */

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#include <ext/continued_fractions.h>
#include <bits/specfun_util.h>
#include <ext/math_constants.h>
#include <ext/complex_util.h> // is_complex

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename _Tp, typename _Tzeta>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp, _Tzeta>>>
    __cyl_hankel_ratio_c_frac(_Tnu __nu, _Tp __z, _Tzeta __zeta)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp, _Tzeta>;
      using _Real = std::__detail::__num_traits_t<_Val>;

      auto __c
	= [__nu](std::size_t __k)
	  -> _Tnu
	  {
	    if (__k % 2 == 0)
	      return _Real(2 * __k - 3) - _Real{2} * __nu;
	    else
	      return _Real(2 * __k + 1) + _Real{2} * __nu;
	  };

      auto __a_H1
	= [__nu, __zeta, __c](std::size_t __k, _Tp)
	  {
	    using __type = decltype(_Tnu{} * _Tzeta{});
	    if (__k == 1)
	      return __type{1};
	    else
	      return __c(__k) * __zeta;
	  };
      using _NumFun = decltype(__a_H1);

      auto __b_H1 = [](std::size_t, _Tp) -> _Real { return _Real{1}; };
      using _DenFun = decltype(__b_H1);

      auto __w_H1
	= [__nu, __zeta, __c](std::size_t __k, _Tp)
	  {
	    return (_Real{-1}
		    + std::sqrt(_Real{1} + _Real{8} * __c(__k) * __zeta))
		 / _Real{2};
	  };
      using _TailFun = decltype(__w_H1);

      _SteedContinuedFraction<_Tp, _NumFun, _DenFun, _TailFun>
      _H1(__a_H1, __b_H1, __w_H1);

      // b_0 is 0 not 1 so subtract 1.
      return _H1(__z) - _Real{1};
    }

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_1_ratio_c_frac(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = _Real{-1} / (_Cmplx{0, 2} * __z);
      return -__cyl_hankel_ratio_c_frac(__nu, __z, __zeta);
    }

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_2_ratio_c_frac(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = _Real{1} / (_Cmplx{0, 2} * __z);
      return __cyl_hankel_ratio_c_frac(__nu, __z, __zeta);
    }

  /**
   * This C-fraction is almost useless as AFAICT.
   */
  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_bessel_k_ratio_c_frac(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;

      const auto __zeta = _Real{1} / (_Real{2} * __z);
      if constexpr (!__gnu_cxx::is_complex_v<_Tp>)
	return std::real(__cyl_hankel_ratio_c_frac(__nu, __z, __zeta));
      else
	return __cyl_hankel_ratio_c_frac(__nu, __z, __zeta);
    }


  /**
   * Compute ratios of Hankel functions using the J-fraction.
   */
  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		 __gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_ratio_j_frac(_Tnu __nu, _Tp __z, _Tp __sgn)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto __zeta = _Cmplx{0, 2} * __z;
      using _Tzeta = decltype(__zeta);

      auto __a_H
	= [__nu](std::size_t __k, _Tp)
	  {
	    const auto __kk = _Tnu(2 * __k - 1) / _Tnu{2};
	    return (__nu - __kk) * (__nu + __kk);
	  };
      using _NumFun = decltype(__a_H);

      auto __b_H
	= [__zeta, __sgn](std::size_t __k, _Tp)
	  { return __sgn * _Tzeta(2 * __k) + __zeta; };
      using _DenFun = decltype(__b_H);

      auto __w_H
	= [__zeta](std::size_t __k, _Tp)
	  { return _Tzeta(__k) + __zeta / _Tzeta{2}; };
      using _TailFun = decltype(__w_H);

      _SteedContinuedFraction<_Tp, _NumFun, _DenFun, _TailFun>
      _H(__a_H, __b_H, __w_H);

      return (_Tzeta(2 * __nu + 1) + __sgn * __zeta) / (_Tp{2} * __z)
	   + __sgn * (_H(__z) - __b_H(0, _Tp{})) / __z;
    }

  /**
   * Return the Hankel function ratio of the first kind from the J-fraction.
   */
  template<typename _Tnu, typename _Tp>
    inline std::complex<std::__detail::__num_traits_t<
			__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_1_ratio_j_frac(_Tnu __nu, _Tp __z)
    { return __cyl_hankel_ratio_j_frac(__nu, __z, _Tp{-1}); }

  /**
   * Return the Hankel function ratio of the second kind from the J-fraction.
   */
  template<typename _Tnu, typename _Tp>
    inline std::complex<std::__detail::__num_traits_t<
			__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_2_ratio_j_frac(_Tnu __nu, _Tp __z)
    { return __cyl_hankel_ratio_j_frac(__nu, __z, _Tp{+1}); }

  /**
   * Return the modified Bessel function ratio of the second kind
   * from the J-fraction ratios of Hankel functions.
   */
  template<typename _Tnu, typename _Tp,
	   typename _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>>
    std::conditional_t<__gnu_cxx::is_complex_v<_Val>,
			std::complex<std::__detail::__num_traits_t<_Val>>,
			_Val>
    __cyl_bessel_k_ratio_j_frac(_Tnu __nu, _Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_i = _Cmplx{0, 1};
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Real>;
      const auto __ph = std::arg(__z);

      _Cmplx _Krat;
      if (__ph > -_S_pi && __ph <= _S_pi / _Real{2})
	_Krat = _S_i * __cyl_hankel_1_ratio_j_frac(__nu, _S_i * __z);
      else
	_Krat = -_S_i * __cyl_hankel_2_ratio_j_frac(__nu, -_S_i * __z);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
	return std::real(_Krat);
      else
	return _Krat;
    }


template<typename _Tp>
  std::complex<_Tp>
  cyl_bessel_j_cf2(_Tp __nu, _Tp __x)
  {
    const auto _S_i = std::complex<_Tp>{0, 1};
    const auto _S_fp_min = __gnu_cxx::__sqrt_min(__nu);
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    constexpr int _S_max_iter = 15000;

    //const int __n = std::max(0, static_cast<int>(__nu - __x + _Tp{1.5L}));
    //const auto __mu = __nu - _Tp(__n);
    const auto __mu = __nu;
    const auto __mu2 = __mu * __mu;
    const auto __xi = _Tp{1} / __x;
    auto __a = _Tp{0.25L} - __mu2;
    auto __pq = std::complex<_Tp>(-__xi / _Tp{2}, _Tp{1});
    auto __b = std::complex<_Tp>(_Tp{2} * __x, _Tp{2});
    auto __fact = __a * __xi / std::norm(__pq);
    auto __c = __b + _S_i * __fact * std::conj(__pq);
    auto __d = std::conj(__b) / std::norm(__b);
    auto __dl = __c * __d;
    __pq *= __dl;
    int __i;
    for (__i = 2; __i <= _S_max_iter; ++__i)
      {
	__a += _Tp(2 * (__i - 1));
	__b += _S_i * _Tp{2};
	__d = __a * __d + __b;
	if (std::abs(__d) < _S_fp_min)
	  __d = _S_fp_min;
	__fact = __a / std::norm(__c);
	__c = __b + __fact * std::conj(__c);
	if (std::abs(__c) < _S_fp_min)
	  __c = _S_fp_min;
	__d = std::conj(__d) / std::norm(__d);
	__dl = __c * __d;
	__pq *= __dl;
	if (std::abs(__dl - _Tp{1}) < _S_eps)
	  break;
      }
    return __pq;
  }


template<typename _Tp>
  _Tp
  cyl_bessel_k_cf2(_Tp __nu, _Tp __x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    constexpr int _S_max_iter = 15000;

    const auto __mu = __nu;
    const auto __mu2 = __mu * __mu;
    auto __b = _Tp{2} * (_Tp{1} + __x);
    auto __d = _Tp{1} / __b;
    auto __delh = __d;
    auto __h = __delh;
    auto __q1 = _Tp{0};
    auto __q2 = _Tp{1};
    const auto __a1 = _Tp{0.25L} - __mu2;
    auto __c = __a1;
    auto __q = __c;
    auto __a = -__a1;
    auto __s = _Tp{1} + __q * __delh;
    int __i;
    for (__i = 2; __i <= _S_max_iter; ++__i)
      {
	__a -= _Tp{2 * (__i - 1)};
	__c = -__a * __c / __i;
	const auto __qnew = (__q1 - __b * __q2) / __a;
	__q1 = __q2;
	__q2 = __qnew;
	__q += __c * __qnew;
	__b += _Tp{2};
	__d = _Tp{1} / (__b + __a * __d);
	__delh = (__b * __d - _Tp{1}) * __delh;
	__h += __delh;
	const auto __dels = __q * __delh;
	__s += __dels;
	if (std::abs(__dels / __s) < _S_eps)
	  break;
      }
    if (__i > _S_max_iter)
      std::__throw_runtime_error(__N("__cyl_bessel_ik_steed: "
				     "Steed's method failed"));
    __h *= __a1;

    return __h;
  }

template<typename _Tp>
  void
  test_cyl_hankel_ratio()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto wr = 6 + std::cout.precision();
    auto wc = 4 + 2 * wr;
    using Ret = decltype(__cyl_hankel_1_ratio_j_frac(_Tp{1}, _Tp{1}));

    std::vector<_Tp> nu_vec{_Tp{0}, _Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3},
			    _Tp{1}, _Tp{2}, _Tp{5}, _Tp{10}, _Tp{20}, _Tp{50}, _Tp{100},
			    _Tp{128},
			    _Tp{200}, _Tp{500}, _Tp{1000}};

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
	    auto z = _Tp(i) / _Tp{10};
	    Ret r{};
	    try
	      {
		r = __cyl_hankel_1_ratio_j_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(wr) << z
			  << ' ' << std::setw(wr) << "FAIL\n";
		continue;
	      }
	    const auto h1h2 = std::__detail::__cyl_hankel_h1h2(nu, z);
	    const auto s1 = nu / z - h1h2.__H1_deriv / h1h2.__H1_value;
	    const auto h1nup1 = __gnu_cxx::cyl_hankel_1(nu + 1, z);
	    const auto h1nu = __gnu_cxx::cyl_hankel_1(nu, z);
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
	    auto z = _Tp(i) / _Tp{10};
	    Ret r{};
	    try
	      {
		r = __cyl_hankel_2_ratio_j_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(wr) << z
			  << ' ' << std::setw(wr) << "FAIL\n";
		continue;
	      }
	    const auto h1h2 = std::__detail::__cyl_hankel_h1h2(nu, z);
	    const auto s1 = nu / z - h1h2.__H2_deriv / h1h2.__H2_value;
	    const auto h2nup1 = __gnu_cxx::cyl_hankel_2(nu + 1, z);
	    const auto h2nu = __gnu_cxx::cyl_hankel_2(nu, z);
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
	      << ' ' << std::setw(wr) << "old_cf2"
	      << ' ' << std::setw(wr) << "delta_r / r"
	      << '\n';
    for (auto nu : nu_vec)
      {
	std::cout << "\n nu = " << std::setw(wr) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    Ret r{};
	    try
	      {
		r = __cyl_bessel_k_ratio_j_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(wr) << z
			  << ' ' << std::setw(wr) << "FAIL\n";
		continue;
	      }
	    const auto ik = std::__detail::__cyl_bessel_ik(nu, z);
	    const auto s = nu / z - ik.__K_deriv / ik.__K_value;
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
