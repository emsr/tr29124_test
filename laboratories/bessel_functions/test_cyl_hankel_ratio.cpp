
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
      using _Cmplx = std::complex<_Real>;

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
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = _Real{1} / (_Real{2} * __z);
      if constexpr (!__gnu_cxx::is_complex_v<_Tp>)
	return std::real(__cyl_hankel_ratio_c_frac(__nu, __z, __zeta));
      else
	return __cyl_hankel_ratio_c_frac(__nu, __z, __zeta);
    }


  /**
   * Compute rations of Hankel functions using the J-fraction.
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
   * from the J-fraction rations of Hankel functions.
   */
  template<typename _Tnu, typename _Tp>
    inline std::complex<std::__detail::__num_traits_t<
			__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_bessel_k_ratio_j_frac(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_i = _Cmplx{0, 1};
      const auto _S_pi = __gnu_cxx::math::__pi_v<_Real>;
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
  void
  test_cyl_hankel_ratio()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto wr = 6 + std::cout.precision();
    auto wc = 4 + 2 * wr;
    using Ret = decltype(__cyl_hankel_1_ratio_j_frac(_Tp{1}, _Tp{1}));

    std::cout << "\n\nRatio H^{(1)}_{\\nu+1}(z) / H^{(1)}_{\\nu}(z)\n";
    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
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
	    std::cout << ' ' << std::setw(wr) << z
		      << ' ' << std::setw(wc) << r
		      << ' ' << std::setw(wc) << s1
		      << ' ' << std::setw(wc) << s2
		      << '\n';
	  }
      }

    std::cout << "\n\nRatio H^{(2)}_{\\nu+1}(z) / H^{(2)}_{\\nu}(z)\n";
    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
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
	    std::cout << ' ' << std::setw(wr) << z
		      << ' ' << std::setw(wc) << r
		      << ' ' << std::setw(wc) << s1
		      << ' ' << std::setw(wc) << s2
		      << '\n';
	  }
      }

    std::cout << "\n\nRatio K_{\\nu+1}(x) / K_{\\nu}(x)\n";
    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
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
	    std::cout << ' ' << std::setw(wr) << z
		      << ' ' << std::setw(wc) << r
		      << ' ' << std::setw(wc) << s
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
