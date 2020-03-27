/**
 *
 */

#include <cmath>
#include <complex>

#include <ext/continued_fractions.h>
#include <ext/fp_type_util.h>
#include <ext/complex_util.h> // is_complex

  /**
   * Compute ratios of Bessel functions using the S-fraction.
   *
   * @param  __nu  The order of the Hankel ratio.
   * @param  __x   The argument of the Hankel ratio.
   * @param  __zeta A variable encapsulating the regular and irregular
   *           Bessel functions; Use @f$ \zeta = (iz)^2 @f$ for @f$ J_\nu(z) @f$
   *           and @f$ \zeta = z^2 @f$ for @f$ I_\nu(z) @f$.
   */
  template<typename _Tnu, typename _Tp, typename _Tzeta>
    std::complex<__gnu_cxx::__num_traits_t<
		 __gnu_cxx::fp_promote_t<_Tnu, _Tp, _Tzeta>>>
    __cyl_bessel_ratio_s_frac(_Tnu __nu, _Tp __z, _Tzeta __zeta)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;

      auto __a_J
	= [__nu, __z, __zeta](std::size_t __k, _Tp)
	  {
	    using __type = decltype(_Tnu{} * _Tzeta{});
	    if (__k == 1)
	      return __type(__z / (_Tnu{2} * __nu + _Tnu{2}));
	    else
	      return __zeta
		   / (_Tnu{4} * (__nu + _Tnu(__k - 1)) * (__nu + _Tnu(__k)));
	  };
      using _AFun = decltype(__a_J);

      auto __b_J = [](std::size_t, _Tp) -> _Real { return _Real{1}; };
      using _BFun = decltype(__b_J);

      auto __w_J = [](std::size_t, _Tp) -> _Real { return _Real{0}; };
      using _WFun = decltype(__w_J);

      //_SteedContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      _J(__a_J, __b_J, __w_J);

      // b_0 is 0 not 1 so subtract 1.
      return _J(__z) - _Real{1};
    }

  /**
   * 
   */
  template<typename _Tnu, typename _Tp,
	   typename _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>>
    std::conditional_t<__gnu_cxx::is_complex_v<_Val>,
			std::complex<__gnu_cxx::__num_traits_t<_Val>>,
			_Val>
    __cyl_bessel_j_ratio_s_frac(_Tnu __nu, _Tp __z)
    {
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __iz = _Cmplx{0, 1} * __z;
      const auto __zeta = __iz * __iz;
      const auto _Jrat = __cyl_bessel_ratio_s_frac(__nu, __z, __zeta);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
	return std::real(_Jrat);
      else
	return _Jrat;
    }

  /**
   * 
   */
  template<typename _Tnu, typename _Tp,
	   typename _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>>
    std::conditional_t<__gnu_cxx::is_complex_v<_Val>,
			std::complex<__gnu_cxx::__num_traits_t<_Val>>,
			_Val>
    __cyl_bessel_i_ratio_s_frac(_Tnu __nu, _Tp __z)
    {
      const auto __zeta = __z * __z;
      const auto _Irat = __cyl_bessel_ratio_s_frac(__nu, __z, __zeta);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
	return std::real(_Irat);
      else
	return _Irat;
    }

template<typename _Tp>
  _Tp
  cyl_bessel_j_cf1(_Tp __nu, _Tp __x)
  {
    const auto _S_fp_min = __gnu_cxx::__sqrt_min(__nu);
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    constexpr int _S_max_iter = 15000;

    int __isign = 1;
    const auto __xi = _Tp{1} / __x;
    const auto __xi2 = _Tp{2} * __xi;
    auto __h = std::max(_S_fp_min, __nu * __xi);
    auto __b = __xi2 * __nu;
    auto __d = _Tp{0};
    auto __c = __h;
    int __i;
    for (__i = 1; __i <= _S_max_iter; ++__i)
      {
	__b += __xi2;
	__d = __b - __d;
	if (std::abs(__d) < _S_fp_min)
	  __d = _S_fp_min;
	__d = _Tp{1} / __d;
	__c = __b - _Tp{1} / __c;
	if (std::abs(__c) < _S_fp_min)
	  __c = _S_fp_min;
	const auto __del = __c * __d;
	__h *= __del;
	if (__d < _Tp{0})
	  __isign = -__isign;
	if (std::abs(__del - _Tp{1}) < _S_eps)
	  break;
      }
    return __h;
  }

template<typename _Tp>
  _Tp
  cyl_bessel_i_cf1(_Tp __nu, _Tp __x)
  {
    const auto _S_fp_min = __gnu_cxx::__sqrt_min(__nu);
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    constexpr int _S_max_iter = 15000;

    const auto __xi = _Tp{1} / __x;
    const auto __xi2 = _Tp{2} * __xi;
    auto __h = std::max(_S_fp_min, __nu * __xi);
    auto __b = __xi2 * __nu;
    auto __d = _Tp{0};
    auto __c = __h;
    int __i;
    for (__i = 1; __i <= _S_max_iter; ++__i)
      {
	__b += __xi2;
	__d = _Tp{1} / (__b + __d);
	__c = __b + _Tp{1} / __c;
	const auto __del = __c * __d;
	__h *= __del;
	if (std::abs(__del - _Tp{1}) < _S_eps)
	  break;
      }
    return __h;
  }

template<typename _Tp>
  void
  test_cyl_bessel_ratio()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 6 + std::cout.precision();
    using Ret = decltype(__cyl_bessel_j_ratio_s_frac(_Tp{1}, _Tp{1}));

    // What's with the NaN at 2?
    __cyl_bessel_j_ratio_s_frac(_Tp{1}, _Tp{2});

    std::cout << "\n\nRatio J_{\\nu+1}(x) / J_{\\nu}(x)\n";
    std::cout << ' ' << std::setw(w) << "z"
	      << ' ' << std::setw(w) << "new ratio"
	      << ' ' << std::setw(w) << "lab w/deriv ratio"
	      << ' ' << std::setw(w) << "nu / z - cf1"
	      << '\n';
    for (auto nu : {_Tp{0}, _Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    Ret r{};
	    try
	      {
		r = __cyl_bessel_j_ratio_s_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    const auto cf1 = cyl_bessel_j_cf1(nu, z);
	    const auto jn = std::__detail::__cyl_bessel_jn(nu, z);
	    const auto s = nu / z - jn.__J_deriv / jn.__J_value;
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
    for (auto nu : {_Tp{0}, _Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    Ret r{};
	    try
	      {
		r = __cyl_bessel_i_ratio_s_frac(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    const auto cf1 = cyl_bessel_i_cf1(nu, z);
	    const auto ik = std::__detail::__cyl_bessel_ik(nu, z);
	    const auto s = -nu / z + ik.__I_deriv / ik.__I_value;
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
