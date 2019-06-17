/**
 *
 */

#include <cmath>
#include <complex>

#include <ext/continued_fractions.h>
#include <bits/specfun_util.h>
#include <bits/complex_util.h> // is_complex

  template<typename _Tnu, typename _Tp, typename _Tzeta>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp, _Tzeta>>>
    __cyl_bessel_j_ratio(_Tnu __nu, _Tp __z, _Tzeta __zeta)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

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

      _SteedContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      _J(__a_J, __b_J, __w_J);

      // b_0 is 0 not 1 so subtract 1.
      return _J(__z) - _Real{1};
    }

  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_j_ratio(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __iz = _Cmplx{0, 1} * __z;
      const auto __zeta = __iz * __iz;
      const auto _Jrat = __cyl_bessel_j_ratio(__nu, __z, __zeta);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
	return std::real(_Jrat);
      else
	return _Jrat;
    }

  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_i_ratio(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = __z * __z;
      const auto _Irat = __cyl_bessel_j_ratio(__nu, __z, __zeta);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
	return std::real(_Irat);
      else
	return _Irat;
    }

template<typename _Tp>
  void
  test_cyl_bessel_ratio()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 6 + std::cout.precision();
    using Ret = decltype(__cyl_hankel_j_ratio(_Tp{1}, _Tp{1}));

    std::cout << "\n\nRatio J_{\\nu+1}(x) / J_{\\nu}(x)\n";
    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    Ret r{};
	    try
	      {
		r = __cyl_hankel_j_ratio(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    const auto jn = std::__detail::__cyl_bessel_jn(nu, z);
	    const auto s = nu / z - jn.__J_deriv / jn.__J_value;
	    std::cout << ' ' << std::setw(w) << z
		      << ' ' << std::setw(w) << r
		      << ' ' << std::setw(w) << s
		      << '\n';
	  }
      }

    std::cout << "\n\nRatio I_{\\nu+1}(x) / I_{\\nu}(x)\n";
    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    Ret r{};
	    try
	      {
		r = __cyl_hankel_i_ratio(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    const auto ik = std::__detail::__cyl_bessel_ik(nu, z);
	    const auto s = -nu / z + ik.__I_deriv / ik.__I_value;
	    std::cout << ' ' << std::setw(w) << z
		      << ' ' << std::setw(w) << r
		      << ' ' << std::setw(w) << s
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
