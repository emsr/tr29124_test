/**
 *
 */

#include <complex>

#include <ext/continued_fractions.h>
#include <bits/specfun_util.h>

  template<typename _Tnu, typename _Tp, typename _Tzeta>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp, _Tzeta>>>
    cyl_bessel_ratio(_Tnu __nu, _Tp __z, _Tzeta __zeta)
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

      return _J(__z);
    }

  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    cyl_hankel_j_ratio(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __iz = _Cmplx{0, 1} * __z;
      const auto __zeta = __iz * __iz;
      return cyl_bessel_ratio(__nu, __z, __zeta);
    }

  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    cyl_hankel_i_ratio(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = __z * __z;
      return cyl_bessel_ratio(__nu, __z, __zeta);
    }

template<typename _Tp>
  void
  test_cyl_bessel_ratio()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 6 + std::cout.precision();
    using ret = decltype(cyl_hankel_j_ratio(_Tp{1}, _Tp{1}));

    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    ret r{};
	    try
	      {
		r = cyl_hankel_j_ratio(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    std::cout << ' ' << std::setw(w) << z
		      << ' ' << std::setw(w) << r
		      << '\n';
	  }
      }

    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    ret r{};
	    try
	      {
		r = cyl_hankel_i_ratio(nu, z);
	      }
	    catch (...)
	      {
		std::cout << ' ' << std::setw(w) << z
			  << ' ' << std::setw(w) << "FAIL\n";
		continue;
	      }
	    std::cout << ' ' << std::setw(w) << z
		      << ' ' << std::setw(w) << r
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
