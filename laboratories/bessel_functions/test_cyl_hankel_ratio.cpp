
/**
 *
 */

#include <complex>
#include <iostream>
#include <iomanip>

#include <ext/continued_fractions.h>
#include <bits/specfun_util.h>

  template<typename _Tnu, typename _Tp, typename _Tzeta>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp, _Tzeta>>>
    cyl_hankel_ratio(_Tnu __nu, _Tp __z, _Tzeta __zeta)
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
      using _AFun = decltype(__a_H1);

      auto __b_H1 = [](std::size_t, _Tp) -> _Real { return _Real{1}; };
      using _BFun = decltype(__b_H1);

      auto __w_H1
	= [__nu, __zeta, __c](std::size_t __k, _Tp)
	  {
	    return (_Real{-1}
		    + std::sqrt(_Real{1} + _Real{8} * __c(__k) * __zeta))
		 / _Real{2};
	  };
      using _TailFun = decltype(__w_H1);

      _SteedContinuedFraction<_Tp, _AFun, _BFun, _TailFun>
      _H1(__a_H1, __b_H1, __w_H1);

      return _H1(__z);
    }

  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    cyl_hankel_h1_ratio(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = _Real{-1} / (_Cmplx{0, 2} * __z);
      return -cyl_hankel_ratio(__nu, __z, __zeta);
    }

  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    cyl_hankel_h2_ratio(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = _Real{1} / (_Cmplx{0, 2} * __z);
      return cyl_hankel_ratio(__nu, __z, __zeta);
    }

  template<typename _Tnu, typename _Tp>
    std::complex<std::__detail::__num_traits_t<
		__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    cyl_k_ratio(_Tnu __nu, _Tp __z)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = std::__detail::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = _Real{1} / (_Real{2} * __z);
      return cyl_hankel_ratio(__nu, __z, __zeta);
    }

template<typename _Tp>
  void
  test_cyl_hankel_ratio()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 6 + std::cout.precision();
    using ret = decltype(cyl_hankel_h1_ratio(_Tp{1}, _Tp{1}));

    for (auto nu : {_Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n nu = " << std::setw(w) << nu << '\n';
	for (int i = 1; i <= 20; ++i)
	  {
	    auto z = _Tp(i) / _Tp{10};
	    ret r{};
	    try
	      {
		r = cyl_hankel_h1_ratio(nu, z);
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
  test_cyl_hankel_ratio<double>();

  return 0;
}
