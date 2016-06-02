// $HOME/bin_specfun/bin/g++ -std=gnu++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -g -Wall -Wextra -o test_struve_new test_struve_new.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_struve_new > test_struve_new.new

// g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -Wall -Wextra -o test_struve_new test_struve_new.cpp

// ./test_struve_new > test_struve_new.txt

#include <cmath> // There are issues with <complex> inclusion if this isn't uphere!
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <string>
#include <ext/math_const.h>
#include <bits/numeric_limits.h>
#include <bits/specfun_util.h>
#include <bits/complex_util.h>
#include <bits/summation.h>

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __struve_series(_Tp __nu, _Tp __x, int __sign)
    {
      using _Val = std::__detail::__num_traits_t<_Tp>;

      constexpr auto _S_eps = std::numeric_limits<_Val>::epsilon();
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Val>::__root_pi;

      auto __x2 = __x / _Tp{2};
      auto __xx4 = _Tp(__sign) * __x2 * __x2;
      auto __term = _Tp{1};
      auto __struve = __term;
      for (int __k = 1; __k < 100; ++__k)
	{
      	  __term *= __xx4 / _Tp(__k - 1 + 1.5L) / (__nu + _Tp(__k - 1 + 1.5L));
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      __struve *= _Tp{2} * std::pow(__x2, __nu + _Tp{1})
		/ std::__detail::__gamma(__nu + _Tp{1.5L}) / _S_sqrt_pi;

      return __struve;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __struve_asymp(_Tp __nu, _Tp __x, int __sign)
    {
      using _Val = std::__detail::__num_traits_t<_Tp>;

      constexpr auto _S_eps = std::numeric_limits<_Val>::epsilon();
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Val>::__root_pi;

      auto __x2 = __x / _Tp{2};
      auto __xx4 = _Tp(__sign) * __x2 * __x2;
      auto __term = _Tp{1};
      auto __struve = __term;
      for (int __k = 1; __k < 100; ++__k)
	{
      	  __term *= _Tp(__k - 1 + 0.5L) / (_Tp(-__k - 1 + 0.5L) + __nu) / __xx4;
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      __struve *= _Tp(__sign) * _Tp{2} * std::pow(__x2, __nu - _Tp{1})
		/ std::__detail::__gamma(__nu + _Tp{0.5L}) / _S_sqrt_pi;

      return __struve;
    }

  template<typename _Tp>
    _Tp
    __struve_h(_Tp __nu, _Tp __x)
    { return __struve_series(__nu, __x, -1); }

  template<typename _Tp>
    _Tp
    __struve_k(_Tp __nu, _Tp __x)
    { return __struve_asymp(__nu, __x, +1); }

  template<typename _Tp>
    _Tp
    __struve_l(_Tp __nu, _Tp __x)
    { return __struve_series(__nu, __x, +1); }

  template<typename _Tp>
    _Tp
    __struve_m(_Tp __nu, _Tp __x)
    { return __struve_asymp(__nu, __x, -1); }


/**
 * 
 */
template<typename _Tp>
  void
  plot_struve(std::string filename)
  {
    using _Val = std::__detail::__num_traits_t<_Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Val>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "t"
	 << std::setw(width) << "H"
	 << std::setw(width) << "L"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    for (int i = 0; i <= +2000; ++i)
      {
	auto t = _Tp(0.01Q * i);
	data << std::setw(width) << t;
	for (int n = 0; n <= 20; ++n)
	  {
	    auto nu = _Tp(1.0Q * n);
	    data << '\t'
		 << std::setw(width) << __struve_h(nu, t)
		 << std::setw(width) << __struve_l(nu, t);
	  }
	data << '\n';
      }
    data << "\n\n";
  }


int
main()
{
  //using cmplx = std::complex<double>;
  plot_struve<double>("plot/struve_double.txt");

  return 0;
}
