// $HOME/bin_specfun/bin/g++ -std=gnu++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -g -Wall -Wextra -o test_struve_new test_struve_new.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_struve_new > test_struve_new.new

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
//#include <complex>
#include <string>
#include <ext/math_const.h>
#include <bits/numeric_limits.h>
#include <bits/specfun_util.h>
#include <bits/complex_util.h>
#include <bits/summation.h>
#include <cmath>

  template<typename _Tp>
    _Tp
    __struve_series(_Tp __nu, _Tp __x, int __sign)
    {
      using _Val = std::__detail::__num_traits_t<_Tp>;

      constexpr auto _S_eps = std::numeric_limits<_Val>::epsilon();
      auto __x2 = __x / _Tp{2};
      auto __xx4 = _Tp(__sign) * __x2 * __x2;
      auto __term = _Tp{1};
      auto __struve = __term;
      for (int __k = 1; __k < 100; ++__k)
	{
      	  __term *= __xx4 / (_Tp(__k) + 1.5) / (__nu + _Tp(__k) + 1.5);
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      __struve *= std::pow(__x2, __nu + _Tp{1}) / std::detail::__gamma(_Tp{1.5L}) / std::detail::__gamma(_Tp{1.5L} + __nu);

      return __struve;
    }

  template<typename _Tp>
    _Tp
    __struve_h(_Tp __nu, _Tp __x)
    { return __struve_series(__nu, __x, -1); }

  template<typename _Tp>
    _Tp
    __struve_l(_Tp __nu, _Tp __x)
    { return __struve_series(__nu, __x, +1); }


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
	for (int n = 0; n <= 5; ++n)
	  {
	    auto nu = _Tp(1.0Q * n);
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::real(__struve_h(nu, t))
		 << std::setw(width) << std::real(__struve_l(nu, t))
		 << '\n';
	  }
      }
    data << "\n\n";
  }


int
main()
{
  using cmplx = std::complex<double>;
  plot_struve<cmplx>("plot/struve_double.txt");

  return 0;
}
