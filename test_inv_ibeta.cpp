/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_inv_ibeta test_inv_ibeta.cpp
./test_inv_ibeta > test_inv_ibeta.txt
*/

//#include <ext/cmath>
#include <iostream>
#include <iomanip>
#include <bits/specfun.h>
#include <ext/root_finding.h>

template<typename _Tp>
  _Tp
  __ibeta_inv(_Tp __a, _Tp __b, _Tp _Ibeta)
  {
    constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    constexpr int _S_max_inner_iter = 12;
    constexpr int _S_max_outer_iter = 35;
    if (__a < _Tp{0} || __b < _Tp{0})
      std::__throw_domain_error(__N("__ibeta_inv: "
				"parameters a and b must be positive"));
    else if (_Ibeta < _Tp{0} || _Ibeta > _Tp{1})
      std::__throw_domain_error(__N("__ibeta_inv: "
				"incomplete beta function out of range"));
    const auto _Beta = std::beta(__a, __b);

    __gnu_cxx::__root_bracket([__a, __b](_Tp __x)->_Tp{ return __gnu_cxx::ibeta(__a, __b, __x); },
		 _Tp{0}, _Tp{1});
/*
    auto _Ix = _Ibeta;
    auto _Iy = _Tp{1} - _Ix;
    auto _Ixy_prev = std::min(_Ix, _Iy);
    _Tp __x_prev, __y_prev;
    if (__a <= __b)
      {
	__x_prev = std::pow(__a * _Ix * _Beta, _Tp{1} / __a);
	__y_prev = _Tp{1} - __x_prev;
      }
    else
      {
	__y_prev = std::pow(__b * _Iy * _Beta, _Tp{1} / __b);
	__x_prev = _Tp{1} - __y_prev;
      }
    auto __xy_prev = std::min(__x_prev, __y_prev);
    auto __xy_lower = _Tp{0};
    auto __xy_upper = _Tp{1};
    int __outer_iter = 0;
    while (__outer_iter++ < _S_max_outer_iter)
      {
	int __iter = 0;
	while (__iter++ < _S_max_inner_iter && _Ixy_prev < _Ixy)
	  {
	    __xy = _Tp{2} * __xy_prev;
	    __xy_lower = __xy;
	  }
	__iter = 0;
	while (__iter++ < _S_max_inner_iter && _Ixy_prev >= _Ixy)
	  {
	    __xy = __xy_prev / _Tp{2};
	    __xy_upper = __xy;
	  }
	__xy = (__xy_lower + __xy_upper) / _Tp{2};
	if (std::abs(__xy_upper - __xy_lower) < _S_eps * __xy)
	  break;
      }
*/
    auto __thing = [__a, __b](_Tp __x){ return __gnu_cxx::ibeta(__a, __b, __x); };
    auto __xy_lower = _Tp{0};
    auto __xy_upper = _Tp{1};
    if (__gnu_cxx::__root_bracket(__thing, __xy_lower, __xy_upper))
      {
	return __gnu_cxx::__root_brent(__thing, __xy_lower, __xy_upper, _Tp{10} * _S_eps);
      }
  }

int
main()
{
  double a = 2.4, b = 5.6;
  double x = 0.8;
  double ibet = __gnu_cxx::ibeta(a, b, x);
  auto y = __ibeta_inv(a, b, ibet);
}
