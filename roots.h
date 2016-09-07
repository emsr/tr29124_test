#ifndef ROOTS_H
#define ROOTS_H 1

#include <vector>
//#include <utility> // For pair?
#include <limits>

// Default __eps = _Tp{5} * std::numeric_limits<_Tp>::epsilon()?

// template<typename _Tp>
//   struct __bracket_t {_Tp __value_lower, _Tp __value_upper}; ?

// All the root methods should have the same API.
// brent, newton, safe don't have max_iter.

// __root_brackets should have an output iterator API?

// in __root_newton what if __df is zero?

namespace __gnu_cxx
{

  template<typename _Tp>
    bool
    __root_bracket(_Tp (*__func)(_Tp), _Tp& __x_lower, _Tp& __x_upper,
    		   std::size_t __max_iter = 50);

  template<typename _Tp>
    std::vector<std::pair<_Tp, _Tp>>
    __root_brackets(_Tp (*__func)(_Tp),
		    _Tp __x_lower, _Tp __x_upper, std::size_t __n);

  template<typename _Tp>
    _Tp
    __root_bisect(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper, _Tp __eps,
		  std::size_t __max_iter = std::numeric_limits<_Tp>::digits);

  template<typename _Tp>
    _Tp
    __root_secant(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper, _Tp __eps,
		  std::size_t __max_iter = 40);

  template<typename _Tp>
    _Tp
    __root_false_position(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
			  _Tp __eps, std::size_t __max_iter = 40);

  template<typename _Tp>
    _Tp
    __root_ridder(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter = 100);

  template<typename _Tp>
    _Tp
    __root_brent(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
		 _Tp __eps, std::size_t __max_iter = 100);

  template<typename _Tp>
    _Tp
    __root_newton(void (*__func)(_Tp, _Tp*, _Tp*), _Tp __x_lower, _Tp __x_upper,
    		  _Tp __eps, std::size_t __max_iter = 40);

  template<typename _Tp>
    _Tp
    __root_safe(void (*__func)(_Tp, _Tp*, _Tp*), _Tp __x_lower, _Tp __x_upper,
    		_Tp __eps, std::size_t __max_iter = 100);

} // namespace __gnu_cxx

#include "roots.tcc"

#endif // ROOTS_H
