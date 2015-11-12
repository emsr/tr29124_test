// $HOME/bin/bin/g++ -std=c++14 -o dilog dilog.cpp

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>

#include "specfun_util.h"

namespace std{
namespace __detail
{
template<typename _Tp>
  _Tp
  __dilog(_Tp __x)
  {
    static constexpr unsigned long long _S_maxit = 100000ULL;
    static constexpr _Tp _S_eps = 10 * std::numeric_limits<_Tp>::epsilon();
    static constexpr _Tp _S_pipio6
      = 1.64493406684822643647241516664602518921894990120679843773555822937000747040320087383362890061975870L;
    if (__isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__x == _Tp(1))
      return _S_pipio6;
    else if (__x == -_Tp(1))
      return -0.5L * _S_pipio6;
    else if (__x > _Tp(0.5))
      return _S_pipio6 - std::log(__x) * std::log(_Tp(1) - __x)
	   - __dilog(_Tp(1) - __x);
    else if (__x < -_Tp(0.5))
      return -0.5L * _S_pipio6 - std::log(_Tp(1) + __x) * std::log(-__x)
	   + __dilog(_Tp(1) + __x) - __dilog(_Tp(1) - __x * __x);
    else
      {
	_Tp __sum = 0;
	_Tp __fact = 1;
	for (auto __i = 1ULL; __i < _S_maxit; ++__i)
	  {
	    __fact *= __x;
	    auto __term = __fact / (__i * __i);
	    __sum += __term;
	    if (std::abs(__term) < _S_eps)
	      {
		std::cout << __i << " - ";
		break;
	      }
	    if (__i + 1 == _S_maxit)
	      std::__throw_runtime_error("__dilog: sum failed");
	  }
	return __sum;
      }
  }
}
}

int
main()
{
  std::cout.precision(15);
  std::cout << '\n';
  for (auto i = -1000; i <= +1000; ++i)
    std::cout << "dilog(" << i * 0.001 << ") = " << std::__detail::__dilog(i * 0.001) << '\n';
}
