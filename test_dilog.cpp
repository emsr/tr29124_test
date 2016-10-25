/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -I. -o test_dilog test_dilog.cpp
./test_dilog > test_dilog.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>

namespace std
{
namespace __detail
{

  template<typename _Tp>
    _Tp
    __dilog(_Tp __x)
    {
      static constexpr unsigned long long _S_maxit = 100000ULL;
      static constexpr _Tp _S_eps = 10 * std::numeric_limits<_Tp>::epsilon();
      static constexpr _Tp _S_pipio6
	= 1.644934066848226436472415166646025189219L;
      if (__isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x > _Tp(+1))
	std::__throw_range_error(__N("dilog: argument greater than one"));
      else if (__x < _Tp(-1))
	{
	  auto __lnfact = std::log(_Tp(1) - __x);
	  return -__dilog(_Tp(1) - _Tp(1) / (_Tp(1) - __x))
		 - _Tp(0.5L) * __lnfact * __lnfact;
	}
      else if (__x == _Tp(1))
	return _S_pipio6;
      else if (__x == -_Tp(1))
	return -_Tp(0.5L) * _S_pipio6;
      else if (__x > _Tp(0.5))
	return _S_pipio6 - std::log(__x) * std::log(_Tp(1) - __x)
	     - __dilog(_Tp(1) - __x);
      else if (__x < -_Tp(0.5))
	return -_Tp(0.5L) * _S_pipio6 - std::log(_Tp(1) + __x) * std::log(-__x)
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

} // namespace __detail
} // namespace std

namespace __gnu_cxx
{

  //  Dilogarithm functions

  inline float
  dilogf(float __x)
  { return std::__detail::__dilog<float>(__x); }

  inline long double
  dilogl(long double __x)
  { return std::__detail::__dilog<long double>(__x); }

  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    dilog(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__dilog<__type>(__x);
    }

}

int
main()
{
  std::cout.precision(15);
  std::cout << '\n';
  for (auto i = -4000; i <= +1000; ++i)
    std::cout << "dilog(" << i * 0.001 << ") = " << __gnu_cxx::dilog(i * 0.001) << '\n';
}
