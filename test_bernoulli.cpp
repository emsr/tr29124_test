
/*
$HOME/bin_tr29124/bin/g++ -o test_bernoulli test_bernoulli.cpp
./test_bernoulli > test_bernoulli.txt
*/ 

#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <limits>
#include "cmath_variable_template"


template<typename _RealTp, typename _IntTp,
	 _IntTp _Num = 1, _IntTp _Den = 1>
  constexpr _RealTp
  __frac = _RealTp(_Num) / _RealTp(_Den);

template<typename _RealTp,
	 unsigned long long _Num = 1,
         unsigned long long _Den = 1>
  constexpr _RealTp
  frac = __frac<_RealTp, unsigned long long, _Num, _Den>;


//
//  Recursion is unstable.
//
template<typename _Tp>
  _Tp
  __bernoulli_series(unsigned int __n)
  {
    static constexpr _Tp
    __num[24]
    {
       frac<_Tp>,
      -frac<_Tp, 1UL, 2UL>,
       frac<_Tp, 1UL, 6UL>,            _Tp(0UL),
      -frac<_Tp, 1UL, 30UL>,           _Tp(0UL),
       frac<_Tp, 1UL, 42UL>,           _Tp(0UL),
      -frac<_Tp, 1UL, 30UL>,           _Tp(0UL),
       frac<_Tp, 5UL, 66UL>,           _Tp(0UL),
      -frac<_Tp, 691UL, 2730UL>,       _Tp(0UL),
       frac<_Tp, 7UL, 6UL>,            _Tp(0UL),
      -frac<_Tp, 3617UL, 510UL>,       _Tp(0UL),
       frac<_Tp, 43867UL, 798UL>,      _Tp(0UL),
      -frac<_Tp, 174611UL, 330UL>,     _Tp(0UL),
       frac<_Tp, 854513UL, 138UL>,     _Tp(0UL)
    };

    if (__n == 0)
      return frac<_Tp>;

    if (__n == 1)
      return -frac<_Tp, 1, 2>;

    //  Take care of the rest of the odd ones.
    if (__n % 2 == 1)
      return _Tp(0);

    //  Take care of some small evens that are painful for the series.
    if (__n < 28)
      return __num[__n];


    _Tp __fact = _Tp(1);
    if ((__n / 2) % 2 == 0)
      __fact *= -1;
    for (unsigned int __k = 1; __k <= __n; ++__k)
      __fact *= __k / __gnu_cxx::__const_pi<_Tp>();
    __fact *= _Tp(2);

    _Tp __sum = _Tp(0);
    for (unsigned int __i = 1; __i < 1000; ++__i)
      {
	_Tp __term = std::pow(_Tp(__i), -_Tp(__n));
	if (__term < __gnu_cxx::__epsilon<_Tp>())
          break;
	__sum += __term;
      }

    return __fact * __sum;
  }


int
main()
{

  std::cout.precision(12);
  std::cout.flags(std::ios::showpoint);

  std::cout << std::endl;
  for (unsigned int n = 0; n <= 200; ++n)
    std::cout << "  " << std::setw(4) << n << "  " << std::setw(14) << __bernoulli_series<double>(n) << std::endl;

  return 0;
}

