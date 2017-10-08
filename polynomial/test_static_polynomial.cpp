/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_static_polynomial test_static_polynomial.cpp -lquadmath
./test_static_polynomial

g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -o test_static_polynomial test_static_polynomial.cpp -lquadmath
*/

//  Get past a bug....
// $HOME/bin/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__=0 -o test_static_polynomial test_static_polynomial.cpp

#include "static_polynomial.h"
#include <initializer_list>
#include <iostream>
#include <complex>
#include <sstream>

/* Nope.
template<typename _CoefTp, typename _ArgTp>
  decltype(_CoefTp{} * _ArgTp{})
  static_polynomial_help(_CoefTp... a, _CoefTp a0, _ArgTp x)
  { return a0 + x * static_polynomial_help(a..., x); }

template<typename _CoefTp, typename _ArgTp>
  decltype(_CoefTp{} * _ArgTp{})
  static_polynomial_help(_CoefTp a0, _ArgTp)
  { return a0; }

template<typename _CoefTp, typename _ArgTp>
  decltype(_CoefTp{} * _ArgTp{})
  static_polynomial(_CoefTp... a, _ArgTp x)
  { return static_polynomial_help(a..., x); }
*/

/*
template<typename _CoefTp, typename _ArgTp>
  constexpr decltype(_CoefTp{} * _ArgTp{})
  static_polynomial(std::initializer_list<_CoefTp> __a, _ArgTp __x)
  {
    using _RetTp = decltype(_CoefTp{} * _ArgTp{});
    constexpr auto _Size = __a.size();
    if (_Size == 0)
      return _RetTp{0};
    else if (_Size == 1)
      return _RetTp{*__a.begin()};
    else
      {
	constexpr const auto __aiter = __a.begin();
      	_RetTp __sum = *__aiter;
	//while(++__aiter != __a.end())
	//  ;
	return __sum;
      }
  }
*/

int
main()
{
  float aa[5]{1.0f, 2.0f, -1.5f, 0.2f, -0.1f};

  __gnu_cxx::_StaticPolynomial<float, 5> p(aa);
  //auto y = p(1.2);
  __gnu_cxx::_StaticPolynomial<float, 5> a{{1.0f, 2.0f, -1.5f, 0.2f, -0.1f}};

  //auto y = static_polynomial({1.0f, 2.0f, -1.5f, 0.2f, -0.1f}, 2.0);
  //return 0;
}
