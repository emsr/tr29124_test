/*
$HOME/bin_tr29124/bin/g++ -std=c++14 -o test_beta_inc test_beta_inc.cpp
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_beta_inc
*/

#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <iomanip>


  //  Evaluates the continued fraction for the incomplete beta function
  //  by the modified Lentz's method
  template<typename _Tp>
    _Tp
    __beta_cont_frac(_Tp __a, _Tp __b, _Tp __x)
    {
      constexpr auto _S_itmax = 100;
      constexpr auto _S_fpmin = std::numeric_limits<_Tp>::min();
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      auto __qab = __a + __b;
      auto __qap = __a + _Tp{1};
      auto __qam = __a - _Tp{1};
      auto __c = _Tp{1};
      //  First step of Lentz's method.
      auto __d = _Tp{1} - __qab * __x / __qap;
      if (std::abs(__d) < _S_fpmin)
	__d = _S_fpmin;
      __d = _Tp{1} / __d;
      auto __h = __d;
      for (int __m = 1; __m <= _S_itmax; ++__m)
	{
	  auto __m2 = 2 * __m;
	  auto __aa = __m * (__b - __m) * __x
		     / ((__qam + __m2) * (__a + __m2));
	  __d = _Tp{1} + __aa * __d;

	  //  Even step of the recurrence.
	  if (std::abs(__d) < _S_fpmin)
	    __d = _S_fpmin;
	  __c = _Tp{1} + __aa / __c;
	  if (std::abs(__c) < _S_fpmin)
	    __c = _S_fpmin;
	  __d = _Tp{1} / __d;
	  __h *= __d * __c;
	  __aa = -(__a + __m) * (__qab + __m) * __x
		/ ((__a + __m2) * (__qap + __m2));
	  __d = _Tp{1} + __aa * __d;

	  //  Odd step of the recurrence.
	  if (std::abs(__d) < _S_fpmin)
	    __d = _S_fpmin;
	  __c = _Tp{1} + __aa / __c;
	  if (std::abs(__c) < _S_fpmin)
	    __c = _S_fpmin;
	  __d = _Tp{1} / __d;
	  auto __del = __d * __c;
	  __h *= __del;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    return __h;
	}
      std::__throw_runtime_error(__N("__beta_cont_frac: "
				     "a or b too big, or _S_itmax too small"));
    }

  ///  Returns the incomplete beta function I_x(a;b).
  template<typename _Tp>
    _Tp
    __ibeta(_Tp __a, _Tp __b, _Tp __x)
    {
      if (__x < _Tp{0} || __x > _Tp{1})
	std::__throw_domain_error(__N("__ibeta: argument out of range"));
      if (__isnan(__x) || __isnan(__a) || __isnan(__b))
	return __gnu_cxx::__math_constants<_Tp>::__NaN;
      else
	{
	  _Tp __fact;
	  if (__x == _Tp{0} || __x == _Tp{1})
	    __fact = _Tp{0};
	  else
	    __fact = std::exp(std::lgamma(__a + __b)
			    - std::lgamma(__a) - std::lgamma(__b)
			    + __a * std::log(__x) + __b * std::log(_Tp{1} - __x));
	  if (__x < (__a + _Tp{1}) / (__a + __b + _Tp{2}))
	    return __fact * __beta_cont_frac(__a, __b, __x) / __a;
	  else
	    return _Tp{1}
	         - __fact * __beta_cont_frac(__b, __a, _Tp{1} - __x) / __b;
	}
    }

int
main()
{
  for (int ia = 0; ia <= 10; ++ia)
    {
      double a = 1.0 * ia;
      std::cout << "a = " << std::setw(6) << a << '\n';
      for (int ib = 10; ib >= 0; --ib)
	{
	  double b = 1.0 * ib;
	  std::cout << "b = " << std::setw(6) << b << '\n';
	  for (int ix = 0; ix <= 100; ++ix)
	    {
	      double x = ix * 0.01;
	      std::cout << ' ' << std::setw(6) << x
		        << ' ' << std::setw(12) << __ibeta(a, b, x) << '\n';
	    }
	}
    }
}
