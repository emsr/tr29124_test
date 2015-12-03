
// $HOME/bin_specfun/bin/g++ -std=c++14 -o test_beta_inc test_beta_inc.cpp

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_beta_inc

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
      constexpr auto _S_maxiter = 100;
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
      int __m = 0;
      for (__m = 1; __m <= _S_maxiter; ++__m)
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
	    break;
	}
      if (__m > _S_maxiter)
	std::__throw_logic_error("betacf: a or b too big, or _S_maxiter too small");
      return __h;
    }

  ///  Returns the incomplete beta function I_x(a;b).
  template<typename _Tp>
    _Tp
    __ibeta(_Tp __a, _Tp __b, _Tp __x)
    {
      _Tp __bt;
      if (__x < _Tp{0} || __x > _Tp{1})
	throw std::domain_error("betai: bad argument");
      if (__x == _Tp{0} || __x == _Tp{1})
	__bt = _Tp{0};
      else // Factors in front of the continued fraction.
	__bt = std::exp(std::lgamma(__a + __b)
	     - std::lgamma(__a) - std::lgamma(__b)
             + __a * std::log(__x) + __b * std::log(_Tp{1} - __x));
      if (__x < (__a + _Tp{1}) / (__a + __b + _Tp{2})) //  Use continued fraction directly.
	return __bt * __beta_cont_frac(__a, __b, __x) / __a;
      else  //  Use continued fraction after making the symmetry transformation.
	return _Tp{1} - __bt * __beta_cont_frac(__b, __a, _Tp{1} - __x) / __b;
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
