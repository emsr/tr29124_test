/**
 *
 */

#include <stdexcept>
#include <cmath>
#include <ext/math_const.h>
#include <bits/numeric_limits.h>

template<typename _Tp>
  _Tp
  cyl_bessel_nk_series(_Tp __nu, _Tp __x, int __max_iter = 100)
  {
    const auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
    const auto _S_pi = __gnu_cxx::math::__pi_v<_Tp>;
    const auto _S_fp_min = __gnu_cxx::__sqrt_min(__nu);
    const int __n = std::nearbyint(__nu);
    const auto __mu = __nu - _Tp(__n);
    const auto __mu2 = __mu * __mu;
    const auto __xi = _Tp{1} / __x;
    const auto __xi2 = _Tp{2} * __xi;
    const auto _Wronski = __xi2 / _S_pi;
    const auto __x2 = __x / _Tp{2};

    int __isign = 1;
    auto __h = std::max(_S_fp_min, __nu * __xi);
    { // CF1
      auto __b = __xi2 * __nu;
      auto __d = _Tp{0};
      auto __c = __h;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
	{
	  __b += __xi2;
	  __d = __b - __d;
	  if (std::abs(__d) < _S_fp_min)
	    __d = _S_fp_min;
	  __d = _Tp{1} / __d;
	  __c = __b - _Tp{1} / __c;
	  if (std::abs(__c) < _S_fp_min)
	    __c = _S_fp_min;
	  const auto __del = __c * __d;
	  __h *= __del;
	  if (__d < _Tp{0})
	    __isign = -__isign;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	}
      if (__i > __max_iter)
	return __cyl_bessel_jn_asymp(__nu, __x);
    }
    auto _Jnuk = __isign * _S_fp_min;
    auto _Jpnuk = __h * _Jnuk;
    auto _Jnuk1 = _Jnuk;
    auto _Jpnu1 = _Jpnuk;
    auto __fact0 = __nu * __xi;
    for (int __k = __n; __k >= 1; --__k)
      {
	const auto _Jnutemp = __fact0 * _Jnuk + _Jpnuk;
	__fact0 -= __xi;
	_Jpnuk = __fact0 * _Jnutemp - _Jnuk;
	_Jnuk = _Jnutemp;
      }
    if (_Jnuk == _Tp{0})
      _Jnuk = _S_eps;

    const auto __f = _Jpnuk / _Jnuk;

    const auto __fact = _Tp{1} / __sinc_pi(__mu);
    auto __d = -std::log(__x2);
    auto __e = __mu * __d;
    const auto __fact2 = __sinhc(__e);
    const auto __gamt = __gamma_temme(__mu);
    auto __ff = (_Tp{2} / _S_pi) * __fact
	      * (__gamt.__gamma_1_value * std::cosh(__e)
	       + __gamt.__gamma_2_value * __fact2 * __d);
    __e = std::exp(__e);
    auto __p = __e / (_S_pi * __gamt.__gamma_plus_value);
    auto __q = _Tp{1} / (__e * _S_pi * __gamt.__gamma_minus_value);
    const auto __fact3 = __sinc_pi(__mu / _Tp{2});
    const auto __pifact3 = _S_pi * __fact3;
    const auto __r = __pifact3 * __pifact3 * __mu2 / _Tp{2};
    auto __c = _Tp{1};
    __d = -__x2 * __x2;
    auto __sum = __ff + __r * __q;
    auto __sum1 = __p;
    int __i;
    for (__i = 1; __i <= __max_iter; ++__i)
      {
	__ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	__c *= __d / _Tp(__i);
	__p /= _Tp(__i) - __mu;
	__q /= _Tp(__i) + __mu;
	const auto __del = __c * (__ff + __r * __q);
	__sum += __del;
	const auto __del1 = __c * __p - _Tp(__i) * __del;
	__sum1 += __del1;
	if (std::abs(__del) < _S_eps * (_Tp{1} + std::abs(__sum)))
	  break;
      }
    if (__i > __max_iter)
      std::__throw_runtime_error(__N("__cyl_bessel_nk_series: "
				     "Series failed to converge"));
    auto _Nmu = -__sum;
    auto _Nmup1 = -__sum1 * __xi2;
    auto _Npmu = __mu * __xi * _Nmu - _Nmup1;
    auto _Jmu = _Wronski / (_Npmu - __f * _Nmu);

    // Recur _Nmu back up...
    for (int __k = 0; __k < __n; ++__k)
      {
	
      }

    return _Nmu;
  }

int
main()
{
}
