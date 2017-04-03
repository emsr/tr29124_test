/* quadrature/double_exp_integrate.tcc
 *
 * Copyright (C) 2017 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef DOUBLE_EXP_INTEGRATE_TCC
#define DOUBLE_EXP_INTEGRATE_TCC 1

#include <ext/cmath>

namespace __gnu_test
{

/**
 * @f[
 *    \int_{-1}^{+1}f(x)dx
 * @f]
 * Making the change of variables:
 * @f[
 *    x = tanh\left[\frac{\pi}{2}sinh(u)\right],
 *   dx = \frac{\pi}{2}\frac{cosh(u)}{cosh^2\left[\frac{\pi}{2}sinh(u)\right]}du
 * @f]
 * gives the following integral:
 * @f[
 *    \int_{-\infty}^{+\infty}f(tanh\left[\frac{\pi}{2}sinh(u)\right])
 *     
 * @f]
 * 
 * @f[
 *    
 * @f]
 */
template<typename _Func, typename _Tp>
  _Tp
  double_exp_integrate1(_Func __func, int __n, _Tp __a, _Tp __b, _Tp /*__tol*/)
  {
    __n /= 2;
    const auto __h = _Tp{5} / __n; // 5.0 is rough limit of K in exp(exp(K)).

    const auto _S_pi_4 = __gnu_cxx::__const_pi_quarter<_Tp>();
    auto __sum = _Tp{0};
    for (int __k = -__n; __k <= __n; ++__k)
      {
	const auto __u = __h * _Tp(__k);
	const auto __eu = std::exp(__u);
	const auto __cosh = __eu + _Tp{1} / __eu;
	const auto __sinh = __eu - _Tp{1} / __eu;
        const auto __s = std::exp(_S_pi_4 * __sinh);
	const auto __w = __s + _Tp{1} / __s;
	const auto __x = (__b * __s + __a / __s) / __w;
	if (__x != __a && __x != __b)
	  {
	    const auto __dxdu = _Tp{2} * (__b - __a) * _S_pi_4 * __cosh
	    	     / (__w * __w);
	    __sum += __h * __func(__x) * __dxdu;
	  }
      }
    return __sum;
  }

/**
 * This uses the symmetry of cosh.
 */
template<typename _Func, typename _Tp>
  _Tp
  double_exp_integrate2(_Func __func, int __n, _Tp __a, _Tp __b, _Tp /*__tol*/)
  {
    __n /= 2;
    const auto __h = _Tp{5} / __n; // 5.0 is rough limit of K in exp(exp(K)).

    const auto _S_pi_4 = __gnu_cxx::__const_pi_quarter<_Tp>();
    auto __sum = __func((__a + __b) / _Tp{2}) / _Tp{2};
    for (int __k = -__n; __k < 0; ++__k)
      {
	const auto __u = __h * _Tp(__k);
	const auto __eu = std::exp(__u);
	const auto __cosh = __eu + _Tp{1} / __eu;
	const auto __sinh = __eu - _Tp{1} / __eu;
        const auto __s = std::exp(_S_pi_4 * __sinh);
	const auto __w = __s + _Tp{1} / __s;
	const auto __dxdu = __cosh / (__w * __w);
	const auto __x1 = (__b * __s + __a / __s) / __w;
	const auto __x2 = (__a * __s + __b / __s) / __w;
	if (__x1 != __a && __x1 != __b) 
	  __sum += __dxdu * __func(__x1);
	if (__x2 != __a && __x2 != __b)
	  __sum += __dxdu * __func(__x2);
      }
    return _Tp{2} * (__b - __a) * _S_pi_4 * __h * __sum;
  }

/**
 * This is like the first except it interlaces another set of points
 * for a higher-order integral..
 */
template<typename _Func, typename _Tp>
  _Tp
  double_exp_integrate3(_Func __func, int __n, _Tp __a, _Tp __b, _Tp /*__tol*/)
  {
    __n /= 2;
    const auto __h = _Tp{5} / __n; // 5.0 is rough limit of K in exp(exp(K)).

    const auto _S_pi_4 = __gnu_cxx::__const_pi_quarter<_Tp>();
    auto __sum = _Tp{0};
    for (int __k = -__n; __k <= __n; ++__k)
      {
	const auto __u = __h * _Tp(__k);
	const auto __eu = std::exp(__u);
	const auto __cosh = __eu + _Tp{1} / __eu;
	const auto __sinh = __eu - _Tp{1} / __eu;
        const auto __s = std::exp(_S_pi_4 * __sinh);
	const auto __w = __s + _Tp{1} / __s;
	const auto __x = (__b * __s + __a / __s) / __w;
	if (__x != __a && __x != __b)
	  {
	    const auto __dxdu = __cosh / (__w * __w);
	    __sum += __func(__x) * __dxdu; 
	  }
      }
    // Interlace values (don't go past the rightmost point).
    for (int __k = -__n; __k < __n; ++__k)
      {
	const auto __u = __h * (_Tp(__k) + 0.5);
	const auto __eu = std::exp(__u);
	const auto __cosh = __eu + _Tp{1} / __eu;
	const auto __sinh = __eu - _Tp{1} / __eu;
        const auto __s = std::exp(_S_pi_4 * __sinh);
	const auto __w = __s + _Tp{1} / __s;
	const auto __x = (__b * __s + __a / __s) / __w;
	if (__x != __a && __x != __b)
	  {
	    const auto __dxdu = __cosh / (__w * __w);
	    __sum += __func(__x) * __dxdu; 
	  }
      }
    return __h * (__b - __a) * _S_pi_4 * __sum;
}

/**
 * Progressive version
 */
template<typename _Func, typename _Tp>
  _Tp
  double_exp_integrate4(_Func __func, int __n, _Tp __a, _Tp __b, _Tp /*__tol*/)
  {
    __n /= 2;
    auto __h = _Tp{5} / __n; // 5.0 is rough limit of K in exp(exp(K)).

    const auto _S_pi_4 = __gnu_cxx::__const_pi_quarter<_Tp>();
    auto __sum = __func((__a + __b) / _Tp{2}) / _Tp{2};
    auto __sum1 = _Tp{0};
    auto __sum2 = _Tp{0};
    for (int __k = -__n; __k < 0; ++__k)
      {
	const auto __u = __h * _Tp(__k);
	const auto __eu = std::exp(__u);
	const auto __cosh = __eu + _Tp{1} / __eu;
	const auto __sinh = __eu - _Tp{1} / __eu;
        const auto __s = std::exp(_S_pi_4 * __sinh);
	const auto __w = __s + _Tp{1} / __s;
	const auto __dxdu = __cosh / (__w * __w);
	const auto __x1 = (__b * __s + __a / __s) / __w;
	if (__x1 != __a && __x1 != __b) 
	  __sum1 += __dxdu * __func(__x1);
	const auto __x2 = (__a * __s + __b / __s) / __w;
	if (__x2 != __a && __x2 != __b)
	  __sum2 += __dxdu * __func(__x2);
      }

    for (int __l = 0; __l < 3; ++__l)
      {
	for (int __k = -__n; __k < 0; ++__k)
	  {
	    const auto __u = __h * (_Tp(__k) + 0.5);
	    const auto __eu = std::exp(__u);
	    const auto __cosh = __eu + _Tp{1} / __eu;
	    const auto __sinh = __eu - _Tp{1} / __eu;
            const auto __s = std::exp(_S_pi_4 * __sinh);
	    const auto __w = __s + _Tp{1} / __s;
	    const auto __dxdu = __cosh / (__w * __w);
	    const auto __x1 = (__b * __s + __a / __s) / __w;
	    if (__x1 != __a && __x1 != __b) 
	      __sum1 += __dxdu * __func(__x1);
	    const auto __x2 = (__a * __s + __b / __s) / __w;
	    if (__x2 != __a && __x2 != __b)
	      __sum2 += __dxdu * __func(__x2);
	  }
	__n *= 2;
	__h /= _Tp{2};
      }
    return _Tp{2} * (__b - __a) * _S_pi_4 * __h * (__sum + __sum1 + __sum2);
  }

} // namespace __gnu_test

#endif // DOUBLE_EXP_INTEGRATE_TCC
