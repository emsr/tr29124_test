/* quadrature/trapezoid_integral.tcc
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

#ifndef TRAPEZOID_INTEGRAL_TCC
#define TRAPEZOID_INTEGRAL_TCC 1

#include <cmath>

namespace __gnu_test
{

/**
 * 
 */
template<typename _Func, typename _Tp>
  _Tp
  trapezoid_integral<_Func, _Tp>::operator()()
  {
    auto __sum_prev = this->_M_step();
    for (std::size_t __j = 1; __j < _S_max_iter; ++__j)
      {
	auto __sum = this->_M_step();
	if (std::abs(__sum - __sum_prev) < this->_M_err * std::abs(__sum))
	  return __sum;
	if (std::abs(__sum) < this->_M_err
		&& std::abs(__sum_prev) < this->_M_err
		&& __j > 6)
	  return __sum;
	__sum_prev = __sum;
      }
    return __sum_prev;
  }

/**
 * Chances are, we stepped on a pole.
 */
template<typename _Func, typename _Tp>
  _Tp
  __downdate_func(_Func __func, _Tp __x)
  {
    auto __y = __func(__x);
    if (__isnan(__y) || std::__detail::__isinf(__y))
      return _Tp{0};
    else
      return __y;
  }

/**
 * 
 */
template<typename _Func, typename _Tp>
  _Tp
  trapezoid_integral<_Func, _Tp>::_M_step()
  {
    if (this->_M_iter == 0)
      {
	this->_M_iter = 1;
        this->_M_sum = (this->_M_upper_lim - this->_M_lower_lim)
		     * (__downdate_func(this->_M_fun, this->_M_lower_lim)
		      + __downdate_func(this->_M_fun, this->_M_upper_lim))
		     / _Tp{2};
        this->_M_pow2 = 1;
      }
    else
      {
	++this->_M_iter;
        const auto __del = (this->_M_upper_lim - this->_M_lower_lim)
			 / this->_M_pow2;
	if (std::abs(__del) < _S_min_delta)
	  return this->_M_sum;
        auto __x = this->_M_lower_lim + __del / _Tp{2};
	auto __sum = _Tp{0};
        for (std::size_t __j = 0; __j < this->_M_pow2; ++__j, __x += __del)
	  __sum += __downdate_func(this->_M_fun, __x);
        this->_M_sum = (this->_M_sum + __del * __sum) / _Tp{2};
        this->_M_pow2 *= 2;
      }
    return this->_M_sum;
  }

} // namespace __gnu_test

#endif // TRAPEZOID_INTEGRAL_TCC
