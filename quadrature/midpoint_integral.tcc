/* quadrature/midpoint_integral.tcc
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

#ifndef MIDPOINT_INTEGRAL_TCC
#define MIDPOINT_INTEGRAL_TCC 1

#include <cmath>

namespace __gnu_ext
{

/**
 * Integrate the function by naive subdivision.
 */
template<typename _Func, typename _Tp>
  typename midpoint_integral<_Func, _Tp>::_AreaTp
  midpoint_integral<_Func, _Tp>::operator()()
  {
    auto __sum_prev = this->_M_step();
    for (std::size_t __j = 1; __j < _S_max_iter; ++__j)
      {
	const auto __sum = this->_M_step();
	this->_M_abs_error = std::abs(__sum - __sum_prev);
	if (this->_M_abs_error < this->_M_rel_tol * std::abs(__sum))
	  return __sum;
	if (__j > 6
	    && std::abs(__sum) < this->_M_rel_tol
	    && std::abs(__sum_prev) < this->_M_rel_tol )
	  return __sum;
	__sum_prev = __sum;
      }
    return __sum_prev;
  }

/**
 * 
 */
template<typename _Func, typename _Tp>
  typename midpoint_integral<_Func, _Tp>::_AreaTp
  midpoint_integral<_Func, _Tp>::_M_step()
  {
    if (this->_M_iter == 0)
      {
	this->_M_iter = 1;
        this->_M_pow3 = 1;
        auto __x = (this->_M_lower_lim + this->_M_upper_lim) / _Tp{2};
        this->_M_result = (this->_M_upper_lim - this->_M_lower_lim)
			* this->_M_fun(__x);
      }
    else
      {
	++this->_M_iter;
        const auto __del = (this->_M_upper_lim - this->_M_lower_lim)
			 / _Tp(3 * this->_M_pow3);
	if (std::abs(__del) < _S_min_delta)
	  return this->_M_result;
        const auto __ddel = _Tp{2} * __del;
        auto __x = this->_M_lower_lim + __del / _Tp{2};
        auto __sum = _Tp{0};
        for (auto __j = 1u; __j <= this->_M_pow3; ++__j)
	  {
            __sum += this->_M_fun(__x);
            __x += __ddel;
            __sum += this->_M_fun(__x);
            __x += __del;
          }
	this->_M_result += (this->_M_upper_lim - this->_M_lower_lim) * __sum
			 / this->_M_pow3;
        this->_M_result /= _Tp{3};
        this->_M_pow3 *= 3;
      }
    return this->_M_result;
  }

} // namespace __gnu_ext

#endif // MIDPOINT_INTEGRAL_TCC
