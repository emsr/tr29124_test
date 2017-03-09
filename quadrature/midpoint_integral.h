/* quadrature/midpoint_integral.h
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

#ifndef MIDPOINT_INTEGRAL_H
#define MIDPOINT_INTEGRAL_H 1

#include <cstddef>
#include <limits>

namespace __gnu_test
{

template<typename _Func, typename _Tp>
  class midpoint_integral
  {
  public:

    midpoint_integral(_Func __fun, _Tp __a, _Tp __b, _Tp __err)
    : _M_fun(__fun), _M_a(__a), _M_b(__b), _M_err(__err), _M_sum()
    { }

    _Tp operator()();

  private:

    _Tp _M_step();

    _Func _M_fun;
    _Tp _M_a;
    _Tp _M_b;
    _Tp _M_err;
    _Tp _M_sum;
    std::size_t _M_iter = 0;
    std::size_t _M_pow3 = 0;
    std::size_t _S_max_iter = std::numeric_limits<std::size_t>::digits / 2;
  };

} // namespace __gnu_test

#include "midpoint_integral.tcc"

#endif // MIDPOINT_INTEGRAL_H
