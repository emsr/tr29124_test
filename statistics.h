// Copyright (C) 2015-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

// statistics.h

#ifndef STATISTICS_H
#define STATISTICS_H 1

#include <emsr/complex_util.h>

/**
 * Incremental computation of statistics.
 *
 * @todo Make an insertion operator (operator<<) that takes an object
 * decomposable into two parts. The first part will be the value
 * and the second part will be the weight.  Then _M_count or something
 * will need to be of type _Tp rather than unsigned integral.
 */
template<typename _Tp>
  struct _Statistics
  {
    _Statistics&
    operator<<(_Tp __diff)
    {
      ++_M_count;
      auto __old_mean = _M_mean;
      _M_mean = (_M_type(__diff) + _M_type(_M_count - 1) * _M_mean)
	      / _M_type(_M_count);
      auto __del_mean = _M_mean - __old_mean;
      auto __del_diff = _M_type(__diff) - _M_mean;
      if (_M_count > 1)
	_M_variance = (_M_type(_M_count - 2) * _M_variance * _M_variance
		    + _M_type(_M_count - 1) * __del_mean * __del_mean
		    + __del_diff * __del_diff) / _M_type(_M_count - 1);
      if (__diff < _M_min)
	{
	  _M_min = __diff;
	  _M_min_index = _M_count - 1;
	}
      if (__diff > _M_max)
	{
	  _M_max = __diff;
	  _M_max_index = _M_count - 1;
	}

      return *this;
    }

    static constexpr bool _M_is_complex = emsr::is_complex_v<_Tp>;

    using _M_type = std::conditional_t<emsr::is_complex_v<_Tp>,
				       std::complex<long double>, long double>;

    _Tp
    count() const
    { return _Tp(_M_count); }

    _Tp
    mean() const
    { return _Tp(_M_mean); }

    _Tp
    variance() const
    { return _Tp(_M_variance); }

    _Tp
    std_deviation() const
    { return _Tp(std::sqrt(_M_variance)); }

    _Tp
    min() const
    { return _Tp(_M_min); }

    std::size_t
    min_index() const
    { return _M_min_index; }

    _Tp
    max() const
    { return _Tp(_M_max); }

    std::size_t
    max_index() const
    { return _M_max_index; }

    std::size_t _M_count = 0;
    std::size_t _M_min_index = -1;
    std::size_t _M_max_index = -1;
    _M_type _M_mean = 0;
    _M_type _M_variance = 0;
    _M_type _M_min = std::numeric_limits<long double>::max();
    _M_type _M_max = std::numeric_limits<long double>::lowest();
  };

#endif // STATISTICS_H
