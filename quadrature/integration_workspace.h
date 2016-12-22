// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.
//
// Ported from GSL by Jason Dick
// Originally written by Brian Gaugh

//
// Implements the integration_workspace class which stores temporary data
// for performing integrals
// Based upon gsl-1.9/integration/workspace.c

#ifndef INTEGRATION_WORKSPACE_H
#define INTEGRATION_WORKSPACE_H 1

#include <vector>
#include <limits>
#include <cmath>

namespace __gnu_test
{

  template<typename _Tp>
    class integration_workspace
    {
    private:

      std::size_t _M_capacity;
      std::size_t _M_size;
      std::size_t _M_nrmax;
      std::size_t _M_ii;
      std::size_t _M_maximum_level;
      std::vector<_Tp> _M_lower_lim, _M_upper_lim, _M_result, _M_abs_error;
      std::vector<std::size_t> _M_order, _M_level;

    public:

      integration_workspace(std::size_t __cap)
      : _M_capacity(__cap),
	_M_size(0),
	_M_nrmax(0),
	_M_ii(0),
	_M_maximum_level(0),
	_M_lower_lim(__cap),
	_M_upper_lim(__cap),
	_M_result(__cap),
	_M_abs_error(__cap),
	_M_order(__cap),
	_M_level(__cap)
      { }

      void
      set_initial(_Tp __a0, _Tp __b0, _Tp __result0, _Tp __error0)
      {
	_M_lower_lim[0] = __a0;
	_M_upper_lim[0] = __b0;
	_M_result[0] = __result0;
	_M_abs_error[0] = __error0;
	_M_order[0] = 0;
	_M_level[0] = 0;
	_M_size = 1;
      }

      void sort_error();

      void append(_Tp __a, _Tp __b, _Tp __area, _Tp __error);

      void update(_Tp __a1, _Tp __b1, _Tp __area1, _Tp __error1,
		  _Tp __a2, _Tp __b2, _Tp __area2, _Tp __error2);

      void
      retrieve(_Tp& __lolim, _Tp& __uplim, _Tp& __res, _Tp& __err) const
      {
	__lolim = _M_lower_lim[_M_ii];
	__uplim = _M_upper_lim[_M_ii];
	__res = _M_result[_M_ii];
	__err = _M_abs_error[_M_ii];
      }

      size_t
      size() const
      { return _M_size; }

      size_t
      capacity() const
      { return _M_capacity; }

      _Tp
      lower_lim(size_t __ii) const
      { return _M_lower_lim[__ii]; }

      _Tp
      upper_lim(size_t __ii) const
      { return _M_upper_lim[__ii]; }

      _Tp
      result(size_t __ii) const
      { return _M_result[__ii]; }

      _Tp
      abs_error(size_t __ii) const
      { return _M_abs_error[__ii]; }

      size_t
      order(size_t __ii) const
      { return _M_order[__ii]; }

      size_t
      level(size_t ii) const
      { return _M_level[ii]; }

      _Tp
      set_abs_error(size_t __ii, _Tp __abserr)
      { return _M_abs_error[__ii] = __abserr; }

      void
      set_level(size_t __ii, size_t __lvl)
      { _M_level[__ii] = __lvl; }

      bool
      increase_nrmax()
      {
	int __k;
	int __id = _M_nrmax;
	int __jupbnd;

	std::size_t __last = _M_size - 1;

	if (__last > (1 + _M_capacity / 2))
	  __jupbnd = _M_capacity + 1 - __last;
	else
	  __jupbnd = __last;

	for (__k = __id; __k <= __jupbnd; ++__k)
	  {
	    std::size_t __i_max = _M_order[_M_nrmax];

	    _M_ii = __i_max;

	    if (_M_level[__i_max] < _M_maximum_level)
	      return true;

	    ++_M_nrmax;

	  }
	return false;
      }

      void
      reset_nrmax()
      {
	_M_nrmax = 0;
	_M_ii = _M_order[0];
      }

      std::size_t
      get_nrmax() const
      {
	return _M_nrmax;
      }

      void
      set_nrmax(std::size_t NRMAX)
      {
	_M_nrmax = NRMAX;
      }

      _Tp
      sum_results() const
      {
	_Tp __result_sum = 0;

	for (std::size_t __kk = 0; __kk < _M_size; ++__kk)
	  __result_sum += _M_result[__kk];

	return __result_sum;
      }

      static bool
      subinterval_too_small(_Tp __a1, _Tp __a2, _Tp __b2)
      {
	const _Tp __e = std::numeric_limits<_Tp>::epsilon();
	const _Tp __u = std::numeric_limits<_Tp>::min();

	_Tp __tmp = (1 + 100 * __e) * (std::abs(__a2) + 1000 * __u);

	return std::abs(__a1) <= __tmp && std::abs(__b2) <= __tmp;
      }

      std::size_t
      current_level() const
      {
	return _M_level[_M_ii];
      }

      std::size_t
      max_level() const
      {
	return _M_maximum_level;
      }

      bool
      large_interval() const
      {
	if (_M_level[_M_ii] < _M_maximum_level)
	  return true;
	else
	  return false;
      }
    };

} // namespace __gnu_test

#include "integration_workspace.tcc"

#endif // INTEGRATION_WORKSPACE_H
