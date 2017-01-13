// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2017 Free Software Foundation, Inc.
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
// Implements the integration_workspace class which stores temporary data
// for performing integrals
// Based on gsl/integration/workspace.c

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

      struct interval
      {
	_Tp _M_lower_lim;
	_Tp _M_upper_lim;
	_Tp _M_result;
	_Tp _M_abs_error;
	std::size_t _M_level;
      };

      std::size_t _M_capacity;
      std::size_t _M_size;
      std::size_t _M_nrmax;
      std::size_t _M_current_index;
      std::size_t _M_maximum_level;
      std::vector<_Tp> _M_lower_lim;
      std::vector<_Tp> _M_upper_lim;
      std::vector<_Tp> _M_result;
      std::vector<_Tp> _M_abs_error;
      std::vector<std::size_t> _M_order;
      std::vector<std::size_t> _M_level;
      bool _M_try_resize;

    public:

      integration_workspace(std::size_t __cap)
      : _M_capacity(__cap),
	_M_size(0),
	_M_nrmax(0),
	_M_current_index(0),
	_M_maximum_level(0),
	_M_lower_lim(__cap),
	_M_upper_lim(__cap),
	_M_result(__cap),
	_M_abs_error(__cap),
	_M_order(__cap),
	_M_level(__cap)
      { }

      void
      resize()
      {
	const auto __new_cap = std::size_t(1 + 1.5 * this->capacity());
	this->_M_capacity = __new_cap;
	this->_M_lower_lim.resize(__new_cap);
	this->_M_upper_lim.resize(__new_cap);
	this->_M_result.resize(__new_cap);
	this->_M_abs_error.resize(__new_cap);
	this->_M_order.resize(__new_cap);
	this->_M_level.resize(__new_cap);
      }

      void
      set_initial_limits(_Tp a0, _Tp b0)
      {
	this->_M_size = 0;
	this->_M_nrmax = 0;
	this->_M_current_index = 0;
	if (this->_M_capacity > 0)
	  {
	    this->_M_lower_lim[0] = a0;
	    this->_M_upper_lim[0] = b0;
	    this->_M_result[0] = _Tp{0};
	    this->_M_abs_error[0] = _Tp{0};
	    this->_M_order[0] = 0;
	    this->_M_level[0] = 0;
	  }
	this->_M_maximum_level = 0;
      }

      void
      set_initial_results(_Tp result, _Tp error)
      {
	this->_M_size = 1;
	if (this->_M_capacity > 0)
	  {
	    this->_M_result[0] = result;
	    this->_M_abs_error[0] = error;
	    this->_M_order[0] = 0;
	    this->_M_level[0] = 0;
	  }
      }

      void
      set_initial(_Tp __a0, _Tp __b0, _Tp __result0, _Tp __error0)
      {
	this->_M_size = 1;
	this->_M_nrmax = 0;
	this->_M_current_index = 0;
	this->_M_lower_lim[0] = __a0;
	this->_M_upper_lim[0] = __b0;
	this->_M_result[0] = __result0;
	this->_M_abs_error[0] = __error0;
	this->_M_order[0] = 0;
	this->_M_level[0] = 0;
	this->_M_maximum_level = 0;
      }

      void sort_error();
      void sort_results();

      void append(_Tp __a, _Tp __b, _Tp __area, _Tp __error);

      void update(_Tp __a1, _Tp __b1, _Tp __area1, _Tp __error1,
		  _Tp __a2, _Tp __b2, _Tp __area2, _Tp __error2);

      void
      retrieve(_Tp& __lolim, _Tp& __uplim, _Tp& __res, _Tp& __err) const
      {
	const auto __ii = this->_M_current_index;
	__lolim = this->_M_lower_lim[__ii];
	__uplim = this->_M_upper_lim[__ii];
	__res = this->_M_result[__ii];
	__err = this->_M_abs_error[__ii];
      }

      std::size_t
      size() const
      { return this->_M_size; }

      std::size_t
      capacity() const
      { return this->_M_capacity; }

      std::size_t
      current_index() const
      { return this->_M_current_index; }

      _Tp
      lower_lim(std::size_t __ii) const
      { return this->_M_lower_lim[__ii]; }

      _Tp
      upper_lim(std::size_t __ii) const
      { return this->_M_upper_lim[__ii]; }

      _Tp
      result(std::size_t __ii) const
      { return this->_M_result[__ii]; }

      _Tp
      abs_error(std::size_t __ii) const
      { return this->_M_abs_error[__ii]; }

      size_t
      order(std::size_t __ii) const
      { return this->_M_order[__ii]; }

      size_t
      level(std::size_t ii) const
      { return this->_M_level[ii]; }

      _Tp
      set_abs_error(std::size_t __ii, _Tp __abserr)
      { return this->_M_abs_error[__ii] = __abserr; }

      void
      set_level(std::size_t __ii, std::size_t __lvl)
      { this->_M_level[__ii] = __lvl; }

      std::size_t
      current_level() const
      { return this->_M_level[this->_M_current_index]; }

      std::size_t
      max_level() const
      { return this->_M_maximum_level; }

      bool
      large_interval() const
      {
	if (this->_M_level[this->_M_current_index] < this->_M_maximum_level)
	  return true;
	else
	  return false;
      }

      bool
      increase_nrmax()
      {
	int __id = this->_M_nrmax;
	int __jupbnd;

	auto __last = this->_M_size - 1;

	if (__last > (1 + this->_M_capacity / 2))
	  __jupbnd = this->_M_capacity + 1 - __last;
	else
	  __jupbnd = __last;

	for (int __k = __id; __k <= __jupbnd; ++__k)
	  {
	    auto __i_max = this->_M_order[this->_M_nrmax];

	    this->_M_current_index = __i_max;

	    if (this->_M_level[__i_max] < this->_M_maximum_level)
	      return true;

	    ++this->_M_nrmax;
	  }
	return false;
      }

      void
      reset_nrmax()
      {
	this->_M_nrmax = 0;
	this->_M_current_index = this->_M_order[0];
      }

      std::size_t
      get_nrmax() const
      { return this->_M_nrmax; }

      void
      set_nrmax(std::size_t nrmax)
      { this->_M_nrmax = nrmax; }

      _Tp
      sum_results() const
      {
	auto __result_sum = _Tp{0};
	for (std::size_t __kk = 0; __kk < this->_M_size; ++__kk)
	  __result_sum += this->_M_result[__kk];

	return __result_sum;
      }

      static bool
      subinterval_too_small(_Tp __a1, _Tp __a2, _Tp __b2)
      {
	const auto __e = std::numeric_limits<_Tp>::epsilon();
	const auto __u = std::numeric_limits<_Tp>::min();

	_Tp __tmp = (1 + 100 * __e) * (std::abs(__a2) + 1000 * __u);

	return std::abs(__a1) <= __tmp && std::abs(__b2) <= __tmp;
      }
    };

} // namespace __gnu_test

#include "integration_workspace.tcc"

#endif // INTEGRATION_WORKSPACE_H
