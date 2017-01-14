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

	bool
	operator<(const interval& __iv)
	{ return this->_M_abs_error < __iv._M_abs_error; }
      };

      /**
       * Comparison of quadrature intervals.
       */
      struct interval_comp
      {
	bool
	operator()(const interval& __ivl,
		   const interval& __ivr)
	{ return __ivl._M_abs_error < __ivr._M_abs_error; }
      };

      std::size_t _M_maximum_level;
      std::vector<interval> _M_ival;

    public:

      integration_workspace(std::size_t __cap)
      : _M_maximum_level{0},
	_M_ival{}
      {
	_M_ival.reserve(__cap);
      }

      void sort_error();

      void
      clear()
      { this->_M_ival.clear(); }

      void append(_Tp __a, _Tp __b, _Tp __area, _Tp __error);

      void split(_Tp __ab,
		 _Tp __area1, _Tp __error1,
		 _Tp __area2, _Tp __error2);

      void
      retrieve(_Tp& __lolim, _Tp& __uplim, _Tp& __res, _Tp& __err) const
      {
	__lolim = this->_M_ival[0]._M_lower_lim;
	__uplim = this->_M_ival[0]._M_upper_lim;
	__res = this->_M_ival[0]._M_result;
	__err = this->_M_ival[0]._M_abs_error;
      }

      std::size_t
      size() const
      { return this->_M_ival.size(); }

      std::size_t
      capacity() const
      { return this->_M_ival.capacity(); }

      _Tp
      lower_lim(std::size_t __ii = 0) const
      { return this->_M_ival[__ii]._M_lower_lim; }

      _Tp
      upper_lim(std::size_t __ii = 0) const
      { return this->_M_ival[__ii]._M_upper_lim; }

      _Tp
      result(std::size_t __ii = 0) const
      { return this->_M_ival[__ii]._M_result; }

      _Tp
      abs_error(std::size_t __ii = 0) const
      { return this->_M_ival[__ii]._M_abs_error; }

      size_t
      level(std::size_t ii = 0) const
      { return this->_M_ival[ii]._M_level; }

      // Only used by qagp
      _Tp
      set_abs_error(std::size_t __ii, _Tp __abserr)
      { return this->_M_ival[__ii]._M_abs_error = __abserr; }

      // Only used by qagp
      void
      set_level(std::size_t __ii, std::size_t __lvl)
      { this->_M_ival[__ii]._M_level = __lvl; }

      std::size_t
      current_level() const
      { return this->_M_ival[0]._M_level; }

      std::size_t
      max_level() const
      { return this->_M_maximum_level; }

      bool
      large_interval() const
      {
	if (this->current_level() < this->max_level())
	  return true;
	else
	  return false;
      }

      _Tp
      sum_results() const
      {
	auto __result_sum = _Tp{0};
	for (const auto& __iv : this->_M_ival)
	  __result_sum += __iv._M_result;
	return __result_sum;
      }

      /*
       * FIXME: I think this should be - 
       *  abs(a2 - a1) < tol || abs(b2 - a2) < tol
       */
      static bool
      subinterval_too_small(_Tp __a1, _Tp __a2, _Tp __b2)
      {
	const auto _S_eps = _Tp{100} * std::numeric_limits<_Tp>::epsilon();
	const auto _S_min = _Tp{1000} * std::numeric_limits<_Tp>::min();

	const auto __tmp = (_Tp{1} + _S_eps) * (std::abs(__a2) + _S_min);

	return std::abs(__a1) <= __tmp
	    && std::abs(__b2) <= __tmp;
      }
    };

} // namespace __gnu_test

#include "integration_workspace.tcc"

#endif // INTEGRATION_WORKSPACE_H
