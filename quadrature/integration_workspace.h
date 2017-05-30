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
	std::size_t _M_depth;

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

      // The start of the heap.
      // This allows to skip the actual max error.
      std::size_t _M_start;

      // The current maximum depth.
      std::size_t _M_max_depth;

      std::vector<interval> _M_ival;

    public:

      integration_workspace(std::size_t __cap)
      : _M_start{0},
	_M_max_depth{0},
	_M_ival{}
      {
	_M_ival.reserve(__cap);
      }

      void sort_error();

      void append(_Tp __a, _Tp __b, _Tp __area, _Tp __error,
		  std::size_t __depth = 0);

      void split(_Tp __ab,
		 _Tp __area1, _Tp __error1,
		 _Tp __area2, _Tp __error2);

      void
      retrieve(_Tp& __lolim, _Tp& __uplim, _Tp& __res, _Tp& __err) const
      {
	__lolim = this->_M_ival[this->_M_start]._M_lower_lim;
	__uplim = this->_M_ival[this->_M_start]._M_upper_lim;
	__res = this->_M_ival[this->_M_start]._M_result;
	__err = this->_M_ival[this->_M_start]._M_abs_error;
      }

      std::size_t
      size() const
      { return this->_M_ival.size(); }

      std::size_t
      capacity() const
      { return this->_M_ival.capacity(); }

      typename std::vector<interval>::iterator
      begin()
      { return this->_M_ival.begin() + this->_M_start; }

      typename std::vector<interval>::iterator
      end()
      { return this->_M_ival.end(); }

      const interval&
      top() const
      { return this->_M_ival[this->_M_start]; }

      interval&
      top()
      { return this->_M_ival[this->_M_start]; }

      void
      clear()
      {
	this->_M_start = 0;
	this->_M_max_depth = 0;
	this->_M_ival.clear();
      }

      void
      push(const interval& __iv)
      {
	this->_M_ival.push_back(__iv);
	std::push_heap(this->begin(), this->end(),
		       interval_comp{});
      }

      void
      pop()
      {
	std::pop_heap(this->begin(), this->end());
	this->_M_ival.pop_back();
      }

      _Tp
      lower_lim(std::size_t __ii = 0) const
      { return this->_M_ival[this->_M_start + __ii]._M_lower_lim; }

      _Tp
      upper_lim(std::size_t __ii = 0) const
      { return this->_M_ival[this->_M_start + __ii]._M_upper_lim; }

      _Tp
      result(std::size_t __ii = 0) const
      { return this->_M_ival[this->_M_start + __ii]._M_result; }

      _Tp
      abs_error(std::size_t __ii = 0) const
      { return this->_M_ival[this->_M_start + __ii]._M_abs_error; }

      size_t
      depth(std::size_t __ii = 0) const
      { return this->_M_ival[this->_M_start + __ii]._M_depth; }

      // Only used by qagp
      _Tp
      set_abs_error(std::size_t __ii, _Tp __abserr)
      { return this->_M_ival[this->_M_start + __ii]._M_abs_error = __abserr; }

      // Only used by qagp
      void
      set_depth(std::size_t __ii, std::size_t __d)
      { this->_M_ival[this->_M_start + __ii]._M_depth = __d; }

      void
      set_start(std::size_t __start)
      {
	this->_M_start = __start;
	this->sort_error();
      }

      bool increment_start();

      void
      reset_start()
      {
	this->_M_start = 0;
	this->sort_error();
      }

      std::size_t
      current_depth() const
      { return this->_M_ival[this->_M_start]._M_depth; }

      std::size_t
      max_depth() const
      { return this->_M_max_depth; }

      std::size_t
      start() const
      { return this->_M_start; }

      bool
      large_interval() const
      {
	if (this->current_depth() < this->max_depth())
	  return true;
	else
	  return false;
      }

      /// Return the total integral:
      /// the sum of the results over all integration segments.
      _Tp
      total_integral() const
      {
	auto __result_sum = _Tp{0};
	for (const auto& __iv : this->_M_ival)
	  __result_sum += __iv._M_result;
	return __result_sum;
      }

      /// Return the sum of the absolute errors over all integration segments.
      _Tp
      total_error() const
      {
	auto __tot_error = _Tp{0};
	for (auto& __iv : this->_M_ival)
	  __tot_error += __iv._M_abs_error;
	return __tot_error;
      }

      /// Return the vector of integration intervals.
      const std::vector<interval>&
      intervals() const
      { return this->_M_ival; }

      /*
       * 
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

  template<typename _Tp>
    std::ostream&
    operator<<(std::ostream& __out, const integration_workspace<_Tp>& __ws)
    {
      auto __w = __out.width();
      __out << std::setw(0);
      __out << ' ' << std::setw(2) << __ws.max_depth() << '\n';
      __out << ' ' << std::setw(2) << __ws.start() << '\n';
      for (const auto& __seg : __ws.intervals())
	__out << ' ' << std::setw(2) << __seg._M_depth
	      << ' ' << std::setw(__w) << __seg._M_lower_lim
	      << ' ' << std::setw(__w) << __seg._M_upper_lim
	      << ' ' << std::setw(__w) << __seg._M_result
	      << ' ' << std::setw(__w) << __seg._M_abs_error
	      << '\n';
      return __out;
    }

} // namespace __gnu_test

#include "integration_workspace.tcc"

#endif // INTEGRATION_WORKSPACE_H
