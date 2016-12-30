/* integration/cquad_workspace.h
 *
 * Copyright (C) 2010 Pedro Gonnet
 * Copyright (C) 2016 Free Software Foundation, Inc.
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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Pedro Gonnet
//
// This file contains constants for use in integration schemes.
// Based upon structs in gsl-2.3/integration/gsl_integration.h

#ifndef CQUAD_WORKSPACE_H
#define CQUAD_WORKSPACE_H 1

#include <vector>

namespace __gnu_test
{

  /**
   * Data of a single interval.
   */
  template<typename _Tp>
    struct cquad_interval
    {
      _Tp a, b;
      _Tp c[64];
      _Tp fx[33];
      _Tp igral;
      _Tp err;
      std::size_t depth, rdepth, ndiv;
    };

  /**
   * Comparison of cquad intervals.
   */
  template<typename _Tp>
    bool
    operator<(const cquad_interval<_Tp>& __ivl,
	      const cquad_interval<_Tp>& __ivr)
    { return __ivl.err < __ivr.err; }

  template<typename _Tp>
    struct cquad_interval_comp
    {
      bool
      operator()(const cquad_interval<_Tp>& __ivl,
		 const cquad_interval<_Tp>& __ivr)
      { return __ivr.err < __ivl.err; }
    };

  /**
   * The workspace is just a collection of intervals.
   */
  template<typename _Tp>
    struct cquad_workspace
    {
      std::vector<cquad_interval<_Tp>> ivals;

      cquad_workspace(std::size_t __len = 200)
      : ivals()
      { this->ivals.reserve(__len); }

      std::size_t size() const
      { return this->ivals.size(); }

      typename std::vector<cquad_interval<_Tp>>::iterator
      begin()
      { return this->ivals.begin(); }

      typename std::vector<cquad_interval<_Tp>>::iterator
      end()
      { return this->ivals.end(); }

      const cquad_interval<_Tp>&
      top() const
      { return this->ivals[0]; }

      cquad_interval<_Tp>&
      top()
      { return this->ivals[0]; }

      void
      clear()
      { this->ivals.clear(); }

      void
      push(const cquad_interval<_Tp>& __iv)
      {
	cquad_interval_comp<_Tp> __cmp;
	this->ivals.push_back(__iv);
	std::push_heap(std::begin(this->ivals), std::end(this->ivals), __cmp);
      }

      void
      pop()
      {
	std::pop_heap(std::begin(this->ivals), std::end(this->ivals));
	this->ivals.pop_back();
      }

      void
      update()
      {
	cquad_interval_comp<_Tp> __cmp;
	std::make_heap(std::begin(this->ivals), std::end(this->ivals), __cmp);
      }

      _Tp
      total_integral() const
      {
	auto __tot_igral = _Tp{0};
	for (auto& __iv : ivals)
	  __tot_igral += __iv.igral;
	return __tot_igral;
      }

      _Tp
      total_error() const
      {
	auto __tot_error = _Tp{0};
	for (auto& __iv : ivals)
	  __tot_error += __iv.err;
	return __tot_error;
      }
    };

} // namespace __gnu_test

#endif // CQUAD_WORKSPACE_H
