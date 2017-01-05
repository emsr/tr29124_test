/* integration/cquad_workspace.h
 *
 * Copyright (C) 2010 Pedro Gonnet
 * Copyright (C) 2016-2017 Free Software Foundation, Inc.
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
      _Tp igral, err;
      std::size_t depth, rdepth, ndiv;
    };

  /**
   * The workspace is just a collection of intervals.
   */
  template<typename _Tp>
    struct cquad_workspace
    {
      std::vector<cquad_interval<_Tp>> ivals;
      std::vector<std::size_t> heap;

      std::size_t size() const
      { return ivals.size(); }

      cquad_workspace(std::size_t __n)
      : ivals(__n), heap(__n)
      {
	if (__n < 3)
	  std::__throw_domain_error("cquad_workspace: "
				    "Workspace size n must be at least 3");

	this->initialize_heap();
      }

      void
      initialize_heap()
      {
	for (std::size_t __i = 0; __i < this->heap.size(); ++__i)
	  this->heap[__i] = __i;
      }
    };

} // namespace __gnu_test

#endif // CQUAD_WORKSPACE_H
