
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file ext/root_finding.h
 */

#ifndef ROOT_SEARCH_H
#define ROOT_SEARCH_H 1

#include <vector>
#include <tuple>
#include <limits>

// Default eps = Tp{5} * std::numeric_limits<Tp>::epsilon()?

// template<typename Tp>
//   struct bracket_t {Tp value_lower, Tp value_upper}; ?

// All the root methods should have the same API.
// brent, newton, safe don't have max_iter.

// root_brackets should have an output iterator API?

// in root_newton what if df is zero?

namespace emsr
{

  // This struct contains both the value and the derivative at a point.
  template<typename Tp>
    struct root_state
    {
      Tp value;
      Tp deriv;
    };

  template<typename Tp, typename Func>
    bool
    root_bracket(Func func, Tp& x_lower, Tp& x_upper,
    		   std::size_t max_iter = 50);

  template<typename Tp, typename Func>
    std::vector<std::pair<Tp, Tp>>
    root_brackets(Func func,
		    Tp x_lower, Tp x_upper, std::size_t n);

  template<typename Tp, typename Func>
    Tp
    root_bisect(Func func, Tp x_lower, Tp x_upper, Tp eps,
		  std::size_t max_iter = std::numeric_limits<Tp>::digits);

  template<typename Tp, typename Func>
    Tp
    root_secant(Func func, Tp x_lower, Tp x_upper, Tp eps,
		  std::size_t max_iter = 40);

  template<typename Tp, typename Func>
    Tp
    root_false_position(Func func, Tp x_lower, Tp x_upper,
			  Tp eps, std::size_t max_iter = 40);

  template<typename Tp, typename Func>
    Tp
    root_ridder(Func func, Tp x_lower, Tp x_upper,
		  Tp eps, std::size_t max_iter = 100);

  template<typename Tp, typename Func>
    Tp
    root_brent(Func func, Tp x_lower, Tp x_upper,
		 Tp eps, std::size_t max_iter = 100);

  template<typename Tp, typename StateFunc>
    Tp
    root_newton(StateFunc func, Tp x_lower, Tp x_upper,
    		  Tp eps, std::size_t max_iter = 40);

  template<typename Tp, typename StateFunc>
    Tp
    root_halley(StateFunc func, Tp x_lower, Tp x_upper,
    		  Tp eps, std::size_t max_iter = 40);

  template<typename Tp, typename Func>
    Tp
    root_steffensen(Func func, Tp x_lower, Tp x_upper,
    		      Tp eps, std::size_t max_iter = 40);

  template<typename Tp, typename StateFunc>
    Tp
    root_safe(StateFunc func, Tp x_lower, Tp x_upper,
    		Tp eps, std::size_t max_iter = 100);

} // namespace emsr

#include <emsr/root_search.tcc>

#endif // ROOT_SEARCH_H
