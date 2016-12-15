// Copyright (C) 2015-2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
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

// specfun_stats.h

//
// These are little tools for special function output statistics.
//

#ifndef _GLIBCXX_SPECFUN_STATS_H
#define _GLIBCXX_SPECFUN_STATS_H

#include <utility> // For pair.

template<typename Tp>
  struct Stats
  {
    Stats(Tp tol)
    : toler(tol)
    { }

    Stats&
    operator<<(std::pair<Tp, Tp> funs)
    {
      ++num;
      auto f = funs.first;
      auto f0 = funs.second;
      auto diff = f - f0;
      sum_diff += diff;
      sum_diff2 += diff * diff;
      auto abs_diff = std::abs(diff);
      if (abs_diff > max_abs_diff)
	max_abs_diff = abs_diff;
      if (std::abs(f0) > Tp(10) * toler
       && std::abs(f) > Tp(10) * toler)
	{
	  const auto frac = diff / f0;
	  if (std::abs(frac) > max_abs_frac)
	    max_abs_frac = std::abs(frac);
	}
      return *this;
    }

    Tp
    mean_diff() const
    { return sum_diff / num; }

    Tp
    variance_diff() const
    { return sum_diff2 / num - mean_diff() * mean_diff(); }

    Tp toler = 10 * std::numeric_limits<Tp>::epsilon();
    unsigned int num = 0;
    Tp sum_diff = Tp{};
    Tp sum_diff2 = Tp{};
    Tp max_abs_diff = Tp{};
    Tp max_abs_frac = Tp{};
  };

#endif // _GLIBCXX_SPECFUN_STATS_H
