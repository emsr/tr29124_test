// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2019 Free Software Foundation, Inc.
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
// Ported from GSL by Jason Dick and Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements an extrapolation table for use in integration schemes
// Based on gsl/integration/qelg.c

#ifndef QUADRATURE_POINT_H
#define QUADRATURE_POINT_H 1

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A structure to store quadrature rules.
   */
  template<typename _Tp>
    struct __quadrature_point_t
    {
      _Tp __point;
      _Tp __weight;

      constexpr __quadrature_point_t() = default;

      constexpr __quadrature_point_t(_Tp __pt, _Tp __wt)
      : __point(__pt),
	__weight(__wt)
      { }
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // QUADRATURE_POINT_H
