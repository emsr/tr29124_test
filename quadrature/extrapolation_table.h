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
// Ported from GSL by Jason Dick and Ed Smith-Rowland
// Originally written by Brian Gaugh
//
//This file implements an extrapolation table for use in integration schemes
//Based upon gsl-2.1/integration/qelg.c

#ifndef EXTRAPOLATION_TABLE_H
#define EXTRAPOLATION_TABLE_H 1

#include <array>
#include <utility>
#include <limits>
#include <cmath>
#include <algorithm>

namespace __gnu_test
{

  template<typename _Tp>
    class extrapolation_table
    {

    private:

      size_t _M_nn;
      std::array<_Tp, 52> _M_rlist2;
      size_t _M_nres;
      std::array<_Tp, 3> _M_res3la;

    public:

      extrapolation_table()
      : _M_nn(0),
	_M_nres(0)
      {}

      explicit extrapolation_table(_Tp __y)
      : _M_nn(0),
	_M_nres(0)
      { this->append(__y); }

      void
      append(_Tp __y)
      {
	_M_rlist2[_M_nn] = __y;
	++_M_nn;
      }

      std::pair<_Tp, _Tp>
      qelg();

      size_t
      get_nn() const
      { return _M_nn; }
    };

} // namespace __gnu_test

#include "extrapolation_table.tcc"

#endif // EXTRAPOLATION_TABLE_H
