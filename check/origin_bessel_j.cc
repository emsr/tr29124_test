// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
//
// Copyright (C) 2016 Free Software Foundation, Inc.
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

#include <testsuite_hooks.h>
#include <cmath>
#include <iostream>

void
test01()
{
  const double inf = std::numeric_limits<double>::infinity();
  double jm1o2 = std::cyl_bessel_j(-0.5, 0.0);
  double jm1 = std::cyl_bessel_j(-1.0, 0.0);
  double jm3o2 = std::cyl_bessel_j(-1.5, 0.0);
  double jm2 = std::cyl_bessel_j(-2.0, 0.0);

  bool test [[gnu::unused]] = true;
  VERIFY(jm1o2 == inf);
  VERIFY(jm1 == 0.0);
  VERIFY(jm3o2 == -inf);
  VERIFY(jm2 == 0.0);
}

int
main()
{
  test01();

  return 0;
}
