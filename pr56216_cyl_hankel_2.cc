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

// PR libstdc++/56216 - Crash of Bessel functions at x==0!

#include <testsuite_hooks.h>
#include <cmath>

void
test01()
{
  const double inf = std::numeric_limits<double>::infinity();
  std::complex<double> h20 = __gnu_cxx::cyl_hankel_2(0.0, 0.0);
  double h20r = std::real(h20);
  double h20i = std::imag(h20);
  std::complex<double> h21 = __gnu_cxx::cyl_hankel_2(1.0, 0.0);
  double h21r = std::real(h21);
  double h21i = std::imag(h21);

  bool test [[gnu::unused]] = true;
  VERIFY(h20r == 1.0);
  VERIFY(h20i == inf);
  VERIFY(h21r == 0.0);
  VERIFY(h21i == inf);
}

int
main()
{
  test01();
  return 0;
}
