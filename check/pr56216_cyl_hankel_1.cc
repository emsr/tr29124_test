// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
//
// Copyright (C) 2017 Free Software Foundation, Inc.
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

#include <cmath>
#include <limits>
#if defined(__TEST_DEBUG)
#  include <iostream>
#  define VERIFY(A) \
  if (!(A)) \
    { \
      std::cout << "line " << __LINE__ \
	<< std::endl; \
    }
#else
#  include <testsuite_hooks.h>
#endif

void
test01()
{
  const double inf = std::numeric_limits<double>::infinity();
  std::complex<double> h10 = __gnu_cxx::cyl_hankel_1(0.0, 0.0);
  double h10r = std::real(h10);
  double h10i = std::imag(h10);
  std::complex<double> h11 = __gnu_cxx::cyl_hankel_1(1.0, 0.0);
  double h11r = std::real(h11);
  double h11i = std::imag(h11);

  bool test [[gnu::unused]] = true;
  VERIFY(h10r == 1.0);
  VERIFY(h10i == -inf);
  VERIFY(h11r == 0.0);
  VERIFY(h11i == -inf);
}

int
main()
{
  test01();
  return 0;
}
