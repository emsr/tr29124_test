// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
// Copyright (C) 2018 Free Software Foundation, Inc.
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

// PR libstdc++/83140 - assoc_legendre returns negated value when m is odd

#include <cmath>
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
  long double P11_ok = 0.866025403784438646787;

  auto P11 = std::assoc_legendre(1, 1, 0.5);
  auto diff_ok = P11 - P11_ok;
  VERIFY(std::abs(diff_ok) < 1.0e-12);
}

int
main()
{
  test01();
  return 0;
}
