// { dg-require-c-std "" }
// { dg-add-options ieee }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }

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

// ellint_rf

#include <cmath>
#include <complex>
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
  using cmplx = std::complex<double>;
  using __gnu_cxx::ellint_rf;

  auto rf1 = ellint_rf(1.0, 2.0, 0.0);
  auto rf2 = ellint_rf(cmplx(0.0, 1.0), cmplx(0.0, -1.0), cmplx(0.0, 0.0));
  auto rf3 = ellint_rf(cmplx(-1.0, 1.0), cmplx(0.0, 1.0), cmplx(0.0, 0.0));
  auto rf4 = ellint_rf(0.5, 1.0, 0.0);
  auto rf5 = ellint_rf(2.0, 3.0, 4.0);
  auto rf6 = ellint_rf(cmplx(0.0, 1.0), cmplx(0.0, -1.0), cmplx(2.0, 0.0));
  auto rf7 = ellint_rf(cmplx(-1.0, 1.0), cmplx(1.0, -1.0), cmplx(0.0, 1.0));

  bool test [[gnu::unused]] = true;
  double eps = 1.0e-12;
  VERIFY(std::abs(rf1 - 1.3110287771461) < eps);
  VERIFY(std::abs(rf2 - 1.8540746773014) < eps);
  VERIFY(std::abs(rf3 - cmplx(0.79612586584234, -1.2138566698365)) < eps);
  VERIFY(std::abs(rf4 - 1.8540746773014) < eps);
  VERIFY(std::abs(rf5 - 0.58408284167715) < eps);
  VERIFY(std::abs(rf6 - 1.0441445654064) < eps);
  VERIFY(std::abs(rf7 - cmplx(0.93912050218619, -0.53296252018635)) < eps);
}

int
main()
{
  test01();
  return 0;
}
