// { dg-require-c-std "" }
// { dg-add-options ieee }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }

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

// ellint_rd

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
  using __gnu_cxx::ellint_rd;

  auto rd1 = ellint_rd(0.0, 2.0, 1.0);
  auto rd2 = ellint_rd(2.0, 3.0, 4.0);
  auto rd3 = ellint_rd(cmplx(0.0, 1.0), cmplx(0.0, -1.0), cmplx(2.0, 0.0));
  auto rd4 = ellint_rd(cmplx(0.0, 0.0), cmplx(0.0, 1.0), cmplx(0.0, -1.0));
  auto rd5 = ellint_rd(cmplx(0.0, 0.0), cmplx(-1.0, 1.0), cmplx(0.0, 1.0));
  auto rd6 = ellint_rd(cmplx(-2.0, -1.0), cmplx(0.0, -1.0), cmplx(-1.0, 1.0));

  bool test [[gnu::unused]] = true;
  double eps = 1.0e-12;
  VERIFY(std::abs(rd1 - 1.7972103521034) < eps);
  VERIFY(std::abs(rd2 - 0.16510527294261) < eps);
  VERIFY(std::abs(rd3 - 0.65933854154220) < eps);
  VERIFY(std::abs(rd4 - cmplx(1.2708196271910, +2.7811120159521)) < eps);
  VERIFY(std::abs(rd5 - cmplx(-1.8577235439239, -0.96193450888839)) < eps);
  VERIFY(std::abs(rd6 - cmplx(1.8249027393704, -1.2218475784827)) < eps);
}

int
main()
{
  test01();
  return 0;
}
