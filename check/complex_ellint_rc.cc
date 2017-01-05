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

// ellint_rc

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
  using __gnu_cxx::ellint_rc;

  auto rc1 = ellint_rc(0.0, 0.25);
  auto rc2 = ellint_rc(2.25, 2.0);  //  ln 2 = 0.69314718055995
  auto rc3 = ellint_rc(cmplx(0.0, 0.0), cmplx(0.0, 1.0));
  //auto rc4 = ellint_rc(cmplx(0.0, -1.0), cmplx(0.0, 1.0));
  auto rc5 = ellint_rc(0.25, -2.0);
  auto rc6 = ellint_rc(cmplx(0.0, 1.0), cmplx(-1.0, 0.0));

  bool test [[gnu::unused]] = true;
  double eps = 1.0e-12;
  VERIFY(std::abs(rc1 - 3.1415926535898) < eps);
  VERIFY(std::abs(rc2 - 0.69314718055995) < eps);
  VERIFY(std::abs(rc3 - cmplx(1.1107207345396, -1.1107207345396)) < eps);
  //VERIFY(std::abs(rc4 - cmplx(1.2260849569072, -0.34471136988768)) < eps);
  VERIFY(std::abs(rc5 - 0.23104906018665) < eps);
  VERIFY(std::abs(rc6 - cmplx(0.77778596920447, 0.19832484993429)) < eps);
}

int
main()
{
  test01();
  return 0;
}
