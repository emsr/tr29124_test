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

// ellint_rj

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
  using __gnu_cxx::ellint_rj;

  auto rj1 = ellint_rj(0.0, 1.0, 2.0, 3.0);
  auto rj2 = ellint_rj(2.0, 3.0, 4.0, 5.0);
  auto rj3 = ellint_rj(cmplx(2.0,0.0), cmplx(3.0,0.0),
		       cmplx(4.0,0.0), cmplx(-1.0,1.0));
  auto rj4 = ellint_rj(cmplx(0.0,1.0), cmplx(0.0,-1.0),
		       cmplx(0.0,0.0), cmplx(2.0,0.0));
  auto rj5 = ellint_rj(cmplx(-1.0,1.0), cmplx(-1.0,-1.0),
				  cmplx(1.0,0.0), cmplx(2.0,0.0));
  auto rj6 = ellint_rj(cmplx(0.0,1.0), cmplx(0.0,-1.0),
		       cmplx(0.0,0.0), cmplx(1.0,-1.0));
  auto rj7 = ellint_rj(cmplx(-1.0,1.0), cmplx(-1.0,-1.0),
		       cmplx(1.0,0.0), cmplx(-3.0,1.0));
  auto rj8 = ellint_rj(cmplx(-1.0,1.0), cmplx(-2.0,-1.0),
		       cmplx(0.0,-1.0), cmplx(-1.0,1.0));

  bool test [[gnu::unused]] = true;
  double eps = 1.0e-12;
  VERIFY(std::abs(rj1 - 0.77688623778582) < eps);
  VERIFY(std::abs(rj2 - 0.14297579667157) < eps);
  VERIFY(std::abs(rj3 - cmplx(0.13613945827771, -0.38207561624427)) < eps);
  VERIFY(std::abs(rj4 - 1.6490011662711) < eps);
  VERIFY(std::abs(rj5 - 0.94148358841220) < eps);
  VERIFY(std::abs(rj6 - cmplx(1.8260115229009, +1.2290661908643)) < eps);
  VERIFY(std::abs(rj7 - cmplx(-0.61127970812028, -1.0684038390007)) < eps);
  VERIFY(std::abs(rj8 - cmplx(1.8249027393704, -1.2218475784827)) < eps);
}

int
main()
{
  test01();
  return 0;
}
