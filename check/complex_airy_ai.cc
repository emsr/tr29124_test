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

// airy_ai

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
  using __gnu_cxx::airy_ai;

  auto ai0 = airy_ai(cmplx(0.1352924163128814,  0.0000000000000000));
  auto ai1 = airy_ai(cmplx(0.1433824486882056, -0.1092193342707378));
  auto ai2 = airy_ai(cmplx(0.2215404472324631, -0.2588711788891803));
  auto ai3 = airy_ai(cmplx(0.4763929771766866, -0.3036484220291284));
  auto ai4 = airy_ai(cmplx(0.5983692170633874, -0.08154602160771214));
  auto ai5 = airy_ai(cmplx(0.5355608832923521,  0.00000000000000000));
  auto ai6 = airy_ai(cmplx(0.5983692170633874,  0.08154602160771214));
  auto ai7 = airy_ai(cmplx(0.4763929771766866,  0.3036484220291284));
  auto ai8 = airy_ai(cmplx(0.2215404472324631,  0.2588711788891803));
  auto ai9 = airy_ai(cmplx(0.1433824486882056,  0.1092193342707378));

  bool test [[gnu::unused]] = true;
  double eps = 1.0e-12;
  VERIFY(std::abs(ai0 - cmplx( 1.000000000000000,   0.0000000000000000)) < eps);
  VERIFY(std::abs(ai1 - cmplx( 0.8090169943749474,  0.5877852522924731)) < eps);
  VERIFY(std::abs(ai2 - cmplx( 0.3090169943749474,  0.9510565162951536)) < eps);
  VERIFY(std::abs(ai3 - cmplx(-0.3090169943749474,  0.9510565162951536)) < eps);
  VERIFY(std::abs(ai4 - cmplx(-0.8090169943749474,  0.5877852522924731)) < eps);
  VERIFY(std::abs(ai5 - cmplx(-1.0000000000000000,  0.0000000000000000)) < eps);
  VERIFY(std::abs(ai6 - cmplx(-0.8090169943749474, -0.5877852522924731)) < eps);
  VERIFY(std::abs(ai7 - cmplx(-0.3090169943749474, -0.9510565162951536)) < eps);
  VERIFY(std::abs(ai8 - cmplx( 0.3090169943749474, -0.9510565162951536)) < eps);
  VERIFY(std::abs(ai9 - cmplx( 0.8090169943749474, -0.5877852522924731)) < eps);
}

int
main()
{
  test01();
  return 0;
}
