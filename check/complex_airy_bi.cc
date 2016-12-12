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

// airy_bi

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
  using __gnu_cxx::airy_bi;

  auto bi0 = airy_bi(cmplx( 1.0000000000000000,  0.0000000000000000));
  auto bi1 = airy_bi(cmplx( 0.8090169943749474,  0.5877852522924731));
  auto bi2 = airy_bi(cmplx( 0.3090169943749474,  0.9510565162951536));
  auto bi3 = airy_bi(cmplx(-0.3090169943749474,  0.9510565162951536));
  auto bi4 = airy_bi(cmplx(-0.8090169943749474,  0.5877852522924731));
  auto bi5 = airy_bi(cmplx(-1.0000000000000000,  0.0000000000000000));
  auto bi6 = airy_bi(cmplx(-0.8090169943749474, -0.5877852522924731));
  auto bi7 = airy_bi(cmplx(-0.3090169943749474, -0.9510565162951536));
  auto bi8 = airy_bi(cmplx( 0.3090169943749474, -0.9510565162951536));
  auto bi9 = airy_bi(cmplx( 0.8090169943749474, -0.5877852522924731));

  bool test [[gnu::unused]] = true;
  double eps = 1.0e-12;
  VERIFY(std::abs(bi0 - cmplx( 1.207423594952871,   0.0000000000000000)) < eps);
  VERIFY(std::abs(bi1 - cmplx( 0.9127160108293936,  0.3800456133135556)) < eps);
  VERIFY(std::abs(bi2 - cmplx( 0.6824453575635721,  0.3343047153635002)) < eps);
  VERIFY(std::abs(bi3 - cmplx( 0.5726265660086474,  0.3988641086982559)) < eps);
  VERIFY(std::abs(bi4 - cmplx( 0.2511841251049547,  0.3401447690712719)) < eps);
  VERIFY(std::abs(bi5 - cmplx( 0.1039973894969446,  0.0000000000000000)) < eps);
  VERIFY(std::abs(bi6 - cmplx( 0.2511841251049547, -0.3401447690712719)) < eps);
  VERIFY(std::abs(bi7 - cmplx( 0.5726265660086474, -0.3988641086982559)) < eps);
  VERIFY(std::abs(bi8 - cmplx( 0.6824453575635721, -0.3343047153635002)) < eps);
  VERIFY(std::abs(bi9 - cmplx( 0.9127160108293936, -0.3800456133135556)) < eps);
}

int
main()
{
  test01();
  return 0;
}
