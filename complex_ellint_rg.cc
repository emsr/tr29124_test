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

// ellint_rg

#include <cmath>
#include <complex>
#include <testsuite_hooks.h>

void
test01()
{
  using cmplx = std::complex<double>;
  using __gnu_cxx::ellint_rg;

  auto rg1 = ellint_rg(0.0, 16.0, 16.0);
  auto rg2 = ellint_rg(2, 3, 4);
  auto rg3 = ellint_rg(cmplx(0.0, 0.0), cmplx(0.0, 1.0), cmplx(0.0, -1.0));
  auto rg4 = ellint_rg(cmplx(-1.0, 1.0), cmplx(0.0, 1.0), cmplx(0.0, 0.0));
  auto rg5 = ellint_rg(cmplx(0.0, -1.0), cmplx(-1.0, 1.0), cmplx(0.0, 1.0));
  auto rg6 = ellint_rg(0, 0.0796, 4);

  bool test [[gnu::unused]] = true;
  double eps = 1.0e-12;
  VERIFY(std::abs(rg1 - 3.1415926535898) < eps);
  VERIFY(std::abs(rg2 - 1.725503020692) < eps);
  VERIFY(std::abs(rg3 - 0.42360654239699) < eps);
  VERIFY(std::abs(rg4 - cmplx(0.44660591677018, +0.70768352357515)) < eps);
  VERIFY(std::abs(rg5 - cmplx(0.36023392184473, +0.40348623401722)) < eps);
  VERIFY(std::abs(rg6 - 1.0284758090288) < eps);
}

int
main()
{
  test01();
  return 0;
}
