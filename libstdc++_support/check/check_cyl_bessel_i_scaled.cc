// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
// { dg-skip-if "no extensions in strict dialects" { *-*-* } { "-std=c++*" } }

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

//  cyl_bessel_i_scaled
//  Compare against values generated by the GNU Scientific Library.
//  The GSL can be found on the web: http://www.gnu.org/software/gsl/
#include <limits>
#include <cmath>
#if defined(__TEST_DEBUG)
#  include <iostream>
#  define VERIFY(A) \
  if (!(A)) \
    { \
      std::cout << "line " << __LINE__ \
	<< "  max_abs_frac = " << max_abs_frac \
	<< '\n'; \
    }
#else
#  include <testsuite_hooks.h>
#endif
#include <specfun_testcase.h>

// Test data for nu=100.00000000000000.
// max(|f - f_GSL|): 1.8149544367407344e-16 at index 9
// max(|f - f_GSL| / |f_GSL|): 1.2043137544512246e-12
// mean(f - f_GSL): -1.9201466884302938e-17
// variance(f - f_GSL): 2.5502292518858383e-34
// stddev(f - f_GSL): 1.5969437222037095e-17
const testcase_cyl_bessel_i_scaled<double>
data001[11] =
{
  { 8.5155875815486127e-05, 100.00000000000000, 1000.0000000000000, 0.0 },
  { 0.00012783770668059490, 100.00000000000000, 1100.0000000000000, 0.0 },
  { 0.00017868879940604636, 100.00000000000000, 1200.0000000000000, 0.0 },
  { 0.00023648188537313321, 100.00000000000000, 1300.0000000000000, 0.0 },
  { 0.00029987385710648917, 100.00000000000000, 1400.0000000000000, 0.0 },
  { 0.00036754116875294015, 100.00000000000000, 1500.0000000000000, 0.0 },
  { 0.00043825961543828500, 100.00000000000000, 1600.0000000000000, 0.0 },
  { 0.00051094422898288565, 100.00000000000000, 1700.0000000000000, 0.0 },
  { 0.00058466302033197875, 100.00000000000000, 1800.0000000000000, 0.0 },
  { 0.00065863487030979218, 100.00000000000000, 1900.0000000000000, 0.0 },
  { 0.00073221866792272470, 100.00000000000000, 2000.0000000000000, 0.0 },
};
const double toler001 = 1.0000000000000006e-10;

template<typename Ret, unsigned int Num>
  void
  test(const testcase_cyl_bessel_i_scaled<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::cyl_bessel_i_scaled(data[i].nu, data[i].x);
	const Ret f0 = data[i].f0;
	const Ret diff = f - f0;
	if (std::abs(diff) > max_abs_diff)
	  max_abs_diff = std::abs(diff);
	if (std::abs(f0) > Ret(10) * eps
	 && std::abs(f) > Ret(10) * eps)
	  {
	    const Ret frac = diff / f0;
	    if (std::abs(frac) > max_abs_frac)
	      max_abs_frac = std::abs(frac);
	  }
      }
    VERIFY(max_abs_frac < toler);
  }

int
main()
{
  test(data001, toler001);
  return 0;
}