// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
//
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

//  sph_bessel_i
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
	<< std::endl; \
    }
#else
#  include <testsuite_hooks.h>
#endif
#include <specfun_testcase.h>


// Failure at n=0 x=0.0000000000000000 f=-nan f_GSL=1.0000000000000000
// Failure at n=1 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=2 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=5 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=10 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=20 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=50 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=100 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000//  sph_bessel_i

// Failure at n=0 x=0.0000000000000000 f=-nan f_GSL=1.0000000000000000
// Failure at n=1 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=2 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=5 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=10 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=20 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=50 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
// Failure at n=100 x=0.0000000000000000 f=-nan f_GSL=0.0000000000000000
template<typename Ret, unsigned int Num>
  void
  test(const testcase_sph_bessel_i<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = -Ret(1);
    Ret max_abs_frac = -Ret(1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::sph_bessel_i(data[i].n, data[i].x);
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
  return 0;
}