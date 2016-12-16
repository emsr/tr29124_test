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

//  conf_hyperg
//  Test relations between special functions.

#include <limits>
#include <ext/cmath>

#if defined(__TEST_DEBUG)
#  include <iostream>
#  define VERIFY(A) \
  if (!(A)) \
    { \
      std::cout << "line " << __LINE__ \
	<< "  max_abs_frac = " << stats.max_abs_frac \
	<< std::endl; \
    } \
  else \
    { \
      std::cout << "Success: " \
	<< "  samples = " << stats.num \
	<< "  mean = " << stats.mean_diff() \
	<< "  variance = " << stats.variance_diff() \
	<< std::endl; \
    }
#else
#  include <testsuite_hooks.h>
#endif

#include <specfun_stats.h>

template<typename Tp>
  void
  test_cyl_bessel_j(Tp nu, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = std::pow(z / Tp{2}, nu)
		 * __gnu_cxx::conf_hyperg_lim(nu + Tp{1}, -z * z / Tp{4});
	auto Jnu = std::cyl_bessel_j(nu, z);
	stats << std::make_pair(hyp, Jnu);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    test_cyl_bessel_j(Tp{0.5L}, toler);
    test_cyl_bessel_j(Tp{1.5L}, toler);
    test_cyl_bessel_j(Tp{2.0L}, toler);
    test_cyl_bessel_j(Tp{3.5L}, toler);
    test_cyl_bessel_j(Tp{5.0L}, toler);
  }

int
main()
{
  test<double>();
}
