// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
// { dg-skip-if "no extensions in strict dialects" { *-*-* } { "-std=c++*" } }

// Copyright (C) 2019 Free Software Foundation, Inc.
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

//  bell

//  Compare against values from The On-Line Encyclopedia of Integer Sequences:
//  https://oeis.org/A000110

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

template<typename Ret>
  void
  test()
  {
    constexpr unsigned int num_bells = 20;
    constexpr uint64_t
    bell[num_bells]
    {
      1ull,
      1ull,
      2ull,
      5ull,
      15ull,
      52ull,
      203ull,
      877ull,
      4140ull,
      21147ull,
      115975ull,
      678570ull,
      4213597ull,
      27644437ull,
      190899322ull,
      1382958545ull,
      10480142147ull,
      82864869804ull,
      682076806159ull,
      5832742205057ull
    };

    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    const Ret toler = 100 * eps;
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    const auto bells = emsr::bell<Ret>(num_bells);
    for (unsigned int i = 0; i < num_bells; ++i)
      {
	const Ret f0 = bell[i];
	const Ret f = bells[i];
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
  test<long double>();
  return 0;
}
