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

//  comp_ellint_1/comp_ellint_2
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
  test_legendre(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_pi2 = __gnu_cxx::__math_constants<Tp>::__pi_half;
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto kp = std::sqrt(Tp{1} - k * k);
	auto K = std::comp_ellint_1(k);
	auto Kp = std::comp_ellint_1(kp);
	auto E = std::comp_ellint_2(k);
	auto Ep = std::comp_ellint_2(kp);
	stats << std::make_pair(K * Ep + E * Kp - K * Kp, _S_pi2);
      }
    VERIFY(stats.max_abs_frac < toler);
  }


template<typename Tp>
  void
  test_theta_3(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
   const auto _S_pi2 = __gnu_cxx::__math_constants<Tp>::__pi_half;
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto q = __gnu_cxx::ellnome(k);
	auto K = std::comp_ellint_1(k);
	auto tht3 = __gnu_cxx::theta_3(Tp{0}, q);
	stats << std::make_pair(_S_pi2 * tht3 * tht3, K);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    test_legendre(toler);

    test_theta_3(toler);
  }

int
main()
{
  test<double>();
}
