// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
// { dg-skip-if "no extensions in strict dialects" { *-*-* } { "-std=c++*" } }

// Copyright (C) 2016-2018 Free Software Foundation, Inc.
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

//  cosint

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

// Test data.
// max(|f - f_GSL|): 8.8817841970012523e-16 at index 23
// max(|f - f_GSL| / |f_GSL|): 6.4503877186605986e-14
// mean(f - f_GSL): 7.8843181983145887e-17
// variance(f - f_GSL): 6.3424623459110902e-35
// stddev(f - f_GSL): 7.9639577760753415e-18
const testcase_cosint<double>
data001[100] =
{
  { -1.7278683866572964, 0.10000000000000001, 0.0 },
  { -1.0422055956727820, 0.20000000000000001, 0.0 },
  { -0.64917293297116130, 0.30000000000000004, 0.0 },
  { -0.37880934642524444, 0.40000000000000002, 0.0 },
  { -0.17778407880661296, 0.50000000000000000, 0.0 },
  { -0.022270706959279574, 0.60000000000000009, 0.0 },
  { 0.10051470700889809, 0.70000000000000007, 0.0 },
  { 0.19827861595246721, 0.80000000000000004, 0.0 },
  { 0.27606783046777283, 0.90000000000000002, 0.0 },
  { 0.33740392290096821, 1.0000000000000000, 0.0 },
  { 0.38487337742465083, 1.1000000000000001, 0.0 },
  { 0.42045918289424067, 1.2000000000000002, 0.0 },
  { 0.44573856752853458, 1.3000000000000000, 0.0 },
  { 0.46200658509467712, 1.4000000000000001, 0.0 },
  { 0.47035631719539994, 1.5000000000000000, 0.0 },
  { 0.47173251693187790, 1.6000000000000001, 0.0 },
  { 0.46696836417695470, 1.7000000000000002, 0.0 },
  { 0.45681112941833696, 1.8000000000000000, 0.0 },
  { 0.44194034968159884, 1.9000000000000001, 0.0 },
  { 0.42298082877486504, 2.0000000000000000, 0.0 },
  { 0.40051198784439640, 2.1000000000000001, 0.0 },
  { 0.37507459904983220, 2.2000000000000002, 0.0 },
  { 0.34717561754031617, 2.3000000000000003, 0.0 },
  { 0.31729161743669793, 2.4000000000000004, 0.0 },
  { 0.28587119636538361, 2.5000000000000000, 0.0 },
  { 0.25333661606258423, 2.6000000000000001, 0.0 },
  { 0.22008487863296156, 2.7000000000000002, 0.0 },
  { 0.18648838964317571, 2.8000000000000003, 0.0 },
  { 0.15289532415958818, 2.9000000000000004, 0.0 },
  { 0.11962978600800017, 3.0000000000000000, 0.0 },
  { 0.086991831195536884, 3.1000000000000001, 0.0 },
  { 0.055257411719942251, 3.2000000000000002, 0.0 },
  { 0.024678284607957846, 3.3000000000000003, 0.0 },
  { -0.0045180779307421037, 3.4000000000000004, 0.0 },
  { -0.032128548512480926, 3.5000000000000000, 0.0 },
  { -0.057974351859800821, 3.6000000000000001, 0.0 },
  { -0.081901001284298447, 3.7000000000000002, 0.0 },
  { -0.10377815035689775, 3.8000000000000003, 0.0 },
  { -0.12349934920781513, 3.9000000000000004, 0.0 },
  { -0.14098169788693049, 4.0000000000000000, 0.0 },
  { -0.15616539182812117, 4.1000000000000005, 0.0 },
  { -0.16901315676715675, 4.2000000000000002, 0.0 },
  { -0.17950957251263325, 4.2999999999999998, 0.0 },
  { -0.18766028680044072, 4.4000000000000004, 0.0 },
  { -0.19349112210173874, 4.5000000000000000, 0.0 },
  { -0.19704707972235622, 4.6000000000000005, 0.0 },
  { -0.19839124684247283, 4.7000000000000002, 0.0 },
  { -0.19760361330993520, 4.8000000000000007, 0.0 },
  { -0.19477980602623723, 4.9000000000000004, 0.0 },
  { -0.19002974965664385, 5.0000000000000000, 0.0 },
  { -0.18347626315929888, 5.1000000000000005, 0.0 },
  { -0.17525360226565947, 5.2000000000000002, 0.0 },
  { -0.16550595855892719, 5.3000000000000007, 0.0 },
  { -0.15438592619072442, 5.4000000000000004, 0.0 },
  { -0.14205294755151926, 5.5000000000000000, 0.0 },
  { -0.12867174936980774, 5.6000000000000005, 0.0 },
  { -0.11441078076167903, 5.7000000000000002, 0.0 },
  { -0.099440664689378469, 5.8000000000000007, 0.0 },
  { -0.083932674118556441, 5.9000000000000004, 0.0 },
  { -0.068057243893247132, 6.0000000000000000, 0.0 },
  { -0.051982528980021883, 6.1000000000000005, 0.0 },
  { -0.035873019273454966, 6.2000000000000002, 0.0 },
  { -0.019888220609842099, 6.3000000000000007, 0.0 },
  { -0.0041814110113350053, 6.4000000000000004, 0.0 },
  { 0.011101519514930106, 6.5000000000000000, 0.0 },
  { 0.025823138061263344, 6.6000000000000005, 0.0 },
  { 0.039855440047043486, 6.7000000000000002, 0.0 },
  { 0.053080716720199275, 6.8000000000000007, 0.0 },
  { 0.065392313975951549, 6.9000000000000004, 0.0 },
  { 0.076695278482184534, 7.0000000000000000, 0.0 },
  { 0.086906888071347888, 7.1000000000000005, 0.0 },
  { 0.095957064345180781, 7.2000000000000002, 0.0 },
  { 0.10378866643202767, 7.3000000000000007, 0.0 },
  { 0.11035766582837825, 7.4000000000000004, 0.0 },
  { 0.11563320323793427, 7.5000000000000000, 0.0 },
  { 0.11959752928456585, 7.6000000000000005, 0.0 },
  { 0.12224583191184631, 7.7000000000000002, 0.0 },
  { 0.12358595418360856, 7.8000000000000007, 0.0 },
  { 0.12363800705971785, 7.9000000000000004, 0.0 },
  { 0.12243388253200956, 8.0000000000000000, 0.0 },
  { 0.12001667326059656, 8.0999999999999996, 0.0 },
  { 0.11644000554456650, 8.2000000000000011, 0.0 },
  { 0.11176729308811123, 8.3000000000000007, 0.0 },
  { 0.10607091957864391, 8.4000000000000004, 0.0 },
  { 0.099431358573421905, 8.5000000000000000, 0.0 },
  { 0.091936239592257427, 8.5999999999999996, 0.0 },
  { 0.083679369633444065, 8.7000000000000011, 0.0 },
  { 0.074759719566173138, 8.8000000000000007, 0.0 },
  { 0.065280385004315378, 8.9000000000000004, 0.0 },
  { 0.055347531333133609, 9.0000000000000000, 0.0 },
  { 0.045069332542612181, 9.0999999999999996, 0.0 },
  { 0.034554913419759768, 9.2000000000000011, 0.0 },
  { 0.023913304469275642, 9.3000000000000007, 0.0 },
  { 0.013252418669884865, 9.4000000000000004, 0.0 },
  { 0.0026780588356506551, 9.5000000000000000, 0.0 },
  { -0.0077070360585336679, 9.6000000000000014, 0.0 },
  { -0.017804097705837581, 9.7000000000000011, 0.0 },
  { -0.027519181109809453, 9.8000000000000007, 0.0 },
  { -0.036763956296836403, 9.9000000000000004, 0.0 },
  { -0.045456433004455371, 10.000000000000000, 0.0 },
};
const double toler001 = 5.0000000000000029e-12;

template<typename Ret, unsigned int Num>
  void
  test(const testcase_cosint<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = -Ret(1);
    Ret max_abs_frac = -Ret(1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::cosint(data[i].x);
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