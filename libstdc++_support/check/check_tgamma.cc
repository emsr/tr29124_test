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

//  tgamma
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

// Test data for a=0.50000000000000000.
// max(|f - f_Boost|): 6.6613381477509392e-16 at index 1
// max(|f - f_Boost| / |f_Boost|): 2.2370451178473489e-15
// mean(f - f_Boost): 1.0672491930611856e-16
// variance(f - f_Boost): 1.3043531529359038e-33
// stddev(f - f_Boost): 3.6115829672539765e-17
const testcase_tgamma<double>
data001[11] =
{
  { 1.7724538509055161, 0.50000000000000000, 0.0000000000000000, 0.0 },
  { 0.56241823159440707, 0.50000000000000000, 0.50000000000000000, 0.0 },
  { 0.27880558528066196, 0.50000000000000000, 1.0000000000000000, 0.0 },
  { 0.14758251320409640, 0.50000000000000000, 1.5000000000000000, 0.0 },
  { 0.080647117960317691, 0.50000000000000000, 2.0000000000000000, 0.0 },
  { 0.044926952600007938, 0.50000000000000000, 2.5000000000000000, 0.0 },
  { 0.025356509323463443, 0.50000000000000000, 3.0000000000000000, 0.0 },
  { 0.014447220989525332, 0.50000000000000000, 3.5000000000000000, 0.0 },
  { 0.0082910693806726669, 0.50000000000000000, 4.0000000000000000, 0.0 },
  { 0.0047852639289850743, 0.50000000000000000, 4.5000000000000000, 0.0 },
  { 0.0027746032604128094, 0.50000000000000000, 5.0000000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

// Test data for a=1.0000000000000000.
// max(|f - f_Boost|): 1.0796918914479647e-14 at index 3
// max(|f - f_Boost| / |f_Boost|): 4.8388433492349763e-14
// mean(f - f_Boost): 3.8824688413702738e-15
// variance(f - f_Boost): 1.6566107023167712e-30
// stddev(f - f_Boost): 1.2870938980186221e-15
const testcase_tgamma<double>
data002[11] =
{
  { 1.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.60653065971263342, 1.0000000000000000, 0.50000000000000000, 0.0 },
  { 0.36787944117144233, 1.0000000000000000, 1.0000000000000000, 0.0 },
  { 0.22313016014842982, 1.0000000000000000, 1.5000000000000000, 0.0 },
  { 0.13533528323661270, 1.0000000000000000, 2.0000000000000000, 0.0 },
  { 0.082084998623898800, 1.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.049787068367863944, 1.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.030197383422318501, 1.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.018315638888734179, 1.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.011108996538242306, 1.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.0067379469990854670, 1.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler002 = 2.5000000000000015e-12;

// Test data for a=1.5000000000000000.
// max(|f - f_Boost|): 4.4408920985006262e-16 at index 1
// max(|f - f_Boost| / |f_Boost|): 2.3851749495636526e-15
// mean(f - f_Boost): 2.1069005126409220e-16
// variance(f - f_Boost): 4.4124039876451316e-33
// stddev(f - f_Boost): 6.6425928579472130e-17
const testcase_tgamma<double>
data003[11] =
{
  { 0.88622692545275805, 1.5000000000000000, 0.0000000000000000, 0.0 },
  { 0.71009105827755692, 1.5000000000000000, 0.50000000000000000, 0.0 },
  { 0.50728223381177329, 1.5000000000000000, 1.0000000000000000, 0.0 },
  { 0.34706877589662155, 1.5000000000000000, 1.5000000000000000, 0.0 },
  { 0.23171655200098068, 1.5000000000000000, 2.0000000000000000, 0.0 },
  { 0.15225125499165762, 1.5000000000000000, 2.5000000000000000, 0.0 },
  { 0.098911986634777363, 1.5000000000000000, 3.0000000000000000, 0.0 },
  { 0.063717741866444161, 1.5000000000000000, 3.5000000000000000, 0.0 },
  { 0.040776812467804693, 1.5000000000000000, 4.0000000000000000, 0.0 },
  { 0.025958372317599586, 1.5000000000000000, 4.5000000000000000, 0.0 },
  { 0.016453809148952222, 1.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler003 = 2.5000000000000020e-13;

// Test data for a=2.0000000000000000.
// max(|f - f_Boost|): 1.0880185641326534e-14 at index 3
// max(|f - f_Boost| / |f_Boost|): 3.6904707012474476e-14
// mean(f - f_Boost): 5.8450719012367617e-15
// variance(f - f_Boost): 3.7760020604631092e-30
// stddev(f - f_Boost): 1.9431937784130304e-15
const testcase_tgamma<double>
data004[11] =
{
  { 1.0000000000000000, 2.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.90979598956895014, 2.0000000000000000, 0.50000000000000000, 0.0 },
  { 0.73575888234288467, 2.0000000000000000, 1.0000000000000000, 0.0 },
  { 0.55782540037107453, 2.0000000000000000, 1.5000000000000000, 0.0 },
  { 0.40600584970983805, 2.0000000000000000, 2.0000000000000000, 0.0 },
  { 0.28729749518364578, 2.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.19914827347145578, 2.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.13588822540043324, 2.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.091578194443670907, 2.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.061099480960332686, 2.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.040427681994512805, 2.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler004 = 2.5000000000000015e-12;

// Test data for a=2.5000000000000000.
// max(|f - f_Boost|): 1.0547118733938987e-15 at index 6
// max(|f - f_Boost| / |f_Boost|): 2.5909892861240649e-15
// mean(f - f_Boost): 3.1161941713910642e-16
// variance(f - f_Boost): 2.0960420838387484e-32
// stddev(f - f_Boost): 1.4477714197478650e-16
const testcase_tgamma<double>
data005[11] =
{
  { 1.3293403881791370, 2.5000000000000000, 0.0000000000000000, 0.0 },
  { 1.2795775586565121, 2.5000000000000000, 0.50000000000000000, 0.0 },
  { 1.1288027918891024, 2.5000000000000000, 1.0000000000000000, 0.0 },
  { 0.93051944278679244, 2.5000000000000000, 1.5000000000000000, 0.0 },
  { 0.73036081404311470, 2.5000000000000000, 2.0000000000000000, 0.0 },
  { 0.55284632921662058, 2.5000000000000000, 2.5000000000000000, 0.0 },
  { 0.40706917587130298, 2.5000000000000000, 3.0000000000000000, 0.0 },
  { 0.29330607260055147, 2.5000000000000000, 3.5000000000000000, 0.0 },
  { 0.20769032981158048, 2.5000000000000000, 4.0000000000000000, 0.0 },
  { 0.14498339006538111, 2.5000000000000000, 4.5000000000000000, 0.0 },
  { 0.10001325131715742, 2.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler005 = 2.5000000000000020e-13;

// Test data for a=3.0000000000000000.
// max(|f - f_Boost|): 2.2093438190040615e-14 at index 6
// max(|f - f_Boost| / |f_Boost|): 3.3218744231521981e-14
// mean(f - f_Boost): 1.5575924388661570e-14
// variance(f - f_Boost): 2.7164705406550178e-29
// stddev(f - f_Boost): 5.2119771110923133e-15
const testcase_tgamma<double>
data006[11] =
{
  { 2.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0 },
  { 1.9712246440660586, 3.0000000000000000, 0.50000000000000000, 0.0 },
  { 1.8393972058572117, 3.0000000000000000, 1.0000000000000000, 0.0 },
  { 1.6176936610761163, 3.0000000000000000, 1.5000000000000000, 0.0 },
  { 1.3533528323661270, 3.0000000000000000, 2.0000000000000000, 0.0 },
  { 1.0876262317666590, 3.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.84638016225368706, 3.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.64169439772426817, 3.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.47620661110708867, 3.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.34715614182007209, 3.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.24930403896616229, 3.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler006 = 2.5000000000000015e-12;

// Test data for a=3.5000000000000000.
// max(|f - f_Boost|): 2.2204460492503131e-15 at index 8
// max(|f - f_Boost| / |f_Boost|): 2.0088602900232915e-15
// mean(f - f_Boost): 8.9827135628762675e-16
// variance(f - f_Boost): 3.5140167596208710e-32
// stddev(f - f_Boost): 1.8745710868411660e-16
const testcase_tgamma<double>
data007[11] =
{
  { 3.3233509704478426, 3.5000000000000000, 0.0000000000000000, 0.0 },
  { 3.3061643822613687, 3.5000000000000000, 0.50000000000000000, 0.0 },
  { 3.1898864208941982, 3.5000000000000000, 1.0000000000000000, 0.0 },
  { 2.9411730253797712, 3.5000000000000000, 1.5000000000000000, 0.0 },
  { 2.5914740071910742, 3.5000000000000000, 2.0000000000000000, 0.0 },
  { 2.1932894398643867, 3.5000000000000000, 2.5000000000000000, 0.0 },
  { 1.7937765274356683, 3.5000000000000000, 3.0000000000000000, 0.0 },
  { 1.4253182908044768, 3.5000000000000000, 3.5000000000000000, 0.0 },
  { 1.1053262689684449, 3.5000000000000000, 4.0000000000000000, 0.0 },
  { 0.83966471731387049, 3.5000000000000000, 4.5000000000000000, 0.0 },
  { 0.62669581626153903, 3.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler007 = 2.5000000000000020e-13;

// Test data for a=4.0000000000000000.
// max(|f - f_Boost|): 6.5281113847959205e-14 at index 6
// max(|f - f_Boost| / |f_Boost|): 3.0056056811349900e-14
// mean(f - f_Boost): 5.7933456012258174e-14
// variance(f - f_Boost): 3.7487315647752206e-28
// stddev(f - f_Boost): 1.9361641368373759e-14
const testcase_tgamma<double>
data008[11] =
{
  { 6.0000000000000000, 4.0000000000000000, 0.0000000000000000, 0.0 },
  { 5.9894902646622548, 4.0000000000000000, 0.50000000000000000, 0.0 },
  { 5.8860710587430773, 4.0000000000000000, 1.0000000000000000, 0.0 },
  { 5.6061452737292994, 4.0000000000000000, 1.5000000000000000, 0.0 },
  { 5.1427407629912825, 4.0000000000000000, 2.0000000000000000, 0.0 },
  { 4.5454567987983960, 4.0000000000000000, 2.5000000000000000, 0.0 },
  { 3.8833913326933875, 4.0000000000000000, 3.0000000000000000, 0.0 },
  { 3.2197960074047103, 4.0000000000000000, 3.5000000000000000, 0.0 },
  { 2.6008207222002535, 4.0000000000000000, 4.0000000000000000, 0.0 },
  { 2.0537757350075463, 4.0000000000000000, 4.5000000000000000, 0.0 },
  { 1.5901554917841703, 4.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler008 = 2.5000000000000015e-12;

// Test data for a=4.5000000000000000.
// max(|f - f_Boost|): 1.2434497875801753e-14 at index 10
// max(|f - f_Boost| / |f_Boost|): 3.0501015627178965e-15
// mean(f - f_Boost): 6.5402229087009219e-15
// variance(f - f_Boost): 3.8216725126570655e-30
// stddev(f - f_Boost): 1.9549098477057875e-15
const testcase_tgamma<double>
data009[11] =
{
  { 11.631728396567448, 4.5000000000000000, 0.0000000000000000, 0.0 },
  { 11.625185580724835, 4.5000000000000000, 0.50000000000000000, 0.0 },
  { 11.532481914301135, 4.5000000000000000, 1.0000000000000000, 0.0 },
  { 11.216417216448384, 4.5000000000000000, 1.5000000000000000, 0.0 },
  { 10.601302969335334, 4.5000000000000000, 2.0000000000000000, 0.0 },
  { 9.7044470815824422, 4.5000000000000000, 2.5000000000000000, 0.0 },
  { 8.6065286092970723, 4.5000000000000000, 3.0000000000000000, 0.0 },
  { 7.4107999003765128, 4.5000000000000000, 3.5000000000000000, 0.0 },
  { 6.2130437191475325, 4.5000000000000000, 4.0000000000000000, 0.0 },
  { 5.0862546002754270, 4.5000000000000000, 4.5000000000000000, 0.0 },
  { 4.0767487967586140, 4.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler009 = 2.5000000000000020e-13;

// Test data for a=5.0000000000000000.
// max(|f - f_Boost|): 2.6467716907063732e-13 at index 10
// max(|f - f_Boost| / |f_Boost|): 2.5036057873868745e-14
// mean(f - f_Boost): 2.5854066362543643e-13
// variance(f - f_Boost): 4.1422368986877857e-30
// stddev(f - f_Boost): 2.0352486085704089e-15
const testcase_tgamma<double>
data010[11] =
{
  { 24.000000000000000, 5.0000000000000000, 0.0000000000000000, 0.0 },
  { 23.995869224881059, 5.0000000000000000, 0.50000000000000000, 0.0 },
  { 23.912163676143752, 5.0000000000000000, 1.0000000000000000, 0.0 },
  { 23.554177530668625, 5.0000000000000000, 1.5000000000000000, 0.0 },
  { 22.736327583750931, 5.0000000000000000, 2.0000000000000000, 0.0 },
  { 21.388272453939631, 5.0000000000000000, 2.5000000000000000, 0.0 },
  { 19.566317868570529, 5.0000000000000000, 3.0000000000000000, 0.0 },
  { 17.410678879430510, 5.0000000000000000, 3.5000000000000000, 0.0 },
  { 15.092086444316964, 5.0000000000000000, 4.0000000000000000, 0.0 },
  { 12.770485832993172, 5.0000000000000000, 4.5000000000000000, 0.0 },
  { 10.571838841565098, 5.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler010 = 2.5000000000000015e-12;

template<typename Ret, unsigned int Num>
  void
  test(const testcase_tgamma<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = -Ret(1);
    Ret max_abs_frac = -Ret(1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::tgamma(data[i].a, data[i].x);
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
  test(data002, toler002);
  test(data003, toler003);
  test(data004, toler004);
  test(data005, toler005);
  test(data006, toler006);
  test(data007, toler007);
  test(data008, toler008);
  test(data009, toler009);
  test(data010, toler010);
  return 0;
}