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

//  rising_factorial
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

// Test data for a=0.25000000000000000.
// max(|f - f_Boost|): 9.5268205270873786e+139 at index 7
// max(|f - f_Boost| / |f_Boost|): 2.0512976309882383e-16
// mean(f - f_Boost): 1.1908525658859223e+139
// variance(f - f_Boost): 1.1345038669416680e+279
// stddev(f - f_Boost): 3.3682396989253424e+139
const testcase_rising_factorial<double>
data001[8] =
{
  { 1.0000000000000000, 0.25000000000000000, 0.0000000000000000, 0.0 },
  { 0.25000000000000000, 0.25000000000000000, 1.0000000000000000, 0.0 },
  { 0.31250000000000000, 0.25000000000000000, 2.0000000000000000, 0.0 },
  { 9.7119140625000000, 0.25000000000000000, 5.0000000000000000, 0.0 },
  { 176310.36293506622, 0.25000000000000000, 10.000000000000000, 0.0 },
  { 70619879085128896., 0.25000000000000000, 20.000000000000000, 0.0 },
  { 4.4529789803918110e+62, 0.25000000000000000, 50.000000000000000, 0.0 },
  { 8.1323318185788312e+155, 0.25000000000000000, 100.00000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

// Test data for a=0.50000000000000000.
// max(|f - f_Boost|): 7.6214564216699029e+140 at index 7
// max(|f - f_Boost| / |f_Boost|): 2.0982602304152644e-16
// mean(f - f_Boost): 9.5268205270873786e+139
// variance(f - f_Boost): 7.2608247484266751e+280
// stddev(f - f_Boost): 2.6945917591402739e+140
const testcase_rising_factorial<double>
data002[8] =
{
  { 1.0000000000000000, 0.50000000000000000, 0.0000000000000000, 0.0 },
  { 0.50000000000000000, 0.50000000000000000, 1.0000000000000000, 0.0 },
  { 0.75000000000000000, 0.50000000000000000, 2.0000000000000000, 0.0 },
  { 29.531249999999996, 0.50000000000000000, 5.0000000000000000, 0.0 },
  { 639383.86230468750, 0.50000000000000000, 10.000000000000000, 0.0 },
  { 3.0501459767616064e+17, 0.50000000000000000, 20.000000000000000, 0.0 },
  { 2.4206344837469458e+63, 0.50000000000000000, 50.000000000000000, 0.0 },
  { 5.2587902919564293e+156, 0.50000000000000000, 100.00000000000000, 0.0 },
};
const double toler002 = 2.5000000000000020e-13;

// Test data for a=0.75000000000000000.
// max(|f - f_Boost|): 3.0485825686679612e+141 at index 7
// max(|f - f_Boost| / |f_Boost|): 3.1375717126919433e-16
// mean(f - f_Boost): 3.8107282108349515e+140
// variance(f - f_Boost): 1.1617319597482680e+282
// stddev(f - f_Boost): 1.0778367036561096e+141
const testcase_rising_factorial<double>
data003[8] =
{
  { 1.0000000000000000, 0.75000000000000000, 0.0000000000000000, 0.0 },
  { 0.75000000000000000, 0.75000000000000000, 1.0000000000000000, 0.0 },
  { 1.3125000000000000, 0.75000000000000000, 2.0000000000000000, 0.0 },
  { 64.291992187500000, 0.75000000000000000, 5.0000000000000000, 0.0 },
  { 1649843.9631700516, 0.75000000000000000, 10.000000000000000, 0.0 },
  { 9.3445029406828723e+17, 0.75000000000000000, 20.000000000000000, 0.0 },
  { 9.3161321630923168e+63, 0.75000000000000000, 50.000000000000000, 0.0 },
  { 2.4060965407256963e+157, 0.75000000000000000, 100.00000000000000, 0.0 },
};
const double toler003 = 2.5000000000000020e-13;

// Test data for a=1.0000000000000000.
// max(|f - f_Boost|): 5.8460065493236117e+48 at index 6
// max(|f - f_Boost| / |f_Boost|): 1.9221373823482188e-16
// mean(f - f_Boost): -7.3075081866545146e+47
// variance(f - f_Boost): 2.0431220958214846e+193
// stddev(f - f_Boost): 4.5200908130495396e+96
const testcase_rising_factorial<double>
data004[8] =
{
  { 1.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0 },
  { 1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 0.0 },
  { 2.0000000000000000, 1.0000000000000000, 2.0000000000000000, 0.0 },
  { 120.00000000000000, 1.0000000000000000, 5.0000000000000000, 0.0 },
  { 3628800.0000000005, 1.0000000000000000, 10.000000000000000, 0.0 },
  { 2.4329020081766400e+18, 1.0000000000000000, 20.000000000000000, 0.0 },
  { 3.0414093201713381e+64, 1.0000000000000000, 50.000000000000000, 0.0 },
  { 9.3326215443944151e+157, 1.0000000000000000, 100.00000000000000, 0.0 },
};
const double toler004 = 2.5000000000000020e-13;

// Test data for a=1.2500000000000000.
// max(|f - f_Boost|): 4.8777321098687379e+142 at index 7
// max(|f - f_Boost| / |f_Boost|): 2.6125979777760648e-16
// mean(f - f_Boost): 6.0971651373359223e+141
// variance(f - f_Boost): 2.9740338169555661e+284
// stddev(f - f_Boost): 1.7245387258497753e+142
const testcase_rising_factorial<double>
data005[8] =
{
  { 1.0000000000000000, 1.2500000000000000, 0.0000000000000000, 0.0 },
  { 1.2500000000000000, 1.2500000000000000, 1.0000000000000000, 0.0 },
  { 2.8125000000000000, 1.2500000000000000, 2.0000000000000000, 0.0 },
  { 203.95019531250000, 1.2500000000000000, 5.0000000000000000, 0.0 },
  { 7228724.8803377151, 1.2500000000000000, 10.000000000000000, 0.0 },
  { 5.7202102058954414e+18, 1.2500000000000000, 20.000000000000000, 0.0 },
  { 8.9504877505875406e+64, 1.2500000000000000, 50.000000000000000, 0.0 },
  { 3.2610650592501111e+158, 1.2500000000000000, 100.00000000000000, 0.0 },
};
const double toler005 = 2.5000000000000020e-13;

// Test data for a=1.5000000000000000.
// max(|f - f_Boost|): 4.6768052394588893e+49 at index 6
// max(|f - f_Boost| / |f_Boost|): 1.9129283120516278e-16
// mean(f - f_Boost): -5.8460065493236117e+48
// variance(f - f_Boost): 8.3686281044848008e+196
// stddev(f - f_Boost): 2.8928581203517054e+98
const testcase_rising_factorial<double>
data006[8] =
{
  { 1.0000000000000000, 1.5000000000000000, 0.0000000000000000, 0.0 },
  { 1.5000000000000000, 1.5000000000000000, 1.0000000000000000, 0.0 },
  { 3.7500000000000000, 1.5000000000000000, 2.0000000000000000, 0.0 },
  { 324.84375000000000, 1.5000000000000000, 5.0000000000000000, 0.0 },
  { 13427061.108398438, 1.5000000000000000, 10.000000000000000, 0.0 },
  { 1.2505598504722584e+19, 1.5000000000000000, 20.000000000000000, 0.0 },
  { 2.4448408285844155e+65, 1.5000000000000000, 50.000000000000000, 0.0 },
  { 1.0570168486832425e+159, 1.5000000000000000, 100.00000000000000, 0.0 },
};
const double toler006 = 2.5000000000000020e-13;

// Test data for a=1.7500000000000000.
// max(|f - f_Boost|): 1.9510928439474951e+144 at index 7
// max(|f - f_Boost| / |f_Boost|): 6.0364428862954360e-16
// mean(f - f_Boost): -2.4388660549343689e+143
// variance(f - f_Boost): 4.7584541071289058e+287
// stddev(f - f_Boost): 6.8981549033991011e+143
const testcase_rising_factorial<double>
data007[8] =
{
  { 1.0000000000000000, 1.7500000000000000, 0.0000000000000000, 0.0 },
  { 1.7500000000000000, 1.7500000000000000, 1.0000000000000000, 0.0 },
  { 4.8125000000000000, 1.7500000000000000, 2.0000000000000000, 0.0 },
  { 492.90527343750006, 1.7500000000000000, 5.0000000000000000, 0.0 },
  { 23647763.472104073, 1.7500000000000000, 10.000000000000000, 0.0 },
  { 2.5853124802555945e+19, 1.7500000000000000, 20.000000000000000, 0.0 },
  { 6.3039160970257999e+65, 1.7500000000000000, 50.000000000000000, 0.0 },
  { 3.2321896863748521e+159, 1.7500000000000000, 100.00000000000000, 0.0 },
};
const double toler007 = 2.5000000000000020e-13;

// Test data for a=2.0000000000000000.
// max(|f - f_Boost|): 1.8707220957835557e+50 at index 6
// max(|f - f_Boost| / |f_Boost|): 1.2060469850028041e-16
// mean(f - f_Boost): -2.3384026197294447e+49
// variance(f - f_Boost): 2.1423687947481090e+199
// stddev(f - f_Boost): 4.6285729925627286e+99
const testcase_rising_factorial<double>
data008[8] =
{
  { 1.0000000000000000, 2.0000000000000000, 0.0000000000000000, 0.0 },
  { 2.0000000000000000, 2.0000000000000000, 1.0000000000000000, 0.0 },
  { 6.0000000000000000, 2.0000000000000000, 2.0000000000000000, 0.0 },
  { 720.00000000000000, 2.0000000000000000, 5.0000000000000000, 0.0 },
  { 39916800.000000000, 2.0000000000000000, 10.000000000000000, 0.0 },
  { 5.1090942171709440e+19, 2.0000000000000000, 20.000000000000000, 0.0 },
  { 1.5511187532873824e+66, 2.0000000000000000, 50.000000000000000, 0.0 },
  { 9.4259477598383599e+159, 2.0000000000000000, 100.00000000000000, 0.0 },
};
const double toler008 = 2.5000000000000020e-13;

// Test data for a=2.2500000000000000.
// max(|f - f_Boost|): 1.8730491301895953e+145 at index 7
// max(|f - f_Boost| / |f_Boost|): 7.0909543095688636e-16
// mean(f - f_Boost): 2.3413114127369942e+144
// variance(f - f_Boost): 4.3853913051299996e+289
// stddev(f - f_Boost): 6.6222287072631362e+144
const testcase_rising_factorial<double>
data009[8] =
{
  { 1.0000000000000000, 2.2500000000000000, 0.0000000000000000, 0.0 },
  { 2.2500000000000000, 2.2500000000000000, 1.0000000000000000, 0.0 },
  { 7.3124999999999991, 2.2500000000000000, 2.0000000000000000, 0.0 },
  { 1019.7509765625000, 2.2500000000000000, 5.0000000000000000, 0.0 },
  { 65058523.923039436, 2.2500000000000000, 10.000000000000000, 0.0 },
  { 9.7243573500222505e+19, 2.2500000000000000, 20.000000000000000, 0.0 },
  { 3.6696999777408920e+66, 2.2500000000000000, 50.000000000000000, 0.0 },
  { 2.6414626979925898e+160, 2.2500000000000000, 100.00000000000000, 0.0 },
};
const double toler009 = 2.5000000000000020e-13;

// Test data for a=2.5000000000000000.
// max(|f - f_Boost|): 1.2486994201263969e+145 at index 7
// max(|f - f_Boost| / |f_Boost|): 3.6561857219714646e-16
// mean(f - f_Boost): 1.5608742751579961e+144
// variance(f - f_Boost): 1.9490628022799998e+289
// stddev(f - f_Boost): 4.4148191381754248e+144
const testcase_rising_factorial<double>
data010[8] =
{
  { 1.0000000000000000, 2.5000000000000000, 0.0000000000000000, 0.0 },
  { 2.5000000000000000, 2.5000000000000000, 1.0000000000000000, 0.0 },
  { 8.7500000000000000, 2.5000000000000000, 2.0000000000000000, 0.0 },
  { 1407.6562500000000, 2.5000000000000000, 5.0000000000000000, 0.0 },
  { 102940801.83105469, 2.5000000000000000, 10.000000000000000, 0.0 },
  { 1.7924691190102373e+20, 2.5000000000000000, 20.000000000000000, 0.0 },
  { 8.3939535114731585e+66, 2.5000000000000000, 50.000000000000000, 0.0 },
  { 7.1524806760899398e+160, 2.5000000000000000, 100.00000000000000, 0.0 },
};
const double toler010 = 2.5000000000000020e-13;

// Test data for a=2.7500000000000000.
// max(|f - f_Boost|): 2.4973988402527938e+145 at index 7
// max(|f - f_Boost| / |f_Boost|): 1.6056338999483010e-16
// mean(f - f_Boost): -3.1217485503159922e+144
// variance(f - f_Boost): 7.7962512091199993e+289
// stddev(f - f_Boost): 8.8296382763508497e+144
const testcase_rising_factorial<double>
data011[8] =
{
  { 1.0000000000000000, 2.7500000000000000, 0.0000000000000000, 0.0 },
  { 2.7500000000000000, 2.7500000000000000, 1.0000000000000000, 0.0 },
  { 10.312500000000000, 2.7500000000000000, 2.0000000000000000, 0.0 },
  { 1901.2060546875000, 2.7500000000000000, 5.0000000000000000, 0.0 },
  { 158777840.45555592, 2.7500000000000000, 10.000000000000000, 0.0 },
  { 3.2131740826033816e+20, 2.7500000000000000, 20.000000000000000, 0.0 },
  { 1.8641580458347725e+67, 2.7500000000000000, 50.000000000000000, 0.0 },
  { 1.8792874319350927e+161, 2.7500000000000000, 100.00000000000000, 0.0 },
};
const double toler011 = 2.5000000000000020e-13;

// Test data for a=3.0000000000000000.
// max(|f - f_Boost|): 65536.000000000000 at index 5
// max(|f - f_Boost| / |f_Boost|): 1.1661202413912816e-16
// mean(f - f_Boost): 8192.0000000000000
// variance(f - f_Boost): 1.5628783020470605e+35
// stddev(f - f_Boost): 3.9533255646949453e+17
const testcase_rising_factorial<double>
data012[8] =
{
  { 1.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0 },
  { 3.0000000000000000, 3.0000000000000000, 1.0000000000000000, 0.0 },
  { 12.000000000000000, 3.0000000000000000, 2.0000000000000000, 0.0 },
  { 2520.0000000000000, 3.0000000000000000, 5.0000000000000000, 0.0 },
  { 239500800.00000000, 3.0000000000000000, 10.000000000000000, 0.0 },
  { 5.6200036388880377e+20, 3.0000000000000000, 20.000000000000000, 0.0 },
  { 4.0329087585471938e+67, 3.0000000000000000, 50.000000000000000, 0.0 },
  { 4.8072333575175635e+161, 3.0000000000000000, 100.00000000000000, 0.0 },
};
const double toler012 = 2.5000000000000020e-13;

// Test data for a=3.2500000000000000.
// max(|f - f_Boost|): 1.1972621413014757e+52 at index 6
// max(|f - f_Boost| / |f_Boost|): 1.4049307365611047e-16
// mean(f - f_Boost): -1.4965776766268446e+51
// variance(f - f_Boost): 3.5942984021148690e+206
// stddev(f - f_Boost): 1.8958634977536936e+103
const testcase_rising_factorial<double>
data013[8] =
{
  { 1.0000000000000000, 3.2500000000000000, 0.0000000000000000, 0.0 },
  { 3.2500000000000000, 3.2500000000000000, 1.0000000000000000, 0.0 },
  { 13.812499999999998, 3.2500000000000000, 2.0000000000000000, 0.0 },
  { 3285.8642578125000, 3.2500000000000000, 5.0000000000000000, 0.0 },
  { 354207519.13654804, 3.2500000000000000, 10.000000000000000, 0.0 },
  { 9.6163089350220028e+20, 3.2500000000000000, 20.000000000000000, 0.0 },
  { 8.5218588371982933e+67, 3.2500000000000000, 50.000000000000000, 0.0 },
  { 1.2003980483099661e+162, 3.2500000000000000, 100.00000000000000, 0.0 },
};
const double toler013 = 2.5000000000000020e-13;

// Test data for a=3.5000000000000000.
// max(|f - f_Boost|): 3.9958381444044701e+146 at index 7
// max(|f - f_Boost| / |f_Boost|): 2.7168357719013416e-16
// mean(f - f_Boost): 4.9947976805055876e+145
// variance(f - f_Boost): 1.9958403095347198e+292
// stddev(f - f_Boost): 1.4127421242161359e+146
const testcase_rising_factorial<double>
data014[8] =
{
  { 1.0000000000000000, 3.5000000000000000, 0.0000000000000000, 0.0 },
  { 3.5000000000000000, 3.5000000000000000, 1.0000000000000000, 0.0 },
  { 15.750000000000000, 3.5000000000000000, 2.0000000000000000, 0.0 },
  { 4222.9687500000000, 3.5000000000000000, 5.0000000000000000, 0.0 },
  { 514704009.15527344, 3.5000000000000000, 10.000000000000000, 0.0 },
  { 1.6132222071092132e+21, 3.5000000000000000, 20.000000000000000, 0.0 },
  { 1.7627302374093634e+68, 3.5000000000000000, 50.000000000000000, 0.0 },
  { 2.9325170771968750e+162, 3.5000000000000000, 100.00000000000000, 0.0 },
};
const double toler014 = 2.5000000000000020e-13;

// Test data for a=3.7500000000000000.
// max(|f - f_Boost|): 4.7950057732853641e+147 at index 7
// max(|f - f_Boost| / |f_Boost|): 6.8288376577574558e-16
// mean(f - f_Boost): -5.9937572166067051e+146
// variance(f - f_Boost): 2.8740100457299965e+294
// stddev(f - f_Boost): 1.6952905490593629e+147
const testcase_rising_factorial<double>
data015[8] =
{
  { 1.0000000000000000, 3.7500000000000000, 0.0000000000000000, 0.0 },
  { 3.7500000000000000, 3.7500000000000000, 1.0000000000000000, 0.0 },
  { 17.812500000000000, 3.7500000000000000, 2.0000000000000000, 0.0 },
  { 5357.9443359375000, 3.7500000000000000, 5.0000000000000000, 0.0 },
  { 736151805.74848652, 3.7500000000000000, 10.000000000000000, 0.0 },
  { 2.6581712865173428e+21, 3.7500000000000000, 20.000000000000000, 0.0 },
  { 3.5757940697376094e+68, 3.7500000000000000, 50.000000000000000, 0.0 },
  { 7.0217012229574833e+162, 3.7500000000000000, 100.00000000000000, 0.0 },
};
const double toler015 = 2.5000000000000020e-13;

// Test data for a=4.0000000000000000.
// max(|f - f_Boost|): 9.5780971304118054e+52 at index 6
// max(|f - f_Boost| / |f_Boost|): 1.3534147347811433e-16
// mean(f - f_Boost): -1.1972621413014757e+52
// variance(f - f_Boost): 1.4722246255062504e+210
// stddev(f - f_Boost): 1.2133526385623639e+105
const testcase_rising_factorial<double>
data016[8] =
{
  { 1.0000000000000000, 4.0000000000000000, 0.0000000000000000, 0.0 },
  { 4.0000000000000000, 4.0000000000000000, 1.0000000000000000, 0.0 },
  { 20.000000000000000, 4.0000000000000000, 2.0000000000000000, 0.0 },
  { 6719.9999999999991, 4.0000000000000000, 5.0000000000000000, 0.0 },
  { 1037836800.0000000, 4.0000000000000000, 10.000000000000000, 0.0 },
  { 4.3086694564808297e+21, 4.0000000000000000, 20.000000000000000, 0.0 },
  { 7.1248054734333768e+68, 4.0000000000000000, 50.000000000000000, 0.0 },
  { 1.6504834527476967e+163, 4.0000000000000000, 100.00000000000000, 0.0 },
};
const double toler016 = 2.5000000000000020e-13;

// Test data for a=4.2500000000000000.
// max(|f - f_Boost|): 6.3933410310471521e+147 at index 7
// max(|f - f_Boost| / |f_Boost|): 2.7439022836122984e-16
// mean(f - f_Boost): -7.9916762888089401e+146
// variance(f - f_Boost): 5.1093511924088827e+294
// stddev(f - f_Boost): 2.2603873987458175e+147
const testcase_rising_factorial<double>
data017[8] =
{
  { 1.0000000000000000, 4.2500000000000000, 0.0000000000000000, 0.0 },
  { 4.2500000000000000, 4.2500000000000000, 1.0000000000000000, 0.0 },
  { 22.312500000000000, 4.2500000000000000, 2.0000000000000000, 0.0 },
  { 8341.0400390625000, 4.2500000000000000, 5.0000000000000000, 0.0 },
  { 1444076808.7874651, 4.2500000000000000, 10.000000000000000, 0.0 },
  { 6.8793594689003559e+21, 4.2500000000000000, 20.000000000000000, 0.0 },
  { 1.3962737940947973e+69, 4.2500000000000000, 50.000000000000000, 0.0 },
  { 3.8135722611693530e+163, 4.2500000000000000, 100.00000000000000, 0.0 },
};
const double toler017 = 2.5000000000000020e-13;

// Test data for a=4.5000000000000000.
// max(|f - f_Boost|): 1.2786682062094304e+148 at index 7
// max(|f - f_Boost| / |f_Boost|): 2.8437907145135538e-16
// mean(f - f_Boost): 1.5983352577617880e+147
// variance(f - f_Boost): 2.0437404769635531e+295
// stddev(f - f_Boost): 4.5207747974916350e+147
const testcase_rising_factorial<double>
data018[8] =
{
  { 1.0000000000000000, 4.5000000000000000, 0.0000000000000000, 0.0 },
  { 4.5000000000000000, 4.5000000000000000, 1.0000000000000000, 0.0 },
  { 24.749999999999996, 4.5000000000000000, 2.0000000000000000, 0.0 },
  { 10255.781250000000, 4.5000000000000000, 5.0000000000000000, 0.0 },
  { 1985286892.4560544, 4.5000000000000000, 10.000000000000000, 0.0 },
  { 1.0831634819161862e+22, 4.5000000000000000, 20.000000000000000, 0.0 },
  { 2.6944590771828837e+69, 4.5000000000000000, 50.000000000000000, 0.0 },
  { 8.6718719282821878e+163, 4.5000000000000000, 100.00000000000000, 0.0 },
};
const double toler018 = 2.5000000000000020e-13;

// Test data for a=4.7500000000000000.
// max(|f - f_Boost|): 2.5573364124188608e+148 at index 7
// max(|f - f_Boost| / |f_Boost|): 2.9900573540754426e-16
// mean(f - f_Boost): 3.1966705155235760e+147
// variance(f - f_Boost): 8.1749619078542123e+295
// stddev(f - f_Boost): 9.0415495949832700e+147
const testcase_rising_factorial<double>
data019[8] =
{
  { 1.0000000000000000, 4.7500000000000000, 0.0000000000000000, 0.0 },
  { 4.7500000000000000, 4.7500000000000000, 1.0000000000000000, 0.0 },
  { 27.312500000000000, 4.7500000000000000, 2.0000000000000000, 0.0 },
  { 12501.870117187498, 4.7500000000000000, 5.0000000000000000, 0.0 },
  { 2699223287.7444506, 4.7500000000000000, 10.000000000000000, 0.0 },
  { 1.6835084814609839e+22, 4.7500000000000000, 20.000000000000000, 0.0 },
  { 5.1253048332905732e+69, 4.7500000000000000, 50.000000000000000, 0.0 },
  { 1.9426706716849034e+164, 4.7500000000000000, 100.00000000000000, 0.0 },
};
const double toler019 = 2.5000000000000020e-13;

// Test data for a=5.0000000000000000.
// max(|f - f_Boost|): 5.1146728248377217e+148 at index 7
// max(|f - f_Boost| / |f_Boost|): 1.5932812290183794e-16
// mean(f - f_Boost): 6.3933410310471521e+147
// variance(f - f_Boost): 3.2699847631416849e+296
// stddev(f - f_Boost): 1.8083099189966540e+148
const testcase_rising_factorial<double>
data020[8] =
{
  { 1.0000000000000000, 5.0000000000000000, 0.0000000000000000, 0.0 },
  { 5.0000000000000000, 5.0000000000000000, 1.0000000000000000, 0.0 },
  { 30.000000000000000, 5.0000000000000000, 2.0000000000000000, 0.0 },
  { 15119.999999999998, 5.0000000000000000, 5.0000000000000000, 0.0 },
  { 3632428800.0000000, 5.0000000000000000, 10.000000000000000, 0.0 },
  { 2.5852016738884974e+22, 5.0000000000000000, 20.000000000000000, 0.0 },
  { 9.6184873891350585e+69, 5.0000000000000000, 50.000000000000000, 0.0 },
  { 4.2912569771440110e+164, 5.0000000000000000, 100.00000000000000, 0.0 },
};
const double toler020 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  void
  test(const testcase_rising_factorial<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = -Ret(1);
    Ret max_abs_frac = -Ret(1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::rising_factorial(data[i].a, data[i].x);
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
  test(data011, toler011);
  test(data012, toler012);
  test(data013, toler013);
  test(data014, toler014);
  test(data015, toler015);
  test(data016, toler016);
  test(data017, toler017);
  test(data018, toler018);
  test(data019, toler019);
  test(data020, toler020);
  return 0;
}