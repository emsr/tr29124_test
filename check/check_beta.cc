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

//  beta
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


// Test data for x=10.000000000000000.
// max(|f - f_GSL|): 7.7707014339346327e-17
// max(|f - f_GSL| / |f_GSL|): 9.7017504401032555e-10
// mean(f - f_GSL): -7.7938051002258226e-18
// variance(f - f_GSL): 7.4991865016938964e-36
// stddev(f - f_GSL): 2.7384642597072353e-18
const testcase_beta<double>
data001[10] =
{
  { 1.0825088224469029e-06, 10.000000000000000, 10.000000000000000, 0.0 },
  { 4.9925087406346778e-09, 10.000000000000000, 20.000000000000000, 0.0 },
  { 1.5729567312509485e-10, 10.000000000000000, 30.000000000000000, 0.0 },
  { 1.2168673582561288e-11, 10.000000000000000, 40.000000000000000, 0.0 },
  { 1.5916380099863291e-12, 10.000000000000000, 50.000000000000000, 0.0 },
  { 2.9408957938463963e-13, 10.000000000000000, 60.000000000000000, 0.0 },
  { 6.9411637980691676e-14, 10.000000000000000, 70.000000000000000, 0.0 },
  { 1.9665612972502651e-14, 10.000000000000000, 80.000000000000000, 0.0 },
  { 6.4187824828154399e-15, 10.000000000000000, 90.000000000000000, 0.0 },
  { 2.3455339739604842e-15, 10.000000000000000, 100.00000000000000, 0.0 },
};
const double toler001 = 5.0000000000000024e-08;

// Test data for x=20.000000000000000.
// max(|f - f_GSL|): 2.5195507868058946e-19
// max(|f - f_GSL| / |f_GSL|): 1.5755574688509410e-10
// mean(f - f_GSL): -2.5184035061512766e-20
// variance(f - f_GSL): 7.8300694071554691e-41
// stddev(f - f_GSL): 8.8487679408805091e-21
const testcase_beta<double>
data002[10] =
{
  { 4.9925087406346778e-09, 20.000000000000000, 10.000000000000000, 0.0 },
  { 7.2544445519248436e-13, 20.000000000000000, 20.000000000000000, 0.0 },
  { 1.7681885473062028e-15, 20.000000000000000, 30.000000000000000, 0.0 },
  { 1.7891885039182335e-17, 20.000000000000000, 40.000000000000000, 0.0 },
  { 4.3240677875623635e-19, 20.000000000000000, 50.000000000000000, 0.0 },
  { 1.8857342309689050e-20, 20.000000000000000, 60.000000000000000, 0.0 },
  { 1.2609804003539998e-21, 20.000000000000000, 70.000000000000000, 0.0 },
  { 1.1660809542079041e-22, 20.000000000000000, 80.000000000000000, 0.0 },
  { 1.3907944279729071e-23, 20.000000000000000, 90.000000000000000, 0.0 },
  { 2.0365059099917614e-24, 20.000000000000000, 100.00000000000000, 0.0 },
};
const double toler002 = 1.0000000000000005e-08;

// Test data for x=30.000000000000000.
// max(|f - f_GSL|): 2.1430776255782837e-20
// max(|f - f_GSL| / |f_GSL|): 1.3624517337320060e-10
// mean(f - f_GSL): 2.1431216846143432e-21
// variance(f - f_GSL): 5.6703340185978029e-43
// stddev(f - f_GSL): 7.5301620291981785e-22
const testcase_beta<double>
data003[10] =
{
  { 1.5729567312509485e-10, 30.000000000000000, 10.000000000000000, 0.0 },
  { 1.7681885473062028e-15, 30.000000000000000, 20.000000000000000, 0.0 },
  { 5.6370779640482451e-19, 30.000000000000000, 30.000000000000000, 0.0 },
  { 1.0539424603796547e-21, 30.000000000000000, 40.000000000000000, 0.0 },
  { 6.0118197777273836e-24, 30.000000000000000, 50.000000000000000, 0.0 },
  { 7.4279528553260165e-26, 30.000000000000000, 60.000000000000000, 0.0 },
  { 1.6212207780604767e-27, 30.000000000000000, 70.000000000000000, 0.0 },
  { 5.4783729715317616e-29, 30.000000000000000, 80.000000000000000, 0.0 },
  { 2.6183005659681346e-30, 30.000000000000000, 90.000000000000000, 0.0 },
  { 1.6587948222122229e-31, 30.000000000000000, 100.00000000000000, 0.0 },
};
const double toler003 = 1.0000000000000005e-08;

// Test data for x=40.000000000000000.
// max(|f - f_GSL|): 5.0244436746634082e-22
// max(|f - f_GSL| / |f_GSL|): 4.1289986460511602e-11
// mean(f - f_GSL): 5.0243379525418262e-23
// variance(f - f_GSL): 3.1165397359694066e-46
// stddev(f - f_GSL): 1.7653724071621280e-23
const testcase_beta<double>
data004[10] =
{
  { 1.2168673582561288e-11, 40.000000000000000, 10.000000000000000, 0.0 },
  { 1.7891885039182335e-17, 40.000000000000000, 20.000000000000000, 0.0 },
  { 1.0539424603796547e-21, 40.000000000000000, 30.000000000000000, 0.0 },
  { 4.6508509140090659e-25, 40.000000000000000, 40.000000000000000, 0.0 },
  { 7.5161712118557719e-28, 40.000000000000000, 50.000000000000000, 0.0 },
  { 3.0311331979886071e-30, 40.000000000000000, 60.000000000000000, 0.0 },
  { 2.4175035070466313e-32, 40.000000000000000, 70.000000000000000, 0.0 },
  { 3.2734839142758369e-34, 40.000000000000000, 80.000000000000000, 0.0 },
  { 6.7690629601315579e-36, 40.000000000000000, 90.000000000000000, 0.0 },
  { 1.9797337118812366e-37, 40.000000000000000, 100.00000000000000, 0.0 },
};
const double toler004 = 2.5000000000000013e-09;

// Test data for x=50.000000000000000.
// max(|f - f_GSL|): 1.1204316897292397e-21
// max(|f - f_GSL| / |f_GSL|): 7.0394881417720306e-10
// mean(f - f_GSL): -1.1204318293025886e-22
// variance(f - f_GSL): 1.5498364001411668e-45
// stddev(f - f_GSL): 3.9367961594946298e-23
const testcase_beta<double>
data005[10] =
{
  { 1.5916380099863291e-12, 50.000000000000000, 10.000000000000000, 0.0 },
  { 4.3240677875623635e-19, 50.000000000000000, 20.000000000000000, 0.0 },
  { 6.0118197777273836e-24, 50.000000000000000, 30.000000000000000, 0.0 },
  { 7.5161712118557719e-28, 50.000000000000000, 40.000000000000000, 0.0 },
  { 3.9646612085674138e-31, 50.000000000000000, 50.000000000000000, 0.0 },
  { 5.8425643906418403e-34, 50.000000000000000, 60.000000000000000, 0.0 },
  { 1.8672362180783552e-36, 50.000000000000000, 70.000000000000000, 0.0 },
  { 1.0939382296458962e-38, 50.000000000000000, 80.000000000000000, 0.0 },
  { 1.0442781609881063e-40, 50.000000000000000, 90.000000000000000, 0.0 },
  { 1.4904121110954370e-42, 50.000000000000000, 100.00000000000000, 0.0 },
};
const double toler005 = 5.0000000000000024e-08;

// Test data for x=60.000000000000000.
// max(|f - f_GSL|): 9.0985727438802444e-23
// max(|f - f_GSL| / |f_GSL|): 3.0938099754905715e-10
// mean(f - f_GSL): 9.0985743683683168e-24
// variance(f - f_GSL): 1.0220253769966532e-47
// stddev(f - f_GSL): 3.1969131627190834e-24
const testcase_beta<double>
data006[10] =
{
  { 2.9408957938463963e-13, 60.000000000000000, 10.000000000000000, 0.0 },
  { 1.8857342309689050e-20, 60.000000000000000, 20.000000000000000, 0.0 },
  { 7.4279528553260165e-26, 60.000000000000000, 30.000000000000000, 0.0 },
  { 3.0311331979886071e-30, 60.000000000000000, 40.000000000000000, 0.0 },
  { 5.8425643906418403e-34, 60.000000000000000, 50.000000000000000, 0.0 },
  { 3.4501231469782229e-37, 60.000000000000000, 60.000000000000000, 0.0 },
  { 4.7706855386086599e-40, 60.000000000000000, 70.000000000000000, 0.0 },
  { 1.2902663809722593e-42, 60.000000000000000, 80.000000000000000, 0.0 },
  { 6.0105571058570508e-45, 60.000000000000000, 90.000000000000000, 0.0 },
  { 4.3922898898347209e-47, 60.000000000000000, 100.00000000000000, 0.0 },
};
const double toler006 = 2.5000000000000012e-08;

// Test data for x=70.000000000000000.
// max(|f - f_GSL|): 3.3337312341737893e-23
// max(|f - f_GSL| / |f_GSL|): 4.8028419025367724e-10
// mean(f - f_GSL): 3.3337311814398949e-24
// variance(f - f_GSL): 1.3720695790252640e-48
// stddev(f - f_GSL): 1.1713537377860130e-24
const testcase_beta<double>
data007[10] =
{
  { 6.9411637980691676e-14, 70.000000000000000, 10.000000000000000, 0.0 },
  { 1.2609804003539998e-21, 70.000000000000000, 20.000000000000000, 0.0 },
  { 1.6212207780604767e-27, 70.000000000000000, 30.000000000000000, 0.0 },
  { 2.4175035070466313e-32, 70.000000000000000, 40.000000000000000, 0.0 },
  { 1.8672362180783552e-36, 70.000000000000000, 50.000000000000000, 0.0 },
  { 4.7706855386086599e-40, 70.000000000000000, 60.000000000000000, 0.0 },
  { 3.0453137143486369e-43, 70.000000000000000, 70.000000000000000, 0.0 },
  { 4.0192274082013779e-46, 70.000000000000000, 80.000000000000000, 0.0 },
  { 9.5865870063501807e-49, 70.000000000000000, 90.000000000000000, 0.0 },
  { 3.7409127305819802e-51, 70.000000000000000, 100.00000000000000, 0.0 },
};
const double toler007 = 2.5000000000000012e-08;

// Test data for x=80.000000000000000.
// max(|f - f_GSL|): 1.9079086931087787e-23
// max(|f - f_GSL| / |f_GSL|): 9.7017504401032555e-10
// mean(f - f_GSL): -1.9079087043730904e-24
// variance(f - f_GSL): 4.4939699064476602e-49
// stddev(f - f_GSL): 6.7037078594220229e-25
const testcase_beta<double>
data008[10] =
{
  { 1.9665612972502651e-14, 80.000000000000000, 10.000000000000000, 0.0 },
  { 1.1660809542079041e-22, 80.000000000000000, 20.000000000000000, 0.0 },
  { 5.4783729715317616e-29, 80.000000000000000, 30.000000000000000, 0.0 },
  { 3.2734839142758369e-34, 80.000000000000000, 40.000000000000000, 0.0 },
  { 1.0939382296458962e-38, 80.000000000000000, 50.000000000000000, 0.0 },
  { 1.2902663809722593e-42, 80.000000000000000, 60.000000000000000, 0.0 },
  { 4.0192274082013779e-46, 80.000000000000000, 70.000000000000000, 0.0 },
  { 2.7160590828669411e-49, 80.000000000000000, 80.000000000000000, 0.0 },
  { 3.4593773902125368e-52, 80.000000000000000, 90.000000000000000, 0.0 },
  { 7.4807039968503468e-55, 80.000000000000000, 100.00000000000000, 0.0 },
};
const double toler008 = 5.0000000000000024e-08;

// Test data for x=90.000000000000000.
// max(|f - f_GSL|): 4.3336941575312027e-25
// max(|f - f_GSL| / |f_GSL|): 6.7515828260787790e-11
// mean(f - f_GSL): -4.3336941085471556e-26
// variance(f - f_GSL): 2.3186302007970768e-52
// stddev(f - f_GSL): 1.5227048961624432e-26
const testcase_beta<double>
data009[10] =
{
  { 6.4187824828154399e-15, 90.000000000000000, 10.000000000000000, 0.0 },
  { 1.3907944279729071e-23, 90.000000000000000, 20.000000000000000, 0.0 },
  { 2.6183005659681346e-30, 90.000000000000000, 30.000000000000000, 0.0 },
  { 6.7690629601315579e-36, 90.000000000000000, 40.000000000000000, 0.0 },
  { 1.0442781609881063e-40, 90.000000000000000, 50.000000000000000, 0.0 },
  { 6.0105571058570508e-45, 90.000000000000000, 60.000000000000000, 0.0 },
  { 9.5865870063501807e-49, 90.000000000000000, 70.000000000000000, 0.0 },
  { 3.4593773902125368e-52, 90.000000000000000, 80.000000000000000, 0.0 },
  { 2.4416737907558032e-55, 90.000000000000000, 90.000000000000000, 0.0 },
  { 3.0238531916564246e-58, 90.000000000000000, 100.00000000000000, 0.0 },
};
const double toler009 = 5.0000000000000026e-09;

// Test data for x=100.00000000000000.
// max(|f - f_GSL|): 8.1625173357209578e-25
// max(|f - f_GSL| / |f_GSL|): 3.4800251995234896e-10
// mean(f - f_GSL): 8.1625173533399292e-26
// variance(f - f_GSL): 8.2255172276019113e-52
// stddev(f - f_GSL): 2.8680162530226204e-26
const testcase_beta<double>
data010[10] =
{
  { 2.3455339739604842e-15, 100.00000000000000, 10.000000000000000, 0.0 },
  { 2.0365059099917614e-24, 100.00000000000000, 20.000000000000000, 0.0 },
  { 1.6587948222122229e-31, 100.00000000000000, 30.000000000000000, 0.0 },
  { 1.9797337118812366e-37, 100.00000000000000, 40.000000000000000, 0.0 },
  { 1.4904121110954370e-42, 100.00000000000000, 50.000000000000000, 0.0 },
  { 4.3922898898347209e-47, 100.00000000000000, 60.000000000000000, 0.0 },
  { 3.7409127305819802e-51, 100.00000000000000, 70.000000000000000, 0.0 },
  { 7.4807039968503468e-55, 100.00000000000000, 80.000000000000000, 0.0 },
  { 3.0238531916564246e-58, 100.00000000000000, 90.000000000000000, 0.0 },
  { 2.2087606931991853e-61, 100.00000000000000, 100.00000000000000, 0.0 },
};
const double toler010 = 2.5000000000000012e-08;

template<typename Ret, unsigned int Num>
  void
  test(const testcase_beta<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = -Ret(1);
    Ret max_abs_frac = -Ret(1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = std::beta(data[i].x, data[i].y);
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