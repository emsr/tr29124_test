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

//  gamma_q
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
// max(|f - f_GSL|): 5.5511151231257827e-16 at index 1
// max(|f - f_GSL| / |f_GSL|): 2.5538896158483858e-15
// mean(f - f_GSL): 6.4657875013680993e-17
// variance(f - f_GSL): 4.4761525423739437e-34
// stddev(f - f_GSL): 2.1156919771965727e-17
const testcase_gamma_q<double>
data001[11] =
{
  { 1.0000000000000000, 0.50000000000000000, 0.0000000000000000, 0.0 },
  { 0.31731050786291382, 0.50000000000000000, 0.50000000000000000, 0.0 },
  { 0.15729920705028497, 0.50000000000000000, 1.0000000000000000, 0.0 },
  { 0.083264516663550406, 0.50000000000000000, 1.5000000000000000, 0.0 },
  { 0.045500263896358389, 0.50000000000000000, 2.0000000000000000, 0.0 },
  { 0.025347318677468242, 0.50000000000000000, 2.5000000000000000, 0.0 },
  { 0.014305878435429610, 0.50000000000000000, 3.0000000000000000, 0.0 },
  { 0.0081509715935026879, 0.50000000000000000, 3.5000000000000000, 0.0 },
  { 0.0046777349810472567, 0.50000000000000000, 4.0000000000000000, 0.0 },
  { 0.0026997960632601883, 0.50000000000000000, 4.5000000000000000, 0.0 },
  { 0.0015654022580025473, 0.50000000000000000, 5.0000000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

// Test data for a=1.0000000000000000.
// max(|f - f_GSL|): 8.4099394115355608e-15 at index 3
// max(|f - f_GSL| / |f_GSL|): 3.7690733542884264e-14
// mean(f - f_GSL): 1.4355625274378977e-15
// variance(f - f_GSL): 2.4971136024944202e-31
// stddev(f - f_GSL): 4.9971127688840684e-16
const testcase_gamma_q<double>
data002[11] =
{
  { 1.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.60653065971263342, 1.0000000000000000, 0.50000000000000000, 0.0 },
  { 0.36787944117144228, 1.0000000000000000, 1.0000000000000000, 0.0 },
  { 0.22313016014842982, 1.0000000000000000, 1.5000000000000000, 0.0 },
  { 0.13533528323661270, 1.0000000000000000, 2.0000000000000000, 0.0 },
  { 0.082084998623898800, 1.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.049787068367863931, 1.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.030197383422318501, 1.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.018315638888734186, 1.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.011108996538242306, 1.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.0067379469990854679, 1.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler002 = 2.5000000000000015e-12;

// Test data for a=1.5000000000000000.
// max(|f - f_GSL|): 5.5511151231257827e-16 at index 4
// max(|f - f_GSL| / |f_GSL|): 2.7355139864007151e-15
// mean(f - f_GSL): 1.5297107015431845e-16
// variance(f - f_GSL): 2.5740163134573157e-33
// stddev(f - f_GSL): 5.0734764348100755e-17
const testcase_gamma_q<double>
data003[11] =
{
  { 1.0000000000000000, 1.5000000000000000, 0.0000000000000000, 0.0 },
  { 0.80125195690120099, 1.5000000000000000, 0.50000000000000000, 0.0 },
  { 0.57240670447087971, 1.5000000000000000, 1.0000000000000000, 0.0 },
  { 0.39162517627108873, 1.5000000000000000, 1.5000000000000000, 0.0 },
  { 0.26146412994911028, 1.5000000000000000, 2.0000000000000000, 0.0 },
  { 0.17179714429673318, 1.5000000000000000, 2.5000000000000000, 0.0 },
  { 0.11161022509471247, 1.5000000000000000, 3.0000000000000000, 0.0 },
  { 0.071897772496465076, 1.5000000000000000, 3.5000000000000000, 0.0 },
  { 0.046011705689231353, 1.5000000000000000, 4.0000000000000000, 0.0 },
  { 0.029290886534888223, 1.5000000000000000, 4.5000000000000000, 0.0 },
  { 0.018566135463043237, 1.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler003 = 2.5000000000000020e-13;

// Test data for a=2.0000000000000000.
// max(|f - f_GSL|): 7.7715611723760958e-15 at index 5
// max(|f - f_GSL| / |f_GSL|): 2.7050570585059849e-14
// mean(f - f_GSL): 1.6861512186494565e-15
// variance(f - f_GSL): 4.7985074195142691e-31
// stddev(f - f_GSL): 6.9271259693427471e-16
const testcase_gamma_q<double>
data004[11] =
{
  { 1.0000000000000000, 2.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.90979598956895025, 2.0000000000000000, 0.50000000000000000, 0.0 },
  { 0.73575888234288489, 2.0000000000000000, 1.0000000000000000, 0.0 },
  { 0.55782540037107398, 2.0000000000000000, 1.5000000000000000, 0.0 },
  { 0.40600584970983766, 2.0000000000000000, 2.0000000000000000, 0.0 },
  { 0.28729749518364556, 2.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.19914827347145547, 2.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.13588822540043313, 2.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.091578194443670810, 2.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.061099480960332637, 2.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.040427681994512764, 2.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler004 = 2.5000000000000015e-12;

// Test data for a=2.5000000000000000.
// max(|f - f_GSL|): 1.1657341758564144e-15 at index 6
// max(|f - f_GSL| / |f_GSL|): 3.8068653037403819e-15
// mean(f - f_GSL): 1.9681226345627775e-16
// variance(f - f_GSL): 3.1438180074015644e-33
// stddev(f - f_GSL): 5.6069760186767021e-17
const testcase_gamma_q<double>
data005[11] =
{
  { 1.0000000000000000, 2.5000000000000000, 0.0000000000000000, 0.0 },
  { 0.96256577324729642, 2.5000000000000000, 0.50000000000000000, 0.0 },
  { 0.84914503608460978, 2.5000000000000000, 1.0000000000000000, 0.0 },
  { 0.69998583587862673, 2.5000000000000000, 1.5000000000000000, 0.0 },
  { 0.54941595135277999, 2.5000000000000000, 2.0000000000000000, 0.0 },
  { 0.41588018699550811, 2.5000000000000000, 2.5000000000000000, 0.0 },
  { 0.30621891841327792, 2.5000000000000000, 3.0000000000000000, 0.0 },
  { 0.22064030793671077, 2.5000000000000000, 3.5000000000000000, 0.0 },
  { 0.15623562757772227, 2.5000000000000000, 4.0000000000000000, 0.0 },
  { 0.10906415794977235, 2.5000000000000000, 4.5000000000000000, 0.0 },
  { 0.075235246146512030, 2.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler005 = 2.5000000000000020e-13;

// Test data for a=3.0000000000000000.
// max(|f - f_GSL|): 7.8270723236073536e-15 at index 7
// max(|f - f_GSL| / |f_GSL|): 2.4395015295024001e-14
// mean(f - f_GSL): 2.2557713273065679e-15
// variance(f - f_GSL): 1.1853930726544654e-30
// stddev(f - f_GSL): 1.0887575821340881e-15
const testcase_gamma_q<double>
data006[11] =
{
  { 1.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.98561232203302929, 3.0000000000000000, 0.50000000000000000, 0.0 },
  { 0.91969860292860595, 3.0000000000000000, 1.0000000000000000, 0.0 },
  { 0.80884683053805861, 3.0000000000000000, 1.5000000000000000, 0.0 },
  { 0.67667641618306196, 3.0000000000000000, 2.0000000000000000, 0.0 },
  { 0.54381311588332815, 3.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.42319008112684248, 3.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.32084719886213348, 3.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.23810330555354370, 3.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.17357807091003563, 3.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.12465201948308077, 3.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler006 = 2.5000000000000015e-12;

// Test data for a=3.5000000000000000.
// max(|f - f_GSL|): 1.5543122344752192e-15 at index 4
// max(|f - f_GSL| / |f_GSL|): 2.3366517280200655e-15
// mean(f - f_GSL): 3.9110129276568017e-16
// variance(f - f_GSL): 1.0084869526973163e-34
// stddev(f - f_GSL): 1.0042345108077676e-17
const testcase_gamma_q<double>
data007[11] =
{
  { 1.0000000000000000, 3.5000000000000000, 0.0000000000000000, 0.0 },
  { 0.99482853651651548, 3.5000000000000000, 0.50000000000000000, 0.0 },
  { 0.95984036873010159, 3.5000000000000000, 1.0000000000000000, 0.0 },
  { 0.88500223164315073, 3.5000000000000000, 1.5000000000000000, 0.0 },
  { 0.77977740847571453, 3.5000000000000000, 2.0000000000000000, 0.0 },
  { 0.65996322969428201, 3.5000000000000000, 2.5000000000000000, 0.0 },
  { 0.53974935039555672, 3.5000000000000000, 3.0000000000000000, 0.0 },
  { 0.42887985755305447, 3.5000000000000000, 3.5000000000000000, 0.0 },
  { 0.33259390259930766, 3.5000000000000000, 4.0000000000000000, 0.0 },
  { 0.25265604649656398, 3.5000000000000000, 4.5000000000000000, 0.0 },
  { 0.18857346751344978, 3.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler007 = 2.5000000000000020e-13;

// Test data for a=4.0000000000000000.
// max(|f - f_GSL|): 6.3837823915946501e-15 at index 9
// max(|f - f_GSL| / |f_GSL|): 1.8649891366754866e-14
// mean(f - f_GSL): 2.1144702150815483e-15
// variance(f - f_GSL): 2.7511552083109249e-30
// stddev(f - f_GSL): 1.6586606670174961e-15
const testcase_gamma_q<double>
data008[11] =
{
  { 1.0000000000000000, 4.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.99824837744370920, 4.0000000000000000, 0.50000000000000000, 0.0 },
  { 0.98101184312384615, 4.0000000000000000, 1.0000000000000000, 0.0 },
  { 0.93435754562154993, 4.0000000000000000, 1.5000000000000000, 0.0 },
  { 0.85712346049854693, 4.0000000000000000, 2.0000000000000000, 0.0 },
  { 0.75757613313306615, 4.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.64723188878223115, 4.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.53663266790078512, 4.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.43347012036670907, 4.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.34229595583459133, 4.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.26502591529736169, 4.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler008 = 1.0000000000000008e-12;

// Test data for a=4.5000000000000000.
// max(|f - f_GSL|): 1.6653345369377348e-15 at index 5
// max(|f - f_GSL| / |f_GSL|): 3.8012092456586305e-15
// mean(f - f_GSL): 5.5006504401882751e-16
// variance(f - f_GSL): 6.7302497329313948e-32
// stddev(f - f_GSL): 2.5942724862533997e-16
const testcase_gamma_q<double>
data009[11] =
{
  { 1.0000000000000000, 4.5000000000000000, 0.0000000000000000, 0.0 },
  { 0.99943750269783249, 4.5000000000000000, 0.50000000000000000, 0.0 },
  { 0.99146760662881350, 4.5000000000000000, 1.0000000000000000, 0.0 },
  { 0.96429497268508912, 4.5000000000000000, 1.5000000000000000, 0.0 },
  { 0.91141252683167928, 4.5000000000000000, 2.0000000000000000, 0.0 },
  { 0.83430826019340587, 4.5000000000000000, 2.5000000000000000, 0.0 },
  { 0.73991829209465287, 4.5000000000000000, 3.0000000000000000, 0.0 },
  { 0.63711940716939830, 4.5000000000000000, 3.5000000000000000, 0.0 },
  { 0.53414621690969011, 4.5000000000000000, 4.0000000000000000, 0.0 },
  { 0.43727418891386710, 4.5000000000000000, 4.5000000000000000, 0.0 },
  { 0.35048521232336094, 4.5000000000000000, 5.0000000000000000, 0.0 },
};
const double toler009 = 2.5000000000000020e-13;

// Test data for a=5.0000000000000000.
// max(|f - f_GSL|): 7.5495165674510645e-15 at index 10
// max(|f - f_GSL| / |f_GSL|): 1.7138777873386710e-14
// mean(f - f_GSL): 2.7149999420378827e-15
// variance(f - f_GSL): 2.5709806101536105e-30
// stddev(f - f_GSL): 1.6034277689230689e-15
const testcase_gamma_q<double>
data010[11] =
{
  { 1.0000000000000000, 5.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.99982788437004411, 5.0000000000000000, 0.50000000000000000, 0.0 },
  { 0.99634015317265634, 5.0000000000000000, 1.0000000000000000, 0.0 },
  { 0.98142406377785940, 5.0000000000000000, 1.5000000000000000, 0.0 },
  { 0.94734698265628903, 5.0000000000000000, 2.0000000000000000, 0.0 },
  { 0.89117801891415149, 5.0000000000000000, 2.5000000000000000, 0.0 },
  { 0.81526324452376875, 5.0000000000000000, 3.0000000000000000, 0.0 },
  { 0.72544495330960301, 5.0000000000000000, 3.5000000000000000, 0.0 },
  { 0.62883693517987171, 5.0000000000000000, 4.0000000000000000, 0.0 },
  { 0.53210357637471473, 5.0000000000000000, 4.5000000000000000, 0.0 },
  { 0.44049328506521113, 5.0000000000000000, 5.0000000000000000, 0.0 },
};
const double toler010 = 1.0000000000000008e-12;

template<typename Ret, unsigned int Num>
  void
  test(const testcase_gamma_q<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = -Ret(1);
    Ret max_abs_frac = -Ret(1);
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::gamma_q(data[i].a, data[i].x);
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