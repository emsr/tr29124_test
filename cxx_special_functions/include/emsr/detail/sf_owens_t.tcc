
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/sf_owens_t.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef SF_OWENS_T_TCC
#define SF_OWENS_T_TCC 1

#include <emsr/math_constants.h>

namespace emsr
{
namespace detail
{

/**
  *
  */
template<typename Tp>
  Tp
  znorm2(Tp x)
  {
    const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
    return Tp{0.5L} * std::erfc(x / s_sqrt_2);
  }

/**
 *
 */
template<typename Tp>
  Tp
  znorm1(Tp x)
  {
    const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
    return Tp{0.5L} * std::erf(x / s_sqrt_2);
  }

/**
 *  The CDF of the normal distribution.
 *  i.e. the integrated lower tail of the normal PDF.
 */
template<typename Tp>
  Tp
  gauss(Tp x)
  {
    const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
    return Tp{0.5L} * (Tp{1} + std::erf(x/s_sqrt_2));
  }

/**
 * Return the Owens T function:
 * @f[
 *  T(h,a) = \frac{1}{2\pi}\int_0^a\frac{\exp[-\frac{1}{2}h^2(1+x^2)]}{1+x^2}dx
 * @f]
 *
 * This implementation is a translation of the Fortran implementation in 
 * @see Patefield, M. and Tandy, D.
 *      "Fast and accurate Calculation of Owen's T-Function",
 *      Journal of Statistical Software, 5 (5), 1 - 25 (2000)
 * @param[in]  h  The scale parameter.
 * @param[in]  a  The integration limit.
 * @return The owens T function.
 */
template<typename Tp>
  Tp
  owens_t(Tp h, Tp a)
  {
    constexpr Tp s_eps = std::numeric_limits<Tp>::epsilon();
    constexpr Tp s_min = std::numeric_limits<Tp>::min();

    constexpr std::size_t s_num_c2 = 21;
    constexpr Tp s_c2[s_num_c2]
    {
       0.99999999999999987510L,
      -0.99999999999988796462L,
       0.99999999998290743652L,
      -0.99999999896282500134L,
       0.99999996660459362918L,
      -0.99999933986272476760L,
       0.99999125611136965852L,
      -0.99991777624463387686L,
       0.99942835555870132569L,
      -0.99697311720723000295L,
       0.98751448037275303682L,
      -0.95915857980572882813L,
       0.89246305511006708555L,
      -0.76893425990463999675L,
       0.58893528468484693250L,
      -0.38380345160440256652L,
       0.20317601701045299653L,
      -0.82813631607004984866e-01L,
       0.24167984735759576523e-01L,
      -0.44676566663971825242e-02L,
       0.39141169402373836468e-03L
    };

    constexpr Tp s_h_range[14]
    {
      0.02L, 0.06L, 0.09L, 0.125L, 0.26L, 0.4L,  0.6L,
      1.6L,  1.7L,  2.33L, 2.4L,   3.36L, 3.4L,  4.8L
    };

    constexpr Tp s_a_range[7]
    {
      0.025L,
      0.09L,
      0.15L,
      0.36L,
      0.5L,
      0.9L,
      0.99999L
    };

    constexpr int s_select[8][15]
    {
      {0,  0,  1, 12, 12, 12, 12, 12, 12, 12, 12, 15, 15, 15,  8},
      {0,  1,  1,  2,  2,  4,  4, 13, 13, 14, 14, 15, 15, 15,  8},
      {1,  1,  2,  2,  2,  4,  4, 14, 14, 14, 14, 15, 15, 15,  9},
      {1,  1,  2,  4,  4,  4,  4,  6,  6, 15, 15, 15, 15, 15,  9},
      {1,  2,  2,  4,  4,  5,  5,  7,  7, 16, 16, 16, 11, 11, 10},
      {1,  2,  4,  4,  4,  5,  5,  7,  7, 16, 16, 16, 11, 11, 11},
      {1,  2,  3,  3,  5,  5,  7,  7, 16, 16, 16, 16, 16, 11, 11},
      {1,  2,  3,  3,  5,  5, 17, 17, 17, 17, 16, 16, 16, 11, 11}
    };

    constexpr int s_method[18]
    {
      1, 1, 1, 1, 1, 1, 1, 1, 2,
      2, 2, 3, 4, 4, 4, 4, 5, 6
    };

    constexpr int s_ord[18]
    {
       2,  3,  4,  5,  7, 10, 12, 18, 10,
      20, 30, 20,  4,  7,  8, 20, 13,  0
    };

    constexpr std::size_t s_num_GJ = 13;

    constexpr Tp s_GJ_pts[s_num_GJ]
    {
      0.35082039676451715489e-02L,
      0.31279042338030753740e-01L,
      0.85266826283219451090e-01L,
      0.16245071730812277011L,
      0.25851196049125434828L,
      0.36807553840697533536L,
      0.48501092905604697475L,
      0.60277514152618576821L,
      0.71477884217753226516L,
      0.81475510988760098605L,
      0.89711029755948965867L,
      0.95723808085944261843L,
      0.99178832974629703586L
    };

    constexpr Tp s_GJ_wts[s_num_GJ]
    {
      0.18831438115323502887e-01L,
      0.18567086243977649478e-01L,
      0.18042093461223385584e-01L,
      0.17263829606398753364e-01L,
      0.16243219975989856730e-01L,
      0.14994592034116704829e-01L,
      0.13535474469662088392e-01L,
      0.11886351605820165233e-01L,
      0.10070377242777431897e-01L,
      0.81130545742299586629e-02L,
      0.60419009528470238773e-02L,
      0.38862217010742057883e-02L,
      0.16793031084546090448e-02L
    };

    constexpr Tp s_1_d_sqrt_2pi = 0.39894228040143267794L;
    constexpr Tp s_1_d_2pi = 0.15915494309189533577L;

    if (std::isnan(a) || std::isnan(h))
      return std::numeric_limits<Tp>::quiet_NaN();
    else if (h < Tp{0})
      return owens_t(-h, a);
    else if (a < Tp{0})
      return -owens_t(h, -a);
    else if (a > Tp{1})
      {
	auto normh = znorm2(h);
	auto normah = znorm2(a * h);
	return Tp{0.5L} * normh + Tp{0.5L} * normah
	     - normh * normah - owens_t(h * a, Tp{1} / a);
      }
    if (h == Tp{0})
      return std::atan(a) * s_1_d_2pi;
    if (a == Tp{0})
      return Tp{0};
    if (a == Tp{1})
      return Tp{0.5L} * znorm2(-h) * znorm2(h);
    else
      {
	//  Determine appropriate method from t1...t6

	auto iaint = 7;
	for (int i = 0; i < 7; ++i)
	  if (a <= s_a_range[i])
	    {
	      iaint = i;
	      break;
	    }

	auto ihint = 14;
	for (int i = 0; i < 14; ++i)
	  if (h <= s_h_range[i])
	    {
	      ihint = i;
	      break;
	    }

	auto icode = s_select[iaint][ihint];
	auto m = s_ord[icode];

	if (s_method[icode] == 1)
	  {
	    //  t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18
	    //  jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)^j / j!
	    //  aj = a^(2j-1) / (2*pi)
	    const auto hh = - Tp{0.5L} * h * h;
	    const auto dhs = std::exp(hh);
	    const auto aa = a * a;
	    auto j = 1;
	    auto jj = Tp{1};
	    auto aj = s_1_d_2pi * a;
	    auto dj = std::expm1(hh);
	    auto gj = hh * dhs;

	    auto value = s_1_d_2pi * std::atan(a);
	    while (true)
	      {
		auto z = dj * aj / jj;
		value += z;

		if (j >= m && std::abs(z) < s_eps * std::abs(value))
        	  return value;

		++j;
		jj += Tp{2};
		aj *= aa;
		dj = gj - dj;
		gj *= hh / Tp(j);
	      }
	  }
	else if (s_method[icode] == 2)
	  {
	    //  t2(h, a, m) ; m = 10, 20 or 30
	    //  z = (-1)^(i-1) * zi ; ii = 2i - 1
	    //  vi = (-1)^(i-1) * a^(2i-1) * exp[-(a*h)^2/2] / sqrt(2*pi)
	    auto maxii = m + m + 1;
	    auto ii = 1;
	    const auto hh = h * h;
	    const auto aa = -a * a;
	    const auto y = Tp{1} / hh;
	    auto ah = a * h;
	    auto vi = s_1_d_sqrt_2pi * a * std::exp(-Tp{0.5} * ah * ah);
	    auto z = znorm1(ah) / h;
	    auto z_prev = std::abs(z) + Tp{1};

	    auto value = Tp{0};
	    while (true)
	      {
		value += z;

		if (std::abs(z) < s_eps * std::abs(value)
		 || (ii >= maxii && std::abs(z) > z_prev)
		 || std::abs(z) < s_eps)
		  {
        	    value *= s_1_d_sqrt_2pi * std::exp(-Tp{0.5} * hh);
        	    return value;
		  }

		z_prev = std::abs(z);
		z = y * (vi - Tp(ii) * z);
		vi *= aa;
		ii += 2;
	      }
	  }
	else if (s_method[icode] == 3)
	  {
	    //  t3(h, a, m) ; m = 20
	    //  ii = 2i - 1
	    //  vi = a^(2i-1) * exp[-(a*h)^2/2] / sqrt(2*pi)
	    const auto hh = h * h;
	    const auto aa = a * a;
	    const auto ah = a * h;
	    const auto y = Tp{1} / hh;
	    auto ii = 1;
	    auto vi = s_1_d_sqrt_2pi * a * std::exp(-Tp{0.5L} * ah * ah);
	    auto zi = znorm1(ah) / h;

	    auto value = Tp{0};
	    for (auto i = 0ull; i < s_num_c2; ++i)
	      {
		value += zi * s_c2[i];
		zi = y  * (Tp(ii) * zi - vi);
		vi *= aa;
		ii += 2;
	      }
            value *= s_1_d_sqrt_2pi * std::exp(-Tp{0.5L} * hh);
            return value;
	  }
	else if (s_method[icode] == 4)
	  {
	    //  t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1
	    //  ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)^i / (2*pi)
	    const auto hh = h * h;
	    const auto aa = -a * a;
	    auto ii = 1;
	    auto ai = s_1_d_2pi * a
		      * std::exp(-Tp{0.5L} * hh * (Tp{1} - aa));
	    auto yi = Tp{1};

	    auto value = ai * yi;
	    while (true)
	      {
		ii += 2;
		yi = (Tp{1} - hh * yi) / Tp(ii);
		ai *= aa;
		auto z = ai * yi;
		value += z;

		//if (maxii <= ii)
		if (std::abs(z) > s_min && std::abs(z) < s_eps * std::abs(value))
        	  return value;
	      }
	  }
	else if (s_method[icode] == 5)
	  {
	    //  t5(h, a, m) ; m = 13
	    //  2m - point Gaussian quadrature
	    const auto aa = a * a;
	    const auto hh = - Tp{0.5L} * h * h;
	    auto value = Tp{0};
	    for (auto i = 0ull; i < s_num_GJ; ++i)
	      {
		auto r = Tp{1} + aa * s_GJ_pts[i];
		value += s_GJ_wts[i] * std::exp(hh * r) / r;
	      }
	    value *= a;
	    return value;
	  }
	else if (s_method[icode] == 6)
	  {
	    //  t6(h, a);  approximation for a near 1, (a<=1)
	    const auto normh = znorm2(h);
	    auto value = Tp{0.5L} * normh * (Tp{1} - normh);
	    const auto y = Tp{1} - a;
	    const auto r = std::atan2(y, Tp{1} + a);

	    if (std::abs(r) > s_eps)
	      value -= s_1_d_2pi * r
		       * std::exp (- Tp{0.5L} * y * h * h / r);
	    return value;
	  }

	return Tp{0};
      }
  }

} // namespace detail
} // namespace emsr

#endif // SF_OWENS_T_TCC
