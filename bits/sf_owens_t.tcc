// TR29124 math special functions -*- C++ -*-

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

#ifndef _GLIBCXX_BITS_SF_OWENS_T_TCC
#define _GLIBCXX_BITS_SF_OWENS_T_TCC 1

#pragma GCC system_header

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

/**
 *
 */
template<typename _Tp>
  _Tp
  __znorm2(_Tp __x)
  {
    constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
    return _Tp{0.5L} * std::erfc(__x / _S_sqrt_2);
  }

/**
 *
 */
template<typename _Tp>
  _Tp
  __znorm1(_Tp __x)
  {
    constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
    return _Tp{0.5L} * std::erf(__x / _S_sqrt_2);
  }

/**
 *  The CDF of the normal distribution.
 *  i.e. the integrated lower tail of the normal PDF.
 */
template<typename _Tp>
  _Tp
  __gauss(_Tp __x)
  {
    constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
    return _Tp{0.5L} * (_Tp{1} + std::erf(__x/_S_sqrt_2));
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
 * @param[in]  __h  The scale parameter.
 * @param[in]  __a  The integration limit.
 * @return The owens T function.
 */
template<typename _Tp>
  _Tp
  __owens_t(_Tp __h, _Tp __a)
  {
    constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();
    constexpr _Tp _S_min = std::numeric_limits<_Tp>::min();

    constexpr std::size_t _S_num_c2 = 21;
    constexpr _Tp _S_c2[_S_num_c2]
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

    constexpr _Tp _S_h_range[14]
    {
      0.02L, 0.06L, 0.09L, 0.125L, 0.26L, 0.4L,  0.6L,
      1.6L,  1.7L,  2.33L, 2.4L,   3.36L, 3.4L,  4.8L
    };

    constexpr _Tp _S_a_range[7]
    {
      0.025L,
      0.09L,
      0.15L,
      0.36L,
      0.5L,
      0.9L,
      0.99999L
    };

    constexpr int _S_select[8][15]
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

    constexpr int _S_method[18]
    {
      1, 1, 1, 1, 1, 1, 1, 1, 2,
      2, 2, 3, 4, 4, 4, 4, 5, 6
    };

    constexpr int _S_ord[18]
    {
       2,  3,  4,  5,  7, 10, 12, 18, 10,
      20, 30, 20,  4,  7,  8, 20, 13,  0
    };

    constexpr std::size_t _S_num_GJ = 13;

    constexpr _Tp _S_GJ_pts[_S_num_GJ]
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

    constexpr _Tp _S_GJ_wts[_S_num_GJ]
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

    constexpr _Tp _S_1_d_sqrt_2pi = 0.39894228040143267794L;
    constexpr _Tp _S_1_d_2pi = 0.15915494309189533577L;

    if (__isnan(__a) || __isnan(__h))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__h < _Tp{0})
      return __owens_t(-__h, __a);
    else if (__a < _Tp{0})
      return -__owens_t(__h, -__a);
    else if (__a > _Tp{1})
      {
	auto __normh = __znorm2(__h);
	auto __normah = __znorm2(__a * __h);
	return _Tp{0.5L} * __normh + _Tp{0.5L} * __normah
	     - __normh * __normah - __owens_t(__h * __a, _Tp{1} / __a);
      }
    if (__h == _Tp{0})
      return std::atan(__a) * _S_1_d_2pi;
    if (__a == _Tp{0})
      return _Tp{0};
    if (__a == _Tp{1})
      return _Tp{0.5L} * __znorm2(-__h) * __znorm2(__h);
    else
      {
	//  Determine appropriate method from t1...t6

	auto __iaint = 7;
	for (int __i = 0; __i < 7; ++__i)
	  if (__a <= _S_a_range[__i])
	    {
	      __iaint = __i;
	      break;
	    }

	auto __ihint = 14;
	for (int __i = 0; __i < 14; ++__i)
	  if (__h <= _S_h_range[__i])
	    {
	      __ihint = __i;
	      break;
	    }

	auto __icode = _S_select[__iaint][__ihint];
	auto __m = _S_ord[__icode];

	if (_S_method[__icode] == 1)
	  {
	    //  t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18
	    //  jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)^j / j!
	    //  aj = a^(2j-1) / (2*pi)
	    const auto __hh = - _Tp{0.5L} * __h * __h;
	    const auto __dhs = std::exp(__hh);
	    const auto __aa = __a * __a;
	    auto __j = 1;
	    auto __jj = _Tp{1};
	    auto __aj = _S_1_d_2pi * __a;
	    auto __dj = std::expm1(__hh);
	    auto __gj = __hh * __dhs;

	    auto __value = _S_1_d_2pi * std::atan(__a);
	    while (true)
	      {
		auto __z = __dj * __aj / __jj;
		__value += __z;

		if (__j >= __m && std::abs(__z) < _S_eps * std::abs(__value))
        	  return __value;

		++__j;
		__jj += _Tp{2};
		__aj *= __aa;
		__dj = __gj - __dj;
		__gj *= __hh / _Tp(__j);
	      }
	  }
	else if (_S_method[__icode] == 2)
	  {
	    //  t2(h, a, m) ; m = 10, 20 or 30
	    //  z = (-1)^(i-1) * zi ; ii = 2i - 1
	    //  vi = (-1)^(i-1) * a^(2i-1) * exp[-(a*h)^2/2] / sqrt(2*pi)
	    auto __maxii = __m + __m + 1;
	    auto __ii = 1;
	    const auto __hh = __h * __h;
	    const auto __aa = -__a * __a;
	    const auto __y = _Tp{1} / __hh;
	    auto __ah = __a * __h;
	    auto __vi = _S_1_d_sqrt_2pi * __a * std::exp(-_Tp{0.5} * __ah * __ah);
	    auto __z = __znorm1(__ah) / __h;
	    auto __z_prev = std::abs(__z) + _Tp{1};

	    auto __value = _Tp{0};
	    while (true)
	      {
		__value += __z;

		if (std::abs(__z) < _S_eps * std::abs(__value)
		 || (__ii >= __maxii && std::abs(__z) > __z_prev)
		 || std::abs(__z) < _S_eps)
		  {
        	    __value *= _S_1_d_sqrt_2pi * std::exp(-_Tp{0.5} * __hh);
        	    return __value;
		  }

		__z_prev = std::abs(__z);
		__z = __y * (__vi - _Tp(__ii) * __z);
		__vi *= __aa;
		__ii += 2;
	      }
	  }
	else if (_S_method[__icode] == 3)
	  {
	    //  t3(h, a, m) ; m = 20
	    //  ii = 2i - 1
	    //  vi = a^(2i-1) * exp[-(a*h)^2/2] / sqrt(2*pi)
	    const auto __hh = __h * __h;
	    const auto __aa = __a * __a;
	    const auto __ah = __a * __h;
	    const auto __y = _Tp{1} / __hh;
	    auto __ii = 1;
	    auto __vi = _S_1_d_sqrt_2pi * __a * std::exp(-_Tp{0.5L} * __ah * __ah);
	    auto __zi = __znorm1(__ah) / __h;

	    auto __value = _Tp{0};
	    for (int __i = 0; __i < _S_num_c2; ++__i)
	      {
		__value += __zi * _S_c2[__i];
		__zi = __y  * (_Tp(__ii) * __zi - __vi);
		__vi *= __aa;
		__ii += 2;
	      }
            __value *= _S_1_d_sqrt_2pi * std::exp(-_Tp{0.5L} * __hh);
            return __value;
	  }
	else if (_S_method[__icode] == 4)
	  {
	    //  t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1
	    //  ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)^i / (2*pi)
	    const auto __maxii = __m + __m + 1;
	    const auto __hh = __h * __h;
	    const auto __aa = -__a * __a;
	    auto __ii = 1;
	    auto __ai = _S_1_d_2pi * __a
		      * std::exp(-_Tp{0.5L} * __hh * (_Tp{1} - __aa));
	    auto __yi = _Tp{1};

	    auto __value = __ai * __yi;
	    while (true)
	      {
		__ii += 2;
		__yi = (_Tp{1} - __hh * __yi) / _Tp(__ii);
		__ai *= __aa;
		auto __z = __ai * __yi;
		__value += __z;

		//if (__maxii <= __ii)
		if (std::abs(__z) > _S_min && std::abs(__z) < _S_eps * std::abs(__value))
        	  return __value;
	      }
	  }
	else if (_S_method[__icode] == 5)
	  {
	    //  t5(h, a, m) ; m = 13
	    //  2m - point Gaussian quadrature
	    const auto __aa = __a * __a;
	    const auto __hh = - _Tp{0.5L} * __h * __h;
	    auto __value = _Tp{0};
	    for (int __i = 0; __i < _S_num_GJ; ++__i)
	      {
		auto __r = _Tp{1} + __aa * _S_GJ_pts[__i];
		__value += _S_GJ_wts[__i] * std::exp(__hh * __r) / __r;
	      }
	    __value *= __a;
	    return __value;
	  }
	else if (_S_method[__icode] == 6)
	  {
	    //  t6(h, a);  approximation for a near 1, (a<=1)
	    const auto __normh = __znorm2(__h);
	    auto __value = _Tp{0.5L} * __normh * (_Tp{1} - __normh);
	    const auto __y = _Tp{1} - __a;
	    const auto __r = std::atan2(__y, _Tp{1} + __a);

	    if (std::abs(__r) > _S_eps)
	      __value -= _S_1_d_2pi * __r
		       * std::exp (- _Tp{0.5L} * __y * __h * __h / __r);
	    return __value;
	  }

	return _Tp{0};
      }
  }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_OWENS_T_TCC
