// Special functions -*- C++ -*-

// Copyright (C) 2015
// Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
// USA.
//
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
#ifndef _GLIBCXX_BITS_SF_HYDROGEN_TCC
#define _GLIBCXX_BITS_SF_HYDROGEN_TCC 1

#pragma GCC system_header

#include <complex>

constexpr long double __PI = 3.1415926535897932384626433832795029L;

template <typename _Tp>
  std::complex<_Tp>
  __hydrogen(const unsigned int __n,
             const unsigned int __l, const unsigned int __m,
             const _Tp _Z, const _Tp __r, const _Tp __theta, const _Tp __phi)
  {
    if(__n < 1)
      __throw_domain_error("__hydrogen: level number less than onw");
    if(__l > __n - 1)
      __throw_domain_error("__hydrogen: angular momentum number too large");
    if(_Z <= _Tp(0))
      __throw_domain_error("__hydrogen: non-positive charge");
    if(__r < _Tp(0))
      __throw_domain_error("__hydrogen: negative radius");

    const _Tp _A = _Tp(2) * _Z / __n;

    const _Tp __pre = std::sqrt(_A * _A * _A / (_Tp(2) * __n));
    const _Tp __ln_a = std::lgamma(__n + __l + 1);
    const _Tp __ln_b = std::lgamma(__n - __l);
    const _Tp __ex = std::exp((__ln_b - __ln_a) / _Tp(2));
    const _Tp __norm = __pre * __ex;

    const _Tp __rho = _A * __r;
    const _Tp __ea = std::exp(-__rho / _Tp(2));
    const _Tp __pp = std::pow(__rho, __l);
    const _Tp __lag = std::assoc_laguerre(__n - __l - 1, 2 * __l + 1,
                                               __rho);
    const std::complex<_Tp> __sphh = std::sph_legendre(__l, __m, __theta)
                                   * std::polar(_Tp(1), _Tp(__m) * __phi);

    const std::complex<_Tp> __psi = __norm * __ea * __pp * __lag * __sphh;

    return __psi;
  }

#endif // _GLIBCXX_BITS_SF_HYDROGEN_TCC
