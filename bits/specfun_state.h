// Special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/specfun_state.h
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SPECFUN_STATE_H
#define _GLIBCXX_BITS_SPECFUN_STATE_H 1

#pragma GCC system_header

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp>
    struct __airy_t
    {
      _Tp x;
      _Tp Ai;
      _Tp Aip;
      _Tp Bi;
      _Tp Bip;
    };

  template<typename _Tp>
    struct __cyl_bessel_t
    {
      _Tp x;
      _Tp nu;
      _Tp J;
      _Tp Jp;
      _Tp N;
      _Tp Np;
    };

  template<typename _Tp>
    struct __cyl_mod_bessel_t
    {
      _Tp x;
      _Tp nu;
      _Tp I;
      _Tp Ip;
      _Tp K;
      _Tp Kp;
    };

  template<typename _Tp>
    struct __cyl_hankel_t
    {
      _Tp x;
      _Tp nu;
      _Tp H1;
      _Tp H1p;
      _Tp H2;
      _Tp H2p;
    };

  template<typename _Tp>
    struct __sph_bessel_t
    {
      _Tp x;
      unsigned int n;
      _Tp j;
      _Tp jp;
      _Tp n;
      _Tp np;
    };

  template<typename _Tp>
    struct __sph_mod_bessel_t
    {
      _Tp x;
      unsigned int n;
      _Tp i;
      _Tp ip;
      _Tp k;
      _Tp kp;
    };

  template<typename _Tp>
    struct __sph_hankel_t
    {
      _Tp x;
      unsigned int n;
      _Tp h1;
      _Tp h1p;
      _Tp h2;
      _Tp h2p;
    };

  template<typename _Tp>
    struct __pqgamma_t
    {
    };

  template<typename _Tp>
    struct __lgamma_t
    {
      _Tp lgamma;
      int sign;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_SPECFUN_STATE_H
