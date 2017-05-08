// -*- C++ -*- header.

// Copyright (C) 2016-2017 Free Software Foundation, Inc.
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

/** @file bits/float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_BITS_FLOAT128_IO_H
#define _GLIBCXX_BITS_FLOAT128_IO_H 1

#pragma GCC system_header

#if _GLIBCXX_HAVE_FLOAT128_MATH
#if __has_include(<quadmath.h>)

#include <iosfwd>
#include <quadmath.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	       __float128 __x);

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_istream<_CharT, _Traits>&
    operator>>(std::basic_istream<_CharT, _Traits>& __is, __float128& __x);

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#include <bits/float128_io.tcc>

#endif // __has_include(<quadmath.h>)
#endif // _GLIBCXX_HAVE_FLOAT128_MATH

#endif // _GLIBCXX_BITS_FLOAT128_IO_H
