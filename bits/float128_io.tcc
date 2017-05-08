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

/** @file bits/float128.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_BITS_FLOAT128_IO_TCC
#define _GLIBCXX_BITS_FLOAT128_IO_TCC 1

#pragma GCC system_header

#if _GLIBCXX_HAVE_FLOAT128_MATH

#include <iostream>
#include <iomanip> // For setw().
#include <sstream>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	       __float128 __x)
    {
      auto __sci = __os.flags() & std::ios::scientific;
      auto __hex = __os.flags() & std::ios::fixed
		&& __os.flags() & std::ios::scientific;
      //auto __hex = __os.flags() & (std::ios::fixed | std::ios::scientific);
      auto __upper = __os.flags() & std::ios::uppercase;
      auto __width = __os.width();
      std::ostringstream __fmt;
      __fmt << '%';

      if (__os.flags() & std::ios::showpos)
	__fmt << '+';
      else
	__fmt << ' '; // Space instead of plus standard?

      if (__os.flags() & std::ios::left)
	__fmt << '-';

      __fmt << __os.width() << '.' << __os.precision() << 'Q';

      if (__hex)
	__fmt << (__upper ? 'A' : 'a');
      else if (__sci)
	__fmt << (__upper ? 'E' : 'e');
      else
	__fmt << (__upper ? 'G' : 'g');

      constexpr int __strlen = 1000;
      char __str[__strlen];
      quadmath_snprintf(__str, __strlen, __fmt.str().c_str(), __x) ;
      __os << __str;
      return __os;
    }

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_istream<_CharT, _Traits>&
    operator>>(std::basic_istream<_CharT, _Traits>& __is, __float128& __x)
    {
      constexpr int __strlen = 160;
      char __str[__strlen];
      __is >> std::setw(__strlen) >> __str;
      __x = strtoflt128(__str, 0);
      return __is;
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_HAVE_FLOAT128_MATH

#endif // _GLIBCXX_BITS_FLOAT128_IO_TCC
