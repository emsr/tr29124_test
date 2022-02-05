// -*- C++ -*- header.

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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

/** @file emsr/float128.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef FLOAT128_IO_TCC
#define FLOAT128_IO_TCC 1

#ifdef _GLIBCXX_USE_FLOAT128

#include <iostream>
#include <iomanip> // For setw().
#include <sstream>

namespace std
{

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& os,
	       __float128 x)
    {
      auto sci = os.flags() & std::ios::scientific;
      auto hex = os.flags() & std::ios::fixed
		&& os.flags() & std::ios::scientific;
      //auto hex = os.flags() & (std::ios::fixed | std::ios::scientific);
      auto upper = os.flags() & std::ios::uppercase;
      //auto width = os.width();
      std::ostringstream fmt;
      fmt << '%';

      if (os.flags() & std::ios::showpos)
	fmt << '+';
      else
	fmt << ' '; // Space instead of plus standard?

      if (os.flags() & std::ios::left)
	fmt << '-';

      fmt << os.width() << '.' << os.precision() << 'Q';

      if (hex)
	fmt << (upper ? 'A' : 'a');
      else if (sci)
	fmt << (upper ? 'E' : 'e');
      else
	fmt << (upper ? 'G' : 'g');

      constexpr int strlen = 1000;
      char str[strlen];
      quadmath_snprintf(str, strlen, fmt.str().c_str(), x) ;
      os << str;
      return os;
    }

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_istream<_CharT, _Traits>&
    operator>>(std::basic_istream<_CharT, _Traits>& is, __float128& x)
    {
      constexpr int strlen = 160;
      char str[strlen];
      is >> std::setw(strlen) >> str;
      x = strtoflt128(str, 0);
      return is;
    }

} // namespace std

#endif // _GLIBCXX_USE_FLOAT128

#endif // EXT_FLOAT128_IO_TCC
