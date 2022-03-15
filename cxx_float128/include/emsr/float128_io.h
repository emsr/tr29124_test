
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

/** @file emsr/float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef FLOAT128_IO_H
#define FLOAT128_IO_H 1

#include <emsr/float128.h>

#ifdef EMSR_HAVE_FLOAT128
#if __has_include(<quadmath.h>)

#include <iosfwd>
#include <quadmath.h>

namespace std
{

  template<typename CharT, typename Traits = std::char_traits<CharT>>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
	       __float128 x);

  template<typename CharT, typename Traits = std::char_traits<CharT>>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& is, __float128& x);

} // namespace std

#include <emsr/float128_io.tcc>

#endif // __has_include(<quadmath.h>)
#endif // EMSR_HAVE_FLOAT128

#endif // FLOAT128_IO_H
