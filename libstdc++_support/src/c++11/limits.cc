// Static data members of -*- C++ -*- numeric_limits classes

// Copyright (C) 1999-2016 Free Software Foundation, Inc.
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

// Written by Gabriel Dos Reis <Gabriel.Dos-Reis@cmla.ens-cachan.fr>

//
// ISO C++ 14882:1998
// 18.2.1
//

#include <limits>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

#define const _GLIBCXX_USE_CONSTEXPR

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  // __float128
  const bool numeric_limits<__float128>::is_specialized;
  const int  numeric_limits<__float128>::digits;
  const int  numeric_limits<__float128>::digits10;
  const int  numeric_limits<__float128>::max_digits10;
  const bool numeric_limits<__float128>::is_signed;
  const bool numeric_limits<__float128>::is_integer;
  const bool numeric_limits<__float128>::is_exact;
  const int  numeric_limits<__float128>::radix;
  const int  numeric_limits<__float128>::min_exponent;
  const int  numeric_limits<__float128>::min_exponent10;
  const int  numeric_limits<__float128>::max_exponent;
  const int  numeric_limits<__float128>::max_exponent10;
  const bool numeric_limits<__float128>::has_infinity;
  const bool numeric_limits<__float128>::has_quiet_NaN;
  const bool numeric_limits<__float128>::has_signaling_NaN;
  const float_denorm_style numeric_limits<__float128>::has_denorm;
  const bool numeric_limits<__float128>::has_denorm_loss;
  const bool numeric_limits<__float128>::is_iec559;
  const bool numeric_limits<__float128>::is_bounded;
  const bool numeric_limits<__float128>::is_modulo;
  const bool numeric_limits<__float128>::traps;
  const bool numeric_limits<__float128>::tinyness_before;
  const float_round_style numeric_limits<__float128>::round_style;
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#undef const

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace
