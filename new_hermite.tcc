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

/** @file tr1/sf_hermite.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

//
// ISO C++ TR29124: Special functions
//

#ifndef _GLIBCXX_BITS_NEW_HERMITE_TCC
#define _GLIBCXX_BITS_NEW_HERMITE_TCC 1



  /**
   *  @brief  Compute the Hermite polynomial by recursion.
   */
  template<typename _Tp>
    _Tp
    __hermite_recur(unsigned int __n, _Tp __x)
    {
      auto __hnm1 = _Tp{0};

      auto __hn = _Tp{1};
      if (__n == 0)
        return __hn;

      _Tp __hnp1;
      for (int __i = 0; __i < __n; ++__i)
        {
          __hnp1 = _Tp{2} *__x * __hn
                 - _Tp(2 * __i) * __hnm1;
          __hnm1 = __hn;
          __hn = __hnp1;
        }

      return __hnp1;
    }


  /**
   *  @brief  Compute the normalized Hermite polynomial by recursion.
   *  @todo  Tabulate sqrt(int) or even sqrt(i)/sqrt(i+1) and add a helper.
   */
  template<typename _Tp>
    _Tp
    __hermite_norm_recur(unsigned int __n, _Tp __x)
    {
      const auto __INV_ROOT4_PI = _Tp{0.7511255444649425L};

      auto __htnm1 = _Tp{0};
      auto __htn = __INV_ROOT4_PI;
      if (__n == 0)
        return __htn;

      _Tp __htnp1;
      for (int __i = 0; __i < __n; ++__i)
        {
          __htnp1 = __x * std::sqrt(_Tp{2} / _Tp(__i + 1)) * __htn
                  - std::sqrt(_Tp(__i) / _Tp(__i + 1)) * __htnm1;
          __htnm1 = __htn;
          __htn = __htnp1;
        }

      return __htnp1;
    }

#endif // _GLIBCXX_BITS_NEW_HERMITE_TCC
