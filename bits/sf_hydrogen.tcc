// Special functions -*- C++ -*-

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

/** @file bits/sf_hydrogen.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_SF_HYDROGEN_TCC
#define _GLIBCXX_BITS_SF_HYDROGEN_TCC 1

#pragma GCC system_header

#include <complex>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the bound-state Coulomb wave-function.
   */
  template <typename _Tp>
    std::complex<_Tp>
    __hydrogen(unsigned int __n,
               unsigned int __l, unsigned int __m,
               _Tp __Z, _Tp __r, _Tp __theta, _Tp __phi)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();

      if (__isnan(__Z) || __isnan(__r) || __isnan(__theta) || __isnan(__phi))
	return std::complex<_Tp>{_S_NaN, _S_NaN};
      else if(__n < 1)
	std::__throw_domain_error(__N("__hydrogen: "
				      "level number less than one"));
      else if(__l > __n - 1)
	std::__throw_domain_error(__N("__hydrogen: "
				      "angular momentum number too large"));
      else if(__Z <= _Tp(0))
	std::__throw_domain_error(__N("__hydrogen: non-positive charge"));
      else if(__r < _Tp(0))
	std::__throw_domain_error(__N("__hydrogen: negative radius"));
      else
	{
	  const auto __A = _Tp(2) * __Z / __n;

	  const auto __pre = std::sqrt(__A * __A * __A / (_Tp(2) * __n));
	  const auto __ln_a = __log_gamma(__n + __l + 1);
	  const auto __ln_b = __log_gamma(__n - __l);
	  const auto __ex = std::exp((__ln_b - __ln_a) / _Tp(2));
	  const auto __norm = __pre * __ex;

	  const auto __rho = __A * __r;
	  const auto __ea = std::exp(-__rho / _Tp(2));
	  const auto __pp = std::pow(__rho, __l);
	  const auto __lag = __assoc_laguerre(__n - __l - 1, 2 * __l + 1,
                                        	__rho);
	  const auto __sphh = __sph_legendre(__l, __m, __theta)
 			    * std::polar(_Tp(1), _Tp(__m) * __phi);

	  const auto __psi = __norm * __ea * __pp * __lag * __sphh;

	  return __psi;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_BITS_SF_HYDROGEN_TCC
