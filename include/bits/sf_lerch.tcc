// Special functions -*- C++ -*-

// Copyright (C) 2006-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_hyperg.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1) 

#ifndef _GLIBCXX_BITS_SF_LERCH_TCC
#define _GLIBCXX_BITS_SF_LERCH_TCC 1

#pragma GCC system_header

#include <ext/summation.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

// Implementation-space details.
namespace __detail
{

  /**
   * A functor for a vanWijnGaarden compressor.
   * vanWijnGaarden requires:
   *   _Tp operator()(int) that returns a term in the original defining series.
   */
  template<typename _Tp>
    class __lerch_term
    {
    public:

      using value_type = _Tp;

      __lerch_term(value_type __z, value_type __s, value_type __a)
      : _M_z{__z}, _M_s{__s}, _M_a{__a}
      { }

      value_type
      operator()(std::size_t __i) const
      {
	return std::pow(_M_z, value_type(__i))
	     / std::pow(_M_a + value_type(__i), _M_s);
      }

    private:

      value_type _M_z;
      value_type _M_s;
      value_type _M_a;
    };

  /**
   * This function blows up on nonpositive integeral parameter a.
   *
   * @param __z The argument.
   * @param __s The order @f$ s != 1 @f$.
   * @param __a The scale parameter @f$ a > -1 @f$.
   */
  template<typename _Tp>
    _Tp
    __lerch_sum(_Tp __z, _Tp __s, _Tp __a)
    {
      const auto _S_nan = __gnu_cxx::__quiet_NaN(__s);
      const auto _S_eps = __gnu_cxx::__epsilon(__s);

      const auto __aint = emsr::fp_is_integer(__a);
      if (__aint && __aint() <= 0)
	return _S_nan;
      else if (std::abs(std::abs(__z) - _Tp{1}) < _S_eps
		&& std::real(__s) <= _Tp{1} + _S_eps)
	return _S_nan;
      else if (std::abs(__z) > _Tp{1} + _S_eps)
	return _S_nan;
      else
	{
	  constexpr auto _S_maxit = 100000u;
	  auto __zpow = _Tp{1};
	  auto __sum = std::pow(__a, -__s);
	  for (auto __k = 1u; __k < _S_maxit; ++__k)
	    {
	      __zpow *= __z;
	      auto __term = __zpow * std::pow(__a + __k, -__s);
	      __sum += __term;
	      if (std::abs(__term / __sum) < _S_eps)
		break;
	    }
	  return __sum;
	}
    }

  /**
   * Try the WenigerDelta<MonotoneVanWijngaarden> composition.
   *
   * @param __z The argument.
   * @param __s The order @f$ s != 1 @f$.
   * @param __a The scale parameter @f$ a > -1 @f$.
   */
  template<typename _Tp>
    _Tp
    __lerch_delta_vanwijngaarden_sum(_Tp __z, _Tp __s, _Tp __a)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__s);
      constexpr auto _S_maxit = 1000u;

      __gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_VanWijngaardenSum<_Tp>> _WDvW;
      if (__z >= _Tp{0})
	{
	  using __lerch_t = __lerch_term<_Tp>;
	  using __lerch_cmp_t = __gnu_cxx::_VanWijngaardenCompressor<__lerch_t>;
	  auto _VwT = __lerch_cmp_t(__lerch_t(__z, __s, __a));
	  for (auto __k = 0u; __k < _S_maxit; ++__k)
	    {
	      auto __term = _VwT[__k];
	      _WDvW += __term;
	      if (std::abs(__term) < _S_eps * std::abs(_WDvW()))
		break;
	    }
	  return _WDvW();
	}
      else
	{
	  auto _LT = __lerch_term<_Tp>(__z, __s, __a);
	  for (auto __k = 0u; __k < _S_maxit; ++__k)
	    {
	      auto __term = _LT(__k);
	      _WDvW += __term;
	      if (std::abs(__term) < _S_eps * std::abs(_WDvW()))
		break;
	    }
	  return _WDvW();
	}
    }

  /**
   * Return the Lerch transcendent @f$ \Phi(z,s,a) @f$.
   *
   * The series is:
   * @f[   *
   *   \Phi(z,s,a) = \sum_{k=0}^{\infty}\frac{z^k}{(a+k^s}
   * @f]
   *
   * This function blows up on nonpositive integeral parameter a.
   *
   * @param __z The argument.
   * @param __s The order @f$ s != 1 @f$.
   * @param __a The scale parameter @f$ a > -1 @f$.
   */
  template<typename _Tp>
    _Tp
    __lerch_phi(_Tp __z, _Tp __s, _Tp __a)
    {
      const auto _S_nan = __gnu_cxx::__quiet_NaN(__s);
      const auto _S_eps = __gnu_cxx::__epsilon(__s);

      if (std::isnan(__z) || std::isnan(__s) || std::isnan(__a))
	return _S_nan;
      else if (std::abs(std::abs(__z) - _Tp{1}) < _S_eps
		&& std::real(__s) <= _Tp{1} + _S_eps)
	return _S_nan;
      else if (std::abs(__z) > _Tp{1} + _S_eps)
	return _S_nan;
      else
	{
	  const auto __aint = emsr::fp_is_integer(__a);

	  const auto __sint = emsr::fp_is_integer(__s);
	  const bool __tinyz = std::abs(__z) < _S_eps; // _S_min?
	  const bool __smallz = !__tinyz && (std::abs(__z) < _Tp{0.5});

	  if (__aint && __aint() <= 0)
	    return _S_nan;
	  else if (__a < _Tp{0})
	    {
	      if (__sint)
		{
		  int __sign = __sint() % 2 == 0 ? +1 : -1;
		  if (__tinyz)
		    return __sign * _Tp{1} / std::pow(std::abs(__a), __s);
		  else
		    {
		      const auto __m = -int(std::floor(__a));
		      const auto __a1 = __a + _Tp(__m);
		      auto __sum1 = _Tp{0};
		      for (int __i = 0; __i < __m; ++__i)
			{
			  __sum1 += __sign * std::pow(std::abs(__z), __i)
				 / std::pow(std::abs(__a + __i), _Tp(__sint()));
			  if (__z < _Tp{0})
			    __sign = -__sign;
			}
		      auto __sum = _Tp{0};
		      if (__smallz)
			__sum = __lerch_sum(__z, __s, __a1);
		      else
			__sum
			  = __lerch_delta_vanwijngaarden_sum(__z, __s, __a1);
		      __sign = 1;
		      if (__z < _Tp{0} && __m % 2 != 0)
			__sign = -1;
		      return __sum1
			   + __sum * __sign * std::pow(std::abs(__z), __m);
		    }
		}
	      else // s is not an integer - Phi is complex.
		return _S_nan;
	    }
	  else if (__tinyz)
	    return _Tp{1} / std::pow(__a, __s);
	  else // a > 0
	    {
	      if (__smallz)
		return __lerch_sum(__z, __s, __a);
	      else
		return __lerch_delta_vanwijngaarden_sum(__z, __s, __a);
	    }
	}
    }

} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_SF_LERCH_TCC
