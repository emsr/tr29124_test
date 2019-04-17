// Special functions -*- C++ -*-

// Copyright (C) 2017-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_stirling.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// C++ Mathematical Special Functions
//

#ifndef _GLIBCXX_BITS_SF_STIRLING_TCC
#define _GLIBCXX_BITS_SF_STIRLING_TCC 1

#pragma GCC system_header

#include <vector>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

// Implementation-space details.
namespace __detail
{
  /**
   * Return the Stirling number of the second kind from lookup
   * or by series expansion.
   *
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \genfrac{\{}{\}}{0pt}{0}{n}{m}
   *      = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   * The Stirling number of the second kind is denoted by other symbols
   * in the literature: 
   * @f$ \sigma_n^{(m)} @f$, @f$ \textit{S}_n^{(m)} @f$ and others.
   *
   * @todo Find a way to predict the maximum Stirling number supported
   *       for a given type.
   */
  template<typename _Tp>
    _Tp
    __stirling_2_series(unsigned int __n, unsigned int __m)
    {
      if (__m > _S_num_factorials<_Tp>)
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      auto __lf1 = __log_factorial<_Tp>(__k);
	      auto __lf2 = __log_factorial<_Tp>(__m - __k);
	      _S2 += ((__m - __k) & 1 ? _Tp{-1} : _Tp{1})
		   * std::exp(__n * std::log(__k) - __lf1 - __lf2);
	    }
	  return _S2;
	}
      else
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      _S2 += ((__m - __k) & 1 ? _Tp{-1} : _Tp{1})
		   * std::pow(__k, __n)
		   / __factorial<_Tp>(__k)
		   / __factorial<_Tp>(__m - __k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(_S2);
	}
    }

  /**
   * Return the Stirling number of the second kind by recursion.
   * The recursion is
   * @f[
   *   \genfrac{\{}{\}}{0pt}{0}{n}{m} = m \genfrac{\{}{\}}{0pt}{0}{n-1}{m}
   *                                    + \genfrac{\{}{\}}{0pt}{0}{n-1}{m-1}
   * @f]
   * with starting values
   * @f[
   *   \genfrac{\{}{\}}{0pt}{0}{0}{0\rightarrow m} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   \genfrac{\{}{\}}{0pt}{0}{0\rightarrow n}{0} = {1, 0, 0, ..., 0}
   * @f]
   *
   * The Stirling number of the second kind is denoted by other symbols
   * in the literature: 
   * @f$ \sigma_n^{(m)} @f$, @f$ \textit{S}_n^{(m)} @f$ and others.
   */
  template<typename _Tp>
    _Tp
    __stirling_2_recur(unsigned int __n, unsigned int __m)
    {
      if (__n == 0)
	return _Tp(__m == 0);
      else if (__m == 0)
	return _Tp(__n == 0);
      else
	{
	  std::vector<_Tp> __sigold(__m + 1), __signew(__m + 1);
	  __sigold[1] = _Tp{1};
	  if (__n == 1)
	    return __sigold[__m];
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      __signew[1] = __sigold[1];
	      for (auto __im = 2u; __im <= __m; ++__im)
		__signew[__im] = __im * __sigold[__im] + __sigold[__im - 1];
	      std::swap(__sigold, __signew);
	    }
	  return __signew[__m];
	}
    }

  /**
   * Return the Stirling number of the second kind from lookup
   * or by series expansion.
   *
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
   */
  template<typename _Tp>
    _Tp
    __stirling_2(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return _Tp{0};
      else if (__m == __n)
	return _Tp{1};
      else if (__m == 0 && __n >= 1)
	return _Tp{0};
      else
	return __stirling_2_recur<_Tp>(__n, __m);
    }

  /**
   * Return a vector of Stirling numbers of the second kind by recursion.
   * The recursion is
   * @f[
   *   \genfrac{\{}{\}}{0pt}{0}{n}{m} = m \genfrac{\{}{\}}{0pt}{0}{n-1}{m}
   *                                    + \genfrac{\{}{\}}{0pt}{0}{n-1}{m-1}
   * @f]
   * with starting values
   * @f[
   *   \genfrac{\{}{\}}{0pt}{0}{0}{0\rightarrow m} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   \genfrac{\{}{\}}{0pt}{0}{0\rightarrow n}{0} = {1, 0, 0, ..., 0}
   * @f]
   *
   * The Stirling number of the second kind is denoted by other symbols
   * in the literature: 
   * @f$ \sigma_n^{(m)} @f$, @f$ \textit{S}_n^{(m)} @f$ and others.
   */
  template<typename _Tp>
    std::vector<_Tp>
    __stirling_2_recur(unsigned int __n)
    {
      if (__n == 0)
	return std::vector<_Tp>(1, _Tp{1});
      else
	{
	  std::vector<_Tp> __sigold(__n + 1), __signew(__n + 1);
	  __sigold[0] = __signew[0] = _Tp{0};
	  __sigold[1] = _Tp{1};
	  if (__n == 1)
	    return __sigold;
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      __signew[1] = __sigold[1];
	      for (auto __im = 2u; __im <= __n; ++__im)
		__signew[__im] = __im * __sigold[__im] + __sigold[__im - 1];
	      std::swap(__sigold, __signew);
	    }
	  return __signew;
	}
    }

  /**
   * Return a vector of Stirling numbers of the second kind.
   * or by series expansion.
   *
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
   */
  template<typename _Tp>
    std::vector<_Tp>
    __stirling_2(unsigned int __n)
    { return __stirling_2_recur<_Tp>(__n); }

  /**
   * Return the Stirling number of the second kind.
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
   */
  template<typename _Tp>
    _Tp
    __log_stirling_2(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return -std::numeric_limits<_Tp>::infinity();
      else if (__m == __n)
	return _Tp{0};
      else if (__m == 0 && __n >= 1)
	return -std::numeric_limits<_Tp>::infinity();
      else
	return std::log(__stirling_2<_Tp>(__n, __m));
    }

  /**
   * Return the Stirling number of the first kind by recursion.
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __stirling_1_recur(unsigned int __n, unsigned int __m)
    {
      if (__n == 0)
	return _Tp(__m == 0);
      else if (__m == 0)
	return _Tp(__n == 0);
      else
	{
	  std::vector<_Tp> _Sold(__m + 1), _Snew(__m + 1);
	  _Sold[1] = _Tp{1};
	  if (__n == 1)
	    return _Sold[__m];
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      for (auto __im = 1u; __im <= __m; ++__im)
		_Snew[__im] = _Sold[__im - 1] - __in * _Sold[__im];
	      std::swap(_Sold, _Snew);
	    }
	  return _Snew[__m];
	}
    }

  /**
   * Return the Stirling number of the first kind.
   *
   * The Stirling numbers of the first kind are the coefficients of
   * the Pocchammer polynomials:
   * @f[
   *   (x)_n = \sum_{k=0}^{n} S_n^{(k)} x^k
   * @f]
   *
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
   */
  template<typename _Tp>
    _Tp
    __stirling_1(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return _Tp{0};
      else if (__m == __n)
	return _Tp{1};
      else if (__m == 0 && __n >= 1)
	return _Tp{0};
      else
        return __stirling_1_recur<_Tp>(__n, __m);
    }

  /**
   * Return a vector of Stirling numbers of the first kind by recursion.
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   */
  template<typename _Tp>
    std::vector<_Tp>
    __stirling_1_recur(unsigned int __n)
    {
      if (__n == 0)
	return std::vector<_Tp>(1, _Tp{1});
      else
	{
	  std::vector<_Tp> _Sold(__n + 1), _Snew(__n + 1);
	  _Sold[0] = _Snew[0] = _Tp{0};
	  _Sold[1] = _Tp{1};
	  if (__n == 1)
	    return _Sold;
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      for (auto __im = 1u; __im <= __n; ++__im)
		_Snew[__im] = _Sold[__im - 1] - __in * _Sold[__im];
	      std::swap(_Sold, _Snew);
	    }
	  return _Snew;
	}
    }

  /**
   * Return a vector of Stirling numbers of the first kind.
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   */
  template<typename _Tp>
    std::vector<_Tp>
    __stirling_1(unsigned int __n)
    { return __stirling_1_recur<_Tp>(__n); }

  /**
   * Return the logarithm of the absolute value of Stirling number
   * of the first kind.
   */
  template<typename _Tp>
    _Tp
    __log_stirling_1(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return -std::numeric_limits<_Tp>::infinity();
      else if (__m == __n)
	return _Tp{0};
      else if (__m == 0 && __n >= 1)
	return -std::numeric_limits<_Tp>::infinity();
      else
	return std::log(std::abs(__stirling_1<_Tp>(__n, __m)));
    }

  /**
   * Return the sign of the exponent of the logarithm of the Stirling number
   * of the first kind.
   */
  template<typename _Tp>
    inline _Tp
    __log_stirling_1_sign(unsigned int __n, unsigned int __m)
    { return (__n + __m) & 1 ? _Tp{-1} : _Tp{+1}; }

  /**
   * Return the Lah number by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename _Tp>
    _Tp
    __lah_recur(unsigned int __n, unsigned int __k)
    {
      if (__k > __n)
	return _Tp{0};
      else if (__n == 0)
	return (__k == 0 ? _Tp{1} : _Tp{0});
      else
	{
	  _Tp _Lnn = 1;
	  for (unsigned int __i = 1u; __i <= __n - __k; ++__i)
	    _Lnn *= _Tp(__n - __i + 1) * _Tp(__n - __i) / _Tp(__i);
	  return _Lnn;
	}
    }

  /**
   * Return the Lah number.
   * Lah numbers are defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename _Tp>
    inline _Tp
    __lah(unsigned int __n, unsigned int __k)
    { return __lah_recur<_Tp>(__n, __k); }

  /**
   * Return a vector of Lah numbers defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename _Tp>
    std::vector<_Tp>
    __lah_recur(unsigned int __n)
    {
      if (__n == 0)
	return std::vector<_Tp>(1, _Tp{1});
      else
	{
	  std::vector<_Tp> _L(__n + 1);
	  _Tp _Lnn = 1;
	  _L[__n] = _Lnn;
	  for (unsigned int __i = 1u; __i <= __n; ++__i)
	    {
	      _Lnn *= _Tp(__n - __i + 1) * _Tp(__n - __i) / _Tp(__i);
	      _L[__n - __i] = _Lnn;
	    }
	  return _L;
	}
    }

  /**
   * Return a vector of Lah numbers.
   * Lah numbers are defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename _Tp>
    inline std::vector<_Tp>
    __lah(unsigned int __n)
    { return __lah_recur<_Tp>(__n); }

  /**
   * Return a vector of the Bell numbers by summation.
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)} = \sum_{k=1}^{n}{n-1 \choose k-1}B(n-k)
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename _Tp>
    std::vector<_Tp>
    __bell_series(unsigned int __n)
    {
      std::vector<_Tp> __bell(__n + 1);
      __bell[0] = _Tp{1};

      /// @todo Test for blowup in Bell number summation.
      for (unsigned int __i = 1; __i <= __n; ++__i)
	for (unsigned int __j = 1; __j <= __i; ++__j)
	  __bell[__i] += __bell[__i - __j]
		       * std::__detail::__binomial<_Tp>(__i - 1, __j - 1);

      return __bell;
    }

  /**
   * Return a vector of the Bell numbers.
   */
  template<typename _Tp>
    inline std::vector<_Tp>
    __bell(unsigned int __n)
    { return __bell_series<_Tp>(__n); }

  /**
   * Evaluate the Bell polynomial
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)}x^k
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename _Tp, typename _Up>
    inline _Up
    __bell(unsigned int __n, _Up __x)
    {
      const auto _Sn = __stirling_2<_Tp>(__n);
      auto __bell = _Sn[__n];
      for (unsigned int __i = 1; __i < __n; ++__i)
	__bell = _Sn[__n - __i] + __x * __bell;
      return __bell;
    }

} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_SF_STIRLING_TCC

