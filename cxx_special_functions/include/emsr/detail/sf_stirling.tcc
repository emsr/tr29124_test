
// Copyright (C) 2017-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

#ifndef SF_STIRLING_TCC
#define SF_STIRLING_TCC 1

#include <vector>

#include <emsr/sf_gamma.h> // factorials.

namespace emsr
{
namespace detail
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
  template<typename Tp>
    Tp
    stirling_2_series(unsigned int n, unsigned int m)
    {
      if (m > s_num_factorials<Tp>)
	{
	  auto S2 = Tp{0};
	  for (auto k = 0u; k <= m; ++k)
	    {
	      auto lf1 = log_factorial<Tp>(k);
	      auto lf2 = log_factorial<Tp>(m - k);
	      S2 += ((m - k) & 1 ? Tp{-1} : Tp{1})
		   * std::exp(n * std::log(k) - lf1 - lf2);
	    }
	  return S2;
	}
      else
	{
	  auto S2 = Tp{0};
	  for (auto k = 0u; k <= m; ++k)
	    {
	      S2 += ((m - k) & 1 ? Tp{-1} : Tp{1})
		   * std::pow(k, n)
		   / factorial<Tp>(k)
		   / factorial<Tp>(m - k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(S2);
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
  template<typename Tp>
    Tp
    stirling_2_recur(unsigned int n, unsigned int m)
    {
      if (n == 0)
	return Tp(m == 0);
      else if (m == 0)
	return Tp(n == 0);
      else
	{
	  std::vector<Tp> sigold(m + 1), signew(m + 1);
	  sigold[1] = Tp{1};
	  if (n == 1)
	    return sigold[m];
	  for (auto in = 1u; in <= n; ++in)
	    {
	      signew[1] = sigold[1];
	      for (auto im = 2u; im <= m; ++im)
		signew[im] = im * sigold[im] + sigold[im - 1];
	      std::swap(sigold, signew);
	    }
	  return signew[m];
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
  template<typename Tp>
    Tp
    stirling_2(unsigned int n, unsigned int m)
    {
      if (m > n)
	return Tp{0};
      else if (m == n)
	return Tp{1};
      else if (m == 0 && n >= 1)
	return Tp{0};
      else
	return stirling_2_recur<Tp>(n, m);
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
  template<typename Tp>
    std::vector<Tp>
    stirling_2_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      else
	{
	  std::vector<Tp> sigold(n + 1), signew(n + 1);
	  sigold[0] = signew[0] = Tp{0};
	  sigold[1] = Tp{1};
	  if (n == 1)
	    return sigold;
	  for (auto in = 1u; in <= n; ++in)
	    {
	      signew[1] = sigold[1];
	      for (auto im = 2u; im <= n; ++im)
		signew[im] = im * sigold[im] + sigold[im - 1];
	      std::swap(sigold, signew);
	    }
	  return signew;
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
  template<typename Tp>
    std::vector<Tp>
    stirling_2(unsigned int n)
    { return stirling_2_recur<Tp>(n); }

  /**
   * Return the Stirling number of the second kind.
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
   */
  template<typename Tp>
    Tp
    log_stirling_2(unsigned int n, unsigned int m)
    {
      if (m > n)
	return -std::numeric_limits<Tp>::infinity();
      else if (m == n)
	return Tp{0};
      else if (m == 0 && n >= 1)
	return -std::numeric_limits<Tp>::infinity();
      else
	return std::log(stirling_2<Tp>(n, m));
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
  template<typename Tp>
    Tp
    stirling_1_recur(unsigned int n, unsigned int m)
    {
      if (n == 0)
	return Tp(m == 0);
      else if (m == 0)
	return Tp(n == 0);
      else
	{
	  std::vector<Tp> Sold(m + 1), Snew(m + 1);
	  Sold[1] = Tp{1};
	  if (n == 1)
	    return Sold[m];
	  for (auto in = 1u; in <= n; ++in)
	    {
	      for (auto im = 1u; im <= m; ++im)
		Snew[im] = Sold[im - 1] - in * Sold[im];
	      std::swap(Sold, Snew);
	    }
	  return Snew[m];
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
  template<typename Tp>
    Tp
    stirling_1(unsigned int n, unsigned int m)
    {
      if (m > n)
	return Tp{0};
      else if (m == n)
	return Tp{1};
      else if (m == 0 && n >= 1)
	return Tp{0};
      else
        return stirling_1_recur<Tp>(n, m);
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
  template<typename Tp>
    std::vector<Tp>
    stirling_1_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      else
	{
	  std::vector<Tp> Sold(n + 1), Snew(n + 1);
	  Sold[0] = Snew[0] = Tp{0};
	  Sold[1] = Tp{1};
	  if (n == 1)
	    return Sold;
	  for (auto in = 1u; in <= n; ++in)
	    {
	      for (auto im = 1u; im <= n; ++im)
		Snew[im] = Sold[im - 1] - in * Sold[im];
	      std::swap(Sold, Snew);
	    }
	  return Snew;
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
  template<typename Tp>
    std::vector<Tp>
    stirling_1(unsigned int n)
    { return stirling_1_recur<Tp>(n); }

  /**
   * Return the logarithm of the absolute value of Stirling number
   * of the first kind.
   */
  template<typename Tp>
    Tp
    log_stirling_1(unsigned int n, unsigned int m)
    {
      if (m > n)
	return -std::numeric_limits<Tp>::infinity();
      else if (m == n)
	return Tp{0};
      else if (m == 0 && n >= 1)
	return -std::numeric_limits<Tp>::infinity();
      else
	return std::log(std::abs(stirling_1<Tp>(n, m)));
    }

  /**
   * Return the sign of the exponent of the logarithm of the Stirling number
   * of the first kind.
   */
  template<typename Tp>
    inline Tp
    log_stirling_1_sign(unsigned int n, unsigned int m)
    { return (n + m) & 1 ? Tp{-1} : Tp{+1}; }

  /**
   * Return the Lah number by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    Tp
    lah_recur(unsigned int n, unsigned int k)
    {
      if (k > n)
	return Tp{0};
      else if (n == 0)
	return (k == 0 ? Tp{1} : Tp{0});
      else
	{
	  Tp Lnn = 1;
	  for (unsigned int i = 1u; i <= n - k; ++i)
	    Lnn *= Tp(n - i + 1) * Tp(n - i) / Tp(i);
	  return Lnn;
	}
    }

  /**
   * Return the Lah number.
   * Lah numbers are defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    inline Tp
    lah(unsigned int n, unsigned int k)
    { return lah_recur<Tp>(n, k); }

  /**
   * Return a vector of Lah numbers defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    std::vector<Tp>
    lah_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      else
	{
	  std::vector<Tp> L(n + 1);
	  Tp Lnn = 1;
	  L[n] = Lnn;
	  for (unsigned int i = 1u; i <= n; ++i)
	    {
	      Lnn *= Tp(n - i + 1) * Tp(n - i) / Tp(i);
	      L[n - i] = Lnn;
	    }
	  return L;
	}
    }

  /**
   * Return a vector of Lah numbers.
   * Lah numbers are defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    inline std::vector<Tp>
    lah(unsigned int n)
    { return lah_recur<Tp>(n); }

  /**
   * Return a vector of the Bell numbers by summation.
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)} = \sum_{k=1}^{n}{n-1 \choose k-1}B(n-k)
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename Tp>
    std::vector<Tp>
    bell_series(unsigned int n)
    {
      std::vector<Tp> bell(n + 1);
      bell[0] = Tp{1};

      /// @todo Test for blowup in Bell number summation.
      for (unsigned int i = 1; i <= n; ++i)
	for (unsigned int j = 1; j <= i; ++j)
	  bell[i] += bell[i - j] * emsr::detail::binomial<Tp>(i - 1, j - 1);

      return bell;
    }

  /**
   * Return a vector of the Bell numbers.
   */
  template<typename Tp>
    inline std::vector<Tp>
    bell(unsigned int n)
    { return bell_series<Tp>(n); }

  /**
   * Evaluate the Bell polynomial
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)}x^k
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename Tp, typename Up>
    inline Up
    bell(unsigned int n, Up x)
    {
      const auto Sn = stirling_2<Tp>(n);
      auto bell = Sn[n];
      for (unsigned int i = 1; i < n; ++i)
	bell = Sn[n - i] + x * bell;
      return bell;
    }

} // namespace detail
} // namespace emsr

#endif // SF_STIRLING_TCC
