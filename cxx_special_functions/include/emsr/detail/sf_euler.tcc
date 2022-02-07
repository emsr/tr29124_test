
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

/** @file bits/sf_euler.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// C++ Mathematical Special Functions
//

#ifndef SF_EULER_TCC
#define SF_EULER_TCC 1

#include <vector>

#include <emsr/sf_gamma.h> // binomial

namespace emsr
{
namespace detail
{

  /**
   * Return the Euler number from lookup or by series expansion.
   *
   * The Euler numbers are given by the recursive sum:
   * @f[
   *   E_n = B_n(1) = B_n
   * @f]
   * where @f$ E_0 = 1 @f$, @f$ E_1 = 0 @f$, @f$ E_2 = -1 @f$
   *
   * @todo Find a way to predict the maximum Euler number for a type.
   */
  template<typename Tp>
    Tp
    euler_series(unsigned int n)
    {
      static constexpr std::size_t s_len = 22;
      static constexpr Tp
      s_num[s_len]
      {
	 1ll, 0,
	-1ll, 0ll,
	 5ll, 0ll,
	-61ll, 0ll,
	 1385ll, 0ll,
	-50521ll, 0ll,
	 2702765ll, 0ll,
	-199360981ll, 0ll,
	 19391512145ll, 0ll,
	-2404879675441ll, 0ll,
	 370371188237525ll, 0ll,
	//-69348874393137901ll, 0ll,
      };

      if (n == 0)
	return Tp{1};
      else if (n & 1)
	return Tp{0};
      else if (n == 2)
        return Tp{-1};
      else if (n < s_len)
	return s_num[n];
      else
	{
	  std::vector<Tp> En(n + 1);
	  En[0] = Tp{1};
	  En[1] = Tp{0};
	  En[2] = Tp{-1};

	  for (auto i = 3u; i <= n; ++i)
	    {
	      En[i] = 0;

	      if (i % 2 == 0)
		{
		  for (auto j = 2u; j <= i; j += 2u)
		    En[i] -= binomial<Tp>(i, j)
			      * En[i - j];
		}
	    }
	  return En[n];
	}
    }

  /**
   * @brief This returns Euler number @f$ E_n @f$.
   *
   * @param n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename Tp>
    inline Tp
    euler(unsigned int n)
    { return euler_series<Tp>(n); }

  /**
   * Return the Euler polynomial @f$ E_n(x) @f$ of order n at argument x.
   *
   * The derivative is proportional to the previous polynomial:
   * @f[
   *   E_n'(x) = n E_{n-1}(x)
   * @f]
   *
   * @f[
   *   E_n(1/2) = \frac{E_n}{2^n}, \mbox{ where } E_n
   *             \mbox{ is the n-th Euler number.}
   * @f]
   */
  template<typename Tp>
    Tp
    euler(unsigned int n, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  auto bx1 = bernoulli(n + 1, x );
	  auto bx2 = bernoulli(n + 1, Tp{0.5L} * x );

	  auto E_n = Tp{2} * (bx1 - bx2 * std::pow(Tp{2}, Tp(n + 1)))
		    / Tp(n + 1);

	  return E_n;
	}
    }

  /**
   * Return the Eulerian number of the first kind.
   * The Eulerian numbers of the first kind are defined by recursion:
   * @f[
   *   \genfrac{\langle}{\rangle}{0pt}{0}{n}{m}
   *     = (n-m)\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}
   *     + (m+1)\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}
   *   \mbox{ for } n > 0
   * @f]
   * Note that @f$ A(n,m) @f$ is a common older notation.
   */
  template<typename Tp>
    Tp
    eulerian_1_recur(unsigned int n, unsigned int m)
    {
      if (m == 0)
	return Tp{1};
      else if (m >= n)
	return Tp{0};
      else if (m == n - 1)
	return Tp{1};
      else if (n - m - 1 < m) // Symmetry.
	return eulerian_1_recur<Tp>(n, n - m - 1);
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(m + 1), Anew(m + 1);
	  Aold[0] = Tp{1};
	  Anew[0] = Tp{1};
	  Anew[1] = Tp{1};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= m; ++im)
		Anew[im] = (in - im) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew[m];
	}
    }

  /**
   * Return the Eulerian number of the first kind.
   * The Eulerian numbers of the first kind are defined by recursion:
   * @f[
   *   \genfrac{\langle}{\rangle}{0pt}{0}{n}{m}
   *     = (n-m)\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}
   *     + (m+1)\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}
   *   \mbox{ for } n > 0
   * @f]
   * Note that @f$ A(n,m) @f$ is a common older notation.
   */
  template<typename Tp>
    inline Tp
    eulerian_1(unsigned int n, unsigned int m)
    { return eulerian_1_recur<Tp>(n, m); }

  /**
   * Return a vector Eulerian numbers of the first kind by recursion.
   * The recursion is
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    std::vector<Tp>
    eulerian_1_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      //else if (m == n - 1)
	//return Tp{1};
      //else if (n - m - 1 < m) // Symmetry.
	//return eulerian_1_recur<Tp>(n, n - m - 1);
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(n + 1), Anew(n + 1);
	  Aold[0] = Anew[0] = Tp{1};
	  Anew[1] = Tp{1};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= n; ++im)
		Anew[im] = (in - im) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew;
	}
    }

  /**
   * Return a vector Eulerian numbers of the first kind.
   * The Eulerian numbers are defined by recursion:
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline std::vector<Tp>
    eulerian_1(unsigned int n)
    { return eulerian_1_recur<Tp>(n); }

  /**
   * Return the Eulerian number of the second kind by recursion:
   * @f[
   *   \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n}{m}\right\rangle
   *   = (2n-m-1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}\right\rangle
   *   + (m+1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}\right\rangle
   *       \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    Tp
    eulerian_2_recur(unsigned int n, unsigned int m)
    {
      if (m == 0)
	return Tp{1};
      else if (m >= n)
	return Tp{0};
      else if (n == 0)
	return Tp{1};
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(m + 1), Anew(m + 1);
	  Aold[0] = Tp{1};
	  Anew[0] = Tp{1};
	  Anew[1] = Tp{2};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= m; ++im)
		Anew[im] = (2 * in - im - 1) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew[m];
	}
    }

  /**
   * Return the Eulerian number of the second kind.
   * The Eulerian numbers of the second kind are defined by recursion:
   * @f[
   *   \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n}{m}\right\rangle
   *   = (2n-m-1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}\right\rangle
   *   + (m+1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}\right\rangle
   *       \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline Tp
    eulerian_2(unsigned int n, unsigned int m)
    { return eulerian_2_recur<Tp>(n, m); }

  /**
   * Return a vector of Eulerian numbers of the second kind.
   * @f[
   *   \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n}{m}\right\rangle
   *   = (2n-m-1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}\right\rangle
   *   + (m+1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}\right\rangle
   *       \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    std::vector<Tp>
    eulerian_2_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      //else if (m >= n)
	//return Tp{0};
      //else if (n == 0)
	//return Tp{1};
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(n + 1), Anew(n + 1);
	  Aold[0] = Anew[0] = Tp{1};
	  Anew[1] = Tp{2};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= n; ++im)
		Anew[im] = (2 * in - im - 1) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew;
	}
    }

  /**
   * Return a vector of Eulerian numbers of the second kind.
   * @f[
   *   \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n}{m}\right\rangle
   *   = (2n-m-1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}\right\rangle
   *   + (m+1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}\right\rangle
   *       \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline std::vector<Tp>
    eulerian_2(unsigned int n)
    { return eulerian_2_recur<Tp>(n); }

} // namespace detail
} // namespace emsr

#endif // SF_EULER_TCC
