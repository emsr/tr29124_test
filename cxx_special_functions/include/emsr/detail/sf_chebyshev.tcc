
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_chebyshev.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_CHEBYSHEV_TCC
#define SF_CHEBYSHEV_TCC 1

#include <tuple>

#include <emsr/specfun_state.h>

namespace emsr
{
namespace detail
{

  /**
   * Return a Chebyshev polynomial of non-negative order @f$ n @f$
   * and real argument @f$ x @f$ by the recursion
   * @f[
   *    C_n(x) = 2xC_{n-1} - C_{n-2}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   * @param C0 The value of the zeroth-order Chebyshev polynomial at @f$ x @f$
   * @param C1 The value of the first-order Chebyshev polynomial at @f$ x @f$
   */
  template<typename Tp>
    std::tuple<Tp, Tp, Tp>
    chebyshev_recur(unsigned int n, Tp x, Tp C0, Tp C1)
    {
      auto Ck = Tp{2} * x * C1 - C0;
      for (unsigned int j = 2; j < n; ++j)
      {
	C0 = C1;
	C1 = Ck;
	Ck = Tp{2} * x * C1 - C0;
      }
      return std::make_tuple(Ck, C1, C0);
    }

  /**
   * Return the Chebyshev polynomial of the first kind @f$ T_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the first kind is defined by:
   * @f[
   *    T_n(x) = \cos(n \theta)
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    emsr::chebyshev_t_t<Tp>
    chebyshev_t(unsigned int n, Tp x)
    {
      auto T0 = Tp{1};
      if (n == 0)
	return {n, x, T0, Tp{0}, Tp{0}};

      auto T1 = x;
      if (n == 1)
	return {n, x, T1, T0, Tp{0}};

      auto Ts = chebyshev_recur(n, x, T0, T1);
      return {n, x, std::get<0>(Ts), std::get<1>(Ts), std::get<2>(Ts)};
    }

  /**
   * Return the Chebyshev polynomial of the second kind @f$ U_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the second kind is defined by:
   * @f[
   *    U_n(x) = \frac{\sin \left[(n+1)\theta \right]}{\sin(\theta)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    emsr::chebyshev_u_t<Tp>
    chebyshev_u(unsigned int n, Tp x)
    {
      auto U0 = Tp{1};
      if (n == 0)
	return {n, x, U0, Tp{0}, Tp{0}};

      auto U1 = Tp{2} * x;
      if (n == 1)
	return {n, x, U1, U0, Tp{0}};

      auto Us = chebyshev_recur(n, x, U0, U1);
      return {n, x, std::get<0>(Us), std::get<1>(Us), std::get<2>(Us)};
    }

  /**
   * Return the Chebyshev polynomial of the third kind @f$ V_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the third kind is defined by:
   * @f[
   *    V_n(x) = \frac{\cos \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\cos \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    emsr::chebyshev_v_t<Tp>
    chebyshev_v(unsigned int n, Tp x)
    {
      auto V0 = Tp{1};
      if (n == 0)
	return {n, x, V0, Tp{0}, Tp{0}};

      auto V1 = Tp{2} * x - Tp{1};
      if (n == 1)
	return {n, x, V1, V0, Tp{0}};

      auto Vs = chebyshev_recur(n, x, V0, V1);
      return {n, x, std::get<0>(Vs), std::get<1>(Vs), std::get<2>(Vs)};
    }

  /**
   * Return the Chebyshev polynomial of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the fourth kind is defined by:
   * @f[
   *    W_n(x) = \frac{\sin \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\sin \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    emsr::chebyshev_w_t<Tp>
    chebyshev_w(unsigned int n, Tp x)
    {
      auto W0 = Tp{1};
      if (n == 0)
	return {n, x, W0, Tp{0}, Tp{0}};

      auto W1 = Tp{2} * x + Tp{1};
      if (n == 1)
	return {n, x, W1, W0, Tp{0}};

      auto Ws = chebyshev_recur(n, x, W0, W1);
      return {n, x, std::get<0>(Ws), std::get<1>(Ws), std::get<2>(Ws)};
    }

} // namespace detail
} // namespace emsr

#endif // SF_CHEBYSHEV_TCC
