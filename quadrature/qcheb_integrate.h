/* quadrature/qaws_integration_table.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
 * Copyright (C) 2016-2017 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef QCHEB_INTEGRATE_H
#define QCHEB_INTEGRATE_H 1

#include <array>

namespace __gnu_test
{

  template<typename _FuncTp, typename _Tp>
    void
    qcheb_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		    std::array<_Tp, 13>& __cheb12,
		    std::array<_Tp, 25>& __cheb24)
    {
      _Tp __fval[25], __v[12];

      /* These are the values of cos(pi*k/24) for k=1..11 needed for the
	 Chebyshev expansion of f(x) */

      constexpr _Tp
      __x[11]
      {
	0.9914448613738104,
	0.9659258262890683,
	0.9238795325112868,
	0.8660254037844386,
	0.7933533402912352,
	0.7071067811865475,
	0.6087614290087206,
	0.5000000000000000,
	0.3826834323650898,
	0.2588190451025208,
	0.1305261922200516
      };

      const _Tp __center = (__b + __a) / _Tp{2};
      const _Tp __half_length = (__b - __a) / _Tp{2};

      __fval[0] = __func(__b) / _Tp{2};
      __fval[12] = __func(__center);
      __fval[24] = __func(__a) / _Tp{2};

      for (int __i = 1; __i < 12; ++__i)
	{
	  const std::size_t __j = 24 - __i;
	  const _Tp __u = __half_length * __x[__i - 1];
	  __fval[__i] = __func(__center + __u);
	  __fval[__j] = __func(__center - __u);
	}

      for (int __i = 0; __i < 12; ++__i)
	{
	  const std::size_t __j = 24 - __i;
	  __v[__i] = __fval[__i] - __fval[__j];
	  __fval[__i] = __fval[__i] + __fval[__j];
	}

      {
	const _Tp __alam1 = __v[0] - __v[8];
	const _Tp __alam2 = __x[5] * (__v[2] - __v[6] - __v[10]);
	__cheb12[3] = __alam1 + __alam2;
	__cheb12[9] = __alam1 - __alam2;
      }

      {
	const _Tp __alam1 = __v[1] - __v[7] - __v[9];
	const _Tp __alam2 = __v[3] - __v[5] - __v[11];
	{
	  const _Tp __alam = __x[2] * __alam1 + __x[8] * __alam2;
	  __cheb24[3] = __cheb12[3] + __alam;
	  __cheb24[21] = __cheb12[3] - __alam;
	}

	{
	  const _Tp __alam = __x[8] * __alam1 - __x[2] * __alam2;
	  __cheb24[9] = __cheb12[9] + __alam;
	  __cheb24[15] = __cheb12[9] - __alam;
	}
      }

      {
	const _Tp __part1 = __x[3] * __v[4];
	const _Tp __part2 = __x[7] * __v[8];
	const _Tp __part3 = __x[5] * __v[6];

	{
	  const _Tp __alam1 = __v[0] + __part1 + __part2;
	  const _Tp __alam2 = __x[1] * __v[2] + __part3 + __x[9] * __v[10];

	  __cheb12[1] = __alam1 + __alam2;
	  __cheb12[11] = __alam1 - __alam2;
	}

	{
	  const _Tp __alam1 = __v[0] - __part1 + __part2;
	  const _Tp __alam2 = __x[9] * __v[2] - __part3 + __x[1] * __v[10];
	  __cheb12[5] = __alam1 + __alam2;
	  __cheb12[7] = __alam1 - __alam2;
	}
      }

      {
	const _Tp __alam = (__x[0] * __v[1] + __x[2] * __v[3] + __x[4] * __v[5]
		    + __x[6] * __v[7] + __x[8] * __v[9] + __x[10] * __v[11]);
	__cheb24[1] = __cheb12[1] + __alam;
	__cheb24[23] = __cheb12[1] - __alam;
      }

      {
	const _Tp __alam = (__x[10] * __v[1] - __x[8] * __v[3] + __x[6] * __v[5]
		    - __x[4] * __v[7] + __x[2] * __v[9] - __x[0] * __v[11]);
	__cheb24[11] = __cheb12[11] + __alam;
	__cheb24[13] = __cheb12[11] - __alam;
      }

      {
	const _Tp __alam = (__x[4] * __v[1] - __x[8] * __v[3] - __x[0] * __v[5]
		    - __x[10] * __v[7] + __x[2] * __v[9] + __x[6] * __v[11]);
	__cheb24[5] = __cheb12[5] + __alam;
	__cheb24[19] = __cheb12[5] - __alam;
      }

      {
	const _Tp __alam = (__x[6] * __v[1] - __x[2] * __v[3] - __x[10] * __v[5]
		    + __x[0] * __v[7] - __x[8] * __v[9] - __x[4] * __v[11]);
	__cheb24[7] = __cheb12[7] + __alam;
	__cheb24[17] = __cheb12[7] - __alam;
      }

      for (int __i = 0; __i < 6; ++__i)
	{
	  const std::size_t __j = 12 - __i;
	  __v[__i] = __fval[__i] - __fval[__j];
	  __fval[__i] = __fval[__i] + __fval[__j];
	}

      {
	const _Tp __alam1 = __v[0] + __x[7] * __v[4];
	const _Tp __alam2 = __x[3] * __v[2];

	__cheb12[2] = __alam1 + __alam2;
	__cheb12[10] = __alam1 - __alam2;
      }

      __cheb12[6] = __v[0] - __v[4];

      {
	const _Tp __alam = __x[1] * __v[1] + __x[5] * __v[3] + __x[9] * __v[5];
	__cheb24[2] = __cheb12[2] + __alam;
	__cheb24[22] = __cheb12[2] - __alam;
      }

      {
	const _Tp __alam = __x[5] * (__v[1] - __v[3] - __v[5]);
	__cheb24[6] = __cheb12[6] + __alam;
	__cheb24[18] = __cheb12[6] - __alam;
      }

      {
	const _Tp __alam = __x[9] * __v[1] - __x[5] * __v[3] + __x[1] * __v[5];
	__cheb24[10] = __cheb12[10] + __alam;
	__cheb24[14] = __cheb12[10] - __alam;
      }

      for (int __i = 0; __i < 3; ++__i)
	{
	  const std::size_t __j = 6 - __i;
	  __v[__i] = __fval[__i] - __fval[__j];
	  __fval[__i] = __fval[__i] + __fval[__j];
	}

      __cheb12[4] = __v[0] + __x[7] * __v[2];
      __cheb12[8] = __fval[0] - __x[7] * __fval[2];

      {
	const _Tp __alam = __x[3] * __v[1];
	__cheb24[4] = __cheb12[4] + __alam;
	__cheb24[20] = __cheb12[4] - __alam;
      }

      {
	const _Tp __alam = __x[7] * __fval[1] - __fval[3];
	__cheb24[8] = __cheb12[8] + __alam;
	__cheb24[16] = __cheb12[8] - __alam;
      }

      __cheb12[0] = __fval[0] + __fval[2];

      {
	const _Tp __alam = __fval[1] + __fval[3];
	__cheb24[0] = __cheb12[0] + __alam;
	__cheb24[24] = __cheb12[0] - __alam;
      }

      __cheb12[12] = __v[0] - __v[2];
      __cheb24[12] = __cheb12[12];

      for (int __i = 1; __i < 12; ++__i)
	__cheb12[__i] *= _Tp{1} / _Tp{6};

      __cheb12[0] *= _Tp{1} / _Tp{12};
      __cheb12[12] *= _Tp{1} / _Tp{12};

      for (int __i = 1; __i < 24; ++__i)
	__cheb24[__i] *= _Tp{1} / _Tp{12};

      __cheb24[0] *= _Tp{1} / _Tp{24};
      __cheb24[24] *= _Tp{1} / _Tp{24};
    }

} // namespace __gnu_test

#endif // QCHEB_INTEGRATE_H
