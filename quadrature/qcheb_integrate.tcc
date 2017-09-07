/* quadrature/qcheb_integrate.tcc
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

#ifndef QCHEB_INTEGRATE_TCC
#define QCHEB_INTEGRATE_TCC 1

#include <array>

namespace __gnu_ext
{

  template<typename _FuncTp, typename _Tp>
    void
    qcheb_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		    std::array<_Tp, 13>& __cheb12,
		    std::array<_Tp, 25>& __cheb24)
    {
      _Tp __fval[25], __v[12];

      // These are the values of cos(pi*k/24) for k=1..11 needed for the
      // Chebyshev expansion of f(x).  These are the zeros of the Chebyshev
      // finction of the second kind of order 23: U_23(x).

      constexpr _Tp
      __x[11]
      {
	9.914448613738104111442846968605486e-01L,
	9.659258262890682867486612158530536e-01L,
	9.238795325112867561257834975394469e-01L,
	8.660254037844386467595427060757126e-01L,
	7.933533402912351645734146973742314e-01L,
	7.071067811865475243919762573395221e-01L,
	6.087614290087206394044894932434070e-01L,
	4.999999999999999999855184455596035e-01L,
	3.826834323650897717110798781478690e-01L,
	2.588190451025207623287087436359508e-01L,
	1.305261922200515915256103766723547e-01L,
      };

      const auto __center = (__b + __a) / _Tp{2};
      const auto __half_length = (__b - __a) / _Tp{2};

      __fval[0] = __func(__b) / _Tp{2};
      __fval[12] = __func(__center);
      __fval[24] = __func(__a) / _Tp{2};

      for (int __i = 1; __i < 12; ++__i)
	{
	  const std::size_t __j = 24 - __i;
	  const auto __u = __half_length * __x[__i - 1];
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
	const auto __alam1 = __v[0] - __v[8];
	const auto __alam2 = __x[5] * (__v[2] - __v[6] - __v[10]);
	__cheb12[3] = __alam1 + __alam2;
	__cheb12[9] = __alam1 - __alam2;
      }

      {
	const auto __alam1 = __v[1] - __v[7] - __v[9];
	const auto __alam2 = __v[3] - __v[5] - __v[11];
	{
	  const auto __alam = __x[2] * __alam1 + __x[8] * __alam2;
	  __cheb24[3] = __cheb12[3] + __alam;
	  __cheb24[21] = __cheb12[3] - __alam;
	}

	{
	  const auto __alam = __x[8] * __alam1 - __x[2] * __alam2;
	  __cheb24[9] = __cheb12[9] + __alam;
	  __cheb24[15] = __cheb12[9] - __alam;
	}
      }

      {
	const auto __part1 = __x[3] * __v[4];
	const auto __part2 = __x[7] * __v[8];
	const auto __part3 = __x[5] * __v[6];

	{
	  const auto __alam1 = __v[0] + __part1 + __part2;
	  const auto __alam2 = __x[1] * __v[2] + __part3 + __x[9] * __v[10];

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
	const auto __alam = (__x[0] * __v[1] + __x[2] * __v[3]
			  + __x[4] * __v[5] + __x[6] * __v[7]
			  + __x[8] * __v[9] + __x[10] * __v[11]);
	__cheb24[1] = __cheb12[1] + __alam;
	__cheb24[23] = __cheb12[1] - __alam;
      }

      {
	const auto __alam = (__x[10] * __v[1] - __x[8] * __v[3]
			  + __x[6] * __v[5] - __x[4] * __v[7]
			  + __x[2] * __v[9] - __x[0] * __v[11]);
	__cheb24[11] = __cheb12[11] + __alam;
	__cheb24[13] = __cheb12[11] - __alam;
      }

      {
	const auto __alam = (__x[4] * __v[1] - __x[8] * __v[3]
			  - __x[0] * __v[5] - __x[10] * __v[7]
			  + __x[2] * __v[9] + __x[6] * __v[11]);
	__cheb24[5] = __cheb12[5] + __alam;
	__cheb24[19] = __cheb12[5] - __alam;
      }

      {
	const auto __alam = (__x[6] * __v[1] - __x[2] * __v[3]
			  - __x[10] * __v[5] + __x[0] * __v[7]
			  - __x[8] * __v[9] - __x[4] * __v[11]);
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
	const auto __alam1 = __v[0] + __x[7] * __v[4];
	const auto __alam2 = __x[3] * __v[2];

	__cheb12[2] = __alam1 + __alam2;
	__cheb12[10] = __alam1 - __alam2;
      }

      __cheb12[6] = __v[0] - __v[4];

      {
	const auto __alam = __x[1] * __v[1] + __x[5] * __v[3] + __x[9] * __v[5];
	__cheb24[2] = __cheb12[2] + __alam;
	__cheb24[22] = __cheb12[2] - __alam;
      }

      {
	const auto __alam = __x[5] * (__v[1] - __v[3] - __v[5]);
	__cheb24[6] = __cheb12[6] + __alam;
	__cheb24[18] = __cheb12[6] - __alam;
      }

      {
	const auto __alam = __x[9] * __v[1] - __x[5] * __v[3] + __x[1] * __v[5];
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
	const auto __alam = __x[3] * __v[1];
	__cheb24[4] = __cheb12[4] + __alam;
	__cheb24[20] = __cheb12[4] - __alam;
      }

      {
	const auto __alam = __x[7] * __fval[1] - __fval[3];
	__cheb24[8] = __cheb12[8] + __alam;
	__cheb24[16] = __cheb12[8] - __alam;
      }

      __cheb12[0] = __fval[0] + __fval[2];

      {
	const auto __alam = __fval[1] + __fval[3];
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

} // namespace __gnu_ext

#endif // QCHEB_INTEGRATE_TCC
