/* integration/glfixed_integrate.tcc
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

#ifndef GLFIXED_INTEGRATE_TCC
#define GLFIXED_INTEGRATE_TCC 1

namespace __gnu_test
{

  template<typename _FuncTp, typename _Tp>
    _Tp
    glfixed_integrate(const gauss_legendre_table<_Tp>& __t,
		      const _FuncTp& __func,
		      _Tp __a,
		      _Tp __b)
    {
      const int __n = __t.order;

      auto __m = (__n + 1) >> 1;
      auto __A = (__b - __a) / _Tp{2};
      auto __B = (__b + __a) / _Tp{2};

      if (__n & 1) // n is odd.
	{
	  auto __sum = __t.wt(0) * __func(__B);
	  for (int __i = 1; __i < __m; ++__i)
	    {
	      auto __Ax = __A * __t.pt(__i);
	      __sum += __t.wt(__i) * (__func(__B + __Ax) + __func(__B - __Ax));
	    }
	  return __A * __sum;
	}
      else // n is even.
	{
	  auto __sum = _Tp{0};
	  for (int __i = 0; __i < __m; ++__i)
	    {
	      auto __Ax = __A * __t.pt(__i);
	      __sum += __t.wt(__i) * (__func(__B + __Ax) + __func(__B - __Ax));
	    }
	  return __A * __sum;
	}
    }

} // namespace __gnu_test

#endif // GLFIXED_INTEGRATE_TCC
