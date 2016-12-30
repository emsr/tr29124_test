/* integration/cquad_integrate.h
 *
 * Copyright (C) 2010 Pedro Gonnet
 * Copyright (C) 2016 Free Software Foundation, Inc.
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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Pedro Gonnet
//
// This file implements the cquad integration scheme.
// Based upon structs in gsl-2.3/integration/cquad.c

#ifndef CQUAD_INTEGRATE_H
#define CQUAD_INTEGRATE_H 1

#include "cquad_const.h"
#include "cquad_workspace.h"

namespace __gnu_test
{

  /**
   * Compute the product of the fx with one of the inverse
   * Vandermonde-like matrices.
   */
  template<typename _Tp>
    void
    _Vinvfx(const _Tp* __fx, _Tp* __c, const int __d)
    {
      switch (__d)
	{
	case 0:
	  for (int __i = 0; __i <= 4; ++__i)
	    {
	      __c[__i] = _Tp{0};
	      for (int __j = 0; __j <= 4; ++__j)
		__c[__i] += V1inv[__i * 5 + __j] * __fx[__j * 8];
	    }
	  break;
	case 1:
	  for (int __i = 0; __i <= 8; ++__i)
	    {
	      __c[__i] = _Tp{0};
	      for (int __j = 0; __j <= 8; ++__j)
		__c[__i] += V2inv[__i * 9 + __j] * __fx[__j * 4];
	    }
	  break;
	case 2:
	  for (int __i = 0; __i <= 16; ++__i)
	    {
	      __c[__i] = _Tp{0};
	      for (int __j = 0; __j <= 16; ++__j)
		__c[__i] += V3inv[__i * 17 + __j] * __fx[__j * 2];
	    }
	  break;
	case 3:
	  for (int __i = 0; __i <= 32; ++__i)
	    {
	      __c[__i] = _Tp{0};
	      for (int __j = 0; __j <= 32; ++__j)
		__c[__i] += V4inv[__i * 33 + __j] * __fx[__j];
	    }
	  break;
	}
    }

  /**
   * Downdate the interpolation given by the n coefficients c
   * by removing the nodes with indices in NaNs.
   */
  template<typename _Tp>
    void
    downdate(_Tp* __c, std::size_t __n, std::size_t __d,
	     std::size_t* __NaN, std::size_t __num_NaNs)
    {
      constexpr std::size_t __bidx[4] = { 0, 6, 16, 34 };
      _Tp __b_new[34], __alpha;

      for (std::size_t __i = 0; __i <= __n + 1; ++__i)
	__b_new[__i] = bee[__bidx[__d] + __i];
      for (std::size_t __i = 0; __i < __num_NaNs; ++__i)
	{
	  __b_new[__n + 1] = __b_new[__n + 1] / Lalpha[__n];
	  __b_new[__n] = (__b_new[__n] + xi[__NaN[__i]] * __b_new[__n + 1])
		       / Lalpha[__n - 1];
	  for (std::size_t __j = __n - 1; __j > 0; --__j)
	    __b_new[__j] = (__b_new[__j] + xi[__NaN[__i]] * __b_new[__j + 1]
			- Lgamma[__j + 1] * __b_new[__j + 2]) / Lalpha[__j - 1];
	  for (std::size_t __j = 0; __j <= __n; ++__j)
	    __b_new[__j] = __b_new[__j + 1];
	  __alpha = __c[__n] / __b_new[__n];
	  for (std::size_t __j = 0; __j < __n; ++__j)
	    __c[__j] -= __alpha * __b_new[__j];
	  __c[__n] = 0;
	  --__n;
	}
    }

  /**
   *
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    cquad_integrate(cquad_workspace<_Tp>& __ws,
		    const _FuncTp& __func,
		    _Tp __a, _Tp __b,
		    _Tp __epsabs, _Tp __epsrel)
    {
      // Some constants that we will need.
      constexpr std::size_t __n[4] = { 4, 8, 16, 32 };
      constexpr std::size_t __skip[4] = { 8, 4, 2, 1 };
      constexpr std::size_t __idx[4] = { 0, 5, 14, 31 };
      constexpr std::size_t __ndiv_max = 20;
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
      const auto _S_inf = std::numeric_limits<_Tp>::infinity();
      constexpr _Tp _S_sqrt2 = M_SQRT2;
      constexpr _Tp __w = _S_sqrt2 / _Tp{2};

      _Tp __result, __abserr;

      // Actual variables (as opposed to constants above).
      bool __split;
      std::size_t __num_NaNs, __NaN[32];
      _Tp __nc, __ncdiff;

      // Check for unreasonable accuracy demands.
      if (__epsabs < _Tp{0} || __epsrel < _Tp{0})
	std::__throw_domain_error("tolerances may not be negative");
      if (__epsabs <= _Tp{0} && __epsrel < _S_eps)
	std::__throw_domain_error("unreasonable accuracy requirement");

      // Create the first interval.
      __ws.initialize_heap();
      auto& __iv = __ws.ivals[0];
      auto __m = (__a + __b) / _Tp{2};
      auto __h = (__b - __a) / _Tp{2};
      __num_NaNs = 0;
      for (std::size_t __i = 0; __i <= __n[3]; ++__i)
	{
	  __iv.fx[__i] = __func(__m + xi[__i] * __h);
	  if (std::isinf(__iv.fx[__i]) || std::isnan(__iv.fx[__i]))
	    {
	      __NaN[__num_NaNs++] = __i;
	      __iv.fx[__i] = _Tp{0};
	    }
	}
      _Vinvfx(__iv.fx, &(__iv.c[__idx[0]]), 0);
      _Vinvfx(__iv.fx, &(__iv.c[__idx[3]]), 3);
      _Vinvfx(__iv.fx, &(__iv.c[__idx[2]]), 2);
      for (std::size_t __i = 0; __i < __num_NaNs; ++__i)
	__iv.fx[__NaN[__i]] = _S_NaN;
      __iv.a = __a;
      __iv.b = __b;
      __iv.depth = 3;
      __iv.rdepth = 1;
      __iv.ndiv = 0;
      __iv.igral = _Tp{2} * __h * __iv.c[__idx[3]] * __w;
      __nc = _Tp{0};
      for (std::size_t __i = __n[2] + 1; __i <= __n[3]; ++__i)
	{
	  auto __temp = __iv.c[__idx[3] + __i];
	  __nc += __temp * __temp;
	}
      __ncdiff = __nc;
      for (std::size_t __i = 0; __i <= __n[2]; ++__i)
	{
	  auto __temp = __iv.c[__idx[2] + __i] - __iv.c[__idx[3] + __i];
	  __ncdiff += __temp * __temp;
	  __nc += __iv.c[__idx[3] + __i] * __iv.c[__idx[3] + __i];
	}
      __ncdiff = std::sqrt(__ncdiff);
      __nc = std::sqrt(__nc);
      __iv.err = __ncdiff * _Tp{2} * __h;
      if ( __ncdiff / __nc > _Tp{0.1} && __iv.err < _Tp{2} * __h * __nc)
	__iv.err = _Tp{2} * __h * __nc;

      auto __igral = __iv.igral;
      auto __igral_final = _Tp{0};
      auto __err = __iv.err;
      auto __err_final = _Tp{0};
      std::size_t __nivals = 1;

#ifdef INTEGRATION_DEBUG
      fprintf(stderr,"\n");
#endif

      // Main loop...
      while (__nivals > 0 && __err > _Tp{0} &&
	     !(__err <= std::abs(__igral) * __epsrel || __err <= __epsabs)
	     && !(__err_final > std::abs(__igral) * __epsrel
		  && __err - __err_final < std::abs(__igral) * __epsrel)
	     && !(__err_final > __epsabs && __err - __err_final < __epsabs))
	{
	  // Put our finger on the interval with the largest error.
	  auto& __iv = __ws.ivals[__ws.heap[0]];
	  __m = (__iv.a + __iv.b) / _Tp{2};
	  __h = (__iv.b - __iv.a) / _Tp{2};

#ifdef INTEGRATION_DEBUG
	  printf
	    ("cquad: processing ival %i (of %i) with [%e,%e] int=%e, err=%e, depth=%i\n",
	     __ws.heap[0], __nivals, __iv.a, __iv.b, __iv.igral, __iv.err, __iv.depth);
#endif
	  // Should we try to increase the degree?
	  if (__iv.depth < 3)
	    {
	      // Keep tabs on some variables.
	      auto __d = ++__iv.depth;

	      // Get the new (missing) function values.
	      for (std::size_t __i = __skip[__d];
			__i <= 32; __i += 2 * __skip[__d])
		__iv.fx[__i] = __func(__m + xi[__i] * __h);
	      __num_NaNs = 0;
	      for (std::size_t __i = 0; __i <= 32; __i += __skip[__d])
		if (std::isinf(__iv.fx[__i]) || std::isnan(__iv.fx[__i]))
		  {
		    __NaN[__num_NaNs++] = __i;
		    __iv.fx[__i] = _Tp{0};
		  }

	      // Compute the new coefficients.
	      _Vinvfx(__iv.fx, &(__iv.c[__idx[__d]]), __d);

	      // Downdate any NaNs.
	      if (__num_NaNs > 0)
		{
		  downdate(&(__iv.c[__idx[__d]]), __n[__d], __d,
			   __NaN, __num_NaNs);
		  for (std::size_t __i = 0; __i < __num_NaNs; ++__i)
		    __iv.fx[__NaN[__i]] = _S_NaN;
		}

	      // Compute the error estimate.
	      __nc = _Tp{0};
	      for (std::size_t __i = __n[__d - 1] + 1; __i <= __n[__d]; ++__i)
		{
		  auto __temp = __iv.c[__idx[__d] + __i];
		  __nc += __temp * __temp;
		}
	      __ncdiff = __nc;
	      for (std::size_t __i = 0; __i <= __n[__d - 1]; ++__i)
		{
		  auto __temp = __iv.c[__idx[__d - 1] + __i]
			      - __iv.c[__idx[__d] + __i];
		  __ncdiff += __temp * __temp;
		  __nc += __iv.c[__idx[__d] + __i] * __iv.c[__idx[__d] + __i];
		}
	      __ncdiff = std::sqrt(__ncdiff);
	      __nc = std::sqrt(__nc);
	      __iv.err = __ncdiff * _Tp{2} * __h;

	      // Compute the local integral.
	      __iv.igral = _Tp{2} * __h * __w * __iv.c[__idx[__d]];

	      // Split the interval prematurely?
	      __split = (__nc > _Tp{0} && __ncdiff / __nc > _Tp{0.1});
	    }
	  else // Maximum degree reached, just split.
	    __split = true;


	  // Should we drop this interval?
	  if ((__m + __h * xi[0]) >= (__m + __h * xi[1])
	      || (__m + __h * xi[31]) >= (__m + __h * xi[32])
	      || __iv.err < std::abs(__iv.igral) * _S_eps * 10)
	    {
#ifdef INTEGRATION_DEBUG
	      printf
		("cquad: dumping ival %i (of %i) with [%e,%e] int=%e, err=%e, depth=%i\n",
		 __ws.heap[0], __nivals, __iv.a, __iv.b, __iv.igral, __iv.err,
		 __iv.depth);
#endif
	      // Keep this interval's contribution.
	      __err_final += __iv.err;
	      __igral_final += __iv.igral;

	      // Swap with the last element on the heap.
	      std::swap(__ws.heap[__nivals - 1], __ws.heap[0]);
	      --__nivals;

	      // Fix up the heap.
	      std::size_t __i = 0;
	      while (2 * __i + 1 < __nivals)
		{
		  // Get the kids.
		  std::size_t __j = 2 * __i + 1;

		  // If the j+1st entry exists and is larger than the jth,
		  // use it instead.
		  if (__j + 1 < __nivals
		      && __ws.ivals[__ws.heap[__j + 1]].err
		       >=__ws.ivals[__ws.heap[__j]].err)
		    ++__j;

		  // Do we need to move the ith entry up?
		  if (__ws.ivals[__ws.heap[__j]].err
			 <= __ws.ivals[__ws.heap[__i]].err)
		    break;
		  else
		    {
		      std::swap(__ws.heap[__i], __ws.heap[__j]);
		      __i = __j;
		    }
		}
	    }
	  else if (__split) // Do we need to split this interval?
	    {
	      // Some values we will need often...
	      auto __d = __iv.depth;

	      // Generate the interval on the left.
	      auto& __ivl = __ws.ivals[__ws.heap[__nivals++]];
	      __ivl.a = __iv.a;
	      __ivl.b = __m;
	      __ivl.depth = 0;
	      __ivl.rdepth = __iv.rdepth + 1;
	      __ivl.fx[0] = __iv.fx[0];
	      __ivl.fx[32] = __iv.fx[16];
	      for (std::size_t __i = __skip[0]; __i < 32; __i += __skip[0])
		__ivl.fx[__i] = __func((__ivl.a + __ivl.b)
			      / _Tp{2} + xi[__i] * __h / _Tp{2});
	      __num_NaNs = 0;
	      for (std::size_t __i = 0; __i <= 32; __i += __skip[0])
		{
		  if (std::isinf(__ivl.fx[__i]) || std::isnan(__ivl.fx[__i]))
		    {
		      __NaN[__num_NaNs++] = __i;
		      __ivl.fx[__i] = _Tp{0};
		    }
		}
	      _Vinvfx(__ivl.fx, __ivl.c, 0);
	      if (__num_NaNs > 0)
		{
		  downdate(__ivl.c, __n[0], 0, __NaN, __num_NaNs);
		  for (std::size_t __i = 0; __i < __num_NaNs; ++__i)
		    __ivl.fx[__NaN[__i]] = _S_NaN;
		}
	      for (std::size_t __i = 0; __i <= __n[__d]; ++__i)
		{
		  __ivl.c[__idx[__d] + __i] = _Tp{0};
		  for (std::size_t __j = __i; __j <= __n[__d]; ++__j)
		    __ivl.c[__idx[__d] + __i] += Tleft[__i * 33 + __j]
						* __iv.c[__idx[__d] + __j];
		}
	      __ncdiff = _Tp{0};
	      for (std::size_t __i = 0; __i <= __n[0]; ++__i)
		{
		  auto __temp = __ivl.c[__i] - __ivl.c[__idx[__d] + __i];
		  __ncdiff += __temp * __temp;
		}
	      for (std::size_t __i = __n[0] + 1; __i <= __n[__d]; ++__i)
		{
		  auto __temp = __ivl.c[__idx[__d] + __i];
		  __ncdiff += __temp * __temp;
		}
	      __ncdiff = std::sqrt(__ncdiff);
	      __ivl.err = __ncdiff * __h;

	      // Check for divergence.
	      __ivl.ndiv = __iv.ndiv + (std::abs (__iv.c[0]) > 0
				      && __ivl.c[0] / __iv.c[0] > 2);
	      if (__ivl.ndiv > __ndiv_max && 2 * __ivl.ndiv > __ivl.rdepth)
		{
		  __result = std::copysign(_S_inf, __igral);
		  return std::make_tuple(__result, __abserr);
		}

	      // Compute the local integral.
	      __ivl.igral = __h * __w * __ivl.c[0];


	      // Generate the interval on the right.
	      auto& __ivr = __ws.ivals[__ws.heap[__nivals++]];
	      __ivr.a = __m;
	      __ivr.b = __iv.b;
	      __ivr.depth = 0;
	      __ivr.rdepth = __iv.rdepth + 1;
	      __ivr.fx[0] = __iv.fx[16];
	      __ivr.fx[32] = __iv.fx[32];
	      for (std::size_t __i = __skip[0]; __i < 32; __i += __skip[0])
		__ivr.fx[__i] = __func((__ivr.a + __ivr.b)
			      / _Tp{2} + xi[__i] * __h / _Tp{2});
	      __num_NaNs = 0;
	      for (std::size_t __i = 0; __i <= 32; __i += __skip[0])
		{
		  if (std::isinf(__ivr.fx[__i]) || std::isnan(__ivr.fx[__i]))
		    {
		      __NaN[__num_NaNs++] = __i;
		      __ivr.fx[__i] = _Tp{0};
		    }
		}
	      _Vinvfx (__ivr.fx, __ivr.c, 0);
	      if (__num_NaNs > 0)
		{
		  downdate(__ivr.c, __n[0], 0, __NaN, __num_NaNs);
		  for (std::size_t __i = 0; __i < __num_NaNs; ++__i)
		    __ivr.fx[__NaN[__i]] = _S_NaN;
		}
	      for (std::size_t __i = 0; __i <= __n[__d]; ++__i)
		{
		  __ivr.c[__idx[__d] + __i] = _Tp{0};
		  for (std::size_t __j = __i; __j <= __n[__d]; ++__j)
		    __ivr.c[__idx[__d] + __i] += Tright[__i * 33 + __j]
						 * __iv.c[__idx[__d] + __j];
		}
	      __ncdiff = _Tp{0};
	      for (std::size_t __i = 0; __i <= __n[0]; ++__i)
		{
		  auto __temp = __ivr.c[__i] - __ivr.c[__idx[__d] + __i];
		  __ncdiff += __temp * __temp;
		}
	      for (std::size_t __i = __n[0] + 1; __i <= __n[__d]; ++__i)
		{
		  auto __temp = __ivr.c[__idx[__d] + __i];
		  __ncdiff += __temp * __temp;
		}
	      __ncdiff = std::sqrt(__ncdiff);
	      __ivr.err = __ncdiff * __h;

	      // Check for divergence.
	      __ivr.ndiv = __iv.ndiv + (std::abs (__iv.c[0]) > 0
				      && __ivr.c[0] / __iv.c[0] > 2);
	      if (__ivr.ndiv > __ndiv_max && 2 * __ivr.ndiv > __ivr.rdepth)
		{
		  __result = std::copysign(_S_inf, __igral);
		  return std::make_tuple(__result, __abserr);
		}

	      // Compute the local integral.
	      __ivr.igral = __h * __w * __ivr.c[0];


	      /* Fix-up the heap: we now have one interval on top
		 that we don't need any more and two new, unsorted
		 ones at the bottom. */

	      // Flip the last interval to the top of the heap and sift down.
	      std::swap(__ws.heap[__nivals - 1], __ws.heap[0]);
	      --__nivals;

	      // Sift this interval back down the heap.
	      std::size_t __i = 0;
	      while (2 * __i + 1 < __nivals - 1)
		{
		  std::size_t __j = 2 * __i + 1;
		  if (__j + 1 < __nivals - 1
		      && __ws.ivals[__ws.heap[__j + 1]].err >=
		      __ws.ivals[__ws.heap[__j]].err)
		    ++__j;
		  if (__ws.ivals[__ws.heap[__j]].err
			 <= __ws.ivals[__ws.heap[__i]].err)
		    break;
		  else
		    {
		      std::swap(__ws.heap[__i], __ws.heap[__j]);
		      __i = __j;
		    }
		}

	      // Now grab the last interval and sift it up the heap.
	      __i = __nivals - 1;
	      while (__i > 0)
		{
		  std::size_t __j = (__i - 1) / 2;
		  if (__ws.ivals[__ws.heap[__j]].err
			 < __ws.ivals[__ws.heap[__i]].err)
		    {
		      std::swap(__ws.heap[__i], __ws.heap[__j]);
		      __i = __j;
		    }
		  else
		    break;
		}
	    }
	  else // Otherwise, just fix-up the heap.
	    {
	      std::size_t __i = 0;
	      while (2 * __i + 1 < __nivals)
		{
		  std::size_t __j = 2 * __i + 1;
		  if (__j + 1 < __nivals
		      && __ws.ivals[__ws.heap[__j + 1]].err
			>= __ws.ivals[__ws.heap[__j]].err)
		    ++__j;
		  if (__ws.ivals[__ws.heap[__j]].err
			<= __ws.ivals[__ws.heap[__i]].err)
		    break;
		  else
		    {
		      std::swap(__ws.heap[__i], __ws.heap[__j]);
		      __i = __j;
		    }
		}
	    }

	  // If the heap is about to overflow, remove the last two intervals.
	  while (__nivals > __ws.size() - 2)
	    {
	      __iv = __ws.ivals[__ws.heap[__nivals - 1]];

#ifdef INTEGRATION_DEBUG
	      printf
		("cquad: dumping ival %i (of %i) with [%e,%e] int=%e, err=%e, depth=%i\n",
		 __ws.heap[0], __nivals, __iv.a, __iv.b, __iv.igral, __iv.err,
		 __iv.depth);
#endif
	      __err_final += __iv.err;
	      __igral_final += __iv.igral;
	      --__nivals;
	    }

	  // Collect the value of the integral and error.
	  __igral = __igral_final;
	  __err = __err_final;
	  for (std::size_t __i = 0; __i < __nivals; ++__i)
	    {
	      __igral += __ws.ivals[__ws.heap[__i]].igral;
	      __err += __ws.ivals[__ws.heap[__i]].err;
	    }
	}

      // Dump the contents of the heap.
#ifdef INTEGRATION_DEBUG
      for (std::size_t __i = 0; __i < __nivals; ++__i)
	{
	  auto __iv = __ws.ivals[__ws.heap[__i]];
	  printf
	    ("cquad: ival %i (%i) with [%e,%e], int=%e, err=%e, depth=%i, rdepth=%i\n",
	     __i, __ws.heap[__i], __iv.a, __iv.b, __iv.igral, __iv.err, __iv.depth,
	     __iv.rdepth);
	}
#endif

      return std::make_tuple(__igral, __err);
    }

} // namespace __gnu_test

#endif // CQUAD_INTEGRATE_H
