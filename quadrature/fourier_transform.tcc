// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2017 Free Software Foundation, Inc.
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
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

#ifndef FOURIER_TRANSFORM_TCC
#define FOURIER_TRANSFORM_TCC 1

#include <numeric>
#include <vector>
#include <complex>
#include <ext/math_const.h>

namespace __gnu_cxx
{

  /**
   * Discrete Fourier transform on complex data.
   */
  template<typename _Tp>
    void
    __discrete_fourier_transform(bool __do_forward,
				 std::vector<std::complex<_Tp>>& __z)
    {
      const auto __sign(__do_forward ? _Tp{+1} : _Tp{-1});
      const auto __len = __z.size();

      std::vector<std::complex<_Tp>> __result;
      __result.reserve(__len);

      // Do matrix multiplication.
      for (std::size_t __i = 0; __i < __len; ++__i)
	{
	  __phase_iterator __coefficient_iter(__sign, __i, __len);
	  __result.push_back(std::inner_product(__z.begin(), __z.end(),
			     __coefficient_iter, std::complex<_Tp>{0}));
	}

      // Rescale if forward.
      if (__do_forward)
	{
	  const auto __norm = _Tp{1} / _Tp(__len);
	  for (std::size_t __i = 0; __i < __len; ++__i)
	    __result[__i] *= __norm;
	}

      // Copy data back (swap is fast!).
      __z.swap(__result);
    }

  /**
   * Fast Fourier Transform on complex data.
   */
  template<typename _Tp>
    void
    fast_fourier_transform(std::vector<std::complex<_Tp>>& __z)
    {
      const auto __len = __z.size();
      if (__len % 2 == 1) // Too bad, we're odd.
	__discrete_fourier_transform(true, __z);
      else // Good, we're even.
	{
	  // Let's divide and conquer!
	  std::vector<std::complex<_Tp>> __odd;
	  std::vector<std::complex<_Tp>> __even;
	  const auto __halflen = __len / 2;
	  __odd.reserve(__halflen);
	  __even.reserve(__halflen);
	  for (auto __run = __z.cbegin(); __run != __z.cend(); ++__run)
	    {
              __even.push_back(*__run);
              ++__run;
              __odd.push_back(*__run);
	    }
	  fast_fourier_transform(__even);
	  fast_fourier_transform(__odd);
	  __phase_iterator __omega_iter(_Tp{1}, 1, __len);
	  for (std::size_t __i = 0, __j = __halflen; __i < __halflen;
		++__i, ++__j, ++__omega_iter)
	    {
              __z[__i] = (__even[__i] + *__omega_iter * __odd[__i]) / _Tp{2};
              // The next line works because omega^(length/2) = -1.
              __z[__j] = (__even[__i] - *__omega_iter * __odd[__i]) / _Tp{2};
	    }
	}
    }

  /**
   * Inverse Fast Fourier Transform on complex data.
   */
  template<typename _Tp>
    void
    inv_fast_fourier_transform(std::vector<std::complex<_Tp>>& __z)
    {
      const std::size_t __len = __z.size();
      if (__len % 2 == 1) // Too bad, we're odd.
	__discrete_fourier_transform(false, __z);
      else // Good, we are even.
	{
	  // Let's divide and conquer!
	  std::vector<std::complex<_Tp>> __odd;
	  std::vector<std::complex<_Tp>> __even;
	  const auto __halflen = __len / 2;
	  __odd.reserve(__halflen);
	  __even.reserve(__halflen);
	  for (auto __run = __z.cbegin(); __run != __z.cend(); ++__run)
	    {
	      __even.push_back(*__run);
	      ++__run;
	      __odd.push_back(*__run);
	    }
	  inv_fast_fourier_transform(__even);
	  inv_fast_fourier_transform(__odd);
	  __phase_iterator __omega_iter(_Tp{-1}, 1, __len);
	  for (std::size_t __i = 0, __j = __halflen; __i < __halflen;
		++__i, ++__j, ++__omega_iter)
	    {
	      __z[__i] = __even[__i] + *__omega_iter * __odd[__i];
	      // The next line works because omega^(length/2) = -1.
	      __z[__j] = __even[__i] - *__omega_iter * __odd[__i];
	    }
	}
    }

  /**
   * Fast Fourier Transform on real data.
   */
  template<typename _Tp>
    void
    fast_fourier_transform(std::vector<_Tp>& __x)
    {
      const auto __len = __x.size();
      if (__len % 2 == 1) // Too bad, we're odd.
	std::__throw_domain_error("fast_fourier_transform: "
				  "Real data must have even length.");
      else // Good, we're even.
	{
	  // @todo I really want to just make a view to a complex vector.
	  const auto __halflen = __len / 2;
	  std::vector<std::complex<_Tp>> __z;
	  __z.reserve(__halflen);
	  for (std::size_t __i = 0; __i < __halflen; ++__i)
	    __z.emplace_back(__x[2 * __i], __x[2 * __i + 1]);
	  fast_fourier_transform(__z);
	  __z.emplace_back(__z[0]); // Use symmetry.  We need N/2 transform.
	  const auto _S_i2 = std::complex<_Tp>{0, 2};
	  __phase_iterator __omega_iter(_Tp{+1}, 1, __halflen);
	  for (std::size_t __i = 0; __i < __halflen; ++__i)
	    {
	      const auto __z1 = __z[__i];
	      const auto __z2 = std::conj(__z[__halflen - __i]);
	      const auto __f = (__z1 + __z2) / _Tp{2}
			     + *__omega_iter * (__z1 - __z2) / _S_i2;
	      __x[2 * __i] = __f.real();
	      __x[2 * __i + 1] = __f.imag();
	    }
	}
    }

  /**
   * Inverse Fast Fourier Transform on real data.
   */
  template<typename _Tp>
    void
    inv_fast_fourier_transform(std::vector<_Tp>& __x)
    {
      const auto __len = __x.size();
      if (__len % 2 == 1) // Too bad, we're odd.
	std::__throw_domain_error("inv_fast_fourier_transform: "
				  "Real data must have even length.");
      else // Good, we're even.
	{
	  // @todo I really want to just make a view to a complex vector.
	  const auto __halflen = __len / 2;
	  std::vector<std::complex<_Tp>> __z;
	  __z.reserve(__halflen);
	  const auto _S_i2 = std::complex<_Tp>{0, 2};
	  __phase_iterator __omega_iter(_Tp{-1}, 1, __halflen);
	  for (std::size_t __i = 0; __i < __halflen; ++__i)
	    {
	      const auto __z1 = std::complex<_Tp>(__x[2 * __i], __x[2 * __i + 1]);
	      const auto __z2 = std::complex<_Tp>(__x[__len - 2 * __i - 2], -__x[__len - 2 * __i - 1]);
	      const auto __ze = (__z1 + __z2) / _Tp{2};
	      const auto __zo = -*__omega_iter * (__z1 - __z2) / _S_i2;
	      __z.emplace_back(__ze + __zo);
	    }
	  inv_fast_fourier_transform(__z);
	  for (std::size_t __i = 0; __i < __halflen; ++__i)
	    {
	      __x[2 * __i] = __z[__i].real();
	      __x[2 * __i + 1] = __z[__i].imag();
	    }
	}
    }

  /**
   * Fast Sine Transform on real data.
   */
  template <typename _Tp>
    void
    fast_sine_transform(std::vector<_Tp>& __x)
    {
      const auto __len = __x.size();
      const auto __halflen = __len / 2;
      const auto __n2 = __len - 1;
      __phase_iterator __omega_iter(_Tp{+1}, 1, __len);
      __x[0] = _Tp{0};
      for (std::size_t __k = 1; __k <= __halflen; ++__k, ++__omega_iter)
	{
	  const auto __y1 = __omega_iter.sin()
			  * (__x[__k] + __x[__n2 - __k]);
	  const auto __y2 = (__x[__k] - __x[__n2 - __k]) / _Tp{2};
	  __x[__k] = __y1 + __y2;
	  __x[__n2 - __k] = __y1 - __y2;
	}
      fast_fourier_transform(__x);
      __x[0] *= _Tp{0.5L};
      __x[1] = _Tp{0};
      auto __sum = _Tp{0};
      for (std::size_t __i = 2; __i < __halflen; __i += 2)
	{
	  __sum += std::exchange(__x[__i], __x[__i + 1]);
	  __x[__i + 1] = __sum;
	}
    }

  /**
   * Fast Sine Transform on real data.
   */
  template <typename _Tp>
    void
    inv_fast_sine_transform(std::vector<_Tp>& __x)
    {
      fast_sine_transform(__x);
      const auto __norm = _Tp{2} / __x.size();
      for (std::size_t __i = 0; __i < __x.size(); ++__i)
	__x[__i] *= __norm;
    }

  /**
   * Fast Fourier Transform on input range.
   */
  template <typename _CmplxIter>
    void
    fast_fourier_transform(const _CmplxIter& __from, const _CmplxIter& __to)
    {
      using _Cmplx = typename _CmplxIter::value_type;
      using _Tp = typename _Cmplx::value_type;
      std::vector<std::complex<_Tp>> __z(__from, __to);
      fast_fourier_transform(__z);
      std::copy(__z.begin(), __z.end(), __from);
    }

  /**
   * Inverse Fast Fourier Transform on input range.
   */
  template <typename _CmplxIter>
    void
    inv_fast_fourier_transform(const _CmplxIter& __from, const _CmplxIter& __to)
    {
      using _Cmplx = typename _CmplxIter::value_type;
      using _Tp = typename _Cmplx::value_type;
      std::vector<std::complex<_Tp>> __z(__from, __to);
      inv_fast_fourier_transform(__z);
      std::copy(__z.begin(), __z.end(), __from);
    }

} // namespace __gnu_cxx

#endif // FOURIER_TRANSFORM_TCC
