
#ifndef __FFT_H
#define __FFT_H 1

// FFT
// ===
/*
 | Discrete Fourier Transform can be regarded as evaluating a
 | polynomial of degree N-1 on the powers
 | omega^0, omega, omega^2, ..., omega^(N-1) where omega is a the
 | Nth root of unity.
 |
 | Given a polynomial of even degree
 |
 | p(t) = a_0 + a_1 t + a_2 t^2 + ...
 |
 | we find:
 |
 | p(t) = a_0 + a_2 t^2 + a_4 t^4 + ... + t(a_1 + a_3 t^2 + ...)
 |      = q_e(t^2) + t q_o (t^2)
 |
 | where q_e and q_o are polynomials formed with the even and odd
 | coefficients of p(t). Thus, we get
 |
 | p(1) = q_e(1) + omega^0 q_o(1)
 | p(omega) = q_e(omega^2) + omega^1 q_o(omega^2)
 | p(omega^2) = q_e(omega^4) + omega^2 q_o(omega^4)
 | ...
 |
 | Note how on the RHS the Fourier transforms of q_e and q_o appear. Thus:
 */

#include <numeric>
#include <vector>
#include <complex>
#include <ext/math_const.h>

namespace __gnu_cxx
{

  /**
   * 
   */
  template<typename _Tp>
    class __phase_iterator
    : public std::iterator<std::input_iterator_tag,
                           std::complex<_Tp>,
                           std::ptrdiff_t,
                           const std::complex<_Tp>*,
                           const std::complex<_Tp>&>
    {

    private:

      std::complex<_Tp> __omega_pow_i;
      std::complex<_Tp> __omega_pow_ik;
      std::size_t __k;

      _Tp
      _M_rational_arg(std::size_t __i, std::size_t __m)
      {
#if REPERIOD
	return _Tp(2 * __i) / _Tp(__m);
#else
	const auto _S_2pi = __gnu_cxx::__const_2_pi<_Tp>();
	return _S_2pi * _Tp(__i) / _Tp(__m);
#endif
      }

    public:

      __phase_iterator(_Tp __sign,
                       std::size_t __i,
                       std::size_t __len,
                       bool __past_end = false)
      :
	__omega_pow_i(std::polar(_Tp{1},
			-__sign * _M_rational_arg(__i, __len))),
	__omega_pow_ik(_Tp{1}),
	__k(__past_end ? __len : 0)
      { }

      std::complex<_Tp>
      operator*() const
      { return __omega_pow_ik; }

      __phase_iterator&
      operator++()
      {
	__omega_pow_ik *= __omega_pow_i;
	++__k;
	return *this;
      }

      __phase_iterator
      operator++(int)
      {
	__phase_iterator __dummy(*this);
	++__k;
	return __dummy;
      }

      bool
      operator==(const __phase_iterator& __other) const
      {
	return (this->__omega_pow_i == __other.__omega_pow_i)
	    && (this->__k == __other.__k);
      }

      bool
      operator!=(const __phase_iterator& __other) const
      { return !(*this == __other); }

    }; // __phase_iterator


  /**
   * Discreet Fourier transform on complex data.
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
   * Do Fast Fourier Transform on complex data.
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
   * Do inverse Fast Fourier Transform on complex data.
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
   * Do Fast Fourier Transform on real data.
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
	  __z.emplace_back(__z[0]);
	  fast_fourier_transform(__z); // Use symmetry.  We need N/2 transform.
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
   * Do inverse Fast Fourier Transform on real data.
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
   * Do Fast Fourier Transform on input range.
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
   * Do inverse Fast Fourier Transform on input range.
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

#endif // __FFT_H
