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


#ifndef FOURIER_TRANSFORM_H
#define FOURIER_TRANSFORM_H 1

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

      std::complex<_Tp> _M_omega_pow_i;
      std::complex<_Tp> _M_omega_pow_ik;
      std::size_t _M_k;

      static _Tp
      _S_rational_arg(std::size_t __i, std::size_t __m)
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
	_M_omega_pow_i(std::polar(_Tp{1},
			-__sign * _S_rational_arg(__i, __len))),
	_M_omega_pow_ik(_Tp{1}),
	_M_k(__past_end ? __len : 0)
      { }

      std::complex<_Tp>
      operator*() const
      { return _M_omega_pow_ik; }

      __phase_iterator&
      operator++()
      {
	++_M_k;
	_M_omega_pow_ik *= _M_omega_pow_i;
	return *this;
      }

      __phase_iterator
      operator++(int)
      {
	__phase_iterator __dummy(*this);
	++_M_k;
	return __dummy;
      }

      _Tp
      cos() const
      { return std::real(this->_M_omega_pow_ik); }

      _Tp
      sin() const
      { return std::imag(this->_M_omega_pow_ik); }

      bool
      operator==(const __phase_iterator& __other) const
      {
	return (this->_M_omega_pow_i == __other._M_omega_pow_i)
	    && (this->_M_k == __other._M_k);
      }

      bool
      operator!=(const __phase_iterator& __other) const
      { return !(*this == __other); }

    }; // __phase_iterator

/**
 * Fast Fourier Transform
 *
 * Discrete Fourier Transform can be regarded as evaluating a
 * polynomial of degree N-1 on the powers
 * @f$ \omega^0, \omega, \omega^2, ..., \omega^(N-1) @f$
 * where @f$ \omega @f$ is the Nth root of unity.
 *
 * Given a polynomial of even degree
 * @f[
 *    p(t) = a_0 + a_1 t + a_2 t^2 + ...
 * @f]
 * we find:
 * @f[
 *    p(t) = a_0 + a_2 t^2 + a_4 t^4 + ... + t(a_1 + a_3 t^2 + ...)
 *         = q_{even}(t^2) + t q_{odd} (t^2)
 * @f]
 * where q_e and q_o are polynomials formed with the even and odd
 * coefficients of p(t). Thus, we get
 *
 * p(1) = q_{even}(1) + \omega^0 q_{odd}(1)
 * p(\omega) = q_{even}(\omega^2) + \omega^1 q_{odd}(\omega^2)
 * p(\omega^2) = q_{even}(\omega^4) + \omega^2 q_{odd}(\omega^4)
 * ...
 *
 * Note how on the RHS the Fourier transforms of q_e and q_o appear. Thus:
 */

  /**
   * Discrete Fourier transform type.
   */
  template<typename _Tp>
    class fourier_transform_t
    {
    private:

      std::vector<std::complex<_Tp>> _M_xform;

    public:

      fourier_transform_t(std::size_t __n)
      : _M_xform{}
      { this->_M_xform.reserve(__n / 2 + 1); }

      fourier_transform_t(std::vector<_Tp> __data);

      std::size_t
      size() const
      { return 2 * this->_M_xform.size() - 2; }

      std::complex<_Tp>
      operator[](std::size_t __k) const
      {
	if (__k < this->_M_xform.size())
	  return this->_M_xform[__k];
	else
	  return std::conj(this->_M_xform[this->_M_xform.size() - __k]);
	// This is real array indexing.
	//if (__k < len / 2)
	//	? std::complex(xform[2 * __k], xform[2 * __k + 1])
	//	: std::complex(xform[2 * len - 2 * i - 2],
	//		      -xform[2 * len - 2 * __k - 1]);
      }
    };

  /**
   * Discrete Fourier transform specialization for transform of complex data.
   */
  template<typename _Tp>
    class fourier_transform_t<std::complex<_Tp>>
    {
    private:

      std::vector<std::complex<_Tp>> _M_xform;

    public:

      fourier_transform_t(std::size_t __n)
      : _M_xform{}
      { this->_M_xform.reserve(__n / 2 + 1); }

      fourier_transform_t(const std::vector<std::complex<_Tp>>& __data);

      std::size_t
      size() const
      { return 2 * this->_M_xform.size() - 2; }

      std::complex<_Tp>
      operator[](std::size_t __k) const
      { return this->_M_xform[__k]; }
    };

  /**
   * Discrete Fourier transform on complex data.
   */
  template<typename _Tp>
    void
    __discrete_fourier_transform(bool __do_forward,
				 std::vector<std::complex<_Tp>>& __z);

  /**
   * Fast Fourier Transform on complex data.
   */
  template<typename _Tp>
    void
    fast_fourier_transform(std::vector<std::complex<_Tp>>& __z);

  /**
   * Inverse Fast Fourier Transform on complex data.
   */
  template<typename _Tp>
    void
    inv_fast_fourier_transform(std::vector<std::complex<_Tp>>& __z);

  /**
   * Fast Fourier Transform on real data.
   */
  template<typename _Tp>
    void
    fast_fourier_transform(std::vector<_Tp>& __x);

  /**
   * Inverse Fast Fourier Transform on real data.
   */
  template<typename _Tp>
    void
    inv_fast_fourier_transform(std::vector<_Tp>& __x);

  /**
   * Fast Fourier Transform on input range.
   */
  template <typename _CmplxIter>
    void
    fast_fourier_transform(const _CmplxIter& __from, const _CmplxIter& __to);

  /**
   * Inverse Fast Fourier Transform on input range.
   */
  template <typename _CmplxIter>
    void
    inv_fast_fourier_transform(const _CmplxIter& __from, const _CmplxIter& __to);

} // namespace __gnu_cxx

#include "fourier_transform.tcc"

#endif // FOURIER_TRANSFORM_H
