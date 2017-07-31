// Special functions -*- C++ -*-

// Copyright (C) 2016-2017 Free Software Foundation, Inc.
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
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/sf_theta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SF_THETA_TCC
#define _GLIBCXX_BITS_SF_THETA_TCC 1

#pragma GCC system_header

#include <vector>
#include <tuple>
#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Compute and return the exponential @f$ \theta_2 @f$ function
   * by series expansion:
   * @f[
   *    \theta_2(\nu, x) = \frac{1}{\sqrt{\pi x}}
   *                       \sum_{k=-\infty}^{\infty}(-1)^k e^{-(\nu+k)^2/x}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_2_sum(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      auto __sum = std::exp(-__nu * __nu / __x);
      auto __sign = _Tp{-1};
      for (auto __k = 1; __k < 20; ++__k)
	{
	  const auto __nup = __nu + _Tp(__k);
	  const auto __termp = __sign * std::exp(-__nup * __nup / __x);
	  const auto __num = __nu - _Tp(__k);
	  const auto __termm = __sign * std::exp(-__num * __num / __x);
	  __sum += __termp + __termm;
	  __sign = -__sign;
	  if (std::abs(__termp) < _S_eps * std::abs(__sum)
	   && std::abs(__termm) < _S_eps * std::abs(__sum))
	    break;
	}
      return __sum / std::sqrt(_S_pi * __x);
    }

  /**
   * Compute and return the exponential @f$ \theta_3 @f$ function
   * by series expansion:
   * @f[
   *    \theta_3(\nu, x) = \frac{1}{\sqrt{\pi x}}
   *                       \sum_{k=-\infty}^{\infty} e^{-(\nu+k)^2/x}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_3_sum(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      auto __sum = std::exp(-__nu * __nu / __x);
      for (auto __k = 1; __k < 20; ++__k)
	{
	  const auto __nup = __nu + _Tp(__k);
	  const auto __termp = std::exp(-__nup * __nup / __x);
	  const auto __num = __nu - _Tp(__k);
	  const auto __termm = std::exp(-__num * __num / __x);
	  __sum += __termp + __termm;
	  if (std::abs(__termp) < _S_eps * std::abs(__sum)
	   && std::abs(__termm) < _S_eps * std::abs(__sum))
	    break;
	}
      return __sum / std::sqrt(_S_pi * __x);
    }

  /**
   * Compute and return the exponential @f$ \theta_2 @f$ function
   * by asymptotic series expansion:
   * @f[
   *    \theta_2(\nu, x) = 2\sum_{k=0}^{\infty} e^{-((k+1/2)\pi)^2 x}
   *                        \cos((2k+1)\nu\pi)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_2_asymp(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      auto __sum = _Tp{0};
      for (auto __k = 0; __k < 20; ++__k)
	{
	  const auto __thing = _Tp(2 * __k + 1) * _S_pi;
	  const auto __cosarg = __nu * __thing;
	  const auto __exparg = __thing * __thing * __x / _Tp{4};
	  const auto __term = std::exp(-__exparg) * std::cos(__cosarg);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Tp{2} * __sum;
    }

  /**
   * Compute and return the exponential @f$ \theta_3 @f$ function
   * by asymptotic series expansion:
   * @f[
   *    \theta_3(\nu, x) = 1 + 2\sum_{k=1}^{\infty} e^{-(k\pi)^2 x}
   *                           \cos(2k\nu\pi)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_3_asymp(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      auto __sum = _Tp{0};
      for (auto __k = 1; __k < 20; ++__k)
	{
	  const auto __thing = _Tp(2 * __k) * _S_pi;
	  const auto __cosarg = __nu * __thing;
	  const auto __exparg = __thing * __thing * __x / _Tp{4};
	  const auto __term = std::exp(-__exparg) * std::cos(__cosarg);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Tp{1} + _Tp{2} * __sum;
    }

  /**
   * Return the exponential theta-2 function of period @c nu and argument @c x.
   *
   * The exponential theta-2 function is defined by
   * @f[
   *    \theta_2(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{k=-\infty}^{+\infty}
   *    (-1)^k \exp\left( \frac{-(\nu + k)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_2(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__x) <= _Real{1} / _S_pi)
	return __theta_2_sum(__nu, __x);
      else
	return __theta_2_asymp(__nu, __x);
    }

  /**
   * Return the exponential theta-1 function of period @c nu and argument @c x.
   *
   * The exponential theta-1 function is defined by
   * @f[
   *    \theta_1(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{k=-\infty}^{+\infty}
   *    (-1)^k \exp\left( \frac{-(\nu + k - 1/2)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_1(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else if (__gnu_cxx::__fp_is_zero(__x))
	return _Tp{0};
      else
	return __theta_2(__nu - _Tp{0.5L}, __x);
    }

  /**
   * Return the exponential theta-3 function of period @c nu and argument @c x.
   *
   * The exponential theta-3 function is defined by
   * @f[
   *    \theta_3(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{k=-\infty}^{+\infty}
   *    \exp\left( \frac{-(\nu+k)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 1) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_3(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__x) <= _Real{1} / _S_pi)
	return __theta_3_sum(__nu, __x);
      else
	return __theta_3_asymp(__nu, __x);
    }

  /**
   * Return the exponential theta-4 function of period @c nu and argument @c x.
   *
   * The exponential theta-4 function is defined by
   * @f[
   *    \theta_4(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{k=-\infty}^{+\infty}
   *    (-1)^k \exp\left( \frac{-(\nu + k)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_4(_Tp __nu, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else
	return __theta_3(__nu + _Tp{0.5L}, __x);
    }

  /**
   * Use MacLaurin series to calculate the elliptic nome
   * given the elliptic argument k.
   * @f[
   *    q(k) = exp\left(-\pi\frac{K(k')}{K(k)}\right)
   * @f]
   * where @f$ k' = \sqrt{1 - k^2} @f$ is the complementary elliptic argument
   * and @f$  @f$ is the Legendre elliptic integral of the first kind.
   */
  template<typename _Tp>
    _Tp
    __ellnome_series(_Tp __k)
    {
      const auto __m = __k * __k; 
      return __m * ((_Tp{1} / _Tp{16})
	   + __m * ((_Tp{1} / _Tp{32})
	   + __m * ((_Tp{21} / _Tp{1024})
	   + __m * ((_Tp{31} / _Tp{2048})
	   + __m * (_Tp{6257} / _Tp{524288})))));
    }

  /**
   * Use the arithmetic-geometric mean to calculate the elliptic nome
   * given the elliptic argument k.
   * @f[
   *    q(k) = exp\left(-\pi\frac{K(k')}{K(k)}\right)
   * @f]
   * where @f$ k' = \sqrt{1 - k^2} @f$ is the complementary elliptic argument
   * and @f$  @f$ is the Legendre elliptic integral of the first kind.
   */
  template<typename _Tp>
    _Tp
    __ellnome_k(_Tp __k)
    {
      const auto _S_pi = _Tp{3.1415926535897932384626433832795029L};
      const auto __kp = std::sqrt(_Tp{1} - __k * __k);
      const auto __K = __comp_ellint_1(__k);
      const auto __Kp = __comp_ellint_1(__kp);
      return std::exp(-_S_pi * __Kp / __K);
    }

  /**
   * Return the elliptic nome given the modulus @c k.
   * @f[
   *    q(k) = exp\left(-\pi\frac{K(k')}{K(k)}\right)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __ellnome(_Tp __k)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      if (__isnan(__k))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__ellnome:"
				      " argument k out of range"));
      else if (__k < std::pow(_Tp{67} * _S_eps, _Tp{0.125L}))
	return __ellnome_series(__k);
      else
	return __ellnome_k(__k);
    }

  /**
   * Return the Neville @f$ \theta_s @f$ function
   * @f[
   *  \theta_s(k,x) = \sqrt{\frac{\pi}{2 k k' K(k)}}
   *                  \theta_1\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_s(_Tp __k, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::abs(__x));

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_s:"
				      " argument k out of range"));
      else
	{
	  const auto __kc = std::sqrt(_Tp{1} - __k * __k);
	  const auto _Kk = __comp_ellint_1(__k);
	  const auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / (__k * __kc * _Kk))
	       * __theta_1(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_c @f$ function
   * @f[
   *    \theta_c(k,x) = \sqrt{\frac{\pi}{2 k K(k)}}
   *                  \theta_1\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_c(_Tp __k, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::abs(__x));

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_c:"
				      " argument k out of range"));
      else
	{
	  const auto _Kk = __comp_ellint_1(__k);
	  const auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / (__k * _Kk))
	       * __theta_2(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_d @f$ function
   * @f[
   *    \theta_d(k,x) = \sqrt{\frac{\pi}{2K(k)}}
   *                  \theta_3\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_d(_Tp __k, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::abs(__x));

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_d:"
				      " argument k out of range"));
      else
	{
	  const auto _Kk = __comp_ellint_1(__k);
	  const auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / _Kk)
	       * __theta_3(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_n @f$ function
   *
   * The Neville theta-n function is defined by
   * @f[
   *  \theta_n(k,x) = \sqrt{\frac{\pi}{2k'K(k)}}
   *                  \theta_4\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __theta_n(_Tp __k, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::abs(__x));

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_n:"
				      " argument k out of range"));
      else
	{
	  const auto __kc = std::sqrt(_Tp{1} - __k * __k);
	  const auto _Kk = __comp_ellint_1(__k);
	  const auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / (__kc * _Kk))
	       * __theta_4(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * A struct representing the Jacobi and Weierstrass lattice.
   * The two types for the frequencies and the subsequent type calculus
   * allow us to treat the rectangulr lattice (real nome, pure imaginary
   * lattice parameter) specially.
   */
  template<typename _Tp_Omega1, typename _Tp_Omega3 = std::complex<_Tp_Omega1>>
    struct __jacobi_lattice_t
    {
      static_assert(__gnu_cxx::is_complex_v<_Tp_Omega1>
		 || __gnu_cxx::is_complex_v<_Tp_Omega3>,
		    "One frequecy type must be complex.");
      using _Real_Omega1 = __num_traits_t<_Tp_Omega1>;
      using _Real_Omega3 = __num_traits_t<_Tp_Omega3>;
      using _Real = __gnu_cxx::__promote_fp_t<_Real_Omega1, _Real_Omega3>;
      using _Cmplx = std::complex<_Real>;
      using _Tp_Nome = std::conditional_t<__gnu_cxx::is_complex_v<_Tp_Omega1>
				       && __gnu_cxx::is_complex_v<_Tp_Omega3>,
					  _Cmplx, _Real>;

      /**
       * A struct representing a complex scalar lattice parameter
       * or half period ratio.
       */
      struct __tau_t
      {
	_Cmplx __val;

	explicit __tau_t(_Cmplx __tau)
	: __val(__tau)
	{ }
      };

      /**
       * A struct representing a complex argument reduced
       * to the 'central' lattice cell.
       */
      struct __arg_t
      {
	int __m;
	int __n;
	_Cmplx __z;
      };

      /// Construct the lattice from two complex lattice frequencies.
      __jacobi_lattice_t(const _Tp_Omega1& __omega1,
			 const _Tp_Omega3& __omega3)
      : _M_omega_1(__omega1),
        _M_omega_3(__omega3)
      {
	if (__isnan(_M_omega_1) || __isnan(_M_omega_3))
	  std::__throw_domain_error("Invalid input");
	else if (std::imag(this->__tau().__val) <= _Real{0})
	  std::__throw_domain_error("__jacobi_lattice_t: "
			"Lattice parameter must have positive imaginary part.");
	else
	  {
	    auto __det = std::real(_M_omega_3) * std::imag(_M_omega_1)
		       - std::imag(_M_omega_3) * std::real(_M_omega_1);
	    if (std::abs(__det) == 0)
	      std::__throw_domain_error("__jacobi_lattice_t: "
			"Lattice frequencies must be linearly independent.");
	  }
      }

      /// Construct the lattice from a single complex lattice parameter
      /// or half period ratio.
      explicit __jacobi_lattice_t(const __tau_t& __tau)
      : _M_omega_1(2 * _S_pi),
        _M_omega_3(2 * _S_pi)
      {
	if (__isnan(__tau.__val))
	  std::__throw_domain_error("Invalid input");
	else if (std::imag(__tau.__val) <= _Real{0})
	  std::__throw_domain_error("__jacobi_lattice_t: "
			"Lattice parameter must have positive imaginary part.");
	else
	  {
	    if constexpr (__gnu_cxx::is_complex_v<_Tp_Omega3>)
	      _M_omega_3 *= __tau.__val;
	    else
	      _M_omega_1 *= __tau.__val;
	  }
      }

      /// Construct the lattice from a single scalar elliptic nome.
      explicit __jacobi_lattice_t(_Tp_Nome __q)
      : __jacobi_lattice_t(__tau_t(_Cmplx{0, -1} * std::log(__q) / _S_pi))
      { }

      /// Return the acalar lattice parameter or half period ratio.
      __tau_t
      __tau() const
      { return __tau_t(this->_M_omega_3 / this->_M_omega_1); }

      /// Return the first lattice frequency.
      _Tp_Omega1
      __omega_1() const
      { return this->_M_omega_1; }

      /// Return the second lattice frequency.
      _Cmplx
      __omega_2() const
      { return -(_Cmplx(this->_M_omega_1) + _Cmplx(this->omega_3)); }

      /// Return the third lattice frequency.
      _Tp_Omega3
      __omega_3() const
      { return this->_M_omega_3; }

      _Tp_Nome
      __ellnome() const;

      __arg_t
      __reduce(const _Cmplx& __z) const;

      static constexpr auto _S_pi = __gnu_cxx::__const_pi<_Real>();
      _Tp_Omega1 _M_omega_1;
      _Tp_Omega3 _M_omega_3;
    };

  /**
   * Return the elliptic nome corresponding to the lattice parameter.
   */
  template<typename _Tp_Omega1, typename _Tp_Omega3>
    typename __jacobi_lattice_t<_Tp_Omega1, _Tp_Omega3>::_Tp_Nome
    __jacobi_lattice_t<_Tp_Omega1, _Tp_Omega3>::__ellnome() const
    {
      const auto _S_i = _Cmplx{0, 1};
      const auto _S_pi = __gnu_cxx::__const_pi<_Real>();
      if constexpr (__gnu_cxx::is_complex_v<_Tp_Nome>)
	return std::exp(_S_i * _S_pi * this->__tau().__val);
      else
	return std::real(std::exp(_S_i * _S_pi * this->__tau().__val));
    }

  /**
   * Reduce the argument to the fundamental lattice parallelogram
   * @f$ (0, 2\pi, 2\pi (1 + \tau), 2\pi \tau) @f$.
   * This is sort of like a 2D lattice remquo.
   *
   * @param __z The argument to be reduced.
   * @return A struct containing the argument reduced to the interior
   *         of the fundamental parallelogram and two integers indicating
   *         the number of periods in the 'real' and 'tau' directions.
   */
  template<typename _Tp1, typename _Tp3>
    typename __jacobi_lattice_t<_Tp1, _Tp3>::__arg_t
    __jacobi_lattice_t<_Tp1, _Tp3>::
    __reduce(const typename __jacobi_lattice_t<_Tp1, _Tp3>::_Cmplx& __z) const
    {
      const auto _S_pi = __gnu_cxx::__const_pi<_Real>();

      const auto __tau = this->__tau().__val;
      const auto __tau_r = std::real(__tau);
      const auto __tau_i = std::imag(__tau);
      const auto __z_r = std::real(__z);
      const auto __z_i = std::imag(__z);

      // Solve z = (z_r, z_i) = pi a (1, 0) + pi b (tau_r, tau_i).
      const auto __b = __z_i / __tau_i / _S_pi;
      const int __n = std::floor(__b);
      const auto __nu = __b - __n;
      const auto __a = (__z_r - __b * __tau_r * _S_pi) / _S_pi;
      const int __m = std::floor(__a);
      const auto __mu = __a - __m;

      return {__m, __n,
      	      _S_pi * _Cmplx(__mu + __nu * __tau_r, __nu * __tau_i)};
    }

  /**
   * A struct for the non-zero theta functions and their derivatives
   * at zero argument.
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>>
    struct __jacobi_theta_0_t
    {
      __jacobi_theta_0_t(const __jacobi_lattice_t<_Tp1, _Tp3>& __lattice);

      using _Type = typename __jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;

      _Type th1p;
      _Type th1ppp;
      _Type th2;
      _Type th2pp;
      _Type th3;
      _Type th3pp;
      _Type th4;
      _Type th4pp;

      _Type
      dedekind_eta() const
      { return std::cbrt(th2 * th3 * th4 / _Type{2}); }
    };

  /**
   * Return a struct of the Jacobi theta functions and up to three non-zero
   * derivatives evaluated at zero argument.
   */
  template<typename _Tp1, typename _Tp3>
    __jacobi_theta_0_t<_Tp1, _Tp3>::
    __jacobi_theta_0_t(const __jacobi_lattice_t<_Tp1, _Tp3>& __lattice)
    {
      using _Real = __num_traits_t<_Type>;
      constexpr std::size_t _S_max_iter = 50;
      const auto __q = __lattice.__ellnome();
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__q));

      const auto __fact = _Real{2} * std::pow(__q, _Real{0.25L});
      this->th1p = __fact;
      this->th2 = __fact;
      this->th3 = _Real{1};
      this->th4 = _Real{1};
      this->th1ppp = _Real{0};
      this->th2pp = _Real{0};
      this->th3pp = _Real{0};
      this->th4pp = _Real{0};
      auto __q2n = _Type{1};
      for (std::size_t __n = 1; __n < _S_max_iter; ++__n)
	{
	  __q2n *= __q;
	  const auto __tp = _Real{1} + __q2n;
	  this->th3 *= __tp * __tp;
	  const auto __tm = _Real{1} - __q2n;
	  this->th4 *= __tm * __tm;

	  this->th3pp += __q2n / __tp / __tp;
	  this->th4pp += __q2n / __tm / __tm;

	  __q2n *= __q;
	  const auto __tm2 = _Real{1} - __q2n;
	  this->th3 *= __tm2;
	  this->th4 *= __tm2;
	  this->th2 *= __tm2;
	  this->th1p *= __tm2 * __tm2 * __tm2;
	  const auto __tp2 = _Real{1} + __q2n;
	  this->th2 *= __tp2 * __tp2;

	  this->th1ppp += __q2n / __tm2 / __tm2;
	  this->th2pp += __q2n / __tp2 / __tp2;

	  if (std::abs(__q2n) < _S_eps)
	    break;
	}
      // Could check th1p =? th2pp * th3pp * th4pp at this point.
      // Could check th1ppp =? th2pp + th3pp + th4pp at this point.
      this->th1ppp = (_Real{-1} + _Real{24} * this->th1ppp) * this->th1p;
      this->th2pp = (_Real{-1} - _Real{8} * this->th2pp) * this->th2;
      this->th3pp = _Real{-8} * this->th3;
      this->th4pp = _Real{8} * this->th4;
    }

  /**
   * A struct of the Weierstrass elliptic function roots.
   * @f[
   *    e_1 = \frac{\pi^2}{12\omega_1^2}(\theta_2^4(q,0) + 2\theta_4^4(q,0))
   * @f]
   * @f[
   *    e_2 = \frac{\pi^2}{12\omega_1^2}(\theta_2^4(q,0) - \theta_4^4(q,0))
   * @f]
   * @f[
   *    e_3 = \frac{\pi^2}{12\omega_1^2}(-2\theta_2^4(q,0) - \theta_4^4(q,0))
   * @f]
   * Note that @f$ e_1 + e_2 + e_3 = 0 @f$
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>>
    struct __weierstrass_roots_t
    {
      using _Type = typename __jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;
      using _Real = __num_traits_t<_Type>;
      using _Cmplx = std::complex<_Real>;

      _Type __e1, __e2, __e3;

      __weierstrass_roots_t(const __jacobi_lattice_t<_Tp1, _Tp3>& __lattice);
    };

  /**
   * Constructor for the Weierstrass roots.
   *
   * @param __lattice The Jacobi latticce.
   */
  template<typename _Tp1, typename _Tp3>
    __weierstrass_roots_t<_Tp1, _Tp3>::
    __weierstrass_roots_t(const __jacobi_lattice_t<_Tp1, _Tp3>& __lattice)
    {
      const auto _S_pi = __gnu_cxx::__const_pi<_Real>();

      const auto __tht0 = __jacobi_theta_0_t(__lattice);

      const auto __th22 = __tht0.th2 * __tht0.th2;
      const auto __th24 = __th22 * __th22;
      const auto __th42 = __tht0.th4 * __tht0.th4;
      const auto __th44 = __th42 * __th42;
      const auto __fr = _S_pi / __lattice.__omega_1();
      const auto __fc = __fr * __fr / _Real{12};

      __e1 = __fc * (__th24 + _Real{2} * __th44);
      __e2 = __fc * (__th24 - __th44);
      __e3 = __fc * (_Real{-2} * __th24 - __th44);
    }

  /**
   * A struct of the Weierstrass elliptic function invariants.
   * @f[
   *    g_2 = 2(e_1 e_2 + e_2 e_3 + e_3 e_1)
   * @f]
   * @f[
   *    g_3 = 4(e_1 e_2 e_3)
   * @f]
   */
  template<typename _Tp1, typename _Tp3>
    struct __weierstrass_invariants_t
    {
      using _Type = typename __jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;
      using _Real = __num_traits_t<_Type>;
      using _Cmplx = std::complex<_Real>;

      _Type __g_2, __g_3;

      __weierstrass_invariants_t(const __jacobi_lattice_t<_Tp1, _Tp3>&);

      _Type
      __klein_j() const
      {
	const auto __g_2p3 = __g_2 * __g_2 * __g_2;
	return _Type{1738} * __g_2p3 / (__g_2p3 - _Type{27} * __g_3 * __g_3);
      }
    };

  /**
   * Constructor for the Weierstrass invariants.
   * @f[
   *    g_2 = 2(e_1 e_2 + e_2 e_3 + e_3 e_1)
   * @f]
   * @f[
   *    g_3 = 4(e_1 e_2 e_3)
   * @f]
   */
  template<typename _Tp1, typename _Tp3>
    __weierstrass_invariants_t<_Tp1, _Tp3>::
    __weierstrass_invariants_t(const __jacobi_lattice_t<_Tp1, _Tp3>& __lattice)
    {
      const auto __roots = __weierstrass_roots_t<_Tp1, _Tp3>(__lattice);
      __g_2 = _Real{2} * (__roots.__e1 * __roots.__e1
        		+ __roots.__e2 * __roots.__e2
        		+ __roots.__e3 * __roots.__e3);
      __g_3 = _Real{4} * __roots.__e1 * __roots.__e2 * __roots.__e3;
    }

  /**
   * Return the Jacobi @f$ \theta_1 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta-1 function is defined by
   * @f[
   *  \theta_1(q,x) = 2\sum_{n=1}^{\infty}(-1)^n
   *                   q^{(n+\frac{1}{2})^2}\sin{(2n+1)x}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_1_sum(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;

      _Tp __sum{};
      _Real __sign{-1};
      for (std::size_t __n = 0; __n < _S_max_iter; ++__n)
	{
	  __sign *= -1;
	  const auto __term = __sign
			    * std::pow(__q, _Real((__n + 0.5L) * (__n + 0.5L)))
			    * std::sin(_Real(2 * __n + 1) * __x);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Real{2} * __sum;
    }

  /**
   * Return the Jacobi @f$ \theta_1 @f$ function by accumulation of the product.
   *
   * The Jacobi or elliptic theta-1 function is defined by
   * @f[
   *  \theta_1(q,x) = 2 q^{1/4} \sin(x) \prod_{n=1}^{\infty}
   *                   (1 - q^{2n})(1 - 2q^{2n}\cos(2x) + q^{4n})
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_1_prod(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;
      const auto __q2 = __q * __q;
      const auto __q4 = __q2 * __q2;
      const auto __cos2x = std::cos(_Real{2} * __x);

      auto __q2n = _Tp{1};
      auto __q4n = _Tp{1};
      auto __prod = _Tp{1};
      for (std::size_t __n = 1; __n < _S_max_iter; ++__n)
	{
	  __q2n *= __q2;
	  __q4n *= __q4;
	  const auto __fact = (_Real{1} - __q2n)
			    * (_Real{1} - _Real{2} * __q2n * __cos2x + __q4n);
	  __prod *= __fact;
	  if (std::abs(__fact) < _S_eps)
	    break;
	}

      return _Real{2} * std::pow(__q, _Tp{0.25L}) * std::sin(__x) * __prod;
    }

  /**
   * Return the Jacobi @f$ \theta_1 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_1(q,x) = 2\sum_{n=1}^{\infty}(-1)^n
   *                   q^{(n+\frac{1}{2})^2}\sin{(2n+1)x}
   * @f]
   *
   * Regarding the nome and the theta function as functions of the lattice
   * parameter @f$ \tau -i log(q)/ \pi @f$ or @f$ q = e^{i\pi\tau} @f$
   * the lattice parameter is transformed to maximize its imaginary part:
   * @f[
   *   \theta_1(\tau+1,x) = -i e^{i\pi/4}\theta_1(\tau,x)
   * @f]
   * and
   * @f[
   *   \sqrt{-i\tau}\theta_1(\tau,x) = e^{(i\tau x^2/\pi)}\theta_1(\tau',\tau' x)
   * @f]
   * where the new lattice parameter is @f$ \tau' = -1/\tau @f$.
   *
   * The argument is reduced with
   * @f[
   *   \theta_1(q, x+(m+n\tau)\pi) = (-1)^{m+n}q^{-n^2}e^{-2inx}\theta_1(q, x)
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __jacobi_theta_1(std::complex<_Tp> __q, std::complex<_Tp> __x)
    {
      using _Real = __num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      const auto _S_i = std::complex<_Real>{0, 1};
      constexpr auto _S_q_min = _Real{0.001L};
      constexpr auto _S_q_max = _Real{0.95e-1L};

      if (__isnan(__q) || __isnan(__x))
	return _Tp{_S_NaN};
      else if (std::abs(__q) >= _Real{1})
	std::__throw_domain_error(__N("__jacobi_theta_1:"
				      " nome q out of range"));
      else if (std::abs(__q) < _S_q_min || std::abs(__q) > _S_q_max)
	return __jacobi_theta_1_prod(__q, __x);
      else if (std::abs(__x) < _S_eps)
	return std::complex<_Tp>{0, 0};
      else
	{
	  const auto __lattice = __jacobi_lattice_t<_Cmplx, _Cmplx>(__q);
	  auto __tau = __lattice.__tau().__val;

	  const auto __x_red = __lattice.__reduce(__x);
	  auto __fact = std::complex<_Tp>{1, 0};
	  if (__x_red.__m != 0)
	    __fact *= __gnu_cxx::__parity<_Tp>(__x_red.__m);
	  if (__x_red.__n != 0)
	    __fact *= __gnu_cxx::__parity<_Tp>(__x_red.__n)
	    	    * std::exp(_S_i * _Real{-2 * __x_red.__n} * __x_red.__z)
		    * std::pow(__q, -__x_red.__n * __x_red.__n);
	  __x = __x_red.__z;

	  // theta_1(tau+1, z) = exp(i tau/4) theta_1(tau, z)
	  const auto __itau = std::floor(std::real(__tau));
	  __tau -= __itau;
	  __fact *= __polar_pi(_Real{1}, __itau / _Real{4});

	  if (std::imag(__tau) < 0.5)
	    {
	      const auto __fact2 = _S_i * std::sqrt(-_S_i * __tau);
	      __tau = _Real{-1} / __tau;
	      const auto __phase = std::exp(_S_i * __tau * __x * __x / _S_pi);
	      __fact *= __phase / __fact2;
	      __q = std::exp(_S_i * _S_pi * __tau);
	      __x *= __tau;
	    }

	  return __fact * __jacobi_theta_1_sum(__q, __x);
	}
    }

  /**
   * Return the Jacobi @f$ \theta_1 @f$ function for real nome and argument.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_1(q,x) = 2\sum_{n=1}^{\infty}(-1)^n
   *                   q^{(n+\frac{1}{2})^2}\sin{(2n+1)x}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_1(_Tp __q, const _Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;

      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto __ret = __jacobi_theta_1(_Cmplx(__q), _Cmplx(__x));

      if (std::abs(__ret) > _S_eps
	  && std::abs(std::imag(__ret)) > _S_eps * std::abs(__ret))
	std::__throw_runtime_error("__jacobi_theta_1: "
				 "Unexpected large imaginary part");
      else
	return std::real(__ret);
    }

  /**
   * Return the Jacobi @f$ \theta_2 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta-2 function is defined by
   * @f[
   *  \theta_2(q,x) = 2\sum_{n=1}^{\infty}
   *                   q^{(n+\frac{1}{2})^2}\cos{(2n+1)x}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_2_sum(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;

      _Tp __sum{};
      for (std::size_t __n = 0; __n < _S_max_iter; ++__n)
	{
	  const auto __term = std::pow(__q, _Real((__n + 0.5L) * (__n + 0.5L)))
			    * std::cos(_Real(2 * __n + 1) * __x);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Real{2} * __sum;
    }

  /**
   * Return the Jacobi @f$ \theta_2 @f$ function by accumulation of the product.
   *
   * The Jacobi or elliptic theta-2 function is defined by
   * @f[
   *  \theta_2(q,x) = 2 q^{1/4} \sin(x) \prod_{n=1}^{\infty}
   *                   (1 - q^{2n})(1 + 2q^{2n}\cos(2x) + q^{4n})
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_2_prod(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;
      const auto __q2 = __q * __q;
      const auto __q4 = __q2 * __q2;
      const auto __cos2x = std::cos(_Real{2} * __x);

      auto __q2n = _Tp{1};
      auto __q4n = _Tp{1};
      auto __prod = _Tp{1};
      for (std::size_t __n = 1; __n < _S_max_iter; ++__n)
	{
	  __q2n *= __q2;
	  __q4n *= __q4;
	  const auto __fact = (_Real{1} - __q2n)
			    * (_Real{1} + _Real{2} * __q2n * __cos2x + __q4n);
	  __prod *= __fact;
	  if (std::abs(__fact) < _S_eps)
	    break;
	}

      return _Real{2} * std::pow(__q, _Tp{0.25L}) * std::cos(__x) * __prod;
    }

  // Pre-declare Jacobi theta_4 sum ...
  template<typename _Tp>
    _Tp
    __jacobi_theta_4_sum(_Tp __q, _Tp __x);

  // ... and product.
  template<typename _Tp>
    _Tp
    __jacobi_theta_4_prod(_Tp __q, _Tp __x);

  /**
   * Return the Jacobi @f$ \theta_2 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_2(q,x) = 2\sum_{n=1}^{\infty}
   *                   q^{(n+\frac{1}{2})^2}\cos{(2n+1)x}
   * @f]
   *
   * Regarding the nome and the theta function as functions of the lattice
   * parameter @f$ \tau -i log(q)/ \pi @f$ or @f$ q = e^{i\pi\tau} @f$
   * the lattice parameter is transformed to maximize its imaginary part:
   * @f[
   *   \theta_2(\tau+1,x) = e^{i\pi/4}\theta_2(\tau,x)
   * @f]
   * and
   * @f[
   *   \sqrt{-i\tau}\theta_2(\tau,x) = e^{(i\tau x^2/\pi)}\theta_4(\tau',\tau' x)
   * @f]
   * where the new lattice parameter is @f$ \tau' = -1/\tau @f$.
   *
   * The argument is reduced with
   * @f[
   *  \theta_2(q, x + (m+n\tau)\pi) = (-1)^{m}q^{-n^2}e^{-2inx}\theta_2(q, x)
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __jacobi_theta_2(std::complex<_Tp> __q, std::complex<_Tp> __x)
    {
      using _Real = __num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      const auto _S_i = std::complex<_Real>{0, 1};
      constexpr auto _S_q_min = _Real{0.001L};
      constexpr auto _S_q_max = _Real{0.95e-1L};

      if (__isnan(__q) || __isnan(__x))
	return _Tp{_S_NaN};
      else if (std::abs(__q) >= _Real{1})
	std::__throw_domain_error(__N("__jacobi_theta_2:"
				      " nome q out of range"));
      else if (std::abs(__q) < _S_q_min || std::abs(__q) > _S_q_max)
	return __jacobi_theta_2_prod(__q, __x);
      else if (std::abs(__x) < _S_eps)
	return __jacobi_theta_0_t(__jacobi_lattice_t<_Cmplx, _Cmplx>(__q)).th2;
      else
	{
	  const auto __lattice = __jacobi_lattice_t<_Cmplx, _Cmplx>(__q);
	  auto __tau = __lattice.__tau().__val;

	  const auto __x_red = __lattice.__reduce(__x);
	  auto __fact = std::complex<_Tp>{1, 0};
	  if (__x_red.__m != 0)
	    __fact *= __gnu_cxx::__parity<_Tp>(__x_red.__m);
	  if (__x_red.__n != 0)
	    __fact *= std::exp(_S_i * _Real{-2 * __x_red.__n} * __x_red.__z)
		    * std::pow(__q, -__x_red.__n * __x_red.__n);
	  __x = __x_red.__z;

	  // theta_2(tau+1, z) = theta_2(tau, z)
	  const auto __itau = std::floor(std::real(__tau));
	  __tau -= __itau;
	  __fact *= __polar_pi(_Real{1}, __itau / _Real{4});

	  if (std::imag(__tau) < 0.5)
	    {
	      const auto __fact2 = std::sqrt(-_S_i * __tau);
	      __tau = _Real{-1} / __tau;
	      const auto __phase = std::exp(_S_i * __tau * __x * __x / _S_pi);
	      __fact *= __phase / __fact2;
	      __q = std::exp(_S_i * _S_pi * __tau);
	      __x *= __tau;
	      return __fact * __jacobi_theta_4_sum(__q, __x);
	    }
	  else
	    return __fact * __jacobi_theta_2_sum(__q, __x);
	}
    }

  /**
   * Return the Jacobi @f$ \theta_2 @f$ function for real nome and argument.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_2(q,x) = 2\sum_{n=1}^{\infty}
   *                   q^{(n+\frac{1}{2})^2}\cos{(2n+1)x}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_2(_Tp __q, const _Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));

      const auto __ret = __jacobi_theta_2(_Cmplx(__q), _Cmplx(__x));

      if (std::abs(__ret) > _S_eps
	  && std::abs(std::imag(__ret)) > _S_eps * std::abs(__ret))
	std::__throw_runtime_error("__jacobi_theta_2: "
				 "Unexpected large imaginary part");
      else
	return std::real(__ret);
    }

  /**
   * Return the Jacobi @f$ \theta_3 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta-3 function is defined by
   * @f[
   *  \theta_3(q,x) = 1 + 2\sum_{n=1}^{\infty} q^{n^2}\cos{2nx}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_3_sum(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;

      _Tp __sum{};
      for (std::size_t __n = 1; __n < _S_max_iter; ++__n)
	{
	  const auto __term = std::pow(__q, _Real(__n * __n))
			    * std::cos(_Real(2 * __n) * __x);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Real{1} + _Real{2} * __sum;
    }

  /**
   * Return the Jacobi @f$ \theta_3 @f$ function by accumulation of the product.
   *
   * The Jacobi or elliptic theta-3 function is defined by
   * @f[
   *  \theta_3(q,x) = \prod_{n=1}^{\infty}
   *                   (1 - q^{2n})(1 + 2q^{2n-1}\cos(2x) + q^{4n-2})
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_3_prod(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;
      const auto __q2 = __q * __q;
      const auto __q4 = __q2 * __q2;
      const auto __cos2x = std::cos(_Real{2} * __x);

      auto __q2nm1 = __q;
      auto __q4nm2 = __q2;
      auto __prod = _Tp{1};
      for (std::size_t __n = 1; __n < _S_max_iter; ++__n)
	{
	  const auto __fact = (_Real{1} - __q2nm1 * __q)
			    * (_Real{1} + _Real{2} * __q2nm1 * __cos2x + __q4nm2);
	  __prod *= __fact;
	  if (std::abs(__fact) < _S_eps)
	    break;
	  __q2nm1 *= __q2;
	  __q4nm2 *= __q4;
	}

      return __prod;
    }

  /**
   * Return the Jacobi @f$ \theta_3 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_3(q,x) = 1 + 2\sum_{n=1}^{\infty} q^{n^2}\cos{2nx}
   * @f]
   *
   * Regarding the nome and the theta function as functions of the lattice
   * parameter @f$ \tau -i log(q)/ \pi @f$ or @f$ q = e^{i\pi\tau} @f$
   * the lattice parameter is transformed to maximize its imaginary part:
   * @f[
   *   \theta_3(\tau+1,x) = \theta_3(\tau,x)
   * @f]
   * and
   * @f[
   *   \sqrt{-i\tau}\theta_3(\tau,x) = e^{(i\tau x^2/\pi)}\theta_3(\tau',\tau' x)
   * @f]
   * where the new lattice parameter is @f$ \tau' = -1/\tau @f$.
   *
   * The argument is reduced with
   * @f[
   *   \theta_3(q, x + (m+n\tau)\pi) = q^{-n^2} e^{-2inx} \theta_3(q, x)
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __jacobi_theta_3(std::complex<_Tp> __q, std::complex<_Tp> __x)
    {
      using _Real = __num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      const auto _S_i = std::complex<_Real>{0, 1};
      constexpr auto _S_q_min = _Real{0.001L};
      constexpr auto _S_q_max = _Real{0.95e-1L};

      if (__isnan(__q) || __isnan(__x))
	return _Tp{_S_NaN};
      else if (std::abs(__q) >= _Real{1})
	std::__throw_domain_error(__N("__jacobi_theta_3:"
				      " nome q out of range"));
      else if (std::abs(__q) < _S_q_min || std::abs(__q) > _S_q_max)
	return __jacobi_theta_3_prod(__q, __x);
      else if (std::abs(__x) < _S_eps)
	return __jacobi_theta_0_t(__jacobi_lattice_t<_Cmplx, _Cmplx>(__q)).th3;
      else
	{
	  const auto __lattice = __jacobi_lattice_t<_Cmplx, _Cmplx>(__q);
	  auto __tau = __lattice.__tau().__val;

	  const auto __x_red = __lattice.__reduce(__x);
	  auto __fact = std::complex<_Tp>{1, 0};
	  if (__x_red.__n != 0)
	    __fact *= std::exp(_S_i * _Real{-2 * __x_red.__n} * __x_red.__z)
		    * std::pow(__q, -__x_red.__n * __x_red.__n);
	  __x = __x_red.__z;

	  // theta_3(tau+1, z) = theta_3(tau, z)
	  const auto __itau = std::floor(std::real(__tau));
	  __tau -= __itau;

	  if (std::imag(__tau) < 0.5)
	    {
	      const auto __fact2 = std::sqrt(-_S_i * __tau);
	      __tau = _Real{-1} / __tau;
	      const auto __phase = std::exp(_S_i * __tau * __x * __x / _S_pi);
	      __fact *= __phase / __fact2;
	      __q = std::exp(_S_i * _S_pi * __tau);
	      __x *= __tau;
	    }

	  return __fact * __jacobi_theta_3_sum(__q, __x);
	}
    }

  /**
   * Return the Jacobi @f$ \theta_3 @f$ function for real nome and argument.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_3(q,x) = 1 + 2\sum_{n=1}^{\infty} q^{n^2}\cos{2nx}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_3(_Tp __q, const _Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));

      const auto __ret = __jacobi_theta_3(_Cmplx(__q), _Cmplx(__x));

      if (std::abs(__ret) > _S_eps
	  && std::abs(std::imag(__ret)) > _S_eps * std::abs(__ret))
	std::__throw_runtime_error("__jacobi_theta_3: "
				 "Unexpected large imaginary part");
      else
	return std::real(__ret);
    }

  /**
   * Return the Jacobi @f$ \theta_4 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_4(q,x) = 1 + 2\sum_{n=1}^{\infty}(-1)^n q^{n^2}\cos{2nx}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_4_sum(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;

      _Tp __sum{};
      _Real __sign{1};
      for (std::size_t __n = 1; __n < _S_max_iter; ++__n)
	{
	  __sign *= -1;
	  const auto __term = __sign * std::pow(__q, _Real(__n * __n))
			    * std::cos(_Real(2 * __n) * __x);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Real{1} + _Real{2} * __sum;
    }

  /**
   * Return the Jacobi @f$ \theta_4 @f$ function by accumulation of the product.
   *
   * The Jacobi or elliptic theta-4 function is defined by
   * @f[
   *  \theta_4(q,x) = \prod_{n=1}^{\infty}
   *                   (1 - q^{2n})(1 - 2q^{2n-1}\cos(2x) + q^{4n-2})
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_4_prod(_Tp __q, _Tp __x)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      constexpr std::size_t _S_max_iter = 50;
      const auto __q2 = __q * __q;
      const auto __q4 = __q2 * __q2;
      const auto __cos2x = std::cos(_Real{2} * __x);

      auto __q2nm1 = __q;
      auto __q4nm2 = __q2;
      auto __prod = _Tp{1};
      for (std::size_t __n = 1; __n < _S_max_iter; ++__n)
	{
	  const auto __fact = (_Real{1} - __q2nm1 * __q)
			    * (_Real{1} - _Real{2} * __q2nm1 * __cos2x + __q4nm2);
	  __prod *= __fact;
	  if (std::abs(__fact) < _S_eps)
	    break;
	  __q2nm1 *= __q2;
	  __q4nm2 *= __q4;
	}

      return __prod;
    }

  /**
   * Return the Jacobi @f$ \theta_4 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta-4 function is defined by
   * @f[
   *  \theta_4(q,x) = 1 + 2\sum_{n=1}^{\infty}(-1)^n q^{n^2}\cos{2nx}
   * @f]
   *
   * Regarding the nome and the theta function as functions of the lattice
   * parameter @f$ \tau -i log(q)/ \pi @f$ or @f$ q = e^{i\pi\tau} @f$
   * the lattice parameter is transformed to maximize its imaginary part:
   * @f[
   *   \theta_4(\tau+1,x) = \theta_4(\tau,x)
   * @f]
   * and
   * @f[
   *   \sqrt{-i\tau}\theta_4(\tau,x) = e^{(i\tau x^2/\pi)}\theta_2(\tau',\tau' x)
   * @f]
   * where the new lattice parameter is @f$ \tau' = -1/\tau @f$.
   *
   * The argument is reduced with
   * @f[
   *   \theta_4(q, z+(m + n\tau)\pi) = (-1)^n q^{-n^2}e^{-2inz}\theta_4(q, z)
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __jacobi_theta_4(std::complex<_Tp> __q, std::complex<_Tp> __x)
    {
      using _Real = __num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__x));
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));
      const auto _S_pi = __gnu_cxx::__const_pi(std::abs(__x));
      const auto _S_i = std::complex<_Real>{0, 1};
      constexpr auto _S_q_min = _Real{0.001L};
      constexpr auto _S_q_max = _Real{0.95e-1L};

      if (__isnan(__q) || __isnan(__x))
	return _Tp{_S_NaN};
      else if (std::abs(__q) >= _Real{1})
	std::__throw_domain_error(__N("__jacobi_theta_4:"
				      " nome q out of range"));
      else if (std::abs(__q) < _S_q_min || std::abs(__q) > _S_q_max)
	return __jacobi_theta_4_prod(__q, __x);
      else if (std::abs(__x) < _S_eps)
	return __jacobi_theta_0_t(__jacobi_lattice_t<_Cmplx, _Cmplx>(__q)).th4;
      else
	{
	  const auto __lattice = __jacobi_lattice_t<_Cmplx, _Cmplx>(__q);
	  auto __tau = __lattice.__tau().__val;

	  const auto __x_red = __lattice.__reduce(__x);
	  auto __fact = std::complex<_Tp>{1, 0};
	  if (__x_red.__n != 0)
	    __fact *= std::exp(_S_i * _Real{-2 * __x_red.__n} * __x_red.__z)
		    * std::pow(__q, -__x_red.__n * __x_red.__n);
	  if (__x_red.__n != 0)
	    __fact *= __gnu_cxx::__parity<_Tp>(__x_red.__n);
	  __x = __x_red.__z;

	  // theta_4(tau+1, z) = theta_4(tau, z)
	  const auto __itau = std::floor(std::real(__tau));
	  __tau -= __itau;

	  if (std::imag(__tau) < 0.5)
	    {
	      const auto __fact2 = std::sqrt(-_S_i * __tau);
	      __tau = _Real{-1} / __tau;
	      const auto __phase = std::exp(_S_i * __tau * __x * __x / _S_pi);
	      __fact *= __phase / __fact2;
	      __q = std::exp(_S_i * _S_pi * __tau);
	      __x *= __tau;
	      return __fact * __jacobi_theta_2_sum(__q, __x);
	    }
	  else
	    return __fact * __jacobi_theta_4_sum(__q, __x);
	}
    }

  /**
   * Return the Jacobi @f$ \theta_4 @f$ function for real nome and argument.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_4(q,x) = 1 + 2\sum_{n=1}^{\infty}(-1)^n q^{n^2}\cos{2nx}
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   * @param __x The argument.
   */
  template<typename _Tp>
    _Tp
    __jacobi_theta_4(_Tp __q, const _Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__x));

      const auto __ret = __jacobi_theta_4(_Cmplx(__q), _Cmplx(__x));

      if (std::abs(__ret) > _S_eps
	  && std::abs(std::imag(__ret)) > _S_eps * std::abs(__ret))
	std::__throw_runtime_error("__jacobi_theta_4: "
				 "Unexpected large imaginary part");
      else
	return std::real(__ret);
    }

  /**
   * Return a structure containing the three primary Jacobi elliptic functions:
   * @f$ sn(k, u), cn(k, u), dn(k, u) @f$.
   *
   * @param __k The elliptic modulus @f$ |k| < 1 @f$.
   * @param __u The argument.
   * @return An object containing the three principal Jacobi elliptic functions,
   *         @f$ sn(k, u), cn(k, u), dn(k, u) @f$ and the means to compute
   *         the remaining nine as well as the amplitude.
   */
  template<typename _Tp>
    __gnu_cxx::__jacobi_ellint_t<_Tp>
    __jacobi_ellint(_Tp __k, _Tp __u)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__u));
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::abs(__u));

      if (__isnan(__k) || __isnan(__u))
	return __gnu_cxx::__jacobi_ellint_t<_Tp>{_S_NaN, _S_NaN, _S_NaN};
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__jacobi_ellint:"
				      " argument k out of range"));
      else if (std::abs(_Tp{1} - __k) < _Tp{2} * _S_eps)
	{
	  auto __sn = std::tanh(__u);
	  auto __cn = _Tp{1} / std::cosh(__u);
	  auto __dn = __cn;
	  return __gnu_cxx::__jacobi_ellint_t<_Tp>{__sn, __cn, __dn};
	}
      else if (std::abs(__k) < _Tp{2} * _S_eps)
	{
	  auto __sn = std::sin(__u);
	  auto __cn = std::cos(__u);
	  auto __dn = _Tp{1};
	  return __gnu_cxx::__jacobi_ellint_t<_Tp>{__sn, __cn, __dn};
	}
      else
	{
	  const auto _S_CA = std::sqrt(_S_eps);
	  const auto _S_N = 100;
	  std::vector<_Tp> __m;
	  std::vector<_Tp> __n;
	  __m.reserve(20);
	  __n.reserve(20);
	  _Tp __c, __d;
	  auto __mc = _Tp{1} - __k * __k;
	  bool __bo = (__mc < _Tp{0});
	  if (__bo)
	    {
	      __d = __k * __k;
	      __mc /= -_Tp{1} / __d;
	      __u *= (__d = __k);
	    }
	  auto __a = _Tp{1};
	  auto __dn = _Tp{1};
	  auto __l = _S_N;
	  for (auto __i = 0; __i < _S_N; ++__i)
	    {
	      __l = __i;
	      __m.push_back(__a);
	      __n.push_back(__mc = std::sqrt(__mc));
	      __c = 0.5 * (__a + __mc);
	      if (std::abs(__a - __mc) <= _S_CA * __a)
		break;
	      __mc *= __a;
	      __a = __c;
	    }
	  __u *= __c;
	  auto __sn = std::sin(__u);
	  auto __cn = std::cos(__u);
	  if (__sn != _Tp{0})
	    {
	      __a = __cn / __sn;
	      __c *= __a;
	      for (auto __ii = __l; __ii >= 0; --__ii)
		{
		  _Tp __b = __m[__ii];
		  __a *= __c;
		  __c *= __dn;
		  __dn = (__n[__ii] + __a) / (__b + __a);
		  __a = __c / __b;
		}
	      __a = _Tp{1} / std::hypot(_Tp{1}, __c);
	      __sn = std::copysign(__a, __sn);
	      __cn = __c * __sn;
	    }
	  if (__bo)
	    {
	      std::swap(__dn, __cn);
	      __sn /= __d;
	    }
	  return __gnu_cxx::__jacobi_ellint_t<_Tp>{__sn, __cn, __dn};
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_THETA_TCC
