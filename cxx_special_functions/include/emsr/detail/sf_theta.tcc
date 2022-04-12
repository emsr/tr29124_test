
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

/** @file bits/sf_theta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_THETA_TCC
#define SF_THETA_TCC 1

#include <stdexcept>
#include <vector>
#include <tuple>
#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/sf_ellint.h>
#include <emsr/sf_trig.h> // reperiodized polar

namespace emsr
{
namespace detail
{

  /**
   * Compute and return the exponential @f$ \theta_2 @f$ function
   * by series expansion:
   * @f[
   *    \theta_2(\nu, x) = \frac{1}{\sqrt{\pi x}}
   *                       \sum_{k=-\infty}^{\infty}(-1)^k e^{-(\nu+k)^2/x}
   * @f]
   */
  template<typename Tp>
    Tp
    theta_2_sum(Tp nu, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      auto sum = std::exp(-nu * nu / x);
      auto sign = Tp{-1};
      for (auto k = 1; k < 20; ++k)
	{
	  const auto nup = nu + Tp(k);
	  const auto termp = sign * std::exp(-nup * nup / x);
	  const auto num = nu - Tp(k);
	  const auto termm = sign * std::exp(-num * num / x);
	  sum += termp + termm;
	  sign = -sign;
	  if (std::abs(termp) < s_eps * std::abs(sum)
	   && std::abs(termm) < s_eps * std::abs(sum))
	    break;
	}
      return sum / std::sqrt(s_pi * x);
    }

  /**
   * Compute and return the exponential @f$ \theta_3 @f$ function
   * by series expansion:
   * @f[
   *    \theta_3(\nu, x) = \frac{1}{\sqrt{\pi x}}
   *                       \sum_{k=-\infty}^{\infty} e^{-(\nu+k)^2/x}
   * @f]
   */
  template<typename Tp>
    Tp
    theta_3_sum(Tp nu, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      auto sum = std::exp(-nu * nu / x);
      for (auto k = 1; k < 20; ++k)
	{
	  const auto nup = nu + Tp(k);
	  const auto termp = std::exp(-nup * nup / x);
	  const auto num = nu - Tp(k);
	  const auto termm = std::exp(-num * num / x);
	  sum += termp + termm;
	  if (std::abs(termp) < s_eps * std::abs(sum)
	   && std::abs(termm) < s_eps * std::abs(sum))
	    break;
	}
      return sum / std::sqrt(s_pi * x);
    }

  /**
   * Compute and return the exponential @f$ \theta_2 @f$ function
   * by asymptotic series expansion:
   * @f[
   *    \theta_2(\nu, x) = 2\sum_{k=0}^{\infty} e^{-((k+1/2)\pi)^2 x}
   *                        \cos((2k+1)\nu\pi)
   * @f]
   */
  template<typename Tp>
    Tp
    theta_2_asymp(Tp nu, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      auto sum = Tp{0};
      for (auto k = 0; k < 20; ++k)
	{
	  const auto thing = Tp(2 * k + 1) * s_pi;
	  const auto cosarg = nu * thing;
	  const auto exparg = thing * thing * x / Tp{4};
	  const auto term = std::exp(-exparg) * std::cos(cosarg);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return Tp{2} * sum;
    }

  /**
   * Compute and return the exponential @f$ \theta_3 @f$ function
   * by asymptotic series expansion:
   * @f[
   *    \theta_3(\nu, x) = 1 + 2\sum_{k=1}^{\infty} e^{-(k\pi)^2 x}
   *                           \cos(2k\nu\pi)
   * @f]
   */
  template<typename Tp>
    Tp
    theta_3_asymp(Tp nu, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      auto sum = Tp{0};
      for (auto k = 1; k < 20; ++k)
	{
	  const auto thing = Tp(2 * k) * s_pi;
	  const auto cosarg = nu * thing;
	  const auto exparg = thing * thing * x / Tp{4};
	  const auto term = std::exp(-exparg) * std::cos(cosarg);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return Tp{1} + Tp{2} * sum;
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
   * @param nu The periodic (period = 2) argument
   * @param x The argument
   */
  template<typename Tp>
    Tp
    theta_2(Tp nu, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;

      if (std::isnan(nu) || std::isnan(x))
	return s_NaN;
      else if (std::abs(x) <= Real{1} / s_pi)
	return theta_2_sum(nu, x);
      else
	return theta_2_asymp(nu, x);
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
   * @param nu The periodic (period = 2) argument
   * @param x The argument
   */
  template<typename Tp>
    Tp
    theta_1(Tp nu, Tp x)
    {
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));

      if (std::isnan(nu) || std::isnan(x))
	return s_NaN;
      else if (emsr::fp_is_zero(x))
	return Tp{0};
      else
	return theta_2(nu - Tp{0.5L}, x);
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
   * @param nu The periodic (period = 1) argument
   * @param x The argument
   */
  template<typename Tp>
    Tp
    theta_3(Tp nu, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;

      if (std::isnan(nu) || std::isnan(x))
	return s_NaN;
      else if (std::abs(x) <= Real{1} / s_pi)
	return theta_3_sum(nu, x);
      else
	return theta_3_asymp(nu, x);
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
   * @param nu The periodic (period = 2) argument
   * @param x The argument
   */
  template<typename Tp>
    Tp
    theta_4(Tp nu, Tp x)
    {
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));

      if (std::isnan(nu) || std::isnan(x))
	return s_NaN;
      else
	return theta_3(nu + Tp{0.5L}, x);
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
  template<typename Tp>
    Tp
    ellnome_series(Tp k)
    {
      const auto m = k * k; 
      return m * ((Tp{1} / Tp{16})
	   + m * ((Tp{1} / Tp{32})
	   + m * ((Tp{21} / Tp{1024})
	   + m * ((Tp{31} / Tp{2048})
	   + m * (Tp{6257} / Tp{524288})))));
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
  template<typename Tp>
    Tp
    ellnome_k(Tp k)
    {
      const auto s_pi = Tp{3.1415926535897932384626433832795029L};
      const auto kp = std::sqrt(Tp{1} - k * k);
      const auto K = emsr::comp_ellint_1(k);
      const auto Kp = emsr::comp_ellint_1(kp);
      return std::exp(-s_pi * Kp / K);
    }

  /**
   * Return the elliptic nome given the modulus @c k.
   * @f[
   *    q(k) = exp\left(-\pi\frac{K(k')}{K(k)}\right)
   * @f]
   */
  template<typename Tp>
    Tp
    ellnome(Tp k)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      if (std::isnan(k))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (std::abs(k) > Tp{1})
	throw std::domain_error("ellnome: argument k out of range");
      else if (k < std::pow(Tp{67} * s_eps, Tp{0.125L}))
	return ellnome_series(k);
      else
	return ellnome_k(k);
    }

  /**
   * Return the Neville @f$ \theta_s @f$ function
   * @f[
   *  \theta_s(k,x) = \sqrt{\frac{\pi}{2 k k' K(k)}}
   *                  \theta_1\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   */
  template<typename Tp>
    Tp
    theta_s(Tp k, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};

      if (std::isnan(k) || std::isnan(x))
	return s_NaN;
      else if (std::abs(k) > Tp{1})
	throw std::domain_error("theta_s: argument k out of range");
      else
	{
	  const auto kc = std::sqrt(Tp{1} - k * k);
	  const auto _Kk = emsr::comp_ellint_1(k);
	  const auto q = ellnome(k);
	  return std::sqrt(s_pi_2 / (k * kc * _Kk))
	       * theta_1(q, s_pi_2 * x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_c @f$ function
   * @f[
   *    \theta_c(k,x) = \sqrt{\frac{\pi}{2 k K(k)}}
   *                  \theta_1\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   */
  template<typename Tp>
    Tp
    theta_c(Tp k, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};

      if (std::isnan(k) || std::isnan(x))
	return s_NaN;
      else if (std::abs(k) > Tp{1})
	throw std::domain_error("theta_c: argument k out of range");
      else
	{
	  const auto _Kk = emsr::comp_ellint_1(k);
	  const auto q = ellnome(k);
	  return std::sqrt(s_pi_2 / (k * _Kk))
	       * theta_2(q, s_pi_2 * x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_d @f$ function
   * @f[
   *    \theta_d(k,x) = \sqrt{\frac{\pi}{2K(k)}}
   *                  \theta_3\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   */
  template<typename Tp>
    Tp
    theta_d(Tp k, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};

      if (std::isnan(k) || std::isnan(x))
	return s_NaN;
      else if (std::abs(k) > Tp{1})
	throw std::domain_error("theta_d: argument k out of range");
      else
	{
	  const auto _Kk = emsr::comp_ellint_1(k);
	  const auto q = ellnome(k);
	  return std::sqrt(s_pi_2 / _Kk)
	       * theta_3(q, s_pi_2 * x / _Kk);
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
  template<typename Tp>
    Tp
    theta_n(Tp k, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};

      if (std::isnan(k) || std::isnan(x))
	return s_NaN;
      else if (std::abs(k) > Tp{1})
	throw std::domain_error("theta_n: argument k out of range");
      else
	{
	  const auto kc = std::sqrt(Tp{1} - k * k);
	  const auto _Kk = emsr::comp_ellint_1(k);
	  const auto q = ellnome(k);
	  return std::sqrt(s_pi_2 / (kc * _Kk))
	       * theta_4(q, s_pi_2 * x / _Kk);
	}
    }

  /**
   * A struct representing the Jacobi and Weierstrass lattice.
   * The two types for the frequencies and the subsequent type calculus
   * allow us to treat the rectangulr lattice (real nome, pure imaginary
   * lattice parameter) specially.
   */
  template<typename Tp_Omega1, typename Tp_Omega3 = std::complex<Tp_Omega1>>
    struct jacobi_lattice_t
    {
      static_assert(emsr::is_complex_v<Tp_Omega1>
		 || emsr::is_complex_v<Tp_Omega3>,
		    "One frequecy type must be complex.");
      using _Real_Omega1 = emsr::num_traits_t<Tp_Omega1>;
      using _Real_Omega3 = emsr::num_traits_t<Tp_Omega3>;
      using Real = emsr::fp_promote_t<_Real_Omega1, _Real_Omega3>;
      using Cmplx = std::complex<Real>;
      using _Tp_Nome = std::conditional_t<emsr::is_complex_v<Tp_Omega1>
				       && emsr::is_complex_v<Tp_Omega3>,
					  Cmplx, Real>;

      /**
       * A struct representing a complex scalar lattice parameter
       * or half period ratio.
       */
      struct tau_t
      {
	Cmplx val;

	explicit tau_t(Cmplx tau)
	: val(tau)
	{ }
      };

      /**
       * A struct representing a complex argument reduced
       * to the 'central' lattice cell.
       */
      struct arg_t
      {
	int m;
	int n;
	Cmplx z;
      };

      /// Construct the lattice from two complex lattice frequencies.
      jacobi_lattice_t(const Tp_Omega1& omega1,
			 const Tp_Omega3& omega3)
      : m_omega_1(omega1),
        m_omega_3(omega3)
      {
	if (std::isnan(m_omega_1) || std::isnan(m_omega_3))
	  throw std::domain_error("Invalid input");
	else if (std::imag(this->tau().val) <= Real{0})
	  throw std::domain_error("jacobi_lattice_t: "
			"Lattice parameter must have positive imaginary part.");
	else
	  {
	    auto det = std::real(m_omega_3) * std::imag(m_omega_1)
		       - std::imag(m_omega_3) * std::real(m_omega_1);
	    if (std::abs(det) == 0)
	      throw std::domain_error("jacobi_lattice_t: "
			"Lattice frequencies must be linearly independent.");
	  }
      }

      /// Construct the lattice from a single complex lattice parameter
      /// or half period ratio.
      explicit jacobi_lattice_t(const tau_t& tau)
      : m_omega_1(2 * s_pi),
        m_omega_3(2 * s_pi)
      {
	if (std::isnan(tau.val))
	  throw std::domain_error("Invalid input");
	else if (std::imag(tau.val) <= Real{0})
	  throw std::domain_error("jacobi_lattice_t: Lattice parameter must have positive imaginary part.");
	else
	  {
	    if constexpr (emsr::is_complex_v<Tp_Omega3>)
	      m_omega_3 *= tau.val;
	    else
	      m_omega_1 *= tau.val;
	  }
      }

      /// Construct the lattice from a single scalar elliptic nome.
      explicit jacobi_lattice_t(_Tp_Nome q)
      : jacobi_lattice_t(tau_t(Cmplx{0, -1} * std::log(q) / s_pi))
      {
	if (std::abs(q) == Real{0})
	  throw std::domain_error("jacobi_lattice_t: Nome must be nonzero.");
      }

      /// Return the acalar lattice parameter or half period ratio.
      tau_t
      tau() const
      { return tau_t(this->m_omega_3 / this->m_omega_1); }

      /// Return the first lattice frequency.
      Tp_Omega1
      omega_1() const
      { return this->m_omega_1; }

      /// Return the second lattice frequency.
      Cmplx
      omega_2() const
      { return -(Cmplx(this->m_omega_1) + Cmplx(this->m_omega_3)); }

      /// Return the third lattice frequency.
      Tp_Omega3
      omega_3() const
      { return this->m_omega_3; }

      _Tp_Nome
      ellnome() const;

      arg_t
      reduce(const Cmplx& z) const;

      static constexpr auto s_pi = emsr::pi_v<Real>;
      Tp_Omega1 m_omega_1;
      Tp_Omega3 m_omega_3;
    };

  /**
   * Return the elliptic nome corresponding to the lattice parameter.
   */
  template<typename Tp_Omega1, typename Tp_Omega3>
    typename jacobi_lattice_t<Tp_Omega1, Tp_Omega3>::_Tp_Nome
    jacobi_lattice_t<Tp_Omega1, Tp_Omega3>::ellnome() const
    {
      const auto s_i = Cmplx{0, 1};
      const auto s_pi = emsr::pi_v<Real>;
      if constexpr (emsr::is_complex_v<_Tp_Nome>)
	return std::exp(s_i * s_pi * this->tau().val);
      else
	return std::real(std::exp(s_i * s_pi * this->tau().val));
    }

  /**
   * Reduce the argument to the fundamental lattice parallelogram
   * @f$ (0, 2\pi, 2\pi (1 + \tau), 2\pi \tau) @f$.
   * This is sort of like a 2D lattice remquo.
   *
   * @param z The argument to be reduced.
   * @return A struct containing the argument reduced to the interior
   *         of the fundamental parallelogram and two integers indicating
   *         the number of periods in the 'real' and 'tau' directions.
   */
  template<typename _Tp1, typename _Tp3>
    typename jacobi_lattice_t<_Tp1, _Tp3>::arg_t
    jacobi_lattice_t<_Tp1, _Tp3>::
    reduce(const typename jacobi_lattice_t<_Tp1, _Tp3>::Cmplx& z) const
    {
      const auto s_pi = emsr::pi_v<Real>;

      const auto tau = this->tau().val;
      const auto tau_r = std::real(tau);
      const auto tau_i = std::imag(tau);
      const auto z_r = std::real(z);
      const auto z_i = std::imag(z);

      // Solve z = (z_r, z_i) = pi a (1, 0) + pi b (tau_r, tau_i).
      const auto b = z_i / tau_i / s_pi;
      const int n = std::floor(b);
      const auto nu = b - n;
      const auto a = (z_r - b * tau_r * s_pi) / s_pi;
      const int m = std::floor(a);
      const auto mu = a - m;

      return {m, n,
      	      s_pi * Cmplx(mu + nu * tau_r, nu * tau_i)};
    }

  /**
   * A struct for the non-zero theta functions and their derivatives
   * at zero argument.
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>>
    struct jacobi_theta_0_t
    {
      jacobi_theta_0_t(const jacobi_lattice_t<_Tp1, _Tp3>& lattice);

      using _Type = typename jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;
      using Real = emsr::num_traits_t<_Type>;
      using Cmplx = std::complex<Real>;

      _Type th1p;
      _Type th1ppp;
      _Type th2;
      _Type th2pp;
      _Type th3;
      _Type th3pp;
      _Type th4;
      _Type th4pp;
      _Type eta_1;
      Cmplx eta_2;
      Cmplx eta_3;

      _Type
      dedekind_eta() const
      { return std::cbrt(th2 * th3 * th4 / _Type{2}); }
    };

  /**
   * Return a struct of the Jacobi theta functions and up to three non-zero
   * derivatives evaluated at zero argument.
   */
  template<typename _Tp1, typename _Tp3>
    jacobi_theta_0_t<_Tp1, _Tp3>::
    jacobi_theta_0_t(const jacobi_lattice_t<_Tp1, _Tp3>& lattice)
    {
      constexpr std::size_t s_max_iter = 50;
      const auto q = lattice.ellnome();
      const auto s_eps = emsr::epsilon(std::abs(q));

      const auto fact = Real{2} * std::pow(q, Real{0.25L});
      this->th1p = fact;
      this->th2 = fact;
      this->th3 = Real{1};
      this->th4 = Real{1};
      this->th1ppp = Real{0};
      this->th2pp = Real{0};
      this->th3pp = Real{0};
      this->th4pp = Real{0};
      auto q2n = _Type{1};
      for (std::size_t n = 1; n < s_max_iter; ++n)
	{
	  q2n *= q;
	  const auto tp = Real{1} + q2n;
	  this->th3 *= tp * tp;
	  const auto tm = Real{1} - q2n;
	  this->th4 *= tm * tm;

	  this->th3pp += q2n / tp / tp;
	  this->th4pp += q2n / tm / tm;

	  q2n *= q;
	  const auto tm2 = Real{1} - q2n;
	  this->th3 *= tm2;
	  this->th4 *= tm2;
	  this->th2 *= tm2;
	  this->th1p *= tm2 * tm2 * tm2;
	  const auto tp2 = Real{1} + q2n;
	  this->th2 *= tp2 * tp2;

	  this->th1ppp += q2n / tm2 / tm2;
	  this->th2pp += q2n / tp2 / tp2;

	  if (std::abs(q2n) < s_eps)
	    break;
	}
      // Could check th1p =? th2pp * th3pp * th4pp at this point.
      // Could check th1ppp =? th2pp + th3pp + th4pp at this point.
      this->th1ppp = (Real{-1} + Real{24} * this->th1ppp) * this->th1p;
      this->th2pp = (Real{-1} - Real{8} * this->th2pp) * this->th2;
      this->th3pp = Real{-8} * this->th3;
      this->th4pp = Real{8} * this->th4;

      const auto s_pi = emsr::pi_v<Real>;
      this->eta_1 = -s_pi * s_pi * this->th1ppp
		  / _Type{12} / lattice.omega_1() / this->th1p;
      const auto s_i = Cmplx{0, 1};
      this->eta_2 = (lattice.omega_2() * this->eta_1
		  + s_i * s_pi / Real{2}) / lattice.omega_1();
      this->eta_3 = (lattice.omega_3() * this->eta_1
		  - s_i * s_pi / Real{2}) / lattice.omega_1();
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
    struct weierstrass_roots_t
    {
      using _Type = typename jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;
      using Real = emsr::num_traits_t<_Type>;
      using Cmplx = std::complex<Real>;

      _Type e1, e2, e3;

      explicit
      weierstrass_roots_t(const jacobi_lattice_t<_Tp1, _Tp3>& lattice);

      weierstrass_roots_t(const jacobi_theta_0_t<_Tp1, _Tp3>& theta0,
			    _Tp1 omega1);

      /// Return the discriminant
      /// @f$ \Delta = 16(e_2 - e_3)^2(e_3 - e_1)^2(e_1 - e_2)^2 @f$.
      _Type
      delta() const
      {
	const auto del1 = e2 - e3;
	const auto del2 = e3 - e1;
	const auto del3 = e1 - e2;
	const auto del = del1 * del2 * del3;
	return _Type{16} * del * del;
      }
    };

  /**
   * Constructor for the Weierstrass roots.
   *
   * @param lattice The Jacobi lattice.
   */
  template<typename _Tp1, typename _Tp3>
    weierstrass_roots_t<_Tp1, _Tp3>::
    weierstrass_roots_t(const jacobi_lattice_t<_Tp1, _Tp3>& lattice)
#if cplusplus > 201403L
    : weierstrass_roots_t(jacobi_theta_0_t(lattice),
			    lattice.omega_1())
#else
    : weierstrass_roots_t(jacobi_theta_0_t<_Tp1, _Tp3>(lattice),
			    lattice.omega_1())
#endif
    { }


  /**
   * Constructor for the Weierstrass roots.
   *
   * @param theta0 Exponential theta functions of argument 0.
   * @param omega_1 The first lattice parameter.
   */
  template<typename _Tp1, typename _Tp3>
    weierstrass_roots_t<_Tp1, _Tp3>::
    weierstrass_roots_t(const jacobi_theta_0_t<_Tp1, _Tp3>& theta0,
			  _Tp1 omega_1)
    {
      const auto s_pi = emsr::pi_v<Real>;

      const auto th22 = theta0.th2 * theta0.th2;
      const auto th24 = th22 * th22;
      const auto th42 = theta0.th4 * theta0.th4;
      const auto th44 = th42 * th42;
      const auto fr = s_pi / omega_1;
      const auto fc = fr * fr / Real{12};

      e1 = fc * (th24 + Real{2} * th44);
      e2 = fc * (th24 - th44);
      e3 = fc * (Real{-2} * th24 - th44);
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
    struct weierstrass_invariants_t
    {
      using _Type = typename jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;
      using Real = emsr::num_traits_t<_Type>;
      using Cmplx = std::complex<Real>;

      _Type g_2, g_3;

      weierstrass_invariants_t(const jacobi_lattice_t<_Tp1, _Tp3>&);

      /// Return the discriminant @f$ \Delta = g_2^3 - 27 g_3^2 @f$.
      _Type
      delta() const
      {
	const auto g_2p3 = g_2 * g_2 * g_2;
	return g_2p3 - _Type{27} * g_3 * g_3;
      }

      /// Return Klein's invariant @f$ J = 1738 g_2^3 / (g_2^3 - 27 g_3^2) @f$.
      _Type
      klein_j() const
      {
	const auto g_2p3 = g_2 * g_2 * g_2;
	return _Type{1738} * g_2p3 / (g_2p3 - _Type{27} * g_3 * g_3);
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
    weierstrass_invariants_t<_Tp1, _Tp3>::
    weierstrass_invariants_t(const jacobi_lattice_t<_Tp1, _Tp3>& lattice)
    {
      const auto roots = weierstrass_roots_t<_Tp1, _Tp3>(lattice);
      g_2 = Real{2} * (roots.e1 * roots.e1
        		+ roots.e2 * roots.e2
        		+ roots.e3 * roots.e3);
      g_3 = Real{4} * roots.e1 * roots.e2 * roots.e3;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_1_sum(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;

      Tp sum{};
      Real sign{-1};
      for (std::size_t n = 0; n < s_max_iter; ++n)
	{
	  sign *= -1;
	  const auto term = sign
			    * std::pow(q, Real((n + 0.5L) * (n + 0.5L)))
			    * std::sin(Real(2 * n + 1) * x);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return Real{2} * sum;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_1_prod(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;
      const auto q2 = q * q;
      const auto q4 = q2 * q2;
      const auto cos2x = std::cos(Real{2} * x);

      auto q2n = Tp{1};
      auto q4n = Tp{1};
      auto prod = Tp{1};
      for (std::size_t n = 1; n < s_max_iter; ++n)
	{
	  q2n *= q2;
	  q4n *= q4;
	  const auto fact = (Real{1} - q2n)
			    * (Real{1} - Real{2} * q2n * cos2x + q4n);
	  prod *= fact;
	  if (std::abs(fact) < s_eps)
	    break;
	}

      return Real{2} * std::pow(q, Tp{0.25L}) * std::sin(x) * prod;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    std::complex<Tp>
    jacobi_theta_1(std::complex<Tp> q, std::complex<Tp> x)
    {
      using Real = emsr::num_traits_t<Tp>;
      using Cmplx = std::complex<Real>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Real>{0, 1};
      constexpr auto s_q_min = Real{0.001L};
      constexpr auto s_q_max = Real{0.95e-1L};

      if (std::isnan(q) || std::isnan(x))
	return Tp{s_NaN};
      else if (std::abs(q) >= Real{1})
	throw std::domain_error("jacobi_theta_1: nome q out of range");
      else if (std::abs(q) < s_q_min || std::abs(q) > s_q_max)
	return jacobi_theta_1_prod(q, x);
      else if (std::abs(x) < s_eps)
	return std::complex<Tp>{0, 0};
      else
	{
	  const auto lattice = jacobi_lattice_t<Cmplx, Cmplx>(q);
	  auto tau = lattice.tau().val;

	  const auto x_red = lattice.reduce(x);
	  auto fact = std::complex<Tp>{1, 0};
	  if (x_red.m != 0)
	    fact *= emsr::parity<Tp>(x_red.m);
	  if (x_red.n != 0)
	    fact *= emsr::parity<Tp>(x_red.n)
	    	    * std::exp(s_i * Real(-2 * x_red.n) * x_red.z)
		    * std::pow(q, -x_red.n * x_red.n);
	  x = x_red.z;

	  // theta_1(tau+1, z) = exp(i tau/4) theta_1(tau, z)
	  const auto itau = std::floor(std::real(tau));
	  tau -= itau;
	  fact *= polar_pi(Real{1}, itau / Real{4});

	  if (std::imag(tau) < 0.5)
	    {
	      const auto fact2 = s_i * std::sqrt(-s_i * tau);
	      tau = Real{-1} / tau;
	      const auto phase = std::exp(s_i * tau * x * x / s_pi);
	      fact *= phase / fact2;
	      q = std::exp(s_i * s_pi * tau);
	      x *= tau;
	    }

	  return fact * jacobi_theta_1_sum(q, x);
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_1(Tp q, const Tp x)
    {
      using Cmplx = std::complex<Tp>;

      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto ret = jacobi_theta_1(Cmplx(q), Cmplx(x));

      if (std::abs(ret) > s_eps
	  && std::abs(std::imag(ret)) > s_eps * std::abs(ret))
	throw std::runtime_error("jacobi_theta_1: "
				 "Unexpected large imaginary part");
      else
	return std::real(ret);
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_2_sum(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;

      Tp sum{};
      for (std::size_t n = 0; n < s_max_iter; ++n)
	{
	  const auto term = std::pow(q, Real((n + 0.5L) * (n + 0.5L)))
			    * std::cos(Real(2 * n + 1) * x);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return Real{2} * sum;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_2_prod(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;
      const auto q2 = q * q;
      const auto q4 = q2 * q2;
      const auto cos2x = std::cos(Real{2} * x);

      auto q2n = Tp{1};
      auto q4n = Tp{1};
      auto prod = Tp{1};
      for (std::size_t n = 1; n < s_max_iter; ++n)
	{
	  q2n *= q2;
	  q4n *= q4;
	  const auto fact = (Real{1} - q2n)
			    * (Real{1} + Real{2} * q2n * cos2x + q4n);
	  prod *= fact;
	  if (std::abs(fact) < s_eps)
	    break;
	}

      return Real{2} * std::pow(q, Tp{0.25L}) * std::cos(x) * prod;
    }

  // Pre-declare Jacobi theta_4 sum ...
  template<typename Tp>
    Tp
    jacobi_theta_4_sum(Tp q, Tp x);

  // ... and product.
  template<typename Tp>
    Tp
    jacobi_theta_4_prod(Tp q, Tp x);

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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    std::complex<Tp>
    jacobi_theta_2(std::complex<Tp> q, std::complex<Tp> x)
    {
      using Real = emsr::num_traits_t<Tp>;
      using Cmplx = std::complex<Real>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Real>{0, 1};
      constexpr auto s_q_min = Real{0.001L};
      constexpr auto s_q_max = Real{0.95e-1L};

      if (std::isnan(q) || std::isnan(x))
	return Tp{s_NaN};
      else if (std::abs(q) >= Real{1})
	throw std::domain_error("jacobi_theta_2: nome q out of range");
      else if (std::abs(q) < s_q_min || std::abs(q) > s_q_max)
	return jacobi_theta_2_prod(q, x);
      else if (std::abs(x) < s_eps)
#if cplusplus > 201403
	return jacobi_theta_0_t(jacobi_lattice_t<Cmplx, Cmplx>(q)).th2;
#else
	return jacobi_theta_0_t<Cmplx, Cmplx>(
				jacobi_lattice_t<Cmplx, Cmplx>(q)).th2;
#endif
      else
	{
	  const auto lattice = jacobi_lattice_t<Cmplx, Cmplx>(q);
	  auto tau = lattice.tau().val;

	  const auto x_red = lattice.reduce(x);
	  auto fact = std::complex<Tp>{1, 0};
	  if (x_red.m != 0)
	    fact *= emsr::parity<Tp>(x_red.m);
	  if (x_red.n != 0)
	    fact *= std::exp(s_i * Real(-2 * x_red.n) * x_red.z)
		    * std::pow(q, -x_red.n * x_red.n);
	  x = x_red.z;

	  // theta_2(tau+1, z) = theta_2(tau, z)
	  const auto itau = std::floor(std::real(tau));
	  tau -= itau;
	  fact *= polar_pi(Real{1}, itau / Real{4});

	  if (std::imag(tau) < 0.5)
	    {
	      const auto fact2 = std::sqrt(-s_i * tau);
	      tau = Real{-1} / tau;
	      const auto phase = std::exp(s_i * tau * x * x / s_pi);
	      fact *= phase / fact2;
	      q = std::exp(s_i * s_pi * tau);
	      x *= tau;
	      return fact * jacobi_theta_4_sum(q, x);
	    }
	  else
	    return fact * jacobi_theta_2_sum(q, x);
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_2(Tp q, const Tp x)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));

      const auto ret = jacobi_theta_2(Cmplx(q), Cmplx(x));

      if (std::abs(ret) > s_eps
	  && std::abs(std::imag(ret)) > s_eps * std::abs(ret))
	throw std::runtime_error("jacobi_theta_2: "
				 "Unexpected large imaginary part");
      else
	return std::real(ret);
    }

  /**
   * Return the Jacobi @f$ \theta_3 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta-3 function is defined by
   * @f[
   *  \theta_3(q,x) = 1 + 2\sum_{n=1}^{\infty} q^{n^2}\cos{2nx}
   * @f]
   *
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_3_sum(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;

      Tp sum{};
      for (std::size_t n = 1; n < s_max_iter; ++n)
	{
	  const auto term = std::pow(q, Real(n * n))
			    * std::cos(Real(2 * n) * x);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return Real{1} + Real{2} * sum;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_3_prod(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;
      const auto q2 = q * q;
      const auto q4 = q2 * q2;
      const auto cos2x = std::cos(Real{2} * x);

      auto q2nm1 = q;
      auto q4nm2 = q2;
      auto prod = Tp{1};
      for (std::size_t n = 1; n < s_max_iter; ++n)
	{
	  const auto fact = (Real{1} - q2nm1 * q)
			    * (Real{1} + Real{2} * q2nm1 * cos2x + q4nm2);
	  prod *= fact;
	  if (std::abs(fact) < s_eps)
	    break;
	  q2nm1 *= q2;
	  q4nm2 *= q4;
	}

      return prod;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    std::complex<Tp>
    jacobi_theta_3(std::complex<Tp> q, std::complex<Tp> x)
    {
      using Real = emsr::num_traits_t<Tp>;
      using Cmplx = std::complex<Real>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Real>{0, 1};
      constexpr auto s_q_min = Real{0.001L};
      constexpr auto s_q_max = Real{0.95e-1L};

      if (std::isnan(q) || std::isnan(x))
	return Tp{s_NaN};
      else if (std::abs(q) >= Real{1})
	throw std::domain_error("jacobi_theta_3: nome q out of range");
      else if (std::abs(q) < s_q_min || std::abs(q) > s_q_max)
	return jacobi_theta_3_prod(q, x);
      else if (std::abs(x) < s_eps)
#if cplusplus > 201403L
	return jacobi_theta_0_t(jacobi_lattice_t<Cmplx, Cmplx>(q)).th3;
#else
	return jacobi_theta_0_t<Cmplx, Cmplx>(
				jacobi_lattice_t<Cmplx, Cmplx>(q)).th3;
#endif
      else
	{
	  const auto lattice = jacobi_lattice_t<Cmplx, Cmplx>(q);
	  auto tau = lattice.tau().val;

	  const auto x_red = lattice.reduce(x);
	  auto fact = std::complex<Tp>{1, 0};
	  if (x_red.n != 0)
	    fact *= std::exp(s_i * Real(-2 * x_red.n) * x_red.z)
		    * std::pow(q, -x_red.n * x_red.n);
	  x = x_red.z;

	  // theta_3(tau+1, z) = theta_3(tau, z)
	  const auto itau = std::floor(std::real(tau));
	  tau -= itau;

	  if (std::imag(tau) < 0.5)
	    {
	      const auto fact2 = std::sqrt(-s_i * tau);
	      tau = Real{-1} / tau;
	      const auto phase = std::exp(s_i * tau * x * x / s_pi);
	      fact *= phase / fact2;
	      q = std::exp(s_i * s_pi * tau);
	      x *= tau;
	    }

	  return fact * jacobi_theta_3_sum(q, x);
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_3(Tp q, const Tp x)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));

      const auto ret = jacobi_theta_3(Cmplx(q), Cmplx(x));

      if (std::abs(ret) > s_eps
	  && std::abs(std::imag(ret)) > s_eps * std::abs(ret))
	throw std::runtime_error("jacobi_theta_3: "
				 "Unexpected large imaginary part");
      else
	return std::real(ret);
    }

  /**
   * Return the Jacobi @f$ \theta_4 @f$ function by summation of the series.
   *
   * The Jacobi or elliptic theta function is defined by
   * @f[
   *  \theta_4(q,x) = 1 + 2\sum_{n=1}^{\infty}(-1)^n q^{n^2}\cos{2nx}
   * @f]
   *
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_4_sum(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;

      Tp sum{};
      Real sign{1};
      for (std::size_t n = 1; n < s_max_iter; ++n)
	{
	  sign *= -1;
	  const auto term = sign * std::pow(q, Real(n * n))
			    * std::cos(Real(2 * n) * x);
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}
      return Real{1} + Real{2} * sum;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_4_prod(Tp q, Tp x)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));
      constexpr std::size_t s_max_iter = 50;
      const auto q2 = q * q;
      const auto q4 = q2 * q2;
      const auto cos2x = std::cos(Real{2} * x);

      auto q2nm1 = q;
      auto q4nm2 = q2;
      auto prod = Tp{1};
      for (std::size_t n = 1; n < s_max_iter; ++n)
	{
	  const auto fact = (Real{1} - q2nm1 * q)
			    * (Real{1} - Real{2} * q2nm1 * cos2x + q4nm2);
	  prod *= fact;
	  if (std::abs(fact) < s_eps)
	    break;
	  q2nm1 *= q2;
	  q4nm2 *= q4;
	}

      return prod;
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    std::complex<Tp>
    jacobi_theta_4(std::complex<Tp> q, std::complex<Tp> x)
    {
      using Real = emsr::num_traits_t<Tp>;
      using Cmplx = std::complex<Real>;
      const auto s_NaN = emsr::quiet_NaN(std::abs(x));
      const auto s_eps = emsr::epsilon(std::abs(x));
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Real>{0, 1};
      constexpr auto s_q_min = Real{0.001L};
      constexpr auto s_q_max = Real{0.95e-1L};

      if (std::isnan(q) || std::isnan(x))
	return Tp{s_NaN};
      else if (std::abs(q) >= Real{1})
	throw std::domain_error("jacobi_theta_4: nome q out of range");
      else if (std::abs(q) < s_q_min || std::abs(q) > s_q_max)
	return jacobi_theta_4_prod(q, x);
      else if (std::abs(x) < s_eps)
#if cplusplus > 201403L
	return jacobi_theta_0_t(jacobi_lattice_t<Cmplx, Cmplx>(q)).th4;
#else
	return jacobi_theta_0_t<Cmplx, Cmplx>(
				jacobi_lattice_t<Cmplx, Cmplx>(q)).th4;
#endif
      else
	{
	  const auto lattice = jacobi_lattice_t<Cmplx, Cmplx>(q);
	  auto tau = lattice.tau().val;

	  const auto x_red = lattice.reduce(x);
	  auto fact = std::complex<Tp>{1, 0};
	  if (x_red.n != 0)
	    fact *= std::exp(s_i * Real(-2 * x_red.n) * x_red.z)
		    * std::pow(q, -x_red.n * x_red.n);
	  if (x_red.n != 0)
	    fact *= emsr::parity<Tp>(x_red.n);
	  x = x_red.z;

	  // theta_4(tau+1, z) = theta_4(tau, z)
	  const auto itau = std::floor(std::real(tau));
	  tau -= itau;

	  if (std::imag(tau) < 0.5)
	    {
	      const auto fact2 = std::sqrt(-s_i * tau);
	      tau = Real{-1} / tau;
	      const auto phase = std::exp(s_i * tau * x * x / s_pi);
	      fact *= phase / fact2;
	      q = std::exp(s_i * s_pi * tau);
	      x *= tau;
	      return fact * jacobi_theta_2_sum(q, x);
	    }
	  else
	    return fact * jacobi_theta_4_sum(q, x);
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
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   * @param x The argument.
   */
  template<typename Tp>
    Tp
    jacobi_theta_4(Tp q, const Tp x)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_eps = emsr::epsilon(std::abs(x));

      const auto ret = jacobi_theta_4(Cmplx(q), Cmplx(x));

      if (std::abs(ret) > s_eps
	  && std::abs(std::imag(ret)) > s_eps * std::abs(ret))
	throw std::runtime_error("jacobi_theta_4: "
				 "Unexpected large imaginary part");
      else
	return std::real(ret);
    }

  /**
   * Return a structure containing the three primary Jacobi elliptic functions:
   * @f$ sn(k, u), cn(k, u), dn(k, u) @f$.
   *
   * @param k The elliptic modulus @f$ |k| < 1 @f$.
   * @param u The argument.
   * @return An object containing the three principal Jacobi elliptic functions,
   *         @f$ sn(k, u), cn(k, u), dn(k, u) @f$ and the means to compute
   *         the remaining nine as well as the amplitude.
   */
  template<typename Tp>
    jacobi_ellint_t<Tp>
    jacobi_ellint(Tp k, Tp u)
    {
      const auto s_eps = emsr::epsilon(std::abs(u));
      const auto s_NaN = emsr::quiet_NaN(std::abs(u));

      if (std::isnan(k) || std::isnan(u))
	return jacobi_ellint_t<Tp>{s_NaN, s_NaN, s_NaN};
      else if (std::abs(k) > Tp{1})
	throw std::domain_error("jacobi_ellint: argument k out of range");
      else if (std::abs(Tp{1} - k) < Tp{2} * s_eps)
	{
	  auto sn = std::tanh(u);
	  auto cn = Tp{1} / std::cosh(u);
	  auto dn = cn;
	  return jacobi_ellint_t<Tp>{sn, cn, dn};
	}
      else if (std::abs(k) < Tp{2} * s_eps)
	{
	  auto sn = std::sin(u);
	  auto cn = std::cos(u);
	  auto dn = Tp{1};
	  return emsr::jacobi_ellint_t<Tp>{sn, cn, dn};
	}
      else
	{
	  const auto s_CA = std::sqrt(s_eps * Tp{0.01});
	  const auto s_N = 100;
	  std::vector<Tp> m;
	  std::vector<Tp> n;
	  m.reserve(20);
	  n.reserve(20);
	  Tp c, d = Tp{1};
	  auto mc = Tp{1} - k * k;
	  const bool bo = (mc < Tp{0});
	  if (bo)
	    {
              mc /= -k * k;
              d = k;
              u *= d;
	    }
	  auto a = Tp{1};
	  auto dn = Tp{1};
	  auto l = s_N;
	  for (auto i = 0; i < s_N; ++i)
	    {
	      l = i;
	      m.push_back(a);
	      n.push_back(mc = std::sqrt(mc));
	      c = (a + mc) / Tp{2};
	      if (std::abs(a - mc) < s_CA * a)
		break;
	      mc *= a;
	      a = c;
	    }
	  u *= c;
	  auto sn = std::sin(u);
	  auto cn = std::cos(u);
	  if (sn != Tp{0})
	    {
	      a = cn / sn;
	      c *= a;
	      for (auto ii = l; ii >= 0; --ii)
		{
		  const auto b = m[ii];
		  a *= c;
		  c *= dn;
		  dn = (n[ii] + a) / (b + a);
		  a = c / b;
		}
	      a = Tp{1} / std::hypot(Tp{1}, c);
	      sn = std::copysign(a, sn);
	      cn = c * sn;
	    }
	  if (bo)
	    {
	      std::swap(dn, cn);
	      sn /= d;
	    }
	  return jacobi_ellint_t<Tp>{sn, cn, dn};
	}
    }
} // namespace detail
} // namespace emsr

#endif // SF_THETA_TCC
