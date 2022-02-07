
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

/** @file bits/sf_polylog.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Florian Goth and Edward Smith-Rowland.
//
// References:
// (1) David C. Wood, "The Computation of Polylogarithms."
//

#ifndef SF_POLYLOG_TCC
#define SF_POLYLOG_TCC 1

#include <stdexcept>
#include <complex>

#include <emsr/complex_util.h>
#include <emsr/fp_type_util.h>
#include <emsr/sf_trig.h>
#include <emsr/sf_zeta.h>

namespace emsr
{
namespace detail
{

  /**
   * This class manages the termination of series.
   * Termination conditions involve both a maximum iteration count
   * and a relative precision.
   */
  template<typename Tp>
    class Terminator
    {
    private:

      using Real = emsr::num_traits_t<Tp>;
      const std::size_t m_max_iter;
      std::size_t m_curr_iter;
      Real m_toler;

    public:

      Terminator(std::size_t max_iter, Real mul = Real{1})
      : m_max_iter(max_iter), m_curr_iter{0},
	m_toler(std::abs(mul) * std::numeric_limits<Real>::epsilon())
      { }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->m_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or the maximum number of terms has been reached.
      bool
      operator()(Tp term, Tp sum)
      {
	if (this->m_curr_iter >= this->m_max_iter
	    || ++this->m_curr_iter == this->m_max_iter)
	  return true;
	else if (std::abs(term) < this->m_toler * std::abs(sum))
	  return true;
	else
	  return false;
      }
    };

  /**
   * This class manages the termination of asymptotic series.
   * In particular, this termination watches for the growth of the sequence
   * of terms to stop the series.
   *
   * Termination conditions involve both a maximum iteration count
   * and a relative precision.
   */
  template<typename Tp>
    class AsympTerminator
    {
    private:

      using Real = emsr::num_traits_t<Tp>;
      const std::size_t m_max_iter;
      std::size_t m_curr_iter;
      Real m_toler;
      Real m_prev_term = std::numeric_limits<Real>::max();
      bool m_stop_asymp = false;

    public:

      AsympTerminator(std::size_t max_iter, Real mul = Real{1})
      : m_max_iter(max_iter), m_curr_iter{0},
	m_toler(std::abs(mul) * std::numeric_limits<Real>::epsilon())
      { }

      /// Filter a term before applying it to the sum.
      Tp
      operator<<(Tp term)
      {
	if (std::abs(term) > this->m_prev_term)
	  {
	    this->m_stop_asymp = true;
	    return Tp{0};
	  }
	else
	  return term;
      }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->m_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or because the terms are starting to grow or
      //  the maximum number of terms has been reached.
      bool
      operator()(Tp term, Tp sum)
      {
	if (this->m_stop_asymp)
	  return true;
	else
	  {
	    const auto aterm = std::abs(term);
	    this->m_stop_asymp = (aterm > this->m_prev_term);
	    this->m_prev_term = aterm;
	    if (this->m_curr_iter >= this->m_max_iter
	    || ++this->m_curr_iter == this->m_max_iter)
	      return true;
	    else if (aterm < this->m_toler * std::abs(sum))
	      return true;
	    else if (this->m_stop_asymp)
	      return true;
	    else
	      return false;
	  }
      }
    };

  template<typename Tp>
    std::complex<Tp>
    clamp_pi(std::complex<Tp> z)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i2pi = std::complex<Tp>{0, Tp{2} * s_pi};
      while (z.imag() > s_pi)
	z -= s_i2pi;
      while (z.imag() <= -s_pi)
	z += s_i2pi;
      return z;
    }

  template<typename Tp>
    std::complex<Tp>
    clamp_0_m2pi(std::complex<Tp> z)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_2pi = emsr::tau_v<Real>;
      while (z.imag() > Tp{0})
	z = std::complex<Tp>(z.real(), z.imag() - s_2pi);
      while (z.imag() <= -s_2pi)
	z = std::complex<Tp>(z.real(), z.imag() + s_2pi);
      return z;
    }

  /**
   * This function treats the cases of positive integer index s
   * for complex argument w.
   *
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) \frac{w^k}{k!}
   *             + \left[H_{s-1} - \log(-w)\right] \frac{w^{s-1}}{(s-1)!}
   * @f]
   * The radius of convergence is @f$ |w| < 2 \pi @f$.
   * Note that this series involves a @f$ \log(-x) @f$.
   * gcc and Mathematica differ in their implementation
   * of @f$ \log(e^{i\pi}) @f$:
   * gcc: @f$ \log(e^{+-i\pi}) = +-i\pi @f$
   * whereas Mathematica doesn't preserve the sign in this case:
   * @f$ \log(e^{+- i\pi}) = +i \pi @f$
   *
   * @param s the positive integer index.
   * @param w the argument.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_pos(unsigned int s, std::complex<Tp> w)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_2pi = emsr::tau_v<Real>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_pipio6 = emsr::pi_sqr_div_6_v<Real>;
      std::complex<Tp> res = emsr::detail::riemann_zeta<Tp>(s);
      auto wk = w;
      auto fact = Tp{1};
      auto harmonicN = Tp{1}; // HarmonicNumber_1
      for (unsigned int k = 1; k <= s - 2; ++k)
	{
	  res += fact * emsr::detail::riemann_zeta<Tp>(s - k) * wk;
	  wk *= w;
	  const auto temp = Tp{1} / Tp(1 + k);
	  fact *= temp;
	  harmonicN += temp;
	}
      // harmonicN now contains H_{s-1}.
      // fact should be 1/(s-1)!
      res += (harmonicN - std::log(-w)) * wk * fact;
      wk *= w;
      fact /= s; // 1/s!
      res -= wk * fact / Tp{2};
      wk *= w;
      // Now comes the remainder of the series.
      const auto pref = wk / s_pi / s_2pi;
      fact /= Tp(s + 1); // 1/(s+1)!
      // Subtract the zeroth order term.
      res -= s_pipio6 * fact * pref;
      fact *= Tp{2} / Tp(s + 2) * Tp{3} / Tp(s + 3);
      const auto wbar = w / s_2pi;
      const auto w2 = -wbar * wbar;
      auto w2k = w2;
      auto rzarg = Tp{2};
      const unsigned int maxit = 200;
      Terminator<std::complex<Tp>> done(maxit);
      while (true)
	{
	  rzarg += Tp{2};
	  const auto rzeta = emsr::detail::riemann_zeta(rzarg);
	  const auto term = pref * fact * rzeta * w2k;
	  res -= term;
	  if (done(term, res))
	    break;
	  w2k *= w2;
	  fact *= Tp(rzarg) / Tp(s + rzarg)
		  * Tp(rzarg + 1) / Tp(s + rzarg + 1);
	}
      return res;
    }

  /**
   * This function treats the cases of positive integer index s
   * for real argument w.
   *
   * This specialization is worthwhile to catch the differing behaviour
   * of log(x).
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) \frac{w^k}{k!}
   *             + \left[H_{s-1} - \log(-w)\right] \frac{w^{s-1}}{(s-1)!}
   * @f]
   * The radius of convergence is @f$ |w| < 2 \pi @f$.
   * Note that this series involves a @f$ \log(-x) @f$.
   * gcc and Mathematica differ in their implementation
   * of @f$ \log(e^{i\pi}) @f$:
   * gcc: @f$ \log(e^{+-i\pi}) = +-i\pi @f$
   * whereas Mathematica doesn't preserve the sign in this case:
   * @f$ \log(e^{+- i\pi}) = +i\pi @f$
   *
   * @param s the positive integer index.
   * @param w the argument.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_pos(unsigned int s, Tp w)
    {
      const auto s_2pi = emsr::tau_v<Tp>;
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pipio6 = emsr::pi_sqr_div_6_v<Tp>;
      auto res = emsr::detail::riemann_zeta<Tp>(s);
      auto wk = w;
      auto fact = Tp{1};
      auto harmonicN = Tp{1}; // HarmonicNumber_1
      for (unsigned int k = 1; k <= s - 2; ++k)
	{
	  res += fact * emsr::detail::riemann_zeta<Tp>(s - k) * wk;
	  wk *= w;
	  const auto temp = Tp{1} / Tp(1 + k);
	  fact *= temp;
	  harmonicN += temp;
	}
      // harmonicN now contains H_{s-1}
      // fact should be 1/(s-1)!
      const auto imagtemp = fact * wk
			    * (harmonicN - std::log(std::complex<Tp>(-w)));
      res += std::real(imagtemp);
      wk *= w;
      fact /= s; // 1/s!
      res -= wk * fact / Tp{2};
      wk *= w;
      // Now comes the remainder of the series.
      const auto pref = wk / s_pi / s_2pi;
      fact /= Tp(s + 1); // 1/(s+1)!
      // Subtract the zeroth order term.
      res -= s_pipio6 * fact * pref;
      fact *= Tp{2} / Tp(s + 2) * Tp{3} / Tp(s + 3);
      const auto wbar = w / s_2pi;
      const auto w2 = -wbar * wbar;
      auto w2k = w2;
      auto rzarg = Tp{2};
      const unsigned int maxit = 200;
      Terminator<Tp> done(maxit);
      while (true)
	{
	  rzarg += Tp{2};
	  const auto rzeta = emsr::detail::riemann_zeta(rzarg);
	  const auto term = pref * fact * rzeta * w2k;
	  res -= term;
	  if (done(term, res))
	    break;
	  w2k *= w2;
	  fact *= Tp(rzarg) / Tp(s + rzarg)
		  * Tp(rzarg + 1) / Tp(s + rzarg + 1);
	}
      return std::complex<Tp>(res, std::imag(imagtemp));
    }

  /**
   * This function treats the cases of negative real index s.
   * Theoretical convergence is present for @f$ |w| < 2\pi @f$.
   * We use an optimized version of
   * @f[
   *   Li_s(e^w) = \Gamma(1-s)(-w)^{s-1} + \frac{(2\pi)^{-s}}{\pi} A_p(w)
   * @f]
   * @f[
   *   A_p(w) = \sum_k \frac{\Gamma(1+k-s)}{k!}
   *         \sin\left(\frac{\pi}{2} (s-k)\right)
   *         \left(\frac{w}{2\pi}\right)^k \zeta(1+k-s)
   * @f]
   * @param s  The negative real index
   * @param w  The complex argument
   * @return  The value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg(Tp s, std::complex<Tp> w)
    {
      const auto s_i = std::complex<Tp>{0, 1};
      const auto s_2pi = emsr::tau_v<Tp>;
      // Basic general loop, but s is a negative quantity here
      // FIXME Large s makes problems.
      // The series should be rearrangeable so that we only need
      // the ratio Gamma(1-s)/(2 pi)^s
      auto ls = log_gamma(Tp{1} - s);
      auto res = std::exp(ls - (Tp{1} - s) * std::log(-w));
      const auto wup = w / s_2pi;
      auto w2k = wup;
      const auto pref = Tp{2} * std::pow(s_2pi, -(Tp{1} - s));
      // Here we factor up the ratio of Gamma(1 - s + k)/k! .
      // This ratio should be well behaved even for large k in the series
      // afterwards
      // Note that we have a problem for large s.
      // Since s is negative we evaluate the Gamma Function
      // on the positive real axis where it is real.
      auto gam = std::exp(ls);

      const auto phase = emsr::polar_pi(Tp{1}, s / Tp{2});
      const auto cp = std::real(phase);
      const auto sp = std::imag(phase);
      // Here we add the expression that would result from ignoring
      // the zeta function in the series.
      const auto p = s_2pi - s_i * w;
      const auto q = s_2pi + s_i * w;
      // This can be optimized for real values of w
      res += s_i * gam * (std::conj(phase) * std::pow(p, s - Tp{1})
	     - phase * std::pow(q, s - Tp{1}));
      // The above expression is the result of
      // sum_k Gamma(1+k-s)/k! * sin(pi (s-k)/2) (w/2/pi)^k
      // Therefore we only need to sample values of zeta(n) on the real axis
      // that really differ from one
      std::complex<Tp> sum = sp * gam * emsr::detail::riemann_zeta_m_1(Tp{1} - s);
      unsigned int j = 1;
      gam *= (Tp{1} - s);
      constexpr unsigned int maxit = 200;
      Terminator<std::complex<Tp>> done(maxit);
      while (true)
	{
	  const auto rzarg = Tp(1 + j) - s;
	  const auto rz = emsr::detail::riemann_zeta_m_1(rzarg);
	  Tp sine;
	  // Save repeated recalculation of the sines.
	  if (j & 1)
	    { // odd
	      sine = cp;
	      if (!((j - 1) / 2 & 1))
		sine = -sine;
	    }
	  else
	    { // even
	      sine = sp;
	      if ((j / 2) & 1)
		sine = -sine;
	    }
	  const auto term =  w2k * (gam * sine * rz);
	  w2k *= wup;
	  ++j;
	  gam  *= rzarg / Tp(j); // == 1/(j+1) we incremented j above.
	  sum += term;
	  if (done(term, sum))
	    break;
	}
      res += pref * sum;
      return res;
    }

  /**
   * Compute the polylogarithm for negative integer order.
   * @f[
   *   Li_{-p}(e^w) = p!(-w)^{-(p+1)}
   *     - \sum_{k=0}^{\infty} \frac{B_{p+2k+q+1}}{(p+2k+q+1)!}
   *                           \frac{(p+2k+q)!}{(2k+q)!}w^{2k+q}
   * @f]
   * where @f$ q = (p+1) mod 2 @f$.
   *
   * @param n the negative integer index @f$ n = -p @f$.
   * @param w the argument w.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg(int n, std::complex<Tp> w)
    {
      const auto s_inf = std::numeric_limits<Tp>::infinity();
      if (emsr::fp_is_zero(w))
	return std::complex<Tp>{0};
      else if (emsr::fp_is_equal(w, Tp{1}))
	return std::complex<Tp>{s_inf, Tp{0}};
      else
	{
	  const int p = -n;
	  const int pp = 1 + p;
	  const int q = p & 1 ? 0 : 1;
	  const auto w2 = w * w;
	  auto wp = p & 1 ? std::complex<Tp>{1} : w;
	  unsigned int __2k = q;
	  auto gam = factorial<Tp>(p + __2k);
	  const auto pfact = factorial<Tp>(p);
	  auto res = pfact * std::pow(-w, Tp(-pp));
	  auto sum = std::complex<Tp>{};
	  constexpr unsigned int maxit = 300;
	  Terminator<std::complex<Tp>> done(maxit);
	  while (true)
	    {
	      const auto id = (p + __2k + 1) / 2;
	      if (id == emsr::detail::Num_Euler_Maclaurin_zeta)
		break;
	      const auto term = gam * wp
			* Tp(emsr::detail::Euler_Maclaurin_zeta[id]);
	      sum += term;
	      if (done(term, sum))
		break;
	      gam *= Tp(p + __2k + 1) / Tp(__2k + 1)
		     * Tp(p + __2k + 2) / Tp(__2k + 2);
	      wp *= w2;
	      __2k += 2;
	    }
	  res -= sum;
	  return res;
	}
    }

  /**
   * This function treats the cases of positive real index s.
   *
   * The defining series is
   * @f[
   *   Li_s(e^w) = A_s(w) + B_s(w) + \Gamma(1-s)(-w)^{s-1}
   * @f]
   * with
   * @f[
   *   A_s(w) = \sum_{k=0}^{m} \zeta(s-k)w^k/k!
   * @f]
   * @f[
   *   B_s(w) = \sum_{k=m+1}^{\infty} \sin(\pi/2(s-k))
   *             \Gamma(1-s+k)\zeta(1-s+k) (w/2/\pi)^k/k!
   * @f]
   *
   * @param s the positive real index s.
   * @param w The complex argument w.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_pos(Tp s, std::complex<Tp> w)
    { // positive s
      const auto s_2pi = emsr::tau_v<Tp>;
      const auto s_pi = emsr::pi_v<Tp>;
      std::complex<Tp> res = emsr::detail::riemann_zeta(s);
      auto wk = w;
      const auto phase = emsr::polar_pi(Tp{1}, s / Tp{2});
      const auto cp = std::real(phase);
      const auto sp = std::imag(phase);
      // This is \Gamma(1-s)(-w)^{s-1}
      res += s_pi / (Tp{2} * sp * cp)
	  * std::exp(-log_gamma(s) + (s - Tp{1}) * std::log(-w));
      auto fact = Tp{1};
      const auto m = static_cast<unsigned int>(std::floor(s));
      for (unsigned int k = 1; k <= m; ++k)
	{
	  res += wk * fact
		 * emsr::detail::riemann_zeta(s - Tp(k));
	  wk *= w;
	  fact /= Tp(1 + k);
	}
      // fac should now be 1/(m+1)!
      const auto pref = Tp{2} * std::pow(s_2pi, s - Tp{1});
      // Factor this out for now so we can compare with sum.
      res /= pref;
      // Now comes the remainder of the series
      unsigned int j = 0;
      constexpr unsigned int maxit = 100;
      Terminator<std::complex<Tp>> done(maxit);
      auto wup = w / s_2pi;
      auto wbark = std::pow(wup, Tp(m + 1));
      // It is 1 < 2 - s + m < 2 => Gamma(2-s+m) will not overflow
      // Here we factor up the ratio of Gamma(1 - s + k) / k!.
      // This ratio should be well behaved even for large k
      auto gam = gamma(Tp(2 + m) - s) * fact;
      std::complex<Tp> sum{};
      while (true)
	{
	  const auto idx = m + 1 + j;
	  const auto zetaarg = Tp(1 + idx) - s;
	  const auto rz = emsr::detail::riemann_zeta(zetaarg);
	  auto sine = cp;
	  if (idx & 1) // Save the repeated calculation of the sines.
	    { // odd
	      sine = cp;
	      if (!((idx - 1) / 2 & 1))
		sine = -sine;
	    }
	  else
	    { // even
	      sine = sp;
	      if ((idx / 2) & 1)
		sine = -sine;
	    }
	  const auto term = wbark * sine * gam * rz;
	  wbark *= wup;
	  gam *= zetaarg / Tp(1 + idx);
	  ++j;
	  sum += term;
	  if (done(term, res + sum))
	    break;
	}
      res += sum;
      return pref * res;
    }

  /**
   * This function implements the asymptotic series for the polylog.
   * It is given by
   * @f[
   *    2 \sum_{k=0}^{\infty} \zeta(2k) w^{s-2k}/\Gamma(s-2k+1)
   *       -i \pi w^{s-1}/\Gamma(s)
   * @f]
   * for @f$ Re(w) >> 1 @f$
   *
   * Don't check this against Mathematica 8.
   * For real w the imaginary part of the polylog is given by
   * @f$ Im(Li_s(e^w)) = -\pi w^{s-1}/\Gamma(s) @f$.
   * Check this relation for any benchmark that you use.
   *
   * @param s the real index s.
   * @param w the large complex argument w.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_asymp(Tp s, std::complex<Tp> w)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      // wgamma = w^{s-1} / Gamma(s)
      auto wgamma = std::pow(w, s - Tp{1}) * gamma_reciprocal(s);
      auto res = std::complex<Tp>(Tp{0}, -s_pi) * wgamma;
      // wgamma = w^s / Gamma(s+1)
      wgamma *= w / s;
      constexpr unsigned int maxiter = 100;
      AsympTerminator<std::complex<Tp>> done(maxiter);
      // zeta(0) w^s / Gamma(s + 1)
      std::complex<Tp> oldterm = -Tp{0.5L} * wgamma;
      res += Tp{2} * oldterm;
      std::complex<Tp> term;
      auto wq = Tp{1} / (w * w);
      int k = 1;
      while (true)
	{
	  wgamma *= wq * (s + Tp(1 - 2 * k)) * (s + Tp(2 - 2 * k));
	  term = emsr::detail::riemann_zeta<Tp>(2 * k) * wgamma;
	  res += done << Tp{2} * term;
	  if (done(Tp{2} * term, res))
	    break;
	  oldterm = term;
	  ++k;
	}
      return res;
    }

  /**
   * Theoretical convergence for Re(w) < 0.
   *
   * Seems to beat the other expansions for @f$ Re(w) < -\pi/2 - \pi/5 @f$.
   * Note that this is an implementation of the basic series:
   * @f[
   *   Li_s(e^z) = \sum_{k=1}^{\infty} e^{kz} k^{-s}
   * @f]
   *
   * @param s is an arbitrary type, integral or float.
   * @param w something with a negative real part.
   * @return the value of the polylogarithm.
   */
  template<typename PowTp, typename Tp>
    Tp
    polylog_exp_sum(PowTp s, Tp w)
    {
      auto ew = std::exp(w);
      const auto up = ew;
      auto res = ew;
      unsigned int maxiter = 500;
      Terminator<Tp> done(maxiter);
      bool terminate = false;
      unsigned int k = 2;
      while (!terminate)
	{
	  ew *= up;
	  Tp temp = std::pow(k, s); // This saves us a type conversion
	  const auto term = ew / temp;
	  res += term;
	  terminate = done(term, res);
	  ++k;
	}
      return res;
  }

  /**
   * Here s is a positive integer and the function descends
   * into the different kernels depending on w.
   *
   * @param s a positive integer.
   * @param w an arbitrary complex number.
   * @return The value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_pos_int(unsigned int s, std::complex<Tp> w)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_2pi = emsr::tau_v<Real>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};
      const auto s_max_asymp = Tp{5};
      const auto rw = w.real();
      const auto iw = w.imag();
      if (emsr::fp_is_real(w)
	  && emsr::fp_is_equal(std::remainder(iw, s_2pi), Tp{0}))
	{
	  if (s == 1)
	    return std::numeric_limits<Tp>::infinity();
	  else
	    return emsr::detail::riemann_zeta<Tp>(s);
	}
      else if (0 == s)
	{
	  const auto t = std::exp(w);
	  return emsr::fp_is_zero(Tp{1} - t)
	       ? std::numeric_limits<Tp>::quiet_NaN()
	       : t / (Tp{1} - t);
	}
      else if (1 == s)
	{
	  const auto t = std::exp(w);
	  return emsr::fp_is_zero(Tp{1} - t)
	       ? std::numeric_limits<Tp>::quiet_NaN()
	       : -std::log(Tp{1} - t);
	}
      else
	{
	  if (rw < -(s_pi_2 + s_pi / Tp{5}))
	    // Choose the exponentially converging series
	    return polylog_exp_sum(s, w);
	  else if (rw < s_max_asymp)
	    // The transition point chosen here, is quite arbitrary
	    // and needs more testing.
	    // The reductions of the imaginary part yield the same results
	    // as Mathematica.
	    // Necessary to improve the speed of convergence
	    return polylog_exp_pos(s, clamp_pi(w));
	  else
	    // Wikipedia says that this is required for Wood's formula.
	    return polylog_exp_asymp(static_cast<Tp>(s),
				       clamp_0_m2pi(w));
	}
    }

  /**
   * Here s is a positive integer and the function descends
   * into the different kernels depending on w.
   *
   * @param s a positive integer
   * @param w an arbitrary real argument w
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_pos_int(unsigned int s, Tp w)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};
      const auto s_max_asymp = Tp{5};
      if (emsr::fp_is_zero(w))
	{
	  if (s == 1)
	    return std::numeric_limits<Tp>::infinity();
	  else
	    return emsr::detail::riemann_zeta<Tp>(s);
	}
      else if (s == 0)
	{
	  const auto t = std::exp(w);
	  return emsr::fp_is_zero(Tp{1} - t)
	       ? std::numeric_limits<Tp>::infinity()
	       : t / (Tp{1} - t);
	}
      else if (s == 1)
	{
	  const auto t = std::exp(w);
	  return emsr::fp_is_zero(Tp{1} - t)
	       ? -std::numeric_limits<Tp>::infinity()
	       : -std::log(Tp{1} - t);
	}
      else
	{
	  if (w < -(s_pi_2 + s_pi / Tp{5}))
	    // Choose the exponentially converging series
	    return polylog_exp_sum(s, w);
	  else if (w < s_max_asymp)
	    return polylog_exp_pos(s, w);
	  else
	    return polylog_exp_asymp(static_cast<Tp>(s),
				       std::complex<Tp>(w));
	}
    }

  /**
   * This treats the case where s is a negative integer.
   *
   * @param s a negative integer.
   * @param w an arbitrary complex number
   * @return the value of the polylogarith,.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg_int(int s, std::complex<Tp> w)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_2pi = emsr::tau_v<Real>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};
      const auto s_max_asymp = Tp{5};
      if ((((-s) & 1) == 0) && emsr::fp_is_imag(w))
	{
	  // Now s is odd and w on the unit-circle.
	  const auto iw = imag(w); // Get imaginary part.
	  const auto rem = std::remainder(iw, s_2pi);
	  if (emsr::fp_is_equal(std::abs(rem), Tp{0.5L}))
	    // Due to: Li_{-n}(-1) + (-1)^n Li_{-n}(1/-1) = 0.
	    return Tp{0};
	  else
	    // No asymptotic expansion available... check the reduction.
	    return polylog_exp_neg(s, std::complex<Tp>(w.real(), rem));
	}
      else
	{
	  if (std::real(w) < -(s_pi_2 + s_pi / Tp{5}))
	    // Choose the exponentially converging series
	    return polylog_exp_sum(s, w);
	  else if (std::real(w) < s_max_asymp)
	    // Arbitrary transition point...
	    // The reductions of the imaginary part yield the same results
	    // as Mathematica.
	    // Necessary to improve the speed of convergence.
	    return polylog_exp_neg(s, clamp_pi(w));
	  else
	    // Wikipedia says that this clamping is required for Wood's formula.
	    return polylog_exp_asymp(Tp(s), clamp_0_m2pi(w));
	}
    }

  /**
   * This treats the case where s is a negative integer and w is a real.
   *
   * @param s a negative integer.
   * @param w the argument.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg_int(int s, Tp w)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_max_asymp = Tp{5};
      if (w < -(s_pi_2 + s_pi / Tp{5}))
	// Choose exponentially converging series.
	return polylog_exp_sum(s, w);
      else if (emsr::fp_is_zero(w))
	return std::numeric_limits<Tp>::infinity();
      else if (w < s_max_asymp)
	// Arbitrary transition point less than 2 pi.
	return polylog_exp_neg(s, std::complex<Tp>(w));
      else
	return polylog_exp_asymp(Tp(s), std::complex<Tp>(w));
    }

  /**
   * Return the polylog where s is a positive real value
   * and for complex argument.
   *
   * @param s A positive real number.
   * @param w the complex argument.
   * @return The value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_pos_real(Tp s, std::complex<Tp> w)
    {
      const auto s_2pi = emsr::tau_v<Tp>;
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_max_asymp = Tp{5};
      const auto rw = w.real();
      const auto iw = w.imag();
      if (emsr::fp_is_real(w)
	  && emsr::fp_is_zero(std::remainder(iw, s_2pi)))
	{
	  if (emsr::fp_is_equal(s, Tp{1}))
	    return std::numeric_limits<Tp>::infinity();
	  else
	    return emsr::detail::riemann_zeta(s);
	}
      if (rw < -(s_pi_2 + s_pi / Tp{5}))
        // Choose exponentially converging series.
	return polylog_exp_sum(s, w);
      if (rw < s_max_asymp)
	// Arbitrary transition point.
	// The reductions of the imaginary part yield the same results
	// as Mathematica then.
	// Branch cuts??
	return polylog_exp_pos(s, clamp_pi(w));
      else
	// Wikipedia says that this is required for Wood's formula
	return polylog_exp_asymp(s, clamp_0_m2pi(w));
    }

  /**
   * Return the polylog where s is a positive real value and the argument
   * is real.
   *
   * @param s  A positive real number tht does not reduce to an integer.
   * @param w  The real argument w.
   * @return  The value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_pos_real(Tp s, Tp w)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_max_asymp = Tp{5};
      if (emsr::fp_is_zero(w))
	{
	  if (emsr::fp_is_equal(s, Tp{1}))
	    return std::numeric_limits<Tp>::infinity();
	  else
	    return emsr::detail::riemann_zeta(s);
	}
      if (w < -(s_pi_2 + s_pi / Tp{5}))
	// Choose exponentially converging series.
	return polylog_exp_sum(s, w);
      if (w < s_max_asymp)
	// Arbitrary transition point
	return polylog_exp_pos(s, std::complex<Tp>(w));
      else
	return polylog_exp_asymp(s, std::complex<Tp>(w));
    }

  /**
   * Return the polylog where s is a negative real value
   * and for complex argument.
   * Now we branch depending on the properties of w in the specific functions
   *
   * @param s  A negative real value that does not reduce
   *             to a negative integer.
   * @param w  The complex argument.
   * @return  The value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg_real(Tp s, std::complex<Tp> w)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_max_asymp = Tp{5};
      const auto rw = w.real();
      if (rw < -(s_pi_2 + s_pi / Tp{5}))
	// Choose exponentially converging series.
	return polylog_exp_sum(s, w);
      else if (rw < s_max_asymp)
	// The reductions of the imaginary part yield the same results
	// as Mathematica then.
	// Necessary to improve the speed of convergence.
	// Branch cuts??
	return polylog_exp_neg(s, clamp_pi(w));
      else
	// Wikipedia says that this is required for Wood's formula
	return polylog_exp_asymp(s, clamp_0_m2pi(w));
    }

  /**
   * Return the polylog where s is a negative real value and for real argument.
   * Now we branch depending on the properties of w in the specific functions.
   *
   * @param s  A negative real value.
   * @param w  A real argument.
   * @return  The value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg_real(Tp s, Tp w)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_max_asymp = Tp{5};
      if (w < -(s_pi_2 + s_pi / Tp{5}))
	// Choose exponentially converging series.
	return polylog_exp_sum(s, w);
      else if (w < s_max_asymp)
	// Arbitrary transition point
	return polylog_exp_neg(s, std::complex<Tp>(w));
      else
	return polylog_exp_asymp(s, std::complex<Tp>(w));
    }

  /**
   * This is the frontend function which calculates @f$ Li_s(e^w) @f$
   * First we branch into different parts depending on the properties of s.
   * This function is the same irrespective of a real or complex w,
   * hence the template parameter ArgType.
   *
   * @note: I *really* wish we could return a variant<Tp, std::complex<Tp>>.
   *
   * @param s  The real order.
   * @param w  The real or complex argument.
   * @return  The real or complex value of Li_s(e^w).
   */
  template<typename Tp, typename ArgType>
    emsr::fp_promote_t<std::complex<Tp>, ArgType>
    polylog_exp(Tp s, ArgType w)
    {
      if (std::isnan(s) || std::isnan(w))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (s > Tp{25})
	// Cutoff chosen by some testing on the real axis.
	return polylog_exp_sum(s, w);
      else
	{
	  const auto p = emsr::fp_is_integer(s, Tp{5});
	  if (p)
	    { // The order s is an integer.
	      if (p() >= 0)
		return polylog_exp_pos_int(p(), w);
	      else
		return polylog_exp_neg_int(p(), w);
	    }
	  else
	    {
	      if (s > Tp{0})
		return polylog_exp_pos_real(s, w);
	      else
		return polylog_exp_neg_real(s, w);
	    }
	}
    }

  /**
   * Return the polylog Li_s(x) for two real arguments.
   * 
   * The polylog is defined by
   * @f[
   *    Li_s(x) = \sum_{k=1}^{\infty} \frac{x^k}{k^s}
   * @f]
   *
   * @param s  The real order.
   * @param x  The real argument.
   * @return The complex value of the polylogarithm.
   */
  template<typename Tp>
    Tp
    polylog(Tp s, Tp x)
    {
      if (std::isnan(s) || std::isnan(x))
	return emsr::quiet_NaN(s);
      else if (emsr::fp_is_zero(x))
	return Tp{0};
      else
	{
	  const auto n = emsr::fp_is_integer(s, Tp{5});
	  if (n && n() == 1)
	    return -std::log(Tp{1} - x);
	  else if (n && n() == 0)
	    return x / (Tp{1} - x);
	  else if (emsr::fp_is_equal(x, Tp{-1}))
	    // Prevent blowups caused by reflecting the branch point.
	    return std::real(polylog_exp(s, Tp{0})
				* (std::pow(Tp{2}, Tp{1} - s) - Tp{1}));
	  else if (x < Tp{0})
	    { // Use the reflection formula to access negative values.
	      const auto y = std::log(-x);
	      return std::real(polylog_exp(s, Tp{2} * y)
				    * std::pow(Tp{2}, Tp{1} - s)
			     - polylog_exp(s, y));
	    }
	  else
	    {
	      const auto y = std::log(x);
	      return std::real(polylog_exp(s, y));
	    }
	}
    }

  /**
   * Return the polylog in those cases where we can calculate it.
   *
   * @param s  The real order.
   * @param w  The complex argument.
   * @return  The complex value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog(Tp s, std::complex<Tp> w)
    {
      if (std::isnan(s) || std::isnan(w))
	return emsr::quiet_NaN(s);
      else if (emsr::fp_is_real(w))
	return polylog(s, std::real(w));
      else
	return polylog_exp(s, std::log(w));
    }

 /**
  * Return the periodic zeta function F(z,s) for two real arguments.
  *
  * The periodic zeta function is defined by
  * @f[
  *    F(z,s) = \sum_{k=1}^{\infty} \frac{e^{i2\pi kz}}{k^s} = Li_s(e^{i2\pi kz})
  * @f]
  *
  * @param s  The real order.
  * @param z  The real or complex argument.
  * @return The complex value of the periodic zeta function.
  */
  template<typename Tp, typename ArgType>
    emsr::fp_promote_t<std::complex<Tp>, ArgType>
    periodic_zeta(ArgType z, Tp s)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_i = Cmplx{0, 1};
      const auto s_pi = emsr::pi_v<Tp>;
      if (std::isnan(s) || std::isnan(z))
	return emsr::quiet_NaN(s);
      else if (emsr::fp_is_zero(z))
	return emsr::detail::riemann_zeta(s);
      else
	return polylog_exp(s, s_i * Tp{2} * s_pi * z);
    }

  /**
   * Return the Hurwitz Zeta function for real s and complex a.
   * This uses Jonquiere's identity:
   * @f[
   *    \frac{(i2\pi)^s}{\Gamma(s)}\zeta(a,1-s) = 
   *          Li_s(e^{i2\pi a}) + (-1)^s Li_s(e^{-i2\pi a})
   * @f]
   * @param s The real argument
   * @param a The complex parameter
   */
  template<typename Tp>
    std::complex<Tp>
    hurwitz_zeta_polylog(Tp s, std::complex<Tp> a)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_2pi = emsr::tau_v<Tp>;
      const auto s_i2pi = Cmplx{0, s_2pi};
      if ((a.imag() >= Tp{0}
		&& (a.real() >= Tp{0} && a.real() <  Tp{1}))
       || (a.imag() <  Tp{0}
		&& (a.real() >  Tp{0} && a.real() <= Tp{1})))
	{
	  const auto t = Tp{1} - s;
	  const auto lpe = polylog_exp(t, s_i2pi * a);
	  /// @todo This hurwitz_zeta_polylog prefactor is prone to overflow.
	  /// positive integer orders s?
	  const auto thing = std::exp(Cmplx(Tp{0}, -s_pi_2 * t));
	  return gamma(t)
	       * std::pow(s_2pi, -t)
	       * (lpe * thing + std::conj(lpe * thing));
	}
      else
	throw std::domain_error("hurwitz_zeta_polylog: Bad argument");
    }

  /**
   * Return the Dirichlet eta function.
   * Currently, s must be real (complex type but negligible imaginary part.)
   * Otherwise std::domain_error is thrown.
   * The Dirichlet eta function, in terms of the polylogarithm, is
   * @f[
   *   \eta(s) = -\Re[Li_s(-1)]
   * @f]
   *
   * @param s  The complex (but on-real-axis) argument.
   * @return  The complex Dirichlet eta function.
   * @throw  std::domain_error if the argument has a significant imaginary part.
   */
  template<typename Tp>
    std::complex<Tp>
    dirichlet_eta(std::complex<Tp> s)
    {
      if (std::isnan(s))
	return emsr::quiet_NaN(std::imag(s));
      else if (emsr::fp_is_real(s))
	return -polylog(std::real(s), Tp{-1});
      else
	throw std::domain_error("dirichlet_eta: Bad argument");
    }

  /**
   * Return the Dirichlet eta function for real argument.
   * The Dirichlet eta function, in terms of the polylogarithm, is
   * @f[
   *   \eta(s) = -\Re[Li_s(-1)]
   * @f]
   *
   * @param s  The real argument.
   * @return  The Dirichlet eta function.
   */
  template<typename Tp>
    Tp
    dirichlet_eta(Tp s)
    {
      if (std::isnan(s))
	return emsr::quiet_NaN(s);
      else if (s < Tp{0})
	{
	  const auto p = emsr::fp_is_integer(s, Tp{5});
	  if (p && (p() % 2 == 0))
	    return Tp{0};
	  else
	    {
	      const auto s_pi = emsr::pi_v<Tp>;
	      const auto sc = Tp{1} - s;
	      const auto p2 = std::pow(Tp{2}, -sc);
	      return Tp{2} * (Tp{1} - p2) / (Tp{1} - Tp{2} * p2)
		   * std::pow(s_pi, -sc) * s * emsr::sin_pi(s / Tp{2})
		   * gamma(-s) * dirichlet_eta(sc);
	    }
	}
      else
	return -std::real(polylog(s, Tp{-1}));
    }

  /**
   * Return the Dirichlet beta function.
   * Currently, s must be real (complex type but negligible imaginary part.)
   * Otherwise std::domain_error is thrown.
   * The Dirichlet beta function, in terms of the polylogarithm, is
   * @f[
   *   \beta(s) = \Im[Li_s(i)]
   * @f]
   *
   * @param s  The complex (but on-real-axis) argument.
   * @return  The Dirichlet Beta function of real argument.
   * @throw  std::domain_error if the argument has a significant imaginary part.
   */
  template<typename Tp>
    Tp
    dirichlet_beta(std::complex<Tp> s)
    {
      const auto s_i = std::complex<Tp>{0, 1};
      if (std::isnan(s))
	return emsr::quiet_NaN(std::imag(s));
      else if (emsr::fp_is_real(s))
	return std::imag(polylog(s.real(), s_i));
      else
	throw std::domain_error("dirichlet_beta: Bad argument.");
    }

  /**
   * Return the Dirichlet beta function for real argument.
   * The Dirichlet beta function, in terms of the polylogarithm, is
   * @f[
   *   \beta(s) = \Im[Li_s(i)]
   * @f]
   *
   * @param s  The real argument.
   * @return  The Dirichlet Beta function of real argument.
   */
  template<typename Tp>
    Tp
    dirichlet_beta(Tp s)
    {
      const auto s_i = std::complex<Tp>{0, 1};
      if (std::isnan(s))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	return std::imag(polylog(s, s_i));
    }

  /**
   * Return the Dirichlet lambda function for real argument.
   * @f[
   *   \lambda(s) = \frac{1}{2}(\zeta(s) + \eta(s))
   * @f]
   *
   * @param s  The real argument.
   * @return  The Dirichlet lambda function.
   */
  template<typename Tp>
    Tp
    dirichlet_lambda(Tp s)
    {
      if (std::isnan(s))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	return (emsr::detail::riemann_zeta(s) + dirichlet_eta(s)) / Tp{2};
    }

  /**
   * Return Clausen's function of integer order m and complex argument @c z.
   * The notation and connection to polylog is from Wikipedia
   *
   * @param m  The non-negative integral order.
   * @param z  The complex argument.
   * @return  The complex Clausen function.
   */
  template<typename Tp>
    std::complex<Tp>
    clausen(unsigned int m, std::complex<Tp> z)
    {
      if (std::isnan(z))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (m == 0)
	throw std::domain_error("clausen: Non-positive order");
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  const auto ple = polylog_exp(Tp(m), s_i * z);
	  if (m & 1)
	    return ple;
	  else
	    return s_i * std::conj(ple);
	}
    }

  /**
   * Return Clausen's function of integer order m and real argument x.
   * The notation and connection to polylog is from Wikipedia
   *
   * @param m  The integer order m >= 1.
   * @param x  The real argument.
   * @return  The Clausen function.
   */
  template<typename Tp>
    Tp
    clausen(unsigned int m, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (m == 0)
	throw std::domain_error("clausen: Non-positive order");
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  const auto ple = polylog_exp(Tp(m), s_i * x);
	  if (m & 1)
	    return std::real(ple);
	  else
	    return std::imag(ple);
	}
    }

  /**
   * Return Clausen's sine sum Sl_m for positive integer order m
   * and complex argument z.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param m  The integer order m >= 1.
   * @param z  The complex argument.
   * @return  The Clausen sine sum Sl_m(w),
   */
  template<typename Tp>
    Tp
    clausen_sl(unsigned int m, std::complex<Tp> z)
    {
      if (std::isnan(z))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (m == 0)
	throw std::domain_error("clausen_sl: Non-positive order");
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  const auto ple = polylog_exp(Tp(m), s_i * z);
	  if (m & 1)
	    return std::imag(ple);
	  else
	    return std::real(ple);
	}
    }

  /**
   * Return Clausen's sine sum Sl_m for positive integer order m
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param m  The integer order m >= 1.
   * @param x  The real argument.
   * @return  The Clausen sine sum Sl_m(w),
   */
  template<typename Tp>
    Tp
    clausen_sl(unsigned int m, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (m == 0)
	throw std::domain_error("clausen_sl: Non-positive order");
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  const auto ple = polylog_exp(Tp(m), s_i * x);
	  if (m & 1)
	    return std::imag(ple);
	  else
	    return std::real(ple);
	}
    }

  /**
   * Return Clausen's cosine sum Cl_m for positive integer order m
   * and complex argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param m  The integer order m >= 1.
   * @param z  The complex argument.
   * @return  The Clausen cosine sum Cl_m(w),
   */
  template<typename Tp>
    Tp
    clausen_cl(unsigned int m, std::complex<Tp> z)
    {
      if (std::isnan(z))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (m == 0)
	throw std::domain_error("clausen_cl: Non-positive order");
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  const auto ple = polylog_exp(Tp(m), s_i * z);
	  if (m & 1)
	    return std::real(ple);
	  else
	    return std::imag(ple);
	}
    }

  /**
   * Return Clausen's cosine sum Cl_m for positive integer order m
   * and real argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param m  The integer order m >= 1.
   * @param x  The real argument.
   * @return  The real Clausen cosine sum Cl_m(w),
   */
  template<typename Tp>
    Tp
    clausen_cl(unsigned int m, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (m == 0)
	throw std::domain_error("clausen_cl: Non-positive order");
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  const auto ple = polylog_exp(Tp(m), s_i * x);
	  if (m & 1)
	    return std::real(ple);
	  else
	    return std::imag(ple);
	}
    }

  /**
   * Return the Fermi-Dirac integral of integer or real order s
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   * @see http://dlmf.nist.gov/25.12.16
   *
   * @f[
   *    F_s(x) = \frac{1}{\Gamma(s+1)}\int_0^\infty \frac{t^s}{e^{t-x} + 1}dt
   *           = -Li_{s+1}(-e^x)
   * @f]
   *
   * @param s  The order s > -1.
   * @param x  The real argument.
   * @return  The real Fermi-Dirac integral F_s(x),
   */
  template<typename Sp, typename Tp>
    Tp
    fermi_dirac(Sp s, Tp x)
    {
      if (std::isnan(s) || std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (s <= Sp{-1})
	throw std::domain_error("fermi_dirac: Order must be greater than -1");
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  const auto s_pi = emsr::pi_v<Tp>;
	  return -std::real(polylog_exp(s + Sp{1}, x + s_i * s_pi));
	}
    }

  /**
   * Return the Bose-Einstein integral of integer or real order s
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   * @see http://dlmf.nist.gov/25.12.16
   *
   * @f[
   *    G_s(x) = \frac{1}{\Gamma(s+1)}\int_0^\infty \frac{t^s}{e^{t-x} - 1}dt
   *           = Li_{s+1}(e^x)
   * @f]
   *
   * @param s  The order s >= 0.
   * @param x  The real argument.
   * @return  The real Bose-Einstein integral G_s(x),
   */
  template<typename Sp, typename Tp>
    Tp
    bose_einstein(Sp s, Tp x)
    {
      if (std::isnan(s) || std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (s <= Sp{0} && x < Tp{0})
	throw std::domain_error("bose_einstein: Order must be greater than 0");
      else
	return std::real(polylog_exp(s + Sp{1}, x));
    }

} // namespace detail
} // namespace emsr

#endif // SF_POLYLOG_TCC
