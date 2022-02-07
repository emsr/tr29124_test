
// Copyright (C) 2006-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_hyperg.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1) Handbook of Mathematical Functions,
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 6, pp. 555-566
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl

#ifndef SF_HYPERG_TCC
#define SF_HYPERG_TCC 1

#include <stdexcept>

#include <emsr/notsospecfun.h>
#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>
#include <emsr/math_util.h>
#include <emsr/sf_gamma.h>

namespace emsr
{
namespace detail
{

/*
  template<typename Tp, typename Iter>
    Tp
    gamma_ratio(Iter a, Iter b)
    {
      auto sign_F = Tp{1};
      auto F = Tp{0};
      for (auto bb : b)
	{
	  bool ok = true;
	  const auto int_b  = std::floor(bb + Tp{0.5L});
	  const auto b_integer = std::abs(bb - intd) < toler;
	  if (b_integer && int_b <= Tp{0})
	    return Tp{0};
	  else
	    {
	      try
		{
		  sign_F *= emsr::detail::log_gamma_sign(bb);
		  F -= emsr::detail::log_gamma(bb);
		}
	      catch (...)
		{
		  ok = false;
		  break;
		}
	    }
	}
      for (auto aa : a)
	{
	  const auto int_a  = std::floor(aa + Tp{0.5L});
	  const auto a_integer = std::abs(aa - intd) < toler;
	  if (a_integer && int_a <= Tp{0})
	    return s_inf;
	  else
	    {
	      try
		{
		  sign_F *= emsr::detail::log_gamma_sign(aa);
		  F += emsr::detail::log_gamma(aa);
		}
	      catch (...)
		{
		  ok = false;
		  break;
		}
	    }
	}
      if (F > s_log_max)
	throw std::runtime_error("gamma_ratio: overflow of gamma function ratios");
    }
*/
  /**
   * @brief This routine returns the confluent hypergeometric limit function
   * 	    by series expansion.
   *
   * @f[
   *   {}_0F_1(-;c;x) = \Gamma(c)
   * 		\sum_{n=0}^{\infty} \frac{1}{\Gamma(c+n)} \frac{x^n}{n!}
   * @f]
   *
   * If a and b are integers and a < 0 and either b > 0 or b < a
   * then the series is a polynomial with a finite number of
   * terms.
   *
   * @param  c  The "denominator" parameter.
   * @param  x  The argument of the confluent hypergeometric limit function.
   * @return  The confluent hypergeometric limit function.
   */
  template<typename Tp>
    Tp
    conf_hyperg_lim_series(Tp c, Tp x)
    {
      const auto eps = emsr::epsilon(x);

      auto term = Tp{1};
      auto Fac = Tp{1};
      const unsigned int max_iter = 100000;
      unsigned int i;
      for (i = 0; i < max_iter; ++i)
	{
	  term *= x / ((c + Tp(i)) * Tp(1 + i));
	  Fac += term;
	  if (std::abs(term) < eps)
	    break;
	}
      if (i == max_iter)
	throw std::runtime_error("conf_hyperg_lim_series: series failed to converge");

      return Fac;
    }


  /**
   * @brief  Return the confluent hypergeometric limit function
   *	     @f$ {}_0F_1(-;c;x) @f$.
   *
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric limit function.
   * @return  The confluent limit hypergeometric function.
   */
  template<typename Tp>
    Tp
    conf_hyperg_lim(Tp c, Tp x)
    {
      const auto c_nint = emsr::fp_is_integer(c);
      if (std::isnan(c) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (c_nint && c_nint() <= 0)
	return emsr::infinity(x);
      //else if (x < Tp{0})
	//return conf_hyperg_lim_luke(c, x);
      else
	return conf_hyperg_lim_series(c, x);
    }


  /**
   * @brief This routine returns the confluent hypergeometric function
   * 	    by series expansion.
   *
   * @f[
   *   {}_1F_1(a;c;x) = \frac{\Gamma(c)}{\Gamma(a)}
   * 			\sum_{n=0}^{\infty}
   * 			\frac{\Gamma(a+n)}{\Gamma(c+n)}
   * 			\frac{x^n}{n!}
   * @f]
   *
   * @param  a  The "numerator" parameter.
   * @param  c  The "denominator" parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The confluent hypergeometric function.
   */
  template<typename Tp>
    Tp
    conf_hyperg_series(Tp a, Tp c, Tp x)
    {
      const auto eps = emsr::epsilon(x);

      auto term = Tp{1};
      auto Fac = Tp{1};
      const unsigned int max_iter = 100000;
      unsigned int i;
      for (i = 0; i < max_iter; ++i)
	{
	  term *= (a + Tp(i)) * x
		  / ((c + Tp(i)) * Tp(1 + i));
	  Fac += term;
	  if (std::abs(term) < eps)
	    break;
	}
      if (i == max_iter)
	throw std::runtime_error("conf_hyperg_series: series failed to converge");

      return Fac;
    }


  /**
   * @brief  Return the hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   *	     by an iterative procedure described in
   *	     Luke, Algorithms for the Computation of Mathematical Functions.
   *
   * Like the case of the 2F1 rational approximations, these are
   * probably guaranteed to converge for x < 0, barring gross
   * numerical instability in the pre-asymptotic regime.
   */
  template<typename Tp>
    Tp
    conf_hyperg_luke(Tp a, Tp c, Tp xin)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto big = emsr::root_max(Val{6});
      const int nmax = 20000;
      const auto eps = emsr::epsilon<Val>();
      const auto x  = -xin;
      const auto x3 = x * x * x;
      const auto t0 = a / c;
      const auto t1 = (a + Val{1}) / (Val{2} * c);
      const auto t2 = (a + Val{2}) / (Val{2} * (c + Val{1}));

      auto F = Tp{1};

      auto Bnm3 = Tp{1};
      auto Bnm2 = Tp{1} + t1 * x;
      auto Bnm1 = Tp{1} + t2 * x * (Val{1} + t1 / Val{3} * x);

      auto Anm3 = Tp{1};
      auto Anm2 = Bnm2 - t0 * x;
      auto Anm1 = Bnm1 - t0 * (Tp{1} + t2 * x) * x
		  + t0 * t1 * (c / (c + Val{1})) * x * x;

      int n = 3;
      while(true)
	{
	  auto npam1 = Val(n - 1) + a;
	  auto npcm1 = Val(n - 1) + c;
	  auto npam2 = Val(n - 2) + a;
	  auto npcm2 = Val(n - 2) + c;
	  auto tnm1  = Val(2 * n - 1);
	  auto tnm3  = Val(2 * n - 3);
	  auto tnm5  = Val(2 * n - 5);
	  auto F1 =  (Val(n - 2) - a) / (Val(2 * tnm3) * npcm1);
	  auto F2 =  (Val(n) + a) * npam1
		    / (Val(4 * tnm1 * tnm3) * npcm2 * npcm1);
	  auto F3 = -npam2 * npam1 * (Val(n - 2) - a)
		    / (Val(8 * tnm3 * tnm3 * tnm5)
		    * (Val(n - 3) + c) * npcm2 * npcm1);
	  auto E  = -npam1 * (Val(n - 1) - c)
		    / (Val(2 * tnm3) * npcm2 * npcm1);

	  auto An = (Val{1} + F1 * x) * Anm1
		    + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
	  auto Bn = (Val{1} + F1 * x) * Bnm1
		    + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
	  auto r = An / Bn;

	  const auto prec = std::abs((F - r) / F);
	  F = r;

	  if (prec < eps || n > nmax)
	    break;

	  if (std::abs(An) > big || std::abs(Bn) > big)
	    {
	      An   /= big;
	      Bn   /= big;
	      Anm1 /= big;
	      Bnm1 /= big;
	      Anm2 /= big;
	      Bnm2 /= big;
	      Anm3 /= big;
	      Bnm3 /= big;
	    }
	  else if (std::abs(An) < Tp{1} / big
		|| std::abs(Bn) < Tp{1} / big)
	    {
	      An   *= big;
	      Bn   *= big;
	      Anm1 *= big;
	      Bnm1 *= big;
	      Anm2 *= big;
	      Bnm2 *= big;
	      Anm3 *= big;
	      Bnm3 *= big;
	    }

	  ++n;
	  Bnm3 = Bnm2;
	  Bnm2 = Bnm1;
	  Bnm1 = Bn;
	  Anm3 = Anm2;
	  Anm2 = Anm1;
	  Anm1 = An;
	}

      if (n >= nmax)
	throw std::runtime_error("conf_hyperg_luke: iteration failed to converge");

      return F;
    }


  /**
   * @brief  Return the confluent hypergeometric function
   * 	     @f$ {}_1F_1(a;c;x) = M(a,c,x) @f$.
   *
   * @param  a  The @a numerator parameter.
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The confluent hypergeometric function.
   */
  template<typename Tp>
    Tp
    conf_hyperg(Tp a, Tp c, Tp x)
    {
      const auto c_nint = emsr::fp_is_integer(c);
      if (std::isnan(a) || std::isnan(c) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (c_nint && c_nint() <= 0)
	return emsr::infinity(x);
      else if (a == Tp{0})
	return Tp{1};
      else if (c == a)
	return std::exp(x);
      else if (x < Tp{0})
	return conf_hyperg_luke(a, c, x);
      else
	return conf_hyperg_series(a, c, x);
    }


  /**
   * @brief  Return the Tricomi confluent hypergeometric function
   * @f[
   *   U(a,c,x) = \frac{\Gamma(1-c)}{\Gamma(a-c+1)} {}_1F_1(a;c;x)
   *       + \frac{\Gamma(c-1)}{\Gamma(a)} x^{1-c} {}_1F_1(a-c+1;2-c;x)
   * @f]
   * @param  a  The @a numerator parameter.
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The Tricomi confluent hypergeometric function.
   */
  template<typename Tp>
    Tp
    tricomi_u_naive(Tp a, Tp c, Tp x)
    {
      auto U1 = Tp{};
      auto b = a - c + Tp{1};
      auto ib = emsr::fp_is_integer(b);
      if (!ib || (ib && ib() > 0))
	U1 = std::tgamma(Tp{1} - c)
	       * conf_hyperg(a, c, x)
	       / std::tgamma(b);

      auto U2 = Tp{};
      auto ia = emsr::fp_is_integer(a);
      if (!ia || (ia && ia() > 0))
	U2 = std::tgamma(c - Tp{1})
	       * std::pow(x, Tp{1} - c)
	       * conf_hyperg(b, Tp{2} - c, x)
	       / std::tgamma(a);

      return U1 + U2;
    }

  /**
   * @brief  Return the Tricomi confluent hypergeometric function
   * @f[
   *   U(a,c,x) = \frac{\Gamma(1-c)}{\Gamma(a-c+1)} {}_1F_1(a;c;x)
   *       + \frac{\Gamma(c-1)}{\Gamma(a)} x^{1-c} {}_1F_1(a-c+1;2-c;x)
   * @f]
   * @param  a  The @a numerator parameter.
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The Tricomi confluent hypergeometric function.
   */
  template<typename Tp>
    Tp
    tricomi_u(Tp a, Tp c, Tp x)
    {
      return tricomi_u_naive(a, c, x);
    }


  /**
   * @brief Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * by series expansion.
   *
   * The hypergeometric function is defined by
   * @f[
   *   {}_2F_1(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
   * 			\sum_{n=0}^{\infty}
   * 			\frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)}
   * 			\frac{x^n}{n!}
   * @f]
   *
   * This works and it's pretty fast.
   *
   * @param  a  The first @a numerator parameter.
   * @param  b  The second @a numerator parameter.
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The confluent hypergeometric function.
   */
  template<typename Tp>
    Tp
    hyperg_series(Tp a, Tp b, Tp c, Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto eps = emsr::epsilon<Val>();

      auto term = Tp{1};
      auto Fabc = Tp{1};
      const unsigned int max_iter = 100000;
      unsigned int i;
      for (i = 0u; i < max_iter; ++i)
	{
	  term *= (a + Tp(i)) * (b + Tp(i)) * x
		  / ((c + Tp(i)) * Tp(1 + i));
	  Fabc += term;
	  if (std::abs(term) < eps)
	    break;
	}
      if (i == max_iter)
	throw std::runtime_error("Series failed to converge in hyperg_series.");

      return Fabc;
    }


  /**
   * @brief Return the hypergeometric polynomial @f$ {}_2F_1(-m,b;c;x) @f$
   * by Holm recursion.
   *
   * The hypergeometric function is defined by
   * @f[
   *   {}_2F_1(-m,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
   * 			\sum_{n=0}^{\infty}
   * 			\frac{\Gamma(n-m)\Gamma(b+n)}{\Gamma(c+n)}
   * 			\frac{x^n}{n!}
   * @f]
   *
   * @f[
   * @f]
   *
   * @param  m  The first @a numerator parameter.
   * @param  b  The second @a numerator parameter.
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The confluent hypergeometric function.
   */
  template<typename Tp>
    Tp
    hyperg_recur(int m, Tp b, Tp c, Tp x)
    {
      auto term = Tp{1};
      auto Fabc = Tp{1};
      for (int i = 0; i < -m; ++i)
	{
	  term *= Tp(m + i) * (b + Tp(i)) * x
		  / ((c + Tp(i)) * Tp(1 + i));
	  Fabc += term;
/// @fixme: go recur!
	}

      return Fabc;
    }


  /**
   * @brief  Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * 	     by an iterative procedure described in
   * 	     Luke, Algorithms for the Computation of Mathematical Functions.
   */
  template<typename Tp>
    Tp
    hyperg_luke(Tp a, Tp b, Tp c, Tp xin)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto big = emsr::root_max(Val{6});
      const int nmax = 20000;
      const auto eps = emsr::epsilon<Val>();
      const auto x  = -xin;
      const auto x3 = x * x * x;
      const auto t0 = a * b / c;
      const auto t1 = (a + Tp{1}) * (b + Tp{1}) / (Tp{2} * c);
      const auto t2 = (a + Tp{2}) * (b + Tp{2})
		     / (Tp{2} * (c + Tp{1}));

      auto F = Tp{1};

      auto Bnm3 = Tp{1};
      auto Bnm2 = Tp{1} + t1 * x;
      auto Bnm1 = Tp{1} + t2 * x * (Tp{1} + t1 / Tp{3} * x);

      auto Anm3 = Tp{1};
      auto Anm2 = Bnm2 - t0 * x;
      auto Anm1 = Bnm1 - t0 * (Tp{1} + t2 * x) * x
		  + t0 * t1 * (c / (c + Tp{1})) * x * x;

      int n = 3;
      while (true)
	{
	  const auto npam1 = Val(n - 1) + a;
	  const auto npbm1 = Val(n - 1) + b;
	  const auto npcm1 = Val(n - 1) + c;
	  const auto npam2 = npam1 - Val{1};
	  const auto npbm2 = npbm1 - Val{1};
	  const auto npcm2 = npcm1 - Val{1};
	  const auto tnm1  = Val(2 * n - 1);
	  const auto tnm3  = tnm1 - Val{2};
	  const auto tnm5  = tnm3 - Val{2};
	  const auto n2 = n * n;
	  const auto F1 = (Val(3 * n2) + (a + b - Val{6}) * Val(n)
			    + Val{2} - a * b - Val{2} * (a + b))
			  / (Val(2 * tnm3) * npcm1);
	  const auto F2 = -(Val(3 * n2) - (a + b + Val{6}) * Val(n)
			    + Val{2} - a * b) * npam1 * npbm1
			  / (Val(4 * tnm1 * tnm3) * npcm2 * npcm1);
	  const auto F3 = (npam2 * npam1 * npbm2 * npbm1
			  * (Val(n - 2) - a) * (Val(n - 2) - b))
			  / (Val(8 * tnm3 * tnm3 * tnm5)
			  * (Val(n - 3) + c) * npcm2 * npcm1);
	  const auto E  = -npam1 * npbm1 * (Val(n - 1) - c)
			  / (Val(2 * tnm3) * npcm2 * npcm1);

	  auto An = (Val{1} + F1 * x) * Anm1
		    + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
	  auto Bn = (Val{1} + F1 * x) * Bnm1
		    + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
	  const auto r = An / Bn;

	  const auto prec = std::abs((F - r) / F);
	  F = r;

	  if (prec < eps || n > nmax)
	    break;

	  if (std::abs(An) > big || std::abs(Bn) > big)
	    {
	      An   /= big;
	      Bn   /= big;
	      Anm1 /= big;
	      Bnm1 /= big;
	      Anm2 /= big;
	      Bnm2 /= big;
	      Anm3 /= big;
	      Bnm3 /= big;
	    }
	  else if (std::abs(An) < Val{1} / big
		|| std::abs(Bn) < Val{1} / big)
	    {
	      An   *= big;
	      Bn   *= big;
	      Anm1 *= big;
	      Bnm1 *= big;
	      Anm2 *= big;
	      Bnm2 *= big;
	      Anm3 *= big;
	      Bnm3 *= big;
	    }

	  ++n;
	  Bnm3 = Bnm2;
	  Bnm2 = Bnm1;
	  Bnm1 = Bn;
	  Anm3 = Anm2;
	  Anm2 = Anm1;
	  Anm1 = An;
	}

      if (n >= nmax)
	throw std::runtime_error("hyperg_luke: iteration failed to converge");

      return F;
    }


  /**
   * @brief  Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * by the reflection formulae in Abramowitz & Stegun formula
   * 15.3.6 for d = c - a - b not integral and formula 15.3.11 for
   * d = c - a - b integral.  This assumes a, b, c != negative
   * integer.
   *
   * The hypergeometric function is defined by
   * @f[
   *   {}_2F_1(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
   *    		\sum_{n=0}^{\infty}
   *    		\frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)}
   *    		\frac{x^n}{n!}
   * @f]
   *
   * The reflection formula for nonintegral @f$ d = c - a - b @f$ is:
   * @f[
   *   {}_2F_1(a,b;c;x) = \frac{\Gamma(c)\Gamma(d)}{\Gamma(c-a)\Gamma(c-b)}
   *    		      {}_2F_1(a,b;1-d;1-x)
   *    	      + \frac{\Gamma(c)\Gamma(-d)}{\Gamma(a)\Gamma(b)}
   *    		      {}_2F_1(c-a,c-b;1+d;1-x)
   * @f]
   *
   * The reflection formula for integral @f$ m = c - a - b @f$ is:
   * @f[
   *   {}_2F_1(a,b;a+b+m;x)
   *        = \frac{\Gamma(m)\Gamma(a+b+m)}{\Gamma(a+m)\Gamma(b+m)}
   *          \sum_{k=0}^{m-1} \frac{(m+a)_k(m+b)_k}{k!(1-m)_k} (1 - x)^k
   *    		 + (-1)^m 
   * @f]
   */
  template<typename Tp>
    Tp
    hyperg_reflect(Tp a, Tp b, Tp c, Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto s_log_max = emsr::log_max<Val>();
      const auto eps = emsr::epsilon<Val>();
      const auto epsfact = Val{1000};
      const auto d = c - a - b;
      const auto d_nint = emsr::fp_is_integer(d, epsfact);

      if (d_nint)
	{
	  const auto ln_omx = emsr::log1p(-x);
	  const auto ad = std::abs(d);
	  Tp F1, F2;

	  Tp d1, d2;
	  if (std::real(d) >= Val{0})
	    {
	      d1 = d;
	      d2 = Tp{0};
	    }
	  else
	    {
	      d1 = Tp{0};
	      d2 = d;
	    }

	  const auto lng_c = emsr::detail::log_gamma(c);

	  if (ad < eps)
	    {
	      // d = c - a - b = 0.
	      F1 = Tp{0};
	    }
	  else
	    {

	      bool ok_d1 = true;
	      Tp lng_ad, lng_ad1, lng_bd1;
	      Tp sgn_ad, sgn_ad1, sgn_bd1;
	      try
		{
		  sgn_ad = emsr::detail::log_gamma_sign(ad);
		  lng_ad = emsr::detail::log_gamma(ad);
		  sgn_ad1 = emsr::detail::log_gamma_sign(a + d1);
		  lng_ad1 = emsr::detail::log_gamma(a + d1);
		  sgn_bd1 = emsr::detail::log_gamma_sign(b + d1);
		  lng_bd1 = emsr::detail::log_gamma(b + d1);
		}
	      catch(...)
		{
		  ok_d1 = false;
		}

	      if (ok_d1)
		{
		  /* Gamma functions in the denominator are ok.
		   * Proceed with evaluation.
		   */
		  auto sum1 = Tp{1};
		  auto term = Tp{1};
		  auto ln_pre1 = lng_ad + lng_c + d2 * ln_omx
				 - lng_ad1 - lng_bd1;

		  if (std::abs(ln_pre1) > s_log_max)
		    throw std::runtime_error("hyperg_reflect: overflow of gamma functions");

		  /* Do F1 sum.
		   */
		  for (int i = 1; i < ad; ++i)
		    {
		      const int j = i - 1;
		      term *= (a + d2 + Val(j)) * (b + d2 + Val(j))
			     / (Val{1} + d2 + Val(j))
			     / Val(i) * (Val{1} - x);
		      sum1 += term;
		    }

		  if (std::abs(ln_pre1) > s_log_max)
		    throw std::runtime_error("hyperg_reflect: overflow of gamma functions");
		  else
		    F1 = sgn_ad * sgn_ad1 * sgn_bd1
			 * std::exp(ln_pre1) * sum1;
		}
	      else
		{
		  // Gamma functions in the denominator were not ok (they diverged).
		  // So the F1 term is zero.
		  F1 = Tp{0};
		}
	    }

	  bool ok_d2 = true;
	  Tp lng_ad2, lng_bd2;
	  Tp sgn_ad2, sgn_bd2;
	  try
	    {
	      sgn_ad2 = emsr::detail::log_gamma_sign(a + d2);
	      lng_ad2 = emsr::detail::log_gamma(a + d2);
	      sgn_bd2 = emsr::detail::log_gamma_sign(b + d2);
	      lng_bd2 = emsr::detail::log_gamma(b + d2);
	    }
	  catch(...)
	    {
	      ok_d2 = false;
	    }

	  if (ok_d2)
	    {
	      // Gamma functions in the denominator are ok.
	      // Proceed with evaluation.
	      const int maxiter = 2000;
	      const auto psi_1 = -emsr::egamma_v<Val>;
	      const auto psi_1pd = emsr::digamma(Val{1} + ad);
	      const auto psi_apd1 = emsr::digamma(a + d1);
	      const auto psi_bpd1 = emsr::digamma(b + d1);

	      auto psi_term = psi_1 + psi_1pd - psi_apd1
			      - psi_bpd1 - ln_omx;
	      auto fact = Tp{1};
	      auto sum2 = psi_term;
	      auto ln_pre2 = lng_c + d1 * ln_omx
			     - lng_ad2 - lng_bd2;

	      if (std::abs(ln_pre2) > s_log_max)
		throw std::runtime_error("hyperg_reflect: overflow of gamma functions");

	      int j;
	      for (j = 1; j < maxiter; ++j)
		{
		  // Values for psi functions use recurrence;
		  // Abramowitz & Stegun 6.3.5
		  const auto term1 = Tp{1} / Tp(j)
				     + Tp{1} / (ad + j);
		  const auto term2 = Tp{1} / (a + d1 + Tp(j - 1))
				     + Tp{1} / (b + d1 + Tp(j - 1));
		  psi_term += term1 - term2;
		  fact *= (a + d1 + Tp(j - 1))
			  * (b + d1 + Tp(j - 1))
			  / ((ad + j) * j) * (Tp{1} - x);
		  const auto delta = fact * psi_term;
		  sum2 += delta;
		  if (std::abs(delta) < eps * std::abs(sum2))
		    break;
		}
	      if (j == maxiter)
		throw std::runtime_error("hyperg_reflect: sum F2 failed to converge");

	      if (sum2 == Tp{0})
		F2 = Tp{0};
	      else
		F2 = sgn_ad2 * sgn_bd2 * std::exp(ln_pre2) * sum2;
	    }
	  else
	    {
	      // Gamma functions in the denominator not ok (they diverge).
	      // So the F2 term is zero.
	      F2 = Tp{0};
	    }

	  const auto sgn_2 = (d_nint() % 2 == 1 ? -Tp{1} : Tp{1});
	  const auto F = F1 + sgn_2 * F2;

	  return F;
	}
      else // d = c - a - b not an integer.
	{
	  // These gamma functions appear in the denominator, so we
	  // catch their harmless domain errors and set the terms to zero.
	  bool ok1 = true;
	  auto sgn_g1ca = Tp{0}, ln_g1ca = Tp{0};
	  auto sgn_g1cb = Tp{0}, ln_g1cb = Tp{0};
	  try
	    {
	      sgn_g1ca = emsr::detail::log_gamma_sign(c - a);
	      ln_g1ca = emsr::detail::log_gamma(c - a);
	      sgn_g1cb = emsr::detail::log_gamma_sign(c - b);
	      ln_g1cb = emsr::detail::log_gamma(c - b);
	    }
	  catch(...)
	    {
	      ok1 = false;
	    }

	  bool ok2 = true;
	  auto sgn_g2a = Tp{0}, ln_g2a = Tp{0};
	  auto sgn_g2b = Tp{0}, ln_g2b = Tp{0};
	  try
	    {
	      sgn_g2a = emsr::detail::log_gamma_sign(a);
	      ln_g2a = emsr::detail::log_gamma(a);
	      sgn_g2b = emsr::detail::log_gamma_sign(b);
	      ln_g2b = emsr::detail::log_gamma(b);
	    }
	  catch(...)
	    {
	      ok2 = false;
	    }

	  const auto sgn_gc = emsr::detail::log_gamma_sign(c);
	  const auto ln_gc = emsr::detail::log_gamma(c);
	  const auto sgn_gd = emsr::detail::log_gamma_sign(d);
	  const auto ln_gd = emsr::detail::log_gamma(d);
	  const auto sgn_gmd = emsr::detail::log_gamma_sign(-d);
	  const auto ln_gmd = emsr::detail::log_gamma(-d);

	  const auto sgn1 = sgn_gc * sgn_gd  * sgn_g1ca * sgn_g1cb;
	  const auto sgn2 = sgn_gc * sgn_gmd * sgn_g2a  * sgn_g2b;

	  Tp pre1, pre2;
	  if (ok1 && ok2)
	    {
	      auto ln_pre1 = ln_gc + ln_gd  - ln_g1ca - ln_g1cb;
	      auto ln_pre2 = ln_gc + ln_gmd - ln_g2a  - ln_g2b
			    + d * std::log(Tp{1} - x);
	      if (std::abs(ln_pre1) < s_log_max
		&& std::abs(ln_pre2) < s_log_max)
		{
		  pre1 = sgn1 * std::exp(ln_pre1);
		  pre2 = sgn2 * std::exp(ln_pre2);
		}
	      else
		throw std::runtime_error("hyperg_reflect: overflow of gamma functions");
	    }
	  else if (ok1 && !ok2)
	    {
	      auto ln_pre1 = ln_gc + ln_gd - ln_g1ca - ln_g1cb;
	      if (std::abs(ln_pre1) < s_log_max)
		{
		  pre1 = sgn1 * std::exp(ln_pre1);
		  pre2 = Tp{0};
		}
	      else
		throw std::runtime_error("hyperg_reflect: overflow of gamma functions");
	    }
	  else if (!ok1 && ok2)
	    {
	      auto ln_pre2 = ln_gc + ln_gmd - ln_g2a - ln_g2b
			     + d * std::log(Tp{1} - x);
	      if (std::abs(ln_pre2) < s_log_max)
		{
		  pre1 = Tp{0};
		  pre2 = sgn2 * std::exp(ln_pre2);
		}
	      else
		throw std::runtime_error("hyperg_reflect: overflow of gamma functions");
	    }
	  else
	    throw std::runtime_error("hyperg_reflect: underflow of gamma functions");

	  const auto F1 = hyperg_series(a, b, Tp{1} - d,
					    Tp{1} - x);
	  const auto F2 = hyperg_series(c - a, c - b, Tp{1} + d,
					    Tp{1} - x);

	  const auto F = pre1 * F1 + pre2 * F2;

	  return F;
	}
    }


  /**
   * @brief Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$.
   *
   * The hypergeometric function is defined by
   * @f[
   *   {}_2F_1(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
   * 			\sum_{n=0}^{\infty}
   * 			\frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)}
   * 			\frac{x^n}{n!}
   * @f]
   *
   * @param  a  The first @a numerator parameter.
   * @param  b  The second @a numerator parameter.
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The confluent hypergeometric function.
   */
  template<typename Tp>
    Tp
    hyperg(Tp a, Tp b, Tp c, Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto s_log_max = emsr::log_max<Val>();
      const auto epsfact = Val{1000};
      const auto toler = epsfact * emsr::epsilon<Val>();
      const auto a_nint = emsr::fp_is_integer(a, epsfact);
      const auto b_nint = emsr::fp_is_integer(b, epsfact);
      const auto c_nint = emsr::fp_is_integer(c, epsfact);
      if (std::abs(x - Tp{1}) < toler && std::abs(c - b - a) > Val{0}
       && !(c_nint && c_nint() <= 0))
	{
	  const auto log_gamc = emsr::detail::log_gamma(c);
	  const auto sign_gamc = emsr::detail::log_gamma_sign(c);
	  const auto log_gamcab = emsr::detail::log_gamma(c - a - b);
	  const auto log_gamca = emsr::detail::log_gamma(c - a);
	  const auto sign_gamca = emsr::detail::log_gamma_sign(c - a);
	  const auto log_gamcb = emsr::detail::log_gamma(c - b);
	  const auto sign_gamcb = emsr::detail::log_gamma_sign(c - b);
	  const auto log_pre = log_gamc + log_gamcab
			       - log_gamca - log_gamcb;
	  const auto sign = sign_gamc * sign_gamca * sign_gamcb;
	  if (sign == Val{0})
	    return emsr::quiet_NaN(x);
	  if (std::abs(log_pre) < s_log_max)
	    return sign * std::exp(log_pre);
	  else
	    throw std::domain_error("hyperg: overflow of gamma functions");
	}
      else if (std::abs(x) >= Val{1})
	throw std::domain_error("hyperg: argument outside unit circle");
      else if (std::isnan(a) || std::isnan(b)
	    || std::isnan(c) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (c_nint && c_nint() <= 0)
	return emsr::infinity(x);
      else if (std::abs(c - b) < toler || std::abs(c - a) < toler)
	return std::pow(Tp{1} - x, c - a - b);
      else if (std::real(a) >= Val{0}
	    && std::real(b) >= Val{0} && std::real(c) >= Val{0}
	    && std::real(x) >= Val{0} && std::abs(x) < Val{0.995L})
	return hyperg_series(a, b, c, x);
      else if (std::abs(a) < Val{10} && std::abs(b) < Val{10})
	{
	  // For non-positive integer a and b the hypergeometric function
	  // is a finite polynomial.
	  if (a_nint && a_nint() < 0)
	    return hyperg_series(Tp(a_nint()), b, c, x);
	  else if (b_nint && b_nint() < 0)
	    return hyperg_series(Tp(b_nint()), a, c, x);
	  else if (std::real(x) < -Val{0.25L})
	    return hyperg_luke(a, b, c, x);
	  else if (std::abs(x) < Val{0.5L})
	    return hyperg_series(a, b, c, x);
	  else if (std::abs(c) > Val{10})
	    return hyperg_series(a, b, c, x);
	  else
	    return hyperg_reflect(a, b, c, x);
	}
      else
	return hyperg_luke(a, b, c, x);
    }

} // namespace detail
} // namespace emsr

#endif // SF_HYPERG_TCC
