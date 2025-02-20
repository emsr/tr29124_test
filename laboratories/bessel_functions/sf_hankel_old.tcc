
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

/** @file emsr/sf_hankel.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef SF_HANKEL_OLD_TCC
#define SF_HANKEL_OLD_TCC 1

#include <complex>
#include <limits>
#include <vector>

#include <emsr/complex_util.h>

namespace emsr
{
namespace detail
{

  /**
   * Compute the Debye region in te complex plane.
   */
  template<typename Tp>
    void
    debye_region(std::complex<Tp> alpha, int& indexr, char& aorb)
    {
      static constexpr Tp
	s_pi(3.141592653589793238462643383279502884195e+0L);

      aorb = ' ';

      auto alphar = std::real(alpha);
      auto alphai = std::imag(alpha);

      auto f1 = Tp{1}
		- alphai * std::cos(alphai) / std::sin(alphai)
		- alphar * std::sinh(alphar) / std::cosh(alphar);

      auto f2 = Tp{1}
		+ (s_pi - alphai) * std::cos(alphai) / std::sin(alphai)
		- alphar * std::sinh(alphar) / std::cosh(alphar);

      if (f1 > Tp{0} && f2 > Tp{0})
	indexr = 1;
      else if (f2 > Tp{0})
	{
	  if (alphar > Tp{0})
	    indexr = 2;
	  else
	    indexr = 3;
	}
      else if (f1 > Tp{0})
	{
	  if (alphar > Tp{0})
	    indexr = 4;
	  else
	    indexr = 5;
	}
      else
	{
	  if (alphar > Tp{0})
            indexr = 6;
	  else
            indexr = 7;
          if (alphai <= (s_pi / Tp{2}))
            aorb = 'A';
          else
            aorb = 'B';
	}
      return;
    }


  /**
   * @brief Compute parameters depending on z and nu that appear
   * in the uniform asymptotic expansions of the Hankel functions
   * and their derivatives, except the arguments to the Airy functions.
   */
  template<typename Tp>
    void
    hankel_params(std::complex<Tp> nu, std::complex<Tp> zhat,
		    std::complex<Tp>& p, std::complex<Tp>& p2,
		    std::complex<Tp>& nup2, std::complex<Tp>& num2,
		    std::complex<Tp>& num1d3, std::complex<Tp>& num2d3,
		    std::complex<Tp>& num4d3, std::complex<Tp>& zeta,
		    std::complex<Tp>& zetaphf, std::complex<Tp>& zetamhf,
		    std::complex<Tp>& zetam3hf, std::complex<Tp>& zetrat)
    {
      using cmplx = std::complex<Tp>;

      static constexpr auto s_inf     = emsr::max<Tp>();
      static constexpr auto s_sqrt_max = emsr::sqrt_max<Tp>();

      static constexpr auto s_1d4   = Tp{0.25L};
      static constexpr auto s_1d3   = Tp{0.3333333333333333333333333333333333333333L};
      static constexpr auto s_1d2   = Tp{0.5L};
      static constexpr auto s_2d3   = Tp{0.6666666666666666666666666666666666666666L};
      static constexpr auto s_2pi   = Tp{6.283185307179586476925286766559005768391L};
      static constexpr auto s_lncon = Tp{0.2703100720721095879853420769762327577152L}; // -(2/3)ln(2/3)
      static constexpr auto s_sqrt2 = Tp{1.414213562373095048801688724209698078569L};
      static constexpr auto s_4d3   = Tp{1.333333333333333333333333333333333333333L};

      static constexpr cmplx zone{Tp{1}, Tp{0}};
      static constexpr cmplx s_j{Tp{0}, Tp{1}};

      // Separate real and imaginary parts of zhat.
      auto rezhat = std::real(zhat);
      auto imzhat = std::imag(zhat);
      auto abs_rezhat = std::abs(rezhat);
      auto abs_imzhat = std::abs(imzhat);

      // If 1 - zhat^2 can be computed without overflow.
      if (abs_rezhat < s_sqrt_max && abs_imzhat < s_sqrt_max)
	{
	  // Find max and min of abs(dx) and abs(dy).
	  auto du = abs_rezhat;
	  auto dv = abs_imzhat;
	  if (du < dv)
	    std::swap(du, dv);
	  if (du >= s_1d2 && dv > s_inf / (Tp{2} * du))
	    throw std::runtime_error("hankel_params: unable to compute 1-zhat^2");
	}
      else
	throw std::runtime_error("hankel_params: unable to compute 1-zhat^2");

      // Compute 1 - zhat^2 and related constants.
      auto w = cmplx{Tp{1} - (rezhat - imzhat) * (rezhat + imzhat),
			-Tp{2} * rezhat * imzhat};
      w = std::sqrt(w);
      p = Tp{1} / w;
      p2 = p * p;

      // If nu^2 can be computed without overflow.
      if (std::abs(nu) <= s_sqrt_max)
	{
	  nup2 = nu * nu;
	  num2 = Tp{1} / nup2;
	  // Compute nu^(-1/3), nu^(-2/3), nu^(-4/3).
	  num4d3 = -std::log(nu);
	  num1d3 = std::exp(s_1d3 * num4d3);
	  num2d3 = std::exp(s_2d3 * num4d3);
	  num4d3 = std::exp(s_4d3 * num4d3);
	}
      else
	throw std::runtime_error("hankel_params: unable to compute nu^2");

      // Compute xi = ln(1+(1-zhat^2)^(1/2)) - ln(zhat) - (1-zhat^2)^(1/2)
      // using default branch of logarithm and square root.
      auto xi = std::log(zone + w) - std::log(zhat) - w;
      zetam3hf = s_2d3 / xi;

      // Compute principal value of ln(xi) and then adjust imaginary part.
      auto lnxi = std::log(xi);

      // Prepare to adjust logarithm of xi to appropriate Riemann sheet.
      auto npi = Tp{0};

      // Find adjustment necessary to get on proper Riemann sheet.
      if (imzhat == Tp{0})  // zhat is real.
	{
	  if (rezhat > Tp{1})
	    npi = s_2pi;
	}
      else // zhat is not real.
	{
	  // zhat is in upper half-plane.
	  if (imzhat > Tp{0})
	    {
	      // xi lies in upper half-plane.
	      if (std::imag(xi) > Tp{0})
		npi = -s_2pi;
	      else
		npi = +s_2pi;
	    }
	}

      // Adjust logarithm of xi.
      lnxi += npi * s_j;

      // Compute ln(zeta), zeta, zeta^(+1/2), zeta^(-1/2).
      auto lnzeta = s_2d3 * lnxi + s_lncon;
      zeta = std::exp(lnzeta);
      zetaphf = std::sqrt(zeta);
      zetamhf = Tp{1} / zetaphf;

      // Compute (4 * zeta / (1 - zhat^2))^(1/4).
      w = std::log(w);
      zetrat = s_sqrt2 * std::exp(s_1d4 * lnzeta - s_1d2 * w);

      return;
    }


  /**
   * @brief Compute the arguments for the Airy function evaluations
   * carefully to prevent premature overflow.  Note that the
   * major work here is in @c safe_div.  A faster, but less safe
   * implementation can be obtained without use of safe_div.
   *
   * @param[in] num2d3 @f$ \nu^{-2/3} @f$ - output from hankel_params
   * @param[in] zeta   zeta in the uniform asymptotic expansions - output
   *    		 from hankel_params
   * @param[out] argp @f$ e^{+i2\pi/3} \nu^{2/3} \zeta @f$
   * @param[out] argm @f$ e^{-i2\pi/3} \nu^{2/3} \zeta @f$
   * @throws std::runtime_error if unable to compute Airy function arguments
   */
  template<typename Tp>
    void
    airy_arg(std::complex<Tp> num2d3, std::complex<Tp> zeta,
	       std::complex<Tp>& argp, std::complex<Tp>& argm)
    {
      using cmplx = std::complex<Tp>;

      // expp and expm are exp(2*pi*i/3) and its reciprocal, respectively.
      static constexpr auto expp = cmplx{-0.5L,  0.8660254037844386467637231707529361834727L};
      static constexpr auto expm = cmplx{-0.5L, -0.8660254037844386467637231707529361834727L};

      try
	{
	  argm = safe_div(num2d3, zeta);
	  argp = expp * argm;
	  argm = expm * argm;
	}
      catch (...)
	{
	   throw std::runtime_error("airy_arg: unable to compute Airy function arguments");
	}
    }


  /**
   * @brief Compute outer factors and associated functions of @c z and @c nu
   * appearing in Olver's uniform asymptotic expansions of the
   * Hankel functions of the first and second kinds and their derivatives.
   * The various functions of z and nu returned by @c hankel_uniform_outer
   * are available for use in computing further terms in the expansions.
   */
  template<typename Tp>
    void
    hankel_uniform_outer(std::complex<Tp> nu, std::complex<Tp> z, Tp eps,
			   std::complex<Tp>& zhat, std::complex<Tp>& __1dnsq,
			   std::complex<Tp>& num1d3, std::complex<Tp>& num2d3,
			   std::complex<Tp>& p, std::complex<Tp>& p2,
			   std::complex<Tp>& etm3h, std::complex<Tp>& etrat,
			   std::complex<Tp>& _Aip, std::complex<Tp>& o4dp,
			   std::complex<Tp>& _Aim, std::complex<Tp>& o4dm,
			   std::complex<Tp>& od2p, std::complex<Tp>& od0dp,
			   std::complex<Tp>& od2m, std::complex<Tp>& od0dm)
    {
      using cmplx = std::complex<Tp>;

      static constexpr cmplx e2pd3{-0.5L,  0.8660254037844386467637231707529361834727L};
      static constexpr cmplx d2pd3{-0.5L, -0.8660254037844386467637231707529361834727L};

      try
	{
	  zhat = safe_div(z, nu);
	  // Try to compute other nu and z dependent parameters except args to Airy functions.
	  cmplx num4d3, nup2, zeta, zetaphf, zetamhf;
	  hankel_params(nu, zhat, p, p2, nup2,
			  __1dnsq, num1d3, num2d3, num4d3,
			  zeta, zetaphf, zetamhf, etm3h, etrat);


	  // Try to compute Airy function arguments.
	  cmplx argp, argm;
	  airy_arg(num2d3, zeta, argp, argm);

	  // Compute Airy functions and derivatives.
	  auto airyp = _Airy<std::complex<Tp>>()(argp);
	  auto airym = _Airy<std::complex<Tp>>()(argm);
	  // Compute partial outer terms in expansions.
	  o4dp = -zetamhf * num4d3 * e2pd3 * airyp.Aip;
	  o4dm = -zetamhf * num4d3 * d2pd3 * airym.Aip;
	  od2p = -zetaphf * num2d3 * airyp.Ai;
	  od0dp = e2pd3 * airyp.Aip;
	  od2m = -zetaphf * num2d3 * airym.Ai;
	  od0dm = d2pd3 * airym.Aip;
	}
      catch (...)
	{
	  throw std::runtime_error("hankel_uniform_outer: unable to compute z/nu");
	}

      return;
    }


  /**
   * @brief Compute the sums in appropriate linear combinations appearing
   * in Olver's uniform asymptotic expansions for the Hankel functions
   * of the first and second kinds and their derivatives, using up to
   * nterms (less than 5) to achieve relative error @c eps.
   *
   * @param[in] p      
   * @param[in] p2     
   * @param[in] num2   
   * @param[in] zetam3hf 
   * @param[in] _Aip     The Airy function value @f$ Ai() @f$.
   * @param[in] o4dp   
   * @param[in] _Aim     The Airy function value @f$ Ai() @f$.
   * @param[in] o4dm   
   * @param[in] od2p   
   * @param[in] od0dp  
   * @param[in] od2m   
   * @param[in] od0dm  
   * @param[in] eps    The error tolerance
   * @param[out] _H1sum  The Hankel function of the first kind.
   * @param[out] _H1psum The derivative of the Hankel function of the first kind.
   * @param[out] _H2sum  The Hankel function of the second kind.
   * @param[out] _H2psum The derivative of the Hankel function of the second kind.
   */
  template<typename Tp>
    void
    hankel_uniform_sum(std::complex<Tp> p, std::complex<Tp> p2,
			 std::complex<Tp> num2, std::complex<Tp> zetam3hf,
			 std::complex<Tp> _Aip, std::complex<Tp> o4dp,
			 std::complex<Tp> _Aim, std::complex<Tp> o4dm,
			 std::complex<Tp> od2p, std::complex<Tp> od0dp,
			 std::complex<Tp> od2m, std::complex<Tp> od0dm,
			 Tp eps,
			 std::complex<Tp>& _H1sum, std::complex<Tp>& _H1psum,
			 std::complex<Tp>& _H2sum, std::complex<Tp>& _H2psum)
    {
      using cmplx = std::complex<Tp>;

      int nterms = 4;

      static constexpr auto zone = cmplx{1, 0};

      // Coefficients for u and v polynomials appearing in Olver's
      // uniform asymptotic expansions for the Hankel functions
      // and their derivatives.

      static constexpr Tp
      s_a[66]
      {
	 0.1000000000000000e+01,
	-0.2083333333333333e+00,
	 0.1250000000000000e+00,
	 0.3342013888888889e+00,
	-0.4010416666666667e+00,
	 0.7031250000000000e-01,
	-0.1025812596450617e+01,
	 0.1846462673611111e+01,
	-0.8912109136581421e+00,
	 0.7324218750000000e-01,
	 0.4669584423426247e+01,
	-0.1120700261622299e+02,
	 0.8789123535156250e+01,
	-0.2364086866378784e+01,
	 0.1121520996093750e+00,
	-0.2821207255820024e+02,
	 0.8463621767460073e+02,
	-0.9181824154324002e+02,
	 0.4253499984741211e+02,
	-0.7368794441223145e+01,
	 0.2271080017089844e+00,
	 0.2125701300392171e+03,
	-0.7652524681411816e+03,
	 0.1059990452528000e+04,
	-0.6995796273761325e+03,
	 0.2181905059814453e+03,
	-0.2649143028259277e+02,
	 0.5725014209747314e+00,
	-0.1919457662318407e+04,
	 0.8061722181737309e+04,
	-0.1358655000643414e+05,
	 0.1165539333686453e+05,
	-0.5305646972656250e+04,
	 0.1200902954101563e+04,
	-0.1080909194946289e+03,
	 0.1727727532386780e+01,
	 0.2020429133096615e+05,
	-0.9698059838863751e+05,
	 0.1925470012325315e+06,
	-0.2034001772804155e+06,
	 0.1222004649830175e+06,
	-0.4119265625000000e+05,
	 0.7109514160156250e+04,
	-0.4939153137207031e+03,
	 0.6074041843414307e+01,
	-0.2429191879005513e+06,
	 0.1311763614662977e+07,
	-0.2998015918538107e+07,
	 0.3763271297656404e+07,
	-0.2813563226586534e+07,
	 0.1268365250000000e+07,
	-0.3316451875000000e+06,
	 0.4521876953125000e+05,
	-0.2499830566406250e+04,
	 0.2438052940368652e+02,
	 0.3284469853072038e+07,
	-0.1970681911843223e+08,
	 0.5095260249266464e+08,
	-0.7410514821153266e+08,
	 0.6634451227472903e+08,
	-0.3756717666076335e+08,
	 0.1328876700000000e+08,
	-0.2785618250000000e+07,
	 0.3081864062500000e+06,
	-0.1388608984375000e+05,
	 0.1100171432495117e+03
      };

      static constexpr Tp
      s_b[66]
      {  0.1000000000000000e+01,
	 0.2916666666666667e+00,
	-0.3750000000000000e+00,
	-0.3949652777777778e+00,
	 0.5156250000000000e+00,
	-0.1171875000000000e+00,
	 0.1146496431327160e+01,
	-0.2130533854166667e+01,
	 0.1089257836341858e+01,
	-0.1025390625000000e+00,
	-0.5075635242854617e+01,
	 0.1238668710214120e+02,
	-0.9961006673177083e+01,
	 0.2793920993804932e+01,
	-0.1441955566406250e+00,
	 0.3015773273462785e+02,
	-0.9140711508856879e+02,
	 0.1005628359759295e+03,
	-0.4753911590576172e+02,
	 0.8502454757690430e+01,
	-0.2775764465332031e+00,
	-0.2247169946128867e+03,
	 0.8146235951180321e+03,
	-0.1138508263826370e+04,
	 0.7604126384523180e+03,
	-0.2411579284667969e+03,
	 0.3002362060546875e+02,
	-0.6765925884246826e+00,
	 0.2013089743407110e+04,
	-0.8497490948317704e+04,
	 0.1440997727955136e+05,
	-0.1245921356699312e+05,
	 0.5730098632812500e+04,
	-0.1315274658203125e+04,
	 0.1208074951171875e+03,
	-0.1993531703948975e+01,
	-0.2106404840887960e+05,
	 0.1014913238950858e+06,
	-0.2024212064239434e+06,
	 0.2150230445535821e+06,
	-0.1300843659496637e+06,
	 0.4424396093750000e+05,
	-0.7727732910156250e+04,
	 0.5459063720703125e+03,
	-0.6883914470672607e+01,
	 0.2520859497081193e+06,
	-0.1365304986690037e+07,
	 0.3131261070473134e+07,
	-0.3946845507298180e+07,
	 0.2965647725320941e+07,
	-0.1345235875000000e+07,
	 0.3545172500000000e+06,
	-0.4883626953125000e+05,
	 0.2737909667968750e+04,
	-0.2724882698059082e+02,
	-0.3395807814193124e+07,
	 0.2042343072273885e+08,
	-0.5295074376688679e+08,
	 0.7725855877372554e+08,
	-0.6943030354332107e+08,
	 0.3949369854080250e+08,
	-0.1404812500000000e+08,
	 0.2965335500000000e+07,
	-0.3310150312500000e+06,
	 0.1509357617187500e+05,
	-0.1215978927612305e+03
      };

      // lambda and mu coefficients appearing in the expansions.
      static constexpr Tp
      s_lambda[21]
      {
	 0.1041666666666667e+00,
	 0.8355034722222222e-01,
	 0.1282265745563272e+00,
	 0.2918490264641405e+00,
	 0.8816272674437577e+00,
	 0.3321408281862768e+01,
	 0.1499576298686255e+02,
	 0.7892301301158652e+02,
	 0.4744515388682643e+03,
	 0.3207490090890662e+04,
	 0.2408654964087401e+05,
	 0.1989231191695098e+06,
	 0.1791902007775344e+07,
	 0.1748437718003412e+08,
	 0.1837073796763307e+09,
	 0.2067904032945155e+10,
	 0.2482751937593589e+11,
	 0.3166945498173489e+12,
	 0.4277112686513472e+13,
	 0.6097113241139256e+14,
	 0.9148694223435640e+15
      };

      static constexpr Tp
      s_mu[21]
      {
	-0.1458333333333333e+00,
	-0.9874131944444445e-01,
	-0.1433120539158951e+00,
	-0.3172272026784136e+00,
	-0.9424291479571203e+00,
	-0.3511203040826354e+01,
	-0.1572726362036805e+02,
	-0.8228143909718595e+02,
	-0.4923553705236705e+03,
	-0.3316218568547973e+04,
	-0.2482767424520859e+05,
	-0.2045265873151298e+06,
	-0.1838444917068210e+07,
	-0.1790568747352892e+08,
	-0.1878356353993943e+09,
	-0.2111438854691369e+10,
	-0.2531915342298413e+11,
	-0.3226140741130003e+12,
	-0.4352813796009286e+13,
	-0.6199585732586975e+14,
	-0.9295073331010611e+15
      };

      std::vector<cmplx> u;
      u.reserve(100);
      std::vector<cmplx> v;
      v.reserve(100);

      auto xtsq = std::real(p2);
      auto ytsq = std::imag(p2);
      auto ytsq2 = ytsq * ytsq;
      auto dr = Tp{2} * xtsq;
      auto ds = std::norm(p2);

      // Compute Debye polynomials u_0,1,2 and v_0,1,2.
      auto pk = p;
      u.push_back(pk * (s_a[1] * p2 + s_a[2]));
      v.push_back(pk * (s_b[1] * p2 + s_b[2]));
      pk *= p;
      u.push_back(pk * cmplx((s_a[3] * xtsq + s_a[4])
				   * xtsq + s_a[5] - s_a[3] * ytsq2,
			(Tp{2} * s_a[3] * xtsq + s_a[4]) * ytsq));
      v.push_back(pk * cmplx((s_b[3] * xtsq + s_b[4])
				   * xtsq + s_b[5] - s_b[3] * ytsq2,
			(Tp{2} * s_b[3] * xtsq + s_b[4]) * ytsq));
      pk *= p;
      u.push_back(pk * cmplx(((s_a[6] * xtsq + s_a[7])
				   * xtsq + s_a[8]) * xtsq
     			+ s_a[9] - (Tp{3} * s_a[6] * xtsq + s_a[7]) * ytsq2,
     			((Tp{3} * s_a[6] * xtsq + Tp{2} * s_a[7]) * xtsq + s_a[8]
     			- s_a[6] * ytsq2) * ytsq));
      v.push_back(pk * cmplx(((s_b[6] * xtsq + s_b[7])
				   * xtsq + s_b[8]) * xtsq
     			+ s_b[9] - (Tp{3} * s_b[6] * xtsq + s_b[7]) * ytsq2,
     			((Tp{3} * s_b[6] * xtsq + Tp{2} * s_b[7]) * xtsq + s_b[8]
     			- s_b[6] * ytsq2) * ytsq));

      // Compute A_0,1, B_0,1, C_0,1, D_0,1 ... note that
      // B_k and C_k are computed up to -zeta^(-1/2) -zeta^(1/2) factors,
      // respectively.  These recurring factors are included as appropriate
      // in the outer factors, thus saving repeated multiplications by them.
      auto _A0 = zone;
      auto _Ak = u[1]
	      + zetam3hf * (s_mu[1] * zetam3hf + s_mu[0] * u[0]);
      auto _B0 = u[0] + s_lambda[0] * zetam3hf;
      auto _Bk = u[2] + zetam3hf * (zetam3hf * (s_lambda[2] * zetam3hf
					 + s_lambda[1] * u[0])
		     + s_lambda[0] * u[1]);
      auto _C0 = v[0] + s_mu[0] * zetam3hf;
      auto _Ck = v[2] + zetam3hf * (zetam3hf * (s_mu[2] * zetam3hf
					 + s_mu[1] * v[0])
		     + s_mu[0] * v[1]);
      auto _D0 = zone;
      auto _Dk = v[1] + zetam3hf * (s_lambda[1] * zetam3hf
			+ s_lambda[0] * v[0]);

      // Compute terms.
      auto _Aterm = _Ak * num2;
      auto _Bterm = _Bk * num2;
      auto _Cterm = _Ck * num2;
      auto _Dterm = _Dk * num2;
      // Compute sum of first two terms to initialize the Kahan summing scheme.
      auto _Asum = _A0 + _Aterm;
      auto _Atemp = _Aterm - (_Asum - _A0);
      auto _Bsum = _B0 + _Bterm;
      auto _Btemp = _Bterm - (_Bsum - _B0);
      auto _Csum = _C0 + _Cterm;
      auto _Ctemp = _Cterm - (_Csum - _C0);
      auto _Dsum = _D0 + _Dterm;
      auto _Dtemp = _Dterm - (_Dsum - _D0);

      // Combine sums in form appearing in expansions.
      _H1sum = _Aip * _Asum + o4dp * _Bsum;
      _H2sum = _Aim * _Asum + o4dm * _Bsum;
      _H1psum = od2p * _Csum + od0dp * _Dsum;
      _H2psum = od2m * _Csum + od0dm * _Dsum;
      auto _H1save = _Aip * _A0 + o4dp * _B0;
      auto _H2save = _Aim * _A0 + o4dm * _B0;
      auto _H1psave = od2p * _C0 + od0dp * _D0;
      auto _H2psave = od2m * _C0 + od0dm * _D0;

      auto converged
	= (l1_norm(_H1sum - _H1save) < eps * l1_norm(_H1sum)
	&& l1_norm(_H2sum - _H2save) < eps * l1_norm(_H2sum)
	&& l1_norm(_H1psum - _H1psave) < eps * l1_norm(_H1psum)
	&& l1_norm(_H2psum - _H2psave) < eps * l1_norm(_H2psum));

      // Save current sums for next convergence test.
      _H1save = _H1sum;
      _H2save = _H2sum;
      _H1psave = _H1psum;
      _H2psave = _H2psum;

      // Maintain index into u_k and v_k coefficients.
      auto index = 10;
      auto indexp = 15;
      // Maintain power of nu^(-2).
      auto num2k = num2;

      for (auto k = 2; k <= nterms; ++k)
	{
	  // Initialize for evaluation of two new u and v polynomials
	  // via Horner's rule modified for complex arguments
          // and real coefficients.
	  auto indexend = indexp;
	  auto ukta = s_a[index];
	  auto vkta = s_b[index];
	  ++index;
	  auto uktb = s_a[index];
	  auto vktb = s_b[index];
	  ++index;
	  auto ukpta = s_a[indexp];
	  auto vkpta = s_b[indexp];
	  ++indexp;
	  auto ukptb = s_a[indexp];
	  auto vkptb = s_b[indexp];
	  ++indexp;

	  // Loop until quantities to evaluate lowest order u and v 
	  // polynomials and partial quantities to evaluate
	  // next highest order polynomials computed.
	  for (; index < indexend; ++index, ++indexp)
	    {
	      auto term = ds * ukta;
	      ukta = uktb + dr * ukta;
	      uktb = s_a[index] - term;
	      term = ds * vkta;
	      vkta = vktb + dr * vkta;
	      vktb = s_b[index] - term;

	      term = ds * ukpta;
	      ukpta = ukptb + dr * ukpta;
	      ukptb = s_a[indexp] - term;
	      term = ds * vkpta;
	      vkpta = vkptb + dr * vkpta;
	      vkptb = s_b[indexp] - term;
	    }

	  // One more iteration for highest order polynomials.
	  auto term = ds * ukpta;
	  ukpta = ukptb + dr * ukpta;
	  ukptb = s_a[indexp] - term;
	  term = ds * vkpta;
	  vkpta = vkptb + dr * vkpta;
	  vkptb = s_b[indexp] - term;
	  ++indexp;

	  // Post multiply and form new polynomials.
	  pk *= p;
	  u.push_back(pk * (ukta * p2 + uktb));
	  v.push_back(pk * (vkta * p2 + vktb));

	  pk *= p;
	  u.push_back(pk * (ukpta * p2 + ukptb));
	  v.push_back(pk * (vkpta * p2 + vkptb));

	  // Update indices in preparation for next iteration.
	  index = indexp;
	  auto i2k = 2 * k - 1;
	  auto i2km1 = i2k - 1;
	  auto i2kp1 = i2k + 1;
	  indexp += i2kp1 + 3;

	  // Start Horner's rule evaluation of A, B, C, and D polynomials.
	  _Ak = s_mu[i2k] * zetam3hf + s_mu[i2km1] * u[0];
	  _Dk = s_lambda[i2k] * zetam3hf + s_lambda[i2km1] * v[0];
	  _Bk = s_lambda[i2kp1] * zetam3hf + s_lambda[i2k] * u[0];
	  _Ck = s_mu[i2kp1] * zetam3hf + s_mu[i2k] * v[0];

	  // Do partial Horner's rule evaluations of A, B, C, and D.
	  for(auto l = 1; l <= i2km1; ++l)
	    {
	      auto i2kl = i2km1 - l;
	      _Ak = _Ak * zetam3hf + s_mu[i2kl] * u[l];
	      _Dk = _Dk * zetam3hf + s_lambda[i2kl] * v[l];
	      i2kl = i2k - l;
	      _Bk = _Bk * zetam3hf + s_lambda[i2kl] * u[l];
	      _Ck = _Ck * zetam3hf + s_mu[i2kl] * v[l];
	    }

	  // Complete the evaluations of A, B, C, and D.
	  _Ak = _Ak * zetam3hf + u[i2k];
	  _Dk = _Dk * zetam3hf + v[i2k];
	  _Bk = zetam3hf
	      * (_Bk * zetam3hf + s_lambda[0] * u[i2k]) + u[i2kp1];
	  _Ck = zetam3hf
	      * (_Ck * zetam3hf + s_mu[0] * v[i2k]) + v[i2kp1];

	  // Evaluate new terms for sums.
	  num2k *= num2;
	  _Aterm = _Ak * num2k + _Atemp;
	  _Bterm = _Bk * num2k + _Btemp;
	  _Cterm = _Ck * num2k + _Ctemp;
	  _Dterm = _Dk * num2k + _Dtemp;

	  // Update sums using Kahan summing scheme.
	  _Atemp = _Asum;
	  _Asum += _Aterm;
	  _Atemp = _Aterm - (_Asum - _Atemp);
	  _Btemp = _Bsum;
	  _Bsum += _Bterm;
	  _Btemp = _Bterm - (_Bsum - _Btemp);
	  _Ctemp = _Csum;
	  _Csum += _Cterm;
	  _Ctemp = _Cterm - (_Csum - _Ctemp);
	  _Dtemp = _Dsum;
	  _Dsum += _Dterm;
	  _Dtemp = _Dterm - (_Dsum - _Dtemp);

	  // Combine sums in form appearing in expansions.
	  _H1sum  = _Aip  * _Asum  + o4dp * _Bsum;
	  _H2sum  = _Aim  * _Asum  + o4dm * _Bsum;
	  _H1psum = od2p * _Csum + od0dp * _Dsum;
	  _H2psum = od2m * _Csum + od0dm * _Dsum;

	  // If convergence criteria met this term, see if it was before.
	  if (l1_norm(_H1sum - _H1save) < eps * l1_norm(_H1sum)
	   && l1_norm(_H2sum - _H2save) < eps * l1_norm(_H2sum)
	   && l1_norm(_H1psum - _H1psave) < eps * l1_norm(_H1psum)
	   && l1_norm(_H2psum - _H2psave) < eps * l1_norm(_H2psum))
	    {
	      if (converged) // Converged twice in a row - done!
		return;
	      else // Converged once...
		converged = true;
	    }
	  else
	    converged = false;
	  // Save combined sums for comparison next iteration.
	  _H1save = _H1sum;
	  _H2save = _H2sum;
	  _H1psave = _H1psum;
	  _H2psave = _H2psum;
	}

      throw std::runtime_error("hankel_uniform_sum: all allowable terms used");

      return;
    }


  /**
   * @brief Compute approximate values for the Hankel functions
   * of the first and second kinds using Olver's uniform asymptotic
   * expansion to of order @c nu along with their derivatives.
   *
   * @param[in] nu The order for which the Hankel functions are evaluated.
   * @param[in] z  The argument at which the Hankel functions are evaluated.
   * @param[out] _H1  The Hankel function of the first kind.
   * @param[out] _H1p The derivative of the Hankel function of the first kind.
   * @param[out] _H2  The Hankel function of the second kind.
   * @param[out] _H2p The derivative of the Hankel function of the second kind.
   */
  template<typename Tp>
    void
    hankel_uniform_olver(std::complex<Tp> nu, std::complex<Tp> z,
			   std::complex<Tp>& _H1, std::complex<Tp>& _H2,
			   std::complex<Tp>& _H1p, std::complex<Tp>& _H2p)
    {
      using namespace std::literals::complex_literals;
      using cmplx = std::complex<Tp>;

      static constexpr Tp
	s_pi(3.141592653589793238462643383279502884195e+0L);
      static constexpr Tp
	s_pi_3(1.047197551196597746154214461093167628063e+0L);
      static constexpr cmplx s_j{1il};
      static constexpr cmplx con1p{ 1.0L, 1.732050807568877293527446341505872366945L}; // 2*exp( pi*j/3) (1,sqrt(3))
      static constexpr cmplx con1m{ 1.0L,-1.732050807568877293527446341505872366945L}; // 2*exp(-pi*j/3)
      static constexpr cmplx con2p{-2.0L, 3.464101615137754587054892683011744733891L}; // 4*exp( 2*pi*j/3) (-2,2sqrt(3))
      static constexpr cmplx con2m{-2.0L,-3.464101615137754587054892683011744733891L}; // 4*exp(-2*pi*j/3)
      static constexpr Tp eps   = 1.0e-06L;
      static constexpr Tp epsai = 1.0e-12L;

      // Extended to accommodate negative real orders.
      bool nuswitch = false;
      if (std::real(nu) < Tp{0})
	{
	  nuswitch = true;
	  nu = -nu;
	}

      // Compute outer factors in the uniform asymptotic expansions
      // for the Hankel functions and their derivatives along with
      // other important functions of nu and z.
      cmplx p, p2,
	    __1dnsq, etm3h, _Aip, o4dp, _Aim, o4dm,
	    od2p, od0dp, od0dm, tmp, zhat, num1d3,
	    num2d3, etrat, od2m, r_factor;
      hankel_uniform_outer(nu, z, epsai, zhat, __1dnsq, num1d3,
			     num2d3, p, p2, etm3h, etrat,
			     _Aip, o4dp, _Aim, o4dm, od2p,
			     od0dp, od2m, od0dm);

      // Compute further terms in the expansions in their appropriate linear combinations.

      hankel_uniform_sum(p, p2, __1dnsq, etm3h,
			   _Aip, o4dp, _Aim, o4dm,
			   od2p, od0dp, od2m, od0dm, eps,
			   _H1, _H1p, _H2, _H2p);

      // Assemble approximations.
      tmp = etrat * num1d3;
      _H1 = con1m * tmp * _H1;
      _H2 = con1p * tmp * _H2;
      tmp = num2d3 / (zhat * etrat);
      _H1p = con2p * tmp * _H1p;
      _H2p = con2m * tmp * _H2p;

      if (nuswitch)
	{
	  r_factor = std::exp(s_j * nu * s_pi);
	  _H1  *= r_factor;
	  _H1p *= r_factor;
	  _H2  /= r_factor;
	  _H2p /= r_factor;
	  nu  = -nu;
	}

      return;
    }


  /**
   * @brief This routine computes the uniform asymptotic approximations
   * of the Hankel functions and their derivatives including a patch
   * for the case when the order equals or nearly equals the argument.
   * At such points, Olver's expressions have zero denominators (and
   * numerators) resulting in numerical problems.  This routine
   * averages results from four surrounding points in the complex plane
   * to obtain the result in such cases.
   *
   * @param[in] nu The order for which the Hankel functions are evaluated.
   * @param[in] z  The argument at which the Hankel functions are evaluated.
   * @param[out] _H1  The Hankel function of the first kind.
   * @param[out] _H1p The derivative of the Hankel function of the first kind.
   * @param[out] _H2  The Hankel function of the second kind.
   * @param[out] _H2p The derivative of the Hankel function of the second kind.
   */
  template<typename Tp>
    void
    hankel_uniform(std::complex<Tp> nu, std::complex<Tp> z,
		     std::complex<Tp>& _H1, std::complex<Tp>& _H2,
		     std::complex<Tp>& _H1p, std::complex<Tp>& _H2p)
    {
      using cmplx = std::complex<Tp>;
      Tp test = std::pow(std::abs(nu), Tp{1} / Tp{3}) / Tp{5};

      if (std::abs(z - nu) > test)
	hankel_uniform_olver(nu, z, _H1, _H2, _H1p, _H2p);
      else
	{
	  Tp r = Tp{2} * test;
	  std::complex<Tp> s_anu[4]{nu + cmplx{r, Tp()},
				      nu + cmplx{Tp(), r},
				      nu - cmplx{r, Tp()},
				      nu - cmplx{Tp(), r}};

	  _H1  = cmplx{};
	  _H2  = cmplx{};
	  _H1p = cmplx{};
	  _H2p = cmplx{};
	  for (auto tnu : s_anu)
	    {
	      std::complex<Tp> th1, th2, th1p, th2p;
	      hankel_uniform_olver(tnu, z, th1, th2, th1p, th2p);
	      _H1  += th1;
	      _H2  += th2;
	      _H1p += th1p;
	      _H2p += th2p;
	    }
	  _H1  /= Tp{4};
	  _H2  /= Tp{4};
	  _H1p /= Tp{4};
	  _H2p /= Tp{4};
	}

      return;
    }


  /**
   *
   * @param[in] nu The order for which the Hankel functions are evaluated.
   * @param[in] z  The argument at which the Hankel functions are evaluated.
   * @param[in] alpha
   * @param[in] indexr
   * @param[out] aorb
   * @param[out] morn
   * @param[out] _H1  The Hankel function of the first kind.
   * @param[out] _H1p The derivative of the Hankel function of the first kind.
   * @param[out] _H2  The Hankel function of the second kind.
   * @param[out] _H2p The derivative of the Hankel function of the second kind.
   */
  template<typename Tp>
    void
    hankel_debye(std::complex<Tp> nu, std::complex<Tp> z,
		   std::complex<Tp> alpha,
		   int indexr, char& aorb, int& morn,
		   std::complex<Tp>& _H1, std::complex<Tp>& _H2,
		   std::complex<Tp>& _H1p, std::complex<Tp>& _H2p)
    {
      using namespace std::literals::complex_literals;
      using cmplx = std::complex<Tp>;

      static constexpr Tp
	s_pi(3.141592653589793238462643383279502884195e+0L);
      static constexpr cmplx s_j{1.0il};
      static constexpr Tp s_toler = Tp{1.0e-8L};
      const auto maxexp
	= std::floor(std::numeric_limits<Tp>::max_exponent
		   * std::log(std::numeric_limits<Tp>::radix));

      auto alphar = std::real(alpha);
      auto alphai = std::imag(alpha);
      auto thalpa = std::sinh(alpha) / std::cosh(alpha);
      auto snhalp = std::sinh(alpha);
      auto denom = std::sqrt(s_pi * z / Tp{2})
		   * std::sqrt(-s_j * std::sinh(alpha));
      if (std::abs(std::real(nu * (thalpa - alpha))) > maxexp)
	throw std::runtime_error("hankel_debye: argument would overflow Hankel function evaluation");
      auto s1 = std::exp(+nu * (thalpa - alpha) - s_j * s_pi / Tp{4})
		/ denom;
      auto s2 = std::exp(-nu * (thalpa - alpha) + s_j * s_pi / Tp{4})
		/ denom;
      auto exparg = nu * (thalpa - alpha) - s_j * s_pi / Tp{4};
      if (indexr == 0)
	{
	  _H1 = Tp{0.5L} * s1 - s2;
	  _H2 = Tp{0.5L} * s1 + s2;
	  _H1p = snhalp * (Tp{0.5L} * s1 + s2);
	  _H2p = snhalp * (Tp{0.5L} * s1 - s2);
	}
      else if (indexr == 1)
	{
	  _H1 = s1;
	  _H2 = s2;
	  _H1p = +snhalp * s1;
	  _H2p = -snhalp * s2;
	}
      else if (indexr == 2)
	{
	  auto jdbye = s1 / Tp{2};
	  _H2 = s2;
	  _H1 = Tp{2} * jdbye - _H2;
	  _H1p = +snhalp * (s1 + s2);
	  _H2p = -snhalp * s2;
	}
      else if (indexr == 3)
	{
	  _H1 = s1;
	  _H2 = s2 - s1;
	  _H1p = +snhalp * s1;
	  _H2p = -snhalp * (s1 + s2);
	}
      else if (indexr == 4)
	{
	  _H1 = s1;
	  _H2 = s2 - std::exp(+Tp{2} * s_j * nu * s_pi) * s1;
	  _H1p = +snhalp * s1;
	  _H2p = -snhalp
		* (s2 + std::exp(+Tp{2} * s_j * nu * s_pi) * s1);
	}
      else if (indexr == 5)
	{
	  _H1 = s1 - std::exp(-Tp{2} * s_j * nu * s_pi) * s2;
	  _H2 = s2;
	  _H1p = +snhalp
		* (s1 + std::exp(-Tp{2} * s_j * nu * s_pi) * s2);
	  _H2p = -snhalp * s2;
	}
      else if (aorb == 'A')
	{
	  cmplx sinrat;
	  if ((std::abs(std::imag(nu)) < s_toler)
	   && (std::abs(std::fmod(std::real(nu), 1)) < s_toler))
	    sinrat = morn;
	  else
	    sinrat = std::sin(Tp(morn) * nu * s_pi)
		    / std::sin(nu * s_pi);
	  if (indexr == 6)
	    {
	      _H2 = s2
		   - std::exp(s_j * Tp(morn + 1) * nu * s_pi)
		   * sinrat * s1;
	      _H1 = s1 - _H2;
	      _H2p = -snhalp
		    * (s2 + std::exp(s_j * Tp(morn + 1) * nu * s_pi)
			     * sinrat * s1);
	      _H1p = +snhalp
		    * ((Tp{1} + std::exp(s_j * Tp(morn + 1) * nu * s_pi)
			  * sinrat) * s1 + s2);
	    }
	  else if (indexr == 7)
	    {
	      _H1 = s1
		   - std::exp(-s_j * Tp(morn + 1) * nu * s_pi)
		    * sinrat * s2;
	      _H2 = s2 - _H1;
	      _H1p = +snhalp
		    * (s1 + std::exp(-s_j * Tp(morn + 1) * nu * s_pi)
			     * sinrat * s2);
	      _H2p = -snhalp
		     * ((Tp{1} + std::exp(-s_j * Tp(morn + 1) * nu * s_pi)
			   * sinrat) * s2 + s1);
	    }
	  else
	    throw std::runtime_error("hankel_debye: unexpected region");
	}
      else
	{
	  cmplx sinrat;
	  if ((std::abs(std::imag(nu)) < s_toler)
	   && (std::abs(std::fmod(std::real(nu), 1)) < s_toler))
	    sinrat = -morn;
	  else
	    sinrat = std::sin(Tp(morn) * nu * s_pi) / std::sin(nu * s_pi);
	  if (indexr == 6)
	    {
	      _H1 = s1 - std::exp(s_j * Tp(morn - 1) * nu * s_pi)
		   * sinrat * s2;
	      _H2 = s2 - std::exp(Tp{2} * s_j * nu * s_pi) * _H2;
	      _H1p = +snhalp
		    * (s1 + std::exp(s_j * Tp(morn - 1) * nu * s_pi)
			    * sinrat * s2);
	      _H2p = -snhalp
		    * ((Tp{1} + std::exp(s_j * Tp(morn + 1) * nu * s_pi)
			  * sinrat) * s2
		      + std::exp(Tp{2} * s_j * nu * s_pi) * s1);
	    }
	  else if (indexr == 7)
	    {
	      _H2 = s2
		   - std::exp(-s_j * Tp(morn - 1) * nu * s_pi)
		   * sinrat * s1;
	      _H1 = s1 - std::exp(-Tp{2} * s_j * nu * s_pi) * _H2;
	      _H2p = -snhalp
		    * (s2 + std::exp(-s_j * Tp(morn - 1) * nu * s_pi)
			    * sinrat * s1);
	      _H1p = +snhalp
		    * ((Tp{1} + std::exp(-s_j * Tp(morn + 1) * nu * s_pi)
				    * sinrat) * s1
				+ std::exp(-Tp{2} * s_j * nu * s_pi) * s2);
	    }
	  else
	    throw std::runtime_error("hankel_debye: unexpected region");
	}

      return;
    }


  /**
   *
   * @param[in] nu The order for which the Hankel functions are evaluated.
   * @param[in] z  The argument at which the Hankel functions are evaluated.
   * @param[out] _H1  The Hankel function of the first kind.
   * @param[out] _H1p The derivative of the Hankel function of the first kind.
   * @param[out] _H2  The Hankel function of the second kind.
   * @param[out] _H2p The derivative of the Hankel function of the second kind.
   */
  template<typename Tp>
    void
    hankel(std::complex<Tp> nu, std::complex<Tp> z,
	     std::complex<Tp>& _H1, std::complex<Tp>& _H2,
	     std::complex<Tp>& _H1p, std::complex<Tp>& _H2p)
    {
      static constexpr Tp
	s_pi(3.141592653589793238462643383279502884195e+0L);

      int indexr;

      auto test = std::abs((nu - z) / std::pow(nu, Tp{1}/Tp{3}));
      if (test < Tp{4})
	hankel_uniform(z, nu, _H1, _H2, _H1p, _H2p);
      else
	{
	  auto sqtrm = std::sqrt((nu / z) * (nu / z) - Tp{1});
	  auto alpha = std::log((nu / z) + sqtrm);
	  if (std::imag(alpha) < Tp{0})
	    alpha = -alpha;
	  auto alphar = std::real(alpha);
	  auto alphai = std::imag(alpha);
	  char aorb;
	  if (std::real(nu) > std::real(z)
	   && std::abs(std::imag(nu / z)) <= Tp{0})
	    {
	      indexr = 0;
	      aorb = ' ';
	    }
	  else
	    debye_region(alpha, indexr, aorb);
	  auto morn = 0;
	  if (aorb == 'A')
	    {
	      auto mfun = ((alphar * std::tanh(alphar) - Tp{1})
			  * std::tan(alphai) + alphai) / s_pi;
	      morn = int(mfun);
	      if (mfun < 0 && std::fmod(mfun, 1) != Tp{0})
		--morn;
	    }
	  else if (aorb == 'B')
	    {
	      auto nfun = ((Tp{1} - alphar * std::tanh(alphar))
			  * std::tan(alphai) - alphai) / s_pi;
	      morn = int(nfun) + 1;
	      if (nfun < Tp{0} && std::fmod(nfun, Tp{1}) != Tp{0})
		--morn;
	    }
	  hankel_debye(nu, z, alpha, indexr, aorb, morn,
			 _H1, _H2, _H1p, _H2p);
	}

      return;
    }


  /**
   * @brief Return the complex cylindrical Hankel function of the first kind.
   *
   * @param[in] nu The order for which the cylindrical Hankel function of the first kind is evaluated.
   * @param[in] z  The argument at which the cylindrical Hankel function of the first kind is evaluated.
   * @return The complex cylindrical Hankel function of the first kind.
   */
  template<typename Tp>
    std::complex<Tp>
    cyl_hankel_1(std::complex<Tp> nu, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      hankel(nu, z, _H1, _H1p, _H2, _H2p);
      return _H1;
    }

  /**
   * @brief Return the complex cylindrical Hankel function of the second kind.
   *
   * @param[in] nu The order for which the cylindrical Hankel function of the second kind is evaluated.
   * @param[in] z  The argument at which the cylindrical Hankel function of the second kind is evaluated.
   * @return The complex cylindrical Hankel function of the second kind.
   */
  template<typename Tp>
    std::complex<Tp>
    cyl_hankel_2(std::complex<Tp> nu, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      hankel(nu, z, _H1, _H1p, _H2, _H2p);
      return _H2;
    }

  /**
   * @brief Return the complex cylindrical Bessel function.
   *
   * @param[in] nu The order for which the cylindrical Bessel function is evaluated.
   * @param[in] z  The argument at which the cylindrical Bessel function is evaluated.
   * @return The complex cylindrical Bessel function.
   */
  template<typename Tp>
    std::complex<Tp>
    cyl_bessel(std::complex<Tp> nu, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      hankel(nu, z, _H1, _H1p, _H2, _H2p);
      return (_H1 + _H2) / Tp{2};
    }

  /**
   * @brief Return the complex cylindrical Neumann function.
   *
   * @param[in] nu The order for which the cylindrical Neumann function is evaluated.
   * @param[in] z  The argument at which the cylindrical Neumann function is evaluated.
   * @return The complex cylindrical Neumann function.
   */
  template<typename Tp>
    std::complex<Tp>
    cyl_neumann(std::complex<Tp> nu, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      hankel(nu, z, _H1, _H1p, _H2, _H2p);
      return (_H1 - _H2) / std::complex<Tp>{0, 2};
    }

  /**
   * @brief Helper to compute complex spherical Hankel functions
   *        and their derivatives.
   *
   * @param[in] n The order for which the spherical Hankel functions are evaluated.
   * @param[in] z The argument at which the spherical Hankel functions are evaluated.
   * @param[out] _H1  The spherical Hankel function of the first kind.
   * @param[out] _H1p The derivative of the spherical Hankel function of the first kind.
   * @param[out] _H2  The spherical Hankel function of the second kind.
   * @param[out] _H2p The derivative of the spherical Hankel function of the second kind.
   */
  template<typename Tp>
    void
    sph_hankel(unsigned int n, std::complex<Tp> z,
		 std::complex<Tp>& _H1, std::complex<Tp>& _H1p,
		 std::complex<Tp>& _H2, std::complex<Tp>& _H2p)
    {
      static constexpr Tp
	s_pi(3.141592653589793238462643383279502884195e+0L);
      std::complex<Tp> nu(n + Tp{0.5});
      hankel(nu, z, _H1, _H1p, _H2, _H2p);
      std::complex<Tp> fact = std::sqrt(s_pi / (Tp{2} * z));
      _H1 *= fact;
      _H1p = fact * _H1p - _H1 / (Tp{2} * z);
      _H2 *= fact;
      _H2p = fact * _H2p - _H2 / (Tp{2} * z);
    }

  /**
   * @brief Return the complex spherical Hankel function of the first kind.
   *
   * @param[in] n The order for which the spherical Hankel function of the first kind is evaluated.
   * @param[in] z The argument at which the spherical Hankel function of the first kind is evaluated.
   * @return The complex spherical Hankel function of the first kind.
   */
  template<typename Tp>
    std::complex<Tp>
    sph_hankel_1(unsigned int n, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      sph_hankel(n, z, _H1, _H1p, _H2, _H2p);
      return _H1;
    }

  /**
   * @brief Return the complex spherical Hankel function of the second kind.
   *
   * @param[in] n The order for which the spherical Hankel function of the second kind is evaluated.
   * @param[in] z The argument at which the spherical Hankel function of the second kind is evaluated.
   * @return The complex spherical Hankel function of the second kind.
   */
  template<typename Tp>
    std::complex<Tp>
    sph_hankel_2(unsigned int n, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      sph_hankel(n, z, _H1, _H1p, _H2, _H2p);
      return _H2;
    }

  /**
   * @brief Return the complex spherical Bessel function.
   *
   * @param[in] n The order for which the spherical Bessel function is evaluated.
   * @param[in] z The argument at which the spherical Bessel function is evaluated.
   * @return The complex spherical Bessel function.
   */
  template<typename Tp>
    std::complex<Tp>
    sph_bessel(unsigned int n, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      sph_hankel(n, z, _H1, _H1p, _H2, _H2p);
      return (_H1 + _H2) / Tp{2};
    }

  /**
   * @brief Return the complex spherical Neumann function.
   *
   * @param[in] n The order for which the spherical Neumann function is evaluated.
   * @param[in] z The argument at which the spherical Neumann function is evaluated.
   * @return The complex spherical Neumann function.
   */
  template<typename Tp>
    std::complex<Tp>
    sph_neumann(unsigned int n, std::complex<Tp> z)
    {
      std::complex<Tp> _H1, _H1p, _H2, _H2p;
      sph_hankel(n, z, _H1, _H1p, _H2, _H2p);
      return (_H1 - _H2) / std::complex<Tp>{0, 2};
    }

} // namespace detail
} // namespace emsr

#endif // SF_HANKEL_OLD_TCC
