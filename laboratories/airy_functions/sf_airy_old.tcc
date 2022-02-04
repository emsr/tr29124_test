// TR29124 math special functions -*- C++ -*-

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

/** @file emsr/sf_airy.tcc
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_SF_AIRY_TCC
#define _GLIBCXX_BITS_SF_AIRY_TCC 1

#pragma GCC system_header

#include <bits/complex_util.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

// Implementation-space details.
namespace __detail
{
  /**
   * @brief This function evaluates @f$ Ai(z) @f$ and @f$ Ai'(z) @f$ from their asymptotic
   * expansions for @f$ |arg(z)| < 2*\pi/3 @f$.  For speed, the number
   * of terms needed to achieve about 16 decimals accuracy is tabled
   * and determined from abs(z).
   *
   * Note that for speed and since this function
   * is called by another, checks for valid arguments are not
   * made.
   *
   * @see Digital Library of Mathematical Functions Sec. 9.7 Asymptotic Expansions
   *     http://dlmf.nist.gov/9.7
   *
   * @param[in]  z Complex input variable set equal to the value at which
   *  	    @f$ Ai(z) @f$ and @f$ Bi(z) @f$ and their derivative are evaluated.
   *  	    This function assumes @f$ |z| > 15 @f$ and @f$ |arg(z)| < 2\pi/3 @f$.
   * @param[inout] _Ai  The value computed for @f$ Ai(z) @f$.
   * @param[inout] _Aip The value computed for @f$ Ai'(z) @f$.
   * @param[in]    sign  The sign of the series terms amd exponent.
   *                    The default (-1) gives the Airy Ai functions
   *                    for @f$ |arg(z)| < \pi @f$.
   *                    The value +1 gives the Airy Bi functions
   *                    for @f$ |arg(z)| < \pi/3 @f$.
   */
  template<typename Tp>
    void
    airy_asymp_absarg_ge_pio3(std::complex<Tp> z,
				std::complex<Tp>& _Ai,
				std::complex<Tp>& _Aip,
				int sign = -1)
    {
      constexpr Tp s_2d3   = Tp{2} / Tp{3};
      // 1/(2 sqrt(pi)))
      constexpr Tp s_pmhd2 = Tp{2.820947917738781434740397257803862929219e-01L};
      constexpr int s_ncoeffs = 15;
      constexpr int s_numnterms = 5;
      constexpr int s_nterms[5]{ s_ncoeffs, 12, 11, 11, 9 };

      // Coefficients for the expansion.
      constexpr Tp
      s_u[s_ncoeffs]
      {
	0.5989251356587907e+05,
	0.9207206599726415e+04,
	0.1533169432012796e+04,
	0.2784650807776026e+03,
	0.5562278536591708e+02,
	0.1234157333234524e+02,
	0.3079453030173167e+01,
	0.8776669695100169e+00,
	0.2915913992307505e+00,
	0.1160990640255154e+00,
	0.5764919041266972e-01,
	0.3799305912780064e-01,
	0.3713348765432099e-01,
	0.6944444444444444e-01,
	0.1000000000000000e+01
      };

      constexpr Tp
      s_v[s_ncoeffs]
      {
	-0.6133570666385206e+05,
	-0.9446354823095932e+04,
	-0.1576357303337100e+04,
	-0.2870332371092211e+03,
	-0.5750830351391427e+02,
	-0.1280729308073563e+02,
	-0.3210493584648621e+01,
	-0.9204799924129446e+00,
	-0.3082537649010791e+00,
	-0.1241058960272751e+00,
	-0.6266216349203231e-01,
	-0.4246283078989483e-01,
	-0.4388503086419753e-01,
	-0.9722222222222222e-01,
	 0.1000000000000000e+01
      };

      // Compute -zeta and z^(1/4).
      auto z1d4 = std::sqrt(z);
      auto zetam = z * z1d4;
      zetam *= s_2d3;
      z1d4 = std::sqrt(z1d4);

      // Compute outer factors in the expansions.
      auto zoutpr = std::exp(Tp(sign) * zetam);
      zoutpr *= s_pmhd2;
      auto zout = zoutpr / z1d4;
      zoutpr *= -z1d4;

      // Determine number of terms to use.
      auto nterm = s_nterms[std::min(s_numnterms - 1,
					(int(std::abs(z)) - 10) / 5)];
      // Initialize for modified Horner's rule evaluation of sums.
      // It is assumed that at least three terms are used.
      zetam = Tp(sign) / zetam;
      auto r = Tp{2} * std::real(zetam);
      auto s = std::norm(zetam);
      auto index = s_ncoeffs - nterm;// + 1;
      auto al = s_u[index];
      auto alpr = s_v[index];
      ++index;
      auto be = s_u[index];
      auto bepr = s_v[index];
      ++index;

      for (int k = index; k < s_ncoeffs; ++k)
	{
	  auto term = s * al;
	  al = be + r * al;
	  be = s_u[k] - term;
	  term = s * alpr;
	  alpr = bepr + r * alpr;
	  bepr = s_v[k] - term;
	}

      _Ai = zout * al * zetam + be;
      _Aip = zoutpr * alpr * zetam + bepr;

      return;
    }


  /**
   * @brief This function evaluates @f$ Ai(z) @f$ and @f$ Ai'(z) @f$ from
   * their asymptotic expansions for @f$ |arg(-z)| < pi/3 @f$.  For speed,
   * the number of terms needed to achieve about 16 decimals accuracy
   * is tabled and determined from @f$ |z| @f$.
   *
   * Note that for speed and since this function is called by another,
   * checks for valid arguments are not made.
   * This function assumes @f$ |z| > 15 @f$ and @f$ |arg(-z)| < pi/3 @f$.
   *
   * @param[in] z  The value at which the Airy function and its derivative
   *              are evaluated.
   * @param[out] _Ai  The computed value of the Airy function @f$ Ai(z) @f$.
   * @param[out] _Aip The computed value of the Airy function derivative @f$ Ai'(z) @f$.
   */
  template<typename Tp>
    void
    airy_asymp_absarg_lt_pio3(std::complex<Tp> z,
				std::complex<Tp>& _Ai,
				std::complex<Tp>& _Aip)
    {
      constexpr Tp s_2d3 {Tp{2} / Tp{3}};
      constexpr Tp s_9d4 {Tp{9} / Tp{4}};
      constexpr Tp s_pimh{5.641895835477562869480794515607725858438e-01};
      constexpr Tp s_pid4 = emsr::math_constants<Tp>::pi_quarter;

      constexpr std::complex<Tp> s_zone{1};
      constexpr int s_ncoeffs = 9;
      constexpr int s_numnterms = 5;
      constexpr int s_nterms[s_numnterms]{ s_ncoeffs, 7, 6, 6, 5 };

      // coefficients for the expansion.
      constexpr Tp
      s_u_cos[s_ncoeffs]
      {
	0.2519891987160237e+08,
	0.4195248751165511e+06,
	0.9207206599726415e+04,
	0.2784650807776026e+03,
	0.1234157333234524e+02,
	0.8776669695100169e+00,
	0.1160990640255154e+00,
	0.3799305912780064e-01,
	0.6944444444444444e-01
      };
      constexpr Tp
      s_u_sin[s_ncoeffs]
      {
	0.3148257417866826e+07,
	0.5989251356587907e+05,
	0.1533169432012796e+04,
	0.5562278536591708e+02,
	0.3079453030173167e+01,
	0.2915913992307505e+00,
	0.5764919041266972e-01,
	0.3713348765432099e-01,
	0.1000000000000000e+01
      };

      constexpr Tp
      s_v_sin[s_ncoeffs]
      {
	-0.2569790838391133e+08,
	-0.4289524004000691e+06,
	-0.9446354823095932e+04,
	-0.2870332371092211e+03,
	-0.1280729308073563e+02,
	-0.9204799924129446e+00,
	-0.1241058960272751e+00,
	-0.4246283078989483e-01,
	-0.9722222222222222e-01
      };
      constexpr Tp
      s_v_cos[s_ncoeffs]
      {
	-0.3214536521400865e+07,
	-0.6133570666385206e+05,
	-0.1576357303337100e+04,
	-0.5750830351391427e+02,
	-0.3210493584648621e+01,
	-0.3082537649010791e+00,
	-0.6266216349203231e-01,
	-0.4388503086419753e-01,
	 0.1000000000000000e+01
      };

      // Set up working value of z.
      z = -z;
      // Compute zeta and z^(1/4).
      auto z1d4 = std::sqrt(z);
      auto zeta = z * z1d4;
      zeta *= s_2d3;
      z1d4 = std::sqrt(z1d4);

      // Compute sine and cosine factors in the expansions.
      auto zetaarg = zeta + s_pid4;
      auto sinzeta = std::sin(zetaarg);
      auto coszeta = std::cos(zetaarg);

      // Determine number of terms to use.
      auto nterm = s_nterms[std::min(s_numnterms - 1,
					(int(std::abs(z)) - 10) / 5)];
      // Initialize for modified Horner's rule evaluation of sums
      // it is assumed that at least three terms are used.
      z = std::pow(s_zone / z, Tp(3));
      z *= s_9d4;
      auto r = -2 * std::real(z);
      auto s = std::norm(z);
      auto index = s_ncoeffs - nterm;

      auto als = s_u_sin[index];
      auto alc = s_u_cos[index];
      auto alprs = s_v_sin[index];
      auto alprc = s_v_cos[index];
      ++index;

      auto bes = s_u_sin[index];
      auto bec = s_u_cos[index];
      auto beprs = s_v_sin[index];
      auto beprc = s_v_cos[index];
      ++index;

      // Loop until components contributing to sums are computed.
      for (int k = index; k < s_ncoeffs; ++k)
	{
	  auto term = s * als;
	  als = bes + r * als;
	  bes = s_u_sin[k] - term;
	  term = s * alc;
	  alc = bec + r * alc;
	  bec = s_u_cos[k] - term;
	  term = s * alprs;
	  alprs = beprs + r * alprs;
	  beprs = s_v_sin[k] - term;
	  term = s * alprc;
	  alprc = beprc + r * alprc;
	  beprc = s_v_cos[k] - term;
	}

      // Complete evaluation of the Airy functions.
      zeta = s_zone / zeta;
      _Ai = sinzeta * als * z + bes
           - zeta * coszeta * alc * z + bec;
      _Ai *= s_pimh / z1d4;
      _Aip = coszeta * alprc * z + beprc
	   + zeta * sinzeta * alprs * z + beprs;
      _Aip *= -s_pimh * z1d4;

      return;
    }


  /**
   * Compute the modified Bessel functions of the first kind of
   * orders +-1/3 and +-2/3 needed to compute the Airy functions
   * and their derivatives from their representation in terms of the
   * modified Bessel functions.  This function is only used for z
   * less than two in modulus and in the closed right half plane.
   * This stems from the fact that the values of the modified 
   * Bessel functions occuring in the representations of the Airy
   * functions and their derivatives are almost equal for z large
   * in the right half plane.  This means that loss of significance
   * occurs if these representations are used for z to large in
   * magnitude.  This algorithm is also not used for z too small,
   * since a low order rational approximation can be used instead.
   *
   * This routine is an implementation of a modified version of
   * Miller's backward recurrence algorithm for computation by
   * from the recurrence relation
   * @f[
   *  I_{\nu-1} = (2\nu/z)I_\nu + I_{\nu+1}
   * @f]
   * satisfied by the modified Bessel functions of the first kind.
   * the normalization relationship used is
   * @f[
   *  \frac{z/2)^\nu e^z}{\Gamma(\nu+1)} = I_\nu(z)
   *    + 2\sum_{k=1}^{\infty}\frac{(k+\nu)\Gamma(2\nu+k)}{k!\Gamma(1+2\nu)}
   *                  I_{\nu+k}(z).
   * @f]
   *
   * This modification of the algorithm is given in part in
   *
   * Olver, F. W. J. and Sookne, D. J., Note on Backward Recurrence
   * Algorithms, Math. of Comp., Vol. 26, no. 120, Oct. 1972.
   *
   * And further elaborated for the Bessel functions in
   *
   * Sookne, D. J., Bessel Functions I and J of Complex Argument 
   * and Integer Order, J. Res. NBS - Series B, Vol 77B, Nos.
   * 3 & 4, July-December, 1973.
   *
   * Insight was also gained from
   *
   * Cody, W. J., Preliminary Report on Software for the Modified
   *  Bessel Functions of the First Kind, Argonne National
   *  Laboratory, Applied Mathematics Division, Tech. Memo. 
   *  No. 357, August, 1980.
   *
   * Cody implements the algorithm of Sookne for fractional order
   * and nonnegative real argument.  Like Cody, we do not change
   * the convergence testing mechanism in any substantial way.
   * However, we do trim the overhead by making the additional
   * assumption that performing the convergence test for the
   * functions of order 2/3 will suffice for order 1/3 as well.
   * This assumption has not been established by rigourous
   * analysis at this time.  For speed the convergence
   * tests are performed in the 1-norm instead of the usual Euclidean
   * norm used in the complex plane using the inequality
   * @f[
   *  |x| + |y| <= \sqrt(2) \sqrt(x^2 + y^2) 
   * @f]
   *
   * @param[in]  z     The argument of the modified Bessel functions.
   * @param[in]  eps   The maximum relative error required in the results.
   * @param[out] _Ip1d3 The value of @f$ I_(+1/3)(z) @f$.
   * @param[out] _Im1d3 The value of @f$ I_(-1/3)(z) @f$.
   * @param[out] _Ip2d3 The value of @f$ I_(+2/3)(z) @f$.
   * @param[out] _Im2d3 The value of @f$ I_(-2/3)(z) @f$.
   */
  template<typename Tp>
    void
    airy_bessel_i(const std::complex<Tp>& z, Tp eps,
		    std::complex<Tp>& _Ip1d3,
		    std::complex<Tp>& _Im1d3,
		    std::complex<Tp>& _Ip2d3,
		    std::complex<Tp>& _Im2d3)
    {
      using cmplx = std::complex<Tp>;

      constexpr cmplx s_zero{0}, s_zone{1};
      constexpr Tp
    	s_1d3  {Tp{1} / Tp{3}}, s_2d3  {Tp{2} / Tp{3}},
    	s_4d3  {Tp{4} / Tp{3}}, s_5d3  {Tp{5} / Tp{3}},
    	s_8d3  {Tp{8} / Tp{3}}, s_10d3 {Tp{10} / Tp{3}},
    	s_14d3 {Tp{14} / Tp{3}}, s_16d3 {Tp{16} / Tp{3}},
    	s_gamma4d3{8.929795115692492112185643136582258813769e-1L},
	s_gamma5d3{9.027452929509336112968586854363425236809e-1L},
    	s_sqrt2 = emsr::math_constants<Tp>::root_2;

      // Compute 1/z for use in recurrence for speed and abs(z).
      auto __1dz = safe_div(Tp{1}, z);

      // Initialize for forward recursion based on order 2/3.
      int n = 0;
      auto d2n = Tp(2 * n) + s_4d3;
      auto plast2 = s_zone;
      auto p2 = d2n * __1dz;

      auto weak_test = 20 * s_sqrt2 / eps;
      bool converged = false;

      while (true)
	{
	  do
	    {
    	      // Update n dependent quantities.
    	      ++n;
    	      d2n += 2;
    	      // Interchange values.
    	      auto pold2 = plast2;
    	      plast2 = p2;
    	      // Recur forward one step.
    	      p2 = __1dz * d2n * plast2 + pold2;
	    }
	  while (l1_norm(p2) < weak_test);

	  // If strong convergence, then weak and strong convergence.
	  if (converged)
    	    break;

	  // Calculate strong convergence test.
	  // See the Olver and Sookne papers cited for details.
	  auto rep2 = std::real(p2);
	  auto imp2 = std::imag(p2);
	  auto replast2 = std::real(plast2);
	  auto implast2 = std::imag(plast2);
	  // Compute scale factor to avoid possible overflow.
	  auto lamn = linf_norm(p2);
	  // Compute the kappa_n of strong convergence lemma.
	  auto kapn = std::sqrt(((rep2 / lamn) * (rep2 / lamn)
    			       + (imp2 / lamn) * (imp2 / lamn))
    			     / ((replast2 / lamn) * (replast2 / lamn)
    			      + (implast2 / lamn) * (implast2 / lamn)));
	  // Compute quantity needed for lambda_n of strong convergence lemma
	  lamn = Tp(n + 1) / std::abs(z);
	  // Determine appropriate value for rho_n of lemma
	  if (kapn + 1 / kapn > 2 * lamn)
    	    kapn = lamn + std::sqrt(lamn * lamn - 1);
	  // Compute test value - sqrt(2) multiple already included.
	  weak_test *= std::sqrt(kapn - 1 / kapn);
	  // Set strong convergence test flag.
	  converged = true;
	}

      // Prepare for backward recurrence for both orders 1/3 and 2/3.
      auto rn = Tp(n);
      ++n;
      d2n = Tp(2 * n);
      auto plast1 = s_zero;
      plast2 = s_zero;
      // Carefully compute 1/p2 to avoid overflow in complex divide.
      auto p1 = safe_div(Tp{1}, p2);
      p2 = p1;
      // Set up n dependent parameters used in normalization sum.
      auto rnpn1 = rn + s_1d3;
      auto rnpn2 = rn + s_2d3;
      auto rnp2n1 = (rn - Tp(1)) + s_2d3;
      auto rnp2n2 = (rn - Tp(1)) + s_4d3;
      // Initialize normalization sum
      auto fact1 = rnpn1 * rnp2n1 / rn;
      auto sum1 = fact1 * p1;
      auto fact2 = rnpn2 * rnp2n2 / rn;
      auto sum2 = fact2 * p2;
      // Set ending loop index to correspond to k=1 term of the
      // normalization relationship.
      auto nend = n - 3;

      // If backward recurrence loop will be nontrivial...
      if (nend > 0)
	{
	  // ...loop until backward recursion to k=1 term of normalization.
	  for (int l = 1; l <= nend; ++l)
	    {
    	      // Update n dependent quantities.
    	      --n;
    	      d2n -= 2;
    	      fact1 = d2n + s_2d3;
    	      fact2 = d2n + s_4d3;
    	      // Interchanges for order 1/3 recurrence.
    	      auto pold1 = plast1;
    	      plast1 = p1;
    	      // Recur back one step for order 1/3.
    	      p1 = __1dz * fact1 * plast1 + pold1;
    	      // Interchanges for order 2/3 recurrence
    	      auto pold2 = plast2;
    	      plast2 = p2;
    	      // Recur back one step for order 2/3.
    	      p2 = __1dz * fact2 * plast2 + pold2;
    	      // Update quantities for computing normalization sums.
    	      rn -= Tp(1);
    	      rnpn1 = rn + s_1d3;
    	      rnp2n1 = rn - s_1d3;
    	      rnpn2 = rn + s_2d3;
    	      rnp2n2 = rnpn1;
    	      fact1 = rnp2n1 / rn;
    	      fact2 = rnp2n2 / rn;
    	      // Update normalization sums.
    	      sum1 += rnpn1 * p1;
    	      sum1 *= fact1;
    	      sum2 += rnpn2 * p2;
    	      sum2 *= fact2;
	    }
	}

      // Perform last two recurrence steps for order 1/3
      auto pold1 = plast1;
      plast1 = p1;
      p1 = s_14d3 * plast1 * __1dz + pold1;
      sum1 += s_4d3 * p1;
      pold1 = plast1;
      plast1 = p1;
      p1 = s_8d3 * plast1 * __1dz + pold1;
      sum1 = Tp(2) * sum1 + p1;

      // Compute scale factor and scale results for order 1/3 case
      auto zd2pow = std::pow(Tp(0.5L) * z, -s_1d3);
      pold1 = zd2pow * std::exp(-z);
      sum1 *= s_gamma4d3 * pold1;
      plast1 /= sum1;
      _Ip1d3 = p1 / sum1;

      // Perform last two recurrence steps for order 2/3
      auto pold2 = plast2;
      plast2 = p2;
      p2 = s_16d3 * plast2 * __1dz + pold2;
      sum2 += s_5d3 * p2;
      pold2 = plast2;
      plast2 = p2;
      p2 = s_10d3 * plast2 * __1dz + pold2;
      sum2 = Tp(2) * sum2 + p2;

      // Compute scale factor and scale results for order 2/3 case
      sum2 *= s_gamma5d3 * zd2pow * pold1;
      plast2 /= sum2;
      _Ip2d3 = p2 / sum2;

      // Recur back one step from order +1/3 to get order -2/3
      _Im2d3 = s_2d3 * _Ip1d3 * __1dz + plast1;

      // Recur back one step from order +2/3 to get order -1/3
      _Im1d3 = s_4d3 * _Ip2d3 * __1dz + plast2;

      return;
    }


  /**
   * @brief Compute approximations to the modified Bessel functions
   * of the second kind of orders 1/3 and 2/3 for moderate arguments.
   *
   * This routine computes
   * @f[
   *   E_\nu(z) = \exp{z} \sqrt{2 z/\pi} K_\nu(z), for \nu = 1/3 and \nu = 2/3
   * @f]
   * using a rational approximation given in
   *
   * Luke, Y. L., Mathematical functions and their approximations,
   *  Academic Press, pp 366-367, 1975.
   *
   * Though the approximation converges in @f$ |\arg(z)| <= pi @f$,
   * The convergence weakens as abs(arg(z)) increases.  Also, in
   * the case of real order between 0 and 1, convergence weakens
   * as the order approaches 1.  For these reasons, we only use
   * this code for @f$ |\arg(z)| <= 3pi/4 @f$ and the convergence test
   * is performed only for order 2/3.
   *
   * The coding of this function is also influenced by the fact
   * that it will only be used for about @f$ 2 <= |z| <= 15 @f$.
   * Hence, certain considerations of overflow, underflow, and
   * loss of significance are unimportant for our purpose.
   *
   * @param[in] z The value for which the quantity @f$ E_\nu @f$ is to be computed.
   *  	     it is recommended that abs(z) not be too small and that
   *   	     @f$ |\arg(z)| <= 3pi/4 @f$.
   * @param[in] eps The maximum relative error allowable in the computed
   *   	     results. The relative error test is based on the
   *   	     comparison of successive iterates.
   * @param[out] _Kp1d3 The value computed for @f$ E_{+1/3}(z) @f$.
   * @param[out] _Kp2d3 The value computed for @f$ E_{+2/3}(z) @f$.
   *
   * @note In the worst case, say, @f$ z=2 @f$ and @f$ \arg(z) = 3pi/4 @f$,
   * 20 iterations should give 7 or 8 decimals of accuracy.
   */
  template<typename Tp>
    void
    airy_bessel_k(const std::complex<Tp>& z, Tp eps,
		    std::complex<Tp>& _Kp1d3,
		    std::complex<Tp>& _Kp2d3)
    {
      using cmplx = std::complex<Tp>;

      constexpr Tp s_an1i{  Tp{437} /  Tp{9}},
		    s_an2i{  Tp{425} /  Tp{9}},
		    s_p12i{  Tp{283} /  Tp{9}},
		    s_p22i{  Tp{295} /  Tp{9}},
		    s_p13i{-  Tp{25} / Tp{27}},
		    s_p23i{   Tp{35} / Tp{27}},
		    s_p11i{-Tp{2135} / Tp{27}},
		    s_p21i{-Tp{2195} / Tp{27}};

      constexpr std::complex<Tp> s_zone{1};
      constexpr int s_iter_max = 100;

      constexpr Tp
      s_f[8]
      { 144, 77, 62208, 95472, 17017, 65, 90288, 13585 };

      constexpr Tp
      s_phi[6]
      { 67, 91152, 12697, 79, 96336, 19633 };


      // Initialize polynomials for recurrence
      auto f10 = s_zone;
      auto f20 = s_zone;
      auto f11 = s_zone + s_f[0] * z / s_f[1];
      auto f12 = z * (s_f[2] * z + s_f[3]);
      f12 = s_zone + f12 / s_f[4];
      auto f21 = s_zone + s_f[0] * z / s_f[5];
      auto f22 = z * (s_f[2] * z + s_f[6]);
      f22 = s_zone + f22 / s_f[7];

      auto phi10 = s_zone;
      auto phi20 = s_zone;
      auto phi11 = cmplx((s_f[0] * std::real(z) + s_phi[0]) / s_f[1],
			     std::imag(f11));
      auto phi12 = z * (s_f[2] * z + s_phi[1]);
      phi12 = (phi12 + s_phi[2]) / s_f[4];
      auto phi21 = cmplx((s_f[0] * std::real(z) + s_phi[3]) / s_f[5],
			      std::imag(f21));
      auto phi22 = z * (s_f[2] * z + s_phi[4]);
      phi22 = (phi22 + s_phi[5]) / s_f[7];

      // Initialize for recursion
      auto ratold = phi22 / f22;
      auto rat1 = phi12 / f12;
      auto delt = Tp{32};
      auto an1 = s_an1i;
      auto an2 = s_an2i;
      auto p11 = s_p11i;
      auto p12 = s_p12i;
      auto p13 = s_p13i;
      auto p21 = s_p21i;
      auto p22 = s_p22i;
      auto p23 = s_p23i;
      auto eta = Tp{24};
      auto gamm = Tp{3};
      auto gam = Tp{5};
      auto q = Tp{16} * gam;

      // Loop until maximum iterations used or convergence.
      for (int i = 1; i < s_iter_max; ++i)
	{
	  // Evaluate next term in recurrence for order 1/3 polynomials.
	  auto qz = q * z;
	  auto fact1 = qz - p11;
	  auto fact2 = qz - p12;
	  auto f13 = fact1 * f12
		     + fact2 * f11
		     - p13 * f10;
	  f13 /= an1;
	  auto phi13 = fact1 * phi12
		       + fact2 * phi11
		       - p13 * phi10;
	  phi13 /= an1;

	  // Evaluate next term in recurrence for order 2/3 polynomials.
	  fact1 = qz - p21;
	  fact2 = qz - p22;
	  auto f23 = fact1 * f22
		     + fact2 * f21
		     - p23 * f20;
	  f23 /= an2;
	  auto phi23 = fact1 * phi22
		       + fact2 * phi21
		       - p23 * phi20;
	  phi23 /= an2;

	  // Check for convergence.
	  auto ratnew = phi23 / f23;
	  rat1 = phi13 / f13;

	  if (std::abs(ratnew - ratold) < eps * std::abs(ratnew))
	    {
	      // Convergence.
	      _Kp2d3 = ratnew;
	      _Kp1d3 = phi13 / f13;
	      return;
	    }

	  // Prepare for next iteration.
	  ratold = ratnew;
	  f20 = f21;
	  f21 = f22;
	  f22 = f23;
	  phi20 = phi21;
	  phi21 = phi22;
	  phi22 = phi23;
	  f10 = f11;
	  f11 = f12;
	  f12 = f13;
	  phi10 = phi11;
	  phi11 = phi12;
	  phi12 = phi13;
	  delt = delt + Tp{24};
	  p12 = p12 + delt;
	  p22 = p22 + delt;
	  eta += Tp{8};
	  an1 += eta;
	  an2 += eta;
	  auto anm1 = an1 - delt - Tp{16};
	  auto anm2 = an2 - delt - Tp{16};
	  gamm = gam;
	  gam += Tp{2};
	  p23 = -gam / gamm;
	  p13 = p23 * anm1;
	  p23 = p23 * anm2;
	  p11 = -an1 - p12 - p13;
	  p21 = -an2 - p22 - p23;
	  q = Tp{16} * gam;
	}

      throw std::runtime_error("airy_bessel_k: maximum iterations exceeded");

      return;
    }


  /**
   * @brief This function computes rational approximations
   * to the hypergeometric functions related to the modified Bessel
   * functions of orders @f$ \nu = +-1/3 @f$ and @f$ \nu +- 2/3 @f$.  That is,
   * @f$ A(z)/B(z) @f$, Where @f$ A(z) @f$ and @f$ B(z) @f$ are
   * cubic polynomials with real coefficients, approximates
   * @f[
   *  \frac{\Gamma(\nu+1)}{(z/2)^nu}I_\nu(z) = _0F_1 (;\nu+1;z^2/4),
   * @f]
   * where the function on the right is a
   * confluent hypergeometric limit function.  For @f$ |z| <= 1/4 @f$ and
   * @f$ |arg(z)| <= pi/2 @f$, the approximations are accurate to
   * about 16 decimals.
   *
   * For further details including the four term recurrence
   * relation satisfied by the numerator and denominator poly-
   * nomials in the higher order approximants, see
   *
   * Luke, Y.L., Mathematical Functions and their Approximations,
   * Academic Press, pp 361-363, 1975.
   *
   * An asymptotic expression for the error is given as well as
   * other useful expressions in the event one wants to extend
   * this function to incorporate higher order approximants.
   *
   * Note also that for speed and since this function is called by another,
   * checks that are not absolutely necessary are not made.
   *
   * @param[in] z The argument at which the hypergeometric given
   *  		above is to be evaluated.  Since the approximation
   *  		is of fixed order, @f$ |z| @f$ must be small to ensure
   *  		sufficient accuracy of the computed results.
   * @param[out] _Ai  The Airy function @f$ Ai(z) @f$.
   * @param[out] _Aip The Airy function derivative @f$ Ai'(z) @f$.
   * @param[out] _Bi  The Airy function @f$ Bi(z) @f$.
   * @param[out] _Bip The Airy function derivative @f$ Bi'(z) @f$.
   */
  template<typename Tp>
    void
    airy_hyperg_rational(const std::complex<Tp>& z,
			   std::complex<Tp>& _Ai,
			   std::complex<Tp>& _Aip,
			   std::complex<Tp>& _Bi,
			   std::complex<Tp>& _Bip)
    {
      using cmplx = std::complex<Tp>;

      constexpr cmplx s_zone{1};

      constexpr Tp s_ap1d3[4]{  81, 32400,  2585520,  37920960};
      constexpr Tp s_bp1d3[4]{ -35,  5040,  -574560,  37920960};
      constexpr Tp s_am1d3[4]{  81, 22680,  1156680,   7711200};
      constexpr Tp s_bm1d3[4]{ -10,  1260,  -128520,   7711200};
      constexpr Tp s_ap2d3[4]{ 162, 75735,  7270560, 139352400};
      constexpr Tp s_bp2d3[4]{-110, 16830, -2019600, 139352400};
      constexpr Tp s_am2d3[4]{ 162, 36855,  1415232,   4481568};
      constexpr Tp s_bm2d3[4]{  -7,	819,   -78624,   4481568};
      constexpr Tp s_Ai0{3.550280538878172392600631860041831763980e-1L};
      constexpr Tp s_Aip0{-2.588194037928067984051835601892039634793e-1L};
      constexpr Tp s_Bi0{6.149266274460007351509223690936135535960e-1L};
      constexpr Tp s_Bip0{4.482883573538263579148237103988283908668e-1L};

      // Check to see if z^3 will underflow and act accordingly.
      auto zzz = z * z * z;

      // The confluent hypergeometric limit functions related to
      // the modified Bessel functions of order +1/3, -1/3, +2/3, and -2/3
      // respectively.
      std::complex<Tp> _Fp1d3, _Fm1d3, _Fp2d3, _Fm2d3;
      if (std::abs(zzz) < Tp{10} * emsr::min<Tp>())
	{
	  _Fp1d3  = s_zone;
	  _Fm1d3 = s_zone;
	  _Fp2d3  = s_zone;
	  _Fm2d3 = s_zone;
	}
      else
	{
	  auto r = Tp{2} * std::real(zzz);
	  auto s = std::norm(zzz);

	  // The following polynomial evaluations are done using
	  // a modified of Horner's rule which exploits the fact that
	  // the polynomial coefficients are all real.
	  // The algorithm is discussed in detail in:
	  // Knuth, D. E., The Art of Computer Programming: Seminumerical
	  // Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
	  //
	  // If n is the degree of the polynomial, n - 3 multiplies are
	  // saved and 4 * n - 6 additions are saved.
	  auto horner
	  {
	    [r, s, zzz](const auto (&s_c)[4])
	    {
	      auto aa = s_c[0];
	      auto t  = s * aa;
	      aa = s_c[1] + r * aa;
	      auto bb = s_c[2] - t;
	      t  = s * aa;
	      aa = bb + r * aa;
	      bb = s_c[3] - t;
	      return aa * zzz + bb;
	    }
	  };

	  _Fp1d3 = horner(s_ap1d3) / horner(s_bp1d3);
	  _Fm1d3 = horner(s_am1d3) / horner(s_bm1d3);
	  _Fp2d3 = horner(s_ap2d3) / horner(s_bp2d3);
	  _Fm2d3 = horner(s_am2d3) / horner(s_bm2d3);
	}

      _Ai = s_Ai0 * _Fm1d3 + s_Aip0 * z * _Fp1d3;
      _Aip = s_Ai0 * z * z * _Fp2d3 / Tp{2} + s_Aip0 * _Fm2d3;
      _Bi = s_Bi0 * _Fm1d3 + s_Bip0 * z * _Fp1d3;
      _Bip = s_Bi0 * z * z * _Fp2d3 / Tp{2} + s_Bip0 * _Fm2d3;

      return;
    }


  /**
   * @brief This function computes the Airy function @f$ Ai(z) @f$
   * and its first derivative in the complex z-plane.
   *
   * The algorithm used exploits numerous representations of the
   * Airy function and its derivative.
   * The representations are recorded here for reference:
   *
   * @f[
   * \mbox{ (1a)    } Ai(z) = \frac{\sqrt{z}}{3}
   *                          \left[I_{-1/3}(\zeta) - I_{1/3}(\zeta)\right]
   * @f]
   * @f[
   * \mbox{ (1b)    } Bi(z) = \sqrt{\frac{z}{3}}
   *                          \left[I_{-1/3}(\zeta) + I_{1/3}(\zeta)\right]
   * @f]
   * @f[
   * \mbox{ (2)     } Ai(z) = \frac{\sqrt{z/3}}{\pi} K_{1/3}(\zeta)
   *    	            = \frac{2^{2/3}3^{-5/6}}{\sqrt(\pi)}
   *  	                      z \exp(-\zeta) U(5/6; 5/3; 2 \zeta)
   * @f]
   * @f[
   * \mbox{ (3a)    } Ai(-z) = \frac{\sqrt{z}}{3}
   *                           \left[J_{-1/3}(\zeta) + J_{1/3}(\zeta)\right]
   * @f]
   * @f[
   * \mbox{ (3b)    } Bi(-z) = \sqrt{\frac{z}{3}}
   *                           \left[J_{-1/3}(\zeta) - J_{1/3}(\zeta)\right]
   * @f]
   *
   * @f[
   * \mbox{ (4a)    } Ai'(z) = \frac{z}{3}
   *                           \left[I_{2/3}(\zeta) - I_{-2/3}(\zeta)\right]
   * @f]
   * @f[
   * \mbox{ (4b)    } Bi'(z) = \frac{z}{\sqrt{3}}
   *                           \left[I_{-2/3}(\zeta) + I_{2/3}(\zeta)\right]
   * @f]
   * @f[
   * \mbox{ (5a)    } Ai'(z) = -\frac{z}{\pi\sqrt{3}} K_{2/3}(\zeta)
   *  	                     = -\frac{4^{2/3}3^{-7/6}}{\sqrt{\pi}}
   *                         	z^2 \exp(-\zeta) U(7/6; 7/3; 2 \zeta)
   * @f]
   * @f[
   * \mbox{ (6a)    } Ai'(-z) = \frac{z}{3}
   *                 \left[J_{2/3}(\zeta) - J_{-2/3}(\zeta)\right]
   * @f]
   * @f[
   * \mbox{ (6b)    } Bi'(-z) = \frac{z}{\sqrt{3}}
   *                 \left[J_{-2/3}(\zeta) + J_{2/3}(\zeta)\right]
   * @f]
   * Where @f$ \zeta = - \frac{2}{3}z^{3/2} @f$ and @f$ U(a;b;z) @f$ is
   * the confluent hypergeometric function defined in
   *
   * @see Stegun, I. A. and Abramowitz, M., Handbook of Mathematical Functions,
   *  Natl. Bureau of Standards, AMS 55, pp 504-515, 1964.
   *
   * The asymptotic expansions derivable from these representations and
   * Hankel's asymptotic expansions for the Bessel functions are used for
   * large modulus of z.  The implementation has taken advantage of the
   * error bounds given in
   *
   * @see Olver, F. W. J., Error Bounds for Asymptotic Expansions,
   * with an Application to Cylinder Functions of Large Argument,
   * in Asymptotic Solutions of Differential Equations (Wilcox, Ed.),
   * Wiley and Sons, pp 163-183, 1964
   *
   * @see Olver, F. W. J., Asymptotics and Special Functions, Academic Press,
   * pp 266-268, 1974.
   *
   * For small modulus of z, a rational approximation is used.  This
   * approximant is derived from
   *
   * Luke, Y. L., Mathematical Functions and their Approximations,
   * Academic Press, pp 361-363, 1975.
   *
   * The identities given below are for Bessel functions of the first
   * kind in terms of modified Bessel functions of the first kind are
   * also used with the rational approximant.
   *
   * For moderate modulus of z, three techniques are used.  Two use
   * a backward recursion algorithm with (1), (3), (4), and (6). The
   * third uses the confluent hypergeometric representations given by
   * (2) and (5).  The backward recursion algorithm generates values of
   * the modified Bessel functions of the first kind of orders + or -
   * 1/3 and + or - 2/3 for z in the right half plane.  Values for
   * the corresponding Bessel functions of the first kind are recovered
   * via the identities
   * @f[
   *      J_\nu(z) = e^{i\nu\pi/2} I_\nu(z e^{-i\pi/2}),
   *  	  0 <= arg(z) <= \pi/2
   * @f]
   * and
   * @f[
   *      J_\nu(z) = e^{-\nu i\pi/2} I_\nu(z e^{i\pi/2}),
   *  	 -\pi/2 <= arg(z) < 0.
   * @f]
   *
   * The particular backward recursion algorithm used is discussed in 
   *
   * @see Olver, F. W. J, Numerical solution of second-order linear
   * difference equations, NBS J. Res., Series B, VOL 71B,
   * pp 111-129, 1967.
   *
   * @see Olver, F. W. J. and Sookne, D. J., Note on backward recurrence
   * algorithms, Math. Comp. Vol 26, No. 120, pp 941-947,
   * Oct. 1972
   *
   * @see Sookne, D. J., Bessel Functions I and J of Complex Argument and
   * Integer Order, NBS J. Res., Series B, Vol 77B, Nos 3& 4,
   * pp 111-113, July-December, 1973. 
   *
   * The following paper was also useful
   *
   * @see Cody, W. J., Preliminary report on software for the modified
   * Bessel functions of the first kind, Applied Mathematics
   * Division, Argonne National Laboratory, Tech. Memo. no. 357.
   *
   * A backward recursion algorithm is also used to compute the
   * confluent hypergeometric function.  The recursion relations
   * and a convergence theorem are given in
   *
   * @see Wimp, J., On the computation of Tricomi's psi function, Computing,
   * Vol 13, pp 195-203, 1974.
   *
   * @param[in] z   The argument at which the Airy function and its derivative
   *  	     are computed.
   * @param[in] eps Relative error required.  Currently, eps is used only
   *		     in the backward recursion algorithms.
   * @param[out] _Ai  The value computed for Ai(z).
   * @param[out] _Aip The value computed for Ai'(z).
   * @param[out] _Bi  The value computed for Bi(z).
   * @param[out] _Bip The value computed for Bi'(z).
   */
  template<typename Tp>
    void
    airy(const std::complex<Tp>& z, Tp eps,
           std::complex<Tp>& _Ai,
           std::complex<Tp>& _Aip,
           std::complex<Tp>& _Bi,
           std::complex<Tp>& _Bip)
    {
      using cmplx = std::complex<Tp>;

      static constexpr Tp
	s_sqrt3 = emsr::math_constants<Tp>::root_3;
      static constexpr cmplx s_j{0, 1},
	s_eppid6{s_sqrt3 / Tp{2},  Tp{0.5L}},
	s_empid6{s_sqrt3 / Tp{2}, -Tp{0.5L}},
	s_eppid3{Tp{0.5L},  s_sqrt3 / Tp{2}},
	s_empid3{Tp{0.5L}, -s_sqrt3 / Tp{2}};
      static constexpr Tp
	s_2d3{Tp{2} / Tp{3}},
	s_rsqpi{2.820947917738781434740397257803862929219e-01L}, // 1/(2sqrt(pi))
	s_pi = emsr::math_constants<Tp>::pi;
      static constexpr Tp s_small{0.25}, s_big{15};

      auto absz = std::abs(z);
      // Check size of abs(z) and select appropriate methods.
      if (absz < s_big)
	{
	  // Moderate or small abs(z)
	  // Check argument for right or left half plane.
	  if (std::real(z) >= Tp{0})
	    {
	      // Argument in closed right half plane. Compute zeta as defined
	      // in the representations in terms of Bessel functions.
	      auto sqrtz = std::sqrt(z);
	      auto zeta = s_2d3 * z * sqrtz;

	      // Check for abs(z) too large for accuracy of
	      // representations (1) and (4).
	      if (absz >= Tp{2})
		{
		  // Use rational approximation for modified Bessel functions
		  // of orders 1/3 and 2/3.
		  airy_bessel_k(zeta, eps, _Ai, _Aip/* , _Bi, _Bip*/); // FIXME
		  // Recover Ai(z) and Ai'(z).
		  auto p1d4c = std::sqrt(sqrtz);
		  zeta = s_rsqpi * std::exp(-zeta);
		  _Ai *= zeta / p1d4c;
		  _Aip *= -zeta * p1d4c;
		  // FIXME: _Bi *= 
		  // FIXME: _Bip *= 
		}
	      else
		{
		  // Check for abs(z) small enough for rational approximation.
		  if (absz <= s_small)
		    {
		      // Use rational approximation along with (1) and (4).
		      airy_hyperg_rational(z, _Ai, _Aip, _Bi, _Bip);
		    }
		  else
		    {
		      // Use Miller's backward recurrence along with (1), (4).
		      cmplx _Ip1d3, _Im1d3, _Ip2d3, _Im2d3;
		      airy_bessel_i(zeta, eps,
				      _Ip1d3, _Im1d3, _Ip2d3, _Im2d3);
		      // Recover Ai(z) and Ai'(z), Bi(z) and Bi'(z).
		      _Ai = sqrtz * (_Im1d3 - _Ip1d3) / Tp{3};
		      _Aip = z * (_Ip2d3 - _Im2d3) / Tp{3};
		      _Bi = sqrtz * (_Im2d3 * _Ip2d3) / s_sqrt3;
		      _Bip = z * (_Im2d3 + _Ip2d3) / s_sqrt3;
		    }
		}
	    }
	  else
	    {
	      // Argument lies in left half plane.  Compute zeta as defined
	      // in the representations in terms of Bessel functions.
	      auto sqrtz = std::sqrt(-z);
	      auto zeta = -s_2d3 * z * sqrtz;
	      // Set up arguments to recover Bessel functions of the first kind
	      // in (3) and (6).
	      cmplx z2zeta, p1d3f, m1d3f, p2d3f, m2d3f;
	      if (std::imag(zeta) >= Tp{0})
		{
		  // Argument lies in the second quadrant.
		  z2zeta = -s_j * zeta;
		  p1d3f = s_eppid6;
		  m1d3f = s_empid6;
		  p2d3f = s_eppid3;
		  m2d3f = s_empid3;
		}
	      else
		{
		  // Argument lies in the third quadrant.
		  z2zeta = s_j * zeta;
		  p1d3f = s_empid6;
		  m1d3f = s_eppid6;
		  p2d3f = s_empid3;
		  m2d3f = s_eppid3;
		}

	      // Use approximation depending on size of z.
	      if (absz <= s_small)
		{
		  // Use rational approximation.
		  zeta = -z;
		  airy_hyperg_rational(z, _Ai, _Aip, _Bi, _Bip);
		}
	      else
		{
		  // Use Miller's backward recurrence.
		  cmplx _Ip1d3, _Im1d3, _Ip2d3, _Im2d3;
		  airy_bessel_i(z2zeta, eps,
				  _Ip1d3, _Im1d3, _Ip2d3, _Im2d3);
		  // Recover Ai(z) and Ai'(z).
		  _Ai = sqrtz * (m1d3f * _Im1d3 + p1d3f * _Ip1d3) / Tp{3};
		  _Aip = z * (m2d3f * _Im2d3 - p2d3f * _Ip2d3) / Tp{3};
		  // FIXME: _Bi = sqrtz * (m2d3f * _Im2d3 + p2d3f * _Ip2d3) / s_sqrt3;
		  // FIXME: _Bip = z * (I_{2/3}(\zeta) + I_{-2/3}(\zeta)) / s_sqrt3
		}
	    }
	}
      else
	{ // abs(z) is large...
	  // Check arg(z) to see which asymptotic form is appropriate.
	  if (std::abs(std::arg(z)) < Tp{2} * s_pi / Tp{3})
	    airy_asymp_absarg_ge_pio3(z, _Ai, _Aip/* , _Bi, _Bip*/); // FIXME
	  else
	    airy_asymp_absarg_lt_pio3(z, _Ai, _Aip/* , _Bi, _Bip*/); // FIXME
	}
      return;
  }


  /**
   * @brief  Return the complex Airy Ai function.
   */
  template<typename Tp>
    std::complex<Tp>
    airy_ai(std::complex<Tp> z)
    {
      std::complex<Tp> _Ai, _Aip, _Bi, _Bip;
      airy(z, _Ai, _Aip, _Bi, _Bip);
      return _Ai;
    }


  /**
   * @brief  Return the complex Airy Bi function.
   */
  template<typename Tp>
    std::complex<Tp>
    airy_bi(std::complex<Tp> z)
    {
      std::complex<Tp> _Ai, _Aip, _Bi, _Bip;
      airy(z, _Ai, _Aip, _Bi, _Bip);
      return _Bi;
    }
} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_SF_AIRY_TCC
