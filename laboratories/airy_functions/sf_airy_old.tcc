// TR29124 math special functions -*- C++ -*-

// Copyright (C) 2016-2018 Free Software Foundation, Inc.
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

/** @file bits/sf_airy.tcc
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
   * @param[in]  __z Complex input variable set equal to the value at which
   *  	    @f$ Ai(z) @f$ and @f$ Bi(z) @f$ and their derivative are evaluated.
   *  	    This function assumes @f$ |z| > 15 @f$ and @f$ |arg(z)| < 2\pi/3 @f$.
   * @param[inout] _Ai  The value computed for @f$ Ai(z) @f$.
   * @param[inout] _Aip The value computed for @f$ Ai'(z) @f$.
   * @param[in]    __sign  The sign of the series terms amd exponent.
   *                    The default (-1) gives the Airy Ai functions
   *                    for @f$ |arg(z)| < \pi @f$.
   *                    The value +1 gives the Airy Bi functions
   *                    for @f$ |arg(z)| < \pi/3 @f$.
   */
  template<typename _Tp>
    void
    __airy_asymp_absarg_ge_pio3(std::complex<_Tp> __z,
				std::complex<_Tp>& _Ai,
				std::complex<_Tp>& _Aip,
				int __sign = -1)
    {
      constexpr _Tp _S_2d3   = _Tp{2} / _Tp{3};
      // 1/(2 sqrt(pi)))
      constexpr _Tp _S_pmhd2 = _Tp{2.820947917738781434740397257803862929219e-01L};
      constexpr int _S_ncoeffs = 15;
      constexpr int _S_numnterms = 5;
      constexpr int _S_nterms[5]{ _S_ncoeffs, 12, 11, 11, 9 };

      // Coefficients for the expansion.
      constexpr _Tp
      _S_u[_S_ncoeffs]
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

      constexpr _Tp
      _S_v[_S_ncoeffs]
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
      auto __z1d4 = std::sqrt(__z);
      auto __zetam = __z * __z1d4;
      __zetam *= _S_2d3;
      __z1d4 = std::sqrt(__z1d4);

      // Compute outer factors in the expansions.
      auto __zoutpr = std::exp(_Tp(__sign) * __zetam);
      __zoutpr *= _S_pmhd2;
      auto __zout = __zoutpr / __z1d4;
      __zoutpr *= -__z1d4;

      // Determine number of terms to use.
      auto __nterm = _S_nterms[std::min(_S_numnterms - 1,
					(int(std::abs(__z)) - 10) / 5)];
      // Initialize for modified Horner's rule evaluation of sums.
      // It is assumed that at least three terms are used.
      __zetam = _Tp(__sign) / __zetam;
      auto __r = _Tp{2} * std::real(__zetam);
      auto __s = std::norm(__zetam);
      auto __index = _S_ncoeffs - __nterm;// + 1;
      auto __al = _S_u[__index];
      auto __alpr = _S_v[__index];
      ++__index;
      auto __be = _S_u[__index];
      auto __bepr = _S_v[__index];
      ++__index;

      for (int __k = __index; __k < _S_ncoeffs; ++__k)
	{
	  auto __term = __s * __al;
	  __al = __be + __r * __al;
	  __be = _S_u[__k] - __term;
	  __term = __s * __alpr;
	  __alpr = __bepr + __r * __alpr;
	  __bepr = _S_v[__k] - __term;
	}

      _Ai = __zout * __al * __zetam + __be;
      _Aip = __zoutpr * __alpr * __zetam + __bepr;

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
   * @param[in] __z  The value at which the Airy function and its derivative
   *              are evaluated.
   * @param[out] _Ai  The computed value of the Airy function @f$ Ai(z) @f$.
   * @param[out] _Aip The computed value of the Airy function derivative @f$ Ai'(z) @f$.
   */
  template<typename _Tp>
    void
    __airy_asymp_absarg_lt_pio3(std::complex<_Tp> __z,
				std::complex<_Tp>& _Ai,
				std::complex<_Tp>& _Aip)
    {
      constexpr _Tp _S_2d3 {_Tp{2} / _Tp{3}};
      constexpr _Tp _S_9d4 {_Tp{9} / _Tp{4}};
      constexpr _Tp _S_pimh{5.641895835477562869480794515607725858438e-01};
      constexpr _Tp _S_pid4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;

      constexpr std::complex<_Tp> _S_zone{1};
      constexpr int _S_ncoeffs = 9;
      constexpr int _S_numnterms = 5;
      constexpr int _S_nterms[_S_numnterms]{ _S_ncoeffs, 7, 6, 6, 5 };

      // coefficients for the expansion.
      constexpr _Tp
      _S_u_cos[_S_ncoeffs]
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
      constexpr _Tp
      _S_u_sin[_S_ncoeffs]
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

      constexpr _Tp
      _S_v_sin[_S_ncoeffs]
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
      constexpr _Tp
      _S_v_cos[_S_ncoeffs]
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
      __z = -__z;
      // Compute zeta and z^(1/4).
      auto __z1d4 = std::sqrt(__z);
      auto __zeta = __z * __z1d4;
      __zeta *= _S_2d3;
      __z1d4 = std::sqrt(__z1d4);

      // Compute sine and cosine factors in the expansions.
      auto __zetaarg = __zeta + _S_pid4;
      auto __sinzeta = std::sin(__zetaarg);
      auto __coszeta = std::cos(__zetaarg);

      // Determine number of terms to use.
      auto __nterm = _S_nterms[std::min(_S_numnterms - 1,
					(int(std::abs(__z)) - 10) / 5)];
      // Initialize for modified Horner's rule evaluation of sums
      // it is assumed that at least three terms are used.
      __z = std::pow(_S_zone / __z, _Tp(3));
      __z *= _S_9d4;
      auto __r = -2 * std::real(__z);
      auto __s = std::norm(__z);
      auto __index = _S_ncoeffs - __nterm;

      auto __als = _S_u_sin[__index];
      auto __alc = _S_u_cos[__index];
      auto __alprs = _S_v_sin[__index];
      auto __alprc = _S_v_cos[__index];
      ++__index;

      auto __bes = _S_u_sin[__index];
      auto __bec = _S_u_cos[__index];
      auto __beprs = _S_v_sin[__index];
      auto __beprc = _S_v_cos[__index];
      ++__index;

      // Loop until components contributing to sums are computed.
      for (int __k = __index; __k < _S_ncoeffs; ++__k)
	{
	  auto __term = __s * __als;
	  __als = __bes + __r * __als;
	  __bes = _S_u_sin[__k] - __term;
	  __term = __s * __alc;
	  __alc = __bec + __r * __alc;
	  __bec = _S_u_cos[__k] - __term;
	  __term = __s * __alprs;
	  __alprs = __beprs + __r * __alprs;
	  __beprs = _S_v_sin[__k] - __term;
	  __term = __s * __alprc;
	  __alprc = __beprc + __r * __alprc;
	  __beprc = _S_v_cos[__k] - __term;
	}

      // Complete evaluation of the Airy functions.
      __zeta = _S_zone / __zeta;
      _Ai = __sinzeta * __als * __z + __bes
           - __zeta * __coszeta * __alc * __z + __bec;
      _Ai *= _S_pimh / __z1d4;
      _Aip = __coszeta * __alprc * __z + __beprc
	   + __zeta * __sinzeta * __alprs * __z + __beprs;
      _Aip *= -_S_pimh * __z1d4;

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
   * @param[in]  __z     The argument of the modified Bessel functions.
   * @param[in]  __eps   The maximum relative error required in the results.
   * @param[out] _Ip1d3 The value of @f$ I_(+1/3)(z) @f$.
   * @param[out] _Im1d3 The value of @f$ I_(-1/3)(z) @f$.
   * @param[out] _Ip2d3 The value of @f$ I_(+2/3)(z) @f$.
   * @param[out] _Im2d3 The value of @f$ I_(-2/3)(z) @f$.
   */
  template<typename _Tp>
    void
    __airy_bessel_i(const std::complex<_Tp>& __z, _Tp __eps,
		    std::complex<_Tp>& _Ip1d3,
		    std::complex<_Tp>& _Im1d3,
		    std::complex<_Tp>& _Ip2d3,
		    std::complex<_Tp>& _Im2d3)
    {
      using __cmplx = std::complex<_Tp>;

      constexpr __cmplx _S_zero{0}, _S_zone{1};
      constexpr _Tp
    	_S_1d3  {_Tp{1} / _Tp{3}}, _S_2d3  {_Tp{2} / _Tp{3}},
    	_S_4d3  {_Tp{4} / _Tp{3}}, _S_5d3  {_Tp{5} / _Tp{3}},
    	_S_8d3  {_Tp{8} / _Tp{3}}, _S_10d3 {_Tp{10} / _Tp{3}},
    	_S_14d3 {_Tp{14} / _Tp{3}}, _S_16d3 {_Tp{16} / _Tp{3}},
    	_S_gamma4d3{8.929795115692492112185643136582258813769e-1L},
	_S_gamma5d3{9.027452929509336112968586854363425236809e-1L},
    	_S_sqrt2 = __gnu_cxx::__math_constants<_Tp>::__root_2;

      // Compute 1/z for use in recurrence for speed and abs(z).
      auto __1dz = __safe_div(_Tp{1}, __z);

      // Initialize for forward recursion based on order 2/3.
      int __n = 0;
      auto __d2n = _Tp(2 * __n) + _S_4d3;
      auto __plast2 = _S_zone;
      auto __p2 = __d2n * __1dz;

      auto __weak_test = 20 * _S_sqrt2 / __eps;
      bool __converged = false;

      while (true)
	{
	  do
	    {
    	      // Update n dependent quantities.
    	      ++__n;
    	      __d2n += 2;
    	      // Interchange values.
    	      auto __pold2 = __plast2;
    	      __plast2 = __p2;
    	      // Recur forward one step.
    	      __p2 = __1dz * __d2n * __plast2 + __pold2;
	    }
	  while (__l1_norm(__p2) < __weak_test);

	  // If strong convergence, then weak and strong convergence.
	  if (__converged)
    	    break;

	  // Calculate strong convergence test.
	  // See the Olver and Sookne papers cited for details.
	  auto __rep2 = std::real(__p2);
	  auto __imp2 = std::imag(__p2);
	  auto __replast2 = std::real(__plast2);
	  auto __implast2 = std::imag(__plast2);
	  // Compute scale factor to avoid possible overflow.
	  auto __lamn = __linf_norm(__p2);
	  // Compute the kappa_n of strong convergence lemma.
	  auto __kapn = std::sqrt(((__rep2 / __lamn) * (__rep2 / __lamn)
    			       + (__imp2 / __lamn) * (__imp2 / __lamn))
    			     / ((__replast2 / __lamn) * (__replast2 / __lamn)
    			      + (__implast2 / __lamn) * (__implast2 / __lamn)));
	  // Compute quantity needed for lambda_n of strong convergence lemma
	  __lamn = _Tp(__n + 1) / std::abs(__z);
	  // Determine appropriate value for rho_n of lemma
	  if (__kapn + 1 / __kapn > 2 * __lamn)
    	    __kapn = __lamn + std::sqrt(__lamn * __lamn - 1);
	  // Compute test value - sqrt(2) multiple already included.
	  __weak_test *= std::sqrt(__kapn - 1 / __kapn);
	  // Set strong convergence test flag.
	  __converged = true;
	}

      // Prepare for backward recurrence for both orders 1/3 and 2/3.
      auto __rn = _Tp(__n);
      ++__n;
      __d2n = _Tp(2 * __n);
      auto __plast1 = _S_zero;
      __plast2 = _S_zero;
      // Carefully compute 1/p2 to avoid overflow in complex divide.
      auto __p1 = __safe_div(_Tp{1}, __p2);
      __p2 = __p1;
      // Set up n dependent parameters used in normalization sum.
      auto __rnpn1 = __rn + _S_1d3;
      auto __rnpn2 = __rn + _S_2d3;
      auto __rnp2n1 = (__rn - _Tp(1)) + _S_2d3;
      auto __rnp2n2 = (__rn - _Tp(1)) + _S_4d3;
      // Initialize normalization sum
      auto __fact1 = __rnpn1 * __rnp2n1 / __rn;
      auto __sum1 = __fact1 * __p1;
      auto __fact2 = __rnpn2 * __rnp2n2 / __rn;
      auto __sum2 = __fact2 * __p2;
      // Set ending loop index to correspond to k=1 term of the
      // normalization relationship.
      auto __nend = __n - 3;

      // If backward recurrence loop will be nontrivial...
      if (__nend > 0)
	{
	  // ...loop until backward recursion to k=1 term of normalization.
	  for (int __l = 1; __l <= __nend; ++__l)
	    {
    	      // Update n dependent quantities.
    	      --__n;
    	      __d2n -= 2;
    	      __fact1 = __d2n + _S_2d3;
    	      __fact2 = __d2n + _S_4d3;
    	      // Interchanges for order 1/3 recurrence.
    	      auto __pold1 = __plast1;
    	      __plast1 = __p1;
    	      // Recur back one step for order 1/3.
    	      __p1 = __1dz * __fact1 * __plast1 + __pold1;
    	      // Interchanges for order 2/3 recurrence
    	      auto __pold2 = __plast2;
    	      __plast2 = __p2;
    	      // Recur back one step for order 2/3.
    	      __p2 = __1dz * __fact2 * __plast2 + __pold2;
    	      // Update quantities for computing normalization sums.
    	      __rn -= _Tp(1);
    	      __rnpn1 = __rn + _S_1d3;
    	      __rnp2n1 = __rn - _S_1d3;
    	      __rnpn2 = __rn + _S_2d3;
    	      __rnp2n2 = __rnpn1;
    	      __fact1 = __rnp2n1 / __rn;
    	      __fact2 = __rnp2n2 / __rn;
    	      // Update normalization sums.
    	      __sum1 += __rnpn1 * __p1;
    	      __sum1 *= __fact1;
    	      __sum2 += __rnpn2 * __p2;
    	      __sum2 *= __fact2;
	    }
	}

      // Perform last two recurrence steps for order 1/3
      auto __pold1 = __plast1;
      __plast1 = __p1;
      __p1 = _S_14d3 * __plast1 * __1dz + __pold1;
      __sum1 += _S_4d3 * __p1;
      __pold1 = __plast1;
      __plast1 = __p1;
      __p1 = _S_8d3 * __plast1 * __1dz + __pold1;
      __sum1 = _Tp(2) * __sum1 + __p1;

      // Compute scale factor and scale results for order 1/3 case
      auto __zd2pow = std::pow(_Tp(0.5L) * __z, -_S_1d3);
      __pold1 = __zd2pow * std::exp(-__z);
      __sum1 *= _S_gamma4d3 * __pold1;
      __plast1 /= __sum1;
      _Ip1d3 = __p1 / __sum1;

      // Perform last two recurrence steps for order 2/3
      auto __pold2 = __plast2;
      __plast2 = __p2;
      __p2 = _S_16d3 * __plast2 * __1dz + __pold2;
      __sum2 += _S_5d3 * __p2;
      __pold2 = __plast2;
      __plast2 = __p2;
      __p2 = _S_10d3 * __plast2 * __1dz + __pold2;
      __sum2 = _Tp(2) * __sum2 + __p2;

      // Compute scale factor and scale results for order 2/3 case
      __sum2 *= _S_gamma5d3 * __zd2pow * __pold1;
      __plast2 /= __sum2;
      _Ip2d3 = __p2 / __sum2;

      // Recur back one step from order +1/3 to get order -2/3
      _Im2d3 = _S_2d3 * _Ip1d3 * __1dz + __plast1;

      // Recur back one step from order +2/3 to get order -1/3
      _Im1d3 = _S_4d3 * _Ip2d3 * __1dz + __plast2;

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
   * @param[in] __z The value for which the quantity @f$ E_\nu @f$ is to be computed.
   *  	     it is recommended that abs(z) not be too small and that
   *   	     @f$ |\arg(z)| <= 3pi/4 @f$.
   * @param[in] __eps The maximum relative error allowable in the computed
   *   	     results. The relative error test is based on the
   *   	     comparison of successive iterates.
   * @param[out] _Kp1d3 The value computed for @f$ E_{+1/3}(z) @f$.
   * @param[out] _Kp2d3 The value computed for @f$ E_{+2/3}(z) @f$.
   *
   * @note In the worst case, say, @f$ z=2 @f$ and @f$ \arg(z) = 3pi/4 @f$,
   * 20 iterations should give 7 or 8 decimals of accuracy.
   */
  template<typename _Tp>
    void
    __airy_bessel_k(const std::complex<_Tp>& __z, _Tp __eps,
		    std::complex<_Tp>& _Kp1d3,
		    std::complex<_Tp>& _Kp2d3)
    {
      using __cmplx = std::complex<_Tp>;

      constexpr _Tp _S_an1i{  _Tp{437} /  _Tp{9}},
		    _S_an2i{  _Tp{425} /  _Tp{9}},
		    _S_p12i{  _Tp{283} /  _Tp{9}},
		    _S_p22i{  _Tp{295} /  _Tp{9}},
		    _S_p13i{-  _Tp{25} / _Tp{27}},
		    _S_p23i{   _Tp{35} / _Tp{27}},
		    _S_p11i{-_Tp{2135} / _Tp{27}},
		    _S_p21i{-_Tp{2195} / _Tp{27}};

      constexpr std::complex<_Tp> _S_zone{1};
      constexpr int _S_iter_max = 100;

      constexpr _Tp
      _S_f[8]
      { 144, 77, 62208, 95472, 17017, 65, 90288, 13585 };

      constexpr _Tp
      _S_phi[6]
      { 67, 91152, 12697, 79, 96336, 19633 };


      // Initialize polynomials for recurrence
      auto __f10 = _S_zone;
      auto __f20 = _S_zone;
      auto __f11 = _S_zone + _S_f[0] * __z / _S_f[1];
      auto __f12 = __z * (_S_f[2] * __z + _S_f[3]);
      __f12 = _S_zone + __f12 / _S_f[4];
      auto __f21 = _S_zone + _S_f[0] * __z / _S_f[5];
      auto __f22 = __z * (_S_f[2] * __z + _S_f[6]);
      __f22 = _S_zone + __f22 / _S_f[7];

      auto __phi10 = _S_zone;
      auto __phi20 = _S_zone;
      auto __phi11 = __cmplx((_S_f[0] * std::real(__z) + _S_phi[0]) / _S_f[1],
			     std::imag(__f11));
      auto __phi12 = __z * (_S_f[2] * __z + _S_phi[1]);
      __phi12 = (__phi12 + _S_phi[2]) / _S_f[4];
      auto __phi21 = __cmplx((_S_f[0] * std::real(__z) + _S_phi[3]) / _S_f[5],
			      std::imag(__f21));
      auto __phi22 = __z * (_S_f[2] * __z + _S_phi[4]);
      __phi22 = (__phi22 + _S_phi[5]) / _S_f[7];

      // Initialize for recursion
      auto __ratold = __phi22 / __f22;
      auto __rat1 = __phi12 / __f12;
      auto __delt = _Tp{32};
      auto __an1 = _S_an1i;
      auto __an2 = _S_an2i;
      auto __p11 = _S_p11i;
      auto __p12 = _S_p12i;
      auto __p13 = _S_p13i;
      auto __p21 = _S_p21i;
      auto __p22 = _S_p22i;
      auto __p23 = _S_p23i;
      auto __eta = _Tp{24};
      auto __gamm = _Tp{3};
      auto __gam = _Tp{5};
      auto __q = _Tp{16} * __gam;

      // Loop until maximum iterations used or convergence.
      for (int __i = 1; __i < _S_iter_max; ++__i)
	{
	  // Evaluate next term in recurrence for order 1/3 polynomials.
	  auto __qz = __q * __z;
	  auto __fact1 = __qz - __p11;
	  auto __fact2 = __qz - __p12;
	  auto __f13 = __fact1 * __f12
		     + __fact2 * __f11
		     - __p13 * __f10;
	  __f13 /= __an1;
	  auto __phi13 = __fact1 * __phi12
		       + __fact2 * __phi11
		       - __p13 * __phi10;
	  __phi13 /= __an1;

	  // Evaluate next term in recurrence for order 2/3 polynomials.
	  __fact1 = __qz - __p21;
	  __fact2 = __qz - __p22;
	  auto __f23 = __fact1 * __f22
		     + __fact2 * __f21
		     - __p23 * __f20;
	  __f23 /= __an2;
	  auto __phi23 = __fact1 * __phi22
		       + __fact2 * __phi21
		       - __p23 * __phi20;
	  __phi23 /= __an2;

	  // Check for convergence.
	  auto __ratnew = __phi23 / __f23;
	  __rat1 = __phi13 / __f13;

	  if (std::abs(__ratnew - __ratold) < __eps * std::abs(__ratnew))
	    {
	      // Convergence.
	      _Kp2d3 = __ratnew;
	      _Kp1d3 = __phi13 / __f13;
	      return;
	    }

	  // Prepare for next iteration.
	  __ratold = __ratnew;
	  __f20 = __f21;
	  __f21 = __f22;
	  __f22 = __f23;
	  __phi20 = __phi21;
	  __phi21 = __phi22;
	  __phi22 = __phi23;
	  __f10 = __f11;
	  __f11 = __f12;
	  __f12 = __f13;
	  __phi10 = __phi11;
	  __phi11 = __phi12;
	  __phi12 = __phi13;
	  __delt = __delt + _Tp{24};
	  __p12 = __p12 + __delt;
	  __p22 = __p22 + __delt;
	  __eta += _Tp{8};
	  __an1 += __eta;
	  __an2 += __eta;
	  auto __anm1 = __an1 - __delt - _Tp{16};
	  auto __anm2 = __an2 - __delt - _Tp{16};
	  __gamm = __gam;
	  __gam += _Tp{2};
	  __p23 = -__gam / __gamm;
	  __p13 = __p23 * __anm1;
	  __p23 = __p23 * __anm2;
	  __p11 = -__an1 - __p12 - __p13;
	  __p21 = -__an2 - __p22 - __p23;
	  __q = _Tp{16} * __gam;
	}

      std::__throw_runtime_error(__N("__airy_bessel_k:"
				     " maximum iterations exceeded"));

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
   * @param[in] __z The argument at which the hypergeometric given
   *  		above is to be evaluated.  Since the approximation
   *  		is of fixed order, @f$ |z| @f$ must be small to ensure
   *  		sufficient accuracy of the computed results.
   * @param[out] _Ai  The Airy function @f$ Ai(z) @f$.
   * @param[out] _Aip The Airy function derivative @f$ Ai'(z) @f$.
   * @param[out] _Bi  The Airy function @f$ Bi(z) @f$.
   * @param[out] _Bip The Airy function derivative @f$ Bi'(z) @f$.
   */
  template<typename _Tp>
    void
    __airy_hyperg_rational(const std::complex<_Tp>& __z,
			   std::complex<_Tp>& _Ai,
			   std::complex<_Tp>& _Aip,
			   std::complex<_Tp>& _Bi,
			   std::complex<_Tp>& _Bip)
    {
      using __cmplx = std::complex<_Tp>;

      constexpr __cmplx _S_zone{1};

      constexpr _Tp _S_ap1d3[4]{  81, 32400,  2585520,  37920960};
      constexpr _Tp _S_bp1d3[4]{ -35,  5040,  -574560,  37920960};
      constexpr _Tp _S_am1d3[4]{  81, 22680,  1156680,   7711200};
      constexpr _Tp _S_bm1d3[4]{ -10,  1260,  -128520,   7711200};
      constexpr _Tp _S_ap2d3[4]{ 162, 75735,  7270560, 139352400};
      constexpr _Tp _S_bp2d3[4]{-110, 16830, -2019600, 139352400};
      constexpr _Tp _S_am2d3[4]{ 162, 36855,  1415232,   4481568};
      constexpr _Tp _S_bm2d3[4]{  -7,	819,   -78624,   4481568};
      constexpr _Tp _S_Ai0{3.550280538878172392600631860041831763980e-1L};
      constexpr _Tp _S_Aip0{-2.588194037928067984051835601892039634793e-1L};
      constexpr _Tp _S_Bi0{6.149266274460007351509223690936135535960e-1L};
      constexpr _Tp _S_Bip0{4.482883573538263579148237103988283908668e-1L};

      // Check to see if z^3 will underflow and act accordingly.
      auto __zzz = __z * __z * __z;

      // The confluent hypergeometric limit functions related to
      // the modified Bessel functions of order +1/3, -1/3, +2/3, and -2/3
      // respectively.
      std::complex<_Tp> _Fp1d3, _Fm1d3, _Fp2d3, _Fm2d3;
      if (std::abs(__zzz) < _Tp{10} * __gnu_cxx::__min<_Tp>())
	{
	  _Fp1d3  = _S_zone;
	  _Fm1d3 = _S_zone;
	  _Fp2d3  = _S_zone;
	  _Fm2d3 = _S_zone;
	}
      else
	{
	  auto __r = _Tp{2} * std::real(__zzz);
	  auto __s = std::norm(__zzz);

	  // The following polynomial evaluations are done using
	  // a modified of Horner's rule which exploits the fact that
	  // the polynomial coefficients are all real.
	  // The algorithm is discussed in detail in:
	  // Knuth, D. E., The Art of Computer Programming: Seminumerical
	  // Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
	  //
	  // If n is the degree of the polynomial, n - 3 multiplies are
	  // saved and 4 * n - 6 additions are saved.
	  auto __horner
	  {
	    [__r, __s, __zzz](const auto (&_S_c)[4])
	    {
	      auto __aa = _S_c[0];
	      auto __t  = __s * __aa;
	      __aa = _S_c[1] + __r * __aa;
	      auto __bb = _S_c[2] - __t;
	      __t  = __s * __aa;
	      __aa = __bb + __r * __aa;
	      __bb = _S_c[3] - __t;
	      return __aa * __zzz + __bb;
	    }
	  };

	  _Fp1d3 = __horner(_S_ap1d3) / __horner(_S_bp1d3);
	  _Fm1d3 = __horner(_S_am1d3) / __horner(_S_bm1d3);
	  _Fp2d3 = __horner(_S_ap2d3) / __horner(_S_bp2d3);
	  _Fm2d3 = __horner(_S_am2d3) / __horner(_S_bm2d3);
	}

      _Ai = _S_Ai0 * _Fm1d3 + _S_Aip0 * __z * _Fp1d3;
      _Aip = _S_Ai0 * __z * __z * _Fp2d3 / _Tp{2} + _S_Aip0 * _Fm2d3;
      _Bi = _S_Bi0 * _Fm1d3 + _S_Bip0 * __z * _Fp1d3;
      _Bip = _S_Bi0 * __z * __z * _Fp2d3 / _Tp{2} + _S_Bip0 * _Fm2d3;

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
   * @param[in] __z   The argument at which the Airy function and its derivative
   *  	     are computed.
   * @param[in] __eps Relative error required.  Currently, eps is used only
   *		     in the backward recursion algorithms.
   * @param[out] _Ai  The value computed for Ai(z).
   * @param[out] _Aip The value computed for Ai'(z).
   * @param[out] _Bi  The value computed for Bi(z).
   * @param[out] _Bip The value computed for Bi'(z).
   */
  template<typename _Tp>
    void
    __airy(const std::complex<_Tp>& __z, _Tp __eps,
           std::complex<_Tp>& _Ai,
           std::complex<_Tp>& _Aip,
           std::complex<_Tp>& _Bi,
           std::complex<_Tp>& _Bip)
    {
      using __cmplx = std::complex<_Tp>;

      static constexpr _Tp
	_S_sqrt3 = __gnu_cxx::__math_constants<_Tp>::__root_3;
      static constexpr __cmplx _S_j{0, 1},
	_S_eppid6{_S_sqrt3 / _Tp{2},  _Tp{0.5L}},
	_S_empid6{_S_sqrt3 / _Tp{2}, -_Tp{0.5L}},
	_S_eppid3{_Tp{0.5L},  _S_sqrt3 / _Tp{2}},
	_S_empid3{_Tp{0.5L}, -_S_sqrt3 / _Tp{2}};
      static constexpr _Tp
	_S_2d3{_Tp{2} / _Tp{3}},
	_S_rsqpi{2.820947917738781434740397257803862929219e-01L}, // 1/(2sqrt(pi))
	_S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      static constexpr _Tp _S_small{0.25}, _S_big{15};

      auto __absz = std::abs(__z);
      // Check size of abs(z) and select appropriate methods.
      if (__absz < _S_big)
	{
	  // Moderate or small abs(z)
	  // Check argument for right or left half plane.
	  if (std::real(__z) >= _Tp{0})
	    {
	      // Argument in closed right half plane. Compute zeta as defined
	      // in the representations in terms of Bessel functions.
	      auto __sqrtz = std::sqrt(__z);
	      auto __zeta = _S_2d3 * __z * __sqrtz;

	      // Check for abs(z) too large for accuracy of
	      // representations (1) and (4).
	      if (__absz >= _Tp{2})
		{
		  // Use rational approximation for modified Bessel functions
		  // of orders 1/3 and 2/3.
		  __airy_bessel_k(__zeta, __eps, _Ai, _Aip/* , _Bi, _Bip*/); // FIXME
		  // Recover Ai(z) and Ai'(z).
		  auto __p1d4c = std::sqrt(__sqrtz);
		  __zeta = _S_rsqpi * std::exp(-__zeta);
		  _Ai *= __zeta / __p1d4c;
		  _Aip *= -__zeta * __p1d4c;
		  // FIXME: _Bi *= 
		  // FIXME: _Bip *= 
		}
	      else
		{
		  // Check for abs(z) small enough for rational approximation.
		  if (__absz <= _S_small)
		    {
		      // Use rational approximation along with (1) and (4).
		      __airy_hyperg_rational(__z, _Ai, _Aip, _Bi, _Bip);
		    }
		  else
		    {
		      // Use Miller's backward recurrence along with (1), (4).
		      __cmplx _Ip1d3, _Im1d3, _Ip2d3, _Im2d3;
		      __airy_bessel_i(__zeta, __eps,
				      _Ip1d3, _Im1d3, _Ip2d3, _Im2d3);
		      // Recover Ai(z) and Ai'(z), Bi(z) and Bi'(z).
		      _Ai = __sqrtz * (_Im1d3 - _Ip1d3) / _Tp{3};
		      _Aip = __z * (_Ip2d3 - _Im2d3) / _Tp{3};
		      _Bi = __sqrtz * (_Im2d3 * _Ip2d3) / _S_sqrt3;
		      _Bip = __z * (_Im2d3 + _Ip2d3) / _S_sqrt3;
		    }
		}
	    }
	  else
	    {
	      // Argument lies in left half plane.  Compute zeta as defined
	      // in the representations in terms of Bessel functions.
	      auto __sqrtz = std::sqrt(-__z);
	      auto __zeta = -_S_2d3 * __z * __sqrtz;
	      // Set up arguments to recover Bessel functions of the first kind
	      // in (3) and (6).
	      __cmplx __z2zeta, __p1d3f, __m1d3f, __p2d3f, __m2d3f;
	      if (std::imag(__zeta) >= _Tp{0})
		{
		  // Argument lies in the second quadrant.
		  __z2zeta = -_S_j * __zeta;
		  __p1d3f = _S_eppid6;
		  __m1d3f = _S_empid6;
		  __p2d3f = _S_eppid3;
		  __m2d3f = _S_empid3;
		}
	      else
		{
		  // Argument lies in the third quadrant.
		  __z2zeta = _S_j * __zeta;
		  __p1d3f = _S_empid6;
		  __m1d3f = _S_eppid6;
		  __p2d3f = _S_empid3;
		  __m2d3f = _S_eppid3;
		}

	      // Use approximation depending on size of z.
	      if (__absz <= _S_small)
		{
		  // Use rational approximation.
		  __zeta = -__z;
		  __airy_hyperg_rational(__z, _Ai, _Aip, _Bi, _Bip);
		}
	      else
		{
		  // Use Miller's backward recurrence.
		  __cmplx _Ip1d3, _Im1d3, _Ip2d3, _Im2d3;
		  __airy_bessel_i(__z2zeta, __eps,
				  _Ip1d3, _Im1d3, _Ip2d3, _Im2d3);
		  // Recover Ai(z) and Ai'(z).
		  _Ai = __sqrtz * (__m1d3f * _Im1d3 + __p1d3f * _Ip1d3) / _Tp{3};
		  _Aip = __z * (__m2d3f * _Im2d3 - __p2d3f * _Ip2d3) / _Tp{3};
		  // FIXME: _Bi = __sqrtz * (__m2d3f * _Im2d3 + __p2d3f * _Ip2d3) / _S_sqrt3;
		  // FIXME: _Bip = __z * (I_{2/3}(\zeta) + I_{-2/3}(\zeta)) / _S_sqrt3
		}
	    }
	}
      else
	{ // abs(z) is large...
	  // Check arg(z) to see which asymptotic form is appropriate.
	  if (std::abs(std::arg(__z)) < _Tp{2} * _S_pi / _Tp{3})
	    __airy_asymp_absarg_ge_pio3(__z, _Ai, _Aip/* , _Bi, _Bip*/); // FIXME
	  else
	    __airy_asymp_absarg_lt_pio3(__z, _Ai, _Aip/* , _Bi, _Bip*/); // FIXME
	}
      return;
  }


  /**
   * @brief  Return the complex Airy Ai function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __airy_ai(std::complex<_Tp> __z)
    {
      std::complex<_Tp> _Ai, _Aip, _Bi, _Bip;
      __airy(__z, _Ai, _Aip, _Bi, _Bip);
      return _Ai;
    }


  /**
   * @brief  Return the complex Airy Bi function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __airy_bi(std::complex<_Tp> __z)
    {
      std::complex<_Tp> _Ai, _Aip, _Bi, _Bip;
      __airy(__z, _Ai, _Aip, _Bi, _Bip);
      return _Bi;
    }
} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_SF_AIRY_TCC