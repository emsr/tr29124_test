// Special functions -*- C++ -*-

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

/** @file emsr/sf_bessel.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
//   (1) Handbook of Mathematical Functions,
//       ed. Milton Abramowitz and Irene A. Stegun,
//       Dover Publications,
//       Section 9, pp. 355-434, Section 10 pp. 435-478
//   (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
//   (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//       W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//       2nd ed, pp. 240-245


#include <cmath>
#include <stdexcept>

  /**
   *   @brief Compute the gamma functions required by the Temme series
   *          expansions of @f$ N_\nu(x) @f$ and @f$ K_\nu(x) @f$.
   *   @f[
   *     \Gamma_1 = \frac{1}{2\mu}
   *                [\frac{1}{\Gamma(1 - \mu)} - \frac{1}{\Gamma(1 + \mu)}]
   *   @f]
   *   and
   *   @f[
   *     \Gamma_2 = \frac{1}{2}
   *                [\frac{1}{\Gamma(1 - \mu)} + \frac{1}{\Gamma(1 + \mu)}]
   *   @f]
   *   where @f$ -1/2 <= \mu <= 1/2 @f$ is @f$ \mu = \nu - N @f$ and @f$ N @f$.
   *   is the nearest integer to @f$ \nu @f$.
   *   The values of \f$ \Gamma(1 + \mu) \f$ and \f$ \Gamma(1 - \mu) \f$
   *   are returned as well.
   *
   *   The accuracy requirements on this are exquisite.
   *
   *   @param mu     The input parameter of the gamma functions.
   *   @param gam1   The output function \f$ \Gamma_1(\mu) \f$
   *   @param gam2   The output function \f$ \Gamma_2(\mu) \f$
   *   @param gampl  The output function \f$ \Gamma(1 + \mu) \f$
   *   @param gammi  The output function \f$ \Gamma(1 - \mu) \f$
   */
  template <typename Tp>
    void
    gamma_temme(const Tp mu,
		   Tp & gam1, Tp & gam2, Tp & gampl, Tp & gammi)
    {
      gampl = Tp{1} / std::::tgamma(Tp{1} + mu);
      gammi = Tp{1} / std::tgamma(Tp{1} - mu);

      if (std::abs(mu) < std::numeric_limits<Tp>::epsilon())
	gam1 = -Tp(GAMMA_E);
      else
	gam1 = (gammi - gampl) / (Tp{2} * mu);

      gam2 = (gammi + gampl) / (Tp{2});

      return;
    }


  /**
   *   @brief  Compute the Bessel @f$ J_\nu(x) @f$ and Neumann
   *           @f$ N_\nu(x) @f$ functions and their first derivatives
   *           @f$ J'_\nu(x) @f$ and @f$ N'_\nu(x) @f$ respectively.
   *           These four functions are computed together for numerical
   *           stability.
   *
   *   @param  nu  The order of the Bessel functions.
   *   @param  x   The argument of the Bessel functions.
   *   @param  _Jnu  The output Bessel function of the first kind.
   *   @param  _Nnu  The output Neumann function (Bessel function of the second kind).
   *   @param  _Jpnu  The output derivative of the Bessel function of the first kind.
   *   @param  _Npnu  The output derivative of the Neumann function.
   */
  template <typename Tp>
    void
    bessel_jn(const Tp nu, const Tp x,
		Tp & _Jnu, Tp & _Nnu, Tp & _Jpnu, Tp & _Npnu)
    {

      //if (std::isnan(nu) || std::isnan(x))
      //  return std::numeric_limits<Tp>::quiet_NaN();

      if (x == Tp{0})
	{
	  if (nu == Tp{0})
	    {
	      _Jnu = Tp{1};
	      _Nnu = -std::numeric_limits<Tp>::infinity();
	      _Jpnu = Tp{0};
	      _Npnu = std::numeric_limits<Tp>::infinity();
	    }
	  else
	    {
	      _Jnu = Tp{0};
	      _Nnu = -std::numeric_limits<Tp>::infinity();
	      //_Jpnu = ???
	      //_Npnu = ???
	    }
	  return;
	}

      const auto eps = std::numeric_limits<Tp>::epsilon();
      const auto fp_min = Tp(10) * std::numeric_limits<Tp>::min();
      const int max_iter = 15000;
      const auto x_min = Tp{2};

      if (x < Tp{0} || nu < Tp{0})
	throw std::runtime_error("Bad arguments in bessel_jn.");

      const int nl = (x < x_min
		    ? static_cast<int>(nu + Tp{0.5L})
		    : std::max(0, static_cast<int>(nu - x + Tp(1.5L))));

      const auto mu = nu - nl;
      const auto mu2 = mu * mu;
      const auto xi = Tp{1} / x;
      const auto xi2 = Tp{2} * xi;
      auto w = xi2 / Tp(PI);
      int isign = 1;
      auto h = nu * xi;
      if (h < fp_min)
	h = fp_min;
      auto b = xi2 * nu;
      auto d = Tp{0};
      auto c = h;
      int i;
      for (i = 1; i <= max_iter; ++i)
	{
	  b += xi2;
	  d = b - d;
	  if (std::abs(d) < fp_min)
	    d = fp_min;
	  c = b - Tp{1} / c;
	  if (std::abs(c) < fp_min)
	    c = fp_min;
	  d = Tp{1} / d;
	  const auto del = c * d;
	  h = del * h;
	  if (d < Tp{0})
	    isign = -isign;
	  if (std::abs(del - Tp{1}) < eps)
	    break;
	}
      if (i > max_iter)
	throw std::runtime_error( "Argument x too large in bessel_jn; "
				  "try asymptotic expansion." );
      auto _Jnul = isign * fp_min;
      auto _Jpnul = h * _Jnul;
      auto _Jnul1 = _Jnul;
      auto _Jpnu1 = _Jpnul;
      auto fact = nu * xi;
      for ( int l = nl; l >= 1; --l )
	{
	  const auto _Jnutemp = fact * _Jnul + _Jpnul;
	  fact -= xi;
	  _Jpnul = fact * _Jnutemp - _Jnul;
	  _Jnul = _Jnutemp;
	}
      if (_Jnul == Tp{0})
	_Jnul = eps;
      auto f = _Jpnul / _Jnul;
      Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (x < x_min)
	{
	  const auto x2 = x / Tp{2};
	  const auto pimu = Tp(PI) * mu;
	  auto fact = (std::abs(pimu) < eps
		      ? Tp{1} : pimu / std::sin(pimu));
	  auto d = -std::log(x2);
	  auto e = mu * d;
	  auto fact2 = (std::abs(e) < eps
		       ? Tp{1} : std::sinh(e) / e);
	  Tp gam1, gam2, gampl, gammi;
	  gamma_temme(mu, gam1, gam2, gampl, gammi);
	  auto ff = (Tp{2} / Tp(PI))
		   * fact * (gam1 * std::cosh(e) + gam2 * fact2 * d);
	  e = std::exp(e);
	  auto p = e / (Tp(PI) * gampl);
	  auto q = Tp{1} / (e * Tp(PI) * gammi);
	  const auto pimu2 = pimu / Tp{2};
	  auto fact3 = (std::abs(pimu2) < eps
		       ? Tp{1} : std::sin(pimu2) / pimu2 );
	  auto r = Tp(PI) * pimu2 * fact3 * fact3;
	  auto c = Tp{1};
	  d = -x2 * x2;
	  auto sum = ff + r * q;
	  auto sum1 = p;
	  for (i = 1; i <= max_iter; ++i)
	    {
	      ff = (i * ff + p + q) / (i * i - mu2);
	      c *= d / Tp(i);
	      p /= Tp(i) - mu;
	      q /= Tp(i) + mu;
	      const auto del = c * (ff + r * q);
	      sum += del; 
	      const auto del1 = c * p - i * del;
	      sum1 += del1;
	      if ( std::abs(del) < eps * (Tp{1} + std::abs(sum)) )
		break;
	    }
	  if ( i > max_iter )
	    throw std::runtime_error("Bessel y series failed to converge "
				     "in bessel_jn." );
	  _Nmu = -sum;
	  _Nnu1 = -sum1 * xi2;
	  _Npmu = mu * xi * _Nmu - _Nnu1;
	  _Jmu = w / (_Npmu - f * _Nmu);
	}
      else
	{
	  auto a = Tp{0.25L} - mu2;
	  auto q = Tp{1};
	  auto p = -xi / Tp{2};
	  auto br = Tp{2} * x;
	  auto bi = Tp{2};
	  auto fact = a * xi / (p * p + q * q);
	  auto cr = br + q * fact;
	  auto ci = bi + p * fact;
	  auto den = br * br + bi * bi;
	  auto dr = br / den;
	  auto di = -bi / den;
	  auto dlr = cr * dr - ci * di;
	  auto dli = cr * di + ci * dr;
	  auto temp = p * dlr - q * dli;
	  q = p * dli + q * dlr;
	  p = temp;
	  int i;
	  for (i = 2; i <= max_iter; ++i)
	    {
	      a += Tp(2 * (i - 1));
	      bi += Tp{2};
	      dr = a * dr + br;
	      di = a * di + bi;
	      if (std::abs(dr) + std::abs(di) < fp_min)
		dr = fp_min;
	      fact = a / (cr * cr + ci * ci);
	      cr = br + cr * fact;
	      ci = bi - ci * fact;
	      if (std::abs(cr) + std::abs(ci) < fp_min)
		cr = fp_min;
	      den = dr * dr + di * di;
	      dr /= den;
	      di /= -den;
	      dlr = cr * dr - ci * di;
	      dli = cr * di + ci * dr;
	      temp = p * dlr - q * dli;
	      q = p * dli + q * dlr;
	      p = temp;
	      if (std::abs(dlr - Tp{1}) + std::abs(dli) < eps)
		break;
	  }
	  if (i > max_iter)
	    throw std::runtime_error("Lentz's method failed in bessel_jn.");
	  const auto gam = (p - f) / q;
	  _Jmu = std::sqrt(w / ((p - f) * gam + q));

	  _Jmu = ::copysign(_Jmu, _Jnul);

	  _Nmu = gam * _Jmu;
	  _Npmu = (p + q / gam) * _Nmu;
	  _Nnu1 = mu * xi * _Nmu - _Npmu;
      }
      fact = _Jmu / _Jnul;
      _Jnu = fact * _Jnul1;
      _Jpnu = fact * _Jpnu1;
      for (i = 1; i <= nl; ++i)
	{
	  const auto _Nnutemp = (mu + i) * xi2 * _Nnu1 - _Nmu;
	  _Nmu = _Nnu1;
	  _Nnu1 = _Nnutemp;
	}
      _Nnu = _Nmu;
      _Npnu = nu * xi * _Nmu - _Nnu1;

      return;
    }



  template <typename Tp>
    void
    bessel_ik(const Tp nu, const Tp x,
		Tp & _Inu, Tp & _Knu, Tp & _Ipnu, Tp & _Kpnu)
    {

      //if (std::isnan(nu) || std::isnan(x))
      //  return std::numeric_limits<Tp>::quiet_NaN();

      if (x == Tp{0})
	{
	  if (nu == Tp{0})
	    {
	      _Inu = Tp{1};
	      _Knu = std::numeric_limits<Tp>::infinity();
	      _Ipnu = Tp{0};
	      _Kpnu = -std::numeric_limits<Tp>::infinity();
	    }
	  else
	    {
	      _Inu = Tp{0};
	      _Knu = std::numeric_limits<Tp>::infinity();
	      //_Ipnu = ???
	      //_Kpnu = ???
	    }
	  return;
	}

      if (x < Tp{0} || nu < Tp{0})
	throw std::runtime_error("Bad arguments in bessel_ik.");

      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp fp_min = Tp(10) * std::numeric_limits<Tp>::epsilon();
      const int max_iter = 15000;
      const Tp x_min = Tp{2};

      const int nl = static_cast<int>(nu + Tp{0.5L});

      const auto mu = nu - nl;
      const auto mu2 = mu * mu;
      const auto xi = Tp{1} / x;
      const auto xi2 = Tp{2} * xi;
      auto h = nu * xi;
      if ( h < fp_min )
	h = fp_min;
      auto b = xi2 * nu;
      auto d = Tp{0};
      auto c = h;
      int i;
      for ( i = 1; i <= max_iter; ++i )
	{
	  b += xi2;
	  d = Tp{1} / (b + d);
	  c = b + Tp{1} / c;
	  const auto del = c * d;
	  h = del * h;
	  if (std::abs(del - Tp{1}) < eps)
	    break;
	}
      if (i > max_iter)
	throw std::runtime_error( "Argument x too large in bessel_jn; "
				  "try asymptotic expansion." );
      auto _Inul = fp_min;
      auto _Ipnul = h * _Inul;
      auto _Inul1 = _Inul;
      auto _Ipnu1 = _Ipnul;
      auto fact = nu * xi;
      for (int l = nl; l >= 1; --l)
	{
	  const Tp _Inutemp = fact * _Inul + _Ipnul;
	  fact -= xi;
	  _Ipnul = fact * _Inutemp + _Inul;
	  _Inul = _Inutemp;
	}
      auto f = _Ipnul / _Inul;
      auto _Kmu, _Knu1;
      if (x < x_min)
	{
	  const auto x2 = x / Tp{2};
	  const auto pimu = Tp(PI) * mu;
	  const auto fact = (std::abs(pimu) < eps
			    ? Tp{1} : pimu / std::sin(pimu));
	  auto d = -std::log(x2);
	  auto e = mu * d;
	  const auto fact2 = (std::abs(e) < eps
			    ? Tp{1} : std::sinh(e) / e);
	  auto gam1, gam2, gampl, gammi;
	  gamma_temme(mu, gam1, gam2, gampl, gammi);
	  auto ff = fact * (gam1 * std::cosh(e) + gam2 * fact2 * d);
	  auto sum = ff;
	  e = std::exp(e);
	  auto p = e / (Tp{2} * gampl);
	  auto q = Tp{1} / (Tp{2} * e * gammi);
	  auto c = Tp{1};
	  d = x2 * x2;
	  auto sum1 = p;
	  int i;
	  for (i = 1; i <= max_iter; ++i)
	    {
	      ff = (i * ff + p + q) / (i * i - mu2);
	      c *= d / i;
	      p /= i - mu;
	      q /= i + mu;
	      const auto del = c * ff;
	      sum += del; 
	      const auto del1 = c * (p - i * ff);
	      sum1 += del1;
	      if (std::abs(del) < eps * std::abs(sum))
		break;
	    }
	  if (i > max_iter)
	    throw std::runtime_error("Bessel k series failed to converge "
				     "in bessel_jn." );
	  _Kmu = sum;
	  _Knu1 = sum1 * xi2;
	}
      else
	{
	  Tp b = Tp{2} * (Tp{1} + x);
	  Tp d = Tp{1} / b;
	  Tp delh = d;
	  Tp h = delh;
	  Tp q1 = Tp{0};
	  Tp q2 = Tp{1};
	  Tp a1 = Tp{0.25L} - mu2;
	  Tp q = c = a1;
	  Tp a = -a1;
	  Tp s = Tp{1} + q * delh;
	  int i;
	  for (i = 2; i <= max_iter; ++i)
	    {
	      a -= 2 * (i - 1);
	      c = -a * c / i;
	      const auto qnew = (q1 - b * q2) / a;
	      q1 = q2;
	      q2 = qnew;
	      q += c * qnew;
	      b += Tp{2};
	      d = Tp{1} / (b + a * d);
	      delh = (b * d - Tp{1}) * delh;
	      h += delh;
	      const auto dels = q * delh;
	      s += dels;
	      if ( std::abs(dels / s) < eps )
		break;
	    }
	  if (i > max_iter)
	    throw std::runtime_error("Steed's method failed in bessel_jn.");
	  h = a1 * h;
	  _Kmu = std::sqrt(Tp(PI) / (Tp{2} * x)) * std::exp(-x) / s;
	  _Knu1 = _Kmu * (mu + x + Tp{0.5L} - h) * xi;
	}

      Tp _Kpmu = mu * xi * _Kmu - _Knu1;
      Tp _Inumu = xi / (f * _Kmu - _Kpmu);
      _Inu = _Inumu * _Inul1 / _Inul;
      _Ipnu = _Inumu * _Ipnu1 / _Inul;
      for ( i = 1; i <= nl; ++i )
	{
	  const Tp _Knutemp = (mu + i) * xi2 * _Knu1 + _Kmu;
	  _Kmu = _Knu1;
	  _Knu1 = _Knutemp;
	}
      _Knu = _Kmu;
      _Kpnu = nu * xi * _Kmu - _Knu1;
  
      return;
    }


  ///
  ///
  ///
  template <typename Tp>
    void
    airy(const Tp x,
	   Tp & _Ai, Tp & _Bi, Tp & _Aip, Tp & _Bip)
    {
      const auto SQRT3 = std::sqrt(Tp{3});
      const auto absx = std::abs(x);
      const auto rootx = std::sqrt(absx);
      const auto z = Tp{2} * absx * rootx / Tp{3};
      if (x > Tp{0})
	{
	  Tp _Inu, _Ipnu, _Knu, _Kpnu;

	  bessel_jn(Tp{1}/Tp{3}, z, _Inu, _Knu, _Ipnu, _Kpnu);
	  _Ai = rootx * _Knu / (Tp(SQRT3) * PI);
	  _Bi = rootx * (_Knu / PI + Tp{2} * _Inu / Tp(SQRT3));

	  bessel_jn(Tp{2}/Tp{3}, z, _Inu, _Knu, _Ipnu, _Kpnu);
	  _Aip = -x * _Knu / (Tp(SQRT3) * PI);
	  _Bip = x * (_Knu / PI + Tp{2} * _Inu / Tp(SQRT3));
	}
      else if (x < Tp{0})
	{
	  Tp _Jnu, _Jpnu, _Nnu, _Npnu;

	  bessel_jn(Tp{1}/Tp{3}, z, _Jnu, _Nnu, _Jpnu, _Npnu);
	  _Ai = rootx * (_Jnu - _Nnu / Tp(SQRT3)) / Tp{2};
	  _Bi = -rootx * (_Nnu + _Jnu / Tp(SQRT3)) / Tp{2};

	  bessel_jn(Tp{2}/Tp{3}, z, _Jnu, _Nnu, _Jpnu, _Npnu);
	  _Aip = absx * (_Nnu / Tp(SQRT3) + _Jnu) / Tp{2};
	  _Bip = absx * (_Jnu / Tp(SQRT3) - _Nnu) / Tp{2};
	}
      else
	{
	  // References : Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
	  // The number is Ai(0) or 3**(-2/3)/Gamma(2/3).
	  _Ai = Tp{0.3550280538878172392600631860041831763979791741991772L};
	  _Bi = _Ai * SQRT3;

	  // References : Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
	  // The number is Ai'(0) or -3**(-1/3)/Gamma(1/3)
	  _Aip = -Tp{0.25881940379280679840518356018920396347909113835493L};
	  _Bip = -_Aip * SQRT3;
	}

      return;
    }


  ///
  ///
  ///
  template <typename Tp>
  void
    sph_bessel_jn(const int n, const Tp x,
		    Tp & jn, Tp & nn, Tp & jpn, Tp & npn)
    {

      if ( n < 0 || x < Tp{0} )
	throw std::runtime_error( "Bad arguments in sph_bessel." );

      const auto nu = Tp(n) + Tp{0.5L};

      Tp _Jnu, _Jpnu, _Nnu, _Npnu;
      bessel_jn( x, nu, _Jnu, _Nnu, _Jpnu, _Npnu );

      const auto factor = Tp(SQRTPIO2) / std::sqrt(x);

      jn = factor * _Jnu;
      nn = factor * _Nnu;
      jpn = factor * _Jpnu - jn / (Tp{2} * x);
      npn = factor * _Npnu - nn / (Tp{2} * x);

      return;
    }


  ///
  ///
  ///
  template <typename Tp>
  void
    sph_bessel_ik(const int n, const Tp x,
		    Tp & in, Tp & kn, Tp & ipn, Tp & kpn)
    {

      if ( n < 0 || x < Tp{0} )
	throw std::runtime_error( "Bad arguments in sph_bessel." );

      const auto order = Tp(n) + Tp{0.5L};

      Tp _Inu, _Ipnu, _Knu, _Kpnu;
      bessel_ik( x, order, _Inu, _Knu, _Ipnu, _Kpnu );

      const auto factor = Tp(SQRTPIO2) / std::sqrt(x);

      in = factor * _Inu;
      kn = factor * _Knu;
      ipn = factor * _Ipnu - in / (Tp{2} * x);
      kpn = factor * _Kpnu - kn / (Tp{2} * x);

      return;
    }


