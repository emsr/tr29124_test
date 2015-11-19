// Special functions -*- C++ -*-

// Copyright (C) 2006-2015 Free Software Foundation, Inc.
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

/** @file bits/sf_bessel.tcc
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
   *   @param __mu     The input parameter of the gamma functions.
   *   @param __gam1   The output function \f$ \Gamma_1(\mu) \f$
   *   @param __gam2   The output function \f$ \Gamma_2(\mu) \f$
   *   @param __gampl  The output function \f$ \Gamma(1 + \mu) \f$
   *   @param __gammi  The output function \f$ \Gamma(1 - \mu) \f$
   */
  template <typename _Tp>
    void
    __gamma_temme(const _Tp __mu,
		   _Tp & __gam1, _Tp & __gam2, _Tp & __gampl, _Tp & __gammi)
    {
      __gampl = _Tp{1} / std::::tgamma(_Tp{1} + __mu);
      __gammi = _Tp{1} / std::tgamma(_Tp{1} - __mu);

      if (std::abs(__mu) < std::numeric_limits<_Tp>::epsilon())
	__gam1 = -_Tp(__GAMMA_E);
      else
	__gam1 = (__gammi - __gampl) / (_Tp{2} * __mu);

      __gam2 = (__gammi + __gampl) / (_Tp{2});

      return;
    }


  /**
   *   @brief  Compute the Bessel @f$ J_\nu(x) @f$ and Neumann
   *           @f$ N_\nu(x) @f$ functions and their first derivatives
   *           @f$ J'_\nu(x) @f$ and @f$ N'_\nu(x) @f$ respectively.
   *           These four functions are computed together for numerical
   *           stability.
   *
   *   @param  __nu  The order of the Bessel functions.
   *   @param  __x   The argument of the Bessel functions.
   *   @param  _Jnu  The output Bessel function of the first kind.
   *   @param  _Nnu  The output Neumann function (Bessel function of the second kind).
   *   @param  _Jpnu  The output derivative of the Bessel function of the first kind.
   *   @param  _Npnu  The output derivative of the Neumann function.
   */
  template <typename _Tp>
    void
    __bessel_jn(const _Tp __nu, const _Tp __x,
		_Tp & _Jnu, _Tp & _Nnu, _Tp & _Jpnu, _Tp & _Npnu)
    {

      //if (std::isnan(__nu) || std::isnan(__x))
      //  return std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp{0})
	{
	  if (__nu == _Tp{0})
	    {
	      _Jnu = _Tp{1};
	      _Nnu = -std::numeric_limits<_Tp>::infinity();
	      _Jpnu = _Tp{0};
	      _Npnu = std::numeric_limits<_Tp>::infinity();
	    }
	  else
	    {
	      _Jnu = _Tp{0};
	      _Nnu = -std::numeric_limits<_Tp>::infinity();
	      //_Jpnu = ???
	      //_Npnu = ???
	    }
	  return;
	}

      const auto __eps = std::numeric_limits<_Tp>::epsilon();
      const auto __fp_min = _Tp(10) * std::numeric_limits<_Tp>::min();
      const int __max_iter = 15000;
      const auto __x_min = _Tp{2};

      if (__x < _Tp{0} || __nu < _Tp{0})
	throw std::runtime_error("Bad arguments in __bessel_jn.");

      const int __nl = (__x < __x_min
		    ? static_cast<int>(__nu + _Tp{0.5L})
		    : std::max(0, static_cast<int>(__nu - __x + _Tp(1.5L))));

      const auto __mu = __nu - __nl;
      const auto __mu2 = __mu * __mu;
      const auto __xi = _Tp{1} / __x;
      const auto __xi2 = _Tp{2} * __xi;
      auto __w = __xi2 / _Tp(__PI);
      int __isign = 1;
      auto __h = __nu * __xi;
      if (__h < __fp_min)
	__h = __fp_min;
      auto __b = __xi2 * __nu;
      auto __d = _Tp{0};
      auto __c = __h;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
	{
	  __b += __xi2;
	  __d = __b - __d;
	  if (std::abs(__d) < __fp_min)
	    __d = __fp_min;
	  __c = __b - _Tp{1} / __c;
	  if (std::abs(__c) < __fp_min)
	    __c = __fp_min;
	  __d = _Tp{1} / __d;
	  const auto __del = __c * __d;
	  __h = __del * __h;
	  if (__d < _Tp{0})
	    __isign = -__isign;
	  if (std::abs(__del - _Tp{1}) < __eps)
	    break;
	}
      if (__i > __max_iter)
	throw std::runtime_error( "Argument x too large in __bessel_jn; "
				  "try asymptotic expansion." );
      auto _Jnul = __isign * __fp_min;
      auto _Jpnul = __h * _Jnul;
      auto _Jnul1 = _Jnul;
      auto _Jpnu1 = _Jpnul;
      auto __fact = __nu * __xi;
      for ( int __l = __nl; __l >= 1; --__l )
	{
	  const auto _Jnutemp = __fact * _Jnul + _Jpnul;
	  __fact -= __xi;
	  _Jpnul = __fact * _Jnutemp - _Jnul;
	  _Jnul = _Jnutemp;
	}
      if (_Jnul == _Tp{0})
	_Jnul = __eps;
      auto __f = _Jpnul / _Jnul;
      _Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (__x < __x_min)
	{
	  const auto __x2 = __x / _Tp{2};
	  const auto __pimu = _Tp(__PI) * __mu;
	  auto __fact = (std::abs(__pimu) < __eps
		      ? _Tp{1} : __pimu / std::sin(__pimu));
	  auto __d = -std::log(__x2);
	  auto __e = __mu * __d;
	  auto __fact2 = (std::abs(__e) < __eps
		       ? _Tp{1} : std::sinh(__e) / __e);
	  _Tp __gam1, __gam2, __gampl, __gammi;
	  __gamma_temme(__mu, __gam1, __gam2, __gampl, __gammi);
	  auto __ff = (_Tp{2} / _Tp(__PI))
		   * __fact * (__gam1 * std::cosh(__e) + __gam2 * __fact2 * __d);
	  __e = std::exp(__e);
	  auto __p = __e / (_Tp(__PI) * __gampl);
	  auto __q = _Tp{1} / (__e * _Tp(__PI) * __gammi);
	  const auto __pimu2 = __pimu / _Tp{2};
	  auto __fact3 = (std::abs(__pimu2) < __eps
		       ? _Tp{1} : std::sin(__pimu2) / __pimu2 );
	  auto __r = _Tp(__PI) * __pimu2 * __fact3 * __fact3;
	  auto __c = _Tp{1};
	  __d = -__x2 * __x2;
	  auto __sum = __ff + __r * __q;
	  auto __sum1 = __p;
	  for (__i = 1; __i <= __max_iter; ++__i)
	    {
	      __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	      __c *= __d / _Tp(__i);
	      __p /= _Tp(__i) - __mu;
	      __q /= _Tp(__i) + __mu;
	      const auto __del = __c * (__ff + __r * __q);
	      __sum += __del; 
	      const auto __del1 = __c * __p - __i * __del;
	      __sum1 += __del1;
	      if ( std::abs(__del) < __eps * (_Tp{1} + std::abs(__sum)) )
		break;
	    }
	  if ( __i > __max_iter )
	    throw std::runtime_error("Bessel y series failed to converge "
				     "in __bessel_jn." );
	  _Nmu = -__sum;
	  _Nnu1 = -__sum1 * __xi2;
	  _Npmu = __mu * __xi * _Nmu - _Nnu1;
	  _Jmu = __w / (_Npmu - __f * _Nmu);
	}
      else
	{
	  auto __a = _Tp{0.25L} - __mu2;
	  auto __q = _Tp{1};
	  auto __p = -__xi / _Tp{2};
	  auto __br = _Tp{2} * __x;
	  auto __bi = _Tp{2};
	  auto __fact = __a * __xi / (__p * __p + __q * __q);
	  auto __cr = __br + __q * __fact;
	  auto __ci = __bi + __p * __fact;
	  auto __den = __br * __br + __bi * __bi;
	  auto __dr = __br / __den;
	  auto __di = -__bi / __den;
	  auto __dlr = __cr * __dr - __ci * __di;
	  auto __dli = __cr * __di + __ci * __dr;
	  auto __temp = __p * __dlr - __q * __dli;
	  __q = __p * __dli + __q * __dlr;
	  __p = __temp;
	  int __i;
	  for (__i = 2; __i <= __max_iter; ++__i)
	    {
	      __a += _Tp(2 * (__i - 1));
	      __bi += _Tp{2};
	      __dr = __a * __dr + __br;
	      __di = __a * __di + __bi;
	      if (std::abs(__dr) + std::abs(__di) < __fp_min)
		__dr = __fp_min;
	      __fact = __a / (__cr * __cr + __ci * __ci);
	      __cr = __br + __cr * __fact;
	      __ci = __bi - __ci * __fact;
	      if (std::abs(__cr) + std::abs(__ci) < __fp_min)
		__cr = __fp_min;
	      __den = __dr * __dr + __di * __di;
	      __dr /= __den;
	      __di /= -__den;
	      __dlr = __cr * __dr - __ci * __di;
	      __dli = __cr * __di + __ci * __dr;
	      __temp = __p * __dlr - __q * __dli;
	      __q = __p * __dli + __q * __dlr;
	      __p = __temp;
	      if (std::abs(__dlr - _Tp{1}) + std::abs(__dli) < __eps)
		break;
	  }
	  if (__i > __max_iter)
	    throw std::runtime_error("Lentz's method failed in __bessel_jn.");
	  const auto __gam = (__p - __f) / __q;
	  _Jmu = std::sqrt(__w / ((__p - __f) * __gam + __q));

	  _Jmu = ::copysign(_Jmu, _Jnul);

	  _Nmu = __gam * _Jmu;
	  _Npmu = (__p + __q / __gam) * _Nmu;
	  _Nnu1 = __mu * __xi * _Nmu - _Npmu;
      }
      __fact = _Jmu / _Jnul;
      _Jnu = __fact * _Jnul1;
      _Jpnu = __fact * _Jpnu1;
      for (__i = 1; __i <= __nl; ++__i)
	{
	  const auto _Nnutemp = (__mu + __i) * __xi2 * _Nnu1 - _Nmu;
	  _Nmu = _Nnu1;
	  _Nnu1 = _Nnutemp;
	}
      _Nnu = _Nmu;
      _Npnu = __nu * __xi * _Nmu - _Nnu1;

      return;
    }



  template <typename _Tp>
    void
    __bessel_ik(const _Tp __nu, const _Tp __x,
		_Tp & _Inu, _Tp & _Knu, _Tp & _Ipnu, _Tp & _Kpnu)
    {

      //if (std::isnan(__nu) || std::isnan(__x))
      //  return std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp{0})
	{
	  if (__nu == _Tp{0})
	    {
	      _Inu = _Tp{1};
	      _Knu = std::numeric_limits<_Tp>::infinity();
	      _Ipnu = _Tp{0};
	      _Kpnu = -std::numeric_limits<_Tp>::infinity();
	    }
	  else
	    {
	      _Inu = _Tp{0};
	      _Knu = std::numeric_limits<_Tp>::infinity();
	      //_Ipnu = ???
	      //_Kpnu = ???
	    }
	  return;
	}

      if (__x < _Tp{0} || __nu < _Tp{0})
	throw std::runtime_error("Bad arguments in __bessel_ik.");

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __fp_min = _Tp(10) * std::numeric_limits<_Tp>::epsilon();
      const int __max_iter = 15000;
      const _Tp __x_min = _Tp{2};

      const int __nl = static_cast<int>(__nu + _Tp{0.5L});

      const auto __mu = __nu - __nl;
      const auto __mu2 = __mu * __mu;
      const auto __xi = _Tp{1} / __x;
      const auto __xi2 = _Tp{2} * __xi;
      auto __h = __nu * __xi;
      if ( __h < __fp_min )
	__h = __fp_min;
      auto __b = __xi2 * __nu;
      auto __d = _Tp{0};
      auto __c = __h;
      int __i;
      for ( __i = 1; __i <= __max_iter; ++__i )
	{
	  __b += __xi2;
	  __d = _Tp{1} / (__b + __d);
	  __c = __b + _Tp{1} / __c;
	  const auto __del = __c * __d;
	  __h = __del * __h;
	  if (std::abs(__del - _Tp{1}) < __eps)
	    break;
	}
      if (__i > __max_iter)
	throw std::runtime_error( "Argument x too large in __bessel_jn; "
				  "try asymptotic expansion." );
      auto _Inul = __fp_min;
      auto _Ipnul = __h * _Inul;
      auto _Inul1 = _Inul;
      auto _Ipnu1 = _Ipnul;
      auto __fact = __nu * __xi;
      for (int __l = __nl; __l >= 1; --__l)
	{
	  const _Tp _Inutemp = __fact * _Inul + _Ipnul;
	  __fact -= __xi;
	  _Ipnul = __fact * _Inutemp + _Inul;
	  _Inul = _Inutemp;
	}
      auto __f = _Ipnul / _Inul;
      auto _Kmu, _Knu1;
      if (__x < __x_min)
	{
	  const auto __x2 = __x / _Tp{2};
	  const auto __pimu = _Tp(__PI) * __mu;
	  const auto __fact = (std::abs(__pimu) < __eps
			    ? _Tp{1} : __pimu / std::sin(__pimu));
	  auto __d = -std::log(__x2);
	  auto __e = __mu * __d;
	  const auto __fact2 = (std::abs(__e) < __eps
			    ? _Tp{1} : std::sinh(__e) / __e);
	  auto __gam1, __gam2, __gampl, __gammi;
	  __gamma_temme(__mu, __gam1, __gam2, __gampl, __gammi);
	  auto __ff = __fact * (__gam1 * std::cosh(__e) + __gam2 * __fact2 * __d);
	  auto __sum = __ff;
	  __e = std::exp(__e);
	  auto __p = __e / (_Tp{2} * __gampl);
	  auto __q = _Tp{1} / (_Tp{2} * __e * __gammi);
	  auto __c = _Tp{1};
	  __d = __x2 * __x2;
	  auto __sum1 = __p;
	  int __i;
	  for (__i = 1; __i <= __max_iter; ++__i)
	    {
	      __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	      __c *= __d / __i;
	      __p /= __i - __mu;
	      __q /= __i + __mu;
	      const auto __del = __c * __ff;
	      __sum += __del; 
	      const auto __del1 = __c * (__p - __i * __ff);
	      __sum1 += __del1;
	      if (std::abs(__del) < __eps * std::abs(__sum))
		break;
	    }
	  if (__i > __max_iter)
	    throw std::runtime_error("Bessel k series failed to converge "
				     "in __bessel_jn." );
	  _Kmu = __sum;
	  _Knu1 = __sum1 * __xi2;
	}
      else
	{
	  _Tp __b = _Tp{2} * (_Tp{1} + __x);
	  _Tp __d = _Tp{1} / __b;
	  _Tp __delh = __d;
	  _Tp __h = __delh;
	  _Tp __q1 = _Tp{0};
	  _Tp __q2 = _Tp{1};
	  _Tp __a1 = _Tp{0.25L} - __mu2;
	  _Tp __q = __c = __a1;
	  _Tp __a = -__a1;
	  _Tp __s = _Tp{1} + __q * __delh;
	  int __i;
	  for (__i = 2; __i <= __max_iter; ++__i)
	    {
	      __a -= 2 * (__i - 1);
	      __c = -__a * __c / __i;
	      const auto __qnew = (__q1 - __b * __q2) / __a;
	      __q1 = __q2;
	      __q2 = __qnew;
	      __q += __c * __qnew;
	      __b += _Tp{2};
	      __d = _Tp{1} / (__b + __a * __d);
	      __delh = (__b * __d - _Tp{1}) * __delh;
	      __h += __delh;
	      const auto __dels = __q * __delh;
	      __s += __dels;
	      if ( std::abs(__dels / __s) < __eps )
		break;
	    }
	  if (__i > __max_iter)
	    throw std::runtime_error("Steed's method failed in __bessel_jn.");
	  __h = __a1 * __h;
	  _Kmu = std::sqrt(_Tp(__PI) / (_Tp{2} * __x)) * std::exp(-__x) / __s;
	  _Knu1 = _Kmu * (__mu + __x + _Tp{0.5L} - __h) * __xi;
	}

      _Tp _Kpmu = __mu * __xi * _Kmu - _Knu1;
      _Tp _Inumu = __xi / (__f * _Kmu - _Kpmu);
      _Inu = _Inumu * _Inul1 / _Inul;
      _Ipnu = _Inumu * _Ipnu1 / _Inul;
      for ( __i = 1; __i <= __nl; ++__i )
	{
	  const _Tp _Knutemp = (__mu + __i) * __xi2 * _Knu1 + _Kmu;
	  _Kmu = _Knu1;
	  _Knu1 = _Knutemp;
	}
      _Knu = _Kmu;
      _Kpnu = __nu * __xi * _Kmu - _Knu1;
  
      return;
    }


  ///
  ///
  ///
  template <typename _Tp>
    void
    __airy(const _Tp __x,
	   _Tp & _Ai, _Tp & _Bi, _Tp & _Aip, _Tp & _Bip)
    {
      const auto __SQRT3 = std::sqrt(_Tp{3});
      const auto __absx = std::abs(__x);
      const auto __rootx = std::sqrt(__absx);
      const auto __z = _Tp{2} * __absx * __rootx / _Tp{3};
      if (__x > _Tp{0})
	{
	  _Tp _Inu, _Ipnu, _Knu, _Kpnu;

	  __bessel_jn(_Tp{1}/_Tp{3}, __z, _Inu, _Knu, _Ipnu, _Kpnu);
	  _Ai = __rootx * _Knu / (_Tp(__SQRT3) * __PI);
	  _Bi = __rootx * (_Knu / __PI + _Tp{2} * _Inu / _Tp(__SQRT3));

	  __bessel_jn(_Tp{2}/_Tp{3}, __z, _Inu, _Knu, _Ipnu, _Kpnu);
	  _Aip = -__x * _Knu / (_Tp(__SQRT3) * __PI);
	  _Bip = __x * (_Knu / __PI + _Tp{2} * _Inu / _Tp(__SQRT3));
	}
      else if (__x < _Tp{0})
	{
	  _Tp _Jnu, _Jpnu, _Nnu, _Npnu;

	  __bessel_jn(_Tp{1}/_Tp{3}, __z, _Jnu, _Nnu, _Jpnu, _Npnu);
	  _Ai = __rootx * (_Jnu - _Nnu / _Tp(__SQRT3)) / _Tp{2};
	  _Bi = -__rootx * (_Nnu + _Jnu / _Tp(__SQRT3)) / _Tp{2};

	  __bessel_jn(_Tp{2}/_Tp{3}, __z, _Jnu, _Nnu, _Jpnu, _Npnu);
	  _Aip = __absx * (_Nnu / _Tp(__SQRT3) + _Jnu) / _Tp{2};
	  _Bip = __absx * (_Jnu / _Tp(__SQRT3) - _Nnu) / _Tp{2};
	}
      else
	{
	  // References : Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
	  // The number is Ai(0) or 3**(-2/3)/Gamma(2/3).
	  _Ai = _Tp{0.3550280538878172392600631860041831763979791741991772L};
	  _Bi = _Ai * __SQRT3;

	  // References : Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
	  // The number is Ai'(0) or -3**(-1/3)/Gamma(1/3)
	  _Aip = -_Tp{0.25881940379280679840518356018920396347909113835493L};
	  _Bip = -_Aip * __SQRT3;
	}

      return;
    }


  ///
  ///
  ///
  template <typename _Tp>
  void
    __sph_bessel_jn(const int __n, const _Tp __x,
		    _Tp & __jn, _Tp & __nn, _Tp & __jpn, _Tp & __npn)
    {

      if ( __n < 0 || __x < _Tp{0} )
	throw std::runtime_error( "Bad arguments in sph_bessel." );

      const auto __nu = _Tp(__n) + _Tp{0.5L};

      _Tp _Jnu, _Jpnu, _Nnu, _Npnu;
      __bessel_jn( __x, __nu, _Jnu, _Nnu, _Jpnu, _Npnu );

      const auto __factor = _Tp(__SQRTPIO2) / std::sqrt(__x);

      __jn = __factor * _Jnu;
      __nn = __factor * _Nnu;
      __jpn = __factor * _Jpnu - __jn / (_Tp{2} * __x);
      __npn = __factor * _Npnu - __nn / (_Tp{2} * __x);

      return;
    }


  ///
  ///
  ///
  template <typename _Tp>
  void
    __sph_bessel_ik(const int __n, const _Tp __x,
		    _Tp & __in, _Tp & __kn, _Tp & __ipn, _Tp & __kpn)
    {

      if ( __n < 0 || __x < _Tp{0} )
	throw std::runtime_error( "Bad arguments in sph_bessel." );

      const auto __order = _Tp(__n) + _Tp{0.5L};

      _Tp _Inu, _Ipnu, _Knu, _Kpnu;
      __bessel_ik( __x, __order, _Inu, _Knu, _Ipnu, _Kpnu );

      const auto __factor = _Tp(__SQRTPIO2) / std::sqrt(__x);

      __in = __factor * _Inu;
      __kn = __factor * _Knu;
      __ipn = __factor * _Ipnu - __in / (_Tp{2} * __x);
      __kpn = __factor * _Kpnu - __kn / (_Tp{2} * __x);

      return;
    }


