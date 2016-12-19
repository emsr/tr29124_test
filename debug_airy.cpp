/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o debug_airy debug_airy.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
./debug_airy

$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -I. -o debug_airy debug_airy.cpp -lquadmath -L. -lwgsl
./debug_airy
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include <bits/specfun.h>
#include "wrap_gsl.h"

double
chebyshev_eval(double a, double b, double * c, int m, double x)
{
  if ((x - a) * (x - b) > 0.0)
  {
    std::cerr << "x not in range in routine chebyshev_eval.";
    return 0.0;
  }

  double d = 0.0, dd = 0.0, y;
  double y2 = 2.0 * (y = (2.0 * x - a - b) / ( b - a ));
  for (int j = m - 1; j >= 1; --j)
  {
    double sv = d;
    d = y2 * d - dd + c[j];
    dd = sv;
  }
  return y * d - dd + 0.5 * c[0];
}

void
bessel_cheb(double x, double& gam1, double& gam2, double& gampl, double& gammi)
{
  double xx;
  static double c1[]
  {
    -1.142022680371168,
     0.0065165112670737,
     0.0003087090173086,
    -0.0000034706269649,
     0.0000000069437664,
     0.0000000000367765,
    -0.0000000000001356
  };
  static double c2[]
  {
    1.843740587300905,
   -0.0786528408447867,
    0.0012719271366546,
   -0.0000049717367042,
   -0.0000000331261198,
    0.0000000002423096,
   -0.0000000000001702,
   -0.00000000000000149
  };

  xx = 8.0 * x * x - 1.0;
  gam1 = chebyshev_eval( -1.0, +1.0, c1, 7, xx );
  gam2 = chebyshev_eval( -1.0, +1.0, c2, 8, xx );
  gampl = gam2 - x * gam1;
  gammi = gam2 + x * gam1;
}

  template<typename _Tp>
    void
    __gamma_temme(_Tp __mu,
		  _Tp & __gam1, _Tp & __gam2, _Tp & __gampl, _Tp & __gammi)
    {
      constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr _Tp _S_gamma_E = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      __gampl = _Tp{1} / std::tgamma(_Tp{1} + __mu);
      __gammi = _Tp{1} / std::tgamma(_Tp{1} - __mu);

      if (std::abs(__mu) < _S_eps)
	__gam1 = -_S_gamma_E;
      else
	__gam1 = (__gammi - __gampl) / (_Tp{2} * __mu);

      __gam2 = (__gammi + __gampl) / _Tp{2};

      return;
    }

void
old_bessel_chunk(double xnu, double x, double& rj, double& ry, double& rjp, double& ryp)
{
    const double  EPS    =  1.0e-16;
    const double  FPMIN  =  1.0e-30;
    const int     MAXIT  =  10000;
    const double  XMIN   =  2.0;
    const double  PI     =  3.14159265358979;

    int i, isign, l, n;
    double b, c, d, del, del1, e, f, fact, fact2,
           fact3, ff, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
           rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
           w, x2, xi, xi2, xmu, xmu2;

    if ( x <= 0.0 || xnu < 0.0 ) std::cout << "Bad arguments in bessel_jy.\n";
std::cout << '\n';
    n = ( x < XMIN ? (int)(xnu + 0.5) : std::max( 0, (int)(xnu - x + 1.5) ) );
    xmu = xnu - n;
    xmu2 = xmu*xmu;
    xi = 1.0/x;
    xi2 = 2.0*xi;
    w = xi2/PI;
std::cout << "n   =" << n << '\n';
std::cout << "xmu  =" << xmu << '\n';
std::cout << "xmu2 =" << xmu2 << '\n';
std::cout << "xi   =" << xi << '\n';
std::cout << "xi2  =" << xi2 << '\n';
std::cout << "w    =" << w << '\n';
    isign = 1;
    h = xnu*xi;
    if ( h < FPMIN ) h = FPMIN;
    b = xi2*xnu;
    d = 0.0;
    c = h;
std::cout << "h   =" << h << '\n';
std::cout << "b   =" << b << '\n';
std::cout << "d   =" << d << '\n';
std::cout << "c   =" << c << '\n';
    for ( i = 1; i <= MAXIT; ++i ) {
        b += xi2;
        d = b - d;
        if ( abs(d) < FPMIN ) d = FPMIN;
        c = b - 1.0/c;
        if ( abs(c) < FPMIN ) c = FPMIN;
        d = 1.0/d;
        del = c*d;
        h = del*h;
        if ( d < 0.0 ) isign = -isign;
        if ( abs(del - 1.0) < EPS ) break;
    }
    if ( i > MAXIT ) std::cout << "Argument x too large in bessel_jy; try asymptotic expansion.\n";
std::cout << "h   =" << h << '\n';
std::cout << "b   =" << b << '\n';
std::cout << "d   =" << d << '\n';
std::cout << "c   =" << c << '\n';
    rjl = isign*FPMIN;
    rjpl = h*rjl;
    rjl1 = rjl;
    rjp1 = rjpl;
    fact = xnu*xi;
    for ( l = n; l >= 1; --l ) {
        rjtemp = fact*rjl + rjpl;
        fact -= xi;
        rjpl = fact*rjtemp - rjl;
        rjl = rjtemp;
    }
    if ( rjl == 0.0 ) rjl = EPS;
std::cout << "Jnul   =" << rjl << '\n';
std::cout << "Jpnul  =" << rjpl << '\n';
std::cout << "Jnul1  =" << rjl1 << '\n';
std::cout << "Jpnu1  =" << rjp1 << '\n';
    f= rjpl/rjl;
std::cout << "f   =" << f << '\n';
    if ( x < XMIN ) {
        x2 = 0.5*x;
        pimu = PI*xmu;
        fact = ( abs(pimu) < EPS ? 1.0 : pimu/sin(pimu) );
        d = -log(x2);
        e = xmu*d;
        fact2 = ( abs(e) < EPS ? 1.0 : sinh(e)/e );
        bessel_cheb( xmu, gam1, gam2, gampl, gammi );
        ff = (2.0/PI)*fact*(gam1*cosh(e) + gam2*fact2*d);
        e = exp(e);
        p = e/(PI*gampl);
        q = 1.0/(e*PI*gammi);
        pimu2 = 0.5*pimu;
        fact3 = (abs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2 );
        r = PI*pimu2*fact3*fact3;
        c = 1.0;
        d = -x2*x2;
        sum = ff + r*q;
        sum1 = p;
std::cout << "x2  =" << x2 << '\n';
std::cout << "pimu  =" << pimu << '\n';
std::cout << "fact  =" << fact << '\n';
std::cout << "fact2 =" << fact2 << '\n';
std::cout << "fact3 =" << fact3 << '\n';
std::cout << "ff    =" << ff << '\n';
std::cout << "gam1  =" << gam1 << '\n';
std::cout << "gam2  =" << gam2 << '\n';
std::cout << "gampl =" << gampl << '\n';
std::cout << "gammi =" << gammi << '\n';
std::cout << "d     =" << d << '\n';
std::cout << "e     =" << e << '\n';
std::cout << "p     =" << p << '\n';
std::cout << "q     =" << q << '\n';
std::cout << "r     =" << r << '\n';
        for ( i = 1; i <= MAXIT; ++i ) {
            ff = (i*ff + p + q)/(i*i - xmu2);
            c *= d/i;
            p /= i - xmu;
            q /= i + xmu;
            del = c*(ff + r*q);
            sum += del; 
            del1 = c*p - i*del;
            sum1 += del1;
std::cout << "  i =" << i <<  "  del = " << del <<  "  del1 = " << del1 << '\n';
            if ( abs(del) < EPS*(1.0 + abs(sum)) ) break;
        }
        if ( i > MAXIT ) std::cout << "Bessel y series failed to converge in bessel_jy.\n";
std::cout << "i     =" << i << '\n';
std::cout << "ff    =" << ff << '\n';
std::cout << "c     =" << c << '\n';
std::cout << "p     =" << p << '\n';
std::cout << "q     =" << q << '\n';
std::cout << "sum   =" << sum << '\n';
std::cout << "sum1  =" << sum1 << '\n';
        rymu = -sum;
        ry1 = -sum1*xi2;
        rymup = xmu*xi*rymu - ry1;
        rjmu = w/(rymup - f*rymu);
std::cout << "Nmu   =" << rymu << '\n';
std::cout << "Nnu1  =" << ry1 << '\n';
std::cout << "Npmu  =" << rymup << '\n';
std::cout << "Jmu   =" << rjmu << '\n';
    fact = rjmu/rjl;
    rj = rjl1*fact;
    rjp = rjp1*fact;
std::cout << "fact  =" << fact << '\n';
std::cout << "Jnu   =" << rj << '\n';
std::cout << "Jpnu  =" << rjp << '\n';
std::cout << "n     =" << n << '\n';
    for ( i = 1; i <= n; ++i ) {
        rytemp = (xmu + i)*xi2*ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
std::cout << "  i =" << i << "  Nmu=" << rymu << "  Nnu1=" << ry1 << '\n';
    }
    ry = rymu;
    ryp = xnu*xi*rymu - ry1;
std::cout << "Nnu   =" << ry << '\n';
std::cout << "Npnu  =" << ryp << '\n';
    }
}

template<typename _Tp>
  void
  new_bessel_chunk(_Tp __nu, _Tp __x, _Tp& _Jnu, _Tp& _Jpnu, _Tp& _Nnu, _Tp& _Npnu)
  {
      constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr _Tp _S_inf = std::numeric_limits<_Tp>::infinity();
      constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();
      if (__x == _Tp{0})
	{
	  if (__nu == _Tp{0})
	    {
	      _Jnu = _Tp{1};
	      _Jpnu = _Tp{0};
	    }
	  else if (__nu == _Tp{1})
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0.5L};
	    }
	  else
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0};
	    }
	  _Nnu = -_S_inf;
	  _Npnu = _S_inf;
	  return;
	}

      constexpr _Tp _S_fp_min = std::sqrt(std::numeric_limits<_Tp>::min());
      constexpr int _S_max_iter = 15000;
      constexpr _Tp _S_x_min = _Tp{2};

std::cout << '\n';
      const int __n = (__x < _S_x_min
		    ? std::nearbyint(__nu)
		    : std::max(0,
			       static_cast<int>(__nu - __x + _Tp{1.5L})));

      const _Tp __mu = __nu - __n;
      const _Tp __mu2 = __mu * __mu;
      const _Tp __xi = _Tp{1} / __x;
      const _Tp __xi2 = _Tp{2} * __xi;
      const _Tp _Wronski = __xi2 / _S_pi;
std::cout << "n    =" << __n << '\n';
std::cout << "xmu  =" << __mu << '\n';
std::cout << "xmu2 =" << __mu2 << '\n';
std::cout << "xi   =" << __xi << '\n';
std::cout << "xi2  =" << __xi2 << '\n';
std::cout << "w    =" << _Wronski << '\n';
      int __isign = 1;
      _Tp __h = __nu * __xi;
      if (__h < _S_fp_min)
	__h = _S_fp_min;
      _Tp __b = __xi2 * __nu;
      _Tp __d = _Tp{0};
      _Tp __c = __h;
std::cout << "h   =" << __h << '\n';
std::cout << "b   =" << __b << '\n';
std::cout << "d   =" << __d << '\n';
std::cout << "c   =" << __c << '\n';
      int __i;
      for (__i = 1; __i <= _S_max_iter; ++__i)
	{
	  __b += __xi2;
	  __d = __b - __d;
	  if (std::abs(__d) < _S_fp_min)
	    __d = _S_fp_min;
	  __c = __b - _Tp{1} / __c;
	  if (std::abs(__c) < _S_fp_min)
	    __c = _S_fp_min;
	  __d = _Tp{1} / __d;
	  const _Tp __del = __c * __d;
	  __h *= __del;
	  if (__d < _Tp{0})
	    __isign = -__isign;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	}
      if (__i > _S_max_iter)
	std::cout << "__bessel_jn: argument x too large;"
				       " try asymptotic expansion\n";
std::cout << "h   =" << __h << '\n';
std::cout << "b   =" << __b << '\n';
std::cout << "d   =" << __d << '\n';
std::cout << "c   =" << __c << '\n';
      _Tp _Jnul = __isign * _S_fp_min;
      _Tp _Jpnul = __h * _Jnul;
      _Tp _Jnul1 = _Jnul;
      _Tp _Jpnu1 = _Jpnul;
      _Tp __fact = __nu * __xi;
      for (int __l = __n; __l >= 1; --__l)
	{
	  const _Tp _Jnutemp = __fact * _Jnul + _Jpnul;
	  __fact -= __xi;
	  _Jpnul = __fact * _Jnutemp - _Jnul;
	  _Jnul = _Jnutemp;
	}
      if (_Jnul == _Tp{0})
	_Jnul = _S_eps;
std::cout << "Jnul   =" << _Jnul << '\n';
std::cout << "Jpnul  =" << _Jpnul << '\n';
std::cout << "Jnul1  =" << _Jnul1 << '\n';
std::cout << "Jpnu1  =" << _Jpnu1 << '\n';

      _Tp __f = _Jpnul / _Jnul;
std::cout << "f   =" << __f << '\n';
      _Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (__x < _S_x_min)
	{
	  const _Tp __x2 = __x / _Tp{2};
	  const _Tp __pimu = _S_pi * __mu;
	  const _Tp __fact = (std::abs(__pimu) < _S_eps
		      ? _Tp{1}
		      : __pimu / std::sin(__pimu));
	  _Tp __d = -std::log(__x2);
	  _Tp __e = __mu * __d;
	  const _Tp __fact2 = (std::abs(__e) < _S_eps
			    ? _Tp{1}
			    : std::sinh(__e) / __e);
	  _Tp __gam1, __gam2, __gampl, __gammi;
	  __gamma_temme(__mu, __gam1, __gam2, __gampl, __gammi);
	  _Tp __ff = (_Tp{2} / _S_pi) * __fact
		   * (__gam1 * std::cosh(__e) + __gam2 * __fact2 * __d);
	  __e = std::exp(__e);
	  _Tp __p = __e / (_S_pi * __gampl);
	  _Tp __q = _Tp{1} / (__e * _S_pi * __gammi);
	  const _Tp __pimu2 = __pimu / _Tp{2};
	  _Tp __fact3 = (std::abs(__pimu2) < _S_eps
		       ? _Tp{1} : std::sin(__pimu2) / __pimu2 );
	  _Tp __r = _S_pi * __pimu2 * __fact3 * __fact3;
	  _Tp __c = _Tp{1};
	  __d = -__x2 * __x2;
std::cout << "x2  =" << __x2 << '\n';
std::cout << "pimu  =" << __pimu << '\n';
std::cout << "fact  =" << __fact << '\n';
std::cout << "fact2 =" << __fact2 << '\n';
std::cout << "fact3 =" << __fact3 << '\n';
std::cout << "ff    =" << __ff << '\n';
std::cout << "gam1  =" << __gam1 << '\n';
std::cout << "gam2  =" << __gam2 << '\n';
std::cout << "gampl =" << __gampl << '\n';
std::cout << "gammi =" << __gammi << '\n';
std::cout << "d     =" << __d << '\n';
std::cout << "e     =" << __e << '\n';
std::cout << "p     =" << __p << '\n';
std::cout << "q     =" << __q << '\n';
std::cout << "r     =" << __r << '\n';
	  _Tp __sum = __ff + __r * __q;
	  _Tp __sum1 = __p;
	  int __i;
	  for (__i = 1; __i <= _S_max_iter; ++__i)
	    {
	      __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	      __c *= __d / _Tp(__i);
	      __p /= _Tp(__i) - __mu;
	      __q /= _Tp(__i) + __mu;
	      const _Tp __del = __c * (__ff + __r * __q);
	      __sum += __del;
	      const _Tp __del1 = __c * __p - _Tp(__i) * __del;
	      __sum1 += __del1;
std::cout << "  i =" << __i <<  "  del = " << __del <<  "  del1 = " << __del1 << '\n';
	      if (std::abs(__del) < _S_eps * (_Tp{1} + std::abs(__sum)))
		break;
	    }
	  if (__i > _S_max_iter)
	    std::cout << "__bessel_jn: Y-series failed to converge\n";
std::cout << "i     =" << __i << '\n';
std::cout << "ff    =" << __ff << '\n';
std::cout << "c     =" << __c << '\n';
std::cout << "p     =" << __p << '\n';
std::cout << "q     =" << __q << '\n';
std::cout << "sum   =" << __sum << '\n';
std::cout << "sum1  =" << __sum1 << '\n';
	  _Nmu = -__sum;
	  _Nnu1 = -__sum1 * __xi2;
	  _Npmu = __mu * __xi * _Nmu - _Nnu1;
	  _Jmu = _Wronski / (_Npmu - __f * _Nmu);
std::cout << "Nmu   =" << _Nmu << '\n';
std::cout << "Nnu1  =" << _Nnu1 << '\n';
std::cout << "Npmu  =" << _Npmu << '\n';
std::cout << "Jmu   =" << _Jmu << '\n';
	}
      __fact = _Jmu / _Jnul;
      _Jnu = __fact * _Jnul1;
      _Jpnu = __fact * _Jpnu1;
std::cout << "fact  =" << __fact << '\n';
std::cout << "Jnu   =" << _Jnu << '\n';
std::cout << "Jpnu  =" << _Jpnu << '\n';
std::cout << "n     =" << __n << '\n';
      for (int __i = 1; __i <= __n; ++__i)
	{
	  const _Tp _Nnutemp = (__mu + __i) * __xi2 * _Nnu1 - _Nmu;
	  _Nmu = _Nnu1;
	  _Nnu1 = _Nnutemp;
std::cout << "  i =" << __i << "  Nmu=" << _Nmu << "  Nnu1=" << _Nnu1 << '\n';
	}
      _Nnu = _Nmu;
      _Npnu = __nu * __xi * _Nmu - _Nnu1;
std::cout << "Nnu   =" << _Nnu << '\n';
std::cout << "Npnu  =" << _Npnu << '\n';
  }

int
main()
{
  std::cout.precision(16);
  std::cout.flags(std::ios::showpoint);

  for (auto i = -200; i <= +200; ++i)
    {
      auto x = i * 0.05;
      auto Ai = __gnu_cxx::airy_ai(x);
      auto Bi = __gnu_cxx::airy_bi(x);
      auto Ai_gsl = gsl::airy_ai(x);
      auto Bi_gsl = gsl::airy_bi(x);
      //double Ai_tr1, Bi_tr1, Aip_tr1, Bip_tr1;
      //std::tr1::__detail::__airy(x, Ai_tr1, Bi_tr1, Aip_tr1, Bip_tr1);
      std::cout << ' ' << std::setw(8) << x
                << ' ' << std::setw(20) << Ai
                << ' ' << std::setw(20) << Bi
                << ' ' << std::setw(20) << Ai_gsl
                << ' ' << std::setw(20) << Bi_gsl
                << ' ' << std::setw(20) << std::abs(Ai - Ai_gsl)
                << ' ' << std::setw(20) << std::abs(Bi - Bi_gsl) << '\n';
    }

  for (auto i = -10; i <= +10; ++i)
    {
      auto x = i * 0.05;
      double gam1, gam2, gampl, gammi;
      bessel_cheb(x, gam1, gam2, gampl, gammi);
      std::cout << '\n';
      std::cout << ' ' << std::setw(8) << x
                << ' ' << std::setw(20) << gam1
                << ' ' << std::setw(20) << gam2
                << ' ' << std::setw(20) << gampl
                << ' ' << std::setw(20) << gammi << '\n';
      __gamma_temme(x, gam1, gam2, gampl, gammi);
      std::cout << ' ' << std::setw(8) << x
                << ' ' << std::setw(20) << gam1
                << ' ' << std::setw(20) << gam2
                << ' ' << std::setw(20) << gampl
                << ' ' << std::setw(20) << gammi << '\n';
    }

  double nu = 0.3333333333333333;
  double x = 1.0;
  double J, Y, Jp, Yp;
  old_bessel_chunk(nu, x, J, Y, Jp, Yp);
  std::cout << "OLD: Jnu=" << J << '\n';
  std::cout << "OLD: Nnu=" << Y << '\n';
  new_bessel_chunk(nu, x, J, Y, Jp, Yp);
  std::cout << "NEW: Jnu=" << J << '\n';
  std::cout << "NEW: Nnu=" << Y << '\n';
  std::cout << "GSL: Jnu=" << gsl::cyl_bessel_j(nu, x) << '\n';
  std::cout << "GSL: Nnu=" << gsl::cyl_neumann(nu, x) << '\n';
}
