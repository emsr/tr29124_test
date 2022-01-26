/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/special_functions.h>

#include <wrap_gsl.h>

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
    gamma_temme(_Tp mu,
		  _Tp & gam1, _Tp & gam2, _Tp & gampl, _Tp & gammi)
    {
      constexpr _Tp s_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr _Tp s_gamma_E = emsr::egamma_v<_Tp>;
      gampl = _Tp{1} / std::tgamma(_Tp{1} + mu);
      gammi = _Tp{1} / std::tgamma(_Tp{1} - mu);

      if (std::abs(mu) < s_eps)
	gam1 = -s_gamma_E;
      else
	gam1 = (gammi - gampl) / (_Tp{2} * mu);

      gam2 = (gammi + gampl) / _Tp{2};

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
    double b, c, d, del, e, f, fact, fact2,
           ff, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
           rjl1, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
           w, xi, xi2, xmu, xmu2;

    if ( x <= 0.0 || xnu < 0.0 ) std::cout << "Bad arguments in bessel_jy.\n";
std::cout << '\n';
    n = ( x < XMIN ? int(xnu + 0.5) : std::max( 0, int(xnu - x + 1.5) ) );
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
        double x2 = 0.5*x;
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
        double fact3 = (abs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2 );
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
            double del1 = c*p - i*del;
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
        double rjmu = w/(rymup - f*rymu);
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
  new_bessel_chunk(_Tp nu, _Tp x, _Tp& _Jnu, _Tp& _Jpnu, _Tp& _Nnu, _Tp& _Npnu)
  {
      constexpr _Tp s_pi = emsr::pi_v<_Tp>;
      constexpr _Tp s_inf = std::numeric_limits<_Tp>::infinity();
      constexpr _Tp s_eps = std::numeric_limits<_Tp>::epsilon();
      if (x == _Tp{0})
	{
	  if (nu == _Tp{0})
	    {
	      _Jnu = _Tp{1};
	      _Jpnu = _Tp{0};
	    }
	  else if (nu == _Tp{1})
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0.5L};
	    }
	  else
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0};
	    }
	  _Nnu = -s_inf;
	  _Npnu = s_inf;
	  return;
	}

      constexpr _Tp s_fp_min = std::sqrt(std::numeric_limits<_Tp>::min());
      constexpr int s_max_iter = 15000;
      constexpr _Tp s_x_min = _Tp{2};

std::cout << '\n';
      const int n = (x < s_x_min
		    ? std::nearbyint(nu)
		    : std::max(0,
			       static_cast<int>(nu - x + _Tp{1.5L})));

      const _Tp mu = nu - n;
      const _Tp mu2 = mu * mu;
      const _Tp xi = _Tp{1} / x;
      const _Tp xi2 = _Tp{2} * xi;
      const _Tp _Wronski = xi2 / s_pi;
std::cout << "n    =" << n << '\n';
std::cout << "xmu  =" << mu << '\n';
std::cout << "xmu2 =" << mu2 << '\n';
std::cout << "xi   =" << xi << '\n';
std::cout << "xi2  =" << xi2 << '\n';
std::cout << "w    =" << _Wronski << '\n';
      int isign = 1;
      _Tp h = nu * xi;
      if (h < s_fp_min)
	h = s_fp_min;
      _Tp b = xi2 * nu;
      _Tp d = _Tp{0};
      _Tp c = h;
std::cout << "h   =" << h << '\n';
std::cout << "b   =" << b << '\n';
std::cout << "d   =" << d << '\n';
std::cout << "c   =" << c << '\n';
      int i;
      for (i = 1; i <= s_max_iter; ++i)
	{
	  b += xi2;
	  d = b - d;
	  if (std::abs(d) < s_fp_min)
	    d = s_fp_min;
	  c = b - _Tp{1} / c;
	  if (std::abs(c) < s_fp_min)
	    c = s_fp_min;
	  d = _Tp{1} / d;
	  const _Tp del = c * d;
	  h *= del;
	  if (d < _Tp{0})
	    isign = -isign;
	  if (std::abs(del - _Tp{1}) < s_eps)
	    break;
	}
      if (i > s_max_iter)
	std::cout << "bessel_jn: argument x too large;"
				       " try asymptotic expansion\n";
std::cout << "h   =" << h << '\n';
std::cout << "b   =" << b << '\n';
std::cout << "d   =" << d << '\n';
std::cout << "c   =" << c << '\n';
      _Tp _Jnul = isign * s_fp_min;
      _Tp _Jpnul = h * _Jnul;
      _Tp _Jnul1 = _Jnul;
      _Tp _Jpnu1 = _Jpnul;
      _Tp fact = nu * xi;
      for (int l = n; l >= 1; --l)
	{
	  const _Tp _Jnutemp = fact * _Jnul + _Jpnul;
	  fact -= xi;
	  _Jpnul = fact * _Jnutemp - _Jnul;
	  _Jnul = _Jnutemp;
	}
      if (_Jnul == _Tp{0})
	_Jnul = s_eps;
std::cout << "Jnul   =" << _Jnul << '\n';
std::cout << "Jpnul  =" << _Jpnul << '\n';
std::cout << "Jnul1  =" << _Jnul1 << '\n';
std::cout << "Jpnu1  =" << _Jpnu1 << '\n';

      _Tp f = _Jpnul / _Jnul;
std::cout << "f   =" << f << '\n';
      _Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (x < s_x_min)
	{
	  const _Tp x2 = x / _Tp{2};
	  const _Tp pimu = s_pi * mu;
	  const _Tp fact = (std::abs(pimu) < s_eps
		      ? _Tp{1}
		      : pimu / std::sin(pimu));
	  _Tp d = -std::log(x2);
	  _Tp e = mu * d;
	  const _Tp fact2 = (std::abs(e) < s_eps
			    ? _Tp{1}
			    : std::sinh(e) / e);
	  _Tp gam1, gam2, gampl, gammi;
	  gamma_temme(mu, gam1, gam2, gampl, gammi);
	  _Tp ff = (_Tp{2} / s_pi) * fact
		   * (gam1 * std::cosh(e) + gam2 * fact2 * d);
	  e = std::exp(e);
	  _Tp p = e / (s_pi * gampl);
	  _Tp q = _Tp{1} / (e * s_pi * gammi);
	  const _Tp pimu2 = pimu / _Tp{2};
	  _Tp fact3 = (std::abs(pimu2) < s_eps
		       ? _Tp{1} : std::sin(pimu2) / pimu2 );
	  _Tp r = s_pi * pimu2 * fact3 * fact3;
	  _Tp c = _Tp{1};
	  d = -x2 * x2;
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
	  _Tp sum = ff + r * q;
	  _Tp sum1 = p;
	  int i;
	  for (i = 1; i <= s_max_iter; ++i)
	    {
	      ff = (i * ff + p + q) / (i * i - mu2);
	      c *= d / _Tp(i);
	      p /= _Tp(i) - mu;
	      q /= _Tp(i) + mu;
	      const _Tp del = c * (ff + r * q);
	      sum += del;
	      const _Tp del1 = c * p - _Tp(i) * del;
	      sum1 += del1;
std::cout << "  i =" << i <<  "  del = " << del <<  "  del1 = " << del1 << '\n';
	      if (std::abs(del) < s_eps * (_Tp{1} + std::abs(sum)))
		break;
	    }
	  if (i > s_max_iter)
	    std::cout << "bessel_jn: Y-series failed to converge\n";
std::cout << "i     =" << i << '\n';
std::cout << "ff    =" << ff << '\n';
std::cout << "c     =" << c << '\n';
std::cout << "p     =" << p << '\n';
std::cout << "q     =" << q << '\n';
std::cout << "sum   =" << sum << '\n';
std::cout << "sum1  =" << sum1 << '\n';
	  _Nmu = -sum;
	  _Nnu1 = -sum1 * xi2;
	  _Npmu = mu * xi * _Nmu - _Nnu1;
	  _Jmu = _Wronski / (_Npmu - f * _Nmu);
std::cout << "Nmu   =" << _Nmu << '\n';
std::cout << "Nnu1  =" << _Nnu1 << '\n';
std::cout << "Npmu  =" << _Npmu << '\n';
std::cout << "Jmu   =" << _Jmu << '\n';
	}
      fact = _Jmu / _Jnul;
      _Jnu = fact * _Jnul1;
      _Jpnu = fact * _Jpnu1;
std::cout << "fact  =" << fact << '\n';
std::cout << "Jnu   =" << _Jnu << '\n';
std::cout << "Jpnu  =" << _Jpnu << '\n';
std::cout << "n     =" << n << '\n';
      for (int i = 1; i <= n; ++i)
	{
	  const _Tp _Nnutemp = (mu + i) * xi2 * _Nnu1 - _Nmu;
	  _Nmu = _Nnu1;
	  _Nnu1 = _Nnutemp;
std::cout << "  i =" << i << "  Nmu=" << _Nmu << "  Nnu1=" << _Nnu1 << '\n';
	}
      _Nnu = _Nmu;
      _Npnu = nu * xi * _Nmu - _Nnu1;
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
      auto Ai = emsr::airy_ai(x);
      auto Bi = emsr::airy_bi(x);
      auto Ai_gsl = gsl::airy_ai(x);
      auto Bi_gsl = gsl::airy_bi(x);
      //double Ai_tr1, Bi_tr1, Aip_tr1, Bip_tr1;
      //std::tr1::detail::airy(x, Ai_tr1, Bi_tr1, Aip_tr1, Bip_tr1);
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
      gamma_temme(x, gam1, gam2, gampl, gammi);
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
