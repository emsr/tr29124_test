

  ///
  ///
  ///
  template <typename _Tp>
  void
  __bessel_cheb( _Tp __mu, _Tp & __gam1, _Tp & __gam2, _Tp & __gampl, _Tp & __gammi )
  {

    const long double __GAMMA_E = 0.5772156649015329L;

    __gampl = _Tp(1) / ::tgamma(_Tp(1) + __mu);
    __gammi = _Tp(1) / ::tgamma(_Tp(1) - __mu);

    if (std::abs(__mu) < std::numeric_limits<_Tp>::epsilon())
      __gam1 = -_Tp(__GAMMA_E);
    else
      __gam1 = (__gammi - __gampl) / (_Tp(2) * __mu);

    __gam2 = (__gammi + __gampl) / (_Tp(2));

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  __bessel_jn( const _Tp __nu, const _Tp __x,
               _Tp & rj, _Tp & ry, _Tp & rjp, _Tp & ryp )
  {

    //if (std::isnan(__nu) || std::isnan(__x))
    //  return std::numeric_limits<_Tp>::quiet_NaN();

    if (__x == _Tp(0))
      {
        if (__nu == _Tp(0))
          {
            rj = _Tp(1);
            ry = -std::numeric_limits<_Tp>::infinity();
            rjp = _Tp(0);
            ryp = std::numeric_limits<_Tp>::infinity();
          }
        else
          {
            rj = _Tp(0);
            ry = -std::numeric_limits<_Tp>::infinity();
            //ry = ???
          }
        return;
      }

    if (__x < _Tp(0) || __nu < _Tp(0))
      throw std::domain_error( "Bad arguments in __bessel_jn." );

    const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp __fp_min = std::numeric_limits<_Tp>::min();
    const int __max_iter = 10000;
    const _Tp __x_min = _Tp(2);

    int i, l, nl;
    _Tp b, c, d, del, del1, e, f, fact, fact2,
           fact3, ff, h, p, pimu, pimu2, q, r, rjl,
           rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
           x2;

    nl = (__x < __x_min
        ? static_cast<int>(__nu + _Tp(0.5L))
        : std::max(0, static_cast<int>(__nu - __x + _Tp(1.5L))));

    const _Tp __mu = __nu - nl;
    const _Tp __mu2 = __mu * __mu;
    const _Tp __xi = _Tp(1) / __x;
    const _Tp __xi2 = _Tp(2) * __xi;
    const _Tp w = __xi2 / M_PI;
    int isign = 1;
    h = __nu * __xi;
    if ( h < __fp_min )
      h = __fp_min;
    b = __xi2 * __nu;
    d = _Tp(0);
    c = h;
    for (i = 1; i <= __max_iter; ++i)
      {
        b += __xi2;
        d = b - d;
        if (std::abs(d) < __fp_min)
          d = __fp_min;
        c = b - _Tp(1) / c;
        if (std::abs(c) < __fp_min)
          c = __fp_min;
        d = _Tp(1) / d;
        del = c * d;
        h = del * h;
        if (d < _Tp(0))
          isign = -isign;
        if (std::abs(del - _Tp(1)) < __eps)
          break;
      }
    if (i > __max_iter)
      throw std::runtime_error( "Argument x too large in __bessel_jn; try asymptotic expansion." );
    rjl = isign * __fp_min;
    rjpl = h * rjl;
    rjl1 = rjl;
    rjp1 = rjpl;
    fact = __nu * __xi;
    for (l = nl; l >= 1; --l)
      {
        rjtemp = fact * rjl + rjpl;
        fact -= __xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
      }
    if (rjl == _Tp(0))
      rjl = __eps;
    f= rjpl / rjl;
    if (__x < __x_min)
      {
        x2 = _Tp(0.5L) * __x;
        pimu = M_PI * __mu;
        fact = ( std::abs(pimu) < __eps ? _Tp(1) : pimu / std::sin(pimu) );
        d = -std::log(x2);
        e = __mu * d;
        fact2 = (std::abs(e) < __eps ? _Tp(1) : sinh(e) / e);
        _Tp gam1, gam2, gampl, gammi;
        __bessel_cheb( __mu, gam1, gam2, gampl, gammi );
        ff = (_Tp(2) / M_PI) * fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        e = std::exp(e);
        p = e / (M_PI * gampl);
        q = _Tp(1) / (e * M_PI * gammi);
        pimu2 = _Tp(0.5L) * pimu;
        fact3 = (std::abs(pimu2) < __eps ? _Tp(1) : std::sin(pimu2) / pimu2);
        r = M_PI * pimu2 * fact3 * fact3;
        c = _Tp(1);
        d = -x2 * x2;
        sum = ff + r * q;
        sum1 = p;
        for (i = 1; i <= __max_iter; ++i)
          {
            ff = (i * ff + p + q) / (i * i - __mu2);
            c *= d / i;
            p /= i - __mu;
            q /= i + __mu;
            del = c * (ff + r * q);
            sum += del; 
            del1 = c * p - i * del;
            sum1 += del1;
            if (std::abs(del) < __eps * (_Tp(1) + std::abs(sum)))
              break;
          }
        if (i > __max_iter)
          throw std::runtime_error("Bessel y series failed to converge in __bessel_jn." );
        rymu = -sum;
        ry1 = -sum1 * __xi2;
        rymup = __mu * __xi * rymu - ry1;
        rjmu = w / (rymup - f * rymu);
      }
    else 
      {
        _Tp a = _Tp(0.25L) - __mu2;
        _Tp q = _Tp(1);
        _Tp p = -_Tp(0.5L) * __xi;
        _Tp br = _Tp(2) * __x;
        _Tp bi = _Tp(2);
        _Tp fact = a * __xi / (p * p + q * q);
        _Tp cr = br + q * fact;
        _Tp ci = bi + p * fact;
        _Tp den = br * br + bi * bi;
        _Tp dr = br / den;
        _Tp di = -bi / den;
        _Tp dlr = cr * dr - ci * di;
        _Tp dli = cr * di + ci * dr;
        _Tp temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        for (i = 2; i <= __max_iter; ++i)
          {
            a += 2 * (i - 1);
            bi += _Tp(2);
            dr = a * dr + br;
            di = a * di + bi;
            if (std::abs(dr) + std::abs(di) < __fp_min)
              dr = __fp_min;
            fact = a / (cr * cr + ci * ci);
            cr = br + cr * fact;
            ci = bi - ci * fact;
            if (std::abs(cr) + std::abs(ci) < __fp_min)
              cr = __fp_min;
            den = dr * dr + di * di;
            dr /= den;
            di /= -den;
            dlr = cr * dr - ci * di;
            dli = cr * di + ci * dr;
            temp = p * dlr - q * dli;
            q = p * dli + q * dlr;
            p = temp;
            if (std::abs(dlr - _Tp(1)) + std::abs(dli) < __eps)
              break;
          }
        if (i > __max_iter)
          throw std::runtime_error( "Lentz's method failed in __bessel_jn." );
        _Tp gam = (p - f) / q;
        rjmu = std::sqrt(w / ((p - f) * gam + q));
        rjmu = ::copysign( rjmu, rjl );
        rymu = rjmu * gam;
        rymup = rymu * (p + q / gam);
        ry1 = __mu * __xi * rymu - rymup;
    }
    fact = rjmu / rjl;
    rj = rjl1 * fact;
    rjp = rjp1 * fact;
    for (i = 1; i <= nl; ++i)
      {
        rytemp = (__mu + i) * __xi2 * ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
      }
    ry = rymu;
    ryp = __nu * __xi * rymu - ry1;

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  __bessel_ik( _Tp __nu, _Tp x,
               _Tp & I_nu, _Tp & K_nu, _Tp & Ip_nu, _Tp & Kp_nu )
  {

    //if (std::isnan(__nu) || std::isnan(__x))
    //  return std::numeric_limits<_Tp>::quiet_NaN();

    if (x == _Tp(0))
      {
        if (__nu == _Tp(0))
          {
            I_nu = _Tp(1);
            K_nu = std::numeric_limits<_Tp>::infinity();
            Ip_nu = _Tp(0);
            Kp_nu = -std::numeric_limits<_Tp>::infinity();
          }
        else
          {
            I_nu = _Tp(0);
            K_nu = std::numeric_limits<_Tp>::infinity();
            //Ip_nu = ???
            //Kp_nu = ???
          }
        return;
      }

    if (x < _Tp(0) || __nu < _Tp(0))
      throw std::domain_error( "Bad arguments in __bessel_ik." );

    const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp __fp_min = std::numeric_limits<_Tp>::min();
    const int __max_iter = 10000;
    const _Tp __x_min = _Tp(2);

    int i, l, nl;
    _Tp a, a1, del, del1, delh, dels, e, f, fact, fact2, ff,
      gam1, gam2, gammi, gampl, p, pimu, q, q1, q2, qnew, ril, ril1, rimu,
      rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2;

    if (x <= _Tp(0) || __nu < _Tp(0))
      throw std::domain_error( "Bad arguments in bessel_ik." );

    nl = static_cast<int>(__nu + _Tp(0.5L));

    const _Tp __mu = __nu - nl;
    const _Tp __mu2 = __mu * __mu;
    const _Tp __xi = _Tp(1) / x;
    const _Tp __xi2 = _Tp(2) * __xi;
    _Tp h = __nu * __xi;
    if (h < __fp_min)
      h = __fp_min;
    _Tp b = __xi2 * __nu;
    _Tp d = _Tp(0);
    _Tp c = h;
    for (i = 1; i <= __max_iter; ++i)
      {
        b += __xi2;
        d = _Tp(1) / (b + d);
        c = b + _Tp(1) / c;
        del = c * d;
        h = del * h;
        if (std::abs(del - _Tp(1)) < __eps)
          break;
      }
    if (i > __max_iter)
      throw std::runtime_error( "Argument x too large in bessel_ik; try asymptotic expansion." );
    ril = __fp_min;
    ripl = h * ril;
    ril1 = ril;
    rip1 = ripl;
    fact = __nu * __xi;
    for (l = nl; l >= 1; --l)
      {
        ritemp = fact * ril + ripl;
        fact -= __xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
      }
    f= ripl/ril;
    if (x < __x_min)
      {
        x2 = _Tp(0.5L) * x;
        pimu = M_PI * __mu;
        fact = (std::abs(pimu) < __eps ? _Tp(1) : pimu / std::sin(pimu));
        d = -std::log(x2);
        e = __mu * d;
        fact2 = (std::abs(e) < __eps ? _Tp(1) : sinh(e) / e);
        __bessel_cheb( __mu, gam1, gam2, gampl, gammi );
        ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        sum = ff;
        e = std::exp(e);
        p = _Tp(0.5L) * e / gampl;
        q = _Tp(0.5L) / (e * gammi);
        c = _Tp(1);
        d = x2 * x2;
        sum1 = p;
        for (i = 1; i <= __max_iter; ++i)
          {
            ff = (i * ff + p + q) / (i * i - __mu2);
            c *= d / i;
            p /= i - __mu;
            q /= i + __mu;
            del = c * ff;
            sum += del; 
            del1 = c * (p - i * ff);
            sum1 += del1;
            if ( std::abs(del) < __eps * std::abs(sum) )
              break;
          }
        if (i > __max_iter)
          throw std::runtime_error("Bessel k series failed to converge in __bessel_jn." );
        rkmu = sum;
        rk1 = sum1 * __xi2;
      }
    else
      {
        b = _Tp(2) * (_Tp(1) + x);
        d = _Tp(1) / b;
        h = delh = d;
        q1 = _Tp(0);
        q2 = _Tp(1);
        a1 = _Tp(0.25L) - __mu2;
        q = c = a1;
        a = -a1;
        s = _Tp(1) + q * delh;
        for (i = 2; i <= __max_iter; ++i)
          {
            a -= 2 * (i - 1);
            c = -a * c / i;
            qnew = (q1 - b * q2) / a;
            q1 = q2;
            q2 = qnew;
            q += c * qnew;
            b += _Tp(2);
            d = _Tp(1) / (b + a * d);
            delh = (b * d - _Tp(1)) * delh;
            h += delh;
            dels = q * delh;
            s += dels;
            if (std::abs(dels / s) < __eps)
              break;
          }
        if (i > __max_iter)
          throw std::runtime_error( "Steed's method failed in bessel_ik." );
        h = a1 * h;
        rkmu = std::sqrt(M_PI / (_Tp(2) * x)) * std::exp(-x) / s;
        rk1 = rkmu * (__mu + x + _Tp(0.5L) - h) * __xi;
      }
    rkmup = __mu * __xi * rkmu - rk1;
    rimu = __xi / (f * rkmu - rkmup);
    I_nu = rimu * ril1 / ril;
    Ip_nu = rimu * rip1 / ril;
    for (i = 1; i <= nl; ++i)
      {
        rktemp = (__mu + i) * __xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
      }
    K_nu = rkmu;
    Kp_nu = __nu * __xi * rkmu - rk1;

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  __airy( _Tp __x, _Tp & __Ai, _Tp & __Bi, _Tp & __Aip, _Tp & __Bip )
  {

    const _Tp __SQRT3 = std::sqrt(_Tp(3));
    _Tp __absx = std::abs(__x);
    _Tp __rootx = std::sqrt(__absx);
    _Tp __z = _Tp(2) * __absx * __rootx / _Tp(3);

    if (__x > _Tp(0))
      {
        _Tp __I_nu, __K_nu, __Ip_nu, __Kp_nu;

        __bessel_ik( _Tp(1) / _Tp(3), __z, __I_nu, __K_nu, __Ip_nu, __Kp_nu );
        __Ai = __rootx * __K_nu / (__SQRT3 * M_PI);
        __Bi = __rootx * (__K_nu / M_PI + _Tp(2) * __I_nu / __SQRT3);

        __bessel_ik( _Tp(2) / _Tp(3), __z, __I_nu, __K_nu, __Ip_nu, __Kp_nu );
        __Aip = -__x * __K_nu / (__SQRT3 * M_PI);
        __Bip = __x * (__K_nu / M_PI + _Tp(2) * __I_nu / __SQRT3);
      }
    else if (__x < _Tp(0))
      {

        _Tp __J_nu, _N_nu, __Jp_nu, __Np_nu;

        __bessel_jn( _Tp(1) / _Tp(3), __z, __J_nu, _N_nu, __Jp_nu, __Np_nu );
        __Ai = _Tp(0.5L) * __rootx * (__J_nu - _N_nu / __SQRT3);
        __Bi = -_Tp(0.5L) * __rootx * (_N_nu + __J_nu / __SQRT3);

        __bessel_jn( _Tp(2) / _Tp(3), __z, __J_nu, _N_nu, __Jp_nu, __Np_nu );
        __Aip = _Tp(0.5L) * __absx * (_N_nu / __SQRT3 + __J_nu);
        __Bip = _Tp(0.5L) * __absx * (__J_nu / __SQRT3 - _N_nu);
      }
    else
      {
        // References : Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
        // The number is Ai(0) or 3**(-2/3)/Gamma(2/3).
        __Ai = 0.35502805388781723926L;
        __Bi = __Ai * __SQRT3;

        // References : Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
        // The number is Ai'(0) or -3**(-1/3)/Gamma(1/3)
        __Aip = -0.25881940379280679840L;
        __Bip = -__Aip * __SQRT3;
      }

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  __sph_bessel( const int __n, const _Tp __x,
                   _Tp & __j_n, _Tp & __n_n, _Tp & __jp_n, _Tp & __np_n )
  {

    if (__n < 0 || __x < _Tp(0))
      throw std::domain_error( "Bad arguments in sph_bessel." );

    _Tp __nu = _Tp(__n) + _Tp(0.5L);

    _Tp __J_nu, __N_nu, __Jp_nu, __Np_nu;
    __bessel_jn( __x, __nu, __J_nu, __N_nu, __Jp_nu, __Np_nu );

    const _Tp SQRTPIO2 = std::sqrt(_Tp(M_PI / 2));
    _Tp __factor = SQRTPIO2 / std::sqrt(__x);
    __j_n = __factor * __J_nu;
    __jp_n = __factor * __Jp_nu - __j_n / (_Tp(2) * __x);
    __n_n = __factor * __N_nu;
    __np_n = __factor * __Np_nu - __n_n / (_Tp(2) * __x);

    return;
  }


