

  ///
  ///
  ///
  template <typename _Tp>
  void
  bessel_cheb( _Tp mu, _Tp & gam1, _Tp & gam2, _Tp & gampl, _Tp & gammi )
  {

    const long double GAMMA_E = 0.5772156649015329L;

    gampl = _Tp(1) / ::tgamma(_Tp(1) + mu);
    gammi = _Tp(1) / ::tgamma(_Tp(1) - mu);

    if (std::abs(mu) < std::numeric_limits<_Tp>::epsilon())
      gam1 = -_Tp(GAMMA_E);
    else
      gam1 = (gammi - gampl) / (_Tp(2) * mu);

    gam2 = (gammi + gampl) / (_Tp(2));

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  bessel_jn( const _Tp nu, const _Tp x,
               _Tp & rj, _Tp & ry, _Tp & rjp, _Tp & ryp )
  {

    //if (std::isnan(nu) || std::isnan(x))
    //  return std::numeric_limits<_Tp>::quiet_NaN();

    if (x == _Tp(0))
      {
        if (nu == _Tp(0))
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

    if (x < _Tp(0) || nu < _Tp(0))
      throw std::domain_error( "Bad arguments in bessel_jn." );

    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp fp_min = std::numeric_limits<_Tp>::min();
    const int max_iter = 10000;
    const _Tp x_min = _Tp(2);

    int i, l, nl;
    _Tp b, c, d, del, del1, e, f, fact, fact2,
           fact3, ff, h, p, pimu, pimu2, q, r, rjl,
           rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
           x2;

    nl = (x < x_min
        ? static_cast<int>(nu + _Tp(0.5L))
        : std::max(0, static_cast<int>(nu - x + _Tp(1.5L))));

    const _Tp mu = nu - nl;
    const _Tp mu2 = mu * mu;
    const _Tp xi = _Tp(1) / x;
    const _Tp xi2 = _Tp(2) * xi;
    const _Tp w = xi2 / M_PI;
    int isign = 1;
    h = nu * xi;
    if ( h < fp_min )
      h = fp_min;
    b = xi2 * nu;
    d = _Tp(0);
    c = h;
    for (i = 1; i <= max_iter; ++i)
      {
        b += xi2;
        d = b - d;
        if (std::abs(d) < fp_min)
          d = fp_min;
        c = b - _Tp(1) / c;
        if (std::abs(c) < fp_min)
          c = fp_min;
        d = _Tp(1) / d;
        del = c * d;
        h = del * h;
        if (d < _Tp(0))
          isign = -isign;
        if (std::abs(del - _Tp(1)) < eps)
          break;
      }
    if (i > max_iter)
      throw std::runtime_error( "Argument x too large in bessel_jn; try asymptotic expansion." );
    rjl = isign * fp_min;
    rjpl = h * rjl;
    rjl1 = rjl;
    rjp1 = rjpl;
    fact = nu * xi;
    for (l = nl; l >= 1; --l)
      {
        rjtemp = fact * rjl + rjpl;
        fact -= xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
      }
    if (rjl == _Tp(0))
      rjl = eps;
    f= rjpl / rjl;
    if (x < x_min)
      {
        x2 = _Tp(0.5L) * x;
        pimu = M_PI * mu;
        fact = ( std::abs(pimu) < eps ? _Tp(1) : pimu / std::sin(pimu) );
        d = -std::log(x2);
        e = mu * d;
        fact2 = (std::abs(e) < eps ? _Tp(1) : sinh(e) / e);
        _Tp gam1, gam2, gampl, gammi;
        bessel_cheb( mu, gam1, gam2, gampl, gammi );
        ff = (_Tp(2) / M_PI) * fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        e = std::exp(e);
        p = e / (M_PI * gampl);
        q = _Tp(1) / (e * M_PI * gammi);
        pimu2 = _Tp(0.5L) * pimu;
        fact3 = (std::abs(pimu2) < eps ? _Tp(1) : std::sin(pimu2) / pimu2);
        r = M_PI * pimu2 * fact3 * fact3;
        c = _Tp(1);
        d = -x2 * x2;
        sum = ff + r * q;
        sum1 = p;
        for (i = 1; i <= max_iter; ++i)
          {
            ff = (i * ff + p + q) / (i * i - mu2);
            c *= d / i;
            p /= i - mu;
            q /= i + mu;
            del = c * (ff + r * q);
            sum += del; 
            del1 = c * p - i * del;
            sum1 += del1;
            if (std::abs(del) < eps * (_Tp(1) + std::abs(sum)))
              break;
          }
        if (i > max_iter)
          throw std::runtime_error("Bessel y series failed to converge in bessel_jn." );
        rymu = -sum;
        ry1 = -sum1 * xi2;
        rymup = mu * xi * rymu - ry1;
        rjmu = w / (rymup - f * rymu);
      }
    else 
      {
        _Tp a = _Tp(0.25L) - mu2;
        _Tp q = _Tp(1);
        _Tp p = -_Tp(0.5L) * xi;
        _Tp br = _Tp(2) * x;
        _Tp bi = _Tp(2);
        _Tp fact = a * xi / (p * p + q * q);
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
        for (i = 2; i <= max_iter; ++i)
          {
            a += 2 * (i - 1);
            bi += _Tp(2);
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
            if (std::abs(dlr - _Tp(1)) + std::abs(dli) < eps)
              break;
          }
        if (i > max_iter)
          throw std::runtime_error( "Lentz's method failed in bessel_jn." );
        _Tp gam = (p - f) / q;
        rjmu = std::sqrt(w / ((p - f) * gam + q));
        rjmu = ::copysign( rjmu, rjl );
        rymu = rjmu * gam;
        rymup = rymu * (p + q / gam);
        ry1 = mu * xi * rymu - rymup;
    }
    fact = rjmu / rjl;
    rj = rjl1 * fact;
    rjp = rjp1 * fact;
    for (i = 1; i <= nl; ++i)
      {
        rytemp = (mu + i) * xi2 * ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
      }
    ry = rymu;
    ryp = nu * xi * rymu - ry1;

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  bessel_ik( _Tp nu, _Tp x,
               _Tp & I_nu, _Tp & K_nu, _Tp & Ip_nu, _Tp & Kp_nu )
  {

    //if (std::isnan(nu) || std::isnan(x))
    //  return std::numeric_limits<_Tp>::quiet_NaN();

    if (x == _Tp(0))
      {
        if (nu == _Tp(0))
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

    if (x < _Tp(0) || nu < _Tp(0))
      throw std::domain_error( "Bad arguments in bessel_ik." );

    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp fp_min = std::numeric_limits<_Tp>::min();
    const int max_iter = 10000;
    const _Tp x_min = _Tp(2);

    int i, l, nl;
    _Tp a, a1, del, del1, delh, dels, e, f, fact, fact2, ff,
      gam1, gam2, gammi, gampl, p, pimu, q, q1, q2, qnew, ril, ril1, rimu,
      rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2;

    if (x <= _Tp(0) || nu < _Tp(0))
      throw std::domain_error( "Bad arguments in bessel_ik." );

    nl = static_cast<int>(nu + _Tp(0.5L));

    const _Tp mu = nu - nl;
    const _Tp mu2 = mu * mu;
    const _Tp xi = _Tp(1) / x;
    const _Tp xi2 = _Tp(2) * xi;
    _Tp h = nu * xi;
    if (h < fp_min)
      h = fp_min;
    _Tp b = xi2 * nu;
    _Tp d = _Tp(0);
    _Tp c = h;
    for (i = 1; i <= max_iter; ++i)
      {
        b += xi2;
        d = _Tp(1) / (b + d);
        c = b + _Tp(1) / c;
        del = c * d;
        h = del * h;
        if (std::abs(del - _Tp(1)) < eps)
          break;
      }
    if (i > max_iter)
      throw std::runtime_error( "Argument x too large in bessel_ik; try asymptotic expansion." );
    ril = fp_min;
    ripl = h * ril;
    ril1 = ril;
    rip1 = ripl;
    fact = nu * xi;
    for (l = nl; l >= 1; --l)
      {
        ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
      }
    f= ripl/ril;
    if (x < x_min)
      {
        x2 = _Tp(0.5L) * x;
        pimu = M_PI * mu;
        fact = (std::abs(pimu) < eps ? _Tp(1) : pimu / std::sin(pimu));
        d = -std::log(x2);
        e = mu * d;
        fact2 = (std::abs(e) < eps ? _Tp(1) : sinh(e) / e);
        bessel_cheb( mu, gam1, gam2, gampl, gammi );
        ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        sum = ff;
        e = std::exp(e);
        p = _Tp(0.5L) * e / gampl;
        q = _Tp(0.5L) / (e * gammi);
        c = _Tp(1);
        d = x2 * x2;
        sum1 = p;
        for (i = 1; i <= max_iter; ++i)
          {
            ff = (i * ff + p + q) / (i * i - mu2);
            c *= d / i;
            p /= i - mu;
            q /= i + mu;
            del = c * ff;
            sum += del; 
            del1 = c * (p - i * ff);
            sum1 += del1;
            if ( std::abs(del) < eps * std::abs(sum) )
              break;
          }
        if (i > max_iter)
          throw std::runtime_error("Bessel k series failed to converge in bessel_jn." );
        rkmu = sum;
        rk1 = sum1 * xi2;
      }
    else
      {
        b = _Tp(2) * (_Tp(1) + x);
        d = _Tp(1) / b;
        h = delh = d;
        q1 = _Tp(0);
        q2 = _Tp(1);
        a1 = _Tp(0.25L) - mu2;
        q = c = a1;
        a = -a1;
        s = _Tp(1) + q * delh;
        for (i = 2; i <= max_iter; ++i)
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
            if (std::abs(dels / s) < eps)
              break;
          }
        if (i > max_iter)
          throw std::runtime_error( "Steed's method failed in bessel_ik." );
        h = a1 * h;
        rkmu = std::sqrt(M_PI / (_Tp(2) * x)) * std::exp(-x) / s;
        rk1 = rkmu * (mu + x + _Tp(0.5L) - h) * xi;
      }
    rkmup = mu * xi * rkmu - rk1;
    rimu = xi / (f * rkmu - rkmup);
    I_nu = rimu * ril1 / ril;
    Ip_nu = rimu * rip1 / ril;
    for (i = 1; i <= nl; ++i)
      {
        rktemp = (mu + i) * xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
      }
    K_nu = rkmu;
    Kp_nu = nu * xi * rkmu - rk1;

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  airy( _Tp x, _Tp & Ai, _Tp & Bi, _Tp & Aip, _Tp & Bip )
  {

    const _Tp SQRT3 = std::sqrt(_Tp(3));
    _Tp absx = std::abs(x);
    _Tp rootx = std::sqrt(absx);
    _Tp z = _Tp(2) * absx * rootx / _Tp(3);

    if (x > _Tp(0))
      {
        _Tp I_nu, K_nu, Ip_nu, Kp_nu;

        bessel_ik( _Tp(1) / _Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu );
        Ai = rootx * K_nu / (SQRT3 * M_PI);
        Bi = rootx * (K_nu / M_PI + _Tp(2) * I_nu / SQRT3);

        bessel_ik( _Tp(2) / _Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu );
        Aip = -x * K_nu / (SQRT3 * M_PI);
        Bip = x * (K_nu / M_PI + _Tp(2) * I_nu / SQRT3);
      }
    else if (x < _Tp(0))
      {

        _Tp J_nu, _N_nu, Jp_nu, Np_nu;

        bessel_jn( _Tp(1) / _Tp(3), z, J_nu, _N_nu, Jp_nu, Np_nu );
        Ai = _Tp(0.5L) * rootx * (J_nu - _N_nu / SQRT3);
        Bi = -_Tp(0.5L) * rootx * (_N_nu + J_nu / SQRT3);

        bessel_jn( _Tp(2) / _Tp(3), z, J_nu, _N_nu, Jp_nu, Np_nu );
        Aip = _Tp(0.5L) * absx * (_N_nu / SQRT3 + J_nu);
        Bip = _Tp(0.5L) * absx * (J_nu / SQRT3 - _N_nu);
      }
    else
      {
        // References : Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
        // The number is Ai(0) or 3**(-2/3)/Gamma(2/3).
        Ai = 0.35502805388781723926L;
        Bi = Ai * SQRT3;

        // References : Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
        // The number is Ai'(0) or -3**(-1/3)/Gamma(1/3)
        Aip = -0.25881940379280679840L;
        Bip = -Aip * SQRT3;
      }

    return;
  }


  ///
  ///
  ///
  template <typename _Tp>
  void
  sph_bessel( const int n, const _Tp x,
                   _Tp & j_n, _Tp & n_n, _Tp & jp_n, _Tp & np_n )
  {

    if (n < 0 || x < _Tp(0))
      throw std::domain_error( "Bad arguments in sph_bessel." );

    _Tp nu = _Tp(n) + _Tp(0.5L);

    _Tp J_nu, N_nu, Jp_nu, Np_nu;
    bessel_jn( x, nu, J_nu, N_nu, Jp_nu, Np_nu );

    const _Tp SQRTPIO2 = std::sqrt(_Tp(M_PI / 2));
    _Tp factor = SQRTPIO2 / std::sqrt(x);
    j_n = factor * J_nu;
    jp_n = factor * Jp_nu - j_n / (_Tp(2) * x);
    n_n = factor * N_nu;
    np_n = factor * Np_nu - n_n / (_Tp(2) * x);

    return;
  }


