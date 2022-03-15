

  ///
  ///
  ///
  template <typename Tp>
  void
  bessel_cheb( Tp mu, Tp & gam1, Tp & gam2, Tp & gampl, Tp & gammi )
  {

    const long double GAMMA_E = 0.5772156649015329L;

    gampl = Tp(1) / ::tgamma(Tp(1) + mu);
    gammi = Tp(1) / ::tgamma(Tp(1) - mu);

    if (std::abs(mu) < std::numeric_limits<Tp>::epsilon())
      gam1 = -Tp(GAMMA_E);
    else
      gam1 = (gammi - gampl) / (Tp(2) * mu);

    gam2 = (gammi + gampl) / (Tp(2));

    return;
  }


  ///
  ///
  ///
  template <typename Tp>
  void
  bessel_jn( const Tp nu, const Tp x,
               Tp & rj, Tp & ry, Tp & rjp, Tp & ryp )
  {

    //if (std::isnan(nu) || std::isnan(x))
    //  return std::numeric_limits<Tp>::quiet_NaN();

    if (x == Tp(0))
      {
        if (nu == Tp(0))
          {
            rj = Tp(1);
            ry = -std::numeric_limits<Tp>::infinity();
            rjp = Tp(0);
            ryp = std::numeric_limits<Tp>::infinity();
          }
        else
          {
            rj = Tp(0);
            ry = -std::numeric_limits<Tp>::infinity();
            //ry = ???
          }
        return;
      }

    if (x < Tp(0) || nu < Tp(0))
      throw std::domain_error( "Bad arguments in bessel_jn." );

    const Tp eps = std::numeric_limits<Tp>::epsilon();
    const Tp fp_min = std::numeric_limits<Tp>::min();
    const int max_iter = 10000;
    const Tp x_min = Tp(2);

    int i, l, nl;
    Tp b, c, d, del, del1, e, f, fact, fact2,
           fact3, ff, h, p, pimu, pimu2, q, r, rjl,
           rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
           x2;

    nl = (x < x_min
        ? static_cast<int>(nu + Tp(0.5L))
        : std::max(0, static_cast<int>(nu - x + Tp(1.5L))));

    const Tp mu = nu - nl;
    const Tp mu2 = mu * mu;
    const Tp xi = Tp(1) / x;
    const Tp xi2 = Tp(2) * xi;
    const Tp w = xi2 / M_PI;
    int isign = 1;
    h = nu * xi;
    if ( h < fp_min )
      h = fp_min;
    b = xi2 * nu;
    d = Tp(0);
    c = h;
    for (i = 1; i <= max_iter; ++i)
      {
        b += xi2;
        d = b - d;
        if (std::abs(d) < fp_min)
          d = fp_min;
        c = b - Tp(1) / c;
        if (std::abs(c) < fp_min)
          c = fp_min;
        d = Tp(1) / d;
        del = c * d;
        h = del * h;
        if (d < Tp(0))
          isign = -isign;
        if (std::abs(del - Tp(1)) < eps)
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
    if (rjl == Tp(0))
      rjl = eps;
    f= rjpl / rjl;
    if (x < x_min)
      {
        x2 = Tp(0.5L) * x;
        pimu = M_PI * mu;
        fact = ( std::abs(pimu) < eps ? Tp(1) : pimu / std::sin(pimu) );
        d = -std::log(x2);
        e = mu * d;
        fact2 = (std::abs(e) < eps ? Tp(1) : sinh(e) / e);
        Tp gam1, gam2, gampl, gammi;
        bessel_cheb( mu, gam1, gam2, gampl, gammi );
        ff = (Tp(2) / M_PI) * fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        e = std::exp(e);
        p = e / (M_PI * gampl);
        q = Tp(1) / (e * M_PI * gammi);
        pimu2 = Tp(0.5L) * pimu;
        fact3 = (std::abs(pimu2) < eps ? Tp(1) : std::sin(pimu2) / pimu2);
        r = M_PI * pimu2 * fact3 * fact3;
        c = Tp(1);
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
            if (std::abs(del) < eps * (Tp(1) + std::abs(sum)))
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
        Tp a = Tp(0.25L) - mu2;
        Tp q = Tp(1);
        Tp p = -Tp(0.5L) * xi;
        Tp br = Tp(2) * x;
        Tp bi = Tp(2);
        Tp fact = a * xi / (p * p + q * q);
        Tp cr = br + q * fact;
        Tp ci = bi + p * fact;
        Tp den = br * br + bi * bi;
        Tp dr = br / den;
        Tp di = -bi / den;
        Tp dlr = cr * dr - ci * di;
        Tp dli = cr * di + ci * dr;
        Tp temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        for (i = 2; i <= max_iter; ++i)
          {
            a += 2 * (i - 1);
            bi += Tp(2);
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
            if (std::abs(dlr - Tp(1)) + std::abs(dli) < eps)
              break;
          }
        if (i > max_iter)
          throw std::runtime_error( "Lentz's method failed in bessel_jn." );
        Tp gam = (p - f) / q;
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
  template <typename Tp>
  void
  bessel_ik( Tp nu, Tp x,
               Tp & I_nu, Tp & K_nu, Tp & Ip_nu, Tp & Kp_nu )
  {

    //if (std::isnan(nu) || std::isnan(x))
    //  return std::numeric_limits<Tp>::quiet_NaN();

    if (x == Tp(0))
      {
        if (nu == Tp(0))
          {
            I_nu = Tp(1);
            K_nu = std::numeric_limits<Tp>::infinity();
            Ip_nu = Tp(0);
            Kp_nu = -std::numeric_limits<Tp>::infinity();
          }
        else
          {
            I_nu = Tp(0);
            K_nu = std::numeric_limits<Tp>::infinity();
            //Ip_nu = ???
            //Kp_nu = ???
          }
        return;
      }

    if (x < Tp(0) || nu < Tp(0))
      throw std::domain_error( "Bad arguments in bessel_ik." );

    const Tp eps = std::numeric_limits<Tp>::epsilon();
    const Tp fp_min = std::numeric_limits<Tp>::min();
    const int max_iter = 10000;
    const Tp x_min = Tp(2);

    int i, l, nl;
    Tp a, a1, del, del1, delh, dels, e, f, fact, fact2, ff,
      gam1, gam2, gammi, gampl, p, pimu, q, q1, q2, qnew, ril, ril1, rimu,
      rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2;

    if (x <= Tp(0) || nu < Tp(0))
      throw std::domain_error( "Bad arguments in bessel_ik." );

    nl = static_cast<int>(nu + Tp(0.5L));

    const Tp mu = nu - nl;
    const Tp mu2 = mu * mu;
    const Tp xi = Tp(1) / x;
    const Tp xi2 = Tp(2) * xi;
    Tp h = nu * xi;
    if (h < fp_min)
      h = fp_min;
    Tp b = xi2 * nu;
    Tp d = Tp(0);
    Tp c = h;
    for (i = 1; i <= max_iter; ++i)
      {
        b += xi2;
        d = Tp(1) / (b + d);
        c = b + Tp(1) / c;
        del = c * d;
        h = del * h;
        if (std::abs(del - Tp(1)) < eps)
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
        x2 = Tp(0.5L) * x;
        pimu = M_PI * mu;
        fact = (std::abs(pimu) < eps ? Tp(1) : pimu / std::sin(pimu));
        d = -std::log(x2);
        e = mu * d;
        fact2 = (std::abs(e) < eps ? Tp(1) : sinh(e) / e);
        bessel_cheb( mu, gam1, gam2, gampl, gammi );
        ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        sum = ff;
        e = std::exp(e);
        p = Tp(0.5L) * e / gampl;
        q = Tp(0.5L) / (e * gammi);
        c = Tp(1);
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
        b = Tp(2) * (Tp(1) + x);
        d = Tp(1) / b;
        h = delh = d;
        q1 = Tp(0);
        q2 = Tp(1);
        a1 = Tp(0.25L) - mu2;
        q = c = a1;
        a = -a1;
        s = Tp(1) + q * delh;
        for (i = 2; i <= max_iter; ++i)
          {
            a -= 2 * (i - 1);
            c = -a * c / i;
            qnew = (q1 - b * q2) / a;
            q1 = q2;
            q2 = qnew;
            q += c * qnew;
            b += Tp(2);
            d = Tp(1) / (b + a * d);
            delh = (b * d - Tp(1)) * delh;
            h += delh;
            dels = q * delh;
            s += dels;
            if (std::abs(dels / s) < eps)
              break;
          }
        if (i > max_iter)
          throw std::runtime_error( "Steed's method failed in bessel_ik." );
        h = a1 * h;
        rkmu = std::sqrt(M_PI / (Tp(2) * x)) * std::exp(-x) / s;
        rk1 = rkmu * (mu + x + Tp(0.5L) - h) * xi;
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
  template <typename Tp>
  void
  airy( Tp x, Tp & Ai, Tp & Bi, Tp & Aip, Tp & Bip )
  {

    const Tp SQRT3 = std::sqrt(Tp(3));
    Tp absx = std::abs(x);
    Tp rootx = std::sqrt(absx);
    Tp z = Tp(2) * absx * rootx / Tp(3);

    if (x > Tp(0))
      {
        Tp I_nu, K_nu, Ip_nu, Kp_nu;

        bessel_ik( Tp(1) / Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu );
        Ai = rootx * K_nu / (SQRT3 * M_PI);
        Bi = rootx * (K_nu / M_PI + Tp(2) * I_nu / SQRT3);

        bessel_ik( Tp(2) / Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu );
        Aip = -x * K_nu / (SQRT3 * M_PI);
        Bip = x * (K_nu / M_PI + Tp(2) * I_nu / SQRT3);
      }
    else if (x < Tp(0))
      {

        Tp J_nu, _N_nu, Jp_nu, Np_nu;

        bessel_jn( Tp(1) / Tp(3), z, J_nu, _N_nu, Jp_nu, Np_nu );
        Ai = Tp(0.5L) * rootx * (J_nu - _N_nu / SQRT3);
        Bi = -Tp(0.5L) * rootx * (_N_nu + J_nu / SQRT3);

        bessel_jn( Tp(2) / Tp(3), z, J_nu, _N_nu, Jp_nu, Np_nu );
        Aip = Tp(0.5L) * absx * (_N_nu / SQRT3 + J_nu);
        Bip = Tp(0.5L) * absx * (J_nu / SQRT3 - _N_nu);
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
  template <typename Tp>
  void
  sph_bessel( const int n, const Tp x,
                   Tp & j_n, Tp & n_n, Tp & jp_n, Tp & np_n )
  {

    if (n < 0 || x < Tp(0))
      throw std::domain_error( "Bad arguments in sph_bessel." );

    Tp nu = Tp(n) + Tp(0.5L);

    Tp J_nu, N_nu, Jp_nu, Np_nu;
    bessel_jn( x, nu, J_nu, N_nu, Jp_nu, Np_nu );

    const Tp SQRTPIO2 = std::sqrt(Tp(M_PI / 2));
    Tp factor = SQRTPIO2 / std::sqrt(x);
    j_n = factor * J_nu;
    jp_n = factor * Jp_nu - j_n / (Tp(2) * x);
    n_n = factor * N_nu;
    np_n = factor * Np_nu - n_n / (Tp(2) * x);

    return;
  }


