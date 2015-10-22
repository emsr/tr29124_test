

#include <cmath>
#include <stdexcept>
//#include <tr1/special_function_util.h>
    //  pi.
    const long double __PI       = 3.1415926535897932384626433832795029L;
    //  pi / 2.
    const long double __PI_2     = 1.5707963267948966192313216916397514L;
    //  pi / 3.
    const long double __PI_3     = 1.0471975511965977461542144610932L;
    //  pi / 4.
    const long double __PI_4     = 0.7853981633974483096156608458198757L;
    //  1 / pi.
    const long double __1_PI     = 0.3183098861837906715377675267450287L;
    //  2 / sqrt(pi).
    const long double __2_SQRTPI = 1.1283791670955125738961589031215452L;
    //  sqrt(2).
    const long double __SQRT2    = 1.4142135623730950488016887242096981L;
    //  sqrt(3).
    const long double __SQRT3    = 1.7320508075688772935274463415058723L;
    //  sqrt(pi/2).
    const long double __SQRTPIO2 = 1.2533141373155002512078826424055226L;
    //  1 / sqrt(2).
    const long double __SQRT1_2  = 0.7071067811865475244008443621048490L;
    //  log(pi).
    const long double __LNPI     = 1.1447298858494001741434273513531L;
    //  Euler's constant.
    const long double __GAMMA_E  = 0.57721566490153286060651209008240L;
    //  Euler-Mascheroni
    const long double __EULER    = 2.7182818284590452353602874713526625L;


    ///
    ///  @brief Compute the gamma functions required by th Temme series
    ///         expansions of @f$ N_nu(x) @f$ and @f$ K_nu(x) @f$.
    ///
    template <typename _Tp>
    void
    __gamma_temme(const _Tp __mu,
                   _Tp & __gam1, _Tp & __gam2, _Tp & __gampl, _Tp & __gammi)
    {
      __gampl = _Tp(1) / ::tgamma(_Tp(1) + __mu);
      __gammi = _Tp(1) / ::tgamma(_Tp(1) - __mu);

      if (std::abs(__mu) < std::numeric_limits<_Tp>::epsilon())
        __gam1 = -_Tp(__GAMMA_E);
      else
        __gam1 = (__gammi - __gampl) / (_Tp(2) * __mu);

      __gam2 = (__gammi + __gampl) / (_Tp(2));

      return;
    }


    template <typename _Tp>
    void
    __bessel_jn(const _Tp __nu, const _Tp __x,
                _Tp & __Jnu, _Tp & __Nnu, _Tp & __Jpnu, _Tp & __Npnu)
    {

      //if (std::isnan(__nu) || std::isnan(__x))
      //  return std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp(0))
        {
          if (__nu == _Tp(0))
            {
              __Jnu = _Tp(1);
              __Nnu = -std::numeric_limits<_Tp>::infinity();
              __Jpnu = _Tp(0);
              __Npnu = std::numeric_limits<_Tp>::infinity();
            }
          else
            {
              __Jnu = _Tp(0);
              __Nnu = -std::numeric_limits<_Tp>::infinity();
              //__Jpnu = ???
              //__Npnu = ???
            }
          return;
        }

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __fp_min = _Tp(10) * std::numeric_limits<_Tp>::min();
      const int __max_iter = 15000;
      const _Tp __x_min = _Tp(2);

      if (__x < _Tp(0) || __nu < _Tp(0))
        throw std::runtime_error("Bad arguments in __bessel_jn.");

      const int __nl = (__x < __x_min
                    ? static_cast<int>(__nu + _Tp(0.5L))
                    : std::max(0, static_cast<int>(__nu - __x + _Tp(1.5L))));

      const _Tp __mu = __nu - __nl;
      const _Tp __mu2 = __mu * __mu;
      const _Tp __xi = _Tp(1) / __x;
      const _Tp __xi2 = _Tp(2) * __xi;
      _Tp __w = __xi2 / _Tp(__PI);
      int __isign = 1;
      _Tp __h = __nu * __xi;
      if (__h < __fp_min)
        __h = __fp_min;
      _Tp __b = __xi2 * __nu;
      _Tp __d = _Tp(0);
      _Tp __c = __h;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
        {
          __b += __xi2;
          __d = __b - __d;
          if (std::abs(__d) < __fp_min)
            __d = __fp_min;
          __c = __b - _Tp(1) / __c;
          if (std::abs(__c) < __fp_min)
            __c = __fp_min;
          __d = _Tp(1) / __d;
          const _Tp __del = __c * __d;
          __h = __del * __h;
          if (__d < _Tp(0))
            __isign = -__isign;
          if (std::abs(__del - _Tp(1)) < __eps)
            break;
        }
      if (__i > __max_iter)
        throw std::runtime_error( "Argument x too large in __bessel_jn; "
                                  "try asymptotic expansion." );
      _Tp __Jnul = __isign * __fp_min;
      _Tp __Jpnul = __h * __Jnul;
      _Tp __Jnul1 = __Jnul;
      _Tp __Jpnu1 = __Jpnul;
      _Tp __fact = __nu * __xi;
      for ( int __l = __nl; __l >= 1; --__l )
        {
          const _Tp __Jnutemp = __fact * __Jnul + __Jpnul;
          __fact -= __xi;
          __Jpnul = __fact * __Jnutemp - __Jnul;
          __Jnul = __Jnutemp;
        }
      if (__Jnul == _Tp(0))
        __Jnul = __eps;
      _Tp __f= __Jpnul / __Jnul;
      _Tp __Nmu, __Nnu1, __Npmu, __Jmu;
      if (__x < __x_min)
        {
          const _Tp __x2 = __x / _Tp(2);
          const _Tp __pimu = _Tp(__PI) * __mu;
          _Tp __fact = (std::abs(__pimu) < __eps
                      ? _Tp(1) : __pimu / std::sin(__pimu));
          _Tp __d = -std::log(__x2);
          _Tp __e = __mu * __d;
          _Tp __fact2 = (std::abs(__e) < __eps
                       ? _Tp(1) : std::sinh(__e) / __e);
          _Tp __gam1, __gam2, __gampl, __gammi;
          __gamma_temme(__mu, __gam1, __gam2, __gampl, __gammi);
          _Tp __ff = (_Tp(2) / _Tp(__PI))
                   * __fact * (__gam1 * std::cosh(__e) + __gam2 * __fact2 * __d);
          __e = std::exp(__e);
          _Tp __p = __e / (_Tp(__PI) * __gampl);
          _Tp __q = _Tp(1) / (__e * _Tp(__PI) * __gammi);
          const _Tp __pimu2 = __pimu / _Tp(2);
          _Tp __fact3 = (std::abs(__pimu2) < __eps
                       ? _Tp(1) : std::sin(__pimu2) / __pimu2 );
          _Tp __r = _Tp(__PI) * __pimu2 * __fact3 * __fact3;
          _Tp __c = _Tp(1);
          __d = -__x2 * __x2;
          _Tp __sum = __ff + __r * __q;
          _Tp __sum1 = __p;
          for (__i = 1; __i <= __max_iter; ++__i)
            {
              __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
              __c *= __d / _Tp(__i);
              __p /= _Tp(__i) - __mu;
              __q /= _Tp(__i) + __mu;
              const _Tp __del = __c * (__ff + __r * __q);
              __sum += __del; 
              const _Tp __del1 = __c * __p - __i * __del;
              __sum1 += __del1;
              if ( std::abs(__del) < __eps * (_Tp(1) + std::abs(__sum)) )
                break;
            }
          if ( __i > __max_iter )
            throw std::runtime_error("Bessel y series failed to converge "
                                     "in __bessel_jn." );
          __Nmu = -__sum;
          __Nnu1 = -__sum1 * __xi2;
          __Npmu = __mu * __xi * __Nmu - __Nnu1;
          __Jmu = __w / (__Npmu - __f * __Nmu);
        }
      else
        {
          _Tp __a = _Tp(0.25L) - __mu2;
          _Tp __q = _Tp(1);
          _Tp __p = -__xi / _Tp(2);
          _Tp __br = _Tp(2) * __x;
          _Tp __bi = _Tp(2);
          _Tp __fact = __a * __xi / (__p * __p + __q * __q);
          _Tp __cr = __br + __q * __fact;
          _Tp __ci = __bi + __p * __fact;
          _Tp __den = __br * __br + __bi * __bi;
          _Tp __dr = __br / __den;
          _Tp __di = -__bi / __den;
          _Tp __dlr = __cr * __dr - __ci * __di;
          _Tp __dli = __cr * __di + __ci * __dr;
          _Tp __temp = __p * __dlr - __q * __dli;
          __q = __p * __dli + __q * __dlr;
          __p = __temp;
          int __i;
          for (__i = 2; __i <= __max_iter; ++__i)
            {
              __a += _Tp(2 * (__i - 1));
              __bi += _Tp(2);
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
              if (std::abs(__dlr - _Tp(1)) + std::abs(__dli) < __eps)
                break;
          }
          if (__i > __max_iter)
            throw std::runtime_error("Lentz's method failed in __bessel_jn.");
          const _Tp __gam = (__p - __f) / __q;
          __Jmu = std::sqrt(__w / ((__p - __f) * __gam + __q));

          __Jmu = ::copysign(__Jmu, __Jnul);

          __Nmu = __gam * __Jmu;
          __Npmu = (__p + __q / __gam) * __Nmu;
          __Nnu1 = __mu * __xi * __Nmu - __Npmu;
      }
      __fact = __Jmu / __Jnul;
      __Jnu = __fact * __Jnul1;
      __Jpnu = __fact * __Jpnu1;
      for (__i = 1; __i <= __nl; ++__i)
        {
          const _Tp __Nnutemp = (__mu + __i) * __xi2 * __Nnu1 - __Nmu;
          __Nmu = __Nnu1;
          __Nnu1 = __Nnutemp;
        }
      __Nnu = __Nmu;
      __Npnu = __nu * __xi * __Nmu - __Nnu1;

      return;
    }



    template <typename _Tp>
    void
    __bessel_ik(const _Tp __nu, const _Tp __x,
                _Tp & __Inu, _Tp & __Knu, _Tp & __Ipnu, _Tp & __Kpnu)
    {

      //if (std::isnan(__nu) || std::isnan(__x))
      //  return std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp(0))
        {
          if (__nu == _Tp(0))
            {
              __Inu = _Tp(1);
              __Knu = std::numeric_limits<_Tp>::infinity();
              __Ipnu = _Tp(0);
              __Kpnu = -std::numeric_limits<_Tp>::infinity();
            }
          else
            {
              __Inu = _Tp(0);
              __Knu = std::numeric_limits<_Tp>::infinity();
              //__Ipnu = ???
              //__Kpnu = ???
            }
          return;
        }

      if (__x < _Tp(0) || __nu < _Tp(0))
        throw std::runtime_error("Bad arguments in __bessel_ik.");

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __fp_min = _Tp(10) * std::numeric_limits<_Tp>::epsilon();
      const int __max_iter = 15000;
      const _Tp __x_min = _Tp(2);

      const int __nl = static_cast<int>(__nu + _Tp(0.5L));

      const _Tp __mu = __nu - __nl;
      const _Tp __mu2 = __mu * __mu;
      const _Tp __xi = _Tp(1) / __x;
      const _Tp __xi2 = _Tp(2) * __xi;
      _Tp __h = __nu * __xi;
      if ( __h < __fp_min )
        __h = __fp_min;
      _Tp __b = __xi2 * __nu;
      _Tp __d = _Tp(0);
      _Tp __c = __h;
      int __i;
      for ( __i = 1; __i <= __max_iter; ++__i )
        {
          __b += __xi2;
          __d = _Tp(1) / (__b + __d);
          __c = __b + _Tp(1) / __c;
          const _Tp __del = __c * __d;
          __h = __del * __h;
          if (std::abs(__del - _Tp(1)) < __eps)
            break;
        }
      if (__i > __max_iter)
        throw std::runtime_error( "Argument x too large in __bessel_jn; "
                                  "try asymptotic expansion." );
      _Tp __Inul = __fp_min;
      _Tp __Ipnul = __h * __Inul;
      _Tp __Inul1 = __Inul;
      _Tp __Ipnu1 = __Ipnul;
      _Tp __fact = __nu * __xi;
      for (int __l = __nl; __l >= 1; --__l)
        {
          const _Tp __Inutemp = __fact * __Inul + __Ipnul;
          __fact -= __xi;
          __Ipnul = __fact * __Inutemp + __Inul;
          __Inul = __Inutemp;
        }
      _Tp __f = __Ipnul / __Inul;
      _Tp __Kmu, __Knu1;
      if (__x < __x_min)
        {
          const _Tp __x2 = __x / _Tp(2);
          const _Tp __pimu = _Tp(__PI) * __mu;
          const _Tp __fact = (std::abs(__pimu) < __eps
                            ? _Tp(1) : __pimu / std::sin(__pimu));
          _Tp __d = -std::log(__x2);
          _Tp __e = __mu * __d;
          const _Tp __fact2 = (std::abs(__e) < __eps
                            ? _Tp(1) : std::sinh(__e) / __e);
          _Tp __gam1, __gam2, __gampl, __gammi;
          __gamma_temme(__mu, __gam1, __gam2, __gampl, __gammi);
          _Tp __ff = __fact * (__gam1 * std::cosh(__e) + __gam2 * __fact2 * __d);
          _Tp __sum = __ff;
          __e = std::exp(__e);
          _Tp __p = __e / (_Tp(2) * __gampl);
          _Tp __q = _Tp(1) / (_Tp(2) * __e * __gammi);
          _Tp __c = _Tp(1);
          __d = __x2 * __x2;
          _Tp __sum1 = __p;
          int __i;
          for (__i = 1; __i <= __max_iter; ++__i)
            {
              __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
              __c *= __d / __i;
              __p /= __i - __mu;
              __q /= __i + __mu;
              const _Tp __del = __c * __ff;
              __sum += __del; 
              const _Tp __del1 = __c * (__p - __i * __ff);
              __sum1 += __del1;
              if (std::abs(__del) < __eps * std::abs(__sum))
                break;
            }
          if (__i > __max_iter)
            throw std::runtime_error("Bessel k series failed to converge "
                                     "in __bessel_jn." );
          __Kmu = __sum;
          __Knu1 = __sum1 * __xi2;
        }
      else
        {
          _Tp __b = _Tp(2) * (_Tp(1) + __x);
          _Tp __d = _Tp(1) / __b;
          _Tp __delh = __d;
          _Tp __h = __delh;
          _Tp __q1 = _Tp(0);
          _Tp __q2 = _Tp(1);
          _Tp __a1 = _Tp(0.25L) - __mu2;
          _Tp __q = __c = __a1;
          _Tp __a = -__a1;
          _Tp __s = _Tp(1) + __q * __delh;
          int __i;
          for (__i = 2; __i <= __max_iter; ++__i)
            {
              __a -= 2 * (__i - 1);
              __c = -__a * __c / __i;
              const _Tp __qnew = (__q1 - __b * __q2) / __a;
              __q1 = __q2;
              __q2 = __qnew;
              __q += __c * __qnew;
              __b += _Tp(2);
              __d = _Tp(1) / (__b + __a * __d);
              __delh = (__b * __d - _Tp(1)) * __delh;
              __h += __delh;
              const _Tp __dels = __q * __delh;
              __s += __dels;
              if ( std::abs(__dels / __s) < __eps )
                break;
            }
          if (__i > __max_iter)
            throw std::runtime_error("Steed's method failed in __bessel_jn.");
          __h = __a1 * __h;
          __Kmu = std::sqrt(_Tp(__PI) / (_Tp(2) * __x)) * std::exp(-__x) / __s;
          __Knu1 = __Kmu * (__mu + __x + _Tp(0.5L) - __h) * __xi;
        }

      _Tp __Kpmu = __mu * __xi * __Kmu - __Knu1;
      _Tp __Inumu = __xi / (__f * __Kmu - __Kpmu);
      __Inu = __Inumu * __Inul1 / __Inul;
      __Ipnu = __Inumu * __Ipnu1 / __Inul;
      for ( __i = 1; __i <= __nl; ++__i )
        {
          const _Tp __Knutemp = (__mu + __i) * __xi2 * __Knu1 + __Kmu;
          __Kmu = __Knu1;
          __Knu1 = __Knutemp;
        }
      __Knu = __Kmu;
      __Kpnu = __nu * __xi * __Kmu - __Knu1;
  
      return;
    }


    ///
    ///
    ///
    template <typename _Tp>
    void
    __airy(const _Tp __x,
           _Tp & __ai, _Tp & __bi, _Tp & __aip, _Tp & __bip)
    {
      const _Tp __SQRT3 = std::sqrt(_Tp(3));
      const _Tp __absx = std::abs(__x);
      const _Tp __rootx = std::sqrt(__absx);
      const _Tp __z = _Tp(2) * __absx * __rootx / _Tp(3);
      if (__x > _Tp(0))
        {
          _Tp __Inu, __Ipnu, __Knu, __Kpnu;

          __bessel_jn(_Tp(1)/_Tp(3), __z, __Inu, __Knu, __Ipnu, __Kpnu);
          __ai = __rootx * __Knu / (_Tp(__SQRT3) * __PI);
          __bi = __rootx * (__Knu / __PI + _Tp(2) * __Inu / _Tp(__SQRT3));

          __bessel_jn(_Tp(2)/_Tp(3), __z, __Inu, __Knu, __Ipnu, __Kpnu);
          __aip = -__x * __Knu / (_Tp(__SQRT3) * __PI);
          __bip = __x * (__Knu / __PI + _Tp(2) * __Inu / _Tp(__SQRT3));
        }
      else if (__x < _Tp(0))
        {
          _Tp __Jnu, __Jpnu, __Nnu, __Npnu;

          __bessel_jn(_Tp(1)/_Tp(3), __z, __Jnu, __Nnu, __Jpnu, __Npnu);
          __ai = __rootx * (__Jnu - __Nnu / _Tp(__SQRT3)) / _Tp(2);
          __bi = -__rootx * (__Nnu + __Jnu / _Tp(__SQRT3)) / _Tp(2);

          __bessel_jn(_Tp(2)/_Tp(3), __z, __Jnu, __Nnu, __Jpnu, __Npnu);
          __aip = __absx * (__Nnu / _Tp(__SQRT3) + __Jnu) / _Tp(2);
          __bip = __absx * (__Jnu / _Tp(__SQRT3) - __Nnu) / _Tp(2);
        }
      else
        {
          // References : Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
          // The number is Ai(0) or 3**(-2/3)/Gamma(2/3).
          __ai = _Tp(0.35502805388781723926L);
          __bi = __ai * __SQRT3;

          // References : Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
          // The number is Ai'(0) or -3**(-1/3)/Gamma(1/3)
          __aip = -_Tp(0.25881940379280679840L);
          __bip = -__aip * __SQRT3;
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

      if ( __n < 0 || __x < _Tp(0) )
        throw std::runtime_error( "Bad arguments in sph_bessel." );

      const _Tp __nu = _Tp(__n) + _Tp(0.5L);

      _Tp __Jnu, __Jpnu, __Nnu, __Npnu;
      __bessel_jn( __x, __nu, __Jnu, __Nnu, __Jpnu, __Npnu );

      const _Tp __factor = _Tp(__SQRTPIO2) / std::sqrt(__x);

      __jn = __factor * __Jnu;
      __nn = __factor * __Nnu;
      __jpn = __factor * __Jpnu - __jn / (_Tp(2) * __x);
      __npn = __factor * __Npnu - __nn / (_Tp(2) * __x);

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

      if ( __n < 0 || __x < _Tp(0) )
        throw std::runtime_error( "Bad arguments in sph_bessel." );

      const _Tp __order = _Tp(__n) + _Tp(0.5L);

      _Tp __Inu, __Ipnu, __Knu, __Kpnu;
      __bessel_ik( __x, __order, __Inu, __Knu, __Ipnu, __Kpnu );

      const _Tp __factor = _Tp(__SQRTPIO2) / std::sqrt(__x);

      __in = __factor * __Inu;
      __kn = __factor * __Knu;
      __ipn = __factor * __Ipnu - __in / (_Tp(2) * __x);
      __kpn = __factor * __Kpnu - __kn / (_Tp(2) * __x);

      return;
    }


