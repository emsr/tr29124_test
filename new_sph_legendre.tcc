
    template <typename _Tp>
    _Tp __sph_legendre(const unsigned int __l, const unsigned int __m,
                       const _Tp __x)
    {

      if (__m < 0 || __l < __m || __x < -_Tp(1) || __x > _Tp(1))
        {
          std::__throw_domain_error(__N("Bad argument "
                                        "in __sph_legendre."));
        }
      else if (__m == 0)
        {
          _Tp __P = __legendre_p(__l, __x);
          _Tp __fact = std::sqrt(_Tp(2 * __l + 1) / _Tp(4 * M_PI));
          __P *= __fact;
          return __P;
        }
      else if (__x == _Tp(1) || __x == -_Tp(1))
        {
          //  m > 0 here
          return _Tp(0);
        }
      else
        {
          // m > 0 and |x| < 1 here

          // Starting value for recursion.
          // Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
          const _Tp __sgn = ( __m % 2 == 1 ? -_Tp(1) : _Tp(1));
          const _Tp __y_mp1m_factor = __x * std::sqrt(_Tp(2 * __m + 3));
#if _GLIBCXX_USE_C99_MATH_TR1
          const _Tp __lncirc = std::tr1::log1p(-__x * __x);
#else
          const _Tp __lncirc = std::log(_Tp(1) - __x * __x)
#endif
          //  Gamma(m+1/2) / Gamma(m)
          const _Tp __lnpoch = std::tr1::lgamma(_Tp(__m + _Tp(0.5L)))
                             - std::tr1::lgamma(_Tp(__m));
          const _Tp __lnpre_val = -0.25L * _TR1_M_LNPI + 0.5L * (__lnpoch + __m * __lncirc);
          const _Tp __sr = std::sqrt((_Tp(2) + _Tp(1) / __m) / _Tp(4 * _TR1_M_PI));
          _Tp __y_mm = __sgn * __sr * std::exp(__lnpre_val);
          _Tp __y_mp1m = __y_mp1m_factor * __y_mm;

          if (__l == __m)
            {
              return __y_mm;
            }
          else if (__l == __m + 1)
            {
              return __y_mp1m;
            }
          else
            {
              _Tp __y_lm = _Tp(0);

              // Compute Y_l^m, l > m+1, upward recursion on l.
              for ( int __ell = __m + 2; __ell <= __l; ++__ell)
                {
                  const _Tp __rat1 = _Tp(__ell - __m) / _Tp(ell + __m);
                  const _Tp __rat2 = _Tp(__ell - __m - 1) / _Tp(ell + m - 1);
                  const _Tp __factor1 = std::sqrt(__rat1 * _Tp(2 * __ell + 1) * _Tp(2 * __ell - 1));
                  const _Tp __factor2 = std::sqrt(__rat1 * __rat2 * _Tp(2 * __ell + 1) / _Tp(2 * __ell - 3));
                  __y_lm = (__x * __y_mp1m * __factor1 - (__ell + __m - 1) * __y_mm * __factor2) / _Tp(__ell - __m);
                  __y_mm = __y_mp1m;
                  __y_mp1m = __y_lm;
                }

              return __y_lm;
            }
        }
    }
