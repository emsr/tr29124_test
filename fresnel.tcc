#include <limits>

    ///
    ///  @brief This routine returns the Fresnel cosine and sine integrals
    ///         as a pair by series expansion for positive argument.
    ///
    template <typename _Tp>
    void
    __fresnel_series(const _Tp __ax, _Tp & __c, _Tp & __s)
    {
      const int __max_iter = 100;
      const _Tp __eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();
      const _Tp __PIO2 = _Tp(M_PI / 2.0);

      //  Evaluate S and C by series expansion.
      _Tp __sum = _Tp(0);
      _Tp __sums = _Tp(0);
      _Tp __sumc = __ax;
      _Tp __sign = _Tp(1);
      _Tp __fact = __PIO2 * __ax * __ax;
      bool __odd = true;
      _Tp __term = __ax;
      int __n = 3;
      int __k = 0;
      for (__k = 1; __k <= __max_iter; ++__k)
        {
          __term *= __fact / __k;
          __sum += __sign * __term / __n;
          _Tp __test = std::abs(__sum) * __eps;
          if ( __odd )
            {
              __sign = -__sign;
              __sums = __sum;
              __sum = __sumc;
            }
          else
            {
              __sumc = __sum;
              __sum = __sums;
            }

          if (__term < __test)
            break;

          __odd = ! __odd;

          __n += 2;
      }
      if (__k > __max_iter)
        throw std::runtime_error("Series evaluation failed"
                                 " in __fresnel_series.");

      __c = __sumc;
      __s = __sums;

      return;
    }


    ///
    ///  @brief This routine computes the Fresnel cosine and sine integrals
    ///         by continued fractions for positive argument.
    ///
    template <typename _Tp>
    void
    __fresnel_cont_frac(const _Tp __ax, _Tp & __c, _Tp & __s)
    {
      const int __max_iter = 100;
      const _Tp __eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();
      const _Tp __fp_min = std::numeric_limits<_Tp>::min();
      const _Tp __PI = _Tp(M_PI);

      //  Evaluate S and C by Lentz's complex continued fraction method.
      const _Tp __pix2 = __PI * __ax * __ax;
      std::complex<_Tp> __b(_Tp(1), -__pix2);
      std::complex<_Tp> __cc(_Tp(1) / __fp_min, _Tp(0));
      std::complex<_Tp> __h = _Tp(1) / __b;
      std::complex<_Tp> __d = __h;
      int __n = -1;
      int __k = 0;
      for (__k = 2; __k <= __max_iter; ++__k)
        {
          __n += 2;
          const _Tp __a = -_Tp(__n * (__n + 1));
          __b += _Tp(4);
          __d = _Tp(1) / (__a * __d + __b);
          __cc = __b + __a / __cc;
          const std::complex<_Tp> __del = __cc * __d;
          __h *= __del;
          if (std::abs(__del.real() - _Tp(1))
            + std::abs(__del.imag()) < __eps)
            break;
        }
      if (__k > __max_iter)
        throw std::runtime_error("Continued fraction evaluation"
                                 " failed in __fresnel_cont_frac.");

      __h = std::complex<_Tp>(__ax, -__ax) * __h;
      std::complex<_Tp> __phase = std::polar(_Tp(1), __pix2/_Tp(2));
      std::complex<_Tp> __cs = std::complex<_Tp>(_Tp(0.5L), _Tp(0.5L))
                             * (_Tp(1) - __phase * __h);
      __c = __cs.real();
      __s = __cs.imag();

      return;
    }


    ///
    ///  @brief This routine returns the Fresnel cosine and sine integrals
    ///         as a pair.
    ///
    ///  The Fresnel cosine integral is defined by:
    ///  @f[
    ///      C(x) = \int_0^x \cos(\frac{\pi}{2}t^2) dt
    ///  @f]
    ///
    ///  The Fresnel sine integral is defined by:
    ///  @f[
    ///      S(x) = \int_0^x \sin(\frac{\pi}{2}t^2) dt
    ///  @f]
    ///
    template <typename _Tp>
    std::pair<_Tp, _Tp>
    __fresnel(const _Tp __x)
    {

      const _Tp __fp_min = std::numeric_limits<_Tp>::min();
      const _Tp __x_min = _Tp(1.5L);

      _Tp __c = _Tp(0);
      _Tp __s = _Tp(0);

      const _Tp __ax = std::abs(__x);
      if (__ax < std::sqrt(__fp_min))
        {
          __c = __ax;
          __s = _Tp(0);
        }
      else if (__ax < __x_min)
        __fresnel_series(__ax, __c, __s);
      else
        __fresnel_cont_frac(__ax, __c, __s);

      if (__x < _Tp(0))
        {
          __c = -__c;
          __s = -__s;
        }

      return std::make_pair(__c, __s);
    }

