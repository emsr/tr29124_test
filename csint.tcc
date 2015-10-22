#include <limits>

    ///
    ///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
    ///         integrals by continued fraction for positive argument.
    ///
    template <typename _Tp>
    void
    __csint_cont_frac(const _Tp __ax, _Tp & __ci, _Tp & __si)
    {
      const int __max_iter = 100;
      const _Tp __fp_min = std::numeric_limits<_Tp>::min();
      const _Tp __eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();

      const _Tp __PIO2 = _Tp(M_PI / 2.0);

      //    Evaluate Ci and Si by Lentz's modified method of continued fracions.
      std::complex<_Tp> __b(_Tp(1), __ax);
      std::complex<_Tp> __c(_Tp(1) / __fp_min, _Tp(0));
      std::complex<_Tp> __h = _Tp(1) / __b;
      std::complex<_Tp> __d = __h;
      int __i = 0;
      for (__i = 2; __i <= __max_iter; ++__i)
        {
          const _Tp __a = -_Tp((__i - 1) * (__i - 1));
          __b += _Tp(2);
          __d = _Tp(1) / (__a * __d + __b);
          __c = __b + __a / __c;
          const std::complex<_Tp> __del = __c * __d;
          __h *= __del;
          if (std::abs(__del.real() - _Tp(1))
            + std::abs(__del.imag()) < __eps)
            break;
        }
      if (__i > __max_iter)
        throw std::runtime_error( "Continued fraction evaluation failed in __csint." );

      __h *= std::polar(_Tp(1), -__ax);
      __ci = -__h.real();
      __si = __PIO2 + __h.imag();

      return;
    }


    ///
    ///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
    ///         integrals by series summation for positive argument.
    ///
    template <typename _Tp>
    void
    __csint_series(const _Tp __ax, _Tp & __ci, _Tp & __si)
    {
      const int __max_iter = 100;
      const _Tp __fp_min = std::numeric_limits<_Tp>::min();
      const _Tp __eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();

      const _Tp __GAMMA = _Tp(0.5772156649015329L);

      //  Evaluate Ci and Si by series simultaneously.
      _Tp __sums = _Tp(0);
      _Tp __sumc = _Tp(0);

      if (__ax < std::sqrt(__fp_min))
        {
          //  Avoid underflow.
          __sumc = _Tp(0);
          __sums = __ax;
        }
      else
        {
          //  Evaluate Si and Ci by series expansion.
          _Tp __sum = _Tp(0);
          _Tp __sign = _Tp(1);
          _Tp __fact = _Tp(1);
          bool __odd = true;
          int __k = 0;
          for ( __k = 1; __k <= __max_iter; ++__k )
            {
              __fact *= __ax / __k;
              const _Tp __term = __fact / __k;
              __sum += __sign * __term;
              const _Tp __err = __term / std::abs(__sum);

              if (__odd)
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

              if (__err < __eps)
                break;

              __odd = ! __odd;
            }
          if (__k > __max_iter)
            throw std::runtime_error("Series evaluation failed in __csint.");

          __ci = __sumc + std::log(__ax) + __GAMMA;
          __si = __sums;
        }

      return;
    }


    ///
    ///  @brief This routine returns the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
    ///         integrals as a pair.
    ///
    ///  The cosine integral is defined by:
    ///  @f[
    ///      Ci(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cos(t) - 1}{t}
    ///  @f]
    ///
    ///  The sine integral is defined by:
    ///  @f[
    ///      Si(x) = \int_0^x dt \frac{\sin(t)}{t}
    ///  @f]
    ///
    template <typename _Tp>
    std::pair<_Tp, _Tp>
    __csint(_Tp __x)
    {

      const _Tp __x_min = _Tp(2);

      _Tp __ci = _Tp(0);
      _Tp __si = _Tp(0);

      const _Tp __ax = std::abs(__x);
      if (__ax == _Tp(0))
        {
          __ci = -std::numeric_limits<_Tp>::quiet_NaN();
          __si = _Tp(0);
          return std::make_pair(__ci, __si);
        }
      if (__ax > __x_min)
        __csint_cont_frac(__ax, __ci, __si);
      else
        __csint_series(__ax, __ci, __si);

      if (__x < _Tp(0))
        __si = -__si;

      return std::make_pair(__ci, __si);
    }


