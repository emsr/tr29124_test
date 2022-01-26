#include <limits>

    ///
    ///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
    ///         integrals by continued fraction for positive argument.
    ///
    template <typename _Tp>
    void
    csint_cont_frac(const _Tp ax, _Tp & ci, _Tp & si)
    {
      const int max_iter = 100;
      const _Tp fp_min = std::numeric_limits<_Tp>::min();
      const _Tp eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();

      const _Tp PIO2 = _Tp(M_PI / 2.0);

      //    Evaluate Ci and Si by Lentz's modified method of continued fracions.
      std::complex<_Tp> b(_Tp(1), ax);
      std::complex<_Tp> c(_Tp(1) / fp_min, _Tp(0));
      std::complex<_Tp> h = _Tp(1) / b;
      std::complex<_Tp> d = h;
      int i = 0;
      for (i = 2; i <= max_iter; ++i)
        {
          const _Tp a = -_Tp((i - 1) * (i - 1));
          b += _Tp(2);
          d = _Tp(1) / (a * d + b);
          c = b + a / c;
          const std::complex<_Tp> del = c * d;
          h *= del;
          if (std::abs(del.real() - _Tp(1))
            + std::abs(del.imag()) < eps)
            break;
        }
      if (i > max_iter)
        throw std::runtime_error( "Continued fraction evaluation failed in csint." );

      h *= std::polar(_Tp(1), -ax);
      ci = -h.real();
      si = PIO2 + h.imag();

      return;
    }


    ///
    ///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
    ///         integrals by series summation for positive argument.
    ///
    template <typename _Tp>
    void
    csint_series(const _Tp ax, _Tp & ci, _Tp & si)
    {
      const int max_iter = 100;
      const _Tp fp_min = std::numeric_limits<_Tp>::min();
      const _Tp eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();

      const _Tp GAMMA = _Tp(0.5772156649015329L);

      //  Evaluate Ci and Si by series simultaneously.
      _Tp sums = _Tp(0);
      _Tp sumc = _Tp(0);

      if (ax < std::sqrt(fp_min))
        {
          //  Avoid underflow.
          sumc = _Tp(0);
          sums = ax;
        }
      else
        {
          //  Evaluate Si and Ci by series expansion.
          _Tp sum = _Tp(0);
          _Tp sign = _Tp(1);
          _Tp fact = _Tp(1);
          bool odd = true;
          int k = 0;
          for ( k = 1; k <= max_iter; ++k )
            {
              fact *= ax / k;
              const _Tp term = fact / k;
              sum += sign * term;
              const _Tp err = term / std::abs(sum);

              if (odd)
                {
                  sign = -sign;
                  sums = sum;
                  sum = sumc;
                }
              else
                {
                  sumc = sum;
                  sum = sums;
                }

              if (err < eps)
                break;

              odd = ! odd;
            }
          if (k > max_iter)
            throw std::runtime_error("Series evaluation failed in csint.");

          ci = sumc + std::log(ax) + GAMMA;
          si = sums;
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
    csint(_Tp x)
    {

      const _Tp x_min = _Tp(2);

      _Tp ci = _Tp(0);
      _Tp si = _Tp(0);

      const _Tp ax = std::abs(x);
      if (ax == _Tp(0))
        {
          ci = -std::numeric_limits<_Tp>::quiet_NaN();
          si = _Tp(0);
          return std::make_pair(ci, si);
        }
      if (ax > x_min)
        csint_cont_frac(ax, ci, si);
      else
        csint_series(ax, ci, si);

      if (x < _Tp(0))
        si = -si;

      return std::make_pair(ci, si);
    }


