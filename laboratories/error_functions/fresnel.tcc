
#include <limits>
#include <cmath>

    ///
    ///  @brief This routine returns the Fresnel cosine and sine integrals
    ///         as a pair by series expansion for positive argument.
    ///
    template <typename _Tp>
    void
    fresnel_series(const _Tp ax, _Tp & c, _Tp & s)
    {
      const int max_iter = 100;
      const _Tp eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();
      const _Tp PIO2 = _Tp(M_PI / 2.0);

      //  Evaluate S and C by series expansion.
      _Tp sum = _Tp(0);
      _Tp sums = _Tp(0);
      _Tp sumc = ax;
      _Tp sign = _Tp(1);
      _Tp fact = PIO2 * ax * ax;
      bool odd = true;
      _Tp term = ax;
      int n = 3;
      int k = 0;
      for (k = 1; k <= max_iter; ++k)
        {
          term *= fact / k;
          sum += sign * term / n;
          _Tp test = std::abs(sum) * eps;
          if ( odd )
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

          if (term < test)
            break;

          odd = ! odd;

          n += 2;
      }
      if (k > max_iter)
        throw std::runtime_error("Series evaluation failed"
                                 " in fresnel_series.");

      c = sumc;
      s = sums;

      return;
    }


    ///
    ///  @brief This routine computes the Fresnel cosine and sine integrals
    ///         by continued fractions for positive argument.
    ///
    template <typename _Tp>
    void
    fresnel_cont_frac(const _Tp ax, _Tp & c, _Tp & s)
    {
      const int max_iter = 100;
      const _Tp eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();
      const _Tp fp_min = std::numeric_limits<_Tp>::min();
      const _Tp PI = _Tp(M_PI);

      //  Evaluate S and C by Lentz's complex continued fraction method.
      const _Tp pix2 = PI * ax * ax;
      std::complex<_Tp> b(_Tp(1), -pix2);
      std::complex<_Tp> cc(_Tp(1) / fp_min, _Tp(0));
      std::complex<_Tp> h = _Tp(1) / b;
      std::complex<_Tp> d = h;
      int n = -1;
      int k = 0;
      for (k = 2; k <= max_iter; ++k)
        {
          n += 2;
          const _Tp a = -_Tp(n * (n + 1));
          b += _Tp(4);
          d = _Tp(1) / (a * d + b);
          cc = b + a / cc;
          const std::complex<_Tp> del = cc * d;
          h *= del;
          if (std::abs(del.real() - _Tp(1))
            + std::abs(del.imag()) < eps)
            break;
        }
      if (k > max_iter)
        throw std::runtime_error("Continued fraction evaluation"
                                 " failed in fresnel_cont_frac.");

      h = std::complex<_Tp>(ax, -ax) * h;
      std::complex<_Tp> phase = std::polar(_Tp(1), pix2/_Tp(2));
      std::complex<_Tp> cs = std::complex<_Tp>(_Tp(0.5L), _Tp(0.5L))
                             * (_Tp(1) - phase * h);
      c = cs.real();
      s = cs.imag();

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
    fresnel(const _Tp x)
    {

      const _Tp fp_min = std::numeric_limits<_Tp>::min();
      const _Tp x_min = _Tp(1.5L);

      _Tp c = _Tp(0);
      _Tp s = _Tp(0);

      const _Tp ax = std::abs(x);
      if (ax < std::sqrt(fp_min))
        {
          c = ax;
          s = _Tp(0);
        }
      else if (ax < x_min)
        fresnel_series(ax, c, s);
      else
        fresnel_cont_frac(ax, c, s);

      if (x < _Tp(0))
        {
          c = -c;
          s = -s;
        }

      return std::make_pair(c, s);
    }

