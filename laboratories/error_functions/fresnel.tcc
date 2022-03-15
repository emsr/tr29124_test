
#include <limits>
#include <cmath>

    ///
    ///  @brief This routine returns the Fresnel cosine and sine integrals
    ///         as a pair by series expansion for positive argument.
    ///
    template <typename Tp>
    void
    fresnel_series(const Tp ax, Tp & c, Tp & s)
    {
      const int max_iter = 100;
      const Tp eps = Tp(5) * std::numeric_limits<Tp>::epsilon();
      const Tp PIO2 = Tp(M_PI / 2.0);

      //  Evaluate S and C by series expansion.
      Tp sum = Tp(0);
      Tp sums = Tp(0);
      Tp sumc = ax;
      Tp sign = Tp(1);
      Tp fact = PIO2 * ax * ax;
      bool odd = true;
      Tp term = ax;
      int n = 3;
      int k = 0;
      for (k = 1; k <= max_iter; ++k)
        {
          term *= fact / k;
          sum += sign * term / n;
          Tp test = std::abs(sum) * eps;
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
    template <typename Tp>
    void
    fresnel_cont_frac(const Tp ax, Tp & c, Tp & s)
    {
      const int max_iter = 100;
      const Tp eps = Tp(5) * std::numeric_limits<Tp>::epsilon();
      const Tp fp_min = std::numeric_limits<Tp>::min();
      const Tp PI = Tp(M_PI);

      //  Evaluate S and C by Lentz's complex continued fraction method.
      const Tp pix2 = PI * ax * ax;
      std::complex<Tp> b(Tp(1), -pix2);
      std::complex<Tp> cc(Tp(1) / fp_min, Tp(0));
      std::complex<Tp> h = Tp(1) / b;
      std::complex<Tp> d = h;
      int n = -1;
      int k = 0;
      for (k = 2; k <= max_iter; ++k)
        {
          n += 2;
          const Tp a = -Tp(n * (n + 1));
          b += Tp(4);
          d = Tp(1) / (a * d + b);
          cc = b + a / cc;
          const std::complex<Tp> del = cc * d;
          h *= del;
          if (std::abs(del.real() - Tp(1))
            + std::abs(del.imag()) < eps)
            break;
        }
      if (k > max_iter)
        throw std::runtime_error("Continued fraction evaluation"
                                 " failed in fresnel_cont_frac.");

      h = std::complex<Tp>(ax, -ax) * h;
      std::complex<Tp> phase = std::polar(Tp(1), pix2/Tp(2));
      std::complex<Tp> cs = std::complex<Tp>(Tp(0.5L), Tp(0.5L))
                             * (Tp(1) - phase * h);
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
    template <typename Tp>
    std::pair<Tp, Tp>
    fresnel(const Tp x)
    {

      const Tp fp_min = std::numeric_limits<Tp>::min();
      const Tp x_min = Tp(1.5L);

      Tp c = Tp(0);
      Tp s = Tp(0);

      const Tp ax = std::abs(x);
      if (ax < std::sqrt(fp_min))
        {
          c = ax;
          s = Tp(0);
        }
      else if (ax < x_min)
        fresnel_series(ax, c, s);
      else
        fresnel_cont_frac(ax, c, s);

      if (x < Tp(0))
        {
          c = -c;
          s = -s;
        }

      return std::make_pair(c, s);
    }

