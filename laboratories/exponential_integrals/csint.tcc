#include <limits>

    ///
    ///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
    ///         integrals by continued fraction for positive argument.
    ///
    template <typename Tp>
    void
    csint_cont_frac(const Tp ax, Tp & ci, Tp & si)
    {
      const int max_iter = 100;
      const Tp fp_min = std::numeric_limits<Tp>::min();
      const Tp eps = Tp(5) * std::numeric_limits<Tp>::epsilon();

      const Tp PIO2 = Tp(M_PI / 2.0);

      //    Evaluate Ci and Si by Lentz's modified method of continued fracions.
      std::complex<Tp> b(Tp(1), ax);
      std::complex<Tp> c(Tp(1) / fp_min, Tp(0));
      std::complex<Tp> h = Tp(1) / b;
      std::complex<Tp> d = h;
      int i = 0;
      for (i = 2; i <= max_iter; ++i)
        {
          const Tp a = -Tp((i - 1) * (i - 1));
          b += Tp(2);
          d = Tp(1) / (a * d + b);
          c = b + a / c;
          const std::complex<Tp> del = c * d;
          h *= del;
          if (std::abs(del.real() - Tp(1))
            + std::abs(del.imag()) < eps)
            break;
        }
      if (i > max_iter)
        throw std::runtime_error( "Continued fraction evaluation failed in csint." );

      h *= std::polar(Tp(1), -ax);
      ci = -h.real();
      si = PIO2 + h.imag();

      return;
    }


    ///
    ///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
    ///         integrals by series summation for positive argument.
    ///
    template <typename Tp>
    void
    csint_series(const Tp ax, Tp & ci, Tp & si)
    {
      const int max_iter = 100;
      const Tp fp_min = std::numeric_limits<Tp>::min();
      const Tp eps = Tp(5) * std::numeric_limits<Tp>::epsilon();

      const Tp GAMMA = Tp(0.5772156649015329L);

      //  Evaluate Ci and Si by series simultaneously.
      Tp sums = Tp(0);
      Tp sumc = Tp(0);

      if (ax < std::sqrt(fp_min))
        {
          //  Avoid underflow.
          sumc = Tp(0);
          sums = ax;
        }
      else
        {
          //  Evaluate Si and Ci by series expansion.
          Tp sum = Tp(0);
          Tp sign = Tp(1);
          Tp fact = Tp(1);
          bool odd = true;
          int k = 0;
          for ( k = 1; k <= max_iter; ++k )
            {
              fact *= ax / k;
              const Tp term = fact / k;
              sum += sign * term;
              const Tp err = term / std::abs(sum);

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
    template <typename Tp>
    std::pair<Tp, Tp>
    csint(Tp x)
    {

      const Tp x_min = Tp(2);

      Tp ci = Tp(0);
      Tp si = Tp(0);

      const Tp ax = std::abs(x);
      if (ax == Tp(0))
        {
          ci = -std::numeric_limits<Tp>::quiet_NaN();
          si = Tp(0);
          return std::make_pair(ci, si);
        }
      if (ax > x_min)
        csint_cont_frac(ax, ci, si);
      else
        csint_series(ax, ci, si);

      if (x < Tp(0))
        si = -si;

      return std::make_pair(ci, si);
    }


