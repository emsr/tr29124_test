//
//  Copyright 1998-2020
//  Alion Science and Technology
//  US Govt Retains rights in accordance
//  with DoD FAR Supp 252.227 - 7013.
//

#include <cmath>
#include <complex>

enum class FockAiryType
{
  Value = 1,
  Deriv = 2
};

template<typename Tp>
  std::complex<Tp>
  airey_fock_series(FockAiryType lm, const std::complex<Tp>& e);

/**
 * Computes Fock-Airy function (kk = 1) or its derivative (kk = 2).
 *
 * @param[in] lm A flag indicating whether t compute the function vale
 *               or the function derivative.
 * @param[in] t  The complex argument.
 * @param[out] w1 The Fock-Airy function w_1(t) if FockAiryType == Value
 *                or the derivative w'_1(t) if FockAiryType == Deriv.
 * @param[out] w2 The Fock-Airy function w_2(t) if FockAiryType == Value
 *                or the derivative w'_2(t) if FockAiryType == Deriv.
 */
template<typename Tp>
  void
  airy_fock(FockAiryType kk, std::complex<Tp> t,
            std::complex<Tp>& w1, std::complex<Tp>& w2)
  {
    using cmplx = std::complex<Tp>;
    constexpr Tp dpi = 3.141592653589793238462643383279502884195L;
    constexpr Tp dpid4 = dpi / Tp{4};
    constexpr Tp phase_lim = Tp{5} * dpi / Tp{9};
    constexpr cmplx j{Tp{0}, Tp{1}};
    constexpr Tp s_small = 1.0e-6;
    constexpr Tp sqrt_pi = 1.772453850905516027298167483341145182798L;
    constexpr Tp sqrt_3 = 1.732050807568877293527446341505872366945L;
    constexpr Tp Ai_value_0 = 3.550280538878172392600631860041831763980e-1L;
    constexpr Tp Ai_deriv_0 = 2.588194037928067984051835601892039634793e-1L;

    if (std::abs(t) < 5.5)
      {
        constexpr Tp eps = Tp{5} * std::numeric_limits<Tp>::epsilon();

        auto dchk1 = Tp{0};
        auto dchk2 = Tp{0};

        const auto t2 = t * t;
        const auto t3 = t * t2;

        Tp d1, d2, da1, da2;
        cmplx dt1, dt2;
        cmplx dy1, dy2;
        if (kk == FockAiryType::Value)
          {
            d1 = Tp{0};
            d2 = Tp{1};
            da1 = Tp{-2};
            da2 = Tp{-1};
            dt1 = Tp{1};
            dt2 = t;
            dy1 = Tp{1};
            dy2 = t;
          }
        else
          {
            d1 = Tp{2};
            d2 = Tp{0};
            da1 = Tp{1};
            da2 = Tp{-1};
            dt1 = t2 / Tp{2};
            dt2 = Tp{1};
            dy1 = dt1;
            dy2 = Tp{1};
          }

        do
          {
            da1 += Tp{3};
            da2 += Tp{3};
            d1 += Tp{3};
            d2 += Tp{3};
            const auto dq1 = da1 / (d1 * (d1 - Tp{1}) * (d1 - Tp{2}));
            const auto dq2 = da2 / (d2 * (d2 - Tp{1}) * (d2 - Tp{2}));
            dt1 *= dq1 * t3;
            dt2 *= dq2 * t3;
            dy1 += dt1;
            dy2 += dt2;
            const auto dtst1 = std::abs(dy1);
            const auto dtst2 = std::abs(dy2);

            bool ok1 = false;
            auto tmp = std::abs(dtst1 - dchk1);
            if (tmp < eps)
              ok1 = true;
            else
              {
                tmp = std::abs(tmp / dtst1);
                if (tmp < eps)
                  ok1 = true;
              }

            bool ok2 = false;
            tmp = std::abs(dtst2 - dchk2);
            if (tmp < eps)
              ok2 = true;
            else
              {
                tmp = std::abs(tmp / dtst2);
                if (tmp < eps)
                  ok2 = true;
              }

            if (!ok1 || !ok2)
              {
                dchk1 = dtst1;
                dchk2 = dtst2;
                continue;
              }
            else
              break;
          }
        while (true);

        const auto du = sqrt_3 * sqrt_pi * (Ai_value_0 * dy1 + Ai_deriv_0 * dy2);
        const auto dv = sqrt_pi * (Ai_value_0 * dy1 - Ai_deriv_0 * dy2);
        w1 = du - j * dv;
        w2 = du + j * dv;
      }
    else
      {
        auto fa = Tp{1};
        auto fb = Tp{1};
        const auto phase = std::arg(t);
        const auto z = Tp{2} * std::pow(t, Tp{3}/Tp{2}) / Tp{3};
        const auto zm = std::abs(z);
        const auto pha = std::arg(z);
        auto e = z;
        auto aa = e;
        auto ab = e.imag();
        auto st = t;
        if (std::abs(phase) > phase_lim)
          {
            e = j * e * std::copysign(Tp{1}, phase);
            aa = -j * e;
            ab = -e.real() - (dpid4 * (Tp{3} - (Tp{2} * Tp(kk))));
            st = -t;
          }

        const auto z4 = std::pow(st, 0.25);
        const auto ee = std::exp(aa.real());
        const auto s = std::polar(ee, ab);

        const auto sum1 = airey_fock_series(kk, aa);
        const auto sum2 = airey_fock_series(kk, -aa);

        if (std::abs(phase) > phase_lim)
          {
            w1 = sum1 * s;
            w2 = sum2 / s;
            if (kk == FockAiryType::Value)
              {
                w1 /= z4;
                w2 /= z4;
              }
            else
              {
                w1 *= z4;
                w2 *= z4;
              }
          }
        else
          {
            auto a = s * fa * sum1;

            if (std::abs(phase) < dpi / Tp{9})
              {
                const auto b = 0.5 * fb * j * sum2 / s;

                if (kk == FockAiryType::Value)
                  {
                    w1 = (a - b) / z4;
                    w2 = (a + b) / z4;
                  }
                else
                  {
                    w1 = (a + b) * z4;
                    w2 = (a - b) * z4;
                  }
              }
            else
              {
                const auto b = fb * sum2 / s;

                if (t.imag() >= Tp{0})
                  {
                    if (kk == FockAiryType::Value)
                      {
                        w2 = (a + j * b) / z4;
                        if (std::abs(fa - Tp{1}) >= s_small)
                          a = s * sum1;
                        w1 = a / z4;
                      }
                    else
                      {
                        w2 = (a - j * b) * z4;
                        if (std::abs(fa - Tp{1}) >= s_small)
                          a = s * sum1;
                        w1 = a * z4;
                      }
                  }
                else
                  {
                    if (kk == FockAiryType::Value)
                      {
                        w1 = (a - j * b) / z4;
                        if (std::abs(fa - Tp{1}) >= s_small)
                          a = s * sum1;
                        w2 = a / z4;
                      }
                    else
                      {
                        w1 = (a + j * b) * z4;
                        if (std::abs(fa - Tp{1}) >= s_small)
                          a = s * sum1;
                        w2 = a * z4;
                      }
                  }
              }
          }
      }

    return;
  }

/**
 * Sum a series used in the Fock-Airy functions.
 *
 * @param[in]  lm  Flag for which series:
 *   lm == 1: Compute value series.
 *   lm == 2: Compute derivative series.
 * @param[in]  e  Argument.
 * @return  Output sum
 */
template<typename Tp>
  std::complex<Tp>
  airey_fock_series(FockAiryType lm, const std::complex<Tp>& e)
  {
    auto s = Tp{1} / e;

    std::complex<Tp> t;
    Tp ss;
    if (lm == FockAiryType::Value)
      {
        t = Tp{5} * s / Tp{72};
        ss = Tp{1};
      }
    else
      {
        t = Tp{7} * s / Tp{72};
        ss = Tp{-1};
      }
    auto a = Tp{1};
    auto f = Tp{1};
    auto sum = Tp{1} + ss * t;
    auto test = std::abs(sum);
    auto term2 = Tp{1000};
    auto term1 = Tp{10000};

    do
      {
        const auto check = test;
        term1 = term2;
        f += Tp{1};
        a += Tp{2};
        const auto b = Tp{6} * f;
        const auto tt = (b - Tp{6} + ss) * (b - Tp{3}) / (Tp{216} * f * a) * (b - ss);
        t *= tt * s;
        term2 = std::abs(t);
        sum += ss * t;
        test = std::abs(sum);
        if (std::abs(test - check) / test <= 5.0e-2)
          break;
      }
    while (term1 >= term2);

    return sum;
  }
