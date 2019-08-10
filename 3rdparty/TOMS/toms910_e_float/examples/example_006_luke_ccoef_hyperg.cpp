#include <deque>
#include <numeric>
#include <algorithm>
#include <functional>

#include <e_float/e_float.h>
#include <examples/examples.h>
#include <functions/functions.h>
#include <utility/util_alternating_sum.h>
#include <utility/util_digit_scale.h>

namespace examples
{
  namespace nr_006
  {
    template<typename T> class HypergPFQ_Base
    {
    private:

      const HypergPFQ_Base& operator=(const HypergPFQ_Base&);

    protected:

      static const INT32 ZERO  = static_cast<INT32>(0);
      static const INT32 ONE   = static_cast<INT32>(1);
      static const INT32 TWO   = static_cast<INT32>(2);
      static const INT32 THREE = static_cast<INT32>(3);
      static const INT32 FOUR  = static_cast<INT32>(4);

      const   T       Z;
      const   e_float W;
              INT32   N;
      mutable std::deque<T> C;


    protected:

      HypergPFQ_Base(const T& z,
                     const e_float& w) : Z(z),
                                         W(w),
                                         N(static_cast<INT32>(Util::DigitScale() * static_cast<double>(400.0))) { }

    public:

      virtual ~HypergPFQ_Base() { }

      virtual void ccoef(void) const = 0;

      virtual T series(void) const
      {
        using ef::chebyshev_t;
        using efz::chebyshev_t;

        // Compute the Chebyshev coefficients.
        // Get the values of the shifted Chebyshev polynomials.
        std::vector<T> Tn_shifted;
        const T z_shifted = ((Z / W) * static_cast<INT32>(2)) - static_cast<INT32>(1);

        chebyshev_t(static_cast<UINT32>(C.size()), z_shifted, &Tn_shifted);

        // Luke: C     ---------- COMPUTE SCALE FACTOR                       ----------
        // Luke: C
        // Luke: C     ---------- SCALE THE COEFFICIENTS                     ----------
        // Luke: C

        // The coefficient scaling is preformed after the Chebyshev summation,
        // and it is carried out with a single division operation.
        const T scale = std::accumulate(C.begin(), C.end(), T(0), Util::alternating_sum<T>());

        // Compute the result of the series expansion using unscaled coefficients.
        const T sum = std::inner_product(C.begin(), C.end(), Tn_shifted.begin(), T(0));

        // Return the properly scaled result.
        return sum / scale;
      }
    };

    template<typename T> class Ccoef3Hyperg1F1 : public HypergPFQ_Base<T>
    {
    private:

      const T AP;
      const T CP;

    public:

      Ccoef3Hyperg1F1(const T& a, const T& c, const T& z, const e_float& w) : HypergPFQ_Base<T>(z, w),
                                                                              AP(a),
                                                                              CP(c) { }

      virtual ~Ccoef3Hyperg1F1() { }

    public:

      virtual void ccoef(void) const
      {
        // See Luke 1977 page 74.
        const INT32 N1 = static_cast<INT32>(HypergPFQ_Base<T>::N + static_cast<INT32>(1));
        const INT32 N2 = static_cast<INT32>(HypergPFQ_Base<T>::N + static_cast<INT32>(2));

        // Luke: C     ---------- START COMPUTING COEFFICIENTS USING         ----------
        // Luke: C     ---------- BACKWARD RECURRENCE SCHEME                 ----------
        // Luke: C
        T A3(HypergPFQ_Base<T>::ZERO);
        T A2(HypergPFQ_Base<T>::ZERO);
        T A1(HypergPFQ_Base<T>::ONE);

        HypergPFQ_Base<T>::C.resize(1u, A1);

        INT32 X  = N1;
        INT32 X1 = N2;

        T XA  =  X + AP;
        T X3A = (X + HypergPFQ_Base<T>::THREE) - AP;

        const T Z1 = HypergPFQ_Base<T>::FOUR / HypergPFQ_Base<T>::W;

        for(INT32 k = static_cast<INT32>(0); k < N1; k++)
        {
          --X;
          --X1;
          --XA;
          --X3A;

          const T X3A_over_X2 = X3A / static_cast<INT32>(X + HypergPFQ_Base<T>::TWO);

          // The terms have been slightly re-arranged resulting in fewer multiplications.
          // Parentheses have been added to avoid reliance on operator precedence.
          const T PART1 =  A1 * (((X + CP) * Z1) - X3A_over_X2);
          const T PART2 =  A2 * (Z1 * ((X + HypergPFQ_Base<T>::THREE) - CP) + (XA / X1));
          const T PART3 =  A3 * X3A_over_X2;

          const T term = (((PART1 + PART2) + PART3) * X1) / XA;

          HypergPFQ_Base<T>::C.push_front(term);

          A3 = A2;
          A2 = A1;
          A1 = HypergPFQ_Base<T>::C.front();
        }

        HypergPFQ_Base<T>::C.front() /= static_cast<INT32>(2);
      }
    };

    template<typename T> class Ccoef6Hyperg1F2 : public HypergPFQ_Base<T>
    {
    private:

      const T AP;
      const T BP;
      const T CP;

    public:

      Ccoef6Hyperg1F2(const T& a,
                      const T& b,
                      const T& c,
                      const T& z,
                      const e_float& w) : HypergPFQ_Base<T>(z, w),
                                          AP(a),
                                          BP(b),
                                          CP(c) { }

      virtual ~Ccoef6Hyperg1F2() { }

    public:

      virtual void ccoef(void) const
      {
        // See Luke 1977 page 85.
        const INT32 N1 = static_cast<INT32>(HypergPFQ_Base<T>::N + static_cast<INT32>(1));

        // Luke: C     ---------- START COMPUTING COEFFICIENTS USING         ----------
        // Luke: C     ---------- BACKWARD RECURRENCE SCHEME                 ----------
        // Luke: C
        T A4(HypergPFQ_Base<T>::ZERO);
        T A3(HypergPFQ_Base<T>::ZERO);
        T A2(HypergPFQ_Base<T>::ZERO);
        T A1(HypergPFQ_Base<T>::ONE);

        HypergPFQ_Base<T>::C.resize(1u, A1);

        INT32 X  = N1;
        T     PP = X + AP;

        const T Z1 = HypergPFQ_Base<T>::FOUR / HypergPFQ_Base<T>::W;

        for(INT32 k = static_cast<INT32>(0); k < N1; k++)
        {
          --X;
          --PP;

          const INT32 TWO_X    = static_cast<INT32>(X * HypergPFQ_Base<T>::TWO);
          const INT32 X_PLUS_1 = static_cast<INT32>(X + HypergPFQ_Base<T>::ONE);
          const INT32 X_PLUS_3 = static_cast<INT32>(X + HypergPFQ_Base<T>::THREE);
          const INT32 X_PLUS_4 = static_cast<INT32>(X + HypergPFQ_Base<T>::FOUR);

          const T QQ = T(TWO_X + HypergPFQ_Base<T>::THREE) / static_cast<INT32>(TWO_X + static_cast<INT32>(5));
          const T SS = (X + BP) * (X + CP);

          // The terms have been slightly re-arranged resulting in fewer multiplications.
          // Parentheses have been added to avoid reliance on operator precedence.
          const T PART1 =   A1 * (((PP - (QQ * (PP + HypergPFQ_Base<T>::ONE))) * HypergPFQ_Base<T>::TWO) + (SS * Z1));
          const T PART2 =  (A2 * (X + HypergPFQ_Base<T>::TWO)) * ((((TWO_X + HypergPFQ_Base<T>::ONE) * PP) / X_PLUS_1) - ((QQ * HypergPFQ_Base<T>::FOUR) * (PP + HypergPFQ_Base<T>::ONE)) + (((TWO_X + HypergPFQ_Base<T>::THREE) * (PP + HypergPFQ_Base<T>::TWO)) / X_PLUS_3) + ((Z1 * HypergPFQ_Base<T>::TWO) * (SS - (QQ * (X_PLUS_1 + BP)) * (X_PLUS_1 + CP))));
          const T PART3 =   A3 * ((((X_PLUS_3 - AP) - (QQ * (X_PLUS_4 - AP))) * HypergPFQ_Base<T>::TWO) + (((QQ * Z1) * (X_PLUS_4 - BP)) * (X_PLUS_4 - CP)));
          const T PART4 = ((A4 * QQ) * (X_PLUS_4 - AP)) / X_PLUS_3;

          const T term = (((PART1 - PART2) + (PART3 - PART4)) * X_PLUS_1) / PP;

          HypergPFQ_Base<T>::C.push_front(term);

          A4 = A3;
          A3 = A2;
          A2 = A1;
          A1 = HypergPFQ_Base<T>::C.front();
        }

        HypergPFQ_Base<T>::C.front() /= static_cast<INT32>(2);
      }
    };

    template<typename T> class Ccoef2Hyperg2F1 : public HypergPFQ_Base<T>
    {
    private:

      const T AP;
      const T BP;
      const T CP;

    public:

      Ccoef2Hyperg2F1(const T& a,
                      const T& b,
                      const T& c,
                      const T& z,
                      const e_float& w) : HypergPFQ_Base<T>(z, w),
                                          AP(a),
                                          BP(b),
                                          CP(c)
      {
        // Set anew the number of terms in the Chebyshev expansion.
        HypergPFQ_Base<T>::N = static_cast<INT32>(Util::DigitScale() * static_cast<double>(1400.0));
      }

      virtual ~Ccoef2Hyperg2F1() { }

    public:

      virtual void ccoef(void) const
      {
        using ef::inv;
        using efz::inv;

        // See Luke 1977 page 59.
        const INT32 N1 = static_cast<INT32>(HypergPFQ_Base<T>::N + static_cast<INT32>(1));
        const INT32 N2 = static_cast<INT32>(HypergPFQ_Base<T>::N + static_cast<INT32>(2));

        // Luke: C     ---------- START COMPUTING COEFFICIENTS USING         ----------
        // Luke: C     ---------- BACKWARD RECURRENCE SCHEME                 ----------
        // Luke: C
        T A3(HypergPFQ_Base<T>::ZERO);
        T A2(HypergPFQ_Base<T>::ZERO);
        T A1(HypergPFQ_Base<T>::ONE);

        HypergPFQ_Base<T>::C.resize(1u, A1);

        INT32 X  = N1;
        INT32 X1 = N2;
        INT32 X3  = static_cast<INT32>((X * HypergPFQ_Base<T>::TWO) + HypergPFQ_Base<T>::THREE);

        T X3A = (X + HypergPFQ_Base<T>::THREE) - AP;
        T X3B = (X + HypergPFQ_Base<T>::THREE) - BP;

        const T Z1 = HypergPFQ_Base<T>::FOUR / HypergPFQ_Base<T>::W;

        for(INT32 k = static_cast<INT32>(0); k < N1; k++)
        {
          --X;
          --X1;
          --X3A;
          --X3B;
          X3 -= HypergPFQ_Base<T>::TWO;

          const INT32 X_PLUS_2 = static_cast<INT32>(X + HypergPFQ_Base<T>::TWO);

          const T XAB = inv((X + AP) * (X + BP));

          // The terms have been slightly re-arranged resulting in fewer multiplications.
          // Parentheses have been added to avoid reliance on operator precedence.
          const T PART1 = (A1 * X1) * (HypergPFQ_Base<T>::TWO - (((AP + X1) * (BP + X1)) * ((T(X3) / X_PLUS_2) * XAB)) + ((CP + X) * (XAB * Z1)));
          const T PART2 = (A2 * XAB) * ((X3A * X3B) - (X3 * ((X3A + X3B) - HypergPFQ_Base<T>::ONE)) + (((HypergPFQ_Base<T>::THREE - CP) + X) * (X1 * Z1)));
          const T PART3 = (A3 * X1) * (X3A / X_PLUS_2) * (X3B * XAB);

          const T term = (PART1 + PART2) - PART3;

          HypergPFQ_Base<T>::C.push_front(term);

          A3 = A2;
          A2 = A1;
          A1 = HypergPFQ_Base<T>::C.front();
        }

        HypergPFQ_Base<T>::C.front() /= static_cast<INT32>(2);
      }
    };
  }
}

e_float examples::nr_006::luke_ccoef3_hyperg_1f1(const e_float& a, const e_float& b, const e_float& x)
{
  const Ccoef3Hyperg1F1<e_float> c3_h1f1(a, b, x, e_float(-10));

  c3_h1f1.ccoef();

  return c3_h1f1.series();
}

ef_complex examples::nr_006::luke_ccoef3_hyperg_1f1(const ef_complex& a, const ef_complex& b, const ef_complex& z)
{
  const Ccoef3Hyperg1F1<ef_complex> c3_h1f1(a, b, z, e_float(-10));

  c3_h1f1.ccoef();

  return c3_h1f1.series();
}

e_float examples::nr_006::luke_ccoef6_hyperg_1f2(const e_float& a, const e_float& b, const e_float& c, const e_float& x)
{
  const Ccoef6Hyperg1F2<e_float> c6_h1f2(a, b, c, x, e_float(10));

  c6_h1f2.ccoef();

  return c6_h1f2.series();
}

ef_complex examples::nr_006::luke_ccoef6_hyperg_1f2(const ef_complex& a, const ef_complex& b, const ef_complex& c, const ef_complex& z)
{
  const Ccoef6Hyperg1F2<ef_complex> c6_h1f2(a, b, c, z, e_float(10));

  c6_h1f2.ccoef();

  return c6_h1f2.series();
}

e_float examples::nr_006::luke_ccoef2_hyperg_2f1(const e_float& a, const e_float& b, const e_float& c, const e_float& x)
{
  const Ccoef2Hyperg2F1<e_float> c2_h2f1(a, b, c, x, e_float(-10));

  c2_h2f1.ccoef();

  return c2_h2f1.series();
}

ef_complex examples::nr_006::luke_ccoef2_hyperg_2f1(const ef_complex& a, const ef_complex& b, const ef_complex& c, const ef_complex& z)
{
  const Ccoef2Hyperg2F1<ef_complex> c2_h2f1(a, b, c, z, e_float(-10));

  c2_h2f1.ccoef();

  return c2_h2f1.series();
}

ef_complex examples::nr_006::luke_ccoef3_hyperg_1f1_test(void)
{
  // N[Hypergeometric1F1[(1/2)+(I/4), (2/3)+((Pi/10) I), (-5-(1/3))-(EulerGamma I)],100]
  static const ef_complex a(ef::half(), ef::quarter());
  static const ef_complex b(ef::two_third(), ef::pi() / 10);
  static const ef_complex z(-5 - ef::third(), -ef::euler_gamma());

  return luke_ccoef3_hyperg_1f1(a, b, z);
}

ef_complex examples::nr_006::luke_ccoef6_hyperg_1f2_test(void)
{
  // N[HypergeometricPFQ[{(1/2)+(I/4)}, {(2/3)+((Pi/10) I), (3/2)+(Catalan I)}, (5+(1/3))+(EulerGamma I)],100]
  static const ef_complex a(ef::half(), ef::quarter());
  static const ef_complex b(ef::two_third(), ef::pi() / 10);
  static const ef_complex c(ef::three_half(), ef::catalan());
  static const ef_complex z(5 + ef::third(), ef::euler_gamma());

  return luke_ccoef6_hyperg_1f2(a, b, c, z);
}

ef_complex examples::nr_006::luke_ccoef2_hyperg_2f1_test(void)
{
  // N[Hypergeometric2F1[(1/2)+(I/4), (2/3)+((Pi/10) I), (3/2)+(Catalan I), (-5-(1/3))-(EulerGamma I)],100]
  static const ef_complex a(ef::half(), ef::quarter());
  static const ef_complex b(ef::two_third(), ef::pi() / 10);
  static const ef_complex c(ef::three_half(), ef::catalan());
  static const ef_complex z(-5 - ef::third(), -ef::euler_gamma());

  return luke_ccoef2_hyperg_2f1(a, b, c, z);
}
