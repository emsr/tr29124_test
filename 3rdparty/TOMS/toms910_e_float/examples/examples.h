#ifndef _EXAMPLES_2010_01_02_H_
  #define _EXAMPLES_2010_01_02_H_

  #include <e_float/e_float.h>
  #include <functions/complex/e_float_complex.h>

  namespace examples
  {
    namespace nr_001
    {
      void basic_usage_real(void);
    }

    namespace nr_002
    {
      void basic_usage_imag(void);
    }

    namespace nr_003
    {
      e_float jahnke_emden_lambda(const e_float& v, const e_float& x);
    }

    namespace nr_004
    {
      e_float bessel_jv_derivative_wrt_v     (const e_float& v, const e_float& x);
      e_float bessel_jv_derivative_wrt_v_test(void);
    }

    namespace nr_005
    {
      e_float recursive_trapezoid_j0     (const e_float& x);
      e_float recursive_trapezoid_j0_test(void);
    }

    namespace nr_006
    {
      e_float    luke_ccoef3_hyperg_1f1         (const e_float& a, const e_float& b, const e_float& x);
      e_float    luke_ccoef6_hyperg_1f2         (const e_float& a, const e_float& b, const e_float& c, const e_float& x);
      e_float    luke_ccoef2_hyperg_2f1         (const e_float& a, const e_float& b, const e_float& c, const e_float& x);
      ef_complex luke_ccoef3_hyperg_1f1         (const ef_complex& a, const ef_complex& b, const ef_complex& z);
      ef_complex luke_ccoef6_hyperg_1f2         (const ef_complex& a, const ef_complex& b, const ef_complex& c, const ef_complex& z);
      ef_complex luke_ccoef2_hyperg_2f1         (const ef_complex& a, const ef_complex& b, const ef_complex& c, const ef_complex& z);
      ef_complex luke_ccoef3_hyperg_1f1_test    (void);
      ef_complex luke_ccoef6_hyperg_1f2_test    (void);
      ef_complex luke_ccoef2_hyperg_2f1_test    (void);
    }

    namespace nr_007
    {
      e_float    poly_gamma_mathematica(const e_float& v,    const e_float& x);
      ef_complex poly_gamma_mathematica(const ef_complex& v, const ef_complex& z);
      ef_complex poly_gamma_mathematica_test(void);
    }
  }

#endif // _EXAMPLES_2010_01_02_H_
