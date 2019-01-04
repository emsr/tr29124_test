
#ifndef _E_FLOAT_MPFR_PROTOS_2009_11_15_H_
  #define _E_FLOAT_MPFR_PROTOS_2009_11_15_H_

  // Wrap the function names of MPFR.
  extern "C"
  {
    #define MPFR_SIGN(x) ((x)->_mpfr_sign)

    void mpfr_init(mpfr_ptr);
    void mpfr_set_default_prec(mp_prec_t);
    void mpfr_clear(mpfr_ptr);

    int mpfr_set4   (mpfr_ptr, mpfr_srcptr, mp_rnd_t, int);
    #define mpfr_set(x, y, rnd) mpfr_set4((x), (y), (rnd), MPFR_SIGN(y))
    int mpfr_set_ui (mpfr_ptr, unsigned long, mp_rnd_t);
    int mpfr_set_d  (mpfr_ptr, double, mp_rnd_t);
    int mpfr_set_ld (mpfr_ptr, long double, mp_rnd_t);

    #define mpfr_init_set(x, y, rnd)     ( mpfr_init(x), static_cast<void>(mpfr_set   ((x), (y),  (rnd))) )
    #define mpfr_init_set_ui(x, i, rnd)  ( mpfr_init(x), static_cast<void>(mpfr_set_ui((x), (i),  (rnd))) )
    int mpfr_init_set_str(mpfr_ptr, const char*, int, mp_rnd_t);

    int mpfr_add(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_sub(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_mul(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_div(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t);

    int mpfr_add_ui(mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t);
    int mpfr_sub_ui(mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t);
    int mpfr_mul_ui(mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t);
    int mpfr_div_ui(mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t);

    int mpfr_cmp3(mpfr_srcptr, mpfr_srcptr, int);
    #define mpfr_cmp(b, c) mpfr_cmp3((b), (c), 1)

    int mpfr_neg(mpfr_ptr, mpfr_srcptr, mp_rnd_t);

    int mpfr_nan_p    (mpfr_srcptr);
    int mpfr_inf_p    (mpfr_srcptr);
    int mpfr_number_p (mpfr_srcptr);
    int mpfr_zero_p   (mpfr_srcptr);
    int mpfr_integer_p(mpfr_srcptr);
    int mpfr_sgn      (mpfr_srcptr);

    int mpfr_sprintf(char*, const char*, ...);

    double        mpfr_get_d     (mpfr_srcptr, mp_rnd_t);
    unsigned long mpfr_get_si    (mpfr_srcptr, mp_rnd_t);
    double        mpfr_get_d_2exp(signed long int*, mpfr_srcptr, mp_rnd_t);
    char*         mpfr_get_str   (char*, mp_exp_t*, int, size_t, mpfr_srcptr, mp_rnd_t);

    int mpfr_floor(mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_sqrt (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_cbrt (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_root (mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t);
    int mpfr_sin  (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_cos  (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_tan  (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_asin (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_acos (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_atan (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_sinh (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_cosh (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_tanh (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_asinh(mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_acosh(mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_atanh(mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_exp  (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_log  (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_gamma(mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_zeta (mpfr_ptr, mpfr_srcptr, mp_rnd_t);
    int mpfr_jn   (mpfr_ptr, long, mpfr_srcptr, mp_rnd_t);
    int mpfr_yn   (mpfr_ptr, long, mpfr_srcptr, mp_rnd_t);
  }

#endif // _E_FLOAT_MPFR_PROTOS_2009_11_15_H_
