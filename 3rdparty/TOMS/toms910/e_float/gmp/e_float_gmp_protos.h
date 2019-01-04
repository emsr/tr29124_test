
#ifndef _E_FLOAT_GMP_PROTOS_2009_05_02_H_
  #define _E_FLOAT_GMP_PROTOS_2009_05_02_H_

  // Wrap the function names of GMP.
  extern "C"
  {
    void   __gmpf_set             (mpf_ptr, mpf_srcptr);
    void   __gmpf_set_default_prec(unsigned long int);
    void   __gmpf_set_prec        (mpf_ptr, unsigned long int);
    void   __gmpf_set_prec_raw    (mpf_ptr, unsigned long int);
    void   __gmpf_init            (mpf_ptr);
    void   __gmpf_init_set        (mpf_ptr, mpf_srcptr);
    int    __gmpf_init_set_str    (mpf_ptr, const char*, int);
    void   __gmpf_init_set_ui     (mpf_ptr, unsigned long int);
    void   __gmpf_clear           (mpf_ptr);
    void   __gmpf_neg             (mpf_ptr, mpf_srcptr);
    void   __gmpf_add_ui          (mpf_ptr, mpf_srcptr, unsigned long int);
    void   __gmpf_sub_ui          (mpf_ptr, mpf_srcptr, unsigned long int);
    void   __gmpf_mul_ui          (mpf_ptr, mpf_srcptr, unsigned long int);
    void   __gmpf_div_ui          (mpf_ptr, mpf_srcptr, unsigned long int);
    void   __gmpf_ui_div          (mpf_ptr, unsigned long int, mpf_srcptr);
    void   __gmpf_add             (mpf_ptr, mpf_srcptr, mpf_srcptr);
    void   __gmpf_mul             (mpf_ptr, mpf_srcptr, mpf_srcptr);
    void   __gmpf_div             (mpf_ptr, mpf_srcptr, mpf_srcptr);
    void   __gmpf_floor           (mpf_ptr, mpf_srcptr);
    void   __gmpf_sqrt            (mpf_ptr, mpf_srcptr);
    int    __gmpf_cmp             (mpf_srcptr, mpf_srcptr);
    int    __gmpf_integer_p       (mpf_srcptr);
    double __gmpf_get_d           (mpf_srcptr);
    double __gmpf_get_d_2exp      (signed long int*, mpf_srcptr);
    char*  __gmpf_get_str         (char*, mp_exp_t*, int, size_t, mpf_srcptr);
    int    __gmp_sprintf          (char*, const char*, ...);
  }

  inline void   mpf_set             (mpf_ptr dst, mpf_srcptr src)                           {        ::__gmpf_set                (dst, src); }
  inline void   mpf_set_default_prec(unsigned long int p)                                   {        ::__gmpf_set_default_prec   (p); }
  inline void   mpf_set_prec_raw    (mpf_ptr dst, unsigned long int p)                      {        ::__gmpf_set_prec_raw       (dst, p); }
  inline void   mpf_init            (mpf_ptr dst)                                           {        ::__gmpf_init               (dst); }
  inline void   mpf_init_set        (mpf_ptr dst, mpf_srcptr src)                           {        ::__gmpf_init_set           (dst, src); }
  inline int    mpf_init_set_str    (mpf_ptr dst, const char* str, int n)                   { return ::__gmpf_init_set_str       (dst, str, n); }
  inline void   mpf_init_set_ui     (mpf_ptr dst, unsigned long int u)                      {        ::__gmpf_init_set_ui        (dst, u); }
  inline void   mpf_clear           (mpf_ptr dst)                                           {        ::__gmpf_clear              (dst); }
  inline void   mpf_neg             (mpf_ptr dst, mpf_srcptr src)                           {        ::__gmpf_neg                (dst, src); }
  inline void   mpf_add_ui          (mpf_ptr dst, mpf_srcptr src, unsigned long int u)      {        ::__gmpf_add_ui             (dst, src, u); }
  inline void   mpf_sub_ui          (mpf_ptr dst, mpf_srcptr src, unsigned long int u)      {        ::__gmpf_sub_ui             (dst, src, u); }
  inline void   mpf_mul_ui          (mpf_ptr dst, mpf_srcptr src, unsigned long int u)      {        ::__gmpf_mul_ui             (dst, src, u); }
  inline void   mpf_div_ui          (mpf_ptr dst, mpf_srcptr src, unsigned long int u)      {        ::__gmpf_div_ui             (dst, src, u); }
  inline void   mpf_ui_div          (mpf_ptr dst, unsigned long int u, mpf_srcptr src)      {        ::__gmpf_ui_div             (dst, u, src); }
  inline void   mpf_add             (mpf_ptr dst, mpf_srcptr u, mpf_srcptr v)               {        ::__gmpf_add                (dst, u, v); }
  inline void   mpf_mul             (mpf_ptr dst, mpf_srcptr u, mpf_srcptr v)               {        ::__gmpf_mul                (dst, u, v); }
  inline void   mpf_floor           (mpf_ptr dst, mpf_srcptr x)                             {        ::__gmpf_floor              (dst, x); }
  inline void   mpf_sqrt            (mpf_ptr dst, mpf_srcptr x)                             {        ::__gmpf_sqrt               (dst, x); }
  inline int    mpf_cmp             (mpf_srcptr x, mpf_srcptr y)                            { return ::__gmpf_cmp                (x, y); }
  inline int    mpf_sgn             (mpf_srcptr x)                                          { return x->_mp_size < 0 ? -1 : (x->_mp_size > 0); }
  inline int    mpf_integer_p       (mpf_srcptr src)                                        { return ::__gmpf_integer_p          (src); }
  inline double mpf_get_d           (mpf_srcptr x)                                          { return ::__gmpf_get_d              (x); }
  inline double mpf_get_d_2exp      (signed long int* e, mpf_srcptr x)                      { return ::__gmpf_get_d_2exp         (e, x); }
  inline char*  mpf_get_str         (char* str, mp_exp_t* e, int i, size_t s, mpf_srcptr x) { return ::__gmpf_get_str            (str, e, i, s, x); }
  #define       gmp_sprintf         __gmp_sprintf

#endif // _E_FLOAT_GMP_PROTOS_2009_05_02_H_
