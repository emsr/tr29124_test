
#ifndef _E_FLOAT_MPFR_TYPES_2009_11_15_H_
  #define _E_FLOAT_MPFR_TYPES_2009_11_15_H_

  // Wrap the types of MPFR.
  extern "C"
  {
    typedef unsigned long     mpfr_prec_t;
    typedef int               mpfr_sign_t;
    typedef long int          mp_exp_t;
    typedef unsigned long int mp_limb_t;

    #define mp_prec_t mpfr_prec_t
    #define mp_rnd_t  mpfr_rnd_t

    typedef struct
    {
      mpfr_prec_t _mpfr_prec;
      mpfr_sign_t _mpfr_sign;
      mp_exp_t    _mpfr_exp;
      mp_limb_t*  _mpfr_d;
    }
    __mpfr_struct;

    #define __gmp_const const

    typedef __mpfr_struct mpfr_t[1];
    typedef __mpfr_struct *mpfr_ptr;
    typedef __gmp_const __mpfr_struct *mpfr_srcptr;

    typedef enum
    {
      GMP_RNDN    =  0,
      GMP_RNDZ    =  1,
      GMP_RNDU    =  2,
      GMP_RNDD    =  3,
      GMP_RND_MAX =  4,
      GMP_RNDNA   = -1
    }
    mpfr_rnd_t;
  }

#endif // _E_FLOAT_MPFR_TYPES_2009_11_15_H_
