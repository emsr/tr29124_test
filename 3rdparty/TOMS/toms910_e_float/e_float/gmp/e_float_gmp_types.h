
#ifndef _E_FLOAT_GMP_TYPES_2009_05_02_H_
  #define _E_FLOAT_GMP_TYPES_2009_05_02_H_

  // Wrap the types of GMP.
  extern "C"
  {
    typedef long int          mp_size_t;
    typedef long int          mp_exp_t;
    typedef unsigned long int mp_limb_t;

    typedef struct struct__mpf_struct
    {
      int _mp_prec;
      int _mp_size;
      mp_exp_t _mp_exp;
      mp_limb_t *_mp_d;
    }
    __mpf_struct;

    typedef       __mpf_struct mpf_t[1];
    typedef       __mpf_struct* mpf_ptr;
    typedef const __mpf_struct* mpf_srcptr;
  }

#endif // _E_FLOAT_GMP_TYPES_2009_05_02_H_
