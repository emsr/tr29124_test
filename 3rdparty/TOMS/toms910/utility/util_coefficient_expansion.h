
#ifndef _UTIL_COEFFICIENT_EXPANSION_2009_11_23_H_
  #define _UTIL_COEFFICIENT_EXPANSION_2009_11_23_H_

  namespace Util
  {
    template<typename T1, typename T2 = T1> struct coefficient_expansion
    {
    private:

      const T1 x_expand;
      T1       x_expand_pow_k;

    private:

      coefficient_expansion& operator=(const coefficient_expansion&);

      static const T2& one_t2(void)
      {
        static const T2 val_t2(1u);
        return val_t2;
      }

    public:

      coefficient_expansion(const T1& expand, const T2& init = one_t2()) : x_expand      (expand),
                                                                           x_expand_pow_k(init) { }

      T1 operator()(const T1& sum, const T2& ck)
      {
        const T1 ck_x_pow_k(ck * x_expand_pow_k);

        x_expand_pow_k *= x_expand;

        return sum + ck_x_pow_k;
      }
    };
  }

#endif // _UTIL_COEFFICIENT_EXPANSION_2009_11_23_H_
