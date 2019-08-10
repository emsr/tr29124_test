
#ifndef _UTIL_ALTERNATING_SUM_2010_01_11_H_
  #define _UTIL_ALTERNATING_SUM_2010_01_11_H_

  namespace Util
  {
    template<typename T1, typename T2 = T1> struct alternating_sum
    {
    private:

      alternating_sum& operator=(const alternating_sum&);

      bool b_neg_term;
      const T2 initial;

    public:

      alternating_sum(const bool b_neg = false, const T2& init = T2(0)) : b_neg_term(b_neg),
                                                                          initial   (init) { }
                                                                          
      T1 operator()(const T1& sum, const T2& ck)
      {
        const T1 the_sum = (!b_neg_term ? (sum + ck) : (sum - ck));
        b_neg_term = !b_neg_term;
        return the_sum + initial;
      }
    };
  }

#endif // _UTIL_ALTERNATING_SUM_2010_01_11_H_
