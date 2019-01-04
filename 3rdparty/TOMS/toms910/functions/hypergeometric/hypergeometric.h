
#ifndef _HYPER_GEOMETRIC_2008_01_06_H_
  #define _HYPER_GEOMETRIC_2008_01_06_H_

  #include <deque>

  #include <e_float/e_float.h>

  namespace ef
  {
    e_float hyperg_0f0    (const e_float& x);
    e_float hyperg_0f1    (const e_float& b, const e_float& x);
    e_float hyperg_0f1_reg(const e_float& a, const e_float& x);
    e_float hyperg_1f0    (const e_float& a, const e_float& x);
    e_float hyperg_1f1    (const e_float& a, const e_float& b, const e_float& x);
    e_float hyperg_1f1_reg(const e_float& a, const e_float& b, const e_float& x);
    e_float hyperg_2f0    (const e_float& a, const e_float& b, const e_float& x);
    e_float hyperg_2f1    (const e_float& a, const e_float& b, const e_float& c, const e_float& x);
    e_float hyperg_2f1_reg(const e_float& a, const e_float& b, const e_float& c, const e_float& x);
    e_float hyperg_pfq    (const std::deque<e_float>& a, const  std::deque<e_float>& b, const e_float& x);

    inline e_float conf_hyperg(const e_float& a, const e_float& c, const e_float& x)                   { return ef::hyperg_1f1(a, c, x); }
    inline e_float      hyperg(const e_float& a, const e_float& b, const e_float& c, const e_float& x) { return ef::hyperg_2f1(a, b, c, x); }
  }
  
#endif // _HYPER_GEOMETRIC_2008_01_06_H_
