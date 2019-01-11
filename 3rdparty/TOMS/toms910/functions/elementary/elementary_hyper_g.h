
#ifndef _ELEMENTARY_HYPER_G_2009_12_09_H_
  #define _ELEMENTARY_HYPER_G_2009_12_09_H_

  #include <deque>

  namespace ef
  {
    e_float hyp0F0(const e_float& x);
    e_float hyp0F1(const e_float& b, const e_float& x);
    e_float hyp1F0(const e_float& a, const e_float& x);
    e_float hyp1F1(const e_float& a, const e_float& b, const e_float& x);
    e_float hyp2F0(const e_float& a, const e_float& b, const e_float& x);
    e_float hyp2F1(const e_float& a, const e_float& b, const e_float& c, const e_float& x);
    e_float hypPFQ(const std::deque<e_float>& a, const  std::deque<e_float>& b, const e_float& x);
  }

#endif // _ELEMENTARY_HYPER_G_2009_12_09_H_
