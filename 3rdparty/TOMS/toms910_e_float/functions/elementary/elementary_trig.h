
#ifndef _ELEMENTARY_TRIG_2009_12_09_H_
  #define _ELEMENTARY_TRIG_2009_12_09_H_

  namespace ef
  {
    void    sincos  (const e_float& x, e_float* const p_sin, e_float* const p_cos);
    e_float sin     (const e_float& x);
    e_float cos     (const e_float& x);
    e_float tan     (const e_float& x);
    e_float csc     (const e_float& x);
    e_float sec     (const e_float& x);
    e_float cot     (const e_float& x);
    e_float asin    (const e_float& x);
    e_float acos    (const e_float& x);
    e_float atan    (const e_float& x);
    e_float atan2   (const e_float& y, const e_float& x);
  }

#endif // _ELEMENTARY_TRIG_2009_12_09_H_
