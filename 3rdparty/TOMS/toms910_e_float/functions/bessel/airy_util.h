
#ifndef _AIRY_UTIL_2008_02_20_H_
  #define _AIRY_UTIL_2008_02_20_H_

  #include <e_float/e_float.h>
  #include <functions/constants/constants.h>
  #include <functions/elementary/elementary.h>

  namespace AiryUtil
  {
    e_float make_zeta(const e_float& sqrt_x);

    const std::tr1::array<e_float, 2u>& Three_Pow_1_3and2_3  (void);
    const std::tr1::array<e_float, 2u>& Ai_and_AiPrime_OfZero(void);
    const std::tr1::array<e_float, 2u>& Bi_and_BiPrime_OfZero(void);
    const std::tr1::array<e_float, 2u>& Gamma_Of_1_3and2_3   (void);

    void BesselJ_Of_1_3and2_3(const e_float& x, e_float& J1_3, e_float& Jm1_3, e_float& J2_3, e_float& Jm2_3);
  }
  
#endif // _AIRY_UTIL_2008_02_20_H_
