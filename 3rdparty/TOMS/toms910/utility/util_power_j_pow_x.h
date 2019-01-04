
#ifndef _UTIL_POWER_J_POW_X_2009_01_25_H_
  #define _UTIL_POWER_J_POW_X_2009_01_25_H_

  #include <map>

  #include <functions/complex/e_float_complex.h>

  namespace Util
  {
    e_float    j_pow_x(const UINT32 j, const e_float&    x, std::map<UINT32, e_float>&    n_pow_x_prime_factor_map);
    ef_complex j_pow_x(const UINT32 j, const ef_complex& x, std::map<UINT32, ef_complex>& n_pow_x_prime_factor_map);
  }

#endif // _UTIL_POWER_J_POW_X_2009_01_25_H_
