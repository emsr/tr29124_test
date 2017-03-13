/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -Izeta -o debug_polylog debug_polylog.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./debug_polylog > debug_polylog.txt

LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH $HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o debug_polylog debug_polylog.cpp -lquadmath -Lwrappers/debug -lwrap_cephes

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -Izeta -o debug_polylog debug_polylog.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
PATH=wrappers/debug:$PATH ./debug_polylog > debug_polylog.txt
*/

#include <ext/cmath>

#include "zeta/riemann_zeta.tcc"
#include "zeta/PolyLog.h"

int
main()
{
  //const auto _S_pi = __gnu_cxx::__const_pi(1.0);
  const auto _S_2pi = __gnu_cxx::__const_2_pi(1.0);

  //  Start with asymp. DONE
  //double s = -12.9;
  //double x = 0.0;
  //std::__detail::__polylog_exp_asymp(s, 100.0 * std::polar(1.0, _S_2pi * x));
  //PolyLog_Exp_asym(s, 100.0 * std::polar(1.0, _S_2pi * x));

  // OK, now on to Test series 6 [PolyLog_Exp_pos(s, exp(i2pix)]: DONE
  //s = 0.145;
  //x = 0.05;
  //std::__detail::__polylog_exp_pos(s, std::polar(1.0, _S_2pi * x));
  //PolyLog_Exp_pos(s, std::polar(1.0, _S_2pi * x));
  // Guess what? my __sincos_pi (and actually __builtin_sincos called sincos...
  // ... which was picked up by cephes/cmath/sincos.c.
  // Which assumes the angle is in degrees. LOL!

  // OK, now on to Test series 6 [PolyLog_Exp_pos(s, exp(i2pix)]:
  int s = 3;
  double x = 0.05;
  // Complex arg seems OK.. But why issue with Clausen?
  //std::__detail::__polylog_exp_pos(s, std::polar(1.0, _S_2pi * x));
  //PolyLog_Exp_pos(s, std::polar(1.0, _S_2pi * x));
  std::__detail::__polylog_exp_pos(s, x);
  PolyLog_Exp_pos(s, x);

  //Test function 3 [PolyLog_Exp_neg(s<0, exp(i2pik)]:
  //s = -7.9;
  //std::__detail::__polylog_exp_neg(s, std::polar(1.0, _S_2pi * x));
  //PolyLog_Exp_neg(s, std::polar(1.0, _S_2pi * x));
  // Is riemann_seta_m_1 broken for s < 0?

}
