// $HOME/bin/bin/g++ -o test_promote_fp test_promote_fp.cpp

#include <type_traits>
#include "promote_fp.h"

int
main()
{
  int i;
  short j;
  float f;
  double d;
  long double ld;

  __promote_fp_t<int, short> x;
  __promote_fp_t<int, long double> y;
  __promote_fp_t<float, float, int, float> z;
  static_assert(std::is_same<decltype(x), double>::value, "");
  static_assert(std::is_same<decltype(y), long double>::value, "");
  static_assert(std::is_same<decltype(z), double>::value, "");

  static_assert(std::is_same<decltype(d), double>::value, "");

  double & rd = d;
  //static_assert(std::is_same<decltype(rd), double>::value, "");
  static_assert(std::is_same<__promote_fp_t<decltype(rd)>, double>::value, "");

  const double cd = d;
  //static_assert(std::is_same<decltype(cd), double>::value, "");
  static_assert(std::is_same<__promote_fp_t<decltype(cd)>, double>::value, "");

  volatile double vd = d;
  //static_assert(std::is_same<decltype(vd), double>::value, "");
  static_assert(std::is_same<__promote_fp_t<decltype(vd)>, double>::value, "");

  const volatile double cvd = d;
  //static_assert(std::is_same<decltype(cvd), double>::value, "");
  static_assert(std::is_same<__promote_fp_t<decltype(cvd)>, double>::value, "");
}