#ifndef ROUNDED_CMATH_H
#define ROUNDED_CMATH_H 1

#include <cmath>

#include "rounding_tools.h"

long double
sqrt(long double __d, fpround __rm);
double
sqrt(double __d, fpround __rm);
float
sqrt(float __f, fpround __rm);

double
cbrt(double __d, fpround __rm);
float
cbrt(float __f, fpround __rm);

double
sin(double __d, fpround __rm);
float
sin(float __f, fpround __rm);

double
cos(double __d, fpround __rm);
float
cos(float __f, fpround __rm);

double
tan(double __d, fpround __rm);
float
tan(float __f, fpround __rm);

double
asin(double __d, fpround __rm);
float
asin(float __f, fpround __rm);

double
acos(double __d, fpround __rm);
float
acos(float __f, fpround __rm);

double
atan(double __d, fpround __rm);
float
atan(float __f, fpround __rm);

double
sinh(double __d, fpround __rm);
float
sinh(float __f, fpround __rm);

double
cosh(double __d, fpround __rm);
float
cosh(float __f, fpround __rm);

double
tanh(double __d, fpround __rm);
float
tanh(float __f, fpround __rm);

double
asinh(double __d, fpround __rm);
float
asinh(float __f, fpround __rm);

double
acosh(double __d, fpround __rm);
float
acosh(float __f, fpround __rm);

double
atanh(double __d, fpround __rm);
float
atanh(float __f, fpround __rm);

double
log(double __d, fpround __rm);
float
log(float __f, fpround __rm);

double
log10(double __d, fpround __rm);
float
log10(float __f, fpround __rm);

double
log2(double __d, fpround __rm);
float
log2(float __f, fpround __rm);

double
exp(double __d, fpround __rm);
float
exp(float __f, fpround __rm);

//double
//exp10(double __d, fpround __rm);
//float
//exp10(float __f, fpround __rm);

double
exp2(double __d, fpround __rm);
float
exp2(float __f, fpround __rm);

// So, ... does glibc maths do rounding all by themselves?
long double
sqrt(long double __d, fpround __rm)
{
  round_sentinel __rs(__rm);
  return std::sqrt(__d);
}
long double
sin(long double __d, fpround __rm)
{
  round_sentinel __rs(__rm);
  return std::sin(__d);
}
long double
cos(long double __d, fpround __rm)
{
  round_sentinel __rs(__rm);
  return std::cos(__d);
}
long double
tan(long double __d, fpround __rm)
{
  round_sentinel __rs(__rm);
  return std::tan(__d);
}

inline double
sqrt(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::sqrt(__res);

  return __res;
}

inline float
sqrt(float __f, fpround __rm)
{
  double __res = __f;

  __res = sqrt(__res, __rm);

  return to_float(__res, __rm);
}


inline double
cbrt(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::cbrt(__res);

  return __res;
}

inline float
cbrt(float __f, fpround __rm)
{
  double __res = __f;

  __res = cbrt(__res, __rm);

  return to_float(__res, __rm);
}

inline double
sin(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::sin(__res);

  return __res;
}

inline float
sin(float __f, fpround __rm)
{
  double __res = __f;

  __res = sin(__res, __rm);

  return to_float(__res, __rm);
}

inline double
cos(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::cos(__res);

  return __res;
}

inline float
cos(float __f, fpround __rm)
{
  double __res = __f;

  __res = cos(__res, __rm);

  return to_float(__res, __rm);
}

// tan() for double interval
//
inline double
tan(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::tan(__res);

  return __res;
}

inline float
tan(float __f, fpround __rm)
{
  double __res = __f;

  __res = tan(__res, __rm);

  return to_float(__res, __rm);
}

inline double
asin(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::asin(__res);

  return __res;
}

inline float
asin(float __f, fpround __rm)
{
  double __res = __f;

  __res = asin(__res, __rm);

  return to_float(__res, __rm);
}

inline double
acos(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::acos(__res);

  return __res;
}

inline float
acos(float __f, fpround __rm)
{
  double __res = __f;

  __res = acos(__res, __rm);

  return to_float(__res, __rm);
}

inline double
atan(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::atan(__res);

  return __res;
}

inline float
atan(float __f, fpround __rm)
{
  double __res = __f;

  __res = atan(__res, __rm);

  return to_float(__res, __rm);
}



inline double
sinh(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::sinh(__res);

  return __res;
}

inline float
sinh(float __f, fpround __rm)
{
  double __res = __f;

  __res = sinh(__res, __rm);

  return to_float(__res, __rm);
}

inline double
cosh(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::cosh(__res);

  return __res;
}

inline float
cosh(float __f, fpround __rm)
{
  double __res = __f;

  __res = cosh(__res, __rm);

  return to_float(__res, __rm);
}

inline double
tanh(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::tanh(__res);

  return __res;
}

inline float
tanh(float __f, fpround __rm)
{
  double __res = __f;

  __res = tanh(__res, __rm);

  return to_float(__res, __rm);
}

inline double
asinh(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::asinh(__res);

  return __res;
}

inline float
asinh(float __f, fpround __rm)
{
  double __res = __f;

  __res = asinh(__res, __rm);

  return to_float(__res, __rm);
}

inline double
acosh(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::acosh(__res);

  return __res;
}

inline float
acosh(float __f, fpround __rm)
{
  double __res = __f;

  __res = acosh(__res, __rm);

  return to_float(__res, __rm);
}

inline double
atanh(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::atanh(__res);

  return __res;
}

inline float
atanh(float __f, fpround __rm)
{
  double __res = __f;

  __res = atanh(__res, __rm);

  return to_float(__res, __rm);
}

inline double
log(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::log(__res);

  return __res;
}

inline float
log(float __f, fpround __rm)
{
  double __res = __f;

  __res = log(__res, __rm);

  return to_float(__res, __rm);
}


inline double
log10(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::log10(__res);

  return __res;
}

inline float
log10(float __f, fpround __rm)
{
  double __res = __f;

  __res = log10(__res, __rm);

  return to_float(__res, __rm);
}


inline double
log2(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::log2(__res);

  return __res;
}


inline float
log2(float __f, fpround __rm)
{
  double __res = __f;

  __res = log2(__res, __rm);

  return to_float(__res, __rm);
}

double
exp(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::exp(__res);

  return __res;
}

float
exp(float __f, fpround __rm)
{
  double __res = __f;

  __res = exp(__res, __rm);

  return to_float(__res, __rm);
}
/*
double
exp10(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::exp10(__res);

  return __res;
}

float
exp10(float __f, fpround __rm)
{
  double __res = __f;

  __res = exp10(__res, __rm);

  return to_float(__res, __rm);
}
*/
double
exp2(double __d, fpround __rm)
{
  long double __res = __d;

  round_sentinel __rs(__rm);

  //FIXME!
  __res = std::exp2(__res);

  return __res;
}

float
exp2(float __f, fpround __rm)
{
  double __res = __f;

  __res = exp2(__res, __rm);

  return to_float(__res, __rm);
}

#endif // ROUNDED_CMATH_H
