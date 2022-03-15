#ifndef ROUNDED_CMATH_H
#define ROUNDED_CMATH_H 1

#include <cmath>

#include "rounding_tools.h"

long double
sqrt(long double d, fpround rm);
double
sqrt(double d, fpround rm);
float
sqrt(float f, fpround rm);

double
cbrt(double d, fpround rm);
float
cbrt(float f, fpround rm);

double
sin(double d, fpround rm);
float
sin(float f, fpround rm);

double
cos(double d, fpround rm);
float
cos(float f, fpround rm);

double
tan(double d, fpround rm);
float
tan(float f, fpround rm);

double
asin(double d, fpround rm);
float
asin(float f, fpround rm);

double
acos(double d, fpround rm);
float
acos(float f, fpround rm);

double
atan(double d, fpround rm);
float
atan(float f, fpround rm);

double
sinh(double d, fpround rm);
float
sinh(float f, fpround rm);

double
cosh(double d, fpround rm);
float
cosh(float f, fpround rm);

double
tanh(double d, fpround rm);
float
tanh(float f, fpround rm);

double
asinh(double d, fpround rm);
float
asinh(float f, fpround rm);

double
acosh(double d, fpround rm);
float
acosh(float f, fpround rm);

double
atanh(double d, fpround rm);
float
atanh(float f, fpround rm);

double
log(double d, fpround rm);
float
log(float f, fpround rm);

double
log10(double d, fpround rm);
float
log10(float f, fpround rm);

double
log2(double d, fpround rm);
float
log2(float f, fpround rm);

double
exp(double d, fpround rm);
float
exp(float f, fpround rm);

//double
//exp10(double d, fpround rm);
//float
//exp10(float f, fpround rm);

double
exp2(double d, fpround rm);
float
exp2(float f, fpround rm);

// So, ... does glibc maths do rounding all by themselves?
long double
sqrt(long double d, fpround rm)
{
  round_sentinel rs(rm);
  return std::sqrt(d);
}
long double
sin(long double d, fpround rm)
{
  round_sentinel rs(rm);
  return std::sin(d);
}
long double
cos(long double d, fpround rm)
{
  round_sentinel rs(rm);
  return std::cos(d);
}
long double
tan(long double d, fpround rm)
{
  round_sentinel rs(rm);
  return std::tan(d);
}

inline double
sqrt(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::sqrt(res);

  return res;
}

inline float
sqrt(float f, fpround rm)
{
  double res = f;

  res = sqrt(res, rm);

  return to_float(res, rm);
}


inline double
cbrt(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::cbrt(res);

  return res;
}

inline float
cbrt(float f, fpround rm)
{
  double res = f;

  res = cbrt(res, rm);

  return to_float(res, rm);
}

inline double
sin(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::sin(res);

  return res;
}

inline float
sin(float f, fpround rm)
{
  double res = f;

  res = sin(res, rm);

  return to_float(res, rm);
}

inline double
cos(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::cos(res);

  return res;
}

inline float
cos(float f, fpround rm)
{
  double res = f;

  res = cos(res, rm);

  return to_float(res, rm);
}

// tan() for double interval
//
inline double
tan(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::tan(res);

  return res;
}

inline float
tan(float f, fpround rm)
{
  double res = f;

  res = tan(res, rm);

  return to_float(res, rm);
}

inline double
asin(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::asin(res);

  return res;
}

inline float
asin(float f, fpround rm)
{
  double res = f;

  res = asin(res, rm);

  return to_float(res, rm);
}

inline double
acos(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::acos(res);

  return res;
}

inline float
acos(float f, fpround rm)
{
  double res = f;

  res = acos(res, rm);

  return to_float(res, rm);
}

inline double
atan(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::atan(res);

  return res;
}

inline float
atan(float f, fpround rm)
{
  double res = f;

  res = atan(res, rm);

  return to_float(res, rm);
}



inline double
sinh(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::sinh(res);

  return res;
}

inline float
sinh(float f, fpround rm)
{
  double res = f;

  res = sinh(res, rm);

  return to_float(res, rm);
}

inline double
cosh(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::cosh(res);

  return res;
}

inline float
cosh(float f, fpround rm)
{
  double res = f;

  res = cosh(res, rm);

  return to_float(res, rm);
}

inline double
tanh(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::tanh(res);

  return res;
}

inline float
tanh(float f, fpround rm)
{
  double res = f;

  res = tanh(res, rm);

  return to_float(res, rm);
}

inline double
asinh(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::asinh(res);

  return res;
}

inline float
asinh(float f, fpround rm)
{
  double res = f;

  res = asinh(res, rm);

  return to_float(res, rm);
}

inline double
acosh(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::acosh(res);

  return res;
}

inline float
acosh(float f, fpround rm)
{
  double res = f;

  res = acosh(res, rm);

  return to_float(res, rm);
}

inline double
atanh(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::atanh(res);

  return res;
}

inline float
atanh(float f, fpround rm)
{
  double res = f;

  res = atanh(res, rm);

  return to_float(res, rm);
}

inline double
log(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::log(res);

  return res;
}

inline float
log(float f, fpround rm)
{
  double res = f;

  res = log(res, rm);

  return to_float(res, rm);
}


inline double
log10(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::log10(res);

  return res;
}

inline float
log10(float f, fpround rm)
{
  double res = f;

  res = log10(res, rm);

  return to_float(res, rm);
}


inline double
log2(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::log2(res);

  return res;
}


inline float
log2(float f, fpround rm)
{
  double res = f;

  res = log2(res, rm);

  return to_float(res, rm);
}

double
exp(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::exp(res);

  return res;
}

float
exp(float f, fpround rm)
{
  double res = f;

  res = exp(res, rm);

  return to_float(res, rm);
}
/*
double
exp10(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::exp10(res);

  return res;
}

float
exp10(float f, fpround rm)
{
  double res = f;

  res = exp10(res, rm);

  return to_float(res, rm);
}
*/
double
exp2(double d, fpround rm)
{
  long double res = d;

  round_sentinel rs(rm);

  //FIXME!
  res = std::exp2(res);

  return res;
}

float
exp2(float f, fpround rm)
{
  double res = f;

  res = exp2(res, rm);

  return to_float(res, rm);
}

#endif // ROUNDED_CMATH_H
