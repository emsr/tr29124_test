#include <functions/complex/e_float_complex.h>
#include <functions/elementary/elementary.h>

#if defined(__GNUC__)
  static inline INT32 _isnan (double x) { return static_cast<INT32>(std::isnan   <double>(x)); }
  static inline INT32 _finite(double x) { return static_cast<INT32>(std::isfinite<double>(x)); }
#endif

e_float ef::floor(const e_float& x)
{
  if(!ef::isfinite(x) || ef::isint(x)) { return x; }

  return ef::isneg(x) ? ef::integer_part(x - ef::one())
                      : ef::integer_part(x);
}

e_float ef::ceil(const e_float& x)
{
  if(!ef::isfinite(x) || ef::isint(x)) { return x; }

  return ef::isneg(x) ? ef::integer_part(x)
                      : ef::integer_part(x + ef::one());
}

INT32 ef::sgn(const e_float& x)
{
  if(ef::iszero(x))
  {
    return static_cast<INT32>(0);
  }
  else
  {
    return (ef::isneg(x) ? static_cast<INT32>(-1) : static_cast<INT32>(1));
  }
}

bool ef::isfinite(const double x) { return static_cast<INT32>(::_finite(x)) == static_cast<INT32>(1); } // NOCOVER_LINE
bool ef::isnan   (const double x) { return static_cast<INT32>(::_isnan (x)) == static_cast<INT32>(1); } // NOCOVER_LINE

double ef::to_double(const e_float& x)    { return x.extract_double(); }
double ef::to_double(const ef_complex& z) { return ef::to_double(z.real()); } // NOCOVER_LINE

INT64 ef::to_int64(const double x)      { return static_cast<INT64>(x); } // NOCOVER_LINE
INT64 ef::to_int64(const e_float& x)    { return x.extract_int64(); }
INT64 ef::to_int64(const ef_complex& z) { return ef::to_int64(z.real()); }

// NOCOVER_BLK_BEG
bool ef::isint(const double x)
{
  static const double delta = std::numeric_limits<double>::min() * static_cast<double>(2.0);

  const double xx = ::fabs(x);

  if((xx - ::floor(xx)) < delta)
  {
    return true;
  }
  else if((::ceil(xx) - xx) < delta)
  {
    return true;
  }
  else
  {
    return false;
  }
}

INT32 ef::to_int32(const double x)
{
  static const INT64 n32_max = static_cast<INT64>(std::numeric_limits<INT32>::max());
  static const INT64 n32_min = static_cast<INT64>(std::numeric_limits<INT32>::min());

  INT64 n64 = ef::to_int64(x);

  if(n64 < n32_min) { n64 = n32_min; }
  if(n64 > n32_max) { n64 = n32_max; }

  return static_cast<INT32>(n64);
}
// NOCOVER_BLK_END

INT32 ef::to_int32(const e_float& x)
{
  static const INT64 n32_max = static_cast<INT64>(std::numeric_limits<INT32>::max());
  static const INT64 n32_min = static_cast<INT64>(std::numeric_limits<INT32>::min());

  INT64 n64 = ef::to_int64(x);

  if(n64 < n32_min) { n64 = n32_min; }
  if(n64 > n32_max) { n64 = n32_max; }

  return static_cast<INT32>(n64);
}

INT32 ef::to_int32(const ef_complex& z)
{
  return ef::to_int32(z.real());
}

void ef::to_parts(const e_float& x, double& mantissa, INT64& exponent)
{
  x.extract_parts(mantissa, exponent);
}

e_float ef::integer_part(const e_float& x)
{
  return x.extract_integer_part();
}

e_float ef::decimal_part(const e_float& x)
{
  return x.extract_decimal_part();
}

// NOCOVER_BLK_BEG
bool ef::small_arg(const double x)
{
  static const double one_sixth = static_cast<double>(1.0) / static_cast<double>(6.0);
  static const double small_tol = ::pow(std::numeric_limits<double>::epsilon(), one_sixth);

  return ::fabs(x) < small_tol;
}
// NOCOVER_BLK_END

bool ef::small_arg(const e_float& x)
{
  static const double lim_d = static_cast<double>(static_cast<INT32>(ef::tol())) / static_cast<double>(10.0);
  static const INT64  lim_n = static_cast<INT64>(lim_d);
  static const INT64  lim   = (lim_n < 6 ? 6 : lim_n);

  return x.order() < -lim;
}

bool ef::small_arg(const ef_complex& z)
{
  return ef::small_arg(efz::abs(z));
}

// NOCOVER_BLK_BEG
bool ef::large_arg(const double x)
{
  static const double one_sixth = static_cast<double>(1.0) / static_cast<double>(6.0);
  static const double small_tol = ::pow(std::numeric_limits<double>::epsilon(), one_sixth);
  static const double large_tol = static_cast<double>(1.0) / small_tol;

  return ::fabs(x) > large_tol;
}
// NOCOVER_BLK_END

bool ef::large_arg(const e_float& x)
{
  static const double lim_d = static_cast<double>(static_cast<INT32>(ef::tol())) / static_cast<double>(10.0);
  static const INT64  lim_n = static_cast<INT64>(lim_d);
  static const INT64  lim   = (lim_n < 6 ? 6 : lim_n);

  return x.order() > lim;
}

bool ef::large_arg(const ef_complex& z) { return ef::large_arg(z.real()); } // NOCOVER_LINE

bool ef::near_one(const double x)      { return ef::small_arg(::fabs(static_cast<double>(1.0) - x)); } // NOCOVER_LINE
bool ef::near_one(const e_float& x)    { return ef::small_arg(ef::fabs(ef::one() - x)); }
bool ef::near_one(const ef_complex& z) { return ef::near_one(z.real()) && ef::iszero(z.imag()); } // NOCOVER_LINE

// NOCOVER_BLK_BEG
bool ef::near_int(const double x)
{
  if(ef::isint(x))
  {
    return true;
  }
  else
  {
    const double xx = ::fabs(x);

    if(ef::small_arg(xx - ::floor(xx)))
    {
      return true;
    }
    else if(ef::small_arg(::ceil(xx) - xx))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}
// NOCOVER_BLK_END

bool ef::near_int(const e_float& x)
{
  if(ef::isint(x))
  {
    return true;
  }
  else
  {
    const e_float xx = ef::fabs(x);

    if(ef::small_arg(xx - ef::floor(xx)))
    {
      return true;
    }
    else if(ef::small_arg(ef::ceil(xx) - xx))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

bool ef::near_int(const ef_complex& z)
{
  if(ef::isint(z))
  {
    return true;
  }
  else
  {
    return ef::iszero(z.imag()) && ef::near_int(z.real());
  }
}
