
#include <functions/complex/e_float_complex.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary_complex.h>
#include <functions/elementary/elementary_math.h>
#include <functions/elementary/elementary_trans.h>
#include <functions/elementary/elementary_trig.h>
#include <utility/util_power_x_pow_n.h>

// NOCOVER_BLK_BEG
std::complex<double> efz::to_double(const ef_complex& z)
{
  const double rx = ef::to_double(z.real());
  const double ix = ef::to_double(z.imag());
  return std::complex<double>(rx, ix);
}
// NOCOVER_BLK_END

// NOCOVER_BLK_BEG
bool operator==(const ef_complex& u, const e_float& v) { return ((u.real() == v) &&  ef::iszero(u.imag())); }
bool operator!=(const ef_complex& u, const e_float& v) { return ((u.real() != v) || !ef::iszero(u.imag())); }

bool operator==(const e_float& u, const ef_complex& v) { return ((u == v.real()) &&  ef::iszero(v.imag())); }
bool operator!=(const e_float& u, const ef_complex& v) { return ((u != v.real()) || !ef::iszero(v.imag())); }

bool operator==(const ef_complex& u, const INT32 n) { return ((u.real() == n) &&  ef::iszero(u.imag())); }
bool operator!=(const ef_complex& u, const INT32 n) { return ((u.real() != n) || !ef::iszero(u.imag())); }

bool operator==(const INT32 n, const ef_complex& u) { return ((n == u.real()) &&  ef::iszero(u.imag())); }
bool operator!=(const INT32 n, const ef_complex& u) { return ((n != u.real()) || !ef::iszero(u.imag())); }
// NOCOVER_BLK_END

e_float efz::abs(const ef_complex& z) { return ef::sqrt(efz::norm(z)); }
e_float efz::arg(const ef_complex& z) { return ef::atan2(z.imag(), z.real()); }

ef_complex efz::polar(const e_float& mod, const e_float& arg)
{
  e_float s, c;
  ef::sincos(arg, &s, &c);
  return ef_complex(c, s) * mod;
}

ef_complex efz::sin(const ef_complex& z)
{
  e_float sin_x, cos_x, sinh_y, cosh_y;

  ef::sincos  (z.real(), &sin_x, &cos_x);
  ef::sinhcosh(z.imag(), &sinh_y, &cosh_y);

  return ef_complex(sin_x * cosh_y, cos_x * sinh_y);
}

ef_complex efz::cos(const ef_complex& z)
{
  e_float sin_x, cos_x, sinh_y, cosh_y;

  ef::sincos  (z.real(), &sin_x,  &cos_x);
  ef::sinhcosh(z.imag(), &sinh_y, &cosh_y);

  return ef_complex(cos_x * cosh_y, -(sin_x * sinh_y));
}

void efz::sincos(const ef_complex& z, ef_complex* const p_sin, ef_complex* const p_cos)
{
  e_float sin_x, cos_x, sinh_y, cosh_y;

  ef::sincos  (z.real(), &sin_x,  &cos_x);
  ef::sinhcosh(z.imag(), &sinh_y, &cosh_y);

  if(p_sin)
  {
    *p_sin = ef_complex(sin_x * cosh_y, cos_x * sinh_y);
  }
  
  if(p_cos)
  {
    *p_cos = ef_complex(cos_x * cosh_y, -(sin_x * sinh_y));
  }
}

ef_complex efz::tan(const ef_complex& z)
{
  ef_complex s, c;
  efz::sincos(z, &s, &c);
  return s * efz::inv(c);
}

ef_complex efz::csc(const ef_complex& z) { return efz::inv(efz::sin(z)); }
ef_complex efz::sec(const ef_complex& z) { return efz::inv(efz::cos(z)); }
ef_complex efz::cot(const ef_complex& z) { return efz::inv(efz::tan(z)); }

ef_complex efz::asin(const ef_complex& z)
{
  return -efz::iz(efz::log(efz::iz(z) + efz::sqrt(ef::one() - (z * z))));
}

ef_complex efz::acos(const ef_complex& z)
{
  return ef_complex(ef::pi_half(), ef::zero()) - efz::asin(z);
}

ef_complex efz::atan(const ef_complex& z)
{
  const ef_complex izz = efz::iz(z);
  return efz::iz(efz::log(ef::one() - izz) - efz::log(ef::one() + izz)) / static_cast<INT32>(2);
}

ef_complex efz::inv(const ef_complex& z)
{
  // Compute inverse 1 / (x + iy) = (x - iy) / (x^2 + y^2)
  return ef_complex(z.real(), -z.imag()) * efz::norm(z).calculate_inv();
}

ef_complex efz::sqrt(const ef_complex& z)
{
  // Equation from MPFUN documentation page 12.
  // See: http://www.nersc.gov/~dhb/mpdist/mpdist.html

  // Pure zero?
  if(ef::iszero(z))
  {
    return ef::zero();
  }
  else
  {
    // sqrt(*this) = (s, I / 2s)     for R >= 0
    //               (|I| / 2s, +-s) for R <  0
    // where s = sqrt{ [ |R| + sqrt(R^2 + I^2) ] / 2 },
    // and the +- sign is the same as the sign of I.
    const e_float s = ef::sqrt((ef::fabs(z.real()) + efz::abs(z)) / static_cast<INT32>(2));
    
    if(ef::iszero(z.real()) || !ef::isneg(z.real()))
    {
      return ef_complex(s, (z.imag() / s) / static_cast<INT32>(2));
    }
    else
    {
      const bool imag_is_pos = ef::iszero(z.imag()) || !ef::isneg(z.imag());

      return ef_complex((ef::fabs(z.imag()) / s) / static_cast<INT32>(2), (imag_is_pos ? s : -s));
    }
  }
}

ef_complex efz::exp(const ef_complex& z)
{
  e_float s, c;
  ef::sincos(z.imag(), &s, &c);
  return ef_complex(c , s) * ef::exp(z.real());
}

ef_complex efz::log(const ef_complex& z)
{
  return ef_complex(ef::log(efz::norm(z)) / static_cast<INT32>(2), ef::atan2(z.imag(), z.real()));
}

ef_complex efz::log10(const ef_complex& z)
{
  return efz::log(z) / ef::ln10();
}

ef_complex efz::loga(const ef_complex& a, const ef_complex& z)
{
  return efz::log(z) / efz::log(a);
}

ef_complex efz::pown(const ef_complex& z, const INT64 p)
{
  return Util::x_pow_n_template<ef_complex>(z, p);
}

ef_complex efz::pow(const ef_complex& z, const ef_complex& a)
{
  return efz::exp(a * efz::log(z));
}

ef_complex efz::rootn(const ef_complex& z, const INT32 p)
{
  if(p < static_cast<INT32>(0))
  {
    return efz::pown(ef::one() / z, static_cast<INT64>(-p));
  }
  else if(p == static_cast<INT32>(0))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  else if(p == static_cast<INT32>(1))
  {
    return z;
  }
  else
  {
    return efz::polar(ef::rootn(efz::norm(z), static_cast<INT32>(2) * p), efz::arg(z) / p);
  }
}

ef_complex efz::sinh(const ef_complex& z)
{
  e_float sin_y, cos_y, sinh_x, cosh_x;

  ef::sincos  (z.imag(), &sin_y,  &cos_y);
  ef::sinhcosh(z.real(), &sinh_x, &cosh_x);

  return ef_complex(cos_y * sinh_x, cosh_x * sin_y);
}

ef_complex efz::cosh(const ef_complex& z)
{
  e_float sin_y, cos_y, sinh_x, cosh_x;
  
  ef::sincos  (z.imag(), &sin_y,  &cos_y);
  ef::sinhcosh(z.real(), &sinh_x, &cosh_x);

  return ef_complex(cos_y * cosh_x, sin_y * sinh_x);
}

void efz::sinhcosh(const ef_complex& z, ef_complex* const p_sinh, ef_complex* const p_cosh)
{
  e_float sin_y, cos_y, sinh_x, cosh_x;

  ef::sincos  (z.imag(), &sin_y,  &cos_y);
  ef::sinhcosh(z.real(), &sinh_x, &cosh_x);

  if(p_sinh)
  {
    *p_sinh = ef_complex(cos_y * sinh_x, cosh_x * sin_y);
  }

  if(p_cosh)
  {
    *p_cosh = ef_complex(cos_y * cosh_x, sin_y  * sinh_x);
  }
}

ef_complex efz::tanh(const ef_complex& z)
{
  ef_complex sh, ch;
  efz::sinhcosh(z, &sh, &ch);
  return sh * efz::inv(ch);
}

ef_complex efz::asinh(const ef_complex& z)
{
  return efz::log(z + efz::sqrt((z * z) + ef::one()));
}

ef_complex efz::acosh(const ef_complex& z)
{
  const ef_complex zp(z.real() + ef::one(), z.imag());
  const ef_complex zm(z.real() - ef::one(), z.imag());

  return efz::log(z + (zp * efz::sqrt(zm / zp)));
}

ef_complex efz::atanh(const ef_complex& z)
{
  return (efz::log(ef::one() + z) - efz::log(ef::one() - z)) / static_cast<INT32>(2);
}
