#include <algorithm>

#include <e_float/e_float.h>
#include <utility/util_lexical_cast.h>

std::ostream& operator<<(std::ostream& os, const e_float_base& f)
{
  std::string str;
  f.wr_string(str, os);
  return (os << str);
}

std::istream& operator>>(std::istream& is, e_float_base& f)
{
  std::string str;
  static_cast<void>(is >> str);
  f.rd_string(str.c_str());
  return is;
}

// NOCOVER_BLK_BEG
e_float e_float_base::get_pos(void) const { return *static_cast<const e_float*>(static_cast<const void*>(this)); }
e_float e_float_base::get_neg(void) const { return (e_float(*static_cast<const e_float*>(static_cast<const void*>(this)))).negate(); }
e_float e_float_base::get_abs(void) const { return (!isneg() ? get_pos() : get_neg()); }
INT64   e_float_base::get_int(void) const { return extract_int64(); }
double  e_float_base::get_flt(void) const { return extract_double(); }
// NOCOVER_BLK_END

// NOCOVER_BLK_BEG
std::string e_float_base::get_str(void) const
{
  const std::streamsize my_prec = static_cast<std::streamsize>(std::numeric_limits<e_float>::digits10);

  std::stringstream ss;

  static_cast<void>(ss.setf(std::ios_base::showpos | std::ios_base::scientific));
  static_cast<void>(ss.precision(my_prec));

  ss << *this;

  return ss.str();  
}
// NOCOVER_BLK_END

const std::string::size_type& e_float_base::width_of_exponent_field(void)
{
  static const std::string::size_type width_of_e_n64 =
      Util::lexical_cast(std::numeric_limits<INT64>::max()).length();

  static const std::string::size_type width_of_e_long =
      Util::lexical_cast(std::numeric_limits<long>::max()).length();

  static const std::string::size_type width_of_e =
      std::max(width_of_e_n64, width_of_e_long);

  return width_of_e;
}

// NOCOVER_BLK_BEG
e_float e_float_base::my_cbrt         (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_rootn        (const e_float& x, const UINT32 p) { static_cast<void>(x); static_cast<void>(p); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_exp          (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_log          (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_sin          (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_cos          (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_tan          (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_asin         (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_acos         (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_atan         (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_sinh         (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_cosh         (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_tanh         (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_asinh        (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_acosh        (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_atanh        (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_gamma        (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_riemann_zeta (const e_float& x) { static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_cyl_bessel_jn(const INT32 n, const e_float& x) { static_cast<void>(n); static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
e_float e_float_base::my_cyl_bessel_yn(const INT32 n, const e_float& x) { static_cast<void>(n); static_cast<void>(x); return std::numeric_limits<e_float>::quiet_NaN(); }
// NOCOVER_BLK_END
