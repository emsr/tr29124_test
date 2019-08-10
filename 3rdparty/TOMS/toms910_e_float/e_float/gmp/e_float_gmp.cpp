
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <e_float/e_float.h>
#include <e_float/gmp/e_float_gmp_protos.h>
#include <functions/constants/constants.h>
#include <utility/util_lexical_cast.h>

void gmp::e_float::init(void)
{
  static bool precision_is_initialized = false;

  if(!precision_is_initialized)
  {
    precision_is_initialized = true;

    ::mpf_set_default_prec(static_cast<unsigned long int>(ef_digits2));
  }
}

const double& gmp::e_float::d_log2(void)
{
  static const double value_log2 = static_cast<double>(0.3010299956639811952137389);
  return value_log2;
}

const INT64& gmp::e_float::max_exp2(void)
{
  static const INT64 val_max_exp2 = static_cast<INT64>(static_cast<double>(ef_max_exp10) / d_log2());
  return val_max_exp2;
}

const INT64& gmp::e_float::min_exp2(void)
{
  static const INT64 val_min_exp2 = static_cast<INT64>(static_cast<double>(ef_min_exp10) / d_log2());
  return val_min_exp2;
}

gmp::e_float::e_float(const INT32 n) : fpclass  (ef_finite),
                                       prec_elem(ef_digits10_tol)
{
  init();

  const bool b_neg = (n < static_cast<INT32>(0));

  from_uint32(UINT32(b_neg ? -n : n));

  if(b_neg)
  {
    ::mpf_neg(rop, rop);
  }
}

gmp::e_float::e_float(const INT64 n) : fpclass  (ef_finite),
                                       prec_elem(ef_digits10_tol)
{
  init();

  const bool b_neg = (n < static_cast<INT64>(0));

  from_uint64(UINT64(b_neg ? -n : n));

  if(b_neg)
  {
    ::mpf_neg(rop, rop);
  }
}

gmp::e_float::e_float(const UINT32 u) : fpclass  (ef_finite),
                                        prec_elem(ef_digits10_tol)
{
  init();
  from_uint32(u);
}

gmp::e_float::e_float(const UINT64 u) : fpclass  (ef_finite),
                                        prec_elem(ef_digits10_tol)
{
  init();
  from_uint64(u);
}

void gmp::e_float::from_uint64(const UINT64 u)
{
  ::mpf_init_set_str(rop, Util::lexical_cast(u).c_str(), 10);
}

void gmp::e_float::from_uint32(const UINT32 u)
{
  ::mpf_init_set_ui(rop, static_cast<unsigned long int>(u));
}

gmp::e_float::e_float() : fpclass  (ef_finite),
                          prec_elem(ef_digits10_tol)
{
  init();
  ::mpf_init(rop);
}

gmp::e_float::e_float(const double d) : fpclass  (ef_finite),
                                        prec_elem(ef_digits10_tol)
{
  init();

  std::stringstream ss;

  ss << std::scientific
     << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<double>::digits10 - 1))
     << d;

  ::mpf_init_set_str(rop, ss.str().c_str(), 10);
}

gmp::e_float::e_float(const long double ld) : fpclass  (ef_finite),
                                              prec_elem(ef_digits10_tol)
{
  init();

  std::stringstream ss;

  ss << std::scientific
     << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<long double>::digits10 - 1))
     << ld;

  ::mpf_init_set_str(rop, ss.str().c_str(), 10);
}

gmp::e_float::e_float(const char* const s) : fpclass  (ef_finite),
                                             prec_elem(ef_digits10_tol)
{
  init();
  static_cast<void>(rd_string(s));
}

gmp::e_float::e_float(const std::string& str) : fpclass  (ef_finite),
                                                prec_elem(ef_digits10_tol)
{
  init();
  static_cast<void>(rd_string(str.c_str()));
}

gmp::e_float::e_float(const e_float& mp) : fpclass  (mp.fpclass),
                                           prec_elem(mp.prec_elem)
{
  init();
  ::mpf_init_set(rop, mp.rop);
}

gmp::e_float::e_float(const double mantissa, const INT64 exponent) : fpclass  (ef_finite),
                                                                     prec_elem(ef_digits10_tol)
{
  init();

  const bool mantissa_is_iszero = (::fabs(mantissa) < (std::numeric_limits<double>::min() * static_cast<double>(2.0)));

  if(mantissa_is_iszero)
  {
    if(exponent == static_cast<INT64>(0))
    {
      ::mpf_init_set(rop, ef::one().rop);
    }
    else
    {
      ::mpf_init_set(rop, ef::zero().rop);
    }
  }
  else
  {
    ::mpf_init_set(rop, ef::zero().rop);
    operator=(e_float(mantissa));
    operator*=(e_float("1E" + Util::lexical_cast(exponent)));
  }
}

gmp::e_float::e_float(const ::mpf_t& op) : fpclass  (ef_finite),
                                           prec_elem(ef_digits10_tol)
{
  init();
  ::mpf_init_set(rop, op);
}

gmp::e_float::~e_float()
{
  ::mpf_set_prec_raw(rop, static_cast<unsigned long int>(ef_digits2));
  ::mpf_clear       (rop);
}

void gmp::e_float::precision(const INT32 prec_digits)
{
  const unsigned long int digits2_request = static_cast<unsigned long int>(static_cast<UINT64>(static_cast<double>(prec_digits) / d_log2()));
  const unsigned long int d2              = static_cast<unsigned long int>(ef_digits2);
  const unsigned long int digits2_set     = std::min(digits2_request, d2);

  prec_elem = static_cast<INT32>(static_cast<INT64>(static_cast<double>(digits2_set) * d_log2()));

  ::mpf_set_prec_raw(rop, digits2_set);
}

gmp::e_float& gmp::e_float::operator=(const e_float& v)
{
  fpclass   = v.fpclass;
  prec_elem = v.prec_elem;

  ::mpf_set         (rop, v.rop);
  ::mpf_set_prec_raw(rop, static_cast<unsigned long>(static_cast<UINT64>(static_cast<double>(prec_elem) / d_log2())));

  return *this;
}

gmp::e_float& gmp::e_float::operator+=(const e_float& v)
{
  if(isnan())
  {
    return *this;
  }

  if(isinf())
  {
    if(v.isinf() && (isneg() != v.isneg()))
    {
      *this = std::numeric_limits<e_float>::quiet_NaN();
    }

    return *this;
  }

  ::mpf_add(rop, rop, v.rop);

  // Check for overflow.
  long u_exp2_signed;
  static_cast<void>(::mpf_get_d_2exp(&u_exp2_signed, rop));

  if(   (u_exp2_signed >= std::numeric_limits<e_float>::max_exponent)
     && (ef::fabs(*this) > std::numeric_limits<e_float>::max())
    )
  {
    const bool b_result_is_neg = isneg();

    *this = (!b_result_is_neg ?  std::numeric_limits<e_float>::infinity()
                              : -std::numeric_limits<e_float>::infinity());
  }

  return *this;
}

gmp::e_float& gmp::e_float::operator-=(const e_float& v)
{
  // Use *this - v = -(-*this + v).
  return (negate().operator+=(v)).negate();
}

gmp::e_float& gmp::e_float::operator*=(const e_float& v)
{
  const bool b_u_is_inf =   isinf();
  const bool b_v_is_inf = v.isinf();

  if(   (isnan() || v.isnan())
     || (b_u_is_inf && v.iszero())
     || (b_v_is_inf &&   iszero())
    )
  {
    return *this = std::numeric_limits<e_float>::quiet_NaN();
  }

  if(b_u_is_inf || b_v_is_inf)
  {
    const bool b_result_is_neg = (isneg() == v.isneg());

    *this = (!b_result_is_neg ?  std::numeric_limits<e_float>::infinity()
                              : -std::numeric_limits<e_float>::infinity());

    return *this;
  }

  // Get the base-2 exponent of *this and v and...
  long u_exp2_signed;
  long v_exp2_signed;
  static_cast<void>(::mpf_get_d_2exp(&u_exp2_signed,   rop));
  static_cast<void>(::mpf_get_d_2exp(&v_exp2_signed, v.rop));

  // Check for overflow or underflow.
  const bool u_exp2_is_neg = (u_exp2_signed < static_cast<long>(0));
  const bool v_exp2_is_neg = (v_exp2_signed < static_cast<long>(0));

  if(u_exp2_is_neg == v_exp2_is_neg)
  {
    // Get the unsigned base-2 exponents of *this and v and...
    const INT64 u_exp2 = !u_exp2_is_neg ? u_exp2_signed : -u_exp2_signed;
    const INT64 v_exp2 = !v_exp2_is_neg ? v_exp2_signed : -v_exp2_signed;

    // Check the range of the upcoming multiplication.
    const bool b_result_is_out_of_range = v_exp2 >= static_cast<long>(static_cast<long>(ef_max_exp) - u_exp2);

    if(b_result_is_out_of_range)
    {
      if(u_exp2_is_neg)
      {
        *this = ef::zero();
      }
      else
      {
        const bool b_result_is_neg = (isneg() == v.isneg());

        *this = (!b_result_is_neg ?  std::numeric_limits<e_float>::infinity()
                                  : -std::numeric_limits<e_float>::infinity());
      }

      return *this;
    }
  }

  // Multiply *this by v.
  ::mpf_mul(rop, rop, v.rop);

  return *this;
}

gmp::e_float& gmp::e_float::operator/=(const e_float& v)
{
  return operator*=(e_float(v).calculate_inv());
}

gmp::e_float& gmp::e_float::mul_by_int(const INT32 n)
{
  // Multiply *this with a constant signed integer.

  const bool b_n_is_neg = (n < static_cast<INT32>(0));

  const bool b_u_is_inf  = isinf();
  const bool b_n_is_zero = (n == static_cast<INT32>(0));

  if(isnan() || (b_u_is_inf && b_n_is_zero))
  {
    return *this = std::numeric_limits<e_float>::quiet_NaN();
  }

  if(b_u_is_inf)
  {
    const bool b_result_is_neg = (isneg() != b_n_is_neg);

    *this = (!b_result_is_neg ?  std::numeric_limits<e_float>::infinity()
                              : -std::numeric_limits<e_float>::infinity());

    return *this;
  }

  const unsigned long nn = static_cast<unsigned long>(!b_n_is_neg ? n : -n);

  ::mpf_mul_ui(rop, rop, static_cast<unsigned long>(nn));

  if(b_n_is_neg)
  {
    negate();
  }

  // Check for overflow.
  long u_exp2_signed;
  static_cast<void>(::mpf_get_d_2exp(&u_exp2_signed, rop));

  if(   (u_exp2_signed >= std::numeric_limits<e_float>::max_exponent)
     && (ef::fabs(*this) > std::numeric_limits<e_float>::max())
    )
  {
    const bool b_result_is_neg = (isneg() != b_n_is_neg);

    *this = (!b_result_is_neg ?  std::numeric_limits<e_float>::infinity()
                              : -std::numeric_limits<e_float>::infinity());
  }

  return *this;
}

gmp::e_float& gmp::e_float::div_by_int(const INT32 n)
{
  const bool b_n_is_neg = n < static_cast<INT32>(0);

  if(isnan())
  {
    return *this;
  }

  if(isinf())
  {
    const bool b_result_is_neg = (isneg() != b_n_is_neg);

    *this = (!b_result_is_neg ?  std::numeric_limits<e_float>::infinity()
                              : -std::numeric_limits<e_float>::infinity());

    return *this;
  }

  if(n == static_cast<INT32>(0))
  {
    // Divide by 0.
    if(iszero())
    {
      return *this = std::numeric_limits<e_float>::quiet_NaN();
    }
    else
    {
      *this = (!isneg() ?  std::numeric_limits<e_float>::infinity()
                        : -std::numeric_limits<e_float>::infinity());

      return *this;
    }
  }

  if(iszero())
  {
    return *this;
  }

  const unsigned long un = static_cast<unsigned long>(!b_n_is_neg ? n : -n);
  
  ::mpf_div_ui(rop, rop, un);

  if(b_n_is_neg)
  {
    negate();
  }

  // Check for underflow.
  long u_exp2_signed;
  static_cast<void>(::mpf_get_d_2exp(&u_exp2_signed, rop));

  if(   (u_exp2_signed <= std::numeric_limits<e_float>::min_exponent)
     && (ef::fabs(*this) < std::numeric_limits<e_float>::min())
    )
  {
    return *this = ef::zero();
  }

  return *this;
}

gmp::e_float& gmp::e_float::calculate_inv(void)
{
  // Compute the inverse of *this.

  bool b_result_is_neg = isneg();

  if(iszero())
  {
    *this = (!b_result_is_neg ?  std::numeric_limits<e_float>::infinity()
                              : -std::numeric_limits<e_float>::infinity());

    return *this;
  }

  if(isnan())
  {
    return *this;
  }

  if(isinf())
  {
    return *this = ef::zero();
  }

  if(isone())
  {
    *this = (!b_result_is_neg ? ef::one() : -ef::one());

    return *this;
  }

  ::mpf_ui_div(rop, static_cast<unsigned long int>(1u), rop);

  return *this;
}

gmp::e_float& gmp::e_float::calculate_sqrt(void)
{
  // Compute the square root of *this.

  if(isneg() || !isfinite())
  {
    return *this = std::numeric_limits<e_float>::quiet_NaN();
  }

  if(iszero() || isone())
  {
    return *this;
  }

  ::mpf_sqrt(rop, rop);

  return *this;
}

INT32 gmp::e_float::cmp_data(const ::mpf_t& v) const
{
  const INT32 result = static_cast<INT32>(::mpf_cmp(rop, v));
  
  if(result > static_cast<INT32>(0))
  {
    return static_cast<INT32>(1);
  }
  else if(result < static_cast<INT32>(0))
  {
    return static_cast<INT32>(-1);
  }
  else
  {
    return static_cast<INT32>(0);
  }
}

INT32 gmp::e_float::cmp(const e_float& v) const
{
  const INT32 result = cmp_data(v.rop);

  if(result == static_cast<INT32>(0))
  {
    return (isfinite() && v.isfinite()) ? static_cast<INT32>(0) : static_cast<INT32>(1);
  }
  else
  {
    return result;
  }
}

bool gmp::e_float::iszero(void) const
{
  // Check if the value of *this is identically 0 or very close to 0.
  return isint() && (cmp(ef::zero()) == static_cast<INT32>(0));
}

bool gmp::e_float::isone(void) const
{
  // Check if the value of *this is identically 1 or very close to 1.
  return isint() && (cmp(ef::one()) == static_cast<INT32>(0));
}

bool gmp::e_float::isint(void) const
{
  // Check if the value of *this is pure integer or very close to pure integer.
  return ::mpf_integer_p(rop) != 0;
}

bool gmp::e_float::isneg(void) const
{
  return ::mpf_sgn(rop) < 0;
}

const gmp::e_float& gmp::e_float::my_value_nan(void) const
{
  static e_float val(0u);

  val.fpclass = ef_NaN;

  static const e_float qnan(val);
  
  return qnan;
}

const gmp::e_float& gmp::e_float::my_value_inf(void) const
{
  static e_float val(0u);

  val.fpclass = ef_inf;

  static const e_float inf(val);

  return inf;
}

const gmp::e_float& gmp::e_float::my_value_max(void) const
{
  static const INT64 exp10_max = std::numeric_limits<e_float>::max_exponent10;

  static const e_float val("1E" + Util::lexical_cast(exp10_max));

  return val;
}

const gmp::e_float& gmp::e_float::my_value_min(void) const
{
  static const INT64 exp10_min = std::numeric_limits<e_float>::min_exponent10;

  static const e_float val("1E" + Util::lexical_cast(exp10_min));

  return val;
}

e_float& gmp::e_float::negate(void)
{
  ::mpf_neg(rop, rop);
  return *this;
}

e_float& gmp::e_float::operator++(void) { ::mpf_add_ui(rop, rop, static_cast<unsigned long>(1u)); return *this; }
e_float& gmp::e_float::operator--(void) { ::mpf_sub_ui(rop, rop, static_cast<unsigned long>(1u)); return *this; }

void gmp::e_float::extract_parts(double& mantissa, INT64& exponent) const
{
  const bool b_neg = isneg();

  long n2;
  const double d2    = ::mpf_get_d_2exp(&n2, (ef::fabs(*this)).rop);
  const double x_exp = static_cast<double>(static_cast<double>(n2) * d_log2());

  const double  x_exp_integer_part = static_cast<double>(static_cast<long>(x_exp));
  const double  x_exp_decimal_part = static_cast<double>(x_exp - x_exp_integer_part);

  double m = d2 * ::pow(static_cast<double>(10.0), x_exp_decimal_part);
  INT64  e = static_cast<INT64>(x_exp_integer_part);

  if(m < static_cast<double>(1.0))
  {
    m *= static_cast<double>(10.0);
    e -= static_cast<INT64>(1);
  }

  mantissa = static_cast<double>(!b_neg ? m : -m);
  exponent = e;
}

double gmp::e_float::extract_double(void) const
{
  const bool b_neg = isneg();

  const e_float xx = ef::fabs(*this);

  static const e_float dbl_max(std::numeric_limits<double>::max());
  static const e_float dbl_min(std::numeric_limits<double>::min());

  if(xx > dbl_max)
  {
    return !b_neg ?  std::numeric_limits<double>::max()
                  : -std::numeric_limits<double>::max();
  }
  else if(xx < dbl_min)
  {
    return !b_neg ?  std::numeric_limits<double>::min()
                  : -std::numeric_limits<double>::min();
  }

  const double dx = ::mpf_get_d(xx.rop);

  return !b_neg ? dx : -dx;
}

INT64 gmp::e_float::extract_int64(void) const
{
  const bool b_neg = isneg();

  // Make a rounded copy.
  e_float xr = *this;

  if(isint())
  {
    !b_neg ? xr += ef::half() : xr -= ef::half();
  }

  const e_float nx = ef::fabs(xr.extract_integer_part());

  static const e_float n64_max(std::numeric_limits<INT64>::max());

  if(nx > n64_max)
  {
    return !b_neg ?  std::numeric_limits<INT64>::max()
                  : -std::numeric_limits<INT64>::max();
  }
  
  if(nx < ef::one())
  {
    return static_cast<INT64>(0);
  }

  if(nx.isone())
  {
    return static_cast<INT64>(!b_neg ? 1 : -1);
  }

  static const char c0 = static_cast<char>('\0');

  std::vector<char> str(32u, c0);

  mp_exp_t p10;

  static_cast<void>(::mpf_get_str(&str[0], &p10, 10, str.size() - 1u, nx.rop));

  std::string str_n64(static_cast<std::size_t>(p10), static_cast<char>('0'));

  std::copy(str.begin(),
            std::find(str.begin(), str.end(), c0),
            str_n64.begin());

  std::stringstream ss;

  ss << str_n64;

  INT64 n;
  ss >> n;

  return !b_neg ? n : -n;
}

gmp::e_float gmp::e_float::extract_integer_part(void) const
{
  const bool b_neg = isneg();

  const e_float xx = ef::fabs(*this);
    
  e_float nx(xx);
  ::mpf_floor(nx.rop, xx.rop);

  return !b_neg ? nx : -nx;
}

gmp::e_float gmp::e_float::extract_decimal_part(void) const
{
  const bool b_neg = isneg();

  const e_float xx = ef::fabs(*this);

  const e_float dx = xx - xx.extract_integer_part();

  return !b_neg ? dx : -dx;
}

INT64 gmp::e_float::order(void) const
{
  const e_float xx = ef::fabs(*this);

  if(xx.iszero() || xx.isone())
  {
    return static_cast<INT64>(0);
  }
  else
  {
    signed long int n2;
    const double d2    = ::mpf_get_d_2exp(&n2, xx.rop);
    const double lg10x = static_cast<double>(::log10(d2) + (static_cast<double>(n2) * d_log2()));

    return lg10x < static_cast<double>(0.0) ? static_cast<INT64>(lg10x - static_cast<double>(0.5))
                                            : static_cast<INT64>(lg10x + static_cast<double>(0.5));
  }
}

void gmp::e_float::wr_string(std::string& str, std::ostream& os) const
{
  if(isnan())
  {
    str = "NaN";
    return;
  }

  if(isinf())
  {
    str = "INF";
    return;
  }

  static const std::streamsize p_min = static_cast<std::streamsize>(10);
  static const std::streamsize p_lim = static_cast<std::streamsize>(ef_digits10_tol);
         const std::streamsize p     = std::max(os.precision(), p_min);

  const std::streamsize my_precision = std::min(p, p_lim);

  const std::ios_base::fmtflags f = os.flags();

  const bool my_uppercase  = ((f & std::ios_base::uppercase)  != static_cast<std::ios_base::fmtflags>(0u));
  const bool my_showpos    = ((f & std::ios_base::showpos)    != static_cast<std::ios_base::fmtflags>(0u));
  const bool my_scientific = ((f & std::ios_base::scientific) != static_cast<std::ios_base::fmtflags>(0u));

  // Create a format string such as "%+.99Fe"
  const std::string str_fmt =   (my_showpos ? "%+." : "%.")
                              + Util::lexical_cast(my_precision - (my_scientific ? 1 : 0))
                              + (my_scientific ? (my_uppercase ? "FE" : "Fe") : (my_uppercase ? "FG" : "Fg"));

  std::tr1::array<char, static_cast<std::size_t>(e_float::ef_digits10_tol + 32)> buf = {{ static_cast<char>(0) }};

  static_cast<void>(gmp_sprintf(buf.data(), str_fmt.c_str(), rop));

  // Set the result string and remove the '\0' padding by using the c_str() representation.
  str = std::string(buf.data());

  const std::size_t pos_E = (my_uppercase ? str.rfind('E') : str.rfind('e'));

  if(pos_E != std::string::npos)
  {
    // Pad the exponent number field with additional zeros such that the width
    // of the exponent number field is equal to the width of ef_max_exp10.
    const std::size_t pos_exp = static_cast<std::string::size_type>(pos_E + 2u);

    const std::string::size_type width_of_exp = str.length() - pos_exp;

    str.insert(pos_exp, std::string(width_of_exponent_field() - width_of_exp, static_cast<char>('0')));
  }
}

namespace gmp
{
  static bool has_exp_or_has_dec_predicate(const char& c)
  {
    return    (c == static_cast<char>('e'))
           || (c == static_cast<char>('E'))
           || (c == static_cast<char>('.'));
  }
}

bool gmp::e_float::rd_string(const char* const s)
{
  std::string str(s);

  // Remove spaces and tabs
  static const char spc = static_cast<char>(' ');
  static const char tab = static_cast<char>('\t');

  str.erase(std::remove(str.begin(), str.end(), spc), str.end());
  str.erase(std::remove(str.begin(), str.end(), tab), str.end());

  // Get possible + sign and remove it
  
  if(   !str.empty()
     &&  str.at(static_cast<std::size_t>(0u)) == static_cast<char>('+')
    )
  {
    str.erase(static_cast<std::size_t>(0u),
              static_cast<std::size_t>(1u));
  }

  // Get possible - sign and remove it
  
  bool b_negate = false;

  if(   !str.empty()
     &&  str.at(static_cast<std::size_t>(0u)) == static_cast<char>('-')
    )
  {
    b_negate = true;
    str.erase(static_cast<std::size_t>(0u),
              static_cast<std::size_t>(1u));
  }

  // Remove leading zeros for all input types
  while(   !str.empty()
        &&  str.at(static_cast<std::size_t>(0u)) == static_cast<char>('0')
       )
  {
    str.erase(static_cast<std::size_t>(0u),
              static_cast<std::size_t>(1u));
  }

  // Scale very long pure integer input strings. Convert these into a string with
  // a decimal point and an exponent.

  const bool is_pure_integer = std::find_if(str.begin(),
                                            str.end(),
                                            gmp::has_exp_or_has_dec_predicate) == str.end();

  bool b_scaled = false;

  if(is_pure_integer && (str.length() > static_cast<std::size_t>(ef::tol())))
  {
    b_scaled = true;

    const std::size_t exp = static_cast<std::size_t>(str.length() - static_cast<std::size_t>(1u));

    const std::string str_exp = "E" + Util::lexical_cast(exp);

    str = str.substr(static_cast<std::size_t>(0u),
                     static_cast<std::size_t>(static_cast<std::size_t>(ef::tol()) - 1u));

    str.insert(static_cast<std::size_t>(1u), ".");

    str += str_exp;
  }

  // Set the e_float value.
  const INT32 n_set_result = static_cast<INT32>(::mpf_init_set_str(rop, str.c_str(), 10));

  if(b_negate)
  {
    negate();
  }

  return (n_set_result == static_cast<INT32>(0));
}
