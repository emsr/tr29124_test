
#include <sstream>
#include <iomanip>

#include <e_float/e_float.h>
#include <e_float/f90/e_float_f90_protos.h>
#include <utility/util_lexical_cast.h>
#include <utility/util_numeric_cast.h>

const double& f90::e_float::d_log2(void)
{
  static const double value_log2 = static_cast<double>(0.3010299956639811952137389);
  return value_log2;
}

f90::e_float::e_float()
{
  static const double dz = static_cast<double>(0.0);
  ::QXSETDBL(&data, &dz);
}

f90::e_float::e_float(const INT32 n)
{
  ::QXSETN32(&data, &n);
}

f90::e_float::e_float(const INT64 n)
{
  ::QXSETN64(&data, &n);
}

f90::e_float::e_float(const UINT32 u)
{
  const INT64 n = static_cast<INT64>(u);
  ::QXSETN64(&data, &n);
}

f90::e_float::e_float(const UINT64 u)
{
  if(u > static_cast<UINT64>(std::numeric_limits<INT64>::max()))
  {
    std::stringstream ss;
    ss << u;
    static_cast<void>(rd_string(ss.str().c_str()));
  }
  else
  {
    const INT64 n = static_cast<INT64>(u);
    ::QXSETN64(&data, &n);
  }
}

f90::e_float::e_float(const double d)
{
  ::QXSETDBL(&data, &d);
}

f90::e_float::e_float(const long double ld)
{
  std::stringstream ss;

  ss << std::scientific
     << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<long double>::digits10 - 1))
     << ld;

  static_cast<void>(rd_string(ss.str().c_str()));
}

f90::e_float::e_float(const e_float& f)
{
  data = f.data;
}

f90::e_float::e_float(const char* s)
{
  static_cast<void>(rd_string(s));
}

f90::e_float::e_float(const std::string& str)
{
  static_cast<void>(rd_string(str.c_str()));
}

f90::e_float::e_float(const double mantissa, const INT64 exponent)
{
  const bool mantissa_is_iszero = (::fabs(mantissa) < (std::numeric_limits<double>::min() * static_cast<double>(2.0)));

  if(mantissa_is_iszero)
  {
    if(exponent == static_cast<INT64>(0))
    {
      static const INT32 n_one = static_cast<INT32>(1);
      ::QXSETN32(&data, &n_one);
    }
    else
    {
      static const INT32 n_zero = static_cast<INT32>(0);
      ::QXSETN32(&data, &n_zero);
    }
  }
  else
  {
    ::QXSETDBL(&data, &mantissa);
    const e_float ef_exp("1E" + Util::lexical_cast(exponent));
    operator*=(ef_exp);
  }
}

const f90::e_float& f90::e_float::my_value_nan(void) const
{
  static const e_float val = e_float(-1).calculate_sqrt();
  return val;
}

const f90::e_float& f90::e_float::my_value_inf(void) const
{
  static const e_float val = e_float(0).calculate_inv();
  return val;
}

const f90::e_float& f90::e_float::my_value_max(void) const
{
  static bool b_is_init = false;

  static e_float my_max;

  if(!b_is_init)
  {
    b_is_init = true;
    ::QXHUGE(&my_max.data, &my_max.data);
  }

  static const e_float val = my_max;

  return val;
}

const f90::e_float& f90::e_float::my_value_min(void) const
{
  static bool b_is_init = false;

  static e_float my_min;

  if(!b_is_init)
  {
    b_is_init = true;
    ::QXTINY(&my_min.data, &my_min.data);
  }

  static const e_float val = my_min;

  return val;
}

f90::e_float& f90::e_float::operator+=(const e_float& v)
{
  ::QXADD(&data, &v.data);
  return *this;
}

f90::e_float& f90::e_float::operator-=(const e_float& v)
{
  ::QXSUB(&data, &v.data);
  return *this;
}

f90::e_float& f90::e_float::operator*=(const e_float& v)
{
  ::QXMUL(&data, &v.data);
  return *this;
}

f90::e_float& f90::e_float::operator/=(const e_float& v)
{
  ::QXDIV(&data, &v.data);
  return *this;
}

f90::e_float& f90::e_float::mul_by_int(const INT32 n)
{
  ::QXMULN(&data, &n);
  return *this;
}

f90::e_float& f90::e_float::div_by_int(const INT32 n)
{
  ::QXDIVN(&data, &n);
  return *this;
}

f90::e_float& f90::e_float::negate(void)
{
  ::QXNEGATE(&data);
  return *this;
}

INT32 f90::e_float::cmp(const e_float& v) const
{
  INT32 n;
  ::QXCMP(&data, &v.data, &n);
  return n;
}

f90::e_float& f90::e_float::calculate_sqrt(void)
{
  ::QXSQRT(&data);
  return *this;
}

f90::e_float& f90::e_float::calculate_inv(void)
{
  ::QXINV(&data);
  return *this;
}

bool f90::e_float::isnan   (void) const { INT32 n; ::QXISNAN   (&data, &n); return (n == static_cast<INT32>(1)); }
bool f90::e_float::isinf   (void) const { INT32 n; ::QXISINF   (&data, &n); return (n == static_cast<INT32>(1)); }
bool f90::e_float::isfinite(void) const { INT32 n; ::QXISFINITE(&data, &n); return (n == static_cast<INT32>(1)); }
bool f90::e_float::iszero  (void) const { return (cmp(ef::zero()) == static_cast<INT32>(0)); }
bool f90::e_float::isone   (void) const { return (cmp(ef::one())  == static_cast<INT32>(0)); }
bool f90::e_float::isneg   (void) const { INT32 n; ::QXISNEG(&data, &n); return (n == static_cast<INT32>(1)); }

bool f90::e_float::isint(void) const
{
  static const e_float delta = std::numeric_limits<e_float>::min() * static_cast<INT32>(2);

  const e_float xx = ef::fabs(*this);

  e_float the_floor;
  ::QXFLOOR(&the_floor.data, &xx.data);

  e_float the_ceil;
  ::QXCEIL(&the_ceil.data, &xx.data);

  if((xx - the_floor) < delta)
  {
    return true;
  }
  else if((the_ceil - xx) < delta)
  {
    return true;
  }
  else
  {
    return false;
  }
}

f90::e_float& f90::e_float::operator++(void) { ::QXADD(&data, &ef::one().data); return *this; }
f90::e_float& f90::e_float::operator--(void) { ::QXSUB(&data, &ef::one().data); return *this; }

void f90::e_float::extract_parts(double& mantissa, INT64& exponent) const
{
  const bool b_neg = isneg();

  const e_float xx = ef::fabs(*this);

  INT32  n2;
  double d2;
  ::QXPARTS2(&xx.data, &d2, &n2);

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

double f90::e_float::extract_double(void) const
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
  else
  {
    double dx;
    ::QX2DBL(&data, &dx);
    return dx;
  }
}

INT64 f90::e_float::extract_int64(void) const
{
  const bool b_neg = isneg();

  const e_float nx = ef::fabs(*this);

  static const e_float n64_max(std::numeric_limits<INT64>::max());

  if(nx > n64_max)
  {
    return (!b_neg ?  std::numeric_limits<INT64>::max()
                   : -std::numeric_limits<INT64>::max());
  }
  else if(nx < ef::one())
  {
    return static_cast<INT64>(0);
  }
  else
  {
    INT64 n;
    ::QX2N64(&data, &n);
    return n;
  }
}

f90::e_float f90::e_float::extract_integer_part(void) const
{
  const bool b_neg = isneg();

  const e_float xx = ef::fabs(*this);

  e_float nx;
  ::QXFLOOR(&nx.data, &xx.data);

  return (!b_neg ? nx : -nx);
}

f90::e_float f90::e_float::extract_decimal_part(void) const
{
  const bool b_neg = isneg();

  const e_float xx = ef::fabs(*this);

  const e_float dx = xx - xx.extract_integer_part();

  return (!b_neg ? dx : -dx);
}

INT64 f90::e_float::order(void) const
{
  INT32 n;
  ::QXORDER2(&data, &n);
  return static_cast<INT32>(static_cast<double>(n) * d_log2());
}

void f90::e_float::wr_string(std::string& str, std::ostream& os) const
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

  // TBD: Support for output formats other than scientific notation needs to be implemented.
  const std::ios_base::fmtflags f = os.flags();

  const bool my_uppercase  = ((f & std::ios_base::uppercase)  != static_cast<std::ios_base::fmtflags>(0u));
  const bool my_showpos    = ((f & std::ios_base::showpos)    != static_cast<std::ios_base::fmtflags>(0u));
  const bool my_scientific = ((f & std::ios_base::scientific) != static_cast<std::ios_base::fmtflags>(0u));

  std::tr1::array<char, 64u> s = {{ static_cast<char>('\0') }};

  ::QX2STR(&data, s.data());

  str.assign(s.begin(), s.begin() + 56u);
}

bool f90::e_float::rd_string(const char* const s)
{
  // Put the input string into the standard e_float input A.BBBE+/-000n, where A
  // has 1 digit, the exponent field is signed and four digits wide, and BBBB
  // fills 39 decimal places.

  std::string str(s);

  // Get possible exponent and remove it
  INT64 exp = static_cast<INT64>(0);

  std::size_t pos;

  if(   ((pos = str.find('e')) != std::string::npos)
     || ((pos = str.find('E')) != std::string::npos)
    )
  {
    exp = Util::numeric_cast<INT64>(std::string(str.c_str() + (pos + 1u)));

    str.erase(pos, str.length() - pos);
  }

  // Get possible +/- sign and remove it
  bool neg = false;

  if((pos = str.find(static_cast<char>('-'))) != std::string::npos)
  {
    neg = true;
    str.erase(pos, static_cast<std::size_t>(1u));
  }

  if((pos = str.find(static_cast<char>('+'))) != std::string::npos)
  {
    str.erase(pos, static_cast<std::size_t>(1u));
  }

  // Remove leading zeros for all input types
  const std::string::iterator fwd_it_leading_zero = std::find_if(str.begin(), str.end(), char_is_nonzero_predicate);

  if(fwd_it_leading_zero != str.begin())
  {
    if(fwd_it_leading_zero == str.end())
    {
      // The string contains nothing but leading zeros. The string is zero.
      operator=(ef::zero());
      return true;
    }
    else
    {
      str.erase(str.begin(), fwd_it_leading_zero);
    }
  }

  // Find possible decimal point
  pos = str.find(static_cast<char>('.'));

  if(pos != std::string::npos)
  {
    // Remove all trailing insignificant zeros.
    const std::string::const_reverse_iterator rev_it_non_zero_elem =
          std::find_if(str.rbegin(), str.rend(), char_is_nonzero_predicate);

    if(rev_it_non_zero_elem != str.rbegin())
    {
      const std::string::size_type ofs = str.length() - std::distance<std::string::const_reverse_iterator>(str.rbegin(), rev_it_non_zero_elem);
      str.erase(str.begin() + ofs, str.end());
    }

    // Check if input is identically zero
    if(str == std::string("."))
    {
      operator=(ef::zero());
      return true;
    }
  
    // Remove leading significant zeros just after the decimal point
    // and adjust the exponent accordingly.
    // Note that the while-loop operates only on strings of the form ".000abcd..."
    // and peels away the zeros just after the decimal point.
    if(str.at(static_cast<std::size_t>(0u)) == static_cast<char>('.'))
    {
      const std::string::iterator fwd_it_first_non_zero_decimal =
            std::find_if(str.begin() + 1u, str.end(), char_is_nonzero_predicate);

      std::string::size_type delta_exp = static_cast<std::size_t>(0u);

      if(str.at(static_cast<std::size_t>(1u)) == static_cast<char>('0'))
      {
        delta_exp = std::distance<std::string::const_iterator>(str.begin() + 1u,
                                                               fwd_it_first_non_zero_decimal);
      }

      // Bring one single digit into the mantissa and adjust exponent accordingly
      str.erase(str.begin(), fwd_it_first_non_zero_decimal);
      str.insert(static_cast<std::size_t>(1u), ".");
      exp -= static_cast<INT64>(delta_exp + 1u);
    }
  }
  else
  {
    // Input string has no decimal point: Append decimal point.
    str.append(".");
  }

  // Shift the decimal point if necessary.
  pos = str.find(static_cast<char>('.'));

  if(pos != static_cast<std::size_t>(1u))
  {
    str.erase(pos, 1u);
    str.insert(1u, ".");
    exp += static_cast<INT64>(pos - 1u);
  }

  if(exp > ef_max_exp10)
  {
    *this = (!neg ?  std::numeric_limits<e_float>::infinity()
                  : -std::numeric_limits<e_float>::infinity());

    return true;
  }
  else if(exp < ef_min_exp10)
  {
    *this = ef::zero();

    return true;
  }

  // Truncate the decimal part if it is too long.
  const std::size_t max_dec = static_cast<std::size_t>(41u);

  if(str.length() > max_dec)
  {
    str.erase(max_dec, max_dec - str.length());
  }
  else if(str.length() < max_dec)
  {
    str.insert(str.end(), max_dec - str.length(), static_cast<char>('0'));
  }

  // Create the exponent string.
  std::stringstream ss;
  ss << std::setw(4)
     << std::setfill(static_cast<char>('0'))
     << ((exp < static_cast<INT64>(0)) ? static_cast<INT64>(-exp) : exp);

  str += ((exp < static_cast<INT64>(0) ? "E-" : "E+") + ss.str());

  // Set the data.
  ::QXSETSTR(&data, str.c_str());

  // Negate if necessary.
  if(neg)
  {
    ::QXNEGATE(&data);
  }

  return true;
}

e_float f90::e_float::my_cbrt(const e_float& x) { e_float res; ::QXCBRT(&res.data, &x.data); return res; }
e_float f90::e_float::my_rootn(const e_float& x, const UINT32 p)
{
  INT32 pp = static_cast<INT32>(p);
  e_float res;
  ::QXROOTN(&res.data, &x.data, &pp);
  return res;
}
e_float f90::e_float::my_exp  (const e_float& x) { e_float res; ::QXEXP  (&res.data, &x.data); return res; }
e_float f90::e_float::my_log  (const e_float& x) { e_float res; ::QXLOG  (&res.data, &x.data); return res; }
e_float f90::e_float::my_sin  (const e_float& x) { e_float res; ::QXSIN  (&res.data, &x.data); return res; }
e_float f90::e_float::my_cos  (const e_float& x) { e_float res; ::QXCOS  (&res.data, &x.data); return res; }
e_float f90::e_float::my_tan  (const e_float& x) { e_float res; ::QXTAN  (&res.data, &x.data); return res; }
e_float f90::e_float::my_asin (const e_float& x) { e_float res; ::QXASIN (&res.data, &x.data); return res; }
e_float f90::e_float::my_acos (const e_float& x) { e_float res; ::QXACOS (&res.data, &x.data); return res; }
e_float f90::e_float::my_atan (const e_float& x) { e_float res; ::QXATAN (&res.data, &x.data); return res; }
e_float f90::e_float::my_sinh (const e_float& x) { e_float res; ::QXSINH (&res.data, &x.data); return res; }
e_float f90::e_float::my_cosh (const e_float& x) { e_float res; ::QXCOSH (&res.data, &x.data); return res; }
e_float f90::e_float::my_tanh (const e_float& x) { e_float res; ::QXTANH (&res.data, &x.data); return res; }
e_float f90::e_float::my_asinh(const e_float& x) { e_float res; ::QXASINH(&res.data, &x.data); return res; }
e_float f90::e_float::my_acosh(const e_float& x) { e_float res; ::QXACOSH(&res.data, &x.data); return res; }
e_float f90::e_float::my_atanh(const e_float& x) { e_float res; ::QXATANH(&res.data, &x.data); return res; }
