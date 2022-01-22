/**
 *
 */

#include <mpreal.h>

#include <emsr/numeric_limits.h>
#include <ext/float128_io.h>
#include <emsr/numeric_limits_mpreal.h>

#include <iostream>
#include <iomanip>
#include <map>
#include <typeinfo>
#include <typeindex>
#include <string>

template<typename _Tp>
  void
  test_numeric_limits(_Tp x)
  {
    std::map<std::float_denorm_style, std::string>
    denorm
    {
      {std::denorm_indeterminate, "denorm_indeterminate"},
      {std::denorm_absent,        "denorm_absent"},
      {std::denorm_present,       "denorm_present"}
    };

    std::map<std::float_round_style, std::string>
    round
    {
      {std::round_indeterminate, "round_indeterminate"},
      {std::round_toward_zero, "round_toward_zero"},
      {std::round_to_nearest, "round_to_nearest"},
      {std::round_toward_infinity, "round_toward_infinity"},
      {std::round_toward_neg_infinity, "round_toward_neg_infinity"}
    };

    auto name{std::type_index{typeid(x)}.name()};

    std::cout << std::setiosflags(std::ios::boolalpha);
    std::cout << std::setprecision(emsr::digits10(x));

    std::cout << '\n';
    std::cout << "type             : " << std::quoted(name) << '\n';
    std::cout << "is_specialized   : " << emsr::is_specialized(x) << '\n';
    std::cout << "min              : " << emsr::lim_min(x) << '\n';
    std::cout << "max              : " << emsr::lim_max(x) << '\n';
    std::cout << "lowest           : " << emsr::lowest(x) << '\n';
    std::cout << "digits           : " << emsr::digits(x) << '\n';
    std::cout << "digits10         : " << emsr::digits10(x) << '\n';
    std::cout << "max_digits10     : " << emsr::max_digits10(x) << '\n';
    std::cout << "is_signed        : " << emsr::is_signed(x) << '\n';
    std::cout << "is_integer       : " << emsr::is_integer(x) << '\n';
    std::cout << "is_exact         : " << emsr::is_exact(x) << '\n';
    std::cout << "radix            : " << emsr::radix(x) << '\n';
    std::cout << "epsilon          : " << emsr::epsilon(x) << '\n';
    std::cout << "round_error      : " << emsr::round_error(x) << '\n';
    std::cout << "min_exponent     : " << emsr::min_exponent(x) << '\n';
    std::cout << "min_exponent10   : " << emsr::min_exponent10(x) << '\n';
    std::cout << "max_exponent     : " << emsr::max_exponent(x) << '\n';
    std::cout << "max_exponent10   : " << emsr::max_exponent10(x) << '\n';
    std::cout << "has_infinity     : " << emsr::has_infinity(x) << '\n';
    std::cout << "has_quiet_NaN    : " << emsr::has_quiet_NaN(x) << '\n';
    std::cout << "has_signaling_NaN: " << emsr::has_signaling_NaN(x) << '\n';
    std::cout << "has_denorm       : " << denorm[emsr::has_denorm(x)] << '\n';
    std::cout << "has_denorm_loss  : " << emsr::has_denorm_loss(x) << '\n';
    std::cout << "infinity         : " << emsr::infinity(x) << '\n';
    std::cout << "quiet_NaN        : " << emsr::quiet_NaN(x) << '\n';
    std::cout << "signaling_NaN    : " << emsr::signaling_NaN(x) << '\n';
    std::cout << "denorm_min       : " << emsr::denorm_min(x) << '\n';
    std::cout << "is_iec559        : " << emsr::is_iec559(x) << '\n';
    std::cout << "is_bounded       : " << emsr::is_bounded(x) << '\n';
    std::cout << "is_modulo        : " << emsr::is_modulo(x) << '\n';
    std::cout << "traps            : " << emsr::traps(x) << '\n';
    std::cout << "tinyness_before  : " << emsr::tinyness_before(x) << '\n';
    std::cout << "round_style      : " << round[emsr::round_style(x)] << '\n';
    std::cout << '\n';
    std::cout << "max_integer      : " << emsr::max_integer<_Tp>(x) << '\n';
    std::cout << "sqrt_max         : " << emsr::sqrt_max<_Tp>(x) << '\n';
    std::cout << "cbrt_max         : " << emsr::cbrt_max<_Tp>(x) << '\n';
    std::cout << "root_max(5)      : " << emsr::root_max(_Tp{5}) << '\n';
    std::cout << "log_max          : " << emsr::log_max<_Tp>(x) << '\n';
    std::cout << "log10_max        : " << emsr::log10_max<_Tp>(x) << '\n';
    std::cout << "sqrt_min         : " << emsr::sqrt_min<_Tp>(x) << '\n';
    std::cout << "cbrt_min         : " << emsr::cbrt_min<_Tp>(x) << '\n';
    std::cout << "root_min(5)      : " << emsr::root_min(_Tp{5}) << '\n';
    std::cout << "log_min          : " << emsr::log_min<_Tp>(x) << '\n';
    std::cout << "log10_min        : " << emsr::log10_min<_Tp>(x) << '\n';
    std::cout << "sqrt_eps         : " << emsr::sqrt_eps<_Tp>(x) << '\n';
    std::cout << "cbrt_eps         : " << emsr::cbrt_eps<_Tp>(x) << '\n';
    std::cout << "root_eps(5)      : " << emsr::root_eps(_Tp{5}) << '\n';
    std::cout << "log_eps          : " << emsr::log_eps<_Tp>(x) << '\n';
    std::cout << "log10_eps        : " << emsr::log10_eps<_Tp>(x) << '\n';
  }

int
main()
{
  test_numeric_limits(1.0F);
  test_numeric_limits(1.0);
  test_numeric_limits(1.0L);
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //test_numeric_limits(1.0Q);
#endif

  const volatile auto x = 6.66666F;
  test_numeric_limits(x);

  long double y = 123.465L;
  mpfr::mpreal b(y, 256);
  test_numeric_limits(b);
}
