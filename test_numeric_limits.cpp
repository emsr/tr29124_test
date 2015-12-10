// $HOME/bin/bin/g++ -o test_numeric_limits test_numeric_limits.cpp -lquadmath

// LD_LIBRARY_PATH=/$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_numeric_limits

#include "numeric_limits.h"
#include "float128.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <typeinfo>
#include <typeindex>
#include <experimental/string_view>

template<typename _Tp>
  void
  test(_Tp __x)
  {
    //  enum float_denorm_style
    //  {
    //    denorm_indeterminate  = -1,
    //    denorm_absent         = 0,
    //    denorm_present        = 1
    //  };
    std::map<int, std::experimental::string_view>
    denorm
    {
      {std::denorm_indeterminate, "denorm_indeterminate"},
      {std::denorm_absent,        "denorm_absent"},
      {std::denorm_present,       "denorm_present"}
    };

    std::map<std::float_round_style, std::experimental::string_view>
    round
    {
      {std::round_indeterminate, "round_indeterminate"},
      {std::round_toward_zero, "round_toward_zero"},
      {std::round_to_nearest, "round_to_nearest"},
      {std::round_toward_infinity, "round_toward_infinity"},
      {std::round_toward_neg_infinity, "round_toward_neg_infinity"}
    };

    auto name{std::type_index{typeid(__x)}.name()};

    std::cout << std::setiosflags(std::ios::boolalpha);
    std::cout << std::setprecision(__gnu_cxx::digits10(__x));

    std::cout << '\n';
    std::cout << "type             : " << std::quoted(name) << '\n';
    std::cout << "is_specialized   : " << __gnu_cxx::is_specialized(__x) << '\n';
    std::cout << "min              : " << __gnu_cxx::min(__x) << '\n';
    std::cout << "max              : " << __gnu_cxx::max(__x) << '\n';
    std::cout << "lowest           : " << __gnu_cxx::lowest(__x) << '\n';
    std::cout << "digits           : " << __gnu_cxx::digits(__x) << '\n';
    std::cout << "digits10         : " << __gnu_cxx::digits10(__x) << '\n';
    std::cout << "max_digits10     : " << __gnu_cxx::max_digits10(__x) << '\n';
    std::cout << "is_signed        : " << __gnu_cxx::xxx_is_signed(__x) << '\n';
    std::cout << "is_integer       : " << __gnu_cxx::is_integer(__x) << '\n';
    std::cout << "is_exact         : " << __gnu_cxx::is_exact(__x) << '\n';
    std::cout << "radix            : " << __gnu_cxx::radix(__x) << '\n';
    std::cout << "epsilon          : " << __gnu_cxx::epsilon(__x) << '\n';
    std::cout << "round_error      : " << __gnu_cxx::round_error(__x) << '\n';
    std::cout << "min_exponent     : " << __gnu_cxx::min_exponent(__x) << '\n';
    std::cout << "min_exponent10   : " << __gnu_cxx::min_exponent10(__x) << '\n';
    std::cout << "max_exponent     : " << __gnu_cxx::max_exponent(__x) << '\n';
    std::cout << "max_exponent10   : " << __gnu_cxx::max_exponent10(__x) << '\n';
    std::cout << "has_infinity     : " << __gnu_cxx::has_infinity(__x) << '\n';
    std::cout << "has_quiet_NaN    : " << __gnu_cxx::has_quiet_NaN(__x) << '\n';
    std::cout << "has_signaling_NaN: " << __gnu_cxx::has_signaling_NaN(__x) << '\n';
    std::cout << "has_denorm       : " << denorm[__gnu_cxx::has_denorm(__x)] << '\n';
    std::cout << "has_denorm_loss  : " << __gnu_cxx::has_denorm_loss(__x) << '\n';
    std::cout << "infinity         : " << __gnu_cxx::infinity(__x) << '\n';
    std::cout << "quiet_NaN        : " << __gnu_cxx::quiet_NaN(__x) << '\n';
    std::cout << "signaling_NaN    : " << __gnu_cxx::signaling_NaN(__x) << '\n';
    std::cout << "denorm_min       : " << __gnu_cxx::denorm_min(__x) << '\n';
    std::cout << "is_iec559        : " << __gnu_cxx::is_iec559(__x) << '\n';
    std::cout << "is_bounded       : " << __gnu_cxx::is_bounded(__x) << '\n';
    std::cout << "is_modulo        : " << __gnu_cxx::is_modulo(__x) << '\n';
    std::cout << "traps            : " << __gnu_cxx::traps(__x) << '\n';
    std::cout << "tinyness_before  : " << __gnu_cxx::tinyness_before(__x) << '\n';
    std::cout << "round_style      : " << round[__gnu_cxx::round_style(__x)] << '\n';
    std::cout << '\n';
    std::cout << "sqrt_max         : " << __gnu_cxx::sqrt_max<_Tp>(__x) << '\n';
    std::cout << "cbrt_max         : " << __gnu_cxx::cbrt_max<_Tp>(__x) << '\n';
    std::cout << "root_max(5)      : " << __gnu_cxx::root_max(_Tp{5}) << '\n';
    std::cout << "log_max          : " << __gnu_cxx::log_max<_Tp>(__x) << '\n';
    std::cout << "log10_max        : " << __gnu_cxx::log10_max<_Tp>(__x) << '\n';
    std::cout << "sqrt_min         : " << __gnu_cxx::sqrt_min<_Tp>(__x) << '\n';
    std::cout << "cbrt_min         : " << __gnu_cxx::cbrt_min<_Tp>(__x) << '\n';
    std::cout << "root_min(5)      : " << __gnu_cxx::root_min(_Tp{5}) << '\n';
    std::cout << "log_min          : " << __gnu_cxx::log_min<_Tp>(__x) << '\n';
    std::cout << "log10_min        : " << __gnu_cxx::log10_min<_Tp>(__x) << '\n';
    std::cout << "sqrt_eps         : " << __gnu_cxx::sqrt_eps<_Tp>(__x) << '\n';
    std::cout << "cbrt_eps         : " << __gnu_cxx::cbrt_eps<_Tp>(__x) << '\n';
    std::cout << "root_eps(5)      : " << __gnu_cxx::root_eps(_Tp{5}) << '\n';
    std::cout << "log_eps          : " << __gnu_cxx::log_eps<_Tp>(__x) << '\n';
    std::cout << "log10_eps        : " << __gnu_cxx::log10_eps<_Tp>(__x) << '\n';
  }

int
main()
{
  test(1.0F);
  test(1.0);
  test(1.0L);
  test(1.0Q);

  const volatile auto x = 6.66666F;
  test(x);
}
