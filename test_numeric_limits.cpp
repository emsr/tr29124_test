// $HOME/bin/bin/g++ -o test_numeric_limits test_numeric_limits.cpp -lquadmath

// LD_LIBRARY_PATH=/$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_numeric_limits

#include "numeric_limits.h"
#include "float128.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <typeinfo>
#include <typeindex>

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
    std::map<int, std::string>
    denorm
    {
      {std::denorm_indeterminate, "denorm_indeterminate"},
      {std::denorm_absent,        "denorm_absent"},
      {std::denorm_present,       "denorm_present"}
    };

    auto name{std::type_index{typeid(__x)}.name()};

    std::cout << std::setiosflags(std::ios::boolalpha);
    std::cout << std::setprecision(std::digits10(__x));

    std::cout << '\n';
    std::cout << "type             : " << std::quoted(name) << '\n';
    std::cout << "is_specialized   : " << std::is_specialized(__x) << '\n';
    std::cout << "min              : " << std::min(__x) << '\n';
    std::cout << "max              : " << std::max(__x) << '\n';
    std::cout << "lowest           : " << std::lowest(__x) << '\n';
    std::cout << "digits           : " << std::digits(__x) << '\n';
    std::cout << "digits10         : " << std::digits10(__x) << '\n';
    std::cout << "max_digits10     : " << std::max_digits10(__x) << '\n';
    std::cout << "is_signed        : " << std::xxx_is_signed(__x) << '\n';
    std::cout << "is_integer       : " << std::is_integer(__x) << '\n';
    std::cout << "is_exact         : " << std::is_exact(__x) << '\n';
    std::cout << "radix            : " << std::radix(__x) << '\n';
    std::cout << "epsilon          : " << std::epsilon(__x) << '\n';
    std::cout << "round_error      : " << std::round_error(__x) << '\n';
    std::cout << "min_exponent     : " << std::min_exponent(__x) << '\n';
    std::cout << "min_exponent10   : " << std::min_exponent10(__x) << '\n';
    std::cout << "max_exponent     : " << std::max_exponent(__x) << '\n';
    std::cout << "max_exponent10   : " << std::max_exponent10(__x) << '\n';
    std::cout << "has_infinity     : " << std::has_infinity(__x) << '\n';
    std::cout << "has_quiet_NaN    : " << std::has_quiet_NaN(__x) << '\n';
    std::cout << "has_signaling_NaN: " << std::has_signaling_NaN(__x) << '\n';
    std::cout << "has_denorm       : " << denorm[std::has_denorm(__x)] << '\n';
    std::cout << "has_denorm_loss  : " << std::has_denorm_loss(__x) << '\n';
    std::cout << "infinity         : " << std::infinity(__x) << '\n';
    std::cout << "quiet_NaN        : " << std::quiet_NaN(__x) << '\n';
    std::cout << "signaling_NaN    : " << std::signaling_NaN(__x) << '\n';
    std::cout << "denorm_min       : " << std::denorm_min(__x) << '\n';
    std::cout << "is_iec559        : " << std::is_iec559(__x) << '\n';
    std::cout << "is_bounded       : " << std::is_bounded(__x) << '\n';
    std::cout << "is_modulo        : " << std::is_modulo(__x) << '\n';
    std::cout << "traps            : " << std::traps(__x) << '\n';
    std::cout << "tinyness_before  : " << std::tinyness_before(__x) << '\n';
    std::cout << "round_style      : " << std::round_style(__x) << '\n';
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
