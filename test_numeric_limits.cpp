// $HOME/bin/bin/g++ -std=gnu++14 -I. -o test_numeric_limits test_numeric_limits.cpp -lquadmath -lmpfr

// LD_LIBRARY_PATH=/$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_numeric_limits > test_numeric_limits.txt

#include <mpreal.h>

#include <bits/numeric_limits.h>
#include <bits/float128.h>
#include <bits/numeric_limits_mpreal.h>

#include <iostream>
#include <iomanip>
#include <map>
#include <typeinfo>
#include <typeindex>
#include <string>

template<typename _Tp>
  void
  test(_Tp __x)
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

    auto name{std::type_index{typeid(__x)}.name()};

    std::cout << std::setiosflags(std::ios::boolalpha);
    std::cout << std::setprecision(__gnu_cxx::__digits10(__x));

    std::cout << '\n';
    std::cout << "type             : " << std::quoted(name) << '\n';
    std::cout << "is_specialized   : " << __gnu_cxx::__is_specialized(__x) << '\n';
    std::cout << "min              : " << __gnu_cxx::__min(__x) << '\n';
    std::cout << "max              : " << __gnu_cxx::__max(__x) << '\n';
    std::cout << "lowest           : " << __gnu_cxx::__lowest(__x) << '\n';
    std::cout << "digits           : " << __gnu_cxx::__digits(__x) << '\n';
    std::cout << "digits10         : " << __gnu_cxx::__digits10(__x) << '\n';
    std::cout << "max_digits10     : " << __gnu_cxx::__max_digits10(__x) << '\n';
    std::cout << "is_signed        : " << __gnu_cxx::__is_signed(__x) << '\n';
    std::cout << "is_integer       : " << __gnu_cxx::__is_integer(__x) << '\n';
    std::cout << "is_exact         : " << __gnu_cxx::__is_exact(__x) << '\n';
    std::cout << "radix            : " << __gnu_cxx::__radix(__x) << '\n';
    std::cout << "epsilon          : " << __gnu_cxx::__epsilon(__x) << '\n';
    std::cout << "round_error      : " << __gnu_cxx::__round_error(__x) << '\n';
    std::cout << "min_exponent     : " << __gnu_cxx::__min_exponent(__x) << '\n';
    std::cout << "min_exponent10   : " << __gnu_cxx::__min_exponent10(__x) << '\n';
    std::cout << "max_exponent     : " << __gnu_cxx::__max_exponent(__x) << '\n';
    std::cout << "max_exponent10   : " << __gnu_cxx::__max_exponent10(__x) << '\n';
    std::cout << "has_infinity     : " << __gnu_cxx::__has_infinity(__x) << '\n';
    std::cout << "has_quiet_NaN    : " << __gnu_cxx::__has_quiet_NaN(__x) << '\n';
    std::cout << "has_signaling_NaN: " << __gnu_cxx::__has_signaling_NaN(__x) << '\n';
    std::cout << "has_denorm       : " << denorm[__gnu_cxx::__has_denorm(__x)] << '\n';
    std::cout << "has_denorm_loss  : " << __gnu_cxx::__has_denorm_loss(__x) << '\n';
    std::cout << "infinity         : " << __gnu_cxx::__infinity(__x) << '\n';
    std::cout << "quiet_NaN        : " << __gnu_cxx::__quiet_NaN(__x) << '\n';
    std::cout << "signaling_NaN    : " << __gnu_cxx::__signaling_NaN(__x) << '\n';
    std::cout << "denorm_min       : " << __gnu_cxx::__denorm_min(__x) << '\n';
    std::cout << "is_iec559        : " << __gnu_cxx::__is_iec559(__x) << '\n';
    std::cout << "is_bounded       : " << __gnu_cxx::__is_bounded(__x) << '\n';
    std::cout << "is_modulo        : " << __gnu_cxx::__is_modulo(__x) << '\n';
    std::cout << "traps            : " << __gnu_cxx::__traps(__x) << '\n';
    std::cout << "tinyness_before  : " << __gnu_cxx::__tinyness_before(__x) << '\n';
    std::cout << "round_style      : " << round[__gnu_cxx::__round_style(__x)] << '\n';
    std::cout << '\n';
    std::cout << "sqrt_max         : " << __gnu_cxx::__sqrt_max<_Tp>(__x) << '\n';
    std::cout << "cbrt_max         : " << __gnu_cxx::__cbrt_max<_Tp>(__x) << '\n';
    std::cout << "root_max(5)      : " << __gnu_cxx::__root_max(_Tp{5}) << '\n';
    std::cout << "log_max          : " << __gnu_cxx::__log_max<_Tp>(__x) << '\n';
    std::cout << "log10_max        : " << __gnu_cxx::__log10_max<_Tp>(__x) << '\n';
    std::cout << "sqrt_min         : " << __gnu_cxx::__sqrt_min<_Tp>(__x) << '\n';
    std::cout << "cbrt_min         : " << __gnu_cxx::__cbrt_min<_Tp>(__x) << '\n';
    std::cout << "root_min(5)      : " << __gnu_cxx::__root_min(_Tp{5}) << '\n';
    std::cout << "log_min          : " << __gnu_cxx::__log_min<_Tp>(__x) << '\n';
    std::cout << "log10_min        : " << __gnu_cxx::__log10_min<_Tp>(__x) << '\n';
    std::cout << "sqrt_eps         : " << __gnu_cxx::__sqrt_eps<_Tp>(__x) << '\n';
    std::cout << "cbrt_eps         : " << __gnu_cxx::__cbrt_eps<_Tp>(__x) << '\n';
    std::cout << "root_eps(5)      : " << __gnu_cxx::__root_eps(_Tp{5}) << '\n';
    std::cout << "log_eps          : " << __gnu_cxx::__log_eps<_Tp>(__x) << '\n';
    std::cout << "log10_eps        : " << __gnu_cxx::__log10_eps<_Tp>(__x) << '\n';
  }

int
main()
{
  test(1.0F);
  test(1.0);
  test(1.0L);
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  test(1.0Q);
#endif

  const volatile auto x = 6.66666F;
  test(x);

  long double y = 123.465L;
  mpfr::mpreal b(y, 256);
  test(b);
}
