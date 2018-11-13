/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_limits test_limits.cpp src/c++11/limits.cc -lquadmath
./test_limits > test_limits.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_limits test_limits.cpp -lquadmath
./test_limits > test_limits.txt
*/

#include <iostream>
#include <limits>
#include <string_view>
#include <map>
#include <bits/float128_io.h>


template<typename Numeric>
  void
  test_limits(Numeric proto [[maybe_unused]] = Numeric{})
  {
    std::map<std::float_denorm_style, std::string_view>
    float_denorm_style
    {
      {std::denorm_indeterminate, "denorm_indeterminate"},
      {std::denorm_absent, "denorm_absent"},
      {std::denorm_present, "denorm_present"}
    };

    std::map<std::float_round_style, std::string_view>
    float_round_style
    {
      {std::round_indeterminate, "round_indeterminate"},
      {std::round_toward_zero, "round_toward_zero"},
      {std::round_to_nearest, "round_to_nearest"},
      {std::round_toward_infinity, "round_toward_infinity"},
      {std::round_toward_neg_infinity, "round_toward_neg_infinity"}
    };

    auto flags = std::cout.flags(std::ios_base::boolalpha);
    auto prec = std::cout.precision(std::numeric_limits<Numeric>::max_digits10);
    std::cout << "  sizeof              " << sizeof(Numeric) << '\n';
    std::cout << "  alignof             " << alignof(Numeric) << '\n';

    std::cout << "  denorm_style        " << float_denorm_style[std::numeric_limits<Numeric>::has_denorm] << '\n';
    std::cout << "  has_denorm_loss     " << std::numeric_limits<Numeric>::has_denorm_loss << '\n';
    std::cout << "  has_infinity        " << std::numeric_limits<Numeric>::has_infinity << '\n';
    std::cout << "  has_quiet_NaN       " << std::numeric_limits<Numeric>::has_quiet_NaN << '\n';
    std::cout << "  has_signaling_NaN   " << std::numeric_limits<Numeric>::has_signaling_NaN << '\n';
    std::cout << "  is_bounded          " << std::numeric_limits<Numeric>::is_bounded << '\n';
    std::cout << "  is_exact            " << std::numeric_limits<Numeric>::is_exact << '\n';
    std::cout << "  is_iec559           " << std::numeric_limits<Numeric>::is_iec559 << '\n';
    std::cout << "  is_integer          " << std::numeric_limits<Numeric>::is_integer << '\n';
    std::cout << "  is_modulo           " << std::numeric_limits<Numeric>::is_modulo << '\n';
    std::cout << "  is_signed           " << std::numeric_limits<Numeric>::is_signed << '\n';
    std::cout << "  is_specialized      " << std::numeric_limits<Numeric>::is_specialized << '\n';
    std::cout << "  tinyness_before     " << std::numeric_limits<Numeric>::tinyness_before << '\n';
    std::cout << "  traps               " << std::numeric_limits<Numeric>::traps << '\n';
    std::cout << "  round_style         " << float_round_style[std::numeric_limits<Numeric>::round_style] << '\n';
    std::cout << "  digits              " << std::numeric_limits<Numeric>::digits << '\n';
    std::cout << "  digits10            " << std::numeric_limits<Numeric>::digits10 << '\n';
    std::cout << "  max_digits10        " << std::numeric_limits<Numeric>::max_digits10 << '\n';
    std::cout << "  max_exponent        " << std::numeric_limits<Numeric>::max_exponent << '\n';
    std::cout << "  max_exponent10      " << std::numeric_limits<Numeric>::max_exponent10 << '\n';
    std::cout << "  min_exponent        " << std::numeric_limits<Numeric>::min_exponent << '\n';
    std::cout << "  min_exponent10      " << std::numeric_limits<Numeric>::min_exponent10 << '\n';
    std::cout << "  radix               " << std::numeric_limits<Numeric>::radix << '\n';

    std::cout << "  denorm_min          " << std::numeric_limits<Numeric>::denorm_min() << '\n';
    std::cout << "  epsilon             " << std::numeric_limits<Numeric>::epsilon() << '\n';
    std::cout << "  infinity            " << std::numeric_limits<Numeric>::infinity() << '\n';
    std::cout << "  max                 " << std::numeric_limits<Numeric>::max() << '\n';
    std::cout << "  min                 " << std::numeric_limits<Numeric>::min() << '\n';
    std::cout << "  lowest              " << std::numeric_limits<Numeric>::lowest() << '\n';
    std::cout << "  quiet_NaN           " << std::numeric_limits<Numeric>::quiet_NaN() << '\n';
    std::cout << "  round_error         " << std::numeric_limits<Numeric>::round_error() << '\n';
    std::cout << "  signaling_NaN       " << std::numeric_limits<Numeric>::signaling_NaN() << '\n';

    std::cout.precision(prec);
    std::cout.flags(flags);

    return;
  }
 
int
main()
{
  std::cout << "\nchar\n" << std::flush;
  test_limits<char>();

  std::cout << "\nsigned char\n" << std::flush;
  test_limits<signed char>();

  std::cout << "\nunsigned char\n" << std::flush;
  test_limits<unsigned char>();

  std::cout << "\nwchar_t\n" << std::flush;
  test_limits<wchar_t>();

  std::cout << "\nchar16_t\n" << std::flush;
  test_limits<char16_t>();

  std::cout << "\nchar32_t\n" << std::flush;
  test_limits<char32_t>();

  std::cout << "\nshort\n" << std::flush;
  test_limits<short>();

  std::cout << "\nsigned short\n" << std::flush;
  test_limits<signed short>();

  std::cout << "\nunsigned short\n" << std::flush;
  test_limits<unsigned short>();

  std::cout << "\nint\n" << std::flush;
  test_limits<int>();

  std::cout << "\nsigned int\n" << std::flush;
  test_limits<signed int>();

  std::cout << "\nunsigned int\n" << std::flush;
  test_limits<unsigned int>();

  std::cout << "\nlong\n" << std::flush;
  test_limits<long>();

  std::cout << "\nsigned long\n" << std::flush;
  test_limits<signed long>();

  std::cout << "\nunsigned long\n" << std::flush;
  test_limits<unsigned long>();

  std::cout << "\nlong long\n" << std::flush;
  test_limits<long long>();

  std::cout << "\nsigned long long\n" << std::flush;
  test_limits<signed long long>();

  std::cout << "\nunsigned long long\n" << std::flush;
  test_limits<unsigned long long>();

  std::cout << "\nfloat\n" << std::flush;
  test_limits<float>();

  std::cout << "\ndouble\n" << std::flush;
  test_limits<double>();

  std::cout << "\nlong double\n" << std::flush;
  test_limits<long double>();

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\n__float128\n" << std::flush;
  test_limits<__float128>();
#endif

  return 0;
}
