/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_limits test_limits.cpp src/c++11/limits.cc -lquadmath
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
    std::cout << "  sizeof              " << sizeof(Numeric) << std::endl;
    std::cout << "  alignof             " << alignof(Numeric) << std::endl;

    std::cout << "  denorm_style        " << float_denorm_style[std::numeric_limits<Numeric>::has_denorm] << std::endl;
    std::cout << "  has_denorm_loss     " << std::numeric_limits<Numeric>::has_denorm_loss << std::endl;
    std::cout << "  has_infinity        " << std::numeric_limits<Numeric>::has_infinity << std::endl;
    std::cout << "  has_quiet_NaN       " << std::numeric_limits<Numeric>::has_quiet_NaN << std::endl;
    std::cout << "  has_signaling_NaN   " << std::numeric_limits<Numeric>::has_signaling_NaN << std::endl;
    std::cout << "  is_bounded          " << std::numeric_limits<Numeric>::is_bounded << std::endl;
    std::cout << "  is_exact            " << std::numeric_limits<Numeric>::is_exact << std::endl;
    std::cout << "  is_iec559           " << std::numeric_limits<Numeric>::is_iec559 << std::endl;
    std::cout << "  is_integer          " << std::numeric_limits<Numeric>::is_integer << std::endl;
    std::cout << "  is_modulo           " << std::numeric_limits<Numeric>::is_modulo << std::endl;
    std::cout << "  is_signed           " << std::numeric_limits<Numeric>::is_signed << std::endl;
    std::cout << "  is_specialized      " << std::numeric_limits<Numeric>::is_specialized << std::endl;
    std::cout << "  tinyness_before     " << std::numeric_limits<Numeric>::tinyness_before << std::endl;
    std::cout << "  traps               " << std::numeric_limits<Numeric>::traps << std::endl;
    std::cout << "  round_style         " << float_round_style[std::numeric_limits<Numeric>::round_style] << std::endl;
    std::cout << "  digits              " << std::numeric_limits<Numeric>::digits << std::endl;
    std::cout << "  digits10            " << std::numeric_limits<Numeric>::digits10 << std::endl;
    std::cout << "  max_digits10        " << std::numeric_limits<Numeric>::max_digits10 << std::endl;
    std::cout << "  max_exponent        " << std::numeric_limits<Numeric>::max_exponent << std::endl;
    std::cout << "  max_exponent10      " << std::numeric_limits<Numeric>::max_exponent10 << std::endl;
    std::cout << "  min_exponent        " << std::numeric_limits<Numeric>::min_exponent << std::endl;
    std::cout << "  min_exponent10      " << std::numeric_limits<Numeric>::min_exponent10 << std::endl;
    std::cout << "  radix               " << std::numeric_limits<Numeric>::radix << std::endl;

    std::cout << "  denorm_min          " << std::numeric_limits<Numeric>::denorm_min() << std::endl;
    std::cout << "  epsilon             " << std::numeric_limits<Numeric>::epsilon() << std::endl;
    std::cout << "  infinity            " << std::numeric_limits<Numeric>::infinity() << std::endl;
    std::cout << "  max                 " << std::numeric_limits<Numeric>::max() << std::endl;
    std::cout << "  min                 " << std::numeric_limits<Numeric>::min() << std::endl;
    std::cout << "  lowest              " << std::numeric_limits<Numeric>::lowest() << std::endl;
    std::cout << "  quiet_NaN           " << std::numeric_limits<Numeric>::quiet_NaN() << std::endl;
    std::cout << "  round_error         " << std::numeric_limits<Numeric>::round_error() << std::endl;
    std::cout << "  signaling_NaN       " << std::numeric_limits<Numeric>::signaling_NaN() << std::endl;

    std::cout.precision(prec);
    std::cout.flags(flags);

    return;
  }
 
int
main()
{
  std::cout << std::endl << "char" << std::endl;
  test_limits<char>();

  std::cout << std::endl << "signed char" << std::endl;
  test_limits<signed char>();

  std::cout << std::endl << "unsigned char" << std::endl;
  test_limits<unsigned char>();

  std::cout << std::endl << "wchar_t" << std::endl;
  test_limits<wchar_t>();

  std::cout << std::endl << "char16_t" << std::endl;
  test_limits<char16_t>();

  std::cout << std::endl << "char32_t" << std::endl;
  test_limits<char32_t>();

  std::cout << std::endl << "short" << std::endl;
  test_limits<short>();

  std::cout << std::endl << "signed short" << std::endl;
  test_limits<signed short>();

  std::cout << std::endl << "unsigned short" << std::endl;
  test_limits<unsigned short>();

  std::cout << std::endl << "int" << std::endl;
  test_limits<int>();

  std::cout << std::endl << "signed int" << std::endl;
  test_limits<signed int>();

  std::cout << std::endl << "unsigned int" << std::endl;
  test_limits<unsigned int>();

  std::cout << std::endl << "long" << std::endl;
  test_limits<long>();

  std::cout << std::endl << "signed long" << std::endl;
  test_limits<signed long>();

  std::cout << std::endl << "unsigned long" << std::endl;
  test_limits<unsigned long>();

  std::cout << std::endl << "long long" << std::endl;
  test_limits<long long>();

  std::cout << std::endl << "signed long long" << std::endl;
  test_limits<signed long long>();

  std::cout << std::endl << "unsigned long long" << std::endl;
  test_limits<unsigned long long>();

  std::cout << std::endl << "float" << std::endl;
  test_limits<float>();

  std::cout << std::endl << "double" << std::endl;
  test_limits<double>();

  std::cout << std::endl << "long double" << std::endl;
  test_limits<long double>();

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << std::endl << "__float128" << std::endl;
  test_limits<__float128>();
#endif

  return 0;
}
