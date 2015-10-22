#include <iostream>
#include <limits>


template<typename Numeric>
void test_limits( void )
{

  std::ios_base::fmtflags prev = std::cout.flags( std::ios_base::boolalpha );

  std::cout << "  denorm_style        " << std::numeric_limits<Numeric>::has_denorm << std::endl;
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
  std::cout << "  round_style         " << std::numeric_limits<Numeric>::round_style << std::endl;
  std::cout << "  digits              " << std::numeric_limits<Numeric>::digits << std::endl;
  std::cout << "  digits10            " << std::numeric_limits<Numeric>::digits10 << std::endl;
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
  std::cout << "  quiet_NaN           " << std::numeric_limits<Numeric>::quiet_NaN() << std::endl;
  std::cout << "  round_error         " << std::numeric_limits<Numeric>::round_error() << std::endl;
  std::cout << "  signaling_NaN       " << std::numeric_limits<Numeric>::signaling_NaN() << std::endl;

  std::cout.flags( prev );

  return;
}
 
int main( int, char ** )
{

  std::cout << std::endl << "char" << std::endl;
  test_limits<char>();

  std::cout << std::endl << "short" << std::endl;
  test_limits<short>();

  std::cout << std::endl << "int" << std::endl;
  test_limits<int>();

  std::cout << std::endl << "long" << std::endl;
  test_limits<long>();

  std::cout << std::endl << "long long" << std::endl;
  test_limits<long long>();

  std::cout << std::endl << "float" << std::endl;
  test_limits<float>();

  std::cout << std::endl << "double" << std::endl;
  test_limits<double>();

  std::cout << std::endl << "long double" << std::endl;
  test_limits<long double>();

  return 0;
}
