
// g++ -std=c++11 -o test_pow test_pow.cpp

#include <iostream>
#include <limits>
#include <ext/cmath>

template<typename _Tp>
  constexpr _Tp
  sqrt_max()
  { return ldexp(_Tp{0.5}, std::numeric_limits<_Tp>::max_exponent / 2); }

template<typename _Tp>
  constexpr _Tp
  cbrt_max()
  { return ldexp(_Tp{0.5}, std::numeric_limits<_Tp>::max_exponent / 3); }

template<typename _Tp>
  constexpr _Tp
  root_max(_Tp __root)
  { return ldexp(_Tp{0.5}, std::numeric_limits<_Tp>::max_exponent / __root); }

template<typename _Tp>
  constexpr _Tp
  log_max()
  { return ldexp(_Tp{0.5}, log2(std::numeric_limits<_Tp>::max_exponent * __gnu_cxx::__math_constants<_Tp>::__ln_2)); }


template<typename _Tp>
  constexpr _Tp
  sqrt_min()
  { return ldexp(_Tp{2}, std::numeric_limits<_Tp>::min_exponent / 2); }

template<typename _Tp>
  constexpr _Tp
  cbrt_min()
  { return ldexp(_Tp{2}, std::numeric_limits<_Tp>::min_exponent / 3); }

template<typename _Tp>
  constexpr _Tp
  root_min(_Tp __root)
  { return ldexp(_Tp{2}, std::numeric_limits<_Tp>::min_exponent / __root); }

template<typename _Tp>
  constexpr _Tp
  log_min()
  { return ldexp(_Tp{2}, -log2(-std::numeric_limits<_Tp>::min_exponent * __gnu_cxx::__math_constants<_Tp>::__ln_2)); }


template<typename _Tp>
  constexpr _Tp
  sqrt_eps()
  { return ldexp(_Tp{1}, -std::numeric_limits<_Tp>::digits / 2); }

template<typename _Tp>
  constexpr _Tp
  cbrt_eps()
  { return ldexp(_Tp{1}, -std::numeric_limits<_Tp>::digits / 3); }

template<typename _Tp>
  constexpr _Tp
  root_eps(_Tp __root)
  { return ldexp(_Tp{1}, -std::numeric_limits<_Tp>::digits / __root); }

template<typename _Tp>
  constexpr _Tp
  log_eps()
  { return ldexp(_Tp{1}, -log2(std::numeric_limits<_Tp>::digits * __gnu_cxx::__math_constants<_Tp>::__ln_2)); }


template<typename _Tp>
  void
  test_limits()
  {
    std::cout << "  tinyness_before     " << std::numeric_limits<_Tp>::tinyness_before << std::endl;
    std::cout << "  digits              " << std::numeric_limits<_Tp>::digits << std::endl;
    std::cout << "  digits10            " << std::numeric_limits<_Tp>::digits10 << std::endl;
    std::cout << "  max_digits10        " << std::numeric_limits<_Tp>::max_digits10 << std::endl;
    std::cout << "  max_exponent        " << std::numeric_limits<_Tp>::max_exponent << std::endl;
    std::cout << "  max_exponent10      " << std::numeric_limits<_Tp>::max_exponent10 << std::endl;
    std::cout << "  min_exponent        " << std::numeric_limits<_Tp>::min_exponent << std::endl;
    std::cout << "  min_exponent10      " << std::numeric_limits<_Tp>::min_exponent10 << std::endl;
    std::cout << "  denorm_min          " << std::numeric_limits<_Tp>::denorm_min() << std::endl;
    std::cout << "  epsilon             " << std::numeric_limits<_Tp>::epsilon() << std::endl;
    std::cout << "  max                 " << std::numeric_limits<_Tp>::max() << std::endl;
    std::cout << "  min                 " << std::numeric_limits<_Tp>::min() << std::endl;
    std::cout << "  sqrt_max            " << sqrt_max<_Tp>() << std::endl;
    std::cout << "  cbrt_max            " << cbrt_max<_Tp>() << std::endl;
    std::cout << "  root_max(5)         " << root_max(_Tp{5}) << std::endl;
    std::cout << "  log_max             " << log_max<_Tp>() << std::endl;
    std::cout << "  sqrt_min            " << sqrt_min<_Tp>() << std::endl;
    std::cout << "  cbrt_min            " << cbrt_min<_Tp>() << std::endl;
    std::cout << "  root_min(5)         " << root_min(_Tp{5}) << std::endl;
    std::cout << "  log_min             " << log_min<_Tp>() << std::endl;
    std::cout << "  sqrt_eps            " << sqrt_eps<_Tp>() << std::endl;
    std::cout << "  cbrt_eps            " << cbrt_eps<_Tp>() << std::endl;
    std::cout << "  root_eps(5)         " << root_eps(_Tp{5}) << std::endl;
    std::cout << "  log_eps             " << log_eps<_Tp>() << std::endl;
  }

int
main()
{
  std::cout << std::endl << "float" << std::endl;
  test_limits<float>();

  std::cout << std::endl << "double" << std::endl;
  test_limits<double>();

  std::cout << std::endl << "long double" << std::endl;
  test_limits<long double>();
}
