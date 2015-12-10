
// g++ -std=c++11 -o test_pow_funs test_pow.cpp

// ./test_pow_funs > ./test_pow_funs.txt

// g++ -DBITS -std=c++11 -o test_pow_bits test_pow.cpp

// ./test_pow_bits > ./test_pow_bits.txt

#include <iostream>
#include <limits>
#include <ext/cmath>

#if BITS

template<typename _Tp>
  constexpr _Tp
  sqrt_max()
  { return ldexp(_Tp{1}, std::numeric_limits<_Tp>::max_exponent / 2); }

template<typename _Tp>
  constexpr _Tp
  cbrt_max()
  { return ldexp(_Tp{1}, std::numeric_limits<_Tp>::max_exponent / 3); }

template<typename _Tp>
  constexpr _Tp
  root_max(_Tp __root)
  { return ldexp(_Tp{1}, std::numeric_limits<_Tp>::max_exponent / __root); }

template<typename _Tp>
  constexpr _Tp
  log_max()
  { return log2(std::numeric_limits<_Tp>::max_exponent * __gnu_cxx::__math_constants<_Tp>::__ln_2); }

template<typename _Tp>
  constexpr _Tp
  log10_max()
  { return log2(std::numeric_limits<_Tp>::max_exponent10 * __gnu_cxx::__math_constants<_Tp>::__ln_2); }


template<typename _Tp>
  constexpr _Tp
  sqrt_min()
  { return ldexp(_Tp{1}, std::numeric_limits<_Tp>::min_exponent / 2); }

template<typename _Tp>
  constexpr _Tp
  cbrt_min()
  { return ldexp(_Tp{1}, std::numeric_limits<_Tp>::min_exponent / 3); }

template<typename _Tp>
  constexpr _Tp
  root_min(_Tp __root)
  { return ldexp(_Tp{1}, std::numeric_limits<_Tp>::min_exponent / __root); }

template<typename _Tp>
  constexpr _Tp
  log_min()
  { return -log2(-std::numeric_limits<_Tp>::min_exponent * __gnu_cxx::__math_constants<_Tp>::__ln_2); }

template<typename _Tp>
  constexpr _Tp
  log10_min()
  { return -log2(-std::numeric_limits<_Tp>::min_exponent10 * __gnu_cxx::__math_constants<_Tp>::__ln_2); }


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
  { return -log2(std::numeric_limits<_Tp>::digits * __gnu_cxx::__math_constants<_Tp>::__ln_2); }

template<typename _Tp>
  constexpr _Tp
  log10_eps()
  { return -log2(std::numeric_limits<_Tp>::digits10 * __gnu_cxx::__math_constants<_Tp>::__ln_2); }

#else

template<typename _Tp>
  constexpr _Tp
  sqrt_max()
  { return std::sqrt(std::numeric_limits<_Tp>::max()); }

template<typename _Tp>
  constexpr _Tp
  cbrt_max()
  { return std::cbrt(std::numeric_limits<_Tp>::max()); }

template<typename _Tp>
  constexpr _Tp
  root_max(_Tp __root)
  { return std::pow(std::numeric_limits<_Tp>::max(), 1 / __root); }

template<typename _Tp>
  constexpr _Tp
  log_max()
  { return std::log(std::numeric_limits<_Tp>::max()); }

template<typename _Tp>
  constexpr _Tp
  log10_max()
  { return std::log10(std::numeric_limits<_Tp>::max()); }


template<typename _Tp>
  constexpr _Tp
  sqrt_min()
  { return std::sqrt(std::numeric_limits<_Tp>::min()); }

template<typename _Tp>
  constexpr _Tp
  cbrt_min()
  { return std::cbrt(std::numeric_limits<_Tp>::min()); }

template<typename _Tp>
  constexpr _Tp
  root_min(_Tp __root)
  { return std::pow(std::numeric_limits<_Tp>::min(), 1 / __root); }

template<typename _Tp>
  constexpr _Tp
  log_min()
  { return std::log(std::numeric_limits<_Tp>::min()); }

template<typename _Tp>
  constexpr _Tp
  log10_min()
  { return std::log10(std::numeric_limits<_Tp>::min()); }


template<typename _Tp>
  constexpr _Tp
  sqrt_eps()
  { return std::sqrt(std::numeric_limits<_Tp>::epsilon()); }

template<typename _Tp>
  constexpr _Tp
  cbrt_eps()
  { return std::cbrt(std::numeric_limits<_Tp>::epsilon()); }

template<typename _Tp>
  constexpr _Tp
  root_eps(_Tp __root)
  { return std::pow(std::numeric_limits<_Tp>::epsilon(), 1 / __root); }

template<typename _Tp>
  constexpr _Tp
  log_eps()
  { return std::log(std::numeric_limits<_Tp>::epsilon()); }

template<typename _Tp>
  constexpr _Tp
  log10_eps()
  { return std::log10(std::numeric_limits<_Tp>::epsilon()); }

#endif

template<typename _Tp>
  void
  test_limits()
  {
    std::cout << "  tinyness_before     " << std::numeric_limits<_Tp>::tinyness_before << '\n';
    std::cout << "  digits              " << std::numeric_limits<_Tp>::digits << '\n';
    std::cout << "  digits10            " << std::numeric_limits<_Tp>::digits10 << '\n';
    std::cout << "  max_digits10        " << std::numeric_limits<_Tp>::max_digits10 << '\n';
    std::cout << "  max_exponent        " << std::numeric_limits<_Tp>::max_exponent << '\n';
    std::cout << "  max_exponent10      " << std::numeric_limits<_Tp>::max_exponent10 << '\n';
    std::cout << "  min_exponent        " << std::numeric_limits<_Tp>::min_exponent << '\n';
    std::cout << "  min_exponent10      " << std::numeric_limits<_Tp>::min_exponent10 << '\n';
    std::cout << "  denorm_min          " << std::numeric_limits<_Tp>::denorm_min() << '\n';
    std::cout << "  epsilon             " << std::numeric_limits<_Tp>::epsilon() << '\n';
    std::cout << "  max                 " << std::numeric_limits<_Tp>::max() << '\n';
    std::cout << "  min                 " << std::numeric_limits<_Tp>::min() << '\n';
    std::cout << "  sqrt_max            " << sqrt_max<_Tp>() << '\n';
    std::cout << "  cbrt_max            " << cbrt_max<_Tp>() << '\n';
    std::cout << "  root_max(5)         " << root_max(_Tp{5}) << '\n';
    std::cout << "  log_max             " << log_max<_Tp>() << '\n';
    std::cout << "  log10_max           " << log10_max<_Tp>() << '\n';
    std::cout << "  sqrt_min            " << sqrt_min<_Tp>() << '\n';
    std::cout << "  cbrt_min            " << cbrt_min<_Tp>() << '\n';
    std::cout << "  root_min(5)         " << root_min(_Tp{5}) << '\n';
    std::cout << "  log_min             " << log_min<_Tp>() << '\n';
    std::cout << "  log10_min           " << log10_min<_Tp>() << '\n';
    std::cout << "  sqrt_eps            " << sqrt_eps<_Tp>() << '\n';
    std::cout << "  cbrt_eps            " << cbrt_eps<_Tp>() << '\n';
    std::cout << "  root_eps(5)         " << root_eps(_Tp{5}) << '\n';
    std::cout << "  log_eps             " << log_eps<_Tp>() << '\n';
    std::cout << "  log10_eps           " << log10_eps<_Tp>() << '\n';
  }

int
main()
{
  std::cout << std::endl << "float" << '\n';
  test_limits<float>();

  std::cout << std::endl << "double" << '\n';
  test_limits<double>();

  std::cout << std::endl << "long double" << '\n';
  test_limits<long double>();
}
