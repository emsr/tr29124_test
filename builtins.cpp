#include <iostream>

template <typename _Tp>
bool __isnan(const _Tp __x)
{
  return __builtin_isnan(__x);
}

template <>
bool __isnan<float>(const float __x)
{
  return __builtin_isnanf(__x);
}

template <>
bool __isnan<double>(const double __x)
{
  return __builtin_isnan(__x);
}

template <>
bool __isnan<long double>(const long double __x)
{
  return __builtin_isnanl(__x);
}

//  Plain old overloads.
bool __oisnan(const float __x)
{
  return __builtin_isnanf(__x);
}
bool __oisnan(const double __x)
{
  return __builtin_isnan(__x);
}
bool __oisnan(const long double __x)
{
  return __builtin_isnanl(__x);
}


int main(int, char **)
{
  double x = 1.5;
  double y = __builtin_lgamma(x);
  if (__builtin_isnan(y))
    return 1;

  std::cout << __isnan(1.0F) << std::endl;
  std::cout << __isnan(1.0) << std::endl;
  std::cout << __isnan(1.0L) << std::endl;

  std::cout << __isnan(__builtin_nanf("0")) << std::endl;
  std::cout << __isnan(__builtin_nan("0")) << std::endl;
  std::cout << __isnan(__builtin_nanl("0")) << std::endl;

  std::cout << __oisnan(__builtin_nanf("0")) << std::endl;
  std::cout << __oisnan(__builtin_nan("0")) << std::endl;
  std::cout << __oisnan(__builtin_nanl("0")) << std::endl;

  return 0;
}

