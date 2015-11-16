#include <type_traits>

  // For complex and cmath
  template<typename _Tp, bool = std::is_integral<_Tp>::value>
    struct __promote_help
    { typedef double __type; };

  // No nested __type member for non-integer non-floating point types,
  // allows this type to be used for SFINAE to constrain overloads in
  // <cmath> and <complex> to only the intended types.
  template<typename _Tp>
    struct __promote_help<_Tp, false>
    { };

  template<>
    struct __promote_help<long double>
    { typedef long double __type; };

  template<>
    struct __promote_help<double>
    { typedef double __type; };

  template<>
    struct __promote_help<float>
    { typedef float __type; };

  template<typename _Tp, typename... _Tps>
    struct __promote_num
    { using __type = decltype(_Tp() + typename __promote_help<_Tps...>::__type()); };

  template<typename... _Tps>
    using __promote_num_t = typename __promote_num<_Tps...>::__type;
