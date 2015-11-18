#include <type_traits>

  // For complex and cmath
  template<typename _Tp, bool = std::is_integral<_Tp>::value>
    struct __promote_help
    { using __type = double; };

  // No nested __type member for non-integer non-floating point types,
  // allows this type to be used for SFINAE to constrain overloads in
  // <cmath> and <complex> to only the intended types.
  template<typename _Tp>
    struct __promote_help<_Tp, false>
    { };

  template<>
    struct __promote_help<long double>
    { using __type = long double; };

  template<>
    struct __promote_help<double>
    { using __type = double; };

  template<>
    struct __promote_help<float>
    { using __type = float; };

  template<typename... _Tps>
    using __promote_help_t = typename __promote_help<_Tps...>::__type;

  // Decay refs and cv...
  // Alternatively we could decay refs and propagate cv to promoted type.
  template<typename _Tp, typename... _Tps>
    struct __promote_num
    { using __type = decltype(__promote_help_t<std::decay_t<_Tp>>()
		   + typename __promote_num<_Tps...>::__type()); };

  template<>
    template<typename _Tp>
      struct __promote_num<_Tp>
      { using __type = decltype(__promote_help_t<std::decay_t<_Tp>>()); };

  template<typename... _Tps>
    using __promote_num_t = typename __promote_num<_Tps...>::__type;
