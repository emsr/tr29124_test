#ifndef TRAITS_UTILS_H
#define TRAITS_UTILS_H 

#include <type_traits>

namespace __gnu_cxx
{

namespace __detail
{
  template<typename, typename = std::void_t<>>
    struct __has_value_type
    : std::false_type
    { };
 
  template<typename _Tp>
    struct __has_value_type<_Tp, std::void_t<typename _Tp::value_type>>
    : std::true_type
    { };

  template<typename _Tp>
    inline constexpr bool
    __has_value_type_v = __has_value_type<_Tp>::value;
}

  template<typename _Tp, bool _HasValueType>
    struct __value_type_impl;

  template<typename _Tp>
    struct __value_type_impl<_Tp, false>
    {
      using __type = std::decay_t<_Tp>;
    };

  template<typename _Tp>
    struct __value_type_impl<_Tp, true>
    {
      using __vtype = typename _Tp::value_type;
      using __type = __value_type_impl<__vtype,
				       __detail::__has_value_type_v<__vtype>>;
    };

  template<typename _Tp>
    struct __value_type
    {
      using __type = __value_type_impl<_Tp, __detail::__has_value_type_v<_Tp>>;
    };

  template<typename _Tp>
    using __value_t = typename __value_type<_Tp>::__type;

}

#endif // TRAITS_UTILS_H
