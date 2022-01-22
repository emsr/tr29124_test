#ifndef TRAITS_UTILS_H
#define TRAITS_UTILS_H 

#include <type_traits>

namespace emsr
{

namespace detail
{
  template<typename, typename = std::void_t<>>
    struct has_value_type
    : std::false_type
    { };
 
  template<typename _Tp>
    struct has_value_type<_Tp, std::void_t<typename _Tp::value_type>>
    : std::true_type
    { };

  template<typename _Tp>
    inline constexpr bool
    has_value_type_v = has_value_type<_Tp>::value;
}

  template<typename _Tp, bool _HasValueType>
    struct value_type_impl;

  template<typename _Tp>
    struct value_type_impl<_Tp, false>
    {
      using type = std::decay_t<_Tp>;
    };

  template<typename _Tp>
    struct value_type_impl<_Tp, true>
    {
      using vtype = typename _Tp::value_type;
      using type = value_type_impl<vtype,
				       detail::has_value_type_v<vtype>>;
    };

  template<typename _Tp>
    struct value_type
    {
      using type = value_type_impl<_Tp, detail::has_value_type_v<_Tp>>;
    };

  template<typename _Tp>
    using value_t = typename value_type<_Tp>::type;

}

#endif // TRAITS_UTILS_H
