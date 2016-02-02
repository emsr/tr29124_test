#include <type_traits>
#include <complex>

  /**
   *  This is a more modern version of __promote_N in ext/type_traits.
   *  This is used for numeric argument promotion of complex and cmath.
   */
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
    struct __promote_help<float>
    { using __type = float; };

  template<>
    struct __promote_help<double>
    { using __type = double; };

  template<>
    struct __promote_help<long double>
    { using __type = long double; };

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    struct __promote_help<__float128>
    { using __type = __float128; };
#endif

  template<typename... _Tps>
    using __promote_help_t = typename __promote_help<_Tps...>::__type;

  // Decay refs and cv...
  // Alternatively we could decay refs and propagate cv to promoted type.
  template<typename _Tp, typename... _Tps>
    struct __promote_num
    { using __type = decltype(__promote_help_t<std::decay_t<_Tp>>{}
		   + typename __promote_num<_Tps...>::__type{}); };

  template<>
    template<typename _Tp>
      struct __promote_num<_Tp>
      { using __type = decltype(__promote_help_t<std::decay_t<_Tp>>{}); };

  template<typename... _Tps>
    using __promote_num_t = typename __promote_num<_Tps...>::__type;

  // Assume complex value_type is floating point.
  template<>
    template<typename _Tp>
      struct __promote_help<std::complex<_Tp>, false>
      {
      private:
	using __vtype = typename std::complex<_Tp>::value_type;
      public:
	using __type = decltype(std::complex<__promote_help_t<__vtype>>{});
      };

  // primary template handles types that have no nested ::value_type member:
  template<typename, typename = std::void_t<>>
    struct __has_value_type
    : std::false_type
    { };

  // Specialization recognizes types that do have a nested ::value_type member:
  template<typename _Tp>
    struct __has_value_type<_Tp, std::void_t<typename _Tp::value_type>>
    : std::true_type
    { };

  // Try to get something like array<Tp, Num> or vector<Tp, Alloc>
  // This should be able to do complex too.
  template<>
    template<template<typename _Arg, typename... _Args> typename _Tp>
      struct __promote_help<typename _Tp<_Arg, _Args...>,
			    typename = std::void_t<typename _Tp<_Arg, _Args...>::value_type>>
      {
      private:
	using __vtype = typename _Tp<_Arg, _Args...>::value_type;
	using __ptype = __promote_num_t<__vtype>;
      public:
	using __type = decltype(_Tp<__ptype, _Args...>>{});
      };
