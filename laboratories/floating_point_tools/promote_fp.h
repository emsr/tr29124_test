#include <type_traits>
#include <complex>

  /**
   *  This is a more modern version of promote_N in ext/type_traits.
   *  This is used for numeric argument promotion of complex and cmath.
   */
  template<typename _Tp, bool = std::is_integral<_Tp>::value>
    struct fp_promote_help
    { using type = double; };

  // No nested type member for non-integer non-floating point types,
  // allows this type to be used for SFINAE to constrain overloads in
  // <cmath> and <complex> to only the intended types.
  template<typename _Tp>
    struct fp_promote_help<_Tp, false>
    { };

  template<>
    struct fp_promote_help<float>
    { using type = float; };

  template<>
    struct fp_promote_help<double>
    { using type = double; };

  template<>
    struct fp_promote_help<long double>
    { using type = long double; };

#if !defined(STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    struct fp_promote_help<__float128>
    { using type = __float128; };
#endif

  template<typename... _Tps>
    using fp_promote_help_t = typename fp_promote_help<_Tps...>::type;

  // Decay refs and cv...
  // Alternatively we could decay refs and propagate cv to promoted type.
  template<typename _Tp, typename... _Tps>
    struct fp_promote
    { using type = decltype(fp_promote_help_t<std::decay_t<_Tp>>{}
		   + typename fp_promote<_Tps...>::type{}); };

  template<>
    template<typename _Tp>
      struct fp_promote<_Tp>
      { using type = decltype(fp_promote_help_t<std::decay_t<_Tp>>{}); };

  template<typename... _Tps>
    using fp_promote_t = typename fp_promote<_Tps...>::type;

  // Assume complex value_type is floating point.
  template<>
    template<typename _Tp>
      struct fp_promote_help<std::complex<_Tp>, false>
      {
      private:
	using vtype = typename std::complex<_Tp>::value_type;
      public:
	using type = decltype(std::complex<fp_promote_help_t<vtype>>{});
      };

  // primary template handles types that have no nested ::value_type member:
  template<typename, typename = std::void_t<>>
    struct has_value_type
    : std::false_type
    { };

  // Specialization recognizes types that do have a nested ::value_type member:
  template<typename _Tp>
    struct has_value_type<_Tp, std::void_t<typename _Tp::value_type>>
    : std::true_type
    { };

  // Try to get something like array<Tp, Num> or vector<Tp, Alloc>
  // This should be able to do complex too.
  template<>
    template<template<typename _Arg, typename... _Args> typename _Tp>
      struct fp_promote_help<typename _Tp<_Arg, _Args...>,
			       typename = std::void_t<typename _Tp<_Arg, _Args...>::value_type>>
      {
      private:
	using vtype = typename _Tp<_Arg, _Args...>::value_type;
	using ptype = fp_promote_t<vtype>;
      public:
	using type = decltype(_Tp<ptype, _Args...>>{});
      };
