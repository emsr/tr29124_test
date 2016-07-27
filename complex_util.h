#ifndef COMPLEX_UTIL_H
#define COMPLEX_UTIL_H 1

/**
 * Introspection class to detect if a type is std::complex.
 */
template<typename _Tp>
  struct is_complex : public std::false_type
  { };

/**
 * Introspection class to detect if a type is std::complex.
 */
template<>
  template<typename _Tp>
    struct is_complex<std::complex<_Tp>> : public std::true_type
    { };

/**
 * Introspection type to detect if a type is std::complex.
 */
template<typename _Tp>
  using is_complex_t = typename is_complex<_Tp>::type;

/**
 * Introspection variable template to detect if a type is std::complex.
 */
template<typename _Tp>
  constexpr bool is_complex_v = is_complex<_Tp>::value;

#endif // STATISTICS_H
