
#include <type_traits>


template<typename _ArgT, typename _Coef0>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner(_ArgT __x, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return __arg_t{__c0};
  }

template<typename _ArgT, typename _Coef0, typename... _Coef>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner(_ArgT __x, _Coef0 __c0, _Coef... __c)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return __arg_t{__c0} + __x * horner(__x, __c...);
  }


template<typename _ArgT, typename _Coef0>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner_big_end(_ArgT, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return __arg_t{__c0};
  }

template<typename _ArgT, typename _Coef1, typename _Coef0>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner_big_end(_ArgT __x, _Coef1 __c1, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return horner_big_end(__x, __x * __arg_t{__c1} + __arg_t{__c0});
  }

template<typename _ArgT, typename _CoefN, typename _CoefNm1, typename... _Coef>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner_big_end(_ArgT __x, _CoefN __cn, _CoefNm1 __cnm1, _Coef... __c)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return horner_big_end(__x, __x * __arg_t{__cn} + __arg_t{__cnm1}, __c...);
  }

/**
 * As I've found out, you can't have a constexpr initializer_list.
 * What Darko did was have a static polynomial templated on an (empty) tag type,
 * the order, and the float type.  The tag type groups all the little monomials into a single polynomial!
 * It would be nice if I could find a way to have a static constexpr Float member
 * and lose the macros.
 */

// Let's try a variadic template polynomial ...
/* ... 'double' is not a valid type for a template non-type parameter
template<typename Tp, Tp... Coef>
  struct Poly
  {
    template<typename Up>
      Up
      eval(Up x)
      { return horner(x, Coef...); }
  };
*/
/*
template<typename Tag, typename Tp, std::ptrdiff_t Order>
  struct Poly
  {
    inline static constexpr Tp Coef;
  };
*/
