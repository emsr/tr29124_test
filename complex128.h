#ifndef COMPLEX128_H
#define COMPLEX128_H 1

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)

#include "float128.h"

namespace std
{

  inline __float128
  abs(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cabsq(__z); }

  inline __float128
  arg(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cargq(__z); }

  inline __float128
  imag(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cimagq(__z); }

  inline __float128
  real(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return crealq(__z); }

  inline __complex128 __z
  acos(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cacosq(__z); }

  inline __complex128 __z
  acosh(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cacoshq(__z); }

  inline __complex128 __z
  asin(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return casinq(__z); }

  inline __complex128 __z
  asinh(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return casinhq(__z); }

  inline __complex128 __z
  atan(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return catanq(__z); }

  inline __complex128 __z
  atanh(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return catanhq(__z); }

  inline __complex128 __z
  cos(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return ccosq(__z); }

  inline __complex128 __z
  cosh(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return ccoshq(__z); }

  inline __complex128 __z
  exp(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cexpq(__z); }

  inline __complex128 __z
  expi(__float128) _GLIBCXX_USE_NOEXCEPT
  { return cexpiq(__z); }

  inline __complex128 __z
  log(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return clogq(__z); }

  inline __complex128 __z
  log10(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return clog10q(__z); }

  inline __complex128 __z
  conj(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return conjq(__z); }

  inline __complex128 __z
  pow(__complex128 __z, __complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cpowq(__z); }

  inline __complex128 __z
  proj(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return cprojq(__z); }

  inline __complex128 __z
  sin(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return csinq(__z); }

  inline __complex128 __z
  sinh(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return csinhq(__z); }

  inline __complex128 __z
  sqrt(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return csqrtq(__z); }

  inline __complex128 __z
  tan(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return ctanq(__z); }

  inline __complex128 __z
  tanh(__complex128 __z) _GLIBCXX_USE_NOEXCEPT
  { return ctanhq(__z); }

  /// complex<__float128> specialization
  template<>
    struct complex<__float128>
    {
      typedef __float128 value_type;
      typedef __complex__ __float128 _ComplexT;

      _GLIBCXX_CONSTEXPR complex(_ComplexT __z) : _M_value(__z) { }

      _GLIBCXX_CONSTEXPR complex(__float128 __r = 0.0Q, 
				 __float128 __i = 0.0Q)
#if __cplusplus >= 201103L
      : _M_value{ __r, __i } { }
#else
      {
	__real__ _M_value = __r;
	__imag__ _M_value = __i;
      }
#endif

      _GLIBCXX_CONSTEXPR complex(const complex<float>& __z)
      : _M_value(__z.__rep()) { }

      _GLIBCXX_CONSTEXPR complex(const complex<double>& __z)
      : _M_value(__z.__rep()) { }

      _GLIBCXX_CONSTEXPR complex(const complex<long double>& __z)
      : _M_value(__z.__rep()) { }

#if __cplusplus >= 201103L
      // _GLIBCXX_RESOLVE_LIB_DEFECTS
      // DR 387. std::complex over-encapsulated.
      __attribute ((__abi_tag__ ("cxx11")))
      constexpr __float128 
      real() const { return __real__ _M_value; }

      __attribute ((__abi_tag__ ("cxx11")))
      constexpr __float128 
      imag() const { return __imag__ _M_value; }
#else
      __float128& 
      real() { return __real__ _M_value; }

      const __float128& 
      real() const { return __real__ _M_value; }

      __float128& 
      imag() { return __imag__ _M_value; }

      const __float128& 
      imag() const { return __imag__ _M_value; }
#endif

      // _GLIBCXX_RESOLVE_LIB_DEFECTS
      // DR 387. std::complex over-encapsulated.
      void 
      real(__float128 __val) { __real__ _M_value = __val; }

      void 
      imag(__float128 __val) { __imag__ _M_value = __val; }

      complex&
      operator=(__float128 __r)
      {
	_M_value = __r;
	return *this;
      }

      complex&
      operator+=(__float128 __r)
      {
	_M_value += __r;
	return *this;
      }

      complex&
      operator-=(__float128 __r)
      {
	_M_value -= __r;
	return *this;
      }

      complex&
      operator*=(__float128 __r)
      {
	_M_value *= __r;
	return *this;
      }

      complex&
      operator/=(__float128 __r)
      {
	_M_value /= __r;
	return *this;
      }

      // The compiler knows how to do this efficiently
      // complex& operator=(const complex&);

      template<typename _Tp>
        complex&
        operator=(const complex<_Tp>& __z)
	{
	  __real__ _M_value = __z.real();
	  __imag__ _M_value = __z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator+=(const complex<_Tp>& __z)
	{
	  __real__ _M_value += __z.real();
	  __imag__ _M_value += __z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator-=(const complex<_Tp>& __z)
	{
	  __real__ _M_value -= __z.real();
	  __imag__ _M_value -= __z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator*=(const complex<_Tp>& __z)
	{
	  _ComplexT __t;
	  __real__ __t = __z.real();
	  __imag__ __t = __z.imag();
	  _M_value *= __t;
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator/=(const complex<_Tp>& __z)
	{
	  _ComplexT __t;
	  __real__ __t = __z.real();
	  __imag__ __t = __z.imag();
	  _M_value /= __t;
	  return *this;
	}

      _GLIBCXX_CONSTEXPR _ComplexT __rep() const { return _M_value; }

    private:
      _ComplexT _M_value;
    };

  // @todo Ctors from larger types are marked explicit in the smaller classes.
  inline _GLIBCXX_CONSTEXPR
  complex<float>::complex(const complex<__float128>& __z)
  : _M_value(__z.__rep()) { }

  inline _GLIBCXX_CONSTEXPR
  complex<double>::complex(const complex<__float128>& __z)
  : _M_value(__z.__rep()) { }

  inline _GLIBCXX_CONSTEXPR
  complex<long double>::complex(const complex<__float128>& __z)
  : _M_value(__z.__rep()) { }

#if __cplusplus > 201103L

inline namespace literals {
inline namespace complex_literals {
/*
  constexpr std::complex<__float128>
  operator""iq(long double __num)
  { return std::complex<__float128>{0.0Q, __num}; }

  constexpr std::complex<__float128>
  operator""iq(unsigned long long __num)
  { return std::complex<__float128>{0.0Q, static_cast<__float128>(__num)}; }
*/
  constexpr std::complex<__float128>
  operator""iq(const char* __str)
  { return strtoflt128(__str, 0); }

} // inline namespace complex_literals
} // inline namespace literals

#endif // C++14

} // namespace std

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // COMPLEX128_H
