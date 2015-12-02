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

} // namespace std

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // COMPLEX128_H
