#ifndef FLOAT128_H
#define FLOAT128_H 1

#if defined(__GNUC__) && !defined(__clang__)
#  if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
#    if __has_include(<quadmath.h>)
#      define EMSR_HAVE_FLOAT128
#    endif
#  endif
#endif

#if defined(__clang__)
// FIXME: I think clang has this
#endif

#if defined (_MSC_VER)
#endif

#endif // FLOAT128_H
