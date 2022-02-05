
#ifdef EMSR_HAVE_FLOAT128

  __float128
  operator""_maxflt(const char* str)
  { return strtoflt128(str, 0); }

#else

  long double
  operator""_maxflt(const char* str)
  { return strtold(str, 0); }

#endif
