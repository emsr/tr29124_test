
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
__float128
operator""_maxflt(const char* __str)
{
  return strtoflt128(__str, 0);
}
#else
long double
operator""_maxflt(const char* __str)
{
  return strtold(__str, 0);
}
#endif
