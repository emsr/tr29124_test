#ifndef FLOAT128_TCC
#define FLOAT128_TCC 1

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)

#include <iostream>
#include <iomanip> // For setw().
#include <sstream>
#include <quadmath.h>

namespace std
{

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	       const __float128& __x)
    {
      std::ostringstream __fmt;
      __fmt << "%." << __os.precision() << "Qg";
std::cerr << __fmt.str() << std::endl;
      constexpr int __strlen = 1000;
      char __str[__strlen];
      quadmath_snprintf(__str, __strlen, __fmt.str().c_str(), __x) ;
      __os.precision(30);
      __os << __str;
      return __os;
    }

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_istream<_CharT, _Traits>&
    operator>>(std::basic_istream<_CharT, _Traits>& __is, __float128& __x)
    {
      constexpr int __strlen = 1000;
      char __str[__strlen];
      __is >> std::setw(__strlen) >> __str;
      __x = strtoflt128(__str, 0);
      return __is;
    }

} // namespace std

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // FLOAT128_TCC
