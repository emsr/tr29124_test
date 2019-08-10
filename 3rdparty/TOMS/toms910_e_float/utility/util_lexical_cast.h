
#ifndef _UTIL_LEXICAL_CAST_2009_11_24_H_
  #define _UTIL_LEXICAL_CAST_2009_11_24_H_

  #include <string>
  #include <sstream>

  namespace Util
  {
    template<typename T> inline std::string lexical_cast(const T& t)
    {
      std::stringstream ss;
      ss << t;
      return ss.str();
    }
  }

#endif // _UTIL_LEXICAL_CAST_2009_11_24_H_
