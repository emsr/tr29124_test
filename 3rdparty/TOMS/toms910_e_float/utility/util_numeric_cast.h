
#ifndef _UTIL_NUMERIC_CAST_2009_11_24_H_
  #define _UTIL_NUMERIC_CAST_2009_11_24_H_

  #include <string>
  #include <sstream>

  namespace Util
  {
    template<typename T> inline T numeric_cast(const std::string& str)
    {
      std::stringstream ss;
      ss << str;
      T t;
      ss >> t;
      return t;
    }
  }

#endif // _UTIL_NUMERIC_CAST_2009_11_24_H_
