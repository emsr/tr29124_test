
#ifndef STL_ENDIAN_H
#define STL_ENDIAN_H 1

namespace std
{

  enum class endian
  {
#ifdef _WIN32
      little = 0,
      big    = 1,
      native = little
#else
      little = __ORDER_LITTLE_ENDIAN__,
      big    = __ORDER_BIG_ENDIAN__,
      native = __BYTE_ORDER__
#endif
  };

} // namespace std

#endif // STL_ENDIAN_H
