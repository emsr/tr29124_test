namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // Tuple interface to class template complex.
#define __cxx_lib_tuple_complex 201705

  /// tuple_size
  template<typename Tp>
    class tuple_size;

  /// Partial specialization for std::complex
  template<typename Tp>
    struct tuple_size<std::complex<Tp>>
    : public integral_constant<std::size_t, 2> { };

  /// tuple_element
  template<std::size_t _Int, typename Tp>
    class tuple_element;

  /// Partial specialization for std::complex
  template<std::size_t _Int, typename Tp>
    struct tuple_element<_Int, std::complex<Tp>>
    {
      static_assert(_Int < 2, "index is out of bounds");
      typedef Tp type;
    };

  template<typename Tp>
    struct is_tuple_like_impl<std::complex<Tp>> : true_type
    { };

  /// Decompose complex as lvalues.
  template<std::size_t _Int, typename Tp>
    constexpr Tp&&
    get(std::complex<Tp>&& z)
    {
      static_assert(_Int < 2, "index is out of bounds");
      return _Int == 0 ? z.real() : z.imag();
    }

  /// Decompose complex as lvalues.
  template<std::size_t _Int, typename Tp>
    constexpr const Tp
    get(const std::complex<Tp>& z)
    {
      static_assert(_Int < 2, "index is out of bounds");
      return _Int == 0 ? z.real() : z.imag();
    }

  /// Decompose complex as rvalues.
  template<std::size_t _Int, typename Tp>
    constexpr Tp&
    get(std::complex<Tp>& z)
    {
      static_assert(_Int < 2, "index is out of bounds");
      auto w = reinterpret_cast<Tp(&)[2]>(z);
      return _Int == 0 ? w[0] : w[1];
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std
