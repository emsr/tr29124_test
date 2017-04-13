namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // Tuple interface to class template complex.

  /// tuple_size
  template<typename _Tp>
    class tuple_size;

  /// Partial specialization for std::complex
  template<typename _Tp>
    struct tuple_size<std::complex<_Tp>>
    : public integral_constant<std::size_t, 2> { };

  /// tuple_element
  template<std::size_t _Int, typename _Tp>
    class tuple_element;

  /// Partial specialization for std::complex
  template<std::size_t _Int, typename _Tp>
    struct tuple_element<_Int, std::complex<_Tp>>
    {
      static_assert(_Int < 2, "index is out of bounds");
      typedef _Tp type;
    };

  template<typename _Tp>
    struct __is_tuple_like_impl<std::complex<_Tp>> : true_type
    { };

  template<std::size_t _Int, typename _Tp>
    constexpr const _Tp
    get(const std::complex<_Tp>& __z)
    {
      static_assert(_Int < 2, "index is within bounds");
      return _Int == 0 ? __z.real() : __z.imag();
    }

  template<std::size_t _Int, typename _Tp>
    constexpr _Tp&
    get(std::complex<_Tp>& __z)
    {
      static_assert(_Int < 2, "index is within bounds");
      auto __w = reinterpret_cast<_Tp(&)[2]>(__z);
      return _Int == 0 ? __w[0] : __w[1];
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std
