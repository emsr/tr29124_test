#ifndef MATRIX_VANDERMONDE_H
#define MATRIX_VANDERMONDE_H 1

namespace emsr
{

template<typename Tp>
  void
  vandermonde_moment(std::size_t n, const Tp* x, const Tp* q, Tp* w);

template<typename Tp>
  void
  vandermonde(std::size_t n, const Tp* x, const Tp* y, Tp* c);

} // namespace emsr

#include <emsr/detail/matrix_vandermonde.tcc>

#endif // MATRIX_VANDERMONDE_H
