#ifndef MATRIX_CHOLESKY_DECOMP_H
#define MATRIX_CHOLESKY_DECOMP_H 1

namespace emsr
{

/**
 * This class represents a Cholesky decomposition of a square matrix.
 */
template<typename HermitianMatrix, typename _Vector>
  class cholesky_decomposition
  {

  public:

    using value_type = decltype(HermitianMatrix{}[0][0]);

    template<typename HermitianMatrix2, typename _Vector2>
      cholesky_decomposition(std::size_t n, const HermitianMatrix2& a, _Vector2& d);

    template<typename Vector2, typename VectorOut>
      void backsubstitute(const Vector2  b, VectorOut& x) const;

    template<typename InVecIter, typename OutVecIter>
      void
      backsubstitution(InVecIter b_begin, InVecIter b_end,
		       OutVecIter x_begin) const;

    template<typename HermitianMatrix2>
      void inverse(HermitianMatrix2& a_inv) const;

  private:

    std::size_t m_n;

    HermitianMatrix m_a;

    std::vector<value_type> m_d;
  };

/**
 * 
 */
template<typename HermitianMatrix, typename Vector>
  void
  cholesky_decomp(std::size_t n, HermitianMatrix& a, Vector& d);

/**
 * Solve the system @f$ Ax = b @f$ with a Cholesky decomposition.
 */
template<typename HermitianMatrix, typename Vector>
  void
  cholesky_backsub(std::size_t n, const HermitianMatrix& a,
		   const Vector& d, const Vector& b, Vector& x);

/**
 * 
 */
template<typename HermitianMatrix, typename Vector>
  void
  cholesky_invert(std::size_t n, const HermitianMatrix& a, const Vector& d,
		  HermitianMatrix& a_inv);

} // namespace emsr

#include <emsr/detail/matrix_cholesky_decomp.tcc>

#endif // MATRIX_CHOLESKY_DECOMP_H

