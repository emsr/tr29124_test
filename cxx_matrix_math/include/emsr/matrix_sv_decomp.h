#ifndef MATRIX_SV_DECOMP_H
#define MATRIX_SV_DECOMP_H

#include <vector>
#include <utility>

namespace emsr
{

/**
 *  This class represents an singular value decomposition of a matrix.
 */
template<typename Matrix>
  class sv_decomposition
  {

  public:

    using NumTp = std::decay_t<decltype(Matrix{}[0][0])>;

    sv_decomposition(std::size_t m_n_rows, std::size_t n_cols, const Matrix& a);

    template<typename Matrix2>
      sv_decomposition(std::size_t m_n_rows, std::size_t n_cols,
		       const Matrix2& a);

    template<typename Vector2, typename VectorOut>
      void
      backsubstitution(const Vector2& b, VectorOut& x) const;

    template<typename Matrix2, typename Vector2, typename VectorOut>
      void
      improve(const Matrix2 a_orig,
	      const Vector2& b, VectorOut& x) const;

  private:

    std::size_t m_n_rows;

    std::size_t m_n_cols;

    Matrix m_a;

    std::vector<NumTp> m_w;

    std::vector<std::vector<NumTp>> m_v;
  };

/**
 *  
 */
template<typename Matrix, typename Vector, typename MatrixV>
  void
  sv_decomp(const std::size_t n_rows, const std::size_t n_cols,
	    Matrix& a, Vector& w, MatrixV& v);

/**
 *  
 */
template<typename Matrix, typename Vector, typename MatrixV>
  void
  sv_backsub(std::size_t n_rows, std::size_t n_cols,
	     const Matrix& u,
	     const Vector& w, const MatrixV& v,
	     const Vector& b, Vector& x);

/**
 *  Improves a solution vector x of the linear set A.x = b.
 *  The Matrix a and the SV decomposition of a -- u, w, v and the
 *  right-hand side Vector are input along with the solution vector x.
 *  The solution vector x is improved and modified on output.
 */
template<typename Matrix, typename Vector, typename MatrixV>
  void
  sv_improve(std::size_t n_rows, std::size_t n_cols,
	     const Matrix& a, const Matrix& u,
	     const Vector& w, const MatrixV& v,
	     const Vector& b, Vector& x);

} // namespace emsr

#include <emsr/detail/matrix_sv_decomp.tcc>

#endif // MATRIX_SV_DECOMP_H
