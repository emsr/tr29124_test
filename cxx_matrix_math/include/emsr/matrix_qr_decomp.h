#ifndef MATRIX_QR_DECOMP_H
#define MATRIX_QR_DECOMP_H 1

#include <vector>

namespace emsr
{

/**
 *  This class represents an QR decomposition.
 */
template<typename NumTp, typename Matrix>
  class qr_decomposition
  {

  public:

    using value_type = std::decay_t<decltype(Matrix{}[0][0])>;

    template<typename Matrix2>
      qr_decomposition(std::size_t m_n_rows, std::size_t n_cols,
		       Matrix2& a);

    template<typename Vector2, typename VectorOut>
      void backsubstitution(Vector2& b, VectorOut& x);

    template<typename InVecIter, typename OutVecIter>
      void
      backsubstitution(InVecIter b_begin, InVecIter b_end,
		       OutVecIter x_begin) const;

    template<typename Matrix2>
      void inverse(Matrix2& a_inv);

    template<typename Matrix2, typename Vector2>
      void update(Matrix2& r, Matrix2& qt,
		  Vector2& u, Vector2& v);

  private:

    std::size_t m_n_rows;

    std::size_t m_n_cols;

    Matrix m_a;

    std::vector<std::size_t> m_c;

    std::vector<std::size_t> m_d;

    bool m_singular;
  };

/**
 * Constructs the QR decomposition of a[0..n_rows - 1][0..n_cols - 1].  The upper triangular matrix R
 * is returned in the upper triangle of a except the diagonal elements of R which are returned
 * in d[0..n_cols - 1].  The orthogonal matrix Q is represented as a product of n - 1 Householder
 * matrices Q_0...Q_n-2 where Q_j = 1 - u_j x u_j/c_j.  The ith component of u_j is zero for
 * i = 
 */
template<typename Matrix, typename Vector>
  void
  qr_decomp(std::size_t n_rows, std::size_t n_cols,
	    Matrix& a,
	    Vector& c, Vector& d, bool& singular);

/**
 * This routine solves the set of equations Rx = b where R is the upper triangular
 * matrix stored in a[0..n_rows - 1][0..n_cols - 1] and d[0..n_cols - 1].
 * Here n_rows >= n_cols.
 */
template<typename Matrix, typename Vector, typename VectorB>
  void
  r_backsub(std::size_t n_rows, std::size_t n_cols,
	    const Matrix& a, const Vector& d,
	    VectorB& b);

/**
 * This routine solves the set of equations Ax = b.
 * The inputs are the QR decomposition of the matrix in a[0..n_rows - 1][0..n_cols - 1],
 * c[0..n_cols - 1], and d[0..n_cols - 1].
 * The vector b[0..n_rows - 1] is input as the "RHS" and output and the solution.
 * Here n_rows >= n_cols.
 */
template<typename Matrix, typename Vector, typename VectorB>
  void
  qr_backsub(const std::size_t n_rows, const std::size_t n_cols,
	     const Matrix& a, const Vector& c, const Vector& d,
	     VectorB& b);

/**
 * Inverts a matrix given the QR decomposed matrix.
 * The inverse matrix is allocated in this routine so make sure the pointer is freed first.
 *
 * The inverse matrix is NOT in QR form.
 */
template<typename Matrix, typename Vector>
  void
  qr_invert(std::size_t n_rows, std::size_t n_cols,
	    const Matrix& a_qr, const Vector& c, const Vector& d,
	    Matrix& a_inv);

/**
 *  Update the QR decomposition.
 */
template<typename Matrix, typename Vector>
  void
  qr_update(std::size_t n_rows, std::size_t n_cols,
	    Matrix& r, Matrix& qt,
	    Vector& u, Vector& v);

/**
 *  Do a Jacobi rotation on rows i and i+1 of the matrices r[0..n_cols - 1][0..n_cols - 1]
 *  and qt[0..n_cols - 1][0..n_rows - 1].
 *  The angle is specified with a and b such that
 *    cos(theta) = a/sqrt(a^2 + b^2)
 *    sin(theta) = b/sqrt(a^2 + b^2).
 */
template<typename NumTp, typename Matrix, typename Vector>
  void
  jacobi_rotate(const int i, const int n_rows, const int n_cols,
		Matrix& r, Matrix& qt,
		NumTp a, NumTp b);


} // namespace emsr

#include <emsr/detail/matrix_qr_decomp.tcc>

#endif // MATRIX_QR_DECOMP_H

