#ifndef MATRIX_LU_DECOMP_H
#define MATRIX_LU_DECOMP_H 1

#include <vector>

namespace emsr
{

/**
 * This class represents an lower-upper decomposition of a square matrix.
 */
template<typename NumTp, typename SquareMatrix>
  class lu_decomposition
  {

  public:

    using value_type = std::decay_t<decltype(SquareMatrix{}[0][0])>;

    template<typename SquareMatrix2>
      lu_decomposition(std::size_t n, const SquareMatrix2& a);

    template<typename Vector, typename VectorOut>
      void backsubstitute(const Vector& b, VectorOut& x) const;

    template<typename SquareMatrix2, typename Vector, typename VectorOut>
      void
      improve(const SquareMatrix2& a_orig,
	      const Vector& b, VectorOut& x) const;

    template<typename SquareMatrix2>
      void
      inverse(SquareMatrix2& a_inv) const;

    NumTp determinant() const;

    NumTp trace() const;

  private:

    std::size_t m_n;

    SquareMatrix m_a;

    std::vector<std::size_t> m_index;

    int m_parity;
  };

/**
 * Given an n*n matrix a[0..n-1][0..n-1], this routine replaces it
 * by the LU (Lower-triangular Upper-triangular)  decomposition
 * of a rowwise permutation of itself.
 * The matrix size n and a[][] are input.  a[][] is output, index[]
 * is an output vector which row permutation effected by the partial pivoting;
 * d is output as the parity of the row permutation
 */
template<typename NumTp, typename SquareMatrix, typename Vector>
  void
  lu_decomp(std::size_t n, SquareMatrix& a,
	    Vector& index, NumTp& parity);

/**
 * Solve the set of n linear equations a.x = b.  Here a[0..n-1][0..n-1] is input, not as the original matrix a but as 
 * its LU decomposition, determined by the routine lu_decomp().  b[0..n-1] is input as the right hand side vector b 
 * and returns with the left-hand solution vector x.  a, n, and index are not modified by this routine and can be left 
 * in place for successive calls with different right hand sides b[0..n-1].  This routine takes into account the 
 * possibility that b will begin with a lot of zeros so that it is efficient for use in matrix inversion.
 */
template<typename SquareMatrix, typename VectorInt, typename Vector>
  void
  lu_backsub(const std::size_t n,
	     const SquareMatrix& a,
	     const VectorInt& index,
	     Vector& b);

/**
 * Improves a solution vector x of the linear set A.x = b.  The matrix a and the
 * LU decomposition of a a_lu (with its row permutation vector index) and the
 * right-hand side vector are input along with the solution vector x.  
 * The solution vector x is improved and modified on output.
 */
template<typename SquareMatrix, typename VectorInt, typename Vector>
  void
  lu_improve(const std::size_t n, const SquareMatrix& a,
	     const SquareMatrix& a_lu,
	     const VectorInt& index, const Vector& b, Vector& x);

/**
 * Inverts a matrix given the LU decomposed matrix.
 *
 * The inverse matrix is NOT in LU form.
 */
template<typename SquareMatrix, typename VectorInt>
  void
  lu_invert(const std::size_t n, const SquareMatrix& a_lu,
	    const VectorInt& index, SquareMatrix& a_inv);

/**
 * Compute determinant of LU decomposed matrix.
 */
template<typename NumTp, typename SquareMatrix>
  NumTp
  lu_determinant(const std::size_t n, const SquareMatrix& a_lu, const NumTp parity);

/**
 * Compute trace of LU decomposed matrix.
 */
template<typename SquareMatrix>
  auto
  lu_trace(const std::size_t n, const SquareMatrix& a_lu);

} // namespace emsr

#include <emsr/detail/matrix_lu_decomp.tcc>

#endif // MATRIX_LU_DECOMP_H
