#ifndef MATRIX_LU_DECOMP_TCC
#define MATRIX_LU_DECOMP_TCC 1

#include <cstdlib>
#include <vector>
#include <cmath>
#include <limits>

#include <emsr/matrix_util.h>

namespace emsr
{

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
	    Vector& index, NumTp& parity)
  {
    const NumTp TINY = NumTp(1.0e-20L);

    std::vector<NumTp> scale(n);
    parity = NumTp{1};

    // Loop over rows to get the implicit scaling information.    
    for (std::size_t i = 0; i < n; ++i)
      {
	NumTp big{0};
	for (std::size_t j = 0; j < n; ++j)
	  if (const auto temp = std::abs(a[i][j]); temp > big)
	    big = temp;
	if (big == NumTp{0})
	  throw std::logic_error("lu_decomp: singular matrix");

	// Save the scaling for the row.
	scale[i] = NumTp{1} / big;
      }

    // This is the loop over columns of Crout's method.
    for (std::size_t j = 0; j < n; ++j)
      {
	// Lower triangle.
	for (std::size_t i = 0; i < j; ++i)
	  {
	    auto sum = a[i][j];
	    for (std::size_t k = 0; k < i; ++k)
	      sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	  }

	// Initialize for the search for the largest pivot point.
	auto imax = std::numeric_limits<std::size_t>::max();
	auto big = NumTp{0};
	for (std::size_t i = j; i < n; ++i)
	  {
	    auto sum = a[i][j];
	    for (std::size_t k = 0; k < j; ++k)
	      sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	    if (const auto dum = scale[i] * std::abs(sum); dum >= big)
	      {
		big = dum;
		imax = i;
	      }
	  }
	// Interchange rows if required.
	if (j != imax)
	  {
	    for (std::size_t k = 0; k < n; ++k)
	      std::swap(a[imax][k], a[j][k]);

	    // Change parity.
	    parity = -parity;

	    // Interchange the scale factor.
	    std::swap(scale[imax], scale[j]);
	  }
	index[j] = imax;
	if (a[j][j] == NumTp{0})
	  a[j][j] = TINY;

	// Now finally divide by the pivot element
	if (j != n - 1)
	  {
	    const auto scale = NumTp{1} / a[j][j];
	    for (std::size_t i = j + 1; i < n; ++i)
	      a[i][j] *= scale;
	  }
      } // Go back for the next column in the reduction.
  }

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
	     Vector& b)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    //  When i_start is set to a non-negative value, it will become the index
    //  of the first nonvanishing element of b[0..n-1].
    //  Do the forward substitution unsrambling the permutation as we go.
    int i_start = -1;
    for (std::size_t i = 0; i < n; ++i)
      {
	auto i_perm = index[i];
	auto sum = b[i_perm];
	b[i_perm] = b[i];
	if (i_start > -1)
	  for (std::size_t j = i_start; j <= i - 1; ++j)
	    sum -= a[i][j] * b[j];
	else if (sum != NumTp{0})
	  i_start = i;
	b[i] = sum;
      }

    //  Now do the backsubstitution.
    for (std::ptrdiff_t i = n - 1; i >= 0; --i)
      {
	auto sum = b[i];
	for (std::size_t j = i + 1; j < n; ++j)
	  sum -= a[i][j] * b[j];
	b[i] = sum / a[i][i];
      }

    return;
  }

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
	     const VectorInt& index, const Vector& b, Vector& x)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    std::vector<NumTp> r(n);

    for (std::size_t i = 0; i < n; ++i)
      {
	r[i] = -b[i];
	for (std::size_t j = 0; j < n; ++j)
	  r[i] += a[i][j] * x[j];
      }

    lu_backsub(a_lu, n, index, r);

    for (std::size_t i = 0; i < n; ++i)
      x[i] -= r[i];

    return;
  }

/**
 * Inverts a matrix given the LU decomposed matrix.
 *
 * The inverse matrix is NOT in LU form.
 */
template<typename SquareMatrix, typename VectorInt>
  void
  lu_invert(const std::size_t n, const SquareMatrix& a_lu,
	    const VectorInt& index, SquareMatrix& a_inv)
  {
    using NumTp = std::decay_t<decltype(a_inv[0][0])>;

    for (std::size_t j = 0; j < n; ++j)
      {
	std::vector<NumTp> col(n);
	col[j] = NumTp{1};

	lu_backsub(n, a_lu, index, col);

	for (std::size_t i = 0; i < n; ++i)
	  a_inv[i][j] = col[i];
      }

    return;
  }

/**
 * Compute determinant of LU decomposed matrix.
 */
template<typename NumTp, typename SquareMatrix>
  NumTp
  lu_determinant(const std::size_t n, const SquareMatrix& a_lu, const NumTp parity)
  {
    auto det = parity;

    for (std::size_t i = 0; i < n; ++i)
      det *= a_lu[i][i];

    return det;
  }

/**
 * Compute trace of LU decomposed matrix.
 */
template<typename SquareMatrix>
  auto
  lu_trace(const std::size_t n, const SquareMatrix& a_lu)
  {
    using NumTp = std::decay_t<decltype(a_lu[0][0])>;

    auto trace = NumTp{0};

    for (std::size_t i = 0; i < n; ++i)
      {
	trace += a_lu[i][i];
	for (std::ptrdiff_t j = i - 1; j >= 0; --j)
	  trace += a_lu[i][j] * a_lu[j][i];
      }

    return trace;
  }

// Implement class methods.

template<typename SquareMatrix>
  lu_decomposition<SquareMatrix>::
  lu_decomposition(std::size_t n, const SquareMatrix& a)
  : m_n(n), m_a(a), m_index(std::vector<std::size_t>(n)), m_parity(+1)
  {
    lu_decomp(this->m_n, this->m_a, this->m_index, this->m_parity);
  }

template<typename SquareMatrix>
  template<typename SquareMatrix2>
  lu_decomposition<SquareMatrix>::
  lu_decomposition(std::size_t n, const SquareMatrix2& a)
  : m_n(n), m_a(n, std::vector<NumTp>(n)),
    m_index(std::vector<std::size_t>(n)), m_parity(+1)
  {
    // Copy a.
    for (std::size_t i_row = 0; i_row < this->m_n; ++i_row)
      for (std::size_t i_col = 0; i_col < this->m_n; ++i_col)
        this->m_a[i_row][i_col] = a[i_row][i_col];
    lu_decomp(this->m_n, this->m_a, this->m_index, this->m_parity);
  }

template<typename SquareMatrix>
  template<typename Vector>
  void
  lu_decomposition<SquareMatrix>::
  backsubstitute(Vector& b) const
  {
    lu_backsub(this->m_n, this->m_a, this->m_index, b);
  }

template<typename SquareMatrix>
  template<typename SquareMatrix2, typename Vector, typename VectorOut>
  void
  lu_decomposition<SquareMatrix>::
  improve(const SquareMatrix2& a_orig,
	  const Vector& b, VectorOut& x) const
  {
    lu_improve(this->m_n, a_orig,
	       this->m_a, this->m_index,
	       b, x);
  }

template<typename SquareMatrix>
  typename lu_decomposition<SquareMatrix>::NumTp
  lu_decomposition<SquareMatrix>::
  determinant() const
  {
    return lu_determinant(this->m_n, this->m_a, this->m_parity);
  }

template<typename SquareMatrix>
  typename lu_decomposition<SquareMatrix>::NumTp
  lu_decomposition<SquareMatrix>::
  trace() const
  {
    return lu_trace(this->m_n, this->m_a);
  }

} // namespace emsr

#endif // MATRIX_LU_DECOMP_TCC
