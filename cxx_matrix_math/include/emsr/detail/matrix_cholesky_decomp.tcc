#ifndef MATRIX_CHOLESKY_DECOMP_TCC
#define MATRIX_CHOLESKY_DECOMP_TCC 1

#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace emsr
{

/**
 * 
 */
template<typename HermitianMatrix, typename Vector>
  void
  cholesky_decomp(std::size_t n, HermitianMatrix& a, Vector& d)
  {
    for (std::size_t i = 0; i < n; ++i)
      {
	for (std::size_t j = i; j < n; ++j)
	  {
	    auto sum = a[i][j];
	    for (std::ptrdiff_t k = i - 1; k >= 0; --k)
	      sum -= a[i][k] * a[j][k];
	    if (i == j)
	      {
		if (sum <= 0)
		  throw std::logic_error("cholesky_decomp: Matrix must be positive definite");
		d[i] = std::sqrt(sum);
	      }
	    else
	      a[j][i] = sum / d[i];
	  }
      }
    //for (std::size_t j = 0; j < n; ++j)
    //  for (std::size_t i = 0; i < j; ++i)
//	a[i][j] = 0;
  }

/**
 * Solve the system @f$ Ax = b @f$ with a Cholesky decomposition.
 */
template<typename HermitianMatrix, typename Vector>
  void
  cholesky_backsub(std::size_t n, const HermitianMatrix& a,
		   const Vector& d, const Vector& b, Vector& x)
  {
    for (std::size_t i = 0; i < n; ++i)
      {
	auto sum = b[i];
	for (std::ptrdiff_t k = i - 1; k >= 0; --k)
	  sum -= a[i][k] * x[k];
	x[i] = sum / d[i];
      }
    for (std::ptrdiff_t i = n - 1; i >= 0; --i)
      {
	auto sum = x[i];
	for (std::size_t k = i + 1; k < n; ++k)
	  sum -= a[k][i] * x[k];
	x[i] = sum / d[i];
      }
  }

/**
 * This inverts the lower triangular matrix (not the original matrix).
 */
template<typename HermitianMatrix, typename Vector>
  void
  cholesky_invert(std::size_t n, const HermitianMatrix& a, const Vector& d,
		  HermitianMatrix& a_inv)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    for (std::size_t i = 0; i < n; ++i)
      {
	a_inv[i][i] = NumTp{1} / d[i];
	for (std::size_t j = i + 1; j < n; ++j)
	  {
	    auto sum = NumTp{0};
	    for (std::size_t k = i; k < j; ++k)
	      sum -= a[j][k] * a_inv[k][i];
	    a_inv[j][i] = sum / d[j];
	  }
        // Zero upper triangle.
        for (std::size_t j = 0; j < i; ++j)
          a_inv[j][i] = NumTp{0};
      }
  }

} // namespace emsr

#endif // MATRIX_CHOLESKY_DECOMP_TCC
