#ifndef MATRIX_GAUSS_JORDAN_TCC
#define MATRIX_GAUSS_JORDAN_TCC 1

#include <cstdlib>
#include <vector>
#include <cmath>
#include <limits>

namespace emsr
{

/**
 * Linear equation solution by Gauss-Jordan elimination. a[1..n][1..n] is an input matrix.
 * b[1..n][1..m] is the input matrix of m right-hand side vectors.
 * On output, a is replaced by its inverse, and b is replaced by the corresponding set of
 *solution vectors.
 */
template<typename NumTp, typename SquareMatrix, typename Matrix>
  void
  gauss_jordan(SquareMatrix& a, std::size_t n, Matrix& b, std::size_t m)
  {
    const auto imax = std::numeric_limits<std::size_t>::max();

    std::vector<std::size_t> index_col(n, imax);
    std::vector<std::size_t> index_row(n, imax);
    std::vector<std::size_t> index_pivot(n, imax);

    for (std::size_t i = 0; i < n; ++i)
      {
	// Loop over columns to be reduced.
	std::size_t irow = imax;
	std::size_t icol = imax;
	auto big = NumTp{};
	for (std::size_t j = 0; j < n; ++j)
	  {
	    // Loop over rows looking for the pivot elements.
	    if (index_pivot[j] != 0)
	      for (std::size_t k = 0; k < n; ++k)
		{
		  if (index_pivot[k] == imax)
		    {
		      if (const auto tmp = std::abs(a[j][k]); tmp >= big)
			{
			  big = tmp;
			  irow = j;
			  icol = k;
			}
		    }
		  else if (index_pivot[k] > 0)
		    std::__throw_runtime_error("gauss_jordan: Singular matrix");
		}
	  }
	++index_pivot[icol];

	// With the pivot elements in hand, we swap rows to put
	// the pivot elements on the diagonal.
	// The columns are not physically moved, only relabeled:
	//   index_pivot[i] is the column of the original ith pivot element,
	//                  it is the ith column that is reduced,
	//   index_row[i] is the row in which that pivot element
	//                was originally located.
	// If index_row[i] != index_col[i] a column interchange is implied.
        if (irow != icol)
	  {
	    for (auto l = 0u; l < n; ++l)
	      std::swap(a[irow][l], a[icol][l]);
	    for (auto l = 0u; l < m; ++l)
	      std::swap(b[irow][l], b[icol][l]);
	  }
	index_row[i] = irow;
	index_col[i] = icol;
	if (a[icol][icol] == NumTp{})
	  std::__throw_runtime_error("gauss_jordan: Singular matrix");
	const auto pivinv = NumTp{1} / a[icol][icol];
	a[icol][icol] = NumTp{1};
	for (auto l = 0u; l < n; ++l)
	  a[icol][l] *= pivinv;
	for (auto l = 0u; l < m; ++l)
	  b[icol][l] *= pivinv;
	for (auto ll = 0u; ll < n; ++ll)
	  {
	    if (ll != icol)
	      {
		const auto dum = a[ll][icol];
		a[ll][icol] = NumTp{};
		for (auto l = 0u; l < n; ++l)
		  a[ll][l] -= dum * a[icol][l];
		for (auto l = 0u; l < m; ++l)
		  b[ll][l] -= dum * b[icol][l];
	      }
	  }
      }
    for (int l = n - 1; l >= 0; --l)
      if (index_row[l] != index_col[l])
        for (auto k = 0u; k < n; ++k)
	  std::swap(a[k][index_row[l]], a[k][index_col[l]]);
  }

} // namespace emsr

#endif // MATRIX_GAUSS_JORDAN_TCC

