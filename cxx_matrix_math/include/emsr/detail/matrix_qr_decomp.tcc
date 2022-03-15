#ifndef MATRIX_QR_DECOMP_TCC
#define MATRIX_QR_DECOMP_TCC 1

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

namespace emsr
{

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
	    Vector& c, Vector& d, bool& singular)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    NumTp sigma, sum, tau;

    //  Allocate the two vectors.
    //c.resize(n);
    //d.resize(n);

    singular = false;
    for (std::size_t k = 0; k < n_cols - 1; ++k)
      {
	//  See if the matrix is singular in this column.
	NumTp scale = NumTp{0};
	for (std::size_t i = k; i < n_rows; ++i)
	  if (scale < std::abs(a[i][k]))
	    scale = std::abs(a[i][k]);

	if (scale == NumTp{0})
	  {
	    singular = true;
	    c[k] = d[k] = NumTp{0};
	  }
	else
	  {
	    sum = NumTp{0};
	    for (std::size_t i = k; i < n_rows; ++i)
	      {
		a[i][k] /= scale;
		sum += a[i][k] * a[i][k];
	      }
	    sigma = std::sqrt(sum);
	    if (a[k][k] < NumTp{0})
	      sigma = -sigma;
	    a[k][k] += sigma;
	    c[k] = sigma * a[k][k];
	    d[k] = -scale * sigma;
	    for (std::size_t j = k + 1; j < n_cols; ++j)
	      {
		sum = NumTp{0};
		for (std::size_t i = k; i < n_rows; ++i)
		  sum += a[i][k] * a[i][j];
		tau = sum/c[k];
		for (std::size_t i = k; i < n_rows; ++i)
		  a[i][j] -= tau * a[i][k];
	      }
	  }
      }
    c[n_cols - 1] = NumTp{0};
    d[n_cols - 1] = a[n_cols - 1][n_cols - 1];

    if (d[n_cols - 1] == NumTp{0})
      singular = true;

    return;
  }

/**
 * This routine solves the set of equations Rx = b where R is the upper triangular
 * matrix stored in a[0..n_rows - 1][0..n_cols - 1] and d[0..n_cols - 1].
 * Here n_rows >= n_cols.
 */
template<typename Matrix, typename Vector, typename VectorB>
  void
  r_backsub(std::size_t n_rows, std::size_t n_cols,
	    const Matrix& a, const Vector& d,
	    VectorB& b)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    b[n_cols - 1] /= d[n_cols - 1];
    for (int i = n_cols - 2; i >= 0; --i)
      {
	NumTp sum = NumTp{0};
	for (std::size_t j = i + 1; j < n_rows; ++j)
	  sum += a[i][j] * b[j];
	b[i] = (b[i] - sum) / d[i];
      }

    return;
  }

/**
 * This routine solves the set of equations Ax = b.
 * The inputs are the QR decomposition of the matrix in a[0..n_rows - 1][0..n_cols - 1],
 * c[0..n_cols - 1], and d[0..n_cols - 1].
 * The vector b[0..n_rows - 1] is input as the "RHS" and output as the solution.
 * Here n_rows >= n_cols.
 */
template<typename Matrix, typename Vector, typename VectorB>
  void
  qr_backsub(const std::size_t n_rows, const std::size_t n_cols,
	     const Matrix& a, const Vector& c, const Vector& d,
	     VectorB& b)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    //  Form Qt.b.
    for (std::size_t j = 0; j < n_cols - 1; ++j)
      {
	NumTp sum = NumTp{0};
	for (std::size_t i = j; i < n_rows; ++i)
	  sum += a[i][j] * b[i];
	NumTp tau = sum/c[j];
	for (std::size_t i = j; i < n_rows; ++i)
	  b[i] -= tau * a[i][j];
      }

    //  Solve R.x = Qt.b
    r_backsub(n_rows, n_cols, a, d, b);

    return;
  }

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
	    Matrix& a_inv)
  {
    using NumTp = std::decay_t<decltype(a_qr[0][0])>;

    std::vector<NumTp> col(n_rows);

    for (std::size_t j = 0; j < n_rows; ++j)
      {
	for (std::size_t i = 0; i < n_rows; ++i)
	  col[i] = NumTp{0};
	col[j] = NumTp{1};

	qr_backsub(n_rows, n_cols, a_qr, c, d, col);

	for (std::size_t i = 0; i < n_cols; ++i)
	  a_inv[i][j] = col[i];
      }

    return;
  }

/**
 *  Update the QR decomposition.
 */
template<typename Matrix, typename Vector>
  void
  qr_update(std::size_t n_rows, std::size_t n_cols,
	    Matrix& r, Matrix& qt,
	    Vector& u, Vector& v)
  {
    using NumTp = std::decay_t<decltype(qt[0][0])>;

    //  Find largest k such that uk != 0.
    std::ptrdiff_t k;
    for (std::ptrdiff_t k = n_cols - 1; k >= 0; --k)
      if (u[k] != NumTp{0})
	break;

    // Transform R + u x v to upper triangular.
    for (std::ptrdiff_t i = k - 1; i >= 0; --i)
      {
	jacobi_rotate(i, n_rows, n_cols, r, qt, u[i], -u[i + 1]);
	if (u[i] == 0)
	  u[i] = std::abs(u[i + 1]);
	else if (std::abs(u[i]) > std::abs(u[i + 1]))
	  u[i] = std::abs(u[i]) * std::sqrt(NumTp{1} + (u[i + 1] / u[i]) * (u[i + 1] / u[i]));
	else
	  u[i] = std::abs(u[i + 1]) * std::sqrt(NumTp{1} + (u[i] / u[i+1]) * (u[i] / u[i + 1]));
      }
    for (std::size_t j = 0; j < n_cols; ++j)
      r[0][j] += u[0] * v[j];

    // Transform upper Hessenberg matrix to upper triangular.
    for (std::size_t i = 0; i < k; ++i)
      jacobi_rotate(i, n_rows, n_cols, r, qt, r[i][i], -r[i + 1][i]);

    return;
  }

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
		NumTp a, NumTp b)
  {
    NumTp c, s;
    if (a == NumTp{0})
      {
	// Avoid underflow or overflow.
	c = NumTp{0};
	s = NumTp{1};
	if (b < NumTp{0})
	  s = -NumTp{1};
      }
    else
      {
	auto fact = std::hypot(a, b);
	c = a / fact;
	s = b / fact;
      }

    // Multiply r by the Jacobi rotation.
    for (std::size_t j = i; j < n_cols; ++j)
      {
	auto y = r[  i  ][j];
	auto w = r[i + 1][j];
	r[  i  ][j] = c * y - s * w;
	r[i + 1][j] = s * y + c * w;
      }

    // Multiply qt by the Jacobi rotation.
    for (std::size_t j = 0; j < n_rows; ++j)
      {
	auto y = qt[  i  ][j];
	auto w = qt[i + 1][j];
	qt[  i  ][j] = c * y - s * w;
	qt[i + 1][j] = s * y + c * w;
      }

    return;
  }

// Implement the class...

template<typename Matrix>
  qr_decomposition<Matrix>::
  qr_decomposition(std::size_t n_rows, std::size_t n_cols, Matrix& a)
  : m_n_rows(n_rows), m_n_cols(n_cols), m_a(a), m_c(n_cols), m_d(n_cols), m_singular(false)
  {
    qr_decomp(this->m_n_rows, this->m_n_cols,
	      this->m_a, this->m_c, this->m_d, this->m_singular);
  }

template<typename Matrix>
  template<typename Matrix2>
  qr_decomposition<Matrix>::
  qr_decomposition(std::size_t n_rows, std::size_t n_cols, Matrix2& a)
  : m_n_rows(n_rows), m_n_cols(n_cols),
    m_a(n_rows, std::vector<NumTp>(n_cols)),
    m_c(n_cols), m_d(n_cols), m_singular(false)
  {
    // Copy a.
    for (std::size_t i_row = 0; i_row < this->m_n_rows; ++i_row)
      for (std::size_t i_col = 0; i_col < this->m_n_cols; ++i_col)
        this->m_a[i_row][i_col] = a[i_row][i_col];
    qr_decomp(this->m_n_rows, this->m_n_cols,
	      this->m_a, this->m_c, this->m_d, this->m_singular);
  }

template<typename Matrix>
  template<typename Vector2>
  void
  qr_decomposition<Matrix>::
  backsubstitution(Vector2& b)
  {
    qr_backsub(this->m_n_rows, this->m_n_cols,
	       this->m_a, this->m_c, this->m_d,
	       b);
  }

template<typename Matrix>
  template<typename Matrix2>
  void
  qr_decomposition<Matrix>::
  inverse(Matrix2& a_inv)
  {
    qr_invert(this->m_n_rows, this->m_n_cols,
	    this->m_a, this->m_c, this->m_d,
	    a_inv);
  }
/*
template<typename Matrix>
  void
  qr_decomposition<Matrix>::
  update()
  {
    qr_update(this->m_n_rows, this->m_n_cols,
	    Matrix& r, Matrix& qt,
	    Vector& u, Vector& v)
  }
*/
} // namespace emsr

#endif // MATRIX_QR_DECOMP_TCC

