#ifndef MATRIX_SV_DECOMP_TCC
#define MATRIX_SV_DECOMP_TCC

#include <cstdlib>
#include <sstream>
#include <vector>
#include <cmath>

namespace emsr
{

/**
 *  
 */
template<typename Matrix, typename Vector>
  void
  sv_decomp(const std::size_t n_rows, const std::size_t n_cols,
	    Matrix& a, Vector& w, Matrix& v)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    int l;
    NumTp c, f, h, s;

    const int max_iters = 30;

    std::vector<NumTp> rv1(n_cols);

    NumTp g = NumTp{0};
    NumTp scale = NumTp{0};
    NumTp anorm = NumTp{0};

    //  Householder reduction to bidiagonal form.
    for (std::size_t i = 0; i < n_cols; ++i)
      {
	auto l = i + 1;
	rv1[i] = scale * g;
	g = s = scale = NumTp{0};
	if (i <= n_rows - 1)
	  {
	    for (std::size_t k = i; k < n_rows; ++k)
	      scale += std::abs(a[k][i]);
	    if (scale)
	      {
		for (std::size_t k = i; k < n_rows; ++k)
		  {
		    a[k][i] /= scale;
		    s += a[k][i] * a[k][i];
		  }
		f = a[i][i];
		g = -std::copysign(std::sqrt(s), f);
		h = f * g - s;
		a[i][i] = f - g;
		for (std::size_t j = l; j < n_cols; ++j)
		  {
		    s = NumTp{0};
		    for (std::size_t k = i; k < n_rows; ++k)
		      s += a[k][i] * a[k][j];
		    f = s / h;
		    for (std::size_t k = i; k < n_rows; ++k)
		      a[k][j] += f * a[k][i];
		  }
		for (std::size_t k = i; k < n_rows; ++k)
		  a[k][i] *= scale;
	      }
	  }
	w[i] = scale * g;
	g = s = scale = NumTp{0};
	if (i <= n_rows - 1 && i != n_cols - 1)
	  {
	    for (std::size_t k = l; k < n_cols; ++k)
	      scale += std::abs(a[i][k]);
	    if (scale)
	      {
		for (std::size_t k = l; k < n_cols; ++k)
		  {
		    a[i][k] /= scale;
		    s += a[i][k] * a[i][k];
		  }
		f = a[i][l];
		g = -std::copysign(std::sqrt(s), f);
		h = f * g - s;
		a[i][l] = f - g;
		for (std::size_t k = l; k < n_cols; ++k)
		  rv1[k] = a[i][k] / h;
		for (std::size_t j = l; j < n_rows; ++j)
		  {
		    s = NumTp{0};
		    for (std::size_t k = l; k < n_cols; ++k)
		      s += a[j][k] * a[i][k];
		    for (std::size_t k = l; k < n_cols; ++k)
		      a[j][k] += s * rv1[k];
		  }
		for (std::size_t k = l; k < n_cols; ++k)
		  a[i][k] *= scale;
	      }
	  }
	anorm = std::max(anorm, std::abs(w[i]) + std::abs(rv1[i]));
      }

    //  Accumulation of right-hand decomposition V.
    for (std::ptrdiff_t i = n_cols - 1; i >= 0; --i)
      {
	if (i < std::ptrdiff_t(n_cols - 1))
	  {
	    if (g)
	      {
		for (std::size_t j = l; j < n_cols; ++j)
		  v[j][i] = (a[i][j]/a[i][l])/g;
		for (std::size_t j = l; j < n_cols; ++j)
		  {
		    s = NumTp{0};
		    for (std::size_t k = l; k < n_cols; ++k)
		      s += a[i][k] * v[k][j];
		    for (std::size_t k = l; k < n_cols; ++k)
		      v[k][j] += s * v[k][i];
		  }
	      }
	    for (std::size_t j = l; j < n_cols; ++j)
	      v[i][j] = v[j][i] = NumTp{0};
	  }
	v[i][i] = NumTp{1};
	g = rv1[i];
	l = i;
      }

    //  Accumulation of left-hand decompositions.
    for (std::ptrdiff_t i = std::min(n_rows, n_cols) - 1; i >= 0; --i)
      {
	l = i + 1;
	g = w[i];
	for (std::size_t j = l; j < n_cols; ++j)
	  a[i][j] = NumTp{0};
	if (g)
	  {
	    g = NumTp{1} / g;
	    for (std::size_t j = l; j < n_cols; ++j)
	      {
		s = NumTp{0};
		for (std::size_t k = l; k < n_rows; ++k)
		  s += a[k][i] * a[k][j];
		f = (s / a[i][i]) * g;
		for (std::size_t k = i; k < n_rows; ++k)
		  a[k][j] += f * a[k][i];
	      }
	    for (std::size_t j = i; j < n_rows; ++j)
	      a[j][i] *= g;
	  }
	else
	  for (std::size_t j = i; j < n_rows; ++j)
	    a[j][i] = NumTp{0};
	++a[i][i];
      }

    //  Diagonalization of the bidiagonal form;
    for (std::ptrdiff_t k = n_cols - 1; k >= 0; --k)
      {
	for (std::size_t iter = 1; iter <= max_iters; ++iter)
	  {
	    bool flag = true;
	    std::ptrdiff_t l, nm;
	    for (l = k; l >= 0; --l)
	      {
		nm = l - 1;
		if (std::abs(rv1[l]) + anorm == anorm)
		  {
		    flag = false;
		    break;
		  }
		if (std::abs(w[nm]) + anorm == anorm)
		  break;
	      }
	    if (flag)
	      {
		auto c = NumTp{0};
		auto s = NumTp{1};
		for (std::ptrdiff_t i = l; i < k; ++i)
		  {
		    const auto f = s * rv1[i];
		    rv1[i] = c * rv1[i];
		    if (std::abs(f) + anorm == anorm)
		      break;
		    auto g = w[i];
		    auto h = std::hypot(f, g);
		    w[i] = h;
		    h = NumTp{1} / h;
		    c = g * h;
		    s = -f * h;
		    for (std::size_t j = 0; j < n_rows; ++j)
		      {
			auto y = a[j][nm];
			auto z = a[j][i];
			a[j][nm] = y * c + z * s;
			a[j][i] = z * c - y * s;
		      }
		  }
	      }
	    auto z = w[k];
	    if (l == k)
	      {
		//  Convergence!!!
		if (z < NumTp{0})
		  {
		    //  Make singular value non negative.
		    w[k] = -z;
		    for (std::size_t j = 0; j < n_cols; ++j)
		      v[j][k] = -v[j][k];
		  }
		break;
	      }
	    if (iter == max_iters)
	      {
		std::ostringstream msg;
		msg << "sv_decomp: No convergence in " << max_iters << " iterations.";
		throw std::logic_error(msg.str().c_str());
	      }

	    //  Shift from bottom 2x2 minor.
	    auto x = w[l];
	    nm = k - 1;
	    auto y = w[nm];
	    g = rv1[nm];
	    h = rv1[k];
	    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
	    g = std::hypot(f, NumTp{1});
	    f = ((x - z) * (x + z) + h * ((y / (f + std::copysign(g, f))) - h)) / x;

	    //  Next QR transformation.
	    c = s = NumTp{1};
	    for (std::ptrdiff_t j = l; j <= nm; ++j)
	      {
		auto i = j + 1;
		g = rv1[i];
		y = w[i];
		h = s * g;
		g = c * g;
		z = std::hypot(f, h);
		rv1[j] = z;
		c = f / z;
		s = h / z;
		f = x * c + g * s;
		g = g * c - x * s;
		h = y * s;
		y *= c;
		for (std::size_t jj = 0; jj < n_cols; ++jj)
		  {
		    x = v[jj][j];
		    z = v[jj][i];
		    v[jj][j] = x * c + z * s;
		    v[jj][i] = z * c - x * s;
		  }
		z = std::hypot(f, h);
		w[j] = z;
		//  Rotation can be arbitrary if z = 0.
		if (z)
		  {
		    z = NumTp{1} / z;
		    c = f * z;
		    s = h * z;
		  }
		f = c * g + s * y;
		x = c * y - s * g;
		for (std::size_t jj = 0; jj < n_rows; ++jj)
		  {
		    auto y = a[jj][j];
		    auto z = a[jj][i];
		    a[jj][j] = y * c + z * s;
		    a[jj][i] = z * c - y * s;
		  }
	      }
	    rv1[l] = NumTp{0};
	    rv1[k] = f;
	    w[k] = x;
	  }
      }

    return;
  }

/**
 *  
 */
template<typename Matrix, typename Vector>
  void
  sv_backsub(std::size_t n_rows, std::size_t n_cols,
	     const Matrix& u,
	     const Vector& w, const Matrix& v,
	     const Vector& b, Vector& x)
  {
    using NumTp = std::decay_t<decltype(u[0][0])>;

    std::vector<NumTp> tmp(n_cols);

    for (std::size_t j = 0; j < n_cols; ++j)
      {
	NumTp s = NumTp{0};
	if (w[j] != NumTp{0})
	  {
	    for (std::size_t i = 0; i < n_rows; ++i)
	      s += u[i][j] * b[i];
	    s /= w[j];
	  }
	tmp[j] = s;
      }
    for (std::size_t j = 0; j < n_cols; ++j)
      {
	NumTp s = NumTp{0};
	for (std::size_t jj = 0; jj < n_cols; ++jj)
	  s += v[j][jj] * tmp[jj];
	x[j] = s;
      }

    return;
  }

/**
 *  Improves a solution vector x of the linear set A.x = b.
 *  The Matrix a and the SV decomposition of a -- u, w, v and the
 *  right-hand side Vector are input along with the solution vector x.
 *  The solution vector x is improved and modified on output.
 */
template<typename Matrix, typename Vector>
  void
  sv_improve(std::size_t n_rows, std::size_t n_cols,
	     const Matrix& a, const Matrix& u,
	     const Vector& w, const Matrix& v,
	     const Vector& b, Vector& x)
  {
    using NumTp = std::decay_t<decltype(a[0][0])>;

    std::vector<NumTp> r(n_rows);
    std::vector<NumTp> dx(n_cols);

    for (std::size_t i = 0; i < n_rows; ++i)
      {
	r[i] = -b[i];
	for (std::size_t j = 0; j < n_cols; ++j)
	  r[i] += a[i][j] * x[j];
      }

    sv_backsub(n_rows, n_cols, u, w, v, r, dx);

    for (std::size_t i = 0; i < n_cols; ++i)
      x[i] -= dx[i];

    return;
  }

// Implement the class...

template<typename Matrix>
  sv_decomposition<Matrix>::
  sv_decomposition(std::size_t n_rows, std::size_t n_cols, const Matrix& a)
  : m_n_rows(n_rows), m_n_cols(n_cols),
    m_a(a), m_w(n_cols), m_v(n_cols, std::vector<NumTp>(n_cols))
  {
    sv_decomp(this->m_n_rows, this->m_n_cols,
	      this->m_a, this->m_w, this->m_v);
  }

template<typename Matrix>
  template<typename Matrix2>
  sv_decomposition<Matrix>::
  sv_decomposition(std::size_t n_rows, std::size_t n_cols, const Matrix2& a)
  : m_n_rows(n_rows), m_n_cols(n_cols),
    m_a(n_rows, std::vector<NumTp>(n_cols)),
    m_w(n_cols), m_v(n_cols, std::vector<NumTp>(n_cols))
  {
    // Copy a.
    for (std::size_t i_row = 0; i_row < this->m_n_rows; ++i_row)
      for (std::size_t i_col = 0; i_col < this->m_n_cols; ++i_col)
        this->m_a[i_row][i_col] = a[i_row][i_col];
    sv_decomp(this->m_n_rows, this->m_n_cols,
	      this->m_a, this->m_w, this->m_v);
  }

template<typename Matrix>
  template<typename Vector2, typename VectorOut>
  void
  sv_decomposition<Matrix>::
  backsubstitution(const Vector2& b, VectorOut& x) const
  {
    sv_backsub(this->m_n_rows, this->m_n_cols,
	       this->m_a, this->m_w, this->m_v,
	       b, x);
  }

template<typename Matrix>
  template<typename Matrix2, typename Vector2, typename VectorOut>
  void
  sv_decomposition<Matrix>::
  improve(const Matrix2 a_orig, const Vector2& b, VectorOut& x) const
  {
    sv_improve(this->m_n_rows, this->m_n_cols,
	       a_orig, this->m_a, this->m_w, this->m_v,
	       b, x);
  }

} // namespace emsr

#endif // MATRIX_SV_DECOMP_TCC
