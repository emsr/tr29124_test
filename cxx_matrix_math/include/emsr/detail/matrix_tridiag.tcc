#ifndef MATRIX_TRIDIAG_TCC
#define MATRIX_TRIDIAG_TCC 1

#include <cmath>

namespace emsr
{

/**
 * Solves a tridiagonal set of equations where a[1..n] is the subdiagonal vector,
 * b[1..n] is the diagonal vector, and c[1..n] is the superdiagonal vector, and
 * r[1..n] is the right hand side vector.  The solution is x[1..n].
 * a, b, c, and r are not modified.
 */
template<typename Tp>
  void
  tridiagonal(const Tp* a, const Tp* b, const Tp* c,
	      const Tp* r, Tp* x, std::size_t n)
  {
    std::vector<Tp> gam(n);

    if (b[0] == Tp{0})
      std::__throw_runtime_error("First diagonal must be non-zero.");
    Tp bet = b[0];
    x[0] = r[0] / bet;

    // Decomposition and forward substitution.
    for (int j = 1; j < n; ++j)
      {
	gam[j] = c[j - 1] / bet;
	bet = b[j] - a[j] * gam[j];
	if (bet == Tp{0})
	  std::__throw_runtime_error("Matrix must be non-singular.");
	x[j] = (r[j] - a[j] * x[j - 1]) / bet;
      }

    // Backsubstitution.
    for (int j = n - 2; j >= 0; --j)
      x[j] -= gam[j + 1] * x[j + 1];
  }

/**
 * Solves for a vector x[1..n] the cyclic set of linear equations.
 * a[[1..n], b[1..n], c[1..n], and r[1..n] are input vectors of the three diagonal rows and the
 * right side respectively.  alpha and beta are the lower and upper corner entries respectively.
 * The input is not modified.
 */
template<typename Tp>
  void
  cyclic(const Tp* a, const Tp* b, const Tp* c,
	 Tp alpha, Tp beta,
	 const Tp* r, Tp* x, std::size_t n)
  {
    if (n <= 2)
      std::__throw_domain_error("cyclic: n must be greater than 2.");

    const auto gamma = -b[0];

    std::vector<Tp> bb(n);
    bb[0] = b[0] - gamma;
    bb[n - 1] = b[n - 1] - alpha * beta / gamma;
    for (unsigned long i = 1; i < n - 1; ++i)
      bb[0] = b[i];
    tridiagonal(a, bb.data(), c, r, x, n);

    std::vector<Tp> u(n);
    std::vector<Tp> z(n);
    u[0] = gamma;
    u[n - 1] = alpha;
    for (unsigned long i = 1; i < n - 1; ++i)
      u[i] = Tp{0};
    tridiagonal(a, bb, c, u.data(), z.data(), n);

    Tp fact = (x[0] + beta * x[n - 1] / gamma)
		/ (1.0 + z[0] + beta * z[n - 1] / gamma);
    for (unsigned long i = 0; i < n; ++i)
      x[i] -= fact * z[i];
  }

} // namespace emsr

#endif // MATRIX_TRIDIAG_TCC

