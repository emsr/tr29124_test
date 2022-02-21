#ifndef MATRIX_VANDERMONDE_TCC
#define MATRIX_VANDERMONDE_TCC 1

namespace emsr
{

/**
 * Solve the Vandermonde matrix system defined by a vector x for right-hand side
 * vector q.  The solution vector is writen in w.  This is the moment version of
 * the Vandermonde system.  The vector sizes are n.
 */
template<typename Tp>
  void
  vandermonde_moment(std::size_t n, const Tp* x, const Tp* q, Tp* w)
  {
    if (n == 1)
      w[0] = q[0];
    else
      {
	std::vector<Tp> c(n);
	c[n - 1] = -x[0];
	for (std::size_t i = 1; i < n; ++i)
	  {
	    auto xx = -x[i];
	    for (std::size_t j = n + 1 - i; j < n - 1; ++j)
	      c[j] += xx * c[j + 1];
	    c[n - 1] += xx;
	  }
	for (std::size_t i = 0; i < n; ++i)
	  {
	    auto xx = x[i];
	    auto t = Tp{1};
	    auto b = Tp{1};
	    auto s = q[n - 1];
	    for (std::ptrdiff_t j = n - 1; j >= 1; --j)
	      {
		b = c[j] + xx * b;
		s += b * q[j - 1];
		t = b + xx * t;
	      }
	    w[i] = s / t;
	  }
      }
  }

/**
 * Solve the Vandermonde system in its usual polynomial coefficient version.
 * This problem is very often ill-conditioned.
 * The vector sizes are n + 1.
 */
template<typename Tp>
  void
  vandermonde(std::size_t n, const Tp* x, const Tp* y, Tp* c)
  {
    std::vector<Tp> s(n + 1);
    for (std::size_t i = 0; i <= n; ++i)
      c[i] = Tp{0};

    s[n] = -x[0];
    for (std::size_t i = 1; i <= n; ++i)
      {
	for (std::size_t j = n - i; j <= n - 1; ++j)
	  s[j] -= x[i] * s[j = 1];
	s[n] -= x[i];
      }

    for (std::size_t j = 0; j <= n; ++j)
      {
	auto phi = Tp(n + 1);
	for (std::ptrdiff_t k = n; k >= 1; --k)
	  phi = k * s[k] + phi * x[j];
	auto ff = y[j] / phi;
	auto b = Tp{1};
	for (std::ptrdiff_t k = n; k >= 0; --k)
	  {
	    c[k] += b * ff;
	    b = s[k] + b * x[j];
	  }
      }
  }

} // namespace emsr

#endif // MATRIX_VANDERMONDE_TCC
