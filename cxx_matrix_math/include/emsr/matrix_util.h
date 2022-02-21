#if ! defined(MATRIX_UTIL_H)
#define MATRIX_UTIL_H

#include <iosfwd>

namespace emsr
{

  template<typename Numeric>
    struct promote
    {
      typedef Numeric type;
    };

  template<>
    struct promote<float>
    {
      typedef double type;
    };

  template<>
    struct promote<double>
    {
      typedef long double type;
    };

  template<typename Numeric>
    using promote_t = typename promote<Numeric>::type;

  template<typename Numeric, std::size_t M, std::size_t N>
    void
    copy_matrix(Numeric (&mat)[M][N], const Numeric (&mat_in)[M][N])
    {
      for (std::size_t i = 0; i < M; ++i)
	for (std::size_t j = 0; j < N; ++j)
	  mat[i][j] = mat_in[i][j];
    }

  template<typename Numeric, std::size_t M, std::size_t K, std::size_t N>
    void
    mul_matrix(Numeric (&c)[M][N], const Numeric (&a)[M][K], const Numeric (&b)[K][N])
    {
      for (std::size_t i = 0; i < M; ++i)
	for (std::size_t j = 0; j < N; ++j)
	  {
	    c[i][j] = Numeric{0};
	    for (std::size_t k = 0; k < K; ++k)
	      c[i][j] += a[i][k] * b[k][j];
	  }
    }

  template<typename Numeric, std::size_t M, std::size_t K>
    void
    mul_matrix(Numeric (&c)[M], const Numeric (&a)[M][K], const Numeric (&b)[K])
    {
      for (int i = 0; i < M; ++i)
	{
	  c[i] = Numeric{0};
	  for (std::size_t k = 0; k < K; ++k)
	    c[i] += a[i][k] * b[k];
	}
    }

  template<typename Numeric, std::size_t M, std::size_t N>
    void
    print_matrix(const Numeric (&mat)[M][N])
    {
      for (auto& row : mat)
	{
	  for (auto& col : row)
	    std::cout << ' ' << std::setw(10) << col;
	  std::cout << '\n';
	}
    }

  template<typename Numeric, std::size_t M>
    void
    print_matrix(const Numeric (&mat)[M])
    {
      for (auto& row : mat)
	std::cout << ' ' << std::setw(10) << row;
      std::cout << '\n';
    }

}

#endif // MATRIX_UTIL_H
