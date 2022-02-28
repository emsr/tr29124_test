
#include <iostream>
#include <iomanip>

#include <emsr/matrix.h>
#include <emsr/matrix_util.h>

int
main()
{
  const std::size_t m = 3, n = 3;
  double a_c[m][n]{{1, 3, 5}, {3, 13, 23}, {5, 23, 42}};
  std::cout << "\n Input matrix for Cholesky decomposition:\n";
  emsr::print_matrix(a_c);

  double d_c[3];
  emsr::cholesky_decomp(3, a_c, d_c);

  std::cout << "\n Output matrix of Cholesky decompostion:\n";
  emsr::print_matrix(a_c);

  std::cout << "\n Output diagonal vector of Cholesky decompostion:\n";
  emsr::print_matrix(d_c);

  decltype(a_c) a_c_l;
  std::cout << "\n Lower triangular matrix (with zero upper triangle):\n";
  for (int i = 0; i < 3; ++i)
    {
      a_c_l[i][i] = d_c[i];
      for (int j = 0; j < i; ++j)
	a_c_l[i][j] = a_c[i][j];
      for (int j = i + 1; j < 3; ++j)
	a_c_l[i][j] = 0.0;
    }
  emsr::print_matrix(a_c_l);

  const double eps = 1.0e-15;
  int num_errors = 0;
  double a_c_l_test[m][n]{{1, 0, 0}, {3, 2, 0}, {5, 4, 1}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (std::abs(a_c_l[i][j] - a_c_l_test[i][j] > eps))
        {
          ++num_errors;
          std::cout << "error: " << a_c_l[i][j] << " versus " << a_c_l_test[i][j] << '\n';
        }
  return num_errors;
}
