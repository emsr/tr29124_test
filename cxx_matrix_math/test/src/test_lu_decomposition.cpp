
#include <iostream>
#include <iomanip>

#include <emsr/matrix.h>
#include <emsr/matrix_util.h>

int
main()
{
  const std::size_t m = 3, n = 3;
  double a_lu[m][n]{{1, 2, 3}, {2, -4, 6}, {3, -9, -3}};
  //double a_lu_inv[m][n];
  int index[m];
  double parity;
  emsr::lu_decomp(3, a_lu, index, parity);

  double a_lu_p[m][n]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  for (int i = 0; i < 3; ++i)
    {
      if (index[i] == i)
        continue;
      for (int j = 0; j < i; ++j)
        {
          auto temp = a_lu_p[i][j];
          a_lu_p[i][j] = a_lu_p[index[i]][j];
          a_lu_p[index[i]][j] = temp;
        }
    }
  std::cout << "\n P\n";
  emsr::print_matrix(a_lu_p);

  // Reconstruct factors.
  double a_lu_l[m][n], a_lu_u[m][n];
  for (int i = 0; i < 3; ++i)
    {
      a_lu_l[i][i] = 1.0;
      for (int j = 0; j < i; ++j)
	a_lu_l[i][j] = a_lu[i][j];
      for (int j = i + 1; j < 3; ++j)
	a_lu_l[i][j] = 0.0;

      for (int j = 0; j < i; ++j)
	a_lu_u[i][j] = 0.0;
      for (int j = i; j < 3; ++j)
	a_lu_u[i][j] = a_lu[i][j];
    }
  std::cout << "\n L\n";
  emsr::print_matrix(a_lu_l);
  std::cout << "\n U\n";
  emsr::print_matrix(a_lu_u);

  const double eps = 1.0e-15;
  int num_errors = 0;
  double a_lu_l_test[m][n]{{1.0, 0.0, 0.0}, {1.0/3.0, 1.0, 0.0}, {2.0/3.0, 0.4, 1.0}};
  double a_lu_u_test[m][n]{{3.0, -9.0, -3.0}, {0.0, 5.0, 4.0}, {0.0, 0.0, 6.4}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      {
        if (std::abs(a_lu_l[i][j] - a_lu_l_test[i][j] > eps))
          {
            ++num_errors;
            std::cout << "error: " << a_lu_l[i][j] << " versus " << a_lu_l_test[i][j] << '\n';
          }
        if (std::abs(a_lu_u[i][j] - a_lu_u_test[i][j] > eps))
          {
            ++num_errors;
            std::cout << "error: " << a_lu_u[i][j] << " versus " << a_lu_u_test[i][j] << '\n';
          }
      }
  return num_errors;
}
