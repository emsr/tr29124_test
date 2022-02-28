
#include <iostream>
#include <iomanip>

#include <emsr/matrix.h>
#include <emsr/matrix_util.h>

int
main()
{
  const std::size_t m = 3, n = 2;
  double a_svd[m][n]{{1, 2}, {3, 4}, {5, 6}};
  double v_svd[n][n];
  double w_svd[n];
  emsr::sv_decomp(m, n, a_svd, w_svd, v_svd);
  std::cout << "\n A\n";
  emsr::print_matrix(a_svd);
  std::cout << "\n W\n";
  emsr::print_matrix(w_svd);
  std::cout << "\n V\n";
  emsr::print_matrix(v_svd);

//{{-0.2298477,   0.88346102,  0.40824829},
// {-0.52474482,  0.24078249, -0.81649658},
// {-0.81964194, -0.40189603,  0.40824829}};

  double w_svd_test[2]{9.52551809, 0.51430058};

  double v_svd_test[2][2]{{-0.61962948, -0.78489445},
                         {-0.78489445, 0.61962948}};

  // Try underdetermined system...

  const std::size_t mm = 4, nn = 5;
  double aa_svd[mm][nn]{{1, 0, 0, 0, 2},
                        {0, 0, 3, 0, 0},
                        {0, 0, 0, 0, 0},
                        {0, 2, 0, 0, 0}};
  double vv_svd[nn][nn];
  double ww_svd[nn];
  emsr::sv_decomp(mm, nn, aa_svd, ww_svd, vv_svd);
  std::cout << "\n A\n";
  emsr::print_matrix(aa_svd);
  std::cout << "\n W\n";
  emsr::print_matrix(ww_svd);
  std::cout << "\n V\n";
  emsr::print_matrix(vv_svd);

  double r[mm][nn], iu[3][3], iv[3][3];
  for (int i = 0; i < mm; ++i)
    for (int j = 0; j < nn; ++j)
      {
        r[i][j] = 0.0;
        iu[i][j] = 0.0;
        iv[i][j] = 0.0;
        for (int k = 0; k < nn; ++k)
	  {
            r[i][j] += aa_svd[i][k] * ww_svd[k] * vv_svd[j][k];
            iu[i][j] += aa_svd[k][i] * aa_svd[k][j];
            iv[i][j] += vv_svd[k][i] * vv_svd[k][j];
	  }
      }

  std::cout << "\n Verify U~.U = I:\n";
  emsr::print_matrix(iu);
  std::cout << "\n Verify V~.V = I:\n";
  emsr::print_matrix(iv);
  std::cout << "\n Reconstruction of input matrix from SV decomposition:\n";
  emsr::print_matrix(r);
}
