/*
gcc -g -o debug_coulomb_gsl debug_coulomb_gsl.c -Lgsl/.libs -Lgsl/cblas/.libs -lgsl -lgslcblas
LD_LIBRARY_PATH=gsl/.libs:gsl/cblas/.libs ./debug_coulomb_gsl 
*/

#include <gsl/gsl_sf_coulomb.h>

int
main()
{
  double lambda = 0.0;
  double eta = -2.0;
  double rho = 0.1;

  const int k_lam_G = 0;
  gsl_sf_result F, Fp, G, Gp;
  double exp_F, exp_G;
  int stat = gsl_sf_coulomb_wave_FG_e(eta, rho, lambda, k_lam_G,
                                      &F, &Fp, &G, &Gp, &exp_F, &exp_G);

  rho = 2.6;
  stat = gsl_sf_coulomb_wave_FG_e(eta, rho, lambda, k_lam_G,
                                  &F, &Fp, &G, &Gp, &exp_F, &exp_G);
}
