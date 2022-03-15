
// See NRiC about rational Chebyshev fits.
// I was hoping that generalized numer and denom polynomials
// would drive toward a tiny b_0 to model a 1/x. But alas, no.
// The denom u and/or bb might be the problem.

// What I actualy solved was p(x)/[q(x) - 1]
// because b_0 appears on both the left and rhs.
// The original set b_0 to 1 and it was the rhs.
// I combobulated the answer but this gets into the loss of precision 1 +/- small.

// Crazy idea: these fitters look at absolute errors. What about relative errors?
// What happens if f(x) is at or near 0? Probably sadness.
// OTOH, that Canadian paper has a weight function with the two most useful choices:
// w(x) = 1 for absolute error
// w(x) = f(x) for relative error.
// They don't mention |f(x)| for some reason.
// If w(x) = f(x) is near 0 then you could eliminate the point - downdate it somehow.
// OTOH, the initial function evaluations are designed to hit _extrema_.

#include <iostream>
#include <vector>
#include <numbers>
#include <cmath>
#include <limits>

#include <emsr/matrix_sv_decomp.h>
#include <emsr/rational_polynomial.h>

template<typename Real, typename Func, typename FuncWt>
void
rational_fit(Func func, FuncWt weight, Real a, Real b, std::size_t deg_numer, std::size_t deg_denom,
             std::vector<Real>& coeff, Real& dev)
{
  const std::size_t num_pts_per_coeff = 8;
  const std::size_t max_iters = 5;
  const Real huge = std::numeric_limits<Real>::max();
  const Real pi_2 = std::numbers::pi_v<Real> / Real{2};

  std::size_t num_coeffs = deg_numer + deg_denom + 2;
  std::size_t num_points = num_pts_per_coeff * num_coeffs;
  std::vector<Real> err(num_points);
  std::vector<Real> fs(num_points);
  std::vector<Real> wts(num_points);
  std::vector<Real> wt(num_points);
  std::vector<Real> xs(num_points);

  dev = huge;

  // Minimax point of the Chebyshev polynomial of the first kind of order num_poimts - 1.
  for (std::size_t i = 0; i < num_points; ++i)
  {
    if (i < num_points / 2)
      {
        auto theta = pi_2 * i / (num_points - 1);
        auto sintheta = std::sin(theta);
        xs[i] = a + (b - a) * sintheta * sintheta;
      }
    else
      {
        auto theta = pi_2 * (num_points - i - 1) / (num_points - 1);
        auto sintheta = std::sin(theta);
        xs[i] = b - (b - a) * sintheta * sintheta;
      }
    fs[i] = func(xs[i]);
    wts[i] = weight(xs[i]);
    err[i] = Real{1} / wts[i];
    wt[i] = std::abs(err[i]);
  }

  using MatTp = std::vector<std::vector<Real>>;

  auto e = Real{0};
  for (std::size_t iter = 1; iter <= max_iters; ++iter)
    {
      std::vector<Real> bb(num_points);
      MatTp u(num_points, std::vector<Real>(num_coeffs));
      for (std::size_t i = 0; i < num_points; ++i)
        {
          auto power = wt[i];
          bb[i] = power * (fs[i] + std::copysign(e, err[i]));
          for (std::size_t j = 0; j < deg_numer + 1; ++j)
            {
              u[i][j] = power;
              power *= xs[i];
            }
          power = -bb[i];
          for (std::size_t j = deg_numer + 1; j < num_coeffs; ++j)
            {
              u[i][j] = power;
              power *= xs[i];
            }
        }
      std::vector<Real> coeff_temp(num_coeffs);
      emsr::sv_decomposition<MatTp> svd(num_points, num_coeffs, u);
      svd.backsubstitution(bb, coeff_temp);
      auto dev_max = Real{0};
      auto sum = Real{0};
      emsr::Polynomial numer(coeff_temp.begin(), coeff_temp.begin() + deg_numer + 1),
                       denom(coeff_temp.begin() + deg_numer + 1, coeff_temp.begin() + deg_numer + deg_denom + 2);
      std::cout << "# numer: " << numer << '\n';
      std::cout << "# denom: " << denom << '\n';
      emsr::RationalPolynomial rat(numer, denom);
      for (std::size_t j = 0; j < num_points; ++j)
        {
          //err[j] = rat(xs[j]) - fs[j];
          //err[j] = rat.numer(xs[j]) / (Real{0} + xs[j] * rat.denom(xs[j])) - fs[j];
          // This thing actually solves for b_0 - 1.
          err[j] = rat.numer(xs[j]) / (Real{1} + rat.denom(xs[j])) - fs[j];
          err[j] /= wts[j];
          wt[j] = std::abs(err[j]);
          sum += wt[j];
          if (wt[j] > dev_max)
            dev_max = wt[j];
        }
      e = sum / num_points;
      if (dev_max <= dev)
        {
          dev = dev_max;
          for (std::size_t j = 0; j < num_coeffs; ++j)
            coeff[j] = coeff_temp[j];
        }

      std::cout << "\n\n";
      for (int i = 1; i <= 100; ++i)
        {
          auto x = Real{0.01} * i * std::numbers::pi_v<Real>;
          //auto r = rat(x);
          //auto r = rat.numer(x) / (Real{0.0} + x * rat.denom(x));
          auto r = rat.numer(x) / (Real{1} + rat.denom(x));
          auto f = func(x);
          std::cout << ' ' << x << ' ' << r << ' ' << f << ' ' << r - f << '\n';
        }
    }
}

int
main()
{
  const auto pi = std::numbers::pi_v<double>;

  auto fun = [](double x) -> double { return std::cos(x)/(1.0 - std::exp(x)); };
  // Possible weight functions
  auto one = [](double) -> double { return 1.0; };
  auto rel = [fun](double x) -> double { return fun(x); };

  double dev;
  std::size_t deg_numer = 4;
  std::size_t deg_denom = 4;
  std::vector<double> coeff(deg_numer + deg_denom + 2);

  rational_fit(fun, rel, 0.01, pi, deg_numer, deg_denom, coeff, dev);

  std::cout << "\n\n dev = " << dev << '\n';
  for (std::size_t k = 0; k < deg_numer + 1; ++k)
    std::cout << ' ' << coeff[k];
  std::cout << '\n';
  for (std::size_t k = deg_numer + 1; k < deg_numer + deg_denom + 2; ++k)
    std::cout << ' ' << coeff[k];
  std::cout << '\n';

  emsr::Polynomial<double> numer(coeff.begin(), coeff.begin() + deg_numer + 1),
                           denom(coeff.begin() + deg_numer + 1, coeff.begin() + deg_numer + deg_denom + 2);
  emsr::RationalPolynomial rat(numer, denom);
  std::cout << "\n\n";
  for (int i = 1; i <= 100; ++i)
    {
      auto x = 0.01 * i * pi;
      //auto r = rat(x);
      //auto r = rat.numer(x) / (0.0 + x * rat.denom(x));
      auto r = rat.numer(x) / (1.0 + rat.denom(x));
      auto f = fun(x);
      std::cout << ' ' << x << ' ' << r << ' ' << f << ' ' << r - f << '\n';
    }

  // Look for extrema..
  for (int i = 2; i < 100; ++i)
    {
      auto xm = 0.01 * (i - 1) * pi;
      auto x = 0.01 * i * pi;
      auto xp = 0.01 * (i + 1) * pi;
      auto delm = rat.numer(xm) / (1.0 + rat.denom(xm)) - fun(xm);
      auto del = rat.numer(x) / (1.0 + rat.denom(x)) - fun(x);
      auto delp = rat.numer(xp) / (1.0 + rat.denom(xp)) - fun(xp);
      if (del > delm && del > delp)
        std::cout << "max: x = " << x << ";  del = " << del << '\n';
      if (del < delm && del < delp)
        std::cout << "min: x = " << x << ";  del = " << del << '\n';
    }
}
