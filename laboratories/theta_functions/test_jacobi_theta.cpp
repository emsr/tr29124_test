/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/numeric_limits.h>
#include <emsr/sf_theta.h>

/**
 * Plot the thetas over x= [0, 2pi] for q = 0, 0.1, ..., 0.9
 */
template<typename Tp>
  void
  test_jacobi_theta(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto s_pi = emsr::pi_v<Tp>;

    const auto del1 = Tp{1} / Tp{10};
    const auto del01 = Tp{1} / Tp{100};
    for (int i = 0; i <= 9; ++i)
      {
	auto q = i * del1;
	std::cout << '\n' << '\n' << " q   = " << std::setw(w) << q << '\n';
	std::cout << ' ' << std::setw(w) << "x"
		  << ' ' << std::setw(w) << "jacobi_theta_1"
		  << ' ' << std::setw(w) << "jacobi_theta_2"
		  << ' ' << std::setw(w) << "jacobi_theta_3"
		  << ' ' << std::setw(w) << "jacobi_theta_4"
		  << '\n';
	std::cout << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << '\n';
	for (int j = 0; j <= 200; ++j)
	  {
	    auto x = j * del01;
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << emsr::detail::jacobi_theta_1(q, s_pi * x)
		      << ' ' << std::setw(w) << emsr::detail::jacobi_theta_2(q, s_pi * x)
		      << ' ' << std::setw(w) << emsr::detail::jacobi_theta_3(q, s_pi * x)
		      << ' ' << std::setw(w) << emsr::detail::jacobi_theta_4(q, s_pi * x)
		      << '\n';
	  }
      }
    std::cout.flush();
  }

/**
 * Reproduce DLMF plots.
 */
template<typename Tp>
  void
  plot_jacobi_theta(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto s_pi = emsr::pi_v<Tp>;

    std::cout << "\n\n\n";
    Tp q0 = Tp{0.15L};
    std::cout << "q = " << q0 << '\n';
    for (int i = 0; i <= 200; ++i)
      {
	auto x = Tp(i) / Tp{100};
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << emsr::detail::jacobi_theta_1(q0, s_pi * x)
		  << ' ' << std::setw(w) << emsr::detail::jacobi_theta_2(q0, s_pi * x)
		  << ' ' << std::setw(w) << emsr::detail::jacobi_theta_3(q0, s_pi * x)
		  << ' ' << std::setw(w) << emsr::detail::jacobi_theta_4(q0, s_pi * x)
		  << '\n';
      }
    std::cout.flush();

    // Vary x for fixed q.
    int nx = 400;
    auto delx = 2 / Tp(nx);

    std::cout << "\n\n\n";
    std::cout << "theta_1; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {Tp{0.05L}, Tp{0.5L}, Tp{0.7L}, Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_1(q, s_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_2; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {Tp{0.05L}, Tp{0.5L}, Tp{0.7L}, Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_2(q, s_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_3; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {Tp{0.05L}, Tp{0.5L}, Tp{0.7L}, Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_3(q, s_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_4; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {Tp{0.05L}, Tp{0.5L}, Tp{0.7L}, Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_4(q, s_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    // Vary q for fixed x.
    int nq = 200;
    auto delq = 1 / Tp(nq);

    std::cout << "\n\n\n";
    std::cout << "theta_1; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {Tp{0.0L}, Tp{0.4L}, Tp{5.0L}, Tp{10.0L}, Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_1(q, x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_2; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {Tp{0.0L}, Tp{0.4L}, Tp{5.0L}, Tp{10.0L}, Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_2(q, x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_3; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {Tp{0.0L}, Tp{0.4L}, Tp{5.0L}, Tp{10.0L}, Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_3(q, x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_4; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {Tp{0.0L}, Tp{0.4L}, Tp{5.0L}, Tp{10.0L}, Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << emsr::detail::jacobi_theta_4(q, x);
	std::cout << '\n';
      }
    std::cout.flush();
  }


int
main()
{
  emsr::detail::jacobi_theta_1(0.205, 40.0);
  emsr::detail::jacobi_theta_1(0.210, 40.0);

  plot_jacobi_theta(1.0);

  test_jacobi_theta(1.0);
}
