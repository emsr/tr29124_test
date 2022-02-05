/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/numeric_limits.h>
#include <emsr/sf_theta.h>

template<typename Tp>
  void
  test_theta(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = 8 + std::cout.precision();
    std::cout << std::showpoint << std::scientific;

    std::cout << "\n\n Theta function values\n";
    std::cout << " =====================\n";
    const auto del1 = Tp{1} / Tp{10};
    const auto del01 = Tp{1} / Tp{100};
    for (int i = 0; i <= 20; ++i)
      {
	auto nu = i * del1;
	std::cout << '\n' << " nu   = " << std::setw(width) << nu << '\n';
	std::cout << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "theta_1"
		  << ' ' << std::setw(width) << "theta_2"
		  << ' ' << std::setw(width) << "theta_3"
		  << ' ' << std::setw(width) << "theta_4"
		  << '\n';
	std::cout << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << '\n';
	for (int j = 0; j <= 100; ++j)
	  {
	    auto x = j * del01;
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << emsr::theta_1(nu, x)
		      << ' ' << std::setw(width) << emsr::theta_2(nu, x)
		      << ' ' << std::setw(width) << emsr::theta_3(nu, x)
		      << ' ' << std::setw(width) << emsr::theta_4(nu, x)
		      << '\n';
	  }
      }

    std::cout << "\n\n Theta function values and compares with Jacobi elliptic functions\n";
    std::cout << " =================================================================\n";
    auto k = Tp{1} / Tp{3};
    std::cout << '\n' << " k    = " << std::setw(width) << k;
    std::cout << '\n' << " q(k) = " << std::setw(width) << emsr::ellnome(k) << '\n';
    std::cout << ' ' << std::setw(width) << "x"
	      << ' ' << std::setw(width) << "theta_s"
	      << ' ' << std::setw(width) << "theta_c"
	      << ' ' << std::setw(width) << "theta_n"
	      << ' ' << std::setw(width) << "theta_d"
	      << ' ' << std::setw(width) << "th_s/th_n - sn"
	      << ' ' << std::setw(width) << "th_c/th_n - cn"
	      << ' ' << std::setw(width) << "th_d/th_n - dn"
	      << '\n';
    std::cout << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << '\n';
    for (int j = -1000; j <= 1000; ++j)
      {
	auto x = j * del01;
	auto s = emsr::theta_s(k, x);
	auto c = emsr::theta_c(k, x);
	auto d = emsr::theta_d(k, x);
	auto n = emsr::theta_n(k, x);
	auto [sn, cn, dn] = emsr::detail::jacobi_ellint(k, x);
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << s
		  << ' ' << std::setw(width) << c
		  << ' ' << std::setw(width) << d
		  << ' ' << std::setw(width) << n
		  << ' ' << std::setw(width) << s / n - sn
		  << ' ' << std::setw(width) << c / n - cn
		  << ' ' << std::setw(width) << d / n - dn
		  << '\n';
      }

    std::cout << "\n\n Theta function values and compares with Jacobi elliptic functions\n";
    std::cout << " =================================================================\n";
    for (int i = -10; i <= 10; ++i)
      {
	auto k = i * del1;
	std::cout << '\n' << " k    = " << std::setw(width) << k;
	std::cout << '\n' << " q(k) = " << std::setw(width) << emsr::ellnome(k) << '\n';
	std::cout << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "theta_s"
		  << ' ' << std::setw(width) << "theta_c"
		  << ' ' << std::setw(width) << "theta_n"
		  << ' ' << std::setw(width) << "theta_d"
		  << ' ' << std::setw(width) << "th_s/th_n - sn"
		  << ' ' << std::setw(width) << "th_c/th_n - cn"
		  << ' ' << std::setw(width) << "th_d/th_n - dn"
		  << '\n';
	std::cout << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << '\n';
	for (int j = 0; j <= 100; ++j)
	  {
	    auto x = j * del01;
	    auto s = emsr::theta_s(k, x);
	    auto c = emsr::theta_c(k, x);
	    auto d = emsr::theta_d(k, x);
	    auto n = emsr::theta_n(k, x);
	    auto [sn, cn, dn] = emsr::detail::jacobi_ellint(k, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << s
		      << ' ' << std::setw(width) << c
		      << ' ' << std::setw(width) << d
		      << ' ' << std::setw(width) << n
		      << ' ' << std::setw(width) << s / n - sn
		      << ' ' << std::setw(width) << c / n - cn
		      << ' ' << std::setw(width) << d / n - dn
		      << '\n';
	  }
      }
  }

int
main()
{
  std::cout << "\nfloat\n===========\n";
  test_theta<float>();

  std::cout << "\ndouble\n===========\n";
  test_theta<double>();

  std::cout << "\nlong double\n===========\n";
  test_theta<long double>();

  //std::cout << "\nfloat128\n===========\n";
  //test_theta<__float128>();
}
