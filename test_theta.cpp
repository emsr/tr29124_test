// $HOME/bin_specfun/bin/g++ -std=gnu++1z -o test_theta test_theta.cpp

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_theta > test_theta.txt

#include <cmath>

template<typename _Tp>
  void
  test_theta()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;

    std::cout << "\n\n Theta function values\n";
    std::cout << " =====================\n";
    std::cout << ' ' << std::setw(width) << "x"
	      << ' ' << std::setw(width) << "theta_1"
	      << ' ' << std::setw(width) << "theta_2"
	      << ' ' << std::setw(width) << "theta_3"
	      << ' ' << std::setw(width) << "theta_4"
	      << '\n';
    for (int i = 0; i <= 20; ++i)
      {
	auto nu = 0.1 * i;
	std::cout << '\n' << ' ' << std::setw(width) << nu << '\n';
	for (int j = 0; j < 100; ++j)
	  {
	    auto k = j * 0.01;
	    auto q = __gnu_cxx::ellnome(k);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << __gnu_cxx::theta_1(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::theta_2(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::theta_3(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::theta_4(nu, x)
		      << '\n';
	  }
      }

    std::cout << "\n\n Theta function values and compares with Jacobi elliptic functions\n";
    std::cout << " =================================================================\n";
    std::cout << ' ' << std::setw(width) << "x"
	      << ' ' << std::setw(width) << "theta_s"
	      << ' ' << std::setw(width) << "theta_c"
	      << ' ' << std::setw(width) << "theta_n"
	      << ' ' << std::setw(width) << "theta_d"
	      << ' ' << std::setw(width) << "th_s/th_n - sn"
	      << ' ' << std::setw(width) << "th_c/th_n - cn"
	      << ' ' << std::setw(width) << "th_d/th_n - dn"
	      << '\n';
    for (int i = 0; i <= 20; ++i)
      {
	auto k = 0.1 * i;
	std::cout << '\n' << ' ' << std::setw(width) << k << '\n';
	for (int j = 0; j < 100; ++j)
	  {
	    auto x = j * 0.01;
	    auto s = __gnu_cxx::theta_s(k, x);
	    auto c = __gnu_cxx::theta_c(k, x);
	    auto d = __gnu_cxx::theta_d(k, x);
	    auto n = __gnu_cxx::theta_n(k, x);
	    auto sncndn = std::__detail::__jacobi_sncndn(k, x);
	    std::cout << ' ' << std::setw() << k
		      << ' ' << x
		      << ' ' << s
		      << ' ' << c
		      << ' ' << d
		      << ' ' << n
		      << ' ' << s / n - std::get<0>(sncndn)
		      << ' ' << c / n - std::get<1>(sncndn)
		      << ' ' << d / n - std::get<2>(sncndn)
		      << '\n';
	  }
      }
  }

int
main()
{
  std::cout << "\nlong double\n===========\n";
  test_theta<long double>();
}
