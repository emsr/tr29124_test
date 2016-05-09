// $HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_jacobi_inv test_jacobi_inv.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_jacobi_inv > test_jacobi_inv.new

// g++ -std=gnu++17 -g -Wall -Wextra -DNO_LOGBQ -I. -o test_jacobi_inv test_jacobi_inv.cpp -lquadmath

// ./test_jacobi_inv > test_jacobi_inv.txt

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

template<typename _Tp>
  _Tp
  jacobi_acn(_Tp __k, _Tp __v)
  { return std::ellint_1(__k, std::acos(__v)); }

template<typename _Tp>
  _Tp
  jacobi_asn(_Tp __k, _Tp __v)
  { return std::ellint_1(__k, std::asin(__v)); }

template<typename _Tp>
  _Tp
  jacobi_adn(_Tp __k, _Tp __v)
  { return std::ellint_1(__k, std::asin(std::sqrt(_Tp{1} - __v * __v) / __k)); }

int
main()
{
  using Tp = double;
  std::cout.precision(std::numeric_limits<Tp>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  auto k = 0.5;

  std::cout << "\n\nJacobi cn inverse\n";
  std::cout << "++++++++++++++++++++++++++++++\n";
  std::cout << std::setw(width) << "y"
	    << std::setw(width) << "cn(k,u)"
	    << std::setw(width) << "acn(k,cn(k,u))"
	      << '\n';
  std::cout << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << '\n';
  for (auto iu = 0; iu < 150; ++iu)
    {
      auto u = Tp(iu * 0.01L);
      auto v = __gnu_cxx::jacobi_cn(k, u);
      auto w = jacobi_acn(k, v);
      std::cout << std::setw(width) << u
		<< std::setw(width) << v
		<< std::setw(width) << std::abs(w - u)
		<< '\n';
    }

  std::cout << "\n\nJacobi sn inverse\n";
  std::cout << "++++++++++++++++++++++++++++++\n";
  std::cout << std::setw(width) << "y"
	    << std::setw(width) << "sn(k,u)"
	    << std::setw(width) << "asn(k,sn(k,u))"
	      << '\n';
  std::cout << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << '\n';
  for (auto iu = 0; iu < 150; ++iu)
    {
      auto u = Tp(iu * 0.01L);
      auto v = __gnu_cxx::jacobi_sn(k, u);
      auto w = jacobi_asn(k, v);
      std::cout << std::setw(width) << u
		<< std::setw(width) << v
		<< std::setw(width) << std::abs(w - u)
		<< '\n';
    }

  std::cout << "\n\nJacobi dn inverse\n";
  std::cout << "++++++++++++++++++++++++++++++\n";
  std::cout << std::setw(width) << "y"
	    << std::setw(width) << "dn(k,u)"
	    << std::setw(width) << "adn(k,dn(k,u))"
	      << '\n';
  std::cout << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << '\n';
  for (auto iu = 0; iu < 150; ++iu)
    {
      auto u = Tp(iu * 0.01L);
      auto v = __gnu_cxx::jacobi_dn(k, u);
      auto w = jacobi_adn(k, v);
      std::cout << std::setw(width) << u
		<< std::setw(width) << v
		<< std::setw(width) << std::abs(w - u)
		<< '\n';
    }
}
