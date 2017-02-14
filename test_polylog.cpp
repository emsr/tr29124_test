/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_polylog test_polylog.cpp -lquadmath -L. -lwrap_cephes
LD_LIBRARY_PATH=.:cephes:$LD_LIBRARY_PATH ./test_polylog > test_polylog.txt

LD_LIBRARY_PATH=.:cephes:$LD_LIBRARY_PATH $HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_polylog test_polylog.cpp -lquadmath -L. -lwrap_cephes

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_polylog test_polylog.cpp -lquadmath -L. -lwrap_cephes
LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH ./test_polylog > test_polylog.txt
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <bits/float128_io.h>

#include "wrap_cephes.h"

  /**
   * Stolen from test_bernoulli.cpp
   * Use for nonpositive integer order.
   */
  template<typename _Tp>
    _Tp
    __stirling_2_recur(unsigned int __n, unsigned int __m)
    {
      if (__n == 0)
	return _Tp(__m == 0);
      else if (__m == 0)
	return _Tp(__n == 0);
      else
	{
	  std::vector<_Tp> __sigold(__m + 1), __signew(__m + 1);
	  __sigold[1] = _Tp{1};
	  if (__n == 1)
	    return __sigold[__m];
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      __signew[1] = __sigold[1];
	      for (auto __im = 2u; __im <= __m; ++__im)
		__signew[__im] = __im * __sigold[__im] + __sigold[__im - 1];
	      std::swap(__sigold, __signew);
	    }
	  return __signew[__m];
	}
    }

  template<typename _Tp>
    _Tp
    __polylog_nonpos_int(int __n, _Tp __x)
    {
      const auto __arg = _Tp{-1} / (_Tp{1} - __x);
      auto __term  =__arg * __stirling_2_recur<_Tp>(1 - __n, 1);
      auto __sum = __term;
      for (int __k = 1; __k <= -__n; ++__k)
	{
	  __term *= __arg * __stirling_2_recur<_Tp>(1 - __n, 1 + __k);
	  __sum += __term;
	}
      __sum *= __gnu_cxx::__parity<_Tp>(1 - __n);
      return __sum;
    }

template<typename Tp>
  void
  test_polylog_nonpos_int(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\nTest negative integer order\n";
    for (auto n : {-0, -1, -2, -3, -4, -5})
      {
	std::cout << "\n\nn = " << n << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = -200; i <= 10; ++i)
	  {
	    auto x = del * i;
	    auto Ls_ceph = __polylog_nonpos_int(n, x);
	    auto Ls_gnu = __gnu_cxx::polylog(Tp(n), x);
	    std::cout << ' ' << n
		      << ' ' << x
		      << ' ' << Ls_ceph
		      << ' ' << Ls_gnu
		      << '\n';
	  }
      }
    std::cout << std::endl;
  }

template<typename Tp>
  void
  test_polylog_0(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against algebraic function\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Li0 = x / (Tp{1} - x);
	auto Lis_gnu = __gnu_cxx::polylog(Tp{0}, x);
	std::cout << ' ' << x
		  << ' ' << Li0
		  << ' ' << Lis_gnu
		  << '\n';
      }
    std::cout << std::endl;
  }

template<typename Tp>
  void
  test_polylog_1(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against algebraic function\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Li1 = -std::log(Tp{1} - x);
	auto Lis_gnu = __gnu_cxx::polylog(Tp{1}, x);
	std::cout << ' ' << x
		  << ' ' << Li1
		  << ' ' << Lis_gnu
		  << '\n';
      }
    std::cout << std::endl;
  }

template<typename Tp>
  void
  test_polylog_dilog(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against local dilog\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Ls_ceph = __gnu_cxx::dilog(x);
	auto Ls_gnu = __gnu_cxx::polylog(Tp(2), x);
	std::cout << ' ' << x
		  << ' ' << Ls_ceph
		  << ' ' << Ls_gnu
		  << '\n';
      }
    std::cout << std::endl;
  }

template<typename Tp>
  void
  test_polylog_cephes(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\nTest against Cephes for integer order\n";
    for (auto n : {0, 1, 2, 3, 4, 5})
      {
	std::cout << "\n\nn = " << n << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = -200; i <= 10; ++i)
	  {
	    auto x = del * i;
	    auto Ls_ceph = cephes::polylog(n, x);
	    auto Ls_gnu = __gnu_cxx::polylog(Tp(n), x);
	    std::cout << ' ' << n
		      << ' ' << x
		      << ' ' << Ls_ceph
		      << ' ' << Ls_gnu
		      << '\n';
	  }
      }
    std::cout << std::endl;
  }

template<typename Tp>
  void
  TestPolyLog(Tp proto = Tp{})
  {
    const auto _S_pi = __gnu_cxx::__const_pi(proto);
    const auto _S_2pi = __gnu_cxx::__const_2_pi(proto);

    std::cout.precision(__gnu_cxx::__digits10(proto) - 1);
    std::cout << std::scientific;

    std::cout << '\n';

    //std::size_t n = 5000;

    // this part of code was for performance testing.
    // the old implementation takes about 2.8s on my core2 and the new one 0.8s
    //     for(std::size_t i = 0; i < n; ++i)
    //     {
    //       Tp x = Tp{10}* static_cast<Tp>(i)/n + 1.1;
    // //      std::cout << std::scientific<<x<<' '<<
    //       std::__detail::__riemann_zeta(x)
    // //      std::tr1::__detail::__riemann_zeta(x)
    //       ;//between 1 and 10 riemann_zeta_glob is called
    // //      <<'\n';
    //     }

    // Something that didn't work in the original implementation
    std::cout << std::riemann_zeta(std::complex<Tp>{5.1, 0.5}) << '\n';
    std::cout << __gnu_cxx::hurwitz_zeta(Tp{5.1}, Tp{0.5}) << '\n';
    std::cout << __gnu_cxx::hurwitz_zeta(Tp{5.1}, std::complex<Tp>{0.5}) << '\n';
    std::cout << std::__detail::__hurwitz_zeta_polylog(Tp{5.1}, std::complex<Tp>{0.5}) << '\n';
    std::cout << std::__detail::__polylog_exp(Tp{2.5}, std::complex<Tp>(Tp{15}, Tp{1})) << '\n';

    for(std::size_t k = 0; k < 32; ++k)
    {
      std::cout << "=======  " << k << "  ==========" << '\n';
      auto w = std::complex<Tp>(Tp{0}, _S_2pi * k / Tp{32});
      std::cout << std::__detail::__polylog_exp(Tp{4}, w) << '\n';
      std::cout << std::__detail::__polylog_exp(-Tp{4}, w) << '\n';
      std::cout << std::__detail::__polylog_exp(Tp{2.6}, w) << '\n';
      std::cout << std::__detail::__polylog_exp(Tp{-2.6}, w) << '\n';
    }
    std::cout << std::endl;

    std::cout << std::__detail::__polylog_exp(Tp{2.6}, std::complex<Tp>(_S_pi, _S_pi)) << '\n';

    for(std::size_t k = 0; k < 10; ++k)
    {
      auto w = std::complex<Tp>(-_S_pi / 2 - _S_pi * 4 / 20, 0); // x == 0 ???
      std::cout << std::__detail::__polylog_exp(-Tp{4}, w) << '\n';
      std::cout << std::__detail::__polylog_exp_negative_real_part(-Tp{4}, w) << '\n';
    }
    std::cout << std::endl;

    std::cout << std::__detail::__polylog_exp_neg(Tp{-50.5}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    std::cout << std::__detail::__polylog_exp_neg(Tp{-5}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    std::cout << std::__detail::__polylog_exp_pos(Tp{2.3}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    //Don't trust Mathematica for small s
    std::cout << std::__detail::__polylog_exp_asymp(Tp{60.4}, std::complex<Tp>(Tp{30}, Tp{0})) << '\n';

    // auto l = 2;
    // auto p = std::atan(l);
    // Tp alpha[] = {0.5, 1, 1.5, 4};
    // std::ofstream data("el20.txt");
    // for(std::size_t a = 0; a < sizeof(alpha) / sizeof(Tp); ++a)
    // {
    //   for(int s = -1 ; s <= 1; s += 2)
    //   {
    //     for(auto k = -_S_pi; k < _S_pi; k += Tp{0.002})
    //       data << k << ' ' << std::sqrt(Tp{1} + l * l) * real(std::exp(std::complex<Tp>(0, -s * p)) / (std::exp(std::complex<Tp>(0, k)) - std::exp(-alpha[a]))) << '\n';
    //     data << "&" << '\n';
    //   }
    // }

    const auto del01 = Tp{1} / Tp{100};
    const auto del05 = Tp{1} / Tp{20};

    std::ofstream test("test.dat");
    for (auto s = Tp{2.5}; s < Tp{3.5}; s += del01)
      test << s << ' ' << std::real(std::__detail::__polylog(s, Tp{2})) - Tp{2} << '\n';
    std::cout << std::endl;

    std::cout << std::__detail::__polylog(Tp{3.1}, Tp{2}) << '\n';
    std::cout << std::__detail::__polylog_exp_pos(Tp{3.1}, std::complex<Tp>(std::log(Tp{2}))) << '\n';

    // Test function 1:
    for (std::size_t k = 3; k < 8; ++k)
      for (Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_pos(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << std::endl;

    // Test function 2
    for (std::size_t k = 3; k < 8; ++k)
      for (Tp x = 0; x < 6.28; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_pos(k, x)
		  << '\n';
    std::cout << std::endl;

    // Test function 3
    for (Tp k = -Tp{8}; k < 0; k += Tp{1}/Tp{13})
      for(Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_neg(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << std::endl;

    // Test function 4 + 5
    for (int k = -40; k < 0; ++k)
      for (Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_neg(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << std::endl;

    // Test series 6
    for (Tp k = Tp{1} / Tp{7}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = Tp{0}; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_pos(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << std::endl;

    // Test series 7
    for (Tp k = -Tp{13}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = Tp{0}; x < Tp{1}; x += del01)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_asymp(k, Tp{100} + std::polar(Tp{1}, _S_2pi * 0)) // x == 0 ???
		  << '\n';
    std::cout << std::endl;

    // Test series 8
    for (Tp k = -Tp{13}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = -Tp{7} / Tp{10} * _S_pi; x > -_S_2pi; x -= del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_negative_real_part(k, x)
		  << '\n';
    std::cout << std::endl;
  }

int
main()
{
  test_polylog_0(1.0);

  test_polylog_1(1.0);

  test_polylog_dilog(1.0);

  test_polylog_cephes(1.0);

  test_polylog_nonpos_int(1.0);

  TestPolyLog<double>();

  TestPolyLog<float>();

  TestPolyLog<long double>();

  TestPolyLog<__float128>();

  return 0;
}
