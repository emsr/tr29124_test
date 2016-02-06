// $HOME/bin_specfun/bin/g++ -std=c++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_polylog test_polylog.cpp

// ./test_polylog > test_polylog.txt

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>

//#include "sf_polylog.tcc"

template<typename Tp>
  void
  TestPolyLog()
  {
    uint n = 5000;
    // this part of code was for performance testing.
    // the old implementation takes about 2.8s on my core2 and the new one 0.8s
    std::cout.precision(14);
    //     for(uint i = 0; i < n; ++i)
    //     {
    //       Tp x = Tp{10}* static_cast<Tp>(i)/n + 1.1;
    // //      std::cout << std::scientific<<x<<" "<<
    //       std::__detail::__riemann_zeta(x)
    // //      std::tr1::__detail::__riemann_zeta(x)
    //       ;//between 1 and 10 riemann_zeta_glob is called
    // //      <<std::endl;
    //     }

    // something that didn't work in the original implementation
    std::cout << __gnu_cxx::hurwitz_zeta(Tp{5.1}, Tp{0.5}) << std::endl;
    std::cout << __gnu_cxx::hurwitz_zeta(Tp{5.1}, std::complex<Tp>{0.5}) << std::endl;
    std::cout << std::__detail::__polylog_exp(Tp{2.5}, std::complex<Tp>(Tp{15}, Tp{1})) << std::endl;
    for(uint k = 0; k < 32; ++k)
    {
      std::cout << "=======  "<<k<<"  ==========" << std::endl;
      std::cout << std::__detail::__polylog_exp(Tp{4}, std::complex<Tp>(Tp{0}, Tp{2}*M_PI * k/Tp{32})) << std::endl;
      std::cout << std::scientific << std::__detail::__polylog_exp(-Tp{4}, std::complex<Tp>(Tp{0}, Tp{2}*M_PI * k/Tp{32})) << std::endl;
      std::cout << std::__detail::__polylog_exp(Tp{2.6}, std::complex<Tp>(Tp{0}, Tp{2}*M_PI * k/Tp{32})) << std::endl;
      std::cout << std::__detail::__polylog_exp(Tp{-2.6}, std::complex<Tp>(Tp{0}, Tp{2}*M_PI * k/Tp{32})) << std::endl;
    }
    std::cout << std::__detail::__polylog_exp(Tp{2.6}, std::complex<Tp>(M_PI, M_PI)) << std::endl;
    for(uint k = 0; k < 10; ++k)
    {
      std::cout << std::__detail::__polylog_exp(-Tp{4}, std::complex<Tp>(-M_PI/2 - M_PI*4/20, 0)) << std::endl;
      std::cout << std::__detail::__polylog_exp_negative_real_part(-Tp{4}, std::complex<Tp>(-M_PI/2 - M_PI*4/20, 0)) << std::endl;
    }
    std::cout << std::__detail::__polylog_exp_neg(Tp{-50.5}, std::complex<Tp>(Tp{1}, Tp{1})) << std::endl;
    std::cout << std::__detail::__polylog_exp_neg(Tp{-5}, std::complex<Tp>(Tp{1}, Tp{1})) << std::endl;
    std::cout << std::__detail::__polylog_exp_pos(Tp{2.3}, std::complex<Tp>(Tp{1}, Tp{1})) << std::endl;
    std::cout << std::__detail::__polylog_exp_asym(Tp{60.4}, std::complex<Tp>(Tp{30}, Tp{0})) << std::endl;//Don't trust Mathematica for small s
    // auto l = 2;
    // auto p = std::atan(l);
    // Tp alpha[] = {0.5, 1, 1.5, 4};
    // std::ofstream data("el20.txt");
    // for(uint a = 0; a < sizeof(alpha)/sizeof(Tp); ++a)
    // {
    //   for(int s = -1 ; s <= 1; s += 2)
    //   {
    // for(auto k = -M_PI; k < M_PI; k+= Tp{0.002})
    //   data<<k<<" "<<std::sqrt(Tp{1} + l*l)*real(std::exp(std::complex<Tp>(0, -s*p)) / (std::exp(std::complex<Tp>(0, k)) - std::exp(-alpha[a]))) << std::endl;
    // data<<"&"<<std::endl;
    //   }
    // }

    std::ofstream test("test.dat");
    for (double s = Tp{2.5}; s < Tp{3.5}; s += Tp{0.01})
      test << s << " " << std::real(std::__detail::__polylog(s, Tp{2})) - Tp{2} << std::endl;
    std::cout << std::__detail::__polylog(Tp{3.1}, Tp{2}) << std::endl;
    std::cout << std::__detail::__polylog_exp_pos(Tp{3.1}, std::complex<Tp>(std::log(Tp{2}))) << std::endl;

    //test function 1:
    for (uint k = 3; k < 8; ++k)
    {
      for (Tp x = 0; x < Tp{1}; x += Tp{0.05})
      {
	std::__detail::__polylog_exp_pos(k, std::polar(Tp{1}, Tp{2}*M_PI* x));
      }
    }

    //test function 2
    for (uint k = 3; k < 8; ++k)
    {
      for (Tp x = 0; x < 6.28; x += Tp{0.05})
      {
	std::__detail::__polylog_exp_pos(k, x);
      }
    }


    //test function 3
    for (Tp k = -Tp{8}; k < 0; k += Tp{1}/Tp{13})
    {
      for(Tp x = 0; x < Tp{1}; x += Tp{0.05})
      {
	std::__detail::__polylog_exp_neg(k, std::polar(Tp{1}, Tp{2}*M_PI* x));
      }
    }

    //test function 4 + 5
    for (int k = -40; k < 0; ++k)
    {
      for (Tp x = 0; x < Tp{1}; x += Tp{0.05})
      {
	std::__detail::__polylog_exp_neg(k, std::polar(Tp{1}, Tp{2}*M_PI* x));
      }
    }

    //test series 6
    for (Tp k = Tp{1}/Tp{7}; k < Tp{13}; k += Tp{1}/Tp{11})
    {
      for (Tp x = Tp{0}; x < Tp{1}; x += Tp{0.05})
      {
	std::__detail::__polylog_exp_pos(k, std::polar(Tp{1}, Tp{2}*M_PI* x));
      }
    }
    //test series 7
    for (Tp k = -Tp{13}; k < Tp{13}; k += Tp{1}/Tp{11})
    {
      for (Tp x = Tp{0}; x < Tp{1}; x += Tp{0.01})
      {
	std::__detail::__polylog_exp_asym(k, Tp{100} + std::polar(Tp{1}, Tp{2}*M_PI* 0) );
      }
    }

    //test series 8
    for (Tp k = -Tp{13}; k < Tp{13}; k += Tp{1}/Tp{11})
    {
      for (Tp x = -Tp{7}/Tp{10}*M_PI; x > -Tp{2}*M_PI; x -= Tp{0.05})
      {
	std::cout << k << " " << x << " " << std::__detail::__polylog_exp_negative_real_part(k, x) << std::endl;
      }
    }
  }

int
main()
{
  TestPolyLog<double>();
  return 0;
}
