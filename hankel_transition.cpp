// $HOME/bin_specfun/bin/g++ -g -std=gnu++14 -o hankel_transition hankel_transition.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./hankel_transition > hankel_transition.txt

// g++ -std=gnu++14 -DNO_CBRT -o hankel_transition hankel_transition.cpp -lquadmath

// ./hankel_transition > hankel_transition.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <tuple>
//#include "float128.h"
#include "polynomial.h"
#include "rational.h"
#include "numeric_limits.h"

int
main()
{
  auto sign = [](int s, int r){return (s + r) % 2 == 1 ? -1 : +1; };
  using rational = double;

  __gnu_cxx::_Polynomial<int> poo({0, -2});
  std::vector<__gnu_cxx::_Polynomial<int>> phi;
  std::vector<__gnu_cxx::_Polynomial<int>> psi;
  phi.push_back(__gnu_cxx::_Polynomial<int>(1));
  phi.push_back(__gnu_cxx::_Polynomial<int>(0));
  phi.push_back(__gnu_cxx::_Polynomial<int>({0, -2}));
  psi.push_back(__gnu_cxx::_Polynomial<int>(0));
  psi.push_back(__gnu_cxx::_Polynomial<int>(-1));
  psi.push_back(__gnu_cxx::_Polynomial<int>(0));
  for (int s = 3; s <= 18; ++s)
    {
      phi.push_back(poo * phi[s - 2] -2 * (s - 2) * phi[s - 3]);
      psi.push_back(poo * psi[s - 2] -2 * (s - 2) * psi[s - 3]);
    }
  std::cout << "phi\n";
  for (const auto& p : phi)
    std::cout << p << '\n';
  std::cout << "psi\n";
  for (const auto& p : psi)
    std::cout << p << '\n';
  std::vector<__gnu_cxx::_Polynomial<rational>> A(6);
  std::vector<__gnu_cxx::_Polynomial<rational>> B(6);
  for (int s = 0; s < 6; ++s)
    {
      for (int r = 0; r <= s; ++r)
	A[s] += phi[2 * s + r] * __gnu_cxx::_Polynomial<rational>(sign(r, s), r);
      std::cout << "A_" << s << ": " << A[s] << '\n';
      for (int r = 0; r <= s; ++r)
	B[s] += psi[2 * s + r] * __gnu_cxx::_Polynomial<rational>(sign(r, s), r);
      std::cout << "B_" << s << ": " << B[s] << '\n';
    }

  std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::rational<int>>>
  P
  {
    {{1, 1}},
    {{}, {-1, 5}},
    {{}, {}, {3, 35}, {}, {}, {-9, 100}},
    {{-1, 225}, {}, {}, {-173, 3150}, {}, {}, {957, 7000}},
    {{}, {947, 346500}, {}, {}, {5903, 138600}, {}, {}, {-23573, 147000}, {}, {}, {27, 20000}}
  };

  std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::rational<int>>>
  Q
  {
    {{}, {}, {3, 10}},
    {{1, 70}, {}, {}, {-17, 70}},
    {{}, {-37, 3150}, {}, {}, {611, 3150}, {}, {}, {-9, 1000}},
    {{}, {}, {79, 12375}, {}, {}, {-110767, 693000}, {}, {}, {549, 28000}}
  };

  std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::rational<int>>>
  R
  {
    {{1, 1}},
    {{}, {-4, 5}},
    {{}, {}, {57, 70}, {}, {}, {-9, 100}},
    {{23, 3150}, {}, {}, {-2617, 3150}, {}, {}, {699, 3500}},
    {{}, {-1159, 115500}, {}, {}, {3889, 4620}, {}, {}, {-46631, 147000}, {}, {}, {27, 20000}}
  };

  std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::rational<int>>>
  S
  {
    {{-1, 5}, {}, {}, {3, 5}},
    {{}, {1, 5}, {}, {}, {-131, 140}},
    {{}, {}, {-593, 3150}, {}, {}, {5437, 4500}, {}, {}, {-9, 500}},
    {{947, 346500}, {}, {}, {31727, 173250}, {}, {}, {-999443, 693000}, {}, {}, {369, 7000}}
  };

  for (int i = -100; i <= +100; ++i)
  {
    auto z = 0.01L * i;
    std::cout << setw(width) << z
	      << setw(width) << P(z)
	      << setw(width) << Q(z)
	      << setw(width) << R(z)
	      << setw(width) << S(z)
	      << '\n';
  }
}
