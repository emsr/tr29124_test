// $HOME/bin_specfun/bin/g++ -g -std=gnu++14 -o hankel_transition hankel_transition.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./hankel_transition > hankel_transition.txt

// g++ -std=gnu++14 -DNO_CBRT -o hankel_transition hankel_transition.cpp -lquadmath

// ./hankel_transition > hankel_transition.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <tuple>
//#include "float128.h"
#include "polynomial"
#include "numeric_limits.h"

int
main()
{
  auto sign = [](int s, int r){return (s + r) % 2 == 1 ? -1 : +1; };
  using rational = double;

  std::polynomial<int> poo({0, -2});
  std::vector<std::polynomial<int>> phi;
  std::vector<std::polynomial<int>> psi;
  phi.push_back(std::polynomial<int>(1));
  phi.push_back(std::polynomial<int>(0));
  phi.push_back(std::polynomial<int>({0, -2}));
  psi.push_back(std::polynomial<int>(0));
  psi.push_back(std::polynomial<int>(-1));
  psi.push_back(std::polynomial<int>(0));
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
  std::vector<std::polynomial<rational>> A(6);
  std::vector<std::polynomial<rational>> B(6);
  for (int s = 0; s < 6; ++s)
    {
      for (int r = 0; r <= s; ++r)
	A[s] += phi[2 * s + r] * std::polynomial<rational>(sign(r, s), r);
      std::cout << "A_" << s << ": " << A[s] << '\n';
      for (int r = 0; r <= s; ++r)
	B[s] += psi[2 * s + r] * std::polynomial<rational>(sign(r, s), r);
      std::cout << "B_" << s << ": " << B[s] << '\n';
    }
}
