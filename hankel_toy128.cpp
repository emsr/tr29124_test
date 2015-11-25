// $HOME/bin/bin/g++ -std=gnu++14 -o hankel_toy128 hankel_toy128.cpp -L$HOME/bin/lib64 -lquadmath

#include <limits>
#include <iostream>
#include <iomanip>
#include "float128.h"
#include "polynomial"

int
main()
{
  auto index = 0;
  auto indexp = 0;
  for (int k = 0; k <= 20; ++k)
    {
      std::cout << '\n';
      std::cout << "k      = " << k << '\n';
      std::cout << "index  = " << index << '\n';
      std::cout << "indexp = " << indexp << '\n';
      auto indexpold = indexp;
      index += 2;
      indexp += 2;
      for (; index < indexpold; ++index, ++indexp)
        /*std::cout << "  index  = " << index << '\n'*/ ;
      ++indexp;
      index = indexp;
      indexp += 2 * k + 3;
    }

  for (int k = 0; k <= 20; ++k)
    {
      std::cout << '\n';
      std::cout << "k      = " << k << '\n';
      std::cout << "index  = " << k * (2 * k + 1) << '\n';
      std::cout << "indexp = " << (k + 1) * (2 * k + 1) << '\n';
    }

  auto prec = std::numeric_limits<__float128>::max_digits10;
  auto width = prec + 6;

  std::cout.precision(prec);
  std::cout << std::scientific;

  std::cout << '\n';
  __float128 lambda = 1.0Q;
  __float128 mu = -1.0Q;
  __float128 denom = 1.0Q;
  for (int s = 1; s <= 50; ++s)
    {
      std::cout << std::setw(width) << lambda << '\t' << mu << '\n';
      denom *= s * 144;
      __float128 numer = 1.0Q;
      for (int m = 2 * s + 1; m <= 6 * s - 1; m += 2)
        numer *= m;
      lambda = numer / denom;
      mu = -(6 * s + 1) * lambda / (6 * s - 1);
    }

  std::polynomial<__float128> upol1{0.0Q, 0.0Q, 0.5Q, 0.0Q, -0.5Q};
  std::polynomial<__float128> upol2{+0.125Q, 0.0Q, -0.625Q};
  std::polynomial<__float128> vpol1{0.0Q, -0.5Q, 0.0Q, +0.5Q};
  std::polynomial<__float128> vpol2{0.0Q, 0.0Q, -1.0Q, 0.0Q, +1.0Q};
  std::polynomial<__float128> u{1.0Q};
  std::vector<std::polynomial<__float128>> uvec;
  std::polynomial<__float128> v{1.0Q};
  std::vector<std::polynomial<__float128>> vvec;
  for (auto k = 1; k <= 20; ++k)
    {
      uvec.push_back(u);
      vvec.push_back(v);
      v = vpol1 * u + vpol2 * u.derivative();
      u = upol1 * u.derivative() + (upol2 * u).integral(0.0);
      v += u;
    }
  std::cout << '\n';
  for (const auto & u : uvec)
    std::cout << u << '\n';
  std::cout << '\n';
  for (const auto & v : vvec)
    std::cout << v << '\n';

  std::vector<std::vector<std::tuple<int, int, __float128>>> uentry;
  auto ku = 0;
  for (const auto & u : uvec)
    {
      uentry.resize(u.degree() + 1);
      for (int i = 0; i <= u.degree(); ++i)
	if (u.coefficient(i) != 0)
	  uentry[i].push_back(std::make_tuple(ku, i, u.coefficient(i)));
      ++ku;
    }
  std::cout << '\n';
  auto iu = 0;
  for (const auto & u : uentry)
    {
      for (const auto & c : u)
	std::cout << ' ' << std::setw(3) << ++iu
		  << ' ' << std::setw(3) << std::get<0>(c)
		  << ' ' << std::setw(3) << std::get<1>(c)
		  << ' ' << std::setw(width) << std::get<2>(c) << '\n';
    }
  std::vector<std::vector<std::tuple<int, int, __float128>>> ventry;
  auto kv = 0;
  for (const auto & v : vvec)
    {
      ventry.resize(v.degree() + 1);
      for (int i = 0; i <= v.degree(); ++i)
	if (v.coefficient(i) != 0)
	  ventry[i].push_back(std::make_tuple(kv, i, v.coefficient(i)));
      ++kv;
    }
  std::cout << '\n';
  auto iv = 0;
  for (const auto & v : ventry)
    {
      for (const auto & c : v)
	std::cout << ' ' << std::setw(3) << ++iv
		  << ' ' << std::setw(3) << std::get<0>(c)
		  << ' ' << std::setw(3) << std::get<1>(c)
		  << ' ' << std::setw(width) << std::get<2>(c) << '\n';
    }

  std::cout << '\n';
  for (const auto& u : uvec)
    for (auto c = u.crbegin(); c != u.crend(); ++c)
      if (*c != 0)
	std::cout << *c << '\n';

  std::cout << '\n';
  for (const auto& v : vvec)
    for (auto c = v.crbegin(); c != v.crend(); ++c)
      if (*c != 0)
	std::cout << *c << '\n';
}
