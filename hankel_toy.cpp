// $HOME/bin/bin/g++ -std=gnu++14 -o hankel_toy hankel_toy.cpp

#include <limits>
#include <iostream>
#include <iomanip>
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

  auto prec = std::numeric_limits<long double>::max_digits10;
  auto width = prec + 6;

  std::cout.precision(prec);
  std::cout << std::scientific;

  std::cout << '\n';
  long double lambda = 1.0L;
  long double mu = -1.0L;
  long double denom = 1.0L;
  for (int s = 1; s <= 50; ++s)
    {
      std::cout << std::setw(width) << lambda << '\t' << mu << '\n';
      denom *= s * 144;
      long double numer = 1.0L;
      for (int m = 2 * s + 1; m <= 6 * s - 1; m += 2)
        numer *= m;
      lambda = numer / denom;
      mu = -(6 * s + 1) * lambda / (6 * s - 1);
    }

  std::vector<long double> ustore;
  std::vector<std::vector<std::pair<int, long double>>> uentry;
  std::polynomial<long double> upol1{0.0L, 0.0L, 0.5L, 0.0L, -0.5L};
  std::polynomial<long double> upol2{0.125L, 0.0L, -0.625L};
  std::polynomial<long double> u{1.0L};
  std::cout << '\n';
  for (auto k = 1; k <= 20; ++k)
    {
      std::cout << u << '\n';
      uentry.resize(u.degree() + 1);
      for (int i = 0; i <= u.degree(); ++i)
	if (u.coefficient(i) != 0)
	  {
	    ustore.push_back(u.coefficient(i));
	    uentry[i].push_back(std::make_pair(k - 1, u.coefficient(i)));
	  }
      u = upol1 * u.derivative() + (upol2 * u).integral(0.0);
    }

  std::cout << '\n';
  auto iu = 0;
  for (const auto & c : ustore)
    std::cout << ' ' << std::setw(3) << ++iu
	      << ' ' << std::setw(width) << c << '\n';

  std::cout << '\n';
  iu = 0;
  for (const auto & p : uentry)
    for (const auto & c : p)
      std::cout << ' ' << std::setw(3) << ++iu
		<< ' ' << std::setw(3) << c.first
		<< ' ' << std::setw(width) << c.second << '\n';

  std::vector<std::vector<std::pair<int, long double>>> ventry;
  std::polynomial<long double> v{1.0L};
  std::cout << '\n';
  for (int k = 1; k <= 20; ++k)
    {
      std::cout << v << '\n';
      
    }
  std::cout << '\n';
  auto iv = 0;
  for (const auto & p : ventry)
    for (const auto & c : p)
      std::cout << ' ' << std::setw(3) << ++iv
		<< ' ' << std::setw(3) << c.first
		<< "  " << std::setw(width) << c.second << '\n';
}
