// $HOME/bin/bin/g++ -o hankhelp hankhelp.cpp

#include <limits>
#include <iostream>
#include <iomanip>

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

  auto width = std::numeric_limits<long double>::digits10;

  std::cout.precision(width);
  std::cout << std::scientific;

  std::cout << '\n';
  long double lambda = 1.0L;
  long double mu = -1.0L;
  long double denom = 1.0L;
  for (int s = 1; s <= 50; ++s)
    {
      std::cout << std::setw(width + 2) << lambda << '\t' << mu << '\n';
      denom *= s * 144;
      long double numer = 1.0L;
      for (int m = 2 * s + 1; m <= 6 * s - 1; m += 2)
        numer *= m;
      lambda = numer / denom;
      mu = -(6 * s + 1) * lambda / (6 * s - 1);
    }

  std::cout << '\n';
  long double u = 1.0L;
  for (int k = 1; k <= 10; ++k)
    {
      
    }
}
