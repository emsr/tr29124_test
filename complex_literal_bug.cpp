// g++ -std=c++14 -o complex_literal_bug complex_literal_bug.cpp

// $HOME/bin/bin/g++ -std=c++14 -o complex_literal_bug complex_literal_bug.cpp

#include <complex>

std::complex<double>
func(std::complex<double> z)
{
  return 2.0 * z;
}

int
main()
{
  using namespace std::literals::complex_literals;
  auto w = func(1.0 + 3.0i);
  w /= 2.0i;
  auto k = w / 4.0i;
  auto q = w * 1.5i;
}
