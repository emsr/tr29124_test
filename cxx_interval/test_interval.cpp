
#include <limits>
#include <iostream>
#include <iomanip>

#include <ext/interval>

template<typename _Tp>
  void
  test_rounded_cmath()
  {
    std::cout.precision(std::numeric_limits<_Tp>::max_digits10);
    const auto w = 6 + std::cout.precision();

    const auto a = _Tp{5};
    interval<_Tp> ivl(a);
    std::cout << ivl << '\n';

    std::cout << std::setw(w) << "sqrt(a)\n";
    std::cout << std::setw(w) << sqrt(a, fpround::upward) << '\n';
    std::cout << std::setw(w) << sqrt(a, fpround::to_nearest) << '\n';
    std::cout << std::setw(w) << sqrt(a, fpround::downward) << '\n';
    std::cout << std::setw(w) << sqrt(a, fpround::toward_zero) << '\n';
    std::cout << std::setw(w) << sqrt(ivl) << '\n';

    std::cout << std::setw(w) << "sin(a)\n";
    std::cout << std::setw(w) << sin(a, fpround::upward) << '\n';
    std::cout << std::setw(w) << sin(a, fpround::to_nearest) << '\n';
    std::cout << std::setw(w) << sin(a, fpround::downward) << '\n';
    std::cout << std::setw(w) << sin(a, fpround::toward_zero) << '\n';
    std::cout << std::setw(w) << sin(ivl) << '\n';

    std::cout << std::setw(w) << "cos(a)\n";
    std::cout << std::setw(w) << cos(a, fpround::upward) << '\n';
    std::cout << std::setw(w) << cos(a, fpround::to_nearest) << '\n';
    std::cout << std::setw(w) << cos(a, fpround::downward) << '\n';
    std::cout << std::setw(w) << cos(a, fpround::toward_zero) << '\n';
    std::cout << std::setw(w) << cos(ivl) << '\n';

    std::cout << std::setw(w) << "tan(a)\n";
    std::cout << std::setw(w) << tan(a, fpround::upward) << '\n';
    std::cout << std::setw(w) << tan(a, fpround::to_nearest) << '\n';
    std::cout << std::setw(w) << tan(a, fpround::downward) << '\n';
    std::cout << std::setw(w) << tan(a, fpround::toward_zero) << '\n';
    std::cout << std::setw(w) << tan(ivl) << '\n';
  }

int
main()
{
  test_rounded_cmath<float>();
  test_rounded_cmath<double>();
  test_rounded_cmath<long double>();

  return 0;
}
