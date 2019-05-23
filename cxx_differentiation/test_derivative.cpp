/*
g++ -std=c++17 -g -Iinclude -o test_derivative test_derivative.cpp
*/

#include <ext/derivative.h>

#include <cmath>
#include <string_view>
#include <iostream>
#include <iomanip>

template<typename Tp>
  Tp
  f1(Tp x)
  { return std::exp(x); }

template<typename Tp>
  Tp
  df1(Tp x)
  { return std::exp(x); }

template<typename Tp>
  Tp
  f2(Tp x)
  {
    if (x >= Tp{0})
      return x * std::sqrt(x);
    else
      return Tp{0};
  }

template<typename Tp>
  Tp
  df2(Tp x)
  {
    if (x >= Tp{0})
      return Tp{1.5L} * std::sqrt(x);
    else
      return Tp{0};
  }

template<typename Tp>
  Tp
  f3(Tp x)
  {
    if (x != Tp{0})
      return std::sin(Tp{1} / x);
    else
      return Tp{0};
  }

template<typename Tp>
  Tp
  df3(Tp x)
  {
    if (x != Tp{0})
      return -std::cos(Tp{1} / x) / (x * x);
    else
      return Tp{0};
  }

template<typename Tp>
  Tp
  f4(Tp x)
  { return std::exp(-x * x); }

template<typename Tp>
  Tp
  df4(Tp x)
  { return Tp{-2} * x * std::exp(-x * x); }

template<typename Tp>
  Tp
  f5(Tp x)
  { return x * x; }

template<typename Tp>
  Tp
  df5(Tp x)
  { return Tp{2} * x; }

template<typename Tp>
  Tp
  f6(Tp x)
  { return Tp{1} / x; }

template<typename Tp>
  Tp
  df6(Tp x)
  { return Tp{-1} / (x * x); }

template<typename Differ, typename Func, typename Deriv, typename Tp>
  void
  test_derivative(Differ differ, Func func, Deriv deriv, Tp x,
		  std::string_view desc)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 6 + std::cout.precision();

    auto [res, err_trunc, err_round] = differ(func, x, Tp{1e-4});
    auto err_abs = err_trunc + err_round;
    auto exp = deriv(x);

    if (std::abs(res - exp) > err_abs)
      std::cout << "FAIL: ";
    else
      std::cout << "PASS: ";
    std::cout << " (" << std::setw(w) << res << " observed vs "
	      << std::setw(w) << exp << " expected)";
    std::cout << " (" << std::setw(w) << err_trunc << " trunc and "
	      << std::setw(w) << err_round << " round)";
    std::cout << ' ' << desc << '\n';
  }

template<typename Tp>
  void
  test_all()
  {
    using Func = Tp(Tp);
    test_derivative(derivative_central<Func,Tp>, f1<Tp>, df1<Tp>, Tp{1}, "exp(x), x=1, central deriv");
    test_derivative(derivative_forward<Func,Tp>, f1<Tp>, df1<Tp>, Tp{1}, "exp(x), x=1, forward deriv");
    test_derivative(derivative_backward<Func,Tp>, f1<Tp>, df1<Tp>, Tp{1}, "exp(x), x=1, backward deriv");
    test_derivative(derivative_ridder<Func,Tp>, f1<Tp>, df1<Tp>, Tp{1}, "exp(x), x=1, Ridder's deriv");

    test_derivative(derivative_central<Func,Tp>, f2<Tp>, df2<Tp>, Tp{0.1L}, "x^(3/2), x=0.1, central deriv");
    test_derivative(derivative_forward<Func,Tp>, f2<Tp>, df2<Tp>, Tp{0.1L}, "x^(3/2), x=0.1, forward deriv");
    test_derivative(derivative_backward<Func,Tp>, f2<Tp>, df2<Tp>, Tp{0.1L}, "x^(3/2), x=0.1, backward deriv");
    test_derivative(derivative_ridder<Func,Tp>, f2<Tp>, df2<Tp>, Tp{0.1L}, "x^(3/2), x=0.1, Ridder's deriv");

    test_derivative(derivative_central<Func,Tp>, f3<Tp>, df3<Tp>, Tp{0.45L}, "sin(1/x), x=0.45, central deriv");
    test_derivative(derivative_forward<Func,Tp>, f3<Tp>, df3<Tp>, Tp{0.45L}, "sin(1/x), x=0.45, forward deriv");
    test_derivative(derivative_backward<Func,Tp>, f3<Tp>, df3<Tp>, Tp{0.45L}, "sin(1/x), x=0.45, backward deriv");
    test_derivative(derivative_ridder<Func,Tp>, f3<Tp>, df3<Tp>, Tp{0.45L}, "sin(1/x), x=0.45, Ridder's deriv");

    test_derivative(derivative_central<Func,Tp>, f4<Tp>, df4<Tp>, Tp{0.5L}, "exp(-x^2), x=0.5, central deriv");
    test_derivative(derivative_forward<Func,Tp>, f4<Tp>, df4<Tp>, Tp{0.5L}, "exp(-x^2), x=0.5, forward deriv");
    test_derivative(derivative_backward<Func,Tp>, f4<Tp>, df4<Tp>, Tp{0.5L}, "exp(-x^2), x=0.5, backward deriv");
    test_derivative(derivative_ridder<Func,Tp>, f4<Tp>, df4<Tp>, Tp{0.5L}, "exp(-x^2), x=0.5, Ridder's deriv");

    test_derivative(derivative_central<Func,Tp>, f5<Tp>, df5<Tp>, Tp{0}, "x^2, x=0, central deriv");
    test_derivative(derivative_forward<Func,Tp>, f5<Tp>, df5<Tp>, Tp{0}, "x^2, x=0, forward deriv");
    test_derivative(derivative_backward<Func,Tp>, f5<Tp>, df5<Tp>, Tp{0}, "x^2, x=0, backward deriv");
    test_derivative(derivative_ridder<Func,Tp>, f5<Tp>, df5<Tp>, Tp{0}, "x^2, x=0, Ridder's deriv");

    test_derivative(derivative_central<Func,Tp>, f6<Tp>, df6<Tp>, Tp{10}, "1/x, x=10, central deriv");
    test_derivative(derivative_forward<Func,Tp>, f6<Tp>, df6<Tp>, Tp{10}, "1/x, x=10, forward deriv");
    test_derivative(derivative_backward<Func,Tp>, f6<Tp>, df6<Tp>, Tp{10}, "1/x, x=10, backward deriv");
    test_derivative(derivative_ridder<Func,Tp>, f6<Tp>, df6<Tp>, Tp{10}, "1/x, x=10, Ridder's deriv");
  }

int
main()
{
  std::cout << "\nfloat\n";
  test_all<float>();

  std::cout << "\ndouble\n";
  test_all<double>();

  std::cout << "\nlong double\n";
  test_all<long double>();
}
