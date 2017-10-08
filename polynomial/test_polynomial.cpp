/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -I.. -o test_polynomial test_polynomial.cpp -lquadmath
./test_polynomial > test_polynomial.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_polynomial test_polynomial.cpp -lquadmath
./test_polynomial > test_polynomial.txt
*/

#include <cmath>
#include <iostream>
#include <complex>
#include <sstream>

#include "polynomial.h"

int
main()
{
  using namespace std::literals::complex_literals;

  std::cout.setf(std::ios_base::boolalpha);

  __gnu_cxx::_Polynomial<double> P({0.0, 1.0, 2.0, 3.0});
  std::cout << "P = " << P << '\n';
  std::cout << "+P = " << +P << '\n';
  std::cout << "-P = " << -P << '\n';
  std::cout << "P = " << P << '\n';
  std::cout << "degree(P) = " << P.degree() << '\n';
  __gnu_cxx::_Polynomial<double> Q({2.0, 1.0});
  std::cout << "Q = " << Q << '\n';
  std::cout << "degree(Q) = " << Q.degree() << '\n';
  std::cout << "P + Q = " << P + Q << '\n';
  std::cout << "P - Q = " << P - Q << '\n';
  std::cout << "P * Q = " << P * Q << '\n';
  std::cout << "P / Q = " << P / Q << '\n';
  std::cout << "P % Q = " << P % Q << '\n';
  double b = 5.0;
  std::cout << "b = " << b << '\n';
  std::cout << "P + b = " << P + b << '\n';
  std::cout << "P - b = " << P - b << '\n';
  std::cout << "P * b = " << P * b << '\n';
  std::cout << "P / b = " << P / b << '\n';
  std::cout << "P % b = " << P % b << '\n';
  double a = 2.0;
  std::cout << "a = " << a << '\n';
  std::cout << "a + Q = " << a + Q << '\n';
  std::cout << "a - Q = " << a - Q << '\n';
  std::cout << "a * Q = " << a * Q << '\n';
  std::cout << "a / Q = " << a / Q << '\n';
  std::cout << "a % Q = " << a % Q << '\n';

  __gnu_cxx::_Polynomial<double> B;// = b;
  B = b;
  std::cout << "B = " << B << '\n';
  std::cout << "P % B = " << P % B << '\n';

  Q = {0.0, -2.0, 4.0, -6.0, 8.0, -12.0};
  std::cout << "Q = " << Q << '\n';

  __gnu_cxx::_Polynomial<double> P2;
  P2 = P;
  std::cout << "P2 = " << P2 << '\n';
  std::cout << "P2 == P = " << (P2 == P) << '\n';

  for (int i = 0; i <= 100; ++i)
  {
    double x = i * 0.1;
    std::cout << "P(" << x << ") = " << P(x) << '\n';
  }

  __gnu_cxx::_Polynomial<std::complex<double>>
  CP({std::complex<double>(0.0, -1.0),
      std::complex<double>(1.0, -2.0), 
      std::complex<double>(2.0, -3.0), 
      std::complex<double>(3.0, -4.0)});
  std::cout << "CP = " << CP << '\n';
  std::cout << "CP * CP = " << CP * CP << '\n';

  __gnu_cxx::_Polynomial<int> IP({0, 1, 2, 3});
  std::cout << "IP = " << IP << '\n';
  std::cout << "IP * IP = " << IP * IP << '\n';

  std::array<double, 10> arr;
  P.eval(1.0, arr);
  std::cout << "P(" << 1.0 << ") =";
  for (unsigned i = 0; i < arr.size(); ++i)
    std::cout << " " << arr[i];
  std::cout << '\n';

  P.eval(1.0, arr.begin(), arr.end());
  std::cout << "P(" << 1.0 << ") =";
  for (auto iarr = arr.cbegin(); iarr != arr.cend(); ++iarr)
    std::cout << " " << *iarr;
  std::cout << '\n';

  std::istringstream is("(-2.0, -1.0, 0.0)");
  __gnu_cxx::_Polynomial<double> R;
  is >> R;
  std::cout << "R = " << R << '\n';

  std::istringstream is2("(5.0)");
  __gnu_cxx::_Polynomial<double> S;
  is2 >> S;
  std::cout << "S = " << S << '\n';

  std::istringstream is3("42.0");
  __gnu_cxx::_Polynomial<double> T;
  is3 >> T;
  std::cout << "T = " << T << '\n';

  __gnu_cxx::_Polynomial<double> u({1.0, 3.0, 3.0, 1.0});
  __gnu_cxx::_Polynomial<double> v({1.0, 1.0});
  __gnu_cxx::_Polynomial<double> q, r;
  __gnu_cxx::divmod(u, v, q, r);
  std::cout << "u = " << u << '\n';
  std::cout << "v = " << v << '\n';
  std::cout << "q = " << q << '\n';
  std::cout << "r = " << r << '\n';

  __gnu_cxx::_Polynomial<double> u1({1.0, -3.0, 3.0, -1.0});
  __gnu_cxx::_Polynomial<double> v1({1.0, -2.0, 1.0});
  __gnu_cxx::_Polynomial<double> q1, r1;
  __gnu_cxx::divmod(u1, v1, q1, r1);
  std::cout << "u1 = " << u1 << '\n';
  std::cout << "v1 = " << v1 << '\n';
  std::cout << "q1 = " << q1 << '\n';
  std::cout << "r1 = " << r1 << '\n';

  __gnu_cxx::_Polynomial<double> u2({1.0, 1.0});
  __gnu_cxx::_Polynomial<double> v2({1.0, 3.0, 3.0, 1.0});
  __gnu_cxx::_Polynomial<double> q2, r2;
  __gnu_cxx::divmod(u2, v2, q2, r2);
  std::cout << "u2 = " << u2 << '\n';
  std::cout << "v2 = " << v2 << '\n';
  std::cout << "q2 = " << q2 << '\n';
  std::cout << "r2 = " << r2 << '\n';

  __gnu_cxx::_Polynomial<double> u3({1.0, 0.0, 0.0, 1.0});
  __gnu_cxx::_Polynomial<double> v3({1.0, 3.0, 3.0, 1.0});
  __gnu_cxx::_Polynomial<double> q3, r3;
  __gnu_cxx::divmod(u3, v3, q3, r3);
  std::cout << "u3 = " << u3 << '\n';
  std::cout << "v3 = " << v3 << '\n';
  std::cout << "q3 = " << q3 << '\n';
  std::cout << "r3 = " << r3 << '\n';

  std::cout << "P = " << P << '\n';
  std::cout << "P' = " << P.derivative() << '\n';
  std::cout << "I = " << P.integral(42.0) << '\n';
  
  __gnu_cxx::_Polynomial<__gnu_cxx::_Polynomial<double>>
  Pp({{9.0, 8.0, 7.0}, {6.0, 5.0, 4.0}, {3.0, 2.0, 1.0}});
  std::cout << "Pp = " << Pp << '\n';
  std::cout << "Pp(-1) = " << Pp(-1.0) << '\n';
  std::cout << "Pp(2) = " << Pp(2.0) << '\n';

  //std::array<double, 5> aaa{{1.1, 2.2, 3.3, 4.4, 5.5}};
  //std::cout << "aaa = " << __gnu_cxx::_Polynomial_eval(aaa, 3.13149) << '\n';

  std::cout << "P = " << P << '\n';
  std::cout << "P(1) = " << P(1.0) << '\n';
  std::cout << "P(i) = " << P(1.0i) << '\n'; // Fucked...
  std::cout << "P(i) = " << P(std::complex<double>{0, 1}) << '\n';

  __gnu_cxx::_Polynomial<double> e([](unsigned int k) -> double { return 1.0 / __gnu_cxx::factorial<double>(k); }, 20);
  std::cout << "e = " << e << '\n';
  std::cout << "e(1) = " << e(1) << '\n';
}

