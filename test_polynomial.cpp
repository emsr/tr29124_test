// $HOME/bin/bin/g++ -o test_polynomial test_polynomial.cpp

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_polynomial

//  Get past a bug....
// $HOME/bin/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__=0 -o test_polynomial test_polynomial.cpp

#include "polynomial"
#include <iostream>
#include <complex>
#include <sstream>

int
main()
{
  std::cout.setf(std::ios_base::boolalpha);

  std::polynomial<double> P({0.0, 1.0, 2.0, 3.0});
  std::cout << "P = " << P << std::endl;
  std::cout << "+P = " << +P << std::endl;
  std::cout << "-P = " << -P << std::endl;
  std::cout << "P = " << P << std::endl;
  std::cout << "degree(P) = " << P.degree() << std::endl;
  std::polynomial<double> Q({2.0, 1.0});
  std::cout << "Q = " << Q << std::endl;
  std::cout << "degree(Q) = " << Q.degree() << std::endl;
  std::cout << "P + Q = " << P + Q << std::endl;
  std::cout << "P - Q = " << P - Q << std::endl;
  std::cout << "P * Q = " << P * Q << std::endl;
  std::cout << "P / Q = " << P / Q << std::endl;
  std::cout << "P % Q = " << P % Q << std::endl;
  double b = 5.0;
  std::cout << "b = " << b << std::endl;
  std::cout << "P + b = " << P + b << std::endl;
  std::cout << "P - b = " << P - b << std::endl;
  std::cout << "P * b = " << P * b << std::endl;
  std::cout << "P / b = " << P / b << std::endl;
  std::cout << "P % b = " << P % b << std::endl;
  double a = 2.0;
  std::cout << "a = " << a << std::endl;
  std::cout << "a + Q = " << a + Q << std::endl;
  std::cout << "a - Q = " << a - Q << std::endl;
  std::cout << "a * Q = " << a * Q << std::endl;
  std::cout << "a / Q = " << a / Q << std::endl;
  std::cout << "a % Q = " << a % Q << std::endl;

  std::polynomial<double> B;// = b;
  B = b;
  std::cout << "B = " << B << std::endl;
  std::cout << "P % B = " << P % B << std::endl;

  Q = {0.0, -2.0, 4.0, -6.0, 8.0, -12.0};
  std::cout << "Q = " << Q << std::endl;

  std::polynomial<double> P2;
  P2 = P;
  std::cout << "P2 = " << P2 << std::endl;
  std::cout << "P2 == P = " << (P2 == P) << std::endl;

  for (int i = 0; i <= 100; ++i)
  {
    double x = i * 0.1;
    std::cout << "P(" << x << ") = " << P(x) << std::endl;
  }

  std::polynomial<std::complex<double>>
  CP({std::complex<double>(0.0, -1.0),
      std::complex<double>(1.0, -2.0), 
      std::complex<double>(2.0, -3.0), 
      std::complex<double>(3.0, -4.0)});
  std::cout << "CP = " << CP << std::endl;
  std::cout << "CP * CP = " << CP * CP << std::endl;

  std::polynomial<int> IP({0, 1, 2, 3});
  std::cout << "IP = " << IP << std::endl;
  std::cout << "IP * IP = " << IP * IP << std::endl;

  std::array<double, 10> arr;
  P.eval(1.0, arr);
  std::cout << "P(" << 1.0 << ") =";
  for (int i = 0; i < arr.size(); ++i)
    std::cout << " " << arr[i];
  std::cout << std::endl;

  P.eval(1.0, arr.begin(), arr.end());
  std::cout << "P(" << 1.0 << ") =";
  for (auto iarr = arr.cbegin(); iarr != arr.cend(); ++iarr)
    std::cout << " " << *iarr;
  std::cout << std::endl;

  std::istringstream is("(-2.0, -1.0, 0.0)");
  std::polynomial<double> R;
  is >> R;
  std::cout << "R = " << R << std::endl;

  std::istringstream is2("(5.0)");
  std::polynomial<double> S;
  is2 >> S;
  std::cout << "S = " << S << std::endl;

  std::istringstream is3("42.0");
  std::polynomial<double> T;
  is3 >> T;
  std::cout << "T = " << T << std::endl;

  std::polynomial<double> u({1.0, 3.0, 3.0, 1.0});
  std::polynomial<double> v({1.0, 1.0});
  std::polynomial<double> q, r;
  std::divmod(u, v, q, r);
  std::cout << "u = " << u << std::endl;
  std::cout << "v = " << v << std::endl;
  std::cout << "q = " << q << std::endl;
  std::cout << "r = " << r << std::endl;

  std::polynomial<double> u1({1.0, -3.0, 3.0, -1.0});
  std::polynomial<double> v1({1.0, -2.0, 1.0});
  std::polynomial<double> q1, r1;
  std::divmod(u1, v1, q1, r1);
  std::cout << "u1 = " << u1 << std::endl;
  std::cout << "v1 = " << v1 << std::endl;
  std::cout << "q1 = " << q1 << std::endl;
  std::cout << "r1 = " << r1 << std::endl;

  std::polynomial<double> u2({1.0, 1.0});
  std::polynomial<double> v2({1.0, 3.0, 3.0, 1.0});
  std::polynomial<double> q2, r2;
  std::divmod(u2, v2, q2, r2);
  std::cout << "u2 = " << u2 << std::endl;
  std::cout << "v2 = " << v2 << std::endl;
  std::cout << "q2 = " << q2 << std::endl;
  std::cout << "r2 = " << r2 << std::endl;

  std::polynomial<double> u3({1.0, 0.0, 0.0, 1.0});
  std::polynomial<double> v3({1.0, 3.0, 3.0, 1.0});
  std::polynomial<double> q3, r3;
  std::divmod(u3, v3, q3, r3);
  std::cout << "u3 = " << u3 << std::endl;
  std::cout << "v3 = " << v3 << std::endl;
  std::cout << "q3 = " << q3 << std::endl;
  std::cout << "r3 = " << r3 << std::endl;

  std::cout << "P = " << P << std::endl;
  std::cout << "P' = " << P.derivative() << std::endl;
  std::cout << "I = " << P.integral(42.0) << std::endl;
  
  std::polynomial<std::polynomial<double>>
  Pp({{9.0, 8.0, 7.0}, {6.0, 5.0, 4.0}, {3.0, 2.0, 1.0}});
  std::cout << "Pp = " << Pp << std::endl;
  //std::cout << "Pp(-1) = " << Pp(-1.0) << std::endl;

  //std::array<double, 5> aaa{{1.1, 2.2, 3.3, 4.4, 5.5}};
  //std::cout << "aaa = " << std::polynomial_eval(aaa, 3.13149) << std::endl;
}

