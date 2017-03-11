
#include <iostream>
#include <iomanip>
#include <cmath>

#include "trapezoid_integral.h"
#include "ext/polynomial.h"

template<typename Tp>
  void
  test_integral(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout.flags(std::ios::showpoint);
    const auto w = 8 + std::cout.precision();

    const auto PI = __gnu_cxx::__const_pi(proto);

    auto sin2 = [](Tp x) ->Tp { Tp s = std::sin(x); return s * s; };
    auto cos2 = [](Tp x) ->Tp { Tp c = std::cos(x); return c * c; };
    auto j1 = [](Tp x) ->Tp { return std::cyl_bessel_j(Tp{1}, x); };
    auto foo = [](Tp x) ->Tp { return (Tp{1} - x) * std::exp(-x / Tp{2}); };
    auto foonum = [](Tp x) ->Tp { return (Tp{1} - x); };
    auto funk1 = [PI](Tp x) ->Tp { return std::cos(x) / std::sqrt(x * (PI - x)); };
    auto funk1num = [](Tp x) ->Tp { return std::cos(x); };
    auto funk2 = [PI](Tp x) ->Tp { return (Tp{2} + std::sin(x)) / std::sqrt(x * (PI - x)); };
    auto funk2num = [](Tp x) ->Tp { return Tp{2} + std::sin(x); };
    auto one = [](Tp) ->Tp { return Tp{1}; };
    auto ex = [](Tp x) ->Tp { return x; };
    __gnu_cxx::_Polynomial<Tp> poly1({1.0l, -0.5l, -3.5l, 2.0l});

    auto fun = [](Tp x) -> Tp { return std::sin(x); };
    using fun_t = decltype(fun);
    __gnu_test::trapezoid_integral<fun_t, Tp> mq(fun, Tp{0}, PI, Tp{0.0000001});
    std::cout << mq() << '\n';

    Tp a = Tp{0};
    Tp b = PI;
    Tp err = Tp{0.0000000001};

    __gnu_test::trapezoid_integral<decltype(one), Tp> t0(one, a, b, err);
    Tp a0 = t0();
    Tp e0 = b - a;
    std::cout << "one     : "
	      << ' ' << std::setw(w) << a0
	      << ' ' << std::setw(w) << e0
	      << ' ' << std::setw(w) << a0 - e0
	      << ' ' << std::setw(w) << t0.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(ex), Tp> t1(ex, a, b, err);
    Tp a1 = t1();
    Tp e1 = (b * b - a * a) / Tp{2};
    std::cout << "ex      : "
	      << ' ' << std::setw(w) << a1
	      << ' ' << std::setw(w) << e1
	      << ' ' << std::setw(w) << a1 - e1
	      << ' ' << std::setw(w) << t1.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(cos2), Tp> t2(cos2, a, b, err);
    Tp a2 = t2();
    Tp e2 = PI / Tp{2};
    std::cout << "cos2    : "
	      << ' ' << std::setw(w) << a2
	      << ' ' << std::setw(w) << e2
	      << ' ' << std::setw(w) << a2 - e2
	      << ' ' << std::setw(w) << t2.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(sin2), Tp> t3(sin2, a, b, err);
    Tp a3 = t3();
    Tp e3 = PI / Tp{2};
    std::cout << "sin2    : "
	      << ' ' << std::setw(w) << a3
	      << ' ' << std::setw(w) << e3
	      << ' ' << std::setw(w) << a3 - e3
	      << ' ' << std::setw(w) << t3.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(j1), Tp> t4(j1, a, b, err);
    Tp a4 = t4();
    Tp e4 = std::cyl_bessel_j(Tp{0}, Tp{0}) - std::cyl_bessel_j(Tp{0}, PI);
    std::cout << "j1      : "
	      << ' ' << std::setw(w) << a4
	      << ' ' << std::setw(w) << e4
	      << ' ' << std::setw(w) << a4 - e4
	      << ' ' << std::setw(w) << t4.abs_error() << '\n';

    a = Tp{0};
    b = Tp{10} * PI;
    __gnu_test::trapezoid_integral<decltype(foo), Tp> t5(foo, a, b, err);
    Tp a5 = t5();
    Tp e5 = Tp{2} * (Tp{1} + b) * std::exp(-b / Tp{2})
	  - Tp{2} * (Tp{1} + a) * std::exp(-a / Tp{2});
    std::cout << "foo     : "
	      << ' ' << std::setw(w) << a5
	      << ' ' << std::setw(w) << e5
	      << ' ' << std::setw(w) << a5 - e5
	      << ' ' << std::setw(w) << t5.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(foonum), Tp> t5n(foonum, a, b, err);
    Tp a5n = t5n();
    Tp e5n = b * (Tp{1} - b / Tp{2})
	   - a * (Tp{1} - a / Tp{2});
    std::cout << "foonum  : "
	      << ' ' << std::setw(w) << a5n
	      << ' ' << std::setw(w) << e5n
	      << ' ' << std::setw(w) << a5n - e5n
	      << ' ' << std::setw(w) << t5n.abs_error() << '\n';

    __gnu_test::trapezoid_integral<__gnu_cxx::_Polynomial<Tp>, Tp> t6(poly1, a, b, err);
    Tp a6 = t6();
    Tp e6 = poly1.integral()(b) - poly1.integral()(a);
    std::cout << "poly1   : "
	      << ' ' << std::setw(w) << a6
	      << ' ' << std::setw(w) << e6
	      << ' ' << std::setw(w) << a6 - e6
	      << ' ' << std::setw(w) << t6.abs_error() << '\n';

    a = Tp{0};
    b = PI;
    __gnu_test::trapezoid_integral<decltype(funk1), Tp> t7(funk1, a, b, err);
    Tp a7 = t7();
    Tp e7 = Tp{0};
    std::cout << "funk1   : "
	      << ' ' << std::setw(w) << a7
	      << ' ' << std::setw(w) << e7
	      << ' ' << std::setw(w) << a7 - e7
	      << ' ' << std::setw(w) << t7.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(funk1num), Tp> t7n(funk1num, a, b, err);
    Tp a7n = t7n();
    Tp e7n = Tp{0};
    std::cout << "funk1num: "
	      << ' ' << std::setw(w) << a7n
	      << ' ' << std::setw(w) << e7n
	      << ' ' << std::setw(w) << a7n - e7n
	      << ' ' << std::setw(w) << t7n.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(funk2), Tp> t8(funk2, a, b, err);
    Tp a8 = t8();
    Tp e8 = Tp{0};
    std::cout << "funk2   : "
	      << ' ' << std::setw(w) << a8
	      << ' ' << std::setw(w) << e8
	      << ' ' << std::setw(w) << a8 - e8
	      << ' ' << std::setw(w) << t8.abs_error() << '\n';

    __gnu_test::trapezoid_integral<decltype(funk2num), Tp> t8n(funk2num, a, b, err);
    Tp a8n = t8n();
    Tp e8n = Tp{0};
    std::cout << "funk2num: "
	      << ' ' << std::setw(w) << a8n
	      << ' ' << std::setw(w) << e8n
	      << ' ' << std::setw(w) << a8n - e8n
	      << ' ' << std::setw(w) << t8n.abs_error() << '\n';
  }

int
main()
{
  test_integral<double>();
}
