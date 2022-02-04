/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>

#include <emsr/float128_io.h>
#include <emsr/polynomial.h>

template<typename Tp>
  void
  run_toy()
  {
    constexpr auto s_1d6 = Tp{1} / Tp{6};
    constexpr auto s_5d6 = Tp{5} / Tp{6};

    auto index = 0;
    auto indexp = 0;
    for (int k = 0; k <= 20; ++k)
      {
	std::cout << '\n';
	std::cout << "k      = " << k << '\n';
	std::cout << "index  = " << index << '\n';
	std::cout << "indexp = " << indexp << '\n';
	//auto indexpold = indexp;
	index += 2;
	indexp += 2;
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

    auto prec = std::numeric_limits<Tp>::max_digits10;
    auto width = prec + 6;

    std::cout.precision(prec);
    std::cout << std::scientific;
    std::cout << std::showpoint;

    bool old = false;

    std::cout << '\n' << std::setw(width) << "lambda\t" << std::setw(width) << "mu\n";
    Tp lambda = Tp{1};
    Tp mu = -Tp{1};
    for (int s = 1; s <= 50; ++s)
      {
	std::cout << std::setw(width) << lambda << '\t' << std::setw(width) << mu << '\n';
	// Turn this into a recursion:
	// for (int m = 2 * s + 1; m <= 6 * s - 1; m += 2)
	//   numer *= m;
	if (old)
	  {
	    lambda *= Tp{1} * (6 * s - 5) * (6 * s - 3) * (6 * s - 1)
		  / ((2 * s - 1) * s * 144);
	    mu = -(6 * s + 1) * lambda / (6 * s - 1);
	  }
	else
	  { // FIXME: Why this no work!?!?!
	    // It's a subtly different thing!  Look at the nenom!
	    // Derive a new one!
	    //lambda *= Tp(s - 1) / Tp{2} + Tp{5} / Tp(72 * s);
	    //mu = -lambda * Tp(s + s_1d6) / Tp(s - s_1d6);
	    lambda *= Tp{3} * Tp(s - s_5d6) * Tp(s - s_1d6)
		   / (s * 4);
	    mu = -lambda * Tp(s + s_1d6) / Tp(s - s_1d6);
	  }
	if (std::isnan(lambda) || std::isinf(lambda)
	 || std::isnan(mu) || std::isinf(mu))
	  break;
      }

    emsr::Polynomial<Tp> upol1{Tp{0}, Tp{0}, Tp{0.5Q}, Tp{0}, -Tp{0.5Q}};
    emsr::Polynomial<Tp> upol2{+Tp{0.125Q}, Tp{0}, -Tp{0.625Q}};
    emsr::Polynomial<Tp> vpol1{Tp{0}, -Tp{0.5Q}, Tp{0}, +Tp{0.5Q}};
    emsr::Polynomial<Tp> vpol2{Tp{0}, Tp{0}, -Tp{1}, Tp{0}, +Tp{1}};
    emsr::Polynomial<Tp> u{Tp{1}};
    std::vector<emsr::Polynomial<Tp>> uvec;
    emsr::Polynomial<Tp> v{Tp{1}};
    std::vector<emsr::Polynomial<Tp>> vvec;
    for (auto k = 1; k <= 20; ++k)
      {
	uvec.push_back(u);
	vvec.push_back(v);
	v = vpol1 * u + vpol2 * u.derivative();
	u = upol1 * u.derivative() + (upol2 * u).integral(0.0);
	v += u;
      }
    std::cout << "\nu\n";
    for (const auto & u : uvec)
      std::cout << u << '\n';
    std::cout << "\nv\n";
    for (const auto & v : vvec)
      std::cout << v << '\n';

    std::vector<std::vector<std::tuple<int, int, Tp>>> uentry;
    auto ku = 0;
    for (const auto & u : uvec)
      {
	uentry.resize(u.degree() + 1);
	for (std::size_t i = 0; i <= u.degree(); ++i)
	  if (u.coefficient(i) != 0)
	    uentry[i].push_back(std::make_tuple(ku, i, u.coefficient(i)));
	++ku;
      }
    std::cout << "\nuentry\n";
    auto iu = 0;
    for (const auto & u : uentry)
      {
	for (const auto & c : u)
	  std::cout << ' ' << std::setw(3) << ++iu
		    << ' ' << std::setw(3) << std::get<0>(c)
		    << ' ' << std::setw(3) << std::get<1>(c)
		    << ' ' << std::setw(width) << std::get<2>(c) << '\n';
      }
    std::vector<std::vector<std::tuple<int, int, Tp>>> ventry;
    auto kv = 0;
    for (const auto & v : vvec)
      {
	ventry.resize(v.degree() + 1);
	for (std::size_t i = 0; i <= v.degree(); ++i)
	  if (v.coefficient(i) != 0)
	    ventry[i].push_back(std::make_tuple(kv, i, v.coefficient(i)));
	++kv;
      }
    std::cout << "\nventry\n";
    auto iv = 0;
    for (const auto & v : ventry)
      {
	for (const auto & c : v)
	  std::cout << ' ' << std::setw(3) << ++iv
		    << ' ' << std::setw(3) << std::get<0>(c)
		    << ' ' << std::setw(3) << std::get<1>(c)
		    << ' ' << std::setw(width) << std::get<2>(c) << '\n';
      }

    std::cout << "\nu\n";
    for (const auto& u : uvec)
      for (auto c = u.crbegin(); c != u.crend(); ++c)
	if (*c != 0)
	  std::cout << std::setw(width) << *c << '\n';

    std::cout << "\nv\n";
    for (const auto& v : vvec)
      for (auto c = v.crbegin(); c != v.crend(); ++c)
	if (*c != 0)
	  std::cout << std::setw(width) << *c << '\n';

    // Try these:  << std::showpos << std::uppercase << std::hexfloat << std::showpos
  }

int
main()
{
  //std::cout << "\nfloat\n=====\n";
  //run_toy<float>();

  //std::cout << "\ndouble\n======\n";
  //run_toy<double>();

  std::cout << "\nlong double\n===========\n";
  run_toy<long double>();

  //std::cout << "\n__float128\n==========\n";
  //run_toy<__float128>();
}
