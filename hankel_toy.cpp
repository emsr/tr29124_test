// $HOME/bin_specfun/bin/g++ -std=gnu++14 -o hankel_toy hankel_toy.cpp /home/ed/tr29124_test/gslextras/Fresnel/fresnel.c -lquadmath -L/usr/local/lib -lgsl -lgslcblas -ljacobi

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./hankel_toy > hankel_toy.txt

// g++ -std=gnu++14 -o hankel_toy hankel_toy.cpp -L$HOME/bin/lib64

// ./hankel_toy > hankel_toy.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include "polynomial.h"
#include <bits/float128.h>

template<typename _Tp>
  void
  run_toy()
  {
    constexpr auto _S_1d6 = _Tp{1} / _Tp{6};
    constexpr auto _S_5d6 = _Tp{5} / _Tp{6};

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

    auto prec = std::numeric_limits<_Tp>::max_digits10;
    auto width = prec + 6;

    std::cout.precision(prec);
    std::cout << std::scientific;
    std::cout << std::showpoint;

    bool old = false;

    std::cout << '\n' << std::setw(width) << "lambda\t" << std::setw(width) << "mu\n";
    _Tp lambda = _Tp{1};
    _Tp mu = -_Tp{1};
    for (int s = 1; s <= 50; ++s)
      {
	std::cout << std::setw(width) << lambda << '\t' << std::setw(width) << mu << '\n';
	// Turn this into a recursion:
	// for (int m = 2 * s + 1; m <= 6 * s - 1; m += 2)
	//   numer *= m;
	if (old)
	  {
	    lambda *= _Tp{1} * (6 * s - 5) * (6 * s - 3) * (6 * s - 1)
		  / ((2 * s - 1) * s * 144);
	    mu = -(6 * s + 1) * lambda / (6 * s - 1);
	  }
	else
	  { // FIXME: Why this no work!?!?!
	    // It's a subtly different thing!  Look at the nenom!
	    // Derive a new one!
	    //lambda *= _Tp(s - 1) / _Tp{2} + _Tp{5} / _Tp(72 * s);
	    //mu = -lambda * _Tp(s + _S_1d6) / _Tp(s - _S_1d6);
	    lambda *= _Tp{3} * _Tp(s - _S_5d6) * _Tp(s - _S_1d6)
		   / (s * 4);
	    mu = -lambda * _Tp(s + _S_1d6) / _Tp(s - _S_1d6);
	  }
	if (std::isnan(lambda) || std::isinf(lambda)
	 || std::isnan(mu) || std::isinf(mu))
	  break;
      }

    __gnu_cxx::_Polynomial<_Tp> upol1{_Tp{0}, _Tp{0}, _Tp{0.5Q}, _Tp{0}, -_Tp{0.5Q}};
    __gnu_cxx::_Polynomial<_Tp> upol2{+_Tp{0.125Q}, _Tp{0}, -_Tp{0.625Q}};
    __gnu_cxx::_Polynomial<_Tp> vpol1{_Tp{0}, -_Tp{0.5Q}, _Tp{0}, +_Tp{0.5Q}};
    __gnu_cxx::_Polynomial<_Tp> vpol2{_Tp{0}, _Tp{0}, -_Tp{1}, _Tp{0}, +_Tp{1}};
    __gnu_cxx::_Polynomial<_Tp> u{_Tp{1}};
    std::vector<__gnu_cxx::_Polynomial<_Tp>> uvec;
    __gnu_cxx::_Polynomial<_Tp> v{_Tp{1}};
    std::vector<__gnu_cxx::_Polynomial<_Tp>> vvec;
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

    std::vector<std::vector<std::tuple<int, int, _Tp>>> uentry;
    auto ku = 0;
    for (const auto & u : uvec)
      {
	uentry.resize(u.degree() + 1);
	for (int i = 0; i <= u.degree(); ++i)
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
    std::vector<std::vector<std::tuple<int, int, _Tp>>> ventry;
    auto kv = 0;
    for (const auto & v : vvec)
      {
	ventry.resize(v.degree() + 1);
	for (int i = 0; i <= v.degree(); ++i)
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
