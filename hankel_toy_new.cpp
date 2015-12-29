// $HOME/bin/bin/g++ -std=gnu++14 -o hankel_toy_new hankel_toy_new.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./hankel_toy_new > hankel_toy_new.txt

// g++ -std=gnu++14 -o hankel_toy_new hankel_toy_new.cpp -lquadmath

// ./hankel_toy_new > hankel_toy_new.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include "float128.h"
#include "complex128.h"
#include "polynomial"

template<typename _Tp>
  void
  run_toy()
  {
    // Figure out the array indexing in the Hankel asymptotic series.
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

    // Figure out a formula for the array indexing in the Hankel asymptotic series.
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

    // Build the lambda_k and mu_k ratios for the asymptotic series.
    std::cout << '\n' << std::setw(width) << "lambda\t" << std::setw(width) << "mu\n";
    _Tp lambda = _Tp{1};
    _Tp mu = -_Tp{1};
    for (int s = 1; s <= 50; ++s)
      {
	std::cout << std::setw(width) << lambda << '\t' << std::setw(width) << mu << '\n';
	lambda *= _Tp{1} * (6 * s - 3) * (6 * s - 3) * (6 * s - 1)
		/ ((2 * s - 1) * s * 144);
	mu = -(6 * s + 1) * lambda / (6 * s - 1);
      }

    // Build the Debye polynomials.
    std::polynomial<_Tp> upol1{_Tp{1}, _Tp{1}, _Tp{0.5L}, _Tp{1}, -_Tp{0.5L}};
    std::polynomial<_Tp> upol2{+_Tp{0.125L}, _Tp{1}, -_Tp{0.625L}};
    std::polynomial<_Tp> vpol1{_Tp{1}, -_Tp{0.5L}, _Tp{1}, +_Tp{0.5L}};
    std::polynomial<_Tp> vpol2{_Tp{1}, _Tp{1}, -_Tp{1}, _Tp{1}, +_Tp{1}};
    std::polynomial<_Tp> u{_Tp{1}};
    std::vector<std::polynomial<_Tp>> uvec;
    std::polynomial<_Tp> v{_Tp{1}};
    std::vector<std::polynomial<_Tp>> vvec;
    for (auto k = 1; k <= 20; ++k)
      {
	uvec.push_back(u);
	vvec.push_back(v);
	v = vpol1 * u + vpol2 * u.derivative();
	u = upol1 * u.derivative() + (upol2 * u).integral(_Tp{0});
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

    // Write the Debye polynomials in reverse as they are stored in the code.
    // This allows application of Horner's rule by traversing the coefficients in order.
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

    // Try to build A, B, C, D functions...
    // These are polynomials in p and zeta^{-3/2}
    for (int i = 0; i <= 100; ++i)
      {
	auto z = std::complex<_Tp>(i * 0.01);
	auto w = std::sqrt(_Tp{1} - z * z);
	auto xi = std::log((_Tp{1} + w) / z) - w;
	auto zeta = std::pow(_Tp{3} * xi / _Tp{2}, _Tp{2} / _Tp{3});
	auto p = _Tp{1} / w;
	auto zetam32 = _Tp{2} / (_Tp{3} * xi);
      }
  }

int
main()
{
  std::cout << "\nfloat\n-----\n";
  run_toy<float>();

  std::cout << "\ndouble\n------\n";
  run_toy<double>();

  std::cout << "\nlong double\n-----------\n";
  run_toy<long double>();

  std::cout << "\n__float128\n----------\n";
  run_toy<__float128>();
}
