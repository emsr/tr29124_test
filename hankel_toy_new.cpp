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
  _Tp
  get_zeta(_Tp __zhat)
  {
    static constexpr auto _S_2d3   = _Tp{0.6666666666666666666666666666666666666666L};
    static constexpr auto _S_lncon = _Tp{0.2703100720721095879853420769762327577152L}; // -(2/3)ln(2/3)
    //if (std::abs(__zhat) <= _Tp{1})
    if (__zhat == _Tp{0})
      return std::numeric_limits<_Tp>::infinity();
    else if (__zhat <= _Tp{1})
      {
	auto __ztemp = std::sqrt((_Tp{1} + __zhat) * (_Tp{1} - __zhat)); // __w
	// Compute xi = ln(1 + (1 - zhat^2)^(1/2)) - ln(zhat) - (1 - zhat^2)^(1/2) = (2/3)(zeta)^(3/2)
	// using default branch of logarithm and square root.
	auto __xi = std::log(_Tp{1} + __ztemp) - std::log(__zhat) - __ztemp;
	//auto __zetam3hf = _S_2d3 / __xi;

	auto __lnxi = std::log(__xi);

	// Compute ln(zeta), zeta.
	auto __lnzeta = _S_2d3 * __lnxi + _S_lncon;
	auto __zeta = std::exp(__lnzeta);
	return __zeta;
      }
    else
      {
	auto __ztemp = std::sqrt((__zhat + _Tp{1}) * (__zhat - _Tp{1})); // __w
	// Compute xi = (zhat^2 - 1)^(1/2) - arcsec(zhat) = (2/3)(-zeta)^(3/2)
	// using default branch of logarithm and square root.
	auto __xi = __ztemp - std::acos(_Tp{1} / __zhat);
	//auto __mzetam3hf = _S_2d3 / __xi;

	auto __lnxi = std::log(__xi);

	// Compute ln(zeta), zeta.
	auto __lnmzeta = _S_2d3 * __lnxi + _S_lncon;
	auto __zeta = -std::exp(__lnmzeta);
	return __zeta;
      }
  }

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
    std::vector<_Tp> lambda;
    std::vector<_Tp> mu;
    lambda.push_back(_Tp{1});
    mu.push_back(-_Tp{1});
    for (int s = 1; s <= 50; ++s)
      {
	std::cout << std::setw(width) << lambda.back() << '\t'
		  << std::setw(width) << mu.back() << '\n';
	lambda.push_back(lambda.back() * (6 * s - 3) * (6 * s - 3) * (6 * s - 1)
			 / ((2 * s - 1) * s * 144));
	mu.push_back(-(6 * s + 1) * lambda.back() / (6 * s - 1));
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
    std::cout << '\n'
	      << std::setw(width) << "z" << ' '
	      << std::setw(width) << "zeta" << ' '
	      << std::setw(width) << "zeta^(3/2)" << ' '
	      << std::setw(width) << "thing" << ' '
	      << std::setw(width) << "p" << ' ';
    std::cout << '\n';
    for (int i = 0; i <= 2000; ++i)
      {
	auto z = i * 0.01L;
	auto zeta = get_zeta<_Tp>(z);
	auto thing = std::sqrt(std::sqrt(4 * zeta / ((_Tp{1} + z) * (_Tp{1} - z))));
	auto p = _Tp{1} / std::sqrt((_Tp{1} + z) * (_Tp{1} - z));
	auto t = _Tp{1.5L} / std::pow(std::abs(zeta), _Tp{1.5L});
	std::cout << std::setw(width) << z << ' '
		  << std::setw(width) << zeta << ' '
		  << std::setw(width) << std::pow(std::abs(zeta), _Tp{-1.5L}) << ' '
		  << std::setw(width) << thing << ' '
		  << std::setw(width) << p << ' ';
	for (int k = 0; k < 6; ++k)
	  {
	    auto tj = _Tp{1};
	    auto A = _Tp{0};
	    for (int j = 0; j <= 2 * k; ++j)
	      {
		A += tj * mu[j] * uvec[2 * k - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << A << ' ';
	  }
	for (int k = 0; k < 6; ++k)
	  {
	    auto tj = _Tp{1};
	    auto B = _Tp{0};
	    for (int j = 0; j <= 2 * k + 1; ++j)
	      {
		B += tj * lambda[j] * uvec[2 * k + 1 - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << B << ' ';
	  }
	for (int k = 0; k < 6; ++k)
	  {
	    auto tj = _Tp{1};
	    auto C = _Tp{0};
	    for (int j = 0; j <= 2 * k + 1; ++j)
	      {
		C += tj * mu[j] * vvec[2 * k + 1 - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << C << ' ';
	  }
	for (int k = 0; k < 6; ++k)
	  {
	    auto tj = _Tp{1};
	    auto D = _Tp{0};
	    for (int j = 0; j <= 2 * k; ++j)
	      {
		D += tj * lambda[j] * vvec[2 * k - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << D << ' ';
	  }
	std::cout << '\n';
      }
  }

int
main()
{
  std::cout << "\nRunning float\n-------------\n";
  run_toy<float>();

  std::cout << "\nRunning double\n--------------\n";
  run_toy<double>();

  std::cout << "\nRunning long double\n-------------------\n";
  run_toy<long double>();

  std::cout << "\nSkipping __float128\n-------------------\n";
  //run_toy<__float128>();
}
