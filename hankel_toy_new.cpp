/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o hankel_toy_new hankel_toy_new.cpp -L$HOME/bin/lib64 -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./hankel_toy_new > hankel_toy_new.txt

$HOME/bin/bin/g++ -std=gnu++17 -DNO_CBRT -o hankel_toy_new hankel_toy_new.cpp -lquadmath
./hankel_toy_new > hankel_toy_new.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <tuple>
//#include <bits/float128.h>
#include <ext/polynomial.h>
#include <bits/numeric_limits.h>

template<typename _Tp>
  struct __hankel_outer_param_t
  {
    __hankel_outer_param_t(std::complex<_Tp> __nu_in, std::complex<_Tp> __z_in);

    std::complex<_Tp> __zhat;
    std::complex<_Tp> __1dnsq;
    std::complex<_Tp> __num1d3;
    std::complex<_Tp> __num2d3;
    std::complex<_Tp> __p;
    std::complex<_Tp> __p2;
    std::complex<_Tp> __etm3h;
    std::complex<_Tp> __etrat;
    std::complex<_Tp> _Aip;
    std::complex<_Tp> __o4dp;
    std::complex<_Tp> _Aim;
    std::complex<_Tp> __o4dm;
    std::complex<_Tp> __od2p;
    std::complex<_Tp> __od0dp;
    std::complex<_Tp> __od2m;
    std::complex<_Tp> __od0dm;
  };


template<typename _Tp>
  struct __hankel_param_t
  {
    __hankel_param_t(std::complex<_Tp> __nu_in, std::complex<_Tp> __zhat_in);

    static constexpr auto _S_2pi   = _Tp{6.283185307179586476925286766559005768391Q};
    static constexpr auto _S_1d3   = _Tp{0.3333333333333333333333333333333333333333Q};
    static constexpr auto _S_2d3   = _Tp{0.6666666666666666666666666666666666666666Q};
    static constexpr auto _S_lncon = _Tp{0.2703100720721095879853420769762327577152Q}; // -(2/3)ln(2/3)
    static constexpr std::complex<_Tp> _S_j{0, 1};

    std::complex<_Tp> __nu;
    std::complex<_Tp> __zhat;
    std::complex<_Tp> __nup2;
    std::complex<_Tp> __num2;
    std::complex<_Tp> __num1d3;
    std::complex<_Tp> __num2d3;
    std::complex<_Tp> __num4d3;
    std::complex<_Tp> __w;
    std::complex<_Tp> __p;
    std::complex<_Tp> __p2;
    std::complex<_Tp> __xi;
    std::complex<_Tp> __zeta;
    std::complex<_Tp> __zetam3hf;
    std::complex<_Tp> __zetaphf;
    std::complex<_Tp> __zetamhf;
    std::complex<_Tp> __thing;
  };

template<typename _Tp>
  __hankel_param_t<_Tp>::__hankel_param_t(std::complex<_Tp> __nu_in,
					  std::complex<_Tp> __zhat_in)
  : __nu(__nu_in), __zhat(__zhat_in)
  {
    //using _Cmplx = std::complex<_Tp>;
    const auto _S_sqrt_max = __gnu_cxx::__sqrt_max<_Tp>();

    // If nu^2 can be computed without overflow.
    if (std::abs(__nu) <= _S_sqrt_max)
      {
	auto __nup2 = __nu * __nu;
	__num2 = _Tp{1} / __nup2;
	// Compute nu^(-1/3), nu^(-2/3), nu^(-4/3).
	__num1d3 = std::pow(__nu, -_S_1d3);
	__num2d3 = __num1d3 * __num1d3;
	__num4d3 = __num2d3 * __num2d3;
      }
    else
      std::__throw_runtime_error(__N("__hankel_param_t: "
				     "unable to compute nu^2"));

    if (std::real(__zhat) <= _Tp{1})
      {
	__w = std::sqrt((_Tp{1} + __zhat) * (_Tp{1} - __zhat));
	__p = _Tp{1} / __w;
	__p2 = __p * __p;
	__xi = std::log(_Tp{1} + __w) - std::log(__zhat) - __w;
	//auto __zetam3hf = _S_2d3 / __xi; // zeta^(-3/2)
	auto __logzeta = _S_2d3 * std::log(__xi) + _S_lncon; // ln(zeta)
	__zeta = std::exp(__logzeta);
	__zetaphf = std::sqrt(__zeta);
	__zetamhf = _Tp{1} / __zetaphf;
      }
    else
      {
	__w = std::sqrt((__zhat + _Tp{1}) * (__zhat - _Tp{1}));
	__p = _S_j / __w;
	__p2 = __p * __p;
	__xi = __w - std::acos(_Tp{1} / __zhat);
	__zetam3hf = _S_j * _S_2d3 / __xi; // i(-zeta)^(-3/2)
	auto __logmzeta = _S_2d3 * std::log(__xi) + _S_lncon; // ln(-zeta)
	auto __mzeta = std::exp(__logmzeta); // -zeta
	__zetaphf = -_S_j * std::sqrt(__mzeta); // -i(-zeta)^(1/2)
	__zetamhf = _Tp{1} / __zetaphf; // i(-zeta)^(-1/2)
	__zeta = -__mzeta;
      }
      __thing = std::pow(_Tp{4} * __zeta / __w, _Tp{0.25L});
  }

template<typename _Tp>
  std::complex<_Tp>
  get_zeta_old(std::complex<_Tp> __zhat)
  {
    using __cmplx = std::complex<_Tp>;

    static constexpr auto _S_2pi   = _Tp{6.283185307179586476925286766559005768391L};
    static constexpr auto _S_2d3   = _Tp{0.6666666666666666666666666666666666666667Q};
    static constexpr auto _S_lncon = _Tp{0.2703100720721095879853420769762327577152Q}; // -(2/3)ln(2/3)
    static constexpr __cmplx _S_j{0, 1};

    auto __rezhat = std::real(__zhat);
    auto __imzhat = std::imag(__zhat);

    // Compute 1 - zhat^2 and related constants.
    auto __w = __cmplx{_Tp{1} - (__rezhat - __imzhat) * (__rezhat + __imzhat),
			-_Tp{2} * __rezhat * __imzhat};
    __w = std::sqrt(__w);

    // Compute xi = ln(1+(1-zhat^2)^(1/2)) - ln(zhat) - (1-zhat^2)^(1/2)
    // using default branch of logarithm and square root.
    auto __xi = std::log(__cmplx{1} + __w) - std::log(__zhat) - __w;
    auto __zetam3hf = _S_2d3 / __xi;

    // Compute principal value of ln(xi) and then adjust imaginary part.
    auto __logxi = std::log(__xi);

    // Prepare to adjust logarithm of xi to appropriate Riemann sheet.
    auto __npi = _Tp{0};

    // Find adjustment necessary to get on proper Riemann sheet.
    if (__imzhat == _Tp{0})  // zhat is real.
      {
	if (__rezhat > _Tp{1})
	  __npi = _S_2pi;
      }
    else // zhat is not real.
      {
	// zhat is in upper half-plane.
	if (__imzhat > _Tp{0})
	  {
	    // xi lies in upper half-plane.
	    if (std::imag(__xi) > _Tp{0})
	      __npi = -_S_2pi;
	    else
	      __npi = +_S_2pi;
	  }
      }

    // Adjust logarithm of xi.
    __logxi += __npi * _S_j;

    // Compute ln(zeta), zeta, zeta^(+1/2), zeta^(-1/2).
    auto __logzeta = _S_2d3 * __logxi + _S_lncon;
    auto __zeta = std::exp(__logzeta);
    return __zeta;
  }

template<typename _Tp>
  std::complex<_Tp>
  get_zeta(std::complex<_Tp> __zhat)
  {
    using __cmplx = std::complex<_Tp>;

    static constexpr auto _S_2d3   = _Tp{0.6666666666666666666666666666666666666667Q};
    static constexpr auto _S_lncon = _Tp{0.2703100720721095879853420769762327577152Q}; // -(2/3)ln(2/3)
    static constexpr __cmplx _S_j{0, 1};

    if (__zhat == _Tp{0})
      return std::numeric_limits<_Tp>::infinity();
    else if (__zhat <= _Tp{1})
      {
	auto __w = std::sqrt((_Tp{1} + __zhat) * (_Tp{1} - __zhat));
	// Compute xi = ln(1 + (1 - zhat^2)^(1/2)) - ln(zhat) - (1 - zhat^2)^(1/2) = (2/3)(zeta)^(3/2)
	// using default branch of logarithm and square root.
	auto __xi = std::log(_Tp{1} + __w) - std::log(__zhat) - __w;

	auto __logxi = std::log(__xi);

	// Compute ln(zeta), zeta.
	auto __logzeta = _S_2d3 * __logxi + _S_lncon;
	auto __zeta = std::exp(__logzeta);
	return __zeta;
      }
    else
      {
	auto __w = std::sqrt((__zhat + _Tp{1}) * (__zhat - _Tp{1}));
	auto __xi = __w - std::acos(_Tp{1} / __zhat);

	auto __logxi = std::log(__xi);

	// Compute ln(-zeta), zeta.
	auto __logmzeta = _S_2d3 * __logxi + _S_lncon;
	auto __zeta = -std::exp(__logmzeta);
	return __zeta;
      }
  }

template<typename _Tp>
  _Tp
  get_zeta(_Tp __zhat)
  {
    static constexpr auto _S_2d3   = _Tp{0.6666666666666666666666666666666666666667Q};
    static constexpr auto _S_lncon = _Tp{0.2703100720721095879853420769762327577152Q}; // -(2/3)ln(2/3)
    if (__zhat == _Tp{0})
      return std::numeric_limits<_Tp>::infinity();
    else if (__zhat <= _Tp{1})
      {
	auto __w = std::sqrt((_Tp{1} + __zhat) * (_Tp{1} - __zhat));
	// Compute xi = ln(1 + (1 - zhat^2)^(1/2)) - ln(zhat) - (1 - zhat^2)^(1/2) = (2/3)(zeta)^(3/2)
	// using default branch of logarithm and square root.
	auto __xi = std::log(_Tp{1} + __w) - std::log(__zhat) - __w;
	//auto __zetam3hf = _S_2d3 / __xi;

	auto __logxi = std::log(__xi);

	// Compute ln(zeta), zeta.
	auto __logzeta = _S_2d3 * __logxi + _S_lncon;
	auto __zeta = std::exp(__logzeta);
	return __zeta;
      }
    else
      {
	auto __w = std::sqrt((__zhat + _Tp{1}) * (__zhat - _Tp{1}));
	// Compute xi = (zhat^2 - 1)^(1/2) - arcsec(zhat) = (2/3)(-zeta)^(3/2)
	// using default branch of logarithm and square root.
	auto __xi = __w - std::acos(_Tp{1} / __zhat);
	//auto __mzetam3hf = _S_2d3 / __xi;

	auto __logxi = std::log(__xi);

	// Compute ln(-zeta), zeta.
	auto __logmzeta = _S_2d3 * __logxi + _S_lncon;
	auto __zeta = -std::exp(__logmzeta);
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
	//auto indexpold = indexp;
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

    auto prec = std::numeric_limits<_Tp>::digits10;
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
    __gnu_cxx::_Polynomial<_Tp> upol1{_Tp{1}, _Tp{1}, _Tp{0.5Q}, _Tp{1}, -_Tp{0.5Q}};
    __gnu_cxx::_Polynomial<_Tp> upol2{+_Tp{0.125Q}, _Tp{1}, -_Tp{0.625Q}};
    __gnu_cxx::_Polynomial<_Tp> vpol1{_Tp{1}, -_Tp{0.5Q}, _Tp{1}, +_Tp{0.5Q}};
    __gnu_cxx::_Polynomial<_Tp> vpol2{_Tp{1}, _Tp{1}, -_Tp{1}, _Tp{1}, +_Tp{1}};
    __gnu_cxx::_Polynomial<_Tp> u{_Tp{1}};
    std::vector<__gnu_cxx::_Polynomial<_Tp>> uvec;
    __gnu_cxx::_Polynomial<_Tp> v{_Tp{1}};
    std::vector<__gnu_cxx::_Polynomial<_Tp>> vvec;
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
	for (std::size_t i = 0; i <= u.degree(); ++i)
	  if (u.coefficient(i) != 0)
	    //uentry[i].push_back(std::make_tuple(ku, i, u.coefficient(i)));
	    uentry[i].emplace_back(ku, i, u.coefficient(i));
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
	for (std::size_t i = 0; i <= v.degree(); ++i)
	  if (v.coefficient(i) != 0)
	    //ventry[i].push_back(std::make_tuple(kv, i, v.coefficient(i)));
	    ventry[i].emplace_back(kv, i, v.coefficient(i));
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

    auto k_max = (uvec.size() - 2) / 2;
    k_max = std::min(k_max, (vvec.size() - 2) / 2);
    k_max = std::min(k_max, lambda.size() - 1);
    k_max = std::min(k_max, mu.size() - 1);
    std::cout << "\nkmax = " << k_max << '\n';
    k_max = 6;
    std::cout << "\nkmax = " << k_max << '\n';
    std::cout << "U_k\n";
    for (std::size_t k = 0; k <= k_max; ++k)
      std::cout << uvec[k] << '\n';
    std::cout << "V_k\n";
    for (std::size_t k = 0; k <= k_max; ++k)
      std::cout << vvec[k] << '\n';

    // Try to build A, B, C, D functions...
    // These are polynomials in p and zeta^{-3/2}
    std::cout << '\n'
	      << std::setw(width) << "z" << ' '
	      << std::setw(width) << "zeta" << ' '
	      << std::setw(width) << "zeta^(3/2)" << ' '
	      << std::setw(width) << "thing" << ' '
	      << std::setw(width) << "p" << ' ';
    for (std::size_t k = 0; k < k_max; ++k)
      std::cout << std::setw(width-1) << "A_" << k << ' ';
    for (std::size_t k = 0; k < k_max; ++k)
      std::cout << std::setw(width-1) << "B_" << k << ' ';
    for (std::size_t k = 0; k < k_max; ++k)
      std::cout << std::setw(width-1) << "C_" << k << ' ';
    for (std::size_t k = 0; k < k_max; ++k)
      std::cout << std::setw(width-1) << "D_" << k << ' ';
    std::cout << '\n';
    for (int i = 0; i <= 2000; ++i)
      {
	auto z = _Tp(i * 0.01Q);
	auto zeta = get_zeta<_Tp>(z);
	auto thing = std::sqrt(std::sqrt(4 * zeta / ((_Tp{1} + z) * (_Tp{1} - z))));
	auto p = std::abs(_Tp{1} / std::sqrt((_Tp{1} + z) * (_Tp{1} - z)));
	if (std::abs(z) > _Tp{1})
	  p = std::abs(_Tp{1} / std::sqrt((z + _Tp{1}) * (z - _Tp{1})));
	auto t = _Tp{1.5L} / std::pow(std::abs(zeta), _Tp{1.5L});
	std::cout << std::setw(width) << z << ' '
		  << std::setw(width) << zeta << ' '
		  << std::setw(width) << std::pow(std::abs(zeta), _Tp{-1.5L}) << ' '
		  << std::setw(width) << thing << ' '
		  << std::setw(width) << p << ' ';
	for (std::size_t k = 0; k < k_max; ++k)
	  {
	    auto tj = _Tp{1};
	    auto A = _Tp{0};
	    for (std::size_t j = 0; j <= 2 * k; ++j)
	      {
//std::cout << "\nuvec[" << 2 * k - j << "]: " << uvec[2 * k - j] << '\n';;
		A += tj * mu[j] * uvec[2 * k - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << A << ' ';
	  }
	for (std::size_t k = 0; k < k_max; ++k)
	  {
	    auto tj = _Tp{1};
	    auto B = _Tp{0};
	    for (std::size_t j = 0; j <= 2 * k + 1; ++j)
	      {
//std::cout << "\nuvec[" << 2 * k + 1 - j << "]: " << uvec[2 * k + 1 - j] << '\n';;
		B += tj * lambda[j] * uvec[2 * k + 1 - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << B << ' ';
	  }
	for (std::size_t k = 0; k < k_max; ++k)
	  {
	    auto tj = _Tp{1};
	    auto C = _Tp{0};
	    for (std::size_t j = 0; j <= 2 * k + 1; ++j)
	      {
//std::cout << "\nvvec[" << 2 * k + 1 - j << "]: " << vvec[2 * k + 1 - j] << '\n';
		C += tj * mu[j] * vvec[2 * k + 1 - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << C << ' ';
	  }
	for (std::size_t k = 0; k < k_max; ++k)
	  {
	    auto tj = _Tp{1};
	    auto D = _Tp{0};
	    for (std::size_t j = 0; j <= 2 * k; ++j)
	      {
//std::cout << "\nvvec[" << 2 * k - j << "]: " << vvec[2 * k - j] << '\n';
		D += tj * lambda[j] * vvec[2 * k - j](p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << D << ' ';
	  }
	std::cout << '\n';
      }

    auto nu = _Tp{1};
    for (int i = 0; i <= 2000; ++i)
      {
	auto zhat = _Tp(i * 0.01Q);
	auto parm = __hankel_param_t<_Tp>(nu, zhat);
	std::cout
	  << ' ' << parm.__zhat
	  << ' ' << parm.__zeta
	  << ' ' << parm.__zetam3hf
	  << ' ' << parm.__thing
	  << ' ' << parm.__p;
	//  << ' ' << parm.
	//  << ' ' << parm.
	//  << ' ' << parm.
	auto t = _Tp{1.5L} / std::pow(parm.__zeta, _Tp{1.5L});
	for (std::size_t k = 0; k < k_max; ++k)
	  {
	    decltype(t) tj = 1;
	    decltype(t) A = 0;
	    for (std::size_t j = 0; j <= 2 * k; ++j)
	      {
//std::cout << "\nuvec[" << 2 * k - j << "]: " << uvec[2 * k - j] << '\n';;
		A += tj * mu[j] * uvec[2 * k - j](parm.__p);
		tj *= t;
	      }
	    std::cout << std::setw(width) << A << ' ';
	  }
	  std::cout << '\n';
      }
  }

int
main()
{
  //std::cout << "\nRunning float\n-------------\n";
  //run_toy<float>();

  std::cout << "\nRunning double\n--------------\n";
  run_toy<double>();

  //std::cout << "\nRunning long double\n-------------------\n";
  //run_toy<long double>();

  //std::cout << "\nSkipping __float128\n-------------------\n";
  //run_toy<__float128>();
}
