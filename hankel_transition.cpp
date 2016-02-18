// $HOME/bin_specfun/bin/g++ -g -std=gnu++1z -o hankel_transition hankel_transition.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./hankel_transition > hankel_transition.txt

// g++ -std=gnu++14 -DNO_CBRT -o hankel_transition hankel_transition.cpp -lquadmath

// ./hankel_transition > hankel_transition.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>
//#include "float128.h"
#include "polynomial.h"
#include "rational.h"
//#include "numeric_limits.h"

template<typename _Tp>
  void
  hankel_transition()
  {
    auto sign = [](int s, int r){return (s + r) % 2 == 1 ? -1 : +1; };
    using rational = _Tp;

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = std::cout.precision() + 6;

    __gnu_cxx::_Polynomial<int> poo({0, -2});
    std::vector<__gnu_cxx::_Polynomial<int>> phi;
    std::vector<__gnu_cxx::_Polynomial<int>> psi;
    phi.push_back(__gnu_cxx::_Polynomial<int>(1));
    phi.push_back(__gnu_cxx::_Polynomial<int>(0));
    phi.push_back(__gnu_cxx::_Polynomial<int>({0, -2}));
    psi.push_back(__gnu_cxx::_Polynomial<int>(0));
    psi.push_back(__gnu_cxx::_Polynomial<int>(-1));
    psi.push_back(__gnu_cxx::_Polynomial<int>(0));
    for (int s = 3; s <= 18; ++s)
      {
	phi.push_back(poo * phi[s - 2] -2 * (s - 2) * phi[s - 3]);
	psi.push_back(poo * psi[s - 2] -2 * (s - 2) * psi[s - 3]);
      }
    std::cout << "\nphi polynomial\n";
    for (const auto& p : phi)
      std::cout << p << '\n';
    std::cout << "\npsi polynomial\n";
    for (const auto& p : psi)
      std::cout << p << '\n';
    std::vector<__gnu_cxx::_Polynomial<rational>> A(6);
    std::vector<__gnu_cxx::_Polynomial<rational>> B(6);
    for (int s = 0; s < 6; ++s)
      {
	for (int r = 0; r <= s; ++r)
	  A[s] += phi[2 * s + r] * __gnu_cxx::_Polynomial<rational>(sign(r, s), r);
	std::cout << "A_" << s << ": " << A[s] << '\n';
	for (int r = 0; r <= s; ++r)
	  B[s] += psi[2 * s + r] * __gnu_cxx::_Polynomial<rational>(sign(r, s), r);
	std::cout << "B_" << s << ": " << B[s] << '\n';
      }

    std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::_Rational<long long>>>
    P
    {
      {{1, 1}},
      {{}, {-1, 5}},
      {{}, {}, {3, 35}, {}, {}, {-9, 100}},
      {{-1, 225}, {}, {}, {-173, 3150}, {}, {}, {957, 7000}},
      {{}, {947, 346500}, {}, {}, {5903, 138600}, {}, {}, {-23573, 147000}, {}, {}, {27, 20000}}
    };

    std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::_Rational<long long>>>
    Q
    {
      {{}, {}, {3, 10}},
      {{1, 70}, {}, {}, {-17, 70}},
      {{}, {-37, 3150}, {}, {}, {611, 3150}, {}, {}, {-9, 1000}},
      {{}, {}, {79, 12375}, {}, {}, {-110767, 693000}, {}, {}, {549, 28000}}
    };

    std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::_Rational<long long>>>
    R
    {
      {{1, 1}},
      {{}, {-4, 5}},
      {{}, {}, {57, 70}, {}, {}, {-9, 100}},
      {{23, 3150}, {}, {}, {-2617, 3150}, {}, {}, {699, 3500}},
      {{}, {-1159, 115500}, {}, {}, {3889, 4620}, {}, {}, {-46631, 147000}, {}, {}, {27, 20000}}
    };

    std::vector<__gnu_cxx::_Polynomial<__gnu_cxx::_Rational<long long>>>
    S
    {
      {{-1, 5}, {}, {}, {3, 5}},
      {{}, {1, 5}, {}, {}, {-131, 140}},
      {{}, {}, {-593, 3150}, {}, {}, {5437, 4500}, {}, {}, {-9, 500}},
      {{947, 346500}, {}, {}, {31727, 173250}, {}, {}, {-999443, 693000}, {}, {}, {369, 7000}}
    };

    std::cout << "\nP polynomial\n";
    for (const auto& p : P)
      std::cout << p << '\n';
    std::cout << "\nP values\n";
    for (int i = -100; i <= +100; ++i)
      {
	auto z = _Tp{0.01Q} * i;
	std::cout << std::setw(width) << z;
	for (const auto& p : P)
	  std::cout << std::setw(width) << p(z);
	std::cout << '\n';
      }

    std::cout << "\nQ polynomial\n";
    for (const auto& q : Q)
      std::cout << q << '\n';
    std::cout << "\nQ values\n";
    for (int i = -100; i <= +100; ++i)
      {
	auto z = _Tp{0.01Q} * i;
	std::cout << std::setw(width) << z;
	for (const auto& q : Q)
	  std::cout << std::setw(width) << q(z);
	std::cout << '\n';
      }

    std::cout << "\nR polynomial\n";
    for (const auto& r : R)
      std::cout << r << '\n';
    std::cout << "\nR values\n";
    for (int i = -100; i <= +100; ++i)
      {
	auto z = _Tp{0.01Q} * i;
	std::cout << std::setw(width) << z;
	for (const auto& r : R)
	  std::cout << std::setw(width) << r(z);
	std::cout << '\n';
      }

    std::cout << "\nS polynomial\n";
    for (const auto& s : S)
      std::cout << s << '\n';
    std::cout << "\nS values\n";
    for (int i = -100; i <= +100; ++i)
      {
	auto z = _Tp{0.01Q} * i;
	std::cout << std::setw(width) << z;
	for (const auto& s : S)
	  std::cout << std::setw(width) << s(z);
	std::cout << '\n';
      }

    auto __nu = _Tp{20};
    const auto _S_2p13 = _Tp{1.259921049894873164767210607278228350570Q};
    const auto _S_2p23 = _Tp{1.587401051968199474751705639272308260393Q};
    const auto _S_2p43 = _Tp{2.519842099789746329534421214556456701140Q};
    const auto _S_2p53 = _Tp{3.174802103936398949503411278544616520785Q};
    const auto _S_pi   = _Tp{3.141592653589793238462643383279502884195Q};
    const auto __nu13 = std::pow(__nu, _Tp{1} / _Tp{3});
    const auto __nu23 = __nu13 * __nu13;
    const auto __nu43 = __nu23 * __nu23;

    std::cout << "\n\nTransition region Bessel functions: J_\\nu(\\nu + a\\nu^{1/3})\n";
    std::cout << "\nnu = " << __nu << "\n"
	      << std::setw(width) << "a"
	      << std::setw(width) << "J_\\nu"
	      << std::setw(width) << "N_\\nu"
	      << std::setw(width) << "J'_\\nu"
	      << std::setw(width) << "N'_\\nu"
	      << '\n';
    for (int __i = -100; __i <= +100; ++__i)
      {
	auto __a = _Tp{0.005Q} * __i;
	const auto __airy_arg = -_S_2p13 * __a;
	_Tp _Ai, _Bi, _Aip, _Bip;
	std::__detail::__airy(__airy_arg, _Ai, _Bi, _Aip, _Bip);

	auto __num2k3 = _Tp{1};

	auto _Jsum1 = _Tp{0};
	auto _Nsum1 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __p : P)
	  {
	    _Jsum1 += __p(__a) * __num2k3;
	    _Nsum1 += __p(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	auto _Jsum2 = _Tp{0};
	auto _Nsum2 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __q : Q)
	  {
	    _Jsum2 += __q(__a) * __num2k3;
	    _Nsum2 += __q(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	const auto _Jt = _S_2p13 * _Ai * _Jsum1 / __nu13
		       + _S_2p23 * _Aip * _Jsum2 / __nu;
	const auto _Nt = -_S_2p13 * _Bi * _Nsum1 / __nu13
			- _S_2p23 * _Bip * _Nsum2 / __nu;

	auto _Jpsum1 = _Tp{0};
	auto _Npsum1 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __r : R)
	  {
	    _Jpsum1 += __r(__a) * __num2k3;
	    _Npsum1 += __r(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	auto _Jpsum2 = _Tp{0};
	auto _Npsum2 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __s : S)
	  {
	    _Jpsum2 += __s(__a) * __num2k3;
	    _Npsum2 += __s(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	const auto _Jtp = -_S_2p23 * _Aip * _Jpsum1 / __nu23
			 + _S_2p13 * _Ai * _Jpsum2 / __nu43;
	const auto _Ntp = _S_2p23 * _Bip * _Npsum1 / __nu23
			- _S_2p13 * _Bi * _Npsum2 / __nu43;

	std::cout << std::setw(width) << __a
		  << std::setw(width) << _Jt
		  << std::setw(width) << _Nt
		  << std::setw(width) << _Jtp
		  << std::setw(width) << _Ntp
		  << std::setw(width) << '\n';
      }

    const auto __mipi3 = std::polar(-_S_pi / _Tp{3});
    const auto __pipi3 = std::polar(+_S_pi / _Tp{3});
    std::cout << "\n\nTransition region Bessel functions: J_\\nu(\\nu + a\\nu^{1/3})\n";
    std::cout << "\nnu = " << __nu << "\n"
	      << std::setw(2*width) << "a"
	      << std::setw(2*width) << "J_\\nu"
	      << std::setw(2*width) << "N_\\nu"
	      << std::setw(2*width) << "J'_\\nu"
	      << std::setw(2*width) << "N'_\\nu"
	      << '\n';

    const auto __eps = std::numeric_limits<_Tp>::epsilon();
    for (int __i = -100; __i <= +100; ++__i)
      {
	auto __a = _Tp{0.005Q} * __i;

	const std::complex<_Tp> __airy_argm = _S_2p13 * __a * __mipi3;
	std::complex<_Tp> _Ami, _Bmi, _Amip, _Bmip;
	std::__detail::__airy(__airy_argm, __eps, _Ami, _Bmi, _Amip, _Bmip);

	const std::complex<_Tp> __airy_argp = _S_2p13 * __a * __pipi3;
	std::complex<_Tp> _Api, _Bpi, _Apip, _Bpip;
	std::__detail::__airy(__airy_argp, __eps, _Api, _Bpi, _Apip, _Bpip);

	auto __num2k3 = _Tp{1};

	auto _H1sum1 = _Tp{0};
	auto _H2sum1 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __p : P)
	  {
	    _H1sum1 += __p(__a) * __num2k3;
	    _H2sum1 += __p(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	auto _H1sum2 = _Tp{0};
	auto _H2sum2 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __q : Q)
	  {
	    _H1sum2 += __q(__a) * __num2k3;
	    _H2sum2 += __q(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	const auto _H1t = _S_2p43 * __mipi3 * _Ami * _H1sum1 / __nu13
			+ _S_2p53 * __mipi3 * _Amip * _H1sum2 / __nu;
	const auto _H2t = _S_2p43 * __pipi3 * _Api * _H2sum1 / __nu13
			+ _S_2p53 * __pipi3 * _Apip * _H2sum2 / __nu;

	auto _H1psum1 = _Tp{0};
	auto _H2psum1 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __r : R)
	  {
	    _H1psum1 += __r(__a) * __num2k3;
	    _H2psum1 += __r(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	auto _H1psum2 = _Tp{0};
	auto _H2psum2 = _Tp{0};
	__num2k3 = _Tp{1};
	for (const auto& __s : S)
	  {
	    _H1psum2 += __s(__a) * __num2k3;
	    _H2psum2 += __s(__a) * __num2k3;
	    __num2k3 /= __nu23;
	  }

	const auto _H1pt = -_S_2p53 * __mipi3 * _Amip * _H1psum1 / __nu23
			  + _S_2p43 * __mipi3 * _Ami * _H1psum2 / __nu43;
	const auto _H2pt = -_S_2p53 * __pipi3 * _Apip * _H2psum1 / __nu23
			  + _S_2p43 * __pipi3 * _Api * _H2psum2 / __nu43;

	std::cout << std::setw(2*width) << __a
		  << std::setw(2*width) << _H1t
		  << std::setw(2*width) << _H2t
		  << std::setw(2*width) << _H1pt
		  << std::setw(2*width) << _H2pt
		  << std::setw(2*width) << '\n';
      }
  }

int
main()
{
  std::cout << "\n\nfloat\n";
  std::cout << "=====\n";
  hankel_transition<float>();

  std::cout << "\n\ndouble\n";
  std::cout << "======\n";
  hankel_transition<double>();

  std::cout << "\n\nlong double\n";
  std::cout << "===========\n";
  hankel_transition<long double>();
}
