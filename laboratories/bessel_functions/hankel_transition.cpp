/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>

#include <emsr/special_functions.h>
#include <emsr/rational_polynomial.h>
#include <emsr/polynomial.h>

template<typename _Tp>
  void
  hankel_transition()
  {
    auto sign = [](int s, int r){return (s + r) % 2 == 1 ? -1 : +1; };
    using rational = _Tp;

    const std::size_t n_AB = 8;
    const std::size_t n_phipsi = 3 * (n_AB - 1) + 1;

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = std::cout.precision() + 6;

    emsr::Polynomial<int> poo({0, -2});
    std::vector<emsr::Polynomial<int>> phi;
    std::vector<emsr::Polynomial<int>> psi;
    phi.push_back(emsr::Polynomial<int>(1));
    phi.push_back(emsr::Polynomial<int>(0));
    phi.push_back(emsr::Polynomial<int>({0, -2}));
    psi.push_back(emsr::Polynomial<int>(0));
    psi.push_back(emsr::Polynomial<int>(-1));
    psi.push_back(emsr::Polynomial<int>(0));
    for (std::size_t s = 3; s < n_phipsi; ++s)
      {
	phi.push_back(poo * phi[s - 2] - int(2 * (s - 2)) * phi[s - 3]);
	psi.push_back(poo * psi[s - 2] - int(2 * (s - 2)) * psi[s - 3]);
      }
    std::cout << "\nphi polynomial\n";
    for (const auto& p : phi)
      std::cout << p << '\n';
    std::cout << "\npsi polynomial\n";
    for (const auto& p : psi)
      std::cout << p << '\n';
    std::vector<emsr::Polynomial<rational>> A(n_AB);
    std::vector<emsr::Polynomial<rational>> B(n_AB);
    for (std::size_t s = 0; s < n_AB; ++s)
      {
	for (std::size_t r = 0; r <= s; ++r)
	  {
	    A[s] += phi[2 * s + r] * emsr::Polynomial<rational>(sign(r, s), r);
	    B[s] += psi[2 * s + r] * emsr::Polynomial<rational>(sign(r, s), r);
	  }
      }
    std::cout << "\nA polynomial\n";
    for (const auto& a : A)
      std::cout << a << '\n';
    std::cout << "\nB polynomial\n";
    for (const auto& b : B)
      std::cout << b << '\n';

    using _RPoly = emsr::RationalPolynomial<long long>;

    std::vector<emsr::Polynomial<_RPoly>>
    P
    {
      {{1, 1}},
      {{}, {-1, 5}},
      {{}, {}, {3, 35}, {}, {}, {-9, 100}},
      {{-1, 225}, {}, {}, {-173, 3150}, {}, {}, {957, 7000}},
      {{}, {947, 346500}, {}, {}, {5903, 138600}, {}, {}, {-23573, 147000}, {}, {}, {27, 20000}}
    };

    std::vector<emsr::Polynomial<_RPoly>>
    Q
    {
      {{}, {}, {3, 10}},
      {{1, 70}, {}, {}, {-17, 70}},
      {{}, {-37, 3150}, {}, {}, {611, 3150}, {}, {}, {-9, 1000}},
      {{}, {}, {79, 12375}, {}, {}, {-110767, 693000}, {}, {}, {549, 28000}}
    };

    std::vector<emsr::Polynomial<_RPoly>>
    R
    {
      {{1, 1}},
      {{}, {-4, 5}},
      {{}, {}, {57, 70}, {}, {}, {-9, 100}},
      {{23, 3150}, {}, {}, {-2617, 3150}, {}, {}, {699, 3500}},
      {{}, {-1159, 115500}, {}, {}, {3889, 4620}, {}, {}, {-46631, 147000}, {}, {}, {27, 20000}}
    };

    std::vector<emsr::Polynomial<_RPoly>>
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

    auto nu = _Tp{20};
    const auto s_2p13 = _Tp{1.259921049894873164767210607278228350570Q};
    const auto s_2p23 = _Tp{1.587401051968199474751705639272308260393Q};
    const auto s_2p43 = _Tp{2.519842099789746329534421214556456701140Q};
    const auto s_2p53 = _Tp{3.174802103936398949503411278544616520785Q};
    const auto s_pi   = _Tp{3.141592653589793238462643383279502884195Q};
    const auto nu13 = std::pow(nu, _Tp{1} / _Tp{3});
    const auto nu23 = nu13 * nu13;
    const auto nu43 = nu23 * nu23;

    std::cout << "\n\nTransition region Bessel functions: J_\\nu(\\nu + a\\nu^{1/3})\n";
    std::cout << "\nnu = " << nu << "\n"
	      << std::setw(width) << "a"
	      << std::setw(width) << "J_\\nu"
	      << std::setw(width) << "N_\\nu"
	      << std::setw(width) << "J'_\\nu"
	      << std::setw(width) << "N'_\\nu"
	      << '\n';
    for (int i = -100; i <= +100; ++i)
      {
	auto a = _Tp{0.005Q} * i;
	const auto airy_arg = -s_2p13 * a;
	auto _Airy = emsr::detail::airy(airy_arg);

	auto num2k3 = _Tp{1};

	auto _Jsum1 = _Tp{0};
	auto _Nsum1 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& p : P)
	  {
	    _Jsum1 += p(a) * num2k3;
	    _Nsum1 += p(a) * num2k3;
	    num2k3 /= nu23;
	  }

	auto _Jsum2 = _Tp{0};
	auto _Nsum2 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& q : Q)
	  {
	    _Jsum2 += q(a) * num2k3;
	    _Nsum2 += q(a) * num2k3;
	    num2k3 /= nu23;
	  }

	const auto _Jt = s_2p13 * _Airy.Ai_value * _Jsum1 / nu13
		       + s_2p23 * _Airy.Ai_deriv * _Jsum2 / nu;
	const auto _Nt = -s_2p13 * _Airy.Bi_value * _Nsum1 / nu13
			- s_2p23 * _Airy.Bi_deriv * _Nsum2 / nu;

	auto _Jpsum1 = _Tp{0};
	auto _Npsum1 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& r : R)
	  {
	    _Jpsum1 += r(a) * num2k3;
	    _Npsum1 += r(a) * num2k3;
	    num2k3 /= nu23;
	  }

	auto _Jpsum2 = _Tp{0};
	auto _Npsum2 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& s : S)
	  {
	    _Jpsum2 += s(a) * num2k3;
	    _Npsum2 += s(a) * num2k3;
	    num2k3 /= nu23;
	  }

	const auto _Jtp = -s_2p23 * _Airy.Ai_deriv * _Jpsum1 / nu23
			 + s_2p13 * _Airy.Ai_value * _Jpsum2 / nu43;
	const auto _Ntp = s_2p23 * _Airy.Bi_deriv * _Npsum1 / nu23
			- s_2p13 * _Airy.Bi_value * _Npsum2 / nu43;

	std::cout << std::setw(width) << a
		  << std::setw(width) << _Jt
		  << std::setw(width) << _Nt
		  << std::setw(width) << _Jtp
		  << std::setw(width) << _Ntp
		  << std::setw(width) << '\n';
      }

    const auto mipi3 = std::polar(-s_pi / _Tp{3});
    const auto pipi3 = std::polar(+s_pi / _Tp{3});
    std::cout << "\n\nTransition region Bessel functions: J_\\nu(\\nu + a\\nu^{1/3})\n";
    std::cout << "\nnu = " << nu << "\n"
	      << std::setw(2*width) << "a"
	      << std::setw(2*width) << "J_\\nu"
	      << std::setw(2*width) << "N_\\nu"
	      << std::setw(2*width) << "J'_\\nu"
	      << std::setw(2*width) << "N'_\\nu"
	      << '\n';

    //const auto eps = std::numeric_limits<_Tp>::epsilon();
    for (int i = -100; i <= +100; ++i)
      {
	auto a = _Tp{0.005Q} * i;

	const std::complex<_Tp> airy_argm = s_2p13 * a * mipi3;
	auto airym = emsr::detail::_Airy<std::complex<_Tp>>()(airy_argm);

	const std::complex<_Tp> airy_argp = s_2p13 * a * pipi3;
	auto airyp = emsr::detail::_Airy<std::complex<_Tp>>()(airy_argp);

	auto num2k3 = _Tp{1};

	auto _H1sum1 = _Tp{0};
	auto _H2sum1 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& p : P)
	  {
	    _H1sum1 += p(a) * num2k3;
	    _H2sum1 += p(a) * num2k3;
	    num2k3 /= nu23;
	  }

	auto _H1sum2 = _Tp{0};
	auto _H2sum2 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& q : Q)
	  {
	    _H1sum2 += q(a) * num2k3;
	    _H2sum2 += q(a) * num2k3;
	    num2k3 /= nu23;
	  }

	const auto _H1t = s_2p43 * mipi3 * airym.Ai_value * _H1sum1 / nu13
			+ s_2p53 * mipi3 * airym.Ai_deriv * _H1sum2 / nu;
	const auto _H2t = s_2p43 * pipi3 * airyp.Ai_value * _H2sum1 / nu13
			+ s_2p53 * pipi3 * airyp.Ai_deriv * _H2sum2 / nu;

	auto _H1psum1 = _Tp{0};
	auto _H2psum1 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& r : R)
	  {
	    _H1psum1 += r(a) * num2k3;
	    _H2psum1 += r(a) * num2k3;
	    num2k3 /= nu23;
	  }

	auto _H1psum2 = _Tp{0};
	auto _H2psum2 = _Tp{0};
	num2k3 = _Tp{1};
	for (const auto& s : S)
	  {
	    _H1psum2 += s(a) * num2k3;
	    _H2psum2 += s(a) * num2k3;
	    num2k3 /= nu23;
	  }

	const auto _H1pt = -s_2p53 * mipi3 * airym.Ai_deriv * _H1psum1 / nu23
			  + s_2p43 * mipi3 * airym.Ai_value * _H1psum2 / nu43;
	const auto _H2pt = -s_2p53 * pipi3 * airyp.Ai_deriv * _H2psum1 / nu23
			  + s_2p43 * pipi3 * airyp.Ai_value * _H2psum2 / nu43;

	std::cout << std::setw(2*width) << a
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
