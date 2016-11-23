/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_chebyshev_trig test_chebyshev_trig.cpp -lquadmath
./test_chebyshev_trig > test_chebyshev_trig.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_chebyshev_trig test_chebyshev_trig.cpp -lquadmath
./test_chebyshev_trig > test_chebyshev_trig.txt

g++ -std=gnu++14 -o test_chebyshev_trig test_chebyshev_trig.cpp -lquadmath
./test_chebyshev_trig > test_chebyshev_trig.txt
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <bits/float128_io.h>

/**
 * Return the Chebyshev polynomial of the first kind
 * of order @f$ n @f$, @f$ T_n(x) @f$ by trigonometric identity:
 * @f[
 *   T_n(x) = \cos(n\theta)
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_t_trig(unsigned int __n, _Tp __x)
  {
    auto __theta = std::acos(__x);
    return std::cos(__n * __theta);
  }

/**
 * Return a vector of zeros of the Chebyshev function of the first kind
 * of order @f$ n @f$, @f$ T_n(x) @f$.
 * The zeros are given by:
 * @f[
 *   x_k = \cos\left(\frac{(k+\frac{1}{2})\pi}{n+1}\right),
 *     k \elem {1, ..., n}
 * @f]
 */
template<typename _Tp>
  std::vector<_Tp>
  __chebyshev_t_zeros(unsigned int __n)
  {
    std::vector<_Tp> __zero;
    for (unsigned int __k = 0; __k < __n; ++__k)
      __zero.push_back(__gnu_cxx::cos_pi(_Tp(__k + 0.5L) / _Tp(__n)));
    return __zero;
  }

/**
 * Return the Chebyshev polynomial of the second kind
 * of order @f$ n @f$, @f$ U_n(x) @f$ by trigonometric identity:
 * @f[
 *   U_n(x) = \frac{\sin((n + 1)\theta)}{\sin(\theta)}
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_u_trig(unsigned int __n, _Tp __x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    if (std::abs(__x + _Tp{1}) < _S_eps)
      return (__n % 2 == 0 ? +1 : -1) * _Tp(__n + 1);
    else if (std::abs(__x - _Tp{1}) < _S_eps)
      return _Tp(__n + 1);
    else
      {
	auto __theta = std::acos(__x);
	return std::sin(_Tp(__n + 1) * __theta)
	     / std::sin(__theta);
      }
  }

/**
 * Return a vector of zeros of the Chebyshev function of the second kind
 * of order @f$ n @f$, @f$ U_n(x) @f$.
 * The zeros are given by:
 * @f[
 *   x_k = \cos\left(\frac{k\pi}{n+1}\right), k \elem {1, ..., n}
 * @f]
 */
template<typename _Tp>
  std::vector<_Tp>
  __chebyshev_u_zeros(unsigned int __n)
  {
    std::vector<_Tp> __zero;
    for (unsigned int __k = 1; __k <= __n; ++__k)
      __zero.push_back(__gnu_cxx::cos_pi(_Tp(__k) / _Tp(__n + 1)));
    return __zero;
  }

/**
 * Return the Chebyshev polynomial of the third kind
 * of order @f$ n @f$, @f$ V_n(x) @f$ by trigonometric identity:
 * @f[
 *   V_n(x) = \frac{\cos((n + \frac{1}{2})\theta)}
 *                 {cos(\frac{1}{2}\theta)}
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_v_trig(unsigned int __n, _Tp __x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    if (std::abs(__x + _Tp{1}) < _S_eps)
      return (__n % 2 == 0 ? +1 : -1) * _Tp(2 * __n + 1);
    else
      {
	auto __theta = std::acos(__x);
	return std::cos(_Tp(__n + 0.5L) * __theta)
	     / std::cos(_Tp{0.5L} * __theta);
      }
  }

/**
 * Return a vector of zeros of the Chebyshev function of the third kind
 * of order @f$ n @f$, @f$ V_n(x) @f$.
 * The zeros are given by:
 * @f[
 *   x_k = \cos\left(\frac{\left(k+\frac{1}{2}\right)\pi}{n+\frac{1}{2}}\right),
 *       k \elem {1, ..., n}
 * @f]
 */
template<typename _Tp>
  std::vector<_Tp>
  __chebyshev_v_zeros(unsigned int __n)
  {
    std::vector<_Tp> __zero;
    for (unsigned int __k = 0; __k < __n; ++__k)
      __zero.push_back(__gnu_cxx::cos_pi(_Tp(__k + 0.5L) / _Tp(__n + 0.5L)));
    return __zero;
  }

/**
 * Return the Chebyshev polynomial of the fourth kind
 * of order @f$ n @f$, @f$ W_n(x) @f$ by trigonometric identity:
 * @f[
 *   W_n(x) = \frac{\sin((n + \frac{1}{2})\theta)}
 *                 {\sin(\frac{1}{2}\theta)}
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_w_trig(unsigned int __n, _Tp __x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    if (std::abs(__x - _Tp{1}) < _S_eps)
      return _Tp(2 * __n + 1);
    else
      {
	auto __theta = std::acos(__x);
	return std::sin(_Tp(__n + 0.5L) * __theta)
	     / std::sin(_Tp{0.5L} * __theta);
      }
  }

/**
 * Return a vector of zeros of the Chebyshev function of the fourth kind
 * of order @f$ n @f$, @f$ W_n(x) @f$.
 * The zeros are given by:
 * @f[
 *   x_k = \cos\left(\frac{k\pi}{n+\frac{1}{2}}\right),
 *       k \elem {1, ..., n}
 * @f]
 */
template<typename _Tp>
  std::vector<_Tp>
  __chebyshev_w_zeros(unsigned int __n)
  {
    std::vector<_Tp> __zero;
    for (unsigned int __k = 1; __k <= __n; ++__k)
      __zero.push_back(__gnu_cxx::cos_pi(_Tp(__k) / _Tp(__n + 0.5L)));
    return __zero;
  }

template<typename Tp>
  void
  test_chebyshev(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<unsigned int> index{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    for (auto n : index)
      {
	std::cout << "\n n = " << std::setw(width) << n << '\n';
	for (int i = -100; i <= 100; ++i)
	  {
	    auto x = Tp{0.01Q} * i;
	    auto Tt = __chebyshev_t_trig(n, x);
	    auto Ut = __chebyshev_u_trig(n, x);
	    auto Vt = __chebyshev_v_trig(n, x);
	    auto Wt = __chebyshev_w_trig(n, x);
	    auto Tg = __gnu_cxx::chebyshev_t(n, x);
	    auto Ug = __gnu_cxx::chebyshev_u(n, x);
	    auto Vg = __gnu_cxx::chebyshev_v(n, x);
	    auto Wg = __gnu_cxx::chebyshev_w(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << Tt
		      << ' ' << std::setw(width) << Tg
		      << ' ' << std::setw(width) << Ut
		      << ' ' << std::setw(width) << Ug
		      << ' ' << std::setw(width) << Vt
		      << ' ' << std::setw(width) << Vg
		      << ' ' << std::setw(width) << Wt
		      << ' ' << std::setw(width) << Wg
		      << ' ' << std::setw(width) << Tt - Tg
		      << ' ' << std::setw(width) << Ut - Ug
		      << ' ' << std::setw(width) << Vt - Vg
		      << ' ' << std::setw(width) << Wt - Wg
		      << '\n';
	  }
      }

    std::cout << "\n\nZeros of Chebyshev T_n(x)\n";
    for (auto n : index)
      {
	std::cout << "\n n = " << std::setw(width) << n << '\n';
	auto tz = __chebyshev_t_zeros<Tp>(n);
	for (auto z : tz)
	  {
	    auto Tt = __chebyshev_t_trig(n, z);
	    auto Tg = __gnu_cxx::chebyshev_t(n, z);
	    std::cout << ' ' << std::setw(width) << z
		      << ' ' << std::setw(width) << Tt
		      << ' ' << std::setw(width) << Tg
		      << '\n';
	  }
      }

    std::cout << "\n\nZeros of Chebyshev U_n(x)\n";
    for (auto n : index)
      {
	std::cout << "\n n = " << std::setw(width) << n << '\n';
	auto uz = __chebyshev_u_zeros<Tp>(n);
	for (auto z : uz)
	  {
	    auto Ut = __chebyshev_u_trig(n, z);
	    auto Ug = __gnu_cxx::chebyshev_u(n, z);
	    std::cout << ' ' << std::setw(width) << z
		      << ' ' << std::setw(width) << Ut
		      << ' ' << std::setw(width) << Ug
		      << '\n';
	  }
      }

    std::cout << "\n\nZeros of Chebyshev V_n(x)\n";
    for (auto n : index)
      {
	std::cout << "\n n = " << std::setw(width) << n << '\n';
	auto vz = __chebyshev_v_zeros<Tp>(n);
	for (auto z : vz)
	  {
	    auto Vt = __chebyshev_v_trig(n, z);
	    auto Vg = __gnu_cxx::chebyshev_v(n, z);
	    std::cout << ' ' << std::setw(width) << z
		      << ' ' << std::setw(width) << Vt
		      << ' ' << std::setw(width) << Vg
		      << '\n';
	  }
      }

    std::cout << "\n\nZeros of Chebyshev W_n(x)\n";
    for (auto n : index)
      {
	std::cout << "\n n = " << std::setw(width) << n << '\n';
	auto wz = __chebyshev_w_zeros<Tp>(n);
	for (auto z : wz)
	  {
	    auto Wt = __chebyshev_w_trig(n, z);
	    auto Wg = __gnu_cxx::chebyshev_w(n, z);
	    std::cout << ' ' << std::setw(width) << z
		      << ' ' << std::setw(width) << Wt
		      << ' ' << std::setw(width) << Wg
		      << '\n';
	  }
      }
  }
int
main()
{
  test_chebyshev(1.0);
}
