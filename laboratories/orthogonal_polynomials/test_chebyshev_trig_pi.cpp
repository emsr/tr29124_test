/**
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <emsr/float128_io.h>

/**
 * @todo return derivatives from these.
 */

  /**
   * Return a reperiodized arccosine:
   * @f[
   *   acos_pi(x) = \frac{1}{\pi}acos(x)
   * @f]
   */
  template<typename _Tp>
    inline _Tp
    acos_pi(_Tp x)
    { return std::acos(x) / emsr::pi_v<_Tp>; }

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
    chebyshev_t_trig(unsigned int n, _Tp x)
    {
      auto theta_pi = acos_pi(x);
      return emsr::detail::cos_pi(n * theta_pi);
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
    chebyshev_t_zeros(unsigned int n)
    {
      std::vector<_Tp> zero;
      for (unsigned int k = 0; k < n; ++k)
	zero.push_back(emsr::cos_pi(_Tp(k + 0.5L) / _Tp(n)));
      return zero;
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
    chebyshev_u_trig(unsigned int n, _Tp x)
    {
      const auto _S_eps = emsr::epsilon(x);
      if (std::abs(x + _Tp{1}) < _S_eps)
	return (n % 2 == 0 ? +1 : -1) * _Tp(n + 1);
      else if (std::abs(x - _Tp{1}) < _S_eps)
	return _Tp(n + 1);
      else
	{
	  const auto theta_pi = acos_pi(x);
	  return emsr::detail::sin_pi(_Tp(n + 1) * theta_pi)
	       / emsr::detail::sin_pi(theta_pi);
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
    chebyshev_u_zeros(unsigned int n)
    {
      std::vector<_Tp> zero;
      for (unsigned int k = 1; k <= n; ++k)
	zero.push_back(emsr::cos_pi(_Tp(k) / _Tp(n + 1)));
      return zero;
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
    chebyshev_v_trig(unsigned int n, _Tp x)
    {
      const auto _S_eps = emsr::epsilon(x);
      if (std::abs(x + _Tp{1}) < _S_eps)
	return (n % 2 == 0 ? +1 : -1) * _Tp(2 * n + 1);
      else
	{
	  const auto theta_pi = acos_pi(x);
	  return emsr::detail::cos_pi(_Tp(n + 0.5L) * theta_pi)
	       / emsr::detail::cos_pi(_Tp{0.5L} * theta_pi);
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
    chebyshev_v_zeros(unsigned int n)
    {
      std::vector<_Tp> zero;
      for (unsigned int k = 0; k < n; ++k)
	zero.push_back(emsr::cos_pi(_Tp(k + 0.5L) / _Tp(n + 0.5L)));
      return zero;
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
    chebyshev_w_trig(unsigned int n, _Tp x)
    {
      const auto _S_eps = emsr::epsilon(x);
      if (std::abs(x - _Tp{1}) < _S_eps)
	return _Tp(2 * n + 1);
      else
	{
	  const auto theta_pi = acos_pi(x);
	  return emsr::detail::sin_pi(_Tp(n + 0.5L) * theta_pi)
	       / emsr::detail::sin_pi(_Tp{0.5L} * theta_pi);
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
    chebyshev_w_zeros(unsigned int n)
    {
      std::vector<_Tp> zero;
      for (unsigned int k = 1; k <= n; ++k)
	zero.push_back(emsr::cos_pi(_Tp(k) / _Tp(n + 0.5L)));
      return zero;
    }

template<typename Tp>
  void
  test_chebyshev(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<unsigned int> index{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    for (auto n : index)
      {
	std::cout << "\n n = " << std::setw(width) << n << '\n';
	std::cout << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "Tt"
		  << ' ' << std::setw(width) << "Tg"
		  << ' ' << std::setw(width) << "Ut"
		  << ' ' << std::setw(width) << "Ug"
		  << ' ' << std::setw(width) << "Vt"
		  << ' ' << std::setw(width) << "Vg"
		  << ' ' << std::setw(width) << "Wt"
		  << ' ' << std::setw(width) << "Wg"
		  << ' ' << std::setw(width) << "Tt - Tg"
		  << ' ' << std::setw(width) << "Ut - Ug"
		  << ' ' << std::setw(width) << "Vt - Vg"
		  << ' ' << std::setw(width) << "Wt - Wg"
		  << '\n';
	const auto del = Tp{1} / Tp{100};
	for (int i = -100; i <= 100; ++i)
	  {
	    auto x = del * i;
	    auto Tt = chebyshev_t_trig(n, x);
	    auto Ut = chebyshev_u_trig(n, x);
	    auto Vt = chebyshev_v_trig(n, x);
	    auto Wt = chebyshev_w_trig(n, x);
	    auto Tg = emsr::chebyshev_t(n, x);
	    auto Ug = emsr::chebyshev_u(n, x);
	    auto Vg = emsr::chebyshev_v(n, x);
	    auto Wg = emsr::chebyshev_w(n, x);
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
	auto tz = chebyshev_t_zeros<Tp>(n);
	for (auto z : tz)
	  {
	    auto Tt = chebyshev_t_trig(n, z);
	    auto Tg = emsr::chebyshev_t(n, z);
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
	auto uz = chebyshev_u_zeros<Tp>(n);
	for (auto z : uz)
	  {
	    auto Ut = chebyshev_u_trig(n, z);
	    auto Ug = emsr::chebyshev_u(n, z);
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
	auto vz = chebyshev_v_zeros<Tp>(n);
	for (auto z : vz)
	  {
	    auto Vt = chebyshev_v_trig(n, z);
	    auto Vg = emsr::chebyshev_v(n, z);
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
	auto wz = chebyshev_w_zeros<Tp>(n);
	for (auto z : wz)
	  {
	    auto Wt = chebyshev_w_trig(n, z);
	    auto Wg = emsr::chebyshev_w(n, z);
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
