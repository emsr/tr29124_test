/**
 *
 */

#include <cmath>
#include <limits>
#include <utility>

/**
 * Return the Marcum Q function for integer order m and a > 0, b >= 0.
 * The Marcum Q function is defined by
 * @f[
 *    Q_M(a,b) = \int_{b}^{\infty}x\left(\frac{x}{a}\right)
 *        \exp\left(-\frac{x^2 + a^2}{2}\right) I_{M-1}(ax)dx
 * @f]
 * where @f$ I_{\nu}(x) @f$ is the modified Bessel function of the first kind.
 * This routine uses the expansion in modified Bessel functions for integer M.
 * @f[
 *    Q_M(a,b) = \exp\left(-\frac{a^2 + b^2}{2}\right)
 *        \sum_{k=1-M}^{\infty}\left(\frac{a}{b}\right)^k I_k(ab)
 * @f]
 */
template<typename _Tp>
  std::pair<_Tp, _Tp>
  marcum_q_bessel_series(unsigned int m, _Tp a, _Tp b)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    if (b == _Tp{0})
      return std::make_pair(_Tp{0}, _Tp{1});

    auto q1 = _Tp{0};
    const auto arg = a * b;
    const auto ab = a / b;
    auto temp = _Tp{1};
    unsigned int k = 0;
    while (true)
      {
	const auto term = temp * emsr::cyl_bessel_i(k, arg);
	q1 += term;
	if (std::abs(term) < _S_eps)
          break;
	temp *= ab;
	++k;
      }

    auto q = q1;

    const _Tp ba = b / a;
    temp = _Tp{1};
    for (unsigned int k = 1; k < m; ++k)
      {
	temp *= ba;
	q += temp * emsr::cyl_bessel_i(k, arg);
      }

    q *= std::exp(-(a * a + b * b) / _Tp{2});

    return std::make_pair(_Tp{1} - q, q);
  }

/**
 * Return the Marcum Q function for ral order m and a > 0, b >= 0.
 * The Marcum Q function is defined by
 * @f[
 *    Q_M(a,b) = \int_{b}^{\infty}x\left(\frac{x}{a}\right)
 *        \exp\left(-\frac{x^2 + a^2}{2}\right) I_{M-1}(ax)dx
 * @f]
 * where @f$ I_{\nu}(x) @f$ is the modified Bessel function of the first kind.
 * This routine uses the expansion in normalized upper gamma functions.
 * @f[
 *    Q_M(a,b) = e^{-a^2/2} \sum_{k=0}^{\infty}
 *       \left(\frac{a^2}{2}\right)^k \frac{Q(M + k, b^2/2)}{k!}
 * @f]
 *
 * @see Recent software developments for special functions
 * in the Santander-Amsterdam project.
 */
template<typename _Tp>
  std::pair<_Tp, _Tp>
  marcum_q_gamma_series(_Tp mu, _Tp a, _Tp b)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    unsigned int _S_max_iter = 100u;
    a *= a / _Tp{2};
    b *= b / _Tp{2};
    auto fact = _Tp{1};
    auto [_MP, _MQ] = emsr::detail::gamma(mu, b);
    for (unsigned int k = 1; k < _S_max_iter; ++k)
      {
	const auto [_Pgam, _Qgam] = emsr::detail::gamma(k + mu, b);
	fact *= a / k;
	const auto termP = fact * _Pgam;
	_MP += termP;
	const auto termQ = fact * _Qgam;
	_MQ += termQ;
	if (std::abs(termP) < _S_eps * _MP
	 && std::abs(termQ) < _S_eps * _MQ)
	  break;
      }
    auto exp = std::exp(-a);
    return std::make_pair(exp * _MP, exp * _MQ);
  }

//
// @brief The Marcum Q function.
//
// Recent software developments for special functions
// in the Santander-Amsterdam project
//
template<typename _Tp>
  _Tp
  marcum_q_integral(unsigned int m, _Tp a, _Tp b)
  {
    //auto rho = [](_Tp theta, _Tp xi)
	//	 -> _Tp
	//	 { return ; };

    //return q;
  }

template<typename _Tp>
  void
  test_marcum_q()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    int m = 2;
    auto a = _Tp{1};
    //auto b = _Tp{2};
    for (int i = 0; i <= 100; ++i)
      {
	auto b = _Tp{i * 0.1};
	auto [p1, q1] = marcum_q_bessel_series(m, a, b);
	auto [p2, q2] = marcum_q_gamma_series(_Tp(m), a, b);
	std::cout << ' ' << std::setw(w) << b
		  << ' ' << std::setw(w) << q1
		  << ' ' << std::setw(w) << q2
		  << '\n';
      }
  }

int
main()
{
  test_marcum_q<double>();
}
