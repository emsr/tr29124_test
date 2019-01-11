/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_marcum_q test_marcum_q.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH ./test_marcum_q > test_marcum_q.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -I. -o test_marcum_q test_marcum_q.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=$HOME/bin/lib64:wrappers/debug:$LD_LIBRARY_PATH ./test_marcum_q > test_marcum_q.txt
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
  marcum_q_bessel_series(unsigned int __m, _Tp __a, _Tp __b)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    if (__b == _Tp{0})
      return std::make_pair(_Tp{0}, _Tp{1});

    auto __q1 = _Tp{0};
    const auto __arg = __a * __b;
    const auto __ab = __a / __b;
    auto __temp = _Tp{1};
    unsigned int __k = 0;
    while (true)
      {
	const auto __term = __temp * std::cyl_bessel_i(__k, __arg);
	__q1 += __term;
	if (std::abs(__term) < _S_eps)
          break;
	__temp *= __ab;
	++__k;
      }

    auto __q = __q1;

    const _Tp __ba = __b / __a;
    __temp = _Tp{1};
    for (unsigned int __k = 1; __k < __m; ++__k)
      {
	__temp *= __ba;
	__q += __temp * std::cyl_bessel_i(__k, __arg);
      }

    __q *= std::exp(-(__a * __a + __b * __b) / _Tp{2});

    return std::make_pair(_Tp{1} - __q, __q);
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
  marcum_q_gamma_series(_Tp __mu, _Tp __a, _Tp __b)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    unsigned int _S_max_iter = 100u;
    __a *= __a / _Tp{2};
    __b *= __b / _Tp{2};
    auto __fact = _Tp{1};
    auto [_MP, _MQ] = std::__detail::__gamma(__mu, __b);
    for (unsigned int __k = 1; __k < _S_max_iter; ++__k)
      {
	const auto [_Pgam, _Qgam] = std::__detail::__gamma(__k + __mu, __b);
	__fact *= __a / __k;
	const auto __termP = __fact * _Pgam;
	_MP += __termP;
	const auto __termQ = __fact * _Qgam;
	_MQ += __termQ;
	if (std::abs(__termP) < _S_eps * _MP
	 && std::abs(__termQ) < _S_eps * _MQ)
	  break;
      }
    auto __exp = std::exp(-__a);
    return std::make_pair(__exp * _MP, __exp * _MQ);
  }

//
// @brief The Marcum Q function.
//
// Recent software developments for special functions
// in the Santander-Amsterdam project
//
template<typename _Tp>
  _Tp
  marcum_q_integral(unsigned int __m, _Tp __a, _Tp __b)
  {
    //auto __rho = [](_Tp __theta, _Tp __xi)
	//	 -> _Tp
	//	 { return ; };

    //return __q;
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
