/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_zeta_trig test_zeta_trig.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_zeta_trig > test_zeta_trig.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_zeta_trig test_zeta_trig.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
PATH=wrappers/debug:$PATH ./test_zeta_trig > test_zeta_trig.txt
*/

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>

/**
 * Return the reperiodized cotangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The cotangent function is given by:
 * @f[
 *    cot(\pi x) = \frac{1}{\pi x}
 *      - \frac{2}{\pi x}\sum_{k=1}^{\infty}\zeta(2k)x^{2k}
 * @f]
 */
template<typename _Tp>
  _Tp
  cot_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = __x * __x;
    auto __xxk = __xx;
    unsigned int __k = 2;
    auto __res = _Tp{0};
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __zeta_even<_Tp>[__k] * __xx;
    return (_Tp{1} - _Tp{2} * __res) / __pi_x;
  }

/**
 * Return the reperiodized cosecant for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The cosecant function is given by:
 * @f[
 *    csc(\pi x) = \frac{1}{\pi x}
 *      + \frac{2}{\pi x}\sum_{k=1}^{\infty}\eta(2k)x^{2k}
 * @f]
 */
template<typename _Tp>
  _Tp
  csc_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = __x * __x;
    auto __xxk = __xx;
    auto __res = _Tp{0};
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __eta_even<_Tp>[__k] * __xx;
    return (_Tp{1} + _Tp{2} * __res) / __pi_x;
  }

/**
 * Return the reperiodized tangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The tangent function is given by:
 * @f[
 *    tan\left(\frac{\pi x}{2}\right) =
 *      -\frac{4}{\pi x}\sum_{k=1}^{\infty}(-1)^k\lambda(2k)x^{2k}
 * @f]
 */
template<typename _Tp>
  _Tp
  tan_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = __x * __x;
    auto __xxk = __xx;
    auto __res = _Tp{0};
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __lambda_even<_Tp>[__k] * __xx;
    return _Tp{4} * __res / __pi_x;
  }

/**
 * Return the reperiodized secant for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The secant function is given by:
 * @f[
 *    sec\left(\frac{\pi x}{2}\right) =
 *      \frac{4}{\pi x}\sum_{k=1}^{\infty}\beta(2k-1)x^{2k-1}
 * @f]
 */
template<typename _Tp>
  _Tp
  sec_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = __x * __x;
    auto __xxk = __x;
    unsigned int __k = 1;
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __beta_odd<_Tp>[__k] * __xx;
    return _Tp{4} * __res / __pi_x;
  }

/**
 * Return the reperiodized hyperbolic cotangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The hyperbolic cotangent function is given by:
 * @f[
 *    coth(\pi x) = \frac{1}{\pi x}
 *      - \frac{2}{\pi x}\sum_{k=1}^{\infty}(-1)^k\zeta(2k)x^{2k}
 * @f]
 */
template<typename _Tp>
  _Tp
  coth_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = -__x * __x;
    auto __xxk = __xx;
    unsigned int __k = 2;
    auto __res = _Tp{0};
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __zeta_even<_Tp>[__k] * __xx;
    return (_Tp{1} - _Tp{2} * __res) / __pi_x;
  }

/**
 * Return the reperiodized csch for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The csch function is given by:
 * @f[
 *    csch(\pi x) = \frac{1}{\pi x}
 *      + \frac{2}{\pi x}\sum_{k=1}^{\infty}(-1)^k\eta(2k)x^{2k}
 * @f]
 */
template<typename _Tp>
  _Tp
  csch_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = -__x * __x;
    auto __xxk = __xx;
    auto __res = _Tp{0};
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __eta_even<_Tp>[__k] * __xx;
    return (_Tp{1} + _Tp{2} * __res) / __pi_x;
  }

/**
 * Return the reperiodized hyperbolic tangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The hyperbolic tangent function is given by:
 * @f[
 *    tan\left(\frac{\pi x}{2}\right) =
 *      -\frac{4}{\pi x}\sum_{k=1}^{\infty}(-1)^k\lambda(2k)x^{2k}
 * @f]
 */
template<typename _Tp>
  _Tp
  tanh_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = -__x * __x;
    auto __xxk = __xx;
    auto __res = _Tp{0};
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __lambda_even<_Tp>[__k] * __xx;
    return -_Tp{4} * __res / __pi_x;
  }

/**
 * Return the reperiodized hyperbolic secant function for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The hyperbolic secant function is given by:
 * @f[
 *    sech\left(\frac{\pi x}{2}\right) =
 *      -\frac{4}{\pi x}\sum_{k=1}^{\infty}(-1)^k\beta(2k-1)x^{2k-1}
 * @f]
 */
template<typename _Tp>
  _Tp
  sech_pi_zeta(_Tp __x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(__x);
    const auto __pi_x = _S_pi * __x;
    const auto __xx = -__x * __x;
    auto __xxk = __x;
    auto __res = _Tp{0};
    for (unsigned int __k = 0; __k < _S_num; ++__k, __xxk *= __xx)
      __res += __beta_odd<_Tp>[__k] * __xx;
    return -_Tp{4} * __res / __pi_x;
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  cot_pi()
  {
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  csc_pi()
  {
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  tan_pi()
  {
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  sec_pi()
  {
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  coth_pi()
  {
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  csch_pi()
  {
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  tanh_pi()
  {
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  sech_pi()
  {
  }

int
main()
{
}
