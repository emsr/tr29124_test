/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_inverse_gamma_distribution test_inverse_gamma_distribution.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_inverse_gamma_distribution > test_inverse_gamma_distribution.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_inverse_gamma_distribution test_inverse_gamma_distribution.cpp -lquadmath -Lwrappers/debug -lwrap_boost
PATH=wrappers/debug:$PATH ./test_inverse_gamma_distribution > test_inverse_gamma_distribution.txt
*/

#include <cmath>

  /**
   * If X has a gamma(a,b) distribution then Y=1/X has an
   * inv_gamma(a, 1/b) distribution.
   * So run x=gamma_distribution(a, 1/b) return 1/x.
   * This might be really general and so is scaling and moving center.
   */

  /**
   * The PDF of the inverse gamma distribution is:
   * @f[
   *   f(x|\alpha,\beta) = \frac{\beta^\alpha}{\Gamma(alpha)}
   *                       x^{-\alpha-1}\exp(-\beta/x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    inverse_gamma_pdf(_Tp __alpha, _Tp __beta, _Tp __x)
    {
      // a>0 b>0
      if (__x == _Tp{0})
	return _Tp{0};
      else
	{
	  auto __bx = __beta / __x;
	  return std::pow(__bx, __alpha)
	       * __gamma_reciprocal(__alpha)
	       * std::exp(-__bx) / __x;
	}
    }

  /**
   * The CDF of the inverse gamma distribution is:
   * @f[
   *   F(x|\alpha,\beta) = Q(\alpha, \beta/x)
   * @f]
   * Where @f$ Q(a,x) @f$ is the incomplete gamma function ratio.
   */
  template<typename _Tp>
    inverse_gamma_cdf(_Tp __alpha, _Tp __beta, _Tp __x)
    {
      if (__x == _Tp{0})
	return _Tp{0};
      else
	return __gnu_cxx::qgamma(__alpha, __beta / __x);
    }

int
main()
{
  ;
}
