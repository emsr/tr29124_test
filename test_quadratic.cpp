/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_quadratic test_quadratic.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_quadratic > test_quadratic.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_quadratic test_quadratic.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
PATH=wrappers/debug:$PATH ./test_quadratic > test_quadratic.txt
*/

#include <cmath>
#include <complex>
#include <experimental/array>
#include <variant>
#include <iostream>

template<typename _Tp>
  using root_t = std::variant<std::monostate, _Tp, std::complex<_Tp>>;

/**
 * This routine solves a quadratic equation:
 * @f[
 *    0 = ax^2 + bx + c
 * @f]
 * placing the (real) roots in r1 and r2
 * or the real and imaginary parts of a complex root in r1 and r2 respectively.
 * The number of unique real roots is returned.
 */
template<typename _Tp>
  std::array<root_t<_Tp>, 2> 
  quadratic(_Tp __c, _Tp __b, _Tp __a)
  {
    const auto __d = __b * __b - _Tp{4} * __a * __c;

    if (__d < _Tp{0})
      {
        // Set the first root to the real part of the true complex root
        // and the second root to the imaginary part of the true complex root.
        const auto __re = -__b / (_Tp{2} * __a);
        const auto __im = -std::copysign(std::sqrt(-__d), __b) / (_Tp{2} * __a);
	return std::experimental::make_array(root_t<_Tp>(std::complex<_Tp>(__re, __im)),
			       root_t<_Tp>(std::complex<_Tp>(__re, -__im)));
      }
    else if (__d == _Tp{0})
      return std::experimental::make_array(root_t<_Tp>(-__b / (_Tp{2} * __a)), root_t<_Tp>());
    else
      {
        const auto __q = _Tp{-0.5} * (__b + std::copysign(std::sqrt(__d), __b));
        return std::experimental::make_array(root_t<_Tp>(__q / __a), root_t<_Tp>(__c / __q));
      }
  }

/**
 * This routine solves a quadratic equation:
 * @f[
 *    0 = ax^2 + bx + c
 * @f]
 * placing the (complex) roots in r1 and r2.
 */
template<typename _Tp>
  std::array<std::complex<_Tp>, 2>
  quad_complex(std::complex<_Tp> __c, std::complex<_Tp> __b, std::complex<_Tp> __a)
  {
    const auto __d = __b * __b - _Tp{4} * __a * __c;

    const auto __sign = std::real(std::conj(__b) * __d) > _Tp{0} ? _Tp{+1} : _Tp{-1};
    const auto __q = _Tp{-0.5} * (__b + __sign * std::sqrt(__d));
    return std::experimental::make_array(__q / __a, __c / __q);
  }

int
main()
{
  auto a1 = 1.0;
  auto b1 = -4.0;
  auto c1 = 1.0;
  auto r1 = quadratic(c1, b1, a1);
  std::cout << ' ' << std::get<1>(r1[0]) << ' ' << std::get<1>(r1[1]) << '\n';
  {
    auto x1 = std::get<1>(r1[0]);
    auto y1 = a1 * x1 * x1 + b1 * x1 + c1;
    std::cout << " y1 = " << y1 << '\n';
    auto x2 = std::get<1>(r1[1]);
    auto y2 = a1 * x2 * x2 + b1 * x2 + c1;
    std::cout << " y2 = " << y2 << '\n';
  }
}
