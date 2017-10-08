#ifndef SOLVER_LOW_DEGREE_TCC
#define SOLVER_LOW_DEGREE_TCC 1

#include <ext/math_const.h>

namespace __gnu_cxx
{

  /**
   * @brief Finds the roots of a quadratic equation of the form:
   * @f[
   *    a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for real coefficients @f$ a_k @f$.
   *
   * For non-degenerate coefficients two roots are returned:
   * Either the roots are real or the roots are a complex conjugate pair.
   *
   * If the quadratic coefficient @f$ a_2 @f$ is zero (degenerate case)
   * at most one valid root is returned.
   * If the linear coefficient @f$ a_1 @f$ is also zero
   * no valid root is returned.
   *
   * @param[in] _CC Array that contains the three coefficients
   *                  of the quadratic equation.
   */
  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 2>
    __quadratic(const _Iter& _CC)
    {
      std::array<solution_t<_Real>, 2> _ZZ;

      if (_CC[2] == _Real{0})
	{
	  // Equation is linear (or completely degenerate).
	  if (_CC[1] == _Real{0})
	    return _ZZ;
	  else
	    {
	      _ZZ[0] = -_CC[0] / _CC[1];
	      return _ZZ;
	    }
	}
      else if (_CC[0] == _Real{0})
	{
	  _ZZ[0] = _Real{0};
	  if (_CC[2] == _Real{0})
	    return _ZZ;
	  else
	    {
	      _ZZ[1] = -_CC[1] / _CC[2];
	      return _ZZ;
	    }
	}
      else
	{
	  // The discriminant of a quadratic equation
	  const auto _QQ = _CC[1] * _CC[1] - _Real{4} * _CC[2] * _CC[0];

	  if (_QQ < _Real{0})
	    {
	      // The roots are complex conjugates.
	      const auto _ReZZ = -_CC[1] / (_Real{2} * _CC[2]);
	      const auto _ImZZ = std::sqrt(std::abs(_QQ)) / (_Real{2} * _CC[2]);
	      _ZZ[0] = std::complex<_Real>(_ReZZ, -_ImZZ);
	      _ZZ[1] = std::complex<_Real>(_ReZZ, _ImZZ);
	    }
	  else
	    {
	      // The roots are real.
	      _Real __temp = -(_CC[1]
			+ std::copysign(std::sqrt(_QQ), _CC[1])) / _Real{2};
	      _ZZ[0] = __temp / _CC[2];
	      _ZZ[1] = _CC[0] / __temp;
	    }
	}

      return _ZZ;
    }


  /**
   * @brief Finds the roots of a cubic equation of the form:
   * @f[
   *    a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for real coefficients @f$ a_k @f$.
   *
   * In the non-degenerate case there are three roots:
   * - All three roots are real
   * - One root is real and the other two are a complex conjugate pair
   *
   * If the cubic coefficient @f$ a_3 @f$ is zero (degenerate case)
   * the problem is referred to the quadratic solver to return, at most,
   * two valid roots.
   *
   * @param[in] _CC Array that contains the four coefficients
   *                  of the cubic equation
   */
  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 3>
    __cubic(const _Iter& _CC)
    {
      using std::experimental::make_array;

      std::array<solution_t<_Real>, 3> _ZZ;

      if (_CC[3] == _Real{0})
	{
	  const auto _ZZ2 = __quadratic<_Real>(_CC);
	  _ZZ[0] = _ZZ2[0];
	  _ZZ[1] = _ZZ2[1];
	  return _ZZ;
	}
      else if (_CC[0] == _Real{0})
	{
	  _ZZ[0] = _Real{0};
	  const auto _ZZ2 = __quadratic<_Real>(make_array(_CC[1], _CC[2],
							  _CC[3]));
	  _ZZ[1] = _ZZ2[0];
	  _ZZ[2] = _ZZ2[1];
	  return _ZZ;
	}
      else
	{
	  //  Normalize cubic equation coefficients.
	  std::array<_Real, 4> _AA3;
	  _AA3[3] = _Real{1};
	  _AA3[2] = _CC[2] / _CC[3];
	  _AA3[1] = _CC[1] / _CC[3];
	  _AA3[0] = _CC[0] / _CC[3];

	  const auto _S_pi = __gnu_cxx::__const_pi(_CC[0]);
	  const auto _S_2pi = _Real{2} * _S_pi;
	  const auto _S_4pi = _Real{4} * _S_pi;
	  const auto _PP = _AA3[2] / _Real{3};
	  const auto _QQ = (_AA3[2] * _AA3[2] - _Real{3} * _AA3[1])
			 / _Real{9};
	  const auto _QQp3 = _QQ * _QQ * _QQ;
	  const auto _RR = (_Real{2} * _AA3[2] * _AA3[2] * _AA3[2]
			  - _Real{9} * _AA3[2] * _AA3[1]
			  + _Real{27} * _AA3[0]) / _Real{54};
	  const auto _RRp2 = _RR * _RR;

	  if (_QQp3 - _RRp2 > _Real{0})
	    {
	      //  Calculate the three real roots.
	      const auto __phi = std::acos(_RR / std::sqrt(_QQp3));
	      const auto __fact = -_Real{2} * std::sqrt(_QQ);
	      _ZZ[0] = __fact * std::cos(__phi / _Real{3}) - _PP;
	      _ZZ[1] = __fact * std::cos((__phi + _S_2pi) / _Real{3}) - _PP;
	      _ZZ[2] = __fact * std::cos((__phi + _S_4pi) / _Real{3}) - _PP;
	    }
	  else
	    {
	      //  Calculate the single real root.
	      const auto __fact = std::cbrt(std::sqrt(_RRp2 - _QQp3)
					     + std::abs(_RR));
	      const auto _BB = -std::copysign(__fact + _QQ / __fact, _RR);
	      _ZZ[0] = _BB - _PP;

	      //  Find the other two roots which are complex conjugates.
	      std::array<_Real, 3> _AA2;
	      _AA2[2] = _Real{1};
	      _AA2[1] = _BB;
	      _AA2[0] = _BB * _BB - _Real{3} * _QQ;
	      const auto _ZZ2 = __quadratic<_Real>(_AA2);
	      _ZZ[1] = std::get<2>(_ZZ2[0]) - _PP;
	      _ZZ[2] = std::get<2>(_ZZ2[1]) - _PP;
	    }

	  return _ZZ;
	}
    }


  /**
   * @brief Finds the roots a quartic equation of the form:
   * @f[
   * 	 a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for real coefficients @f$ a_k @f$.
   *
   * In the non-degenerate case there are four roots:
   * - All four roots are real
   * - Two roots real and two complex roots are a complex conjugate pair
   * - Four complex roots in two complex conjugate pairs
   *
   * If the qartic coefficient @f$ a_4 @f$ is zero (degenerate case)
   * the problem is referred to the cubic solver to return, at most,
   * three valid roots.
   *
   * @param[in] _CC Array that contains the five(5) coefficients
   *                  of the quartic equation.
   */
  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 4>
    __quartic(const _Iter& _CC)
    {
      using std::experimental::make_array;

      const auto _S_pi = __gnu_cxx::__const_pi(_CC[0]);

      std::array<solution_t<_Real>, 4> _ZZ;

      if (_CC[4] == _Real{0})
	{
	  const auto _ZZ3 = __cubic<_Real>(_CC);
	  _ZZ[0] = _ZZ3[0];
	  _ZZ[1] = _ZZ3[1];
	  _ZZ[2] = _ZZ3[2];
	  return _ZZ;
	}
      else if (_CC[0] == _Real{0})
	{
	  _ZZ[0] = _Real{0};
	  const auto _ZZ3 = __cubic<_Real>(make_array(_CC[1], _CC[2],
						      _CC[3], _CC[4]));
	  _ZZ[1] = _ZZ3[0];
	  _ZZ[2] = _ZZ3[1];
	  return _ZZ;
	}
      else
	{
	  //  Normalize quartic equation coefficients.
	  std::array<_Real, 5> _AA4;
	  _AA4[4] = _Real{1};
	  _AA4[3] = _CC[3] / _CC[4];
	  _AA4[2] = _CC[2] / _CC[4];
	  _AA4[1] = _CC[1] / _CC[4];
	  _AA4[0] = _CC[0] / _CC[4];

	  //  Calculate the coefficients of the resolvent cubic equation.
	  std::array<_Real, 4> _AA3;
	  _AA3[3] = _Real{1};
	  _AA3[2] = -_AA4[2];
	  _AA3[1] = _AA4[3] * _AA4[1] - _Real{4} * _AA4[0];
	  _AA3[0] = _AA4[0] * (_Real{4} * _AA4[2] - _AA4[3] * _AA4[3])
		  - _AA4[1] * _AA4[1];

	  // Find the algebraically largest real root of the cubic equation
	  // Note: A cubic equation has either three real roots or one
	  //       real root and two complex roots that are complex
	  //       conjugates. If there is only a single real root then
	  //       subroutine cubic always returns that single real root
	  //       (and therefore the algebraically largest real root of
	  //       the cubic equation) as root[0].
	  _Real _Z3max;
	  auto _ZZ3 = __cubic<_Real>(_AA3);
	  if (_ZZ3[1].index() == 1 && _ZZ3[2].index() == 1)
            {
	      // There is some horrible big with swap and this variant.
	      if (_ZZ3[0] < _ZZ3[1])
		//std::swap(_ZZ3[0], _ZZ3[1]);
		{
		  const auto __tmp = _ZZ3[0];
		  _ZZ3[0] = _ZZ3[1];
		  _ZZ3[1] = __tmp;
		}
	      if (_ZZ3[0] < _ZZ3[2])
		//std::swap(_ZZ3[0], _ZZ3[2]);
		{
		  const auto __tmp = _ZZ3[0];
		  _ZZ3[0] = _ZZ3[2];
		  _ZZ3[2] = __tmp;
		}
	      _Z3max = std::get<1>(_ZZ3[0]);
            }
	  else
	    _Z3max = std::get<1>(_ZZ3[0]);

	  //  Calculate the coefficients for the two quadratic equations
	  const auto __capa = _Real{0.5L} * _AA4[3];
	  const auto __capb = _Real{0.5L} * _Z3max;
	  const auto __capc = std::sqrt(__capa * __capa - _AA4[2] + _Z3max);
	  const auto __capd = std::sqrt(__capb * __capb - _AA4[0]);
	  const auto __cp = __capa + __capc;
	  const auto __cm = __capa - __capc;
	  auto __dp = __capb + __capd;
	  auto __dm = __capb - __capd;
	  const auto __t1 = __cp * __dm + __cm * __dp;
	  const auto __t2 = __cp * __dp + __cm * __dm;
	  if (std::abs(__t2 - _AA4[1]) < std::abs(__t1 - _AA4[1]))
	    std::swap(__dp, __dm);

	  //  Coefficients for the first quadratic equation and find the roots.
	  std::array<_Real, 3> _AA2;
	  _AA2[2] = _Real{1};
	  _AA2[1] = __cp;
	  _AA2[0] = __dp;
	  const auto _ZZ2p = __quadratic<_Real>(_AA2);
	  _ZZ[0] = _ZZ2p[0];
	  _ZZ[1] = _ZZ2p[1];

	  //  Coefficients for the second quadratic equation and find the roots.
	  _AA2[2] = _Real{1};
	  _AA2[1] = __cm;
	  _AA2[0] = __dm;
	  const auto _ZZ2m = __quadratic<_Real>(_AA2);
	  _ZZ[2] = _ZZ2m[0];
	  _ZZ[3] = _ZZ2m[1];

	  return _ZZ;
	}
    }

} // namespace __gnu_cxx

#endif // SOLVER_LOW_DEGREE_TCC
