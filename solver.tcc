#ifndef SOLVER_TCC
#define SOLVER_TCC 1

#include <ext/cmath>

  /**
   * @brief Finds the roots of a quadratic equation of the form:
   * @f[
   *    a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   *
   * @note The coefficients @f$a_k@f$ are real.
   *
   * @note There are two roots: Both are real roots or the roots
   *       are complex conjugates.
   *
   * @param[in] _Coef Array that contains the three coefficients of the quadratic equation
   */
  template<typename _Real>
    std::array<solution_t<_Real>, 2>
    quadratic(const std::array<_Real, 3>& _Coef)
    {
      std::array<solution_t<_Real>, 2> __ret;

      if (_Coef[2] == _Real{0})
	{
	  // Equation is linear.
	  __ret[0] = -_Coef[0] / _Coef[1];
	  return __ret;
	}
      else
	{
	  // The discriminant of a quadratic equation
	  auto __discrm = _Coef[1] * _Coef[1] - _Real{4} * _Coef[2] * _Coef[0];

	  if (__discrm < _Real{0})
	    {
	      // The roots are complex conjugates.
	      auto __xreal = -_Coef[1] / (_Real{2} * _Coef[2]);
	      auto __ximag = std::sqrt(std::abs(__discrm)) / (_Real{2} * _Coef[2]);
	      __ret[0] = std::complex<_Real>(__xreal, -__ximag);
	      __ret[1] = std::complex<_Real>(__xreal, __ximag);
	    }
	  else
	    {
	      // The roots are real.
	      _Real __temp = -(_Coef[1]
			+ std::copysign(std::sqrt(__discrm), _Coef[1])) / _Real{2};
	      __ret[0] = __temp / _Coef[2];
	      __ret[1] = _Coef[0] / __temp;
	    }
	}

      return __ret;
    }


  /**
   * @brief  Subroutine cubic solves a cubic equation of the form:
   * @f[
   *    a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for its roots.
   *
   * @note The coefficients @f$a_k@f$ are real.
   *
   * @note There are three roots:
   * 		 Case 1.  All three roots are real
   * 		 Case 2.  One root is real and the other two
   * 			  are complex conjugates
   *
   * @param[in] _Coef  Array that contains the four coefficients of
   *                  the cubic equation
   */
  template<typename _Real>
    std::array<solution_t<_Real>, 3>
    cubic(const std::array<_Real, 4>& _Coef)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(_Coef[0]);

      std::array<solution_t<_Real>, 3> __ret;

      if (_Coef[3] == _Real{0})
	{
	  const auto __quad = quadratic(std::experimental::make_array(_Coef[0], _Coef[1], _Coef[2]));
	  __ret[0] = __quad[0];
	  __ret[1] = __quad[1];
	  return __ret;
	}

      //  Normalize cubic equation coefficients.
      std::array<_Real, 4> __acube;
      __acube[3] = _Real{1};
      __acube[2] = _Coef[2] / _Coef[3];
      __acube[1] = _Coef[1] / _Coef[3];
      __acube[0] = _Coef[0] / _Coef[3];

      if (__acube[0] == _Real{0})
	{
	  __ret[0] = _Real{0};
	  std::array<_Real, 3> __aquad;
	  __aquad[2] = __acube[3];
	  __aquad[1] = __acube[2];
	  __aquad[0] = __acube[1];
	  const auto __quad = quadratic(__aquad);
	  __ret[1] = __quad[0];
	  __ret[2] = __quad[1];
	}
      else
	{
	  const auto __lambda = __acube[2] / _Real{3};
	  const auto _Q = (__acube[2] * __acube[2] - _Real{3} * __acube[1])
			/ _Real{9};
	  const auto _Qp3 = _Q * _Q * _Q;
	  const auto _R = (_Real{2} * __acube[2] * __acube[2] * __acube[2]
			 - _Real{9} * __acube[2] * __acube[1]
			 + _Real{27} * __acube[0]) / _Real{54};
	  const auto _Rp2 = _R * _R;
	  if (_Qp3 - _Rp2 > _Real{0})
	    {
	      //  Calculate the three real roots.
	      const auto __theta = std::acos(_R / std::sqrt(_Qp3));
	      const auto __fact = -_Real{2} * std::sqrt(_Q);
	      __ret[0] = __fact * std::cos(__theta / _Real{3}) - __lambda;
	      __ret[1] = __fact * std::cos((__theta + _Real{2} * _S_pi) / _Real{3}) - __lambda;
	      __ret[2] = __fact * std::cos((__theta + _Real{4} * _S_pi) / _Real{3}) - __lambda;
	    }
	  else
	    {
	      //  Calculate the single real root.
	      const auto __fact = std::cbrt((std::sqrt(_Rp2 - _Qp3) + std::abs(_R)));
	      const auto _BB = -std::copysign(__fact + _Q / __fact, _R);
	      __ret[0] = _BB - __lambda;

	      //  Find the other two roots which are complex conjugates.
	      std::array<_Real, 3> __aquad;
	      __aquad[2] = _BB * _BB - _Real{3} * _Q;
	      __aquad[1] = _BB;
	      __aquad[0] = _Real{1};
	      auto __quad = quadratic(__aquad);
	      __ret[1] = std::get<2>(__quad[0]) - __lambda;
	      __ret[2] = std::get<2>(__quad[1]) - __lambda;
	    }
	}

      return __ret;
    }


  /**
   * @brief  Solves a quartic equation of the form:
   * @f[
   * 	 a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for its roots using Brown's method of biquadratic equations.
   *
   * @note The coefficients @f$ a_k @f$ are real.
   *
   * @note There are four roots:
   * 	   * Four real roots
   * 	   * Two real roots and two complex roots (conjugate pair)
   * 	   * Four complex roots in conjugate pairs
   *
   * @param[in]  _Coef  Array that contains the five(5) coefficients
   *                   of the quartic equation.
   */
  template<typename _Real>
    std::array<solution_t<_Real>, 4>
    quartic(const std::array<_Real, 5>& _Coef)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(_Coef[0]);

      std::array<solution_t<_Real>, 4> __ret;

      if (_Coef[4] == _Real{0})
	{
	  const auto cube = cubic(std::experimental::make_array(_Coef[0], _Coef[1], _Coef[2], _Coef[3]));
	  __ret[0] = cube[0];
	  __ret[1] = cube[1];
	  __ret[2] = cube[2];
	  return __ret;
	}

      //  Normalize quartic equation coefficients.
      std::array<_Real, 5> aquart;
      aquart[4] = _Real{1};
      aquart[3] = _Coef[3] / _Coef[4];
      aquart[2] = _Coef[2] / _Coef[4];
      aquart[1] = _Coef[1] / _Coef[4];
      aquart[0] = _Coef[0] / _Coef[4];

      //  Calculate the coefficients of the resolvent cubic equation.
      std::array<_Real, 4> acube;
      acube[3] = _Real{1};
      acube[2] = -aquart[2];
      acube[1] = aquart[3] * aquart[1] - _Real{4} * aquart[0];
      acube[0] = aquart[0] * (_Real{4} * aquart[2] - aquart[3] * aquart[3])
	       - aquart[1] * aquart[1];

      // Find the algebraically largest real root of the cubic equation
      // Note: A cubic equation has either three real roots or one
      //       real root and two complex roots that are complex
      //       conjugates. If there is only a single real root then
      //       subroutine cubic always returns that single real root
      //       (and therefore the algebraically largest real root of
      //       the cubic equation) as root[0].
      _Real cube0max;
      auto cube = cubic(acube);
      if (cube[1].index() == 1 && cube[2].index() == 1)
        {
	  if (cube[0] < cube[1])
	    std::swap(cube[0], cube[1]);
	  if (cube[0] < cube[2])
	    std::swap(cube[0], cube[2]);
	  cube0max = std::get<1>(cube[0]);
        }
      else
	cube0max = std::get<1>(cube[0]);
      //  Calculate the coefficients for the two quadratic equations
      const auto capa = _Real{0.5L} * aquart[3];
      const auto capb = _Real{0.5L} * cube0max;
      const auto capc = std::sqrt(capa * capa - aquart[2] + cube0max);
      const auto capd = std::sqrt(capb * capb - aquart[0]);
      const auto c1 = capa + capc;
      const auto c2 = capa - capc;
      auto d1 = capb + capd;
      auto d2 = capb - capd;
      const auto t1 = c1 * d2 + c2 * d1;
      const auto t2 = c1 * d1 + c2 * d2;
      if (std::abs(t2 - aquart[1]) < std::abs(t1 - aquart[1]))
	std::swap(d1, d2);

      //  Coefficients for the first quadratic equation and find the roots...
      std::array<_Real, 3> aquad;
      aquad[2] = _Real{1};
      aquad[1] = c1;
      aquad[0] = d1;
      const auto quad1 = quadratic(aquad);
      __ret[0] = quad1[0];
      __ret[1] = quad1[1];

      //  Coefficients for the second quadratic equation and find the roots...
      aquad[2] = _Real{1};
      aquad[1] = c2;
      aquad[0] = d2;
      const auto quad2 = quadratic(aquad);
      __ret[2] = quad2[0];
      __ret[3] = quad2[1];

      return __ret;
    }


#endif  //  SOLVER_TCC
