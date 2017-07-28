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
   * @param[in] asub Array that contains the three coefficients of the quadratic equation
   *
   * @param[out] xreal1 Real part of first root
   * @param[out] xreal2 Real part of second root
   * @param[out] ximag  Imaginary part of roots
   */
  template<typename _Real>
    std::array<solution_t<_Real>, 2>
    quadratic(const std::array<_Real, 3>& asub)
    {
      std::array<solution_t<_Real>, 2> __ret;

      if (asub[2] == _Real{0})
        {
          // Equation is linear.
          __ret[0] = -asub[0] / asub[1];
        }
      else
        {
          // The discriminant of a quadratic equation
          auto discrm = asub[1] * asub[1] - 4 * asub[2] * asub[0];

          if (discrm < _Real{0})
            {
              // The roots are complex conjugates.
              auto xreal = -asub[1] / (2 * asub[2]);
              auto ximag = std::sqrt(std::abs(discrm)) / (2 * asub[2]);
	      __ret[0] = std::complex<_Real>(xreal, ximag);
	      __ret[1] = std::complex<_Real>(xreal, -ximag);
            }
          else
            {
              // The roots are real.
              _Real signb = (asub[1] < 0 ? -_Real{1} : _Real{1});
              _Real temp = -_Real{0.5L} * (asub[1] + signb * std::sqrt(discrm));
              __ret[0] = temp / asub[2];
              __ret[1] = asub[0] / temp;
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
   * @param[in] asub  Array that contains the four coefficients of
   *                  the cubic equation
   *
   * @param[out] xreal Array containing the real parts of the three
   *                   roots of the cubic equation. Note that xreal(1)
   *                   contains a real root only.
   * @param[out] ximag The imaginary part of the complex conjugate roots
   *                   of the cubic equation. Note that there are either
   *                   no complex roots or there are two complex roots
   *                   that are complex conjugates.
   */
  template<typename _Real>
    std::array<solution_t<_Real>, 3>
    cubic(const std::array<_Real, 4>& asub)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(asub[0]);

      if (asub[3] == _Real{0})
	return quadratic(std::experimental::make_array(asub[0], asub[1], asub[2]));

      std::array<solution_t<_Real>, 3> __ret;

      //  Normalize cubic equation coefficients.
      std::array<_Real, 4> acube;
      acube[3] = _Real{1};
      acube[2] = asub[2] / asub[3];
      acube[1] = asub[1] / asub[3];
      acube[0] = asub[0] / asub[3];

      if (acube[0] == _Real{0})
        {
          __ret[0] = _Real{0};
          std::array<_Real, 3> aq;
          aq[2] = acube[3];
          aq[1] = acube[2];
          aq[0] = acube[1];
          const auto __quad = quadratic(aq);
	  __ret[1] = __quad[0];
	  __ret[2] = __quad[1];
        }
      else
        {
          auto q = (acube[2] * acube[2] - _Real{3} * acube[1]) / _Real{9};
          auto qcube = q * q * q;
          auto r = (_Real{2} * acube[2] * acube[2] * acube[2]
                  - _Real{9} * acube[2] * acube[1]
                  + _Real{27} * acube[0]) / _Real{54};
          auto rsqed = r * r;
          if (qcube - rsqed > _Real{0})
            {
              //  Calculate the three real roots.
              auto theta = std::acos(r / std::sqrt(qcube));
              auto term1 = -2 * std::sqrt(q);
              __ret[0] = term1 * std::cos(theta / _Real{3}) - acube[2] / _Real{3};
              __ret[1] = term1 * std::cos((theta + _Real{2} * _S_pi) / _Real{3}) - acube[2] / _Real{3};
              __ret[2] = term1 * std::cos((theta + _Real{4} * _S_pi) / _Real{3}) - acube[2] / _Real{3};
            }
          else
            {
              //  Calculate the single real root.
              auto term1 = std::cbrt((std::sqrt(rsqed - qcube) + std::abs(r)));
              __ret[0] = std::copysign(term1 + q / term1, r) - (acube[2] / _Real{3});

              //  Find the other two roots which are complex conjugates.
              std::array<_Real, 3> aq;
              aq[2] = acube[3];
              aq[1] = acube[2] + acube[3] * __ret[0];
              aq[0] = aq[2] * __ret[0] + acube[1];
              auto quad = quadratic(aq);
              __ret[1] = quad[0];
              __ret[2] = quad[1];
            }
        }

      return __ret;
    }


  /**
   * @brief  Solves a quartic equation of the form:
   * 	     @f[
   * 		 a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
   * 	     @f]
   * 	     for its roots using brown's method of biquadratic equations.
   *
   * 	 Note 1. The coefficients @f$a_k@f$ are real.
   *
   * 	 Note 2. There are four roots:
   * 		 Case 1.  Four real roots
   * 		 Case 2.  Two real roots and two complex roots (conjugate pair)
   * 		 Case 3.  Four complex roots in conjugate pairs
   *
   * @param[in]  asub  Array that contains the five(5) coefficients of the quartic equation
   *
   * @param[out]  xreal  Array containing the real parts of the four roots
   * @param[out]  ximag  Array containing the imaginary parts of the four roots
   */
  template<typename _Real>
    std::array<solution_t<_Real>, 4>
    quartic(const std::array<_Real, 5>& asub)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(asub[0]);

      if (asub[4] == _Real{0})
	return cubic(std::experimental::make_array(asub[0], asub[1], asub[2], asub[3]));

      std::array<solution_t<_Real>, 4> __ret;

      //  The three coefficients of the quadratic equation.
      _Real aq[3];

      //  Normalize quartic equation coefficients.
      std::array<_Real, 5> aquart;
      aquart[4] = _Real{1};
      aquart[3] = asub[3] / asub[4];
      aquart[2] = asub[2] / asub[4];
      aquart[1] = asub[1] / asub[4];
      aquart[0] = asub[0] / asub[4];

      //  Calculate the coefficients of the resolvent cubic equation.
      std::array<_Real, 4> acube;
      acube[3] = _Real{1};
      acube[2] = -aquart[2];
      acube[1] = aquart[3] * aquart[1] - 4 * aquart[0];
      acube[0] = aquart[0] * (4 * aquart[2] - aquart[3] * aquart[3])
               - aquart[1] * aquart[1];

      //  Find the algebraically largest real root of the cubic equation
      //  Note: A cubic equation has either three real roots or one
      //        real root and two complex roots that are complex
      //        conjugates. If there is only a single real root then
      //        subroutine cubic always returns that single real root
      //        (and therefore the algebraically largest real root of
      //        the cubic equation) as rt(1).
      //  The imaginary part of the complex conjugate root of the resolvent cubic equation.
      _Real rti = 0;
      auto cube = cubic(acube, rti);
      if (rti == 0)
        {
          if (cube[0] < cube[1])
            std::swap(cube[0], cube[1]);
          if (cube[0] < cube[2])
            std::swap(cube[0], cube[2]);
        }
      //  Calculate the coefficients for the two quadratic equations
      auto capa = _Real{0.5L} * aquart[3];
      auto capb = _Real{0.5L} * cube[0];
      auto capc = std::sqrt(capa * capa - aquart[2] + cube[0]);
      auto capd = std::sqrt(capb * capb - aquart[0]);
      auto c1 = capa + capc;
      auto c2 = capa - capc;
      auto d1 = capb + capd;
      auto d2 = capb - capd;
      auto t1 = c1 * d2 + c2 * d1;
      auto t2 = c1 * d1 + c2 * d2;
      if (std::abs(t2 - aquart[1]) < std::abs(t1 - aquart[1]))
        std::swap(d1, d2);

      //  Coefficients for the first quadratic equation and find the roots...
      aq[2] = _Real{1};
      aq[1] = c1;
      aq[0] = d1;
      auto quad1 = quadratic(aq);
      __ret[0] = quad1[0];
      __ret[1] = quad1[1];

      //  Coefficients for the second quadratic equation and find the roots...
      aq[2] = _Real{1};
      aq[1] = c2;
      aq[0] = d2;
      auto quad2 = quadratic(aq);
      __ret[2] = quad2[0];
      __ret[3] = quad2[1];

      return __ret;
    }


#endif  //  SOLVER_TCC
