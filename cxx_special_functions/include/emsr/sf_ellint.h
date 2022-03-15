#ifndef SF_ELLINT_H
#define SF_ELLINT_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_ellint.tcc>

namespace emsr
{

  // Complete elliptic integrals of the first kind

  /**
   * Return the complete elliptic integral of the first kind
   * @f$ K(k) @f$ for real modulus @c k.
   *
   * The complete elliptic integral of the first kind is defined as
   * @f[
   *   K(k) = F(k,\pi/2) = \int_0^{\pi/2}\frac{d\theta}
   * 					     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   * where @f$ F(k,\phi) @f$ is the incomplete elliptic integral of the
   * first kind and the modulus @f$ |k| <= 1 @f$.
   * @see ellint_1 for details of the incomplete elliptic function
   * of the first kind.
   *
   * @tparam Tp The floating-point type of the modulus @c k.
   * @param  k  The modulus, <tt> abs(k) <= 1 </tt>
   * @throw std::domain_error if <tt> abs(k) > 1 </tt>.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    comp_ellint_1(Tp k)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::comp_ellint_1<type>(k);
    }

  // Complete elliptic integrals of the second kind

  /**
   * Return the complete elliptic integral of the second kind @f$ E(k) @f$
   * for real modulus @c k.
   *
   * The complete elliptic integral of the second kind is defined as
   * @f[
   *   E(k) = E(k,\pi/2) = \int_0^{\pi/2}\sqrt{1 - k^2 sin^2\theta}
   * @f]
   * where @f$ E(k,\phi) @f$ is the incomplete elliptic integral of the
   * second kind and the modulus @f$ |k| <= 1 @f$.
   * @see ellint_2 for details of the incomplete elliptic function
   * of the second kind.
   *
   * @tparam Tp The floating-point type of the modulus @c k.
   * @param  k  The modulus, @c abs(k) <= 1
   * @throw std::domain_error if @c abs(k) > 1.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    comp_ellint_2(Tp k)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::comp_ellint_2<type>(k);
    }

  // Complete elliptic integrals of the third kind

  /**
   * Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) = \Pi(k,\nu,\pi/2) @f$ for real modulus @c k.
   *
   * The complete elliptic integral of the third kind is defined as
   * @f[
   *   \Pi(k,\nu) = \Pi(k,\nu,\pi/2) = \int_0^{\pi/2}
   * 		     \frac{d\theta}
   * 		   {(1 - \nu \sin^2\theta)\sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   * where @f$ \Pi(k,\nu,\phi) @f$ is the incomplete elliptic integral of the
   * second kind and the modulus @f$ |k| <= 1 @f$.
   * @see ellint_3 for details of the incomplete elliptic function
   * of the third kind.
   *
   * @tparam Tp The floating-point type of the modulus @c k.
   * @tparam _Tpn The floating-point type of the argument @c nu.
   * @param  k  The modulus, @c abs(k) <= 1.
   * @param  nu  The characteristic.
   * @throw std::domain_error if @c abs(k) > 1.
   */
  template<typename Tp, typename _Tpn>
    inline emsr::fp_promote_t<Tp, _Tpn>
    comp_ellint_3(Tp k, _Tpn nu)
    {
      using type = emsr::fp_promote_t<Tp, _Tpn>;
      return emsr::detail::comp_ellint_3<type>(k, nu);
    }

  // Incomplete elliptic integrals of the first kind

  /**
   * Return the incomplete elliptic integral of the first kind @f$ F(k,\phi) @f$
   * for @c real modulus @c k and angle @f$ \phi @f$.
   *
   * The incomplete elliptic integral of the first kind is defined as
   * @f[
   *   F(k,\phi) = \int_0^{\phi}\frac{d\theta}
   * 				     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   * For  @f$ \phi= \pi/2 @f$ this becomes the complete elliptic integral of
   * the first kind, @f$ K(k) @f$.  @see comp_ellint_1.
   *
   * @tparam Tp The floating-point type of the modulus @c k.
   * @tparam _Tpp The floating-point type of the angle @c phi.
   * @param  k  The modulus, <tt> abs(k) <= 1 </tt>.
   * @param  phi  The integral limit argument in radians.
   * @throw std::domain_error if <tt> abs(k) > 1 </tt>.
   */
  template<typename Tp, typename _Tpp>
    inline emsr::fp_promote_t<Tp, _Tpp>
    ellint_1(Tp k, _Tpp phi)
    {
      using type = emsr::fp_promote_t<Tp, _Tpp>;
      return emsr::detail::ellint_1<type>(k, phi);
    }

  // Incomplete elliptic integrals of the second kind

  /**
   * Return the incomplete elliptic integral of the second kind
   * @f$ E(k,\phi) @f$.
   *
   * The incomplete elliptic integral of the second kind is defined as
   * @f[
   *   E(k,\phi) = \int_0^{\phi} \sqrt{1 - k^2 sin^2\theta}
   * @f]
   * For  @f$ \phi= \pi/2 @f$ this becomes the complete elliptic integral of
   * the second kind, @f$ E(k) @f$.  @see comp_ellint_2.
   *
   * @tparam Tp The floating-point type of the modulus @c k.
   * @tparam _Tpp The floating-point type of the angle @c phi.
   * @param  k  The modulus, <tt> abs(k) <= 1 </tt>
   * @param  phi  The integral limit argument in radians
   * @return  The elliptic function of the second kind.
   * @throw std::domain_error if <tt> abs(k) > 1 </tt>.
   */
  template<typename Tp, typename _Tpp>
    inline emsr::fp_promote_t<Tp, _Tpp>
    ellint_2(Tp k, _Tpp phi)
    {
      using type = emsr::fp_promote_t<Tp, _Tpp>;
      return emsr::detail::ellint_2<type>(k, phi);
    }

  // Incomplete elliptic integrals of the third kind

  /**
   * @brief Return the incomplete elliptic integral of the third kind
   * @f$ \Pi(k,\nu,\phi) @f$.
   *
   * The incomplete elliptic integral of the third kind is defined by:
   * @f[
   *   \Pi(k,\nu,\phi) = \int_0^{\phi}
   * 			 \frac{d\theta}
   * 			 {(1 - \nu \sin^2\theta)
   * 			  \sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   * For  @f$ \phi= \pi/2 @f$ this becomes the complete elliptic integral of
   * the third kind, @f$ \Pi(k,\nu) @f$.  @see comp_ellint_3.
   *
   * @tparam Tp The floating-point type of the modulus @c k.
   * @tparam _Tpn The floating-point type of the argument @c nu.
   * @tparam _Tpp The floating-point type of the angle @c phi.
   * @param  k  The modulus, <tt> abs(k) <= 1 </tt>.
   * @param  nu  The characteristic.
   * @param  phi  The integral limit argument in radians.
   * @return  The elliptic function of the third kind.
   * @throw std::domain_error if <tt> abs(k) > 1 </tt>.
   */
  template<typename Tp, typename _Tpn, typename _Tpp>
    inline emsr::fp_promote_t<Tp, _Tpn, _Tpp>
    ellint_3(Tp k, _Tpn nu, _Tpp phi)
    {
      using type = emsr::fp_promote_t<Tp, _Tpn, _Tpp>;
      return emsr::detail::ellint_3<type>(k, nu, phi);
    }

  // Complete Carlson elliptic R_F functions

  /**
   * Return the complete Carlson elliptic function @f$ R_F(x,y) @f$
   * for real arguments.
   *
   * The complete Carlson elliptic function of the first kind is defined by:
   * @f[
   *    R_F(x,y) = R_F(x,y,y) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)}
   * @f]
   *
   * @param  x  The first argument.
   * @param  y  The second argument.
   */
  template<typename _Tx, typename _Ty>
    inline emsr::fp_promote_t<_Tx, _Ty>
    comp_ellint_rf(_Tx x, _Ty y)
    {
      using type = emsr::fp_promote_t<_Tx, _Ty>;
      return emsr::detail::comp_ellint_rf<type>(x, y);
    }

  // Carlson elliptic R_F functions

  /**
   * Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   * of the first kind for real arguments.
   *
   * The Carlson elliptic function of the first kind is defined by:
   * @f[
   *    R_F(x,y,z) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}}
   * @f]
   *
   * @param  x  The first of three symmetric arguments.
   * @param  y  The second of three symmetric arguments.
   * @param  z  The third of three symmetric arguments.
   */
  template<typename Tp, typename _Up, typename _Vp>
    inline emsr::fp_promote_t<Tp, _Up, _Vp>
    ellint_rf(Tp x, _Up y, _Vp z)
    {
      using type = emsr::fp_promote_t<Tp, _Up, _Vp>;
      return emsr::detail::ellint_rf<type>(x, y, z);
    }

  // Carlson elliptic R_C functions

  /**
   * Return the Carlson elliptic function @f$ R_C(x,y) = R_F(x,y,y) @f$
   * where @f$ R_F(x,y,z) @f$ is the Carlson elliptic function
   * of the first kind.
   *
   * The Carlson elliptic function is defined by:
   * @f[
   *    R_C(x,y) = \frac{1}{2} \int_0^\infty
   *               \frac{dt}{(t + x)^{1/2}(t + y)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *    by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first argument.
   * @param  y  The second argument.
   */
  template<typename Tp, typename _Up>
    inline emsr::fp_promote_t<Tp, _Up>
    ellint_rc(Tp x, _Up y)
    {
      using type = emsr::fp_promote_t<Tp, _Up>;
      return emsr::detail::ellint_rc<type>(x, y);
    }

  // Carlson elliptic R_J functions

  /**
   * Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$
   * 	     of the third kind.
   *
   * The Carlson elliptic function of the third kind is defined by:
   * @f[
   *    R_J(x,y,z,p) = \frac{3}{2} \int_0^\infty
   *                   \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}(t + p)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first of three symmetric arguments.
   * @param  y  The second of three symmetric arguments.
   * @param  z  The third of three symmetric arguments.
   * @param  p  The fourth argument.
   */
  template<typename Tp, typename _Up, typename _Vp, typename _Wp>
    inline emsr::fp_promote_t<Tp, _Up, _Vp, _Wp>
    ellint_rj(Tp x, _Up y, _Vp z, _Wp p)
    {
      using type = emsr::fp_promote_t<Tp, _Up, _Vp, _Wp>;
      return emsr::detail::ellint_rj<type>(x, y, z, p);
    }

  // Carlson elliptic R_D functions

  /**
   * Return the Carlson elliptic function of the second kind
   * @f$ R_D(x,y,z) = R_J(x,y,z,z) @f$ where
   * @f$ R_J(x,y,z,p) @f$ is the Carlson elliptic function
   * of the third kind.
   *
   * The Carlson elliptic function of the second kind is defined by:
   * @f[
   *    R_D(x,y,z) = \frac{3}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{3/2}}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *    by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first of two symmetric arguments.
   * @param  y  The second of two symmetric arguments.
   * @param  z  The third argument.
   */
  template<typename Tp, typename _Up, typename _Vp>
    inline emsr::fp_promote_t<Tp, _Up, _Vp>
    ellint_rd(Tp x, _Up y, _Vp z)
    {
      using type = emsr::fp_promote_t<Tp, _Up, _Vp>;
      return emsr::detail::ellint_rd<type>(x, y, z);
    }

  // Complete Carlson elliptic R_G functions

  /**
   * Return the complete Carlson elliptic function @f$ R_G(x,y) @f$
   * for real arguments.
   *
   * The complete Carlson elliptic function is defined by:
   * @f[
   *    R_G(x,y) = R_G(x,y,y) = \frac{1}{4} \int_0^\infty
   *            dt t (t + x)^{-1/2}(t + y)^{-1}
   *         \left(\frac{x}{t + x} + \frac{2y}{t + y}\right)
   * @f]
   *
   * @param  x  The first argument.
   * @param  y  The second argument.
   */
  template<typename _Tx, typename _Ty>
    inline emsr::fp_promote_t<_Tx, _Ty>
    comp_ellint_rg(_Tx x, _Ty y)
    {
      using type = emsr::fp_promote_t<_Tx, _Ty>;
      return emsr::detail::comp_ellint_rg<type>(x, y);
    }

  // Carlson elliptic R_G functions

  /**
   * Return the symmetric Carlson elliptic function of the second kind
   * @f$ R_G(x,y,z) @f$.
   *
   * The Carlson symmetric elliptic function of the second kind is defined by:
   * @f[
   *    R_G(x,y,z) = \frac{1}{4} \int_0^\infty
   *            dt t [(t + x)(t + y)(t + z)]^{-1/2}
   *         (\frac{x}{t + x} + \frac{y}{t + y} + \frac{z}{t + z})
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *    by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first of three symmetric arguments.
   * @param  y  The second of three symmetric arguments.
   * @param  z  The third of three symmetric arguments.
   */
  template<typename Tp, typename _Up, typename _Vp>
    inline emsr::fp_promote_t<Tp, _Up, _Vp>
    ellint_rg(Tp x, _Up y, _Vp z)
    {
      using type = emsr::fp_promote_t<Tp, _Up, _Vp>;
      return emsr::detail::ellint_rg<type>(x, y, z);
    }

  // Jacobi zeta functions.

  /**
   * Return the Jacobi zeta function of @c k and @f$ \phi @f$.
   *
   * The Jacobi zeta function is defined by
   * @f[
   *    Z(m,\phi) = E(m,\phi) - \frac{E(m)F(m,\phi)}{K(m)}
   * @f]
   * where @f$ E(m,\phi) @f$ is the elliptic function of the second kind,
   * @f$ E(m) @f$ is the complete ellitic function of the second kind,
   * and @f$ F(m,\phi) @f$ is the elliptic function of the first kind.
   *
   * @tparam _Tk the real type of the modulus
   * @tparam _Tphi the real type of the angle limit
   * @param k The modulus
   * @param phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline emsr::fp_promote_t<_Tk, _Tphi>
    jacobi_zeta(_Tk k, _Tphi phi)
    {
      using type = emsr::fp_promote_t<_Tk, _Tphi>;
      return emsr::detail::jacobi_zeta<type>(k, phi);
    }

  // Heuman lambda functions.

  /**
   * Return the Heuman lambda function @f$ \Lambda(k,\phi) @f$
   * of modulus @c k and angular limit @f$ \phi @f$.
   *
   * The complete Heuman lambda function is defined by
   * @f[
   *    \Lambda(k,\phi) = \frac{F(1-m,\phi)}{K(1-m)}
   *    + \frac{2}{\pi} K(m) Z(1-m,\phi)
   * @f]
   * where @f$ m = k^2 @f$, @f$ K(k) @f$ is the complete elliptic function
   * of the first kind, and @f$ Z(k,\phi) @f$ is the Jacobi zeta function.
   *
   * @tparam _Tk the floating-point type of the modulus
   * @tparam _Tphi the floating-point type of the angular limit argument
   * @param k The modulus
   * @param phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline emsr::fp_promote_t<_Tk, _Tphi>
    heuman_lambda(_Tk k, _Tphi phi)
    {
      using type = emsr::fp_promote_t<_Tk, _Tphi>;
      return emsr::detail::heuman_lambda<type>(k, phi);
    }

  // Complete Legendre elliptic integral D.

  /**
   * Return the complete Legendre elliptic integral @f$ D(k) @f$
   * of real modulus @c k.
   *
   * The complete Legendre elliptic integral D is defined by
   * @f[
   *    D(k) = \int_0^{\pi/2} \frac{\sin^2\theta d\theta}{\sqrt{1-k^2sin2\theta}}
   * @f]
   *
   * @tparam _Tk The type of the modulus @c k
   * @param k The modulus <tt>-1 <= k <= +1</tt>
   */
  template<typename _Tk>
    inline emsr::fp_promote_t<_Tk>
    comp_ellint_d(_Tk k)
    {
      using type = emsr::fp_promote_t<_Tk>;
      return emsr::detail::comp_ellint_d<type>(k);
    }

  // Legendre elliptic integrals D.

  /**
   * Return the incomplete Legendre elliptic integral @f$ D(k,\phi) @f$
   * of real modulus @c k and angular limit @f$ \phi @f$.
   *
   * The Legendre elliptic integral D is defined by
   * @f[
   *    D(k,\phi) = \int_0^\phi
   *         \frac{\sin^2\theta d\theta}{\sqrt{1-k^2sin^2\theta}}
   * @f]
   *
   * @param k The modulus <tt>-1 <= k <= +1</tt>
   * @param phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline emsr::fp_promote_t<_Tk, _Tphi>
    ellint_d(_Tk k, _Tphi phi)
    {
      using type = emsr::fp_promote_t<_Tk, _Tphi>;
      return emsr::detail::ellint_d<type>(k, phi);
    }

  // Bulirsch elliptic integrals of the first kind.

  /**
   * Return the Bulirsch elliptic integral @f$ el1(x,k_c) @f$
   * of the first kind of real tangent limit @c x
   * and complementary modulus @f$ k_c @f$.
   *
   * The Bulirsch elliptic integral of the first kind is defined by
   * @f[
   *    el1(x,k_c) = el2(x,k_c,1,1) = \int_0^{\arctan x} \frac{1+1\tan^2\theta}
   *           {\sqrt{(1+\tan^2\theta)(1+k_c^2\tan^2\theta)}}d\theta
   * @f]
   *
   * @param x The tangent of the angular integration limit
   * @param k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   */
  template<typename Tp, typename _Tk>
    inline emsr::fp_promote_t<Tp, _Tk>
    ellint_el1(Tp x, _Tk k_c)
    {
      using type = emsr::fp_promote_t<Tp, _Tk>;
      return emsr::detail::ellint_el1<type>(x, k_c);
    }

  // Bulirsch elliptic integrals of the second kind.

  /**
   * Return the Bulirsch elliptic integral of the second kind
   * @f$ el2(x,k_c,a,b) @f$.
   *
   * The Bulirsch elliptic integral of the second kind is defined by
   * @f[
   *    el2(x,k_c,a,b) = \int_0^{\arctan x} \frac{a+b\tan^2\theta}
   *           {\sqrt{(1+\tan^2\theta)(1+k_c^2\tan^2\theta)}}d\theta
   * @f]
   *
   * @param x The tangent of the angular integration limit
   * @param k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   * @param a The  parameter
   * @param b The  parameter
   */
  template<typename Tp, typename _Tk, typename _Ta, typename _Tb>
    inline emsr::fp_promote_t<Tp, _Tk, _Ta, _Tb>
    ellint_el2(Tp x, _Tk k_c, _Ta a, _Tb b)
    {
      using type = emsr::fp_promote_t<Tp, _Tk, _Ta, _Tb>;
      return emsr::detail::ellint_el2<type>(x, k_c, a, b);
    }

  // Bulirsch elliptic integrals of the third kind.

  /**
   * Return the Bulirsch elliptic integral of the third kind
   * @f$ el3(x,k_c,p) @f$ of real tangent limit @c x,
   * complementary modulus @f$ k_c @f$, and parameter @c p.
   *
   * The Bulirsch elliptic integral of the third kind is defined by
   * @f[
   *    el3(x,k_c,p) = \int_0^{\arctan x} \frac{d\theta}
   *          {(cos^2\theta+p\sin^2\theta)\sqrt{cos^2\theta+k_c^2\sin^2\theta}}
   * @f]
   *
   * @param x The tangent of the angular integration limit
   * @param k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   * @param p The paramenter
   */
  template<typename _Tx, typename _Tk, typename Tp>
    inline emsr::fp_promote_t<_Tx, _Tk, Tp>
    ellint_el3(_Tx x, _Tk k_c, Tp p)
    {
      using type = emsr::fp_promote_t<_Tx, _Tk, Tp>;
      return emsr::detail::ellint_el3<type>(x, k_c, p);
    }

  // Bulirsch complete elliptic integrals.

  /**
   * Return the Bulirsch complete elliptic integral @f$ cel(k_c,p,a,b) @f$
   * of real complementary modulus @f$ k_c @f$, and parameters @c p,
   * @c a, and @c b.
   *
   * The Bulirsch complete elliptic integral is defined by
   * @f[
   *    cel(k_c,p,a,b)=\int_0^{\pi/2}
   *        \frac{a\cos^2\theta + b\sin^2\theta}{cos^2\theta + p\sin^2\theta}
   *        \frac{d\theta}{\sqrt{cos^2\theta + k_c^2\sin^2\theta}}
   * @f]
   *
   * @param k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   * @param p The  parameter
   * @param a The  parameter
   * @param b The  parameter
   */
  template<typename _Tk, typename Tp, typename _Ta, typename _Tb>
    inline emsr::fp_promote_t<_Tk, Tp, _Ta, _Tb>
    ellint_cel(_Tk k_c, Tp p, _Ta a, _Tb b)
    {
      using type = emsr::fp_promote_t<_Tk, Tp, _Ta, _Tb>;
      return emsr::detail::ellint_cel<type>(k_c, p, a, b);
    }

} // namespace emsr

#endif // SF_ELLINT_H
