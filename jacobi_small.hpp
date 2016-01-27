/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-07-26

  Copyright (C) 2007 Université Joseph Fourier

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file jacobi.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-07-26
 */
#ifndef __Jacobi_H
#define __Jacobi_H 1


namespace Life
{
  /**
   * \class Jacobi
   * \brief 1D Jacobi polynomial
   *
   * (excerpt from Karniadakis/Sherwin Appendix A)
   * Jacobi polynomials \f$ P^{\alpha,beta}_n(x)\f$ are a family of
   * polynomial to the singular Sturm-Liouville problem. A significant
   * feature of these polynomials is that they are orthogonal in the
   * interval \f$[-1,1]\f$ with respect to the function
   * \f$(1-x)^\alpha(1+x)^\beta (\alpha,\beta > -1)\f$
   *
   * Several functions related to the one-dimensional jacobi
   * polynomials: Evaluation, evaluation of derivatives, plus
   * computation of the roots via Newton's method.
   *
   * \ingroup Polynomial
   * @author Christophe Prud'homme
   * @see Karniadakis and Sherwin "Spectral/hp element methods for CFD"
   */
  template<typename T = double>
  class Jacobi
  {
  public:

    /** @name Static values
     */
    //@{

    //@}

    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef Jacobi<T> self_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{


    /**
     * default values for a and b give the special case of Legendre
     * polynomials
     */
    Jacobi(int N, value_type a = value_type(0.0), value_type b = value_type(0.0))
      :
      _M_degree(N),
      _M_a(a),
      _M_b(b)
    {}

    Jacobi(Jacobi const & p)
      :
      _M_degree(p._M_degree),
      _M_a(p._M_a),
      _M_b(p._M_b)
    {}

    ~Jacobi()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=(self_type const& p)
    {
      if (this != &p)
	{
	  _M_degree = p._M_degree;
	  _M_a = p._M_a;
	  _M_b = p._M_b;
	}
      return *this;
    }

    /**
     * Evaluates the nth jacobi polynomial with weight parameters a,b
     * at a point \p x. Recurrence relations implemented from the
     * pseudocode given in Karniadakis and Sherwin, Appendix B
     *
     * a and b are defaulted to 0 and the Jacobi polynomial is then
     * the Legendre polynomial
     *
     * \param x point for polynomial evaluation
     * \return the value of the jacobi polynomial at \c x
     */
    value_type operator()(value_type const& x) const;

    //@}

    /** @name Accessors
     */
    //@{

    int degree() const { return _M_degree; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setDegree(int N) { _M_degree = N; }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Evaluates the nth jacobi polynomial with weight parameters a,b
     * at a point x. Recurrence relations implemented from the
     * pseudocode given in Karniadakis and Sherwin, Appendix B
     *
     * a and b are defaulted to 0 and the Jacobi polynomial is then
     * the Legendre polynomial
     *
     * \param x point for polynomial evaluation
     * \return the value of the jacobi polynomial at \c x
     */
    value_type value(value_type const& x) const
    {
      return this->operator()(x);
    }

    //@}



  protected:

  private:
    int _M_degree;
    value_type _M_a;
    value_type _M_b;
  };

  template<typename T>
  typename Jacobi<T>::value_type
  Jacobi<T>::operator()(value_type const& x) const
  {
    const int N = this->_M_degree;
    const value_type one = 1.0;
    const value_type two = 2.0;
    if (N == 0)
      return one;
    else if (N == 1)
      return 0.5 * (_M_a - _M_b + (_M_a + _M_b + two) * x);
    else  // N >= 2
      {
        value_type apb = _M_a + _M_b;
        value_type pn2 = one;
        value_type pn1 = 0.5 * (_M_a - _M_b + (apb + two) * x);
        value_type p = 0.0;
        for (int k = 2; k < N+1; ++k)
	  {
            value_type kv = value_type(k);
            value_type a1 = two * kv * (kv + apb) * (two * kv + apb - two);
            value_type a2 = (two * kv + apb - one) * (_M_a * _M_a - _M_b * _M_b);
            value_type a3 = ((two * kv + apb - two)
                              * (two * kv + apb - one)
                              * (two * kv + apb));
            value_type a4 = (two * (kv + _M_a - one) * (kv + _M_b - one)
                              * (two * kv + apb));

            p = ((a2 + a3 * x) * pn1 - a4 * pn2)/a1;
            pn2 = pn1;
            pn1 = p;
	  }
        return p;
      }
  }

} // Life
#endif /* __Jacobi_H */

