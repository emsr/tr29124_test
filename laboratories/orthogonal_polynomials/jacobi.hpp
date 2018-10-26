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

#include <ublas.hpp>


namespace Life
{
  namespace ublas=boost::numeric::ublas;
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
    Jacobi( int N, value_type a = value_type( 0.0 ), value_type b = value_type( 0.0 ) )
      :
      _M_degree( N ),
      _M_a( a ),
      _M_b( b )
    {}

    Jacobi( Jacobi const & p )
      :
      _M_degree( p._M_degree ),
      _M_a( p._M_a ),
      _M_b( p._M_b )
    {}

    ~Jacobi()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& p )
    {
      if ( this != &p )
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
    value_type operator()( value_type const& x ) const;

    //@}

    /** @name Accessors
     */
    //@{

    int degree() const { return _M_degree; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setDegree( int N ) { _M_degree = N; }

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
    value_type value( value_type const& x ) const
    {
      return this->operator()( x );
    }

    value_type derivate( value_type const& x ) const;

    //@}



  protected:

  private:
    int _M_degree;
    value_type _M_a;
    value_type _M_b;
  };
  template<typename T>
  typename Jacobi<T>::value_type
  Jacobi<T>::operator()( value_type const& x ) const
  {
    const int N = this->_M_degree;
    const value_type one = 1.0;
    const value_type two = 2.0;
    if ( N == 0 )
      return one;
    else if ( N == 1 )
      return 0.5 * ( _M_a - _M_b + ( _M_a + _M_b + two ) * x );
    else  // N >= 2
      {
        value_type apb = _M_a + _M_b;
        value_type pn2 = one;
        value_type pn1 = 0.5 * ( _M_a - _M_b + ( apb + two ) * x );
        value_type p = 0.0;
        for ( int k = 2; k < N+1; ++k )
	  {
            value_type kv = value_type(k);
            value_type a1 = two * kv * ( kv + apb ) * ( two * kv + apb - two );
            value_type a2 = ( two * kv + apb - one ) * ( _M_a * _M_a - _M_b * _M_b );
            value_type a3 = ( ( two * kv + apb - two )
                              * ( two * kv + apb - one )
                              * ( two * kv + apb ) );
            value_type a4 = ( two * ( kv + _M_a - one ) * ( kv + _M_b - one )
                              * ( two * kv + apb ) );

            p = ( ( a2 + a3 * x ) * pn1 - a4 * pn2 )/a1;
            pn2 = pn1;
            pn1 = p;
	  }
        return p;
      }
  }
  template<typename T>
  typename Jacobi<T>::value_type
  Jacobi<T>::derivate( value_type const& x ) const
  {
    const int N = this->_M_degree;
    if (  N == 0 )
      return 0.0;

    Jacobi<T> dp( N-1, _M_a + 1.0, _M_b + 1.0 );
    value_type Nv = value_type( N );
    return 0.5 * ( _M_a  + _M_b + Nv + 1.0 ) * dp( x );
  }
  template<typename T>
  ublas::matrix<T>
  JacobiBatchEvaluation( int N, T a, T b, ublas::vector<T> const& __pts )
  {
    typedef T value_type;
    ublas::matrix<T> res( N+1, __pts.size() );
    ublas::row( res, 0 ) = ublas::scalar_vector<value_type>( res.size2(), 1.0 );
    if ( N > 0 )
      {
        ublas::row( res, 1 ) = 0.5 * ( ublas::scalar_vector<value_type>( res.size2(), a - b) + ( a + b + 2.0 ) * __pts );
        value_type apb = a + b;
        for ( int k = 2; k < N+1; ++k )
	  {
            value_type kv = value_type(k);
            value_type a1 = 2.0 * kv * ( kv + apb ) * ( 2.0 * kv + apb - 2.0 );
            value_type a2 = ( 2.0 * kv + apb - 1.0 ) * ( a * a - b * b );
            value_type a3 = ( 2.0 * kv + apb - 2.0 ) * ( 2.0 * kv + apb - 1.0 ) * ( 2.0 * kv + apb );
            value_type a4 = 2.0 * ( kv + a - 1.0 ) * ( kv + b - 1.0 ) * ( 2.0 * kv + apb );
            a2 = a2 / a1;
            a3 = a3 / a1;
            a4 = a4 / a1;

            ublas::row( res, k ) = a2* ublas::row( res, k-1 ) +
	      a3 * ublas::element_prod(__pts, ublas::row( res, k-1 ) ) -
	      a4 * ublas::row( res, k-2 );
	  }
      }
    return res;
  }

  template<typename T>
  ublas::matrix<T>
  JacobiBatchDerivation( int N, T a, T b, ublas::vector<T> const& __pts )
  {
    typedef T value_type;
    ublas::matrix<T> res( N+1, __pts.size() );
    ublas::row( res, 0 ) = ublas::scalar_vector<value_type>( res.size2(), 0.0 );
    if ( N > 0 )
      {
        ublas::subrange( res, 1, N+1, 0, __pts.size() ) = JacobiBatchEvaluation<T>( N-1, a+1.0, b+1.0, __pts );
        for ( int i = 1;i < N+1; ++i )
	  ublas::row( res, i ) *= 0.5*(a+b+value_type( i )+1.0);
      }
    return res;
  }

  /**
   * Computes the m roots of \f$P_{m}^{a,b}\f$ on \f$[-1,1]\f$ by
   * Newton's method.  The initial guesses are the Chebyshev points as
   * we know an explicit formula for these polynomials.
   */
  template<typename JacobiP, typename Vector>
  void
  roots( JacobiP const& p, Vector& xr )
  {
    const int N = p.degree();
    typedef typename JacobiP::value_type value_type;

    if ( N != 0 )
      {
	value_type eps = 1e-16;
	value_type r;
	int max_iter = 30;
	for ( int k = 0;k < N; ++k )
	  {
	    value_type pi = 4.0*std::atan( value_type( 1.0 ) );
	    // use k-th checbychev point to  initiliaze newton
	    r = -std::cos(( 2.0*value_type( k ) + 1.0) * pi / ( 2.0 * value_type( N ) ) );
	    // use average of r and xr[k-1] as starting point (see KS)
	    if ( k > 0 )
	      r = 0.5 * ( r + xr[k-1] );
	    int j = 0;
	    value_type jf = 2.0*eps;
	    value_type delta  = 2.0*eps;
	    do
	      {
		// use deflation as proposed in KS
		value_type s = 0.0;
		for ( int i = 0;i < k; ++i )
		  {
		    s +=  value_type( 1.0 ) / ( r - xr[i] );
		  }
		jf = p( r );
		value_type jdf = p.derivate( r );
		delta = jf / (jdf - jf * s);

		// newton step done
		r = r - delta;

		++j;
	      }
	    while ( std::abs( jf ) > eps && j < max_iter );
	    // store k-th root
	    xr[k] = r;
	  }
      }
  }

  template<typename T>
  T fact( T  n )
  {
    if ( n == 0.0 )
      return 1;
    return n*fact( n-1.0 );
  }
  template<typename T, typename VectorW,  typename VectorN>
  void
  gaussjacobi( int N, VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ) )
  {
    typedef T value_type;

    Jacobi<T> p( N, a, b );

    roots( p, xr );

    const value_type two = 2.0;
    const value_type power = a+b+1.0;
    value_type a1 = std::pow( two,power);
    value_type a2 = fact( a+value_type( N ) );//gamma(a + m + 1);
    value_type a3 = fact( b+value_type( N ) );//gamma(b + m + 1);
    value_type a4 = fact( a+b+value_type( N ) );//gamma(a + b + m + 1);
    value_type a5 = fact(value_type(N));
    value_type a6 = a1 * a2 * a3;// / a4 / a5;

    for ( int k = 0; k < N; ++k )
      {
        value_type fp = p.derivate( xr[k] );
        value_type dn = fp*fp*( 1.0  - xr[k]*xr[k])*a4*a5;

        wr[k] =a6 / dn;
        //wr[k] = a6 / ( 1.0  - xr[k]*xr[k]) /  ( fp*fp );
      }
  }
  /**
   * \brief Gauss-Lobatto-Jacobi-Bouzitat Quadrature
   *
   * @see Karniadakis et Sherwin appendix B for theoretical details.
   **/
  template<typename T, typename VectorW,  typename VectorN>
  void
  gausslobattojacobi( int N, VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ), bool interior = false  )
  {
    wr.resize( N - 2*interior );
    xr.resize( N - 2*interior );

    typedef T value_type;

    Life::Jacobi<T> p( N-2, a+ 1.0, b+ 1.0 );

    Life::Jacobi<T> q( N-1, a, b );

    VectorN prexr(N-2);

    roots( p, prexr );

    if ( !interior )
      {
	xr[0]= -1.0;
	xr[N-1]= 1.0;
      }
    for(int i = 1-interior; i < N-(1+interior);++i )
      xr[i] = prexr[i-(1-interior)];


    const value_type two = 2.0;

    value_type a1 = std::pow( two,int(a+b+1.0));
    value_type a2 = fact( a+value_type( N ) -1.0 );//gamma(a + Q);
    value_type a3 = fact( b+value_type( N ) -1.0 );//gamma(b + Q);

    value_type a4 = fact( a+b+value_type( N ) );//gamma(a + b + m + 1);
    value_type a5 = fact(value_type(N)-1.0)*(value_type(N)-1.0); // (Q-1)!(Q-1)

    value_type a6 = a1 * a2 * a3;// Numérateur

    for ( int k = 1-interior; k < int(N-1-interior); ++k )
      {
        value_type fq = q.value( xr[k] );

        value_type dn = fq*fq*a4*a5;

        wr[k] =a6 /dn ;
      }
    if ( !interior )
      {
	wr[0] = wr[0]*(b+1.0);
	wr[N-1] = wr[N-1]*(a+1.0);
      }
  }
  /**
   * Integrate the function \p f using the quadrature \f$(w_q,
   * x_q)_q=0...N\f$ which integrates exactely polynomials of degree
   * \f$N\f$
   */
  template<typename T>
  T integrate( int N, boost::function<T( T const&)> const& f )
  {
    int quadrature_degree = (N-1)/2+2;
    std::cout << "integrate: quadrature degree is = " << (quadrature_degree-1) << "\n";
    std::cout << "integrate: integrates exactely polynomials of degree " << 2*(quadrature_degree-1)+1 << "\n";
    typedef T value_type;
    ublas::vector<T> xr( quadrature_degree );
    ublas::vector<T> wr( quadrature_degree );

    // get weights and nodes for Legendre polynomials
    gaussjacobi<T, ublas::vector<T> >( quadrature_degree, wr, xr );

    value_type res = 0.0;
    for ( int k = 0;k < quadrature_degree; ++k)
      res += wr[k]*f( xr[k] );
    return res;
  }
} // Life
#endif /* __Jacobi_H */

