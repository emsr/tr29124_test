
#include <iostream>
#include <utility>


///
///  @brief Evaluate a continued fraction via the Lentz method.
///         The coefficients are to be provided by a function.
///
template<typename _Tp>
  inline _Tp
  continued_fraction_lentz( std::pair<_Tp,_Tp> Coeff(const unsigned int, const _Tp), const _Tp x )
  {
    const _Tp tiny = std::numeric_limits<_Tp>::epsilon();
    const _Tp tol = 0.00001;

    std::pair<_Tp, _Tp> c = Coeff( 0, x );

    _Tp f = c.second;
    if ( f == 0 )
      f = tiny;

    _Tp C = f;
    _Tp D = 0;
    for ( unsigned int i = 1; i < 100; ++i )
      {
	c = Coeff( i, x );

	D = c.second + c.first * D;
	if ( D == 0 )
	  D = tiny;

	C = c.second + c.first / C;
	if ( C == 0 )
	  C = tiny;

	D = 1 / D;

	auto R = C * D;
	f *= R;

	if ( std::abs(1 - R) < tol )
	  break;
      }

    return f;
  }


template<typename _Tp>
  class BesselTerm
  {
  public:
    BesselTerm(const _Tp n)
    { nu = n; }
    std::pair<_Tp,_Tp>
    operator()(const unsigned int n, const _Tp x) const
    { return std::make_pair( _Tp(1), _Tp(1) ); }
  protected:
    _Tp nu;
  };


int
main()
{


  return 0;
}


