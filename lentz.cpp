
#include <iostream>
#include <utility>


///
///  @brief Evaluate a continued fraction via the Lentz method.
///         The coefficients are to be provided by a function.
///
template <class Type>
inline Type continued_fraction_lentz( std::pair<Type,Type> Coeff(const unsigned int, const Type), const Type x )
{

  const Type tiny = std::numeric_limits<Type>::epsilon();
  const Type tol = 0.00001;

  std::pair<Type,Type> c = Coeff( 0, x );

  Type f = c.second;
  if ( f == 0 )
    f = tiny;

  Type C = f;
  Type D = 0;
  for ( unsigned int i = 1; i < 100; ++i )
    {

      c = Coeff( i, x );

      D = c.second + c.first * D;
      if ( D == 0 )
        {
          D = tiny;
        }

      C = c.second + c.first / C;
      if ( C == 0 )
        {
          C = tiny;
        }

      D = 1 / D;

      Type R = C * D;
      f *= R;

      if ( std::abs(1 - R) < tol )
        break;
    }

  return f;
}


template <class Type>
class BesselTerm
{
public:
  BesselTerm(const Type n)
    {
      nu = n;
    }
  std::pair<Type,Type> operator()(const unsigned int n, const Type x) const
    {
      return std::make_pair( Type(1), Type(1) );
    }
protected:
  Type nu;
};


int main(int, char **)
{


  return 0;
}


