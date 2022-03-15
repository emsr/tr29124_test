/**
 *
 */

#include <iostream>
#include <functional>
#include <tuple>
#include <complex>
/*
template<typename Tp>
  std::tuple<Tp, Tp>
  complex_tie(std::complex<Tp>& z)
  {
    auto carr = static_cast<Tp[2]>(z);
    return std::make_tuple(std::ref(carr[0]), std::ref(carr[1]));
  }
*/
/*
template<typename Tp>
  std::complex<Tp&>
  complex_tie(Tp& rez, Tp& imz)
  {
    //auto carr = static_cast<Tp[2]>(z);
    //return std::make_tuple();
    return std::complex<Tp&>(std::ref(rez), std::ref(imz));
  }
*/

template<typename Tp>
  struct _ComplexTie
  {
    Tp& m_rez;
    Tp& m_imz;

    _ComplexTie(Tp& rez, Tp& imz)
    : m_rez(rez),
      m_imz(imz)
    { }

    explicit _ComplexTie(std::complex<Tp>& z)
    : m_rez(reinterpret_cast<Tp(&)[2]>(z)[0]),
      m_imz(reinterpret_cast<Tp(&)[2]>(z)[1])
    { }

    _ComplexTie&
    operator=(std::complex<Tp>& z)
    {
      //auto zrep = z.rep();
      Tp (&zarr)[2] = reinterpret_cast<Tp(&)[2]>(z);
      this->m_rez = zarr[0];
      this->m_imz = zarr[1];
      return *this;
    }
  };

template<typename Tp>
  _ComplexTie<Tp>
  complex_tie(Tp& rez, Tp& imz)
  { return _ComplexTie<Tp>(rez, imz); }

template<typename Tp>
  _ComplexTie<Tp>
  complex_tie(std::complex<Tp>& z)
  { return _ComplexTie<Tp>(z); }


//  Idea: complex_tie returns a impldef class that unpacks a complex into a pair of refs.
//  Overload operator=, construction.
//  This class should be decomposeable.
//  We may have to friend this in complex. Do we have public access to m_rep?

int
main()
{
  std::complex<double> z(3.0, 5.0);

  double rez, imz;
  complex_tie(rez, imz) = z;
  std::cout << rez << '\n';
  std::cout << imz << '\n';
  rez = 1.5;
  imz = -6.0;
  std::cout << z << '\n';
  auto [rez2, imz2] = complex_tie(z);
  std::cout << rez2 << '\n';
  std::cout << imz2 << '\n';
}
