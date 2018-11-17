/*
/home/ESmith-rowland/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o complex_tie complex_tie.cpp 
./complex_tie
*/
#include <iostream>
#include <functional>
#include <tuple>
#include <complex>
/*
template<typename _Tp>
  std::tuple<_Tp, _Tp>
  complex_tie(std::complex<_Tp>& __z)
  {
    auto __carr = static_cast<_Tp[2]>(__z);
    return std::make_tuple(std::ref(__carr[0]), std::ref(__carr[1]));
  }
*/
/*
template<typename _Tp>
  std::complex<_Tp&>
  complex_tie(_Tp& __rez, _Tp& __imz)
  {
    //auto __carr = static_cast<_Tp[2]>(__z);
    //return std::make_tuple();
    return std::complex<_Tp&>(std::ref(__rez), std::ref(__imz));
  }
*/

template<typename _Tp>
  struct _ComplexTie
  {
    _Tp& _M_rez;
    _Tp& _M_imz;

    _ComplexTie(_Tp& __rez, _Tp& __imz)
    : _M_rez(__rez),
      _M_imz(__imz)
    { }

    explicit _ComplexTie(std::complex<_Tp>& __z)
    : _M_rez(reinterpret_cast<_Tp(&)[2]>(__z)[0]),
      _M_imz(reinterpret_cast<_Tp(&)[2]>(__z)[1])
    { }

    _ComplexTie&
    operator=(std::complex<_Tp>& __z)
    {
      //auto __zrep = __z.__rep();
      _Tp (&__zarr)[2] = reinterpret_cast<_Tp(&)[2]>(__z);
      this->_M_rez = __zarr[0];
      this->_M_imz = __zarr[1];
      return *this;
    }
  };

template<typename _Tp>
  _ComplexTie<_Tp>
  complex_tie(_Tp& __rez, _Tp& __imz)
  { return _ComplexTie<_Tp>(__rez, __imz); }

template<typename _Tp>
  _ComplexTie<_Tp>
  complex_tie(std::complex<_Tp>& __z)
  { return _ComplexTie<_Tp>(__z); }


//  Idea: complex_tie returns a impldef class that unpacks a complex into a pair of refs.
//  Overload operator=, construction.
//  This class should be decomposeable.
//  We may have to friend this in complex. Do we have public access to _M_rep?

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
