//  I want 24.6.5 [iterator.range] on steroids.
//
//  template<class C>
//    auto
//    begin(C& c) -> decltype(c.begin()) { return c.begin(); }
//  template<class C>
//    auto
//    begin(const C& c) -> decltype(c.begin()) { return c.begin(); }
//  template<class C>
//    auto
//    end(C& c) -> decltype(c.end()) { return   c.end(); }
//  template<class C>
//    auto
//    end(const C& c) -> decltype(c.end()) { return   c.end(); }
//
//  template<class T, size_t N>
//    T*
//    begin(T (&array)[N]) { return array; }
//  template<class T, size_t N>
//    T*
//    end(T (&array)[N]) { return array + N; }
//
//  You should be able to hit the result with [i][j].
//  It would be nice to get the number of rows and columns.
//
//  Algorithms like SVD would be templatized on iterators.
//  We require [][].

/*

 Starting with this:

template<class Numeric>
  bool
  lu_decomp(const unsigned long n, std::vector<std::vector<NUMERIC>> & a, std::vector<int> & index, double & parity)

 we would want this:

template<class Matrix>
  bool
  lu_decomp(const unsigned long n, Matrix & a, std::vector<int > & index, Matrix::value_type & parity)

 where the numeric type might be Matrix::value_type

*/

#include <cstdint>
#include <valarray>
#include <vector>
#include <array>

namespace emsr
{

template<typename Tp>
  class matrix
  {
    //operator[](size_t);
  };

template<typename Tp, size_t M, size_t N>
  class matrix<Tp (&)[M][N]>
  {
    matrix(Tp (&mat)[M][N]);
    matrix<Tp[N]> operator[](size_t);
    Tp (&m_mat)[M][N];
  };

template<typename Tp>
  class matrix<Tp*>
  {
    Tp*& m_mat;
  };

template<typename Tp>
  class matrix<Tp**>
  {
    Tp **& m_mat;
  };

template<typename Tp, size_t M>
  class matrix<Tp(*)[M]>
  {
    Tp (*&m_mat)[M];
  };

template<typename Tp, size_t M, size_t N>
  class matrix<std::array<Tp, M * N>>
  {
    std::array<Tp, M * N> & m_mat;
  };

template<typename Tp, size_t M, size_t N>
  class matrix<std::array<std::array<Tp, N>, M>>
  {
    std::array<std::array<Tp, N>, M> & m_mat;
  };

template<typename Tp>
  class matrix<std::vector<std::vector<Tp>>>
  {
    std::vector<std::vector<Tp>> & m_mat;
  };

template<typename Tp>
  class matrix<std::valarray<Tp>>
  {
    std::valarray<Tp> & m_mat;
  };



///
///
///
template<typename Tp>
  class banded_diagonal_matrix
  {
    size_t start_col;
    size_t width;
    std::vector<std::vector<Tp>> & stuff;
  };

}  //  namespace emsr

#include "matrix_lu_decomp.h"
#include "matrix_qr_decomp.h"
#include "matrix_sv_decomp.h"
#include "matrix_cholesky_decomp.h"
#include "matrix_gauss_jordan.h"
#include "matrix_tridiag.h"
#include "matrix_vandermonde.h"
