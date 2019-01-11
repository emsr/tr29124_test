
#ifndef _POLYNOMIALS_2007_02_16_H_
  #define _POLYNOMIALS_2007_02_16_H_

  #include <vector>

  namespace ef
  {
    e_float chebyshev_t(const INT32 n, const e_float& x);
    e_float chebyshev_u(const INT32 n, const e_float& x);
    e_float hermite    (const INT32 n, const e_float& x);
    e_float laguerre   (const INT32 n, const e_float& x);
    e_float legendre_p (const INT32 n, const e_float& x);
    e_float legendre_q (const INT32 n, const e_float& x);

    // NOCOVER_BLK_BEG
    e_float chebyshev_t(const UINT32 n, const e_float& x, std::vector<e_float>* vp);
    e_float chebyshev_u(const UINT32 n, const e_float& x, std::vector<e_float>* vp);
    e_float hermite    (const UINT32 n, const e_float& x, std::vector<e_float>* vp);
    e_float laguerre   (const UINT32 n, const e_float& x, std::vector<e_float>* vp);
    // NOCOVER_BLK_END
  }

  namespace efz
  {
    ef_complex chebyshev_t(const INT32 n, const ef_complex& z);
    ef_complex chebyshev_u(const INT32 n, const ef_complex& z);
    ef_complex hermite    (const INT32 n, const ef_complex& z);
    ef_complex laguerre   (const INT32 n, const ef_complex& z);

    // NOCOVER_BLK_BEG
    ef_complex chebyshev_t(const UINT32 n, const ef_complex& z, std::vector<ef_complex>* vp);
    ef_complex chebyshev_u(const UINT32 n, const ef_complex& z, std::vector<ef_complex>* vp);
    ef_complex hermite    (const UINT32 n, const ef_complex& z, std::vector<ef_complex>* vp);
    ef_complex laguerre   (const UINT32 n, const ef_complex& z, std::vector<ef_complex>* vp);
    // NOCOVER_BLK_END
  }

#endif // _POLYNOMIALS_2007_02_16_H_
