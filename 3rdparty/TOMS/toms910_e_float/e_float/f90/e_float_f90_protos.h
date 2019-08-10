
#ifndef _E_FLOAT_F90_PROTOS_2010_01_13_H_
  #define _E_FLOAT_F90_PROTOS_2010_01_13_H_

  // Wrap the function names of the F90 interface.
  extern "C"
  {
    void QXSETDBL  (ef_f90_t* DST, const double* X);
    void QXSETN32  (ef_f90_t* DST, const INT32* N);
    void QXSETN64  (ef_f90_t* DST, const INT64* N);
    void QXSETSTR  (ef_f90_t* DST, const char* s);
    void QXADD     (ef_f90_t* DST, const ef_f90_t* X);
    void QXSUB     (ef_f90_t* DST, const ef_f90_t* X);
    void QXMUL     (ef_f90_t* DST, const ef_f90_t* X);
    void QXDIV     (ef_f90_t* DST, const ef_f90_t* X);
    void QXMULN    (ef_f90_t* DST, const INT32* N);
    void QXDIVN    (ef_f90_t* DST, const INT32* N);
    void QXNEGATE  (ef_f90_t* DST);
    void QXINV     (ef_f90_t* DST);
    void QXSQRT    (ef_f90_t* DST);
    void QXISNEG   (const ef_f90_t* SRC, INT32* N);
    void QXISNAN   (const ef_f90_t* SRC, INT32* N);
    void QXISINF   (const ef_f90_t* SRC, INT32* N);
    void QXISFINITE(const ef_f90_t* SRC, INT32* N);
    void QXORDER2  (const ef_f90_t* SRC, INT32* N);
    void QXPARTS2  (const ef_f90_t* SRC, double* M, INT32* E);
    void QXHUGE    (ef_f90_t* DST, const ef_f90_t* X);
    void QXTINY    (ef_f90_t* DST, const ef_f90_t* X);
    void QXCMP     (const ef_f90_t* X, const ef_f90_t* Y, INT32* N);
    void QX2N64    (const ef_f90_t* SRC, INT64* N);
    void QX2DBL    (const ef_f90_t* SRC, double* X);
    void QX2STR    (const ef_f90_t* SRC, char* STR);
    void QXFLOOR   (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXCEIL    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXCBRT    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXROOTN   (ef_f90_t* DST, const ef_f90_t* X, const INT32* P);
    void QXEXP     (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXLOG     (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXPOW     (ef_f90_t* DST, const ef_f90_t* X, const ef_f90_t* A);
    void QXSIN     (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXCOS     (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXTAN     (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXASIN    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXACOS    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXATAN    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXATAN2   (ef_f90_t* DST, const ef_f90_t* Y, const ef_f90_t* X);
    void QXSINH    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXCOSH    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXTANH    (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXASINH   (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXACOSH   (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXATANH   (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXERF     (ef_f90_t* DST, const ef_f90_t* SRC);
    void QXERFC    (ef_f90_t* DST, const ef_f90_t* SRC);
  }

#endif // _E_FLOAT_F90_PROTOS_2010_01_13_H_
