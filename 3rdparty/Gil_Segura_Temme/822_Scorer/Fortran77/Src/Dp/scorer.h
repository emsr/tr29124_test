#include <stddef.h>
#ifdef __cplusplus
#include <complex>
#define __GFORTRAN_FLOAT_COMPLEX std::complex<float>
#define __GFORTRAN_DOUBLE_COMPLEX std::complex<double>
#define __GFORTRAN_LONG_DOUBLE_COMPLEX std::complex<long double>
extern "C" {
#else
#define __GFORTRAN_FLOAT_COMPLEX float _Complex
#define __GFORTRAN_DOUBLE_COMPLEX double _Complex
#define __GFORTRAN_LONG_DOUBLE_COMPLEX long double _Complex
#endif

void scorer_gi (int ifacg, double x, double *y, double *reg, double *img, double *regp, double *imgp, int *ierrog);
void scorer_hi (int ifach, double x, double *y, double *reh, double *imh, double *rehp, double *imhp, int *ierroh);

#ifdef __cplusplus
}
#endif
