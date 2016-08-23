#ifndef BURKHARDT_SPECIAL_FUNCTIONS_H
#define BURKHARDT_SPECIAL_FUNCTIONS_H 1

#ifdef __cplusplus
extern "C" {
#endif

void airya_(const double * x, double * ai, double * bi, double * ad, double * bd);

void stvhv_(const double * nu, const double * x, double * sh);

void stvlv_(const double * nu, const double * x, double * sl);

#ifdef __cplusplus
}
#endif

#endif // BURKHARDT_SPECIAL_FUNCTIONS_H
