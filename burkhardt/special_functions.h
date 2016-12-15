#ifndef BURKHARDT_SPECIAL_FUNCTIONS_H
#define BURKHARDT_SPECIAL_FUNCTIONS_H 1

#ifdef __cplusplus
extern "C" {
#endif

// airy_ai airy_bi
void airya_(const double * x, double * ai, double * bi, double * ad, double * bd);

// struve_h
void stvhv_(const double * nu, const double * x, double * sh);

// struve_l
void stvlv_(const double * nu, const double * x, double * sl);

// conf_hyperg
void cchg_(const double * a, const double * b, const double * z, double * chg);

// tricomi_u
void chgu_(const double * a, const double * b, const double * x, double * hu, int * md);

// hyperg
void hygfx_(const double * a, const double * b, const double * c, const double * x, double * hf);

#ifdef __cplusplus
}
#endif

#endif // BURKHARDT_SPECIAL_FUNCTIONS_H
