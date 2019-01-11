
#ifndef _TABLES_2008_01_10_H_
  #define _TABLES_2008_01_10_H_

  #include <vector>

  #include <e_float/e_float.h>

  namespace Tables
  {
    typedef const e_float&                           (*pfn_efloat)              (void);
    typedef const std::vector<e_float>&              (*pfn_vector_efloat)       (void);
    typedef const std::vector<std::vector<e_float> >&(*pfn_vector_vector_efloat)(void);

    const std::vector<pfn_efloat>&        A000142(void);
    const std::vector<pfn_efloat>&        A000364(void);
    const std::vector<pfn_efloat>&        A000367(void);
    const std::vector<pfn_efloat>&        A002445(void);
    const std::vector<pfn_efloat>&        A006882(void);
    const std::vector<pfn_vector_efloat>& A007318(void);
    const std::vector<pfn_vector_efloat>& A008277(void);
    const std::vector<pfn_vector_efloat>& A144617(void);
    const std::vector<pfn_efloat>&        A144618(void);
    const std::vector<pfn_vector_efloat>& A144622(void);
    const std::vector<pfn_vector_efloat>& A158503(void);
    const std::vector<pfn_efloat>&        A001164(void);
  }

#endif // _TABLES_2008_01_10_H_
