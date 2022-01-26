#ifndef SF_COULOMB_H
#define SF_COULOMB_H 1

#include <emsr/detail/sf_coulomb.tcc>

namespace emsr
{
  template<typename Tp>
    Tp
    coulomb_norm(unsigned int l, Tp eta)
    {
      return emsr::detail::coulomb_norm(l, eta);
    }

} // namespace emsr

#endif // SF_COULOMB_H
