
#include <LambertW.h>

#include <wrap_lambert.h>

namespace lambert
{

  template<>
    double
    lambert<0>(double x)
    { return utl::LambertW<0>(x); }

  template<>
    double
    lambert<-1>(double x)
    { return utl::LambertW<-1>(x); }

} // namespace lambert

