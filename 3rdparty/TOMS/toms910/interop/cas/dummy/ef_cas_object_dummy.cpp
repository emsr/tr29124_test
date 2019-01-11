#include <interop/cas/dummy/ef_cas_dummy.h>
#include <interop/cas/dummy/ef_cas_object_dummy.h>

const ef_cas::ComputerAlgebraSystemObject& ef_cas::cas_dummy(void)
{
  static const ef_cas::Dummy d;
  return d;
}
