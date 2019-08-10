#include <interop/cas/mathematica/ef_cas_mathematica.h>
#include <interop/cas/mathematica/ef_cas_object_mathematica.h>

const ef_cas::ComputerAlgebraSystemObject& ef_cas::cas_mathematica(char* obj_path,
                                                                   char* obj_name)
{
  static const ef_cas::Mathematica m(obj_path, obj_name);
  return m;
}
