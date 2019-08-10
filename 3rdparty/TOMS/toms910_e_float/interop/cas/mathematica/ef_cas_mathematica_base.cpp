#include <array>
#include <vector>
#include <algorithm>

#include <interop/cas/mathematica/7-0/mathlink.h>
#include <interop/cas/mathematica/ef_cas_mathematica_base.h>

namespace
{
  inline MLINK MLP(void* p) { return static_cast<MLINK>(p); }
}

void ef_cas::MathematicaBase::InitAndOpen(void)
{
  const MLENV ep = ::MLInitialize(static_cast<MLParametersPointer>(0));

  if(ep == static_cast<MLENV>(0))
  {
    o_stat = 1;
  }
  else
  {
    std::vector<char> v_path(my_object_path.length() + 1u);
    std::vector<char> v_name(my_object_name.length() + 1u);

    std::copy(my_object_path.begin(), my_object_path.end(), v_path.begin());
    std::copy(my_object_name.begin(), my_object_name.end(), v_name.begin());

    v_path.back() = static_cast<char>('\0');
    v_name.back() = static_cast<char>('\0');

    std::tr1::array<char*, 4u> args =
    {
      &v_path[0u],
      &v_name[0u],
      "-linkmode",
      "launch"
    };

    int err;
    const MLINK lp = ::MLOpenArgcArgv(ep,
                                      static_cast<int>(args.size()),
                                      args.data(),
                                      &err);

    if(lp == static_cast<MLINK>(0u))
    {
      ::MLDeinitialize(ep);
      o_stat = 2;
    }
    else
    {
      lnk_ptr = static_cast<void*>(lp);
    }
  }
}

int ef_cas::MathematicaBase::ReturnPacketId(void) const
{
  return RETURNPKT;
}

void ef_cas::MathematicaBase::Close(void) const
{
  ::MLClose(::MLP(lnk_ptr));
}

bool ef_cas::MathematicaBase::CheckFunction(const std::string& str, long& count) const
{
  long cnt = 0L;
  const bool b_check = ::MLCheckFunction(::MLP(lnk_ptr), str.c_str(), &cnt) != 0;
  count = cnt;
  return b_check;
}

bool ef_cas::MathematicaBase::EndPacket(void) const { return ::MLEndPacket(::MLP(lnk_ptr)) != 0; }
int  ef_cas::MathematicaBase::Error    (void) const { return ::MLError    (::MLP(lnk_ptr)); }

bool ef_cas::MathematicaBase::GetFunction(std::string& str, int& count) const
{
  const char* s = static_cast<const char*>(0u);
  int cnt = 0;
  const bool b_get = ::MLGetFunction(::MLP(lnk_ptr), &s, &cnt) != 0;
  str = b_get ? std::string(s) : "";
  ::MLDisownSymbol(::MLP(lnk_ptr), s);
  count = cnt;
  return b_get;
}

bool ef_cas::MathematicaBase::GetInteger(int& i) const
{
  int n = 0;
  const bool b_get = ::MLGetInteger(::MLP(lnk_ptr), &n) != 0;
  i = n;
  return b_get;
}

bool ef_cas::MathematicaBase::GetString(std::string& str) const
{
  const char* s = static_cast<const char*>(0u);
  const bool b_get = (::MLGetString(::MLP(lnk_ptr), &s) != 0);
  str = (b_get ? std::string(s) : "");
  ::MLDisownString(::MLP(lnk_ptr), s);
  return b_get;
}

int  ef_cas::MathematicaBase::GetType    (void) const                             { return ::MLGetType    (::MLP(lnk_ptr)); }
int  ef_cas::MathematicaBase::GetNext    (void) const                             { return ::MLGetNext    (::MLP(lnk_ptr)); }
int  ef_cas::MathematicaBase::NextPacket (void) const                             { return ::MLNextPacket (::MLP(lnk_ptr)); }
int  ef_cas::MathematicaBase::NewPacket  (void) const                             { return ::MLNewPacket  (::MLP(lnk_ptr)); }
bool ef_cas::MathematicaBase::PutFunction(const std::string& str, int argc) const { return ::MLPutFunction(::MLP(lnk_ptr), str.c_str(), argc) != 0; }
bool ef_cas::MathematicaBase::PutInteger (const int n) const                      { return ::MLPutInteger (::MLP(lnk_ptr), n) != 0; }
bool ef_cas::MathematicaBase::PutString  (const std::string& str) const           { return ::MLPutString  (::MLP(lnk_ptr), str.c_str()) != 0; }
bool ef_cas::MathematicaBase::PutSymbol  (const std::string& str) const           { return ::MLPutSymbol  (::MLP(lnk_ptr), str.c_str()) != 0; }
