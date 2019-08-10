#include <array>

#include <MathLink/7-0/mathlink.h>
#include <MathLink/MathLinkMathematica.h>

namespace
{
  inline MLINK MLP(void* p) { return static_cast<MLINK>(p); }
}

void MathLinkMathematica::InitAndOpen(char* prg_path, char* lnk_name)
{
  const MLENV ep = ::MLInitialize(static_cast<MLParametersPointer>(0));

  if(ep == static_cast<MLENV>(0))
  {
    o_stat = 1;
  }
  else
  {
    std::tr1::array<char*, 4u> args =
    {
      prg_path,
      lnk_name,
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

int MathLinkMathematica::ReturnPacketId(void) const
{
  return RETURNPKT;
}

void MathLinkMathematica::Close(void) const
{
  ::MLClose(::MLP(lnk_ptr));
}

bool MathLinkMathematica::CheckFunction(const std::string& str, long& count) const
{
  const char* s = static_cast<const char*>(0u);
  long cnt = 0L;
  const bool b_check = ::MLCheckFunction(::MLP(lnk_ptr), str.c_str(), &cnt) != 0;
  count = cnt;
  return b_check;
}

bool MathLinkMathematica::EndPacket(void) const { return ::MLEndPacket(::MLP(lnk_ptr)) != 0; }
int  MathLinkMathematica::Error    (void) const { return ::MLError    (::MLP(lnk_ptr)); }

bool MathLinkMathematica::GetFunction(std::string& str, int& count) const
{
  const char* s = static_cast<const char*>(0u);
  int cnt = 0;
  const bool b_get = ::MLGetFunction(::MLP(lnk_ptr), &s, &cnt) != 0;
  str = b_get ? std::string(s) : "";
  ::MLDisownSymbol(::MLP(lnk_ptr), s);
  count = cnt;
  return b_get;
}

bool MathLinkMathematica::GetInteger(int& i) const
{
  int n = 0;
  const bool b_get = ::MLGetInteger(::MLP(lnk_ptr), &n) != 0;
  i = n;
  return b_get;
}

bool MathLinkMathematica::GetString(std::string& str) const
{
  const char* s = static_cast<const char*>(0u);
  const bool b_get = (::MLGetString(::MLP(lnk_ptr), &s) != 0);
  str = (b_get ? std::string(s) : "");
  ::MLDisownString(::MLP(lnk_ptr), s);
  return b_get;
}

int  MathLinkMathematica::GetType    (void) const                             { return ::MLGetType    (::MLP(lnk_ptr)); }
int  MathLinkMathematica::GetNext    (void) const                             { return ::MLGetNext    (::MLP(lnk_ptr)); }
int  MathLinkMathematica::NextPacket (void) const                             { return ::MLNextPacket (::MLP(lnk_ptr)); }
int  MathLinkMathematica::NewPacket  (void) const                             { return ::MLNewPacket  (::MLP(lnk_ptr)); }
bool MathLinkMathematica::PutFunction(const std::string& str, int argc) const { return ::MLPutFunction(::MLP(lnk_ptr), str.c_str(), argc) != 0; }
bool MathLinkMathematica::PutInteger (const int n) const                      { return ::MLPutInteger (::MLP(lnk_ptr), n) != 0; }
bool MathLinkMathematica::PutString  (const std::string& str) const           { return ::MLPutString  (::MLP(lnk_ptr), str.c_str()) != 0; }
bool MathLinkMathematica::PutSymbol  (const std::string& str) const           { return ::MLPutSymbol  (::MLP(lnk_ptr), str.c_str()) != 0; }
