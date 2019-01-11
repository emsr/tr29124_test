#ifndef _MATH_LINK_BASE_2009_10_17_H_
  #define _MATH_LINK_BASE_2009_10_17_H_

  #include <string>

  class MathLinkBase
  {
  protected:

    MathLinkBase(char* prg_path = "MyProgramPath",
                 char* lnk_name = "MyMathLinkName")
    {
      static_cast<void>(prg_path);
      static_cast<void>(lnk_name);
    }

  public:

    virtual ~MathLinkBase() { }

    virtual int  ReturnPacketId(void) const                                = 0;
    virtual bool IsOpen        (void) const                                = 0;
    virtual bool CheckFunction (const std::string& str, long& count) const = 0;
    virtual bool EndPacket     (void) const                                = 0;
    virtual int  Error         (void) const                                = 0;
    virtual bool GetFunction   (std::string& str, int& count) const        = 0;
    virtual bool GetInteger    (int& i) const                              = 0;
    virtual bool GetString     (std::string& str) const                    = 0;
    virtual int  GetType       (void) const                                = 0;
    virtual int  GetNext       (void) const                                = 0;
    virtual int  NextPacket    (void) const                                = 0;
    virtual int  NewPacket     (void) const                                = 0;
    virtual bool PutFunction   (const std::string& str, int argc) const    = 0;
    virtual bool PutInteger    (const int n) const                         = 0;
    virtual bool PutString     (const std::string& str) const              = 0;
    virtual bool PutSymbol     (const std::string& str) const              = 0;
  };

#endif // _MATH_LINK_BASE_2009_10_17_H_
