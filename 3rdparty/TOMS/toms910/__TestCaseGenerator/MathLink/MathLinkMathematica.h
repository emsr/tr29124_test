#ifndef _MATH_LINK_MATHEMATICA_2009_10_17_H_
  #define _MATH_LINK_MATHEMATICA_2009_10_17_H_

  #include <MathLink/MathLinkBase.h>

  class MathLinkMathematica : public MathLinkBase
  {
  private:

            int   o_stat;
    mutable void* lnk_ptr;

  public:

    MathLinkMathematica(char* prg_path = "MyProgramPath",
                        char* lnk_name = "MyMathLinkName") : MathLinkBase(prg_path, lnk_name),
                                                             o_stat (0),
                                                             lnk_ptr(static_cast<void*>(0u))
    {
      InitAndOpen(prg_path, lnk_name);
    }

    virtual ~MathLinkMathematica()
    {
      if(o_stat > 1)
      {
        static_cast<void>(PutFunction("Exit", 0));
        Close();
      }
    }

  private:

    void InitAndOpen(char* prg_path, char* lnk_name);
    void Close(void) const;
    int  GetOpenStatus(void) const { return o_stat; }

  public:

    virtual bool IsOpen(void) const { return o_stat == 0; }

    virtual int  ReturnPacketId(void) const;
    virtual bool CheckFunction (const std::string& str, long& count) const;
    virtual bool EndPacket     (void) const;
    virtual int  Error         (void) const;
    virtual bool GetFunction   (std::string& str, int& count) const;
    virtual bool GetInteger    (int& i) const;
    virtual bool GetString     (std::string& str) const;
    virtual int  GetType       (void) const;
    virtual int  GetNext       (void) const;
    virtual int  NextPacket    (void) const;
    virtual int  NewPacket     (void) const;
    virtual bool PutFunction   (const std::string& str, int argc) const;
    virtual bool PutInteger    (const int n) const;
    virtual bool PutString     (const std::string& str) const;
    virtual bool PutSymbol     (const std::string& str) const;
  };

#endif // _MATH_LINK_MATHEMATICA_2009_10_17_H_
