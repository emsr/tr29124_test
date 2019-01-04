#ifndef _EF_CAS_MATHEMATICA_BASE_2010_01_15_H_
  #define _EF_CAS_MATHEMATICA_BASE_2010_01_15_H_

  #include <interop/cas/ef_cas_object.h>

  namespace ef_cas
  {
    class MathematicaBase : public ComputerAlgebraSystemObject
    {
    private:

      int o_stat;
      mutable void* lnk_ptr;

    protected:

      MathematicaBase(const char* obj_path = "",
                      const char* obj_name = "") : ComputerAlgebraSystemObject(obj_path, obj_name),
                                                   o_stat (0),
                                                   lnk_ptr(static_cast<void*>(0u))
      {
        InitAndOpen();
      }

    public:

      virtual ~MathematicaBase()
      {
        if(o_stat > 1)
        {
          static_cast<void>(PutFunction("Exit", 0));
          Close();
        }
      }

    private:

      virtual bool get_values(const std::string& str_cmd, std::vector<e_float>&    values) const { static_cast<void>(str_cmd.size()); static_cast<void>(values.size()); return false; }
      virtual bool get_values(const std::string& str_cmd, std::vector<ef_complex>& values) const { static_cast<void>(str_cmd.size()); static_cast<void>(values.size()); return false; }

    protected:

      void InitAndOpen(void);
      void Close(void) const;
      int  GetOpenStatus(void) const { return o_stat; }

      bool IsOpen(void) const { return o_stat == 0; }

      int  ReturnPacketId(void) const;
      bool CheckFunction (const std::string& str, long& count) const;
      bool EndPacket     (void) const;
      int  Error         (void) const;
      bool GetFunction   (std::string& str, int& count) const;
      bool GetInteger    (int& i) const;
      bool GetString     (std::string& str) const;
      int  GetType       (void) const;
      int  GetNext       (void) const;
      int  NextPacket    (void) const;
      int  NewPacket     (void) const;
      bool PutFunction   (const std::string& str, int argc) const;
      bool PutInteger    (const int n) const;
      bool PutString     (const std::string& str) const;
      bool PutSymbol     (const std::string& str) const;
    };
  }

#endif // _EF_CAS_MATHEMATICA_BASE_2010_01_15_H_
