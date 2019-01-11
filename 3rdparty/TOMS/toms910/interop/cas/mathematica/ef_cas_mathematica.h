#ifndef _EF_CAS_MATHEMATICA_2010_01_15_H_
  #define _EF_CAS_MATHEMATICA_2010_01_15_H_

  #include <interop/cas/mathematica/ef_cas_mathematica_base.h>

  namespace ef_cas
  {
    class Mathematica : public MathematicaBase
    {
    private:

              int   o_stat;
      mutable void* lnk_ptr;

    public:

      Mathematica(const char* obj_path = "",
                  const char* obj_name = "") : MathematicaBase(obj_path, obj_name) { }

      virtual ~Mathematica() { }

    private:

      static void Replace_StarHat_With_e(std::string& str);

      bool SendCommand(const std::string& str_cmd) const;

      virtual void create_mp_string(std::string& str, const e_float&    x) const;
      virtual void create_mp_string(std::string& str, const ef_complex& z) const;

      virtual bool get_values(const std::string& str_cmd, std::vector<e_float>&    values) const;
      virtual bool get_values(const std::string& str_cmd, std::vector<ef_complex>& values) const;
    };
  }

#endif // _EF_CAS_MATHEMATICA_2010_01_15_H_
