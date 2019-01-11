#ifndef _EF_CAS_OBJECT_2010_01_15_H_
  #define _EF_CAS_OBJECT_2010_01_15_H_

  #include <string>
  #include <vector>

  #include <e_float/e_float.h>
  #include <functions/complex/e_float_complex.h>

  namespace ef_cas
  {
    class ComputerAlgebraSystemObject
    {
    private:

      ComputerAlgebraSystemObject(const ComputerAlgebraSystemObject&);
      const ComputerAlgebraSystemObject& operator=(const ComputerAlgebraSystemObject&);

    protected:

      const std::string my_object_path;
      const std::string my_object_name;

    protected:

      ComputerAlgebraSystemObject(const char* obj_path = "",
                                  const char* obj_name = "") : my_object_path(obj_path),
                                                               my_object_name(obj_name) { }

    public:

      INT32 my_digits(void) const { return static_cast<INT32>(ef::tol()); }

      virtual bool set_precision      (const UINT32 p) const { static_cast<void>(p); return false; }
      virtual bool set_precision_extra(const UINT32 p) const { static_cast<void>(p); return false; }

      virtual void create_mp_string(std::string& str, const e_float&    x) const { static_cast<void>(str.size()); static_cast<void>(x); }
      virtual void create_mp_string(std::string& str, const ef_complex& z) const { static_cast<void>(str.size()); static_cast<void>(z); }

      virtual bool get_values(const std::string& str_cmd, std::vector<e_float>&    values) const = 0;
      virtual bool get_values(const std::string& str_cmd, std::vector<ef_complex>& values) const = 0;
    };
  }

#endif // _EF_CAS_OBJECT_2010_01_15_H_
