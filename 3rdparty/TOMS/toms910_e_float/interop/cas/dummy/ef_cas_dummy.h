#ifndef _EF_CAS_DUMMY_2010_01_15_H_
  #define _EF_CAS_DUMMY_2010_01_15_H_

  #include <interop/cas/ef_cas_object.h>

  namespace ef_cas
  {
    class Dummy : public ComputerAlgebraSystemObject
    {
    public:

      Dummy() { }

      virtual ~Dummy() { }

    public:

      virtual bool get_values(const std::string& str_cmd, std::vector<e_float>&    values) const { static_cast<void>(str_cmd.size()); static_cast<void>(values.size()); return false; }
      virtual bool get_values(const std::string& str_cmd, std::vector<ef_complex>& values) const { static_cast<void>(str_cmd.size()); static_cast<void>(values.size()); return false; }
    };
  }

#endif // _EF_CAS_DUMMY_2010_01_15_H_
