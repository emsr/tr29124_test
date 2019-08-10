#ifndef _TEST_CASE_DESCRIPTION_2009_10_28_H_
  #define _TEST_CASE_DESCRIPTION_2009_10_28_H_

  #include <string>
  #include <vector>
  #include <array>

  struct TestCaseDescription
  {
  public:

    typedef enum enum_real_complex
    {
      is_real,
      is_complex
    }
    t_real_complex;

    const std::vector<std::string> cpp_seq;
    const std::string              case_name;
    const std::string              math_cmd;
    const t_real_complex           real_complex;

    template<typename FwdIt> TestCaseDescription(FwdIt seq_begin,
                                                 FwdIt seq_end,
                                                 const std::string& cn,
                                                 const std::string& mc,
                                                 const t_real_complex rc = is_real) : cpp_seq     (seq_begin, seq_end),
                                                                                      case_name   (cn),
                                                                                      math_cmd    (mc),
                                                                                      real_complex(rc) { }
  };

#endif // _TEST_CASE_DESCRIPTION_2009_10_28_H_
