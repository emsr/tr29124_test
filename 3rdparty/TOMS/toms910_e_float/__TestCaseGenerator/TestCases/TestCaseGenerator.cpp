#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <MathLink/MathLinkMathematica.h>
#include <TableEvaluator/MathLinkTableEvaluator.h>
#include <TableEvaluator/MathLinkTableEvaluatorComplex.h>
#include <TestCases/TestCaseGenerator.h>

namespace
{
  static const MathLinkMathematica& TheML(void)
  {
    static const MathLinkMathematica the_ml;
    return the_ml;
  }
}

namespace
{
  static const std::string& MarkerOfCppSequence(void)
  {
    static const std::string str("MARKER_OF_CPP_SEQUENCE");
    return str;
  }
  
  static const std::string& MarkerOfDataSequence(void)
  {
    static const std::string str("MARKER_OF_DATA_SEQUENCE");
    return str;
  }

  void TestCaseFrameWork(std::vector<std::string>& frame,
                         const std::string& str_case_name,
                         const TestCaseDescription::t_real_complex rc)
  {
    const bool b_is_complex = (rc == TestCaseDescription::is_complex);

    const std::string str_class_name   = "TestCase_" + str_case_name;
    const std::string str_derived_from = (!b_is_complex ? std::string("Real") : std::string("Imag"));
    const std::string str_include_path = (!b_is_complex ? std::string("#include <test/real/test_case_real.h>")
                                                        : std::string("#include <test/imag/test_case_imag.h>"));

    const std::tr1::array<std::string, 41u> frame_data =
    {
      std::string(""),
      std::string("// Automatically generated file"),
      std::string("#include <functions/functions.h>"),
      str_include_path,
      std::string(""),
      std::string("namespace test"),
      std::string("{"),
      std::string("  namespace " + (!b_is_complex ? std::string("real") : std::string("imag"))),
      std::string("  {"),
      std::string("    class " + str_class_name + " : public TestCase" + str_derived_from),
      std::string("    {"),
      std::string("    public:"),
      std::string("      " + str_class_name + "() { }"),
      std::string("      virtual ~" + str_class_name + "() { }"),
      std::string("    private:"),
      std::string("      virtual const std::string& name(void) const"),
      std::string("      {"),
      std::string("        static const std::string str(\"" + str_class_name + "\");"),
      std::string("        return str;"),
      std::string("      }"),
      std::string("      virtual void e_float_test(std::vector<" + (!b_is_complex ? std::string("e_float") : std::string("ef_complex")) + ">& data) const"),
      std::string("      {"),
      MarkerOfCppSequence(),
      std::string("      }"),
      std::string("      virtual const std::vector<" + (!b_is_complex ? std::string("e_float") : std::string("ef_complex")) + ">& control_data(void) const"),
      std::string("      {"),
      MarkerOfDataSequence(),
      std::string("      }"),
      std::string("    };"),
      std::string(""),
      std::string("    bool test_" + str_case_name + "(const bool b_write_output)"),
      std::string("    {"),
      std::string("      return " + str_class_name + "().execute(b_write_output);"),
      std::string("    }"),
      std::string("  }"),
      std::string("}")
    };

    frame = std::vector<std::string>(frame_data.begin(), frame_data.end());
  }

  std::string ToString(const std::size_t sz)
  {
    std::stringstream ss;
    ss << sz;
    std::string str;
    ss >> str;
    return str;
  }

  static void TransformStringsToCppDataCode(std::vector<std::string>& dst, const std::vector<std::vector<std::string> >& src)
  {
    dst.push_back("        static const std::tr1::array<e_float, " + ::ToString(src.size()) + "u> a =");
    dst.push_back("        {{");

    for(std::vector<std::string>::size_type i = 0u; i < src.size(); i++)
    {
      dst.push_back("           e_float(\"" + (src[i].at(0u) + "\"),"));
    }

    dst.push_back("        }};");
    dst.push_back("        static const std::vector<e_float> v(a.begin(), a.end());");
    dst.push_back("        return v;");
  }

  static void TransformStringsToCppDataCodeComplex(std::vector<std::string>& dst, const std::vector<std::vector<std::string> >& src)
  {
    dst.push_back("        static const std::tr1::array<ef_complex, " + ::ToString(src.size()) + "u> a =");
    dst.push_back("        {{");

    for(std::vector<std::string>::size_type i = 0u; i < src.size(); i++)
    {
      const std::string str_real = "e_float(\"" + (src[i].at(0u) + "\")");
      const std::string str_imag = "e_float(\"" + (src[i].at(1u) + "\")");
      dst.push_back("           ef_complex(" + str_real + ", " + str_imag + "),");
    }

    dst.push_back("        }};");
    dst.push_back("        static const std::vector<ef_complex> v(a.begin(), a.end());");
    dst.push_back("        return v;");
  }
}

bool TestCaseGenerator(const TestCaseDescription& desc)
{
  std::vector<std::string> test_case_code;
  
  ::TestCaseFrameWork(test_case_code, desc.case_name, desc.real_complex);

  const std::vector<std::string>::const_iterator it_cpp_seq = std::find(test_case_code.begin(),
                                                                        test_case_code.end(),
                                                                        ::MarkerOfCppSequence());

  if(it_cpp_seq != test_case_code.end())
  {
    test_case_code.insert(it_cpp_seq, desc.cpp_seq.begin(), desc.cpp_seq.end());
    test_case_code.erase(std::find(test_case_code.begin(),
                                   test_case_code.end(),
                                   ::MarkerOfCppSequence()));
  }

  const bool b_is_complex = (desc.real_complex == TestCaseDescription::is_complex);

  const MathLinkTableEvaluatorBase* p_eval = !b_is_complex ? static_cast<const MathLinkTableEvaluatorBase*>(new MathLinkTableEvaluator       (TheML(), desc.math_cmd))
                                                           : static_cast<const MathLinkTableEvaluatorBase*>(new MathLinkTableEvaluatorComplex(TheML(), desc.math_cmd));

  if(!p_eval->Evaluate())
  {
    delete p_eval;
    return false;
  }
  else
  {
    std::cout << "MathLink eval succeeded: " << desc.case_name << std::endl;
  }

  std::vector<std::string> cpp_data_strings;

  !b_is_complex ? ::TransformStringsToCppDataCode       (cpp_data_strings, p_eval->StringValues())
                : ::TransformStringsToCppDataCodeComplex(cpp_data_strings, p_eval->StringValues());

  delete p_eval;

  const std::vector<std::string>::const_iterator it_data_seq = std::find(test_case_code.begin(),
                                                                         test_case_code.end(),
                                                                         ::MarkerOfDataSequence());

  if(it_data_seq != test_case_code.end())
  {
    test_case_code.insert(it_data_seq, cpp_data_strings.begin(), cpp_data_strings.end());
    test_case_code.erase(std::find(test_case_code.begin(),
                                   test_case_code.end(),
                                   ::MarkerOfDataSequence()));
  }

  const std::string str_outfile = ("test_" + desc.case_name) + ".cpp";

  std::ofstream out(str_outfile.c_str());

  if(!out.is_open())
  {
    return false;
  }

  std::copy(test_case_code.begin(),
            test_case_code.end(),
            std::ostream_iterator<std::string>(out, "\n"));

  out.close();

  std::cout << "File creation succeeded: " << str_outfile << std::endl;

  return true;
}
