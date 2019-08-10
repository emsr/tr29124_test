
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00011_various_elem_math : public TestCaseReal
    {
    public:
      TestCase_case_00011_various_elem_math() { }
      virtual ~TestCase_case_00011_various_elem_math() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00011_various_elem_math");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.clear();
        data.push_back(ef::floor(ef::pi()));
        data.push_back(ef::ceil (ef::pi()));
        data.push_back(ef::floor(-100 - ef::euler_gamma()));
        data.push_back(ef::ceil (-100 - ef::euler_gamma()));
        data.push_back(e_float(ef::to_int32(e_float("1e9"))));
        data.push_back(e_float(ef::to_int64(e_float("1e18"))));
        data.push_back(e_float(ef::to_int32(e_float("1e29"))));
        data.push_back(e_float(ef::to_int64(e_float("1e29"))));
        data.push_back(e_float(ef::to_int32(ef_complex(ef::pi(), ef::euler_gamma()))));
        data.push_back(e_float(ef::to_int64(ef_complex(ef::pi(), ef::euler_gamma()))));
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 10u> a =
        {{
           e_float("3"),
           e_float("4"),
           e_float("-101"),
           e_float("-100"),
           e_float("1000000000"),
           e_float("1000000000000000000"),
           e_float("2147483647"),
           e_float("9223372036854775807"),
           e_float("3"),
           e_float("3"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00011_various_elem_math(const bool b_write_output)
    {
      return TestCase_case_00011_various_elem_math().execute(b_write_output);
    }
  }
}
