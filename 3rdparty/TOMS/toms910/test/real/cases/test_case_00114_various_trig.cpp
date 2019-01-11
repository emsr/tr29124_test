
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00114_various_trig : public TestCaseReal
    {
    public:
      TestCase_case_00114_various_trig() { }
      virtual ~TestCase_case_00114_various_trig() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00114_various_trig");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.clear();
        data.push_back(ef::sec(ef::euler_gamma()));
        data.push_back(ef::sec(ef::catalan()));
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 2u> a =
        {{
           e_float("1.1933383729559773203137729267291947387173342208011247547153095890743890561595544691386633932476257421048110900061166081912943949195747181880197698938173812351027488623304591726956668351876851737682699217189957312877162495237230373808708191601490391769224697378179329863218317580732300824551143583689671403028593009944218611445507878903766448405383719847063127742468033101525589647537043486864948893731"),
           e_float("1.6419687273531364590567262071166440461685175366204953964739895651486200745161209040937844868647740537903996967626064986560711138637337304241821135325361595693541390424196494707138586320032467383542843973174977685986944349019525021804348293697747650763421720272451832097977353331515326976664205296758687867984902307858016265514207983161159140822122404062357504219549141106987355182424111888131203820022"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00114_various_trig(const bool b_write_output)
    {
      return TestCase_case_00114_various_trig().execute(b_write_output);
    }
  }
}
