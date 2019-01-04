
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00773_legendre_pnm_zero_x : public TestCaseReal
    {
    public:
      TestCase_case_00773_legendre_pnm_zero_x() { }
      virtual ~TestCase_case_00773_legendre_pnm_zero_x() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00773_legendre_pnm_zero_x");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const INT32 n = static_cast<INT32>(static_cast<INT32>(2 * k) + static_cast<INT32>(1));
          const INT32 m = static_cast<INT32>(k + 1);
          data[k] = ef::legendre_p(n * n, m * m, ef::zero());
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-1."),
           e_float("0"),
           e_float("-6.13515689984310150146484375e11"),
           e_float("0"),
           e_float("-3.637442993320801013516632292440245261074038220730503530519865662551382001765887252986431121826171875e46"),
           e_float("0"),
           e_float("-5.2727073829369430063489537621654234910701204960546821166041712532966465833466465078110049274360598264277612587881725405069521508283428589504540978824523381847349205341598002657848404638174510949966133921407163143157958984375e107"),
           e_float("0"),
           e_float("-4.0118221411403337104500294038139065008364436934471179510446378343438585008574159883095971025509450702758202816309380001763322782145830173190430758241948531124587805835584261356171406452763957565982350363111466780113922609453810389730914217385045487751887863865303595082956250545985581055772756386844708249810218506766617057904164900476688995489816636704801775014164744170841458981158211827278137207031e197"),
           e_float("0"),
           e_float("-8.9402633922462602768337989949671156793995227216334274026735352772120225752217217542626791161391948259453704678318499382823243232815007394388755731675107743076437740027683359313604497145916354787625472754988793607645391160851271677469696670135642537482431508990618025242965285774617427766927811149127495826079187444107141598958997959104924233309149425848758454399618091281513158613926637486042308805623e317"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00773_legendre_pnm_zero_x(const bool b_write_output)
    {
      return TestCase_case_00773_legendre_pnm_zero_x().execute(b_write_output);
    }
  }
}
