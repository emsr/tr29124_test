
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00783_legendre_qnm_zero_x : public TestCaseReal
    {
    public:
      TestCase_case_00783_legendre_qnm_zero_x() { }
      virtual ~TestCase_case_00783_legendre_qnm_zero_x() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00783_legendre_qnm_zero_x");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const INT32 n = static_cast<INT32>(static_cast<INT32>(2 * k) + static_cast<INT32>(1));
          const INT32 m = static_cast<INT32>(k + 1);
          data[k] = ef::legendre_q(n * n, m * m, ef::zero());
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("0"),
           e_float("-3072."),
           e_float("0"),
           e_float("-1.78461872619598450262016e26"),
           e_float("0"),
           e_float("-7.5073288706651227102734047857230093505795933957042463990223044486662949519597376849168827207564434955009913070001525087692542321183468049412841238371206344364800976056123227085557419551624218392557572060393472624675918865334756748513039499771236846118651822479792588073814244319048345279853591581515937166387067256367241116364190940979106298612170199786487723044075034314473082202226628031111788927863e73"),
           e_float("0"),
           e_float("-1.4108492440794321911997080377799423116846203247004013418392556320863657791058250298911460058357872455428243547426229367737051349832260385385405664172477021625589587421153608242461165294536798946195202447780410905735937874662794638006558253261362217655533451875992790054646448547040873437262462903517580257636362040719148734989778738927674917186694497108576375510532244462103482023121763656587143434259e149"),
           e_float("0"),
           e_float("-1.1787257109757566740290169793853881106481575082659887530618922827883552503297162365877630192258495478036848576280815202994715558710686539850906078995691617070128871172268779824131132145236710495225392726912180158332520388885332158712645698523349845062351240071695965766012362610744186806709757117922894447862498099577192520727427497520763388242869548856762602635307007093439403677729549559241808260389e254"),
           e_float("0"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00783_legendre_qnm_zero_x(const bool b_write_output)
    {
      return TestCase_case_00783_legendre_qnm_zero_x().execute(b_write_output);
    }
  }
}
