
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00612_legendre_pvu_vary_12 : public TestCaseReal
    {
    public:
      TestCase_case_00612_legendre_pvu_vary_12() { }
      virtual ~TestCase_case_00612_legendre_pvu_vary_12() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00612_legendre_pvu_vary_12");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        static const e_float sqrt_1_3 = ef::sqrt(ef::third());
        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const e_float v  = -127 - sqrt_1_3;
          const e_float u  = +17  + sqrt_1_5;
          data[k] = ef::legendre_p(v, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-7.9750241822048537642199573394890705335211095097173794292564421430366910297132104182421264371362954583967827482072326357915608025151530683315739975442609754008922695657557372271590650896056791652002777650884524311901003073859288891911439369577141502569443480636983885541373830813469732302530055601628071991282579926191399971810368596879625163344103444690022450037565871481742513077006491837752743229207e36"),
           e_float("-1.8131557337424523016940940231451658828298571408810494945341283332278707592250124638451112155506550079872241833524857934019604910392770069377482096555534936305717945511790018556198569839573908394983107716656478224013068930287780608041572458192750986339012511377049021634950189625445002039821300684992829457769876235146613277587886504840713789703051725587306030934544308693441248798923517972087119354652e35"),
           e_float("-3.8975520340086221806741087817390030442175238685989334738489314097600343164074587622882802901313238472751912684383708659905010545308088492754942586177790081896290772118010457523092154451120434676120639468956491534500157113684528067402817969695692904799873569933054222299219870411804927673452006884656939499972264004946121929046340764591943891834860903000326510267074637943554443983607243195655311351282e34"),
           e_float("2.6317805238819575989755116994265165636369207675976166371732377370697453822622814084329232751815423136429943649539643206368967037813408646712388006384202931350967524387387824672263379871160408881826582428395057983423223405524216340649824521831819334011456095018139776254085220041533388645283649536625656429708822699037186032497199096109000324703366366472984721202109840194534661393493701479427054013901e35"),
           e_float("3.0726079638906047175968711494710842397941930191262856503388656686248349187847880086225895058594844515446063712265631569570624731760912353711810474839129823963193110541926697003045355606061717408717807748677774927124327767923026188952158413473103788520558202789026129014911909409581215198564112105200053501908156605304563792716427275811391542304981159212480298331445614546885896503775792495225194091114e35"),
           e_float("2.6674202332196598495101835991890774367133498912406670265917456682281116179248191938714619240334388476238364243524090855860746160011016137270878344606746387619364367139677300841902202136107013570760598384425488213549755372496593027750271378282320080628196762888046422773645727345102161749886730733114573221524081563023531077044601660656986782830525511222336084226837983283234765429377622667359499219466e35"),
           e_float("2.936384918494461208822839483248967208711415387437790109460171550242834334200797037025014178754500443531501327589155268850965826519450043754902530513700285790951817155053468371754948573978515101532260098920120817871504914561859872920804132858021190183047566524774098019322972215873751508237755215786820050511165224127331999104115528610144742887124751310362254380523690980100710643786387060692129970611e35"),
           e_float("2.8138754048821921261586151630846888062339034328690480820457136682043951554051904438458053844941596617125727822090446226084235597445262840985172389116175201590552775739380323400321390130408926124107240105449860084840002499504261015524668309836677998738598460388572482568120554370063291118308492106957651887106292891367791675415691148760697028288918685638175028477321461274149928939042343231601695550987e35"),
           e_float("-6.8466989522012223419021436420419551127399334930585533026637407046314510085147223982980246121537702549092892567912963797660116928363360070957534096709018116725028674808308718230416675572203353612485183466907852742099873158169615613865886269984010544187232305361749467315468975431928002466542013586376433541070114135633717898399514559744060026248466580394568683485043892380464008658236212360437499305624e34"),
           e_float("-1.4930615008852121322028766393988685786197324003954976243087271617421151207266147373657723039706949573937516814313460616270888500743952100013195890610998915029045636272775581135831653252257921561956485358841482484985131441427322470079731072119268733189131598829348041563376343718422862658824667016599232579113485999724308412630234260581462478537797335310587760921258458792296791481726338236561848105949e35"),
           e_float("-8.0968691307567293746289597373662152751697104698939939426129380832116232539821305869186729763649459007628792836656013292165155337231665713667321273948639254084437426272305421926726818841203159878435178166097459353784001422070132853549510414471473267105226209431782485075165495256311630897808382390894971518337601669925364254272300829029156400730386760239478565862362253108769247513471580655005112670351e36"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00612_legendre_pvu_vary_12(const bool b_write_output)
    {
      return TestCase_case_00612_legendre_pvu_vary_12().execute(b_write_output);
    }
  }
}