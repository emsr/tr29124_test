
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00604_legendre_pvu_vary_04 : public TestCaseReal
    {
    public:
      TestCase_case_00604_legendre_pvu_vary_04() { }
      virtual ~TestCase_case_00604_legendre_pvu_vary_04() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00604_legendre_pvu_vary_04");
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
          const e_float v  = -((10 - k) * 13) - sqrt_1_3;
          const e_float u  = +((10 - k) * 11) + sqrt_1_5;
          data[k] = ef::legendre_p(v, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("8.66297480271492325085532766021459678367687191114110747988815753060000985741520000087006224675039328110257160079811142712034191684925941772819912781927275941959601840028584414925731283837121851935341527789601921549299100372686445095255226645231046534349053599436672510870612112388646563926649497696021315601386160973047056304858665772222770189000608643234866211115264172217651890855564576358638920141e316"),
           e_float("-7.3011371614425896805038496593782045185403546071358473381123375470835672912793136818569518360904016297388811680529843716014683522217465560585199195242115666810909828927564351766406835096145501897638311804932502326947661867930053098490739749304850233973512787433525413879142026885668142131152831720169252983507524204713238837312130220682062061175607457798602691788272146044151533135717046870955655002121e207"),
           e_float("1.4561866635917319121832161084597428382957184424895594302822221176697923218133447528792970265291505675063703830209052124424408062182747524046947961174059162361753430085761988951374510620564075114835758246008860974668711144075817976651230039138794589121514837987922410902386910907955288612792219802821024612441194466891405398419400655890042274410997533144874279045569997652946152881921676218278110978216e172"),
           e_float("1.4112837057832067778019533396184505677778534899640640374984719996550581379974464354741770283001443216837770489967742068750403833000288354356561125142623412395068980334427433686919003152105400154360551269670138032767454921894205936819878568461844966894614582704541643969173048571439863664078999826752242831936346957089941773638146580042368012736704921632726864821810849947032850828310758116554036863101e145"),
           e_float("-8.9156807294914844500075305811676642528474765799880633907461704392540870329279529920629249045854105521881360073287567709226515376508485135530063999538102333256114062272939235160845929454177454910668566834393569261150594214673958042402305061054893551045191780187665029222173850802866172129613878478212718450686518830126843263137153091325511472350254681902508585394370192709356542399940063624332763276887e118"),
           e_float("5.3054587092931983545442251709652935352915524583968047950472265920154178526882419513798590877991223716058239059617824082858282328557357654366799893510217549920283274260138997355571468549437527850350777322158633428653666541559370987776161755118399177807349590192386868653666129356069958660219201243510148702640478792005856346729902878078096043239724945758282278225801635772257334985538005976954363817214e95"),
           e_float("9.3011536543998252496897332935611143338898408973129452531829821046824950045519686048680129608440888534354659010712699140038604904208091988048459966863684212960921573029235528957415985080448786388217011593558472602200538591144890423373498193795183160825226792653968063559974141646981477622557083437849383905969506538119446177140883084407312737881564211393793738098224875587155878848607183368316395613214e71"),
           e_float("5.4482706166790082530884705463793230250387325476397285507512023916566598448351811183774167433889399650904447546716374124903862064149689125367202171692821997295009966001159841609418289957789623016263407334378149822298922124614731107607827494395453798697997507364259351243297788009182971771773549231569265171808863440226107599962781866587245319789684521487459453324786278186829524894492951754579259505306e49"),
           e_float("8.301395164109077460365767380002679554451806906937079228297465505077061482184417253225900368580961597867858557176599401723907812476376075923833762590419991562143287011754645067812203885769776607691293164867804557416342644415321256663184671191309354520392179470792894790022993647355211994101005052437881677889014505033916304777752942006066954891226100351024000995343808532664765544004130947008877312988e29"),
           e_float("-5.0210291096404005087864108209306435570472495979006753103868323599987180468477752882631235064064330842125661051225746216534684349600474469194795115216926063502382510903686820698870698657101591842237301893615096259216424484674674222765582466071065386749512612404787518448240092140321199177885105697470250891946808321492911317750871910549606747733991408633277654952437551563061550043700075340832020871761e12"),
           e_float("2.3003011413959304992068583965883982510580883230383625787605986376243247640021414867445828733621975518877375559350820008168560907793987103286588746633298463797178402832708393332290679431556256818719211055999200244022288588545362779138961059307982384795991652386686442170751360874952813055595806759982361123255080159039029030619823978337957574016911080061251007205145993961056637558609322013774691341184"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00604_legendre_pvu_vary_04(const bool b_write_output)
    {
      return TestCase_case_00604_legendre_pvu_vary_04().execute(b_write_output);
    }
  }
}