
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00802_conf_hyperg_vary2_mpx : public TestCaseReal
    {
    public:
      TestCase_case_00802_conf_hyperg_vary2_mpx() { }
      virtual ~TestCase_case_00802_conf_hyperg_vary2_mpx() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00802_conf_hyperg_vary2_mpx");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float ak = -ef::third()   - static_cast<INT32>(100 - (10 * k));
          const e_float bk = -ef::quarter() - static_cast<INT32>(0   + (10 * k));
          data[k] = ef::conf_hyperg(ak, bk, ef::euler_gamma() + (100 * k));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("12.726486343212659795175114085473376941333730696323478848523784039281599878210535545770015387782028137431386017036234204611278202301744835885553413635139456942985513858850091202852696071509364402814014106182351457630264055484894494452472453282973825919879207315264156297691587550287030067544769982334277676758369826302446995509249422546140322273847622651399127413292103782079122621187023051958135253694"),
           e_float("1.7102906846150624071727907644902833279101181740740414395691200141114128295357326229598184018071586974305788868226732005242813320896759331529820049866478352308131616580385299934511466118280591665807797223955525453215604691011751876602867042142279614988764022092484606403624052582624780127296945674634697490506196701799688366808202230402365398158057174489323861985675130681287298635642061814742310222288e36"),
           e_float("-1.6738511240954476689023327972067488225369238527928239541550276691757069937610188837578878392627480535008756330531432808809545621250582908680164326861883513381903576053217033369876889195936938549627908208313612153582914496251704593144245085848861950407359760052977066002157144971561121553241467307849469301272459641517668593147033585279429847348040011936369552387162157219068645165293131121878776959224e68"),
           e_float("3.3196010811393005917762637971726072871814249597961523078542816308193719411685839221266493494678548670627330693833243759638909673978764892852196225848662566988632754043082332209771820908999309538628824249869603845943324140061099184967214272399340149991551000548291108636337181020314553180195457240153910676112759473415869309062912630376132034703298269198095601999612433383513195740486426619656193632088e104"),
           e_float("3.0757758819085446582212220539053157808723863266358213031142285220930584116809964296879685041969827776898382974303669793563116102658700055804474921596196176870509797555920750713171063471517578444551600349415605279394097395523032217865490154263602197514571407783290035630240085205867618489850947239115033338082217684897566932557967242612554307373387747851523665097665299142166775337984369149163208191574e157"),
           e_float("2.5531888224411563483468133017647254028423531214601291911670971127274421790943835357394016289801158553071867406044370159371923830736314943470452562665751165134504986023351808954631484855946520850509796750876217202883472481948660152088490492759226840870971590003053656021524542351132053063933226473960284718689009553649482914368256107500413377522236770668978837474224846382054013617513919658169398818854e217"),
           e_float("5.3859010884536087100624035443063169690692306707658459342546042512107903346457778489236722653212717286049869510895260691557959371240012233884644620693524981285124188090262433927568434314394251910849318164974629914246252309140401905774682488919116388711159831511519766812955386261075593688268029966694399527726045836394665687490001125652797642673513273678569322197453300502920978733056627454474386774794e281"),
           e_float("3.4200867949595716797432116500048703125375475214008769890367845280134618010117130015256588023272529730288651899699054285741846552471234480100775256799619311746236245922108536996307569212219967311396692054281019744190199102835516694146647136227999392579766675952094351253461314141066962847900247922167195503296179034128209110248523432117254953145970893038507651858531667797607206833420999955564810149776e349"),
           e_float("3.5692483912585663253184934092935565972275795552022549815607691154263665357292889660322480320954259534169772602166500488081651565843070489701184146214879994695659717157729810010558558396725264629586159244655485566280267713071022937056892458575717850553663755048312225849250415525564063173968456621680589747085114004893802046176145726642041919325508225194383205433392916501158024004189602272777265728422e420"),
           e_float("1.378638941134802569568546918450868328962211219210953276017114536961323940880364201718822535416756071819622213733931306847733234901691131420756685541111258858915752527267995527527696472868540167154661520146510301567693866919584800794547489141619338825196437763091901415144792170072798426908669051319107196738880540053832869316158780669228039308313590569100530968620154354688467656640047567090538094619e495"),
           e_float("6.8114289945048426455566215856411316598289488075313235223518497594838433854968893893577203290272802962029519061369414645085251310180851342388019242214056715460900727479447896366175902069105063978486419969126620553903848602152873767145880303787078402304991240995560285323642923026777006416401134701499912898767423180315043485029063645605898716469440559202735138218916220617031420149499254461359785919321e575"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00802_conf_hyperg_vary2_mpx(const bool b_write_output)
    {
      return TestCase_case_00802_conf_hyperg_vary2_mpx().execute(b_write_output);
    }
  }
}