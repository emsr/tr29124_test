
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00956_polylog_x_small_neg_n : public TestCaseReal
    {
    public:
      TestCase_case_00956_polylog_x_small_neg_n() { }
      virtual ~TestCase_case_00956_polylog_x_small_neg_n() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00956_polylog_x_small_neg_n");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(21u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-((static_cast<INT32>(2) * k) + static_cast<INT32>(1)), (k + ef::euler_gamma()) / static_cast<INT32>(20));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 21u> a =
        {{
           e_float("0.030601670076875596489099911133176414782524577777504945884411231475457691301201108661519639936858766539308620200659062473156736382440243945259677114177334792552634615842122078982802111286815941207641505027398428389983207017922008119252557836062002248173723283104574814906284441245182220003832093654133183802977387037025112262152693159145766836881107855998152901061119992782162878998345792785174851822441"),
           e_float("0.14477062563044256756741488188736501586005306071674410654083575127783610745588408431077112519205251021219474709156557927315261072315530049477798312881874462468306973070553567514328514449915078267577238222669798789755123351641950458667349502356281432347101622567826591196821122770718348646330266938555112940363483167786400511006479630789110333164275473624948533335642367570660411517751737085416456698409"),
           e_float("1.622309891112157604354220214039709395587400777785431358939026157629345600090505258633045209967363328889627083849944824248029087349789820381357579041348173432027505452399307572355287953999081813526038287082316098460517036518558041305723131430096465307516483031377053727073853570244943337015483275795525034528295802486210857437380708969788977617912412179870810414087284512540580528901336181796738346486"),
           e_float("65.444736424266153870830171350621693220118691079731893816091396631290269454213898503155059193279777615245541086986545312155973279902479857369644601571817283038543957558080343290962325236459796763170617297704755051623082724773244032822376026675252735284849184308817225753161243597370649550764892445989502906311794093414500682025722456965040847131583688515262969033086621506682752418080491770313583864167"),
           e_float("7462.7527479798460479913414854772259361924104208051174700713412682439230544891297126609780426173630056931838565808014106888857313941122163705137626661106503608204799116669669925103272165866179555629822735616672227131262838270549958613913480116530064968402129454587162604540921725477451570554731376797006212717698451546762047086879073830036750221639342368990427553671254233268316343858523575435602670847"),
           e_float("2.121734660250871941974180011902910205061045767031822852756059369023659320501968639884617697655280675932126403026273388308112899381533149356698476964790665150791170454076348147429717887497313520412900029267301687681808870592901496913941055329511813220344182969548333544986551100736002768796365783573186598490056461479865611034592999725161925477244699163393190426060810881827601737066617884737127156364e6"),
           e_float("1.4065429484012582617848897462058639487733676840407808308357032869078804294585429989426140014489528836962233189208621628483552458458503114834282698658257075242200854945737343670470735744141921332490073217664561016575883918710678222985445417744411180008786753955156038305350557447446331010554978450639652908466069029378988409551338499781865622437354587281025699229046248750344292291533748618941872137744e9"),
           e_float("2.10839461806596306635712520155371675243199264168649402035626642177039557862696180733769618850070898297285600754383959679102924400030756782314067553057887863713336437746454255580968201358651994531265036356311331351762377608309236652350099650341349507213361349020625564371238747737984072986059103725683346268882603442145834108056727127376278669699579810975375916295093825573768797860345372706474622382e12"),
           e_float("7.1227558545440552813649880277918784152268686811922627788384685616886694318789292241705416184681328580570376117101150413608455955837068233043439556424243067728079569374526006715159604655373776226809256053928223578083085023139291804489023774891511682196115225954849623503494010798049843450282487557222382141436882725728133689986528808698365478508094320845119811581439949780300950961141079002030472803601e15"),
           e_float("5.5394096253118504736438837288120904563847934626928439863812808151129313982203991102903047944406054207941094608204284262167695863878012321317622753848820437139603256066702727821004103876310165654683863632808743916730660398316837836047486611681894376986245565213197672335935619890473954504093828651435147739767987242624322844852358486673082702698182009875373271427050627592054320426934374222611248908868e19"),
           e_float("1.0395028013392962133912581166337602519378693188372706148350697900676813799605633879446221757822625363692331251674223314539379235752811677549888195137451627813301398425666331060342100436580435952102749619238578721385804279308234633109502077929244637099131502189350176019829143598886138238850204343624706613655766195016244739091276249906777318353464625879573254975699834836036999108769764591323440800451e24"),
           e_float("5.0891458129961303582103779720747130003646393753872472300681048617518447160386470382852975027935194452042272860395777942197547922881421755697637457258406251215211082384149109047537070447071479464085665886106019245831858910212680142066968232112196486139309577201677899408157872920456248897997924817664108923140322317628065406170958463925410349959119792150578433884015484945491063108829041257826400476656e28"),
           e_float("7.3272605113213423892292918583751853012310750900240329398480954276635111103101803352124580535752454429619657626457580149264395713530136609118202001600644308299448207438874786122627012278707556238615153943234944152472595011647221232607006521573441174748090523153727334382707471985184689175632382881166055994515468052992375715578636123279609053882093928925735017964092027840603109179052442930100017816415e33"),
           e_float("3.7189643520022782605822508910178790749409469241969083067241660530418969308042859591804205418869387267430414737586736301761996771856431710998789826169860018655512003713285383104270396176645962863830741460111128531736506394258905886706539605633034910225484927934909187525343481464766192728551276924678542334752086048425989983449632941526981799571643073702649369134465696796214113658825342782614024957317e39"),
           e_float("8.8042921609472495954682228841416656921820681099207578416654764479429660821293909940534977576956998122022658450084387864806785397093846483296483746091773694620261474972309298631585723593470050041840934704484356772783924778181925347843731183450398560568459768313752705409261449417375943488195694958156353932440576607235492864143356354821555989117613569382597597288635990498039402283697762024910678357196e45"),
           e_float("1.5318798847413772468020912536418610288537883275306715793816287924360780715526961381505844962557269792975787170366231675558498992371972946886894993661286503821810898853991073487584663188933870940271200016060534271929710188606829400632127988476318522036253436333976377895702667179956432233701248579792961257948343066880872056637805048795429235454609392249059344291127083781542731945007337938486674856766e53"),
           e_float("4.3717851006288660255178559771114586526429364959577233010659952625478541733016613236918765896534500213390687902888720458689514922354374927012380753471775168634204506291300201295417518773739211451439747375682454150821086341690630877067945327637419160198732826776488939153796655727333658280272394228546301096075466640021769498378578163383129270421008450588331158317046014588224968487951606789388120087539e61"),
           e_float("1.0408869649974536581140623698336106925959742303848971693770902576432636595568081113730640687092629995157719508039206837477018977638652747586456150353956445093050074946290697241122452456253746017321880049115670522649985812127690747706085205482853345222738993057666381568636307376837553377265987713250496046627460669131433316803969494917317951955324174570188874881535897907562440689135777777781103517487e72"),
           e_float("1.4236341425774448049769951173916161666541896554423434734057025708282786733002174979395380781939322045315394414193760507563552686488370031189451023144435226020087271478649109446050195905469295970760331119785069602350299163466199706963045759290321330654521549587047876240053619475465420431078496840458825473968092715241572512903122668418003119203117894013602797330099665528903025297750669826666750682874e86"),
           e_float("1.320621652596931870741119459826100247705254338481153031850737484041790773500352610860485519428309410124137016411902881083450104831493265539437919277078507121420058001311738215785655134007628939281643381705982215516757937254603222178950982514379979694821283650281737498370396814212891345027739825766168305341200812193980762412571753304129745549622873960164327952970825462574879616053307538573473980029e113"),
           e_float("2.8287352889582044980280234551692002985941988282291353443604489731395823952012887544607023575901932041049510163849910812006268296686589458305386643606666255984868488981177772930136488247702248128508944930673085139911086585504617651755154487948328233229689959857387909613234008057588727763289822893614657043758874004291151007115385991038998294833036074303595130387447004828951283349023362013267291571252e114"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00956_polylog_x_small_neg_n(const bool b_write_output)
    {
      return TestCase_case_00956_polylog_x_small_neg_n().execute(b_write_output);
    }
  }
}