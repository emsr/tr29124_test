{Part 7A of regression test for SPECFUNX unit  (c) 2018  W.Ehrhardt}

unit t_sfx7a;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}


interface

procedure test_besselpolyx;
procedure test_debyex;
procedure test_einsteinx;
procedure test_eulerpolyx;
procedure test_expnx;
procedure test_transportx;
procedure test_keplerx;
procedure test_langevinx;
procedure test_langevin_invx;
procedure test_rrcfx;


implementation

uses
  AMath,SpecFunX,t_sfx0;


{---------------------------------------------------------------------------}
procedure test_debyex;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
  f1 : THexExtW = ($A3D7,$1C70,$01C7,$FFE0,$3FFE);  {+9.9951182471380889183E-1}
  f72: THexExtW = ($6829,$6D07,$8498,$C5B7,$3FFB);  {.9654143896262086549358612680186751027932e-1}

begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','debyex');

  {Test values debye(n,x), n=1,3,4 from MISCFUN [22]}
  x := 1/512;
  y := debyex(1,x);
  {f := 0.99951182471380889183;}
  f := extended(f1);
  testrel( 1, NE, y, f, cnt,failed);

  x := 1/32;
  y := debyex(1,x);
  f := 0.99221462647120597836;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1/8;
  y := debyex(1,x);
  f := 0.96918395997895308324;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1/2;
  y := debyex(1,x);
  f := 0.88192715679060552968;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1.0;
  y := debyex(1,x);
  f := 0.77750463411224827642;
  testrel( 5, NE, y, f, cnt,failed);

  x := 2.0;
  y := debyex(1,x);
  f := 0.60694728460981007205;
  testrel( 6, NE, y, f, cnt,failed);

  x := 3.0;
  y := debyex(1,x);
  f := 0.48043521957304283829;
  testrel( 7, NE, y, f, cnt,failed);

  x := 4.0;
  y := debyex(1,x);
  f := 0.38814802129793784501;
  testrel( 8, NE, y, f, cnt,failed);

  x := 17/4;
  y := debyex(1,x);
  f := 0.36930802829242526815;
  testrel( 9, NE, y, f, cnt,failed);

  x := 5.0;
  y := debyex(1,x);
  f := 0.32087619770014612104;
  testrel(10, NE, y, f, cnt,failed);

  x := 11/2;
  y := debyex(1,x);
  f := 0.29423996623154246701;
  testrel(11, NE, y, f, cnt,failed);

  x := 6.0;
  y := debyex(1,x);
  f := 0.27126046678502189985;
  testrel(12, NE, y, f, cnt,failed);

  x := 8.0;
  y := debyex(1,x);
  f := 0.20523930310221503723;
  testrel(13, NE, y, f, cnt,failed);

  x := 10.0;
  y := debyex(1,x);
  f := 0.16444346567994602563;
  testrel(14, NE, y, f, cnt,failed);

  x := 20.0;
  y := debyex(1,x);
  f := 0.82246701178200016086e-1;
  testrel(15, NE, y, f, cnt,failed);

  x := 50.0;
  y := debyex(1,x);
  f := 0.32898681336964528729e-1;
  testrel(16, NE, y, f, cnt,failed);

  x := 1/512;
  y := debyex(3,x);
  f := 0.99926776885985461940;
  testrel(17, NE, y, f, cnt,failed);

  x := 1/32;
  y := debyex(3,x);
  f := 0.98833007755734698212;
  testrel(18, NE, y, f, cnt,failed);

  x := 1/8;
  y := debyex(3,x);
  f := 0.95390610472023510237;
  testrel(19, NE, y, f, cnt,failed);

  x := 1/2;
  y := debyex(3,x);
  f := 0.82496296897623372315;
  testrel(20, NE, y, f, cnt,failed);

  x := 1.0;
  y := debyex(3,x);
  f := 0.67441556407781468010;
  testrel(21, NE, y, f, cnt,failed);

  x := 2.0;
  y := debyex(3,x);
  f := 0.44112847372762418113;
  testrel(22, NE, y, f, cnt,failed);

  x := 3.0;
  y := debyex(3,x);
  f := 0.28357982814342246206;
  testrel(23, NE, y, f, cnt,failed);

  x := 4.0;
  y := debyex(3,x);
  f := 0.18173691382177474795;
  testrel(24, NE, y, f, cnt,failed);

  x := 17/4;
  y := debyex(3,x);
  f := 0.16277924385112436877;
  testrel(25, NE, y, f, cnt,failed);

  x := 5.0;
  y := debyex(3,x);
  f := 0.11759741179993396450;
  testrel(26, NE, y, f, cnt,failed);

  x := 11/2;
  y := debyex(3,x);
  f := 0.95240802723158889887e-1;
  testrel(27, NE, y, f, cnt,failed);

  x := 6.0;
  y := debyex(3,x);
  f := 0.77581324733763020269e-1;
  testrel(28, NE, y, f, cnt,failed);

  x := 8.0;
  y := debyex(3,x);
  f := 0.36560295673194845002e-1;
  testrel(29, NE, y, f, cnt,failed);

  x := 10.0;
  y := debyex(3,x);
  f := 0.19295765690345489563e-1;
  testrel(30, NE, y, f, cnt,failed);

  x := 20.0;
  y := debyex(3,x);
  f := 0.24352200674805479827e-2;
  testrel(31, NE, y, f, cnt,failed);

  x := 50.0;
  y := debyex(3,x);
  f := 0.15585454565440389896e-3;
  testrel(32, NE, y, f, cnt,failed);

  x := 1/512;
  y := debyex(4,x);
  f := 0.99921896192761576256;
  testrel(33, NE, y, f, cnt,failed);

  x := 1/32;
  y := debyex(4,x);
  f := 0.98755425280996071022;
  testrel(34, NE, y, f, cnt,failed);

  x := 1/8;
  y := debyex(4,x);
  f := 0.95086788606389739976;
  testrel(35, NE, y, f, cnt,failed);

  x := 1/2;
  y := debyex(4,x);
  f := 0.81384569172034042516;
  testrel(36, NE, y, f, cnt,failed);

  x := 1.0;
  y := debyex(4,x);
  f := 0.65487406888673697092;
  testrel(37, NE, y, f, cnt,failed);

  x := 2.0;
  y := debyex(4,x);
  f := 0.41189273671788528876;
  testrel(38, NE, y, f, cnt,failed);

  x := 3.0;
  y := debyex(4,x);
  f := 0.25187863642883314410;
  testrel(39, NE, y, f, cnt,failed);

  x := 4.0;
  y := debyex(4,x);
  f := 0.15185461258672022043;
  testrel(40, NE, y, f, cnt,failed);

  x := 17/4;
  y := debyex(4,x);
  f := 0.13372661145921413299;
  testrel(41, NE, y, f, cnt,failed);

  x := 5.0;
  y := debyex(4,x);
  f := 0.91471377664481164749e-1;
  testrel(42, NE, y, f, cnt,failed);

  x := 11/2;
  y := debyex(4,x);
  f := 0.71227828197462523663e-1;
  testrel(43, NE, y, f, cnt,failed);

  x := 6.0;
  y := debyex(4,x);
  f := 0.55676547822738862783e-1;
  testrel(44, NE, y, f, cnt,failed);

  x := 8.0;
  y := debyex(4,x);
  f := 0.21967566525574960096e-1;
  testrel(45, NE, y, f, cnt,failed);

  x := 10.0;
  y := debyex(4,x);
  f := 0.96736755602711590082e-2;
  testrel(46, NE, y, f, cnt,failed);

  x := 20.0;
  y := debyex(4,x);
  f := 0.62214648623965450200e-3;
  testrel(47, NE, y, f, cnt,failed);

  x := 50.0;
  y := debyex(4,x);
  f := 0.15927210319002161231e-4;
  testrel(48, NE, y, f, cnt,failed);

  {Rest of test values for n=6,8,12 calculated with Maple }
  {f := x->n*int(t^n/(exp(t)-1),t=0..x)/x^n; and Digits:=50}
  x := 1/512;
  y := debyex(6,x);
  f := 0.9991631848471384035;
  testrel(49, NE, y, f, cnt,failed);

  x := 1/32;
  y := debyex(6,x);
  f := 0.9866681772186796587;
  testrel(50, NE, y, f, cnt,failed);

  x := 1/8;
  y := debyex(6,x);
  f := 0.9474049305411031823;
  testrel(51, NE, y, f, cnt,failed);

  x := 1/2;
  y := debyex(6,x);
  f := 0.8012874593544054948;
  testrel(52, NE, y, f, cnt,failed);

  x := 1.0;
  y := debyex(6,x);
  f := 0.6331114258349510759;
  testrel(53, NE, y, f, cnt,failed);

  x := 2.0;
  y := debyex(6,x);
  f := 0.3804986630746610429;
  testrel(54, NE, y, f, cnt,failed);

  x := 3.0;
  y := debyex(6,x);
  f := 0.2193992525257245836;
  testrel(55, NE, y, f, cnt,failed);

  x := 4.0;
  y := debyex(6,x);
  f := 0.1229278562814578228;
  testrel(56, NE, y, f, cnt,failed);

  x := 17/4;
  y := debyex(6,x);
  f := 0.1060375248597196031;
  testrel(57, NE, y, f, cnt,failed);

  x := 5.0;
  y := debyex(6,x);
  f := 0.6777784974890353731e-1;
  testrel(58, NE, y, f, cnt,failed);

  x := 11/2;
  y := debyex(6,x);
  f := 0.5020600934448088116e-1;
  testrel(59, NE, y, f, cnt,failed);

  x := 6.0;
  y := debyex(6,x);
  f := 0.3719333613705515670e-1;
  testrel(60, NE, y, f, cnt,failed);

  x := 8.0;
  y := debyex(6,x);
  f := 0.1145231921902748610e-1;
  testrel(61, NE, y, f, cnt,failed);

  x := 10.0;
  y := debyex(6,x);
  f := 0.3793849329461595528e-2;
  testrel(62, NE, y, f, cnt,failed);

  x := 20.0;
  y := debyex(6,x);
  f := 0.6804635545479456894e-4;
  testrel(63, NE, y, f, cnt,failed);

  x := 50.0;
  y := debyex(6,x);
  f := 0.2787884082105527120e-6;
  testrel(64, NE, y, f, cnt,failed);

  x := 1/512;
  y := debyex(12,x);
  f := 0.9990988301706686501;
  testrel(65, NE, y, f, cnt,failed);

  x := 1/32;
  y := debyex(12,x);
  f := 0.9856466765478185763;
  testrel(66, NE, y, f, cnt,failed);

  x := 1/8;
  y := debyex(12,x);
  f := 0.9434235095071814036;
  testrel(67, NE, y, f, cnt,failed);

  x := 1/2;
  y := debyex(12,x);
  f := 0.7870231504611680153;
  testrel(68, NE, y, f, cnt,failed);

  x := 1.0;
  y := debyex(12,x);
  f := 0.6088700041762235049;
  testrel(69, NE, y, f, cnt,failed);

  x := 2.0;
  y := debyex(12,x);
  f := 0.3472653175019342084;
  testrel(70, NE, y, f, cnt,failed);

  x := 3.0;
  y := debyex(12,x);
  f := 0.1872401059320712096;
  testrel(71, NE, y, f, cnt,failed);

  x := 4.0;
  y := debyex(12,x);
  {f := 0.9654143896262086549e-1;}
  f := extended(f72);
  testrel(72, NE, y, f, cnt,failed);

  x := 17/4;
  y := debyex(12,x);
  f := 0.8134774441706960165e-1;
  testrel(73, NE, y, f, cnt,failed);

  x := 5.0;
  y := debyex(12,x);
  f := 0.4814185645191148541e-1;
  testrel(74, NE, y, f, cnt,failed);

  x := 11/2;
  y := debyex(12,x);
  f := 0.3367880055948328374e-1;
  testrel(75, NE, y, f, cnt,failed);

  x := 6.0;
  y := debyex(12,x);
  f := 0.2344811348723500784e-1;
  testrel(76, NE, y, f, cnt,failed);

  x := 8.0;
  y := debyex(12,x);
  f := 0.5344588786833453221e-2;
  testrel(77, NE, y, f, cnt,failed);

  x := 10.0;
  y := debyex(12,x);
  f := 0.1198815360618837356e-2;
  testrel(78, NE, y, f, cnt,failed);

  x := 20.0;
  y := debyex(12,x);
  f := 0.1348750701799345211e-5;
  testrel(79, NE, y, f, cnt,failed);

  x := 50.0;
  y := debyex(12,x);
  f := 0.2354677578932315411e-10;
  testrel(80, NE, y, f, cnt,failed);

  {D(n,20) n=20..200}
  x := 20;
  f := 0.2045985597891880435e-6;
  y := debyex(20,x);
  testrel(81, NE, y, f, cnt,failed);

  f := 0.6521968411709099410e-7;
  y := debyex(50,x);
  testrel(82, NE, y, f, cnt,failed);

  f := 0.5964824507515076368e-7;
  y := debyex(60,x);
  testrel(83, NE, y, f, cnt,failed);

  f := 0.5378148028857494088e-7;
  y := debyex(80,x);
  testrel(84, NE, y, f, cnt,failed);

  f := 0.5074077974146328750e-7;
  y := debyex(100,x);
  testrel(85, NE, y, f, cnt,failed);

  f := 0.4714758392337898263e-7;
  y := debyex(150,x);
  testrel(86, NE, y, f, cnt,failed);

  f := 0.4552275137444157216e-7;
  y := debyex(200,x);
  testrel(87, NE, y, f, cnt,failed);


  {D(n,5) n=100..10000}
  x := 5;
  f := 0.3532560165653401053e-1;
  y := debyex(100,x);
  testrel(88, NE, y, f, cnt,failed);

  f := 0.34193472837210103295e-1;
  y := debyex(500,x);
  testrel(89, 20, y, f, cnt,failed);

  f := 0.3410139487334168171e-1;
  y := debyex(750,x);
  testrel(90, 40, y, f, cnt,failed);

  f := 0.3405548547450571961e-1;
  y := debyex(1000,x);
  testrel(91, 40, y, f, cnt,failed);

  f := 0.3398678310287432440e-1;
  y := debyex(2000,x);
  testrel(92, 100, y, f, cnt,failed);

  f := 0.3393196075652582305e-1;
  y := debyex(10000,x);
  testrel(93, 600, y, f, cnt,failed);

  {Test after fix for large x}
  y := debyex(1, 100000.0);
  f := 0.1644934066848226436e-4;
  testrel(94, NE, y, f, cnt,failed);

  x := 20000.0;
  y := debyex(2, x);
  f := 0.1202056903159594285e-7;
  testrel(95, NE, y, f, cnt,failed);

  y := debyex(7, x);
  f := 0.2767488213020584085e-25;
  testrel(96, NE, y, f, cnt,failed);

  y := debyex(10, x);
  f := 0.3545501280865848353e-35;
  testrel(97, NE, y, f, cnt,failed);

  y := debyex(8, 30000.0);
  f := 0.4926197640450862355e-30;
  testrel(98, NE, y, f, cnt,failed);

  y := debyex(20, 12000);
  f := 0.12691995186452421609e-61;
  testrel(99, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_einsteinx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','einsteinx');

  y := einsteinx(1,0);
  f := 1;
  testrel(1, NE, y, f, cnt,failed);

  y := einsteinx(2,0);
  f := 1;
  testrel(2, NE, y, f, cnt,failed);

  y := einsteinx(1,-1);
  f := 0.9206735942077923189;
  testrel(3, NE, y, f, cnt,failed);

  y := einsteinx(2,-1);
  f := 1.581976706869326424;
  testrel(4, NE1, y, f, cnt,failed);

  y := einsteinx(1,1);
  f := 0.9206735942077923189;
  testrel(5, NE, y, f, cnt,failed);

  y := einsteinx(2,1);
  f := 0.5819767068693264244;
  testrel(6, NE, y, f, cnt,failed);

  y := einsteinx(3,1);
  f := -0.458675145387081891;
  testrel(7, NE, y, f, cnt,failed);

  y := einsteinx(4,1);
  f := 1+0.0406518522564083154;
  testrel(8, NE, y, f, cnt,failed);

  y := einsteinx(1,10);
  f := 0.4540405235047541210e-2;
  testrel(9, NE, y, f, cnt,failed);

  y := einsteinx(2,10);
  f := 0.4540199100968776833e-3;
  testrel(10, NE, y, f, cnt,failed);

  y := einsteinx(3,10);
  f := -0.4540096037048920950e-4;
  testrel(11, NE, y, f, cnt,failed);

  y := einsteinx(4,10);
  f := 0.4994208704673668928e-3;
  testrel(12, NE, y, f, cnt,failed);

  y := einsteinx(2,-10);
  f := 10+0.0004540199100968776833;
  testrel(13, NE, y, f, cnt,failed);

  y := einsteinx(1,45);
  f := 0.5796600125612522130e-16;
  testrel(14, NE, y, f, cnt,failed);

  y := einsteinx(1,44);
  f := 0.1506427201883503006e-15;
  testrel(15, NE, y, f, cnt,failed);

  y := einsteinx(1,38);
  f := 0.4532907751717355068e-13;
  testrel(16, NE, y, f, cnt,failed);

  y := einsteinx(4,45);
  f := 0.1316758547052721076e-17;
  testrel(17, NE1, y, f, cnt,failed);

  y := einsteinx(4,44);
  f := 0.3501509508510208432e-17;
  testrel(18, NE, y, f, cnt,failed);

  y := einsteinx(4,MinDouble);
  f := 709.3964185322641062;
  testrel(19, NE, y, f, cnt,failed);

  y := einsteinx(4,100);
  f := 0.3757276735781044323e-41;
  testrel(20, NE, y, f, cnt,failed);

  y := einsteinx(2,-10000);
  f := 10000.0;
  testrel(21, NE, y, f, cnt,failed);

  y := einsteinx(2,-43);
  f := 43.00000000000000001;
  testrel(22, NE, y, f, cnt,failed);

  y := einsteinx(2,-38);
  f := 38.00000000000000119;
  testrel(23, NE, y, f, cnt,failed);

  y := einsteinx(3,44);
  f := -0.7781132241133796516e-19;
  testrel(24, NE, y, f, cnt,failed);

  y := einsteinx(3,38);
  f := -0.3139132792048029678e-16;
  testrel(25, NE, y, f, cnt,failed);

  {----------------------------------}
  y := einsteinx(3,MinExtended);
  f := -11355.13711193302406;
  testrel(26, NE, y, f, cnt,failed);

  y := einsteinx(4,MinExtended);
  f := 11356.13711193302406;
  testrel(27, NE, y, f, cnt,failed);

  y := einsteinx(4,10000);
  f := 0.1135597413701267572e-4338;
  testrel(28, NE, y, f, cnt,failed);

  y := einsteinx(1,10000);
  f := 0.1135483865314736099e-4334;
  testrel(29, NE1, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_keplerx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 4;
  NE2 = 16;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','keplerx');

  {elliptic, Maple: kepler := (m,e) -> fsolve(y-e*sin(y) = m, y=Pi);}
  y := keplerx(1, 0.5);
  f := 1.498701133517848314;
  testrel(1, NE, y, f, cnt,failed);

  y := keplerx(-1000, 0.5);
  f := -1000.497514775673146;
  testrel(2, NE, y, f, cnt,failed);

  y := keplerx(Pi_2, 0.75);
  f := 2.184106679498448926;
  testrel(3, NE, y, f, cnt,failed);

  y := keplerx(0.5, 0.125);
  f := 0.5671542538034771510;
  testrel(4, NE, y, f, cnt,failed);

  y := keplerx(0.125, 0.9375);
  f := 0.7928034322756140260;
  testrel(5, NE, y, f, cnt,failed);

  y := keplerx(1.5, 0.9375);
  f := 2.237023829054169401;
  testrel(6, NE, y, f, cnt,failed);

  y := keplerx(10.125, 0.9375);
  f := 9.790088287071411171;
  testrel(7, NE, y, f, cnt,failed);

  y := keplerx(0.125, 0.99609375);
  f := 0.9136395753855342618;
  testrel(8, NE, y, f, cnt,failed);

  {in difficult region m near 0, e near 1}
  y := 1-ldexp(1,-20);
  y := keplerx(1/1024, y);
  f := 0.1803684433746817368;
  testrel(9, NE2, y, f, cnt,failed);

  {arguments from literature, somewhat inexact with binary}
  y := keplerx(0.2, 0.99);
  f := 1.066997365281563186;
  testrel(10, NE, y, f, cnt,failed);

  y := keplerx(0.06, 0.6);
  f := 0.1491710835982268287;
  testrel(11, NE, y, f, cnt,failed);

  y := keplerx(1, 0.9);
  f := 1.862086686874532255;
  testrel(12, NE, y, f, cnt,failed);

  {hyperbolic, Maple: kepler_hyp := (m,e) -> fsolve(e*sinh(x) - x - m, x = signum(m)*ln(2*abs(m)/e + 1.8));}
  y := keplerx(-1000, 10);
  f := -5.303631719539061703;
  testrel(13, NE, y, f, cnt,failed);

  y := keplerx(1, 2);
  f := 0.8140967963021331692;
  testrel(14, NE, y, f, cnt,failed);

  y := keplerx(6,2);
  f := 2.107689797681256377;
  testrel(15, NE, y, f, cnt,failed);

  y := keplerx(0.5,1.5);
  f := 0.7673431749540970103;
  testrel(16, NE, y, f, cnt,failed);

  y := keplerx(MaxExtended,1.5);
  f := 11356.811088366595730;
  testrel(17, NE, y, f, cnt,failed);

  y := keplerx(10,6);
  f := 1+0.3978298998186000144;
  testrel(18, NE, y, f, cnt,failed);

  y := keplerx(10000,20);
  f := 6.9084468837654158448;
  testrel(19, NE, y, f, cnt,failed);

  y := keplerx(1,20);
  f := 0.5260603476886937829e-1;
  testrel(20, NE, y, f, cnt,failed);

  y := keplerx(0,2);
  f := 0;
  testrel(21, NE, y, f, cnt,failed);

  y := keplerx(1e-6,1.5);
  f := 0.1999999999996000000e-5;
  testrel(22, NE, y, f, cnt,failed);

  {parabolic, Maple: kepler_para := m -> fsolve(x + x^3/3 = m, x = m^(1/3));}
  y := keplerx(2, 1);
  f := 1.287909750704127236;
  testrel(23, NE, y, f, cnt,failed);

  y := keplerx(1, 1);
  f := 0.8177316738868235061;
  testrel(24, NE, y, f, cnt,failed);

  y := keplerx(0.5,1.0);
  f := 0.4662205239107734274;
  testrel(25, NE, y, f, cnt,failed);

  y := keplerx(0.125, 1);
{$ifdef BIT16}
  f := 0.2487178477269259601*0.5;
{$else}
  f := 0.124358923863462980055;
{$endif}
  testrel(26, NE, y, f, cnt,failed);

  y := keplerx(1/1024, 1);
{$ifdef BIT16}
  f := 0.1953124379118875708e-2*0.5;
{$else}
  f := 0.9765621895594378539e-3;
{$endif}
  testrel(27, NE, y, f, cnt,failed);

  y := keplerx(-1000, 1);
  {f:= -14.35316011237345298;}
  f := -(14 + 0.3531601123734529825);
  testrel(28, NE, y, f, cnt,failed);

  y := keplerx(sqrt_epsh, 1);
  f := 0.2328306436538696289e-9;
  testrel(29, NE, y, f, cnt,failed);

  y := keplerx(1.25e30, 1);
  f := 0.1553616252976929433e11;
  testrel(30, NE1, y, f, cnt,failed);

  y := keplerx(MaxExtended, 1);
  f := 0.1528234751400654562e1645;
  testrel(31, NE, y, f, cnt,failed);

  y := keplerx(0.25*MaxExtended, 1);
{$ifdef BIT16}
  f := 0.1925455132470543184e1645*0.5;
{$else}
  f := 0.9627275662352715918e1644;
{$endif}
  testrel(32, NE1, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_langevinx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LangevinLx');

  y := LangevinLx(0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := LangevinLx(1e-10);
  f := 0.3333333333333333333e-10;
  testrel(2, NE, y, f, cnt,failed);

  y := LangevinLx(-1/1024);
  f := -0.3255208126372779994e-3;
  testrel(3, NE, y, f, cnt,failed);

  y := LangevinLx(0.25);
  f := 0.8298816507359656826e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := LangevinLx(0.5+1/1024);
  f := 0.1642632531019772181;
  testrel(5, NE, y, f, cnt,failed);

  y := LangevinLx(1-1/1024);
  f := 0.3127657674731261456;
  testrel(6, NE, y, f, cnt,failed);

  y := LangevinLx(-1);
  f := -0.3130352854993313036;
  testrel(7, NE, y, f, cnt,failed);

  y := LangevinLx(1.125);
  f := 0.3467452117192722997;
  testrel(8, NE, y, f, cnt,failed);

  y := LangevinLx(2);
  f := 0.5373147207275480959;
  testrel(9, NE, y, f, cnt,failed);

  y := LangevinLx(-2.5);
  f := -0.6135673098126084622;
  testrel(10, NE, y, f, cnt,failed);

  y := LangevinLx(10);
  f := 0.9000000041223072534;
  testrel(11, NE, y, f, cnt,failed);

  y := LangevinLx(23);
{$ifdef BIT16}
      {1.234567890123456789}
  f := 1.913043478260869565*0.5;
{$else}
  f := 0.9565217391304347826;
{$endif}
  testrel(12, NE, y, f, cnt,failed);

  y := LangevinLx(25);
  f := 0.96;
  testrel(13, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_langevin_invx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LangevinL_invx');

  y := LangevinL_invx(0.0);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := LangevinL_invx(1e-8);
  f := 0.3000000000000000180e-7;
  testrel(2, NE, y, f, cnt,failed);

  y := LangevinL_invx(2e-5);
  f := 0.60000000014400000005e-4;
  testrel(3, NE, y, f, cnt,failed);

  y := LangevinL_invx(5e-5);
  f := 0.15000000022500000053e-3;
  testrel(4, NE, y, f, cnt,failed);

  y := LangevinL_invx(-1/1024);
  f := -0.2929689176382141675e-2;
  testrel(5, NE, y, f, cnt,failed);

  y := LangevinL_invx(-0.125);
  f := -0.37856827056108226381;
  testrel(6, NE, y, f, cnt,failed);

  y := LangevinL_invx(0.25);
  f := 0.7798973686506122298;
  testrel(7, NE, y, f, cnt,failed);

  y := LangevinL_invx(-0.375);
  f := -1.234664583210388140;
  testrel(8, NE1, y, f, cnt,failed);

  y := LangevinL_invx(0.5);
  f := 1.796755984723713041;
  testrel(9, NE, y, f, cnt,failed);

  y := LangevinL_invx(-0.75);
  f := -3.989053868488535163;
  testrel(10, NE, y, f, cnt,failed);

  y := LangevinL_invx(0.9375);
  f := 15.99999999999351595;
  testrel(11, NE1, y, f, cnt,failed);

  y := LangevinL_invx(0.96875);
  f := 32.0;
  testrel(12, NE, y, f, cnt,failed);

  y := LangevinL_invx(0.984375);
  f := 64.0;
  testrel(13, NE, y, f, cnt,failed);

  y := LangevinL_invx(0.9990234375);
  f := 1024;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_transportx;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE1 = 16;
  NE2 = 220;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','transportx');

  {Tests for MacLeod x, function values for n=2,5,9 from MacLeod}
  {others with Maple transport := (n,x) -> int(t^n*exp(t)/(exp(t)-1)^2, t=0..x);}

  {---------------------}
  x := 0.001953125;
  f := 0.19531247930394515480e-02;
  y := transportx(2,x);
  testrel(1, NE, y, f, cnt,failed);

  f := 0.36379780361036116971e-11;
  y := transportx(5,x);
  testrel(2, NE, y, f, cnt,failed);

  f := 0.26469772870084897671e-22;
  y := transportx(9,x);
  testrel(3, NE, y, f, cnt,failed);

  f := 0.61184093997606970767e-270;
  y := transportx(100,x);
  testrel(4, NE, y, f, cnt,failed);

  {---------------------}
  x := 0.03125;
  f := 0.31249152314331109004e-01;
  y := transportx(2,x);
  testrel(5, NE, y, f, cnt,failed);

  f := 0.23840564453948442379e-06;
  y := transportx(5,x);
  testrel(6, NE, y, f, cnt,failed);

  y := transportx(9,x);
  f := 0.11367943653594246210e-12;
  testrel(7, NE, y, f, cnt,failed);

  y := transportx(100,x);
  f := 0.98737541391723215213e-151;
  testrel(8, NE, y, f, cnt,failed);

  {---------------------}
  x := 0.125;
  f := 0.12494577194783451032e+00;
  y := transportx(2,x);
  testrel(9, NE, y, f, cnt,failed);

  f := 0.60982205372226969189e-04;
  y := transportx(5,x);
  testrel(10, NE, y, f, cnt,failed);

  f := 0.74428246255329800255e-08;
  y := transportx(9,x);
  testrel(11, NE, y, f, cnt,failed);

  y := transportx(100,x);
  f := 0.39618850817575413053e-91;
  testrel(12, NE, y, f, cnt,failed);

  {---------------------}
  x := 0.5;
  f := 0.49655363615640595865e+00;
  y := transportx(2,x);
  testrel(13, NE, y, f, cnt,failed);

  f := 0.15410004586376649337e-01;
  y := transportx(5,x);
  testrel(14, NE, y, f, cnt,failed);

  f := 0.48022728485415366194e-03;
  y := transportx(9,x);
  testrel(15, NE, y, f, cnt,failed);

  f := 0.15615096769038751704e-31;
  y := transportx(100,x);
  testrel(16, NE, y, f, cnt,failed);

  {---------------------}
  x := 1.0;
  y := transportx(2,x);
  f := 0.97303256135517012845e+00;
  testrel(17, NE, y, f, cnt,failed);

  y := transportx(5,x);
  f := 0.23661587923909478926e+00;
  testrel(18, NE, y, f, cnt,failed);

  y := transportx(9,x);
  f := 0.11700243014358676725e+00;
  testrel(19, NE, y, f, cnt,failed);

  y := transportx(100,x);
  f := 0.93148583747959048513e-2;
  testrel(20, NE, y, f, cnt,failed);

  {---------------------}
  x := 1.5;
  f := 0.14121978695932525805e+01;
  y := transportx(2,x);
  testrel(21, NE, y, f, cnt,failed);

  f := 0.11198756851307629651e+01;
  y := transportx(5,x);
  testrel(22, NE, y, f, cnt,failed);

  f := 0.27648973910899914391e+01;
  y := transportx(9,x);
  testrel(23, NE, y, f, cnt,failed);

  f := 0.2285612236102933935e16;
  y := transportx(100,x);
  testrel(24, NE, y, f, cnt,failed);


  {---------------------}
  x := 2.0;
  y := transportx(2,x);
  f := 0.18017185674405776809e+01;
  testrel(25, NE, y, f, cnt,failed);

  y := transportx(5,x);
  f := 0.32292901663684049171e+01;
  testrel(26, NE, y, f, cnt,failed);

  y := transportx(9,x);
  f := 0.24716631405829192997e+02;
  testrel(27, NE, y, f, cnt,failed);

  y := transportx(100,x);
  f := 0.4664586943216141974e28;
  testrel(28, NE, y, f, cnt,failed);

  {---------------------}
  x := 2.5;
  f := 0.21350385339277043015e+01;
  y := transportx(2,x);
  testrel(29, NE, y, f, cnt,failed);

  f := 0.70362973105160654056e+01;
  y := transportx(5,x);
  testrel(30, NE, y, f, cnt,failed);

  f := 0.12827119828849828583e+03;
  y := transportx(9,x);
  testrel(31, NE, y, f, cnt,failed);

  y := transportx(100,x);
  f := 0.15454773358862196301e38;
  testrel(32, NE, y, f, cnt,failed);

  {---------------------}
  x := 3.0;
  f := 0.24110500490169534620e+01;
  y := transportx(2,x);
  testrel(33, NE, y, f, cnt,failed);

  f := 0.12770557691044159511e+02;
  y := transportx(5,x);
  testrel(34, NE, y, f, cnt,failed);

  f := 0.46842894800662208986e+03;
  y := transportx(9,x);
  testrel(35, NE, y, f, cnt,failed);

  f := 0.87254293066564685718e45;
  y := transportx(100,x);
  testrel(36, NE, y, f, cnt,failed);

  {---------------------}
  x := 4.0;
  f := 0.28066664045631179931e+01;
  y := transportx(2,x);
  testrel(37, NE, y, f, cnt,failed);

  f := 0.29488339015245845447e+02;
  y := transportx(5,x);
  testrel(38, NE, y, f, cnt,failed);

  f := 0.31673967371627895718e+04;
  y := transportx(9,x);
  testrel(39, NE, y, f, cnt,failed);

  f := 0.12608777929089504884e58;
  y := transportx(100,x);
  testrel(40, NE1, y, f, cnt,failed);

  {---------------------}
  x := 4.25;
  f := 0.28777421863296234131e+01;
  y := transportx(2,x);
  testrel(41, NE, y, f, cnt,failed);

  f := 0.34471340540362254586e+02;
  y := transportx(5,x);
  testrel(42, NE, y, f, cnt,failed);

  f := 0.46140886546630195390e+04;
  y := transportx(9,x);
  testrel(43, NE, y, f, cnt,failed);

  f := 0.44538412510635393439e60;
  y := transportx(100,x);
  testrel(44, NE, y, f, cnt,failed);

  {---------------------}
  x := 5.0;
  f := 0.30391706043438554330e+01;
  y := transportx(2,x);
  testrel(45, NE, y, f, cnt,failed);

  f := 0.50263092218175187785e+02;
  y := transportx(5,x);
  testrel(46, NE, y, f, cnt,failed);

  f := 0.11952718545392302185e+05;
  y := transportx(9,x);
  testrel(47, NE, y, f, cnt,failed);

  f := 0.2806636405928563965e67;
  y := transportx(100,x);
  testrel(48, NE, y, f, cnt,failed);

  {---------------------}
  x := 5.5;
  f := 0.31125074928667355940e+01;
  y := transportx(2,x);
  testrel(49, NE, y, f, cnt,failed);

  f := 0.60819909101127165207e+02;
  y := transportx(5,x);
  testrel(50, NE, y, f, cnt,failed);

  f := 0.20001612666477027728e+05;
  y := transportx(9,x);
  testrel(51, NE, y, f, cnt,failed);

  f := 0.2579419206218227689e71;
  y := transportx(100,x);
  testrel(52, NE1, y, f, cnt,failed);

  {---------------------}
  x := 6.0;
  f := 0.31656687817738577185e+01;
  y := transportx(2,x);
  testrel(53, NE, y, f, cnt,failed);

  f := 0.70873334429213460498e+02;
  y := transportx(5,x);
  testrel(54, NE, y, f, cnt,failed);

  f := 0.31011073271851366554e+05;
  y := transportx(9,x);
  testrel(55, NE, y, f, cnt,failed);

  f := 0.1027542649495980112e75;
  y := transportx(100,x);
  testrel(56, NE1, y, f, cnt,failed);

  {---------------------}
  x := 8.0;
  f := 0.32623520367816009184e+01;
  y := transportx(2,x);
  testrel(57, NE, y, f, cnt,failed);

  f := 0.10147781242977788097e+03;
  y := transportx(5,x);
  testrel(58, NE, y, f, cnt,failed);

  f := 0.10352949905541130133e+06;
  y := transportx(9,x);
  testrel(59, NE, y, f, cnt,failed);

  f := 0.58772179997748625915e86;
  y := transportx(100,x);
  testrel(60, NE+1, y, f, cnt,failed);  {+1 for FPC <= 2.4.0}

  {---------------------}
  x := 10.0;
  f := 0.32843291144979517358e+01;
  y := transportx(2,x);
  testrel(61, NE, y, f, cnt,failed);  

  f := 0.11638074540242071077e+03;
  y := transportx(5,x);
  testrel(62, NE, y, f, cnt,failed);

  f := 0.19743173017140591390e+06;
  y := transportx(9,x);
  testrel(63, NE, y, f, cnt,failed);

  f := 0.49835724908069047812e95;
  y := transportx(100,x);
  testrel(64, NE, y, f, cnt,failed);

  {---------------------}
  x := 15.0;
  f := 0.32897895167775788137e+01;
  y := transportx(2,x);
  testrel(65, NE, y, f, cnt,failed);

  f := 0.12409623901262967878e+03;
  y := transportx(5,x);
  testrel(66, NE, y, f, cnt,failed);

  f := 0.33826030414658460679e+06;
  y := transportx(9,x);
  testrel(67, NE, y, f, cnt,failed);

  f := 0.2164887690008280379e111;
  y := transportx(100,x);
  testrel(68, NE, y, f, cnt,failed);

  {---------------------}
  x := 20.0;
  f := 0.32898672226665499687e+01;
  y := transportx(2,x);
  testrel(69, NE, y, f, cnt,failed);

  f := 0.12442270155632550228e+03;
  y := transportx(5,x);
  testrel(70, NE1, y, f, cnt,failed);

  f := 0.36179607036750755227e+06;
  y := transportx(9,x);
  testrel(71, NE, y, f, cnt,failed);

  f := 0.6432158007007351265e121;
  y := transportx(100,x);
  testrel(72, NE, y, f, cnt,failed);

  {---------------------}
  x := 30.0;
  f := 0.32898681336064325400e+01;
  y := transportx(2,x);
  testrel(73, NE, y, f, cnt,failed);

  f := 0.12443132790838589548e+03;
  y := transportx(5,x);
  testrel(74, NE1, y, f, cnt,failed);

  f := 0.36360622124777561525e+06;
  y := transportx(9,x);
  testrel(75, NE, y, f, cnt,failed);

  f := 0.2026006520981163383e135;
  y := transportx(100,x);
  testrel(76, NE, y, f, cnt,failed);

  {---------------------}
  x := 50.0;
  f := 0.32898681336964528724e+01;
  y := transportx(2,x);
  testrel(77, NE, y, f, cnt,failed);

  f := 0.12443133061720432435e+03;
  y := transportx(5,x);
  testrel(78, NE, y, f, cnt,failed);

  f := 0.36360880558827162725e+06;
  y := transportx(9,x);
  testrel(79, NE, y, f, cnt,failed);

  f := 0.1464984508159559704e149;
  y := transportx(100,x);
  testrel(80, NE, y, f, cnt,failed);

  {------  other ---------------}
  y := transportx(2,4);
  f := 2.806666404563117993;
  testrel(81, NE, y, f, cnt,failed);

  y := transportx(3,4);
  f := 4.579217437229156370;
  testrel(82, NE, y, f, cnt,failed);

  y := transportx(4,4);
  f := 10.73193239299862222;
  testrel(83, NE, y, f, cnt,failed);

  y := transportx(500,4);
  f := 0.1639463787743363691e298;
  testrel(84, NE1, y, f, cnt,failed);

  y := transportx(4,1e-6);
  f := 0.3333333333333166667e-18;
  testrel(85, NE, y, f, cnt,failed);

  {-------------------------------------------------------------------}

  f := 0.5018512214525259702e1200;
  y := transportx(2000,4);
  testrel(86, NE1, y, f, cnt,failed);

  {logarithm/overflow branch in integrand function}
  y := transportx(4100,16);
  f := 3.436725002515484781e4927; {Wolfram Alpha}
  testrel(87, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_rrcfx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','rrcfx');

  {RRR[q_]:=Sign[q]Abs[q]^(1/5)QPochhammer[q, q^5]QPochhammer[q^4, q^5]/QPochhammer[q^2, q^5]/QPochhammer[q^3, q^5]}

  y := rrcfx(9/10);
  f := 0.6180339887498948482;
  testrel(1, NE, y, f, cnt,failed);

  y := rrcfx(8/10);
  f := 0.6180339887498942546;
  testrel(2, NE, y, f, cnt,failed);

  y := rrcfx(3/5);
  f := 0.6180337209974445856;
  testrel(3, NE, y, f, cnt,failed);

  y := rrcfx(5/10);
  f := 0.6180183781955219566;
  testrel(4, NE, y, f, cnt,failed);

  y := rrcfx(1/4);
  f := 0.6133988985778487247;
  testrel(5, NE, y, f, cnt,failed);

  y := rrcfx(1/8);
  f := 0.5874502218013974770;
  testrel(6, NE, y, f, cnt,failed);

  y := rrcfx(1/16);
  f := 0.5406876571469069971;
  testrel(7, NE, y, f, cnt,failed);

  y := rrcfx(1/100);
  f := 0.3941659056233625130;
  testrel(8, NE, y, f, cnt,failed);

  y := rrcfx(1/1024);
  f := 0.2497560977933519497;
  testrel(9, NE, y, f, cnt,failed);

  y := rrcfx(1/65536.0);
  f := 0.1088171599939251764;
  testrel(10, NE1, y, f, cnt,failed);

  y := rrcfx(1e-12);
  f := 0.003981071705530991436;
  testrel(11, NE, y, f, cnt,failed);

  y := rrcfx(0);
  f := 0;
  testrel(12, NE, y, f, cnt,failed);

  y := rrcfx(-1e-10);
  f := -0.01000000000100000000;
  testrel(13, NE, y, f, cnt,failed);

  y := rrcfx(-1/65536.0);
  f := -0.1088204808807785704;
  testrel(14, NE, y, f, cnt,failed);

  y := rrcfx(-1/1024);
  f := -0.2502443790433515056;
  testrel(15, NE, y, f, cnt,failed);

  y := rrcfx(-1/100);
  f := -0.4021280489548030626;
  testrel(16, NE, y, f, cnt,failed);

  y := rrcfx(-1/32);
  f := -0.5161127890314825420;
  testrel(17, NE, y, f, cnt,failed);

  y := rrcfx(-1/16);
  f := -0.6124802045877800335;
  testrel(18, NE, y, f, cnt,failed);

  y := rrcfx(-1/8);
  f := -0.7523478235636344042;
  testrel(19, NE, y, f, cnt,failed);

  y := rrcfx(-1/4);
{$ifdef BIT16}
  f := -0.4953810982405444272*2;
{$else}
  f := -0.9907621964810888544;
{$endif}
  testrel(20, NE, y, f, cnt,failed);

  y := rrcfx(-1/2);
  f := -1.426271071713779222;
  testrel(21, NE, y, f, cnt,failed);

  y := rrcfx(-3/4);
  f := -1.614250954566703945;
  testrel(22, NE1, y, f, cnt,failed);

  y := rrcfx(-8/10);
  f := -1.617513254575678072;
  testrel(23, NE, y, f, cnt,failed);

  y := rrcfx(-9/10);
  f := -1.618033962325803401;
  testrel(24, NE, y, f, cnt,failed);

  y := rrcfx(-99/100);
  f := -1.618033988749894848;
  testrel(25, NE, y, f, cnt,failed);

  y := rrcfx(-1 + 1/512);
  f := -1.618033988749894848;
  testrel(26, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_eulerpolyx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 8;
  NE2 = 80;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','eulerpolyx');

  y := eulerpolyx(1,-1);
  f := -3/2;
  testrel(1, NE, y, f, cnt,failed);

  y := eulerpolyx(7,-1);
  f := -4.125;
  testrel(2, NE, y, f, cnt,failed);

  y := eulerpolyx(8,-1);
  f := 2;
  testrel(3, NE, y, f, cnt,failed);

  y := eulerpolyx(1,1);
  f := 0.5;
  testrel(4, NE, y, f, cnt,failed);

  y := eulerpolyx(7,1);
  f := -17/8;
  testrel(5, NE, y, f, cnt,failed);

  y := eulerpolyx(8,1);
  f := 0;
  testrel(6, NE, y, f, cnt,failed);

  y := eulerpolyx(10,-0.75);
  f := 34.99948596954345703;
  testrel(7, NE, y, f, cnt,failed);

  y := eulerpolyx(12,1/3);
  f := 571.4513219717710903;
  testrel(8, NE, y, f, cnt,failed);

  y := eulerpolyx(11,-0.125);
  f := 159.6000123716657981;
  testrel(9, NE, y, f, cnt,failed);

  y := eulerpolyx(20, 0.375);
  f := 326326713.8134914185;
  testrel(10, NE, y, f, cnt,failed);

  y := eulerpolyx(15,-1.25);
  f := -41138.37547717522830;
  testrel(11, NE1, y, f, cnt,failed);

  y := eulerpolyx(12,2.5);
  f := 919.346923828125;
  testrel(12, NE, y, f, cnt,failed);

  y := eulerpolyx(7,1e-5);
  f := 2.12499999895;
  testrel(13, NE, y, f, cnt,failed);

  y := eulerpolyx(6,1e-5);
  f := -0.29999999995e-4;
  testrel(14, NE, y, f, cnt,failed);

  y := eulerpolyx(6,0.5e-5);
  f := -0.14999999999375e-4;
  testrel(15, NE, y, f, cnt,failed);

  y := eulerpolyx(2,1e-11);
  f := -0.99999999999e-11;
  testrel(16, NE, y, f, cnt,failed);

  y := eulerpolyx(3,1e-11);
  f := 0.25;
  testrel(17, NE, y, f, cnt,failed);

  y := eulerpolyx(4,1e-11);
  f := 1e-11;
  testrel(18, NE, y, f, cnt,failed);

  y := eulerpolyx(5,1e-11);
  f := -0.5;
  testrel(19, NE, y, f, cnt,failed);

  y := eulerpolyx(15,1e-12);
  f := 58098.0625;
  testrel(20, NE, y, f, cnt,failed);

  y := eulerpolyx(15,1e-10);
  f := 58098.0625;
  testrel(21, NE, y, f, cnt,failed);

  y := eulerpolyx(20,1e-10);
  f := 0.1109652905;
  testrel(22, NE, y, f, cnt,failed);

  y := eulerpolyx(51,2e-10);
  f := 0.8727938146243366230e41;
  testrel(23, NE, y, f, cnt,failed);

  y := eulerpolyx(51,1e-5);
  f := 0.8727938141936301395e41;
  testrel(24, NE, y, f, cnt,failed);

  y := eulerpolyx(101,1e-5);
  f := -0.7363732519859706166e110;
  testrel(25, NE, y, f, cnt,failed);

  y := eulerpolyx(16,3e5);
  f := 0.4304557308744022321e88;
  testrel(26, NE1, y, f, cnt,failed);

  y := eulerpolyx(15,10000+1/3);
  f := 0.9997497667214722518e60;
  testrel(27, NE1, y, f, cnt,failed);

  y := eulerpolyx(10,10000+1/3);
  f := 0.9998332333477834809e40;
  testrel(28, NE1, y, f, cnt,failed);

  y := eulerpolyx(15,sqrt(1e21));
  f := 0.3162277659418379332e158;
  testrel(29, NE1, y, f, cnt,failed);

  y := eulerpolyx(10,sqrt(1e21));
  f := 0.9999999998418861170e105;
  testrel(30, NE1, y, f, cnt,failed);

  y := eulerpolyx(100,0.5);
  f := 0.2290479999881941142e109;
  testrel(31, NE+1, y, f, cnt,failed);    {+1 FPC 2.4/6}

  y := eulerpolyx(70,-987654321.0/65536.0);
  f := 0.2949579273885267560e293;
  testrel(32, NE1, y, f, cnt,failed);

  y := eulerpolyx(2,-100);
  f := 10100.00;
  testrel(33, NE, y, f, cnt,failed);

  y := eulerpolyx(15,-Pi);
  f := -57185905.81470278369;
  testrel(34, NE1, y, f, cnt,failed);

  y := eulerpolyx(20,-3.25);
  f := 0.3479747288397440296e11;
  testrel(35, NE1, y, f, cnt,failed);

  y := eulerpolyx(15,-3.25);
  f := -0.9499879158876228612e8;
  testrel(36, NE1, y, f, cnt,failed);

  y := eulerpolyx(68,0.25);
  f := 0.3488684394193107058e63;
  testrel(37, NE+1, y, f, cnt,failed);    {+1 FPC 2.4/6}

  y := eulerpolyx(70,2/3);
  f := -0.2091003199868267386e66;
  testrel(38, NE1, y, f, cnt,failed);

  y := eulerpolyx(75,-1.75);
  f := 0.1155490536153012231e73;
  testrel(39, NE1, y, f, cnt,failed);

  y := eulerpolyx(120,4.75);
  f := 0.1323744731796364147e140;
  testrel(40, NE2, y, f, cnt,failed);

  y := eulerpolyx(200,4.75);
  f := 2.638209148556641986e275;
  testrel(41, NE2, y, f, cnt,failed);

  y := eulerpolyx(10,-1/3);
  f := 42.72729428102084709;
  testrel(42, NE, y, f, cnt,failed);

  y := eulerpolyx(50,-1/3);
  f := 0.4656096664921116345e40;
  testrel(43, NE, y, f, cnt,failed);

  y := eulerpolyx(100,-1/3);
  f := -0.1983613866757939061e109;
  testrel(44, NE, y, f, cnt,failed);

  y := eulerpolyx(200,-1/3);
  f := -0.3231133124353118137e276;
  testrel(45, NE1, y, f, cnt,failed);

  {extended only}

  y := eulerpolyx(256,0.5);
  f := 5.860541561523165577e379;
  testrel(46, NE, y, f, cnt,failed);

  y := eulerpolyx(500,-100);
  f := 1.986940598312024748e1000;
  testrel(47, NE1, y, f, cnt,failed);

  y := eulerpolyx(500,-100.25);
  f := 6.923721368817791394e1000;
  testrel(48, NE2, y, f, cnt,failed);

  y := eulerpolyx(1000,0.5);
  f := 3.628120031122876848e2070;
  testrel(49, NE, y, f, cnt,failed);

  y := eulerpolyx(2000,1.25);
  f := -1.497270546998412557e4741;
  testrel(50, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_besselpolyx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 5;
  NE2 = 10;
  NE3 = 80;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','besselpolyx');

  y := besselpolyx(5, 3.25);
  f := 463655.2041015625;
  testrel(1, NE, y, f, cnt,failed);

  y := besselpolyx(5, -3.25);
  f := -250574.5712890625;
  testrel(2, NE, y, f, cnt,failed);

  y := besselpolyx(5, 1.2345);
  f := 5874.004716610139213;
  testrel(3, NE, y, f, cnt,failed);

  y := besselpolyx(5, -1.2345);
  f := -1162.354290950471088;
  testrel(4, NE1, y, f, cnt,failed);

  y := besselpolyx(5, 2/3);
  f := 493.2222222222222222;
  testrel(5, NE, y, f, cnt,failed);

  y := besselpolyx(5, -2/3);
  f := -24.55555555555555556;
  testrel(6, NE, y, f, cnt,failed);

  y := besselpolyx(6, Pi);
  f := 13676180.74468108377;
  testrel(7, NE, y, f, cnt,failed);

  y := besselpolyx(6, -Pi);
  f := 7235763.733869276038;
  testrel(8, NE, y, f, cnt,failed);

  y := besselpolyx(10, 1);
  f := 1733584106.0;
  testrel(9, NE, y, f, cnt,failed);

  y := besselpolyx(10, -1);
  f := 234615096.0;
  testrel(10, NE, y, f, cnt,failed);

  y := besselpolyx(10, 1/3);
  f := 176298.9355281207133;
  testrel(11, NE1, y, f, cnt,failed);

  y := besselpolyx(10, -1/3);
  f := 437.0013717421124829;
  testrel(12, NE3, y, f, cnt,failed);

  {negative n}
  {Maple: yu := (n,x) -> (2/x)^(n+1)*KummerU(n+1, 2*n+2, 2/x);}
  y := besselpolyx(-4, -0.5);
  f := -0.125;
  testrel(13, NE, y, f, cnt,failed);

  y := besselpolyx(-10, -0.5);
  f := -8104.888671875;
  testrel(14, NE, y, f, cnt,failed);

  f := 28875761731.0;
  y := besselpolyx(-10,2);
  testrel(15, NE, y, f, cnt,failed);

  y := besselpolyx(-15, -1.5);
  f := 0.3173125422610212804e17;
  testrel(16, NE, y, f, cnt,failed);

  {near zero}
  y := besselpolyx(3, -0.4306640625);
  f := -0.5022343248128890991e-4;
  testrel(17, NE, y, f, cnt,failed);

  y := besselpolyx(11,-0.130859375);
  f := 0.2513213762353303456e-4;
  testabs(18, 1, y, f, cnt,failed);

  {large double}
  y := besselpolyx(100,2);
  f := 0.1392385298391907837e218;
  testrel(19, NE2, y, f, cnt,failed);

  y := besselpolyx(123,-3);
  f := -0.2281475467071922978e300;
  testrel(20, NE, y, f, cnt,failed);

  y := besselpolyx(500,1/256);
  f := 0.1775614659407953041e177;
  testrel(21, NE1, y, f, cnt,failed);

  {extended}
  y := besselpolyx(500,1);
  f := 0.2737253201635514604e1284;
  testrel(22, NE1, y, f, cnt,failed);  {FPC < 2.6}

  f := 0.2519836903971480732e1633;
  y := besselpolyx(500,-5);
  testrel(23, NE1, y, f, cnt,failed);

  f := 0.3759154927372568952e1633;
  y := besselpolyx(500,5);
  testrel(24, NE, y, f, cnt,failed);

  f := 0.1007584050999813734e3284;
  y := besselpolyx(500,10000);
  testrel(25, NE2, y, f, cnt,failed);

  f := 0.1007382554339951413e3284;
  y := besselpolyx(500,-10000);
  testrel(26, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_expnx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 5;
  NE2 = 40;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','expnx');

  y := expnx(0,0);
  f := 1;
  testrel(1, NE, y, f, cnt,failed);

  y := expnx(0,-2);
  f := 1;
  testrel(2, NE, y, f, cnt,failed);

  y := expnx(1,-1);
  f := 0;
  testrel(3, NE, y, f, cnt,failed);

  y := expnx(1,1);
  f := 2;
  testrel(4, NE, y, f, cnt,failed);

  y := expnx(2,-1);
  f := 0.5;
  testrel(5, NE, y, f, cnt,failed);

  y := expnx(2,1);
  f := 2.5;
  testrel(6, NE, y, f, cnt,failed);

  y := expnx(3,-1.596071637983321523);
  f := 0;
  testrel(7, NE, y, f, cnt,failed);

  y := expnx(3,-2);
  f := -1/3;
  testrel(8, NE, y, f, cnt,failed);

  y := expnx(4,-2);
  f := 1/3;
  testrel(9, NE, y, f, cnt,failed);

  y := expnx(5, -0.1);
  f := 0.9048374166666666667;
  testrel(10, NE, y, f, cnt,failed);

  y := expnx(5, 0.1);
  f := 1.105170916666666667;
  testrel(11, NE1, y, f, cnt,failed);

  y := expnx(5,-3);
  f := -0.65;
  testrel(12, NE, y, f, cnt,failed);

  y := expnx(5,3);
  f := 18.4;
  testrel(13, NE, y, f, cnt,failed);

  y := expnx(5,-100);
  f := -7.932843233333333333e7;
  testrel(14, NE, y, f, cnt,failed);

  y := expnx(5,100);
  f := 8.767176766666666667e7;
  testrel(15, NE, y, f, cnt,failed);

  y := expnx(10,-20);
  f := 1859623.680776014109;
  testrel(16, NE, y, f, cnt,failed);

  y := expnx(10,20);
  f := 5245469.677248677249;
  testrel(17, NE, y, f, cnt,failed);

  y := expnx(10,-500);
  f := 2.638275091805486254e20;
  testrel(18, NE, y, f, cnt,failed);

  y := expnx(10,500);
  f := 2.745951877217058133e20;
  testrel(19, NE, y, f, cnt,failed);

  y := expnx(42,-100);
  f := 0.5001735594488707970e33;
  testrel(20, NE, y, f, cnt,failed);

  y := expnx(42,100);
  f := 0.1212808667545658566e34;
  testrel(21, NE1, y, f, cnt,failed);

  y := expnx(170,500);
  f := 0.1392847236069440217e153;
  testrel(22, NE2, y, f, cnt,failed);

  y := expnx(170,-40);
  f := 4.248354255291589001e-18;
  testrel(23, NE, y, f, cnt,failed);


  {AMath/double}
  y := expnx(123,1234.5);
  f := 0.1638697957253044953e176;
  testrel(24, NE2, y, f, cnt,failed);

  y := expnx(123,-200);
  f := -0.5414509219188112089e78;
  testrel(25, NE, y, f, cnt,failed);

  y := expnx(1600,-100);
  f := 0.3720075976020835963e-43;
  testrel(26, NE, y, f, cnt,failed);

  y := expnx(200,-75);
  f := 0.3547670769893360447;
  testrel(27, NE1, y, f, cnt,failed);

  y := expnx(1754,-75);
  f := 0.2678636961808077944e-32;
  testrel(28, NE+1, y, f, cnt,failed);   {+1 FPC<3}

  {AMath/extended}
  y := expnx(1700,1234.5);
  f := 1.369423920185182535e536;
  testrel(29, NE, y, f, cnt,failed);

  y := expnx(1700,-500);
  f := 1.344371686239870847e-168;
  testrel(30, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.
