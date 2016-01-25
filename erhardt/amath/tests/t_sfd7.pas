{Part 7 of regression test for SPECFUN unit  (c) 2010  W.Ehrhardt}

unit t_sfd7;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_agm;
procedure test_bernpoly;
procedure test_lambertw;
procedure test_lambertw1;
procedure test_debye;
procedure test_li_inv;
procedure test_RiemannR;
procedure test_cosint;
procedure test_sinint;
procedure test_fibpoly;
procedure test_lucpoly;
procedure test_catalan;
procedure test_kepler;


implementation

uses
  AMath, SpecFun, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_agm;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','agm');

  y := agm(MinDouble, MaxDouble);
  f := 1.9892072050015465487e305;
  testrel( 1, 2, y, f, cnt,failed);

  y := agm(1, MaxDouble);
  f := 3.9706592162214106567e305;
  testrel( 2, 2, y, f, cnt,failed);

  y := agm(MinDouble, 1);
  f := 0.0022130664755015594768;
  testrel( 3, NE, y, f, cnt,failed);

  y := agm(0, 1);
  f := 0;
  testabs( 4, 0, y, f, cnt,failed);

  y := agm(1, 0);
  f := 0;
  testabs( 5, 0, y, f, cnt,failed);

  y := agm(1, 1e-100);
  f := 0.0067810557455754508824;
  testrel( 5, NE, y, f, cnt,failed);

  y := agm(1, 1e-10);
  f := 0.064344870476013322929;
  testrel( 7, NE, y, f, cnt,failed);

  y := agm(1, 1e-5);
  f := 0.12177452186538904490;
  testrel( 8, NE, y, f, cnt,failed);

  y := agm(1, 0.125);
  f := 0.45196952219967034359;
  testrel( 9, NE, y, f, cnt,failed);

  y := agm(1, 0.5);
  f := 0.72839551552345343459;
  testrel(10, NE, y, f, cnt,failed);

  y := agm(1, 1);
  f := 1;
  testrel(11, 0, y, f, cnt,failed);

  y := agm(1, 1000);
  f := 189.38830240995087556;
  testrel(12, NE, y, f, cnt,failed);

  y := agm(1, 1e20);
  f := 3311261967046375735.6;
  testrel(13, NE, y, f, cnt,failed);

  y := agm(1, 1e200);
  f := 3.4007037462646779626e197;
  testrel(14, NE, y, f, cnt,failed);

  y := agm(1, 1e300);
  f := 2.2694061941578213057e297;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lambertw;
var
  x,y,f: double;
  cnt, failed: integer;
  i: integer;
const
  NE = 1;
  NA = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LambertW');

  x := -exp(-1);
  y := LambertW(x);
  f := -1;
  testrel( 1, NE, y, f, cnt,failed);

  x := -6027/16384.0;
  y := LambertW(x);
  f := -0.9894660903569092419;
  testrel( 2, NE, y, f, cnt,failed);

  x := -3013/8192.0;
  y := LambertW(x);
  f := -0.9790854111351907084;
  testrel( 3, NE, y, f, cnt,failed);

  x := -753/2048;
  y := LambertW(x);
  f := -0.9670887700916448631;
  testrel( 4, NE, y, f, cnt,failed);

  x := -0.25;
  y := LambertW(x);
  f := -0.3574029561813889031;
  testrel( 5, NE, y, f, cnt,failed);

  x := -0.125;
  y := LambertW(x);
  f := -0.1444213531375097292;
  testrel( 6, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := LambertW(x);
  f := -0.977517573730222695e-3;
  testrel( 7, NE, y, f, cnt,failed);

  x := -1e-10;
  y := LambertW(x);
  f := -0.1000000000100000000e-9;
  testrel( 8, NE, y, f, cnt,failed);

  x := 0.0;
  y := LambertW(x);
  f := 0;
  testrel( 9, 1, y, f, cnt,failed);

  x := 1e-10;
  y := LambertW(x);
  f := 0.9999999999000000000e-10;
  testrel(10, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := LambertW(x);
  f := 0.9756102202467530500e-3;
  testrel(11, NE, y, f, cnt,failed);

  x := 0.125;
  y := LambertW(x);
  f := 0.1117801089327885068;
  testrel(12, NE, y, f, cnt,failed);

  x := 0.25;
  y := LambertW(x);
  f := 0.2038883547022401644;
  testrel(13, NE, y, f, cnt,failed);

  x := 1;
  y := LambertW(x);
  f := +0.5671432904097838730;
  testrel(14, NE, y, f, cnt,failed);

  x := 3.125;
  y := LambertW(x);
  f := 1.070918030310010008;
  testrel(15, NE, y, f, cnt,failed);

  x := 10;
  y := LambertW(x);
  f := 1.745528002740699383;
  testrel(16, NE, y, f, cnt,failed);

  x := 100;
  y := LambertW(x);
  f := 3.385630140290050185;
  testrel(17, NE, y, f, cnt,failed);

  x := 1e4;
  y := LambertW(x);
  f := 7.231846038093372706;
  testrel(18, NE, y, f, cnt,failed);

  x := 1e10;
  y := LambertW(x);
  f := 20.02868541330495078;
  testrel(19, NE, y, f, cnt,failed);

  x := MaxDouble;  {2^1024}
  y := LambertW(x);
  f := 703.2270331047701870;
  testrel(20, NE, y, f, cnt,failed);

  x := 1e-300;
  y := LambertW(x);
  f := x;
  testrel(21, NE, y, f, cnt,failed);

  x := MaxDouble;
  i := 1000;
  while x>1e-200 do begin
    y := LambertW(x);
    f := ln(x/y);
    if y>=0.5 then testrel(i, 2, y, f, cnt,failed)
    else testabs(i, 1, y, f, cnt,failed);
    x := 0.125*x;
    inc(i);
  end;

  x := -exp(-1);
  i := 10000;
  while x<-MinDouble do begin
    y := LambertW(x);
    f := ln(x/y);
    if y>=0.5 then testrel(i, 2, y, f, cnt,failed)
    else testabs(i, NA, y, f, cnt,failed);
    if x < -1e-7 then x := 0.95*x
    else x := 0.01*x;
    inc(i);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lambertw1;
var
  x,y,f: double;
  cnt, failed: integer;
  i: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LambertW1');

  x := -exp(-1);
  y := LambertW1(x);
  f := -1;
  testrel( 1, NE, y, f, cnt,failed);

  x := -6027/16384.0;
  y := LambertW1(x);
  f := -1.010608408692175856;
  testrel( 2, NE, y, f, cnt,failed);

  x := -3013/8192.0;
  y := LambertW1(x);
  f := -1.021210331605726089434;
  testrel( 3, NE, y, f, cnt,failed);

  x := -753/2048;
  y := LambertW1(x);
  f := -1.033649565301978476;
  testrel( 4, NE, y, f, cnt,failed);

  x := -0.25;
  y := LambertW1(x);
  f := -2.153292364110349649;
  testrel( 5, NE, y, f, cnt,failed);

  x := -0.125;
  y := LambertW1(x);
  f := -3.261685684576488777;
  testrel( 6, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := LambertW1(x);
  f := -9.144639686625083192;
  testrel( 7, NE, y, f, cnt,failed);

  x := -1e-10;
  y := LambertW1(x);
  f := -26.29523881924692569;
  testrel( 8, NE, y, f, cnt,failed);

  x := -1e-100;
  y := LambertW1(x);
  f := -235.7211588756853137;
  testrel( 9, NE, y, f, cnt,failed);

  x := -1e-299;
  y := LambertW1(x);
  f := -695.0168789367294442;
  testrel(10, NE, y, f, cnt,failed);

  x := -exp(-1);
  i := 1000;
  while x<-MinDouble do begin
    y := LambertW1(x);
    f := ln(x/y);
    testrel(i, 2, y, f, cnt,failed);
    if x < -1e-7 then x := 0.95*x
    else
      x := 0.01*x;
    inc(i);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_debye;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','debye');

  {Test values debye(n,x), n=1,3,4 from MISCFUN [22]}
  x := 1/512;
  y := debye(1,x);
  f := 0.99951182471380889183;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1/32;
  y := debye(1,x);
  f := 0.99221462647120597836;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1/8;
  y := debye(1,x);
  f := 0.96918395997895308324;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1/2;
  y := debye(1,x);
  f := 0.88192715679060552968;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1.0;
  y := debye(1,x);
  f := 0.77750463411224827642;
  testrel( 5, NE, y, f, cnt,failed);

  x := 2.0;
  y := debye(1,x);
  f := 0.60694728460981007205;
  testrel( 6, NE, y, f, cnt,failed);

  x := 3.0;
  y := debye(1,x);
  f := 0.48043521957304283829;
  testrel( 7, NE, y, f, cnt,failed);

  x := 4.0;
  y := debye(1,x);
  f := 0.38814802129793784501;
  testrel( 8, NE, y, f, cnt,failed);

  x := 17/4;
  y := debye(1,x);
  f := 0.36930802829242526815;
  testrel( 9, NE, y, f, cnt,failed);

  x := 5.0;
  y := debye(1,x);
  f := 0.32087619770014612104;
  testrel(10, NE, y, f, cnt,failed);

  x := 11/2;
  y := debye(1,x);
  f := 0.29423996623154246701;
  testrel(11, NE, y, f, cnt,failed);

  x := 6.0;
  y := debye(1,x);
  f := 0.27126046678502189985;
  testrel(12, NE, y, f, cnt,failed);

  x := 8.0;
  y := debye(1,x);
  f := 0.20523930310221503723;
  testrel(13, NE, y, f, cnt,failed);

  x := 10.0;
  y := debye(1,x);
  f := 0.16444346567994602563;
  testrel(14, NE, y, f, cnt,failed);

  x := 20.0;
  y := debye(1,x);
  f := 0.82246701178200016086e-1;
  testrel(15, NE, y, f, cnt,failed);

  x := 50.0;
  y := debye(1,x);
  f := 0.32898681336964528729e-1;
  testrel(16, NE, y, f, cnt,failed);

  x := 1/512;
  y := debye(3,x);
  f := 0.99926776885985461940;
  testrel(17, NE, y, f, cnt,failed);

  x := 1/32;
  y := debye(3,x);
  f := 0.98833007755734698212;
  testrel(18, NE, y, f, cnt,failed);

  x := 1/8;
  y := debye(3,x);
  f := 0.95390610472023510237;
  testrel(19, NE, y, f, cnt,failed);

  x := 1/2;
  y := debye(3,x);
  f := 0.82496296897623372315;
  testrel(20, NE, y, f, cnt,failed);

  x := 1.0;
  y := debye(3,x);
  f := 0.67441556407781468010;
  testrel(21, NE, y, f, cnt,failed);

  x := 2.0;
  y := debye(3,x);
  f := 0.44112847372762418113;
  testrel(22, NE, y, f, cnt,failed);

  x := 3.0;
  y := debye(3,x);
  f := 0.28357982814342246206;
  testrel(23, NE, y, f, cnt,failed);

  x := 4.0;
  y := debye(3,x);
  f := 0.18173691382177474795;
  testrel(24, NE, y, f, cnt,failed);

  x := 17/4;
  y := debye(3,x);
  f := 0.16277924385112436877;
  testrel(25, NE, y, f, cnt,failed);

  x := 5.0;
  y := debye(3,x);
  f := 0.11759741179993396450;
  testrel(26, NE, y, f, cnt,failed);

  x := 11/2;
  y := debye(3,x);
  f := 0.95240802723158889887e-1;
  testrel(27, NE, y, f, cnt,failed);

  x := 6.0;
  y := debye(3,x);
  f := 0.77581324733763020269e-1;
  testrel(28, NE, y, f, cnt,failed);

  x := 8.0;
  y := debye(3,x);
  f := 0.36560295673194845002e-1;
  testrel(29, NE, y, f, cnt,failed);

  x := 10.0;
  y := debye(3,x);
  f := 0.19295765690345489563e-1;
  testrel(30, NE, y, f, cnt,failed);

  x := 20.0;
  y := debye(3,x);
  f := 0.24352200674805479827e-2;
  testrel(31, NE, y, f, cnt,failed);

  x := 50.0;
  y := debye(3,x);
  f := 0.15585454565440389896e-3;
  testrel(32, NE, y, f, cnt,failed);

  x := 1/512;
  y := debye(4,x);
  f := 0.99921896192761576256;
  testrel(33, NE, y, f, cnt,failed);

  x := 1/32;
  y := debye(4,x);
  f := 0.98755425280996071022;
  testrel(34, NE, y, f, cnt,failed);

  x := 1/8;
  y := debye(4,x);
  f := 0.95086788606389739976;
  testrel(35, NE, y, f, cnt,failed);

  x := 1/2;
  y := debye(4,x);
  f := 0.81384569172034042516;
  testrel(36, NE, y, f, cnt,failed);

  x := 1.0;
  y := debye(4,x);
  f := 0.65487406888673697092;
  testrel(37, NE, y, f, cnt,failed);

  x := 2.0;
  y := debye(4,x);
  f := 0.41189273671788528876;
  testrel(38, NE, y, f, cnt,failed);

  x := 3.0;
  y := debye(4,x);
  f := 0.25187863642883314410;
  testrel(39, NE, y, f, cnt,failed);

  x := 4.0;
  y := debye(4,x);
  f := 0.15185461258672022043;
  testrel(40, NE, y, f, cnt,failed);

  x := 17/4;
  y := debye(4,x);
  f := 0.13372661145921413299;
  testrel(41, NE, y, f, cnt,failed);

  x := 5.0;
  y := debye(4,x);
  f := 0.91471377664481164749e-1;
  testrel(42, NE, y, f, cnt,failed);

  x := 11/2;
  y := debye(4,x);
  f := 0.71227828197462523663e-1;
  testrel(43, NE, y, f, cnt,failed);

  x := 6.0;
  y := debye(4,x);
  f := 0.55676547822738862783e-1;
  testrel(44, NE, y, f, cnt,failed);

  x := 8.0;
  y := debye(4,x);
  f := 0.21967566525574960096e-1;
  testrel(45, NE, y, f, cnt,failed);

  x := 10.0;
  y := debye(4,x);
  f := 0.96736755602711590082e-2;
  testrel(46, NE, y, f, cnt,failed);

  x := 20.0;
  y := debye(4,x);
  f := 0.62214648623965450200e-3;
  testrel(47, NE, y, f, cnt,failed);

  x := 50.0;
  y := debye(4,x);
  f := 0.15927210319002161231e-4;
  testrel(48, NE, y, f, cnt,failed);

  {Rest of test values for n=6,8,12 calculated with Maple }
  {f := x->n*int(t^n/(exp(t)-1),t=0..x)/x^n; and Digits:=50}
  x := 1/512;
  y := debye(6,x);
  f := 0.9991631848471384035;
  testrel(49, NE, y, f, cnt,failed);

  x := 1/32;
  y := debye(6,x);
  f := 0.9866681772186796587;
  testrel(50, NE, y, f, cnt,failed);

  x := 1/8;
  y := debye(6,x);
  f := 0.9474049305411031823;
  testrel(51, NE, y, f, cnt,failed);

  x := 1/2;
  y := debye(6,x);
  f := 0.8012874593544054948;
  testrel(52, NE, y, f, cnt,failed);

  x := 1.0;
  y := debye(6,x);
  f := 0.6331114258349510759;
  testrel(53, NE, y, f, cnt,failed);

  x := 2.0;
  y := debye(6,x);
  f := 0.3804986630746610429;
  testrel(54, NE, y, f, cnt,failed);

  x := 3.0;
  y := debye(6,x);
  f := 0.2193992525257245836;
  testrel(55, NE, y, f, cnt,failed);

  x := 4.0;
  y := debye(6,x);
  f := 0.1229278562814578228;
  testrel(56, NE, y, f, cnt,failed);

  x := 17/4;
  y := debye(6,x);
  f := 0.1060375248597196031;
  testrel(57, NE, y, f, cnt,failed);

  x := 5.0;
  y := debye(6,x);
  f := 0.6777784974890353731e-1;
  testrel(58, NE, y, f, cnt,failed);

  x := 11/2;
  y := debye(6,x);
  f := 0.5020600934448088116e-1;
  testrel(59, NE, y, f, cnt,failed);

  x := 6.0;
  y := debye(6,x);
  f := 0.3719333613705515670e-1;
  testrel(60, NE, y, f, cnt,failed);

  x := 8.0;
  y := debye(6,x);
  f := 0.1145231921902748610e-1;
  testrel(61, NE, y, f, cnt,failed);

  x := 10.0;
  y := debye(6,x);
  f := 0.3793849329461595528e-2;
  testrel(62, NE, y, f, cnt,failed);

  x := 20.0;
  y := debye(6,x);
  f := 0.6804635545479456894e-4;
  testrel(63, NE, y, f, cnt,failed);

  x := 50.0;
  y := debye(6,x);
  f := 0.2787884082105527120e-6;
  testrel(64, NE, y, f, cnt,failed);

  x := 1/512;
  y := debye(12,x);
  f := 0.9990988301706686501;
  testrel(65, NE, y, f, cnt,failed);

  x := 1/32;
  y := debye(12,x);
  f := 0.9856466765478185763;
  testrel(66, NE, y, f, cnt,failed);

  x := 1/8;
  y := debye(12,x);
  f := 0.9434235095071814036;
  testrel(67, NE, y, f, cnt,failed);

  x := 1/2;
  y := debye(12,x);
  f := 0.7870231504611680153;
  testrel(68, NE, y, f, cnt,failed);

  x := 1.0;
  y := debye(12,x);
  f := 0.6088700041762235049;
  testrel(69, NE, y, f, cnt,failed);

  x := 2.0;
  y := debye(12,x);
  f := 0.3472653175019342084;
  testrel(70, NE, y, f, cnt,failed);

  x := 3.0;
  y := debye(12,x);
  f := 0.1872401059320712096;
  testrel(71, NE, y, f, cnt,failed);

  x := 4.0;
  y := debye(12,x);
  f := 0.9654143896262086549e-1;
  testrel(72, NE, y, f, cnt,failed);

  x := 17/4;
  y := debye(12,x);
  f := 0.8134774441706960165e-1;
  testrel(73, NE, y, f, cnt,failed);

  x := 5.0;
  y := debye(12,x);
  f := 0.4814185645191148541e-1;
  testrel(74, NE, y, f, cnt,failed);

  x := 11/2;
  y := debye(12,x);
  f := 0.3367880055948328374e-1;
  testrel(75, NE, y, f, cnt,failed);

  x := 6.0;
  y := debye(12,x);
  f := 0.2344811348723500784e-1;
  testrel(76, NE, y, f, cnt,failed);

  x := 8.0;
  y := debye(12,x);
  f := 0.5344588786833453221e-2;
  testrel(77, NE, y, f, cnt,failed);

  x := 10.0;
  y := debye(12,x);
  f := 0.1198815360618837356e-2;
  testrel(78, NE, y, f, cnt,failed);

  x := 20.0;
  y := debye(12,x);
  f := 0.1348750701799345211e-5;
  testrel(79, NE, y, f, cnt,failed);

  x := 50.0;
  y := debye(12,x);
  f := 0.2354677578932315411e-10;
  testrel(80, NE, y, f, cnt,failed);

  {D(n,20) n=20..200}
  x := 20;
  f := 0.2045985597891880435e-6;
  y := debye(20,x);
  testrel(81, NE, y, f, cnt,failed);

  f := 0.6521968411709099410e-7;
  y := debye(50,x);
  testrel(82, NE, y, f, cnt,failed);

  f := 0.5964824507515076368e-7;
  y := debye(60,x);
  testrel(83, NE, y, f, cnt,failed);

  f := 0.5378148028857494088e-7;
  y := debye(80,x);
  testrel(84, NE, y, f, cnt,failed);

  f := 0.5074077974146328750e-7;
  y := debye(100,x);
  testrel(85, NE, y, f, cnt,failed);

  f := 0.4714758392337898263e-7;
  y := debye(150,x);
  testrel(86, NE, y, f, cnt,failed);

  f := 0.4552275137444157216e-7;
  y := debye(200,x);
  testrel(87, NE, y, f, cnt,failed);


  {D(n,5) n=100..10000}
  x := 5;
  f := 0.3532560165653401053e-1;
  y := debye(100,x);
  testrel(88, NE, y, f, cnt,failed);

  f := 0.34193472837210103295e-1;
  y := debye(500,x);
  testrel(89, NE, y, f, cnt,failed);

  f := 0.3410139487334168171e-1;
  y := debye(750,x);
  testrel(90, NE, y, f, cnt,failed);

  f := 0.3405548547450571961e-1;
  y := debye(1000,x);
  testrel(91, NE, y, f, cnt,failed);

  f := 0.3398678310287432440e-1;
  y := debye(2000,x);
  testrel(92, NE, y, f, cnt,failed);

  f := 0.3393196075652582305e-1;
  y := debye(10000,x);
  testrel(93, NE, y, f, cnt,failed);

  {Test after fix for large x}
  y := debye(1, 100000.0);
  f := 0.1644934066848226436e-4;
  testrel(94, NE, y, f, cnt,failed);

  x := 20000.0;
  y := debye(2, x);
  f := 0.1202056903159594285e-7;
  testrel(95, NE, y, f, cnt,failed);

  y := debye(7, x);
  f := 0.2767488213020584085e-25;
  testrel(96, NE, y, f, cnt,failed);

  y := debye(10, x);
  f := 0.3545501280865848353e-35;
  testrel(97, NE, y, f, cnt,failed);

  y := debye(8, 30000.0);
  f := 0.4926197640450862355e-30;
  testrel(98, NE, y, f, cnt,failed);

  y := debye(20, 12000);
  f := 0.12691995186452421609e-61;
  testrel(99, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_RiemannR;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','RiemannR');

  y := RiemannR(ldexpd(1,1023));
  f := 1.269399252238138694e305;
  testrel( 1, NE, y, f, cnt,failed);

  y := RiemannR(1e100);
  f := 4.361971987140703159e97;
  testrel( 2, NE, y, f, cnt,failed);

  y := RiemannR(1e30);
  f := 0.1469239889772043272e29;
  testrel( 3, NE, y, f, cnt,failed);

  y := RiemannR(1e24);
  f := 0.18435599767347541878e23;
  testrel( 4, NE, y, f, cnt,failed);

  y := RiemannR(1e20);
  f := 0.2220819602556027015e19;
  testrel( 5, NE, y, f, cnt,failed);

  y := RiemannR(1e19);
  f := 0.2340576673002289402e18;
  testrel( 6, NE, y, f, cnt,failed);

  y := RiemannR(1e18);
  f := 0.2473995428423949440e17;
  testrel( 7, NE, y, f, cnt,failed);

  y := RiemannR(1e16);
  f := 0.2792383413609771872e15;
  testrel( 8, NE, y, f, cnt,failed);

  y := RiemannR(1e6);
  f := 78527.39942912770486;
  testrel( 9, NE, y, f, cnt,failed);

  y := RiemannR(1000);
  f := 168.3594462811673481;
  testrel(10, NE, y, f, cnt,failed);

  y := RiemannR(100);
  f := 25.66163326692418259;
  testrel(11, NE, y, f, cnt,failed);

  y := RiemannR(10);
  f := 4.564583141005090240;
  testrel(12, NE, y, f, cnt,failed);

  y := RiemannR(8);
  f := 3.901186044934149947;
  testrel(13, NE, y, f, cnt,failed);

  y := RiemannR(4);
  f := 2.426657752706807136;
  testrel(14, NE, y, f, cnt,failed);

  y := RiemannR(2);
  f := 1.541009016187131883;
  testrel(15, NE, y, f, cnt,failed);

  y := RiemannR(1);
  f := 1;
  testrel(16, NE, y, f, cnt,failed);

  y := RiemannR(0.5);
  f := 0.6635262381124574212;
  testrel(17, NE, y, f, cnt,failed);

  y := RiemannR(0.1);
  f := 0.2790883528560020161;
  testrel(18, NE, y, f, cnt,failed);

  y := RiemannR(0.125);
  f := 0.3124612249259569001;
  testrel(19, NE, y, f, cnt,failed);

  y := RiemannR(0.0625);
  f := 0.2216077332920197402;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_li_inv;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','li_inv');

  y := li_inv(0);
  f := 1.0 + 0.4513692348833810503;
  testrel(1, NE, y, f, cnt,failed);

  y := li_inv(0.25);
  f := 1.552800837242485188;
  testrel(2, NE, y, f, cnt,failed);

  y := li_inv(0.5);
  f := 1.671930573009875373;
  testrel(3, NE, y, f, cnt,failed);

  y := li_inv(0.75);
  f := 1.810255009236505581;
  testrel(4, NE, y, f, cnt,failed);

  y := li_inv(1);
  f := 1.969047489224750850;
  testrel(5, NE, y, f, cnt,failed);

  y := succd(1);
  y := li_inv(y);
  f := 1.969047489224751001;
  testrel(6, NE, y, f, cnt,failed);

  y := li_inv(2);
  f := 2.825187152005826843;
  testrel(7, NE, y, f, cnt,failed);

  y := li_inv(3);
  f := 4.045118486231030200;
  testrel(8, NE, y, f, cnt,failed);

  y := li_inv(3.5);
  f := 4.786319700881971309;
  testrel(9, NE, y, f, cnt,failed);

  y := li_inv(4);
  f := 5.609276693050890355;
  testrel(10, NE, y, f, cnt,failed);

  y := li_inv(5);
  f := 7.480870261577641432;
  testrel(11, NE, y, f, cnt,failed);

  y := li_inv(8);
  f := 14.58290311807629198;
  testrel(12, NE, y, f, cnt,failed);

  y := li_inv(10);
  f := 20.284365456596612497;
  testrel(13, NE, y, f, cnt,failed);

  y := li_inv(20);
  f := 56.07960987414566197;
  testrel(14, NE, y, f, cnt,failed);

  y := li_inv(100);
  f := 488.8719098528075319;
  testrel(15, NE, y, f, cnt,failed);

  y := li_inv(1000);
  f := 7762.986220174737687;
  testrel(16, NE, y, f, cnt,failed);

  y := li_inv(-0.25);
  f := 1.365970426374257461;
  testrel(17, NE, y, f, cnt,failed);

  y := li_inv(-0.5);
  f := 1.294838891062147533;
  testrel(18, NE, y, f, cnt,failed);

  y := li_inv(-0.75);
  f := 1.236183126594032207;
  testrel(19, NE, y, f, cnt,failed);

  y := li_inv(-1);
  f := 1.188256066274325355;
  testrel(20, NE, y, f, cnt,failed);

  y := li_inv(-10);
  f := 1.000025489896249632;
  testrel(21, NE, y, f, cnt,failed);

  y := li_inv(-15);
  f := 1.0 + 0.1717517441415356666e-6;
  testrel(22, NE, y, f, cnt,failed);

  y := li_inv(-40);
  f := 1.0 + 0.238528e-17;
  testrel(23, NE, y, f, cnt,failed);

  y := li_inv(-43.5);
  f := 1 + 0.72e-19;
  testrel(24, NE, y, f, cnt,failed);

  y := li_inv(1e10);
  f := 0.2520971600140342078e12;
  testrel(25, NE, y, f, cnt,failed);

  y := li_inv(1e100);
  f := 0.2347125735865764178e103;
  testrel(26, NE, y, f, cnt,failed);

  y := li_inv(1e300);
  f := 0.6963198968074983689e303;
  testrel(27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_bernpoly;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
  NE1 = 4;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bernpoly');

  y := bernpoly(1,-1);
  f := -3/2;
  testrel(1, NE, y, f, cnt,failed);

  y := bernpoly(7,-1);
  f := -7;
  testrel(2, NE, y, f, cnt,failed);

  y := bernpoly(8,-1);
  f := 239/30;
  testrel(3, NE, y, f, cnt,failed);

  y := bernpoly(1,1);
  f := 0.5;
  testrel(4, NE, y, f, cnt,failed);

  y := bernpoly(7,1);
  f := 0;
  testrel(5, NE, y, f, cnt,failed);

  y := bernpoly(8,1);
  f := -0.3333333333333333333e-1;
  testrel(6, NE, y, f, cnt,failed);

  {some 'normal' case}
  y := bernpoly(10,-0.75);
  f := 0.7507730252815015388;
  testrel(7, NE, y, f, cnt,failed);

  y := bernpoly(12,1/3);
  f := 0.1265560621400686431;
  testrel(8, NE, y, f, cnt,failed);

  y := bernpoly(11,-0.125);
  f := -0.9375500085297971964e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := bernpoly(20, 0.375);
  f := 374.14698287953266697;
  testrel(10, NE, y, f, cnt,failed);

  y := bernpoly(15,-1.25);
  f := -343.8455539522692561;
  testrel(11, NE, y, f, cnt,failed);

  y := bernpoly(12,2.5);
  f := 1038.229552462511447;
  testrel(12, NE, y, f, cnt,failed);


  x := 10000+1/3;
  y := bernpoly(15,x);
  f := 0.9997499416835207086e60;
  testrel(13, NE2, y, f, cnt,failed);

  y := bernpoly(10,x);
  f := 0.9998333083377781148e40;
  testrel(14, NE1, y, f, cnt,failed);

  x := sqrt(1e21);
  y := bernpoly(15,x);
  f := 0.3162277659418379332e158;
  testrel(15, NE1, y, f, cnt,failed);

  x := sqrt(1e21);
  y := bernpoly(10,x);
  f := 0.9999999998418861170e105;
  testrel(16, NE1, y, f, cnt,failed);

  y := bernpoly(256,0.5);
  f := 0.7950212504588525285e303;
  testrel(17, NE, y, f, cnt,failed);

  y := bernpoly(100,0.5);
  f := 2.838224957069370696e78;
  testrel(18, NE, y, f, cnt,failed);

  x := -987654321.0/65536.0;  {=15070.4089508056641}
  y := bernpoly(70,x);
  f := 0.2949584500818898822e293;
  testrel(19, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(7,x);
  f := 0.1666666665500000000e-5;
  testrel(20, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(6,x);
  f := 0.2380952375952380955e-1;
  testrel(21, NE, y, f, cnt,failed);

  x := 0.5e-5;
  y := bernpoly(6,x);
  f := 0.2380952379702380953e-1;
  testrel(22, NE, y, f, cnt,failed);

  x := 1e-11;
  y := bernpoly(2,x);
  f := 0.1666666666566666667;
  testrel(23, NE, y, f, cnt,failed);

  y := bernpoly(3,x);
  f := 0.4999999999850000000e-11;
  testrel(24, NE, y, f, cnt,failed);

  y := bernpoly(4,x);
  f := -0.3333333333333333333e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := bernpoly(5,x);
  f := -0.1666666666666666667e-11;
  testrel(26, NE, y, f, cnt,failed);

  x := 1e-12;
  y := bernpoly(15,x);
  f := 0.175e-10;
  testrel(27, NE, y, f, cnt,failed);

  x := 1e-10;
  y := bernpoly(15,x);
  f := 0.175e-8;
  testrel(28, NE, y, f, cnt,failed);

  x := 1e-10;
  y := bernpoly(20,x);
  f := -529.1242424242424241;
  testrel(29, NE, y, f, cnt,failed);

  x := 2e-10;
  y := bernpoly(51,x);
  f := 0.7650884080998503652e17;
  testrel(30, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(51,x);
  f := 0.3825442037982211854e22;
  testrel(31, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(101,x);
  f := -0.2866607204753912463e76;
  testrel(32, NE, y, f, cnt,failed);

  x := 3e5;
  y := bernpoly(16,x);
  f := 0.4304557309700593800e88;
  testrel(33, NE, y, f, cnt,failed);

  x := -100;
  y := bernpoly(2,x);
  f := 10100.16666666666667;
  testrel(34, NE, y, f, cnt,failed);

  x := -Pi;
  y := bernpoly(15,x);
  f := -0.1374730009236393778e9;
  testrel(35, NE1, y, f, cnt,failed);

  x := sqrt(10);
  y := bernpoly(20,x);
  f := 46168783.47767854148;
  testrel(36, NE2, y, f, cnt,failed);
  y := bernpoly(15,x);
  f := 732699.8879814299995;
  testrel(37, NE2, y, f, cnt,failed);

  x := 1/4;
  y := bernpoly(68,x);
  f := 0.8896458292761226510e22;
  testrel(38, NE, y, f, cnt,failed);

  x := 2/3;
  y := bernpoly(70,x);
  f := -0.1606254105135901626e45;
  testrel(39, NE1, y, f, cnt,failed);

  x := -1.75;
  y := bernpoly(75,x);
  f := 0.6794407537645821705e50;
  testrel(40, NE, y, f, cnt,failed);

  x := 4.75;
  y := bernpoly(120,x);
  f := 0.2450175593271593322e71;
  testrel(41, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cosint;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','cosint');

  {special case}
  y := cosint(10000,0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  f := -1234.5;
  y := cosint(0,f);
  testrel(2, NE, y, f, cnt,failed);

  y := cosint(1,-1234.5);
  f := -0.1453956505229364208;
  testrel(3, NE, y, f, cnt,failed);

  {------------------------------}
  y := cosint(10,4);
  f := 1.158627877632916986;
  testrel(4, NE, y, f, cnt,failed);

  y := cosint(10,5);
  f := 1.159689565748002596;
  testrel(5, NE, y, f, cnt,failed);

  y := cosint(10,Pi);
  f := 0.7731263170943631798;
  testrel(6, NE, y, f, cnt,failed);

  y := cosint(10,Pi_2);
  f := 0.3865631585471815899;
  testrel(7, NE, y, f, cnt,failed);

  y := cosint(10,19*Pi_2);
  f := 7.344700012396450208;
  testrel(8, NE, y, f, cnt,failed);

  y := cosint(10,-19*Pi_2);
  f := -7.344700012396450208;
  testrel(9, NE, y, f, cnt,failed);

  y := cosint(10,100);
  f := 24.38341351832059336;
  testrel(10, NE, y, f, cnt,failed);

  y := cosint(5,1);
  f := 0.5286328129112155881;
  testrel(11, NE, y, f, cnt,failed);

  y := cosint(17,Pi_2);
  f := 32768/109395;
  testrel(12, NE, y, f, cnt,failed);

  y := cosint(20,1.25);
  f := 0.2767696820752703606;
  testrel(13, NE, y, f, cnt,failed);

  y := cosint(20,10);
  f := 1.935371243364068614;
  testrel(14, NE, y, f, cnt,failed);

  y := cosint(99,TwoPi);
  f := 0;
  testabs(15, 2, y, f, cnt,failed);

  y := cosint(99,4);
  f := -0.1256451290185490101;
  testrel(16, NE, y, f, cnt,failed);

  y := cosint(9,30);
  f := -0.4063492055789358166;
  testrel(17, NE, y, f, cnt,failed);

  y := cosint(200,10000);
  f := 563.5558003428517485;
  testrel(18, NE, y, f, cnt,failed);

  y := cosint(600,100000);
  f := 3255.962631924894514;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_sinint;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
  NE2 = 4;
  NE3 = 12;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sinint');

  {special case}
  y := sinint(10000,0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  f := -1234.5;
  y := sinint(0,f);
  testrel(2, NE, y, f, cnt,failed);

  y := sinint(1,-1234.5);
  f := 1.989373592132422007;
  testrel(3, NE, y, f, cnt,failed);

  {------------------------------}
  y := sinint(6,0.5);
  f := 0.9185547761246106100e-3;
  testrel(4, NE, y, f, cnt,failed);

  y := sinint(6,1.5);
  f := 0.4204309461889264874;
  testrel(5, NE, y, f, cnt,failed);

  y := sinint(6,2.5);
  f := 0.9771099714670822183;
  testrel(6, NE, y, f, cnt,failed);

  y := sinint(6,3.0);
  f := 0.9817475437693720085;
  testrel(7, NE, y, f, cnt,failed);

  y := sinint(6,5.0);
  f := 1.737945254534703918;
  testrel(18, NE, y, f, cnt,failed);

  y := sinint(6,8.0);
  f := 2.597326791887639688;
  testrel(9, NE, y, f, cnt,failed);

  y := sinint(6,11.0);
  f := 3.440542590614796164;
  testrel(10, NE, y, f, cnt,failed);

  y := sinint(6,-4.0);
  f := -1.009340246947754459;
  testrel(11, NE, y, f, cnt,failed);

  y := sinint(6,-20.0);
  f := -6.025751555775555279;
  testrel(12, NE, y, f, cnt,failed);

  y := sinint(5,0.5);
  f := 0.2226985853239443664e-2;
  testrel(13, NE, y, f, cnt,failed);

  y := sinint(5,1.5);
  f := 0.4628317450440416392;
  testrel(14, NE, y, f, cnt,failed);

  y := sinint(5,2.5);
  f := 1.057683460168249835;
  testrel(15, NE, y, f, cnt,failed);

  y := sinint(5,3.5);
  f := 1.066340666683688580;
  testrel(16, NE, y, f, cnt,failed);

  y := sinint(5,4.5);
  f := 0.7379679184068533471;
  testrel(17, NE, y, f, cnt,failed);

  y := sinint(5,5.5);
  f := 0.2618445135450898748e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := sinint(5,6.5);
  f := 0.16811939182727524784e-4;
  testrel(19, NE, y, f, cnt,failed);

  y := sinint(10,Pi_2);
  f := 0.3865631585471815899;
  testrel(20, NE, y, f, cnt,failed);

  y := sinint(10,Pi);
  f := 0.7731263170943631798;
  testrel(21, NE, y, f, cnt,failed);

  y := sinint(10,-13*Pi_2);
  f := -5.025321061113360669;
  testrel(22, NE, y, f, cnt,failed);

  y := sinint(9,Pi_2);
  f := 0.4063492063492063492;
  testrel(23, NE, y, f, cnt,failed);

  y := sinint(9,Pi);
  f := 0.8126984126984126984;
  testrel(24, NE, y, f, cnt,failed);

  y := sinint(9,-TwoPi);
  f := 0;
  testrel(25, NE, y, f, cnt,failed);

  y := sinint(9,-13*Pi_2);
  f := 0.4063492063492063492;
  testrel(26, NE3, y, f, cnt,failed);

  y := sinint(99,0.5);
  f := 0.1341426041012494184e-33;
  testrel(27, NE, y, f, cnt,failed);

  y := sinint(99,2);
  f := 0.2512885482248272477;
  testrel(28, NE, y, f, cnt,failed);

  y := sinint(99,Pi_2);
  f := 0.1256451290185490101;
  testrel(29, NE2, y, f, cnt,failed);

  y := sinint(100,0.5);
  f := 0.6367642770299571293e-34;
  testrel(30, NE, y, f, cnt,failed);

  y := sinint(100,2);
  f := 0.2500354235665634526;
  testrel(31, NE, y, f, cnt,failed);

  y := sinint(100,Pi_2);
  f := 0.1250184817401874538;
  testrel(32, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_fibpoly;
var
  x,y,f,r: double;
  cnt, failed, k,n: integer;
const
  NE  = 2;
  NE1 = 75;

{$ifdef BIT16}
const
  tol=5;
{$else}
const
  tol=2;
{$endif}

{Values directly from Maple, note that the 16-bit compilers are inaccurate}
{for some literals of length > 19, therefore the tolerance is increased !}
const
  fp15: array[-10..10] of double = (
          -450359962737049.599999999999998,
          14073748835532.8000000000000113,
          -439804651110.399999999999636200,
          13743895347.2000000000116415321,
          -429496729.599999999627470970153,
          13421772.8000000119209289550781,
          -419430.0 - 0.3999996185302734375,  {16 Bit!}
          13107.20001220703125,
          -409.599609375,
          12.8125,
          0,
          12.8125,
          409.599609375,
          13107.20001220703125,
          419430.0 + 0.3999996185302734375,   {16 Bit!}
          13421772.8000000119209289550781,
          429496729.599999999627470970153,
          13743895347.2000000000116415321,
          439804651110.399999999999636200,
          14073748835532.8000000000000113,
          450359962737049.599999999999998);

  fm25: array[-10..10] of double = (
          17491671035336584642638.2191390,
          92899200757966079853.9563808855,
          493392625783670135.494300316212,
          2620434634437155.42014947329880,
          13917268549466.5442822966724634,
          73915357907.6294983029365539548,
          392568420.6778049468994140625,
          2084951.0 + 0.88653564453125,     {16 Bit!}
          11073.291015625,
          58.8125,
          0,
          58.8125,
          -11073.291015625,
          2084951.0 + 0.88653564453125,     {16 Bit!}
          -392568420.6778049468994140625,
          73915357907.6294983029365539548,
          -13917268549466.5442822966724634,
          2620434634437155.42014947329880,
          -493392625783670135.494300316212,
          92899200757966079853.9563808855,
          -17491671035336584642638.2191390);

begin
  {pari: fib(n,x) = ([x,1;1,0]^(n-1))[1,1]}
  cnt := 0;
  failed := 0;
  writeln('Function: ','fibpoly');

  y := fibpoly(10,1.5);
  f := 409.599609375;
  testrel(1, NE, y, f, cnt,failed);

  y := fibpoly(15,-Pi);
  f := 0.2909849191767565688e8; {Fib(15, asd(Pi))}
  testrel(2, NE, y, f, cnt,failed);

  y := fibpoly(123,1.5);
  f := 0.4253529586511730793e37;
  testrel(3, NE, y, f, cnt,failed);

  y := fibpoly(123,-1.5);
  f := 0.4253529586511730793e37;
  testrel(4, NE, y, f, cnt,failed);

  y := fibpoly(-123,1.5);
  f := 0.4253529586511730793e37;
  testrel(5, NE, y, f, cnt,failed);

  y := fibpoly(-123,-1.5);
  f := 0.4253529586511730793e37;
  testrel(6, NE, y, f, cnt,failed);

  y := fibpoly(234,1.5);
  f := 0.1104279415486490206e71;
  testrel(7, NE, y, f, cnt,failed);

  y := fibpoly(234,-1.5);
  f := -0.1104279415486490206e71;
  testrel(8, NE, y, f, cnt,failed);

  y := fibpoly(-234,1.5);
  f := -0.1104279415486490206e71;
  testrel(9, NE, y, f, cnt,failed);

  y := fibpoly(-234,-1.5);
  f := 0.1104279415486490206e71;
  testrel(10, NE, y, f, cnt,failed);

  {Max n for F(n) = fibpoly(n,1) is 23601, 1476 for double}
  y := fibpoly(1476, 1);
  f := 0.1306989223763399318e309;
  testrel(11, NE, y, f, cnt,failed);

  y := fibpoly(32000, 1/32);
  f := 6.875799035044984665e216;
  testrel(12, NE1, y, f, cnt,failed);

  x := 1.5;
  for k:=-10 to 10 do begin
    inc(cnt);
    n := 5*k;
    y := fibpoly(n,x);
    f := fp15[k];
    if f=0.0 then r := y-f
    else r := 1.0-y/f;
    if abs(r) > tol*eps_d then begin
      inc(failed);
      writeln('Test failed: n,x,err= ',n:4, x:8:2, r:30);
    end;
  end;

  x := -2.5;
  for k:=-10 to 10 do begin
    inc(cnt);
    n := 5*k;
    y := fibpoly(n,x);
    f := fm25[k];
    if f=0.0 then r := y-f
    else r := 1.0-y/f;
    if abs(r) > tol*eps_d then begin
      inc(failed);
      writeln('Test failed: n,x,err= ',n:4, x:8:2, r:30);
    end;
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_lucpoly;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  {Pari:  luc(n,x) = trace([x,1;1,0]^n)    }
  {Maple: luc := (n,x) -> fibonacci(n+1,x) + fibonacci(n-1,x); }
  cnt := 0;
  failed := 0;
  writeln('Function: ','lucpoly');

  y := lucpoly(0,0);
  f := 2;
  testrel(1, NE, y, f, cnt,failed);

  y := lucpoly(1,0);
  f := 0;
  testabs(2, 0, y, f, cnt,failed);

  y := lucpoly(2,0);
  f := 2;
  testrel(3, NE, y, f, cnt,failed);

  y := lucpoly(3,0);
  f := 0;
  testabs(4, 0, y, f, cnt,failed);

  y := lucpoly(-1,0);
  f := 0;
  testabs(5, 0, y, f, cnt,failed);

  y := lucpoly(-2,0);
  f := 2;
  testrel(6, NE, y, f, cnt,failed);

  y := lucpoly(-3,0);
  f := 0;
  testabs(7, 0, y, f, cnt,failed);

  y := lucpoly(1,1.5);
  f := 1.5;
  testrel(8, NE, y, f, cnt,failed);

  y := lucpoly(2,1.5);
  f := 4.25;
  testrel(9, NE, y, f, cnt,failed);

  y := lucpoly(9,1.5);
  f := 511.998046875;
  testrel(10, NE, y, f, cnt,failed);

  y := lucpoly(9,-1.5);
  f := -511.998046875;
  testrel(11, NE, y, f, cnt,failed);

  y := lucpoly(-9,1.5);
  f := -511.998046875;
  testrel(12, NE, y, f, cnt,failed);

  y := lucpoly(-9,-1.5);
  f := 511.998046875;
  testrel(13, NE, y, f, cnt,failed);

  y := lucpoly(10,1.5);
  f := 1024.0009765625;
  testrel(14, NE, y, f, cnt,failed);

  y := lucpoly(10,-1.5);
  f := 1024.0009765625;
  testrel(15, NE, y, f, cnt,failed);

  y := lucpoly(-10,1.5);
  f := 1024.0009765625;
  testrel(16, NE, y, f, cnt,failed);

  y := lucpoly(-10,-1.5);
  f := 1024.0009765625;
  testrel(17, NE, y, f, cnt,failed);

  y := lucpoly(125,1.5);
  f := 0.4253529586511730793e38;
  testrel(18, NE, y, f, cnt,failed);

  y := lucpoly(125,-1.5);
  f := -0.4253529586511730793e38;
  testrel(19, NE, y, f, cnt,failed);

  y := lucpoly(-125,1.5);
  f := -0.4253529586511730793e38;
  testrel(20, NE, y, f, cnt,failed);

  y := lucpoly(-125,-1.5);
  f := 0.4253529586511730793e38;
  testrel(21, NE, y, f, cnt,failed);

  y := lucpoly(234,1.5);
  f := 0.2760698538716225515e71;
  testrel(22, NE, y, f, cnt,failed);

  y := lucpoly(234,-1.5);
  f := 0.2760698538716225515e71;
  testrel(23, NE, y, f, cnt,failed);

  y := lucpoly(-234,1.5);
  f := 0.2760698538716225515e71;
  testrel(24, NE, y, f, cnt,failed);

  y := lucpoly(-234,-1.5);
  f := 0.2760698538716225515e71;
  testrel(25, NE, y, f, cnt,failed);

  y := lucpoly(-123,4.5);
  f := -6.407383962038300961e82; {Wolfram Alpha}
  testrel(26, NE, y, f, cnt,failed);

  {Max n for L(n) = lucpoly(n,1) is 23599, 1474 for double}
  y := lucpoly(1474,1);
  f := 1.116302065883468286e308;
  testrel(27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_catalan;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','catalan');

  y := catalan(0);
  f := 1;
  testrel(1, NE, y, f, cnt,failed);

  y := catalan(1e-9);
  f := 1.999999998000000005*0.5;
  testrel(2, NE, y, f, cnt,failed);

  y := catalan(0.125);
  f := 0.4542281453407636314 * 2;
  testrel(3, NE, y, f, cnt,failed);

  y := catalan(0.5);
  f := 0.8488263631567751241;
  testrel(4, NE, y, f, cnt,failed);

  y := catalan(2.5);
  f := 3.104279270973349025;
  testrel(5, NE, y, f, cnt,failed);

  y := catalan(21);
  f := 24466267020.0;
  testrel(6, NE, y, f, cnt,failed);

  y := catalan(35);
{$ifdef BIT16}
  f := 0.38953568686341265775e18 * 8;
{$else}
  f := 3116285494907301262.0;
{$endif}
  testrel(7, NE, y, f, cnt,failed);

  y := catalan(100);
  f := 0.8965199470901314967e57;
  testrel(8, NE, y, f, cnt,failed);

  y := catalan(500);
  f := 0.5394974869170390609e297;
  testrel(9, NE, y, f, cnt,failed);

  y := catalan(-1);
  f := -0.5;
  testrel(10, NE, y, f, cnt,failed);

  y := catalan(-1.25);
  f := -0.3934468663386987420;
  testrel(11, NE, y, f, cnt,failed);

  y := catalan(-12.375);
  f := -0.1218624678667657878e-8;
  testrel(12, NE, y, f, cnt,failed);

  y := catalan(-45.625);
  f := 0.1538947487306520572e-29;
  testrel(13, NE, y, f, cnt,failed);

  y := catalan(-123.875);
  f := 0.4497289298048582343e-78;
  testrel(14, NE, y, f, cnt,failed);

  y := catalan(-456.75);
  f := 0.5916664229531300294e-279;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_kepler;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','kepler');

  {elliptic, Maple: kepler := (m,e) -> fsolve(y-e*sin(y) = m, y=Pi);}
  y := kepler(1, 0.5);
  f := 1.498701133517848314;
  testrel(1, NE, y, f, cnt,failed);

  y := kepler(-1000, 0.5);
  f := -1000.497514775673146;
  testrel(2, NE, y, f, cnt,failed);

  y := kepler(Pi_2, 0.75);
  f := 2.184106679498448926;
  testrel(3, NE, y, f, cnt,failed);

  y := kepler(0.5, 0.125);
  f := 0.5671542538034771510;
  testrel(4, NE, y, f, cnt,failed);

  y := kepler(0.125, 0.9375);
  f := 0.7928034322756140260;
  testrel(5, NE, y, f, cnt,failed);

  y := kepler(1.5, 0.9375);
  f := 2.237023829054169401;
  testrel(6, NE, y, f, cnt,failed);

  y := kepler(10.125, 0.9375);
  f := 9.790088287071411171;
  testrel(7, NE, y, f, cnt,failed);

  y := kepler(0.125, 0.99609375);
  f := 0.9136395753855342618;
  testrel(8, NE, y, f, cnt,failed);

  {in difficult region m near 0, e near 1}
  y := 1-ldexpd(1,-15);
  y := kepler(1/512, y);
  f := 0.2270684928814436695;
  testrel(9, NE, y, f, cnt,failed);

  {arguments from literature, somewhat inexact with binary}
  y := kepler(0.2, 0.99);
  f := 1.066997365281563186;
  testrel(10, NE, y, f, cnt,failed);

  y := kepler(0.06, 0.6);
  f := 0.1491710835982268287;
  testrel(11, NE, y, f, cnt,failed);

  y := kepler(1, 0.9);
  f := 1.862086686874532255;
  testrel(12, NE, y, f, cnt,failed);

  {hyperbolic, Maple: kepler_hyp := (m,e) -> fsolve(e*sinh(x) - x - m, x = signum(m)*ln(2*abs(m)/e + 1.8));}
  y := kepler(-1000, 10);
  f := -5.303631719539061703;
  testrel(13, NE, y, f, cnt,failed);

  y := kepler(1, 2);
  f := 0.8140967963021331692;
  testrel(14, NE, y, f, cnt,failed);

  y := kepler(6,2);
  f := 2.107689797681256377;
  testrel(15, NE, y, f, cnt,failed);

  y := kepler(0.5,1.5);
  f := 0.7673431749540970103;
  testrel(16, NE, y, f, cnt,failed);

  y := kepler(Maxdouble,1.5);
  f := 710.07039496583577777;
  testrel(17, NE, y, f, cnt,failed);

  y := kepler(10,6);
  f := 1+0.3978298998186000144;
  testrel(18, NE, y, f, cnt,failed);

  y := kepler(10000,20);
  f := 6.9084468837654158448;
  testrel(19, NE, y, f, cnt,failed);

  y := kepler(1,20);
  f := 0.5260603476886937829e-1;
  testrel(20, NE, y, f, cnt,failed);

  y := kepler(0,2);
  f := 0;
  testrel(21, NE, y, f, cnt,failed);

  y := kepler(1e-6,1.5);
  f := 0.1999999999996000000e-5;
  testrel(22, NE, y, f, cnt,failed);

  {parabolic, Maple: kepler_para := m -> fsolve(x + x^3/3 = m, x = m^(1/3));}
  y := kepler(2, 1);
  f := 1.287909750704127236;
  testrel(23, NE, y, f, cnt,failed);

  y := kepler(1, 1);
  f := 0.8177316738868235061;
  testrel(24, NE, y, f, cnt,failed);

  y := kepler(0.5,1.0);
  f := 0.4662205239107734274;
  testrel(25, NE, y, f, cnt,failed);

  y := kepler(0.125, 1);
  f := 0.124358923863462980055;
  testrel(26, NE, y, f, cnt,failed);

  y := kepler(1/1024, 1);
  f := 0.9765621895594378539e-3;
  testrel(27, NE, y, f, cnt,failed);

  y := kepler(-1000, 1);
  f := -14.35316011237345298;
  testrel(28, NE, y, f, cnt,failed);

  y := kepler(sqrt_epsh,1.0);
  f := 0.2328306436538696289e-9;
  testrel(29, NE, y, f, cnt,failed);

  y := kepler(1.25e25,1.0);
  f := 0.33471647504108475795e9;
  testrel(30, NE, y, f, cnt,failed);

  y := kepler(Maxdouble, 1.0);
  f := 0.8139772587397598764e103;
  testrel(31, NE, y, f, cnt,failed);

  y := kepler(0.25*Maxdouble, 1.0);
  f := 0.5127735412109745435e103;
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



end.
