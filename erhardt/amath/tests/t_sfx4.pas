{Part 4 of regression test for SPECFUNX unit  (c) 2010  W.Ehrhardt}

unit t_sfx4;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_erfx;
procedure test_erfcx;
procedure test_erfcex;
procedure test_erfgx;
procedure test_inerfcx;
procedure test_erfix;
procedure test_dawsonx;
procedure test_dawson2x;
procedure test_erfinvx;
procedure test_erfcinvx;
procedure test_erf_zx;
procedure test_erf_px;
procedure test_erf_qx;
procedure test_fresnelx;
procedure test_gsix;
procedure test_expint3x;


implementation

uses
  amath, specfunx, t_sfx0;


{---------------------------------------------------------------------------}
procedure test_erfx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfx');
  {-Inf .. -7 .. -1 .. 1 .. 7 .. Inf}

  y := erfx(NegInf_x);
  f := -1;
  testabs( 1, 0, y, f, cnt,failed);

  y := erfx(-100);
  f := -1;
  testrel( 2, NE, y, f, cnt,failed);

  y := erfx(-6.0);
  f := -0.99999999999999997848;
  testrel( 3, NE, y, f, cnt,failed);

  y := erfx(-2.0);
  f := -0.99532226501895273416;
  testrel( 4, NE, y, f, cnt,failed);

  y := erfx(-1.1);
  f := -0.88020506957408169977;
  testrel( 5, NE, y, f, cnt,failed);

  y := erfx(-1.0);
  f := -0.84270079294971486934;
  testrel( 6, NE, y, f, cnt,failed);

  y := erfx(-0.9);
  f := -0.79690821242283212852;
  testrel( 7, NE, y, f, cnt,failed);

  y := erfx(-0.5);
  f := -0.52049987781304653768;
  testrel( 8, NE, y, f, cnt,failed);

  y := erfx(-0.125);
  f := -0.14031620480133381739;
  testrel( 9, NE, y, f, cnt,failed);

  y := erfx(-0.0078125);
  f := -0.88152828951791887128e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := erfx(0.0);
  f := 0;
  testrel(11, NE, y, f, cnt,failed);

  y := erfx(0.0078125);
  f := 0.88152828951791887128e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := erfx(0.125);
  f := 0.14031620480133381739;
  testrel(13, NE, y, f, cnt,failed);

  y := erfx(0.5);
  f := 0.52049987781304653768;
  testrel(14, NE, y, f, cnt,failed);

  y := erfx(0.9);
  f := 0.79690821242283212852;
  testrel(15, NE, y, f, cnt,failed);

  y := erfx(1.0);
  f := 0.84270079294971486934;
  testrel(16, NE, y, f, cnt,failed);

  y := erfx(1.1);
  f := 0.88020506957408169977;
  testrel(17, NE, y, f, cnt,failed);

  y := erfx(2.0);
  f := 0.99532226501895273416;
  testrel(18, NE, y, f, cnt,failed);

  y := erfx(6.0);
  f := 0.99999999999999997848;
  testrel(19, NE, y, f, cnt,failed);

  y := erfx(100);
  f := 1;
  testrel(20, NE, y, f, cnt,failed);

  y := erfx(PosInf_x);
  f := 1;
  testabs(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfcx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfcx');
  {-Inf .. -7 .. -1 .. 1 .. 7 .. 107 .. Inf}

  y := erfcx(NegInf_x);
  f := 2;
  testabs( 1, 0, y, f, cnt,failed);

  y := erfcx(-100);
  f := 2;
  testrel( 2, NE, y, f, cnt,failed);

  y := erfcx(-6.0);
  f := 1.9999999999999999785;
  testrel( 3, NE, y, f, cnt,failed);

  y := erfcx(-2.0);
  f := 1.9953222650189527342;
  testrel( 4, NE, y, f, cnt,failed);

  y := erfcx(-1.1);
  f := 1.8802050695740816998;
  testrel( 5, NE, y, f, cnt,failed);

  y := erfcx(-1.0);
  f := 1.8427007929497148693;
  testrel( 6, NE, y, f, cnt,failed);

  y := erfcx(-0.9);
  f := 1.7969082124228321285;
  testrel( 7, NE, y, f, cnt,failed);

  y := erfcx(-0.5);
  f := 1.5204998778130465377;
  testrel( 8, NE, y, f, cnt,failed);

  y := erfcx(-0.125);
  f := 1.1403162048013338174;
  testrel( 9, NE, y, f, cnt,failed);

  y := erfcx(-0.0078125);
  f := 1.0088152828951791887;
  testrel(10, NE, y, f, cnt,failed);

  y := erfcx(0.0);
  f := 1;
  testrel(11, NE, y, f, cnt,failed);

  y := erfcx(0.0078125);
  f := 0.99118471710482081129;
  testrel(12, NE, y, f, cnt,failed);

  y := erfcx(0.125);
  f := 0.85968379519866618261;
  testrel(13, NE, y, f, cnt,failed);

  y := erfcx(0.5);
  f := 0.47950012218695346232;
  testrel(14, NE, y, f, cnt,failed);

  y := erfcx(0.9);
  f := 0.20309178757716787148;
  testrel(15, NE, y, f, cnt,failed);

  y := erfcx(1.0);
  f := 0.15729920705028513066;
  testrel(16, NE, y, f, cnt,failed);

  y := erfcx(1.1);
  f := 0.11979493042591830023;
  testrel(17, NE, y, f, cnt,failed);

  y := erfcx(2.0);
  f := 0.46777349810472658379e-2;
  testrel(18, NE, y, f, cnt,failed);

  y := erfcx(4.0);
  f := 0.15417257900280018852e-7;
  testrel(19, NE, y, f, cnt,failed);

  y := erfcx(6.0);
  f := 0.21519736712498913117e-16;
  testrel(20, NE, y, f, cnt,failed);

  y := erfcx(10.0);
  f := 0.20884875837625447570e-44;
  testrel(21, NE, y, f, cnt,failed);

  y := erfcx(50.0);
  f := 0.20709207788416560484e-1087;
  testrel(22, NE, y, f, cnt,failed);

  y := erfcx(100);
  f := 0.64059614249217320390e-4345;
  testrel(23, NE, y, f, cnt,failed);

  y := erfcx(PosInf_x);
  f := 0;
  testabs(24, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfix;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfix');

  {erfi := x -> erf(I*x)/I;}
  x := 1e-20;
  y := erfix(x);
  f := 0.1128379167095512574e-19;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1e-5;
  y := erfix(x);
  f := 0.1128379167133125213e-4;
  testrel( 2, NE, y, f, cnt,failed);

  x := -0.125;
  y := erfix(x);
  f := -0.14178547413026538933;
  testrel( 3, NE, y, f, cnt,failed);

  x := 0.5;
  y := erfix(x);
  f := 0.6149520946965109808;
  testrel( 4, NE, y, f, cnt,failed);

  x := 0.75;
  y := erfix(x);
  f := 1.035757284411962968;
  testrel( 5, NE, y, f, cnt,failed);

  x := -1;
  y := erfix(x);
  f := -1.650425758797542876;
  testrel( 6, NE, y, f, cnt,failed);

  x := 2;
  y := erfix(x);
  f := 18.56480241457555260;
  testrel( 7, NE, y, f, cnt,failed);

  x := 4;
  y := erfix(x);
  f := 1296959.730717639232;
  testrel( 8, NE, y, f, cnt,failed);

  x := -10;
  y := erfix(x);
  f := -0.1524307422708669699e43;
  testrel( 9, NE, y, f, cnt,failed);

  x := 26.5;
  y := erfix(x);
  f := 2.0501652832248793153e303;
  testrel(10, NE, y, f, cnt,failed);

  x := 99.762451171875;
  y := erfix(x);
  f := 0.1226243604626303240235e4321;
  testrel(11, NE, y, f, cnt,failed);

  x := 100;
  y := erfix(x);
  f := 4.9689635801475924641e4340;
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfcex;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfcex');

  x := 1e-6;
  f := 0.999998871621832904;
  y := erfcex(x);
  testrel( 1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.8732218450821508096;
  y := erfcex(x);
  testrel( 2, NE, y, f, cnt,failed);

  x := 1/3.0;
  f := 0.7122528886000577139;
  y := erfcex(x);
  testrel( 3, NE, y, f, cnt,failed);

  x := 1;
  f := 0.4275835761558070044;
  y := erfcex(x);
  testrel( 4, NE, y, f, cnt,failed);

  x := 2;
  f := 0.2553956763105057439;
  y := erfcex(x);
  testrel( 5, NE, y, f, cnt,failed);

  x := 10.0;
  f := 0.5614099274382258586e-1;
  y := erfcex(x);
  testrel( 6, NE, y, f, cnt,failed);

  x := 100.0;
  f := 0.5641613782989432904e-2;
  y := erfcex(x);
  testrel( 7, NE, y, f, cnt,failed);

  x := 130.0;
  f := 0.4339791484842787238e-2;
  y := erfcex(x);
  testrel( 8, NE, y, f, cnt,failed);

  x := 5000.0;
  f := 0.1128379144527930586e-3;
  y := erfcex(x);
  testrel( 9, NE, y, f, cnt,failed);

  x := 70000.0;
  f := 0.8059851192716941733e-5;
  y := erfcex(x);
  testrel(10, NE, y, f, cnt,failed);

  x := -1e-15;
  f := 1.000000000000001128;
  y := erfcex(x);
  testrel(11, NE, y, f, cnt,failed);

  x := -0.25;
  f := 1.358642370104722115;
  y := erfcex(x);
  testrel(12, NE, y, f, cnt,failed);

  x := -1;
  f := 5.008980080762283466;
  y := erfcex(x);
  testrel(13, NE, y, f, cnt,failed);

  x := -8/7;
  f := 6.992173861723201774;
  y := erfcex(x);
  testrel(14, NE, y, f, cnt,failed);

  x := -2;
  f := 108.9409043899779724;
  y := erfcex(x);
  testrel(15, NE, y, f, cnt,failed);

  x := -6.75;
  f := 0.1226231102139051807e21;
  y := erfcex(x);
  testrel(16, NE, y, f, cnt,failed);

  x := -26.5;
  f := 0.1924553162418568809e306;
  y := erfcex(x);
  testrel(17, NE, y, f, cnt,failed);

  x := -106.5637380121098418;
  f := PosInf_x;
  y := erfcex(x);
  testabs(18, 0, y, f, cnt,failed);

  x := -106.5625;
  f := 0.9138167535613739207e4932;
  y := erfcex(x);
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfgx;
var
  x,y,f,p: extended;
  cnt,failed,i: integer;
const
  NE  = 4;
  NE2 = 2*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfgx');

  {Maple: gerf := (p,x) -> int(exp(-t^p), t=0..x);}
  {       gerf := (p,x) -> (GAMMA(1/p)-GAMMA(1/p,x^p))/p;}
  {Pari:  gerf(p,x) = intnum(t=0,x,exp(-t^p))}
  {       gerf(p,x) = incgamc(1/p, x^p)/p}

  p := 0.5;
  y := erfgx(p,1/1024);
  f := 0.9564538925403311883e-3;
  testrel( 1, NE, y, f, cnt,failed);

  y := erfgx(p,0.125);
  f := 0.9910074638765149504e-1;
  testrel( 2, NE, y, f, cnt,failed);

  y := erfgx(p,0.5);
  f := 0.3165581866568181266;
  testrel( 3, NE, y, f, cnt,failed);

  y := erfgx(p,1);
  f := 0.5284822353142307136;
  testrel( 4, NE, y, f, cnt,failed);

  y := erfgx(p,2);
  f := 0.8261285649781240111;
  testrel( 5, NE, y, f, cnt,failed);

  y := erfgx(p,10);
  f := 1.647628069579945710;
  testrel( 6, NE, y, f, cnt,failed);

  y := erfgx(p,1000);
  f := 1.999999999998795093;
  testrel( 7, NE, y, f, cnt,failed);


  p := 5;
  y := erfgx(p,1/1024);
  f := 2*0.4882812499999999277e-3;
  testrel( 8, NE, y, f, cnt,failed);

  y := erfgx(p,0.125);
  f := 0.12499936422241396436;
  testrel( 9, NE, y, f, cnt,failed);

  y := erfgx(p,0.5);
  f := 0.4974178699312369020;
  testrel(10, NE, y, f, cnt,failed);

  y := erfgx(p,1);
  f := 0.8700746676858926573;
  testrel(11, NE, y, f, cnt,failed);

  y := erfgx(p,2);
  f := 0.9181687423997604561;
  testrel(12, NE, y, f, cnt,failed);

  y := erfgx(p,10);
  f := 0.9181687423997606106;
  testrel(13, NE, y, f, cnt,failed);


  p := 1.5;
  y := erfgx(p,1/1024);
  f := 2*0.4882752895923654593e-3;
  testrel(14, NE, y, f, cnt,failed);

  y := erfgx(p,0.125);
  f := 0.12282048474733161436;
  testrel(15, NE, y, f, cnt,failed);

  y := erfgx(p,0.5);
  f := 0.4364761380885268077;
  testrel(16, NE, y, f, cnt,failed);

  y := erfgx(p,1);
  f := 0.6997923277614944777;
  testrel(17, NE, y, f, cnt,failed);

  y := erfgx(p,2);
  f := 0.8772524660847196675;
  testrel(18, NE, y, f, cnt,failed);

  y := erfgx(p,10);
  f := 0.9027452929509297575;
  testrel(19, NE, y, f, cnt,failed);

  y := erfgx(p,1000);
  f := 0.9027452929509336113;
  testrel(20, NE, y, f, cnt,failed);

  {pari}
  p := 0.01;
  y := erfgx(p,1/1024);
  f := 0.0003877209165992133333;
  testrel(21, NE, y, f, cnt,failed);

  y := erfgx(p,0.125);
  f := 0.04740070277221238905;
  testrel(22, NE, y, f, cnt,failed);

  y := erfgx(p,0.5);
  f := 0.1870537286914434240;
  testrel(23, NE, y, f, cnt,failed);

  y := erfgx(p,1);
  f := 0.3715578714528098103;
  testrel(24, NE, y, f, cnt,failed);

  y := erfgx(p,2);
  f := 0.7380162212946553285;
  testrel(25, NE, y, f, cnt,failed);

  y := erfgx(p,10);
  f := 3.630877526710845930;
  testrel(26, NE, y, f, cnt,failed);

  y := erfgx(p,1000);
  f := 346.1598375906238578;
  testrel(27, NE, y, f, cnt,failed);

  y := erfgx(p,1e100);
  f := 5.038299480618184746e95;
  testrel(28, NE, y, f, cnt,failed);

  {selected values}
  y := erfgx(20,0);
  f := 0;
  testrel(29, NE, y, f, cnt,failed);

  y := erfgx(0,20);
  f := 7.357588823428846432;
  testrel(30, NE, y, f, cnt,failed);

  y := erfgx(1,20);
  f := 2*0.4999999989694231888;
  testrel(31, NE, y, f, cnt,failed);

  y := erfgx(0.25,10);
  f := 2.525854435942086756;
  testrel(32, NE, y, f, cnt,failed);

  y := erfgx(10,0.125);
  f := 0.1249999999894167889;
  testrel(33, NE, y, f, cnt,failed);

  y := erfgx(20,0.25);
  f := 0.2499999999999891727;
  testrel(34, NE, y, f, cnt,failed);

  y := erfgx(40,0.5);
  f := 0.4999999999999889086;
  testrel(35, NE, y, f, cnt,failed);

  y := erfgx(100,1-1/1024);
  f := 0.991745400651105407;
  testrel(36, NE, y, f, cnt,failed);

  y := erfgx(10000,1-1/1024);
  f := 2*0.4995117158972298877;
  testrel(37, NE, y, f, cnt,failed);

  {expint3 tests}
  for i:=0 to 63 do begin
    x := i/16;
    y := erfgx(3,x);
    f := expint3x(x);
    testrel(100+i, NE2, y, f, cnt,failed);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_inerfcx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
  NE2 = 16;
  NE3 = 24;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','inerfcx');

  y := inerfcx(3,1e-40);
  f := 0.9403159725795938116e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := inerfcx(250,1e-40);
  f := 0.2935791617973459071e-284;
  testrel(2, NE, y, f, cnt,failed);

  y := inerfcx(3,0.5);
  f := 0.2775132137764753636e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := inerfcx(10,20);
  f := 0.2480339545677166807e-17;
  testrel(4, NE, y, f, cnt,failed);

  y := inerfcx(10,2);
  f := 0.4519702487098470796e-8;
  testrel(5, NE, y, f, cnt,failed);

  y := inerfcx(10,1);
  f := 0.1305302483913995620e-6;
  testrel(6, NE, y, f, cnt,failed);

  y := inerfcx(10,0.5);
  f := 0.9245065129805797129e-6;
  testrel(7, NE, y, f, cnt,failed);

  y := inerfcx(10,0.25);
  f := 0.26649953573332452744e-5;
  testrel(8, NE, y, f, cnt,failed);

  y := inerfcx(10,0.125);
  f := 0.4622606908901843597e-5;
  testrel(9, NE, y, f, cnt,failed);

  y := inerfcx(10,1/1024);
  f := 0.8101666489669277473e-5;
  testrel(10, NE, y, f, cnt,failed);

  y := inerfcx(100,8);
  f := 0.3327723273692355744e-132;
  testrel(11, NE, y, f, cnt,failed);

  y := inerfcx(100,2);
  f := 0.8417913154612625873e-106;
  testrel(12, NE, y, f, cnt,failed);

  y := inerfcx(100,0.5);
  f := 0.2448087104431131554e-97;
  testrel(13, NE, y, f, cnt,failed);

  y := inerfcx(100,0.25);
  f := 0.7728184538695565628e-96;
  testrel(14, NE, y, f, cnt,failed);

  y := inerfcx(100,1/1024);
  f := 0.2558072521643705004e-94;
  testrel(15, NE, y, f, cnt,failed);

  y := inerfcx(100,1e-10);
  f := 0.2593734749448854344e-94;
  testrel(16, NE, y, f, cnt,failed);

  y := inerfcx(5,0.125);
  f := 0.6252579526727030871e-2;
  testrel(17, NE, y, f, cnt,failed);

  y := inerfcx(20,0.125);
  f := 0.1189350337542980588e-12;
  testrel(18, NE, y, f, cnt,failed);

  y := inerfcx(100,0.125);
  f := 0.4442699570076624848e-95;
  testrel(19, NE, y, f, cnt,failed);

  y := inerfcx(250,0.125);
  f := 0.1803034248451498754e-285;
  testrel(20, NE, y, f, cnt,failed);

  y := inerfcx(5,0.75);
  f := 0.9967894387400569614e-3;
  testrel(21, NE, y, f, cnt,failed);

  y := inerfcx(20,0.75);
  f := 0.2815855137500476563e-14;
  testrel(22, NE, y, f, cnt,failed);

  y := inerfcx(50,0.75);
  f := 0.4007073297196590972e-43;
  testrel(23, NE, y, f, cnt,failed);

  y := inerfcx(125,0.75);
  f := 0.8592345382450687318e-129;
  testrel(24, NE, y, f, cnt,failed);

  y := inerfcx(250,0.75);
  f := 0.1984987852267995133e-291;
  testrel(25, NE, y, f, cnt,failed);

  y := inerfcx(5,50);     {AE}
  f := 0.1123656973372330445e-11;
  testrel(26, NE, y, f, cnt,failed);

  y := inerfcx(20,50);    {AE}
  f := 0.1077656273122095687e-41;
  testrel(27, NE, y, f, cnt,failed);

  y := inerfcx(100,50);   {AE}
  f := 0.4110694394268175203e-202;
  testrel(28, NE, y, f, cnt,failed);

  y := inerfcx(150,50);   {CF}
  f := 0.1215088638985536402e-302;
  testrel(29, NE, y, f, cnt,failed);

  y := inerfcx(10,1000);
  f := 0.5509482087061094431e-36;
  testrel(30, NE, y, f, cnt,failed);

  y := inerfcx(90,1000);
  f := 0.4547958756280478593e-300;
  testrel(31, NE, y, f, cnt,failed);

  {x<0}
  y := inerfcx(-1,-2);
  f := 0.2066698535409205386e-1;
  testrel(32, NE, y, f, cnt,failed);

  y := inerfcx(0,-2);
  f := 1.9953222650189527342;
  testrel(33, NE, y, f, cnt,failed);

  y := inerfcx(1,-2);
  f := 4.000978022714951495;
  testrel(34, NE, y, f, cnt,failed);

  y := inerfcx(5,-0.5);
  f := 0.4374437542410532658e-1;
  testrel(35, NE, y, f, cnt,failed);

  y := inerfcx(5,-2);
  f := 1.325001048378169994;
  testrel(36, NE, y, f, cnt,failed);

  y := inerfcx(5,-10);
  f := 1750.625;
  testrel(37, NE, y, f, cnt,failed);

  y := inerfcx(5,-300);
  f := 40502250018.75;
  testrel(38, NE, y, f, cnt,failed);

  y := inerfcx(20,-0.5);
  f := 0.5711213876890956573e-11;
  testrel(39, NE, y, f, cnt,failed);

  y := inerfcx(20,-2);
  f := 0.1557830391267366734e-7;
  testrel(40, NE, y, f, cnt,failed);

  y := inerfcx(20,-10);
  f := 196.8746175638803787;
  testrel(41, NE, y, f, cnt,failed);

  y := inerfcx(20,-300);
  f := 0.2869385160990788617e32;
  testrel(42, NE, y, f, cnt,failed);

  y := inerfcx(1000,-500);
  f := 0.1257148984018811968e133;
  testrel(43, NE, y, f, cnt,failed);

  y := inerfcx(100,-1000);
  f := 0.2148330859419028936e143;
  testrel(44, NE, y, f, cnt,failed);

  y := inerfcx(5000,-2000);
  f := 0.3183024088162901648e181;
  testrel(45, NE3, y, f, cnt,failed);

  y := inerfcx(2000,-1000);
  f := 0.1636906036218412122e266;
  testrel(46, NE3, y, f, cnt,failed);

  y := inerfcx(30000,-11000);
  f := 0.2806686227497745515e-44;
  testrel(47, NE3, y, f, cnt,failed);

  {--------- extended only ---------}

  y := inerfcx(500,8);
  f := 0.6053717215571535792e-740;
  testrel(48, NE, y, f, cnt,failed);

  y := inerfcx(1000,8);
  f := 0.3437870361645832591e-1577;
  testrel(49, NE, y, f, cnt,failed);

  y := inerfcx(2000,8);
  f := 0.7612760696923994492e-3376;
  testrel(50, NE2, y, f, cnt,failed);

  y := inerfcx(2800,8);
  f := 0.9088747209936031381e-4888;
  testrel(51, NE, y, f, cnt,failed);

  y := inerfcx(2000,2);
  f := 0.1792400092894991409e-3223;
  testrel(52, NE, y, f, cnt,failed);

  y := inerfcx(2000,0.5);
  f := 0.4510111198053716206e-3183;
  testrel(53, NE, y, f, cnt,failed);

  y := inerfcx(2000,0.125);
  f := 0.8033895338468148036e-3173;
  testrel(54, NE, y, f, cnt,failed);

  y := inerfcx(2000,0);
  f := 0.2164534188917535707e-3169;
  testrel(55, NE, y, f, cnt,failed);

  y := inerfcx(1400,400);   {CF}
  f := 0.3135483200306032641e-4068;
  testrel(56, NE, y, f, cnt,failed);

  y := inerfcx(1200,400);   {AF}
  f := 0.2920331348555060292e-3487;
  testrel(57, NE2, y, f, cnt,failed);

  {x<0}
  y := inerfcx(1000,-10);
  f := 0.1098327100833423869e-1260;
  testrel(58, NE2, y, f, cnt,failed);

  y := inerfcx(1000,-20);
  f := 0.1345511944678638231e-1120;
  testrel(59, NE2, y, f, cnt,failed);

  y := inerfcx(1000,-50);
  f := 0.2790843160428145087e-831;
  testrel(60, NE, y, f, cnt,failed);

  y := inerfcx(100,-1e20);
  f := 0.2143020576250933846e1843;
  testrel(61, NE, y, f, cnt,failed);

  y := inerfcx(500,-1e10);
  f := 0.1639160426434480284e3867;
  testrel(62, NE2, y, f, cnt,failed);

  y := inerfcx(10000,-1200);
  f := 0.1486872156150231891e-4859;
  testrel(63, NE2, y, f, cnt,failed);

  y := inerfcx(10000,-2600);
  f := 0.1531112204308584074e-1507;
  testrel(64, NE2, y, f, cnt,failed);

  y := inerfcx(10000,-4000);
  f := 0.1333441997807055193e363;
  testrel(65, NE3, y, f, cnt,failed);

  y := inerfcx(10000,-6000);
  f := 0.4579105761886398102e2123;
  testrel(66, NE2, y, f, cnt,failed);

  y := inerfcx(10000,-10000);
  f := 0.9022207499862962510e4341;
  testrel(67, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_dawsonx;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 7;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','dawsonx');

  x := 0;
  y := dawsonx(x);
  f := 0;
  testrel( 1, NE, y, f, cnt,failed);

  x := succx(0); {x=0.5^16445}
  y := dawsonx(x);
  {$ifdef HAS_DENORM_LIT}
    f := 0.36451995318824746025e-4950;
  {$else}
    f := 0.36451995318824746025e-950;
    f := f*1e-4000;
  {$endif}
  testrel( 2, NE, y, f, cnt,failed);

  x := 1e-100;
  y := dawsonx(x);
  f := 1e-100;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1e-6;
  y := dawsonx(x);
  f := 0.99999999999933333333e-6;
  testrel( 4, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := dawsonx(x);
  f := 0.97656187911852043720e-3;
  testrel( 5, NE, y, f, cnt,failed);

  x := 0.125;
  y := dawsonx(x);
  f := 0.12370601848283973394;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.5;
  y := dawsonx(x);
  f := 0.42443638350202229593;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.92413887300459176701;
  y := dawsonx(x);
  f := 0.54104422463518169847;
  testrel( 8, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := dawsonx(x);
  f := 0.53815344009396028924;
  testrel( 9, NE, y, f, cnt,failed);

  x := 1.0;
  y := dawsonx(x);
  f := 0.53807950691276841914;
  testrel(10, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := dawsonx(x);
  f := 0.53800469268824945169;
  testrel(11, NE, y, f, cnt,failed);

  x := 3.9990234375;
  y := dawsonx(x);
  f := 0.12938197933297510231;
  testrel(12, NE, y, f, cnt,failed);

  x := 4.0;
  y := dawsonx(x);
  f := 0.12934800123600511559;
  testrel(13, NE, y, f, cnt,failed);

  x := 4.0009765625;
  y := dawsonx(x);
  f := 0.12931404180823832109;
  testrel(14, NE, y, f, cnt,failed);

  x := 3.0e9;
  y := dawsonx(x);
  f := 0.16666666666666666667e-9;
  testrel(15, NE, y, f, cnt,failed);

  x := 4.0e9;
  y := dawsonx(x);
  f := 0.12500000000000000000e-9;
  testrel(16, NE, y, f, cnt,failed);

  x := 1.0e307;
  y := dawsonx(x);
  f := 5e-308;
  testrel(17, NE, y, f, cnt,failed);

  x := 1.4871643691965395E+4931;
  y := dawsonx(x);
  f := 0.33621031431120939728e-4931;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_dawson2x;
var
  p,x,y,f: extended;
  cnt, failed: integer;
const
  NE = 5;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','dawson2x');

  {Pari: gd(p,x)=intnum(t=0,x,exp(t^p))*exp(-x^p)}

  {Standard Dawson}
  p := 2;
  y := dawson2x(p,1);
  f := 0.5380795069127684191;
  testrel(1, NE, y, f, cnt,failed);

  y := dawson2x(p,3);
  f := 0.1782710306105582873;
  testrel(2, NE, y, f, cnt,failed);

  y := dawson2x(p,5);
  f := 0.1021340744242768355;
  testrel(3, NE, y, f, cnt,failed);

  y := dawson2x(p,100);
  f := 0.5000250037509378283e-2;
  testrel(4, NE, y, f, cnt,failed);

  y := dawson2x(p,1e300);
  f := 0.5e-300;
  testrel(5, NE, y, f, cnt,failed);

  y := dawson2x(p,1/16);
  f := 0.6233749361289893996e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := dawson2x(p,1e-5);
  f := 0.9999999999333333333e-5;
  testrel(7, NE, y, f, cnt,failed);

  y := dawson2x(p,1e-10);
  f := 1e-10;
  testrel(8, NE, y, f, cnt,failed);

  y := dawson2x(p,1e-15);
  f := 1e-15;
  testrel(9, NE, y, f, cnt,failed);

  {Special case x=1}
  x := 1;
  y := dawson2x(1+1/1024,x);
  f := 0.6319767998075155365;
  testrel(10, NE, y, f, cnt,failed);

  y := dawson2x(10,x);
  f := 0.4125034088017899048;
  testrel(11, NE, y, f, cnt,failed);

  y := dawson2x(100,x);
  f := 0.3726859445420274430;
  testrel(12, NE, y, f, cnt,failed);

  y := dawson2x(1000,x);
  f := 0.3683638488980294743;
  testrel(13, NE, y, f, cnt,failed);

  y := dawson2x(1e10,x);
  f := 0.3678794412199252323;
  testrel(14, NE, y, f, cnt,failed);

  y := dawson2x(1e300,x);
  f := 0.3678794411714423216;
  testrel(15, NE, y, f, cnt,failed);

  y := dawson2x(1/1024,x);
  f := 0.9990253402056163500;
  testrel(16, NE, y, f, cnt,failed);

  y := dawson2x(1e-9,x);
  f := 0.999999999000000002;
  testrel(17, NE, y, f, cnt,failed);

  y := dawson2x(1e-10,x);
  f := 0.9999999999;
  testrel(18, NE, y, f, cnt,failed);

  y := dawson2x(1e-11,x);
  f := 0.99999999999;
  testrel(19, NE, y, f, cnt,failed);

  y := dawson2x(1e-16,x);
  f := 0.9999999999999999;
  testrel(20, NE, y, f, cnt,failed);

  {Other cases}
  y := dawson2x(1/80,80);
  f := 78.97000764508719285;
  testrel(21, NE, y, f, cnt,failed);

  y := dawson2x(1.5,100);
  f := 0.6668891858788577765e-1;
  testrel(22, NE, y, f, cnt,failed);

  y := dawson2x(1,4);
  f := 0.981684361111265820;
  testrel(23, NE, y, f, cnt,failed);

  y := dawson2x(1,1e-10);
  f := 0.99999999995e-10;
  testrel(24, NE, y, f, cnt,failed);

  y := dawson2x(1/2,100);
  f := 18.00009079985952497;
  testrel(25, NE, y, f, cnt,failed);

  y := dawson2x(1/2,1e10);
  f := 199998.0;
  testrel(26, NE, y, f, cnt,failed);

  y := dawson2x(0.1,10);
  f := 8.964926863611165808;
  testrel(27, NE, y, f, cnt,failed);

  y := dawson2x(0.01,100);
  f := 98.9737752892581009;
  testrel(28, NE, y, f, cnt,failed);

  y := dawson2x(0.01,1e10);
  f := 9876873773.46623806;
  testrel(29, NE, y, f, cnt,failed);

  y := dawson2x(0.01,1e100);
  f := 9.09837270131716558e99;
  testrel(30, NE, y, f, cnt,failed);

  y := dawson2x(1/1024,1/2);
  f := 0.49951299954466997501;
  testrel(31, NE, y, f, cnt,failed);

  y := dawson2x(1e-5,1e20);
  f := 9.99989995593902475e19;
  testrel(32, NE, y, f, cnt,failed);

  y := dawson2x(1e-10,1e10);
  f := 9999999999.0;
  testrel(33, NE, y, f, cnt,failed);

  y := dawson2x(25,2);
  f := 0.2384185859227731617e-8;
  testrel(34, NE, y, f, cnt,failed);

  y := dawson2x(5,50);
  f := 0.3200000008192000047e-7;
  testrel(35, NE, y, f, cnt,failed);

  y := dawson2x(4,100);
  f := 0.2500000018750000328e-6;
  testrel(36, NE, y, f, cnt,failed);

  y := dawson2x(40,1.25);
  f := 0.4154375964424212495e-5;
  testrel(37, NE, y, f, cnt,failed);

  y := dawson2x(40,1.5);
  f := 0.3391415055475334289e-8;
  testrel(38, NE, y, f, cnt,failed);

  y := dawson2x(300,1+1/16);
  f := 4.4722500002707202404e-11;
  testrel(39, NE, y, f, cnt,failed);

  y := dawson2x(300,1-1/16);
  f := 0.937499996353225120;
  testrel(40, NE, y, f, cnt,failed);

  y := dawson2x(500,3);
  f := 0.1650151773721882811e-240;
  testrel(41, NE, y, f, cnt,failed);

  y := dawson2x(200,10);
  f := 5.0e-202;
  testrel(42, NE, y, f, cnt,failed);

  y := dawson2x(200,8);
  f := 0.9639679460411536471e-182;
  testrel(43, NE, y, f, cnt,failed);

  y := dawson2x(1000,1-1/16);
  f := 0.9375;
  testrel(44, NE, y, f, cnt,failed);

  y := dawson2x(1000,1+1/16);
  f := 0.4981845059124721260e-29;
  testrel(45, NE, y, f, cnt,failed);

  y := dawson2x(10000,1+1/16);
  f := 0.545684622222586877e-267;
  testrel(46, NE2, y, f, cnt,failed);

  y := dawson2x(1000,2);
  f := 1.866527237006437758e-304;
  testrel(47, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfinvx;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erf_invx');

  x := 0;
  y := erf_invx(x);
  f := 0;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1e-4900;
  y := erf_invx(x);
  {$ifdef HAS_DENORM_LIT}
    f := 0.886226925452758013649e-4900;
  {$else}
    f := 0.886226925452758013649e-900;
    f := f*1e-4000;
  {$endif}
  testrel( 2, NE, y, f, cnt,failed);

  x := 1e-800;
  y := erf_invx(x);
  f := 0.886226925452758013649e-800;
  testrel( 3, NE, y, f, cnt,failed);

  x := -1e-100;
  y := erf_invx(x);
  f := -0.886226925452758013649e-100;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1e-6;
  y := erf_invx(x);
  f := 0.8862269254529900273156e-6;
  testrel( 5, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := erf_invx(x);
  f := 0.86545619796713755345e-3;
  testrel( 6, NE, y, f, cnt,failed);

  x := -0.125;
  y := erf_invx(x);
  f := -0.1112354518409499591471840;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.5;
  y := erf_invx(x);
  f := 0.47693627620446987338;
  testrel( 8, NE, y, f, cnt,failed);

  x := -0.75;
  y := erf_invx(x);
  f := -0.8134198475976185417;
  testrel( 9, NE, y, f, cnt,failed);

  x := 0.9375;
  y := erf_invx(x);
  f := 1.31715033498613074888;
  testrel(10, NE, y, f, cnt,failed);

  x := -0.9990234375;
  y := erf_invx(x);
  f := -2.3314677736219476723;
  testrel(11, NE, y, f, cnt,failed);

  x := 1.0-ldexp(1,-20);
  y := erf_invx(x);
  f := 3.465505025803330717256;
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erf_zx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','erf_zx');
  {Maple: statevalf[pdf,normald](x)}

  y := erf_zx(NegInf_x);
  f := 0;
  testabs( 1, 0, y, f, cnt,failed);

  y := erf_zx(-100);
  f := 0.13443134677817230434e-2171;
  testrel( 2, NE, y, f, cnt,failed);

  y := erf_zx(-6.0);
  f := 0.60758828498232854869e-8;
  testrel( 3, NE, y, f, cnt,failed);

  y := erf_zx(-2.0);
  f := 0.53990966513188051949e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := erf_zx(-1.1);
  f := 0.21785217703255053138;
  testrel( 5, NE, y, f, cnt,failed);

  y := erf_zx(-1.0);
  f := 0.24197072451914334980;
  testrel( 6, NE, y, f, cnt,failed);

  y := erf_zx(-0.9);
  f := 0.26608524989875482182;
  testrel( 7, NE, y, f, cnt,failed);

  y := erf_zx(-0.5);
  f := 0.35206532676429947777;
  testrel( 8, NE, y, f, cnt,failed);

  y := erf_zx(-0.125);
  f := 0.39583768694474948414;
  testrel( 9, NE, y, f, cnt,failed);

  y := erf_zx(-0.0078125);
  f := 0.39893010583499324766;
  testrel(10, NE, y, f, cnt,failed);

  y := erf_zx(0.0);
  f := 0.39894228040143267794;
  testrel(11, NE, y, f, cnt,failed);

  y := erf_zx(0.0078125);
  f := 0.39893010583499324766;
  testrel(12, NE, y, f, cnt,failed);

  y := erf_zx(0.125);
  f := 0.39583768694474948414;
  testrel(13, NE, y, f, cnt,failed);

  y := erf_zx(0.5);
  f := 0.35206532676429947777;
  testrel(14, NE, y, f, cnt,failed);

  y := erf_zx(0.9);
  f := 0.26608524989875482182;
  testrel(15, NE, y, f, cnt,failed);

  y := erf_zx(1.0);
  f := 0.24197072451914334980;
  testrel(16, NE, y, f, cnt,failed);

  y := erf_zx(1.1);
  f := 0.21785217703255053138;
  testrel(17, NE, y, f, cnt,failed);

  y := erf_zx(2.0);
  f := 0.53990966513188051949e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := erf_zx(4.0);
  f := 0.5*0.26766045152977070355e-3; {0.133830225764885351774074488441e-3}
  testrel(19, NE, y, f, cnt,failed);

  y := erf_zx(6.0);
  f := 0.60758828498232854869e-8;
  testrel(20, NE, y, f, cnt,failed);

  y := erf_zx(10.0);
  f := 0.76945986267064193463e-22;
  testrel(21, NE, y, f, cnt,failed);

  y := erf_zx(50.0);
  f := 0.54051492041927084473e-543;
  testrel(22, NE, y, f, cnt,failed);

  y := erf_zx(100);
  f := 0.13443134677817230434e-2171;
  testrel(23, NE, y, f, cnt,failed);

  y := erf_zx(150);
  f := 0.61374597109796063650e-4886;
  testrel(24, NE, y, f, cnt,failed);

  y := erf_zx(PosInf_x);
  f := 0;
  testabs(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erf_px;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 10;
  NE1 = 25;
  NE2 = 80;
  NE3 = 260;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','erf_px');
  {Maple: statevalf[cdf,normald](x)}

  y := erf_px(NegInf_x);
  f := 0;
  testabs( 1, 0, y, f, cnt,failed);

  y := erf_px(-35);
  f := 0.1124910706472406243978835e-267;
  testrel( 2, NE, y, f, cnt,failed);

  y := erf_px(-30);
  f := 0.4906713927148187059531905e-197;
  testrel( 3, NE, y, f, cnt,failed);

  y := erf_px(-20);
  f := 0.27536241186062336951e-88;
  testrel( 4, NE, y, f, cnt,failed);

  y := erf_px(-10);
  f := 0.76198530241605260660e-23;
  testrel( 5, NE, y, f, cnt,failed);

  y := erf_px(-6.0);
  f := 0.98658764503769814070e-9;
  testrel( 6, NE, y, f, cnt,failed);

  y := erf_px(-2.0);
  f := 0.227501319481792072003e-1;
  testrel( 7, NE, y, f, cnt,failed);

  y := erf_px(-1.1);
  f := 0.13566606094638267518;
  testrel( 8, NE, y, f, cnt,failed);

  y := erf_px(-1.0);
  f := 0.15865525393145705142;
  testrel( 9, NE, y, f, cnt,failed);

  y := erf_px(-0.9);
  f := 0.18406012534675948856;
  testrel(10, NE, y, f, cnt,failed);

  y := erf_px(-0.5);
  f := 0.30853753872598689637;
  testrel(11, NE, y, f, cnt,failed);

  y := erf_px(-0.125);
  f := 0.45026177516988710702;
  testrel(12, NE, y, f, cnt,failed);

  y := erf_px(-0.0078125);
  f := 0.49688329513915741955;
  testrel(13, NE, y, f, cnt,failed);

  y := erf_px(0.0);
  f := 0.5;
  testrel(14, NE, y, f, cnt,failed);

  y := erf_px(0.0078125);
  f := 0.50311670486084258045;
  testrel(15, NE, y, f, cnt,failed);

  y := erf_px(0.125);
  f := 0.54973822483011289298;
  testrel(16, NE, y, f, cnt,failed);

  y := erf_px(0.5);
  f := 0.69146246127401310364;
  testrel(17, NE, y, f, cnt,failed);

  y := erf_px(0.9);
  f := 0.81593987465324051145;
  testrel(18, NE, y, f, cnt,failed);

  y := erf_px(1.0);
  f := 0.84134474606854294859;
  testrel(19, NE, y, f, cnt,failed);

  y := erf_px(1.1);
  f := 0.86433393905361732483;
  testrel(20, NE, y, f, cnt,failed);

  y := erf_px(2.0);
  f := 0.97724986805182079280;
  testrel(21, NE, y, f, cnt,failed);

  y := erf_px(4.0);
  f := 0.99996832875816688008;
  testrel(22, NE, y, f, cnt,failed);

  y := erf_px(6.0);
  f := 0.99999999901341235496;
  testrel(23, NE, y, f, cnt,failed);

  y := erf_px(10.0);
  f := 1;
  testrel(24, NE, y, f, cnt,failed);

  y := erf_px(PosInf_x);
  f := 1;
  testabs(25, NE, y, f, cnt,failed);

  {precision for -150 < x < -40 is increased but still suboptimal}
  y := erf_px(-150);
  f := 0.40914579809030025090e-4888;
  testrel(26, NE3, y, f, cnt,failed);

  y := erf_px(-107);
  f := 0.2836133751848002929e-2488;
  testrel(27, NE1, y, f, cnt,failed);

  y := erf_px(-80);
  f := 0.90242300162047361580e-1392;
  testrel(28, NE2, y, f, cnt,failed);

  y := erf_px(-40);
  f := 0.36558935409150297108e-349;
  testrel(29, NE1, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erf_qx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','erf_qx');
  {Maple: statevalf[cdf,normald](-x)}
  y := erf_qx(0);
  f := 0.5;
  testrel(1, NE, y, f, cnt,failed);

  y := erf_qx(0.5);
  f := 0.3085375387259868964;
  testrel(2, NE, y, f, cnt,failed);

  y := erf_qx(1);
  f := 0.1586552539314570514;
  testrel(3, NE, y, f, cnt,failed);

  y := erf_qx(2);
  f := 0.2275013194817920720e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := erf_qx(5);
  f := 0.2866515718791939117e-6;
  testrel(5, NE, y, f, cnt,failed);

  y := erf_qx(10);
  f := 0.7619853024160526066e-23;
  testrel(6, NE, y, f, cnt,failed);

  y := erf_qx(20);
  f := 0.27536241186062336951e-88;
  testrel(7, NE, y, f, cnt,failed);

  y := erf_qx(35);
  f := 0.1124910706472406244e-267;
  testrel(8, NE, y, f, cnt,failed);

  y := erf_qx(-0.25);
  f := 0.5987063256829237242;
  testrel(9, NE, y, f, cnt,failed);

  y := erf_qx(-1.5);
  f := 0.9331927987311419340;
  testrel(10, NE, y, f, cnt,failed);

  y := erf_qx(-4);
  f := 0.9999683287581668801;
  testrel(11, NE, y, f, cnt,failed);

  y := erf_qx(-8);
  f := 2*0.4999999999999996890; {0.999999999999999377903942572822;}
  testrel(12, NE, y, f, cnt,failed);

  y := erf_qx(-10);
  f := 1.0;
  testrel(13, NE, y, f, cnt,failed);

  y := erf_qx(100);
  f := 0.1344179076744198305e-2173;
{$ifdef FPC}
  testrel(14, NE+1, y, f, cnt,failed);   {FPC271 only}
{$else}
  testrel(14, NE, y, f, cnt,failed);
{$endif}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfcinvx;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfc_invx');

  x := ldexp(1,-16000);
  y := erfc_invx(x);
  f := 105.2859240315651474724078;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1e-4200;
  y := erfc_invx(x);
  f := 98.31427569110739080218294;;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1e-2000;
  y := erfc_invx(x);
  f := 67.82610681154402285494;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1e-800;
  y := erfc_invx(x);
  f := 42.86883824391192220198;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1e-100;
  y := erfc_invx(x);
  f := 15.06557470259264570440;
  testrel( 5, NE, y, f, cnt,failed);

  x := 1e-10;
  y := erfc_invx(x);
  f := 4.572824967389485278741;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := erfc_invx(x);
  f := 2.33146777362194767232;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.125;
  y := erfc_invx(x);
  f := 1.08478704006928314131;
  testrel( 8, NE, y, f, cnt,failed);

  x := 0.2;
  y := erfc_invx(x);
  f := 0.90619380243682322007;
  testrel( 9, NE, y, f, cnt,failed);

  x := 0.5;
  y := erfc_invx(x);
  f := 0.47693627620446987338;
  testrel(10, NE, y, f, cnt,failed);

  x := 0.75;
  y := erfc_invx(x);
  f := 0.225312055012178104725;
  testrel(11, NE, y, f, cnt,failed);

  x := 0.9375;
  y := erfc_invx(x);
  f := 0.55445948772782020299e-1;
  testrel(12, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := erfc_invx(x);
  f := 0.86545619796713755345e-3;
  testrel(13, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := erfc_invx(x);
  f := -0.86545619796713755345e-3;
  testrel(14, NE, y, f, cnt,failed);

  x := 1.125;
  y := erfc_invx(x);
  f := -0.111235451840949959147;
  testrel(15, NE, y, f, cnt,failed);

  x := 1.25;
  y := erfc_invx(x);
  f := -0.225312055012178104725;
  testrel(16, NE, y, f, cnt,failed);

  x := 1.5;
  y := erfc_invx(x);
  f := -0.47693627620446987338;
  testrel(17, NE, y, f, cnt,failed);

  x := 1.75;
  y := erfc_invx(x);
  f := -0.81341984759761854169;
  testrel(18, NE, y, f, cnt,failed);

  x := 1.875;
  y := erfc_invx(x);
  f := -1.08478704006928314131;
  testrel(19, NE, y, f, cnt,failed);

  x := 1.9375;
  y := erfc_invx(x);
  f := -1.31715033498613074888;
  testrel(20, NE, y, f, cnt,failed);

  x := 1.9990234375;
  y := erfc_invx(x);
  f := -2.3314677736219476723;
  testrel(21, NE, y, f, cnt,failed);

  x := 2.0-ldexp(1,-20);
  y := erfc_invx(x);
  f := -3.4655050258033307173;
  testrel(22, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_fresnelx;
var
  x,yc,ys,fs,fc: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Fresnel(CS)x');

  x  := 0;
  fs := 0;
  fc := 0;
  Fresnelx(x,ys,yc);
  testrel( 1, NE, ys, fs, cnt,failed);
  testrel( 2, NE, yc, fc, cnt,failed);

  x  := 1e-300;
  fs := 0.5235987755982988730e-900;
  fc := 1e-300;;
  Fresnelx(x,ys,yc);
  testrel( 3, NE, ys, fs, cnt,failed);
  testrel( 4, NE, yc, fc, cnt,failed);

  x  := 1e-10;
  fs := 0.5235987755982988731e-30;
  fc := 1e-10;
  Fresnelx(x,ys,yc);
  testrel( 5, NE, ys, fs, cnt,failed);
  testrel( 6, NE, yc, fc, cnt,failed);

  x  := 0.00048828125;
  fs := 0.6095491996946437605e-10;
  fc := 0.4882812499999931516e-3;
  Fresnelx(x,ys,yc);
  testrel( 7, NE, ys, fs, cnt,failed);
  testrel( 8, NE, yc, fc, cnt,failed);

  x  := -0.125;
  fs := -0.1022609856621743054e-2;
  fc := -0.1249924702994110995;
  Fresnelx(x,ys,yc);
  testrel( 9, NE, ys, fs, cnt,failed);
  testrel(10, NE, yc, fc, cnt,failed);

  x  := 0.5;
  fs := 0.6473243285999927761e-1;
  fc := 0.4923442258714463929;
  Fresnelx(x,ys,yc);
  testrel(11, NE, ys, fs, cnt,failed);
  testrel(12, NE, yc, fc, cnt,failed);

  x  := 1.0;
  fs := 0.4382591473903547661;
  fc := 0.7798934003768228295;
  Fresnelx(x,ys,yc);
  testrel(13, NE, ys, fs, cnt,failed);
  testrel(14, NE, yc, fc, cnt,failed);

  x  := -1.5;
  fs := -0.6975049600820930131;
  fc := -0.4452611760398215351;
  Fresnelx(x,ys,yc);
  testrel(15, NE, ys, fs, cnt,failed);
  testrel(16, NE, yc, fc, cnt,failed);

  x  := 2.0;
  fs := 0.3434156783636982422;
  fc := 0.4882534060753407545;
  Fresnelx(x,ys,yc);
  testrel(17, NE, ys, fs, cnt,failed);
  testrel(18, NE, yc, fc, cnt,failed);

  x  := 10.0;
  fs := 0.4681699785848822404;
  fc := 0.4998986942055157236;
  Fresnelx(x,ys,yc);
  testrel(19, NE, ys, fs, cnt,failed);
  testrel(20, NE, yc, fc, cnt,failed);

  x  := -100.0;
  fs := -0.4968169011478375533;
  fc := -0.4999998986788178976;
  Fresnelx(x,ys,yc);
  testrel(21, NE, ys, fs, cnt,failed);
  testrel(22, NE, yc, fc, cnt,failed);

  x  := 1000.0;
  fs := 0.4996816901138163061;
  fc := 0.4999999998986788164;
  Fresnelx(x,ys,yc);
  testrel(23, NE, ys, fs, cnt,failed);
  testrel(24, NE, yc, fc, cnt,failed);

  x := 5000.10009765625;
  fs := 0.4999986583294711585;
  fc := 0.5000636465631322710;
  Fresnelx(x,ys,yc);
  testrel(25, NE, ys, fs, cnt,failed);
  testrel(26, NE, yc, fc, cnt,failed);

  x := 6000;
  fs := 0.4999469483523027016;
  fc := 0.4999999999995309204;
  Fresnelx(x,ys,yc);
  testrel(27, NE, ys, fs, cnt,failed);
  testrel(28, NE, yc, fc, cnt,failed);

  x := 1e5;
  fs := 0.4999968169011381621;
  fc := 0.4999999999999998987;
  Fresnelx(x,ys,yc);
  testrel(29, NE, ys, fs, cnt,failed);
  testrel(20, NE, yc, fc, cnt,failed);

  x := 2e6;
  fs := 0.4999998408450569081;
  fc := 0.5;
  Fresnelx(x,ys,yc);
  testrel(31, NE, ys, fs, cnt,failed);
  testrel(32, NE, yc, fc, cnt,failed);

  x  := 1e7;
  fs := 0.4999999681690113816;
  fc := 0.5;
  Fresnelx(x,ys,yc);
  testrel(33, NE, ys, fs, cnt,failed);
  testrel(34, NE, yc, fc, cnt,failed);

  x  := 2e19;
  fs := 0.5;
  fc := 0.5;
  Fresnelx(x,ys,yc);
  testrel(35, NE, ys, fs, cnt,failed);
  testrel(36, NE, yc, fc, cnt,failed);

  x  := 1e-2000;
  fs := 0;
  fc := 1e-2000;
  Fresnelx(x,ys,yc);
  testrel(37, NE, ys, fs, cnt,failed);
  testrel(38, NE, yc, fc, cnt,failed);

  {one test for separate functions}
  x  := 0.75;
  fs := 0.2088771112333835702;
  fc := 0.6935259907871358975;
  ys := FresnelSx(x);
  yc := FresnelCx(x);
  testrel(39, NE, ys, fs, cnt,failed);
  testrel(40, NE, yc, fc, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_gsix;
var
  x,y,t: extended;
  cnt, failed: integer;
const
  NE = 4;
{$ifdef BIT16}
  t5 : THexExtW = ($DD79,$61D8,$E600,$F6AB,$3FFE);  {9.6356046208697728562E-1}
{$endif}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gsix');

  {Test values from MISCFUN [22]}
  x := 1.0/512.0;
  t := 0.59531540040441651584e1;
  y := gsix(x);
  testrel( 1, NE, y, t, cnt,failed);

  x := 1.0/128.0;
  t := 0.45769601268624494109e1;
  y := gsix(x);
  testrel( 2, NE, y, t, cnt,failed);

  x := 1.0/32.0;
  t := 0.32288921331902217638e1;
  y := gsix(x);
  testrel( 2, NE, y, t, cnt,failed);

  x := 1.0/8.0;
  t := 0.19746110873568719362e1;
  y := gsix(x);
  testrel( 4, NE, y, t, cnt,failed);

  x := 1.0/2.0;
  {$ifdef BIT16}
   t := extended(t5);
  {$else}
   t := 0.96356046208697728563;
  {$endif}
  y := gsix(x);
  testrel( 5, NE, y, t, cnt,failed);

  x := 1.0;
  t := 0.60513365250334458174;
  y := gsix(x);
  testrel( 6, NE, y, t, cnt,failed);

  x := 5.0/4.0;
  t := 0.51305506459532198016;
  y := gsix(x);
  testrel( 7, NE, y, t, cnt,failed);

  x := 3.0/2.0;
  t := 0.44598602820946133091;
  y := gsix(x);
  testrel( 8, NE, y, t, cnt,failed);

  x := 15.0/8.0;
  t := 0.37344458206879749357;
  y := gsix(x);
  testrel( 9, NE, y, t, cnt,failed);

  x := 2.0;
  t := 0.35433592884953063055;
  y := gsix(x);
  testrel(10, NE, y, t, cnt,failed);

  x := 17.0/8.0;
  t := 0.33712156518881920994;
  y := gsix(x);
  testrel(11, NE, y, t, cnt,failed);

  x := 5.0/2.0;
  t := 0.29436170729362979176;
  y := gsix(x);
  testrel(12, NE, y, t, cnt,failed);

  x := 3.0;
  t := 0.25193499644897222840;
  y := gsix(x);
  testrel(13, NE, y, t, cnt,failed);

  x := 7.0/2.0;
  t := 0.22028778222123939276;
  y := gsix(x);
  testrel(14, NE, y, t, cnt,failed);

  x := 4.0;
  t := 0.19575258237698917033;
  y := gsix(x);
  testrel(15, NE, y, t, cnt,failed);

  x := 9.0/2.0;
  t := 0.17616303166670699424;
  y := gsix(x);
  testrel(16, NE, y, t, cnt,failed);

  x := 5.0;
  t := 0.16015469479664778673;
  y := gsix(x);
  testrel(17, NE, y, t, cnt,failed);

  x := 23.0/4.0;
  t := 0.14096116876193391066;
  y := gsix(x);
  testrel(18, NE, y, t, cnt,failed);

  x := 6.0;
  t := 0.13554987191049066274;
  y := gsix(x);
  testrel(19, NE, y, t, cnt,failed);

  x := 7.0;
  t := 0.11751605060085098084;
  y := gsix(x);
  testrel(20, NE, y, t, cnt,failed);

  {Test values calculated with Maple}
  x := 100;
  t := 0.88127074334736483777e-2;
  y := gsix(x);
  testrel(21, NE, y, t, cnt,failed);

  x := 1000;
  t := 0.88572736806688441188e-3;
  y := gsix(x);
  testrel(22, NE, y, t, cnt,failed);

  x := 1e18;
  t := 0.88622692545275801315e-18;
  y := gsix(x);
  testrel(23, NE, y, t, cnt,failed);

  x := 1e19;
  t := 0.886226925452758013599e-19;
  y := gsix(x);
  testrel(24, NE, y, t, cnt,failed);

  x := 1e20;
  t := 0.88622692545275801364e-20;
  y := gsix(x);
  testrel(25, NE, y, t, cnt,failed);

  x := 1e-20;
  t := 45.7630940274301472501;
  y := gsix(x);
  testrel(26, NE, y, t, cnt,failed);

  x := 1e-19;
  t := 43.4605089344361015662;
  y := gsix(x);
  testrel(27, NE, y, t, cnt,failed);

  x := 1e-18;
  t := 41.1579238414420558838;
  y := gsix(x);
  testrel(28, NE, y, t, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_expint3x;
var
  x,y,t: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','expint3x');

  {Test values from MISCFUN [22]}
  x := 1.0/512.0;
  t := 0.19531249963620212007e-2;
  y := expint3x(x);
  testrel( 1, NE, y, t, cnt,failed);

  x := 1.0/128.0;
  t := 0.78124990686775522671e-2;
  y := expint3x(x);
  testrel( 2, NE, y, t, cnt,failed);

  x := 1.0/32.0;
  t := 0.31249761583499728667e-1;
  y := expint3x(x);
  testrel( 3, NE, y, t, cnt,failed);

  x := 1.0/8.0;
  t := 0.12493899888803079984;
  y := expint3x(x);
  testrel( 4, NE, y, t, cnt,failed);

  x := 1.0/2.0;
  t := 0.48491714311363971332;
  y := expint3x(x);
  testrel( 5, NE, y, t, cnt,failed);

  x := 1.0;
  t := 0.80751118213967145286;
  y := expint3x(x);
  testrel( 6, NE, y, t, cnt,failed);

  x := 5.0/4.0;
  t := 0.86889265412623270696;
  y := expint3x(x);
  testrel( 7, NE, y, t, cnt,failed);

  x := 3.0/2.0;
  t := 0.88861722235357162648;
  y := expint3x(x);
  testrel( 8, NE, y, t, cnt,failed);

  x := 15.0/8.0;
  t := 0.89286018500218176869;
  y := expint3x(x);
  testrel( 9, NE, y, t, cnt,failed);

  x := 2.0;
  t := 0.89295351429387631138;
  y := expint3x(x);
  testrel(10, NE, y, t, cnt,failed);

  x := 17.0/8.0;
  t := 0.89297479112737843939;
  y := expint3x(x);
  testrel(11, NE, y, t, cnt,failed);

  x := 18.0/8.0;
  t := 0.89297880579798112220;
  y := expint3x(x);
  testrel(12, NE, y, t, cnt,failed);

  x := 5.0/2.0;
  t := 0.89297950317496621294;
  y := expint3x(x);
  testrel(13, NE, y, t, cnt,failed);

  x := 11.0/4.0;
  t := 0.89297951152951902903;
  y := expint3x(x);
  testrel(14, NE, y, t, cnt,failed);

  x := 3.0;
  t := 0.89297951156918122102;
  y := expint3x(x);
  testrel(15, NE, y, t, cnt,failed);

  x := 25.0/8.0;
  t := 0.89297951156924734716;
  y := expint3x(x);
  testrel(16, NE, y, t, cnt,failed);

  x := 13.0/4.0;
  t := 0.89297951156924917298;
  y := expint3x(x);
  testrel(17, NE, y, t, cnt,failed);

  x := 7.0/2.0;
  t := 0.89297951156924921121;
  y := expint3x(x);
  testrel(18, NE, y, t, cnt,failed);

  x := 15.0/4.0;
  t := 0.89297951156924921122;
  y := expint3x(x);
  testrel(19, NE, y, t, cnt,failed);

  x := 4.0;
  t := 0.89297951156924921122;
  y := expint3x(x);
  testrel(20, NE, y, t, cnt,failed);

  {Test values calculated with Maple}
  x := 1e-10;
  t := 1e-10;
  y := expint3x(x);
  testrel(21, NE, y, t, cnt,failed);

  x := 1e-5;
  t := 0.99999999999999975e-5;
  y := expint3x(x);
  testrel(22, NE, y, t, cnt,failed);

  x := 1e-3;
  t := 0.99999999975e-3;
  y := expint3x(x);
  testrel(23, NE, y, t, cnt,failed);

  x := PosInf_x;
  t := 0.89297951156924921122;
  y := expint3x(x);
  testrel(24, NE, y, t, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.
