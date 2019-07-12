{Part 1 of regression test for SPECFUNX unit  (c) 2010-2018  W.Ehrhardt}

unit t_sfx1;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_gammax;
procedure test_gamma1pm1x;
procedure test_gammastarx;
procedure test_inv_gammax;
procedure test_facx;
procedure test_dfacx;
procedure test_lnfacx;
procedure test_binomialx;
procedure test_lnbinomialx;
procedure test_lngammax;
procedure test_lngamma1px;
procedure test_lngamma_invx;
procedure test_rgammax;
procedure test_gamma_ratiox;
procedure test_pochhammerx;
procedure test_poch1x;
procedure test_lnBarnesGx;


implementation


uses
  amath, specfunx, t_sfx0;

{---------------------------------------------------------------------------}
procedure test_gammax;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gammax');

  x := PosInf_x;
  y := gammax(x);
  f := PosInf_d;
  testabs( 1, 0, y, f, cnt,failed);

  x := MaxExtended;
  y := gammax(x);
  f := PosInf_d;
  testabs( 2, 0, y, f, cnt,failed);

  x := -12.75;
  y := gammax(x);
  f := -0.13645255983702753328e-8;
  testrel( 3, NE, y, f, cnt,failed);

  x := +12.75;
  y := gammax(x);
  f := 255371835.69921110046;
  testrel( 4, NE, y, f, cnt,failed);

  x := -13.125;
  y := gammax(x);
  f := 0.95164579566876706501e-9;
  testrel( 5, NE, y, f, cnt,failed);

  x := +13.125;
  y := gammax(x);
  f := 657257524.55014748218;
  testrel( 6, NE, y, f, cnt,failed);

  x := -171.5;
  y := gammax(x);
  f := 0.19316265431711996005e-309;
  testrel( 7, NE, y, f, cnt,failed);

  x := +171.5;
  y := gammax(x);
  f := 0.94833675668247993363e308;
  testrel( 8, NE, y, f, cnt,failed);

  x := -7.875;
  y := gammax(x);
  f := 0.26582609556918071863e-3;
  testrel( 9, NE, y, f, cnt,failed);

  x := +7.875;
  y := gammax(x);
  f := 3921.5886522261462915;
  testrel(10, NE, y, f, cnt,failed);

  x := +0.03;
  y := gammax(x);
  f := 32.784998351794135982;
  testrel(11, NE, y, f, cnt,failed);

  x := -0.03;
  y := gammax(x);
  f := -33.941064736219641591;
  testrel(12, NE, y, f, cnt,failed);

  x := -1.015625;
  y := gammax(x);
  f := 63.599158175314632272;
  testrel(13, NE, y, f, cnt,failed);

  x := -0.984375;
  y := gammax(x);
  f := -64.444973175762601638;
  testrel(14, NE, y, f, cnt,failed);

  x := -1.984375;
  y := gammax(x);
  f := 32.476206954715011849;
  testrel(15, NE, y, f, cnt,failed);

  x := -2.015625;
  y := gammax(x);
  f := -31.553070722636716786;
  testrel(16, NE, y, f, cnt,failed);

  x := 500;
  y := gammax(x);
  f := 0.24402736519822201374e1132;
  testrel(17, 2*NE, y, f, cnt,failed);

  x := 1000.5;
  y := gammax(x);
  f := 0.12723011956950554642e2567;
  testrel(18, 2*NE, y, f, cnt,failed);

  x := 1500.0;
  y := gammax(x);
  f := 0.32079985311865165734e4112;
  testrel(19, 2*NE, y, f, cnt,failed);

  x := 1500.125;
  y := gammax(x);
  f := 0.80025746340459573619e4112;
  testrel(20, 3*NE, y, f, cnt,failed);

  x := 1755.125;
  y := gammax(x);
  f := 0.50352886147710913730e4931;
  testrel(21, 3*NE, y, f, cnt,failed);

  x := 1755.4375;
  y := gammax(x);
  f := 0.51980930952591760081e4932;
  testrel(22, 3*NE, y, f, cnt,failed);

  x := 1755.4375;
  y := gammax(x);
  f := 0.51980930952591760081e4932;
  testrel(23, 3*NE, y, f, cnt,failed);

  {extreme case x < -MAXGAMX without underflow}
  x := -1756.0 - ldexp(1,-24);
  y := gammax(x);
  f := -2.750520002448301157e-4930;
  testrel(25, 6500, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_gamma1pm1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gamma1pm1x');

  x := 0.0;
  y := gamma1pm1x(x);
  f := 0.0;
  testabs( 1, 0, y, f, cnt,failed);

  x := 1e-20;
  y := gamma1pm1x(x);
  f := -0.57721566490153286060e-20;
  testrel( 2, NE, y, f, cnt,failed);

  x := -1e-20;
  y := gamma1pm1x(x);
  f := 0.5772156649015328606e-20;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1e-10;
  y := gamma1pm1x(x);
  f := -0.57721566480262726108279e-10;
  testrel( 4, NE, y, f, cnt,failed);

  x := -1e-10;
  y := gamma1pm1x(x);
  f := 0.5772156650004384601484e-10;
  testrel( 5, NE, y, f, cnt,failed);

  x := 1e-5;
  y := gamma1pm1x(x);
  f := -0.57720577443232650677e-5;
  testrel( 6, NE, y, f, cnt,failed);

  x := -1e-5;
  y := gamma1pm1x(x);
  f := 0.5772255555522350297e-5;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.125;
  y := gamma1pm1x(x);
  f := -0.582573001502985119125963e-1;
  testrel( 8, NE, y, f, cnt,failed);

  x := -0.125;
  y := gamma1pm1x(x);
  f := 0.89652357422896951252377e-1;
  testrel( 9, NE, y, f, cnt,failed);

  x := 0.5;
  y := gamma1pm1x(x);
  f := -0.11377307454724198635;
  testrel(10, NE, y, f, cnt,failed);

  x := -0.5;
  y := gamma1pm1x(x);
  f := 0.772453850905516027298167;
  testrel(11, NE, y, f, cnt,failed);

  x := 0.625;
  y := gamma1pm1x(x);
  f := -0.10342571994340201523;
  testrel(12, NE, y, f, cnt,failed);

  x := -0.625;
  y := gamma1pm1x(x);
  f := 1.370436184416600908646474;
  testrel(13, NE, y, f, cnt,failed);

  x := 1.5;
  y := gamma1pm1x(x);
  f := 0.329340388179137020473626;
  testrel(13, NE, y, f, cnt,failed);

  x := -1.5;
  y := gamma1pm1x(x);
  f := -4.5449077018110320546;
  testrel(15, NE, y, f, cnt,failed);

  x := 2.5;
  y := gamma1pm1x(x);
  f := 2.323350970447842551184064;
  testrel(16, NE, y, f, cnt,failed);

  x := -2.5;
  y := gamma1pm1x(x);
  f := 1.3632718012073547031;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_gammastarx;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gammastarx');

  x := MinDouble;
  y := gammastarx(x);
  f := 0.26744707353778560993e154;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1e-20;
  y := gammastarx(x);
  f := 3989422804.0143267812;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1e-10;
  y := gammastarx(x);
  f := 39894.22813368978815;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1e-4;
  y := gammastarx(x);
  f := 39.93267749222754497;
  testrel( 4, NE, y, f, cnt,failed);

  x := 0.125;
  y := gammastarx(x);
  f := 1.561566128548903816;
  testrel( 5, NE, y, f, cnt,failed);

  x := 0.5;
  y := gammastarx(x);
  f := 1.165821990798562102;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.625;
  y := gammastarx(x);
  f := 1.133875505693637387;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.75;
  y := gammastarx(x);
  f := 1.112114253317550701;
  testrel( 8, NE, y, f, cnt,failed);

  x := 1.0;
  y := gammastarx(x);
  f := 1.084437551419227547;
  testrel( 9, NE, y, f, cnt,failed);

  x := 1.5;
  y := gammastarx(x);
  f := 1.056344244268559867;
  testrel(10, NE, y, f, cnt,failed);

  x := 2.0;
  y := gammastarx(x);
  f := 1.042207120816673058;
  testrel(11, NE, y, f, cnt,failed);

  x := 4.0;
  y := gammastarx(x);
  f := 1.021008303746349370;
  testrel(12, NE, y, f, cnt,failed);

  x := 10;
  y := gammastarx(x);
  f := 1.008365359132400246;
  testrel(13, NE, y, f, cnt,failed);

  x := 1000;
  y := gammastarx(x);
  f := 1.000083336802874000;
  testrel(14, NE, y, f, cnt,failed);

  x := 1e8;
  y := gammastarx(x);
  f := 1.000000000833333334;
  testrel(15, NE, y, f, cnt,failed);

  x := 2e8;
  y := gammastarx(x);
  f := 1.000000000416666667;
  testrel(16, NE, y, f, cnt,failed);

  x := 3e8;
  y := gammastarx(x);
  f := 1+1/(12*x);
  testrel(17, NE, y, f, cnt,failed);

  x := MinExtended;
  y := gammastarx(x);
  f := 0.2175727742838104741e2466;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_facx;
var
  y,f: extended;
  n, cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','facx');
  for n:=0 to trunc(MAXGAMX)-1 do begin
    y := facx(n);
    f := gammax(n+1);
    testrel(n, NE, y, f, cnt,failed);
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_dfacx;
var
  y,f: extended;
  n, cnt, failed: integer;
const
  NE = 150;  {needed for large n with n/2 multiplications}

  function dfbrute(k: integer): extended;
    {-dfac brute force from definition}
  var
    x: extended;
  begin
    if k<0 then x := 0
    else begin
      x := 1;
      while k>0 do begin
        x := x*k;
        dec(k,2);
      end;
    end;
    dfbrute := x;
  end;

begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','dfacx');
  for n:=0 to MAXDFACX do begin
    y := dfacx(n);
    f := dfbrute(n);
    testrel(n, NE, y, f, cnt,failed);
  end;

  n := -1;
  y := dfacx(n);
  f := 1;
  testrel(n, NE, y, f, cnt,failed);
  n := -3;
  y := dfacx(n);
  f := -1;
  testrel(n, NE, y, f, cnt,failed);
  n := -5;
  y := dfacx(n);
  f := 1/3;
  testrel(n, NE, y, f, cnt,failed);
  n := -7;
  y := dfacx(n);
  f := -1/15;
  testrel(n, NE, y, f, cnt,failed);
  n := -9;
  y := dfacx(n);
  f := 1/105;
  testrel(n, NE, y, f, cnt,failed);
  n := -97;
  y := dfacx(n);
  f := 3.523529645491056559e-75;
  testrel(n, NE, y, f, cnt,failed);
  n := -99;
  y := dfacx(n);
  f := -3.632504789166037690e-77;
  testrel(n, NE, y, f, cnt,failed);
  n := -3209;
  y := dfacx(n);
  f := 2.87033402181645533e-4928;
  testrel(n, NE, y, f, cnt,failed);
  n := -3211;
  y := dfacx(n);
  f := -8.944637026539281199e-4932;
  testrel(n, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lnfacx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lnfacx');

  y := lnfacx(-1);
  f := PosInf_x;
  testabs(1, 0, y, f, cnt,failed);

  y := lnfacx(0);
  f := 0;
  testabs(2, 0, y, f, cnt,failed);

  f := 0;
  y := lnfacx(1);
  testabs(3, 0, y, f, cnt,failed);

  f := 0.6931471805599453094;
  y := lnfacx(2);
  testrel(4, NE, y, f, cnt,failed);

  f := 15.10441257307551530;
  y := lnfacx(10);
  testrel(5, NE, y, f, cnt,failed);

  f := 58.00360522298051994;
  y := lnfacx(25);
  testrel(6, NE, y, f, cnt,failed);

  f := 61.26170176100200198;
  y := lnfacx(26);
  testrel(7, NE, y, f, cnt,failed);

  f := 363.7393755555634901;
  y := lnfacx(100);
  testrel(8, NE, y, f, cnt,failed);

  f := 1051299.221899121865;
  y := lnfacx(100000);
  testrel(9, NE, y, f, cnt,failed);

  f := 4.399670565537852434e10;
  y := lnfacx(2147483647);
  testrel(10, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lngammax;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lngammax');

  x := MaxExtended;
  y := lngammax(x);
  f := PosInf_x;
  testabs( 1, 0, y, f, cnt,failed);

  x := -8.0625;
  y := lngammax(x);
  f := -7.9595994278861808680;
  testrel( 2, NE, y, f, cnt,failed);

  x := 8.0625;
  y := lngammax(x);
  f := 8.6513982694450757101;
  testrel( 3, NE, y, f, cnt,failed);

  x := -800.1;
  y := lngammax(x);
  f := -4550.3001717121540297;
  testrel( 4, NE, y, f, cnt,failed);

  x := +800.1;
  y := lngammax(x);
  f := 4545.9345238837669004;
  testrel( 5, NE, y, f, cnt,failed);

  x := -800000000.25;
  y := lngammax(x);
  f := -15600097843.3084878921;
  testrel( 6, NE, y, f, cnt,failed);

  x := +800000000.25;
  y := lngammax(x);
  f := 15600097824.299669082;
  testrel( 7, NE, y, f, cnt,failed);

  x := -7.9;
  y := lngammax(x);
  f := -8.0720397360128358026;
  testrel( 8, NE, y, f, cnt,failed);

  x := +7.9;
  y := lngammax(x);
  f := 8.3242658680088089235;
  testrel( 9, NE, y, f, cnt,failed);

  x := +0.03;
  y := lngammax(x);
  f := 3.4899710434424119167;
  testrel(10, NE, y, f, cnt,failed);

  x := -0.03;
  y := lngammax(x);
  f := 3.5246256304460036695;
  testrel(11, NE, y, f, cnt,failed);

  x := -1.015625;
  y := lngammax(x);
  f := 4.1526002340212853834;
  testrel(12, NE, y, f, cnt,failed);

  x := -0.984375;
  y := lngammax(x);
  f := 4.1658117306220779559;
  testrel(13, NE, y, f, cnt,failed);

  x := -1.984375;
  y := lngammax(x);
  f := 3.4805077275231585394;
  testrel(14, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := lngammax(x);
  f := 0.5644719118551233843e-3;
  testrel(15, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := lngammax(x);
  f := -0.5629031799912046317e-3;
  testrel(16, NE, y, f, cnt,failed);

  x := 1.9990234375;
  y := lngammax(x);
  f := -0.4125677359714894017e-3;
  testrel(17, NE, y, f, cnt,failed);

  x := 2.0009765625;
  y := lngammax(x);
  f := 0.4131827930642542643e-3;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lngamma1px;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lngamma1x');

  x := -0.001;
  y := lngamma1px(x);
  f := 0.5780385328913797240e-3;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-8;
  y := lngamma1px(x);
  f := -0.5772156566768625664e-8;
  testrel(2, NE, y, f, cnt,failed);

  x := ldexp(1,-40);
  y := lngamma1px(x);
  f := -0.5249745890076017815e-12;
  testrel(3, NE, y, f, cnt,failed);

  x := -1e-12;
  y := lngamma1px(x);
  f := 0.5772156649023553276e-12;
  testrel(4, NE, y, f, cnt,failed);

  x := 1e-20;
  y := lngamma1px(x);
  f := -0.5772156649015328606e-20;
  testrel(5, NE, y, f, cnt,failed);

  x := -0.99;
  y := lngamma1px(x);
  f := 4.599479878042021723;
  testrel(6, NE, y, f, cnt,failed);

  x := 7.5;
  y := lngamma1px(x);
  f := 9.549267257300997712;
  testrel(7, NE, y, f, cnt,failed);

  x := 15;
  y := lngamma1px(x);
  f := 27.89927138384089157;
  testrel(8, NE, y, f, cnt,failed);

  x := -1.25;
  y := lngamma1px(x);
  f := 1.589575312551185990;
  testrel(9, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lngamma_invx;
var
  x,y,f,r: extended;
  fc, cnt, failed: integer;
const
  NE  = 3;
  NE1 = 5;
  NEF = NE1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lngamma_invx');

  x := -0.12142;
  y := lngamma_invx(x);
  f := 1.4733581679374343597;
  testrel(1, NE1, y, f, cnt,failed);   {AMD}

  x := -0.12;
  y := lngamma_invx(x);
  f := 1.517523824137325410;
  testrel(2, NE, y, f, cnt,failed);

  x := 0;
  y := lngamma_invx(x);
  f := 2;
  testrel(3, NE, y, f, cnt,failed);

  x := 1e-10;
  y := lngamma_invx(x);
  f := 2.000000000236527212;
  testrel(4, NE, y, f, cnt,failed);

  x := 0.15;
  y := lngamma_invx(x);
  f := 2.292995833067628689;
  testrel(5, NE, y, f, cnt,failed);

  x := 0.5;
  y := lngamma_invx(x);
  f := 2.780026731309318306;
  testrel(6, NE, y, f, cnt,failed);

  x := 1;
  y := lngamma_invx(x);
  f := 3.312440882539160370;
  testrel(7, NE, y, f, cnt,failed);

  x := 2;
  y := lngamma_invx(x);
  f := 4.162830474448703075;
  testrel(8, NE, y, f, cnt,failed);

  x := 3;
  y := lngamma_invx(x);
  f := 4.880725029962703745;
  testrel(9, NE, y, f, cnt,failed);

  x := 10;
  y := lngamma_invx(x);
  f := 8.715310138772789661;
  testrel(10, NE, y, f, cnt,failed);

  x := 20;
  y := lngamma_invx(x);
  f := 13.00506116612975588;
  testrel(11, NE, y, f, cnt,failed);

  x := 1000;
  y := lngamma_invx(x);
  f := 226.5074413075941081;
  testrel(12, NE, y, f, cnt,failed);

  x := 1e20;
  y := lngamma_invx(x);
  f := 0.2419543491507313264e19;
  testrel(13, NE, y, f, cnt,failed);

  x := 1e300;
  y := lngamma_invx(x);
  f := 0.1463595972213524998e298;
  testrel(14, NE1, y, f, cnt,failed);   {AMD}

  x := 1e4000;
  y := lngamma_invx(x);
  f := 0.1086931310837727672e3997;
  testrel(15, NE1, y, f, cnt,failed);

  x := 0.005;
  r := 0.0;
  fc:= 0;
  repeat
    if x<2 then x := 1.01*x
    else x := 1.25*x;
    inc(fc);
    y := lngamma_invx(x);
    f := lngammax(y);
    if (x < 0.5) or (f=0.0) then f := abs(x-f)
    else f := abs(1.0 - x/f);
    if f>r then r := f;
  until x>1e300;
  inc(cnt);
  if r > NEF*eps_x then begin
    inc(failed);
    writeln('Test ', cnt:2, ' failed: Max. rel. error = ', r:20, ' eps > ',NEF, ' eps ');
  end
  else writeln('  [',fc,' tests OK: lngamma(lngamma_inv(x) = x]');

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_inv_gammax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
  NE1 = 8;  {AMD}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','inv_gammax');

  y := inv_gammax(0.8857421875);
  f := 1.479691445431194499;
  testrel(1, NE1, y, f, cnt,failed);

  y := inv_gammax(0.9375);
  f := 1.822325750344712207;
  testrel(2, NE, y, f, cnt,failed);

  y := inv_gammax(1);
  f := 2;
  testrel(3, NE, y, f, cnt,failed);

  y := inv_gammax(2);
  f := 3;
  testrel(4, NE, y, f, cnt,failed);

  y := inv_gammax(3);
  f := 3.405869986309566925;
  testrel(5, NE, y, f, cnt,failed);

  y := inv_gammax(10);
  f := 4.390077650833141892;
  testrel(6, NE, y, f, cnt,failed);

  y := inv_gammax(50);
  f := 5.471527458207833634;
  testrel(7, NE, y, f, cnt,failed);

  y := inv_gammax(100);
  f := 5.892518696343772809;
  testrel(8, NE, y, f, cnt,failed);

  y := inv_gammax(1e4);
  f := 8.336251344938071494;
  testrel(9, NE, y, f, cnt,failed);

  y := inv_gammax(1e8);
  f := 12.37351687346275523;
  testrel(10, NE, y, f, cnt,failed);

  y := inv_gammax(1e20);
  f := 22.21852188335517397;
  testrel(11, NE, y, f, cnt,failed);

  y := inv_gammax(1e300);
  f := 167.9203488839770710;
  testrel(12, NE, y, f, cnt,failed);

  y := inv_gammax(1e4900);
  f := 1745.657863337063577;
  testrel(13, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_rgammax;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','rgammax');

  x := PosInf_d;
  y := rgammax(x);
  f := 0;
  testabs( 1, 0, y, f, cnt,failed);

  x := 1755.5;
  y := rgammax(x);
  f := 0;
  testabs( 2, 0, y, f, cnt,failed);

  x := 10.0;
  y := rgammax(x);
  f := 0.27557319223985890653e-5;
  testrel( 3, NE, y, f, cnt,failed);

  x := 8.00390625;
  y := rgammax(x);
  f := 0.19685641040785775363e-3;
  testrel( 4, NE, y, f, cnt,failed);

  x := 8.0;
  y := rgammax(x);
  f := 0.19841269841269841270e-3;
  testrel( 5, NE, y, f, cnt,failed);

  x := 7.99609375;
  y := rgammax(x);
  f := 0.19998088370232038881e-3;
  testrel( 6, NE, y, f, cnt,failed);

  x := 3.125;
  y := rgammax(x);
  f := 0.44417721917199933367;
  testrel( 7, NE, y, f, cnt,failed);

  x := 2.875;
  y := rgammax(x);
  f := 0.55937456141092137994;
  testrel( 8, NE, y, f, cnt,failed);

  x := 0.125;
  y := rgammax(x);
  f := 0.13273264557288261338;
  testrel( 9, NE, y, f, cnt,failed);

  x := 0.015625;
  y := rgammax(x);
  f := 0.15763417467835646516e-1;
  testrel(10, NE, y, f, cnt,failed);

  x := 0;
  y := rgammax(x);
  f := 0;
  testrel(11, NE, y, f, cnt,failed);

  x := -0.015625;
  y := rgammax(x);
  f := -0.15481578889790291018e-1;
  testrel(12, NE, y, f, cnt,failed);

  x := -0.125;
  y := rgammax(x);
  f := -0.11471548622684911112;
  testrel(13, NE, y, f, cnt,failed);

  x := -0.5;
  y := rgammax(x);
  f := -0.28209479177387814347;
  testrel(14, NE, y, f, cnt,failed);

  x := -0.515625;
  y := rgammax(x);
  f := -0.28194771040097295871;
  testrel(15, NE, y, f, cnt,failed);

  x := -1;
  y := rgammax(x);
  f := 0;
  testrel(16, NE, y, f, cnt,failed);

  x := -5.5;
  y := rgammax(x);
  f := 91.636730015295728167;
  testrel(17, NE, y, f, cnt,failed);

  x := -7.99609375;
  y := rgammax(x);
  f := 156.18471453821053044;
  testrel(18, NE, y, f, cnt,failed);

  x := -17.99609375;
  y := rgammax(x);
  f := 24725224142466.918422;
  testrel(19, NE, y, f, cnt,failed);

  x := 2;
  y := rgammax(x);
  f := 1;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_gamma_ratiox;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 8;
  NE1 = 400;
  NE2 = 7000;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gamma_delta/ratiox');

  y := gamma_delta_ratiox(-4, -2);
  f := 30.0;
  testrel(1, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(-5, -2);
  f := 42.0;
  testrel(2, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(-5, -3);
  f := -336;
  testrel(3, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(-4, -3);
  f := -210;
  testrel(4, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(-4.5, -2.25);
  f := 41.30860527437135181;
  testrel(5, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(-4.5, 2.25);
  f := 0.3443831153885897958e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox( 4.5, 2.25);
  f := 0.2567649439468911646e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox( 4.5, -2.25);
  f := 10.26628120819269988;
  testrel(8, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(-5.5, -11);
  f := -0.2974621527196655273e12;
  testrel(9, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(20,7);
  f := 0.3016307364133451090e-9;
  testrel(10, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(20,10);
  f := 0.1375801570942095918e-13;
  testrel(11, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(4,1/1024);
  f := 0.99877393949161112781;
  testrel(12, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(4,-1/1024);
  f := 1.001227294571164676;
  testrel(13, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(4,0);
  f := 1.0;
  testrel(14, NE, y, f, cnt,failed);

  y := gamma_ratiox(100,90);
  f := 0.565340858599765248e20;
  testrel(15, NE, y, f, cnt,failed);

  y := gamma_ratiox(1-ldexp(1,-20), 1+ldexp(1,-20));
  f := 1.000001100952115336;
  testrel(16, NE, y, f, cnt,failed);

  y := gamma_ratiox(8-ldexp(1,-20), 8+ldexp(1,-20));
  f := 0.9999961554763729811;
  testrel(17, NE, y, f, cnt,failed);

  y := gamma_ratiox(-4.5, -2.25);
  f := 0.3443831153885897958e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := gamma_ratiox(-4.5, 2.25);
  f := -0.5297390756961110946e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := gamma_ratiox( 4.5, 2.25);
  f := 10.26628120819269988;
  testrel(20, NE, y, f, cnt,failed);

  y := gamma_ratiox( 4.5, -2.25);
  f := -6.674104418834535643;
  testrel(21, NE, y, f, cnt,failed);

  y := gamma_ratiox(2020, 2012);
  f := 0.2723078585154254231e27;
  testrel(22, NE, y, f, cnt,failed);

  y := gamma_ratiox(2000, 2012);
  f := 0.2362303670217364563e-39;
  testrel(23, NE1, y, f, cnt,failed);

  y := gamma_delta_ratiox(1e-40, 10);
  f := 0.2755731922398589065e35;
  testrel(24, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(1e-50, 170);
  f := 0.2342431645246009997e-254;
  testrel(25, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(1e-60, 190);
  f := 0.1962744490797566428e-289;
  testrel(26, NE, y, f, cnt,failed);

  y := gamma_delta_ratiox(1e-1000, 1770);
  f := 0.1031401503968751336e-3978;
  testrel(27, NE2, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_pochhammerx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','pochhammerx');

  y := pochhammerx(4,0.1);
  f := 1.135437143836113145;
  testrel( 1, NE, y, f, cnt,failed);

  y := pochhammerx(0.1,4);
  f := 0.7161;
  testrel( 2, NE, y, f, cnt,failed);

  y := pochhammerx(-4,-0.1);
  f := 0;
  testrel( 3, NE, y, f, cnt,failed);

  y := pochhammerx(4-1/32, 1/8);
  f := 1.171297529007599674;
  testrel( 4, NE, y, f, cnt,failed);

  y := pochhammerx(-2,-6);
  f := 0.4960317460317460318e-4;
  testrel( 5, NE, y, f, cnt,failed);

  y := pochhammerx(-2,-5);
  f := -0.3968253968253968254e-3;
  testrel( 6, NE, y, f, cnt,failed);

  y := pochhammerx(-2,-6);
  f := 0.4960317460317460318e-4;
  testrel( 7, NE, y, f, cnt,failed);

  y := pochhammerx(1500,90);
  f := 0.9668996796991433869e287;
  testrel( 8, NE, y, f, cnt,failed);

  y := pochhammerx(2000,20);
  f := 0.1152721650262556377e67;
  testrel( 9, NE, y, f, cnt,failed);

  y := pochhammerx(-4, 2);
  f := 12;
  testrel(10, NE, y, f, cnt,failed);

  y := pochhammerx(-5, 2);
  f := 20;
  testrel(11, NE, y, f, cnt,failed);

  y := pochhammerx(-5, 3);
  f := -60;
  testrel(12, NE, y, f, cnt,failed);

  y := pochhammerx(-4, 3);
  f := -24;
  testrel(13, NE, y, f, cnt,failed);

  y := pochhammerx(-4, -2);
  f := 1/30;   {Fix311}
  testrel(14, NE, y, f, cnt,failed);

  y := pochhammerx(-5, -2);
  f := 1/42;   {Fix311}
  testrel(15, NE, y, f, cnt,failed);

  y := pochhammerx(-5, -3);
  f := -1/336;
  testrel(16, NE, y, f, cnt,failed);

  y := pochhammerx(-4, -3);
  f := -1/210;
  testrel(17, NE, y, f, cnt,failed);

  y := pochhammerx(-5.5,-11.0);
  f := -0.3361772214909036314e-11;
  testrel(18, NE, y, f, cnt,failed);

  y := pochhammerx(-5.5,-11.25);
  f := -0.2337013293452401502e-11;
  testrel(19, NE, y, f, cnt,failed);

  y := pochhammerx(-5.5, 11.25);
  f := 7219.552220411654733;
  testrel(20, NE, y, f, cnt,failed);

  y := pochhammerx( 5.5, 11.25);
  f := 0.1987017505300614113e12;
  testrel(21, NE, y, f, cnt,failed);

  y := pochhammerx( 5.5,-11.25);
  f := 0.1873697945363471336e-3;
  testrel(22, NE, y, f, cnt,failed);

  y := pochhammerx(-5.5, 6);
  f := 162.421875;
  testrel(23, NE, y, f, cnt,failed);

  y := pochhammerx(-5.5, 1);
  f := -5.5;
  testrel(24, NE, y, f, cnt,failed);

  y := pochhammerx(-5.5, 2);
  f := 24.75;
  testrel(25, NE, y, f, cnt,failed);

  y := pochhammerx(-5.5, 100);
  f := 0.1026370005605834561e148;
  testrel(26, NE, y, f, cnt,failed);

  {-------------------------}
  {Values from Wolfram Alpha}
  y := pochhammerx(0,-2);
  f := 0.5;
  testrel(27, NE, y, f, cnt,failed);

  y := pochhammerx(0,0.1);
  f := 0;
  testabs(28, 0, y, f, cnt,failed);

  y := pochhammerx(0,1);
  f := 0;
  testabs(29, 0, y, f, cnt,failed);

  y := pochhammerx(0,-0.1);
  f := 0;
  testabs(30, 0, y, f, cnt,failed);

  y:= pochhammerx(0,-1);
  f := -1;
  testrel(31, NE, y, f, cnt,failed);

  y := pochhammerx(0,-3);
  f := -1/6;
  testrel(32, NE, y, f, cnt,failed);

  y := pochhammerx(0,-4);
  f := 1/24;
  testrel(33, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_poch1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','poch1x');

  y := poch1x(5,2);
  f := 29/2;
  testrel( 1, NE, y, f, cnt,failed);

  y := poch1x(-5,2);
  f := 19/2;
  testrel( 2, NE, y, f, cnt,failed);

  y := poch1x(-5,-2);
  f := 41/84;
  testrel( 3, NE, y, f, cnt,failed);

  y := poch1x(5,-2);
  f := 11/24;
  testrel( 4, NE, y, f, cnt,failed);

  y := poch1x(5,0.01);
  f := 1.518639366136827533;
  testrel( 5, NE, y, f, cnt,failed);

  y := poch1x(6.25,-1/8);
  f := 1.563402663432698836;
  testrel( 6, NE, y, f, cnt,failed);

  y := poch1x(6.25, 1/128);
  f := 1.763164043235102186;
  testrel( 7, NE, y, f, cnt,failed);

  y := poch1x(6.25,-1/256);
  f := 1.744146003990379491;
  testrel( 8, NE, y, f, cnt,failed);

  x := ldexp(1,-30);
  y := poch1x(6.25, x);
  f := 1.750453528391345545;
  testrel( 9, NE, y, f, cnt,failed);

  y := poch1x(6.25, 0);
  f := 1.750453526883736028;
  testrel(10, NE, y, f, cnt,failed);

  y := poch1x(-5.5,-1/8);
  f := 1.088331930355213549;
  testrel(11, NE, y, f, cnt,failed);

  y := poch1x(-5.5, 1/128);
  f := 1.843974509298089515;
  testrel(12, NE, y, f, cnt,failed);

  y := poch1x(-5.5,-1/256);
  f := 1.767826803772617745;
  testrel(13, NE, y, f, cnt,failed);

  x := ldexp(1,-30);
  y := poch1x(-5.5, x);
  f := 1.792911336415276031;
  testrel(14, NE, y, f, cnt,failed);

  y := poch1x(-5.5, 0);
  f := 1.792911330399932942;
  testrel(15, NE, y, f, cnt,failed);

  y := poch1x(1e4, 1/128);
  f := 9.549748713641611024;
  testrel(16, NE, y, f, cnt,failed);

  y := poch1x(1e8, 1/128);
  f := 19.81209403452184574;
  testrel(17, NE, y, f, cnt,failed);

  y := poch1x(-10000.5, 1/128);
  f := 9.591349003185080560;
  testrel(18, NE, y, f, cnt,failed);

  y := poch1x(1e-10, 1e-20);
  f := -0.9999999999577215665e10;
  testrel(19, NE, y, f, cnt,failed);

  y := poch1x(-1e-10, 1e-20);
  f := 0.10000000000422784335e11;
  testrel(20, NE, y, f, cnt,failed);

  y := poch1x(0, 0.5);
  f := -2;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_binomialx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 10;
  NE2 = 1000*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','binomialx');

  y := binomialx(0,0);
  f := 1;
  testrel( 1, NE, y, f, cnt,failed);

  y := binomialx(-2,-4);
  f := 3;
  testrel( 2, NE, y, f, cnt,failed);

  y := binomialx(-2,-3);
  f := -2;
  testrel( 3, NE, y, f, cnt,failed);

  y := binomialx(-2,-2);
  f := 1;
  testrel( 4, NE, y, f, cnt,failed);

  y := binomialx(-2,-1);
  f := 0;
  testrel( 5, NE, y, f, cnt,failed);

  y := binomialx(-2,0);
  f := 1;
  testrel( 6, NE, y, f, cnt,failed);

  y := binomialx(-2,1);
  f := -2;
  testrel( 7, NE, y, f, cnt,failed);

  y := binomialx(-2,2);
  f := 3;
  testrel( 8, NE, y, f, cnt,failed);

  y := binomialx(-2,3);
  f := -4;
  testrel( 9, NE, y, f, cnt,failed);

  y := binomialx(-3,-12);
  f := -55;
  testrel(10, NE, y, f, cnt,failed);

  y := binomialx(-20,-15);
  f := 0;
  testrel(11, NE, y, f, cnt,failed);

  y := binomialx(-4,-2000);
  f := 1329336999;
  testrel(12, NE, y, f, cnt,failed);

  y := binomialx(-200,-300);
  f := 0.27721676421723764965e82;
  testrel(13, NE, y, f, cnt,failed);

  y := binomialx(1700,6);
  f := 33229582324814800.0;
  testrel(14, NE, y, f, cnt,failed);

  y := binomialx(-1700,6);
  f := 33821192478773700.0;
  testrel(15, NE, y, f, cnt,failed);

  y := binomialx(1800,6);
  f := 46846777478732700.0;
  testrel(16, NE, y, f, cnt,failed);

  y := binomialx(10000,100);
  f := 0.65208469245472575695e242;
  testrel(17, NE, y, f, cnt,failed);

  y := binomialx(4170,4000);
  f := 0.11106511568601554897e308;
  testrel(18, NE, y, f, cnt,failed);

  y := binomialx(1307,1000);
  f := 0.6855892205560626147e308;
  testrel(19, NE, y, f, cnt,failed);

  y := binomialx(1756,1509);
  f := 0.1528133755387334589e309;
  testrel(20, NE, y, f, cnt,failed);

  y := binomialx(4171,4000);
  f := 0.2709079517698075174e309;
  testrel(21, NE, y, f, cnt,failed);

  y := binomialx(-2001,-4000);
  f := -0.83144893751691753477e1202;
  testrel(22, NE, y, f, cnt,failed);

  y := binomialx(-5000,-10000);
  f := 0.79589513176621947417e3008;
  testrel(23, NE2, y, f, cnt,failed);

  y := binomialx(30000,3000);
  f := 0.21750093296621173746e4234;
  testrel(24, NE2, y, f, cnt,failed);

  y := binomialx(16955,4972);
  f := 0.92955668247937447804e4453;
  testrel(25, NE2, y, f, cnt,failed);

  y := binomialx(22012,4652);         {t=11347.050411}
  f := 0.99214609960720754934e4928;
  testrel(26, NE2, y, f, cnt,failed);

  {d1}
  y := binomialx(16865,6772);         {0.25 .. 1 MinExt}
  f := 0.53696452467021885443e4932;
  testrel(27, NE2, y, f, cnt,failed);

  y := binomialx(30000,3777);
  f := 0.30697830915154703045e4930;
  testrel(28, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lnbinomialx;
var
  y,f: extended;
  cnt, failed: integer;
  n: longint;
const
  NE = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lnbinomialx');

  n := 1073741824;
  y := lnbinomialx(n, n div 2);
  f := 744261107.3318939566;
  testrel( 1, NE, y, f, cnt,failed);

  n := MaxLongint-1;
  y := lnbinomialx(n, n div 2);
  f := 1488522223.553919024;
  testrel( 2, NE, y, f, cnt,failed);

  y := lnbinomialx(1600, 800);
  f := 1105.120661839164004;
  testrel( 3, NE, y, f, cnt,failed);

  y := lnbinomialx(1600, 1600);
  f := 0;
  testrel( 4, NE, y, f, cnt,failed);

  y := lnbinomialx(1600, 0);
  f := 0;
  testrel( 5, NE, y, f, cnt,failed);

  y := lnbinomialx(2000, 0);
  f := 0;
  testrel( 6, NE, y, f, cnt,failed);

  y := lnbinomialx(0, 0);
  f := 0;
  testrel( 7, NE, y, f, cnt,failed);

  y := lnbinomialx(1600, 1599);
  f := 7.377758908227872606;
  testrel( 8, NE, y, f, cnt,failed);

  y := lnbinomialx(1600, 1);
  f := 7.377758908227872606;
  testrel( 9, NE, y, f, cnt,failed);

  y := lnbinomialx(170258529,62455558);
  f := 0.1119016866895984359e9;
  testrel(10, NE, y, f, cnt,failed);

  y := lnbinomialx(74434281,13439222);
  f := 35149997.81626606341;
  testrel(11, NE, y, f, cnt,failed);

  y := lnbinomialx(12345678,3456789);
  f := 7320401.650278493808;
  testrel(12, NE, y, f, cnt,failed);

  y := lnbinomialx(123,45);
  f := 78.17880048930058032;
  testrel(13, NE, y, f, cnt,failed);

  y := lnbinomialx(150, 100);
{$ifdef BIT16}
  f := 185.6059266841743233*0.5;
{$else}
  f := 92.80296334208716165;
{$endif}
  testrel(14, NE, y, f, cnt,failed);

  y := lnbinomialx(170, 85);
{$ifdef BIT16}
  f := 230.0797190885321981*0.5;
{$else}
  f := 115.03985954426609906;
{$endif}
  testrel(15, NE, y, f, cnt,failed);

  y := lnbinomialx(9300, 3011);
  f := 5851.237546930416409;
  testrel(16, NE, y, f, cnt,failed);

  y := lnbinomialx(6538, 4913);
  f := 3661.600748327122391;
  testrel(17, NE, y, f, cnt,failed);

  y := lnbinomialx(6220, 3844);
  f := 4131.931897687751260;
  testrel(18, NE, y, f, cnt,failed);

  y := lnbinomialx(99, 57);
  f := 64.96647556283297552;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lnBarnesGx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 6;
  NE2 = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lnBarnesGx');

  y := lnBarnesGx(1);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := lnBarnesGx(2);
  f := 0;
  testrel(2, NE, y, f, cnt,failed);

  y := lnBarnesGx(3);
  f := 0;
  testrel(3, NE, y, f, cnt,failed);

  y := lnBarnesGx(4);
  f := 0.6931471805599453094;
  testrel(4, NE, y, f, cnt,failed);

  y := lnBarnesGx(5);
  f := 2.484906649788000310;
  testrel(5, NE, y, f, cnt,failed);

  y := lnBarnesGx(8);
  f := 17.02970343492809292;
  testrel(6, NE, y, f, cnt,failed);

  y := lnBarnesGx(9);
  f := 25.55486479599350722;
  testrel(7, NE, y, f, cnt,failed);

  y := lnBarnesGx(10);
  f := 36.15946769873875745;
  testrel(8, NE, y, f, cnt,failed);

  y := lnBarnesGx(1/1024);
  f := -6.930500534322783107;
  testrel(9, NE, y, f, cnt,failed);

  y := lnBarnesGx(-1/1024);
  f := -6.932446149758604683;
  testrel(10, NE, y, f, cnt,failed);

  y := lnBarnesGx(1e-4);
  f := -9.210240772666171834;
  testrel(11, NE, y, f, cnt,failed);

  y := lnBarnesGx(1e-30);
  f := -69.07755278982137052;
  testrel(12, NE, y, f, cnt,failed);

  y := lnBarnesGx(0.5);
  f := -0.5054330544896953828;
  testrel(13, NE, y, f, cnt,failed);

  y := lnBarnesGx(-0.5);
  f := -1.770945177974340779;
  testrel(14, NE, y, f, cnt,failed);

  y := lnBarnesGx(-1.5);
  f := -2.630992193350821794;
  testrel(15, NE, y, f, cnt,failed);

  y := lnBarnesGx(-2.5);
  f := -2.574748476853147743;
  testrel(16, NE, y, f, cnt,failed);

  y := lnBarnesGx(-1.23);
  f := -3.242224688030460964;
  testrel(18, NE, y, f, cnt,failed);

  y := lnBarnesGx(1.23);
  f := 0.6058671879802079415e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := lnBarnesGx(-1.25);
  f := -3.120035872301288668;
  testrel(19, NE, y, f, cnt,failed);

  y := lnBarnesGx(1.25);
  f := 0.6301661850380737394e-1;
  testrel(20, NE1, y, f, cnt,failed);

  y := lnBarnesGx(3.5);
  f := 0.2308325212726786416;
  testrel(21, NE, y, f, cnt,failed);

  y := lnBarnesGx(4.5);
  f := 1.431806123619752866;
  testrel(22, NE1, y, f, cnt,failed);

  y := lnBarnesGx(Pi);
  f := 0.4482028387476615430e-1;
  testrel(23, NE1, y, f, cnt,failed);

  y := lnBarnesGx(-Pi);
  f := -5.228834647463673243;
  testrel(24, NE2, y, f, cnt,failed);

  y := lnBarnesGx(10.1);
  f := 37.33894856006103116;
  testrel(25, NE, y, f, cnt,failed);

  y := lnBarnesGx(-12.34);
  f := 92.00323773138391428;
  testrel(26, NE, y, f, cnt,failed);

  y := lnBarnesGx(20);
  f := 277.7702652773851256;
  testrel(27, NE, y, f, cnt,failed);

  y := lnBarnesGx(28);
  f := 678.9537981255094512;
  testrel(28, NE, y, f, cnt,failed);

  y := lnBarnesGx(29);
  f := 743.5113367525157823;
  testrel(29, NE, y, f, cnt,failed);

  y := lnBarnesGx(30);
  f := 811.4010798896973173;
  testrel(30, NE, y, f, cnt,failed);

  y := lnBarnesGx(50);
  f := 2915.918514771743880;
  testrel(31, NE, y, f, cnt,failed);

  y := lnBarnesGx(-200.5);
  f := 77215.70142728548527;
  testrel(32, NE, y, f, cnt,failed);

  y := lnBarnesGx(-200.1);
  f := 76633.15311963577404;
  testrel(33, NE1, y, f, cnt,failed);

  y := lnBarnesGx(-200.9);
  f := 77325.83198224738520;
  testrel(34, NE, y, f, cnt,failed);

  y := lnBarnesGx(-20000.5);
  f := 1680960109.118795368;
  testrel(35, NE, y, f, cnt,failed);

  y := lnBarnesGx(-1000000.25);
  f := 6157770726014.298561;
  testrel(36, NE, y, f, cnt,failed);

  y := lnBarnesGx(1e7);
  f := 7.309046405563504434e14;
  testrel(37, NE, y, f, cnt,failed);

  y := lnBarnesGx(1e5*pi);
  f := 5.506048481182415331e11;
  testrel(38, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



end.
