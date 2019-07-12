{Part 6a of regression test for SPECFUNX unit  (c) 2011  W.Ehrhardt}

unit t_sfx6a;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_sph_bess_jnx;
procedure test_sph_bess_ynx;
procedure test_sph_bessel_inx;
procedure test_sph_bessel_inex;
procedure test_sph_bessel_knx;
procedure test_sph_bessel_knex;

procedure test_airy_aix;
procedure test_airy_bix;
procedure test_airy_aipx;
procedure test_airy_bipx;
procedure test_airy_aisx;
procedure test_airy_bisx;
procedure test_airy_scorerx;

procedure test_berx;
procedure test_beix;
procedure test_kerx;
procedure test_keix;
procedure test_kelvin_derx;

procedure test_struve_h0x;
procedure test_struve_h1x;
procedure test_struve_l0x;
procedure test_struve_l1x;
procedure test_struve_hx;
procedure test_struve_lx;


implementation


uses
  amath, specfunx, t_sfx0;


{---------------------------------------------------------------------------}
procedure test_airy_aix;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 8*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','airy_aix');


  y := airy_aix(50);
  f := 0.45849417240748284783e-103;
  testrel(1, NE2, y, f, cnt,failed);

  y := airy_aix(5);
  f := 0.10834442813607441735e-3;
  testrel(2, NE, y, f, cnt,failed);

  y := airy_aix(1);
  f := 0.13529241631288141552;
  testrel(3, NE, y, f, cnt,failed);

  y := airy_aix(0.5);
  f := 0.23169360648083348977;
  testrel(4, NE, y, f, cnt,failed);

  y := airy_aix(0);
  f := 0.35502805388781723926;
  testrel(5, NE, y, f, cnt,failed);

  y := airy_aix(-0.25);
  f := 0.41872461427545292423;
  testrel(6, NE, y, f, cnt,failed);

  y := airy_aix(-0.5);
  f := 0.47572809161053958880;
  testrel(7, NE, y, f, cnt,failed);

  y := airy_aix(-1.0);
  f := 0.53556088329235211880;
  testrel(8, NE, y, f, cnt,failed);

  y := airy_aix(-2.0);
  f := 0.22740742820168557599;
  testrel(9, NE, y, f, cnt,failed);

  y := airy_aix(-5);
  f := 0.35076100902411431979;
  testrel(10, NE, y, f, cnt,failed);

  y := airy_aix(-20);
  f := -0.17640612707798468959;
  testrel(11, NE2, y, f, cnt,failed);

  y := airy_aix(-100);
  f := 0.17675339323955287809;
  testrel(12, NE2, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_airy_bix;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 8*NE;

const
  BiM1: THexExtW = ($0B56,$11EC,$9556,$D4FC,$3FFB);  {+1.0399738949694461189E-1}

begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','airy_bix');

  y := airy_bix(50);
  f := 0.4909099699444219329e102;
  testrel(1, NE2, y, f, cnt,failed);

  y := airy_bix(5);
  f := 657.7920441711711824;
  testrel(2, NE, y, f, cnt,failed);

  y := airy_bix(1);
  f := 1.207423594952871259;
  testrel(3, NE, y, f, cnt,failed);

  y := airy_bix(0.5);
  f := 0.8542770431031554933;
  testrel(4, NE, y, f, cnt,failed);

  y := airy_bix(0);
  f := 0.6149266274460007352;
  testrel(5, NE, y, f, cnt,failed);

  y := airy_bix(-0.25);
  f := 0.5013998734692333890;
  testrel(6, NE, y, f, cnt,failed);

  y := airy_bix(-0.5);
  f := 0.3803526597510538502;
  testrel(7, NE, y, f, cnt,failed);

  y := airy_bix(-1.0);
  {f := 0.10399738949694461189;}
  f := extended(BiM1);
  testrel(8, NE, y, f, cnt,failed);

  y := airy_bix(-2.0);
  f := -0.4123025879563984881;
  testrel(9, NE, y, f, cnt,failed);

  y := airy_bix(-5);
  f := -0.1383691349016005769;
  testrel(10, NE, y, f, cnt,failed);

  y := airy_bix(-20);
  f := -0.2001393093226513493;
  testrel(11, NE2, y, f, cnt,failed);

  y := airy_bix(-100);
  f := 0.2427388768016013161e-1;
  testrel(12,2100, y, f, cnt,failed);   {!!!!!!!}


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_airy_aipx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 8*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','airy_aipx');


  y := airy_aipx(50);
  f := -0.3244331819828799296e-102;
  testrel(1, NE2, y, f, cnt,failed);

  y := airy_aipx(5);
  f := -0.2474138908684624760e-3;
  testrel(2, NE, y, f, cnt,failed);

  y := airy_aipx(1);
  f := -0.1591474412967932128;
  testrel(3, NE, y, f, cnt,failed);

  y := airy_aipx(0.5);
  f := -0.2249105326646838931;
  testrel(4, NE, y, f, cnt,failed);

  y := airy_aipx(0);
  f := -0.2588194037928067984;
  testrel(5, NE, y, f, cnt,failed);

  y := airy_aipx(-0.25);
  f := -0.2463891899201759730;
  testrel(6, NE, y, f, cnt,failed);

  y := airy_aipx(-0.5);
  f := -0.2040816703395473861;
  testrel(7, NE, y, f, cnt,failed);

  y := airy_aipx(-1.0);
  f := -0.1016056711664520939e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := airy_aipx(-2.0);
  f := 0.6182590207416910414;
  testrel(9, NE, y, f, cnt,failed);

  y := airy_aipx(-5);
  f := 0.3271928185544431368;
  testrel(10, NE2, y, f, cnt,failed);

  y := airy_aipx(-20);
  f := 0.8928628567364712384;
  testrel(11, NE2, y, f, cnt,failed);

  y := airy_aipx(-100);
  f := -0.2422970316605838054;
  testrel(12, 600, y, f, cnt,failed);   {!!!!!!!}


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_airy_bipx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 8*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','airy_bipx');


  y := airy_bipx(50);
  f := 0.3468798779545976724e103;
  testrel(1, NE2, y, f, cnt,failed);

  y := airy_bipx(5);
  f := 1435.819080217982519;
  testrel(2, NE, y, f, cnt,failed);

  y := airy_bipx(1);
  f := 0.9324359333927756330;
  testrel(3, NE, y, f, cnt,failed);

  y := airy_bipx(0.5);
  f := 0.5445725641405923018;
  testrel(4, NE, y, f, cnt,failed);

  y := airy_bipx(0);
  f := 0.4482883573538263579;
  testrel(5, NE, y, f, cnt,failed);

  y := airy_bipx(-0.25);
  f := 0.46515148833715370327;
  testrel(6, NE, y, f, cnt,failed);

  y := airy_bipx(-0.5);
  f := 0.50593371362384716657;
  testrel(7, NE, y, f, cnt,failed);

  y := airy_bipx(-1.0);
  f := 0.59237562642279235082;
  testrel(8, NE, y, f, cnt,failed);

  y := airy_bipx(-2.0);
  f := 0.27879516692116952269;
  testrel(9, NE, y, f, cnt,failed);

  y := airy_bipx(-5);
  f := 0.77841177300189924609;
  testrel(10, NE, y, f, cnt,failed);

  y := airy_bipx(-20);
  f := -0.79142903383953647936;
  testrel(11, NE2, y, f, cnt,failed);

  y := airy_bipx(-100);
  f := 1.7675948932340609324;
  testrel(12, NE2, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_airy_aisx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NE1 = 12;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','airy_aisx');

  y := airy_aisx(1e-20);
  f := 0.3550280538878172393;
  testrel(1, NE, y, f, cnt,failed);

  y := airy_aisx(0.125);
  f := 0.3324375951095047298;
  testrel(2, NE, y, f, cnt,failed);

  y := airy_aisx(1);
  f := 0.2635136447491400686;
  testrel(3, NE, y, f, cnt,failed);

  y := airy_aisx(1.125);
  f := 0.2579914489991170531;
  testrel(4, NE, y, f, cnt,failed);

  y := airy_aisx(1.5);
  f := 0.2441848976714084334;
  testrel(5, NE, y, f, cnt,failed);

  y := airy_aisx(1.9375);
  f := 0.2317137238870651524;
  testrel(6, NE1, y, f, cnt,failed);

  y := airy_aisx(2);
  f := 0.2301649186525116059;
  testrel(7, NE1, y, f, cnt,failed);

  y := airy_aisx(3);
  f := 0.2105720427859769851;
  testrel(8, NE, y, f, cnt,failed);

  y := airy_aisx(4);
  f := 0.1970948026430665127;
  testrel(9, NE, y, f, cnt,failed);

  y := airy_aisx(10);
  f := 0.1581236668543461503;
  testrel(10, NE, y, f, cnt,failed);

  y := airy_aisx(100);
  f := 0.8919692093633041318e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := airy_aisx(500);
  f := 0.05965522950749524825;
  testrel(12, NE, y, f, cnt,failed);

  y := airy_aisx(1000);
  f := 0.50164170749970862190e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := airy_aisx(17179869184.0);
  f := 0.7791841414090481631e-3;
  testrel(14, NE, y, f, cnt,failed);

  y := airy_aisx(2e12);
  f := 0.2372124991643971727e-3;
  testrel(15, NE, y, f, cnt,failed);

  y := airy_aisx(3e12);
  f := 0.2143456895262479282e-3;
  testrel(16, NE, y, f, cnt,failed);

  y := airy_aisx(maxdouble);
  f := 0.2436218170273481498e-77;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_airy_bisx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','airy_bisx');


  y := airy_bisx(1e-20);
  f := 0.6149266274460007352;
  testrel(1, NE, y, f, cnt,failed);

  y := airy_bisx(0.125);
  f := 0.6516858507482041529;
  testrel(2, NE, y, f, cnt,failed);

  y := airy_bisx(1);
  f := 0.6199119435726784851;
  testrel(3, NE, y, f, cnt,failed);

  y := airy_bisx(2);
  f := 0.5004372543040949650;
  testrel(4, NE, y, f, cnt,failed);

  y := airy_bisx(2.125);
  f := 0.4902398661293438879;
  testrel(5, NE, y, f, cnt,failed);

  y := airy_bisx(4);
  f := 0.4048094678892980736;
  testrel(6, NE, y, f, cnt,failed);

  y := airy_bisx(10);
  f := 0.3183401053367344453;
  testrel(7, NE, y, f, cnt,failed);

  y := airy_bisx(12);
  f := 0.3039054138807329093;
  testrel(8, NE, y, f, cnt,failed);

  y := airy_bisx(100);
  f := 0.1784310111708354151;
  testrel(9, NE, y, f, cnt,failed);

  y := airy_bisx(500);
  f := 0.1193126822548645875;
  testrel(10, NE, y, f, cnt,failed);

  y := airy_bisx(1000);
  f := 0.1003290024731051856;
  testrel(11, NE, y, f, cnt,failed);

  y := airy_bisx(17179869184.0);
  f := 0.1558368282818096470e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := airy_bisx(2e12);
  f := 0.4744249983287943454e-3;
  testrel(13, NE, y, f, cnt,failed);

  y := airy_bisx(3e12);
  f := 0.4286913790524958564e-3;
  testrel(14, NE, y, f, cnt,failed);

  y := airy_bisx(maxdouble);
  f := 0.4872436340546962995e-77;
  testrel(15, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_airy_scorerx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 8;
  NE1 = 10;
  NE2 = 320;
const
  gm8: THexExtW = ($9F2F,$C445,$1E8A,$8180,$BFFD); {-0.252930597724240196942532167187}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','airy_gix/hix');

  {MISCFUN test values}
  f := 0.20468308070040542435;
  y := airy_gix(-0.001953125);
  testrel(1, NE, y, f, cnt,failed);

  f := 0.18374662832557904078;
  y := airy_gix(-0.125);
  testrel(2, NE, y, f, cnt,failed);

  f :=-0.11667221729601528265;
  y := airy_gix(-1);
  testrel(3, NE, y, f, cnt,failed);

  f := 0.31466934902729557596;
  y := airy_gix(-4);
  testrel(4, NE, y, f, cnt,failed);

  f :=-0.37089040722426257729;
  y := airy_gix(-8);
  testrel(5, NE, y, f, cnt,failed);

  {f := -0.25293059772424019694;}
  f := extended(gm8);
  y := airy_gix(-8.25);
  testrel(6, NE1, y, f, cnt,failed);

  f := 0.28967410658692701936;
  y := airy_gix(-9);
  testrel(7, NE, y, f, cnt,failed);

  f :=-0.34644836492634090590;
  y := airy_gix(-10);
  testrel(8, NE, y, f, cnt,failed);

  f := 0.28076035913873049496;
  y := airy_gix(-11);
  testrel(9, NE, y, f, cnt,failed);

  f := 0.21814994508094865815;
  y := airy_gix(-13);
  testrel(10, NE, y, f, cnt,failed);

  f := 0.20526679000810503329;
  y := airy_gix(0.001953125);
  testrel(11, NE, y, f, cnt,failed);

  f := 0.22123695363784773258;
  y := airy_gix(0.125);
  testrel(12, NE, y, f, cnt,failed);

  f := 0.23521843981043793760;
  y := airy_gix(1);
  testrel(13, NE, y, f, cnt,failed);

  f := 0.82834303363768729338e-1;
  y := airy_gix(4.0);
  testrel(15, NE, y, f, cnt,failed);

  y := airy_gix(7.0);
  f := 0.45757385490989281893e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := airy_gix(7.25);
  f := 0.44150012014605159922e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := airy_gix(8.0);
  f := 0.39951133719508907541e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := airy_gix(9.0);
  f := 0.35467706833949671483e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := airy_gix(10.0);
  f := 0.31896005100679587981e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := airy_gix(12.0);
  f := 0.26556892713512410405e-1;
  testrel(20, NE, y, f, cnt,failed);

  f := 0.40936798278458884024;
  y := airy_hix(-0.001953125);
  testrel(21, NE, y, f, cnt,failed);

  f := 0.37495291608048868619;
  y := airy_hix(-0.125);
  testrel(22, NE, y, f, cnt,failed);

  f := 0.22066960679295989454;
  y := airy_hix(-1);
  testrel(23, NE, y, f, cnt,failed);

  f := 0.77565356679703713590e-1;
  y := airy_hix(-4);
  testrel(24, NE, y, f, cnt,failed);

  f := 0.39638826473124717315e-1;
  y := airy_hix(-8);
  testrel(25, NE, y, f, cnt,failed);

  f := 0.38450072575004151871e-1;
  y := airy_hix(-8.25);
  testrel(26, NE, y, f, cnt,failed);

  f := 0.35273216868317898556e-1;
  y := airy_hix(-9);
  testrel(27, NE, y, f, cnt,failed);

  f := 0.31768535282502272742e-1;
  y := airy_hix(-10);
  testrel(28, NE, y, f, cnt,failed);

  f := 0.28894408288051391369e-1;
  y := airy_hix(-11);
  testrel(29, NE, y, f, cnt,failed);

  f := 0.24463284011678541180e-1;
  y := airy_hix(-13);
  testrel(30, NE, y, f, cnt,failed);

  f := 0.41053540139998941517;
  y := airy_hix(0.001953125);
  testrel(31, NE, y, f, cnt,failed);

  f := 0.44993502381204990817;
  y := airy_hix(0.125);
  testrel(32, NE, y, f, cnt,failed);

  f := 0.97220515514243332184;
  y := airy_hix(1);
  testrel(33, NE1, y, f, cnt,failed);

  f := 0.83764237105104371193e2;
  y := airy_hix(4);
  testrel(34, NE, y, f, cnt,failed);

  f := 0.80327744952044756016e5;
  y := airy_hix(7);
  testrel(35, NE, y, f, cnt,failed);

  f := 0.15514138847749108298e6;
  y := airy_hix(7.25);
  testrel(36, NE, y, f, cnt,failed);

  f := 0.11995859641733262114e7;
  y := airy_hix(8);
  testrel(37, NE, y, f, cnt,failed);

  y := airy_hix(9.0);
  f := 0.21472868855967642259e8;
  testrel(38, NE, y, f, cnt,failed);

  y := airy_hix(10.0);
  f := 0.45564115351632913590e9;
  testrel(39, NE, y, f, cnt,failed);

  y := airy_hix(12.0);
  f := 0.32980722582904761929e12;
  testrel(40, NE1, y, f, cnt,failed);

  {---------------------------------------}
  f := 0.392203077804138178e36;
  y := airy_hix(25);
  testrel(41, NE, y, f, cnt,failed);

  f := 0.6041223996670201399e289;
  y := airy_hix(100);
  testrel(42, NE2, y, f, cnt,failed);

  f := 0.3183098855471709119e-3;
  y := airy_hix(-1000);
  testrel(43, NE, y, f, cnt,failed);

  f := 0.3183098861831540518e-4;
  y := airy_hix(-10000);
  testrel(44, NE, y, f, cnt,failed);

  f := 2.122065907891937809e-7; {Alpha}
  y := airy_hix(-1.5e6);
  testrel(45, NE, y, f, cnt,failed);

  f := 1.061032953945968905e-7; {Alpha}
  y := airy_hix(-3e6);
  testrel(46, NE, y, f, cnt,failed);

  f := 3.183098861837906715e-8; {Alpha}
  y := airy_hix(-1e7);
  testrel(47, NE, y, f, cnt,failed);

  y := airy_gix(20.0);
  f := 0.1591948320056103765e-1;
  testrel(48, NE, y, f, cnt,failed);

  y := airy_gix(50.0);
  f := 0.6366299599144166116e-2;
  testrel(49, NE, y, f, cnt,failed);

  y := airy_gix(100.0);
  f := 0.3183105228162961477e-2;
  testrel(50, NE, y, f, cnt,failed);

  f := -0.1435162480096224794;
  y := airy_gix(-50);
  testrel(51, NE, y, f, cnt,failed);

  f := -0.8358288400262780392e-1;
  y := airy_gix(-1000);
  testrel(52, NE2, y, f, cnt,failed);

  f := -0.4953937439675591109e-1;
  y := airy_gix(-10000);
  testrele(53, 2e-14, y, f, cnt,failed);      {!!!!!!!!!!!!!}

  f := 0.8443789579269917059e-2;
  y := airy_gix(-1e7);
  testrele(54, 2e-10, y, f, cnt,failed);      {!!!!!!!!!!!!!}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_sph_bess_jnx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 8;
  NE2 = 4*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sph_bessel_jnx');

  {Test values from GSL, some corrected with Maple}
  {jn := (n,x) -> sqrt(Pi/2/x)*BesselJ(n+1/2,x); }

  y := sph_bessel_jnx(0,-10.0);
  f := -0.05440211108893698134;
  testrel(1, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(0,0.001);
  f := 0.9999998333333416667;
  testrel(2, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(0,1.0);
  f := 0.84147098480789650670;
  testrel(3, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(0,10.0);
  f := -0.05440211108893698134;
  testrel(4, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(0,100.0);
  f := -0.005063656411097587937;
  testrel(5, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(0,1048576.0);
  f := 3.1518281938718287624e-07;
  testrel(6, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(1,-10.0);
  f := -0.07846694179875154709;
  testrel(7, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(1,0.01);
  f := 0.003333300000119047399;
  testrel(8, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(1,1.0);
  f := 0.30116867893975678925;
  testrel(9, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(1,10.0);
  f := 0.07846694179875154709;
  testrel(10, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(1,100.0);
  f := -0.008673825286987815220;
  testrel(11, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(1,1048576.0);
  f := -9.000855242905546158e-07;
  testrel(12, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(2,-10.0);
  f := 0.07794219362856244547;
  testrel(13, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(2,0.01);
  f := 6.666619047751322551e-06;
  testrel(14, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(2,1.0);
  f := 0.062035052011373861103;
  testrel(15, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(2,10.0);
  f := 0.077942193628562445467;
  testrel(16, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(2,100.0);
  f := 0.48034416524879534799e-2;
  testrel(17, NE2, y, f, cnt,failed);

  y := sph_bessel_jnx(2,1048576.0);
  f := -3.1518539455252413111e-07;
  testrel(18, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(0,0.0);
  f := 1.0;
  testrel(19, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(1,10.0);
  f := 0.07846694179875154709000;
  testrel(20, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(5,1.0);
  f := 0.00009256115861125816357;
  testrel(21, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(10,10.0);
  f := 0.06460515449256426427;
  testrel(22, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(100,100.0);
  f := 0.010880477011438336540;
  testrel(23, NE2, y, f, cnt,failed);

  y := sph_bessel_jnx(2,900.0);
  f := -0.0011089115568832940086;
  testrel(24, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(2,15000.0);
  f := -0.59555920330757505537e-4;
  testrel(25, NE, y, f, cnt,failed);

  y := sph_bessel_jnx(100,1000.0);
  f := -0.00025326311230945818285;
  testrel(26, NE2, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sph_bess_ynx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 8;
  NE2 = 4*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sph_bessel_ynx');

  {Test values from GSL, some corrected with Maple}
  {yn := (n,x) -> sqrt(Pi/2/x)*BesselY(n+1/2,x);}

  y := sph_bessel_ynx(0,0.001);
  f := -999.99950000004166670;
  testrel(1, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(0,1.0);
  f := -0.5403023058681397174;
  testrel(2, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(0,10.0);
  f := 0.08390715290764524523;
  testrel(3, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(0,100.0);
  f := -0.008623188722876839341;
  testrel(4, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(0,65536.0);
  f := 0.000011014324202158573930;
  testrel(5, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(0,4294967296.0);
  f := 2.0649445131370357007e-10;
  testrel(6, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(1,0.01);
  f := -10000.499987500069444;
  testrel(7, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(1,1.0);
  f := -1.3817732906760362241;
  testrel(8, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(1,10.0);
  f := 0.06279282637970150586;
  testrel(9, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(1,100.0);
  f := 0.49774245238688195431e-2;
  testrel(10, NE2, y, f, cnt,failed);

  y := sph_bessel_ynx(1,4294967296.0);
  f := 1.0756463271573404688e-10;
  testrel(11, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(2,0.01);
  f := -3.0000500012499791668e+06;
  testrel(12, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(2,1.0);
  f := -3.605017566159968955;
  testrel(13, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(2,10.0);
  f := -0.06506930499373479347;
  testrel(4, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(2,100.0);
  f := 0.008772511458592903927;
  testrel(15, NE2, y, f, cnt,failed);

  y := sph_bessel_ynx(2,4294967296.0);
  f := -2.0649445123857054207e-10;
  testrel(16, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(0,0.0078125);
  f := -127.9960937698681745;
  testrel(17, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(5,1.0);
  f := -999.44034339223640949;
  testrel(20, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(10,0.01);
  f := -0.65473079797378378406e31;
  testrel(21, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(10,10.0);
  f := -0.172453672088057849;
  testrel(22, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(100,1.0);
  f := -0.66830794632586775138e187;
  testrel(23, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(100,100.0);
  f := -0.22983850491562281089e-1;
  testrel(24, NE, y, f, cnt,failed);

  y := sph_bessel_ynx(2000,1048576.0);
  f := 5.9545201447146155e-07;
  testrel(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_sph_bessel_inex;
var
  f,y: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sph_bessel_inex');

  {GSL arguments, recalculated with Maple}
  y := sph_bessel_inex(0,0);
  f := 1.0;
  testrel(1, NE, y, f, cnt,failed);

  y := sph_bessel_inex(1,0);
  f := 0.0;
  testrel(2, NE, y, f, cnt,failed);

  y := sph_bessel_inex(2,0);
  f := 0.0;
  testrel(3, NE, y, f, cnt,failed);

  y := sph_bessel_inex(0,0.1);
  f := 0.9063462346100907067;
  testrel(4, NE, y, f, cnt,failed);

  y := sph_bessel_inex(1,0.1);
  f := 0.3019141928900222685e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := sph_bessel_inex(2,0.1);
  f := 0.6036559400239012567e-3;
  testrel(6, NE, y, f, cnt,failed);

  y := sph_bessel_inex(0,2);
  f := 0.2454210902778164549;
  testrel(7, NE, y, f, cnt,failed);

  y := sph_bessel_inex(1,2);
  f := 0.1318683645832753176;
  testrel(8, NE, y, f, cnt,failed);

  y := sph_bessel_inex(2,2);
  f := 0.4761854340290347851e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := sph_bessel_inex(0,100);
  f := 0.005;
  testrel(10, NE, y, f, cnt,failed);

  y := sph_bessel_inex(1,100);
  f := 0.00495;
  testrel(11, NE, y, f, cnt,failed);

  y := sph_bessel_inex(2,100);
  f := 0.48515e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := sph_bessel_inex(4,0.001);
  f := 0.1057143434119036501e-14;
  testrel(13, NE, y, f, cnt,failed);

  y := sph_bessel_inex(4,0.1);
  f := 0.957935224205713493e-7;
  testrel(14, NE, y, f, cnt,failed);

  y := sph_bessel_inex(5,2);
  f := 0.4851564602127540059e-3;
  testrel(15, NE, y, f, cnt,failed);

  y := sph_bessel_inex(5,100);
  f := 0.43004467775e-2;
  testrel(16, NE, y, f, cnt,failed);

  y := sph_bessel_inex(100,100);
  f := 0.1389816196429913279e-22;
  testrel(17, NE, y, f, cnt,failed);

  {Other arguments, recalculated with Maple}
  y := sph_bessel_inex(0,1e-5);
  f := 0.9999900000666663333;
  testrel(18, NE, y, f, cnt,failed);

  y := sph_bessel_inex(0,20000);
  f := 0.25e-4;
  testrel(19, NE, y, f, cnt,failed);

  y := sph_bessel_inex(5,20000);
  f := 0.2498125656118764765e-4;
  testrel(20, NE, y, f, cnt,failed);

  y := sph_bessel_inex(0,-8);
  f := 0.6249999296655158005e-1;
  testrel(21, NE, y, f, cnt,failed);

  y := sph_bessel_inex(1,-8);
  f := -0.5468750791262947245e-1;
  testrel(22, NE, y, f, cnt,failed);

  y := sph_bessel_inex(2,-8);
  f := 0.4199217749931552788e-1;
  testrel(23, NE, y, f, cnt,failed);

  y := sph_bessel_inex(3,-8);
  f := -0.2844239697555726752e-1;
  testrel(24, NE, y, f, cnt,failed);

  y := sph_bessel_inex(-1,8);
  f := 0.6250000703344841995e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := sph_bessel_inex(-2,8);
  f := 0.5468749208737052755e-1;
  testrel(26, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sph_bessel_inx;
var
  f,y: extended;
  cnt, failed: integer;
const
  NE = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sph_bessel_inx');

  y := sph_bessel_inx(0,0);
  f := 1.0;
  testrel(1, NE, y, f, cnt,failed);

  y := sph_bessel_inx(1,0);
  f := 0.0;
  testrel(2, NE, y, f, cnt,failed);

  y := sph_bessel_inx(2,0);
  f := 0.0;
  testrel(3, NE, y, f, cnt,failed);

  y := sph_bessel_inx(0,0.1);
  f := 1.001667500198440258;
  testrel(4, NE, y, f, cnt,failed);

  y := sph_bessel_inx(1,0.1);
  f := 0.3336667857363340751e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := sph_bessel_inx(2,0.1);
  f := 0.6671429894380330319e-3;
  testrel(6, NE, y, f, cnt,failed);

  y := sph_bessel_inx(0,2);
  f := 1.813430203923509384;
  testrel(7, NE, y, f, cnt,failed);

  y := sph_bessel_inx(1,2);
  f := 0.974382743580061038;
  testrel(8, NE, y, f, cnt,failed);

  y := sph_bessel_inx(2,2);
  f := 0.3518560885534178270;
  testrel(9, NE, y, f, cnt,failed);

  y := sph_bessel_inx(0,100);
  f := 0.1344058570908067724e42;
  testrel(10, NE, y, f, cnt,failed);

  y := sph_bessel_inx(1,100);
  f := 0.1330617985198987047e42;
  testrel(11, NE, y, f, cnt,failed);

  y := sph_bessel_inx(2,100);
  f := 0.1304140031352098113e42;
  testrel(12, NE, y, f, cnt,failed);

  y := sph_bessel_inx(4,0.001);
  f := 0.1058201106301107226e-14;
  testrel(13, NE, y, f, cnt,failed);

  y := sph_bessel_inx(4,0.1);
  f := 0.1058682151192429726e-6;
  testrel(14, NE, y, f, cnt,failed);

  y := sph_bessel_inx(5,2);
  f := 0.3584848301270655334e-2;
  testrel(15, NE, y, f, cnt,failed);

  y := sph_bessel_inx(5,100);
  f := 0.1156010470006571019e42;
  testrel(16, NE, y, f, cnt,failed);

  y := sph_bessel_inx(100,100);
  f := 0.3735988741596951155e21;
  testrel(17, NE, y, f, cnt,failed);

  y := sph_bessel_inx(0,1e-5);
  f := 1.000000000016666667;
  testrel(18, NE, y, f, cnt,failed);

  y := sph_bessel_inx(0,10000);
  f := 0.4403409112831460794e4339;
  testrel(19, NE, y, f, cnt,failed);

  y := sph_bessel_inx(5,10000);
  f := 0.4396808620892766329e4339;
  testrel(20, NE, y, f, cnt,failed);

  y := sph_bessel_inx(0,-8);
  f := 186.3098532236937733;
  testrel(21, NE, y, f, cnt,failed);

  y := sph_bessel_inx(1,-8);
  f := -163.0211635035605394;
  testrel(22, NE, y, f, cnt,failed);

  y := sph_bessel_inx(2,-8);
  f := 125.1769169098585710;
  testrel(23, NE, y, f, cnt,failed);

  y := sph_bessel_inx(3,-8);
  f := -84.78559043489893256;
  testrel(24, NE, y, f, cnt,failed);

  y := sph_bessel_inx(-1,8);
  f := 186.3098951565222611;
  testrel(25, NE, y, f, cnt,failed);

  y := sph_bessel_inx(-2,8);
  f := 163.0211163291284906;
  testrel(26, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_sph_bessel_knex;
var
  f,y: extended;
  cnt, failed: integer;
const
  NE  = 5;
  NE2 = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sph_bessel_knex');

  {GSL arguments, recalculated with Maple}
  y := sph_bessel_knex(0,0.1);
  f := 15.70796326794896619;
  testrel(1, NE, y, f, cnt,failed);

  y := sph_bessel_knex(0,2);
  f := 0.7853981633974483096;
  testrel(2, NE, y, f, cnt,failed);

  y := sph_bessel_knex(0,100);
  f := 0.1570796326794896619e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := sph_bessel_knex(1,0.1);
  f := 172.7875959474386281;
  testrel(4, NE, y, f, cnt,failed);

  y := sph_bessel_knex(1,2);
  f := 1.178097245096172464;
  testrel(5, NE, y, f, cnt,failed);

  y := sph_bessel_knex(1,100);
  f := 0.1586504290062845585e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := sph_bessel_knex(2,0.1);
  f := 5199.335841691107810;
  testrel(7, NE, y, f, cnt,failed);

  y := sph_bessel_knex(2,2);
  f := 2.552544031041707006;
  testrel(8, NE, y, f, cnt,failed);

  y := sph_bessel_knex(2,100);
  f := 0.1618391455496781987e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := sph_bessel_knex(4,1/256);
  f := 0.1820559981696195444e15;
  testrel(10, NE, y, f, cnt,failed);

  y := sph_bessel_knex(4,1/8);
  f := 0.6117321781440659753e7;
  testrel(11, NE, y, f, cnt,failed);

  y := sph_bessel_knex(5,2);
  f := 138.1073582949200512;
  testrel(12, NE, y, f, cnt,failed);

  y := sph_bessel_knex(100,100);
  f := 0.3985930768060258219e19;
  testrel(13, NE, y, f, cnt,failed);

  {Other arguments, calculated with Maple}
  y := sph_bessel_knex(25,1e-10);
  f := 0.9179080510564203585e292;
  testrel(14, NE, y, f, cnt,failed);

  y := sph_bessel_knex(-25,1e-10);
  f := 0.1873281736849837466e281;
  testrel(15, NE, y, f, cnt,failed);

  y := sph_bessel_knex(100,2);
  f := 0.3021343835088123453e158;
  testrel(16, NE2, y, f, cnt,failed);

  y := sph_bessel_knex(100,10);
  f := 0.1794607980003062273e91;
  testrel(17, NE2, y, f, cnt,failed);

  y := sph_bessel_knex(100,1000);
  f := 0.243431420384866238982038265809;
  testrel(18, NE, y, f, cnt,failed);

  y := sph_bessel_knex(-10,100);
  f := 0.2457213721229154669e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := sph_bessel_knex(-21,100);
  f := 0.1260553345815578451;
  testrel(20, NE, y, f, cnt,failed);

  y := sph_bessel_knex(-4,1/8);
  f := 109189.1942681668538;
  testrel(21, NE, y, f, cnt,failed);

  y := sph_bessel_knex(-5,1/8);
  f := 6117321.781440659753;
  testrel(22, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sph_bessel_knx;
var
  f,y: extended;
  cnt, failed: integer;
const
  NE = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sph_bessel_knx');

  y := sph_bessel_knx(0,0.1);
  f := 14.21315292597463638;
  testrel(1, NE, y, f, cnt,failed);

  y := sph_bessel_knx(0,2);
  f := 0.1062920828969090821;
  testrel(2, NE, y, f, cnt,failed);

  y := sph_bessel_knx(0,100);
  f := 0.5843481678531469047e-45;
  testrel(3, NE, y, f, cnt,failed);

  y := sph_bessel_knx(1,0.1);
  f := 156.3446821857210002;
  testrel(4, NE, y, f, cnt,failed);

  y := sph_bessel_knx(1,2);
  f := 0.15943812434536362316;
  testrel(5, NE, y, f, cnt,failed);

  y := sph_bessel_knx(1,100);
  f := 0.5901916495316783737e-45;
  testrel(6, NE, y, f, cnt,failed);

  y := sph_bessel_knx(2,0.1);
  f := 4704.553618497604642;
  testrel(7, NE, y, f, cnt,failed);

  y := sph_bessel_knx(2,2);
  f := 0.3454492694149545169;
  testrel(8, NE, y, f, cnt,failed);

  y := sph_bessel_knx(2,100);
  f := 0.6020539173390972559e-45;
  testrel(9, NE, y, f, cnt,failed);

  y := sph_bessel_knx(4,1/256);
  f := 0.1813462290970072313e15;
  testrel(10, NE, y, f, cnt,failed);

  y := sph_bessel_knx(4,1/8);
  f := 0.5398517524234661520e7;
  testrel(11, NE, y, f, cnt,failed);

  y := sph_bessel_knx(5,2);
  f := 18.69079845190335641;
  testrel(12, NE, y, f, cnt,failed);

  y := sph_bessel_knx(100,100);
  f := 0.1482796529234324543e-24;
  testrel(13, NE, y, f, cnt,failed);

  y := sph_bessel_knx(25,1e-10);
  f := 0.9179080509646295534e292;
  testrel(14, NE, y, f, cnt,failed);

  y := sph_bessel_knx(-25,1e-10);
  f := 0.1873281736662509293e281;
  testrel(15, NE, y, f, cnt,failed);

  y := sph_bessel_knx(100,2);
  f := 0.4088944236768448154e157;
  testrel(16, NE, y, f, cnt,failed);

  y := sph_bessel_knx(100,10);
  f := 0.8147507624333384617e86;
  testrel(17, NE, y, f, cnt,failed);

  y := sph_bessel_knx(100,1000);
  f := 0.1235647884245663991e-434;
  testrel(18, NE, y, f, cnt,failed);

  y := sph_bessel_knx(-10,100);
  f := 0.9141021732293337889e-45;
  testrel(19, NE, y, f, cnt,failed);

  y := sph_bessel_knx(-21,100);
  f := 0.4689354218261218366e-44;
  testrel(20, NE, y, f, cnt,failed);

  y := sph_bessel_knx(-4,1/8);
  f := 96359.12573736490671;
  testrel(21, NE, y, f, cnt,failed);

  y := sph_bessel_knx(-5,1/8);
  f := 5398517.524234661520;
  testrel(22, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_berx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 3*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','kelvin_berx');


  y := kelvin_berx(0);
  f := 1.0;
  testrel(1, NE, y, f, cnt,failed);

  y := kelvin_berx(1e-5);
  f := 1.0;
  testrel(2, NE, y, f, cnt,failed);

  y := kelvin_berx(0.03125);
  f := 0.9999999850988388123;
  testrel(3, NE, y, f, cnt,failed);

  y := kelvin_berx(0.5);
  f := 0.9990234639908382556;
  testrel(4, NE, y, f, cnt,failed);

  y := kelvin_berx(-1);
  f := 0.9843817812130868840;
  testrel(5, NE, y, f, cnt,failed);

  y := kelvin_berx(2.75);
  f := 0.1284781977714383538;
  testrel(6, NE, y, f, cnt,failed);

  y := kelvin_berx(-3.0);
  f := -0.2213802495986938889;
  testrel(7, NE, y, f, cnt,failed);

  y := kelvin_berx(10);
  f := 138.8404659416326472;
  testrel(8, NE, y, f, cnt,failed);

  y := kelvin_berx(19.5);
  f := 59956.93311223375741;
  testrel(9, NE2, y, f, cnt,failed);

  y := kelvin_berx(-20);
  f := 47489.37026506176015;
  testrel(10, NE2, y, f, cnt,failed);

  y := kelvin_berx(100);
  f := 0.7368706878094957313e29;
  testrele(11, 2e-17, y, f, cnt,failed);

  y := kelvin_berx(1000);
  f := -0.1545186630003373009e306;
  testrele(12, 5e-17, y, f, cnt,failed);

  y := kelvin_berx(15000);
  f := 0.782333794271490522655583e4604;
  testrele(13, 4e-16, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_beix;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 2*NE;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','kelvin_beix');


  y := kelvin_beix(0);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := kelvin_beix(1e-5);
  f := 0.25e-10;
  testrel(2, NE, y, f, cnt,failed);

  y := kelvin_beix(0.03125);
  f := 0.2441406245957801326e-3;
  testrel(3, NE, y, f, cnt,failed);

  y := kelvin_beix(0.5);
  f := 0.6249321838219945865e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := kelvin_beix(-1);
  f := 0.2495660400366597214;
  testrel(5, NE, y, f, cnt,failed);

  y := kelvin_beix(2.75);
  f := 1.704577752311522268;
  testrel(6, NE, y, f, cnt,failed);

  y := kelvin_beix(-3.0);
  f := 1.937586785266042767;
  testrel(7, NE, y, f, cnt,failed);

  y := kelvin_beix(10);
  f := 56.37045855390663823;
  testrel(8, NE, y, f, cnt,failed);

  y := kelvin_beix(19.5);
  f := 64879.42349727954849;
  testrel(9, NE2, y, f, cnt,failed);

  y := kelvin_beix(-20);
  f := 114775.1973600662216;
  testrel(10, NE2, y, f, cnt,failed);

  y := kelvin_beix(100);
  f := 0.1906911409362379763e30;
  testrele(11, 8e-18, y, f, cnt,failed);

  y := kelvin_beix(1000);
  f := 0.2246152918745784947e305;
  testrele(12, 4e-16, y, f, cnt,failed);

  y := kelvin_beix(15000);
  f := 0.1522554751444777512e4604;
  testrele(13, 5e-16, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_kerx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 5;
  NE2 = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','kelvin_kerx');


  y := kelvin_kerx(1e-10);
  f := 23.14178244559886929;
  testrel(1, NE, y, f, cnt,failed);

  y := kelvin_kerx(1e-5);
  f := 11.62885698064827582;
  testrel(2, NE, y, f, cnt,failed);

  y := kelvin_kerx(0.03125);
  f := 3.581859090333561926;
  testrel(3, NE, y, f, cnt,failed);

  y := kelvin_kerx(0.5);
  f := 0.8559058721186342137;
  testrel(4, NE, y, f, cnt,failed);

  y := kelvin_kerx(1);
  f := 0.2867062087283160460;
  testrel(5, NE, y, f, cnt,failed);

  y := kelvin_kerx(2.75);
  f := -0.7072857588003313910e-1;
  testrele(6, 3e-18, y, f, cnt,failed);

  y := kelvin_kerx(3.0);
  f := -0.6702923330379869775e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := kelvin_kerx(10);
  f := 0.1294663302148061222e-3;
  testrel(8, NE2, y, f, cnt,failed);

  y := kelvin_kerx(19.5);
  f := -0.1153142853119002640e-7;
  testrele(9, 3e-17, y, f, cnt,failed);

  y := kelvin_kerx(20);
  f := -0.7715233109860961461e-7;
  testrele(10, 3e-17, y, f, cnt,failed);

  y := kelvin_kerx(100);
  f := -0.9898417996730774015e-32;
  testrele(11, 3e-17, y, f, cnt,failed);

  y := kelvin_kerx(1000);
  f := -0.2566470946629444788e-308;
  testrele(12, 8e-17, y, f, cnt,failed);

  y := kelvin_kerx(15000);
  f := 0.2337915022377877664e-4608;
  testrele(13, 2e-15, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_keix;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','kelvin_keix');


  y := kelvin_keix(1e-10);
  f := -0.7853981633974483096;
  testrel(1, NE, y, f, cnt,failed);

  y := kelvin_keix(1e-5);
  f := -0.7853981630817268851;
  testrel(2, NE, y, f, cnt,failed);

  y := kelvin_keix(0.03125);
  f := -0.7842795805492080246;
  testrel(3, NE, y, f, cnt,failed);

  y := kelvin_keix(0.5);
  f := -0.6715816950943676032;
  testrel(4, NE, y, f, cnt,failed);

  y := kelvin_keix(1);
  f := -0.4949946365187199003;
  testrel(5, NE, y, f, cnt,failed);

  y := kelvin_keix(2.75);
  f := -0.7735398824302068709e-1;
  testrele(6, 3e-18 , y, f, cnt,failed);

  y := kelvin_keix(3.0);
  f := -0.5112188404598678140e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := kelvin_keix(10);
  f := -0.3075245690881441990e-3;
  testrel(8, NE, y, f, cnt,failed);

  y := kelvin_keix(19.5);
  f := -0.2900202390640604669e-6;
  testrel(9, NE, y, f, cnt,failed);

  y := kelvin_keix(20);
  f := -0.1858941511119437205e-6;
  testrel(10, NE, y, f, cnt,failed);

  y := kelvin_keix(100);
  f := -0.2236535526041445723e-31;
  testrele(11, 7e-18, y, f, cnt,failed);

  y := kelvin_keix(1000);
  f := 0.1915021570632197478e-308;
  testrele(12, 7e-18, y, f, cnt,failed);

  y := kelvin_keix(15000);
  f := -0.3467807459820270330e-4608;
  testrele(13, 3e-16, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_kelvin_derx;
var
  x, fker, fkei, fber, fbei: extended;
  kerp, keip, berp, beip: extended;
var
  cnt, failed: integer;
const
  NE  = 16;
  NE1 = 128;
  NE2 = 800;
  EPS1 = 5e-14;
  EPS2 = 1e-12;
  EPS3 = 2e-11;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Derivatives of Kelvin functions');

  {----------------------------------}
  x := 0.03125;
  kelvin_derx(x, berp, beip, kerp, keip);

  fber := -0.1907348631233516143e-5;
  testrel(1, NE, berp, fber, cnt,failed);

  fbei := 0.1562499992238978547e-1;
  testrel(2, NE, beip, fbei, cnt,failed);

  fker := -31.98773736943210584;
  testrel(3, NE, kerp, fker, cnt,failed);

  fkei := 0.6377755103419601698e-1;
  testrel(4, NE, keip, fkei, cnt,failed);

  {----------------------------------}
  x := 1.0;
  kelvin_derx(x, berp, beip, kerp, keip);

  fber := -0.6244575217903096024e-1;
  testrel(5, NE, berp, fber, cnt,failed);

  fbei := 0.4973965114680973269;
  testrel(6, NE, beip, fbei, cnt,failed);

  fker := -0.6946038911006905212;
  testrel(7, NE, kerp, fker, cnt,failed);

  fkei := 0.3523699133361705344;
  testrel(8, NE, keip, fkei, cnt,failed);

  {----------------------------------}
  x := 8.0;
  kelvin_derx(x, berp, beip, kerp, keip);

  fber := 38.31132570089801858;
  testrel(9, NE, berp, fber, cnt,failed);

  fbei := -7.660318413649826509;
  testrel(10, NE, beip, fbei, cnt,failed);

  fker := -0.8797240991414723762e-3;
  testrele(11, EPS1, kerp, fker, cnt,failed);

  fkei := -0.1336312914858890032e-2;
  testrele(12, EPS2, keip, fkei, cnt,failed);

  {----------------------------------}
  x := 12.0;
  kelvin_derx(x, berp, beip, kerp, keip);
  fber := -472.5688163611958810;
  testrel(13, NE, berp, fber, cnt,failed);

  fbei := 272.6700215595577026;
  testrel(14, NE, beip, fbei, cnt,failed);

  fker := 0.1959385209274561835e-4;
  testrele(15, EPS2, kerp, fker, cnt,failed);

  fkei := 0.7381496295976999158e-4;
  testrele(16, EPS2, keip, fkei, cnt,failed);

  {----------------------------------}
  x := 16.0;
  kelvin_derx(x, berp, beip, kerp, keip);
  fber := 5349.301904593333702;
  testrel(17, NE, berp, fber, cnt,failed);

  fbei := -5999.523615372117048;
  testrel(18, NE, beip, fbei, cnt,failed);

  fker := 0.2280693973119430301e-6;
  testrele(19, EPS1, kerp, fker, cnt,failed);

  fkei := -0.3881116596470087699e-5;
  testrele(20, EPS1, keip, fkei, cnt,failed);

  {----------------------------------}
  x := 20.0;
  kelvin_derx(x, berp, beip, kerp, keip);
  fber := -48803.19784717097646;
  testrel(21, NE, berp, fber, cnt,failed);

  fbei := 111855.0252234971238;
  testrel(22, NE, beip, fbei, cnt,failed);

  fker := -0.7501859210700241327e-7;
  testrel(23, NE, kerp, fker, cnt,failed);

  fkei := 0.1906242756745311635e-6;
  testrel(24, NE, keip, fkei, cnt,failed);

  {----------------------------------}
  x := 100.0;
  kelvin_derx(x, berp, beip, kerp, keip);
  fber := -0.8310516898144792882e29;
  testrel(25, NE1, berp, fber, cnt,failed);

  fbei := 0.1859891445884795168e30;
  testrel(26, NE1, beip, fbei, cnt,failed);

  fker := -0.8766246185882256455e-32;
  testrel(27, NE2, kerp, fker, cnt,failed);

  fkei := 0.2292564824625330302e-31;
  testrel(28, NE1, keip, fkei, cnt,failed);

  {----------------------------------}
  x := 500.0;
  kelvin_derx(x, berp, beip, kerp, keip);
  fber := -0.3103065329961790330e152;
  testrel(29, NE2, berp, fber, cnt,failed);

  fbei := 0.5451844726091405450e152;
  testrel(30, NE1, beip, fbei, cnt,failed);

  fker := -0.4220529386541760937e-155;
  testrel(31, NE2, kerp, fker, cnt,failed);

  fkei := 0.1537225557166260039e-154;
  testrel(32, NE2, keip, fkei, cnt,failed);

  {--------------------------------------------------}
  x := 5.0;
  fber := -3.845339473262154524;
  berp := kelvin_berpx(x);
  testrel(33, NE2, berp, fber, cnt,failed);

  fbei := -4.354140514843110922;
  beip := kelvin_beipx(x);
  testrel(34, NE, beip, fbei, cnt,failed);

  fker := 0.1719340382839311280e-1;
  kerp := kelvin_kerpx(x);
  testrel(35, NE2, kerp, fker, cnt,failed);

  fkei := -0.8199865436307894681e-3;
  keip := kelvin_keipx(x);
  testrele(36, EPS1, keip, fkei, cnt,failed);

  {--------------------------------------------------}
  x := 15.0;
  fber := 91.05533316976338726;
  berp := kelvin_berpx(x);
  testrel(37, NE2, berp, fber, cnt,failed);

  fbei := -4087.755236846416413;
  beip := kelvin_beipx(x);
  testrel(38, NE, beip, fbei, cnt,failed);

  fker := 0.5644678075955970080e-5;
  kerp := kelvin_kerpx(x);
  testrele(39, EPS1, kerp, fker, cnt,failed);

  fkei := -0.5882222803057288351e-5;
  keip := kelvin_keipx(x);
  testrele(40, EPS1, keip, fkei, cnt,failed);

  {--------------------------------------------------}
  x := 13.0;
  fber := -1047.339627111197596;
  berp := kelvin_berpx(x);
  testrel(41, NE, berp, fber, cnt,failed);

  fbei := -192.6057411997334099;
  beip := kelvin_beipx(x);
  testrel(42, NE1, beip, fbei, cnt,failed);

  fker := 0.2969170834567421002e-4;
  kerp := kelvin_kerpx(x);
  testrele(43, EPS3, kerp, fker, cnt,failed);

  fkei := 0.2056476542120472207e-4;
  keip := kelvin_keipx(x);
  testrele(44, EPS3, keip, fkei, cnt,failed);

  {--------------------------------------------------}
  x := 1e-6;
  fber := -0.6249999999999999961e-19;
  berp := kelvin_berpx(x);
  testrel(45, NE1, berp, fber, cnt,failed);

  fbei := 0.5e-6;
  beip := kelvin_beipx(x);
  testrel(46, NE, beip, fbei, cnt,failed);

  fker := -999999.9999996073009;
  kerp := kelvin_kerpx(x);
  testrel(47, NE, kerp, fker, cnt,failed);

  fkei := 0.7215721036811392364e-5;
  keip := kelvin_keipx(x);
  testrel(48, NE, keip, fkei, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_struve_h0x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE  = 10;
  NE2 = 25;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','struve_h0x');

  x := 1e-10;
  f := 0.6366197723675813431e-10;
  y := struve_h0x(x);
  testrel(1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.7943940253295662488e-1;
  y := struve_h0x(x);
  testrel(2, NE, y, f, cnt,failed);

  x := 1.5;
  f := 0.7367234656043998716;
  y := struve_h0x(x);
  testrel(3, NE, y, f, cnt,failed);

  x := 5.0;
  f := -0.1852168157766848901;
  y := struve_h0x(x);
  testrel(4, NE, y, f, cnt,failed);

  x := 10.9375;
  f := -0.1005132970627699984;
  y := struve_h0x(x);
  testrel(5, NE, y, f, cnt,failed);

  x := -11.0625;
  f := 0.12160763382511235722;
  y := struve_h0x(x);
  testrel(6, NE, y, f, cnt,failed);

  x := 50;
  f := -0.8533767482611899895e-1;
  y := struve_h0x(x);
  testrel(7, NE, y, f, cnt,failed);

  x := -200.0;
  f := 0.5108275594755779573e-1;
  y := struve_h0x(x);
  testrel(8, NE2, y, f, cnt,failed);  {!!!!!!!!}

  x := 10000.0;
  f := 0.3711467535586744306e-2;
  y := struve_h0x(x);
  testrel(9, NE, y, f, cnt,failed);

  x := 1e10;
  f := -7.676444513815699932354438350541348934573911340922678e-6;  {Alpha}
  y := struve_h0x(x);
  testrel(10, NE, y, f, cnt,failed);

  x := ldexp(1,200);
  f := 1.775696965839498604282508984743015823013271680925546e-31;  {Alpha}
  y := struve_h0x(x);
  testrel(11, NE, y, f, cnt,failed);

  x := ldexp(1,9000);
  f := 7.68501483154872238119647707669764258531260567653646e-1356; {Alpha}
  y := struve_h0x(x);
  testrel(12, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_struve_h1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','struve_h1x');

  x := 1e-10;
  f := 0.21220659078919378103e-20;
  y := struve_h1x(x);
  testrel(1, NE, y, f, cnt,failed);

  x := 9e-10;
  f := 0.17188733853924696262111255e-18;
  y := struve_h1x(x);
  testrel(2, NE, y, f, cnt,failed);

  x := 10e-10;
  f := 2.122065907891937810e-19;
  y := struve_h1x(x);
  testrel(3, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.3312275639297361862e-2;
  y := struve_h1x(x);
  testrel(4, NE, y, f, cnt,failed);

  x := -1.5;
  f := 0.4102884759694156390;
  y := struve_h1x(x);
  testrel(5, NE, y, f, cnt,failed);

  x := 4.0;
  f := 1.069726661308919359;
  y := struve_h1x(x);
  testrel(6, NE, y, f, cnt,failed);

  x := 8.9375;
  f := 0.73350546732887753298;
  y := struve_h1x(x);
  testrel(7, NE, y, f, cnt,failed);

  x := 9;
  f := 0.7485424374510771033;
  y := struve_h1x(x);
  testrel(8, NE2, y, f, cnt,failed);

  x := 9.0625;
  f := 0.7630757591588199560;
  y := struve_h1x(x);
  testrel(9, NE, y, f, cnt,failed);

  x := 50.0;
  f := 0.580078447945441898996;
  y := struve_h1x(x);
 {$ifdef FPC271or3}
  testrel(10, NE+1, y, f, cnt,failed);
 {$else}
  testrel(10, NE, y, f, cnt,failed);
 {$endif}

  x := -100;
  f := 0.6163111032720133845;
  y := struve_h1x(x);
  testrel(11, NE, y, f, cnt,failed);

  x := 10000;
  f := 0.6437161214863153709;
  y := struve_h1x(x);
  testrel(12, NE, y, f, cnt,failed);

  x := 1e40;
  f := 0.6366197723675813431;  {Alpha}
  y := struve_h1x(x);
  testrel(13, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_struve_l0x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','struve_l0x');

  x := -1e-10;
  f := -0.6366197723675813431e-10;
  y := struve_l0x(x);
  testrel(1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.79715713253115014945e-1;
  y := struve_l0x(x);
  testrel(2, NE, y, f, cnt,failed);

  x := 1.0;
  f := 0.71024318593789088874;
  y := struve_l0x(x);
  testrel(3, NE, y, f, cnt,failed);

  x := -5;
  f := -27.10591712655814655;
  y := struve_l0x(x);
  testrel(4, NE, y, f, cnt,failed);

  x := 10.0;
  f := 2815.652249374594856;
  y := struve_l0x(x);
  testrel(5, NE, y, f, cnt,failed);

  x := 21.9375;
  f := 288526427.5058929094;
  y := struve_l0x(x);
  testrel(6, NE, y, f, cnt,failed);

  x := 22;
  f := 306692993.6113665191;
  y := struve_l0x(x);
  testrel(7, NE, y, f, cnt,failed);

  x := 22.0625;
  f := 326004733.7583088458;
  y := struve_l0x(x);
  testrel(8, NE, y, f, cnt,failed);

  x := -50;
  f := -0.2932553783849336327e21;
  y := struve_l0x(x);
  testrel(9, NE, y, f, cnt,failed);

  x := 100.0;
  f := 0.10737517071310738235e43;
  y := struve_l0x(x);
  testrel(10, NE, y, f, cnt,failed);

  x := -700;
  f := -0.1529593347671873736e303;
  y := struve_l0x(x);
  testrel(11, NE, y, f, cnt,failed);

  x := 11356.0;
  f := 0.26389954891827928230e4930;
  y := struve_l0x(x);
  testrel(12, NE, y, f, cnt,failed);

  x := 11600;
  f := PosInf_x;
  y := struve_l0x(x);
  testabs(13, 0, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_struve_l1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','struve_l1x');

  x := -1e-10;
  f := 0.2122065907891937810e-20;
  y := struve_l1x(x);
  testrel(1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.3319183406689451674e-2;
  y := struve_l1x(x);
  testrel(2, NE, y, f, cnt,failed);

  x := 1.0;
  f := 0.2267643810558086368;
  y := struve_l1x(x);
  testrel(3, NE, y, f, cnt,failed);

  x := -5;
  f := 23+0.7282157804082824458;
  y := struve_l1x(x);
  testrel(4, NE, y, f, cnt,failed);

  x := 10.0;
  f := 2670.358285208482969;
  y := struve_l1x(x);
  testrel(5, NE, y, f, cnt,failed);

  x := 21.9375;
  f := 281871698.9854743394;
  y := struve_l1x(x);
  testrel(6, NE, y, f, cnt,failed);

  x := 22;
  f := 299639606.2420829644;
  y := struve_l1x(x);
  testrel(7, NE, y, f, cnt,failed);

  x := 22.0625;
  f := 318528712.6350923864;
  y := struve_l1x(x);
  testrel(8, NE, y, f, cnt,failed);

  x := -50;
  f := 0.2903078590103556797e21;
  y := struve_l1x(x);
  testrel(9, NE, y, f, cnt,failed);

  x := 100.0;
  f := 0.1068369390338162481e43;
  y := struve_l1x(x);
  testrel(10, NE, y, f, cnt,failed);

  x := -700;
  f := 0.1528500390233900688e303;
  y := struve_l1x(x);
  testrel(11, NE, y, f, cnt,failed);

  x := 11356.0;
  f := 0.2638879292740769163e4930;
  y := struve_l1x(x);
  testrel(12, NE, y, f, cnt,failed);

  x := 11600;
  f := PosInf_x;
  y := struve_l1x(x);
  testabs(13, 0, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_struve_lx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 6; {AMD}
  NA  = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','struve_lx');

  y := struve_lx(2,0.125);
  f := 0.8295489743456552627e-4;
  testrel(1, NE, y, f, cnt,failed);

  y := struve_lx(2,1);
  f := 0.4450783303707983406e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := struve_lx(2,10);
  f := 0.2279458526425006324e4;
  testrel(3, NE, y, f, cnt,failed);

  y := struve_lx(2,20);
  f := 0.3931278100849803299e8;
  testrel(4, NE, y, f, cnt,failed);

  y := struve_lx(2,29);
  f := 0.2726977277708515203e12;
  testrel(5, NE, y, f, cnt,failed);

  y := struve_lx(2,30);
  f := 0.7304368285550353082e12;
  testrel(6, NE, y, f, cnt,failed);

  y := struve_lx(2,40);
{$ifdef BIT16}
  f := 0.2831880997051384763e17*0.5;
{$else}
  f := 0.141594049852569238149e17;
{$endif}
  testrel(7, NE, y, f, cnt,failed);

  y := struve_lx(3,0.125);
  f := 0.1481092567124268210e-5;
  testrel(8, NE, y, f, cnt,failed);

  y := struve_lx(3,1);
  f := 0.6291730749650544377e-2;
  testrel(9, NE, y, f, cnt,failed);

  y := struve_lx(3,10);
  f := 0.1754330742822696564e4;
  testrel(10, NE, y, f, cnt,failed);

  y := struve_lx(3,20);
  f := 0.3459239957188510972e8;
  testrel(11, NE, y, f, cnt,failed);

  y := struve_lx(3,29);
  f := 0.2498186283973945593e12;
  testrel(12, NE, y, f, cnt,failed);

  y := struve_lx(3,30);
  f := 0.6711404617594525287e12;
  testrel(13, NE, y, f, cnt,failed);

  y := struve_lx(3,40);
  f := 0.2658291132946718363e17*0.5;
  testrel(14, NE, y, f, cnt,failed);

  y := struve_lx(22,22);
  f := 0.8808824715288026159e4;
  testrel(15, NE, y, f, cnt,failed);

  y := struve_lx(100,10);
  f := 0.5612962056656603369e-88;
  testrel(16, NE, y, f, cnt,failed);

  y := struve_lx(100,125);
  f := 0.1586799477332954879e37;
  testrel(17, NE, y, f, cnt,failed);

  y := struve_lx(100,200);
  f := 0.4352750449727021914e75;
  testrel(18, NE, y, f, cnt,failed);

  y := struve_lx(0,-5);
  f := -0.2710591712655814655e2;
  testrel(19, NE, y, f, cnt,failed);

  y := struve_lx(2,-30);
  f := -0.7304368285550353082e12;
  testrel(20, NE, y, f, cnt,failed);

  y := struve_lx(3,-30);
  f := 0.6711404617594525287e12;
  testrel(21, NE, y, f, cnt,failed);

  y := struve_lx(4,-30);
  f := -0.5962087360394425753e12;
  testrel(22, NE, y, f, cnt,failed);

  y := struve_lx(10,-5);
  f := -0.3267717256563183426e-2;
  testrel(23, NE, y, f, cnt,failed);

  y := struve_lx(1200,3000);
  f := 0.1341124989174625607e1199*0.5;
  testrel(24, NE, y, f, cnt,failed);

  {non-integer nu}

  y := struve_lx(0.5, 1);
  f := 0.4333156537901020906;
  testrel(25, NE, y, f, cnt,failed);

  y := struve_lx(4.125, 30);
  f := 0.5860632942110845909e12;
  testrel(26, NE, y, f, cnt,failed);

  y := struve_lx(0.125, 1);
  f := 0.6385292223521091581;
  testrel(27, NE, y, f, cnt,failed);

  y := struve_lx(0.125, 20);
  f := 0.43540820615210402417e8;
  testrel(28, NE, y, f, cnt,failed);

  y := struve_lx(0.125, 30);
  f := 0.7814652427883059899e12;
  testrel(29, NE, y, f, cnt,failed);

  y := struve_lx(0.125, 600);
  f := 0.6146225307628329318e259;
  testrel(30, NE, y, f, cnt,failed);

  y := struve_lx(2.5, 0.125);
  f := 0.1148594175782397002e-4;
  testrel(31, NE, y, f, cnt,failed);

  y := struve_lx(2.5, 20.0);
  f := 0.3711237359585533815e8;
  testrel(32, NE, y, f, cnt,failed);

  y := struve_lx(2.5, 30.0);
  f := 0.7031240155028873763e12;
  testrel(33, NE, y, f, cnt,failed);

  y := struve_lx(34.5, 1.0);
  f := 0.2257714211993053785e-50;
  testrel(34, NE, y, f, cnt,failed);

  y := struve_lx(34.5, 20.0);
  f := 0.2699008185941967816e-3;
  testrel(35, NE, y, f, cnt,failed);

  y := struve_lx(34.5, 60.0);
  f := 0.3464808963345829804e21;
  testrel(26, NE, y, f, cnt,failed);

  y := struve_lx(34.5, 600.0);
  f := 0.2278283882763869035e259;
  testrel(37, NE1, y, f, cnt,failed);   {AMD FPC264}

  y := struve_lx(0.0078125, 10.5);
  f := 0.4527364872181681461e4;
  testrel(38, NE, y, f, cnt,failed);

  y := struve_lx(0.0078125, 10.5);
  f := 0.4527364872181681461e4;
  testrel(39, NE, y, f, cnt,failed);

  {negative nu}

  y := struve_lx(-4.125, 30);           {AE}
  f := 0.5860632944022376469e12;
  testrel(40, NE, y, f, cnt,failed);

  y := struve_lx(-0.125, 1);
  f := 0.7796445191058644578;
  testrel(41, NE, y, f, cnt,failed);

  y := struve_lx(-0.125, 20);
  f := 0.4354082064985214698e8;
  testrel(42, NE, y, f, cnt,failed);

  y := struve_lx(-0.125, 30);           {AE}
  f := 0.7814652427883314765e12;
  testrel(43, NE, y, f, cnt,failed);

  y := struve_lx(-0.125, 600);          {AE}
  f := 0.6146225307628329318e259;
  testrel(44, NE, y, f, cnt,failed);

  y := struve_lx(-31.25, 1.0);
  f := 0.3643699813187873390e41;
  testrel(45, NE, y, f, cnt,failed);

  y := struve_lx(-31.25, 10);
  f := 0.1222514141568695836e11;
  testrel(46, NE, y, f, cnt,failed);

  y := struve_lx(-31.25, 20.0);
  f := 3.108364052786170322;
  testrel(47, NE, y, f, cnt,failed);

  y := struve_lx(-31.25, 40.0);
  f := 0.1120269549800670326e12;
  testrel(48, NE, y, f, cnt,failed);

  y := struve_lx(-31.25, 47);           {AE}
  f := 0.5932663640494889554e15;
  testrel(49, NE, y, f, cnt,failed);

  y := struve_lx(-30, 20);
{$ifdef BIT16}
  f := -0.4953916489769681009*2;
{$else}
  f := -0.990783297953936201875;
{$endif}
  testrel(50, NE, y, f, cnt,failed);

  y := struve_lx(-30, 20.77734375);
  f := 0.1569850541827846859e-3;
  testabs(51, NA, y, f, cnt,failed);

  y := struve_lx(-20, 30);              {AE}
  f := 0.1126985104448377122e10;
  testrel(52, NE, y, f, cnt,failed);

  y := struve_lx(-2, -1);
  f := 0.3799053485413077280;
  testrel(54, NE, y, f, cnt,failed);

  y := struve_lx(-3, -1);
  f := 1.746385775221039549;
  testrel(55, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_struve_hx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 6;
  NA  = 2;
  NA1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','struve_hx');

  y := struve_hx(3,-2);
  f := 0.8363766505550048094e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := struve_hx(3,-20);
  f := 0.1734324046607079089e2;
  testrel(2, NE, y, f, cnt,failed);

  y := struve_hx(3,-40);
  f := 0.6811268391613181297e2;
  testrel(3, NE, y, f, cnt,failed);

  y := struve_hx(2,-2);
  f := -0.2803180603538537866;
  testrel(4, NE, y, f, cnt,failed);

  y := struve_hx(2,-20);
  f := -0.4197006936131656458e1;
  testrel(5, NE, y, f, cnt,failed);

  y := struve_hx(2,-38);
  f := -0.8135239017003814325e1;
  testrel(6, NE, y, f, cnt,failed);

  y := struve_hx(2,-40);
  f := -0.8377982783015674911e1;
  testrel(7, NE, y, f, cnt,failed);

  y := struve_hx(3,12);
  f := 0.6466466749652133069e1;
  testrel(8, NE, y, f, cnt,failed);

  y := struve_hx(4,12);
  f := 0.1089385313566373143e2;
  testrel(9, NE, y, f, cnt,failed);

  y := struve_hx(12,12);
  f := 0.1916310028811865144e1;
  testrel(10, NE, y, f, cnt,failed);

  y := struve_hx(15,12);
  f := 0.1724346086687928455;
  testrel(11, NE, y, f, cnt,failed);

  y := struve_hx(20,12);
  f := 0.8174274710118640982e-3;
  testrel(12, NE, y, f, cnt,failed);

  y := struve_hx(50, 0.125);
  f := 0.2025435664804135908e-126;
  testrel(13, NE, y, f, cnt,failed);

  y := struve_hx(50,1);
  f := 0.2305285238913316546e-80;
  testrel(14, NE, y, f, cnt,failed);

  y := struve_hx(50,12);
  f := 0.1617549730830877160e-25;
  testrel(15, NE, y, f, cnt,failed);

  y := struve_hx(50,50);
  f := 0.4337948946341474352e5;
  testrel(16, NE, y, f, cnt,failed);

  y := struve_hx(50,100);
  f := 0.2359715627416222343e20;
  testrel(17, NE, y, f, cnt,failed);

  y := struve_hx(2,0.125);
  f := 0.8283154445038754206e-4;
  testrel(18, NE, y, f, cnt,failed);

  y := struve_hx(2,1);
  f := 0.4046463614479462791e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := struve_hx(2,6);
  f := 0.1617186635063330261e1;
  testrel(20, NE, y, f, cnt,failed);

  y := struve_hx(2,8);
  f := 0.3035387252267535533e1*0.5;
  testrel(21, NE, y, f, cnt,failed);

  y := struve_hx(2,10);
  f := 0.2181688722623384164e1;
  testrel(22, NE, y, f, cnt,failed);

  y := struve_hx(2,12);
  f := 0.2816322778697388064e1;
  testrel(23, NE, y, f, cnt,failed);

  y := struve_hx(2,20);
  f := 0.4197006936131656458e1;
  testrel(24, NE, y, f, cnt,failed);

  y := struve_hx(2,30);
  f := 0.6510412837119774873e1;
  testrel(25, NE, y, f, cnt,failed);

  y := struve_hx(2,40);
  f := 0.8377982783015674911e1;
  testrel(26, NE, y, f, cnt,failed);

  y := struve_hx(2,50);
  f := 0.10718870352203625726e2;
  testrel(27, NE, y, f, cnt,failed);

  y := struve_hx(2,100);
  f := 0.2130386405267446571e2;
  testrel(28, NE, y, f, cnt,failed);

  y := struve_hx(10,0.125);
  f := 0.5389034941017067878e-20;
  testrel(29, NE, y, f, cnt,failed);

  y := struve_hx(10,1);
  f := 0.4563623880470109001e-10;
  testrel(30, NE, y, f, cnt,failed);

  y := struve_hx(10,4);
  f := 0.1544763075114888737e-3;
  testrel(31, NE, y, f, cnt,failed);

  y := struve_hx(10,6);
  f := 0.1013711880650545025e-1;
  testrel(32, NE, y, f, cnt,failed);

  y := struve_hx(10,8);
  f := 0.16687353850333171213;
  testrel(33, NE, y, f, cnt,failed);

  y := struve_hx(10,9);
  f := 0.4966834740182485912;
  testrel(34, NE, y, f, cnt,failed);

  y := struve_hx(10,20);
  f := 0.5251927521745221582e3;
  testrel(35, NE, y, f, cnt,failed);

  y := struve_hx(10,30);
  f := 0.1956771527019647543e5;
  testrel(36, NE, y, f, cnt,failed);

  y := struve_hx(20,40);
  f := 0.5615220059442179315e7;
  testrel(37, NE, y, f, cnt,failed);

  y := struve_hx(30,4);
  f := 0.1514561328076532737e-23;
  testrel(38, NE, y, f, cnt,failed);

  y := struve_hx(30,5);
  f := 0.1459461294786301802e-20;
  testrel(39, NE, y, f, cnt,failed);

  y := struve_hx(30,10);
  f := 0.2157359166313841896e-11;
  testrel(40, NE, y, f, cnt,failed);

  y := struve_hx(30,20);
  f := 0.1460043012228651387e-2;
  testrel(41, NE, y, f, cnt,failed);

  y := struve_hx(30,30);
  f := 0.1623595353884071012e3;
  testrel(42, NE, y, f, cnt,failed);

  y := struve_hx(30,40);
  f := 0.6542690683966296798e6;
  testrel(43, NE, y, f, cnt,failed);

  {non-integer nu}

  y := struve_hx(0.5, Pi);
  f := 0.9003163161571060696;
  testrel(44, NE, y, f, cnt,failed);

  y := struve_hx(0.5, 60.0);
  f := 0.2011111376078884702;
  testrel(45, NE, y, f, cnt,failed);

  y := struve_hx(2.5, 20.0);
  f := 0.9058993618079528111e1;
  testrel(46, NE, y, f, cnt,failed);

  y := struve_hx(34.5, 20.0);
  f := 0.7692959060558401110e-5;
  testrel(47, NE, y, f, cnt,failed);

  y := struve_hx(34.5, 40.0);
  f := 0.7702896071503372942e5;
  testrel(48, NE1, y, f, cnt,failed);

  y := struve_hx(34.5, 60.0);
  f := 0.5935220289490881169e11;
  testrel(49, NE, y, f, cnt,failed);

  y := struve_hx(0.125, 0.5);
  f := 0.2578567948941973427;
  testrel(50, NE, y, f, cnt,failed);

  y := struve_hx(0.125, 2);
  f := 0.8149844854995485319;
  testrel(51, NE, y, f, cnt,failed);

  y := struve_hx(0.125, 5);          {int + Yv: abs}
  f := -0.9618309457354139700e-1;
  testabs(52, NA1, y, f, cnt,failed);

  y := struve_hx(0.125, 30);         {int + Yv: abs}
  f := -0.6146886542830816681e-1;
  testabs(53, NA, y, f, cnt,failed);

  y := struve_hx(0.125, 40);         {AE  + Yv: abs}
  f := 0.1506693495161465447;
  testabs(54, NA, y, f, cnt,failed);

  y := struve_hx(0.25, 5);           {int + Yv: abs}
  f := 0.8861990084230123199e-2;
  testabs(55, NA, y, f, cnt,failed);

  y := struve_hx(0.25, 4);
  f := 0.4008926305910176599;
  testrel(56, NE1, y, f, cnt,failed);

  y := struve_hx(0.625, 5);
  f := 0.3910731128047716249;
  testrel(57, NE, y, f, cnt,failed);

  y := struve_hx(0.0078125, 10.5);   {int + Yv: abs}
  f := -0.2779771777082488473e-2;
  testabs(58, NA, y, f, cnt,failed);

  {negative nu}

  y := struve_hx(-0.125, 0.5);
  f := 0.3660759119554399380;
  testrel(59, NE, y, f, cnt,failed);

  y := struve_hx(-0.125, 2);
  f := 0.7492415174287243325;
  testrel(60, NE, y, f, cnt,failed);

  y := struve_hx(-0.125, 5);         {int}
  f := -0.2558580251280904802;
  testabs(61, NA1, y, f, cnt,failed);   {NA1 for FPC 2.4.2+}

  y := struve_hx(-0.125, 30);        {int}
  f := -0.1206128654483063351;
  testabs(62, NA, y, f, cnt,failed);

  y := struve_hx(-0.125, 40);        {AE}
  f := 0.1331278346986515631;
  testabs(63, NA, y, f, cnt,failed);

  y := struve_hx(-31.25, 0.125);
  f := 0.7629289376723737283e68;
  testrel(64, NE, y, f, cnt,failed);

  y := struve_hx(-31.25, 20);
  f := 267.5978503479011890;
  testrel(65, NE, y, f, cnt,failed);

  y := struve_hx(-31.25, 40);        {int}
  f := 0.1077952026367204399;
  testabs(66, NA, y, f, cnt,failed);

  y := struve_hx(-31.25, 48);        {AE}
  f := 0.9127208388797704561e-2;
  testabs(67, NA1, y, f, cnt,failed);

  y := struve_hx(-10, 30);           {int}
  f := 0.7505672513983422036e-1;
  testabs(68, NA, y, f, cnt,failed);

  y := struve_hx(-3, 10);            {int}
  f := 0.2504621870054056667;
  testabs(69, NA, y, f, cnt,failed);

  y := struve_hx(-90, 125);          {int}
  f := -0.8043851578747538147e-1;
  testabs(70, NA, y, f, cnt,failed);

  y := struve_hx(-90, 140);          {AE}
  f := 0.2684228605590296111e-1;
  testabs(71, NA, y, f, cnt,failed);

  y := struve_hx(-2, -1);
  f := 0.8083617270119804962;
  testrel(72, NE, y, f, cnt,failed);

  y := struve_hx(-3, -1);
  f := 2.158664699514703698;
  testrel(73, NE, y, f, cnt,failed);

  {extended only}
  y := struve_hx(1200,3000);
  f := 0.4185705650992682416e634;
  testrel(74, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.
