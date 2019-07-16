{Part 1a of regression test for SPECFUNX unit  (c) 2010  W.Ehrhardt}

unit t_sfx1a;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_psix;
procedure test_psistarx;
procedure test_psi_invx;
procedure test_trigammax;
procedure test_tetrapentagammax;
procedure test_polygammax;
procedure test_BatemanGx;

procedure test_lnbetax;
procedure test_betax;
procedure test_ibetax;
procedure test_ibeta_invx;
procedure test_beta3x;


implementation

uses
  amath, specfunx, t_sfx0;


{---------------------------------------------------------------------------}
procedure test_psix;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 3;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','psix');

  x := PosInf_x;
  y := psix(x);
  f := PosInf_x;
  testabs( 1, 0, y, f, cnt,failed);

  x := MaxDouble;
  y := psix(x);
  f := 709.7827128933839968;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1.46163214496836234126;
  y := psix(x);
  f := 0;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1.0;
  y := psix(x);
  f := -EulerGamma;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1e12;
  y := psix(x);
  f := 27.63102111592804821;
  testrel( 5, NE, y, f, cnt,failed);

  x := 12.125;
  y := psix(x);
  f := 2.453465874567354846;
  testrel( 6, NE, y, f, cnt,failed);

  x := 12.0;
  y := psix(x);
  f := 2.442661679975812017;
  testrel( 7, NE, y, f, cnt,failed);

  x := 11.875;
  y := psix(x);
  f := 2.431739553364003224;
  testrel( 8, NE, y, f, cnt,failed);

  x := 1.46875;
  y := psix(x);
  f := 0.6865411470735776728e-2;
  testrel( 9, NE, y, f, cnt,failed);

  x := 1.453125;
  y := psix(x);
  f := -0.8264331530773935209e-2;
  testrel(10, NE, y, f, cnt,failed);

  x := 1.461632145;
  y := psix(x);
  f := 0.3061488427062461183e-10;
  testabs(11, 1, y, f, cnt,failed);

  x := 0.001953125;
  y := psix(x);
  f := -512.5740074804865255;
  testrel(12, NE, y, f, cnt,failed);

  x := -0.001953125;
  y := psix(x);
  f := 511.4195669796869178;
  testrel(13, NE, y, f, cnt,failed);

  x := -0.5040830082644554093;
  y := psix(x);
  f := 0;
  testabs(14, NE, y, f, cnt,failed);

  x := -0.984375;
  y := psix(x);
  f := -63.53592995729484080;
  testrel(15, NE, y, f, cnt,failed);

  x := -1.015625;
  y := psix(x);
  f := 64.38139996249215472;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_psistarx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 2;
  NE1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','psistarx');

  y := psistarx(1e-18);
  f := -999999999999999959.1;
  testrel( 1, NE, y, f, cnt,failed);

  y := psistarx(1e-9);
  f := -999999979.8539498263;
  testrel( 2, NE1, y, f, cnt,failed);

  y := psistarx(1/1024);
  f := -1017.644138623741169;
  testrel( 3, NE1, y, f, cnt,failed);

  y := psistarx(1/32);
  f := -29.06121738580639156;
  testrel( 4, NE, y, f, cnt,failed);

  y := psistarx(0.5);
  f := -1.270362845461478170;
  testrel( 5, NE1, y, f, cnt,failed);

  y := psistarx(1.4616241455078125);
  f := -0.3795559867445265848;
  testrel( 6, NE, y, f, cnt,failed);

  y := psistarx(2);
  f := -0.2703628454614781700;
  testrel( 7, NE, y, f, cnt,failed);

  y := psistarx(3);
  f := -0.1758279535696425520;
  testrel( 8, NE, y, f, cnt,failed);

  y := psistarx(4);
  f := -0.1301766926880901461;
  testrel( 9, NE, y, f, cnt,failed);

  y := psistarx(5);
  f := -0.1033202440022999019;
  testrel(10, NE, y, f, cnt,failed);

  y := psistarx(6);
  f := -0.8564180079625452809e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := psistarx(7);
  f := -0.7312581395684616571e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := psistarx(8);
  f := -0.6380006372422593172e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := psistarx(8.5);
  f := -0.5997534591685064943e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := psistarx(11);
  f := -0.4614268373164943641e-1;
  testrel(15, NE1, y, f, cnt,failed);

  y := psistarx(12);
  f := -0.4224496981218829349e-1;
  testrel(16, NE1, y, f, cnt,failed);

  y := psistarx(13);
  f := -0.3895434415239138598e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := psistarx(14);
  f := -0.3613923938303634137e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := psistarx(20);
  f := -0.2520828131184194256e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := psistarx(100);
  f := -0.5008333250003967837e-2;
  testrel(20, NE, y, f, cnt,failed);

  y := psistarx(12345);
  f := -0.4050277443266634188e-4;
  testrel(21, NE, y, f, cnt,failed);

  y := psistarx(1e9);
  f := -0.5000000000833333333e-9;
  testrel(22, NE, y, f, cnt,failed);

  y := psistarx(1e20);
  f := -0.5e-20;
  testrel(23, NE, y, f, cnt,failed);

  y := psistarx(MaxDouble);
  f := -0.2781342323134002038e-308;
  testrel(24, NE, y, f, cnt,failed);

  y := psistarx(PosInf_x);
  f := 0;
  testrel(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_psi_invx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 3;
  NE1 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','psi_invx');

  y := psi_invx(654);
  f := 0.1068048968017990679e285;
  testrel(1, NE, y, f, cnt,failed);

  y := psi_invx(100);
  f := 0.2688117141816135448e44;
  testrel(2, NE, y, f, cnt,failed);

  y := psi_invx(49);  {AMath y>=ym}
  f := 0.1907346572495099691e22;
  testrel(3, NE, y, f, cnt,failed);

  y := psi_invx(41);  {DAMath y>=ym}
  f := 0.6398434935300549497e18;
  testrel(4, NE, y, f, cnt,failed);

  y := psi_invx(40);
  f := 0.2353852668370199859e18;
  testrel(5, NE1, y, f, cnt,failed);

  y := psi_invx(30);
  f := 0.1068647458152496215e14;
  testrel(6, NE1, y, f, cnt,failed);

  y := psi_invx(15);
  f := 0.3269017872472097893e7;
  testrel(7, NE1, y, f, cnt,failed);

  y := psi_invx(5);
  f := 148.9128783562188555;
  testrel(8, NE1, y, f, cnt,failed);

  y := psi_invx(2);
  f := 7.883428631186041039;
  testrel(9, NE, y, f, cnt,failed);

  y := psi_invx(1);
  f := 3.203171468376931069;
  testrel(10, NE, y, f, cnt,failed);

  y := psi_invx(0.875);
  f := 2.881817722916420026;
  testrel(11, NE, y, f, cnt,failed);

  y := psi_invx(0.5);
  f := 2.124337554428001093;
  testrel(12, NE, y, f, cnt,failed);

  y := psi_invx(0);
  f := 1.4616321449683623413;
  testrel(13, NE, y, f, cnt,failed);

  y := psi_invx(-Eulergamma);
  f := 1;
  testrel(14, NE, y, f, cnt,failed);

  y := psi_invx(-0.875);
  f := 0.8406646236282619084;
  testrel(15, NE, y, f, cnt,failed);

  y := psi_invx(-1);
  f := 0.7850033253745206691;
  testrel(16, NE, y, f, cnt,failed);

  y := psi_invx(-2);
  f := 0.4926978052047480760;
  testrel(17, NE, y, f, cnt,failed);

  y := psi_invx(-4);
  f := 0.2639162080054398413;
  testrel(18, NE, y, f, cnt,failed);

  y := psi_invx(-8);
  f := 0.1312314676887234664;
  testrel(19, NE, y, f, cnt,failed);

  y := psi_invx(-10);
  f := 0.1043571987701165077;
  testrel(20, NE1, y, f, cnt,failed);

  y := psi_invx(-100);
  f := 0.10056395666750782055e-1;
  testrel(21, NE1, y, f, cnt,failed);

  y := psi_invx(-1e9);
  f := 0.1000000000577215664e-8;
  testrel(22, NE1, y, f, cnt,failed);

  y := psi_invx(-1e10);
  f:= 0.1000000000057721566e-9;
  testrel(23, NE1, y, f, cnt,failed);

  y := psi_invx(-1e14);
  f := 0.1000000000000005772e-13;
  testrel(24, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_trigammax;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','trigammax');

  x := 1.0;
  y := trigammax(x);
  f := 1.6449340668482264365;
  testrel( 1, NE, y, f, cnt,failed);

  x := 2.0;
  y := trigammax(x);
  f := 0.64493406684822643647;
  testrel( 2, NE, y, f, cnt,failed);

  x := 3.0;
  y := trigammax(x);
  f := 0.39493406684822643647;
  testrel( 3, NE, y, f, cnt,failed);

  x := 4.0;
  y := trigammax(x);
  f := 0.28382295573711532536;
  testrel( 4, NE, y, f, cnt,failed);

  x := 5.0;
  y := trigammax(x);
  f := 0.22132295573711532536;
  testrel( 5, NE, y, f, cnt,failed);

  x := 6.0;
  y := trigammax(x);
  f := 0.18132295573711532536;
  testrel( 6, NE, y, f, cnt,failed);

  x := 7.0;
  y := trigammax(x);
  f := 0.15354517795933754758;
  testrel( 7, NE, y, f, cnt,failed);

  x := 8.0;
  y := trigammax(x);
  f := 0.13313701469403142513;
  testrel( 8, NE, y, f, cnt,failed);

  x := 9.0;
  y := trigammax(x);
  f := 0.11751201469403142513;
  testrel( 9, NE, y, f, cnt,failed);

  x := 10.0;
  y := trigammax(x);
  f := 0.10516633568168574612;
  testrel(10, NE, y, f, cnt,failed);

  x := 11.0;
  y := trigammax(x);
  f := 0.95166335681685746122e-1;
  testrel(11, NE, y, f, cnt,failed);

  x := 11.875;
  y := trigammax(x);
  f := 0.87855620835225893327e-1;
  testrel(12, NE, y, f, cnt,failed);

  x := 12.0;
  y := trigammax(x);
  f := 0.86901872871768390750e-1;
  testrel(13, NE, y, f, cnt,failed);

  x := 12.125;
  y := trigammax(x);
  f := 0.85968597520282084583e-1;
  testrel(14, NE, y, f, cnt,failed);

  x := 0.001953125;
  y := trigammax(x);
  f := 262145.64025088744769;
  testrel(15, NE, y, f, cnt,failed);

  x := -0.001953125;
  y := trigammax(x);
  f := 262145.64964201880968;
  testrel(16, NE, y, f, cnt,failed);

  x := -0.5;
  y := trigammax(x);
  f := 8.9348022005446793094;
  testrel(17, NE, y, f, cnt,failed);

  x := -0.984375;
  y := trigammax(x);
  f := 4098.6401449656114839;
  testrel(18, NE, y, f, cnt,failed);

  x := -1.015625;
  y := trigammax(x);
  f := 4098.6527746490979288;
  testrel(19, NE, y, f, cnt,failed);

  x := 1e-8;
  y := trigammax(x);
  f := 10000000000000001.0 + 0.64493404;
  testrel(20, NE, y, f, cnt,failed);

  x := 1e-7;
  y := trigammax(x);
  f := 100000000000001.0 + 0.6449338264;
  testrel(21, NE, y, f, cnt,failed);

  x := 1e-6;
  y := trigammax(x);
  f := 1000000000001.0 + 0.644931662738;
  testrel(22, NE, y, f, cnt,failed);

  x := 1e-5;
  y := trigammax(x);
  f := 10000000001.0 + 0.64491002603486;
  testrel(23, NE, y, f, cnt,failed);

  x := 1e-4;
  y := trigammax(x);
  f := 100000001.0 + 0.6446936879331443;
  testrel(24, NE, y, f, cnt,failed);

  x := 1e-3;
  y := trigammax(x);
  f := 1000001.0 + 0.642533195868978033;
  testrel(25, NE, y, f, cnt,failed);

  x := -1e-4;
  y := trigammax(x);
  f := 100000001.0 + 0.6451745107027036;
  testrel(26, NE, y, f, cnt,failed);

  x := 1e-3;
  y := trigammax(x);
  f := 1000001.0 + 0.642533195868978033;
  testrel(27, NE, y, f, cnt,failed);

  x := -0.1;
  y := trigammax(x);
  f := 101.92253995947720352;
  testrel(28, NE, y, f, cnt,failed);

  x := -1000000000.5;
  y := trigammax(x);
  f := 9.86960440008935862;
  testrel(29, NE, y, f, cnt,failed);

  x := -12345689.875;
  y := trigammax(x);
  f := 67.393874321913027163;
  testrel(30, NE, y, f, cnt,failed);

  x := 100.125;
  y := trigammax(x);
  f := 0.10037556879391693663e-1;
  testrel(31, NE, y, f, cnt,failed);

  x := 10000.375;
  y := trigammax(x);
  f := 0.10000124993228873720e-3;
  testrel(32, NE, y, f, cnt,failed);

  x := 1000000.5;
  y := trigammax(x);
  f := 0.99999999999991666667e-6;
  testrel(33, NE, y, f, cnt,failed);

  x := 1e18;
  y := trigammax(x);
  f := 0.10000000000000000005e-17;
  testrel(34, NE, y, f, cnt,failed);

  x := 1e20;
  y := trigammax(x);
  f := 1e-20;
  testrel(35, NE, y, f, cnt,failed);

  x := 1e300;
  y := trigammax(x);
  f := 1e-300;
  testrel(36, NE, y, f, cnt,failed);

  x := PosInf_x;
  y := trigammax(x);
  f := 0;
  testabs(37, 0, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_tetrapentagammax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','tetra/pentagammax');

  y := tetragammax(1/1024);
  f := -2147483650.3977839164;
  testrel( 1, NE, y, f, cnt,failed);

  y := tetragammax(0.5);
  f := -16.828796644234319996;
  testrel( 2, NE, y, f, cnt,failed);

  y := tetragammax(1);
  f := -2.4041138063191885708;
  testrel( 3, NE, y, f, cnt,failed);

  y := tetragammax(10);
  f := -0.11049834970802067462e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := tetragammax(1000);
  f := -0.1001000499999833333e-5;
  testrel( 5, NE, y, f, cnt,failed);

  y := tetragammax(1e6);
  f := -0.10000010000005000000e-11;
  testrel( 6, NE, y, f, cnt,failed);

  y := tetragammax(1e10);
  f := -0.10000000001000000000e-19;
  testrel( 7, NE, y, f, cnt,failed);

  y := tetragammax(-1/1024);
  f := 2147483645.5895325703;
  testrel( 8, NE, y, f, cnt,failed);

  y := tetragammax(-0.5);
  f := -0.82879664423431999560;
  testrel( 9, NE, y, f, cnt,failed);

  y := tetragammax(-4.0625);
  f := 8191.1307641388450279;
  testrel(10, NE, y, f, cnt,failed);

  y := tetragammax(-2.25);
  f := 123.89694977406016558;
  testrel(11, NE, y, f, cnt,failed);

  y := pentagammax(1/1024);
  f := 6597069766662.4696945;
  testrel(12, NE, y, f, cnt,failed);

  y := pentagammax(0.5);
  f := 97.409091034002437236;
  testrel(13, NE, y, f, cnt,failed);

  y := pentagammax(1);
  f := 6.4939394022668291491;
  testrel(14, NE, y, f, cnt,failed);

  y := pentagammax(10);
  f := 0.23199013042898683856e-2;
  testrel(15, NE, y, f, cnt,failed);

  y := pentagammax(1000);
  f := 0.20030019999990000013e-8;
  testrel(16, NE, y, f, cnt,failed);

  y := pentagammax(1e6);
  f := 0.20000030000020000000e-17;
  testrel(17, NE, y, f, cnt,failed);

  y := pentagammax(1e10);
  f := 0.20000000003000000000e-29;
  testrel(18, NE, y, f, cnt,failed);

  y := pentagammax(-1/1024);
  f := 6597069766662.5183007;
  testrel(19, NE, y, f, cnt,failed);

  y := pentagammax(-0.5);
  f := 193.40909103400243724;
  testrel(20, NE, y, f, cnt,failed);

  y := pentagammax(-4.0625);
  f := 393229.45068149631385;
  testrel(21, NE, y, f, cnt,failed);

  y := pentagammax(-2.25);
  f := 1558.4550231887143400;
  testrel(22, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_polygammax;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NEN = 6;
  NEN1 = 8;
  NEN2 = 24;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','polygammax');

  y := polygammax(5,1/1024);
  f := 0.13835058055282163724e21;
  testrel( 1, NE, y, f, cnt,failed);

  y := polygammax(10,1/1024);
  f := -0.47104517100627956246e40;
  testrel( 2, NE, y, f, cnt,failed);

  y := polygammax(51,1/1024);
  f := 0.53240581940759651815e223;
  testrel( 3, NE, y, f, cnt,failed);

  y := polygammax(500,1/1024);
  f := -0.17647647144055629201e2643;
  testrel( 4, NE, y, f, cnt,failed);

  y := polygammax(5,0.5);
  f := 7691.1135486024354962;
  testrel( 5, NE, y, f, cnt,failed);

  y := polygammax(10,0.5);
  f := -7431824508.8587689755;
  testrel( 6, NE, y, f, cnt,failed);

  y := polygammax(51,0.5);
  f := 0.69856178393124431546e82;
  testrel( 7, NE, y, f, cnt,failed);

  y := polygammax(500,0.5);
  f := -0.79879688530950177226e1285;
  testrel( 8, NE, y, f, cnt,failed);

  y := polygammax(5,5);
  f := 0.12261509635954378911e-1;
  testrel( 9, NE, y, f, cnt,failed);

  y := polygammax(10,5);
  f := -0.86751075791965813173e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := polygammax(51,5);
  f := 0.69861510803247881420e30;
  testrel(11, NE, y, f, cnt,failed);

  y := polygammax(500,5);
  f := -0.79879688530950177226e784;
  testrel(12, NE, y, f, cnt,failed);

  y := polygammax(5,500);
  f := 0.77184767997132824576e-12;
  testrel(13, NE, y, f, cnt,failed);

  y := polygammax(10,500);
  f := -0.37531863599270351615e-21;
  testrel(14, NE, y, f, cnt,failed);

  y := polygammax(51,500);
  f := 0.72039788797923852198e-73;
  testrel(15, NE, y, f, cnt,failed);

  y := polygammax(500,500);
  f := -0.12637986018904568831e-217;
  testrel(16, NE, y, f, cnt,failed);

  y := polygammax(5,1e6);
  f := 0.24000060000060000000e-28;
  testrel(17, NE, y, f, cnt,failed);

  y := polygammax(10,1e6);
  f := -0.36288181440332640000e-54;
  testrel(18, NE, y, f, cnt,failed);

  y := polygammax(51,1e6);
  f := 0.30414868767811536332e-241;
  testrel(19, NE, y, f, cnt,failed);

  y := polygammax(500,1e6);
  f := -0.24408837713359279632e-1868;
  testrel(20, NE, y, f, cnt,failed);

  y := polygammax(2001,32);
  f := 0.32481871300239838228e2726;
{$ifdef FPC}
  testrel(21, 10, y, f, cnt,failed);  {!!!???}
{$else}
  testrel(21, NE, y, f, cnt,failed);
{$endif}

  y := polygammax(2001,1024);
  f := 0.18525166336666044149e-287;
  testrele(22, 3e-16, y, f, cnt,failed);   {!!! exp/ln inaccuracy}

  y := polygammax(1800,980);
  f := -0.46177041402704490422e-307;
  testrel(23, 2*NE, y, f, cnt,failed);

  x := -1.25;
  y := polygammax(4,x);
  f := 24481.06371418458253;
  testrel(24, NEN, y, f, cnt,failed);

  y := polygammax(10,x);
  f := 15220204740174.33759;
  testrel(25, NEN, y, f, cnt,failed);

  x := -2.75;
  y := polygammax(5,x);
  f := 492231.2077521004215;
  testrel(26, NEN1, y, f, cnt,failed);

  y := polygammax(6,x);
  f := -11791224.84887192732;
  testrel(27, NEN1, y, f, cnt,failed);

  x := -ldexp(1,-30);
  y := polygammax(5,x);
  f := 0.1838994649039066630e57;
  testrel(28, NEN, y, f, cnt,failed);

  x := -1 + ldexp(1,-20);
  y := polygammax(6,x);
  f := -0.1003533533933878041e46;
  testrel(29, NEN, y, f, cnt,failed);

  x := -5 + 1/32;
  y := polygammax(1,x);
  f := 1027.113866143815380;
  testrel(30, NEN, y, f, cnt,failed);

  y := polygammax(2,x);
  f := -65535.62604688783524;
  testrel(31, NEN, y, f, cnt,failed);

  y := polygammax(3,x);
  f := 6291469.095472251473;
  testrel(32, NEN, y, f, cnt,failed);

  y := polygammax(4,x);
  f := -805306360.3247959621;
  testrel(33, NEN, y, f, cnt,failed);

  y := polygammax(5,x);
  f := 128849019129.1285377;
  testrel(34, NEN, y, f, cnt,failed);

  y := polygammax(6,x);
  f := -24739011624640.00489;
  testrel(35, NEN, y, f, cnt,failed);

  y := polygammax(7,x);
  f := 5541538604001519.012;
  testrel(36, NEN, y, f, cnt,failed);

  y := polygammax(8,x);
  f := -1418633882621683128.2;
  testrel(37, NEN, y, f, cnt,failed);

  y := polygammax(9,x);
  f := 0.4085665581950521631e21;
  testrel(38, NEN, y, f, cnt,failed);

  y := polygammax(10,x);
  f := -0.1307412986224164445e24;
  testrel(39, NEN, y, f, cnt,failed);

  x := -4 - 1/16;
  y := polygammax(1,x);
  f := 259.0970762766438360;
  testrel(40, NEN, y, f, cnt,failed);

  y := polygammax(2,x);
  f := 8191.130764138845028;
  testrel(41, NEN, y, f, cnt,failed);

  y := polygammax(3,x);
  f := 393229.4506814963139;
  testrel(42, NEN, y, f, cnt,failed);

  y := polygammax(4,x);
  f := 25165808.30887021410;
  testrel(43, NEN, y, f, cnt,failed);

  y := polygammax(5,x);
  f := 2013266184.387213715;
  testrel(44, NEN, y, f, cnt,failed);

  y := polygammax(6,x);
  f := 193273527657.2144316;
  testrel(45, NEN, y, f, cnt,failed);

  y := polygammax(7,x);
  f := 21646635183431.80948;
  testrel(46, NEN, y, f, cnt,failed);

  y := polygammax(8,x);
  f := 2770769301946764.682;
  testrel(47, NEN, y, f, cnt,failed);

  y := polygammax(9,x);
  f := 398990779488245464.4;
  testrel(48, NEN, y, f, cnt,failed);

  y := polygammax(10,x);
  f := 0.6383852471797126186e20;
  testrel(49, NEN, y, f, cnt,failed);

  y := polygammax(11,-4.96875);
  f := 0.4602093711509058946e26;
  testrel(50, NEN, y, f, cnt,failed);

  y := polygammax(12,-4.96875);
  f := -0.1767203985219478632e29;
  testrel(51, NEN, y, f, cnt,failed);

  y := polygammax(11,-4.0625);
  f := 0.1123558035036401932e23;
  testrel(52, NEN, y, f, cnt,failed);

  y := polygammax(12,-4.0625);
  f := 0.2157231427269870486e25;
  testrel(53, NEN, y, f, cnt,failed);

  {Recursive reflection coefficient calculation}
  y := polygammax(20,-0.25);
  f := 0.1070001618789627526e32;
  testrel(54, NEN, y, f, cnt,failed);

  y := polygammax(20,-0.5);
  f := -0.4877729494626098736e15;
  testrel(55, NEN, y, f, cnt,failed);

  y := polygammax(20,-5.25);
  f := 0.1070001618789629770e32;
  testrel(56, NEN, y, f, cnt,failed);

  y := polygammax(20,-5.5);
  f := -0.2175607673650773455e2;
  testrel(57, NEN, y, f, cnt,failed);

  y := polygammax(25,-5.5);
  f := 0.2081879370547485528e34;
  testrel(58, NEN, y, f, cnt,failed);

  y := polygammax(25,-5.125);
  f := 0.4687975578712144987e49;
  testrel(59, NEN, y, f, cnt,failed);

  y := polygammax(43,-15.0625);
  f := 0.5786632577803752727e106;
  testrel(60, NEN, y, f, cnt,failed);

  y := polygammax(50,-15.0625);
  f := 0.7819770151595511607e126;
  testrel(61, NEN, y, f, cnt,failed);

  y := polygammax(120,-10.25);
  f := 0.4727731434981616198e272;
  testrel(62, NEN, y, f, cnt,failed);

  {------}
  f := -0.5998777844943402075e118;
  y := polygammax(100,-1.5);
  testrel(63, NEN, y, f, cnt,failed);

  f := -156.2295965932329948;
  y := polygammax(10,-1.5);
  testrel(64, NEN, y, f, cnt,failed);

  f := 0.9559046668383224511e191;
  y := polygammax(101, -1.5);
  testrel(65, NEN, y, f, cnt,failed);

  f := 0.2366100660490897156e187;
  y := polygammax(99, -1.5);
  testrel(66, NEN, y, f, cnt,failed);

  f := 0.4285818862427966147e27;
  y := polygammax(21, -1.5);
  testrel(67, NEN, y, f, cnt,failed);

  f := 0.15421777401019709997e96;
  y := polygammax(50, -1.25);
  testrel(68, NEN, y, f, cnt,failed);

  x := -1.0078125;
  f := 0.8928599125084502451e172;
  y := polygammax(50, x);
  testrel(69, NEN, y, f, cnt,failed);

  x := -0.5-ldexp(1,-30);
  f := -0.1301172718277658206e74;
  y := polygammax(50, x);
  testrel(70, NEN2, y, f, cnt,failed);

  x := -0.5-ldexp(1,-12);
  f := -0.3411319768828644233e79;
  y := polygammax(50, x);
  testrel(71, NEN2, y, f, cnt,failed);

  x := -0.5-ldexp(1,-10);
  f := -0.1366770526924152279e80;
  y := polygammax(50, x);
  testrel(72, NEN, y, f, cnt,failed);

  x := -0.25;
  f := 0.1522020442846279112e14;
  y := polygammax(10, x);
  testrel(73, NEN, y, f, cnt,failed);

  x := -0.1;
  f := 0.3628799999884331908e18;
  y := polygammax(10, x);
  testrel(74, NEN, y, f, cnt,failed);

  x := -0.125;
  f := 0.2243956035062425046e38;
  y := polygammax(20, x);
  testrel(75, NEN, y, f, cnt,failed);

  x := -1.0/THREE;
  f := 0.2544900217877993057e29;
  y := polygammax(20, x);
  testrel(76, NEN, y, f, cnt,failed);

  {Adjust for double}
  x := -0.5-ldexp(1,-60);
  f := -0.1211811529185720110e65;
  y := polygammax(50, x); {cot poly}
  testrel(77, NEN2, y, f, cnt,failed);

  x := -0.5-ldexp(1,-60);
  f := -0.69671205975479281090e106;
  y := polygammax(70, x);  {cot series}
  testrel(78, NEN, y, f, cnt,failed);

  x := -0.5+ldexp(1,-60);
  f := 0.6967120597547928101e106;
  y := polygammax(70, x);
  testrel(79, NEN1, y, f, cnt,failed);  { <= FPC264}

  x := -0.5+ldexp(1,-60);
  f := 0.8032543162102312049e124;
  y := polygammax(71, x);
  testrel(80, NEN, y, f, cnt,failed);


  {Extended only}
  {-------- Large parameters: no overflow with reasonable acc. ----------}
  y := polygammax(3000,50);
  f := -0.1020934928135008429e4033;
  testrele(81, 8e-16, y, f, cnt,failed);

  y := polygammax(3000,200);
  f := -0.1686414097933211677e2226;
  testrele(82, 2e-15, y, f, cnt,failed);

  {ASE overflow, direct sum}
  y := polygammax(30000,7560);
  f := -0.8255438888141577822e4928;
  testrele(83, 2e-15, y, f, cnt,failed);

  y := polygammax(1400,-10.25);
  f := 0.1059748811873522956e4643;
  testrel(84, NEN, y, f, cnt,failed);

  x := -ldexp(1,-30);
  f := 0.11454291687128789153e526;
  y := polygammax(50, x);
  testrel(85, NEN, y, f, cnt,failed);

  x := -ldexp(1,-30);
  y := polygammax(437, x);
  f := 0.2276438170663825936e4922;
  testrel(86, NEN, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_BatemanGx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE2 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','BatemanGx');

  {Maple: bg:= x->Psi((x+1)/2)-Psi(x/2);}
  y := BatemanGx(0.5);
  f := Pi;
  testrel(1, NE, y, f, cnt,failed);

  y := BatemanGx(1e19);
  f := 0.100000000000000000005e-18;
  testrel(2, NE, y, f, cnt,failed);

  y := BatemanGx(1e6);
  f := 0.5000002499999999999e-6*2;  {f=0.100000049999999999975e-5}
  testrel(3, NE, y, f, cnt,failed);

  y := BatemanGx(1e5);
  f := 0.1000004999999999750e-4;
  testrel(4, NE, y, f, cnt,failed);

  y := BatemanGx(1e4);
  f := 0.1000049999999750000e-3;
  testrel(5, NE, y, f, cnt,failed);

  y := BatemanGx(1000);
  f := 0.1000499999750000500e-2;
  testrel(6, NE, y, f, cnt,failed);

  y := BatemanGx(100);
  f := 0.5024998750249893827e-2*2; {0.10049997500499787655e-1;}
  testrel(7, NE, y, f, cnt,failed);

  y := BatemanGx(20);
  f := 0.5124844523096526763e-1;
  testrel(8, NE2, y, f, cnt,failed);

  y := BatemanGx(19);
  f := 0.5401471266377157448e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := BatemanGx(18);
  f := 0.5709639844733953663e-1;
  testrel(10, NE2, y, f, cnt,failed);

  y := BatemanGx(10);
  f := 0.1049754801499506510;
  testrel(11, NE, y, f, cnt,failed);

  y := BatemanGx(9);
  f := 0.1172467420722715712;
  testrel(12, NE, y, f, cnt,failed);

  y := BatemanGx(8);
  f := 0.1327532579277284288;
  testrel(13, NE, y, f, cnt,failed);

  y := BatemanGx(1);
  f := 1.386294361119890619;
  testrel(14, NE, y, f, cnt,failed);

  y := BatemanGx(2);
  f := 0.6137056388801093812;
  testrel(15, NE, y, f, cnt,failed);

  y := BatemanGx(-1.5);
  f := 5.808259320256459905;
  testrel(16, NE, y, f, cnt,failed);

  y := BatemanGx(-42.125);
  f := -16.44221163230795132;
  testrel(17, NE, y, f, cnt,failed);

  y := BatemanGx(-1234567.875);
  f := 16.41875363763280506;
  testrel(18, NE, y, f, cnt,failed);

  {nearly worst case, too 'big' for small x asymptotic,}
  {uses reflection, 10 double recursion steps, 23 loops}
  y := BatemanGx(-2e-8);
  f := -0.1000000013862943940e9;
  testrel(19, NE2, y, f, cnt,failed);

  y := BatemanGx(sqrt_epsh);  {extended: 0.5^32}
  f := 0.8589934590613705639e10;
  testrel(20, NE, y, f, cnt,failed);

  y := BatemanGx(1e-10);
  f := 0.1999999999861370564e11;
  testrel(21, NE, y, f, cnt,failed);

  y := BatemanGx(-1e-10);
  f := -0.2000000000138629436e11;
  testrel(22, NE, y, f, cnt,failed);

  {Check no overlow}
  y := BatemanGx(1e4000);
  f := 1e-4000;
  testrel(23, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lnbetax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lnbetax');

  y := lnbetax(9,9);
  f := -12.295867644646388427;
  testrel( 1, NE, y, f, cnt,failed);

  y := lnbetax(9,8);
  f := -11.542095842270008275;
  testrel( 2, NE, y, f, cnt,failed);

  y := lnbetax(9,7);
  f := -10.715417269085540343;
  testrel( 3, NE, y, f, cnt,failed);

  y := lnbetax(8,9);
  f := -11.542095842270008275;
  testrel( 4, NE, y, f, cnt,failed);

  y := lnbetax(8,8);
  f := -10.848948661710062966;
  testrel( 5, NE, y, f, cnt,failed);

  y := lnbetax(8,7);
  f := -10.086808609663166205;
  testrel( 6, NE, y, f, cnt,failed);

  y := lnbetax(7,9);
  f := -10.715417269085540343;
  testrel( 7, NE, y, f, cnt,failed);

  y := lnbetax(7,8);
  f := -10.086808609663166205;
  testrel( 8, NE, y, f, cnt,failed);

  y := lnbetax(7,7);
  f := -9.3936614291032208955;
  testrel( 9, NE, y, f, cnt,failed);

  y := lnbetax(1.5,-0.504083);
  f := 1.1422848743102066695;
  testrel(10, NE, y, f, cnt,failed);

  y := lnbetax(1,-4.5);
  f := -1.5040773967762740734;
  testrel(11, NE, y, f, cnt,failed);

  y := lnbetax(6.7,1);
  f := -1.9021075263969203756;
  testrel( 12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_betax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 12;
  NE1 = 20;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','betax');

  y := betax(9,9);
  f := 1/218790;  {Fix311}
  testrel( 1, NE, y, f, cnt,failed);

  y := betax(9,8);
  f := 1/102960;
  testrel( 2, NE, y, f, cnt,failed);

  y := betax(9,7);
  f := 1/45045;
  testrel( 3, NE, y, f, cnt,failed);

  y := betax(8,9);
  f := 1/102960;
  testrel( 4, NE, y, f, cnt,failed);

  y := betax(8,8);
  f := 1/51480;
  testrel( 5, NE, y, f, cnt,failed);

  y := betax(8,7);
  f := 1/24024;
  testrel( 6, NE, y, f, cnt,failed);

  y := betax(7,9);
  f := 1/45045;
  testrel( 7, NE, y, f, cnt,failed);

  y := betax(7,8);
  f := 1/24024;
  testrel( 8, NE, y, f, cnt,failed);

  y := betax(7,7);
  f := 1/12012;
  testrel( 9, NE, y, f, cnt,failed);

  y := betax(1.5,-0.504083);
  f := -3.1339208060001728471;
  testrel(10, NE, y, f, cnt,failed);

  y := betax(1,-5);
  f := -0.2;
  testrel(11, NE, y, f, cnt,failed);

  y := betax(2,1);
  f := 0.5;
  testrel( 12, NE, y, f, cnt,failed);

  {negative x+y}
  f := betax(-4,3);
  y := -1/12;
  testrel( 13, NE, y, f, cnt,failed);

  y := betax(5,-8);
  f := -1/280;
  testrel( 14, NE, y, f, cnt,failed);

  y := betax(1,-8);
  f := -1/8;
  testrel( 15, NE, y, f, cnt,failed);

  y := betax(-7,4);
  f := 1/140;
  testrel( 16, NE, y, f, cnt,failed);

  y := betax(10.5,8.5);
  f := 0.2484217763355545533e-5;
  testrel( 17, NE, y, f, cnt,failed);

  y := betax(-1.5,2);
  f := 4/3;
  testrel( 18, NE, y, f, cnt,failed);

  y := betax(10.25,4.25);
  f := 0.2293445172607700050e-3;
  testrel( 19, NE, y, f, cnt,failed);

  y := betax(100,4.25);
  f := 0.1223606228656168081e-7*2;
  testrel( 20, NE1, y, f, cnt,failed);

  y := betax(1000,4.25);
  f := 0.1463192181497841789e-11;
  testrel( 21, NE1, y, f, cnt,failed);

  y := betax(1000,1e-10);
  f := 0.9999999992515529142e10;
  testrel( 22, NE1, y, f, cnt,failed);

  y := betax(100,-1.5);
  f := 0.2319080638275973289e4;
  testrel( 23, NE, y, f, cnt,failed);

  y := betax(100,-4.25);
  f := -0.3562396916700877082e8;
  testrel( 24, NE1, y, f, cnt,failed);

  y := betax(1e-5,1e-10);
  f := 0.1000009999999998355e11;
  testrel( 25, NE, y, f, cnt,failed);

  y := betax(1e-5,0.5e-20);
  f := 0.2000000000000001e21;
  testrel( 26, NE1, y, f, cnt,failed);

  y := betax(2,0.5e-20);
  f := 2e20;
  testrel( 27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ibetax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 10;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','ibetax');

  y := ibetax(1.2, 1.3, 0);
  f := 0;
  testabs( 1, NE, y, f, cnt,failed);

  y := ibetax(1.2, 1.3, 1e-10);
  f := 1.3443494465428975103e-12;
  testrel( 2, NE, y, f, cnt,failed);

  y := ibetax(1.2, 1.3, 0.001);
  f := 3.3763004250453581348e-4;
  testrel( 3, NE, y, f, cnt,failed);

  y := ibetax(1.2, 1.3, 0.1);
  f := 8.3399782830674834581e-2;
  testrel( 4, NE, y, f, cnt,failed);

  y := ibetax(1.2, 1.3, 0.5);
  f := 5.2978142945129908088e-1;
  testrel( 5, NE, y, f, cnt,failed);

  y := ibetax(1.2, 1.3, 0.9);
  f := 9.3852939722443065921e-1;
  testrel( 6, NE, y, f, cnt,failed);

  y := ibetax(1.2, 1.3, 0.99);
  f := 9.9688643834125438042e-1;
  testrel( 7, NE, y, f, cnt,failed);

  y := ibetax(1.2, 1.3, 0.999);
  f := 9.9984379283306763424e-1;
  testrel( 8, NE, y, f, cnt,failed);

  y := ibetax(1.5, 4.3, 0.9);
  f := 0.99987810769284942730980; {qcalc}
  testrel( 9, NE, y, f, cnt,failed);

  y := ibetax(5, 20, 0.2);
  f := 5.4012267024242087821312e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := ibetax(8, 10, 0.2);
  f := 0.1093431523409920000e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := ibetax(8, 10, 0.8);
  f := 0.99950675026247680000;
  testrel(12, NE, y, f, cnt,failed);

  y := ibetax(1, 0.5, 0.9);
  f := 6.8377223398316206680e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := ibetax(0.1, 0.5, 0.95);
  f := 9.5989434282661978298e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := ibetax(0.1, 0.5, 0.99);
  f := 9.82283678979602092297725e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := ibetax(0.1, 0.5, 1);
  f := 1;
  testrel(16, NE, y, f, cnt,failed);

  y := ibetax(1e-4, 1e4, 1e-4);
  f := 9.99978059362107134278786e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := ibetax(1e-4, 1e4, 1e-3);
  f := 9.99999999586219825547e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := ibetax(1e-4, 1e4, 1e-10);
  f := 9.9867703308020004431e-1;
  testrel(19, NE, y, f, cnt,failed);

  {Extreme cases, reduced accuracy}
  y := ibetax(1.2, 1.3, 1e-100);
  f := 1.34434944656489595584e-120;
  testrel(20, 100, y, f, cnt,failed);

  y := ibetax(800, 0.5, 0.9);
  f := 1.5538647445875693462265e-38;
  testrel(21, 200, y, f, cnt,failed);

  y := ibetax(900, 900, 0.515625);
  f := 9.075743403780043342e-1;
  testrel(22, NE, y, f, cnt,failed);

  {added after changing threshold to 0.95}
  y := ibetax(0.1, 0.5, 0.95);
  f := 9.5989434282661978298e-1;
  testrel(23, NE, y, f, cnt,failed);

  {added after addding basym and ccdf_binom}
  y := ibetax(2400,340,0.875);
  f := 0.4368082445466639340;
  testrel(24, NE, y, f, cnt,failed);

  y := ibetax(4500,40,0.875);
  f := 0.1212321666933000790e-199;
  testrel(25, 12, y, f, cnt,failed);    {!!!!}

  y := ibetax(100,102,0.5);
  f := 0.5560695261428739647;
  testrel(26, NE, y, f, cnt,failed);

  y := ibetax(100000000,100000002,0.5);
  f := 0.5000564189580021571;
  testrel(27, NE, y, f, cnt,failed);

  y := ibetax(100,300,0.25);
  f := 0.5076806528982708500;
  testrel(28, NE, y, f, cnt,failed);

  y := ibetax(1e-23,2e-23,0.99);
  f := 2/3;
  testrel(29, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ibeta_invx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','ibeta_invx');

  y := ibeta_invx(2, 1, 0.1);
  f := 3.162277660168379332e-1;
  testrel( 1, NE, y, f, cnt,failed);

  y := ibeta_invx(1, 2, 0.01);
  f := 5.012562893380045266e-3;
  testrel( 2, NE, y, f, cnt,failed);

  y := ibeta_invx(2, 1, 0.01);
  f := 10e-2;
  testrel( 3, NE, y, f, cnt,failed);

  y := ibeta_invx(2, 1, 0.05);
  f := 2.236067977499789696e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := ibeta_invx(2, 1, 0.99);
  f := 1.989974874213239909*0.5;
  testrel( 5, NE, y, f, cnt,failed);

  y := ibeta_invx(5, 10, 0.2);
  f := 0.22912109609099720490;
  testrel( 6, NE, y, f, cnt,failed);

  y := ibeta_invx(1, 1, 0);
  f := 0;
  testabs( 7, 1, y, f, cnt,failed);

  y := ibeta_invx(1, 1, 1);
  f := 1;
  testabs( 8, 1, y, f, cnt,failed);

  y := ibeta_invx(4, 7, 0.8);
  f := 4.836572374148516672e-1;
  testrel( 9, NE, y, f, cnt,failed);

  y := ibeta_invx(1.5, 4.25, 0.99);
  f := 0.7194720681939346506;
  testrel(10, NE, y, f, cnt,failed);

  y := ibeta_invx(5, 20, 0.5);
  f := 1.919240597664992704e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := ibeta_invx(8, 10, 0.2);
  f := 3.449977502489844803e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := ibeta_invx(8, 10, 0.8);
  f := 5.427825084949223437e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := ibeta_invx(1, 0.5, 0.9);
  f := 0.99;
  testrel(14, NE, y, f, cnt,failed);

  y := ibeta_invx(1800, 0.5, 0.9);
  f := 9.999956130742314815e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := ibeta_invx(0.1, 0.5, 0.98);
  f := 9.872767968062899295e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := ibeta_invx(0.5, 0.125, 0.9);
  f := 9.999999663947204350e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := ibeta_invx(900, 900, 0.9);
  f := 5.151018817560778887e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := ibeta_invx(1800, 0.5, 1e-4);
  f := 9.958036054195801287e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := ibeta_invx(4500,40,0.125);
  f := 0.9895767941483560347;
  testrel(20, NE, y, f, cnt,failed);

  y := ibeta_invx(100,102,0.5);
  f := 0.4950331343781250235;
  testrel(21, NE, y, f, cnt,failed);

  y := ibeta_invx(100,102,0.875);
  f := 0.5355056613331838346;
  testrel(22, NE, y, f, cnt,failed);

  y := ibeta_invx(100,300,0.5);
  f := 0.2495830001077997259;
  testrel(23, NE, y, f, cnt,failed);

  y := ibeta_invx(100,300,0.625);
  f := 0.2565236334554090243;
  testrel(24, NE, y, f, cnt,failed);

  y := ibeta_invx(1000,200,0.125);
  f := 0.8209035891017542587;
  testrel(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_beta3x;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE1  = 8;
  NE2 = 50;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','beta3x');

  {special case}
  y := beta3x(2.5,1.5,1);
  f := 0.1963495408493620774;
  testrel(1, NE, y, f, cnt,failed);

  y := beta3x(1.5,2.5,1);
  f := 0.1963495408493620774;
  testrel(2, NE, y, f, cnt,failed);

  y := beta3x(2.5,1,0.5);
  f := 0.07071067811865475244;
  testrel(3, NE, y, f, cnt,failed);

  y := beta3x(1,0.015625,0.9375);
  f := 2.713390035291286596;
  testrel(4, NE, y, f, cnt,failed);

  y := beta3x(1,2/3,1e-5);
  f := 0.1000001666674074117e-4;
  testrel(5, NE, y, f, cnt,failed);

  {'normal' cases}
  y := beta3x(3.25,3,0.75);
  f := 0.02429951329179160258;
  testrel(6, NE, y, f, cnt,failed);

  y := beta3x(3,3.25,0.25);
  f := 0.003280749582588918860;
  testrel(7, NE, y, f, cnt,failed);

  y := beta3x(50,20,0.25);
  f := 0.7606402961983154568e-34;
  testrel(8, NE, y, f, cnt,failed);

  y := beta3x(1.5,2.5,0.75);
  f := 0.1850262816361021337;
  testrel(9, NE, y, f, cnt,failed);

  y := beta3x(1.5,2.5,0.5);
  f := 0.1398414370913477054;
  testrel(10, NE, y, f, cnt,failed);

  y := beta3x(1/16,2.5,0.75);
  f := 14.78652498220471553;
  testrel(11, NE, y, f, cnt,failed);

  y := beta3x(1/32,0.125,0.875);
  f := 33.51189552677704245;
  testrel(12, NE, y, f, cnt,failed);

  y := beta3x(1e-5,0.25,0.875);
  f := 1.000012078282605376e5;
  testrel(13, NE, y, f, cnt,failed);

  y := beta3x(1e-10,0.25,0.875);
  f := 1.0000000001207838737245e10;
{$ifdef FPC271or3}
  testrel(14, NE+1, y, f, cnt,failed);
{$else}
  testrel(14, NE, y, f, cnt,failed);
{$endif}

  y := beta3x(1/128,0.25,0.875);
  f := 129.1997065384640326;
  testrel(15, NE, y, f, cnt,failed);

  y := beta3x(1/2,2/3,1e-5);
  f := 0.6324562347648557218e-2;
  testrel(16, NE, y, f, cnt,failed);

  y := beta3x(1/2,1/3,1/4);
  f := 1.063845731927502405;
  testrel(17, NE, y, f, cnt,failed);

  y := beta3x(1/2,1e-9,0.75);
  f := 2.633915792528420978;
  testrel(18, NE, y, f, cnt,failed);

  y := beta3x(1/2,1e5,0.75);
  f := 0.005604998222641328069;
  testrel(19, NE, y, f, cnt,failed);

  y := beta3x(1/100,1e8,0.75);
  f := 82.70442197328590292;
  testrel(20, NE, y, f, cnt,failed);

  y := beta3x(2.75, 3.25, 0.625);
  f := 0.2709555897182153629e-1;
  testrel(21, NE, y, f, cnt,failed);

  y := beta3x(8.25, 3.50, 0.875);
  f := 0.1185537392578522183e-2;
  testrel(22, NE1, y, f, cnt,failed);      {NE1 for AMD}

  y := beta3x(1.25, 9.50, 0.125);
  f := 0.3369236106120093626e-1;
  testrel(23, NE, y, f, cnt,failed);

  y := beta3x(3.50, 2.50, 0.625);
  f := 0.2065970176350567637e-1;
  testrel(24, NE, y, f, cnt,failed);

  y := beta3x(4.25, 1.50, 0.750);
  f := 0.4297971957078045810e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := beta3x(8.25, 1.25, 0.875);
  f := 0.2724467366536068364e-1;
  testrel(26, NE, y, f, cnt,failed);

  y := beta3x(1.25, 8.25, 0.125);
  f := 0.3640135574010781261e-1;
  testrel(27, NE, y, f, cnt,failed);

  y := beta3x(2.50, 3.50, 0.375);
  f := 0.1615583714574971314e-1;
  testrel(28, NE, y, f, cnt,failed);

  y := beta3x(2.25, 2.25, 0.500);
  f := 0.5518079396730273567e-1;
  testrel(29, NE, y, f, cnt,failed);

  y := beta3x(2,3.125,0.375);
  f := 0.03878878296398575327;
  testrel(30, NE, y, f, cnt,failed);

  y := beta3x(3,3.25,0.5);
  f := 0.01496500475931052035;
  testrel(31, NE, y, f, cnt,failed);

  {large/extreme parameters}
  y := beta3x(1000,25,0.75);
  f := 4.406883674164145361e-143;
  testrel(32, NE, y, f, cnt,failed);

  y := beta3x(4500,40,0.875);
  f := 1.547547465068642999e-300;
  testrel(33, NE, y, f, cnt,failed);

  {extreme}
  y := beta3x(1e-20,1e-20,0.75);
  f := 1.0e20;
  testrel(34, NE2, y, f, cnt,failed);

  y := beta3x(100,300,0.25);
  f := 3.019408481529601946e-99;
  testrel(35, NE2, y, f, cnt,failed);

  {BASYM}
  y := beta3x(100,102,0.5);
  f := 3.085837704355118753e-62;
  testrel(36, 120, y, f, cnt,failed);   {!!!!!}

  {a<=0 or b<=0}
  y := beta3x(-1.5, 0.25, 0.75);
  f := -1.113130970686713684;
  testrel(37, NE, y, f, cnt,failed);

  y := beta3x(-1.5, 0.5, 0.75);
  f := -1.283000598199168366;
  testrel(38, NE, y, f, cnt,failed);

  y := beta3x(-1.5, 0.5+ldexp(1,-20), 0.75);
  f := -1.283000830999050474;
  testrel(39, NE, y, f, cnt,failed);

  y := beta3x(-1.25, 0.25, 0.75);
  f := -3.241965592113461535;
  testrel(40, NE, y, f, cnt,failed);

  y := beta3x(-1.25, -0.25, 0.75);
  f := -3.776649532917987619;
  testrel(41, NE, y, f, cnt,failed);

  y := beta3x(1.25, -0.25, 0.75);
  f := 1.337062175798913071;
  testrel(42, NE, y, f, cnt,failed);

(* MPMath bug:
  >>> betainc(mpf('-2.00000000000000000000000000'),2,0,0.75)
  mpf('-0.055555555555555555555555555555555555555555555555555555555555573')
  >>> betainc(mpf('-2.000000000000000000000000001'),2,0,0.75)
  mpf('0.44444444444444444444444444368341425442301374552854178131548544')
*)
  y := beta3x(-2, 2, 0.75);
  f := 0.4444444444444444444;
  testrel(43, NE, y, f, cnt,failed);

  y := beta3x(-2, 1, 0.75);
  f := -0.8888888888888888889;
  testrel(44, NE, y, f, cnt,failed);

  y := beta3x(1, -2, 0.75);
  f := 7.5;
  testrel(45, NE, y, f, cnt,failed);

  y := beta3x(-1.5, 0, 0.75);
  f := -0.7018857614682043332;
  testrel(46, NE, y, f, cnt,failed);

  y := beta3x(-1.5, 0.5, 0.75);
  f := -1.2830005981991683656;
  testrel(47, NE, y, f, cnt,failed);

  y := beta3x(-2+ldexp(1,-30), -0.25, 0.75);
  f := 1.509949439278383020e09;
  testrel(48, NE, y, f, cnt,failed);

  {extended}
  y := beta3x(900,700, 1/65536.0);
  f := 0.1618798884281406400e-4337;
  testrel(49, NE, y, f, cnt,failed);

  y := beta3x(5450,25,0.125);
  f := 1.075575511106366976e-4927;
  testrel(50, NE, y, f, cnt,failed);

  y := beta3x(-5454.5, -0.25, 0.125);
  f := -1.738134006073741041e4922;
  testrel(51, NE, y, f, cnt,failed);

  y := beta3x(5454.5, -0.25, 0.125);
  f := 2.700142752066199782e-4930;
  testrel(52, NE, y, f, cnt,failed);

  {extended, extreme, BASYM}
  y := beta3x(8000,8002,0.5);
  f := 1.661565524082334348e-4819;
  testrel(53, 1400, y, f, cnt,failed);   {!!!!!!!!!}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;

end.
