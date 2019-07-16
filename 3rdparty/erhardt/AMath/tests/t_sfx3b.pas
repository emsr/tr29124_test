{Part 3b of regression test for SPECFUNX unit  (c) 2011+  W.Ehrhardt}

unit t_sfx3b;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_laplace_cdfx;
procedure test_laplace_invx;
procedure test_laplace_pdfx;

procedure test_logistic_cdfx;
procedure test_logistic_invx;
procedure test_logistic_pdfx;

procedure test_lognormal_cdfx;
procedure test_lognormal_invx;
procedure test_lognormal_pdfx;

procedure test_moyal_cdfx;
procedure test_moyal_invx;
procedure test_moyal_pdfx;

procedure test_pareto_cdfx;
procedure test_pareto_invx;
procedure test_pareto_pdfx;

procedure test_triangular_pdfx;
procedure test_triangular_invx;
procedure test_triangular_cdfx;

procedure test_uniform_cdfx;
procedure test_uniform_invx;
procedure test_uniform_pdfx;

procedure test_weibull_cdfx;
procedure test_weibull_invx;
procedure test_weibull_pdfx;


implementation

uses
  amath, specfunx, t_sfx0;


{---------------------------------------------------------------------------}
procedure test_laplace_cdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 5;
{$ifdef BIT16}
  NE2 = 7;
{$else}
  NE2 = NE;
{$endif}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','laplace_cdfx');

  {statevalf[cdf,laplace[0,2]](x);}
  a := 0.0;
  b := 2.0;
  y := laplace_cdfx(a, b, -100);
  f := 0.9643749239819588915e-22;
  testrel(1, NE2, y, f, cnt,failed);

  y := laplace_cdfx(a, b, -2);
  f := 0.1839397205857211608;
  testrel(2, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, -1);
  f := 0.3032653298563167118;
  testrel(3, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, -0.5);
  f := 0.3894003915357024341;
  testrel(4, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 0);
  f := 0.5;
  testrel(5, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 0.5);
  f := 0.6105996084642975659;
  testrel(6, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 1);
  f := 0.6967346701436832882;
  testrel(7, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 5);
  f := 0.9589575006880506024;
  testrel(8, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 10);
  f := 0.9966310265004572665;
  testrel(9, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 50);
  f := 0.9999999999930560281;
  testrel(10, NE, y, f, cnt,failed);

  {statevalf[cdf,laplace[1,0.5]](x);}
  a := 1.0;
  b := 0.5;
  y := laplace_cdfx(a, b, -10);
  f := 0.1394734046434462404e-9;
  testrel(11, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, -1);
  f := 0.9157819444367090147e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 0);
  f := 0.6766764161830634595e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 0.5);
  f := 0.1839397205857211608;
  testrel(14, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 0.75);
  f := 0.3032653298563167118;
  testrel(15, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 1);
  f := 0.5;
  testrel(16, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 1.5);
  f := 0.8160602794142788392;
  testrel(17, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 2);
  f := 0.9323323583816936541;
  testrel(18, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 5);
  f := 0.9998322686860487441;
  testrel(19, NE, y, f, cnt,failed);

  y := laplace_cdfx(a, b, 10);
  f := 0.9999999923850101276;
  testrel(20, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_laplace_invx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','laplace_invx');
  a := 0.0;
  b := 2.0;
  y := laplace_invx(a, b, ldexp(1,-20));
  f := -26.33959286127792176;
  testrel(1, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.125);
  f := -2.772588722239781238;
  testrel(2, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.5);
  f := 0.0;
  testrel(3, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.75);
  f := 1.386294361119890619;
  testrel(4, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.875);
  f := 2.772588722239781238;
  testrel(5, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.9);
  f := 3.218875824868200749;
  testrel(6, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.99);
  f := 7.824046010856292117;
  testrel(7, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.9990234375);
  f := 12.47664925007901557;
  testrel(8, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 1-ldexp(1,-20));
  f := 26.33959286127792176;
  testrel(9, NE, y, f, cnt,failed);

  a := 1.0;
  b := 0.5;
  y := laplace_invx(a, b, ldexp(1,-20));
  f := -5.584898215319480439;
  testrel(10, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.125);
  f := 0.3068528194400546906;
  testrel(11, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.5);
  f := 1.0;
  testrel(12, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.75);
  f := 1.346573590279972655;
  testrel(13, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.875);
  f := 1.693147180559945309;
  testrel(14, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.9);
  f := 1.804718956217050187;
  testrel(15, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.99);
  f := 2.956011502714073029;
  testrel(16, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 0.9990234375);
  f := 4.119162312519753892;
  testrel(17, NE, y, f, cnt,failed);

  y := laplace_invx(a, b, 1-ldexp(1,-20));
  f := 7.584898215319480439;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_laplace_pdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 5; {for FPC}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','laplace_pdfx');

  {statevalf[pdf,laplace[0,2]](x);}
  a := 0.0;
  b := 2.0;
  y := laplace_pdfx(a, b, 0);
  f := 0.25;
  testrel(1, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, ldexp(1,-20));
  f := 0.2499998807907388709;
  testrel(2, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, -0.125);
  f := 0.2348532657033689465;
  testrel(3, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 0.5);
  f := 0.1947001957678512171;
  testrel(4, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, -1);
  f := 0.1516326649281583559;
  testrel(5, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 5);
  f := 0.2052124965597469879e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, -10);
  f := 0.1684486749771366774e-2;
  testrel(7, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 100);
  f := 0.4821874619909794458e-22;
  testrel(8, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, -1000);
  f := 0.1781144101685321383e-217;
  testrel(9, NE, y, f, cnt,failed);

  {statevalf[pdf,laplace[1,0.5]](x);}
  a := 1.0;
  b := 0.5;
  y := laplace_pdfx(a, b, -100);
  f := 0.1872900284160809227e-87;
  testrel(10, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, -1);
  f := 0.1831563888873418029e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 0);
  f := 0.1353352832366126919;
  testrel(12, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 0.125);
  f := 0.1737739434504451267;
  testrel(13, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 0.5);
  f := 0.3678794411714423216;
  testrel(14, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 0.9990234375);
  f := 0.998048781107475473;
  testrel(15, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 1);
  f := 1.0;
  testrel(16, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 1.25);
  f := 0.6065306597126334236;
  testrel(17, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 3);
  f := 0.1831563888873418029e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 10);
  f := 0.1522997974471262844e-7;
  testrel(19, NE, y, f, cnt,failed);

  y := laplace_pdfx(a, b, 200);
  f := 0.1415129558908617756e-172;
  testrel(20, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_logistic_cdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','logistic_cdfx');

  {statevalf[cdf,logistic[a,b]](x);}
  a := 0;
  b := 1;
  y := logistic_cdfx(a, b, -100);
  f := 0.3720075976020835963e-43;
  testrel(1, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, -10);
  f := 0.4539786870243439450e-4;
  testrel(2, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, -1);
  f := 0.2689414213699951207;
  testrel(3, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 0);
  f := 0.5;
  testrel(4, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 1);
  f := 0.7310585786300048793;
  testrel(5, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 2);
  f := 0.8807970779778824441;
  testrel(6, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 5);
  f := 0.9933071490757151444;
  testrel(7, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 10);
  f := 0.9999546021312975656;
  testrel(8, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 20);
  f := 0.9999999979388463818;
  testrel(9, NE2, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 50);
  f := 1.0;
  testrel(10, NE, y, f, cnt,failed);

  a := 2.0;
  b := 0.5;
  y := logistic_cdfx(a, b, -100);
  f := 0.2534694904308355121e-88;
  testrel(11, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 0);
  f := 0.1798620996209155803e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 1);
  f := 0.1192029220221175559;
  testrel(13, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 2);
  f := 0.5;
  testrel(1, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 3);
  f := 0.8807970779778824441;
  testrel(14, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 4);
  f := 0.9820137900379084420;
  testrel(15, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 5);
  f := 0.9975273768433652257;
  testrel(16, NE2, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 10);
  f := 0.9999998874648379449;
  testrel(17, NE2, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 20);
  f := 0.9999999999999997680;
  testrel(18, NE, y, f, cnt,failed);

  y := logistic_cdfx(a, b, 25);
  f := 1.0;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_logistic_invx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','logistic_invx');

  a := 0.0;
  b := 1.0;
  y := logistic_invx(a, b, ldexp(1,-20));
  f := -13.86294265752413503;
  testrel(1, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.125);
  f := -1.945910149055313305;
  testrel(2, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.5);
  f := 0;
  testrel(3, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.75);
  f := 1.098612288668109691;
  testrel(4, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.875);
  f := 1.945910149055313305;
  testrel(5, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.9);
  f := 2.197224577336219383;
  testrel(6, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.99);
  f := 4.595119850134589927;
  testrel(7, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.9990234375);
  f := 6.930494765951626481;
  testrel(8, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 1-ldexp(1,-20));
  f := 13.86294265752413503;
  testrel(9, NE, y, f, cnt,failed);

  a := 2.0;
  b := 0.5;
  y := logistic_invx(a, b, ldexp(1,-20));
  f := -4.931471328762067517;
  testrel(10, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.125);
  f := 1.027044925472343347;
  testrel(11, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.5);
  f := 2.0;
  testrel(12, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.75);
  f := 2.549306144334054846;
  testrel(13, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.875);
  f := 2.972955074527656653;
  testrel(14, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.9);
  f := 3.098612288668109691;
  testrel(15, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.99);
  f := 4.297559925067294963;
  testrel(16, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 0.9990234375);
  f := 5.465247382975813241;
  testrel(17, NE, y, f, cnt,failed);

  y := logistic_invx(a, b, 1-ldexp(1,-20));
  f := 8.931471328762067517;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_logistic_pdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 3;
  NE2 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','logistic_pdfx');

  a := 0.0;
  b := 1.0;
  y := logistic_pdfx(a, b, 0);
  f := 0.25;
  testrel(1, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 0.5);
  f := 0.2350037122015944891;
  testrel(2, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, -2);
  f := 0.1049935854035065174;
  testrel(3, NE2, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 5);
  f := 0.6648056670790154914e-2;
  testrel(4, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, -10);
  f := 0.4539580773595167103e-4;
  testrel(5, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 50);
  f := 0.1928749847963917783e-21;
  testrel(6, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, -100);
  f := 0.3720075976020835963e-43;
  testrel(7, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 500);
  f := 0.7124576406741285532e-217;
  testrel(8, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 10000);
  f := 0.1135483865314736099e-4342;
  testrel(9, NE2, y, f, cnt,failed);

  a := 2.0;
  b := 0.5;
  y := logistic_pdfx(a, b, 0);
  f := 0.3532541242658223284e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 0.5);
  f := 0.9035331946182426530e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, -1);
  f := 0.4933018582720095637e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 2);
  f := 0.5;
  testrel(13, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, -5);
  f := 0.1663054672450542702e-5;
  testrel(14, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 5);
  f := 0.4933018582720095637e-2;
  testrel(15, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, -10);
  f := 0.7550269087988129870e-10;
  testrel(16, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 50);
  f := 0.4062185325469621851e-41;
  testrel(17, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 300);
  f := 0.2894134972965153792e-258;
  testrel(18, NE, y, f, cnt,failed);

  y := logistic_pdfx(a, b, 5000);
  {f := 0.1239906368773370148e-4340;}
  f := 0.12399063687733701483976265325e-4340; {FPC271}
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lognormal_cdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 10;
  NE2 = 32;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lognormal_cdfx');
  {statevalf[cdf,lognormal[a,b]](x);}

  a := 2;
  b := 1.5;
  y := lognormal_cdfx(a, b, 0.0009765625);
  f := 0.1305820500266527804e-8;
  testrel(1, NE2, y, f, cnt,failed);  {!!!!}

  y := lognormal_cdfx(a, b, 0.125);
  f := 0.3267772762571824551e-2;
  testrel(2, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 1);
  f := 0.9121121972586786984e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 5);
  f := 0.3972873686243666920;
  testrel(4, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 10);
  f := 0.5799335141280535526;
  testrel(4, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 50);
  f := 0.8987890904238633564;
  testrel(6, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 100);
  f := 0.9587870056652280963;
  testrel(7, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 1000);
  f := 0.9994657439207799539;
  testrel(8, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 10000);
  f := 0.9999992335323244909;
  testrel(9, NE, y, f, cnt,failed);


  a := 0.25;
  b := 0.5;
  y := lognormal_cdfx(a, b, 0.0009765625);
  f := 0.4419579950671215973e-46;
  testrel(10, NE2, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 0.125);
  f := 0.1589648510521725010e-5;
  testrel(11, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 0.5);
  f := 0.2962764936665734492e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 1);
  f := 0.3085375387259868964;
  testrel(13, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 2);
  f := 0.8122705365749475296;
  testrel(14, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 5);
  f := 0.9967247902977106223;
  testrel(15, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 10);
  f := 0.9999797991334357154;
  testrel(16, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 50);
  f := 0.9999999999998796979;
  testrel(17, NE, y, f, cnt,failed);

  y := lognormal_cdfx(a, b, 150);
  f := 1.0;
  testrel(18, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lognormal_invx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lognormal_invx');

  a := 2;
  b := 1.5;

  y := lognormal_invx(a, b, ldexp(1,-20));
  f := 0.5831380283049865408e-2;
  testrel(1, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.125);
  f := 1.315840900655271756;
  testrel(2, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.5);
  f := 7.389056098930650227;
  testrel(3, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.75);
  f := 20.32262150159026255;
  testrel(4, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.875);
  f := 41.49297229319674307;
  testrel(5, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.9);
  f := 50.51788077261473522;
  testrel(6, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.99);
  f := 242.141389834721616;
  testrel(7, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.9990234375);
  f := 769.6251801436839498;
  testrel(8, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 1-ldexp(1,-20));
  f := 9362.817614869889005;
  testrel(9, NE, y, f, cnt,failed);

  a := 0.25;
  b := 0.5;
  y := lognormal_invx(a, b, ldexp(1,-20));
  f := 0.1186591101741005831;
  testrel(10, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.125);
  f := 0.7224011462267458664;
  testrel(11, NE, y, f, cnt,failed);


  y := lognormal_invx(a, b, 0.5);
  f := 1.284025416687741484; {Wolfram alpha: InverseCDF[logNormalDistribution[0.25, 0.5], 0.5]}
  testrel(12, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.75);
  f := 1.799025042487527964;
  testrel(13, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.875);
  f := 2.282279422328921280;
  testrel(14, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.9);
  f := 2.437019515889078871;
  testrel(15, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.99);
  f := 4.108976361480573370;
  testrel(16, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 0.9990234375);
  f := 6.041392544505611605;
  testrel(17, NE, y, f, cnt,failed);

  y := lognormal_invx(a, b, 1-ldexp(1,-20));
  f := 13.89460335814982601;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lognormal_pdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 32;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lognormal_pdfx');

  a := 2.0;
  b := 1.5;
  y := lognormal_pdfx(a, b, 0.0009765625);
  f := 0.5450176094613075364e-5;
  testrel(1, NE2, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 0.125);
  f := 0.5269949016070822002e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 1);
  f := 0.1093400497839957465;
  testrel(3, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 2);
  f := 0.9098358075364976163e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 5);
  f := 0.5141943565196794751e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 10);
  f := 0.2606049016368246096e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 100);
  f := 0.5885925220711086559e-3;
  testrel(7, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 1000);
  f := 0.1259724999118361604e-5;
  testrel(8, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 1e9);
  f := 0.3910637796213916798e-43;
  testrel(9, NE2, y, f, cnt,failed);

  a := 0.25;
  b := 0.5;
  y := lognormal_pdfx(a, b, 0.0009765625);
  f := 0.13062752536277960258e-41;
  testrel(10, NE2, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 0.125);
  f := 0.1235399956173157858e-3;
  testrel(11, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 0.25);
  f := 0.1507955442895558942e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 0.5);
  f := 0.2693624575570867716;
  testrel(13, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 1);
  f := 0.7041306535285989555;
  testrel(14, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 2);
  f := 0.2693624575570867716;
  testrel(15, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 5);
  f := 0.3960550926237422989e-2;
  testrel(16, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 10);
  f := 0.1747765109291960413e-4;
  testrel(17, NE, y, f, cnt,failed);

  y := lognormal_pdfx(a, b, 100);
  f := 0.2672837201236921963e-18;
  testrel(18, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_pareto_pdfx;
var
  y,f,a,k: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','pareto_pdfx');

  {pareto_pdf := (k,a,x) -> a*k^a/x^(a+1);}
  k := 0.8;
  a := 0.5;
  y := pareto_pdfx(k, a, 0.8);
  f := 0.625;
  testrel(1, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 1);
  f := 0.4472135954999579393;
  testrel(2, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 2);
  f := 0.1581138830084189666;
  testrel(3, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 5);
  f := 0.04;
  testrel(4, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 10);
  f := 0.1414213562373095049e-1;
  {WA:0.014142135623730950488016887242096980785697}  {WA: Wolfram Alpha}
  testrel(5, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 100);
  f := 0.4472135954999579393e-3;
  testrel(6, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 10000);
  f := 0.4472135954999579393e-6;
  testrel(7, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 1e9);
  f := 0.1414213562373095049e-13;
  testrel(8, NE, y, f, cnt,failed);

  k := 2;
  a := 2.5;
  y := pareto_pdfx(k, a, 2);
  f := 1.25;
  testrel(9, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 2.5);
  f := 0.5724334022399461623;
  testrel(10, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 3);
  f := 0.3024061410843429751;
  testrel(11, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 5);
  f := 0.5059644256269406931e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 10);
  f := 0.4472135954999579393e-2;
 {WA:0.0044721359549995793928183473374625524708812}
  testrel(13, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 50);
  f := 1.6e-5;
  testrel(14, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 500);
  f := 0.5059644256269406931e-8;
  testrel(15, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 10000);
  f := 0.1414213562373095049e-12;
  testrel(16, NE, y, f, cnt,failed);

  y := pareto_pdfx(k, a, 1e9);
  f := 0.4472135954999579393e-30;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_pareto_invx;
var
  y,f,a,k: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','pareto_invx');

  {pareto_invx := (k,a,x) -> 1-(k/x)^a}
  k := 0.8;
  a := 0.5;
  y := pareto_invx(k, a, 0);
  f := 0.8;
  testrel(1, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, ldexp(1,-20));
  f := 0.80000152588108904006;
  testrel(2, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.125);
  f := 1.044897959183673469;
  {WA: 1.0448979591836734693877551020408163265306}
  testrel(3, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.5);
  f := 3.2;
  testrel(4, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.75);
  f := 12.8;
  testrel(5, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.875);
  f := 51.2;
  testrel(6, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.9);
  f := 80.0;
  testrel(7, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.99);           {not exact}
  f := 8000.0;
  testrel(8, 16, y, f, cnt,failed);       {!!!!!!!!!!!}

  y := pareto_invx(k, a, 0.9990234375);
  f := 838860.8;
  testrel(9, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 1-ldexp(1,-20));
  f := 879609302220.8;
  testrel(10, NE, y, f, cnt,failed);

  k := 2.0;
  a := 2.5;
  y := pareto_invx(k, a, 0);
  f := 2.0;
  testrel(11, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, ldexp(1,-20));
  f := 2.000000762939962442;
  testrel(12, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.125);
  f := 2.109729494498063450;
  {WA  2.1097294944980634503461707435061465224367}
  testrel(13, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.5);
  f := 2.639015821545788519;
  testrel(14, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.75);
  f := 3.482202253184496557;
  testrel(14, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.875);
  f := 4.594793419988140027;
  testrel(16, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.9);
  f := 5.023772863019160222;
  testrel(17, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.99);
  f := 12.61914688960386499;
  testrel(18, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 0.9990234375);
  f := 32.0;
  testrel(19, NE, y, f, cnt,failed);

  y := pareto_invx(k, a, 1-ldexp(1,-20));
  f := 512.0;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_pareto_cdfx;
var
  y,f,a,k: extended;
  cnt, failed: integer;
const
  NE = 4;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','pareto_cdfx');

  {pareto_cdfx := (k,a,x) -> 1-(k/x)^a}
  k := 0.8;
  a := 0.5;
  y := pareto_cdfx(k, a, 0.80078125);
  f := 0.4879239129211806457e-3;
  testrel(1, NE2, y, f, cnt,failed);     {!!!!}

  y := pareto_cdfx(k, a, 1);
  f := 0.1055728090000841214;
  testrel(2, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 2);
  f := 0.3675444679663241336;
  testrel(3, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 5);
  f := 0.6;
  testrel(4, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 10);
  f := 0.7171572875253809902;
  {WA  0.71715728752538099023966225515806038428607}
  testrel(5, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 50);
  f := 0.8735088935932648267;
  testrel(6, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 100);
  f := 0.9105572809000084121;
  testrel(7, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 1000);
  f := 0.9717157287525380990;
  testrel(8, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 50000);
  f := 0.996;
  testrel(9, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 1e12);
  f := 0.9999991055728090001;
  testrel(10, NE, y, f, cnt,failed);

  k := 2;
  a := 2.5;
  y := pareto_cdfx(k, a, 2.1);
  f := 0.1148298658063191116;
  testrel(11, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 2.25);
  f := 0.2550644609721968467;
  testrel(12, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 2.5);
  f := 0.4275665977600538377;
  testrel(13, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 3);
  f := 0.6371126306987884299;
  testrel(14, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 5);
  f := 0.8988071148746118614;
  testrel(15, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 10);
  f := 0.9821114561800016824;
  {WA  0.98211145618000168242872661065014979011648}
  testrel(16, NE2, y, f, cnt,failed);  {!!!!}

  y := pareto_cdfx(k, a, 100);
  f := 0.9999434314575050762;
  testrel(17, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 1000);
  f := 0.9999998211145618000;
  testrel(18, NE, y, f, cnt,failed);

  y := pareto_cdfx(k, a, 1e6);
  f := 0.9999999999999943431;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_uniform_pdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','uniform_pdfx');

  a := -1.0;
  b := 3.0;

  y := uniform_pdfx(a, b, -2);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := uniform_pdfx(a, b, -1);
  f := 0.25;
  testrel(2, NE, y, f, cnt,failed);

  y := uniform_pdfx(a, b, 1);
  f := 0.25;
  testrel(3, NE, y, f, cnt,failed);

  y := uniform_pdfx(a, b, 3);
  f := 0.25;
  testrel(4, NE, y, f, cnt,failed);

  y := uniform_pdfx(a, b, 4);
  f := 0.0;
  testrel(55, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_uniform_invx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','uniform_invx');

  a := -1.0;
  b := 3.0;
  y := uniform_invx(a, b, ldexp(1,-20));
  f := -0.999996185302734375;
  testrel(1, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 0.125);
  f := -0.5;
  testrel(2, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 0.5);
  f := 1.0;
  testrel(3, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 0.75);
  f := 2.0;
  testrel(4, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 0.875);
  f := 2.5;
  testrel(5, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 0.9);
  f := 2.6;
  testrel(6, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 0.99);
  f := 2.96;
  testrel(7, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 0.9990234375);
  f := 2.99609375;
  testrel(8, NE, y, f, cnt,failed);

  y := uniform_invx(a, b, 1-ldexp(1,-20));
  f := 2.999996185302734375;
  testrel(9, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_triangular_cdfx;
var
  y,f,a,b,c: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','triangular_cdfx');

  a := 0.0;
  b := 1.0;
  c := 0.25;

  y := triangular_cdfx(a,b,c,-1);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0);
  f := 0;
  testrel(2, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.1);
  f := 0.04;
  testrel(3, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.25);
  f := 0.25;
  testrel(4, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.5);
  f := 2/3;
  testrel(5, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.8);
  f := 71/75;
  testrel(6, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,1);
  f := 1;
  testrel(7, NE, y, f, cnt,failed);

  a := 0.0;
  b := 1.0;
  c := 1.0;

  y := triangular_cdfx(a,b,c,0);
  f := 0;
  testrel(8, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.25);
  f := 0.0625;
  testrel(9, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.5);
  f := 0.25;
  testrel(10, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.75);
  f := 0.5625;
  testrel(11, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.9);
  f := 0.81;
  testrel(12, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,1.25);
  f := 1;
  testrel(13, NE, y, f, cnt,failed);


  a := -2;
  b := 2;
  c := -1;

  y := triangular_cdfx(a,b,c,-2);
  f := 0;
  testrel(14, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,-1.5);
  f := 0.0625;
  testrel(15, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,-1);
  f := 0.25;
  testrel(16, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,-0.5);
  f := 23/48;
  testrel(17, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0);
  f := 2/3;
  testrel(18, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,0.5);
  f := 0.8125;
  testrel(19, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,1);
  f := 11/12;
  testrel(20, NE, y, f, cnt,failed);

  y := triangular_cdfx(a,b,c,2);
  f := 1;
  testrel(21, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_triangular_pdfx;
var
  y,f,a,b,c: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','triangular_pdfx');

  a := 0.0;
  b := 1.0;
  c := 0.25;

  y := triangular_pdfx(a,b,c,-1);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0);
  f := 0;
  testrel(2, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.1);
  f := 0.8;
  testrel(3, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.25);
  f := 2.0;
  testrel(4, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.5);
  f := 4/3;
  testrel(5, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.8);
  f := 8/15;
  testrel(6, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,1);
  f := 0;
  testrel(7, NE, y, f, cnt,failed);

  a := 0.0;
  b := 1.0;
  c := 1.0;

  y := triangular_pdfx(a,b,c,0);
  f := 0;
  testrel(8, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.25);
  f := 0.5;
  testrel(9, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.5);
  f := 1.0;
  testrel(10, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.75);
  f := 1.5;
  testrel(11, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.9);
  f := 1.8;
  testrel(12, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,1.25);
  f := 0;
  testrel(13, NE, y, f, cnt,failed);


  a := -2;
  b := 2;
  c := -1;

  y := triangular_pdfx(a,b,c,-2);
  f := 0;
  testrel(14, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,-1.5);
  f := 0.25;
  testrel(15, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,-1);
  f := 0.5;
  testrel(16, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,-0.5);
  f := 5/12;
  testrel(17, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0);
  f := 1/3;
  testrel(18, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,0.5);
  f := 0.25;
  testrel(19, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,1);
  f := 1/6;
  testrel(20, NE, y, f, cnt,failed);

  y := triangular_pdfx(a,b,c,2);
  f := 0;
  testrel(21, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_triangular_invx;
var
  y,f,a,b,c: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','triangular_invx');

  a := 0.0;
  b := 1.0;
  c := 0.25;

  y := triangular_invx(a,b,c,0.01);
  f := 0.05;
  testrel(1, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.1);
  f := 0.1581138830084189666;
  testrel(2, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,1/3);
  f := 0.2928932188134524756;
  testrel(3, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.5);
  f := 0.3876275643042054755;
  testrel(4, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,2/3);
  f := 0.5;
  testrel(5, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.9);
  f := 0.7261387212474169433;
  testrel(6, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.99);
  f := 0.9133974596215561353;
  testrel(7, NE, y, f, cnt,failed);

  a := 0.0;
  b := 1.0;
  c := 1.0;

  y := triangular_invx(a,b,c,0.01);
  f := 0.1;
  testrel(8, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.1);
  f := sqrt(0.1);
  testrel(9, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.5);
  f := sqrt(0.5);
  testrel(10, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,2/3);
  f := sqrt(2/3);
  testrel(11, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.9);
  f := sqrt(0.9);
  testrel(12, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.99);
  f := sqrt(0.99);
  testrel(13, NE, y, f, cnt,failed);

  a := -2;
  b := 2;
  c := -1;

  y := triangular_invx(a,b,c,0.01);
  f := -1.8;
  testrel(14, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.1);
  f := -1.367544467966324134;
  testrel(15, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,1/3);
  f := -0.8284271247461900976;
  testrel(16, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.5);
  f := -0.4494897427831780982;
  testrel(17, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,2/3);
  f := 0;
  testrel(18, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.9);
  f := 0.9045548849896677731;
  testrel(19, NE, y, f, cnt,failed);

  y := triangular_invx(a,b,c,0.99);
  f := 1.653589838486224541;
  testrel(20, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_uniform_cdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','uniform_cdfx');

  a := -1.0;
  b := 3.0;
  y := uniform_cdfx(a, b, -2);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := uniform_cdfx(a, b, -1);
  f := 0.0;
  testrel(2, NE, y, f, cnt,failed);

  y := uniform_cdfx(a, b, 0);
  f := 0.25;
  testrel(3, NE, y, f, cnt,failed);

  y := uniform_cdfx(a, b, 0.0009765625);
  f := 0.250244140625;
  testrel(4, NE, y, f, cnt,failed);

  y := uniform_cdfx(a, b, 1);
  f := 0.5;
  testrel(5, NE, y, f, cnt,failed);

  y := uniform_cdfx(a, b, 2.5);
  f := 0.875;
  testrel(6, NE, y, f, cnt,failed);

  y := uniform_cdfx(a, b, 2.99);
  f := 0.9975;
  testrel(7, NE, y, f, cnt,failed);

  y := uniform_cdfx(a, b, 3);
  f := 1.0;
  testrel(8, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_weibull_cdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','weibull_cdfx');

  {statevalf[cdf,weibull[0.5,1]](x);}
  a := 0.5;
  b := 1.0;
  y := weibull_cdfx(a, b, ldexp(1,-20));
  f := 0.9760858180243377653e-3;
  testrel(1, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 0.015625);
  f := 0.1175030974154045971;
  testrel(2, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 0.125);
  f := 0.2978114986734404038;
  testrel(3, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 0.5);
  f := 0.5069313086047602122;
  testrel(4, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 1);
  f := 0.6321205588285576784;
  testrel(5, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 5);
  f := 0.8931220743396142490;
  testrel(6, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 10);
  f := 0.9576707803767950023;
  testrel(7, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 100);
  f := 0.9999546000702375151;
  testrel(8, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 1000);
  f := 0.9999999999999815327;
  testrel(9, 8, y, f, cnt,failed);   {!!!!}


  {statevalf[cdf,weibull[3,2]](x);}
  a := 3.0;
  b := 2.0;
  y := weibull_cdfx(a, b, ldexp(1,-20));
  f := 0.1084202172485504434e-18;
  testrel(10, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 0.125);
  f := 0.2441108251027834869e-3;
  testrel(11, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 0.5);
  f := 0.1550356299459159401e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 1);
  f := 0.1175030974154045971;
  testrel(13, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 2);
  f := 0.6321205588285576784;
  testrel(14, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 3);
  f := 0.9657818816883339653;
  testrel(15, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 4);
  f := 0.9996645373720974882;
  testrel(16, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 6);
  f := 0.9999999999981204712;
  testrel(17, NE, y, f, cnt,failed);

  y := weibull_cdfx(a, b, 8);
  f := 1.0;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_weibull_invx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE  = 3;
  NE2 = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','weibull_invx');

  a := 0.5;
  b := 1.0;
  y := weibull_invx(a, b, ldexp(1,-20));
  f := 0.9094955691354244759e-12;
  testrel(1, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.125);
  f := 0.1783063281624441480e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.5);
  f := 0.48045301391820142467;
  testrel(3, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.75);
  f := 1.921812055672805699;
  testrel(4, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.875);
  f := 4.324077125263812822;
  testrel(5, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.9);
  f := 5.301898110478398011;
  testrel(6, NE2, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.99);
  f := 21.20759244191359204;
  testrel(7, NE2, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.9990234375);
  f := 48.04530139182014247;
  testrel(8, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 1-ldexp(1,-20));
  f := 192.1812055672805699;
  testrel(9, NE, y, f, cnt,failed);

  a := 3.0;
  b := 2.0;
  y := weibull_invx(a, b, ldexp(1,-20));
  f := 0.1968626953365666125e-1;
  testrel(10, NE2, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.125);
  f := 1.022251575238047668;
  testrel(11, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.5);
  f := 1.769994089001035437;
  testrel(12, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.75);
  f := 2.230052810921904143;
  testrel(13, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.875);
  f := 2.552773214308396109;
  testrel(14, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.9);
  f := 2.641000956907370444;
  testrel(15, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.99);
  f := 3.327452698400098850;
  testrel(16, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 0.9990234375);
  f := 3.813336666495212192;
  testrel(17, NE, y, f, cnt,failed);

  y := weibull_invx(a, b, 1-ldexp(1,-20));
  f := 4.804503136453263549;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_weibull_pdfx;
var
  y,f,a,b: extended;
  cnt, failed: integer;
const
  NE = 2;
  NE2 = 3; {for FPC}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','weibull_pdfx');

  a := 0.5;
  b := 1.0;
  y := weibull_pdfx(a, b, ldexp(1,-20));
  f := 511.50024406117153906;
  testrel(1, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 0.125);
  f := 0.9930445019184586250;
  testrel(2, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 0.5);
  f := 0.3486522152763511488;
  testrel(3, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 1);
  f := 0.1839397205857211608;
  testrel(4, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 5);
  f := 0.2389863070707916413e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 10);
  f := 0.6692837279341107373e-2;
  testrel(6, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 100);
  f := 0.2269996488124242577e-5;
  testrel(7, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 10000);
  f := 0.1860037988010417981e-45;
  testrel(8, NE2, y, f, cnt,failed);

  {statevalf[pdf,weibull[3,2]](x);}
  a := 3.0;
  b := 2.0;
  y := weibull_pdfx(a, b, ldexp(1,-20));
  f := 0.3410605131648480892e-12;
  testrel(9, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 0.125);
  f := 0.5857944663134163378e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 0.5);
  f := 0.9229654096925703806e-1;
  testrel(11, 8, y, f, cnt,failed);    {!!!!!!!}

  y := weibull_pdfx(a, b, 1);
  f := 0.3309363384692232761;
  testrel(12, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 2);
  f := 0.5518191617571634824;
  testrel(13, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 3);
  f := 0.1154861493018728671;
  testrel(14, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 4);
  f := 0.2012775767415071033e-2;
  testrel(15, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 5);
  f := 0.1535041059928886882e-5;
  testrel(16, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 6);
  f := 0.2537363902327762448e-10;
  testrel(17, NE, y, f, cnt,failed);

  y := weibull_pdfx(a, b, 10);
  f := 0.1937407737314197868e-52;
  testrel(18, 3, y, f, cnt,failed);     {!!!!!!}

  y := weibull_pdfx(a, b, 15);
  f := 0.5107745300725181818e-181;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_moyal_pdfx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','moyal_pdfx');
  {moyal_pdf := (a,b,x) -> exp(-((x-a)/b + exp(-(x-a)/b))/2)/(b*sqrt(2*Pi));}

  y := moyal_pdfx(0, 2, -8);
  f := 0.2054146402188781323e-11;
  testrel(1, NE2, y, f, cnt,failed);

  y := moyal_pdfx(0, 2, -4);
  f := 0.1347911587940801694e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := moyal_pdfx(0, 2, 0);
  f := 0.1209853622595716749;
  testrel(3, NE, y, f, cnt,failed);

  y := moyal_pdfx(0, 2, 2);
  f := 0.1006581220324439644;
  testrel(4, NE, y, f, cnt,failed);

  y := moyal_pdfx(0, 2, 8);
  f := 0.02674939204444159961; {alpha}
  testrel(5, NE, y, f, cnt,failed);


  y := moyal_pdfx(0, 2, 10);
  f := 0.1631851889962222587e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := moyal_pdfx(0, 2, 100);
  f := 0.2770243997787916495e-11;
  testrel(6, NE, y, f, cnt,failed);

  y := moyal_pdfx(1, 1, -2);
  f := 0.7777371830827616734e-4;
  testrel(8, NE, y, f, cnt,failed);

  y := moyal_pdfx(1, 1, -1);
  f := 0.2695823175881603388e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := moyal_pdfx(1, 1, -0.5);
  f := 0.08983478046364198945; {alpha}
  testrel(10, NE, y, f, cnt,failed);

  y := moyal_pdfx(1, 1, 0);
  f := 0.1689623369069972882;
  testrel(11, NE, y, f, cnt,failed);

  y := moyal_pdfx(1, 1, 2);
  f := 0.2013162440648879288;
  testrel(12, NE, y, f, cnt,failed);

  y := moyal_pdfx(1, 1, 10);
  f := 0.4431574953602824269e-2;
  testrel(13, NE, y, f, cnt,failed);

  y := moyal_pdfx(1, 1, 100);
  f := 0.1268624842535086870e-21;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_moyal_cdfx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 6;
  NE3 = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','moyal_cdfx');
  {moyal_cdf := (a,b,x) -> erfc(exp(-(x-a)/(2*b))/sqrt(2));}

  y := moyal_cdfx(0, 2, -6);
  f := 0.74054583906901796616e-5;
  testrel(1, NE3, y, f, cnt,failed);

  y := moyal_cdfx(0, 2, -4);
  f := 0.6562191672591341096e-2;
  testrel(2, NE2, y, f, cnt,failed);

  y := moyal_cdfx(0, 2, 0);
  f := 0.3173105078629141028;
  testrel(3, NE, y, f, cnt,failed);

  y := moyal_cdfx(0, 2, 2);
  f := 0.544162429362303046517801261567;
  testrel(4, NE, y, f, cnt,failed);

  y := moyal_cdfx(0, 2, 8);
  f := 0.8923467896957689980;  {alpha}
  testrel(5, NE, y, f, cnt,failed);

  y := moyal_cdfx(0, 2, 10);
  f := 0.9345791222280399186;
  testrel(6, NE3, y, f, cnt,failed);

  y := moyal_cdfx(0, 2, 100);
  f := 0.4999999999944595120*2;
  testrel(7, NE, y, f, cnt,failed);

  y := moyal_cdfx(1, 1, -1.5);
  f := 0.4824010414200695478e-3;
  testrel(8, NE, y, f, cnt,failed);

  y := moyal_cdfx(1, 1, -1);
  f := 0.6562191672591341096e-2;
  testrel(9, NE2, y, f, cnt,failed);

  y := moyal_cdfx(1, 1, 0);
  f := 0.4960237520555734989e-1*2;
  testrel(10, NE, y, f, cnt,failed);

  y := moyal_cdfx(1, 1, 2);
  f := 0.5441624293623030465;
  testrel(11, NE, y, f, cnt,failed);

  y := moyal_cdfx(1, 1, 4);
  f := 0.8234342056094611968;  {alpha}
  testrel(12, NE, y, f, cnt,failed);

  y := moyal_cdfx(1, 1, 10);
  f := 0.4955682427419653034*2;
  testrel(13, NE, y, f, cnt,failed);

  y := moyal_cdfx(1, 1, 50);
  f := 0.4999999999908652796*2;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_moyal_invx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','moyal_invx');
  {moyal_inv := (a,b,x) -> a-ln((statevalf[icdf,normald[0,1]](x/2))^2)*b;}

  y := moyal_invx(0, 2, 0.125);
  f := -1.711829126295260056;
  testrel(1, NE, y, f, cnt,failed);

  y := moyal_invx(0, 2, 0.25);
  f := -0.5602628210286337878;
  testrel(2, NE, y, f, cnt,failed);

  y := moyal_invx(0, 2, 0.5);
  f := 1.575195198403564366;
  testrel(3, NE, y, f, cnt,failed);

  y := moyal_invx(0, 2, 0.75);
  f := 4.574781344821882738;
  testrel(4, NE, y, f, cnt,failed);

  y := moyal_invx(0, 2, 0.875);
  f := 7.398130184703767302; {alpha}
  testrel(5, NE, y, f, cnt,failed);

  y := moyal_invx(0, 2, 0.95);
  f := 11.07714294407235304;
  testrel(6, NE, y, f, cnt,failed);

  y := moyal_invx(0, 2, 0.99);
  f := 17.51741060923158137;
  testrel(7, NE, y, f, cnt,failed);

  y := moyal_invx(1, 1, 0.125);
  f := 0.1440854368523699718;
  testrel(8, NE, y, f, cnt,failed);

  y := moyal_invx(1, 1, 0.25);
  f := 0.7198685894856831061;
  testrel(9, NE, y, f, cnt,failed);

  y := moyal_invx(1, 1, 0.5);
  f := 1.787597599201782183;
  testrel(10, NE, y, f, cnt,failed);

  y := moyal_invx(1, 1, 0.75);
  f := 3.287390672410941369;
  testrel(11, NE, y, f, cnt,failed);

  y := moyal_invx(1, 1, 0.875);
  f := 4.699065092351883651;
  testrel(12, NE, y, f, cnt,failed);

  y := moyal_invx(1, 1, 0.95);
  f := 6.538571472036176521;
  testrel(13, NE, y, f, cnt,failed);

  y := moyal_invx(1, 1, 0.99);
  f := 2*4.879352652307895341; {alpha}
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.
