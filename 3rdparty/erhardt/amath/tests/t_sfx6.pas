{Part 6 of regression test for SPECFUNX unit  (c) 2010  W.Ehrhardt}

unit t_sfx6;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_bessel_j0x;
procedure test_bessel_j1x;
procedure test_bessel_y0x;
procedure test_bessel_y1x;
procedure test_jnynx;
procedure test_bessel_i0x;
procedure test_bessel_i0ex;
procedure test_bessel_i1x;
procedure test_bessel_i1ex;
procedure test_bessel_k0x;
procedure test_bessel_k0ex;
procedure test_bessel_k1x;
procedure test_bessel_k1ex;
procedure test_inknx;

procedure test_bessel_jvx;
procedure test_bessel_yvx;
procedure test_bessel_ivx;
procedure test_bessel_kvx;

procedure test_bessel_ivex;
procedure test_bessel_kvex;


implementation

uses
  amath, specfunx, t_sfx0;


{---------------------------------------------------------------------------}
procedure test_bessel_j0x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NA = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_j0x');

  x := 0.0;
  y := bessel_j0x(x);
  f := 1.0;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_j0x(x);
  f := 0.99999999999975000000;
  testrel(2, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := bessel_j0x(x);
  f := 0.99999976158143510929;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_j0x(x);
  f := -0.39714980986384737229;
  testrel(4, NE, y, f, cnt,failed);

  x := 8.0;
  y := bessel_j0x(x);
  f := 0.17165080713755390609;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_j0x(x);
  f := -0.24593576445134833520;
  testrel(6, NE, y, f, cnt,failed);

  x := 11.791534439014281614;
  y := bessel_j0x(x);
  f := 0.59739493394979385704e-19;
  testabs(7, NA, y, f, cnt,failed);

  x := 100.0;
  y := bessel_j0x(x);
  f := 0.19985850304223122424e-1;
  testabs(8, NA, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_j0x(x);
  f := 0.24786686152420174561e-1;
  testabs(9, NA, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_j1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NA = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_j1x');

  x := 1e-20;
  y := bessel_j1x(x);
  f := 0.5e-20;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_j1x(x);
  f := 0.49999999999993750000e-6;
  testrel(2, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := bessel_j1x(x);
  f := -0.48828119179234139950e-3;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_j1x(x);
  f := -0.66043328023549136143e-1;
  testrel(4, NE, y, f, cnt,failed);

  x := -8.0;
  y := bessel_j1x(x);
  f := -0.23463634685391462438;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_j1x(x);
  f := 0.43472746168861436670e-1;
  testrel(6, NE, y, f, cnt,failed);

  x := 13.323691936314223032;
  y := bessel_j1x(x);
  f := -0.85919366866229194255e-19;
  testabs(7, NA, y, f, cnt,failed);

  x := 100.0;
  y := bessel_j1x(x);
  f := -0.77145352014112158033e-1;
  testrel(8, NE, y, f, cnt,failed);

  x := -1000.0;
  y := bessel_j1x(x);
  f := -0.47283119070895239176e-2;
  testabs(9, NA, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_y0x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NA = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_y0x');

  x := 1e-20;
  y := bessel_y0x(x);
  f := -29.391228250285796860;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_y0x(x);
  f := -8.8690314816594437029;
  testrel(2, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := bessel_y0x(x);
  f := -4.4865150767109739412;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_y0x(x);
  f := -0.1694073932506499190e-1;
  testabs(4, NA, y, f, cnt,failed);

  x := 8.0;
  y := bessel_y0x(x);
  f := 0.22352148938756622053;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_y0x(x);
  f := 0.55671167283599391424e-1;
  testrel(6, NE, y, f, cnt,failed);

  x := 13.36109747387276348;
  y := bessel_y0x(x);
  f := 0.3782632594504609474e-18;
  testabs(7, NA, y, f, cnt,failed);

  x := 100.0;
  y := bessel_y0x(x);
  f := -0.7724431336508315225e-1;
  testabs(8, NA, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_y0x(x);
  f := 0.471591797762281334e-2;
  testabs(9, 10, y, f, cnt,failed);

  x := 2000.0;
  y := bessel_y0x(x);
  f := 0.163683664259955773e-1;
  testabs(10, NA, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_y1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NA = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_y1x');

  x := 1e-20;
  y := bessel_y1x(x);
  f := -63661977236758134308.0;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_y1x(x);
  f := -636619.77237217501376;
  testrel(2, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := bessel_y1x(x);
  f := -651.90099301063115047;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_y1x(x);
  f := 0.39792571055710000525;
  testrel(4, NE, y, f, cnt,failed);

  x := 8.0;
  y := bessel_y1x(x);
  f := -0.15806046173124749426;
  testabs(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_y1x(x);
  f := 0.24901542420695388392;
  testrel(6, NE, y, f, cnt,failed);

  x := 14.89744212833672538;
  y := bessel_y1x(x);
  f := 0.32052154822108592568e-19;
  testabs(7, NA, y, f, cnt,failed);

  x := 100.0;
  y := bessel_y1x(x);
  f := -0.20372312002759793305e-1;
  testabs(8, NA, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_y1x(x);
  f := -0.24784331292351778915e-1;
  testabs(9, NA, y, f, cnt,failed);

  x := 2000.0;
  y := bessel_y1x(x);
  f := -0.70942499636719690210e-2;
  testabs(10, NA, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_jnynx;
var
  cnt, failed: longint;
  n,nmax: integer;
  xt,delta,tol,emax,rms,ssq,xmax: extended;

  procedure onetest(x: extended);
  var
    w1,w2,err,Jn,Jnp1,Jmn,Jmnp1,Yn,Ynp1: extended;
    fail1: boolean;
  begin
    Jnp1  := bessel_jnx(n+1, x);
    Jmn   := bessel_jnx(-n, x);
    Jn    := bessel_jnx(n, x);
    Jmnp1 := bessel_jnx(-(n+1), x);
    {Wronskian HMF[1], 9.1.15: J(n+1,x)*J(-n,x) + J(n,x)*J(-n-1,x) = 0}
    err   := abs(Jnp1*Jmn + Jn*Jmnp1);
    fail1 := err > abs(Jnp1*Jmn)*eps_x;
    if fail1 then writeln('Id test failed: ',n:3,x:23:18,err:28);
    if x<0.0 then begin
      x := -x;
      Jnp1 := bessel_jnx(n+1, x);
      Jn   := bessel_jnx(n, x);
    end;
    Yn   := bessel_ynx(n, x);
    Ynp1 := bessel_ynx(n+1, x);
    {Wronskian HMF[1], 9.1.16: J(n+1,x)*Y(n,x) - J(n,x)*Y(n+1,x) = 2/(Pi*x)}
    w1 := Jnp1 * Yn - Jn * Ynp1;
    w2 := 2.0/(Pi*x);
    err:= abs(w1 - w2);
    ssq:= ssq + sqr(err);
    if err>emax then begin
      emax := err;
      nmax := n;
      xmax := x;
    end;
    if err>tol then begin
      writeln('Wronskian test failed: ',n:3,x:23:18,err:28);
      fail1 := true;
    end;
    if fail1 then inc(failed);
  end;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_jn/yn [Wronskian]');
  delta := 0.25;
  emax  := 0.0;
  rms   := 0.0;
  for n:=-30 to 30 do begin
    tol := (abs(n)+4)*40*eps_x;
    xt  := -4.0;
    ssq := 0.0;
    while xt<40.0 do begin
      if abs(xt)>1e-6 then begin
        inc(cnt);
        onetest(xt);
      end;
      xt := xt+delta;
    end;
    rms := rms + ssq;
    delta := delta + 0.00390625;
  end;
  rms := sqrt(rms/cnt);
  writeln('Peak abs error: ',emax/eps_x:1:3, ' eps for n=',nmax,', x=',xmax:20:18);
  writeln('RMS  abs error: ', rms:20, ' = ',rms/eps_x:1:3,' eps');
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_i0x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_i0x');

  x := 0.0;
  y := bessel_i0x(x);
  f := 1.0;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_i0x(x);
  f := 1.0000000000002500000;
  testrel(2, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := bessel_i0x(x);
  f := 1.0000002384185933124;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_i0x(x);
  f := 11.301921952136330496;
  testrel(4, NE, y, f, cnt,failed);

  x := -8.0;
  y := bessel_i0x(x);
  f := 427.56411572180478518;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_i0x(x);
  f := 2815.7166284662544715;
  testrel(6, NE, y, f, cnt,failed);

  x := 50.0;
  y := bessel_i0x(x);
  f := 0.2932553783849336326710e21;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_i0x(x);
  f := 0.1073751707131073823510e43;
  testrel(8, NE, y, f, cnt,failed);

  x := 500.0;
  y := bessel_i0x(x);
  f := 0.2504809476570078096610e216;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_i0x(x);
  f := 0.2485686096075864174610e433;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_i0x(x);
  f := 0.3513456066149315694310e4341;
  testrel(11, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_i0ex;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_i0ex');

  x := 0.0;
  y := bessel_i0ex(x);
  f := 1.0;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_i0ex(x);
  f := 0.9999990000007499996;
  testrel(2, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := bessel_i0ex(x);
  f := 0.9990241523678519670;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_i0ex(x);
  f := 0.2070019212239866979;
  testrel(4, NE, y, f, cnt,failed);

  x := -8.0;
  y := bessel_i0ex(x);
  f := 0.1434317818568503107;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_i0ex(x);
  f := 0.1278333371634286073;
  testrel(6, NE, y, f, cnt,failed);

  x := 50.0;
  y := bessel_i0ex(x);
  f := 0.5656162664745419253e-1;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_i0ex(x);
  f := 0.3994437929909668265e-1;
  testrel(8, NE, y, f, cnt,failed);

  x := 500.0;
  y := bessel_i0ex(x);
  f := 0.1784570650015316724e-1;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_i0ex(x);
  f := 0.1261724045589125659e-1;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_i0ex(x);
  f := 0.3989472674604732106e-2;
  testrel(11, NE, y, f, cnt,failed);

  x := 1e9;
  y := bessel_i0ex(x);
  f := 0.1261566261167775807e-4;
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_i1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_i1x');

  x := 1e-100;
  y := bessel_i1x(x);
  f := 0.5e-100;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_i1x(x);
  f := 0.50000000000006250000e-6;
  testrel(2, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := bessel_i1x(x);
  f := -0.4882813082076632264e-3;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_i1x(x);
  f := 9.75946515370444991;
  testrel(4, NE, y, f, cnt,failed);

  x := -8.0;
  y := bessel_i1x(x);
  f := -399.8731367825600982;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_i1x(x);
  f := 2670.988303701254654;
  testrel(6, NE, y, f, cnt,failed);

  x := -50.0;
  y := bessel_i1x(x);
  f := -0.2903078590103556797e21;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_i1x(x);
  f := 0.1068369390338162481e43;
  testrel(8, NE, y, f, cnt,failed);

  x := -500.0;
  y := bessel_i1x(x);
  f := -0.2502303412176099996e216;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_i1x(x);
  f := 0.24844429420058669730e433;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_i1x(x);
  f := 0.35132803889537488952e4341;
  testrel(11, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_i1ex;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_i1ex');

  x := 1e-100;
  y := bessel_i1ex(x);
  f := 0.5e-100;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_i1ex(x);
  f := 0.49999950000031249985e-6;
  testrel(2, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := bessel_i1ex(x);
  f := -0.48780470374751535568e-3;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_i1ex(x);
  f := 0.17875083950243532702;
  testrel(4, NE, y, f, cnt,failed);

  x := -8.0;
  y := bessel_i1ex(x);
  f := -0.13414249329269817831;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_i1ex(x);
  f := 0.12126268138445551872;
  testrel(6, NE, y, f, cnt,failed);

  x := -50.0;
  y := bessel_i1ex(x);
  f := -0.55993123892895399644e-1;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_i1ex(x);
  f := 0.39744153025130252674e-1;
  testrel(8, NE, y, f, cnt,failed);

  x := -500.0;
  y := bessel_i1ex(x);
  f := -0.17827851852898056461e-1;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_i1ex(x);
  f := 0.12610930256928629470e-1;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_i1ex(x);
  f := 0.39892731959836622643e-2;
  testrel(11, NE, y, f, cnt,failed);

  x := 1e9;
  y := bessel_i1ex(x);
  f := 0.12615662605369926760e-4;
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_k0x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_k0x');

  x := 1e-100;
  y := bessel_k0x(x);
  f := 230.37444081506298085;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_k0x(x);
  f := 13.931442073626419413;
  testrel(2, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := bessel_k0x(x);
  f := 7.0474052399084523204;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_k0x(x);
  f := 0.11159676085853024270e-1;
  testrel(4, NE, y, f, cnt,failed);

  x := 8.0;
  y := bessel_k0x(x);
  f := 0.14647070522281538710e-3;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_k0x(x);
  f := 0.17780062316167651811e-4;
  testrel(6, NE, y, f, cnt,failed);

  x := 50.0;
  y := bessel_k0x(x);
  f := 0.34101677497894955139e-22;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_k0x(x);
  f := 0.46566282291759020189e-44;
  testrel(8, NE, y, f, cnt,failed);

  x := 500.0;
  y := bessel_k0x(x);
  f := 0.39923216091177928774e-218;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_k0x(x);
  f := 0.20115173162429969967e-435;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_k0x(x);
  f := 0.14231001931183701080e-4344;
  testrel(11, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_k0ex;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_k0ex');

  x := 1e-100;
  y := bessel_k0ex(x);
  f := 230.3744408150629809;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_k0ex(x);
  f := 13.93145600507545876;
  testrel(2, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := bessel_k0ex(x);
  f := 7.054290833146906105;
  testrel(3, NE, y, f, cnt,failed);

  x := 4.0;
  y := bessel_k0ex(x);
  f := 0.6092976692566952693;
  testrel(4, NE, y, f, cnt,failed);

  x := 8.0;
  y := bessel_k0ex(x);
  f := 0.4366230186015861126;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_k0ex(x);
  f := 0.3916319344365986657;
  testrel(6, NE, y, f, cnt,failed);

  x := 50.0;
  y := bessel_k0ex(x);
  f := 0.1768071558574293381;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_k0ex(x);
  f := 0.1251756216591265789;
  testrel(8, NE, y, f, cnt,failed);

  x := 500.0;
  y := bessel_k0ex(x);
  f := 0.56035915417234515427e-1;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_k0ex(x);
  f := 0.3962832160075421711e-1;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_k0ex(x);
  f := 0.1253298471769928529e-1;
  testrel(11, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_k1x;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_k1x');

  x := 1e-100;
  y := bessel_k1x(x);
  f := 1e100;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_k1x(x);
  f := 999999.9999927842790;
  testrel(2, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := bessel_k1x(x);
  f := 1023.9963147439890696;
  testrel(3, NE, y, f, cnt,failed);

  x := 2.0;
  y := bessel_k1x(x);
  f := 0.1398658818165224273;
  testrel(4, NE, y, f, cnt,failed);

  x := 6.0;
  y := bessel_k1x(x);
  f := 0.13439197177355090057e-2;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_k1x(x);
  f := 0.18648773453825584597e-4;
  testrel(6, NE, y, f, cnt,failed);

  x := 50.0;
  y := bessel_k1x(x);
  f := 0.34441022267175556126e-22;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_k1x(x);
  f := 0.46798537356369092866e-44;
  testrel(8, NE, y, f, cnt,failed);

  x := 500.0;
  y := bessel_k1x(x);
  f := 0.39963119385460033495e-218;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_k1x(x);
  f := 0.2012522823712501570e-435;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_k1x(x);
  f := 0.1423171346349328645e-4344;
  testrel(11, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_k1ex;
var
  x,y,f: extended;
  cnt, failed: integer;
{$ifdef VER5X}
const
  NE = 8;
{$else}
const
  NE = 6;
{$endif}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_k1ex');

  x := 1e-100;
  y := bessel_k1ex(x);
  f := 1e100;
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-6;
  y := bessel_k1ex(x);
  f := 1000000.999993284272;
  testrel(2, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := bessel_k1ex(x);
  f := 1024.996799583582939;
  testrel(3, NE, y, f, cnt,failed);

  x := 2.0;
  y := bessel_k1ex(x);
  f := 1.0334768470686885731;
  testrel(4, NE, y, f, cnt,failed);

  x := 6.0;
  y := bessel_k1ex(x);
  f := 0.54217591027713353831;
  testrel(5, NE, y, f, cnt,failed);

  x := 10.0;
  y := bessel_k1ex(x);
  f := 0.41076657059578875113;
  testrel(6, NE, y, f, cnt,failed);

  x := 50.0;
  y := bessel_k1ex(x);
  f := 0.17856655855881557460;
  testrel(7, NE, y, f, cnt,failed);

  x := 100.0;
  y := bessel_k1ex(x);
  f := 0.12579995047957852933;
  testrel(8, NE, y, f, cnt,failed);

  x := 500.0;
  y := bessel_k1ex(x);
  f := 0.56091923370555569238e-1;
  testrel(9, NE, y, f, cnt,failed);

  x := 1000.0;
  y := bessel_k1ex(x);
  f := 0.39648130812960210480e-1;
  testrel(10, NE, y, f, cnt,failed);

  x := 10000.0;
  y := bessel_k1ex(x);
  f := 0.12533611351270505734e-1;
  testrel(11, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_inknx;
var
  cnt, failed: longint;
  n,nmax: integer;
  xt,delta,tol,emax,rms,ssq,xmax: extended;

  procedure onetest(x: extended);
  var
    w1,w2,err,Inx,Inp1,Imn,Imnp1,Kn,Knp1: extended;
    fail1: boolean;
  begin
    Inp1  := bessel_inx(n+1, x);
    Imn   := bessel_inx(-n, x);
    Inx   := bessel_inx(n, x);
    Imnp1 := bessel_inx(-(n+1), x);
    {Wronskian HMF[1], 9.6.14: I(n,x)*I(-n-1,x) - I(n+1,x)*I(-n,x) = 0}
    err   := abs(Inx*Imnp1 - Inp1*Imn);
    fail1 := err<>0;
    if fail1 then writeln('Id test failed: ',n:3,x:23:18,err:28);
    Kn   := bessel_knx(n, x);
    Knp1 := bessel_knx(n+1, x);
    {Wronskian HMF[1], 9.6.15: I(n,x)*K(n+1,x) + I(n+1,x)*K(n,x) = 1/x}
    w1 := Inx*Knp1 + Inp1*Kn;
    w2 := 1.0/x;
    err:= abs(w1 - w2);
    ssq:= ssq + sqr(err);
    if err>emax then begin
      emax := err;
      nmax := n;
      xmax := x;
    end;
    if err>tol then begin
      writeln('Wronskian test failed: ',n:3,x:23:18,err:28);
      fail1 := true;
    end;
    if fail1 then inc(failed);
  end;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_in/kn [Wronskian]');
  delta := 0.25;
  emax  := 0.0;
  rms   := 0.0;
  for n:=0 to 30 do begin
    tol := (abs(n)+4)*40*eps_x;
    xt  := 0.03125;
    ssq := 0.0;
    while xt<30.0 do begin
      if abs(xt)>1e-6 then begin
        inc(cnt);
        onetest(xt);
      end;
      xt := xt+delta;
    end;
    rms := rms + ssq;
    delta := delta + 0.00390625;
  end;
  rms := sqrt(rms/cnt);
  writeln('Peak abs error: ',emax/eps_x:1:3, ' eps for n=',nmax,', x=',xmax:20:18);
  writeln('RMS  abs error: ', rms:20, ' = ',rms/eps_x:1:3,' eps');
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_jvx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 10;
  NE2 = 20;
  NE3 = 30;

begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_jvx');


  y := bessel_jvx(2457/1024, 1/1024);
  f := 0.3807399201186033356e-8;
  testrel(1, NE, y, f, cnt,failed);

  y := bessel_jvx(5.5, 3217/1024);
  f := 0.2819330762575060916e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := bessel_jvx(-5.5, 3217/1024);
  f :=  -0.2558200644706479118e1;
  testrel(3, NE, y, f, cnt,failed);

  y := bessel_jvx(-5.5, 1e+04);
  f := 0.2449843111985605522e-02;
  testrel(4, NE, y, f, cnt,failed);

  y := bessel_jvx(5.5, 1e+04);
  f := 0.7593435027226703614e-2;
  testrel(5, NE, y, f, cnt,failed);

  y := bessel_jvx(5.5, 1e+06);
  f := -0.7474242485956301774e-3;
  testrel(6, NE, y, f, cnt,failed);

  y := bessel_jvx(5.125, 1e+06);
  f := -0.7766001248357042806e-3;
  testrel(7, NE, y, f, cnt,failed);

  y := bessel_jvx(5.875, 1e+06);
  f := -0.4663227211151930716e-3;
  testrel(8, NE, y, f, cnt,failed);

  y := bessel_jvx(0.5, 101);
  f := 0.358874487875643822e-1;
  testrel(9, NE2,  y, f, cnt,failed);    {!!}

  y := bessel_jvx(-5.5, 1e+04);
  f := 0.2449843111985605522e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := bessel_jvx(-5.5, 1e+06);
  f := 0.2792432004335795111e-3;
  testrel(11, NE, y, f, cnt,failed);

  y := bessel_jvx(-0.5, 101);
  f := 0.7081847980975942685e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := bessel_jvx(-10486074/1048576, 0.0009765625);
  f := 0.1414740131604946958e36;
  testrel(13, NE, y, f, cnt,failed);

  y := bessel_jvx(-10486074/1048576, 15);
  f := -0.9022392888854233096e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := bessel_jvx(-10486074/1048576, 100);
  f := -0.5476136603168065513e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := bessel_jvx(10486074/1048576, 100);
  f := -0.5470649146151378076e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := bessel_jvx(10486074/1048576, 10000);
  f := 0.7112613531849121437e-2;
  testrel(17, NE, y, f, cnt,failed);

  y := bessel_jvx(5.5, 3.125);
  f := 0.27499449389648972995e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := bessel_jvx(10, 1600);
  f := 0.19642250872681627003e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := bessel_jvx(180, 3);
  f := 0.24437475117168956889e-297;
  testrel(20, NE, y, f, cnt,failed);

  y := bessel_jvx(-180, 3);
  f := 0.24437475117168956889e-297;
  testrel(21, NE, y, f, cnt,failed);

  y := bessel_jvx(181, 3);
  f := 0.20253434642874977521e-299;
  testrel(22, NE, y, f, cnt,failed);

  y := bessel_jvx(-181, 3);
  f := -0.20253434642874977521e-299;
  testrel(23, NE, y, f, cnt,failed);

  y := bessel_jvx(0.25, 0.75);
  f := 0.76921463928572487957;
  testrel(24, NE, y, f, cnt,failed);

  y := bessel_jvx(8.25, 0.75);
  f := 0.43614753464277517659e-8;
  testrel(25, NE, y, f, cnt,failed);

  y := bessel_jvx(2.5, 2.5);
  f := 0.32809141153443809388;
  testrel(26, NE, y, f, cnt,failed);

  y := bessel_jvx(-2.5, 2.5);
  f := 0.57263060443914839171;
  testrel(27, NE, y, f, cnt,failed);

  y := bessel_jvx(2.25, 2.5);
  f := 0.38844121553062244020;
  testrel(28, NE, y, f, cnt,failed);

  y := bessel_jvx(-2.25, 2.5);
  f := 0.61507674032895119303;
  testrel(29, NE, y, f, cnt,failed);

  y := bessel_jvx(3.75, 2000);
  f := 0.23139400409732231894e-3;
  testabs(30, 1, y, f, cnt,failed);      {!!!!!!!!!!!}

  y := bessel_jvx(2.25, ldexp(1,-450));
  f := 0.13286123502124834232e-305;
  testrel(31, NE, y, f, cnt,failed);

  y := bessel_jvx(141.75, 100);
  f := 0.58807854686058285162e-12;
  testrel(32, NE, y, f, cnt,failed);

  y := bessel_jvx(-141.75, 100);
  f := -0.38100980344476687750e10;
  testrel(33, NE, y, f, cnt,failed);

  y := bessel_jvx(-141.75, 20000);
  f := -0.56335566400200594373e-2;
  testrel(34, NE, y, f, cnt,failed);

  y := bessel_jvx(-3, -1200);
  f := 0.17606940067282534608e-1;
  testrel(35, NE, y, f, cnt,failed);

  y := bessel_jvx(3, -1200);
  f := -0.17606940067282534608e-1;
  testrel(36, NE, y, f, cnt,failed);

  y := bessel_jvx(-3, 1200);
  f := -0.17606940067282534608e-1;
  testrel(37, NE, y, f, cnt,failed);

  y := bessel_jvx(3, 1200);
  f := 0.17606940067282534608e-1;
  testrel(38, NE, y, f, cnt,failed);

  y := bessel_jvx(185.75, 3);
  f := 0.24357953735629609598e-309;
  testrel(39, NE, y, f, cnt,failed);

  y := bessel_jvx(500.5, 3);
  f := 0.49602963937019509756e-1047;     {!!!!!!}
  testrel(40, NE3, y, f, cnt,failed);

  y := bessel_jvx(-4001, 200);
  f := -0.11223600055016688175e-4675;    {!!!!!!!}
  testrel(41, NE3, y, f, cnt,failed);

  y := bessel_jvx(5000, 400);
  f := 0.11151663553767544274e-4823;
  testrel(42, NE2, y, f, cnt,failed);    {!!!!!!!}

  y := bessel_jvx(1800, 10);
  f := 0.22531557315715245301e-3821;
  testrel(43, NE2, y, f, cnt,failed);

  {Test case for fixed typo (V1.19.00), old vers. rel.err. ~ 16.6 eps_x}
  y := bessel_jvx(1.25, 380);
  f := 0.3953692150771290695e-1;
  testrel(44, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{Avoid crash with numbers like -0.116480033881283810419007840427e4933 !!!}

{$ifdef VER50}
{$define bugspec}
{$endif}

{$ifdef FPC}
{$ifdef VER1_0}
{$define bugspec}
{$endif}
{$ifdef VER2_0}
{$define bugspec}
{$endif}
{$endif}


{---------------------------------------------------------------------------}
procedure test_bessel_yvx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NA = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_yvx');

  y := bessel_yvx(5.5, 3.125);
  f := -2.6148944032841746878;
  testrel(1, NE, y, f, cnt,failed);

  y := bessel_yvx(10, 1600);
  f := 0.34752108447889821839e-2;
  testabs(2, NA, y, f, cnt,failed);

  y := bessel_yvx(180, 3);
  f := -0.72373840437520354195e295;
  testrel(3, NE, y, f, cnt,failed);

  y := bessel_yvx(-180, 3);
  f := -0.72373840437520354195e295;
  testrel(4, NE, y, f, cnt,failed);

  y := bessel_yvx(181, 3);
  f := -0.86842543249718030285e297;
  testrel(5, NE, y, f, cnt,failed);

  y := bessel_yvx(-181, 3);
  f := 0.86842543249718030285e297;
  testrel(5, NE, y, f, cnt,failed);

  y := bessel_yvx(0.25, 0.75);
  f := -0.43995599531033785111;
  testrel(7, NE, y, f, cnt,failed);

  y := bessel_yvx(8.25, 0.75);
  f := -8883668.6901696101586;
  testrel(7, NE, y, f, cnt,failed);

  y := bessel_yvx(2.5, 2.5);
  f := -0.57263060443914839171;
  testrel(9, NE, y, f, cnt,failed);

  y := bessel_yvx(2.5, 2.5);
  f := -0.57263060443914839171;
  testrel(10, NE, y, f, cnt,failed);

  y := bessel_yvx(2.25, 2.5);
  f := -0.48140865254281476467;
  testrel(11, NE, y, f, cnt,failed);

  y := bessel_yvx(-2.25, 2.5);
  f := -0.65737905140854417532e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := bessel_yvx(3.75, 2000);
  f := 0.17839755956340562681e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := bessel_yvx(5.5, 3.125);
  f := -2.6148944032841746878;
  testrel(14, NE, y, f, cnt,failed);

  y := bessel_yvx(-5.5, 3.125);
  f := -0.27499449389648972995e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := bessel_yvx(2.25, ldexp(1,-450));
  f := -0.10648031421919560336e306;
  testrel(16, NE, y, f, cnt,failed);

  y := bessel_yvx(-10486074/1048576, 0.0009765625);
  f := -0.150382374389531766117869e39;
  testrel(17, NE, y, f, cnt,failed);

  y := bessel_yvx(-10486074/1048576, 100);
  f :=  0.58304189131902600995578e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := bessel_yvx(141.75, 100);
  f := -0.53882923142869650729e10;
  testrel(19, NE, y, f, cnt,failed);

  y := bessel_yvx(-141.75, 100);
  f := -0.38100980344476687750e10;
  testrel(20, NE, y, f, cnt,failed);

  y := bessel_yvx(-141.75, 20000);
  f := 0.30794106544903641394e-3;
  testabs(21, NA, y, f, cnt,failed);

  y := bessel_yvx(185.75, 3);
  f := -0.70361834614225829376e307;
  testrel(22, NE, y, f, cnt,failed);

  y := bessel_yvx(-1700.5, 1.5);
  f := 0.0;
  testrel(23, NE, y, f, cnt,failed);

  y := bessel_yvx(-1688.75, 1.5);
  f := 0.1837802416599446436e4927;
  testrel(24, NE, y, f, cnt,failed);

{$ifndef WINCRT}
  {** Underflow under 16 bit windows}
  y := bessel_yvx(-1688.5, 1.5);
  f := 0.49959862241878760459e-4929;
  testrel(25, NE, y, f, cnt,failed);
{$endif}

  y := bessel_yvx(500.5, 3);
  f := -0.12821717964012682244e1045;
  testrel(26, 20, y, f, cnt,failed);    {!!!!!!!!!!!}

  y := bessel_yvx(-4001, 200);
  f := 0.70972919299081125989e4672;
  testrel(27, 30, y, f, cnt,failed);    {!!!!!!!!!!!}

  y := bessel_yvx(5000, 400);
  f := -0.57270987738165510091e4820;
  testrel(28, NE, y, f, cnt,failed);

  {near overflow}
  y := bessel_yvx(8.125, ldexp(1,-2014));  {Yv_series}
  f := -0.5503169807480362694e4932;
  testrel(29, NE, y, f, cnt,failed);

{$ifdef bugspec}
  f := -0.1164800338812838104e4911;
  f := f*1e22;
{$else}
  f := -0.1164800338812838104e4933;
{$endif}
  y := bessel_yvx(1856.0625, 3);
  testrel(30, 2*NE, y, f, cnt,failed);

  y := bessel_yvx(1724.8, 1.75);
  f := -0.9207826447251883707e4932;
  testrel(31, 3500, y, f, cnt,failed);             {!!!!!!!!!!!}

{$ifdef bugspec}
  f := -0.1186508712139644935e4911;
  f := f*1e22;
{$else}
  f := -0.1186508712139644935e4933;
{$endif}
  y := bessel_yvx(1600.8, 0.975);  {Yv_series}
  testrel(32, 3500, y, f, cnt,failed);             {!!!!!!!!!!!}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_bessel_ivx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_ivx');

  y := bessel_ivx(0, -10);
  f := 2815.7166284662544715;
  testrel(1, NE, y, f, cnt,failed);

  y := bessel_ivx(1, -10);
  f := -2670.9883037012546543;
  testrel(2, NE, y, f, cnt,failed);

  y := bessel_ivx(2, -10);
  f := 2281.5189677260035406;
  testrel(3, NE, y, f, cnt,failed);

  y := bessel_ivx(3, -10);
  f := -1758.3807166108532381;
  testrel(4, NE, y, f, cnt,failed);

  y := bessel_ivx(0.5, 10);
  f := 2778.7846038745710240;
  testrel(5, NE, y, f, cnt,failed);

  y := bessel_ivx(-0.5, 10);
  f := 2778.7846153295749521;
  testrel(6, NE, y, f, cnt,failed);

  y := bessel_ivx(0.25, 100);
  f := 0.10734145166453237066e43;
  testrel(7, NE, y, f, cnt,failed);

  y := bessel_ivx(4.75, 100);
  f := 0.95867597530423643956e42;
  testrel(8, NE, y, f, cnt,failed);

  y := bessel_ivx(5, 100);
  f := 0.94700938730355812462e42;
  testrel(9, NE, y, f, cnt,failed);

  y := bessel_ivx(4.5, 100);
  f := 0.96987756405711455547e42;
  testrel(10, NE, y, f, cnt,failed);

  y := bessel_ivx(3.1, 4.5);
  f := 5.5330623502911278334;
  testrel(11, NE, y, f, cnt,failed);

  y := bessel_ivx(3.1, 1.5);
  f := 0.68893840471849860622e-1;
  testrel(12, NE, y, f, cnt,failed);

  {gsl}
  y := bessel_ivx(0.0001,10.0);
  f := 2815.7166269770030352;
  testrel(13, NE, y, f, cnt,failed);

  y := bessel_ivx( 1.0, 0.001);
  f := 0.0005000000625000026042;
  testrel(14, NE, y, f, cnt,failed);

  y := bessel_ivx( 1.0,   1.0);
  f := 0.5651591039924850272;
  testrel(15, NE, y, f, cnt,failed);

  y := bessel_ivx(30.0,   1.0);
  f := 3.539500588106447747e-42;
  testrel(16, NE, y, f, cnt,failed);

  y := bessel_ivx(30.0, 100.0);
  f := 1.2061548704498434006e+40;
  testrel(17, NE, y, f, cnt,failed);

  y := bessel_ivx(10.0,   1.0);
  f := 2.7529480398368736252e-10;
  testrel(18, NE, y, f, cnt,failed);

  y := bessel_ivx(10.0, 100.0);
  f := 6.498975524720147799e+41;
  testrel(19, NE, y, f, cnt,failed);

  y := bessel_ivx(10.2, 100.0);
  f := 6.368587361287030443e+41;
  testrel(20, NE, y, f, cnt,failed);

  {boost}
  y := bessel_ivx(2.25, 1/(1024*1024));
  f := 2.343792121334813472e-15;
  testrel(21, NE, y, f, cnt,failed);

  y := bessel_ivx(5.5, 3.125);
  f := 0.05835140459893715005;
  testrel(22, NE, y, f, cnt,failed);

  y := bessel_ivx(-5 + 1/1024,2.125);
  f := 0.02679209380095710237;
  testrel(23, NE, y, f, cnt,failed);

  y := bessel_ivx(-5.5, 10);
  f := 597.5776069613691696;
  testrel(24, NE, y, f, cnt,failed);

  y := bessel_ivx(-5.5, 100);
  f := 9.22362906144706872e41;
  testrel(25, NE, y, f, cnt,failed);

  y := bessel_ivx(-10486074/(1024*1024),1/1024);
  f := 1.414740056651813504e35;
  testrel(26, NE, y, f, cnt,failed);

  y := bessel_ivx(-10486074/(1024*1024), 50);
  f := 1.071532772029006715e20;
  testrel(27, NE, y, f, cnt,failed);

  y := bessel_ivx(144794/1024, 100);
  f := 2066.276947573926604;
  testrel(28, NE, y, f, cnt,failed);

  y := bessel_ivx(144794/1024, 200);
  f := 2.236997394722469288e64;
  testrel(29, NE, y, f, cnt,failed);

  y := bessel_ivx(-144794/1024, 100);
  f := 2066.276946727631909;
  testrel(30, NE, y, f, cnt,failed);

  {extended only}
  y := bessel_ivx(0.25, 1000);
  f := 0.24856083807196953877e433;
  testrel(31, NE, y, f, cnt,failed);

  y := bessel_ivx(0.25, 10000);
  f := 0.35134450860672298023e4341;
  testrel(32, NE, y, f, cnt,failed);

  y := bessel_ivx(0.5, 1000);
  f := 0.24853752492344490307e433;
  testrel(33, NE, y, f, cnt,failed);

  y := bessel_ivx(0.5, 10000);
  f := 0.35134121460268650739e4341;
  testrel(34, NE, y, f, cnt,failed);

  y := bessel_ivx(0.5, 11350);
  f := 0.65430725674015483388e4927;
  testrel(35, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_ivex;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_ivex');

  y := bessel_ivex(0, -10);
  f := 0.12783333716342860733;
  testrel(1, NE, y, f, cnt,failed);

  y := bessel_ivex(3, -10);
  f := -0.79830361029840517288e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := bessel_ivex(0.5, 10);
  f := 0.12615662584097981553;
  testrel(3, NE, y, f, cnt,failed);


  y := bessel_ivex(-0.5, 10);
  f := 0.12615662636103618930;
  testrel(4, NE, y, f, cnt,failed);

  y := bessel_ivex(0.25, 100);
  f := 0.39931835556842864588e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := bessel_ivex(4.75, 100);
  f := 0.35663474645176342070e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := bessel_ivex(-3.1, 4.5);
  f := 0.61430696231288476674e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := bessel_ivex(3.1, 1.5);
  f := 0.15372293657724235852e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := bessel_ivex(2.25, 1/(1024*1024));
  f := 0.23437898861215301925e-14;
  testrel(9, NE, y, f, cnt,failed);

  y := bessel_ivex(144794/1024, 100);
  f := 0.76867072324754287139e-40;
  testrel(10, NE, y, f, cnt,failed);

  y := bessel_ivex(-144794/1024, 100);
  f := 0.76867072293271481512e-40;
  testrel(11, NE, y, f, cnt,failed);

  y := bessel_ivex(0.5, 10000);
  f := 0.39894228040143267792e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := bessel_ivex(0.5, 1e8);
  f := 0.39894228040143267794e-4;
  testrel(13, NE, y, f, cnt,failed);


  {gsl}
  y := bessel_ivex(0,1e-10);
  f := 0.99999999990000000001;
  testrel(14, NE, y, f, cnt,failed);

  y := bessel_ivex(0,65536.0);
  f := 0.0015583712551952223537;
  testrel(15, NE, y, f, cnt,failed);

  y := bessel_ivex(1,0.1);
  f := 0.04529844680880932501;
  testrel(16, NE, y, f, cnt,failed);

  y := bessel_ivex(1,65536.0);
  f := 0.0015583593657207350452;
  testrel(17, NE, y, f, cnt,failed);

  y := bessel_ivex(-4,0.1);
  f := 2.3575258620054605307e-07;
  testrel(18, NE, y, f, cnt,failed);

  y := bessel_ivex(4,0.1);
  f := 2.3575258620054605307e-07;
  testrel(19, NE, y, f, cnt,failed);

  y := bessel_ivex(5,2.0);
  f := 0.1329761094188157814e-2;
  testrel(20, NE, y, f, cnt,failed);

  y := bessel_ivex(100, 100.0);
  f := 1.7266862628167695785e-22;
  testrel(21, NE, y, f, cnt,failed);

  y := bessel_ivex(0.0001,10.0);
  f := 0.1278333370958166967;
  testrel(22, NE, y, f, cnt,failed);

  y := bessel_ivex(1.0, 0.001);
  f := 0.0004995003123542213370;
  testrel(23, NE, y, f, cnt,failed);

  y := bessel_ivex(1.0, 1.0);
  f := 0.2079104153497084489;
  testrel(24, NE, y, f, cnt,failed);

  y := bessel_ivex(30.0, 1.0);
  f := 1.3021094983785914437e-42;
  testrel(25, NE, y, f, cnt,failed);

  y := bessel_ivex(30.0, 100.0);
  f := 0.0004486987756920986146;
  testrel(26, NE, y, f, cnt,failed);

  y := bessel_ivex(10.0, 1.0);
  f := 1.0127529864692066036e-10;
  testrel(27, NE, y, f, cnt,failed);

  y := bessel_ivex(10.0, 100.0);
  f := 0.024176682718258828365;
  testrel(28, NE, y, f, cnt,failed);

  y := bessel_ivex(10.2, 100.0);
  f := 0.023691628843913810043;
  testrel(29, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_bessel_kvx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_kvx');


  y := bessel_kvx(0.5, 10);
  f := 0.17993478093705179608e-4;
  testrel(1, NE, y, f, cnt,failed);

  y := bessel_kvx(-0.5, 10);
  f := 0.17993478093705179608e-4;
  testrel(2, NE, y, f, cnt,failed);

  y := bessel_kvx(0.25, 100);
  f := 0.46580764515098397833e-44;
  testrel(3, NE, y, f, cnt,failed);

  y := bessel_kvx(4.75, 100);
  f := 0.52097168282519795500e-44;
  testrel(4, NE, y, f, cnt,failed);

  y := bessel_kvx(5, 100);
  f := 0.52732561132929498946e-44;
  testrel(5, NE, y, f, cnt,failed);

  y := bessel_kvx(4.5, 100);
  f := 0.51501415511003247155e-44;
  testrel(6, NE, y, f, cnt,failed);

  y := bessel_kvx(3.1, 4.5);
  f := 0.1650669628003692240e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := bessel_kvx(3.1, 1.5);
  f := 2.0950762832264046499;
  testrel(8, NE, y, f, cnt,failed);


  {GSL}
  y := bessel_kvx(0.0001,0.001);
  f := 7.023689431812884141;
  testrel(9, NE, y, f, cnt,failed);

  y := bessel_kvx(0.0001,10.0);
  f := 0.000017780062324654874306;
  testrel(10, NE, y, f, cnt,failed);

  y := bessel_kvx( 1.0, 0.001);
  f := 999.9962381560855743;
  testrel(11, NE, y, f, cnt,failed);

  y := bessel_kvx( 1.0,   1.0);
  f := 0.6019072301972345747;
  testrel(12, NE, y, f, cnt,failed);

  y := bessel_kvx(10.0, 0.001);
  f := 1.8579455483904008064e+38;
  testrel(13, NE, y, f, cnt,failed);

  y := bessel_kvx(10.0,   1.0);
  f := 1.8071328990102945469e+08;
  testrel(14, NE, y, f, cnt,failed);

  y := bessel_kvx(10.0, 100.0);
  f := 7.655427977388100611e-45;
  testrel(15, NE, y, f, cnt,failed);

  y := bessel_kvx(10.2, 100.0);
  f := 7.810600225948217841e-45;
  testrel(16, NE, y, f, cnt,failed);

  y := bessel_kvx(30.0,   1.0);
  f := 4.706145526783626883e+39;
  testrel(17, NE, y, f, cnt,failed);

  y := bessel_kvx(30.0, 100.0);
  f := 3.970602055959398739e-43;
  testrel(18, NE, y, f, cnt,failed);

  {Boost}
  y := bessel_kvx(0.5, 0.875);
  f := 0.5585322316466086461;
  testrel(19, NE, y, f, cnt,failed);

  y := bessel_kvx(0.5, 1.125);
  f := 0.3836210106501895471;
  testrel(20, NE, y, f, cnt,failed);

  y := bessel_kvx(2.25, ldexp(1.0, -30));
  f := 5.623973927192832713e20;
  testrel(21, NE, y, f, cnt,failed);

  y := bessel_kvx(5.5, 3217/1024);
  f := 1.306232887750125963;
  testrel(22, NE, y, f, cnt,failed);

  y := bessel_kvx(-5.5, 10);
  f := 0.7330453007985021646e-4;
  testrel(23, NE, y, f, cnt,failed);

  y := bessel_kvx(-5.5, 100);
  f := 5.412745553067922673e-45;
  testrel(24, NE, y, f, cnt,failed);

  y := bessel_kvx(10240/1024, 1/1024);
  f := 2.355225792639220762e38;
  testrel(25, NE, y, f, cnt,failed);

  y := bessel_kvx(10240/1024, 10);
  f := 0.1614255300390670023e-2;
  testrel(26, NE, y, f, cnt,failed);

  y := bessel_kvx(144793/1024, 100);
  f := 1.395652458603025281e-6;
  testrel(27, NE, y, f, cnt,failed);

  y := bessel_kvx(144793/1024, 200);
  f := 9.119504120432254322e-68;
  testrel(28, NE, y, f, cnt,failed);

  y := bessel_kvx(-144793/1024, 50);
  f := 1.3018522971752502517e42;
  testrel(29, NE, y, f, cnt,failed);

  {extended only}
  y := bessel_kvx(0.5, 1000);
  f := 0.20117686460183875401e-435;
  testrel(30, NE, y, f, cnt,failed);

  y := bessel_kvx(0.5, 10000);
  f := 0.14231179810926081515e-4344;
  testrel(31, NE, y, f, cnt,failed);

  y := bessel_kvx(0.5, 11350);
  f := 0.67327487174145877013e-4931;
  testrel(32, NE, y, f, cnt,failed);

  y := bessel_kvx(0.25, 1000);
  f := 0.20115801457440548596e-435;
  testrel(33, NE, y, f, cnt,failed);

  y := bessel_kvx(0.25, 10000);
  f := 0.14231046400910860951e-4344;
  testrel(34, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_bessel_kvex;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bessel_kvex');

  y := bessel_kvex(0.5, 10);
  f := 0.3963327297606011013;
  testrel(1, NE, y, f, cnt,failed);

  y := bessel_kvex(-0.5, 10);
  f := 0.3963327297606011013;
  testrel(2, NE, y, f, cnt,failed);

  y := bessel_kvex(0.25, 100);
  f := 0.1252145515719367696;
  testrel(3, NE, y, f, cnt,failed);

  y := bessel_kvex(4.75, 100);
  f := 0.1400432911003213388;
  testrel(4, NE, y, f, cnt,failed);

  y := bessel_kvex(5, 100);
  f := 0.1417513015132950781;
  testrel(5, NE, y, f, cnt,failed);

  y := bessel_kvex(-4.25, 100);
  f := 0.13694375389121990656;
  testrel(6, NE, y, f, cnt,failed);

  {GSL}
  y := bessel_kvex(100, 100.0);
  f := 2.0475736731166756813e+19;
  testrel(7, NE, y, f, cnt,failed);

  y := bessel_kvex(0.0001,0.001);
  f := 7.0307166342603205496;
  testrel(8, NE, y, f, cnt,failed);

  y := bessel_kvex(0.0001,10.0);
  f := 0.3916319346235421817;
  testrel(9, NE, y, f, cnt,failed);

  y := bessel_kvex( 1.0, 0.001);
  f := 1000.9967345590684524;
  testrel(10, NE, y, f, cnt,failed);

  y := bessel_kvex( 1.0, 1.0);
  f := 1.6361534862632582465;
  testrel(11, NE, y, f, cnt,failed);

  y := bessel_kvex(30.0, 1.0);
  f := 1.2792629867539753925e+40;
  testrel(12, NE, y, f, cnt,failed);

  y := bessel_kvex(30.0, 100.0);
  f := 10.673443449954850040;
  testrel(13, NE, y, f, cnt,failed);

  y := bessel_kvex(10.0, 1.0);
  f := 4.912296520990198599e+08;
  testrel(14, NE, y, f, cnt,failed);

  y := bessel_kvex(10.0, 100.0);
  f := 0.20578687173955779807;
  testrel(14, NE, y, f, cnt,failed);

  y := bessel_kvex(10.0, 1000.0);
  f := 0.04165905142800565788;
  testrel(16, NE, y, f, cnt,failed);

  y := bessel_kvex(10.0, 1.0e+8);
  f := 0.1253314762406078994e-3;
  testrel(17, NE, y, f, cnt,failed);

  y := bessel_kvex(10.2, 100.0);
  f := 0.20995808355244385075;
  testrel(18, NE, y, f, cnt,failed);

  y := bessel_kvex(-144793/1024, 50);
  f := 0.67497208025683004008e64;
  testrel(19, NE, y, f, cnt,failed);

  y := bessel_kvex(0.25, 10000);
  f := 0.12533023881379649469e-1;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.
