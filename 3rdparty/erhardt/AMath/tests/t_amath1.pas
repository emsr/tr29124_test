{Part 1 of regression test for AMath unit,   (c) 2011  W.Ehrhardt}
{Tests for basic/statistic functions}

unit t_amath1;

{$i STD.INC}

{$ifdef BIT16}
{$N+}
{$endif}

interface


procedure test_exp;
procedure test_exp2;
procedure test_exp2m1;
procedure test_exp3;
procedure test_exp10;
procedure test_exp10m1;
procedure test_expm1;
procedure test_exprel;
procedure test_expx2;
procedure test_expmx2h;
procedure test_ln;
procedure test_ln1mexp;
procedure test_ln1pexp;
procedure test_log_add_sub_exp;
procedure test_ln1p;
procedure test_ln1pmx;
procedure test_lnxp1;
procedure test_lncosh;
procedure test_lnsinh;
procedure test_log10;
procedure test_log2;
procedure test_log2p1;
procedure test_log10p1;
procedure test_logistic;
procedure test_logit;
procedure test_sqrt1pmx;
procedure test_trig_degree;

procedure test_compound;
procedure test_comprel;
procedure test_pow1p;
procedure test_powpi2k;

procedure test_polevl;
procedure test_polevale;
procedure test_polevalche;

implementation

uses
  amath, t_amath0;

{---------------------------------------------------------------------------}
procedure test_exp;
const
  NT=13;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx: 1e-15;   ty: 1.0000000000000010000),
           (tx: 1e-9;    ty: 1.0000000010000000005),
           (tx: -1e-10;  ty: 0.99999999990000000001),
           (tx: 1e-5;    ty: 1.0000100000500001667),
           (tx: 0.125;   ty: 1.1331484530668263168),
           (tx: 0.5;     ty: 1.6487212707001281468),
           (tx: 1.0;     ty: 2.7182818284590452354),
           (tx: 2.0;     ty: 7.3890560989306502272),
           (tx: -6.0;    ty: 0.24787521766663584230e-2),
           (tx: 10.0;    ty: 22026.465794806716517),
           (tx: 100.0;   ty: 0.26881171418161354484e44),
           (tx: 1e4;     ty: 0.88068182256629215873e4343)
         );
const
  name = 'exp';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := exp(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
    x := -tx;
    y := 1.0/ty;
    z := exp(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_exp2;
const
  NT=13;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx: 1e-15;   ty: 1.0000000000000006931),
           (tx: 1e-9;    ty: 1.0000000006931471808),
           (tx: -1e-10;  ty: 0.99999999993068528195),
           (tx: 1e-5;    ty: 1.0000069314958283057),
           (tx: 0.125;   ty: 1.0905077326652576592),
           (tx: 0.5;     ty: 1.4142135623730950488),
           (tx: 1.0;     ty: 2.0),
           (tx: 2.0;     ty: 4.0),
           (tx: -6.7;    ty: 0.96183157292571584727e-2),
           (tx: 12.345;  ty: 5202.5384244788822035),
           (tx: 100.0;   ty: 0.12676506002282294015e31),
           (tx: 1.333e4; ty: 0.53683670350017507020e4013)
         );
const
  name = 'exp2';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := exp2(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
    x := -tx;
    y := 1/ty;
    z := exp2(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_exp2m1;
const
  NT=16;
  TData: array[1..NT] of TPair = (
           (tx:  0.0;     ty:  0.0),
           (tx: -1e-15;   ty: -0.69314718055994506919e-15),
           (tx:  1e-9;    ty:  0.69314718080017181643e-9),
           (tx: -1e-10;   ty: -0.69314718053592265872e-10),
           (tx:  1e-5;    ty:  0.69314958283056532091e-5),
           (tx:  0.125;   ty:  0.90507732665257659207e-1),
           (tx: -0.125;   ty: -0.82995956795328768256e-1),
           (tx:  0.5;     ty:  0.41421356237309504880),
           (tx: -0.75;    ty: -0.40539644249863946664),
           (tx:  1.0;     ty:  1.0),
           (tx: -1.0;     ty: -0.5),
           (tx:  2.0;     ty:  3.0),
           (tx: -6.7;     ty: -0.99038168427074284153),
           (tx:  12.345;  ty:  5201.5384244788822035),
           (tx:  100.0;   ty:  0.12676506002282294015e31),
           (tx:  1.333e4; ty:  0.53683670350017507020e4013)
          );
const
  name = 'exp2m1';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := exp2m1(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_exp3;
const
  NT=13;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx: 1e-15;   ty: 1.0000000000000010986),
           (tx: 1e-8;    ty: 1.0000000109861229470),
           (tx: -1e-10;  ty: 0.99999999989013877114),
           (tx: 1e-5;    ty: 1.0000109861832343501),
           (tx: 0.125;   ty: 1.1472026904398770895),
           (tx: 0.5;     ty: 1.7320508075688772935),
           (tx: 1.0;     ty: 3.0),
           (tx: 2.0;     ty: 9.0),
           (tx: -6.7;    ty: 0.63575179255414235962e-3),
           (tx: 100.25;  ty: 0.67827496189528453224e48),
           (tx: 234.5;   ty: 0.76724529126378281597e112),
           (tx: 8000.5;  ty: 0.16165844041233134703e3818)
         );
const
  name = 'exp3';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := exp3(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
    x := -tx;
    y := 1/ty;
    z := exp3(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_exp10;
const
  NT=13;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx: 1e-15;   ty: 1.0000000000000023026),
           (tx: 1e-8;    ty: 1.0000000230258511950),
           (tx: -1e-10;  ty: 0.99999999976974149073),
           (tx: 1e-5;    ty: 1.0000230261160268807),
           (tx: 0.125;   ty: 1.3335214321633240257),
           (tx: 0.5;     ty: 3.1622776601683793320),
           (tx: 1.0;     ty: 10.0),
           (tx: 2.0;     ty: 100.0),
           (tx: -6.7;    ty: 0.19952623149688796014e-6),
           (tx: 100.25;  ty: 0.17782794100389228012e101),
           (tx: 234.5;   ty: 0.31622776601683793320e235),
           (tx: 2000.5;  ty: 0.31622776601683793320e2001)
         );
const
  name = 'exp10';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := exp10(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
    x := -tx;
    y := 1/ty;
    z := exp10(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_exp10m1;
const
  NT=16;
  TData: array[1..NT] of TPair = (
           (tx:  0.0;     ty:  0.0),
           (tx: -1e-15;   ty: -0.2302585092994043033e-14),
           (tx:  1e-9;    ty:  0.2302585095644994741e-8),
           (tx: -2e-10;   ty: -0.4605170184927711746e-9),
           (tx:  5e-6;    ty:  0.1151299173895094496e-4),
           (tx:  0.125;   ty:  0.3335214321633240257),
           (tx: -0.125;   ty: -0.2501057906675441727),
           (tx:  0.5;     ty:  2.162277660168379332),
           (tx: -0.75;    ty: -0.8221720589961077199),
           (tx:  1.0;     ty:  9.0),
           (tx: -1.0;     ty: -0.9),
           (tx:  2.0;     ty:  99),
           (tx: -6.7;     ty: -0.9999998004737685031),
           (tx:  12.345;  ty:  0.2213094709604637722e13),
           (tx:  100.0;   ty:  1e100),
           (tx:  4568.75; ty:  0.5623413251903490803e4569)
          );
const
  name = 'exp10m1';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := exp10m1(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_expm1;
const
  NT=15;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0),
           (tx:  1e-16;  ty: 0.100000000000000005e-15),
           (tx: 1e-9;    ty: 0.10000000005e-8),
           (tx: 1e-5;    ty: 0.10000050000166667083e-4),
           (tx: 0.125;   ty: 0.1331484530668263168),
           (tx: 0.5;     ty: 0.6487212707001281468),
           (tx: 1.0;     ty: 1.7182818284590452354),
           (tx: 2.0;     ty: 6.3890560989306502272),
           (tx: -1e-16;  ty: -0.99999999999999995e-16),
           (tx: -1e-9;   ty: -0.99999999950e-9),
           (tx: -1e-5;   ty: -0.999995000016666625e-5),
           (tx: -0.125;  ty: -0.11750309741540459714),
           (tx: -0.5;    ty: -0.39346934028736657640),
           (tx: -1.0;    ty: -0.63212055882855767840),
           (tx: -2.0;    ty: -0.86466471676338730811)
        );
const
  name = 'expm1';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := expm1(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_exprel;
const
  NT=19;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx:  1e-16;  ty: 1.0000000000000000500),
           (tx: 1e-9;    ty: 1.0000000005000000002),
           (tx: 1e-5;    ty: 1.0000050000166667083),
           (tx: 0.125;   ty: 1.0651876245346105346),
           (tx: 0.5;     ty: 1.2974425414002562937),
           (tx: 1.0;     ty: 1.7182818284590452354),
           (tx: 2.0;     ty: 3.1945280494653251136),
           (tx: -1e-16;  ty: 0.99999999999999995000),
           (tx: -1e-9;   ty: 0.99999999950000000017),
           (tx: -1e-5;   ty: 0.99999500001666662500),
           (tx: -0.125;  ty: 0.94002477932323677708),
           (tx: -0.5;    ty: 0.78693868057473315279),
           (tx: -1.0;    ty: 0.63212055882855767840),
           (tx: -2.0;    ty: 0.43233235838169365405),
           (tx: -10.0;   ty: 0.099995460007023751515),
           (tx: -40.0;   ty: 0.024999999999999999894),
           (tx: -50.0;   ty: 0.020000000000000000000),
           (tx: -1000;   ty: 0.001000000000000000000)
        );
const
  name = 'exprel';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := exprel(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_expx2;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 301/3;   ty: 8.831600902148364246E+4371),
           (tx:-301/3;   ty: 1.132297542744194027E-4372),
           (tx: 151/6;   ty: 1.162078878445555997E+275),
           (tx:-151/6;   ty: 8.605267839801379908E-276),
           (tx: 31/3;    ty: 2.360476487288662176E+46),
           (tx:-31/3;    ty: 4.236432793908657177E-47),
           (tx: 4/3;     ty: 5.916693590664329887),
           (tx:-4/3;     ty: 1.690133154060660767E-1),
           (tx: 0.1;     ty: 1.010050167084168058),
           (tx:-0.1;     ty: 9.900498337491680536E-1)
        );
const
  name = 'expx2';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := expx2(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_expmx2h;
const
  NT=7;
  TData: array[1..NT] of TPair = (
           (tx:-901/6;   ty: 2.107097699916602611E-4897),
           (tx: 301/3;   ty: 1.064094705721344160E-2186),
           (tx:-151/6;   ty: 2.933473681457084405E-138),
           (tx: 31/3;    ty: 6.508788515467880200E-24),
           (tx:-4/3;     ty: 4.111122905071874359E-1),
           (tx: 0.1;     ty: 9.950124791926823134E-1),
           (tx: -0.01;   ty: 9.999500012499791669E-1)
        );
const
  name = 'expmx2h';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := expmx2h(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ln;
const
  NT=13;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;       ty: -11398.805384308300613),  {succx(0)=0.36451995318824746025e-4950}
           (tx: 1e-4900;   ty: -11282.666955670823852),
           (tx: 1e-100;    ty: -230.25850929940456840),
           (tx: 1e-16;     ty: -36.841361487904730944),
           (tx: 0.0078125; ty: -4.8520302639196171659),
           (tx: 0.125;     ty: -2.0794415416798359283),
           (tx: 0.5;       ty: -0.69314718055994530942),
           (tx: 1.0;       ty: 0.0),
           (tx: 2.5;       ty: 0.91629073187415506518),
           (tx: 100.0;     ty: 4.6051701859880913680),
           (tx: 1e22;      ty: 50.656872045869005048),
           (tx: 1e100;     ty: 230.25850929940456840),
           (tx: 1e4931;    ty: 11354.047093553639268)
        );
const
  name = 'ln';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=0.0 then x := succx(0);
    y := ty;
    z := ln(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ln1p;
const
  NT=15;
  TData: array[1..NT] of TPair = (
           (tx: 0.0                      ; ty: 0.0   ),
           (tx: 0.100000000000000005e-15 ; ty: 1e-16 ),
           (tx: 0.10000000005e-8         ; ty: 1e-9  ),
           (tx: 0.10000050000166667083e-4; ty: 1e-5  ),
           (tx: 0.1331484530668263168    ; ty: 0.125 ),
           (tx: 0.6487212707001281468    ; ty: 0.5   ),
           (tx: 1.7182818284590452354    ; ty: 1.0   ),
           (tx: 6.3890560989306502272    ; ty: 2.0   ),
           (tx: -0.99999999999999995e-16 ; ty: -1e-16),
           (tx: -0.99999999950e-9        ; ty: -1e-9 ),
           (tx: -0.999995000016666625e-5 ; ty: -1e-5 ),
           (tx: -0.11750309741540459714  ; ty: -0.125),
           (tx: -0.39346934028736657640  ; ty: -0.5  ),
           (tx: -0.63212055882855767840  ; ty: -1.0  ),
           (tx: -0.86466471676338730811  ; ty: -2.0  )
        );
const
  name = 'ln1p';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := ln1p(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ln1pmx;
const
  NT=18;
  TData: array[1..NT] of TPair = (
           (tx: 0.0   ; ty:  0.0                      ),
           (tx: 1e-20 ; ty: -0.5e-40                  ),
           (tx: 1e-10 ; ty: -0.4999999999666666667e-20),
           (tx: 1e-5  ; ty: -0.4999966666916664667e-10),
           (tx: 0.125 ; ty: -0.7216964343616545461e-2 ),
           (tx: 0.5   ; ty: -0.9453489189183561802e-1 ),
           (tx: 0.75  ; ty: -0.1903842120645773137    ),
           (tx: 1.0   ; ty: -0.3068528194400546906    ),
           (tx: 2.0   ; ty: -0.9013877113318903086    ),
           (tx:10.0   ; ty: -7.6021047272016294559    ),
           (tx: -1e-20; ty: -0.5e-40                  ),
           (tx: -1e-10; ty: -0.5000000000333333333e-20),
           (tx: -1e-5 ; ty: -0.5000033333583335333e-10),
           (tx: -0.125; ty: -0.8531392624522623146e-2 ),
           (tx: -0.5  ; ty: -0.1931471805599453094    ),
           (tx: -0.75 ; ty: -0.6362943611198906188    ),
           (tx: -0.875; ty: -1.2044415416798359283    ),
           (tx: -0.9990234375; ty: -5.9324483680994530942)
        );
const
  name = 'ln1pmx';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := ln1pmx(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ln1mexp;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx: -1e-100; ty: -230.2585092994045684      ),
           (tx: -1e-20;  ty: -46.05170185988091368      ),
           (tx: -1e-10;  ty: -23.02585092999045684      ),
           (tx: -0.25;   ty: -1.508691549446032134      ),
           (tx: -1;      ty: -0.4586751453870818911     ),
           (tx: -10;     ty: -0.4540096037048920950e-4  ),
           (tx: -25;     ty: -0.1388794386506045809e-10 ),
           (tx: -44.5;   ty: -0.4719495271526123416e-19 ),
           (tx: -50;     ty: -0.1928749847963917783e-21 ),
           (tx: -100;    ty: -0.3720075976020835963e-43 ),
           (tx: -702;    ty: -0.1334362117671115081e-304)
        );
const
  name = 'ln1mexp';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := ln1mexp(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ln1pexp;
const
  NT=12;
  TData: array[1..NT] of TPair = (
           (tx: -500;   ty: 0.7124576406741285532e-217),
           (tx: -100;   ty: 0.3720075976020835963e-43),
           (tx: -20;    ty: 0.2061153620314380703e-8),
           (tx: -10;    ty: 0.4539889921686464677e-4),
           (tx: -1;     ty: 0.3132616875182228340),
           (tx: -1e-6;  ty: 0.6931466805600703094),
           (tx: 0;      ty: 0.6931471805599453094),
           (tx: 1;      ty: 1.3132616875182228340),
           (tx: 20;     ty: 20.000000002061153620),
           (tx: 40;     ty: 40.000000000000000004),
           (tx: 50;     ty: 50.0),
           (tx: 700;    ty: 700.0)
        );
const
  name = 'ln1pexp';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := ln1pexp(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_log_add_sub_exp;
const
  name = 'logadd/subexp';
var
  cnt, failed: integer;
  x,y,r,f: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;

  x := -50;
  y := -52;
  r := logaddexp(x,y);
  f := -49.87307198895702750;
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;
  r := logsubexp(x,y);
  f := -50.14541345786885906;
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;

  x := -200;
  y := -210;
  r := logaddexp(x,y);
  f := -199.9999546011007831;
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;
  r := logsubexp(x,y);
  f := -200.0000454009603705;
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;

  x := ln(2.5e-50);
  y := ln(1e-50);
  r := logaddexp(x,y);
  f := ln(3.5e-50);
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;
  r := logsubexp(x,y);
  f := ln(1.5e-50);
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;

  x := -15000;
  y := -15001;
  r := logaddexp(x,y);
  f := -14999.68673831248178;
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;
  r := logsubexp(x,y);
  f := -15000.45867514538708;
  inc(cnt);
  if reldevx(r,f) > EPS then begin
    inc(failed);
    writeln(cnt:3,' ',reldevx(r,f):18);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lncosh;
const
  NT=13;
  TData: array[1..NT] of TPair = (
          (tx: 12345.0;    ty: 12344.3068528194400546906),
          (tx: 1000.0;     ty: 999.3068528194400546906),
          (tx: -100.0;     ty: 99.3068528194400546906),
          (tx: 20.0;       ty: 19.3068528194400546948),
          (tx: -10.0;      ty: 9.3068528215012083109),
          (tx: 5.0;        ty: 4.3068982183392715552),
          (tx: 1.0;        ty: 0.4337808304830271870),
          (tx: 1.125;      ty: 0.5320593783568019023),
          (tx: 0.125;      ty: 0.7792239318898252791e-2),
          (tx: 0.00390625; ty: 0.7629375128775311005e-5),
          (tx: 1e-6;       ty: 0.4999999999999166667e-12),
          (tx: 1e-9;       ty: 0.4999999999999999999e-18),
          (tx: 1e-10;      ty: 0.5e-20));
const
  name = 'lncosh';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := lncosh(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lnsinh;
const
  NT=13;
  TData: array[1..NT] of TPair = (
          (tx: 12345.0;    ty:  12344.306852819440055),
          (tx: 1000.0;     ty:  999.3068528194400547),
          (tx: 100.0;      ty:  99.30685281944005469),
          (tx: 20.0;       ty:  19.30685281944005469),
          (tx: 10.0;       ty:  9.306852817378901066),
          (tx: 5.0;        ty:  4.306807418479684201),
          (tx: 1.0;        ty:  0.1614393615711956336),
          (tx: 1.125;      ty:  0.3204750982550105431),
          (tx: 0.125;      ty: -2.076838730005977444),
          (tx: 0.00390625; ty: -5.545174901349345561),
          (tx: 1e-6;       ty: -13.81551055796410744),
          (tx: 1e-9;       ty: -20.72326583694641116),
          (tx: 1e-10;      ty: -23.02585092994045684));
const
  name = 'lnsinh';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := lnsinh(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lnxp1;
const
  NT=15;
  TData: array[1..NT] of TPair = (
           (tx: 0.0                      ; ty: 0.0   ),
           (tx: 0.100000000000000005e-15 ; ty:  1e-16),
           (tx: 0.10000000005e-8         ; ty: 1e-9  ),
           (tx: 0.10000050000166667083e-4; ty: 1e-5  ),
           (tx: 0.1331484530668263168    ; ty: 0.125 ),
           (tx: 0.6487212707001281468    ; ty: 0.5   ),
           (tx: 1.7182818284590452354    ; ty: 1.0   ),
           (tx: 6.3890560989306502272    ; ty: 2.0   ),
           (tx: -0.99999999999999995e-16 ; ty: -1e-16),
           (tx: -0.99999999950e-9        ; ty: -1e-9 ),
           (tx: -0.999995000016666625e-5 ; ty: -1e-5 ),
           (tx: -0.11750309741540459714  ; ty: -0.125),
           (tx: -0.39346934028736657640  ; ty: -0.5  ),
           (tx: -0.63212055882855767840  ; ty: -1.0  ),
           (tx: -0.86466471676338730811  ; ty: -2.0  )
        );
const
  name = 'lnxp1';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := lnxp1(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_log10;
const
  NT=9;
  TData: array[1..NT] of TPair = (
           (tx: 3e-4000; ty: -3999.5228787452803376),
           (tx: 1.2e-5;  ty: -4.9208187539523751723),
           (tx: 0.125;   ty: -0.90308998699194358564),
           (tx: 1.0;     ty: 0.0),
           (tx: 2.0;     ty: 0.30102999566398119521),
           (tx: 12.0;    ty: 1.0791812460476248277),
           (tx: 123.0;   ty: 2.0899051114393979318),
           (tx: 1234e567;ty: 570.09131515969722288),
           (tx: 8.9e4090;ty: 4090.9493900066449128)
         );
const
  name = 'log10';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := log10(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
    x := 1.0/tx;
    y := -ty;
    z := log10(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
 end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_log2;
const
  NT=9;
  TData: array[1..NT] of TPair = (
           (tx: 3e-4000; ty: -13286.127417048728235),
           (tx: 1.2e-5;  ty: -16.346606068603017906),
           (tx: 0.125;   ty: -3.0),
           (tx: 1.0;     ty: 0.0),
           (tx: 2.4;     ty: 1.2630344058337938336),
           (tx: 12.3;    ty: 3.6205864104518775268),
           (tx: 123.4;   ty: 6.9471985842620555400),
           (tx: 1234e567;ty: 1893.8023564802838691),
           (tx: 8.9e4090;ty: 13589.839713425391038)
         );
const
  name = 'log2';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := log2(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
    x := 1.0/tx;
    y := -ty;
    z := log2(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_log2p1;
const
  NT=15;
  TData: array[1..NT] of TPair = (
           (tx: 0.0   ; ty:  0.0),
           (tx: 1e-20 ; ty:  0.1442695040888963407e-19),
           (tx: 1e-10 ; ty:  0.1442695040816828655e-9),
           (tx: 2e-5  ; ty:  0.2885361228261821942e-4),
           (tx: 0.125 ; ty:  0.1699250014423123629),
           (tx: 0.5   ; ty:  0.5849625007211561815),
           (tx: 1.0   ; ty:  1.0),
           (tx: 2.0   ; ty:  1.584962500721156181),
           (tx:10.0   ; ty:  3.459431618637297256),
           (tx: -2e-10; ty: -0.28853900820664658229e-9),
           (tx: -1e-5 ; ty: -0.14427022544122580475e-4),
           (tx: -0.125; ty: -0.19264507794239589256),
           (tx: -0.5  ; ty: -1.0),
           (tx: -0.75 ; ty: -2.0),
           (tx: -0.99 ; ty: -6.643856189774724696)
        );
const
  name = 'log2p1';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := log2p1(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_log10p1;
const
  NT=15;
  TData: array[1..NT] of TPair = (
           (tx: 0.0   ; ty:  0.0),
           (tx: 1e-20 ; ty:  0.4342944819032518276e-20),
           (tx: 2e-10 ; ty:  0.8685889637196447589e-10),
           (tx: 1.5e-5; ty:  0.6514368370908139095e-5),
           (tx: 0.125 ; ty:  0.5115252244738128895e-1),
           (tx: 0.5   ; ty:  0.1760912590556812421),
           (tx: 1.0   ; ty:  0.3010299956639811952),
           (tx: 2.0   ; ty:  0.4771212547196624373),
           (tx:10.0   ; ty:  1.041392685158225041),
           (tx: -1e-10; ty: -0.4342944819249665517e-10),
           (tx: -2e-5 ; ty: -0.8685976498119553194e-5),
           (tx: -0.125; ty: -0.5799194697768675493e-1),
           (tx: -0.5  ; ty: -0.3010299956639811952),
           (tx: -0.75 ; ty: -0.6020599913279623904),
           (tx: -0.9990234375; ty: -3.010299956639811952)
        );
const
  name = 'log10p1';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := log10p1(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_logit;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx: 1e-20;      ty: -46.05170185988091368),
           (tx: 1e-9;       ty: -20.72326583594641116),
           (tx: 1e-3;       ty: -6.906754778648553519),
           (tx: 0.125;      ty: -1.945910149055313305),
           (tx: 0.25;       ty: -1.098612288668109691),
           (tx: 1/3;        ty: -0.6931471805599453094),  {Fix311}
           (tx: 0.5;        ty: 0.0),
           (tx: 0.501953125;ty: 0.7812539736793652106e-2),
           (tx: 0.6875;     ty: 0.7884573603642701695),
           (tx: 0.99;       ty: 4.595119850134589927),
           (tx: 0.999996185302734375;   ty:12.47664543537447397));
const
  name = 'logit';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := logit(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_logistic;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx: -50;  ty: 0.1928749847963917783e-21),
           (tx: -20;  ty: 0.2061153618190203581e-8),
           (tx: -5;   ty: 0.6692850924284855559e-2),
           (tx: -1;   ty: 0.2689414213699951207),
           (tx: -0.5; ty: 0.3775406687981454354),
           (tx: 0;    ty: 0.5),
           (tx: 0.5;  ty: 0.6224593312018545646),
           (tx: 1;    ty: 0.7310585786300048793),
           (tx: 2;    ty: 0.8807970779778824441),
           (tx: 8;    ty: 0.9996646498695335219),
           (tx: 20;   ty: 0.9999999979388463818));
const
  name = 'logistic';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := logistic(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_polevl;
var
  i,f,cnt,failed: integer;
  ad,ax,a,x: extended;

  {arccosh rational approximations from Cephes, 1 <= x <= 1.5}
  function acosh0x(x: extended): extended;
  const
    p1: array[0..5] of THexExtW = (
         ($0000,$0000,$0000,$8000,$3fff),  {9.9999999999999999999715E-1}
         ($2364,$c41b,$1891,$8cab,$3fff),  {1.0989714347599256302467E0 }
         ($0fd5,$937f,$0515,$d4ed,$3ffd),  {4.1587081802731351459504E-1}
         ($7e8b,$8645,$341f,$8118,$3ffb),  {6.3034445964862182128388E-2}
         ($23a5,$f9aa,$289c,$d7a7,$3ff6),  {3.2906030801088967279449E-3}
         ($4536,$4dba,$9f55,$f3df,$3fef)); {2.9071989653343333587238E-5}
    q1: array[0..5] of THexExtW = (
         ($0000,$0000,$0000,$8000,$3fff),  {1.0000000000000000000028E0 }
         ($cdf2,$6ec5,$c33c,$9755,$3fff),  {1.1823047680932589605190E0 }
         ($e1e0,$c97c,$573a,$fdc5,$3ffd),  {4.9564621536841869854584E-1}
         ($2f83,$9e5e,$80af,$b3b6,$3ffb),  {8.7750439986662958343370E-2}
         ($c319,$c272,$a90a,$c4e3,$3ff7),  {6.0085845375571145826908E-3}
         ($1e7c,$4f16,$e98c,$db03,$3ff1)); {1.0443462486787584738322E-4}
  var
    w,z: extended;
    p: array[0..5] of extended absolute p1;
    q: array[0..5] of extended absolute q1;
  begin
    x := x-1.0;
    w := PolEvalX(x,p,6);
    z := PolEvalX(x,q,6);
    acosh0x := sqrt(2.0*x)*w/z;
  end;

  function acosh0d(x: extended): extended;
  const
    p: array[0..4] of double = (
         1.10855947270161294369E5,
         1.08102874834699867335E5,
         3.43989375926195455866E4,
         3.94726656571334401102E3,
         1.18801130533544501356E2);
    q: array[0..5] of double = (
         7.83869920495893927727E4,
         8.29725251988426222434E4,
         2.97683430363289370382E4,
         4.15352677227719831579E3,
         1.86145380837903397292E2,
         1.00000000000000000000E0);
  var
    w,z: extended;
  begin
    x := x-1.0;
    w := PolEval(x,p,5);
    z := PolEval(x,q,6);
    acosh0d := sqrt(x)*w/z;
  end;

begin
  writeln('Test PolEval/x');
  cnt := 0;
  failed := 0;
  i := 1;
  f := 0;
  while (i<=1000) and (f=0) do begin
    x  := 1.0 + 0.5*random;
    a  := arccosh(x);
    ad := acosh0d(x);
    ax := acosh0x(x);
    inc(cnt,2);
    if abs(a-ax)>4*EPSX*a then begin
      inc(failed);
      f := 1;
      writeln('PolEvalX failed. x=',x:24:17,'   diff=', abs(a-ax));
    end;
    if abs(a-ad)>4*EPSD*a then begin
      inc(failed);
      f := 1;
      writeln('PolEval  failed. x=',x:24:17,'   diff=', abs(a-ad));
    end;
    inc(i);
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_polevale;
const
  d08: array[0..08] of double = (1,-8,28,-56,70,-56,28,-8,1);
  d21: array[0..21] of double = (1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,
                                 352716,293930,203490,116280,54264,20349,5985,1330,210,21,1);
const
  x08: array[0..08] of extended = (1,-8,28,-56,70,-56,28,-8,1);
  x21: array[0..21] of extended = (1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,
                                   352716,293930,203490,116280,54264,20349,5985,1330,210,21,1);
const
  xa08: array[0..08] of extended = (1,8,28,56,70,56,28,8,1);
var
  i,cnt,failed,n,f: integer;
  yx,ex,y,x,ax,ps: extended;
  yd,ed,ad: double;
begin
  writeln('Test PolEvalEE/X');
  cnt := 0;
  failed := 0;
  i := 1;
  while i<=1000 do begin
    if odd(i) then begin
      x  := 1 + (random(128)-256.0)/1024.0;
      y  := power(x-1,8);
      yx := PolEvalEEX(x,x08,9,ex);
      yd := PolEvalEE(x,d08,9,ed);
      ps := PolEvalX(abs(x),xa08,9);
      n  := 8;
    end
    else begin
      x  := 1 - (random(128)-256.0)/1024.0;
      y  := power(x+1,21);
      yx := PolEvalEEX(x,x21,22,ex);
      yd := PolEvalEE(x,d21,22,ed);
      ps := PolEvalX(abs(x),x21,22);
      n  := 21;
    end;
    {ax, ad are the approx.  a priori bounds}
    {ex, ed are the estimated running bounds}
    ax := ps*2*n/(1-2*n*eps_x);
    ad := eps_d*ax;
    ax := eps_x*ax;
    f  := 0;
    inc(cnt,2);
    if abs(y-yx) > ax then begin
      inc(failed);
      inc(f);
      writeln('PolEvalEEX failed. i=',i:3,'   ax=',ax:24:17,'   diff=', abs(y-yx));
    end;
    if abs(y-yd) > ad then begin
      inc(failed);
      inc(f);
      writeln('PolEvalEE  failed. i=',i:3,'   ad=',ad:24:17,'   diff=', abs(y-yd));
    end;
    if f=0 then begin
      if abs(y-yx)>ex then begin
        inc(failed);
        writeln('PolEvalEEX failed. i=',i:3,'   ex=',ex:24,'   diff=', abs(y-yx));
      end;
      if abs(y-yd)>ed then begin
        inc(failed);
        writeln('PolEvalEE  failed. i=',i:3,'   ed=',ed:24,'   diff=', abs(y-yd));
      end;
    end;
    inc(i);
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_polevalche;
const
  d08: array[0..08] of double = (1,-8,28,-56,70,-56,28,-8,1);
  d21: array[0..21] of double = (1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,
                                   352716,293930,203490,116280,54264,20349,5985,1330,210,21,1);
const
  x08: array[0..08] of extended = (1,-8,28,-56,70,-56,28,-8,1);
  x21: array[0..21] of extended = (1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,
                                   352716,293930,203490,116280,54264,20349,5985,1330,210,21,1);
const
  xa08: array[0..08] of extended = (1,8,28,56,70,56,28,8,1);
  da08: array[0..08] of double   = (1,8,28,56,70,56,28,8,1);

const fx08: array[-4..4] of THexExtW = {  (x-1)^8}
  { 0.5555555556}   (($EBAA,$1463,$920F,$C78C,$3FF5),
  { 0.6666666667}    ($A6E2,$A3CC,$CD5A,$9FD1,$3FF2),
  { 0.7777777778}    ($EBB8,$1463,$920F,$C78C,$3FED),
  { 0.8888888889}    ($EBB8,$1463,$920F,$C78C,$3FE5),
  { 1.0000000000}    ($0000,$0000,$0000,$0000,$0000),
  { 1.1111111111}    ($EBB8,$1463,$920F,$C78C,$3FE5),
  { 1.2222222222}    ($EBB8,$1463,$920F,$C78C,$3FED),
  { 1.3333333333}    ($A6E2,$A3CC,$CD5A,$9FD1,$3FF2),
  { 1.4444444444}    ($EBB8,$1463,$920F,$C78C,$3FF5));
const fd08: array[-4..4] of THexDblW = {(d-1)^8}
  { 0.5555555556}   (($8C7A,$41E2,$F192,$3F58),
  { 0.6666666667}    ($799A,$AB54,$FA39,$3F23),
  { 0.7777777778}    ($8C7A,$41E2,$F192,$3ED8),
  { 0.8888888889}    ($8C96,$41E2,$F192,$3E58),
  { 1.0000000000}    ($0000,$0000,$0000,$0000),
  { 1.1111111111}    ($8C96,$41E2,$F192,$3E58),
  { 1.2222222222}    ($8C96,$41E2,$F192,$3ED8),
  { 1.3333333333}    ($798B,$AB54,$FA39,$3F23),
  { 1.4444444444}    ($8C7A,$41E2,$F192,$3F58));
const fx21: array[-4..4] of THexExtW = {(x+1)^21}
  {-1.4444444444}  (($74AA,$FC50,$5EB0,$ACA2,$BFE6),
  {-1.3333333333}   ($74BA,$A5CF,$7D6D,$D239,$BFDD),
  {-1.2222222222}   ($74AA,$FC50,$5EB0,$ACA2,$BFD1),
  {-1.1111111111}   ($74AA,$FC50,$5EB0,$ACA2,$BFBC),
  {-1.0000000000}   ($0000,$0000,$0000,$0000,$0000),
  {-0.8888888889}   ($74AA,$FC50,$5EB0,$ACA2,$3FBC),
  {-0.7777777778}   ($74AA,$FC50,$5EB0,$ACA2,$3FD1),
  {-0.6666666667}   ($74BA,$A5CF,$7D6D,$D239,$3FDD),
  {-0.5555555556}   ($748A,$FC50,$5EB0,$ACA2,$3FE6));
const fd21: array[-4..4] of THexDblW = {(d+1)^21}
  {-1.4444444444}  (($8A07,$D61F,$944B,$BE65),
  {-1.3333333333}   ($B9CC,$ADB4,$472F,$BDDA),
  {-1.2222222222}   ($8A47,$D61F,$944B,$BD15),
  {-1.1111111111}   ($8A47,$D61F,$944B,$BBC5),
  {-1.0000000000}   ($0000,$0000,$0000,$0000),
  {-0.8888888889}   ($8A47,$D61F,$944B,$3BC5),
  {-0.7777777778}   ($8A07,$D61F,$944B,$3D15),
  {-0.6666666667}   ($BA00,$ADB4,$472F,$3DDA),
  {-0.5555555556}   ($8A07,$D61F,$944B,$3E65));

var
  i,cnt,failed,f: integer;
  yx,ex,x,rx,px,cx,erx,x9: extended;
  yd,ed,d,rd,pd,cd,erd: double;
begin
  writeln('Test PolEvalCHE/X');
  cnt := 0;
  failed := 0;

  {cx,cd are the values p~(x) = |a[0]| + |a[1]*x| +...+ |a[n-1]*x^(n-1)| }
  {used for the a priori bounds |p(x)-y| = eps*|p(x)| + (2n*eps)^2*p~(x) }
  {Failure, if the a priori bounds are violated! These calculations are  }
  {very sensitive due to the nearly doubled precision! Therefore be sure }
  {to calculate with the multi-precision values of the actual used double}
  {or extended arguments and compare the results to the hex MP values!   }
  {Note: With round to nearest, u=eps/2 could be used in the formulas.   }

  x9 := 9.0;
  {Note D2/3 are too clever and get a wrong correctly rounded result for }
  {x := +-1.0 + i/9 if i = +-3! Therefore the x9 variable is used!  This }
  {shows again the sensibility and/or the tightness of the error bounds. }
  {Compare the same false D2/3 compile time optimisation in sfc_ellint_3!}
  for i:=-4 to 4 do begin
    x  := 1.0 + i/x9;
    yx := PolEvalCHEX(x,x08,9,erx);
    cx := PolEvalX(abs(x),xa08,9);
    px := extended(fx08[i]);
    rx := abs(yx-px);
    ex := eps_x*abs(px) + sqr(16*eps_x)*cx;
    f  := 0;
    inc(cnt);
    if rx > ex then begin
      writeln('Failed a priori x08, i=',i,rx:30, ex:30);
      inc(failed);
    end;
    if rx > erx then begin
      writeln('Failed dynamic  x08, i=',i,rx:30, erx:30);
      if f=0 then inc(failed);
    end;
  end;
  for i:=-4 to 4 do begin
    d  := 1 + i/9;
    yd := PolEvalCHE(d,d08,9, erd);
    cd := PolEval(abs(d),da08,9);
    pd := double(fd08[i]);
    rd := abs(yd-pd);
    ed := eps_d*abs(pd) + sqr(16*eps_d)*cd;
    f  := 0;
    inc(cnt);
    if rd > ed then begin
      writeln('Failed a priori d08, i=',i,rd:30, ed:30);
      inc(failed);
    end;
    if rd > erd then begin
      writeln('Failed dynamic  d08, i=',i,rd:30, erd:30);
      if f=0 then inc(failed);
    end;
  end;
  for i:=-4 to 4 do begin
    x  := -1.0 + i/x9;
    yx := PolEvalCHEX(x,x21,22,erx);
    cx := PolEvalX(abs(x),x21,22);
    px := extended(fx21[i]);
    rx := abs(yx-px);
    ex := eps_x*abs(px) + sqr(42*eps_x)*abs(cx);
    f  := 0;
    inc(cnt);
    if rx > ex then begin
      writeln('Failed a priori x21, i=',i,rx:30, ex:30);
      inc(failed);
    end;
    if rx > erx then begin
      writeln('Failed dynamic  x21, i=',i,rx:30, erx:30);
      if f=0 then inc(failed);
    end;
  end;
  for i:=-4 to 4 do begin
    d  := -1 + i/9;
    yd := PolEvalCHE(d,d21,22,erd);
    cd := PolEval(abs(d),d21,22);
    pd := double(fd21[i]);
    rd := abs(yd-pd);
    ed := eps_d*abs(pd) + sqr(42*eps_d)*abs(cd);
    f  := 0;
    inc(cnt);
    if rd > ed then begin
      writeln('Failed a priori d21, i=',i,rd:30, ed:30);
      inc(failed);
    end;
    if rd > erd then begin
      writeln('Failed dynamic  d21, i=',i,rd:30, erd:30);
      if f=0 then inc(failed);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sqrt1pmx;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx:      0; ty: 1.0),
           (tx:  0.125; ty: 0.8827822185373187065),
           (tx: -0.125; ty: 1.132782218537318707),
           (tx:      1; ty: 0.4142135623730950488),
           (tx:     -1; ty: 2.414213562373095049),
           (tx:     10; ty: 0.4987562112089027022e-1),
           (tx:    -10; ty: 20.04987562112089027),
           (tx:    100; ty: 0.4999875006249609402e-2),
           (tx:   -100; ty: 200.0049998750062496),
           (tx:    1e6; ty: 0.49999999999987500e-6),
           (tx:   -1e6; ty: 2000000.000000500000));
const
  name = 'sqrt1pmx';
var
  i, cnt, failed: integer;
  x,y,z: extended;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := sqrt1pmx(x);
    inc(cnt);
    if reldevx(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldevx(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_trig_degree;
var
  i,cnt,failed: integer;
  x,y,t: extended;

  procedure test1(n: integer; f,v: extended);
  begin
    inc(cnt);
    if f<>v then begin
      inc(failed);
      writeln('Test ',n, ' failed: ', f:30, v:30);
    end;
  end;
begin
  writeln('Degrees versions of trig/invtrig functions');
  cnt := 0;
  failed := 0;

  test1(01, sind(0), 0);
  test1(02, sind(90), 1);
  test1(03, sind(180), 0);
  test1(04, sind(270), -1);
  test1(05, sind(360), 0);

  test1(06, cosd(0), 1);
  test1(07, cosd(90), 0);
  test1(08, cosd(180), -1);
  test1(09, cosd(270), 0);
  test1(10, cosd(360), 1);

  test1(11, tand(0), 0);
  test1(12, tand(45), 1);
  test1(13, tand(90), PosInf_x);
  test1(14, tand(135), -1);
  test1(15, tand(180), 0);

  test1(16, cotd(0), PosInf_x);
  test1(17, cotd(45), 1);
  test1(18, cotd(90), 0);
  test1(19, cotd(135), -1);
  test1(20, cotd(180), PosInf_x);

  test1(21, arccosd(1) , 0);
  test1(22, arccosd(0) , 90);
  test1(23, arccosd(-1) , 180);

  test1(24, arcsind(1) , 90);
  test1(25, arcsind(0) , 0);
  test1(26, arcsind(-1) , -90);

  test1(27, arctand(1) , 45);
  test1(28, arctand(0) , 0);
  test1(29, arctand(-1) , -45);
  test1(30, arctand(PosInf_x) , 90);
  test1(31, arctand(NegInf_x) , -90);

  test1(32, arccotd(1) , 45);
  test1(33, arccotd(0) , 90);
  test1(34, arccotd(-1) , -45);
  test1(35, arccotd(PosInf_x) , 0);
  test1(36, arccotd(NegInf_x) , 0);

  test1(37, arccotcd(1) , 45);
  test1(38, arccotcd(0) , 90);
  test1(39, arccotcd(-1) , 135);
  test1(40, arccotcd(PosInf_x) , 0);
  test1(41, arccotcd(NegInf_x) , 180);

  for i:=-20 to 20 do begin
    x := i*90.0;
    y := sind(x);
    case abs(i) and 3 of
        0: t := 0;
        1: t := 1;
        2: t := 0;
      else t := -1;
    end;
    if i<0 then t := -t;
    test1(100+i, y , t);
  end;

  for i:=-20 to 20 do begin
    x := i*90.0;
    y := cosd(x);
    case abs(i) and 3 of
        0: t := 1;
        1: t := 0;
        2: t := -1;
      else t := 0;
    end;
    test1(200+i, y , t);
  end;

  for i:=-20 to 20 do begin
    x := i*45.0;
    y := tand(x);
    case abs(i) and 3 of
        0: t := 0;
        1: t := 1;
        2: t := PosInf_x;
      else t := -1;
    end;
    if (i<0) and odd(i) then t := -t;
    test1(300+i, y , t);
  end;

  for i:=-20 to 20 do begin
    x := i*45.0;
    y := cotd(x);
    case abs(i) and 3 of
        0: t := PosInf_x;
        1: t := 1;
        2: t := 0
      else t := -1;
    end;
    if (i<0) and odd(i) then t := -t;
    test1(400+i, y , t);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_compound;
var
  a,b,x,sm: extended;
  y: longint;
  cnt, failed: integer;

  procedure itest;
  begin
    a := compound(x,y);
    inc(cnt);
    if reldevx(a,b) > SEPS then begin
      inc(failed);
      writeln(' compound: ',x:5:2,' ',y:5, ' ',a:19,' ',b:19,' ',reldevx(a,b):16);
    end;
  end;

begin
  writeln('Function compound');
  cnt :=0;
  failed := 0;
  sm  := ldexp(1,-30);

  x := -1.5;
  y := -2;  b :=  4.000000000000000000;      itest;
  y := -1;  b := -2.000000000000000000;      itest;
  y :=  0;  b :=  1.0;                       itest;
  y :=  1;  b := -0.5;                       itest;
  y :=  2;  b :=  0.25;                      itest;

  x := -1.0 + sm;
  y := -2;   b := 0.1152921504606846976e19;  itest;
  y := -5;   b := 0.1427247692705959881e46;  itest;
  y := -1;   b := 1073741824.000000000;      itest;
  y :=-30;   b := 0.8452712498170643942e271; itest;
  y :=  0;   b := 1.0;                       itest;
  y := 30;   b := 0.1183052186166774711e-270;itest;
  y :=  1;   b := 0.9313225746154785156e-9;  itest;
  y :=  5;   b := 0.7006492321624085355e-45; itest;
  y :=  2;   b := 0.8673617379884035472e-18; itest;

  x := -sm;
  y := -2;   b := 1.000000001862645152;   itest;
  y := -5;   b := 1.000000004656612886;   itest;
  y := -1;   b := 1.000000000931322575;   itest;
  y := -200; b := 1.000000186264532357;   itest;
  y :=  0;   b := 1.0;                    itest;
  y := 200;  b := 0.9999998137355023374;  itest;
  y :=  1;   b := 0.9999999990686774254;  itest;
  y :=  5;   b := 0.9999999953433871356;  itest;
  y :=  2;   b := 0.9999999981373548516;  itest;

  x :=  0.0;
  y := -2;   b := 1.0;   itest;
  y := -5;   b := 1.0;   itest;
  y := -1;   b := 1.0;   itest;
  y := -200; b := 1.0;   itest;
  y :=  0;   b := 1.0;   itest;
  y :=  200; b := 1.0;   itest;
  y :=  1;   b := 1.0;   itest;
  y :=  5;   b := 1.0;   itest;
  y :=  2;   b := 1.0;   itest;

  x := +sm;
  y := -2;   b := 0.9999999981373548534;  itest;
  y := -5;   b := 0.9999999953433871399;  itest;
  y := -1;   b := 0.9999999990686774263;  itest;
  y := -200; b := 0.9999998137355025109;  itest;
  y :=  0;   b := 1.0;                    itest;
  y := 200;  b := 1.000000186264532184;   itest;
  y :=  1;   b := 1.000000000931322575;   itest;
  y :=  5;   b := 1.000000004656612882;   itest;
  y :=  2;   b := 1.000000001862645150;   itest;

  x := 1.0-sm;
  y := -2;   b := 0.2500000002328306438;  itest;
  y := -5;   b := 0.3125000007275957624e-1;  itest;
  y := -1;   b := 0.5000000002328306438;  itest;
  y := -200; b := 0.6223015857424629875e-60;  itest;
  y :=  0;   b := 1.0;                    itest;
  y := 200;  b := 0.1606937894601229547e61;   itest;
  y :=  1;   b := 1.999999999068677425;   itest;
  y :=  5;   b := 31.99999992549419410;   itest;
  y :=  2;   b := 3.999999996274709702;   itest;

  x := 1e-9;
  y := 10000;
  b := 1.000010000049995167;
  itest;

  x := 1e-15;
  y := 10000;
  b := 1.000000000010000000;
  itest;

  x := 1e-9;
  y := 2000000000;
  b := 7.389056091541594137;
  itest;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_comprel;
var
  cnt, failed: integer;

  procedure crtest(x,y,f: extended);
  var
    a: extended;
  begin
    a := comprel(x,y);
    inc(cnt);
    if reldevx(a,f) > SEPS then begin
      inc(failed);
      writeln(' comprel: ',x:12:10,' ',y:8:0, ' ',a:24,' ',reldevx(a,f):18);
    end;
  end;

begin
  writeln('Function comprel');
  cnt :=0;
  failed := 0;

  crtest(1/1024, 100, 104.9919149837248596);
  crtest(1/1024, 10000, 17756707.13845046449);
  crtest(-1, 10, 1);
  crtest(0, 10, 10);
  crtest(2, 10, 29524);
  crtest(-2.5, 100, -162624471014086094.6);
  crtest(0.125, -100, -7.999938646726101872);
  crtest(-0.125, -100, -5038299.671174422405);
  crtest(1e-9, 1000000, 1000500.1662078414183);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_pow1p;
var
  a,b,x,y,sm: extended;
  cnt, failed: integer;

  procedure itest;
  begin
    a := pow1p(x,y);
    inc(cnt);
    if reldevx(a,b) > SEPS then begin
      inc(failed);
      writeln(' pow1pm1: ',x:5:2,' ',y:5:2, ' ',a:19,' ',b:19,' ',reldevx(a,b):16);
    end;
  end;

begin
  writeln('Function pow1p');
  cnt :=0;
  failed := 0;
  sm  := ldexp(1,-30);

  x := -0.5;
  y := -2.0;  b := 4.000000000000000000;      itest;
  y := -1.0;  b := 2.000000000000000000;      itest;
  y :=  0.0;  b := 1.0;                       itest;
  y :=  1.0;  b := 0.5;                       itest;
  y :=  2.0;  b := 0.25;                      itest;

  x := -1.0 + sm;
  y := -2.0;  b := 0.1152921504606846976e19;  itest;
  y := -1.5;  b := 0.3518437208883200000e14;  itest;
  y := -1.0;  b := 1073741824.000000000;      itest;
  y := -sm;   b := 1.000000019366308691;      itest;
  y :=  0.0;  b := 1.0;                       itest;
  y :=  sm;   b := 0.9999999806336916839;     itest;
  y :=  1.0;  b := 0.9313225746154785156e-9;  itest;
  y :=  1.5;  b := 0.2842170943040400743e-13; itest;
  y :=  2.0;  b := 0.8673617379884035472e-18; itest;

  x := -sm;
  y := -2.0;  b := 1.000000001862645152;   itest;
  y := -1.5;  b := 1.000000001396983864;   itest;
  y := -1.0;  b := 1.000000000931322575;   itest;
  y := -sm;   b := 1.000000000000000001;   itest;
  y :=  0.0;  b := 1.0;                    itest;
  y :=  sm;   b := 0.9999999999999999991;  itest;
  y :=  1.0;  b := 0.9999999990686774254;  itest;
  y :=  1.5;  b := 0.9999999986030161384;  itest;
  y :=  2.0;  b := 0.9999999981373548516;  itest;

  x :=  0.0;
  y := -2.0;  b := 1.0;   itest;
  y := -1.5;  b := 1.0;   itest;
  y := -1.0;  b := 1.0;   itest;
  y :=  -sm;  b := 1.0;   itest;
  y :=  0.0;  b := 1.0;   itest;
  y :=   sm;  b := 1.0;   itest;
  y :=  1.0;  b := 1.0;   itest;
  y :=  1.5;  b := 1.0;   itest;
  y :=  2.0;  b := 1.0;   itest;

  x := +sm;
  y := -2.0;  b := 0.9999999981373548534;  itest;
  y := -1.5;  b := 0.9999999986030161397;  itest;
  y := -1.0;  b := 0.9999999990686774263;  itest;
  y := -sm;   b := 0.9999999999999999991;  itest;
  y :=  0.0;  b := 1.0;                    itest;
  y :=  sm;   b := 1.000000000000000001;   itest;
  y :=  1.0;  b := 1.000000000931322575;   itest;
  y :=  1.5;  b := 1.000000001396983862;   itest;
  y :=  2.0;  b := 1.000000001862645150;   itest;

  x := 1.0-sm;
  y := -2.0;  b := 0.2500000002328306438;  itest;
  y := -1.5;  b := 0.3535533908402279528;  itest;
  y := -1.0;  b := 0.5000000002328306438;  itest;
  y :=  -sm;  b := 0.9999999993544563839;  itest;
  y :=  0.0;  b := 1.0;                    itest;
  y :=   sm;  b := 1.000000000645543617;   itest;
  y :=  1.0;  b := 1.999999999068677425;   itest;
  y :=  1.5;  b := 2.828427122770556574;   itest;
  y :=  2.0;  b := 3.999999996274709702;   itest;

  x := 1e-9;
  y := 10000;
  b := 1.000010000049995167;
  itest;

  x := 1e-15;
  y := 10000;
  b := 1.000000000010000000;
  itest;

  x := 1e-9;
  y := 1e10;
  b := 22026.46568467438789;
  itest;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_powpi2k;
var
  cnt, failed: integer;

  procedure test1(k,n: longint; f: extended);
  var
    y: extended;
  {$ifdef FPC}
  const NE = 4;  {needed for FPC264/AMD, see below}
  {$else}
  const NE = 2;
  {$endif}
  begin
    y := powpi2k(k,n);
    inc(cnt);
    if reldevx(y,f) > NE*eps_x then begin
      inc(failed);
      writeln(' failed: ',k:5,' ',n:5, ' ',y:19,' ',f:19,' ',reldevx(y,f)/eps_x:6:2);
    end;
  end;

begin
  writeln('Function powpi2k');
  cnt :=0;
  failed := 0;

  test1(0,1,Pi);
  test1(0,2,PiSqr);
  test1(123,0,1);
  test1( -2, 47000, 0.1676579725613674408e-4930);
  test1(  1,  5432, 0.5164697875989590794e4336);
  test1( 10,  -321, 0.1284117556795506586e-1125);
  test1(123,    45, 0.3739127474315378958e1689);    {FPC264/AMD: 3.21 eps_x}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;

end.
