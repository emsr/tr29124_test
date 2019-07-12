{T_DAMATM - regression test main unit for DAMATH  (c) 2013-2016  W.Ehrhardt}

{Many test cases were calculated with Maple V R4, Digits:=20}
{others with Pari/GP 2.3.4 and \p 20, some with Cephes qcalc}

unit t_damatm;

interface

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}


procedure test_damath_main;
  {-DAMath regression test via single procedure call}

implementation

uses
  damath, t_damat0, t_damat1, t_damat2;

const
  Pi_1   = 3.1415926535897932385;
  Pi_2   = 1.5707963267948966192;
  Pi_3   = 1.0471975511965977462;
  Pi_4   = 0.78539816339744830962;
  Pi_6   = 0.52359877559829887308;


{---------------------------------------------------------------------------}
procedure test_arccos;
const
  NT=7;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: Pi_2   ),
           (tx: 1.0;     ty: 0.0    ),
           (tx: -1.0;    ty: Pi_1   ),
           (tx: sqrt_5;  ty: pi_4   ),
           (tx: -sqrt_5; ty: 2.3561944901923449288),
           (tx: 0.5;     ty: Pi_3   ),
           (tx: -0.5;    ty: 2.0943951023931954923)
         );
const
  name = 'arccos';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccos(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccos1m;
const
  NT=8;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0    ),
           (tx: 1e-20;   ty: 0.1414213562373095049e-9),
           (tx: 2e-9;    ty: 0.6324555321390851218e-4),
           (tx: 0.125;   ty: 0.5053605102841573070),
           (tx: 0.5;     ty: Pi_3 ),
           (tx: 0.75;    ty: 1.318116071652817966),
           (tx: 1;       ty: Pi_2 ),
           (tx: 2;       ty: Pi   )
         );
const
  name = 'arccos1m';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccos1m(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccosh;
const
  NT=6;
  TData: array[1..NT] of TPair = (
           (tx: 1.0;       ty: 0.0   ),
           (tx: 1.01;      ty: 0.1413037694856485773511516471),
           (tx: 2.0;       ty: 1.3169578969248167086),
           (tx: 10.0;      ty: 2.9932228461263808979),
           (tx: 1e3;       ty: 7.6009022095419886114),
           (tx: 1e300;     ty: 691.468675078773650514814668527)
         );
const
  name = 'arccosh';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccosh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccosh1p;
const
  NT=6;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;       ty: 0.0   ),
           (tx: 0.0000001; ty: 0.0004472135917731780606347295503),
           (tx: 0.0001;    ty: 0.01414201777525232424406347593),
           (tx: 0.01;      ty: 0.1413037694856485773511516471),
           (tx: 1.0;       ty: 1.3169578969248167086),
           (tx: 9.0;       ty: 2.9932228461263808979)
         );
const
  name = 'arccosh1p';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccosh1p(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccot;
const
  NT=7;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: Pi_2   ),
           (tx: 1.0;     ty: Pi_4   ),
           (tx: 2.0;     ty: 0.46364760900080611622),
           (tx: 10.0;    ty: 0.99668652491162027378e-1),
           (tx: 100.0;   ty: 0.99996666866652382063e-2),
           (tx: 1e4;     ty: 0.99999999666666668667e-4),
           (tx: 1e10;    ty: 1e-10)
         );
const
  name = 'arccot';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccot(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    if x<>0.0 then begin
      x := -tx;
      y := -ty;
      z := arccot(x);
      inc(cnt);
      if reldev(z,y) > EPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccotc;
const
  NT=7;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: Pi_2   ),
           (tx: 1.0;     ty: Pi_4   ),
           (tx: 2.0;     ty: 0.46364760900080611622),
           (tx: 10.0;    ty: 0.99668652491162027378e-1),
           (tx: 100.0;   ty: 0.99996666866652382063e-2),
           (tx: 1e4;     ty: 0.99999999666666668667e-4),
           (tx: 1e10;    ty: 1e-10)
         );
const
  name = 'arccotc';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccotc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    if x<>0.0 then begin
      x := -tx;
      y := Pi-ty;
      z := arccotc(x);
      inc(cnt);
      if reldev(z,y) > EPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccoth;
const
  NT=5;
  TData: array[1..NT] of TPair = (
           (tx: 1.00006103515625; ty: 5.1986191127558264138),
           (tx: 1.0009765625;     ty: 3.8125535741194498774),
           (tx: 1.01;             ty: 2.6516524540295378755),
           (tx: 1e4;              ty: 0.10000000033333333533e-3),
           (tx: 1e10;             ty: 1e-10)
         );
const
  name = 'arccoth';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccoth(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := arccoth(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccsc;
const
  NT=7;
  TData: array[1..NT] of TPair = (
           (tx: 1.0;     ty: Pi_2),
           (tx: 2.0;     ty: Pi_6),
           (tx: sqrt2;   ty: Pi_4),
           (tx: 10.0;    ty: 0.10016742116155979635),
           (tx: 100.0;   ty: 0.10000166674167113126e-1),
           (tx: 1e4;     ty: 0.10000000016666666742e-3),
           (tx: 1e10;    ty: 1e-10)
         );
const
  name = 'arccsc';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccsc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := arccsc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arccsch;
const
  NT=9;
  TData: array[1..NT] of TPair = (
           (tx: 1e-300;  ty: 691.468675078773650514814668527),
           (tx: 1e-5;    ty: 12.206072645555173730),
           (tx: 0.125;   ty: 2.7764722807237176735),
           (tx: 1.0;     ty: 0.88137358701954302523),
           (tx: 2.0;     ty: 0.48121182505960344750),
           (tx: 10.0;    ty: 0.99834078899207563327e-1),
           (tx: 100.0;   ty: 0.99998333408328869351e-2),
           (tx: 1e4;     ty: 0.99999999833333334083e-4),
           (tx: 1e10;    ty: 1e-10)
         );
const
  name = 'arccsch';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arccsch(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := arccsch(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arcsec;
const
  NT=6;
  TData: array[1..NT] of TPair = (
           (tx: 1.0;     ty: 0),
           (tx: 2.0;     ty: Pi_3),
           (tx: sqrt2;   ty: Pi_4),
           (tx: 100.0;   ty: 1.5607961601207295061),
           (tx: 1e4;     ty: 1.5706963267947299526),
           (tx: 1e10;    ty: 1.5707963266948966192)
         );
const
  name = 'arcsec';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arcsec(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := Pi-ty;
    z := arcsec(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arcsech;
const
  NT=8;
  TData: array[1..NT] of TPair = (
           (tx: 1e-305;  ty: 702.9816005437438789349),
           (tx: 1e-5;    ty: 12.206072645505173730),
           (tx: 0.125;   ty: 2.7686593833135738327),
           (tx: 0.5;     ty: 1.3169578969248167086),
           (tx: 0.75;    ty: 0.79536546122390563053),
           (tx: 0.99;    ty: 0.14201444074607706875),
           (tx: 0.999755859375; ty: 0.22099335098031397862e-1),
           (tx: 1.0;     ty: 0.0)
         );
const
  name = 'arcsech';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arcsech(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arcsin;
const
  NT=9;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0),
           (tx: 0.5;     ty: Pi_6),
           (tx: sqrt_5;  ty: Pi_4),
           (tx: 0.125;   ty: 0.12532783116806539687),
           (tx: 0.5;     ty: 0.52359877559829887308),
           (tx: 0.75;    ty: 0.84806207898148100805),
           (tx: 0.99;    ty: 1.4292568534704694005),
           (tx: 0.9999;  ty: 1.5566540733173837416),
           (tx: 1.0;     ty: Pi_2)
         );
const
  name = 'arcsin';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arcsin(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := arcsin(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arcsinh;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0),
           (tx: 1e-10;   ty: 1e-10),
           (tx: 1e-5;    ty: 0.99999999998333333333e-5),
           (tx: 0.125;   ty: 0.12467674692144274393),
           (tx: 1.0;     ty: 0.88137358701954302523),
           (tx: 2.0;     ty: 1.4436354751788103425),
           (tx: 10.0;    ty: 2.9982229502979697388),
           (tx: 100.0;   ty: 5.2983423656105887574),
           (tx: 1e10;    ty: 23.718998110500402150),
           (tx: 1e300;   ty: 691.4686750787736505148)
         );
const
  name = 'arcsinh';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arcsinh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := arcsinh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arctanh;
const
  NT=8;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0),
           (tx: 1e-100;  ty: 1e-100),
           (tx: 1e-5;    ty: 0.10000000000333333333e-4),
           (tx: 0.125;   ty: 0.12565721414045303884),
           (tx: 0.5;     ty: 0.54930614433405484570),
           (tx: 0.75;    ty: 0.97295507452765665255),
           (tx: 0.99;    ty: 2.6466524123622461977),
           (tx: 0.999755859375; ty: 4.5053956347578010201)
         );
const
  name = 'arctanh';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arctanh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := arctanh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cosh;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx: 1e-10;   ty: 1.0),
           (tx: 1e-5;    ty: 1.0000000000500000000),
           (tx: 0.125;   ty: 1.0078226778257108598),
           (tx: 0.5;     ty: 1.1276259652063807852),
           (tx: 1.0;     ty: 1.5430806348152437785),
           (tx: 2.0;     ty: 3.7621956910836314596),
           (tx: 10.0;    ty: 11013.232920103323140),
           (tx: 100.0;   ty: 0.13440585709080677242e44),
           (tx: 700;     ty: 0.50711602736750225473e304)
         );
const
  name = 'cosh';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := cosh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := ty;
    z := cosh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_coshm1;
const
  NT=12;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0),
           (tx: 1e-12;   ty: 0.5e-24),
           (tx: 1e-7;    ty: 0.50000000000000041667e-14),
           (tx: 1e-5;    ty: 0.50000000000416666667e-10),
           (tx: 0.03125; ty: 0.48832098772337639343e-3),
           (tx: 0.125;   ty: 0.78226778257108598469e-2),
           (tx: 0.5;     ty: 0.12762596520638078523),
           (tx: 1.0;     ty: 0.54308063481524377848),
           (tx: 2.0;     ty: 2.76219569108363145956),
           (tx: 10.0;    ty: 11012.2329201033231397),
           (tx: 100.0;   ty: 0.13440585709080677242e44),
           (tx: 700;     ty: 0.50711602736750225473e304)
         );
const
  name = 'coshm1';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := coshm1(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := ty;
    z := coshm1(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cos;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 1e22;      ty:  0.5232147853951389455),
           (tx: -256.0;    ty: -0.3979075993115770952e-1),
           (tx: 5000.0;    ty:  0.1546684061807471215),
           (tx: -40000.0;  ty:  0.3225874736128876322),
           (tx: 1e6;       ty:  0.9367521275331447870),
           (tx: -1e15;     ty: -0.5131937377869702522),
           (tx: 1e17;      ty: -0.8855573282976306850),
           (tx: 1e20;      ty:  0.7639704044417283004),
           (tx: 0.0;       ty:  0.9031822026143886528),     {2^500}
           (tx: 1.0;       ty:  0.9872460775989134842)      {2^1000}
         );
const
  name = 'cos';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=0.0 then x := ldexpd(1.0,500)
    else if x=1.0 then x := ldexpd(1.0,1000);
    y := ty;
    z := cos(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    if x=0.0 then x := -ldexpd(1.0,500)
    else if x=-1.0 then x := -ldexpd(1.0,1000);
    y := ty;
    z := cos(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_vers;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 1e22;         ty:  0.4767852146048610545),
           (tx: 1e-12;        ty:  0.5e-24),
           (tx: 1e-5;         ty:  0.4999999999958333333e-10),
           (tx: 6.28125;      ty:  0.1872706355174341344e-5),  {1.8727063551743413443418944803655789637335052043489084e-6}
           (tx: 6.283203125;  ty:  0.1587373621400802049e-9),
           (tx: 3.140625;     ty:  1.9999995318233016117),
           (tx: 3.1416015625; ty:  1.9999999999603156595),
           (tx: 0.0;          ty:  0.0),
           (tx: 102950.0;     ty:  0.3820983643652536913e-4),  {0.0000382098364365253691269912458198904181918287693493267864.}
           (tx: 1.0;          ty:  0.4596976941318602826)
         );
const
  name = 'vers/hav';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := vers(x);
    if z<>2*hav(x) then writeln('vers <> 2*hav for x=',x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := ty;
    z := vers(x);
    if z<>2*hav(x) then writeln('vers <> 2*hav for x=',x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_versint;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: -1e-10;  ty:  -0.1666666666666666667e-30),
           (tx: 1e-5;    ty:   0.1666666666658333333e-15),
           (tx: 0.125;   ty:   0.3252666147723100426e-3),
           (tx: -0.5;    ty:  -0.2057446139579699973e-1),
           (tx: 1;       ty:   0.1585290151921034933),
           (tx: Pi_3;    ty:   0.1811721474121590994),
           (tx: 0.0;     ty:   0.0),
           (tx: 2.0;     ty:   1.090702573174318305),
           (tx: -40.0;   ty:  -39.25488683952065121),
           (tx: 10000.0; ty:   10000.30561438888825)
         );
const
  name = 'versint';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := versint(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    if x<>0 then begin
      x := -x;
      y := -y;
      z := versint(x);
      inc(cnt);
      if reldev(z,y) > EPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cos_sc;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 1e22;      ty:  0.5232147853951389455),
           (tx: -256.0;    ty: -0.3979075993115770952e-1),
           (tx: 5000.0;    ty:  0.1546684061807471215),
           (tx: -40000.0;  ty:  0.3225874736128876322),
           (tx: 1e6;       ty:  0.9367521275331447870),
           (tx: -1e15;     ty: -0.5131937377869702522),
           (tx: 1e17;      ty: -0.8855573282976306850),
           (tx: 1e20;      ty:  0.7639704044417283004),
           (tx: 0.0;       ty:  0.9031822026143886528),     {2^500}
           (tx: 1.0;       ty:  0.9872460775989134842)      {2^1000}
         );
const
  name = 'cos (sincos)';
var
  i, cnt, failed: integer;
  x,y,z,t: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=0.0 then x := ldexpd(1.0,500)
    else if x=1.0 then x := ldexpd(1.0,1000);
    y := ty;
    sincos(x,t,z);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    if x=0.0 then x := -ldexpd(1.0,500)
    else if x=-1.0 then x := -ldexpd(1.0,1000);
    y := ty;
    sincos(x,t,z);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cosPi;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;         ty:  1.0),
           (tx: 0.25;        ty:  0.70710678118654752440),
           (tx: 0.5;         ty:  0.0),
           (tx: 1.0;         ty: -1.0),
           (tx: 1.75;        ty: +0.70710678118654752440),
           (tx: 1000.125;    ty:  0.92387953251128675613),
           (tx: 100000.0625; ty:  0.98078528040323044913),
           (tx: 10000000.375;ty:  0.38268343236508977173),
           (tx: -1;          ty:  0.92387953251128675613), {2^40+0.125}
           (tx: -2;          ty:  1.0)                     {2^55}
         );
const
  name = 'cosPi/sincosPi';
var
  i, cnt, failed: integer;
  x,y,z,t: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=-1.0 then x := ldexpd(1.0,40)+0.125
    else if x=-2.0 then x := ldexpd(1.0,55);
    y := ty;
    z := cosPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -x;
    y := ty;
    z := cosPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    sincosPi(x,t,z);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cot;
const
  NT=12;
  TData: array[1..NT] of TPair = (
           (tx: pi_6;    ty: sqrt3 ),
           (tx: pi_4;    ty: 1.0   ),
           (tx: pi_3;    ty: 0.57735026918962576449),
           (tx: 0.5;     ty: 1.8304877217124519193),
           (tx: 1.375;   ty: 0.19833732800536536529),
           (tx: Pi_2;    ty: 0.0),
           (tx: 1e6;     ty: -2.6764843396283451088),
           (tx: -1e15;   ty:  0.5979377907242829412),
           (tx: 1e19;    ty:  0.4044008901234283005),
           (tx: 1e22;    ty: -0.6139571270529418504),
           (tx: 1.0;     ty:  2.1040574227173702189),   {2^500}
           (tx: 2.0;     ty: -6.2012281179184650510));  {2^1000}
const
  name = 'cot';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=1.0 then x := ldexpd(1.0,500)
    else if x=2.0 then x := ldexpd(1.0,1000);
    y := ty;
    z := cot(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -x;
    y := -ty;
    z := cot(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    if abs(x)<1e4 then begin
      {no +-Pi test huge args}
      x := tx+Pi;
      y := ty;
      z := cot(x);
      inc(cnt);
      {y > Pi, Use SEPS here}
      if reldev(z,y) > SEPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
      x := tx-Pi;
      y := ty;
      z := cot(x);
      inc(cnt);
      if reldev(z,y) > EPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_coth;
const
  NT=9;
  TData: array[1..NT] of TPair = (
           (tx: 1e-300; ty: 1e300),
           (tx: 1e-5;    ty: 100000.00000333333333),
           (tx: 0.015625;ty: 64.005208248564253991),
           (tx: 0.125;   ty: 8.0416233283755969285),
           (tx: 0.5;     ty: 2.1639534137386528488),
           (tx: 1.0;     ty: 1.3130352854993313036),
           (tx: 2.0;     ty: 1.0373147207275480959),
           (tx: 30.0;    ty: 1.0),
           (tx: 1e300;   ty: 1.0)
         );
const
  name = 'coth';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := coth(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := coth(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_csc;
const
  NT=7;
  TData: array[1..NT] of TPair = (
           (tx: pi_6;    ty: 2.0 ),
           (tx: pi_4;    ty: sqrt2   ),
           (tx: pi_3;    ty: 1.1547005383792515290),
           (tx: 0.5;     ty: 2.0858296429334881858),
           (tx: 1.0;     ty: 1.1883951057781212163),
           (tx: 1.375;   ty: 1.0194791295952594872),
           (tx: Pi_2;    ty: 1.0)
         );
const
  name = 'csc';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := csc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := tx+TwoPi;
    y := ty;
    z := csc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := Pi-tx;
    y := ty;
    z := csc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := csc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_csch;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx: 1e-300;  ty: 1e300),
           (tx: 1e-5;    ty: 99999.999998333333333),
           (tx: 0.015625;ty: 63.997395907506092977),
           (tx: 0.125;   ty: 7.9792045816280844140),
           (tx: 0.5;     ty: 1.9190347513349437195),
           (tx: 1.0;     ty: 0.85091812823932154513),
           (tx: 2.0;     ty: 0.27572056477178320776),
           (tx: 10.0;    ty: 0.90799859712122162834e-4),
           (tx: 100.0;   ty: 0.74401519520416719259e-43),
           (tx: 1e4;     ty: 0.22709677306294721971e-4342),
           (tx: 1e5;     ty: 0.0)
         );
const
  name = 'csch';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := csch(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := csch(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sec;
const
  NT=6;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0    ),
           (tx: pi_6;    ty: 1.1547005383792515290),
           (tx: pi_4;    ty: sqrt2  ),
           (tx: pi_3;    ty: 2.0  ),
           (tx: 2*pi_3;  ty: -2.0 ),
           (tx: 3*pi_4;  ty: -sqrt2)
         );
const
  name = 'sec';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := sec(x);
    inc(cnt);
    if reldev(z,y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := ty;
    z := sec(x);
    inc(cnt);
    if reldev(z,y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := twopi+tx;
    y := ty;
    z := sec(x);
    inc(cnt);
    if reldev(z,y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := pi-tx;
    y := -ty;
    z := sec(x);
    inc(cnt);
    if reldev(z,y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sech;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx: 1e-10;   ty: 1.0),
           (tx: 1e-5;    ty: 0.99999999995000000000),
           (tx: 0.125;   ty: 0.99223804147512576124),
           (tx: 0.5;     ty: 0.88681888397007390866),
           (tx: 1.0;     ty: 0.64805427366388539957),
           (tx: 2.0;     ty: 0.26580222883407969212),
           (tx: 10.0;    ty: 0.90799859337817244080e-4),
           (tx: 100.0;   ty: 0.74401519520416719259e-43),
           (tx: 1e4;     ty: 0.22709677306294721971e-4342),
           (tx: 1e5;     ty: 0.0)
         );
const
  name = 'sech';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := sech(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := ty;
    z := sech(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sinh;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0),
           (tx: 1e-10;   ty: 1.0e-10),
           (tx: 1e-5;    ty: 0.10000000000166666667e-4),
           (tx: 0.125;   ty: 0.12532577524111545698),
           (tx: 0.5;     ty: 0.52109530549374736162),
           (tx: 1.0;     ty: 1.1752011936438014569),
           (tx: 2.0;     ty: 3.6268604078470187677),
           (tx: 10.0;    ty: 11013.232874703393377),
           (tx: 100.0;   ty: 0.13440585709080677242e44),
           (tx: 700;     ty: 0.50711602736750225473e304)
         );
const
  name = 'sinh';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := sinh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := sinh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sinhmx;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0),
           (tx: 1e-10;   ty: 0.1666666666666666667e-30),
           (tx: 1e-5;    ty: 0.1666666666675000000e-15),
           (tx: 0.125;   ty: 0.3257752411154569821e-3),
           (tx: 0.5;     ty: 0.2109530549374736162e-1),
           (tx: 1.0;     ty: 0.1752011936438014569),
           (tx: 2.0;     ty: 1.626860407847018768),
           (tx: 10.0;    ty: 11003.23287470339338),
           (tx: 100.0;   ty: 0.1344058570908067724e44),
           (tx: 700;     ty: 0.5071160273675022547e304)
         );
const
  name = 'sinhmx';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := sinhmx(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := sinhmx(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sinhc;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 1.0),
           (tx: 1e-8;    ty: 1.0),
           (tx: 1e-5;    ty: 1.0000000000166666667),
           (tx: 0.125;   ty: 1.00260620192892365585646033833),
           (tx: 0.5;     ty: 1.042190610987494723244851252829),
           (tx: 1.0;     ty: 1.1752011936438014569),
           (tx: 10.0;    ty: 1101.3232874703393377),
           (tx: 100.0;   ty: 0.13440585709080677242e42),
           (tx: 700;     ty: 0.724451467667860363896663996594e301),
           (tx: 717;     ty: 0.170841892352061574278887710536e309)
         );
const
  name = 'sinhc';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT-1 do with TData[i] do begin
    x := tx;
    y := ty;
    z := sinhc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := ty;
    z := sinhc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;

  with TData[NT] do begin
    {Special case: sinh(x) overflow, sinhc computed with exp(ln..)}
    x := tx;
    y := ty;
    z := sinhc(x);
    inc(cnt);
    if reldev(z,y) > 10*EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sin;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 1e22;     ty: -0.8522008497671888018),
           (tx: -256.0;   ty:  0.9992080341070626991),
           (tx: 5000.0;   ty: -0.9879664387667768473),
           (tx: -40000.0; ty: -0.9465396567857336878),
           (tx: 1e6;      ty: -0.3499935021712929521),
           (tx: -1e15;    ty: -0.8582727931702358355),
           (tx: 1e19;     ty: -0.9270631660486503852),
           (tx: 1e17;     ty: -0.4645301048353726962),
           (tx: 0;        ty:  0.4292573923424282777),  {2^500}
           (tx: 1;        ty: -0.1592017030862424382)   {2^1000}
         );
const
  name = 'sin';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=0.0 then x := ldexpd(1.0,500)
    else if x=1.0 then x := ldexpd(1.0,1000);
    y := ty;
    z := sin(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    if x=0.0 then x := -ldexpd(1.0,500)
    else if x=-1.0 then x := -ldexpd(1.0,1000);
    y := -ty;
    z := sin(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sinPi;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;         ty:  0.0),
           (tx: 0.25;        ty:  0.70710678118654752440),
           (tx: 0.5;         ty:  1.0),
           (tx: 1.0;         ty:  0.0),
           (tx: 1.75;        ty: -0.70710678118654752440),
           (tx: 1000.125;    ty:  0.38268343236508977173),
           (tx: 100000.0625; ty:  0.19509032201612826785),
           (tx: 10000000.375;ty:  0.92387953251128675613),
           (tx: -1  ;        ty:  0.3826834323650897717),  {2^40+0.125}
           (tx: -2;          ty:  0.0)                     {2^55}
         );
const
  name = 'sinPi/sincosPi';
var
  i, cnt, failed: integer;
  x,y,z,t: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=-1.0 then x := ldexpd(1.0,40)+0.125
    else if x=-2.0 then x := ldexpd(1.0,55);
    y := ty;
    z := sinPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -x;
    y := -ty;
    z := sinPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    sincosPi(x,z,t);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sinc;
const
  NT=14;
  TData: array[1..NT] of TPair = (
           (tx: 1e22;     ty: -0.8522008497671888018e-22),
           (tx: 256.0;    ty: -0.3903156383230713668e-2),
           (tx: 5000.0;   ty: -0.1975932877533553695e-3),
           (tx: 40000.0;  ty:  0.2366349141964334219e-4),
           (tx: 1e6;      ty: -0.3499935021712929521e-6),
           (tx: 0.125;    ty:  0.99739786708182151968),
           (tx: 1e-2;     ty:  0.99998333341666646825),
           (tx: 1e-3;     ty:  0.9999998333333416667),
           (tx: 5e-4;     ty:  0.99999995833333385416),
           (tx: 1e-4;     ty:  0.99999999833333333417),
           (tx: 1e-7;     ty:  0.99999999999999833333),
           (tx: 2e-10;    ty:  1.0000000000000000000),
           (tx: 1e-300;   ty:  1.0000000000000000000),
           (tx: 0;        ty:  1)
         );
const
  name = 'sinc';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := sinc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    z := sinc(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sincPi;
const
  NT=16;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;         ty:  1.0),
           (tx: 0.25;        ty:  0.90031631615710606956),
           (tx: 0.5;         ty:  0.63661977236758134306),
           (tx: 1.0;         ty:  0.0),
           (tx: 1.75;        ty: -0.12861661659387229564),
           (tx: 1000.125;    ty:  0.12179669521365237409e-3),
           (tx: 100000.0625; ty:  0.62099139384550745050e-6),
           (tx: 10000000.375;ty:  0.29407997781320225276e-7),
           (tx: -1;          ty:  0.11078729567138019068e-12),  {2^40+0.125}
           (tx: 0.31e-3;     ty:  0.99999984192184367250),
           (tx: 5e-4;        ty:  0.99999958876653402184),
           (tx: 0.31e-4;     ty:  0.99999999841921836251),
           (tx: 1e-10;       ty:  1.0000000000000000000),
           (tx: 0.7e-10;     ty:  1.0000000000000000000),
           (tx: 0.31e-20;    ty:  1.0000000000000000000),
           (tx: 1e-300;      ty:  1.0000000000000000000)
         );
const
  name = 'sincPi';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=-1.0 then x := ldexpd(1.0,40)+0.125;
    y := ty;
    z := sincPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -x;
    z := sincPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_sin_sc;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 1e22;     ty: -0.8522008497671888018),
           (tx: -256.0;   ty:  0.9992080341070626991),
           (tx: 5000.0;   ty: -0.9879664387667768473),
           (tx: -40000.0; ty: -0.9465396567857336878),
           (tx: 1e6;      ty: -0.3499935021712929521),
           (tx: -1e15;    ty: -0.8582727931702358355),
           (tx: 1e19;     ty: -0.9270631660486503852),
           (tx: 1e17;     ty: -0.4645301048353726962),
           (tx: 0;        ty:  0.4292573923424282777),  {2^500}
           (tx: 1;        ty: -0.1592017030862424382)   {2^1000}
         );
const
  name = 'sin (sincos)';
var
  i, cnt, failed: integer;
  x,y,z,t: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=0.0 then x := ldexpd(1.0,500)
    else if x=1.0 then x := ldexpd(1.0,1000);
    y := ty;
    sincos(x,z,t);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    if x=0.0 then x := -ldexpd(1.0,500)
    else if x=-1.0 then x := -ldexpd(1.0,1000);
    y := -ty;
    sincos(x,z,t);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_tan;
const
  NT=12;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;     ty: 0.0    ),
           (tx: pi_6;    ty: sqrt3/3),
           (tx: pi_4;    ty: 1.0    ),
           (tx: pi_3;    ty: sqrt3  ),
           (tx: 2*pi_3;  ty: -sqrt3 ),
           (tx: 3*pi_4;  ty: -1.0   ),
           (tx: 1e6;     ty: -0.3736244539875990292),
           (tx: -1e15;   ty:  1.6724147821275830430),
           (tx: 1e19;    ty:  2.4727937658465273603),
           (tx: 1e22;    ty: -1.6287782256068988493),
           (tx: 1.0;     ty:  0.4752721998948628874),   {2^500}
           (tx: 2.0;     ty: -0.1612583799506580567));  {2^1000}
const
  name = 'tan';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=1.0 then x := ldexpd(1.0,500)
    else if x=2.0 then x := ldexpd(1.0,1000);
    y := ty;
    z := tan(x);
    inc(cnt);
    if reldev(z,y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -x;
    z := tan(x);
    inc(cnt);
    if reldev(z,-y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',-y:19, ' ',z:19,' ',reldev(z,-y):18);
    end;
    if abs(x)<1e4 then begin
      {no +-Pi test huge args}
      x := tx+Pi;
      z := tan(x);
      inc(cnt);
      if reldev(z,y) > SEPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
      x := tx-Pi;
      z := tan(x);
      inc(cnt);
      if reldev(z,y) > SEPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_tanh;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx: 1e-300;  ty: 1e-300),
           (tx: 1e-5;    ty: 0.99999999996666666667e-5),
           (tx: 0.015625;ty: 0.15623728558408865202e-1),
           (tx: 0.125;   ty: 0.12435300177159620805),
           (tx: 0.5;     ty: 0.46211715726000975850),
           (tx: 1.0;     ty: 0.76159415595576488812),
           (tx: 2.0;     ty: 0.9640275800758168840),
           (tx: 10.0;    ty: 0.9999999958776927636),
           (tx: 18.0;    ty: 0.9999999999999995361),
           (tx: 20.0;    ty: 1.0),
           (tx: 1e300;   ty: 1.0)
         );
const
  name = 'tanh';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := tanh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    y := -ty;
    z := tanh(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_tanPi;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;         ty:  0.0),
           (tx: 0.25;        ty:  1.0),
           (tx: 0.28125;     ty:  1.2185035255879763445),
           (tx: 1.0;         ty:  0.0),
           (tx: 1.75;        ty:  -1.0),
           (tx: 1000.125;    ty:  0.4142135623730950488),
           (tx: 100000.0625; ty:  0.1989123673796580069),
           (tx: 10000000.375;ty:  2.4142135623730950488),
           (tx: -1  ;        ty:  0.4142135623730950488),  {2^34+0.125}
           (tx: -2;          ty:  0.0)                     {2^55}
         );
const
  name = 'tanPi';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    if x=-1.0 then x := ldexpd(1.0,34)+0.125
    else if x=-2.0 then x := ldexpd(1.0,55);
    y := ty;
    if x=0.5 then y := PosInf_d;
    z := tanPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -x;
    y := -ty;
    z := tanPi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  inc(cnt);
  z := tanPi(0.5);
  if not isinfd(z) then begin
    inc(failed);
    writeln(x:19,' ',PosInf_d:19, ' ',z:19);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_covers;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 1e22;     ty: 1.8522008497671888018),
           (tx: -256.0;   ty: 0.7919658929373009324e-3),
           (tx: 5000.0;   ty: 1.9879664387667768472),
           (tx: -10000.0; ty: 0.6943856111117478586),
           (tx: 1e6;      ty: 1.3499935021712929521),
           (tx: 51/32;    ty: 0.2634239906243735828e-3), {0.0002634239906243735827877883107694945728082794656929790107}
           (tx: -51/32;   ty: 1.9997365760093756264),
           (tx: 1/1024;   ty: 0.9990234376552204217),
           (tx: -151/32;  ty: 0.2023121705578823636e-4), {0.0000202312170557882363638059164907190091530431528101282650}
           (tx: 0;        ty: 1)
         );
const
  name = 'covers';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := covers(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_archav;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0;            ty: 0),
           (tx: 1e-10;        ty: 0.2000000000033333333e-4),
           (tx: 1/1048576;    ty: 0.1953125310440991432e-2),
           (tx: 1/1024;       ty: 0.6251017699899030937e-1),
           (tx: 0.25;         ty: 1.0471975511965977462),   {1.0471975511965977461542144610931676280657231331250352736583}
           (tx: 0.5;          ty: 1.5707963267948966192),
           (tx: 0.75;         ty: 2.0943951023931954923),
           (tx: 0.9990234375; ty: 3.0790824765908029291),
           (tx: 1-1/1048576;  ty: 3.1396395282793522470),   {3.1396395282793522470306452950615808530523109415655700564779}
           (tx: 1;            ty: 3.1415926535897932385)
         );
const
  name = 'archav';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := archav(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_gd;
const
  NT=10;
  TData: array[1..NT] of TPair = (
           (tx: 0;            ty:  0),
           (tx: -1e-8;        ty: -0.9999999999999999833e-8),
           (tx: 1/1048576;    ty:  0.9536743164061054397e-6),
           (tx: -1/1024;      ty: -0.9765623447796079048e-3),{-0.000976562344779607904844149783145849353952794402901433691}
           (tx: 0.25;         ty:  0.2474357989824314828),
           (tx: -0.5;         ty: -0.4803810791337294486),   {-0.480381079133729448604886289963232005609691792423076398821}
           (tx: 1;            ty:  0.8657694832396586243),   {0.8657694832396586242896018461918444413796791992487600996118}
           (tx: -10;          ty: -1.5707055269354340337),
           (tx: 40;           ty:  1.5707963267948966107),
           (tx: -50;          ty: -1.5707963267948966192)
          );

const
  name = 'gd';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := gd(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_arcgd;
const
  NT=11;
  TData: array[1..NT] of TPair = (
           (tx: 0;            ty:  0),
           (tx: -1e-8;        ty: -0.1000000000000000017e-7),
           (tx: 1/1048576;    ty:  0.9536743164063945603e-6), {9.5367431640639456028966476679373871125762481962914326e-7}
           (tx: -1/1024;      ty: -0.9765626552204661100e-3),
           (tx: 0.25;         ty:  0.2526456103578675383),    {0.2526456103578675382941100556784265475844701350423196409701}
           (tx: -0.5;         ty: -0.5222381032784403302),
           (tx: 1;            ty:  1.2261911708835170708),
           (tx: -1.5;         ty: -3.3406775427983110033),
           (tx: 25/16;        ty:  5.4850838617513340825),
           (tx: -201/128;     ty: -8.326930738131876848),
           (tx: 102943/65536; ty: 12.12871039236958305)       {12.128710392369583049995127447130112365303060526568156991862}
          );

const
  name = 'arcgd';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := arcgd(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_power;
var
  a,b,x,y: double;
  cnt, failed: integer;

  procedure itest;
  begin
    a := power(x,y);
    inc(cnt);
    if reldev(a,b) > SEPS then begin
      inc(failed);
      writeln('  power: ',x:5:2,' ',y:5:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
    end;
    if x>0.0 then begin
      a := logbase(x,b);
      inc(cnt);
      if reldev(a,y) > SEPS then begin
        inc(failed);
        writeln('logbase: ',x:5:2,' ',y:5:2, ' ',a:19,' ',b:19,' ',reldev(a,y):16);
      end;
      a := logN(x,b);
      inc(cnt);
      if reldev(a,y) > SEPS then begin
        inc(failed);
        writeln('   logN: ',x:5:2,' ',y:5:2, ' ',a:19,' ',b:19,' ',reldev(a,y):16);
      end;
    end;
  end;
begin
  writeln('Function power/logbase/logN');
  cnt :=0;
  failed := 0;

  x := +2.3; y := +1.2;  b := 2.7168984324991494779;      itest;
  x := +2.3; y := -1.2;  b := 0.36806675878573280046;     itest;

  x := +2.3; y := +0.2;  b := 1.1812601880431084686;      itest;
  x := +2.3; y := -0.2;  b := 0.84655354520718544105;     itest;

  x := +0.2; y := +2.3;  b := 0.24681354508800385463e-1;  itest;
  x := +0.2; y := -2.3;  b := 40.516414917319060879;      itest;

  x := +2.3; y := +3;    b := 12.167;                     itest;
  x := +2.3; y := -3;    b := 0.82189529053998520587e-1;  itest;
  x := -2.3; y := +3;    b := -12.167;                    itest;
  x := -2.3; y := -3;    b := -0.82189529053998520587e-1; itest;

  x := +2.3; y := +4;    b := 27.9841;                    itest;
  x := +2.3; y := -4;    b := 0.35734577849564574168e-1;  itest;
  x := -2.3; y := +4;    b := 27.9841;                    itest;
  x := -2.3; y := -4;    b := 0.35734577849564574168e-1;  itest;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_powm1;
var
  a,b,x,y,sm: double;
  cnt, failed: integer;

  procedure itest;
  begin
    a := powm1(x,y);
    inc(cnt);
    if reldev(a,b) > SEPS then begin
      inc(failed);
      writeln('  powm1: ',x:5:2,' ',y:5:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
    end;
  end;

begin
  writeln('Function powm1');
  cnt :=0;
  failed := 0;
  sm  := ldexpd(1,-30);

  x := -1.5;
  y := -2.0;  b := -0.55555555555555555555;  itest;
  y := -1.0;  b := -1.6666666666666666667;   itest;
  y :=  0.0;  b :=  0.0;   itest;
  y :=  1.0;  b := -2.5;   itest;
  y :=  2.0;  b :=  1.25;  itest;

  x := -1.0;
  y := -2.0;  b :=  0.0;   itest;
  y := -1.0;  b := -2.0;   itest;
  y :=  0.0;  b :=  0.0;   itest;
  y :=  1.0;  b := -2.0;   itest;
  y :=  2.0;  b :=  0.0;   itest;

  x := -sm;
  y := -2.0;  b :=  1152921504606846975.0;   itest;
  y := -1.0;  b := -1073741825.0000000000;   itest;
  y :=  0.0;  b :=  0.0;                     itest;
  y :=  1.0;  b := -1.0000000009313225746;   itest;
  y :=  2.0;  b := -0.99999999999999999913;  itest;

  x :=  0.0;
  y :=  0.0;  b :=  0.0;   itest;
  y :=   sm;  b := -1.0;   itest;
  y :=  1.0;  b := -1.0;   itest;
  y :=  1.5;  b := -1.0;   itest;
  y :=  2.0;  b := -1.0;   itest;

  x := +sm;
  y := -2.0;  b :=  1152921504606846975.0;      itest;
  y := -1.5;  b :=  35184372088831.000000;      itest;
  y := -1.0;  b :=  1073741823.0000000000;      itest;
  y := -sm;   b :=  0.19366308691123400479e-7;  itest;
  y :=  0.0;  b :=  0.0;                        itest;
  y :=  sm;   b := -0.19366308316069495422e-7;  itest;
  y :=  1.0;  b := -0.99999999906867742538;     itest;
  y :=  1.5;  b := -0.99999999999997157829;     itest;
  y :=  2.0;  b := -0.99999999999999999913;     itest;

  x :=  1.0;
  y := -2.0;  b :=  0.0;   itest;
  y := -1.5;  b :=  0.0;   itest;
  y := -1.0;  b :=  0.0;   itest;
  y :=  -sm;  b :=  0.0;   itest;
  y :=  0.0;  b :=  0.0;   itest;
  y :=   sm;  b :=  0.0;   itest;
  y :=  1.0;  b :=  0.0;   itest;
  y :=  1.5;  b :=  0.0;   itest;
  y :=  2.0;  b :=  0.0;   itest;

  x := 1.0+sm;
  y := -2.0;  b := -0.18626451466288718205e-8;  itest;
  y := -1.5;  b := -0.13969838602969145165e-8;  itest;
  y := -1.0;  b := -0.93132257374811677845e-9;  itest;
  y :=  -sm;  b := -0.86736173758450676361e-18; itest;
  y :=  0.0;  b :=  0.0;                        itest;
  y :=   sm;  b :=  0.86736173758450676436e-18; itest;
  y :=  1.0;  b :=  0.93132257461547851563e-9;  itest;
  y :=    x;  b :=  0.93132257548284025402e-9;  itest;
  y :=  2.0;  b :=  0.18626451500983187692e-8;  itest;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_pow1pm1;
var
  a,b,x,y,sm: double;
  cnt, failed: integer;

  procedure itest;
  begin
    a := pow1pm1(x,y);
    inc(cnt);
    if reldev(a,b) > SEPS then begin
      inc(failed);
      writeln(' pow1pm1: ',x:5:2,' ',y:5:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
    end;
  end;

begin
  writeln('Function pow1pm1');
  cnt :=0;
  failed := 0;
  sm  := ldexpd(1,-30);

  x := -1.5;
  y := -2.0;  b :=  3.0000000000000000000;  itest;
  y := -1.0;  b := -3.0000000000000000000;  itest;
  y :=  0.0;  b :=  0.0;   itest;
  y :=  1.0;  b := -1.5;   itest;
  y :=  2.0;  b := -0.75;  itest;

  x := -1.0;
  y :=  0.0;  b :=  0.0;   itest;
  y :=  1.0;  b := -1.0;   itest;
  y :=  2.0;  b := -1.0;   itest;

  x := -1.0 + sm;
  y := -2.0;  b :=  1152921504606846975.00;    itest;
  y := -1.5;  b :=  35184372088831.0000000;    itest;
  y := -1.0;  b :=  1073741823.00000000000;    itest;
  y := -sm;   b :=  0.19366308691123400479e-7; itest;
  y :=  0.0;  b :=  0.0;                       itest;
  y :=  sm;   b := -0.19366308316069495422e-7; itest;
  y :=  1.0;  b := -0.99999999906867742538;    itest;
  y :=  1.5;  b := -0.99999999999997157829;    itest;
  y :=  2.0;  b := -0.99999999999999999913;    itest;

  x := -sm;
  y := -2.0;  b :=  0.18626451518330422484e-8;  itest;
  y := -1.5;  b :=  0.13969838635495210339e-8;  itest;
  y := -1.0;  b :=  0.93132257548284025442e-9;  itest;
  y := -sm;   b :=  0.86736173839230033131e-18; itest;
  y :=  0.0;  b :=  0.0;                        itest;
  y :=  sm;   b := -0.86736173839230033055e-18; itest;
  y :=  1.0;  b := -0.93132257461547851563e-9;  itest;
  y :=  1.5;  b := -0.13969838615979571216e-8;  itest;
  y :=  2.0;  b := -0.18626451483635952933e-8;  itest;

  x :=  0.0;
  y := -2.0;  b :=  0.0;   itest;
  y := -1.5;  b :=  0.0;   itest;
  y := -1.0;  b :=  0.0;   itest;
  y :=  -sm;  b :=  0.0;   itest;
  y :=  0.0;  b :=  0.0;   itest;
  y :=   sm;  b :=  0.0;   itest;
  y :=  1.0;  b :=  0.0;   itest;
  y :=  1.5;  b :=  0.0;   itest;
  y :=  2.0;  b :=  0.0;   itest;

  x := +sm;
  y := -2.0;  b := -0.18626451466288718205e-8;  itest;
  y := -1.5;  b := -0.13969838602969145165e-8;  itest;
  y := -1.0;  b := -0.93132257374811677845e-9;  itest;
  y := -sm;   b := -0.86736173758450676361e-18; itest;
  y :=  0.0;  b :=  0.0;                        itest;
  y :=  sm;   b :=  0.86736173758450676436e-18; itest;
  y :=  1.0;  b :=  0.93132257461547851563e-9;  itest;
  y :=  1.5;  b :=  0.13969838622484784251e-8;  itest;
  y :=  2.0;  b :=  0.18626451500983187692e-8;  itest;

  x := 1.0-sm;
  y := -2.0;  b := -0.74999999976716935618;     itest;
  y := -1.5;  b := -0.64644660915977204716;     itest;
  y := -1.0;  b := -0.49999999976716935624;     itest;
  y :=  -sm;  b := -0.64554361614450407531e-9;  itest;
  y :=  0.0;  b :=  0.0;                        itest;
  y :=   sm;  b :=  0.64554361656123063592e-9;  itest;
  y :=  1.0;  b :=  0.99999999906867742538;     itest;
  y :=  1.5;  b :=  1.82842712277055657389;     itest;
  y :=  2.0;  b :=  2.99999999627470970241;     itest;

  x := 1e-9;
  y := 10000;
  b := 0.10000049995166617086e-4;
  itest;

  x := 1e-15;
  y := 10000;
  b := 0.1000000000004999500e-10;
  itest;

  x := 1e-9;
  y := 1e10;
  b := 22025.465684674387892;
  itest;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_intpower;
var
  a,b,x: double;
  cnt, failed: integer;
  n: longint;
const
  vs: THexDblW = ($0000,$0000,$1000,$0000);  { 8.69169475979376E-0311}

  procedure itest;
  begin
    a := intpower(x,n);
    inc(cnt);
    if reldev(a,b) > SEPS then begin
      inc(failed);
      writeln('  !!! ',x:6:2,' ',n:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
    end;
  end;

begin
  writeln('Function intpower');
  cnt :=0;
  failed := 0;

  x := 10.0;
  n :=  5;    b := 1e5;   itest;
  n := -5;    b := 1e-5;  itest;
  x := -10.0;
  n :=  6;    b := 1e6;   itest;
  n :=  -6;   b := 1e-6;  itest;

  x := -1.0;
  n :=  17;   b :=  -1;   itest;
  n :=  -17;  b :=  -1;   itest;
  n :=  18;   b :=  1;    itest;
  n :=  -18;  b :=  1;    itest;

  x := 12.5;
  n := 0;     b :=  1;    itest;

  x := 2;
  n := -1030; b := double(vs); itest;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_hypot;
var
  cnt, failed: integer;

  procedure itest(tx,ty,b: double);
  var
    a,x,y: double;
  begin
    x := tx;
    y := ty;
    a := hypot(x,y);
    inc(cnt);
    if reldev(a,b) > EPSD then begin
      inc(failed);
      writeln(x:8:2,' ',y:8:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
    end;
    if (x=0.0) and (y=0.0) then exit;
    if x<>0.0 then begin
      x := -tx;
      y := ty;
      a := hypot(x,y);
      inc(cnt);
      if reldev(a,b) > EPSD then begin
        inc(failed);
        writeln(x:8:2,' ',y:8:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
      end;
    end;
    if y<>0.0 then begin
      x := tx;
      y := -ty;
      a := hypot(x,y);
      inc(cnt);
      if reldev(a,b) > EPSD then begin
        inc(failed);
        writeln(x:8:2,' ',y:8:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
      end;
    end;
    if (y<>0.0) and (x<>0.0) then begin
      x := -tx;
      y := -ty;
      a := hypot(x,y);
      inc(cnt);
      if reldev(a,b) > EPSD then begin
        inc(failed);
        writeln(x:8:2,' ',y:8:2, ' ',a:19,' ',b:19,' ',reldev(a,b):16);
      end;
    end;
  end;

begin
  writeln('Function hypot');
  cnt :=0;
  failed := 0;

  itest(0.0,0.0,0.0);
  itest(1.0,0.0,1.0);
  itest(0.0,1.0,1.0);
  itest(1.0,1.0,sqrt(2));
  itest(1.0,2.0,sqrt(5));
  itest(2.0,1.0,sqrt(5));
  itest(1e20,1.0,1e20);
  itest(1,1e20,1e20);
  itest(1e-20,1.0,1.0);
  itest(1,1e-20,1.0);
  itest(0,MinDouble,Mindouble);
  itest(Mindouble,0,Mindouble);
  itest(Mindouble,Mindouble,Mindouble*sqrt(2));

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_hypot3;
var
  cnt, failed: integer;
  t: double;

  procedure itest(n: integer; x,y,z,b: double);
  var
    a: double;
  begin
    a := hypot3(x,y,z);
    inc(cnt);
    if reldev(a,b) > EPSD then begin
      inc(failed);
      writeln(' Test ',n, ' failed; rel. dev. = ', reldev(a,b):16);
    end;
  end;

begin
  writeln('Function hypot3');
  cnt :=0;
  failed := 0;

  itest( 1, 0.0,0.0,0.0, 0.0);
  itest( 2, 1.0,0.0,0.0, 1.0);
  itest( 3, 0.0,1.0,0.0, 1.0);
  itest( 4, 0.0,0.0,1.0, 1.0);
  itest( 5, 1.0,1.0,1.0, sqrt(3));
  itest( 6, Mindouble,0.0,0.0, Mindouble);
  itest( 7, 0.0,Mindouble,0.0, Mindouble);
  itest( 8, 0.0,0.0,Mindouble, Mindouble);
  itest( 9, Maxdouble,0.0,0.0, Maxdouble);
  itest(10, 0.0,Maxdouble,0.0, Maxdouble);
  itest(11, 0.0,0.0,Maxdouble, Maxdouble);
  itest(12, Mindouble,Mindouble,Mindouble, sqrt(3)*Mindouble);

  t := succd(0);
  itest(14, t,0,0,t);
  itest(15, t,t,0,sqrt(2)*t);
  itest(16, t,t,t,sqrt(3)*t);

  t := 5e-10;
  itest(18, 1,t,0,1.0+1.25e-19);
  itest(19, 1,t,t,1.0+2.5e-19);

  t := 0.5*Maxdouble;
  itest(20, t,0,0,t);
  itest(21, 0,t,t,sqrt(2)*t);
  itest(22, t,t,t,sqrt(3)*t);

  t := 0.25*Maxdouble;
  itest(22, t,2*t,3*t,sqrt(14)*t);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_dot_norm_sumsqr;
var
  xd,yd: array[0..20] of double;
  i,j,f: integer;
  cnt, failed: integer;
  dd,e,ee,s,x,nd,rd,sd,q,sa: double;
begin
  writeln('Test dot, norm, rms, sum, sumsqr');
  cnt :=0;
  failed := 0;
  i := 1;
  randseed := 0;

  inc(cnt,6);
  while (i<=1000) and (failed=0) do begin
    s  := 0.0;
    sa := 0.0;
    for j:=0 to 20 do begin
      xd[j] := random - 0.5;
      yd[j] := random - 0.5;
      x  := xd[j]*yd[j];
      s  := s  + x;
      sa := sa + abs(x)
    end;
    dd := dot2(xd,yd,21);
    if abs(s-dd) > 2*EPSD*sa  then begin
      inc(failed);
      writeln('random dot2 failed.  i=',i,':', dd:24,' ',s:24);
      writeln(reldev(dd,s)/EPSD);
    end;
    dd := dot2(xd,xd,21);
    sd := sumsqr(xd,21);
    rd := rms(xd,21);
    x  := sqrt(dd/21);
    if abs(sd-dd)>4*EPSD*abs(dd) then begin
      inc(failed);
      writeln('random sumsqr failed.  i=',i,':', sd:24,' ',dd:24);
    end;
    if abs(x-rd)>4*EPSD*abs(x) then begin
      inc(failed);
      writeln('random rms failed.  i=',i,':', rd:24,' ',x:24);
    end;

    dd := sqrt(dd);
    nd := norm2(xd,21);
    if abs(nd-dd)>4*EPSD*abs(dd) then begin
      inc(failed);
      writeln('random norm2 failed.  i=',i,':', nd:24,' ',dd:24);
    end;

    inc(i);
  end;

  i := 1;
  f := 0;
  e := 0.5;
  inc(cnt,2);

  while (e>EPSD) and (f=0) do begin
    {dot = 1*1 + e*e - 1*1 = e*e}
    ee := e*e;
    xd[0] := 1;
    xd[1] := e;
    xd[2] := 1;
    yd[0] := 1;
    yd[1] := e;
    yd[2] := -1;
    dd := dot2(xd,yd,3);
    if abs(dd-ee)>EPSD*ee then begin
      inc(failed);
      f := 1;
      writeln('eps dot2 failed.  i=',i,':', dd:24,' ',e:24);
    end;
    e := e/1024;
    inc(i);
  end;

  f := 0;
  q := 0.1;
  inc(cnt,4);
  while (abs(q)<=0.99) and (f=0) do begin
    x := 1.0;
    for i:=0 to 20 do begin
      xd[i] := x;
      x := x*q;
    end;
    sd := sum2(xd,21);
    x := (1.0-intpower(q,21))/(1.0-q);
    if abs(sd-x)>4*EPSD*abs(x) then begin
      inc(failed);
      f := 1;
      writeln('sum2 failed. q=',q:15:12,' ', sd:20,' ',x:20);
    end;
    sd := sumsqr(xd,21);
    x := (1.0-intpower(q,2*21))/((1.0-q)*(1.0+q));
    if abs(sd-x)>8*EPSD*abs(x) then begin
      inc(failed);
      f := 1;
      writeln('sumsqr failed. q=',q:15:12,' ', sd:20,' ',x:20);
    end;
    q := -q*1.125;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cbrt;
const
  NT=14;
  TData: array[1..NT] of TPair = (
           (tx: 0.0;       ty: 0.4291754585626620336847e-107),  {2^(-1070)}
           (tx: 0.0;       ty: 0.0),
           (tx: 1e-307;    ty: 0.464158883361277889241e-102),
           (tx: 1e-100;    ty: 0.46415888336127788924e-33),
           (tx: 1e-16;     ty: 0.46415888336127788924e-5),
           (tx: 0.0078125; ty: 0.19842513149602493434),
           (tx: 0.125;     ty: 0.5),
           (tx: 0.5;       ty: 0.79370052598409973738),
           (tx: 1.0;       ty: 1.0),
           (tx: 2.5;       ty: 1.3572088082974532858),
           (tx: 100.0;     ty: 4.6415888336127788924),
           (tx: 1e22;      ty: 21544346.900318837218),
           (tx: 1e100;     ty: 0.21544346900318837218e34),
           (tx: 1e308;     ty: 0.46415888336127788924e103)
        );
const
  name = 'cbrt';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    if i=1 then x := ldexpd(1,-1070) else x := tx;
    y := ty;
    z := cbrt(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    if x<>0.0 then begin
      x := -x;
      y := -ty;
      z := cbrt(x);
      inc(cnt);
      if reldev(z,y) > EPS then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_nroot;
const
  name = 'nroot';
var
  n, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;

  {nroot(x,) ~ Maple surd(x,n)}
  n := 6;
  x := 1;
  y := nroot(x,n);
  z := 1;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := 5;
  x := -1;
  y := nroot(x,n);
  z := -1;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := 7;
  x := -1e22;
  y := nroot(x,n);
  z := -1389.4954943731376371;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := -5;
  x := -100;
  y := nroot(x,n);
  z := -0.39810717055349725077;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := -4;
  x := 2;
  y := nroot(x,n);
  z := 0.84089641525371454303;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := 6;
  x := 1e308;
  y := nroot(x,n);
  z := 0.2154434690031883721759e52;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := -9;
  x := -4e-308;
  y := nroot(x,n);
  z := -0.14299691483087287440366e35;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := 500;
  x := 10;
  y := nroot(x,n);
  z := 1.0046157902783951424;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := -32103;
  x := 0.125;
  y := nroot(x,n);
  z := 1.0000647761545670534;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := -1;
  x := 3;
  y := nroot(x,n);
  z := 0.33333333333333333333;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := -2;
  x := 5;
  y := nroot(x,n);
  z := 0.44721359549995793928;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := 1;
  x := -17;
  y := nroot(x,n);
  z := -17;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := 2;
  x := 2;
  y := nroot(x,n);
  z := 1.4142135623730950488;
  inc(cnt);
  if reldev(z,y) > EPS then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18,' ',reldev(z,y):18);
  end;

  n := 2;
  x := PosInf_d;
  y := nroot(x,n);
  z := PosInf_d;
  inc(cnt);
  if y<>z then begin
    inc(failed);
    writeln(n:4,' ',x:16,' ',y:18, ' ',z:18);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_rem_2pi;
const
  NS=10;
  SData: array[1..NS] of TPair = (
           (tx: 1e22;     ty: -0.8522008497671888018),
           (tx: -256.0;   ty:  0.9992080341070626991),
           (tx: 5000.0;   ty: -0.9879664387667768473),
           (tx: -40000.0; ty: -0.9465396567857336878),
           (tx: 1e6;      ty: -0.3499935021712929521),
           (tx: -1e15;    ty: -0.8582727931702358355),
           (tx: 1e19;     ty: -0.9270631660486503852),
           (tx: 1e12;     ty: -0.6112387023768894982),
           (tx: 0;        ty:  0.4292573923424282777),  {2^500}
           (tx: 1;        ty: -0.1592017030862424382)   {2^1000}
         );
const
  NT=8;
  TData: array[1..NT] of TPair = (
           (tx: 1e20;     ty:  5.581833149464241096),
           (tx: -512.0;   ty:  3.221195188726091108),
           (tx: 9999.0;   ty:  2.452176277277915212),
           (tx: -20000.0; ty:  5.662018059803342530),
           (tx: 1e8;      ty:  1.942695134504014460),
           (tx: -1e12;    ty:  0.657624759136786467),
           (tx: -1e19;    ty:  1.955090807925410481),
           (tx: 1e22;     ty:  5.263007914620499494)
         );
const
  name = 'rem_2pi';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NS do with SData[i] do begin
    x := tx;
    if x=0.0 then x := ldexpd(1.0,500)
    else if x=1.0 then x := ldexpd(1.0,1000);
    y := ty;
    z := sin(rem_2pi(x));
    inc(cnt);
    if reldev(z,y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    if x=0.0 then x := -ldexpd(1.0,500)
    else if x=-1.0 then x := -ldexpd(1.0,1000);
    y := -ty;
    z := sin(rem_2pi(x));
    inc(cnt);
    if reldev(z,y) > SEPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := rem_2pi(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_rem_2pi_sym;
const
  NS=10;
  SData: array[1..NS] of TPair = (
           (tx: 1e22;     ty: -0.8522008497671888018),
           (tx: -256.0;   ty:  0.9992080341070626991),
           (tx: 5000.0;   ty: -0.9879664387667768473),
           (tx: -40000.0; ty: -0.9465396567857336878),
           (tx: 1e6;      ty: -0.3499935021712929521),
           (tx: -1e15;    ty: -0.8582727931702358355),
           (tx: 1e19;     ty: -0.9270631660486503852),
           (tx: 1e12;     ty: -0.6112387023768894982),
           (tx: 0;        ty:  0.4292573923424282777),  {2^500}
           (tx: 1;        ty: -0.1592017030862424382)   {2^1000}
         );
const
  NT=8;
  TData: array[1..NT] of TPair = (
           (tx: 1e20;     ty:  -0.701352157715345382),
           (tx: -512.0;   ty:  -3.061990118453495369),
           (tx: 9999.0;   ty:   2.452176277277915212),
           (tx: -20000.0; ty:  -0.621167247376243947),
           (tx: 1e8;      ty:   1.942695134504014460),
           (tx: -1e12;    ty:   0.657624759136786467),
           (tx: -1e19;    ty:   1.955090807925410481),
           (tx: 1e21;     ty:  -0.7303362699738673424)
         );
const
  name = 'rem_2pi_sym';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NS do with SData[i] do begin
    x := tx;
    if x=0.0 then x := ldexpd(1.0,500)
    else if x=1.0 then x := ldexpd(1.0,1000);
    y := ty;
    z := sin(rem_2pi_sym(x));
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
    x := -tx;
    if x=0.0 then x := -ldexpd(1.0,500)
    else if x=-1.0 then x := -ldexpd(1.0,1000);
    y := -ty;
    z := sin(rem_2pi(x));
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := rem_2pi_sym(x);
    inc(cnt);
    if reldev(z,y) > EPS then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ',reldev(z,y):18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  if failed>0 then failed := 1;
  inc(total_cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_randg;
{$ifdef BIT16}
const
  NMAX=160000;
{$else}
const
  NMAX=640000;
{$endif}
var
  mx,sx,x,d: double;
  n: longint;
  fail: boolean;
begin
  writeln('Test randg');
  mx := 0;
  sx := 0;
  n := 1;
  while n<NMAX do begin
    x := randg(TwoPi,Pi_4);
    d := x-mx;
    mx := mx+d/n;
    sx := sx + sqr(d)*(n-1)/n;
    inc(n);
  end;
  fail := false;
  x := 0.75*sqrt(n);
  d := (mx-TwoPi)*x;
  if abs(d)>1 then begin
    writeln(' Mean error: ', d);
    fail := true;
  end;
  d := (sqrt(sx/n)-Pi_4)*x;
  if abs(d)>1 then begin
    writeln(' StdDev error: ', d);
    fail := true;
  end;
  inc(total_cnt);
  if fail then begin
    writeln('*** failed');
    inc(total_failed);
  end
  else  writeln(' OK.')
end;


{---------------------------------------------------------------------------}
procedure test_damath_main;
  {-DAMath regression test via single procedure call}
label
  done;
begin
  total_cnt    := 0;
  total_failed := 0;
  writeln('T_DAMATH - main test program for DAMATH unit V',DAMath_Version,'   (c) 2009-2018  W.Ehrhardt');

  {Values for rmNearest}
  EPS    := 3E-15;
  SEPS   := 3E-14;
  EPSD   := 2*eps_d;

{
  test_misc2d;
  goto done;
}

  writeln('Test of basic functions:');
  {First do frexp/ldexp/ilogb, does halt if error}
  test_frexp_ldexp_ilogb;
  Test_IsInfNan;
  Test_Consts;
  test_ceil;
  test_floor;
  test_rint;
  test_succ_pred;
  test_ulp;
  test_dminmax;
  test_sminmax;
  test_hex2float;


  writeln;
  writeln('Tests polynomial, vector, statistic, and other functions:');
  test_polevl;
  test_polevale;
  test_polevalche;
  test_dot_norm_sumsqr;
  test_meansdev;
  test_stat1;
  test_stat_moment;
  test_randg;
  test_angle2;
  test_geometry2;

  writeln;
  writeln('Tests double-double functions:');
  test_ddivrem;
  test_fma_d;
  test_misc2d;


  writeln;
  writeln('Tests for accurate new functions using representative values and symmetries:');
  test_arccos;
  test_arccos1m;
  test_arccosh1p;
  test_arccosh;
  test_arccot;
  test_arccotc;
  test_arccoth;
  test_arccsc;
  test_arccsch;
  test_arcgd;
  test_archav;
  test_arcsec;
  test_arcsech;
  test_arcsin;
  test_arcsinh;
  test_arctanh;
  test_cbrt;
  test_cos;
  test_cosh;
  test_coshm1;
  test_cosPi;
  test_cos_sc;
  test_cot;
  test_coth;
  test_covers;
  test_csc;
  test_csch;
  test_exp10;
  test_exp10m1;
  test_exp2;
  test_exp2m1;
  test_exp3;
  test_exp;
  test_expm1;
  test_exprel;
  test_expx2;
  test_expmx2h;
  test_gd;
  test_hypot;
  test_hypot3;
  test_intpower;
  test_ln;
  test_lncosh;
  test_lnsinh;
  test_ln1mexp;
  test_ln1pexp;
  test_ln1p;
  test_ln1pmx;
  test_lnxp1;
  test_log2p1;
  test_log10p1;
  test_log10;
  test_log2;
  test_logistic;
  test_logit;
  test_log_add_sub_exp;
  test_nroot;
  test_pow1pm1;
  test_power;
  test_powm1;
  test_compound;
  test_comprel;
  test_pow1p;
  test_powpi2k;
  test_rem_2pi;
  test_rem_2pi_sym;
  test_sec;
  test_sech;
  test_sin;
  test_sinc;
  test_sincPi;
  test_sinh;
  test_sinhc;
  test_sinhmx;
  test_sinPi;
  test_sin_sc;
  test_sqrt1pmx;
  test_tan;
  test_tanh;
  test_tanPi;
  test_vers;
  test_versint;
  test_trig_degree;

done:

  writeln;

  if total_failed>0 then begin
    writeln('***** Total number of failed tests: ',total_failed,' of ',total_cnt);
  end
  else writeln('Passed. All ',total_cnt, ' test were OK.');

end;

end.
