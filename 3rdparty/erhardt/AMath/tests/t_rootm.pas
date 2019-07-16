{Main test unit for quadsolve/cubsolve/polyroots,  (c) W.Ehrhardt 2010-2018}
unit t_rootm;

interface

{$i STD.INC}


{$ifdef BIT16}
{$N+}
{$endif}

procedure test_quadsolve;
  {-Quadsolve regression test via single procedure call}

procedure test_cubsolve;
  {-Self test cubsolve}

procedure test_PolyRoots;
  {-Simple tests for PolyRoots/A}

implementation

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  AMath, AMTools;


var
  failed: integer;
  totalfail: integer;
  totalcnt : integer;


{---------------------------------------------------------------------------}
function squad_selftest: integer;
  {-Selftest for squad, Result=0 if OK, else number of the failing test case}
var
  x,y,z,x1,x2,y1,y2: double;
  x0,xz,y0,yz,i: integer;
const
  e: double = 1.0/8192.0;
begin
  {test arithmetik}
  y := (1.0/2147483648.0);
  z := -y*y;
  x := (1.0+y)*(1.0-y); x := x - 1.0;
  y := (1.0+y)*(1.0-y); y := -1.0 + y;
  if x=0 then x0:=1 else x0:=0;
  if y=0 then y0:=1 else y0:=0;
  if x=z then xz:=1 else xz:=0;
  if y=z then yz:=1 else yz:=0;
  i := x0*y0 + 2*xz*yz + 3*x0*yz + 4*xz*y0;
  x := 3.0;
  y := 4.0/x;
  y := y - 1.0;
  y := y * x;
  z := 1.0 - y;
  {z should be 2^-52}
  if (i<>1) or ((z-2.22045E-16) > 1e-20) then begin
    squad_selftest := 1;
    exit;
  end;

  {0*x^2 + 0*x + 1 = 0}
  if squad(0,0,1,x1,y1,x2,y2) <> 0 then begin
    squad_selftest := 2;
    exit;
  end;

  {0*x^2 + 2x + 1 = 0}
  if (squad(0,2,1,x1,y1,x2,y2) <> 1) or (x1<>-0.5) then begin
    squad_selftest := 3;
    exit;
  end;

  {4*x^2 + 3x + 0 = 0}
  if (squad(4,3,0,x1,y1,x2,y2) <> 2) or (x1<>-0.75) or (x2<>0.0) then begin
    squad_selftest := 4;
    exit;
  end;

  {x^2 - 2*x + 1 = 0}
  if (squad(1,-2,1,x1,y1,x2,y2) <> 1) or (x1<>1.0) then begin
    squad_selftest := 5;
    exit;
  end;

  {x^2 - 2*x + (1+e*e) = 0}
  if (squad(1,-2, 1.0+e*e, x1,y1,x2,y2) <> -2) or (x1<>1.0) or (y2<>e) then begin
    squad_selftest := 6;
    exit;
  end;

  {x^2 - 2*x + (1-e*e) = 0}
  if (squad(1,-2, 1.0-e*e, x1,y1,x2,y2) <> 2) or (x1<>1.0-e) or (x2<>1+e) then begin
    squad_selftest := 7;
    exit;
  end;

  {Test cases from Kahan [15]}
  i := squad(94906266.375, -2*94906267.375, 94906268.375, x1,y1,x2,y2);
  if (i <> 2) or (x1<>1.0) or (abs(x2-1.0000000210734241) > 1E-16) then begin
    squad_selftest := 8;
    exit;
  end;

  i := squad(2017810.0*8264001469.0,-2*39213.0*229699315399.0, 45077.0*107932908389.0, x1,y1,x2,y2);
  if (i <> -2) or (abs(x1-0.54015588795707844) > 2e-16) or (abs(y2-0.59969350369679513e-16) > 1E-31) then begin
    squad_selftest := 9;
    exit;
  end;
  squad_selftest := 0;
end;



type
  solvefunc = function(a,b,c: double; var x1,y1,x2,y2: double): integer;

var
  verbose: boolean;

{---------------------------------------------------------------------------}
procedure testnaninf(a,b,c: double);
var
  x1,x2,y1,y2: double;
  n: integer;
begin
  n := squadx(a,b,c,x1,y1,x2,y2);
  inc(totalcnt);
  if n<>0 then begin
    writeln(' *** failed for:');
    writeln(' a = ', a);
    writeln(' b = ', b);
    writeln(' c = ', c);
    writeln(' n = ', n);
    writeln('x1 = ', x1);
    writeln('y1 = ', y1);
    writeln('x2 = ', x2);
    writeln('y2 = ', y2);
    writeln(' n = ', n);
    writeln('x1 = ', x1);
    writeln('y1 = ', y1);
    writeln('x2 = ', x2);
    writeln('y2 = ', y2);
    inc(failed);
    inc(totalfail);
  end;
end;

{---------------------------------------------------------------------------}
procedure test_nan_inf;
begin
  failed := 0;
  writeln('** Test Inf/Nan arguments');
  testnaninf(Nan_d,0,0);
  testnaninf(0,Nan_d,0);
  testnaninf(0,0,Nan_d);
  testnaninf(PosInf_d,0,0);
  testnaninf(0,PosInf_d,0);
  testnaninf(0,0,PosInf_d);
  testnaninf(NegInf_d,0,0);
  testnaninf(0,NegInf_d,0);
  testnaninf(0,0,NegInf_d);
  if failed=0 then writeln('OK.')
  else writeln(' *** number of failed tests: ', failed);
end;

{---------------------------------------------------------------------------}
procedure test(fsolve: solvefunc; const msg: string);
const
  NMax=78;
var
  i,n,nums: integer;
  a,b,c,x1,x2,y1,y2: double;
  maxr1,maxr2,r1,r2: double;
  maxi1,maxi2,i1,i2: double;
  F: array[0..NMax] of double;
begin
  writeln('** Fibonacci test for ', msg);

  failed := 0;

  {Fibonacci test cases from Kahan}
  F[0]  := 0;
  F[1]  := 1;
  maxr1 := 0;
  maxr2 := 0;
  maxi1 := 0;
  maxi2 := 0;
  for n:=2 to NMax do F[n] := F[n-1] + F[n-2];
  {real test cases}
  for i:=1 to NMax div 2 do begin
    inc(totalcnt);
    n := 2*i;
    a := F[n];
    b := -2*F[n-1];
    c := F[n-2];
    nums := fsolve(a,b,c,x1,y1,x2,y2);
    {double real}
    if nums <> +2 then begin
      inc(failed);
      inc(totalfail);
      writeln(n:3,'*** wrong number of solutions: ',nums, ' for n=',n);
    end;
    r1 := abs(x1-(F[n-1]-1)/F[n]);
    r2 := abs(x2-(F[n-1]+1)/F[n]);
    if r1>maxr1 then maxr1 := r1;
    if r2>maxr2 then maxr2 := r2;
    if verbose then writeln(n:3, x1:21:17,' ',r1:13,x2:22:17,' ',r2:13);
  end;
  inc(totalcnt);
  if (maxr1>eps_d) or (maxr2>eps_d) then begin
    inc(failed);
    inc(totalfail);
    writeln('*** Max abs r1:', maxr1:22,'   Max abs r2:', maxr2:22);
  end;
  for i:=2 to NMax div 2 do begin
    inc(totalcnt);
    n := 2*i-1;
    a := F[n];
    b := -2*F[n-1];
    c := F[n-2];
    nums := fsolve(a,b,c,x1,y1,x2,y2);
    {complex case, y1<=0, y2 >=0}
    if nums <> -2 then begin
      inc(failed);
      inc(totalfail);
      writeln(n:3, '**+ wrong number of solutions: ',nums, ' for n=',n);
    end;
    i1 := abs(x1-F[n-1]/F[n]);
    i2 := abs(y2-1/F[n])*F[n];  {relative error}
    if i1>maxi1 then maxi1 := i1;
    if i2>maxi2 then maxi2 := i2;
    if verbose then writeln(n:3,x1:21:17,' ',i1:13,'  ',y2:20,' ',i2:13);
  end;
  inc(totalcnt);
  if (maxi1>eps_d) or (maxi2>eps_d) then begin
    inc(failed);
    inc(totalfail);
    writeln('***  Max abs i1:',maxi1:22,'   Max rel i2:', maxi2:22);
  end;
  if failed=0 then writeln('OK.')
  else writeln(' *** number of failed tests: ', failed);
end;



{---------------------------------------------------------------------------}
procedure singlespec(a,b,c: double; fsolve: solvefunc; en: integer; ex1,ey1,ex2,ey2: double);
var
  n: integer;
  x1,x2,y1,y2: double;
begin
  n := fsolve(a,b,c,x1,y1,x2,y2);
  inc(totalcnt);
  if (n<>en) or (x1<>ex1) or (x2<>ex2) or (y1<>ey1) or (y2<>ey2) then begin
    writeln(' *** failed: a =', a:17, ', b =', b:17, ', c =', c:17);
    writeln(' n  en = ', n ,'   ', en );
    writeln('x1 ex1 = ', x1,'   ', ex1);
    writeln('y1 ey1 = ', y1,'   ', ey1);
    writeln('x2 ex2 = ', x2,'   ', ex2);
    writeln('y2,ey2 = ', y2,'   ', ey2);
    inc(failed);
    inc(totalfail);
  end;
end;


{---------------------------------------------------------------------------}
procedure specialtest(fsolve: solvefunc; all: boolean; const msg: string);
var
  a,b,c: double;
begin
  failed := 0;
  writeln('** Special value test for ',msg);

{Test case H2}
  a := 1e300;
  b := 0;
  c := 1e-100;
  singlespec(a,b,c,fsolve,-2,0,-1e-200,0,1e-200);

{Test case H1}
  a := 1e300;
  b := 0;
  c := -1e-100;
  singlespec(a,b,c,fsolve,2,-1e-200,0,1e-200,0);

{Test case H4}
  a := 1e300;
  b := 1e-100;
  c := 1e-100;
  singlespec(a,b,c,fsolve,-2,0,-1e-200,0,1e-200);

{Test case H5}
  a := -1e300;
  b := 1e-100;
  c := 1e-100;
  singlespec(a,b,c,fsolve,2,-1e-200,0,1e-200,0);

{Test case H3}
  a := 1e-200;
  b := 1e100;
  c := 1e-100;
  singlespec(a,b,c,fsolve,2,-1e300,0,-1e-200,0);

{Test case H3}
  a := 1e-200;
  b := 1e100;
  c := 10;
  singlespec(a,b,c,fsolve,2,-1e300,0,-1e-99,0);

  if all then begin
    a := 1e-300;
    b := 1e100;
    c := 1e-100;
    singlespec(a,b,c,fsolve,2,NegInf_d,0,-1e-200,0);
    a := -1e-300;
    b := 1e100;
    c := 1e-100;
    singlespec(a,b,c,fsolve,2,-1e-200,0,PosInf_d,0);
    a := 1e100;
    b := 1e-200;
    c := 1e300;
    singlespec(a,b,c,fsolve,-2,-5e-301,-1e100,-5e-301,1e100);
  end;
  if failed=0 then writeln('OK.')
  else writeln(' *** number of failed tests: ', failed);
end;


{---------------------------------------------------------------------------}
procedure test_quadsolve;
  {-Quadsolve regression test via single procedure call}
var
  n: integer;
begin
  writeln('Test program for AMTools/quadsolve   (c) W.Ehrhardt 2010-2012');
  failed    := 0;
  totalfail := 0;
  totalcnt  := 0;
  verbose := paramcount<0;
  n := squad_selftest;
  write('squad self test: ');
  if n=0 then writeln('OK')
  else writeln(' failure in test case ',n);

  specialtest({$ifdef FPC_ProcVar}@{$endif}squadx, true, 'squadx');
  specialtest({$ifdef FPC_ProcVar}@{$endif}squad, false, 'squad');
  test_nan_inf;
  test({$ifdef FPC_ProcVar}@{$endif}squadx, 'squadx');
  test({$ifdef FPC_ProcVar}@{$endif}squad, 'squad');

  if totalfail=0 then writeln('All ',totalcnt,' tests passed.')
  else writeln(' *** total number of failed tests: ', totalfail);

end;


{---------------------------------------------------------------------------}
function cubsolve_selftest: integer;
  {-Selftest for cubsove, Result=0 if OK, else number of the failing test case}
const
  tol1 = 1.0e-15;  {~4.5*eps_d}
  tol2 = 1.5e-8;   {~sqrt(eps_d)}
  tol3 = 6.1e-6;   {~cbrt(eps_d)}

  function ccheck(ax,ax1,ay1,ax2,ay2,tx,tx1,ty1,tx2,ty2,tol: double): boolean;
  begin
    ccheck := (abs(ax-tx)>tol) or (abs(ax1-tx1)>tol) or (abs(ax2-tx2)>tol)
               or (abs(ay1-ty1)>tol) or (abs(ay2-ty2)>tol);
  end;

var
  a,b,c,d: double;
  x,y,x1,x2,y1,y2: double;
  tx,tx1,ty1,tx2,ty2: double;
begin

  a := 1; b := -6; c := 11; d := -6;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  if ccheck(x,x1,y1,x2,y2,1,2,0,3,0,tol1) then begin
    cubsolve_selftest := 1;
    exit;
  end;

  a := 1; b := 0; c := 0; d := 1;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  ty1 := -0.8660254037844386468;
  if ccheck(x,x1,y1,x2,y2,-1,0.5,ty1,0.5,-ty1,tol1) then begin
    cubsolve_selftest := 2;
    exit;
  end;

  a := 1; b := -30; c := 299; d := -299e-18;
  y := 8.602325267042626771;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  if ccheck(x,x1,y1,x2,y2,1e-17,15,-y,15,y,tol1) then begin
    cubsolve_selftest := 3;
    exit;
  end;

  {Triple root 1/3 with small inccuracy for d=double(1/9)}
  {will result is errors of order cqrt(eps_d)}
  a := -3; b := 3; c := -1; d := 1/9;
  {d = 0.111111111111111104943205418749130330979824066162109375}
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  tx  :=  0.3333320617675781250;
  tx1 :=  0.3333339691162109375;
  ty1 := -0.1101208246592761544e-5;
  tx2 := tx1;
  ty2 := -ty1;
  if ccheck(x,x1,y1,x2,y2,tx,tx1,ty1,tx2,ty2,tol3) then begin
    cubsolve_selftest := 4;
    exit;
  end;

  {Similar to test 4}
  a := -27; b := 27; c := -9; d := 1;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  if ccheck(x,x1,y1,x2,y2,tx,tx1,ty1,tx2,ty2,tol3) then begin
    cubsolve_selftest := 5;
    exit;
  end;

  {Prevent overflow, large 1/a^3 > MaxDouble; x=-b/a}
  a := 1e-110; b := 1; c := 1; d := -2;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  if ccheck(x,x1,y1,x2,y2,-b/a,-2,0,1,0,tol1) then begin
    cubsolve_selftest := 6;
    exit;
  end;

  {Nothing special, but small a, large x=-b/a}
  a := 1e-60; b := 1; c := 1; d := -2;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  if ccheck(x,x1,y1,x2,y2,-b/a,-2,0,1,0,tol1) then begin
    cubsolve_selftest := 7;
    exit;
  end;

  {similar to test 7}
  a := 1e-10; b := 1; c := 1; d := -2;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
           {1234567890123456789}
  tx  := -0.9999999999000000000e10;
  tx1 := -2.000000000266666667;
  tx2 :=  0.9999999999666666667;
  if ccheck(x,x1,y1,x2,y2,tx,tx1,0,tx2,0,tol1) then begin
    cubsolve_selftest := 8;
    exit;
  end;

  {similar to test 7/8}
  a := 1e-8; b := 1; c := 1; d := -2;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  tx  := -0.9999999899999997000e8;
  tx1 := -2.000000026666667496;
  tx2 :=  0.9999999966666666963;
  if ccheck(x,x1,y1,x2,y2,tx,tx1,0,tx2,0,tol1) then begin
    cubsolve_selftest := 9;
    exit;
  end;

  {nearly quadratic, one root with small absolute value}
  a := 1; b := 2; c := -3; d := 1e-10;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  tx  := -3.000000000008333333;
  tx1 :=  0.3333333333407407407e-10;
  tx2 :=  0.9999999999750000000;
  if ccheck(x,x1,y1,x2,y2,tx,tx1,0,tx2,0,tol1) then begin
    cubsolve_selftest := 10;
    exit;
  end;

  {Test case from Kahan, clustered roots}
  a := 658; b := -190125; c := 18311811; d := -587898164;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  tx  := 96.22963934659218218;
  tx1 := 96.35706482518415207;
  ty1 := -0.6974975204368962543e-1;
  tx2 := tx1;
  ty2 := -ty1;
  if ccheck(x,x1,y1,x2,y2,tx,tx1,ty1,tx2,ty2,tol2) then begin
    cubsolve_selftest := 11;
    exit;
  end;

  {double root x=10, other x=-1}
  a := 1; b := -19; c := 80; d := 100;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  if ccheck(x,x1,y1,x2,y2,-1,10,0,10,0,tol1) then begin
    cubsolve_selftest := 12;
    exit;
  end;

  {clustered roots + double root}
  a := 1; b := -24.125; c := 194; d := -520;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  tx  := 8;
  tx1 := 8;
  tx2 := 8.125;
  if ccheck(x,x1,y1,x2,y2,tx,tx1,0,tx2,0,tol3) then begin
    cubsolve_selftest := 13;
    exit;
  end;

  {From Kahan, two nearly identical roots far from third}
  a := 100000000000000.0; b := -a+2; c := -a; d := -b;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  if ccheck(x,x1,y1,x2,y2,-1,0.99999999999998,0,1,0,tol1) then begin
    cubsolve_selftest := 14;
    exit;
  end;

  {Finally a 'normal' case}
  a := 8; b := -9; c := 10; d := -11;
  cubsolve(a,b,c,d,x,x1,y1,x2,y2);
  tx  := 1.11243727016969651302117199180;
  tx1 := 0.6281364915151743489e-2;
  ty1 := -1.11174875559925879997894189650;
  tx2 := tx1;
  ty2 := -ty1;
  if ccheck(x,x1,y1,x2,y2,tx,tx1,ty1,tx2,ty2,tol1) then begin
    cubsolve_selftest := 15;
    exit;
  end;

  cubsolve_selftest := 0;

end;


{---------------------------------------------------------------------------}
procedure test_cubsolve;
  {-Self test cubsolve}
var
  n: integer;
begin
  writeln('Test program for AMTools/cubsolve   (c) W.Ehrhardt 2013');
  n := cubsolve_selftest;
  write('cubsolve self test: ');
  if n=0 then writeln('passed.')
  else writeln(' failure in test case ',n);
end;


{---------------------------------------------------------------------------}
procedure root2coef(const x,y: array of double; n: integer; var cr,ci: TPolyvec);
  {-Reconstruct coefficients from roots}
var
  i,j: integer;
  zr,zi,tr,ti: double;
  qr,qi: TPolyvec;
begin
  if n<0 then exit;
  cr[0] := 1.0;
  ci[0] := 0.0;
  if n=0 then exit;
  {p[0]=1, p[i]=p[i-1]*(z-z[i])}
  qr[0] := 0;
  qi[0] := 0;
  for i:=1 to n do begin
    {p[i] = p[i-1]*(z-z[i]) = p[i-1]*z - p[i-1]*z[i]}
    {Compute q = p[i-1]*z}
    for j:=0 to i-1 do begin
      qr[j+1] := cr[j];
      qi[j+1] := ci[j];
    end;
    cr[i] := 0.0;
    ci[i] := 0.0;
    zr := x[i];
    zi := y[i];
    for j:=0 to i do begin
      {complex product c[j]*z[i]}
      tr := cr[j]*zr - ci[j]*zi;
      ti := cr[j]*zi + ci[j]*zr;
      cr[j] := qr[j] - tr;
      ci[j] := qi[j] - ti;
    end;
  end;
end;



{$define verbose}
{---------------------------------------------------------------------------}
procedure testroot(cnt,n,maxmult: integer; c,rr,ri: array of double; var OK: boolean);
var
  x,y,coef,cr,ci: TPolyVec;
  i,ierr: integer;
  mr,mc,d,a,h: extended;
  eps: double;
begin
  OK := false;
  eps := 4*nroot(eps_d, maxmult);
  writeln('Test ',cnt:2 );
  if high(c) > 20  then begin
    writeln(' * Invalid coefficient vector for test ');
    exit;
  end;
  {Test both PolyRootsA routines}
  if odd(cnt) then begin
    coef[0] := c[0];  {Avoid warning variable "xx" does not seem to be initialized}
    for i:=1 to n do coef[i] := c[i];
    PolyRootsA(coef,n,x,y,ierr);
  end
  else begin
    PolyRootsOA(c,n,x,y,ierr);
  end;
  if ierr<>0  then begin
    writeln(' * ierr = ' ,ierr);
    exit;
  end;
  mr := -1;
  mc := -1;
  {$ifdef verbose}
   writeln(' Computed roots');
  {$endif}
  for i:=1 to n do begin
   {$ifdef verbose}
    writeln(x[i]:18:12, y[i]:18:12, rr[i]:18:12, ri[i]:18:12);
   {$endif}
    {error absolute if abs. root > 1, else relative}
    a := hypot(rr[i], ri[i]);
    d := hypot(rr[i]-x[i], ri[i]-y[i]);
    if a>1 then d := d/a;
    mr := maxx(mr,d);
  end;

  root2coef(x, y, n, cr, ci);
  {Note: root2coef returna highest coefficient = 1}
  h := c[n];

  {$ifdef verbose}
   writeln(' Reconstructed coefficients');
  {$endif}
  for i:=0 to n do begin
    cr[i] := h*cr[i];
    ci[i] := h*ci[i];
    d := hypot(cr[i]-c[i], ci[i]);
    a := abs(c[i]);
    if a>1 then d := d/a;
    mc := maxx(mc,d);
   {$ifdef verbose}
    writeln(i:3, c[i]:11:2, cr[i]:20:10, '  ', ci[i]:20, '  ', d:20);
   {$endif}
  end;

  writeln(' max-r: ',mr:18, ', max-c: ',mc:18, ', eps: ',eps:18);
  OK := (mr <= eps) and (mc <= eps);
  if OK then writeln(' OK')
  else writeln('  **failed**');
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_PolyRoots;
  {-Simple tests for PolyRoots/A}
const
  coef01: array[0..6] of double = (-0.25,1,1,-1,-1,-1,1);
  rre_01: array[0..6] of double = (0,
                                  -0.72469312227085959948,
                                  -0.51805192664198894946,
                                  -0.51805192664198894946,
                                   0.21597354849329094598,
                                   0.90685904275987088995,
                                   1.6379643843016756625);
  rim_01: array[0..6] of double = (0,
                                   0,
                                   -0.89830421912249224775,
                                   +0.89830421912249224775,
                                   0,
                                   0,
                                   0);

  coef02: array[0..3] of double = (-1,27,-243,729);
  rre_02: array[0..3] of double = (0, 1/9, 1/9, 1/9);
  rim_02: array[0..3] of double = (0,0,0,0);


  coef03: array[0..5] of double = (-120,274,-225,85,-15,1);
  rre_03: array[0..5] of double = (0,1,2,3,4,5);
  rim_03: array[0..5] of double = (0,0,0,0,0,0);

  coef04: array[0..1] of double = (1,2);
  rre_04: array[0..1] of double = (0,-0.5);
  rim_04: array[0..1] of double = (0,0);

  coef05: array[0..4] of double = (2,-2,7,-4,1);
  rre_05: array[0..4] of double = (0,
                                   0.60435090833358811872e-1,
                                   0.60435090833358811872e-1,
                                   1.9395649091666411881,
                                   1.9395649091666411881);
  rim_05: array[0..4] of double = (0,
                                   -0.56432242226560213841,
                                   +0.56432242226560213841,
                                   -1.5643224222656021384,
                                   +1.5643224222656021384);

  coef06: array[0..6] of double = (1,-2,-1,4,-1,-2,1);
  rre_06: array[0..6] of double = (0,-1,-1, 1,1,1,1);
  rim_06: array[0..6] of double = (0,0,0,0,0,0,0);

  coef07: array[0..6] of double = (-1,0,0,0,0,1,1);
  rre_07: array[0..6] of double = ( 0,
                                   -1.28519903324534936790,
                                   -0.67136892027499877778,
                                   -0.67136892027499877778,
                                    0.37333270608088866453,
                                    0.37333270608088866453,
                                    0.88127146163356959440);

  rim_07: array[0..6] of double = (0,
                                    0.0,
                                   -0.7848511587200250234,
                                   +0.7848511587200250234,
                                   -0.8296446029891385578,
                                   +0.8296446029891385578,
                                    0);

  coef08: array[0..8] of double = (0,0,0,-120,274,-225,85,-15,1);
  rre_08: array[0..8] of double = (0,0,0,0,1,2,3,4,5);
  rim_08: array[0..8] of double = (0,0,0,0,0,0,0,0,0);

  coef09: array[0..8] of double =(144,420,340,-63,-183,-46,18,9,1);
  rre_09: array[0..8] of double =(0,-4,-3,-3,-1,-1,-1,2,2);
  rim_09: array[0..8] of double =(0,0,0,0,0,0,0,0,0);

  coef10: array[0..11] of double =(+241920,-213984,-254760,+278476,-6270,
                                    -65857,20790,1023,-1650,+341,-30,1);
  rre_10: array[0..11] of double =(0,-3,-2,-1,1,2,3,4,5,6,7,8);
  rim_10: array[0..11] of double =(0,0,0,0,0,0,0,0,0,0,0,0);

  coef11: array[0..4] of double = (-1, 100000, -100000, -1, 1);
  rre_11: array[0..4] of double = (0,
                                  -316.22776603259947911,
                                   0.10000100002000060002e-4,
                                   0.99998999979999699996,
                                   316.22776603269948011);

  rim_11: array[0..4] of double = (0,0,0,0,0);


var
  OK: boolean;
begin
  writeln;
  writeln('Test program for AMTools/PolyRoots(A)   (c) W.Ehrhardt 2017');
  failed := 0;
  testroot(01,  6, 1, coef01, rre_01, rim_01, OK); if not OK then inc(failed);
  testroot(02,  3, 3, coef02, rre_02, rim_02, OK); if not OK then inc(failed);
  testroot(03,  5, 1, coef03, rre_03, rim_03, OK); if not OK then inc(failed);
  testroot(04,  1, 1, coef04, rre_04, rim_04, OK); if not OK then inc(failed);
  testroot(05,  4, 1, coef05, rre_05, rim_05, OK); if not OK then inc(failed);
  testroot(06,  6, 4, coef06, rre_06, rim_06, OK); if not OK then inc(failed);
  testroot(07,  6, 1, coef07, rre_07, rim_07, OK); if not OK then inc(failed);
  testroot(08,  8, 1, coef08, rre_08, rim_08, OK); if not OK then inc(failed);
  testroot(09,  8, 3, coef09, rre_09, rim_09, OK); if not OK then inc(failed);
  testroot(10, 11, 1, coef10, rre_10, rim_10, OK); if not OK then inc(failed);
  testroot(11,  4, 1, coef11, rre_11, rim_11, OK); if not OK then inc(failed);
  if failed=0 then writeln('All tests passed')
  else writeln(' *** number of failed tests: ', failed);
end;

end.
