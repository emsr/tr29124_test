{Part 2 of regression test for SPECFUN unit  (c) 2010-2018  W.Ehrhardt}

unit t_sfd2;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_carlson;
procedure test_bulirsch;
procedure test_maple;
procedure test_legendre;
procedure test_comp_ellint;
procedure test_heuman_lambda;
procedure test_jacobi_zeta;
procedure test_EllipticKim;
procedure test_EllipticECim;
procedure test_EllipticPiCim;

implementation

uses
  amath, specfun, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_carlson;
var
  y,f: double;
  cnt, failed: integer;
const
  eps = 5e-14;
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Elliptic integrals (Carlson)');

  {The 14 digit literal values are from Carlson's paper [12], because of}
  {the low precision, eps must be choosen as 5e-14 for these test cases.}

  y := ell_rf(1,2,0);
  f := 1.3110287771461;
  testrele( 1, eps, y, f, cnt,failed);

  y := ell_rf(2,3,4);
  f := 0.58408284167715;
  testrele( 2, eps, y, f, cnt,failed);

  y := ell_rc(0,0.25);
  f := pi;
  testrele( 3, eps_d, y, f, cnt,failed);

  y := ell_rc(2.25,2);
  f := ln2;
  testrele( 4, eps_d, y, f, cnt,failed);

  y := ell_rc(0.25,-2);
  f := ln(2)/3.0;
  testrele( 5, eps_d, y, f, cnt,failed);

  y := ell_rj(0,1,2,3);
  f := 0.77688623778582;
  testrele( 6, eps, y, f, cnt,failed);

  y := ell_rj(2,3,4,5);
  f := 0.14297579667157;
  testrele( 7, eps, y, f, cnt,failed);

  y := ell_rj(2,3,4,-0.5);
  f := 0.24723819703052;
  testrele( 8, eps, y, f, cnt,failed);

  y := ell_rj(2,3,4,-5);
  f := -0.12711230042964;
  testrele( 9, eps, y, f, cnt,failed);

  y := ell_rd(0,2,1);
  f := 1.7972103521034;
  testrele(10, eps, y, f, cnt,failed);

  y := ell_rd(2,3,4);
  f := 0.16510527294261;
  testrele(11, eps, y, f, cnt,failed);

  y := ell_rf(0.5,1,0);
  f := 1.8540746773014;
  testrele(12, eps, y, f, cnt,failed);

  y := ell_rf(5e-11,1e-10,1);
  f := 12.36441982979439208156;
  testrele(13, eps_d, y, f, cnt,failed);

  y := ell_rd(5e-11,1e-10,1);
  f := 34.09325949193373672376;
  testrele(14, eps_d, y, f, cnt,failed);

  {------------------------------------------------------------}
  {Carlson R_G, test values computed with MPMath elliprg(x,y,z)}

  y := ell_rg(0,1,1);
  f := Pi_4;
  testrel(15, NE, y, f, cnt,failed);

  y := ell_rg(0,0.5,1);
  f := 0.6753219405238377513;
  testrel(16, NE, y, f, cnt,failed);

  y := ell_rg(4,4,4);
  f := 2;
  testrel(17, NE, y, f, cnt,failed);

  y := ell_rg(0,0,4);
  f := 1;
  testrel(18, NE, y, f, cnt,failed);

  y := ell_rg(0,4,4);
  f := Pi_2;
  testrel(19, NE, y, f, cnt,failed);

  y := ell_rg(1e-9,0,4);
  f := 1 + 0.0000000014926343076522288;
  testrel(20, NE, y, f, cnt,failed);

  y := ell_rg(1,2,3);
  f := 1.401847099990895099;
  testrel(21, NE, y, f, cnt,failed);

  y := ell_rg(0,2,3);
  f := 1.238689348131614792;
  testrel(22, NE, y, f, cnt,failed);

  y := ell_rg(4,2,3);
  f := 1.725503028069227760;
  testrel(23, NE, y, f, cnt,failed);

  y := ell_rg(4,4,1);
  f := 1.709199576156145234;
  testrel(24, NE, y, f, cnt,failed);

  y := ell_rg(1e22, 1e16, 0);
  f := 5.000019485130860306e10;
  testrel(25, NE, y, f, cnt,failed);

  y := ell_rg(1e22, 1e16, 1e10);
  f := 5.000019485151592103e10;
  testrel(26, NE, y, f, cnt,failed);

  {asymptotic}
  y := ell_rg(1e26, 1, 1);
  f := 5e12;
  testrel(27, NE, y, f, cnt,failed);

  y := ell_rg(1e20, 1, 1);
  f := 5.000000000000000001e9;
  testrel(28, NE, y, f, cnt,failed);

  y := ell_rg(1e20, 1, 2);
  f := 5.000000000000000002e9;
  testrel(29, NE, y, f, cnt,failed);

  {double}
  y := ell_rg(1e200, 1e190, 2);
  f := 5.000000003099804957e+99;
  testrel(30, NE, y, f, cnt,failed);

  {asymptotic}
  y := ell_rg(1e300, 1, 2);
  f := 0.5e150;
  testrel(31, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
function bulirsch_el2(x,kc,a,b: extended): extended;
  {-Return the general elliptic integral of the 2nd kind, kc<>0}
  { = integral((a+b*t^2)/[(1+t^2)*sqrt((1+t^2)*(1+kc^2*t^2))], t=0..x)}
const
  ca=2e-10;
  cb=1e-22;
var
  c,d,e,f,g,i,m,p,y,z: extended;
  l: integer;
  done: boolean;
begin
  if x=0.0 then bulirsch_el2 := 0.0
  else if kc=0.0 then begin
    {sfc_el2 := x*(a-b)/sqrt(1+x*x) + b*arcsinh(x)}
    y := x*(a-b)/sqrt(1.0+x*x);
    z := b*arcsinh(x);
    bulirsch_el2 := y+z;
  end
  else begin
    {This is my Pascal translation of the Algol procedure el2 from [10]}
    c := x*x;
    d := 1.0+c;
    p := sqrt((1.0+c*kc*kc)/d);
    d := x/d;
    c := 0.5*d/p;
    z := a-b;
    i := a;
    a := 0.5*(b+a);
    y := abs(1.0/x);
    f := 0.0;
    l := 0;
    m := 1.0;
    kc:= abs(kc);
    done := false;
    repeat
      b := i*kc+b;
      e := m*kc;
      g := e/p;
      d := f*g+d;
      f := c;
      i := a;
      p := g+p;
      c := 0.5*(d/p+c);
      g := m;
      m := kc+m;
      a := 0.5*(b/m+a);
      y := y-e/y;
      if y=0.0 then y := sqrt(e)*cb;
      if abs(g-kc) <= ca*g then done := true
      else begin
        kc := 2.0*sqrt(e);
        l  := l+l;
        if y<0.0 then inc(l);
      end;
    until done;
    if y<0.0 then l := l+1;
    e := (arctan(m/y)+pi*l)*a/m;
    if x<0.0 then e := -e;
    bulirsch_el2 := e+c*z;
  end;
end;



{---------------------------------------------------------------------------}
procedure test_bulirsch;
var
  x,y,f,kc,a,b,p: double;
  cnt, failed, i: integer;
const
  NE1 = 1;
  eps = 5e-16;
  NE2 = 10;      {NE2, NABS used for random tests. Most test case will}
  NABS = 16;     {will pass with smaller values}

  function randx: double;
  var
    r: integer;
  begin
    r := random(8);
    randx := ldexpd(random-0.5,r-4);
  end;

begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Elliptic integrals (Bulirsch/Carlson)');

  {cel1(kc) = EllipticCK(k), test case with Maple V}
  y := cel1(1e-10);
  f := 24.412145291060347459;
  testrel(1, NE1, y, f, cnt,failed);

  y := cel1(0.125);
  f := 3.4754474574968195118;
  testrel(2, NE1, y, f, cnt,failed);

  y := cel1(0.5);
  f := 2.1565156474996432354;
  testrel(3, NE1, y, f, cnt,failed);

  y := cel1(1);
  f := 1.5707963267948966192;
  testrel(4, NE1, y, f, cnt,failed);

  y := cel1(10);
  f := 0.36956373629898746778;
  testrel(5, NE1, y, f, cnt,failed);

  y := cel1(-10);
  f := 0.36956373629898746778;
  testrel(6, NE1, y, f, cnt,failed);

  y := cel1(100);
  f := 0.59915893405069964024e-1;
  testrel(7, NE1, y, f, cnt,failed);

  {test cases from Bulirsch [11], 4.2}
  y := cel(0.1, 4.1, 1.2, 1.1);
  f := 1.5464442694017956;
  testrele(8, eps, y, f, cnt,failed);

  y := cel(0.1,-4.1, 1.2, 1.1);
  f := -6.7687378198360556e-1;
  testrele(9, eps, y, f, cnt,failed);

  {Note: el3 is computed using Carlson functions, all other}
  {Bulirsch functions via AGM and Landen/Bartky transformations.}
  {The following el3 test cases are from Bulirsch [11], 4.1}
  y := el3(1.3, 0.11, 4.21);
  f := 6.6220785847015254e-1;
  testrele(10, eps, y, f, cnt,failed);

  y := el3(1.3, 0.11, 0.82);
  f := 1.1307046442074609;
  testrele(11, eps, y, f, cnt,failed);

  y := el3(1.3, 0.92, 0.71);
  f := 1.0058286266977115;
  testrele(12, eps, y, f, cnt,failed);

  y := el3(1.3, 0.92, 0.23);
  f := 1.1884070823345123;
  testrele(13, eps, y, f, cnt,failed);

  y := el3(1.3, 0.12,-0.11);
  f := 1.7259650355348878;
  testrele(14, eps, y, f, cnt,failed);

  y := el3(1.3, 0.12,-2.11);
  f := 2.4416814520721179e-1;
  testrele(15, eps, y, f, cnt,failed);

  y := el3(1.3, 0.40, 0.1600001);
  f := 1.4004165258366944;
  testrele(16, eps, y, f, cnt,failed);

  y := el3(1.3, 1e-10, 0.82);
  f := 1.1341505395282723;
  testrele(17, eps, y, f, cnt,failed);

  y := el3(1.3e-10, 1e-10, 1e-10);
  f := 1.3000000000000000e-10;
  testrele(18, eps, y, f, cnt,failed);

  y := el3(1.6, 1.90, 9.81);
  f := 3.8572324379967252e-1;
  testrele(19, eps, y, f, cnt,failed);

  y := el3(1.6, 1.90, 1.22);
  f := 7.6656179311956402e-1;
  testrele(20, eps, y, f, cnt,failed);

  y := el3(1.6, 1.90, 0.87);
  f := 8.3210591112618096e-1;
  testrele(21, eps, y, f, cnt,failed);

  y := el3(1.6, 1.90, 0.21);
  f := 1.0521272221906806;
  testrele(22, eps, y, f, cnt,failed);

  y := el3(1.6, 1.90,-0.21);
  f := 1.4730439819954361;
  testrele(23, eps, y, f, cnt,failed);

  y := el3(1.6, 1.90,-4.30);
  f := 2.5476951397193611e-1;
  testrele(24, eps, y, f, cnt,failed);

  y := el3(1.6, 10.1,-1e-5);
  f := 3.9501709821649139e-1;
  testrele(25, eps, y, f, cnt,failed);

  y := el3(1.6, 1.50,2.24999);
  f := 7.0057431688357934e-1;
  testrele(26, eps, y, f, cnt,failed);

  y := el3(1.6, 1e10,1.2);
  f := 2.3734774669772208e-9;
  testrele(27, eps, y, f, cnt,failed);

  y := el3(-1.6, 1e10,1.2);
  f :=-2.3734774669772208e-9;
  testrele(28, eps, y, f, cnt,failed);

  y := el3(1.0, 0.31, 9.9e-2);
  f := 1.0903577921777398;
  testrele(29, eps, y, f, cnt,failed);

  {Tests for kc=0, |kc|=1}
  y := el1(-2,0);
  f := -1.443635475178810342;
  testrel(30, NE1, y, f, cnt,failed);

  y := el1(0.5,0);
  f := 0.4812118250596034475;
  testrel(31, NE1, y, f, cnt,failed);

  y := el1(-2,-1);
  f := -1.107148717794090503;
  testrel(32, NE1, y, f, cnt,failed);

  y := el1(0.5,1);
  f := 0.4636476090008061162;
  testrel(33, NE1, y, f, cnt,failed);

  {Special tests for el2 with small x}
  y := el2(0.75e-4, 1,0,1);
  f := 0.1406249990507812557e-12;
  testrel(34, NE1+1, y, f, cnt,failed);

  y := el2(1e-4, 1,0,1);
  f := 0.3333333293333333762e-12;
  testrel(35, NE1, y, f, cnt,failed);

  y := el2(1e-4, 0,1,1);
  f := 0.9999999983333333408e-4;
  testrel(36, NE1, y, f, cnt,failed);

  y := el2(1e4, 0,1,1);
  f := 9.903487555036128036;
  testrel(37, NE1, y, f, cnt,failed);

  {random tests Carlson [12] formula (4.16)}
  randseed := 1000;
  for i:=1 to 500 do begin
    x  := (random-0.5)*2048.0;
    kc := randx;
    y  := el1(x,kc);
    f  := x*ell_rf(1,1+sqr(kc*x), 1+sqr(x));
    if abs(f)>=1.0 then testrel(1000+i, NE2, y, f, cnt,failed)
    else testabs(1000+i,NABS, y, f, cnt,failed)
  end;

  {random tests Carlson [12] formula (4.17) vs. bulirsch_el2}
  randseed := 2000;
  for i:=1 to 500 do begin
    x  := (random-0.5)*2048.0;
    kc := randx; if kc=0 then kc := 0.5;
    a  := randx;
    b  := randx;
    y  := el2(x,kc,a,b); {This returns the Carlson form}
    f  := bulirsch_el2(x,kc,a,b);
    if abs(f)>=1.0 then testrel(2000+i, NE2, y, f, cnt,failed)
    else testabs(2000+i, NE2, y, f, cnt,failed)
  end;

  {random tests Carlson [12] formula (4.19)}
  randseed := 3000;
  for i:=1 to 500 do begin
    kc := randx; if kc=0 then kc := 0.5;
    a  := randx;
    b  := randx;
    p  := randx;
    y  := cel(kc,p,a,b);
    f  := a*ell_rf(0,sqr(kc), 1) +
          (b-p*a)/3.0*ell_rj(0,sqr(kc),1,p);
    if abs(f)>=1.0 then testrel(3000+i, NE2, y, f, cnt,failed)
    else testabs(3000+i, NABS, y, f, cnt,failed)
  end;

  {random tests Carlson [12] formula (4.19) and cel2(kc,a,b)=cel(kc,1,a,b)}
  randseed := 4000;
  for i:=1 to 500 do begin
    kc := randx; if kc=0 then kc := 0.5;
    a  := randx;
    b  := randx;
    y  := cel2(kc,a,b);
    f  := a*ell_rf(0,sqr(kc), 1) +
          (b-a)/3.0*ell_rj(0,sqr(kc),1,1);
    if abs(f)>=1.0 then testrel(4000+i, NE2, y, f, cnt,failed)
    else testabs(4000+i, NABS, y, f, cnt,failed)
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_maple;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Elliptic integrals (Maple)');
  {---------------------------------------------------------------------}

  y := EllipticCE(-1.0);
  f := 1.5707963267948966192;
  testrel( 1, NE, y, f, cnt,failed);

  y := EllipticCE(0);
  f := 1.0;
  testrel( 2, NE, y, f, cnt,failed);

  y := EllipticCE(1e-10);
  f := 1.0000000000000000001;
  testrel( 3, NE, y, f, cnt,failed);

  y := EllipticCE(0.0009765625);
  f := 1.0000037278026361254;
  testrel( 4, NE, y, f, cnt,failed);

  y := EllipticCE(0.25);
  f := 1.0723027218946042779;
  testrel( 5, NE, y, f, cnt,failed);

  y := EllipticCE(0.5);
  f := 1.2110560275684595248;
  testrel( 6, NE, y, f, cnt,failed);

  y := EllipticCE(1);
  f := 1.5707963267948966192;
  testrel( 7, NE, y, f, cnt,failed);

  y := EllipticCE(5);
  f := 5.2525111349222502362;
  testrel( 8, NE, y, f, cnt,failed);

  y := EllipticCE(100);
  f := 100.02745824306629652;
  testrel( 9, NE, y, f, cnt,failed);

  y := EllipticCE(1e10);
  f := 10000000000.000000001;
  testrel(10, NE, y, f, cnt,failed);


  y := EllipticEC(0);
  f := 1.5707963267948966192;
  testrel(11, NE, y, f, cnt,failed);

  y := EllipticEC(1e-10);
  f := 1.5707963267948966192;
  testrel(12, NE, y, f, cnt,failed);

  y := EllipticEC(0.0009765625);
  f := 1.5707959522878013597;
  testrel(13, NE, y, f, cnt,failed);

  y := EllipticEC(0.25);
  f := 1.5459572561054650350;
  testrel(14, NE, y, f, cnt,failed);

  y := EllipticEC(-0.5);
  f := 1.4674622093394271555;
  testrel(15, NE, y, f, cnt,failed);

  y := EllipticEC(0.9990234375);
  f := 1.0039123555568176367;
  testrel(16, NE, y, f, cnt,failed);

  y := EllipticEC(1);
  f := 1;
  testrel(17, NE, y, f, cnt,failed);


  y := EllipticK(0);
  f := 1.5707963267948966192;
  testrel(21, NE, y, f, cnt,failed);

  y := EllipticK(2e-5);
  f := 1.570796326951976252;
  testrel(22, NE, y, f, cnt,failed);

  y := EllipticK(0.0009765625);
  f := 1.5707967013021258130;
  testrel(23, NE, y, f, cnt,failed);

  y := EllipticK(0.25);
  f := 1.5962422221317835101;
  testrel(24, NE, y, f, cnt,failed);

  y := EllipticK(-0.5);
  f := 1.6857503548125960429;
  testrel(25, NE, y, f, cnt,failed);

  y := EllipticK(0.875);
  f := 2.1854884692782236869;
  testrel(26, NE, y, f, cnt,failed);

  y := EllipticK(0.9990234375);
  f := 4.5074135978990422666;
  testrel(27, NE, y, f, cnt,failed);

  y := EllipticK(1-ldexp(1,-20));
  f := 7.9711961389836741937;
  testrel(28, NE, y, f, cnt,failed);

  y := EllipticK(1);
  f := PosInf_x;
  testabs(29, NE, y, f, cnt,failed);


  y := EllipticCK(0.0);
  f := PosInf_x;
  testabs(31, NE, y, f, cnt,failed);

  y := EllipticCK(1e-300);
  f := 692.16182225933359582;
  testrel(32, NE, y, f, cnt,failed);

  y := EllipticCK(1e-10);
  f := 24.412145291060347459;
  testrel(33, NE, y, f, cnt,failed);

  y := EllipticCK(0.0009765625);
  f := 8.3177679114116699894;
  testrel(34, NE, y, f, cnt,failed);

  y := EllipticCK(0.25);
  f := 2.8012060846652040464;
  testrel(35, NE, y, f, cnt,failed);

  y := EllipticCK(-0.5);
  f := 2.1565156474996432354;
  testrel(36, NE, y, f, cnt,failed);

  y := EllipticCK(1);
  f := 1.5707963267948966192;
  testrel(37, NE, y, f, cnt,failed);

  y := EllipticCK(5);
  f := 0.60322249849552949821;
  testrel(38, NE, y, f, cnt,failed);

  y := EllipticCK(100);
  f := 0.59915893405069964024e-1;
  testrel(39, NE, y, f, cnt,failed);

  y := EllipticCK(1e10);
  f := 0.24412145291060347459e-8;
  testrel(40, NE, y, f, cnt,failed);


  y := EllipticF(0.0,42);
  f := 0;
  testabs(41, NE, y, f, cnt,failed);

  y := EllipticF(0.2, 0.3);
  f := 0.20147959013912327699;
  testrel(42, NE, y, f, cnt,failed);

  y := EllipticF(0.2, 4.0);
  f := 0.23360519703721422990;
  testrel(43, NE, y, f, cnt,failed);

  y := EllipticF(0.2, 0.0);
  f := 0.20135792079033079146;
  testrel(44, NE, y, f, cnt,failed);

  y := EllipticF(1,0.9990234375);
  f := 4.5074135978990422666;
  testrel(45, NE, y, f, cnt,failed);

  y := EllipticF(-0.9990234375,0);
  f := -1.5265985556491813013;
  testrel(46, NE, y, f, cnt,failed);

  y := EllipticF(-0.9990234375,1);
  f := -3.8120652928306447646;
  testrel(47, NE, y, f, cnt,failed);

  y := EllipticF(-0.9990234375,1.0009765625);
  f := -4.4722647220371728115;
  testrel(48, NE, y, f, cnt,failed);

  y := EllipticF(0.0009765625,1023.75);
  f := 0.15127708468293349287e-2;
  testrel(49, NE, y, f, cnt,failed);

  y := EllipticF(0.5,0.5);
  f := 0.52942862705190581774;
  testrel(50, NE, y, f, cnt,failed);


  y := EllipticE(0.0,1.0);
  f := 0;
  testabs(51, NE, y, f, cnt,failed);

  y := EllipticE(0.2, 0.3);
  f := 0.20123638327684510722;
  testrel(52, NE, y, f, cnt,failed);

  y := EllipticE(0.2, -1.0);
  f := 0.2;
  testrel(53, NE, y, f, cnt,failed);

  y := EllipticE(0.2, 0.0);
  f := 0.20135792079033079146;
  testrel(54, NE, y, f, cnt,failed);

  y := EllipticE(1,0.9990234375);
  f := 1.0039123555568176367;
  testrel(55, NE, y, f, cnt,failed);

  y := EllipticE(-0.9990234375,0);
  f := -1.5265985556491813013;
  testrel(56, NE, y, f, cnt,failed);

  y := EllipticE(-0.9990234375,1);
  f := -0.9990234375;
  testrel(57, NE, y, f, cnt,failed);

  y := EllipticE(-0.9990234375,0.0009765625);
  f := -1.5265982022035129847;
  testrel(58, NE, y, f, cnt,failed);

  y := EllipticE(0.0009765625,0.9990234375);
  f := 0.97656250030301704416e-3;
  testrel(59, NE, y, f, cnt,failed);

  {Test Reciprocal-Modulus Transformation}
  y := EllipticE(0.0009765625,1023.75);
  f := 0.76717427171564878748867e-3;
  testrel(60, NE, y, f, cnt,failed);


  y := EllipticCPi(0, 0.25);
  f := 2.8012060846652040464;
  testrel(61, NE, y, f, cnt,failed);

  y := EllipticCPi(-0.2, 1e-10);
  f := 20.5001783026448952759;
  testrel(62, NE, y, f, cnt,failed);

  y := EllipticCPi(-0.2, -0.5);
  f := 1.9386337279430486673;
  testrel(63, NE, y, f, cnt,failed);

  y := EllipticCPi(-0.2, 1);
  f := 1.4339343023863691105;
  testrel(64, NE, y, f, cnt,failed);

  y := EllipticCPi(-0.2, 100);
  f := 0.58146550022385530833e-1;
  testrel(65, NE, y, f, cnt,failed);

  y := EllipticCPi(0.875, 1e-10);
  f := 182.57521136476065425;
  testrel(66, NE, y, f, cnt,failed);

  y := EllipticCPi(0.875, 0.0009765625);
  f := 53.820352555031399325;
  testrel(67, NE, y, f, cnt,failed);

  y := EllipticCPi(0.875, -0.5);
  f := 7.2213865762078484896;
  testrel(68, NE, y, f, cnt,failed);

  y := EllipticCPi(0.875, 1);
  f := 4.44288293815836624702;
  testrel(69, NE, y, f, cnt,failed);

  y := EllipticCPi(0.875, 100);
  f := 0.91912179770862224788e-1;
  testrel(70, NE, y, f, cnt,failed);

  {See also comp_ellint_3}
  y := EllipticPiC(0, 0.25);
  f := 1.5962422221317835101;
  testrel(71, NE, y, f, cnt,failed);

  y := EllipticPiC(-0.2, 1e-10);
  f := 1.4339343023863691105;
  testrel(72, NE, y, f, cnt,failed);

  y := EllipticPiC(-0.2, -0.5);
  f := 1.5338490483665980685;
  testrel(73, NE, y, f, cnt,failed);

  y := EllipticPiC(-0.2, 0.9990234375);
  f := 3.9124854945489819052;
  testrel(74, NE, y, f, cnt,failed);

  y := EllipticPiC(1.2, 0.3);
  f := -0.64689132067447156769e-1;
  testrel(75, NE, y, f, cnt,failed);

  y := EllipticPiC(0.875, 1e-10);
  f := 4.4428829381583662470;
  testrel(76, NE, y, f, cnt,failed);

  y := EllipticPiC(0.875, 0.0009765625);
  f := 4.4428845033223313898;
  testrel(77, NE, y, f, cnt,failed);

  y := EllipticPiC(0.875, -0.5);
  f := 4.9359533239463205273;
  testrel(78, NE, y, f, cnt,failed);

  y := EllipticPiC(0.875, 0.9990234375);
  f := 23.458474987230783737;
  testrel(79, NE, y, f, cnt,failed);

  y := EllipticPiC(1.875, 0.9990234375);
  f := -3.6890037635802376621;
  testrel(80, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, 0, 0.25);
  f := 0.25284646724525372005;
  testrel(81, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, -0.2, 1e-10);
  f := 0.25162629824102647439;
  testrel(82, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, -0.2, -0.5);
  f := 0.25228970720670885550;
  testrel(83, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, -0.2, 1);
  f := 0.25433835602752515312;
  testrel(84, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, -0.2, 3.9990234375);
  f := 0.39103199305948838131;
  testrel(85, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, 0.875, 1e-10);
  f := 0.25748523549627356454;
  testrel(86, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, 0.875, 0.0009765625);
  f := 0.25748523811478012385;
  testrel(87, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, 0.875, -0.5);
  f := 0.25817658256838196729;
  testrel(88, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, 0.875, 1);
  f := 0.26031199088878030704;
  testrel(89, NE, y, f, cnt,failed);

  y := EllipticPi(0.25, 0.875, 3.9990234375);
  f := 0.40460663764275244084;
  testrel(90, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, 0, 0.25);
  f := -0.85365499569974914999;
  testrel(91, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, -0.2, 1e-10);
  f := -0.81516418328369180560;
  testrel(92, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, -0.2, -0.5);
  f := -0.83720106817561389230;
  testrel(93, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, -0.2, 1);
  f := -0.93140141390649194717;
  testrel(94, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, -0.2, 1.25);
  f := -1.0836375159548058208;
  testrel(95, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, 0.875, 1e-10);
  f := -1.0784085177177721527;
  testrel(96, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, 0.875, 0.0009765625);
  f := -1.0784086432466533076;
  testrel(97, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, 0.875, -0.5);
  f := -1.1138851933974080516;
  testrel(98, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, 0.875, 1);
  f := -1.2704072161399052427;
  testrel(99, NE, y, f, cnt,failed);

  y := EllipticPi(-0.75, 0.875, 1.25);
  f := -1.5415123259785248943;
  testrel(100, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_legendre;
var
  x,y,f,nu,k,r: double;
  i,lf,cnt,failed: integer;
const
  NE1  = 1;
  NE1a = 2;    {for k=0.99 and phi ~ n*pi_2}
  NE2  = 1;
  NE3  = 1;
  NE3a = 2;
  NEd  = 1;
  NEda = 8;
  NEb  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Elliptic integrals (Legendre)');

(*
e1 := (x,k)   -> signum(x)*int(1/sqrt(1-k^2*sin(t)^2),t=0..abs(x));
e2 := (x,k)   -> signum(x)*int(sqrt(1-k^2*sin(t)^2),t=0..abs(x));
e3 := (x,n,k) -> signum(x)*int(1/(1+n*sin(t)^2)/sqrt(1-k^2*sin(t)^2),t=0..abs(x));

Maple V will give wrong answers, if the normal integrals are evaluated with
x < 0. Therefore signum(x) and abs(x) are used.
*)

  y := ellint_2(10,1);
  f := 6.544021110889369813;
  testrel(1, NE2, y, f, cnt,failed);

  y := ellint_2(-2,1);
  f := -1.090702573174318305;
  testrel(2, NE2, y, f, cnt,failed);

  y := ellint_2(-100,1);
  f := -63.49363435889024121;
  testrel(3, NE2, y, f, cnt,failed);

  x := Pi_2;
  y := ellint_2(x,0.25);
  f := 1.545957256105465035;
  testrel(4, NE2, y, f, cnt,failed);

  x := -Pi_2;
  y := ellint_2(x,0.99);
  f := -1.028475809028804001;
  testrel(5, NE2, y, f, cnt,failed);

  x := 18998197087967.0/1099511627776.0;
  y := ellint_2(x,0.99);
  f := 11.313233899316843665;
  testrel(6, NE2, y, f, cnt,failed);

  y := ellint_2(x,0.25);
  f := 17.00552981716011301;
  testrel(7, NE2, y, f, cnt,failed);

  x := 1+8*Pi;
  y := ellint_2(x,0.99);
  f := 17.30089027458527057;
  testrel(8, NE2, y, f, cnt,failed);

  x := 1+10*Pi;
  y := ellint_2(x,0.25);
  f := 31.91056242887723082;
  testrel(9, NE2, y, f, cnt,failed);

  y := ellint_2(1,0.99);
  f := 0.8452773301244065537;
  testrel(10, NE2, y, f, cnt,failed);

  y := ellint_2(1,0.25);
  f := 0.991417306767930117;
  testrel(11, NE2, y, f, cnt,failed);

  x := 1+Pi_2;
  y := ellint_2(x,0.99);
  f := 1.516010463169232204;
  testrel(12, NE2, y, f, cnt,failed);

  x := 1+Pi_2;
  y := ellint_2(-x,0.25);
  f := -2.522938091155351336;
  testrel(13, NE2, y, f, cnt,failed);

  x := succd(3*pi_2);
  y := ellint_2(x,0.99);
  f := 3.085427427086412003;
  testrel(14, NE2, y, f, cnt,failed);

  x := succd(pi_2);
  y := ellint_2(x,0.99);
  f := 1.028475809028804001;
  testrel(15, NE2, y, f, cnt,failed);

  lf := 0;
  k  := 0.25;
  for i:=-500 to 500 do begin
    x := i*Pi_2;
    y := ellint_2(x,k);
    f := i*EllipticEC(k);
    if f=0 then r := y-f
    else r := 1-y/f;
    if abs(r) > eps_d then inc(lf);
  end;
  inc(cnt);
  if lf<>0 then begin
    inc(failed);
    writeln(' Failures for ellint_2 n*Pi/2 test: ',lf);
  end;

  {-----------------------------------------------------}
  x := 1;
  y := ellint_1(x,0.99);
  f := 1.218113986212965775;
  testrel(16, NE1, y, f, cnt,failed);

  x := Pi_2;
  y := ellint_1(x,0.99);
  f := 3.356600523361192376;
  {If calculated with the 'real' double values the relative error is zero}
  {x = 1.5707963267948965579989817342720925807952880859375}
  {k = 0.9899999999999999911182158029987476766109466552734}
  {f = 3.3566005233611915084158314894008396005010485965489}
  testrel(17, NE1a, y, f, cnt,failed);

  x := succd(3*Pi_2);
  y := ellint_1(x,0.99);
  f := 10.06980157008357713;
  testrel(18, NE1a, y, f, cnt,failed);

  x := 1+Pi_2;
  y := ellint_1(x,0.99);
  f := 6.109404184682513357;
  testrel(19, NE1, y, f, cnt,failed);

  y := ellint_1(1,1);
  f := 1.226191170883517071;
  testrel(20, NE1, y, f, cnt,failed);

  y := ellint_1(-1.5,1);
  f := -3.340677542798311003;
  testrel(21, NE1, y, f, cnt,failed);

  x := 18998197087967.0/1099511627776.0;
  y := ellint_1(x,0.99);
  f := 36.92260575697309877;
  testrel(22, NE1, y, f, cnt,failed);

  x := 12345*Pi_2;
  y := ellint_1(x,0.99);
  f := 41437.23346089391988;
  testrel(23, NE1a, y, f, cnt,failed);

  x := -10;
  y := ellint_1(x,0.99);
  f := -20.74864531838926736;
  testrel(24, NE1, y, f, cnt,failed);

  x := -100;
  y := ellint_1(x,0.5);
  f := -107.3509831176952950;
  testrel(25, NE1, y, f, cnt,failed);

  x := 10;
  y := ellint_1(x,0.1);
  f := 10.02399268596734540;
  testrel(26, NE1, y, f, cnt,failed);

  x := 100;
  y := ellint_1(x,0.01);
  f := 100.0025110575700321;
  testrel(27, NE1, y, f, cnt,failed);

  x := 1000;
  y := ellint_1(x,0.001);
  f := 1000.000249883885592;
  testrel(28, NE1, y, f, cnt,failed);

  x := 123;
  y := ellint_1(x,0.0);
  f := x;
  testrel(29, NE1, y, f, cnt,failed);

  lf := 0;
  k  := 0.25;
  for i:=-500 to 500 do begin
    x := i*Pi_2;
    y := ellint_1(x,k);
    f := i*EllipticK(k);
    if f=0 then r := y-f
    else r := 1-y/f;
    if abs(r) > eps_d then inc(lf);
  end;
  inc(cnt);
  if lf<>0 then begin
    inc(failed);
    writeln(' Failures for ellint_1 n*Pi/2 test: ',lf);
  end;

  x := -123;
  nu:= 0;
  k := 0.0;
  f := x;
  y := ellint_3(x,nu,k);
  testrel(30, NE3, y, f, cnt,failed);

  x := 100;
  nu:= 1;
  k := 0.0;
  f := -0.5872139151569290767;
  y := ellint_3(x,nu,k);
  testrel(31, NE3, y, f, cnt,failed);

  x := 50;
  nu:= 2;
  k := 0.0;
  f := -0.2789150223828393740;
  y := ellint_3(x,nu,k);
  testrel(32, NE3, y, f, cnt,failed);

  x := 50;
  nu:= 0.5;
  k := 0.0;
  f := 70.81750426633380247;
  y := ellint_3(x,nu,k);
  testrel(33, NE3, y, f, cnt,failed);

  x := pi/3;
  nu:= 0;
  k := 0.5;
  f := 1.089550670051885409;
  y := ellint_3(x,nu,k);
  testrel(34, NE3, y, f, cnt,failed);

  x := pi;
  nu:= -0.5;
  k := 0.5;
  f := 2.732947906009193789;
  y := ellint_3(x,nu,k);
  testrel(35, NE3, y, f, cnt,failed);

  x := 7*pi_2;
  nu:= -0.5;
  k := 0.5;
  f := 9.565317671032178262;
  y := ellint_3(x,nu,k);
  testrel(36, NE3, y, f, cnt,failed);

  x := pi_2;
  nu:= -0.5;
  k := 0.5;
  f := 1.366473953004596895;
  y := ellint_3(x,nu,k);
  testrel(37, NE3, y, f, cnt,failed);

  x := pi/3;
  nu:= -0.5;
  k := 0.5;
  f := 0.9570574331323584890;
  y := ellint_3(x,nu,k);
  testrel(38, NE3, y, f, cnt,failed);

  x := -7*pi/3;
  nu:= 0.1;
  k := 0.01;
  f := -7.702732080131895647;
  y := ellint_3(x,nu,k);
  testrel(39, NE3, y, f, cnt,failed);

  x := 0.5;
  nu:= 2.0;
  k := 0.99;
  y := ellint_3(x,nu,k);
  f := 0.6443632408880116653;
  testrel(40, NE3, y, f, cnt,failed);

  x := 100;
  nu:= 2.0;
  k := 0.99;
  y := ellint_3(x,nu,k);
  f := -132.9522530560008481;
  testrel(41, NE3, y, f, cnt,failed);

  x := Pi/6;
  nu:= 0.25;
  k := 0.5;
  y := ellint_3(x,nu,k);
  f := 0.5414380479799314085;
  testrel(42, NE3, y, f, cnt,failed);

  {Note: extended test case x=Pi/6 nu=4 with nu*sin(x)^2 - 1 = 0 does not}
  {work because with double x,nu the LHS is not zero but of order eps_d.}
  x := 0.5234375;
  nu:= 4.001953125;
  k := 0.5;
  y := ellint_3(x,nu,k);
  f := 3.15877925212406674;
  testabs(43, 0, y, f, cnt,failed);

  x := Pi/6;
  nu:= 5;
  k := 0.5;
  y := ellint_3(x,nu,k);
  f := 0.6617540054437447339;
  testrel(44, NE3a, y, f, cnt,failed);

  x := 100.0;
  nu:= 5;
  k := 0.5;
  y := ellint_3(x,nu,k);
  f := -3.537900009163700912;
  testrel(45, NE3, y, f, cnt,failed);

  x  := 1e20;
  nu := -1e10;
  k  := 0.875;
  f  := 1.000006266655673326e15;
  y  := ellint_3(x,nu,k);
  testrel(46, NE3, y, f, cnt,failed);

  x  := 0.1;
  nu := 50;
  k  := 0.25;
  f  := 0.1245737703427495254;
  y  := ellint_3(x,nu,k);
  testrel(47, NE3, y, f, cnt,failed);

  x  := 10;
  nu := 25;
  k  := 0.75;
  f  := -0.843500077550270004e-1;
  y  := ellint_3(x,nu,k);
  testrel(48, NE3, y, f, cnt,failed);

  x  := 10;
  nu := -1000;
  k  := 0.75;
  f  := 0.3496229795809743825;
  y  := ellint_3(x,nu,k);
  testrel(49, NE3, y, f, cnt,failed);

  x  := 10;
  nu := -128000;
  k  := 0.75;
  f  := 0.3075042579740085930e-1;
  y  := ellint_3(x,nu,k);
  testrel(50, NE3, y, f, cnt,failed);

  x  := 10;
  nu := 1000;
  k  := 0.75;
  f  := -0.2180671013195822018e-2;
  y  := ellint_3(x,nu,k);
  testrel(51, NE3, y, f, cnt,failed); {!!!}

  x  := 10;
  nu := 100;
  k  := 0.75;
  f  := -0.2165904609207162138e-1;
  y  := ellint_3(x,nu,k);
  testrel(52, NE3, y, f, cnt,failed);   {!!!}

  lf := 0;
  nu := 2;
  k  := 0.75;
  for i:=-500 to 500 do begin
    x := i*Pi_2;
    y := ellint_3(x,nu,k);
    f := i*comp_ellint_3(nu,k);
    if f=0 then r := y-f
    else r := 1-y/f;
    if abs(r) > 4*eps_d then inc(lf);
  end;
  inc(cnt);
  if lf<>0 then begin
    inc(failed);
    writeln(' Failures for ellint_3 n*Pi/2 test: ',lf);
  end;

  {------ k > 1 ----------}
  y := ellint_1(0.515625,2);
  f := 0.7468461959637924012;
  testrel(53, NE1, y, f, cnt,failed);

  y := ellint_1(0.5,2);
  f := 0.6774175382039303866;
  testrel(54, NE1, y, f, cnt,failed);

  y := ellint_1(0.25,3);
  f := 0.2820207819635217412;
  testrel(55, NE1, y, f, cnt,failed);

  y := ellint_1(0.125,5);
  f := 0.1349839123459701849;
  testrel(56, NE1, y, f, cnt,failed);

  f := 0.52548500853235841526;
  y := ellint_2(0.625,1.5);
  testrel(57, NE2, y, f, cnt,failed);

  y := ellint_2(0.515625,2);
  f := 0.4054166288690934782;
  testrel(58, NE2, y, f, cnt,failed);

  y := ellint_2(0.5,2);
  f := 0.4018194805534948653;
  testrel(59, NE2, y, f, cnt,failed);

  y := ellint_2(0.25,3);
  f := 0.2244039958745628351;
  testrel(60, NE2, y, f, cnt,failed);

  y := ellint_2(0.125,5);
  f := 0.1163322348977956058;
  testrel(61, NE2, y, f, cnt,failed);

  y := ellint_2(0.125,8);
  f := 0.983518903720969875e-1;
  testrel(62, NE2, y, f, cnt,failed);

  y := ellint_2(-1/128,100);
  f := -0.6921694843063938017e-2;
  testrel(63, NE2, y, f, cnt,failed);

  {Complete case sin(Pi/6)*k=1, arg. x=Pi/6 crashes due to }
  {rounding, predd(Pi/6) is near enough and does not crash!}
  x := predd(Pi/6);
  k := 2;
  f := 0.4062988864599602466;
  y := ellint_2(x,k);
  testrel(64, NE2, y, f, cnt,failed);

  y := ellint_3(0.5234375,0.5,2);
  f := 0.8876398747412474324;
  testrel(65, NE3, y, f, cnt,failed);

  y := ellint_3(0.5,0.5,2);
  f := 0.7148674793327006503;
  testrel(66, NE3, y, f, cnt,failed);

  y := ellint_3(0.125,100,5);
  f := 0.10567011217042676985;
  testrel(67, NE3, y, f, cnt,failed);

  y := ellint_3(0.625,1,1.5);
  f := 0.9161817695423524898;
  testrel(68, NE3, y, f, cnt,failed);

  y := ellint_3(0.7,2,1.25);
  f := 1.573815248832902577;
  testrel(69, NE3, y, f, cnt,failed);

  y := ellint_3(0.125,100,6);
  f := 0.1000924463506501474;
  testrel(70, NE3, y, f, cnt,failed);

  y := ellint_3(0.125,0,6);
  f := 0.1412592900613043378;
  testrel(71, NE3, y, f, cnt,failed);

  y := ellint_3(0.2,6,4);
  f := 0.2558498621248563402;
  testrel(72, NE3, y, f, cnt,failed);

  y := ellint_3(-0.0625,2,10);
  f := -0.6769580807404474460e-1;
  testrel(73, NE3, y, f, cnt,failed);

  y := ellint_d(1,0.5);
  f := 0.2899186629341992247;
  testrel(74, NEd, y, f, cnt,failed);

  y := ellint_d(0.5,2);
  f := 0.6889951441260888032e-1;
  testrel(75, NEd, y, f, cnt,failed);

  y := ellint_d(-0.25,-4);
  f := -0.1026022707357443431e-1;
  testrel(76, NEd, y, f, cnt,failed);

  y := ellint_d(0.1,10);
  f := 0.7304072463625569414e-3;
  testrel(77, NEda, y, f, cnt,failed);

  y := ellint_d(1,1);
  f := 0.3847201860756205642;
  testrel(78, NEd, y, f, cnt,failed);

  y := ellint_d(1.5,1);
  f := 2.343182556194256572;
  testrel(79, NEd, y, f, cnt,failed);

  y := ellint_d(1e-10,1);
  f := 0.3333333333333333333e-30;
  testrel(80, NEd, y, f, cnt,failed);

  y := ellint_d(1e-5,1);
  f := 0.3333333333366666667e-15;
  testrel(81, NEd, y, f, cnt,failed);

  y := ellint_d(-10,0.5);
  f := -5.299691450157785580;
  testrel(82, NEd, y, f, cnt,failed);

  y := ellint_d(2.5*Pi,0.5);
  f := 4.365762909463377748;
  testrel(83, NEd, y, f, cnt,failed);

  y := ellint_d(pi,0.5);
  f := 1.746305163785351099;
  testrel(84, NEd, y, f, cnt,failed);

  y := ellint_d(100,0.25);
  f := 51.44441447208787025;
  testrel(85, NEd, y, f, cnt,failed);

  y := ellint_d(-1e-5,0);
  f := -0.3333333333266666667e-15;
  testrel(86, NEda, y, f, cnt,failed);

  y := ellint_d(-1/1024,0);
  f := -0.3104407989932702365e-9;
  testrel(87, NEd, y, f, cnt,failed);

  y := ellint_d(0.125,0);
  f := 0.6490101863692676008e-3;
  testrel(88, NEd, y, f, cnt,failed);

  y := ellint_d(0.25,0);
  f := 0.5143615348949249932e-2;
  testrel(89, NEd, y, f, cnt,failed);

  lf := 0;
  k  := 0.25;
  for i:=-500 to 500 do begin
    x := i*Pi_2;
    y := ellint_d(x,k);
    f := i*comp_ellint_d(k);
    if f=0 then r := y-f
    else r := 1-y/f;
    if abs(r) > 2*eps_d then inc(lf);
  end;
  inc(cnt);
  if lf<>0 then begin
    inc(failed);
    writeln(' Failures for ellint_d n*Pi/2 test: ',lf);
  end;

  y := ellint_b(10, 1);
  f := 6.544021110889369813;
  testrel(90, NEb, y, f, cnt,failed);

  y := ellint_b(10, 0);
  f := 5.228236312681906914;
  testrel(91, NEb, y, f, cnt,failed);

  y := ellint_b(-20, 1);
  f := -12.91294525072762765;
  testrel(92, NEb, y, f, cnt,failed);

  y := ellint_b(-20, 0);
  f := -10.18627829011983720;
  testrel(93, NEb, y, f, cnt,failed);

  y := ellint_b(1e-6, 0.5);
  f := 0.9999999999997083333e-6;
  testrel(94, NEb, y, f, cnt,failed);

  y := ellint_b(0.5, 0.75);
  f := 0.4705492282917460816;
  testrel(95, NEb, y, f, cnt,failed);

  y := ellint_b(100, 0.75);
  f := 54.39222624285895632;
  testrel(96, NEb, y, f, cnt,failed);

  y := ellint_b(0.5, 1.25);
  f := 0.4924309497019432622;
  testrel(97, NEb, y, f, cnt,failed);

  y := ellint_b(0.5, -2);
  f := 0.6085180237913215063;
  testrel(98, NEb, y, f, cnt,failed);

  y := ellint_b(0.125, 5);
  f := 0.1342378452480432018;
  testrel(99, NEb, y, f, cnt,failed);

  y := 1-ldexp(1,-30);
  y := ellint_b(1e-5, y);
  f := 0.9999999999833333333e-5;
  testrel(100, NEb, y, f, cnt,failed);

  lf := 0;
  k  := 0.25;
  for i:=-500 to 500 do begin
    x := i*Pi_2;
    y := ellint_b(x,k);
    f := i*comp_ellint_b(k);
    if f=0 then r := y-f
    else r := 1-y/f;
    if abs(r) > 2*eps_d then inc(lf);
  end;
  inc(cnt);
  if lf<>0 then begin
    inc(failed);
    writeln(' Failures for ellint_b n*Pi/2 test: ',lf);
  end;

  {cnt includes n*Pi/2 tests}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_comp_ellint;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Complete elliptic integrals (Legendre)');

  y := comp_ellint_1(0);
  f := 1.570796326794896619;
  testrel(1, NE, y, f, cnt,failed);

  y := comp_ellint_1(2e-5);
  f := 1.570796326951976252;
  testrel(2, NE, y, f, cnt,failed);

  y := comp_ellint_1(0.0009765625);
  f := 1.570796701302125813;
  testrel(3, NE, y, f, cnt,failed);

  y := comp_ellint_1(0.25);
  f := 1.596242222131783510;
  testrel(4, NE, y, f, cnt,failed);

  y := comp_ellint_1(-0.5);
  f := 1.685750354812596043;
  testrel(5, NE, y, f, cnt,failed);

  y := comp_ellint_1(0.875);
  f := 2.185488469278223687;
  testrel(6, NE, y, f, cnt,failed);

  y := comp_ellint_1(0.9990234375);
  f := 4.507413597899042267;
  testrel(7, NE, y, f, cnt,failed);

  y := comp_ellint_1(1-ldexp(1,-20));
  f := 7.971196138983674194;
  testrel(8, NE, y, f, cnt,failed);

  y := comp_ellint_1(1);
  f := PosInf_d;
  testabs(9, NE, y, f, cnt,failed);

  y := comp_ellint_2(1);
  f := 1;
  testrel(10, NE, y, f, cnt,failed);

  y := comp_ellint_2(0);
  f := 1.570796326794896619;
  testrel(11, NE, y, f, cnt,failed);

  y := comp_ellint_2(1e-10);
  f := 1.570796326794896619;
  testrel(12, NE, y, f, cnt,failed);

  y := comp_ellint_2(0.0009765625);
  f := 1.570795952287801360;
  testrel(13, NE, y, f, cnt,failed);

  y := comp_ellint_2(0.25);
  f := 1.545957256105465035;
  testrel(14, NE, y, f, cnt,failed);

  y := comp_ellint_2(-0.5);
  f := 1.467462209339427156;
  testrel(15, NE, y, f, cnt,failed);

  y := comp_ellint_2(0.9990234375);
  f := 1.003912355556817637;
  testrel(16, NE, y, f, cnt,failed);

  y := comp_ellint_3(0, 0.25);
  f := 1.5962422221317835101;
  testrel(17, NE, y, f, cnt,failed);

  y := comp_ellint_3(-0.2, 1e-10);
  f := 1.4339343023863691105;
  testrel(18, NE, y, f, cnt,failed);

  y := comp_ellint_3(-0.2, -0.5);
  f := 1.533849048366598069;
  testrel(19, NE, y, f, cnt,failed);

  y := comp_ellint_3(-0.2, 0.9990234375);
  f := 3.912485494548981905;
  testrel(20, NE, y, f, cnt,failed);

  y := comp_ellint_3(1.2, 0.3);
  f := -0.6468913206744715677e-1;
  testrel(21, NE, y, f, cnt,failed);

  y := comp_ellint_3(0.875, 1e-10);
  f := 4.442882938158366247;
  testrel(22, NE, y, f, cnt,failed);

  y := comp_ellint_3(0.875, 0.0009765625);
  f := 4.442884503322331390;
  testrel(23, NE, y, f, cnt,failed);

  y := comp_ellint_3(0.875, -0.5);
  f := 4.935953323946320527;
  testrel(24, NE, y, f, cnt,failed);

  y := comp_ellint_3(0.875, 0.9990234375);
  f := 23.45847498723078374;
  testrel(25, NE, y, f, cnt,failed);

  y := comp_ellint_3(1.875, 0.9990234375);
  f := -3.689003763580237662;
  testrel(26, NE, y, f, cnt,failed);

  y := comp_ellint_d(0.5);
  f := 0.8731525818926755496;
  testrel(27, NE, y, f, cnt,failed);

  y := comp_ellint_d(2e-5);
  f := 0.7853981635152580342;
  testrel(28, NE, y, f, cnt,failed);

  y := comp_ellint_d(1e-9);
  f := 0.785398163397448309910;
  testrel(29, NE, y, f, cnt,failed);

  y := comp_ellint_d(0);
  f := 0.7853981633974483096;
  testrel(20, NE, y, f, cnt,failed);

  y := comp_ellint_d(0.998046875);
  f := 3.167672944303718706;
  testrel(31, NE, y, f, cnt,failed);

  y := comp_ellint_b(0.75);
  f := 0.8576250291834590647;
  testrel(32, NE, y, f, cnt,failed);

  y := 1-ldexp(1,-30);
  y := comp_ellint_b(y);
  {f:= 0.999999990745514169;}
  f := 0.4999999953727570847*2;
  testrel(33, NE, y, f, cnt,failed);

  {extended}
  (*
  y := 1-ldexp(1,-50);
  y := comp_ellint_b(y);
  {f:= 0.999999999999985018;}
  f := 0.4999999999999925089*2;
  testrel(34, NE, y, f, cnt,failed);
  *)

  {double}
  y := 1-ldexpd(1,-40);
  y := comp_ellint_b(y);
  f := 0.999999999987810348;
  testrel(34, NE, y, f, cnt,failed);

  y := comp_ellint_b(0.998046875);
  f := 0.9947873958437291932;
  testrel(35, NE, y, f, cnt,failed);

  y := comp_ellint_b(2e-5);
  f := 0.7853981634367182178;
  testrel(36, NE, y, f, cnt,failed);

  y := comp_ellint_b(0);
  f := Pi_4;
  testrel(37, NE, y, f, cnt,failed);

  y := comp_ellint_b(1);
  f := 1;
  testrel(38, NE, y, f, cnt,failed);

  {E/K for modulus > 1}
  y := comp_ellint_1(10);
  f := 0.1574745561517355953;
  testrel(39, NE, y, f, cnt,failed);

  y := comp_ellint_1(1+ldexp(1,-30));
  f := 11.43692847414619341;
  testrel(40, NE, y, f, cnt,failed);

  y := comp_ellint_1(1+ldexp(1,-40));
  f := 14.902664382032274579;
  testrel(41, NE, y, f, cnt,failed);

  y := comp_ellint_1(1+ldexp(1,-50));
  f := 18.36840028483854276;
  testrel(42, NE, y, f, cnt,failed);

  y := comp_ellint_1(1e10);
  f := 0.1570796326794896619e-9;
  testrel(43, NE, y, f, cnt,failed);

  y := comp_ellint_1(5e4);
  f := 0.3141592653903952504e-4;
  testrel(44, NE, y, f, cnt,failed);

  y := comp_ellint_2(10);
  f := 0.7863836119485898078e-1;
  testrel(45, NE, y, f, cnt,failed);

  y := comp_ellint_2(1+ldexp(1,-30));
  f := 0.9999999898141916125;
  testrel(46, NE, y, f, cnt,failed);

  y := comp_ellint_2(1+ldexp(1,-50));
  f := 0.9999999999999841297;
  testrel(47, NE, y, f, cnt,failed);

  y := comp_ellint_2(1e10);
  f := 0.7853981633974483096e-10;
  testrel(48, NE, y, f, cnt,failed);

  y := comp_ellint_2(5e4);
  f := 0.1570796326873436436e-4;
  testrel(49, NE, y, f, cnt,failed);

  {B/D for modulus > 1}
  y := comp_ellint_b(8);
  f := 0.1955802886578196118;
  testrel(50, NE, y, f, cnt,failed);

  y := comp_ellint_d(8);
  f := 0.1543057744157030961e-2;
  testrel(51, NE, y, f, cnt,failed);

  y := comp_ellint_b(1+ldexp(1,-30));
  f := 1.000000009254485800;
  testrel(52, NE, y, f, cnt,failed);

  y := comp_ellint_d(1+ldexp(1,-30));
  f := 10.43692846489170761;
  testrel(53, NE, y, f, cnt,failed);

  y := comp_ellint_b(1e9);
  f := 0.1570796326794896619e-8;
  testrel(54, NE, y, f, cnt,failed);

  y := comp_ellint_d(1e9);
  f := 0.7853981633974483196e-27;
  testrel(55, NE, y, f, cnt,failed);

  {Pi for modulus > 1}
  y := comp_ellint_3(-0.2, 1.0078125);
  f := 3.037395489367734302;
  testrel(56, NE, y, f, cnt,failed);

  y := comp_ellint_3(-0.2, 2);
  f := 0.8218416369272520393;
  testrel(57, NE, y, f, cnt,failed);

  y := comp_ellint_3(-0.2, 10);
  f := 0.1573171199757702665;
  testrel(58, NE, y, f, cnt,failed);

  y := comp_ellint_3(0.875, 1.0078125);
  f := 14.45363308160301378;
  testrel(59, NE, y, f, cnt,failed);

  y := comp_ellint_3(0.875, 2);
  f := 0.9578412088487968247;
  testrel(60, NE, y, f, cnt,failed);

  y := comp_ellint_3(0.875, 10);
  f := 0.1581689349611240975;
  testrel(61, NE, y, f, cnt,failed);

  y := comp_ellint_3(2.5, 1.0078125);
  f := -1.538741380875526709;
  testrel(62, NE, y, f, cnt,failed);

  {Mathematic epc[n_, k_] := EllipticPi[n, k^2], Maple V/7 are buggy for m>n>1}
  y := comp_ellint_3(2.5, 2);
  f :=- 0.2552700909026052488;
  testrel(63, NE, y, f, cnt,failed);

  y := comp_ellint_3(2.5, 10);
  f := -0.04588911787903200023;
  testrel(64, NE, y, f, cnt,failed);

  y := comp_ellint_3(4, 1.5);
  f := -0.2013855571558890710;
  testrel(65, NE, y, f, cnt,failed);

  y := comp_ellint_3(4, 2-1/1024);
  f := -0.1355154851796886485;
  testrel(66, NE, y, f, cnt,failed);

  y := comp_ellint_3(4, 2);
  f := -0.1354329621533200822;
  testrel(67, NE, y, f, cnt,failed);

  y := comp_ellint_3(4, 2+1/1024);
  f := -0.1353505502040747876;
  testrel(68, NE, y, f, cnt,failed);  { k^2 > nu and k^2 ~ nu}

  y := comp_ellint_3(4, 1000);
  f := -0.0002430031307463734391;
  testrel(69, NE, y, f, cnt,failed);

  y := comp_ellint_3(2, 2);
  f := -0.3639605746942992989;
  testrel(70, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_heuman_lambda;
var
  k,x,y,f: double;
  lf,i,cnt,failed: integer;
const
  NE = 1;
begin

  {Maple hl := (x,k) -> 2/Pi*(e1(x,sqrt(1-k^2))*(EllipticE(k)-EllipticK(k)) + EllipticK(k)*e2(x,sqrt(1-k^2)));}
  cnt := 0;
  failed := 0;
  writeln('Function: ','heuman_lambda');

  k := 0.99;
  x := 1e-8;
  f := 0.6547480354294812625e-8;
  y := heuman_lambda(x,k);
  testrel(1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.8182013215904385693e-1;
  y := heuman_lambda(x,k);
  testrel(2, NE, y, f, cnt,failed);

  x := 0.5;
  f := 0.3259450823104420013;
  y := heuman_lambda(x,k);
  testrel(3, NE, y, f, cnt,failed);

  x := 1;
  f := 0.6448884938344582100;
  y := heuman_lambda(x,k);
  testrel(4, NE, y, f, cnt,failed);

  x := 1.5625;
  f := 0.9948694761638905011;
  y := heuman_lambda(x,k);
  testrel(5, NE, y, f, cnt,failed);

  x := 2;
  f := 1.266353779072924062;
  y := heuman_lambda(x,k);
  testrel(6, NE, y, f, cnt,failed);

  x := 10;
  f := 6.374483952053941207;
  y := heuman_lambda(x,k);
  testrel(7, NE, y, f, cnt,failed);

  x := 100;
  f := 63.65405229968510604;
  y := heuman_lambda(x,k);
  testrel(8, NE, y, f, cnt,failed);

  {-------------------------------------------}
  k := 0.1;
  x := 1e-8;
  f := 0.9974952928612608871e-8;
  y := heuman_lambda(x,k);
  testrel(9, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.1243624720236512411;
  y := heuman_lambda(x,k);
  testrel(10, NE, y, f, cnt,failed);

  x := 0.5;
  f := 0.4782256132133362774;
  y := heuman_lambda(x,k);
  testrel(11, NE, y, f, cnt,failed);

  x := 1;
  f := 0.8393759819123124977;
  y := heuman_lambda(x,k);
  testrel(12, NE, y, f, cnt,failed);

  x := 1.5625;
  f := 0.999583248170516774;
  y := heuman_lambda(x,k);
  testrel(13, NE, y, f, cnt,failed);

  x := 2;
  f := 1.092953526145575069;
  y := heuman_lambda(x,k);
  testrel(14, NE, y, f, cnt,failed);

  x := 10;
  f := 6.542659927234079251;
  y := heuman_lambda(x,k);
  testrel(15, NE, y, f, cnt,failed);

  x := 100;
  f := 63.49490156438618514;
  y := heuman_lambda(x,k);
  testrel(16, NE, y, f, cnt,failed);

  {-------------------------------------------}
  k := 0.5;
  x := 1e-8;
  f := 0.9342154576676941010e-8;
  y := heuman_lambda(x,k);
  testrel(17, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.1164812830322840199;
  y := heuman_lambda(x,k);
  testrel(18, NE, y, f, cnt,failed);

  x := 0.5;
  f := 0.4484638788008964076;
  y := heuman_lambda(x,k);
  testrel(19, NE, y, f, cnt,failed);

  x := 1;
  f := 0.7927451830080710360;
  y := heuman_lambda(x,k);
  testrel(20, NE, y, f, cnt,failed);

  x := 1.5625;
  f := 0.9978538569524656181;
  y := heuman_lambda(x,k);
  testrel(21, NE, y, f, cnt,failed);

  x := 2;
  f := 1+0.1387579786297624343;
  y := heuman_lambda(x,k);
  testrel(22, NE, y, f, cnt,failed);

  x := 10;
  f := 6.509142167436577482;
  y := heuman_lambda(x,k);
  testrel(23, NE, y, f, cnt,failed);

  x := 100;
  f := 63.52624498968529151;
  y := heuman_lambda(x,k);
  testrel(24, NE, y, f, cnt,failed);

  {-------------------------------------------}
  f := 31.73762514629607121;
  y := heuman_lambda(50,0);
  testrel(25, NE, y, f, cnt,failed);

  f := 1.273239544735162686;
  y := heuman_lambda(2,1);
  testrel(26, NE, y, f, cnt,failed);

  {HMF Table 17.8}
  {lambda(30ø\60ø}
  x := DegToRad(30);
  k := sin(DegToRad(60));
  f := 0.392328;
  y := heuman_lambda(x,k);
  testrele(27, 5e-7, y, f, cnt,failed);

  {lambda(55ø\15ø}
  x := DegToRad(55);
  k := sin(DegToRad(15));
  y := heuman_lambda(x,k);
  f := 0.805703;
  testrele(28, 5e-7, y, f, cnt,failed);

  lf := 0;
  k  := 0.75;
  for i:=-500 to 500 do begin
    x := i*Pi_2;
    y := heuman_lambda(x,k);
    if i=0 then f := y
    else f := 1-y/i;
    if abs(f) > NE*eps_d then inc(lf);
  end;
  inc(cnt);
  if lf<>0 then begin
    inc(failed);
    writeln(' Failures for n*Pi/2 test: ',lf);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_jacobi_zeta;
var
  k,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
  NE2 = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_zeta');

  {-------------------------------------------}
  k := 1;
  y := jacobi_zeta(-1.5,k);
  f := -0.997494986604054431;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_zeta(-1.6,k);
  f := +0.999573603041505164;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_zeta(1.5,k);
  f := +0.997494986604054431;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_zeta(1.6,k);
  f := -0.999573603041505164;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_zeta(100,k);
  f := -0.5063656411097587937;
  testrel(5, NE, y, f, cnt,failed);


  {-------------------------------------------}
  k := 0.99;
  y := jacobi_zeta(1e4,k);
         {1234567890123456789}
  f := +0.2090153270661946387;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_zeta(100,k);
  f := -0.3361430018581131097;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_zeta(10,k);
  f := +0.3580623763102360697;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_zeta(8,k);
  f := -0.2545044311256799924;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_zeta(-5,k);
  f := +0.3903337931352334941;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_zeta(2,k);
  f := -0.4540275136517358294;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_zeta(1.5,k);
  f := +0.1376815936485508684;
  testrel(12, NE2, y, f, cnt,failed);   {!!!!!}

  y := jacobi_zeta(1,k);
  f := +0.4720423387274658268;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_zeta(0.5,k);
  f := +0.3199877677658655313;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_zeta(0.1,k);
  f := +0.6914616499994754968e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := jacobi_zeta(1e-8,k);
  f := +0.6935960052824751092e-8;
  testrel(16, NE, y, f, cnt,failed);

  {-------------------------------------------}
  k := 0.5;
  y := jacobi_zeta(1e4,k);
  f := +0.03791141217191336334;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_zeta(100,k);
  f := -0.05751289792171800898;
  testrel(18, NE, y, f, cnt,failed);

  y := jacobi_zeta(10,k);
  f := +0.06028685902587111446;
  testrel(19, NE, y, f, cnt,failed);

  y := jacobi_zeta(8,k);
  f := -0.01999729146882044197;
  testrel(20, NE, y, f, cnt,failed);

  y := jacobi_zeta(-5,k);
  f := +0.03760932996814525948;
  testrel(21, NE, y, f, cnt,failed);

  y := jacobi_zeta(2,k);
  f := -0.05194253745767273272;
  testrel(22, NE, y, f, cnt,failed);

  y := jacobi_zeta(1.5,k);
  f := +0.00981467090913325006;
  testrel(23, NE, y, f, cnt,failed);

  y := jacobi_zeta(1,k);
  f := +0.06184778256509866925;
  testrel(24, NE, y, f, cnt,failed);

  y := jacobi_zeta(0.5,k);
  f := +0.05531701425512965148;
  testrel(25, NE, y, f, cnt,failed);

  y := jacobi_zeta(0.1,k);
  f := +0.01287118127889838976;
  testrel(26, NE, y, f, cnt,failed);

  y := jacobi_zeta(1e-8,k);
  f := +1.294901969618365286e-9;
  testrel(27, NE, y, f, cnt,failed);


  {-------------------------------------------}
  k := 0.1;
  y := jacobi_zeta(1e4,k);
  f := +0.001457130507386052059;
  testrel(28, NE, y, f, cnt,failed);

  y := jacobi_zeta(100,k);
  f := -0.002187390845340071609;
  testrel(29, NE, y, f, cnt,failed);

  y := jacobi_zeta(10,k);
  f := +0.002286925910628055168;
  testrel(30, NE, y, f, cnt,failed);

  y := jacobi_zeta(8,k);
  f := -0.0007224369491555876416;
  testrel(31, NE, y, f, cnt,failed);

  y := jacobi_zeta(-5,k);
  f := +0.001364910307891728524;
  testrel(32, NE, y, f, cnt,failed);

  y := jacobi_zeta(2,k);
  f := -0.001898320166322237758;
  testrel(33, NE, y, f, cnt,failed);

  y := jacobi_zeta(1.5,k);
  f := +0.0003541274403777301120;
  testrel(34, NE, y, f, cnt,failed);

  y := jacobi_zeta(1,k);
  f := +0.002280147954670897267;
  testrel(35, NE, y, f, cnt,failed);

  y := jacobi_zeta(0.5,k);
  f := +0.002107533546992718746;
  testrel(36, NE, y, f, cnt,failed);

  y := jacobi_zeta(0.1,k);
  f := +0.0004973097001265187990;
  testrel(37, NE, y, f, cnt,failed);

  y := jacobi_zeta(1e-8,k);
  f := +5.006281451646925084e-11;
  testrel(38, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_EllipticKim;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','EllipticKim');

  y := EllipticKim(1e-11);
  f := 1.570796326794896619;
  testrel(1, NE, y, f, cnt,failed);

  y := EllipticKim(5e-10);
  f := 1.570796326794896619;
  testrel(2, NE, y, f, cnt,failed);

  y := EllipticKim(1e-6);
  f := 1.570796326794503920;
  testrel(3, NE, y, f, cnt,failed);

  y := EllipticKim(0.125);
  f := 1.564713754401038368;
  testrel(4, NE, y, f, cnt,failed);

  y := EllipticKim(0.5);
  f := 1.484412473422386453;
  testrel(5, NE, y, f, cnt,failed);

  y := EllipticKim(0.9921875);
  f := 1.313811559022087676;
  testrel(6, NE, y, f, cnt,failed);

  y := EllipticKim(1);
  f := 1.311028777146059905;
  testrel(7, NE, y, f, cnt,failed);

  y := EllipticKim(10);
  f := 0.3682192486091410329;
  testrel(8, NE, y, f, cnt,failed);

  y := EllipticKim(2e4);
  f := 0.5644890953612452489e-3;
  testrel(9, NE, y, f, cnt,failed);

  y := EllipticKim(1e5);
  f := 0.1289921982579263854e-3;
  testrel(10, NE, y, f, cnt,failed);

  y := EllipticKim(2e5);
  f := 0.6796183503285681027e-4;
  testrel(11, NE, y, f, cnt,failed);

  y := EllipticKim(1e10);
  f := 0.2441214529106034746e-8;
  testrel(12, NE, y, f, cnt,failed);

  y := EllipticKim(1e11);
  f := 0.2671473038405439314e-9;
  testrel(13, NE, y, f, cnt,failed);

  y := EllipticKim(MaxDouble);
  f := 0.39560089175562875807e-305;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_EllipticECim;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','EllipticECim');

  y := EllipticECim(5e-10);
  f := 1.570796326794896619;
  testrel(1, NE, y, f, cnt,failed);

  y := EllipticECim(6.6e-10);
  f := 1.570796326794896619;
  testrel(2, NE, y, f, cnt,failed);

  y := EllipticECim(1e-6);
  f := 1.570796326795289318;
  testrel(3, NE, y, f, cnt,failed);

  y := EllipticECim(0.125);
  f := 1.576914389652225438;
  testrel(4, NE, y, f, cnt,failed);

  f := 1.664791805391337856;
  y := EllipticECim(0.5);
  testrel(5, NE, y, f, cnt,failed);

  f := 1.910098894513856009;
  y := EllipticECim(1);
  testrel(6, NE, y, f, cnt,failed);

  f := 10.20926091981457201;
  y := EllipticECim(10);
  testrel(7, NE, y, f, cnt,failed);

  y := EllipticECim(2e4);
  f := 0.2000000029474454776e5;
  testrel(8, NE, y, f, cnt,failed);

  y := EllipticECim(1.5e10);
  f := 0.1500000000000000000e11;
  testrel(9, NE, y, f, cnt,failed);

  y := EllipticECim(2e10);
  f := 0.2000000000000000000e11;
  testrel(10, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_EllipticPiCim;
var
  n,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','EllipticPiCim');

  {Mathematica epc[n_, k_] := EllipticPi[n, -k^2]}
  y := EllipticPiCim(0.0, 0.0);
  f := Pi_2;
  testrel(1, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.0, 1e-5);
  f := 1.570796326755626711;
  testrel(2, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.0, 0.75);
  f := 1.400603042332602023;
  testrel(3, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.0, 2);
  f := 1.009452909989211608;
  testrel(4, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.0, -10);
  f := 0.3682192486091410329;
  testrel(5, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.0, 100);
  f := 0.5991339767278711681e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.0, 1e8);
  f := 1.980697510507225562e-7;
  testrel(7, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.5, 0);
  f := 2.221441469079183124;
  testrel(8, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.5, 1e-5);
  f := 2.221441469014118609;
  testrel(9, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.5, 0.75);
  f := 1.943465368727491449;
  testrel(10, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.5, 2);
  f := 1.336657241354634668;
  testrel(11, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.5, -10);
  f := 0.4457715465792184757;
  testrel(12, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.5, 100);
  f := 0.06776581019515492406;
  testrel(13, NE, y, f, cnt,failed);

  y := EllipticPiCim(0.5, 1e8);
  f := 2.059237326846970343e-7;
  testrel(14, NE, y, f, cnt,failed);

  y := EllipticPiCim(2.5, 0);
  f := 0;
  testrel(15, NE, y, f, cnt,failed);

  y := EllipticPiCim(2.5, 1e-5);
  f := 3.141592653377735734e-11;
  testrel(16, NE, y, f, cnt,failed);

  y := EllipticPiCim(2.5, 0.75);
  f := 0.1279781718389826573;
  testrel(17, NE, y, f, cnt,failed);

  y := EllipticPiCim(2.5, 2);
  f := 0.3220741020959794423;
  testrel(18, NE, y, f, cnt,failed);

  y := EllipticPiCim(2.5, -10);
  f := 0.2327564217137776279;
  testrel(19, NE, y, f, cnt,failed);

  y := EllipticPiCim(2.5, 100);
  f := 0.04658877044959606140;
  testrel(20, NE, y, f, cnt,failed);

  y := EllipticPiCim(2.5, 1e8);
  f := 1.847503220444232799e-7;
  testrel(21, NE, y, f, cnt,failed);

  y := EllipticPiCim(-1, 0);
  f := 1.110720734539591562;
  testrel(22, NE, y, f, cnt,failed);

  y := EllipticPiCim(-1, 1e-5);
  f := 1.110720734516587782;
  testrel(23, NE, y, f, cnt,failed);

  y := EllipticPiCim(-1, 0.75);
  f := 1.009560283379708377;
  testrel(24, NE, y, f, cnt,failed);

  y := EllipticPiCim(-1, 2);
  f := 0.7639884329779299318;
  testrel(25, NE, y, f, cnt,failed);

  y := EllipticPiCim(-1, -10);
  f := 0.3071840408785169673;
  testrel(26, NE, y, f, cnt,failed);

  y := EllipticPiCim(-1, 100);
  f := 0.05368357947845765681;
  testrel(27, NE, y, f, cnt,failed);

  y := EllipticPiCim(-1, 1e8);
  f := 1.918374986493202604e-7;
  testrel(28, NE, y, f, cnt,failed);

  {near n=1}
  n := 1.0 - ldexp(1,-20);
  y := EllipticPiCim(n, 0);
  f := 1608.495438637974138;
  testrel(29, NE, y, f, cnt,failed);

  y := EllipticPiCim(n, 0.5);
  f := 1438.834498866204789;
  testrel(30, NE, y, f, cnt,failed);

  y := EllipticPiCim(n, 1);
  f := 1137.733740534597468;
  testrel(31, NE, y, f, cnt,failed);

  y := EllipticPiCim(n, 100);
  f := 16.13405368559726568;
  testrel(32, NE, y, f, cnt,failed);

  n := 1.0 + ldexp(1,-20);
  y := EllipticPiCim(n, 0);
  f := 0;
  testrel(33, NE, y, f, cnt,failed);

  y := EllipticPiCim(n, 0.5);
  f := 0.1525789056308727865;
  testrel(34, NE, y, f, cnt,failed);

  y := EllipticPiCim(n, 1);
  f := 0.3559791215068463688;
  testrel(35, NE, y, f, cnt,failed);

  y := EllipticPiCim(n, 100);
  f := 0.04991114901692350757;
  testrel(36, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.
