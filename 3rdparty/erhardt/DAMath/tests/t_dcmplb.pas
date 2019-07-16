{Test unit for DAMCmplx, branch points related, WE 2014}
unit T_DCMPLB;

{$i STD.INC}

{$ifdef BIT16}
{$N+,F+}
{$endif}

interface

{$ifdef J_OPT}
{$J+}
{$endif}

const
  chkres_err: double = 1e-19;

var
  fail, numtest: integer;

procedure test_arc_sin_cos_tanh;
procedure test_arc_tan_cotc_sinh;
procedure test_arc_sec_csc;
procedure test_arc_cot_csch;
procedure test_arccosh;
procedure test_arcsech;
procedure test_arccoth;
procedure test_arccothc;
procedure test_misc_branches;
procedure test_branch_series;
procedure test_near_branch_points;


implementation


uses
  DAMath, DAMCmplx;

var
  z,w,u: complex;
  d: double;

{---------------------------------------------------------------------------}
procedure check_res(const u,v: complex);
var
  w: complex;
begin
  if v.re=0 then w.re := u.re-v.re else w.re := 1.0-u.re/v.re;
  if v.im=0 then w.im := u.im-v.im else w.im := 1.0-u.im/v.im;
  {writeln(w.re:30, w.im:30);}
  if (abs(w.re) > chkres_err) or (abs(w.im) > chkres_err) then begin
    writeln(w.re:30, w.im:30);
    inc(fail);
  end;
end;


{---------------------------------------------------------------------------}
procedure test_arc_sin_cos_tanh;
begin
  {arcsin, arccos, arctanh; branch cuts at y=0, |x| > 1,}
  {with f(-x + i0) = f(-x) and f(+x - i0) = f(+x)       }

  d := ldexpd(1,-38);

  write('arcsin ');
  z.re := -2;
  z.im := d;
  carcsin(z,w);
  u.re := -1.570796326792796231;
  u.im := +1.316957896924816709;
  check_res(w,u);


  z.re := -2;
  z.im := 0;
  carcsin(z,w);
  u.re := -1.570796326794896619;
  u.im := +1.316957896924816709;
  check_res(w,u);

  z.re := -2;
  z.im := -d;
  carcsin(z,w);
  u.re := -1.570796326792796231;
  u.im := -1.316957896924816709;
  check_res(w,u);

  z.re := 3;
  z.im := d;
  carcsin(z,w);
  u.re := 1.570796326793610399;
  u.im := 1.762747174039086050;
  check_res(w,u);

  z.re := 3;
  z.im := 0;
  carcsin(z,w);
  u.re :=  1.570796326794896619;
  u.im := -1.762747174039086050;
  check_res(w,u);

  z.re := 3;
  z.im := -d;
  carcsin(z,w);
  u.re :=  1.570796326793610399;
  u.im := -1.762747174039086050;
  check_res(w,u);

  z.re := 1;
  z.im := d;
  carcsin(z,w);
  u.re := 1.570794419446263807;
  u.im := 1.907348632813078241e-6;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carcsin(z,w);
  u.re := 1.570796326794896619;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carcsin(z,w);
  u.re :=  1.570794419446263807;
  u.im := -1.907348632813078241e-6;
  check_res(w,u);

  z.re := -1;
  z.im := d;
  carcsin(z,w);
  u.re := -1.570794419446263807;
  u.im :=  1.907348632813078241e-6;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carcsin(z,w);
  u.re := -1.570796326794896619;
  u.im := 0;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carcsin(z,w);
  u.re := -1.570794419446263807;
  u.im := -1.907348632813078241e-6;
  check_res(w,u);


  write('arccos ');
  z.re := -2;
  z.im := d;
  carccos(z,w);
  u.re :=  3.141592653587692850;
  u.im := -1.316957896924816709;
  check_res(w,u);

  z.re := -2;
  z.im := 0;
  carccos(z,w);
  u.re :=  3.141592653589793238;
  u.im := -1.316957896924816709;
  check_res(w,u);

  z.re := -2;
  z.im := -d;
  carccos(z,w);
  u.re := 3.141592653587692850;
  u.im := 1.316957896924816709;
  check_res(w,u);

  z.re := 3;
  z.im := d;
  carccos(z,w);
  u.re :=  1.286219742153748529e-12;
  u.im := -1.762747174039086050;
  check_res(w,u);

  z.re := 3;
  z.im := 0;
  carccos(z,w);
  u.re := 0;
  u.im := 1.762747174039086050;
  check_res(w,u);

  z.re := 3;
  z.im := -d;
  carccos(z,w);
  u.re := 1.286219742153748529e-12;
  u.im := 1.762747174039086050;
  check_res(w,u);

  z.re := 1;
  z.im := d;
  carccos(z,w);
  u.re :=  1.907348632811921759e-6;
  u.im := -1.907348632813078241e-6;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carccos(z,w);
  u.re := 0;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carccos(z,w);
  u.re := 1.907348632811921759e-6;
  u.im := 1.907348632813078241e-6;
  check_res(w,u);

  z.re := -1;
  z.im := d;
  carccos(z,w);
  u.re :=  3.141590746241160427;
  u.im := -1.907348632813078241e-6;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carccos(z,w);
  u.re := 3.141592653589793238;
  u.im := 0;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carccos(z,w);
  u.re := 3.141590746241160427;
  u.im := 1.907348632813078241e-6;
  check_res(w,u);

  write('arctanh ');
  z.re := -2;
  z.im := d;
  carctanh(z,w);
  u.re := -5.493061443340548457e-1;
  u.im :=  1.570796326793683960;
  check_res(w,u);

  z.re := -2;
  z.im := 0;
  carctanh(z,w);
  u.re := -5.493061443340548457e-1;
  u.im := +1.570796326794896619;
  check_res(w,u);

  z.re := -2;
  z.im := -d;
  carctanh(z,w);
  u.re := -5.493061443340548457e-1;
  u.im := -1.570796326793683960;
  check_res(w,u);

  z.re := 3;
  z.im := d;
  carctanh(z,w);
  u.re := 3.465735902799726547e-1;
  u.im := 1.570796326794441872;
  check_res(w,u);

  z.re := 3;
  z.im := 0;
  carctanh(z,w);
  u.re :=  3.465735902799726547e-1;
  u.im := -1.570796326794896619;
  check_res(w,u);

  z.re := 3;
  z.im := -d;
  carctanh(z,w);
  u.re :=  3.465735902799726547e-1;
  u.im := -1.570796326794441872;
  check_res(w,u);

  {poles at z=1,-1}
  z.re := 1;
  z.im := d;
  carctanh(z,w);
  u.re := 13.51637002091893353;
  u.im := 7.853981633983578043e-1;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carctanh(z,w);
  u.re :=  13.51637002091893353;
  u.im := -7.853981633983578043e-1;
  check_res(w,u);

  z.re := -1;
  z.im := d;
  carctanh(z,w);
  u.re := -13.51637002091893353;
  u.im :=  7.853981633983578043e-1;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carctanh(z,w);
  u.re := -13.51637002091893353;
  u.im := -7.853981633983578043e-1;
  check_res(w,u);
end;


{---------------------------------------------------------------------------}
procedure test_arc_tan_cotc_sinh;
begin
  {arctan, arccotc, arcsinh; branch cuts at x=0, |y| > 1,}
  {with f(+0 + iy) = f(iy) and f(-0 -iy) = f(-iy)        }

  d := ldexpd(1,-38);

  write('arctan ');
  z.re := -d;
  z.im := 2;
  carctan(z,w);
  u.re := -1.570796326793683960;
  u.im :=  5.493061443340548457e-1;
  check_res(w,u);

  z.re := 0;
  z.im := 2;
  carctan(z,w);
  u.re := 1.570796326794896619;
  u.im := 5.493061443340548457e-1;
  check_res(w,u);

  z.re := d;
  z.im := 2;
  carctan(z,w);
  u.re := 1.570796326793683960;
  u.im := 5.493061443340548457e-1;
  check_res(w,u);

  z.re := -d;
  z.im := -3;
  carctan(z,w);
  u.re := -1.570796326794441872;
  u.im := -3.465735902799726547e-1;
  check_res(w,u);

  z.re := 0;
  z.im := -3;
  carctan(z,w);
  u.re := -1.570796326794896619;
  u.im := -3.465735902799726547e-1;
  check_res(w,u);

  z.re := d;
  z.im := -3;
  carctan(z,w);
  u.re :=  1.570796326794441872;
  u.im := -3.465735902799726547e-1;
  check_res(w,u);

  {poles at z=i,-i}
  z.re := -d;
  z.im := 1;
  carctan(z,w);
  u.re := -7.853981633983578043e-1;
  u.im := +13.51637002091893353;
  check_res(w,u);

  z.re := d;
  z.im := 1;
  carctan(z,w);
  u.re := 7.853981633983578043e-1;
  u.im := 13.51637002091893353;
  check_res(w,u);


  z.re := -d;
  z.im := -1;
  carctan(z,w);
  u.re := -7.853981633983578043e-1;
  u.im := -13.51637002091893353;
  check_res(w,u);

  z.re := d;
  z.im := -1;
  carctan(z,w);
  u.re :=  7.853981633983578043e-1;
  u.im := -13.51637002091893353;
  check_res(w,u);

  z.re := d;
  z.im := 0.5;
  carctan(z,w);
  u.re := 4.850638409455617269e-12;
  u.im := 5.493061443340548457e-1;
  check_res(w,u);

  z.re := 1e-10;
  z.im := -1;
  carctan(z,w);
  u.re :=  7.853981634224483096e-1;
  u.im := -11.85949905525020107;
  check_res(w,u);


  write('arccotc ');

  z.re := 1e-30;
  z.im := 1e-30;
  carccotc(z,w);
  u.re :=  1.570796326794896619;
  u.im := -1e-30;
  check_res(w,u);

  z.re := 1e-30;
  z.im := 0;
  carccotc(z,w);
  u.re := 1.570796326794896619;
  u.im := 0;
  check_res(w,u);

  z.re := 0;
  z.im := 1e-30;
  carccotc(z,w);
  u.re := 1.570796326794896619;
  u.im := -1e-30;
  check_res(w,u);


  z.re := -d;
  z.im := 2;
  carccotc(z,w);
  u.re :=  3.141592653588580579;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := 0;
  z.im := 2;
  carccotc(z,w);
  u.re := 0;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := d;
  z.im := 2;
  carccotc(z,w);
  u.re :=  1.212659602363904317e-12;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);


  z.re := -d;
  z.im := -3;
  carccotc(z,w);
  u.re := 3.141592653589338491;
  u.im := 3.465735902799726547e-1;
  check_res(w,u);

  z.re := 0;
  z.im := -3;
  carccotc(z,w);
  u.re := 3.141592653589793238;
  u.im := 3.465735902799726547e-1;
  check_res(w,u);

  z.re := d;
  z.im := -3;
  carccotc(z,w);
  u.re := 4.547473508864641190e-13;
  u.im := 3.465735902799726547e-1;
  check_res(w,u);


  {poles at z=i,-i}
  z.re := -d;
  z.im := 1;
  carccotc(z,w);
  u.re :=  2.356194490193254424;
  u.im := -13.51637002091893353;
  check_res(w,u);

  z.re := d;
  z.im := 1;
  carccotc(z,w);
  u.re :=  7.853981633965388149e-1;
  u.im := -13.51637002091893353;
  check_res(w,u);

  z.re := d;
  z.im := 1.0-d;
  carccotc(z,w);
  u.re :=  1.178097245095262970;
  u.im := -13.34308322577803771;
  check_res(w,u);

  z.re := d;
  z.im := 0.5;
  carccotc(z,w);
  u.re :=  1.570796326790045981;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := -d;
  z.im := -1;
  carccotc(z,w);
  u.re := 2.356194490193254424;
  u.im := 13.51637002091893353;
  check_res(w,u);

  z.re := d;
  z.im := -1;
  carccotc(z,w);
  u.re := 7.853981633965388149e-1;
  u.im := 13.51637002091893353;
  check_res(w,u);

  z.re := 1e-10;
  z.im := -1;
  carccotc(z,w);
  u.re := 7.853981633724483096e-1;
  u.im := 11.85949905525020107;
  check_res(w,u);

  z.re := -d;
  z.im := d;
  carccotc(z,w);
  u.re :=  1.570796326798534598;
  u.im := -3.637978807091712952e-12;
  check_res(w,u);


  write('arcsinh ');
  z.re := -d;
  z.im := 2;
  carcsinh(z,w);
  u.re := -1.316957896924816709;
  u.im :=  1.570796326792796231;
  check_res(w,u);

  z.re := 0;
  z.im := 2;
  carcsinh(z,w);
  u.re := 1.316957896924816709;
  u.im := 1.570796326794896619;
  check_res(w,u);

  z.re := d;
  z.im := 2;
  carcsinh(z,w);
  u.re := 1.316957896924816709;
  u.im := 1.570796326792796231;
  check_res(w,u);


  z.re := -d;
  z.im := -3;
  carcsinh(z,w);
  u.re := -1.762747174039086050;
  u.im := -1.570796326793610399;
  check_res(w,u);

  z.re := 0;
  z.im := -3;
  carcsinh(z,w);
  u.re := -1.762747174039086050;
  u.im := -1.570796326794896619;
  check_res(w,u);

  z.re := d;
  z.im := -3;
  carcsinh(z,w);
  u.re :=  1.762747174039086050;
  u.im := -1.570796326793610399;
  check_res(w,u);

  z.re := -d;
  z.im := 1;
  carcsinh(z,w);
  u.re := -1.907348632813078241e-6;
  u.im :=  1.570794419446263807;
  check_res(w,u);

  z.re := 0;
  z.im := 1;
  carcsinh(z,w);
  u.re := 0;
  u.im := 1.570796326794896619;
  check_res(w,u);

  z.re := d;
  z.im := 1;
  carcsinh(z,w);
  u.re := 1.907348632813078241e-6;
  u.im := 1.570794419446263807;
  check_res(w,u);

  z.re := d;
  z.im := 0.5;
  carcsinh(z,w);
  u.re := 4.200776087161108186e-12;
  u.im := 5.235987755982988731e-1;
  check_res(w,u);

  z.re := -d;
  z.im := -1;
  carcsinh(z,w);
  u.re := -1.907348632813078241e-6;
  u.im := -1.570794419446263807;
  check_res(w,u);

  z.re := 0;
  z.im := -1;
  carcsinh(z,w);
  u.re := 0;
  u.im := -1.570796326794896619;
  check_res(w,u);

  z.re := d;
  z.im := -1;
  carcsinh(z,w);
  u.re :=  1.907348632813078241e-6;
  u.im := -1.570794419446263807;
  check_res(w,u);

  z.re := 1e-10;
  z.im := -1;
  carcsinh(z,w);
  u.re :=  1.000000000008333333e-5;
  u.im := -1.570786326794896703;
  check_res(w,u);

end;


{---------------------------------------------------------------------------}
procedure test_arc_cot_csch;
begin
  {arccot, arccsch; branch cuts at x=0, |y| < 1,}
  {with f(-0 + iy) = f(iy) and f(+0 -iy) = f(-iy), y>0}

  d := ldexpd(1,-38);

  write('arccsch ');
  z.re := -d;
  z.im := 0.5;
  carccsch(z,w);
  u.re := -1.316957896924816709;
  u.im := -1.570796326786495067;
  check_res(w,u);

  z.re := 0;
  z.im := 0.5;
  carccsch(z,w);
  u.re := -1.316957896924816709;
  u.im := -1.570796326794896619;
  check_res(w,u);

  z.re := +d;
  z.im := 0.5;
  carccsch(z,w);
  u.re :=  1.316957896924816709;
  u.im := -1.570796326786495067;
  check_res(w,u);

  z.re := -d;
  z.im := -0.75;
  carccsch(z,w);
  u.re := -7.953654612239056305e-1;
  u.im :=  1.570796326787563143;
  check_res(w,u);

  z.re := 0;
  z.im := -0.75;
  carccsch(z,w);
  u.re := 7.953654612239056305e-1;
  u.im := 1.570796326794896619;
  check_res(w,u);

  z.re := +d;
  z.im := -0.75;
  carccsch(z,w);
  u.re := 7.953654612239056305e-1;
  u.im := 1.570796326787563143;
  check_res(w,u);


  z.re := -d;
  z.im := 1;
  carccsch(z,w);
  u.re := -1.907348632809608794e-6;
  u.im := -1.570794419446263804;
  check_res(w,u);

  z.re := 0;
  z.im := 1;
  carccsch(z,w);
  u.re := 0;
  u.im := -1.570796326794896619;
  check_res(w,u);

  z.re := +d;
  z.im := 1;
  carccsch(z,w);
  u.re :=  1.907348632809608794e-6;
  u.im := -1.570794419446263804;
  check_res(w,u);

  z.re := -d;
  z.im := -1;
  carccsch(z,w);
  u.re := -1.907348632809608794e-6;
  u.im :=  1.570794419446263804;
  check_res(w,u);

  z.re := 0;
  z.im := -1;
  carccsch(z,w);
  u.re := 0;
  u.im := 1.570796326794896619;
  check_res(w,u);

  z.re := +d;
  z.im := -1;
  carccsch(z,w);
  u.re := 1.907348632809608794e-6;
  u.im := 1.570794419446263804;
  check_res(w,u);

  z.re := -d;
  z.im := 2;
  carccsch(z,w);
  u.re := -1.050194021790277046e-12;
  u.im := -5.235987755982988731e-1;
  check_res(w,u);

  z.re := 0;
  z.im := 2;
  carccsch(z,w);
  u.re := 0;
  u.im := -5.235987755982988731e-1;
  check_res(w,u);

  z.re := +d;
  z.im := 2;
  carccsch(z,w);
  u.re :=  1.050194021790277046e-12;
  u.im := -5.235987755982988731e-1;
  check_res(w,u);

  z.re := -d;
  z.im := -3;
  carccsch(z,w);
  u.re := -4.287399140512495096e-13;
  u.im :=  3.398369094541219371e-1;
  check_res(w,u);

  z.re := 0;
  z.im := -3;
  carccsch(z,w);
  u.re := 0;
  u.im := 3.398369094541219371e-1;
  check_res(w,u);

  z.re := +d;
  z.im := -3;
  carccsch(z,w);
  u.re := 4.287399140512495096e-13;
  u.im := 3.398369094541219371e-1;
  check_res(w,u);


  {acot := z -> arctan(1/z)}
  write('arccot ');
  z.re := -d;
  z.im := 0.5;
  carccot(z,w);
  u.re := -1.570796326790045981;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := 0;
  z.im := 0.5;
  carccot(z,w);
  u.re := -1.570796326794896619;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := +d;
  z.im := 0.5;
  carccot(z,w);
  u.re :=  1.570796326790045981;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := -d;
  z.im := -0.25;
  carccot(z,w);
  u.re := -1.570796326791016109;
  u.im :=  2.554128118829953416e-1;
  check_res(w,u);

  z.re := 0;
  z.im := -0.25;
  carccot(z,w);
  u.re := 1.570796326794896619;
  u.im := 2.554128118829953416e-1;
  check_res(w,u);

  z.re := +d;
  z.im := -0.25;
  carccot(z,w);
  u.re := 1.570796326791016109;
  u.im := 2.554128118829953416e-1;
  check_res(w,u);


  z.re := -d;
  z.im := 1;
  carccot(z,w);
  u.re := -7.853981633965388149e-1;
  u.im := -13.51637002091893353;
  check_res(w,u);

  z.re := +d;
  z.im := 1;
  carccot(z,w);
  u.re :=  7.853981633965388149e-1;
  u.im := -13.51637002091893353;
  check_res(w,u);

  z.re := d;
  z.im := 1.0-d;
  carccot(z,w);
  u.re :=  1.178097245095262970;
  u.im := -13.34308322577803771;
  check_res(w,u);

  z.re := -d;
  z.im := -1;
  carccot(z,w);
  u.re := -7.853981633965388149e-1;
  u.im :=  13.51637002091893353;
  check_res(w,u);

  z.re := +d;
  z.im := -1;
  carccot(z,w);
  u.re := 7.853981633965388149e-1;
  u.im := 13.51637002091893353;
  check_res(w,u);

  z.re := -d;
  z.im := 2;
  carccot(z,w);
  u.re := -1.212659602363904317e-12;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := 0;
  z.im := 2;
  carccot(z,w);
  u.re := 0;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := d;
  z.im := 2;
  carccot(z,w);
  u.re :=  1.212659602363904317e-12;
  u.im := -5.493061443340548457e-1;
  check_res(w,u);

  z.re := -d;
  z.im := -3;
  carccot(z,w);
  u.re := -4.547473508864641190e-13;
  u.im :=  3.465735902799726547e-1;
  check_res(w,u);

  z.re := 0;
  z.im := -3;
  carccot(z,w);
  u.re := 0;
  u.im := 3.465735902799726547e-1;
  check_res(w,u);

  z.re := +d;
  z.im := -3;
  carccot(z,w);
  u.re := 4.547473508864641190e-13;
  u.im := 3.465735902799726547e-1;
  check_res(w,u);

  z.re := -d;
  z.im := d;
  carccot(z,w);
  u.re := -1.570796326791258640;
  u.im := -3.637978807091712952e-12;
  check_res(w,u);

end;


{---------------------------------------------------------------------------}
procedure test_arc_sec_csc;
begin
  {arcsec, arccsc; branch cuts at y=0, |x| < 1,}
  {with f(x + i0) = f(x) and f(-x -i0) = f(-x), x>0}

  d := ldexpd(1,-38);

  write('arcsec ');
  z.re := 0.5;
  z.im := -d;
  carcsec(z,w);
  u.re :=  8.401552174322216371e-12;
  u.im := -1.316957896924816709;
  check_res(w,u);

  z.re := 0.5;
  z.im := 0;
  carcsec(z,w);
  u.re := 0;
  u.im := 1.316957896924816709;
  check_res(w,u);

  z.re := 0.5;
  z.im := +d;
  carcsec(z,w);
  u.re := 8.401552174322216371e-12;
  u.im := 1.316957896924816709;
  check_res(w,u);

  z.re := -0.75;
  z.im := -d;
  carcsec(z,w);
  u.re :=  3.141592653582459763;
  u.im := -7.953654612239056305e-1;
  check_res(w,u);

  z.re := -0.75;
  z.im := 0;
  carcsec(z,w);
  u.re :=  3.141592653589793238;
  u.im := -7.953654612239056305e-1;
  check_res(w,u);

  z.re := -0.75;
  z.im := +d;
  carcsec(z,w);
  u.re := 3.141592653582459763;
  u.im := 7.953654612239056305e-1;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carcsec(z,w);
  u.re :=  1.907348632815391206e-6;
  u.im := -1.907348632809608794e-6;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carcsec(z,w);
  u.re :=  2.963588665063205022e-6;
  u.im := -1.227558618359092331e-6;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carcsec(z,w);
  u.re := 0;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carcsec(z,w);
  u.re := 1.907348632815391206e-6;
  u.im := 1.907348632809608794e-6;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carcsec(z,w);
  u.re :=  3.141590746241160423;
  u.im := -1.907348632809608794e-6;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carcsec(z,w);
  u.re := 3.141592653589793238;
  u.im := 0;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carcsec(z,w);
  u.re := 3.141590746241160423;
  u.im := 1.907348632809608794e-6;
  check_res(w,u);

  z.re := -1-d;
  z.im := -d;
  carcsec(z,w);
  u.re :=  3.141589690001128175;
  u.im := -1.227558618359092331e-6;
  check_res(w,u);

  z.re := 1;
  z.im := 1e-12;
  carcsec(z,w);
  u.re := 1.000000000000416667e-6;
  u.im := 0.999999999999583333e-6;
  check_res(w,u);


  write('arccsc ');
  z.re := 0.5;
  z.im := -d;
  carccsc(z,w);
  u.re := 1.570796326786495067;
  u.im := 1.316957896924816709;
  check_res(w,u);

  z.re := 0.5;
  z.im := 0;
  carccsc(z,w);
  u.re :=  1.570796326794896619;
  u.im := -1.316957896924816709;
  check_res(w,u);

  z.re := 0.5;
  z.im := +d;
  carccsc(z,w);
  u.re :=  1.570796326786495067;
  u.im := -1.316957896924816709;
  check_res(w,u);

  z.re := -0.75;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.570796326787563143;
  u.im :=  7.953654612239056305e-1;

  z.re := -0.75;
  z.im := 0;
  carccsc(z,w);
  u.re := -1.570796326794896619;
  u.im :=  7.953654612239056305e-1;
  check_res(w,u);

  z.re := -0.75;
  z.im := +d;
  carccsc(z,w);
  u.re := -1.570796326787563143;
  u.im := -7.953654612239056305e-1;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carccsc(z,w);
  u.re := 1.570794419446263804;
  u.im := 1.907348632809608794e-6;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carccsc(z,w);
  u.re := 1.570793363206231556;
  u.im := 1.227558618359092331e-6;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carccsc(z,w);
  u.re := 1.570796326794896619;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carccsc(z,w);
  u.re :=  1.570794419446263804;
  u.im := -1.907348632809608794e-6;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.570794419446263804;
  u.im :=  1.907348632809608794e-6;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carccsc(z,w);
  u.re := -1.570796326794896619;
  u.im := 0;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carccsc(z,w);
  u.re := -1.570794419446263804;
  u.im := -1.907348632809608794e-6;
  check_res(w,u);

  z.re := -1-d;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.570793363206231556;
  u.im :=  1.227558618359092331e-6;
  check_res(w,u);

  z.re := -1+d;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.570795099236278247;
  u.im :=  2.963588665068468056e-6;
  check_res(w,u);
end;


{---------------------------------------------------------------------------}
procedure test_arccosh;
begin
  {arccosh; branch cut at y=0, x < 1,}
  {with f(x + i0) = f(x)}

  d := ldexpd(1,-38);
  write('arccosh ');

  z.re := -2;
  z.im := -d;
  carccosh(z,w);
  u.re :=  1.316957896924816709;
  u.im := -3.141592653587692850;
  check_res(w,u);

  z.re := -2;
  z.im := 0;
  carccosh(z,w);
  u.re := 1.316957896924816709;
  u.im := 3.141592653589793238;
  check_res(w,u);

  z.re := -2;
  z.im := +d;
  carccosh(z,w);
  u.re := 1.316957896924816709;
  u.im := 3.141592653587692850;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carccosh(z,w);
  u.re :=  1.907348632813078241e-6;
  u.im := -3.141590746241160427;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carccosh(z,w);
  u.re := 0;
  u.im := 3.141592653589793238;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carccosh(z,w);
  u.re := 1.907348632813078241e-6;
  u.im := 3.141590746241160427;
  check_res(w,u);

  z.re := 0.5;
  z.im := -d;
  carccosh(z,w);
  u.re :=  4.200776087161108186e-12;
  u.im := -1.047197551196597746;
  check_res(w,u);

  z.re := 0.5;
  z.im := 0;
  carccosh(z,w);
  u.re := 0;
  u.im := 1.047197551196597746;
  check_res(w,u);

  z.re := 0.5;
  z.im := +d;
  carccosh(z,w);
  u.re := 4.200776087161108186e-12;
  u.im := 1.047197551196597746;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carccosh(z,w);
  u.re :=  1.907348632813078241e-6;
  u.im := -1.907348632811921759e-6;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carccosh(z,w);
  u.re := 0;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carccosh(z,w);
  u.re := 1.907348632813078241e-6;
  u.im := 1.907348632811921759e-6;
  check_res(w,u);


  z.re := 2.5;
  z.im := -d;
  carccosh(z,w);
  u.re :=  1.566799236972411079;
  u.im := -1.587744120013611837e-12;
  check_res(w,u);

  z.re := 2.5;
  z.im := 0;
  carccosh(z,w);
  u.re := 1.566799236972411079;
  u.im := 0;
  check_res(w,u);

  z.re := 2.5;
  z.im := +d;
  carccosh(z,w);
  u.re := 1.566799236972411079;
  u.im := 1.587744120013611837e-12;
  check_res(w,u);
end;


{---------------------------------------------------------------------------}
procedure test_arcsech;
begin
  {arccosh; branch cut at y=0, x < 0 or x > 1}
  {with f(x - i0) = f(x)}

  d := ldexpd(1,-38);
  write('arcsech ');

  z.re := -2;
  z.im := -d;
  carcsech(z,w);
  u.re := 1.050194021790277046e-12;
  u.im := 2.094395102393195492;
  check_res(w,u);

  z.re := -2;
  z.im := 0;
  carcsech(z,w);
  u.re := 0;
  u.im := 2.094395102393195492;
  check_res(w,u);

  z.re := -2;
  z.im := +d;
  carcsech(z,w);
  u.re :=  1.050194021790277046e-12;
  u.im := -2.094395102393195492;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carcsech(z,w);
  u.re := 1.907348632809608794e-6;
  u.im := 3.141590746241160423;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carcsech(z,w);
  u.re := 0;
  u.im := 3.141592653589793238;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carcsech(z,w);
  u.re :=  1.907348632809608794e-6;
  u.im := -3.141590746241160423;
  check_res(w,u);

  z.re := -1+d;
  z.im := +d;
  carcsech(z,w);
  u.re :=  2.963588665068468056e-6;
  u.im := -3.141591426031174867;
  check_res(w,u);

  z.re := -1+d;
  z.im := -d;
  carcsech(z,w);
  u.re := 2.963588665068468056e-6;
  u.im := 3.141591426031174867;
  check_res(w,u);

  z.re := 0.5;
  z.im := -d;
  carcsech(z,w);
  u.re := 1.316957896924816709;
  u.im := 8.401552174322216371e-12;
  check_res(w,u);

  z.re := 0.5;
  z.im := 0;
  carcsech(z,w);
  u.re := 1.316957896924816709;
  u.im := 0;
  check_res(w,u);

  z.re := 0.5;
  z.im := +d;
  carcsech(z,w);
  u.re :=  1.316957896924816709;
  u.im := -8.401552174322216371e-12;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carcsech(z,w);
  u.re := 1.907348632809608794e-6;
  u.im := 1.907348632815391206e-6;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carcsech(z,w);
  u.re := 0;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carcsech(z,w);
  u.re :=  1.907348632809608794e-6;
  u.im := -1.907348632815391206e-6;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carcsech(z,w);
  u.re := 1.227558618359092331e-6;
  u.im := 2.963588665063205022e-6;
  check_res(w,u);

  z.re := 2.5;
  z.im := -d;
  carcsech(z,w);
  u.re := 6.350976480054447348e-13;
  u.im := 1.159279480727408600;
  check_res(w,u);

  z.re := 2.5;
  z.im := 0;
  carcsech(z,w);
  u.re := 0;
  u.im := 1.159279480727408600;
  check_res(w,u);

  z.re := 2.5;
  z.im := +d;
  carcsech(z,w);
  u.re :=  6.350976480054447349e-13;
  u.im := -1.159279480727408600;
  check_res(w,u);

end;


{---------------------------------------------------------------------------}
procedure test_arccothc;
begin
  {arccothc; branch cut at y=0, |x| > 1,}
  {with f(-x + i0) = f(-x) and f(+x - i0) = f(+x)}
  d := ldexpd(1,-38);

  write('arccothc ');
  z.re := -2;
  z.im := -d;
  carccothc(z,w);
  u.re := -5.493061443340548457e-1;
  u.im :=  1.212659602363904317e-12;
  check_res(w,u);

  z.re := -2;
  z.im := 0;
  carccothc(z,w);
  u.re := -5.493061443340548457e-1;
  u.im :=  3.141592653589793238;
  check_res(w,u);

  z.re := -2;
  z.im := d;
  carccothc(z,w);
  u.re := -5.493061443340548457e-1;
  u.im :=  3.141592653588580579;
  check_res(w,u);

  z.re := 3;
  z.im := -d;
  carccothc(z,w);
  u.re := 3.465735902799726547e-1;
  u.im := 4.547473508864641190e-13;
  check_res(w,u);

  z.re := 3;
  z.im := 0;
  carccothc(z,w);
  u.re := 3.465735902799726547e-1;
  u.im := 0;
  check_res(w,u);

  z.re := 3;
  z.im := d;
  carccothc(z,w);
  u.re := 3.465735902799726547e-1;
  u.im := 3.141592653589338491;
  check_res(w,u);

  {poles at z=1,-1}
  z.re := 1;
  z.im := -d;
  carccothc(z,w);
  u.re := 13.51637002091893353;
  u.im := 7.853981633965388149e-1;
  check_res(w,u);

  z.re := 1;
  z.im := d;
  carccothc(z,w);
  u.re := 13.51637002091893353;
  u.im := 2.356194490193254424;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carccothc(z,w);
  u.re := -13.51637002091893353;
  u.im :=  7.853981633965388149e-1;
  check_res(w,u);

  z.re := -1;
  z.im := d;
  carccothc(z,w);
  u.re := -13.51637002091893353;
  u.im :=  2.356194490193254424;
  check_res(w,u);

  z.re := 0.5;
  z.im := -d;
  carccothc(z,w);
  u.re := 5.493061443340548457e-1;
  u.im := 1.570796326790045981;
  check_res(w,u);

  z.re := 0.5;
  z.im := 0;
  carccothc(z,w);
  u.re := 5.493061443340548457e-1;
  u.im := 1.570796326794896619;
  check_res(w,u);

  z.re := 0.5;
  z.im := +d;
  carccothc(z,w);
         {1.234567890123456789}
  u.re := 5.493061443340548457e-1;
  u.im := 1.570796326799747258;
  check_res(w,u);

end;


{---------------------------------------------------------------------------}
procedure test_arccoth;
begin
  {arccoth; branch cut at y=0, |x| < 1,}
  {with f(-x - i0) = f(-x) and f(+x + i0) = f(+x)}
  d := ldexpd(1,-38);

  write('arccoth ');
  z.re := 0.5;
  z.im := -d;
  carccoth(z,w);
  u.re := 5.493061443340548457e-1;
  u.im := 1.570796326790045981;
  check_res(w,u);

  z.re := 0.5;
  z.im := 0;
  carccoth(z,w);
  u.re :=  5.493061443340548457e-1;
  u.im := -1.570796326794896619;
  check_res(w,u);

  z.re := 0.5;
  z.im := +d;
  carccoth(z,w);
  u.re :=  5.493061443340548457e-1;
  u.im := -1.570796326790045981;
  check_res(w,u);

  z.re := -0.25;
  z.im := -d;
  carccoth(z,w);
  u.re := -2.554128118829953416e-1;
  u.im :=  1.570796326791016109;
  check_res(w,u);

  z.re := -0.25;
  z.im := 0;
  carccoth(z,w);
  u.re := -2.554128118829953416e-1;
  u.im :=  1.570796326794896619;
  check_res(w,u);

  z.re := -0.25;
  z.im := +d;
  carccoth(z,w);
  u.re := -2.554128118829953416e-1;
  u.im := -1.570796326791016109;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carccoth(z,w);
  u.re := 13.51637002091893353;
  u.im := 7.853981633965388149e-1;
  check_res(w,u);

  {poles at z=1, z=-1}

  z.re := 1;
  z.im := +d;
  carccoth(z,w);
  u.re :=  13.51637002091893353;
  u.im := -7.853981633965388149e-1;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carccoth(z,w);
  u.re := -13.51637002091893353;
  u.im :=  7.853981633965388149e-1;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carccoth(z,w);
  u.re := -13.51637002091893353;
  u.im := -7.853981633965388149e-1;
  check_res(w,u);

  z.re := 3;
  z.im := -d;
  carccoth(z,w);
  u.re := 3.465735902799726547e-1;
  u.im := 4.547473508864641190e-13;
  check_res(w,u);

  z.re := 3;
  z.im := 0;
  carccoth(z,w);
  u.re := 3.465735902799726547e-1;
  u.im := 0;
  check_res(w,u);

  z.re := 3;
  z.im := +d;
  carccoth(z,w);
  u.re :=  3.465735902799726547e-1;
  u.im := -4.547473508864641190e-13;
  check_res(w,u);
end;


{---------------------------------------------------------------------------}
procedure test_misc_branches;
begin
  d := ldexpd(1,-38);
  write('ln ');
  z.re := -2;
  z.im := 0;
  cln(z,w);
  u.re := 0.69314718055994530942;
  u.im := 3.1415926535897932385;
  check_res(w,u);
  z.re := -2;
  z.im := d;
  cln(z,w);
  u.re := 0.69314718055994530942;
  u.im := 3.1415926535879742491;
  check_res(w,u);
  z.re := -2;
  z.im := -d;
  cln(z,w);
  u.re := 0.69314718055994530942;
  u.im := -3.1415926535879742491;
  check_res(w,u);

  write('ln1p ');
  z.re := -4;
  z.im := 0;
  cln1p(z,w);
  u.re := 1.0986122886681096914;
  u.im := 3.1415926535897932385;
  check_res(w,u);
  z.re := -4;
  z.im := d;
  cln1p(z,w);
  u.re := 1.0986122886681096914;
  u.im := 3.1415926535885805789;
  check_res(w,u);
  z.re := -4;
  z.im := -d;
  cln1p(z,w);
  u.re := 1.0986122886681096914;
  u.im := -3.1415926535885805789;
  check_res(w,u);

  write('E1 ');
  z.re := -2;
  z.im := 0;
  ce1(z,w);
  u.re := -4.9542343560018901634;
  u.im := -3.1415926535897932385;
  check_res(w,u);
  z.re := -2;
  z.im := d;
  ce1(z,w);
  u.re := -4.9542343560018901634;
  u.im := -3.1415926535763526237;
  check_res(w,u);
  z.re := -2;
  z.im := -d;
  ce1(z,w);
  u.re := -4.9542343560018901634;
  u.im :=  3.1415926535763526237;
  check_res(w,u);

  write('Ei ');
  z.re := -2;
  z.im := 0;
  cei(z,w);
  u.re := -0.048900510708061119567;
  u.im := 0;
  check_res(w,u);
  z.re := -2;
  z.im := d;
  cei(z,w);
  u.re := -0.048900510708061119567;
  u.im := +3.1415926535895470650;
  check_res(w,u);
  z.re := -2;
  z.im := -d;
  cei(z,w);
  u.re := -0.048900510708061119567;
  u.im := -3.1415926535895470650;
  check_res(w,u);
end;



{---------------------------------------------------------------------------}
procedure test_branch_series;
begin
  d := ldexpd(1,-21);

  write('arcsec ');

  z.re := 1;
  z.im := -d;
  carcsec(z,w);
  u.re :=  6.905341031992181169e-4;
  u.im := -6.905338288056731242e-4;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carcsec(z,w);
  u.re :=  1.072933579684633726e-3;
  u.im := -4.444233904836025431e-4;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carcsec(z,w);
  u.re := 0;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carcsec(z,w);
  u.re := 6.905341031992181169e-4;
  u.im := 6.905338288056731242e-4;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carcsec(z,w);
  u.re :=  3.140902119486594020;
  u.im := -6.905338288056731242e-4;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carcsec(z,w);
  u.re := 3.141592653589793238;
  u.im := 0;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carcsec(z,w);
  u.re := 3.140902119486594020;
  u.im := 6.905338288056731242e-4;
  check_res(w,u);

  z.re := -1-d;
  z.im := -d;
  carcsec(z,w);
  u.re :=  3.140519720010108605;
  u.im := -4.444233904836025431e-4;
  check_res(w,u);


  write('arccsc ');
  z.re := 1;
  z.im := -d;
  carccsc(z,w);
  u.re := 1.570105792691697401;
  u.im := 6.905338288056731242e-4;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carccsc(z,w);
  u.re := 1.569723393215211986;
  u.im := 4.444233904836025431e-4;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carccsc(z,w);
  u.re := 1.570796326794896619;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carccsc(z,w);
  u.re :=  1.570105792691697401;
  u.im := -6.905338288056731242e-4;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.570105792691697401;
  u.im :=  6.905338288056731242e-4;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carccsc(z,w);
  u.re := -1.570796326794896619;
  u.im := 0;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carccsc(z,w);
  u.re := -1.570105792691697401;
  u.im := -6.905338288056731242e-4;
  check_res(w,u);

  z.re := -1-d;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.569723393215211986;
  u.im :=  4.444233904836025431e-4;
  check_res(w,u);

  z.re := -1+d;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.570351902801469359;
  u.im :=  1.072933829432073985e-3;
  check_res(w,u);

  write('arcsech ');
  z.re := -1;
  z.im := -d;
  carcsech(z,w);
  u.re := 6.905338288056731242e-4;
  u.im := 3.140902119486594020;
  check_res(w,u);

  z.re := -1;
  z.im := 0;
  carcsech(z,w);
  u.re := 0;
  u.im := 3.141592653589793238;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carcsech(z,w);
  u.re :=  6.905338288056731242e-4;
  u.im := -3.140902119486594020;
  check_res(w,u);

  z.re := -1+d;
  z.im := +d;
  carcsech(z,w);
  u.re :=  1.072933829432073985e-3;
  u.im := -3.141148229596365978;
  check_res(w,u);

  z.re := -1+d;
  z.im := -d;
  carcsech(z,w);
  u.re := 1.072933829432073985e-3;
  u.im := 3.141148229596365978;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carcsech(z,w);
  u.re := 6.905338288056731242e-4;
  u.im := 6.905341031992181169e-4;
  check_res(w,u);

  z.re := 1;
  z.im := 0;
  carcsech(z,w);
  u.re := 0;
  u.im := 0;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carcsech(z,w);
  u.re :=  6.905338288056731242e-4;
  u.im := -6.905341031992181169e-4;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carcsech(z,w);
  u.re := 4.444233904836025431e-4;
  u.im := 1.072933579684633726e-3;
  check_res(w,u);
end;


{---------------------------------------------------------------------------}
procedure test_near_branch_points;
begin

  write('arcsec ');

  d := ldexpd(1,-10);
  z.re := 1+d;
  z.im := -d;
  carcsec(z,w);
  u.re :=  4.854385287045888337e-2;
  u.im := -2.008440430044694158e-2;
  check_res(w,u);

  z.re := 1+d;
  z.im := -1/16;
  carcsec(z,w);
  u.re :=  2.580346175413484997e-1;
  u.im := -2.411532936877794409e-1;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carcsec(z,w);
  u.re := 3.126270764243306024e-2;
  u.im := 3.123727633882749806e-2;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carcsec(z,w);
  u.re :=  3.110329945947360178;
  u.im := -3.123727633882749806e-2;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carcsec(z,w);
  u.re := 3.110329945947360178;
  u.im := 3.123727633882749806e-2;
  check_res(w,u);

  z.re := -1-d;
  z.im := -0.5*d;
  carcsec(z,w);
  u.re :=  3.096129384918737878;
  u.im := -1.072266569112678745e-2;
  check_res(w,u);


  write('arccsc ');
  z.re := 1;
  z.im := -0.5*d;
  carccsc(z,w);
  u.re := 1.548694745635354515;
  u.im := 2.209258983287612288e-2;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carccsc(z,w);
  u.re := 1.522252473924437736;
  u.im := 2.008440430044694158e-2;
  check_res(w,u);

  z.re := 1;
  z.im := +d;
  carccsc(z,w);
  u.re :=  1.539533619152463559;
  u.im := -3.123727633882749806e-2;
  check_res(w,u);

  z.re := -1;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.539533619152463559;
  u.im :=  3.123727633882749806e-2;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carccsc(z,w);
  u.re := -1.539533619152463559;
  u.im := -3.123727633882749806e-2;
  check_res(w,u);

  z.re := -1-d;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.522252473924437736;
  u.im :=  2.008440430044694158e-2;
  check_res(w,u);

  z.re := -1+d;
  z.im := -d;
  carccsc(z,w);
  u.re := -1.550656040509739050;
  u.im :=  4.856699988721482855e-2;
  check_res(w,u);

  write('arcsech ');
  z.re := -1;
  z.im := -0.5*d;
  carcsech(z,w);
  u.re := 2.209258983287612288e-2;
  u.im := 3.119491072430251135;
  check_res(w,u);

  z.re := -1;
  z.im := +d;
  carcsech(z,w);
  u.re :=  3.123727633882749806e-2;
  u.im := -3.110329945947360178;
  check_res(w,u);

  z.re := -1+d;
  z.im := +d;
  carcsech(z,w);
  u.re :=  4.856699988721482855e-2;
  u.im := -3.121452367304635669;
  check_res(w,u);

  z.re := -1+d;
  z.im := -d;
  carcsech(z,w);
  u.re := 4.856699988721482855e-2;
  u.im := 3.121452367304635669;
  check_res(w,u);

  z.re := 1;
  z.im := -d;
  carcsech(z,w);
  u.re := 3.123727633882749806e-2;
  u.im := 3.126270764243306024e-2;
  check_res(w,u);

  z.re := 1;
  z.im := 2.0*d;
  carcsech(z,w);
  u.re :=  4.415816327071361598e-2;
  u.im := -4.423009376208059029e-2;
  check_res(w,u);

  z.re := 1+d;
  z.im := -d;
  carcsech(z,w);
  u.re := 2.008440430044694158e-2;
  u.im := 4.854385287045888337e-2;
  check_res(w,u);
end;

end.