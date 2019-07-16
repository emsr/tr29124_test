{Part 3 of regression test for SPECFUN unit  (c) 2010  W.Ehrhardt}

unit t_sfd3;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_normstd_pdf;
procedure test_normstd_cdf;
procedure test_normstd_inv;

procedure test_normal_pdf;
procedure test_normal_cdf;
procedure test_normal_inv;

procedure test_beta_cdf;
procedure test_beta_inv;
procedure test_beta_pdf;

procedure test_t_cdf;
procedure test_t_inv;
procedure test_t_pdf;

procedure test_f_pdf;
procedure test_f_cdf;
procedure test_f_inv;


implementation

uses
  amath, specfun, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_normstd_pdf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','normstd_pdf');
  {Maple: statevalf[pdf,normald](x)}

  y := normstd_pdf(NegInf_x);
  f := 0;
  testabs( 1, 0, y, f, cnt,failed);

  y := normstd_pdf(-37.5);
  f := 0.17282337322841052208e-305;
  testrel( 2, NE, y, f, cnt,failed);

  y := normstd_pdf(-6.0);
  f := 0.60758828498232854869e-8;
  testrel( 3, NE, y, f, cnt,failed);

  y := normstd_pdf(-2.0);
  f := 0.53990966513188051949e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := normstd_pdf(-1.1);
  f := 0.21785217703255053138;
  testrel( 5, NE, y, f, cnt,failed);

  y := normstd_pdf(-1.0);
  f := 0.24197072451914334980;
  testrel( 6, NE, y, f, cnt,failed);

  y := normstd_pdf(-0.9);
  f := 0.26608524989875482182;
  testrel( 7, NE, y, f, cnt,failed);

  y := normstd_pdf(-0.5);
  f := 0.35206532676429947777;
  testrel( 8, NE, y, f, cnt,failed);

  y := normstd_pdf(-0.125);
  f := 0.39583768694474948414;
  testrel( 9, NE, y, f, cnt,failed);

  y := normstd_pdf(-0.0078125);
  f := 0.39893010583499324766;
  testrel(10, NE, y, f, cnt,failed);

  y := normstd_pdf(0.0);
  f := 0.39894228040143267794;
  testrel(11, NE, y, f, cnt,failed);

  y := normstd_pdf(0.0078125);
  f := 0.39893010583499324766;
  testrel(12, NE, y, f, cnt,failed);

  y := normstd_pdf(0.125);
  f := 0.39583768694474948414;
  testrel(13, NE, y, f, cnt,failed);

  y := normstd_pdf(0.5);
  f := 0.35206532676429947777;
  testrel(14, NE, y, f, cnt,failed);

  y := normstd_pdf(0.9);
  f := 0.26608524989875482182;
  testrel(15, NE, y, f, cnt,failed);

  y := normstd_pdf(1.0);
  f := 0.24197072451914334980;
  testrel(16, NE, y, f, cnt,failed);

  y := normstd_pdf(1.1);
  f := 0.21785217703255053138;
  testrel(17, NE, y, f, cnt,failed);

  y := normstd_pdf(2.0);
  f := 0.53990966513188051949e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := normstd_pdf(4.0);
  f := 0.13383022576488535177e-3;
  testrel(19, NE, y, f, cnt,failed);

  y := normstd_pdf(6.0);
  f := 0.60758828498232854869e-8;
  testrel(20, NE, y, f, cnt,failed);

  y := normstd_pdf(10.0);
  f := 0.76945986267064193463e-22;
  testrel(21, NE, y, f, cnt,failed);

  y := normstd_pdf(37.5);
  f := 0.17282337322841052208e-305;
  testrel(22, NE, y, f, cnt,failed);

  y := normstd_pdf(38.5);
  f := 0.54251551813365901832e-322;
  testrel(23, NE, y, f, cnt,failed);

  y := normstd_pdf(40);
  f := 0.14632702508383031788e-347;
  testrel(24, NE, y, f, cnt,failed);

  y := normstd_pdf(PosInf_x);
  f := 0;
  testabs(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_normstd_cdf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','normstd_cdf');
  {Maple: statevalf[cdf,normald](x)}

  y := normstd_cdf(NegInf_x);
  f := 0;
  testabs( 1, 0, y, f, cnt,failed);

  y := normstd_cdf(-35);
  f := 0.1124910706472406243978835e-267;
  testrel( 2, NE, y, f, cnt,failed);

  y := normstd_cdf(-30);
  f := 0.4906713927148187059531905e-197;
  testrel( 3, NE, y, f, cnt,failed);

  y := normstd_cdf(-20);
  f := 0.27536241186062336951e-88;
  testrel( 4, NE, y, f, cnt,failed);

  y := normstd_cdf(-10);
  f := 0.76198530241605260660e-23;
  testrel( 5, NE, y, f, cnt,failed);

  y := normstd_cdf(-6.0);
  f := 0.98658764503769814070e-9;
  testrel( 6, NE, y, f, cnt,failed);

  y := normstd_cdf(-2.0);
  f := 0.227501319481792072003e-1;
  testrel( 7, NE, y, f, cnt,failed);

  y := normstd_cdf(-1.1);
  f := 0.13566606094638267518;
  testrel( 8, NE, y, f, cnt,failed);

  y := normstd_cdf(-1.0);
  f := 0.15865525393145705142;
  testrel( 9, NE, y, f, cnt,failed);

  y := normstd_cdf(-0.9);
  f := 0.18406012534675948856;
  testrel(10, NE, y, f, cnt,failed);

  y := normstd_cdf(-0.5);
  f := 0.30853753872598689637;
  testrel(11, NE, y, f, cnt,failed);

  y := normstd_cdf(-0.125);
  f := 0.45026177516988710702;
  testrel(12, NE, y, f, cnt,failed);

  y := normstd_cdf(-0.0078125);
  f := 0.49688329513915741955;
  testrel(13, NE, y, f, cnt,failed);

  y := normstd_cdf(0.0);
  f := 0.5;
  testrel(14, NE, y, f, cnt,failed);

  y := normstd_cdf(0.0078125);
  f := 0.50311670486084258045;
  testrel(15, NE, y, f, cnt,failed);

  y := normstd_cdf(0.125);
  f := 0.54973822483011289298;
  testrel(16, NE, y, f, cnt,failed);

  y := normstd_cdf(0.5);
  f := 0.69146246127401310364;
  testrel(17, NE, y, f, cnt,failed);

  y := normstd_cdf(0.9);
  f := 0.81593987465324051145;
  testrel(18, NE, y, f, cnt,failed);

  y := normstd_cdf(1.0);
  f := 0.84134474606854294859;
  testrel(19, NE, y, f, cnt,failed);

  y := normstd_cdf(1.1);
  f := 0.86433393905361732483;
  testrel(20, NE, y, f, cnt,failed);

  y := normstd_cdf(2.0);
  f := 0.97724986805182079280;
  testrel(21, NE, y, f, cnt,failed);

  y := normstd_cdf(4.0);
  f := 0.99996832875816688008;
  testrel(22, NE, y, f, cnt,failed);

  y := normstd_cdf(6.0);
  f := 0.99999999901341235496;
  testrel(23, NE, y, f, cnt,failed);

  y := normstd_cdf(10.0);
  f := 1;
  testrel(24, NE, y, f, cnt,failed);

  y := normstd_cdf(PosInf_x);
  f := 1;
  testabs(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_normstd_inv;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','normstd_inv');
  {Maple: statevalf[icdf,normald](x)}

  y := normstd_inv(1e-305);
  f := -37.35634609306710368891;
  testrel( 1, NE, y, f, cnt,failed);

  y := normstd_inv(1e-100);
  f := -21.2734535609653242951;
  testrel( 2, NE, y, f, cnt,failed);

  y := normstd_inv(ldexp(1,-24));
  f := -5.2947040848545980574;
  testrel( 3, NE, y, f, cnt,failed);

  y := normstd_inv(0.0009765625);
  f := -3.0972690781987844624;
  testrel( 4, NE, y, f, cnt,failed);

  y := normstd_inv(0.125);
  f := -1.1503493803760081783;
  testrel( 5, NE, y, f, cnt,failed);

  y := normstd_inv(0.250);
  f := -0.67448975019608174320;
  testrel( 6, NE, y, f, cnt,failed);

  y := normstd_inv(0.3125);
  f := -0.48877641111466949891;
  testrel( 7, NE, y, f, cnt,failed);

  y := normstd_inv(0.4375);
  f := -0.15731068461017069552;
  testrel( 8, NE, y, f, cnt,failed);

  y := normstd_inv(0.4990234375);
  f := -0.24478816191106774544e-2;
  testrel( 9, NE, y, f, cnt,failed);

  y := normstd_inv(0.5);
  f := 0;
  testrel(10, NE, y, f, cnt,failed);

  y := normstd_inv(0.5009765625);
  f := 0.24478816191106774544e-2;
  testrel( 11, NE, y, f, cnt,failed);

  y := normstd_inv(0.625);
  f := 0.31863936396437516302;
  testrel(12, NE, y, f, cnt,failed);

  y := normstd_inv(0.750);
  f := 0.67448975019608174320;
  testrel(13, NE, y, f, cnt,failed);

  y := normstd_inv(0.9375);
  f := 1.5341205443525463117;
  testrel(14, NE, y, f, cnt,failed);

  y := normstd_inv(0.9990234375);
  f := 3.0972690781987844624;
  testrel(15, NE, y, f, cnt,failed);

  y := normstd_inv(1-ldexp(1,-20));
  f := 4.7630010342678139570;
  testrel(16, NE, y, f, cnt,failed);

  y := normstd_inv(1-ldexp(1,-50));
  f := 7.9560381254815309622;
  testrel(17, NE, y, f, cnt,failed);

  y := normstd_cdf(-10);
  y := normstd_inv(y);
  f := -10;
  testrel(18, NE, y, f, cnt,failed);

  y := normstd_cdf(2);
  y := normstd_inv(y);
  f := 2;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_normal_pdf;
var
  y,f,mu,sd: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','normal_pdf');
  {Maple: statevalf[pdf,normal[mu,sd]](x)}

  mu := 5;
  sd := 2.0;
  y := normal_pdf(mu, sd, -70);
  f := 0.8641168661420526104e-306;
  testrel(1, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, -1);
  f := 0.2215924205969003588e-2;
  testrel(2, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 0);
  f := 0.8764150246784268681e-2;
  testrel(3, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 2);
  f := 0.6475879783294586381e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 5);
  f := 0.1994711402007163390;
  testrel(5, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 10);
  f := 0.8764150246784268681e-2;
  testrel(6, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 80);
  f := 0.8641168661420526105e-306;
  testrel(7, NE, y, f, cnt,failed);


  mu := -1.5;
  sd := 0.5;
  y := normal_pdf(mu, sd, -20);
  f := 0.4240013103049211254e-297;
  testrel(8, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, -1);
  f := 0.4839414490382866996;
  testrel(9, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 0);
  f := 0.8863696823876014351e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 2);
  f := 0.1826944081672918669e-10;
  testrel(11, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 5);
  f := 0.1599765551401362300e-36;
  testrel(12, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 10);
  f := 0.1074112073004118265e-114;
  testrel(13, NE, y, f, cnt,failed);

  y := normal_pdf(mu, sd, 17);
  f := 0.4240013103049211254e-297;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_normal_cdf;
var
  y,f,mu,sd: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','normal_cdf');

  {statevalf[cdf,normald[mu,sd]](x);}
  mu := 5;
  sd := 2.0;

  y := normal_cdf(mu, sd, -70);
  f := 0.4605353009581954842e-307;
  testrel(1, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, -1);
  f := 0.1349898031630094527e-2;
  testrel(2, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 0);
  f := 0.62096653257761351670e-2;
  testrel(3, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 2);
  f := 0.668072012688580660e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 5);
  f := 0.5;
  testrel(5, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 10);
  f := 0.993790334674223865;
  testrel(6, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 20);
  f := 0.9999999999999680911;
  testrel(7, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 24);
  f := 1.0;
  testrel(8, NE, y, f, cnt,failed);


  mu := -1.5;
  sd := 0.5;
  y := normal_cdf(mu, sd, -20);
  f := 0.5725571222524576820e-299;
  testrel(9, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, -3);
  f := 0.1349898031630094527e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, -1);
  f := 0.8413447460685429486;
  testrel(11, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 0);
  f := 0.9986501019683699055;
  testrel(12, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 1);
  f := 0.9999997133484281208;
  testrel(13, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 2);
  f := 0.9999999999987201875;
  testrel(14, NE, y, f, cnt,failed);

  y := normal_cdf(mu, sd, 4);
  f := 1.0;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_normal_inv;
var
  y,f,mu,sd: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','normal_inv');

  mu := 5;
  sd := 2.0;
  y := normal_inv(mu, sd, ldexp(1,-20));;
  f := -4.526002068535627914;
  testrel(1, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.125);
  f := 2.699301239247983643;
  testrel(2, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.5);
  f := 5.0;
  testrel(3, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.75);
  f := 6.348979500392163486;
  testrel(4, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.875);
  f := 7.300698760752016357;
  testrel(5, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.9);
  f := 7.563103131089200934;
  testrel(6, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.99);
  f := 9.652695748081682202;
  testrel(7, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.9990234375);
  f := 11.19453815639756892;
  testrel(8, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 1-ldexp(1,-20));
  f := 14.52600206853562791;
  testrel(9, NE, y, f, cnt,failed);


  mu := -1.5;
  sd := 0.5;
  y := normal_inv(mu, sd, ldexp(1,-20));;
  f := -3.881500517133906978;
  testrel(10, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.125);
  f := -2.075174690188004089;
  testrel(11, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.5);
  f := -1.5;
  testrel(12, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.75);
  f := -1.162755124901959128;
  testrel(13, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.875);
  f := -0.9248253098119959110;
  testrel(14, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.9);
  f := -0.8592242172276997665;
  testrel(15, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 0.99);
  f := -0.3368260629795794496;
  testrel(16, 3, y, f, cnt,failed);         {!!!????}

  y := normal_inv(mu, sd, 0.9990234375);
  f := 0.4863453909939223118e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := normal_inv(mu, sd, 1-ldexp(1,-20));
  f := 0.8815005171339069785;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_beta_cdf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 2;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','beta_cdf');

  y := beta_cdf(1.2, 1.3, 0);
  f := 0;
  testabs(1, NE, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 1e-10);
  f := 1.3443494465428975103e-12;
  testrel( 2, 10, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 0.001);
  f := 3.3763004250453581348e-4;
  testrel( 3, NE, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 0.1);
  f := 8.3399782830674834581e-2;
  testrel( 4, NE, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 0.5);
  f := 5.2978142945129908088e-1;
  testrel( 5, NE, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 0.9);
  f := 9.3852939722443065921e-1;
  testrel( 6, NE, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 0.99);
  f := 9.9688643834125438042e-1;
  testrel( 7, NE, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 0.999);
  f := 9.9984379283306763424e-1;
  testrel( 8, NE, y, f, cnt,failed);

  y := beta_cdf(1.5, 4.3, 0.9);
  f := 0.99987810769284942730980; {qcalc}
  testrel( 9, NE, y, f, cnt,failed);

  y := beta_cdf(5, 20, 0.2);
  f := 5.4012267024242087821312e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := beta_cdf(8, 10, 0.2);
  f := 0.1093431523409920000e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := beta_cdf(8, 10, 0.8);
  f := 0.99950675026247680000;
  testrel(12, NE, y, f, cnt,failed);

  y := beta_cdf(1, 0.5, 0.9);
  f := 6.8377223398316206680e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := beta_cdf(0.1, 0.5, 0.99);
  f := 9.82283678979602092297725e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := beta_cdf(0.1, 0.5, 1);
  f := 1;
  testrel(15, NE, y, f, cnt,failed);

  y := beta_cdf(1e-4, 1e4, 1e-4);
  f := 9.99978059362107134278786e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := beta_cdf(1e-4, 1e4, 1e-3);
  f := 9.99999999586219825547e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := beta_cdf(1e-4, 1e4, 1e-10);
  f := 9.9867703308020004431e-1;
  testrel(18, NE, y, f, cnt,failed);

  {Extreme cases, reduced accuracy}
  y := beta_cdf(1.2, 1.3, 1e-100);
  f := 1.3443494465648959558e-120;
  testrel(19, 50, y, f, cnt,failed);

  y := beta_cdf(800, 0.5, 0.9);
  f := 1.5538647445875693462265e-38;
  testrel(20, 100, y, f, cnt,failed);

  y := beta_cdf(900, 900, 0.515625);
  f := 9.0757434037800433420e-1;
  testrel(21, 1, y, f, cnt,failed);

  {Case < 0 and >1, not allowed for ibetax}
  y := beta_cdf(1.2, 1.3, -1);
  f := 0;
  testabs(22, 1, y, f, cnt,failed);

  y := beta_cdf(3, 4, -5);
  f := 0;
  testabs(23, 1, y, f, cnt,failed);

  y := beta_cdf(1.2, 1.3, 2);
  f := 1;
  testabs(24, 1, y, f, cnt,failed);

  y := beta_cdf(3, 4, 5);
  f := 1;
  testabs(25, 1, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_beta_inv;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','beta_inv');

  y := beta_inv(2, 1, 0.1);
  f := 3.162277660168379332e-1;
  testrel( 1, NE, y, f, cnt,failed);

  y := beta_inv(1, 2, 0.01);
  f := 5.012562893380045266e-3;
  testrel( 2, NE, y, f, cnt,failed);

  y := beta_inv(2, 1, 0.01);
  f := 10e-2;
  testrel( 3, NE, y, f, cnt,failed);

  y := beta_inv(2, 1, 0.05);
  f := 2.236067977499789696e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := beta_inv(2, 1, 0.99);
  f := 1.989974874213239909*0.5;
  testrel( 5, NE, y, f, cnt,failed);

  y := beta_inv(5, 10, 0.2);
  f := 0.22912109609099720490;
  testrel( 6, NE, y, f, cnt,failed);

  y := beta_inv(1, 1, 0);
  f := 0;
  testabs( 7, 1, y, f, cnt,failed);

  y := beta_inv(1, 1, 1);
  f := 1;
  testabs( 8, 1, y, f, cnt,failed);

  y := beta_inv(4, 7, 0.8);
  f := 4.836572374148516672e-1;
  testrel( 9, NE, y, f, cnt,failed);

  y := beta_inv(1.5, 4.25, 0.99);
  f := 0.7194720681939346506;
  testrel(10, NE, y, f, cnt,failed);

  y := beta_inv(5, 20, 0.5);
  f := 1.919240597664992704e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := beta_inv(8, 10, 0.2);
  f := 3.449977502489844803e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := beta_inv(8, 10, 0.8);
  f := 5.427825084949223437e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := beta_inv(1, 0.5, 0.9);
  f := 0.99;
  testrel(14, NE, y, f, cnt,failed);

  y := beta_inv(1800, 0.5, 0.9);
  f := 9.999956130742314815e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := beta_inv(0.1, 0.5, 0.98);
  f := 9.872767968062899295e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := beta_inv(0.5, 0.125, 0.9);
  f := 9.999999663947204350e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := beta_inv(900, 900, 0.9);
  f := 5.151018817560778887e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := beta_inv(1800, 0.5, 1e-4);
  f := 9.958036054195801287e-1;
  testrel(19, NE, y, f, cnt,failed);

  y := beta_inv(4500,40,0.125);
  f := 0.9895767941483560347;
  testrel(20, NE, y, f, cnt,failed);

  y := beta_inv(100,102,0.5);
  f := 0.4950331343781250235;
  testrel(21, NE, y, f, cnt,failed);

  y := beta_inv(100,102,0.875);
  f := 0.5355056613331838346;
  testrel(22, NE, y, f, cnt,failed);

  y := beta_inv(100,300,0.5);
  f := 0.2495830001077997259;
  testrel(23, NE, y, f, cnt,failed);

  y := beta_inv(100,300,0.625);
  f := 0.2565236334554090243;
  testrel(24, NE, y, f, cnt,failed);

  y := beta_inv(1000,200,0.125);
  f := 0.8209035891017542587;
  testrel(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_beta_pdf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 25;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','beta_pdf');

  y := beta_pdf(2, 3, 0.1);
  f := 0.972;
  testrel( 1, NE, y, f, cnt,failed);

  y := beta_pdf(3, 2, 0.01);
  f := 0.1188e-2;
  testrel( 2, NE, y, f, cnt,failed);

  y := beta_pdf(2, 3, 0.01);
  f := 0.117612;
  testrel( 3, NE, y, f, cnt,failed);

  y := beta_pdf(2, 3, 0.05);
  f := 0.541500;
  testrel( 4, NE, y, f, cnt,failed);

  y := beta_pdf(2, 3, 0.99);
  f := 0.1188e-2;
  testrel( 5, NE, y, f, cnt,failed);

  y := beta_pdf(5, 10, 0.2);
  f := 2.1496311316480;
  testrel( 6, NE, y, f, cnt,failed);

  y := beta_pdf(4, 1, 0);
  f := 0;
  testabs( 7, 1, y, f, cnt,failed);

  y := beta_pdf(4, 1, 1);
  f := 4;
  testabs( 8, 1, y, f, cnt,failed);

  y := beta_pdf(4, 7, 0.8);
  f := 0.27525120e-1;
  testrel( 9, NE, y, f, cnt,failed);

  y := beta_pdf(1.5, 4.3, 0.9);
  f := 0.51874622029711676062e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := beta_pdf(5, 20, 0.5);
  f := 0.25334358215332031250e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := beta_pdf(8, 10, 0.2);
  f := 0.3341140958904320;
  testrel(12, NE, y, f, cnt,failed);

  y := beta_pdf(8, 10, 0.8);
  f := 0.208821309931520e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := beta_pdf(1, 0.5, 0.9);
  f := 1.5811388300841896660;
  testrel(14, NE, y, f, cnt,failed);

  y := beta_pdf(0.1, 0.5, 0.98);
  f := 0.63594080282711715084;
  testrel(15, NE, y, f, cnt,failed);

  y := beta_pdf(0.5, 0.125, 0.9);
  f := 0.84915665654733553012;
  testrel(16, NE, y, f, cnt,failed);

  {extreme cases}
  y := beta_pdf(50, 50, 0.9);
  f := 0.14443647504225557881e-20;
  testrel(17, 200, y, f, cnt,failed);

  y := beta_pdf(180, 0.5, 0.9);
  f := 0.15423001926734076804e-6;
  testrel(18, 200, y, f, cnt,failed);

  y := beta_pdf(180, 0.5, 0.1);
  f := 0.79733066736098512657e-178;
  testrel(19, 200, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_t_cdf;
var
  y,f: double;
  cnt, failed: integer;
  nu: longint;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','t_cdf');

  {statevalf[cdf,studentst[nu]](t);}

  nu := 1;

  y := t_cdf(nu, -1e5);
  f := 0.31830988617318034200e-5;
  testrel( 1, NE, y, f, cnt,failed);

  y := t_cdf(nu, -10);
  f := 0.31725517430553569513e-1;
  testrel( 2, NE, y, f, cnt,failed);

  y := t_cdf(nu, -2.0);
  f := 0.14758361765043327417;
  testrel( 3, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.5);
  f := 0.18716704181099881618;
  testrel( 4, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.0);
  f := 0.25000000000000000000;
  testrel( 5, NE, y, f, cnt,failed);

  y := t_cdf(nu, -0.5);
  f := 0.35241638234956672583;
  testrel( 6, NE, y, f, cnt,failed);

  y := t_cdf(nu, 0);
  f := 0.5;
  testrel( 7, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1.0);
  f := 0.75000000000000000000;
  testrel( 8, NE, y, f, cnt,failed);

  y := t_cdf(nu, 2.0);
  f := 0.85241638234956672583;
  testrel( 9, NE, y, f, cnt,failed);

  y := t_cdf(nu, 4.0);
  f := 0.92202086962263067454;
  testrel(10, NE, y, f, cnt,failed);

  y := t_cdf(nu, 10.0);
  f := 0.9682744825694464304;
  testrel(11, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1e5);
  f := 0.99999681690113826820;
  testrel(12, NE, y, f, cnt,failed);

  nu := 2;

  y := t_cdf(nu, -1e5);
  f := 0.49999999992500000001e-10;
  testrel(13, NE, y, f, cnt,failed);

  y := t_cdf(nu, -10);
  f := 0.49262285116628454233e-2;
  testrel(14, NE, y, f, cnt,failed);

  y := t_cdf(nu, -2.0);
  f := 0.91751709536136983634e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.5);
  f := 0.13619656244550053972;
  testrel(16, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.0);
  f := 0.21132486540518711775;
  testrel(17, NE, y, f, cnt,failed);

  y := t_cdf(nu, -0.5);
  f := 0.33333333333333333333;
  testrel(18, NE, y, f, cnt,failed);

  y := t_cdf(nu, 0);
  f := 0.5;
  testrel(19, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1.0);
  f := 0.78867513459481288225;
  testrel(20, NE, y, f, cnt,failed);

  y := t_cdf(nu, 2.0);
  f := 0.90824829046386301637;
  testrel(21, NE, y, f, cnt,failed);

  y := t_cdf(nu, 4.0);
  f := 0.97140452079103168293;
  testrel(22, NE, y, f, cnt,failed);

  y := t_cdf(nu, 10.0);
  f := 0.99507377148833715458;
  testrel(23, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1e5);
  f := 0.99999999995000000001;
  testrel(24, NE, y, f, cnt,failed);

  nu := 5;

  y := t_cdf(nu, -1e5);
  f := 0.94901672353943244519e-24;
  testrel(25, NE, y, f, cnt,failed);

  y := t_cdf(nu, -10);
  f := 0.85473787871481795358e-4;
  testrel(26, NE, y, f, cnt,failed);

  y := t_cdf(nu, -2.0);
  f := 0.5096973941492917812e-1;
  testrel(27, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.5);
  f := 0.9695184012123671605e-1;
  testrel(28, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.0);
  f := 0.18160873382456131279;
  testrel(29, NE, y, f, cnt,failed);

  y := t_cdf(nu, -0.5);
  f := 0.31914943582046450336;
  testrel(30, NE, y, f, cnt,failed);

  y := t_cdf(nu, 0);
  f := 0.5;
  testrel(31, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1.0);
  f := 0.81839126617543868721;
  testrel(32, NE, y, f, cnt,failed);

  y := t_cdf(nu, 2.0);
  f := 0.94903026058507082188;
  testrel(33, NE, y, f, cnt,failed);

  y := t_cdf(nu, 4.0);
  f := 0.99483829225958427310;
  testrel(34, NE, y, f, cnt,failed);

  y := t_cdf(nu, 10.0);
  f := 0.99991452621212851820;
  testrel(35, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1e5);
  f := 1.0;
  testrel(36, NE, y, f, cnt,failed);

  nu := 20;

  y := t_cdf(nu, -1e5);
  f := 0.90212888902753924405e-88;
  testrel(37, NE, y, f, cnt,failed);

  y := t_cdf(nu, -10);
  f := 0.15818908793571940812e-8;
  testrel(38, NE, y, f, cnt,failed);

  y := t_cdf(nu, -2.0);
  f := 0.29632767723285236484e-1;
  testrel(39, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.5);
  f := 0.74617885584626265112e-1;
  testrel(40, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.0);
  f := 0.16462828858585453212;
  testrel(41, NE, y, f, cnt,failed);

  y := t_cdf(nu, -0.5);
  f := 0.31126592114051179869;
  testrel(42, NE, y, f, cnt,failed);

  y := t_cdf(nu, 0);
  f := 0.5;
  testrel(43, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1.0);
  f := 0.83537171141414546788;
  testrel(44, NE, y, f, cnt,failed);

  y := t_cdf(nu, 2.0);
  f := 0.97036723227671476353;
  testrel(45, NE, y, f, cnt,failed);

  y := t_cdf(nu, 4.0);
  f := 0.99964823835343584085;
  testrel(46, NE, y, f, cnt,failed);

  y := t_cdf(nu, 10.0);
  f := 0.99999999841810912065;
  testrel(47, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1e5);
  f := 1.0;
  testrel(48, NE, y, f, cnt,failed);

  y := t_cdf(120, -1.0);
  f := 0.15966136193221061885;
  testrel(49, NE, y, f, cnt,failed);

  y := t_cdf(550, -1.0);
  f := 0.15887512729919423702;
  testrel(50, NE, y, f, cnt,failed);

  {extreme cases}
  y := t_cdf(120000, -1.0);
  f := 0.15865626214070873891;
  testrele(51, 1e-14, y, f, cnt,failed);

  y := t_cdf(120000000, -1.0);
  f := 0.15865525493965227849;
  testrele(52, 1e-11, y, f, cnt,failed);


  nu := 25000;

  y := t_cdf(nu, -30);
  f := 0.1361571814548456598e-193;
  testrel(53, NE, y, f, cnt,failed);

  y := t_cdf(nu, -10);
  f := 0.8435744407871687785e-23;
  testrel(54, NE, y, f, cnt,failed);

  y := t_cdf(nu, -2.0);
  f := 0.2275553114739848001e-1;
  testrel(55, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.5);
  f := 0.668135152014062988e-1;
  testrel(56, NE, y, f, cnt,failed);

  y := t_cdf(nu, -1.0);
  f := 0.1586600932975529668;
  testrel(57, NE, y, f, cnt,failed);

  y := t_cdf(nu, -0.5);
  f := 0.308539739120663701;
  testrel(58, 2*NE, y, f, cnt,failed);

  y := t_cdf(nu, 0);
  f := 0.5;
  testrel(59, NE, y, f, cnt,failed);

  y := t_cdf(nu, 1.0);
  f := 0.8413399067024470332;
  testrel(60, NE, y, f, cnt,failed);

  y := t_cdf(nu, 2.0);
  f := 0.9772444688526015200;
  testrel(61, NE, y, f, cnt,failed);

  y := t_cdf(nu, 4.0);
  f := 0.9999682376606739350;
  testrel(62, NE, y, f, cnt,failed);

  y := t_cdf(nu, 8.0);
  f := 0.9999999999999993511;
  testrel(63, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_t_inv;
var
  y,f: double;
  cnt, failed: integer;
  nu: longint;
const
  NE  = 1;
  NE2 = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','t_inv');

  {statevalf[icdf,studentst[nu]](p);}

  nu := 1;

  y := t_inv(nu, 1e-6);
  f := -318309.88618274347399;
  testrel( 1, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.001953125);
  f := -162.97261641324996316;
  testrel( 2, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.03125);
  f := -10.153170387608860462;
  testrel( 3, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.125);
  f := -2.4142135623730950488;
  testrel( 4, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.25);
  f := -1;
  testrel( 5, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.5);
  f := 0;
  testrel( 6, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.75);
  f := 1;
  testrel( 7, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.9);
  f := 3.0776835371752534026;
  testrel( 8, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.95);                   {0.95 not exact}
  f := 6.3137515146750430990;
  testrel( 9, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.99);                   {0.95 not exact}
  f := 31.820515953773958040;
  testrel(10, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.9990234375);
  f := 325.94830079770134898;
  testrel(11, NE, y, f, cnt,failed);

  y := t_inv(nu, 1-ldexp(1,-20)); {0.99999904632568359375}
  f := 333772.10721405580178;
  testrel(12, NE, y, f, cnt,failed);

  nu := 2;

  y := t_inv(nu, 1e-6);
  f := -707.10572052593380254;
  testrel(13, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.001953125);
  f := -15.953086800791288032;
  testrel(14, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.03125);
  f := -3.8100038100057150095;
  testrel(15, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.125);
  f := -1.6035674514745463081;
  testrel(16, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.25);
  f := -0.81649658092772603273;
  testrel(17, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.5);
  f := 0;
  testrel(18, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.75);
  f := 0.81649658092772603273;
  testrel(19, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.9);
  f := 1.8856180831641267317;
  testrel(20, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.95);                         {0.95 not exact}
  f := 2.9199855803537256870;
  testrel(21, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.99);                         {0.95 not exact}
  f := 6.9645567342832741871;
  testrel(22, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.9990234375);
  f := 22.594257871383013732;
  testrel(23, NE, y, f, cnt,failed);

  y := t_inv(nu, 1-ldexp(1,-20)); {0.99999904632568359375}
  f := 724.07630813366407066;
  testrel(24, NE, y, f, cnt,failed);

  nu := 5;

  y := t_inv(nu, 1e-6);
  f := -24.771029720515944171;
  testrel(25, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.001953125);
  f := -5.0581169399535194651;
  testrel(26, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.03125);
  f := -2.3885469271281262464;
  testrel(27, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.125);
  f := -1.3009490369230305560;
  testrel(28, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.25);
  f := -0.72668684380042265302;
  testrel(29, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.5);
  f := 0;
  testrel(30, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.75);
  f := 0.72668684380042265302;
  testrel(31, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.9);
  f := 1.4758840488244810785;
  testrel(32, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.95);                     {0.95 not exact}
  f := 2.0150483733330242378;
  testrel(33, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.99);                     {0.95 not exact}
  f := 3.3649299989072185928;
  testrel(34, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.9990234375);
  f := 5.9248458020727425071;
  testrel(35, NE, y, f, cnt,failed);

  y := t_inv(nu, 1-ldexp(1,-20)); {0.99999904632568359375}
  f := 25.008780822627182395;
  testrel(36, NE, y, f, cnt,failed);

  nu := 20;

  y := t_inv(nu, 1e-6);
  f := -6.5965301483587757450;
  testrel(37, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.001953125);
  f := -3.2615195635002924611;
  testrel(38, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.03125);
  f := -1.9728136094297307620;
  testrel(39, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.125);
  f := -1.1847614343569054294;
  testrel(40, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.25);
  f := -0.68695449644880341996;
  testrel(41, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.5);
  f := 0;
  testrel(42, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.75);
  f := 0.68695449644880341996;
  testrel(43, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.9);
  f := 1.3253407069850463431;
  testrel(44, NE, y, f, cnt,failed);

  y := t_inv(nu, 0.95);                      {0.95 not exact}
  f := 1.7247182429207872728;
  testrel(45, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.99);                      {0.99 not exact}
  f := 2.5279770027415734818;
  testrel(46, NE2, y, f, cnt,failed);

  y := t_inv(nu, 0.9990234375);
  f := 3.5620321200529290161;
  testrel(47, NE, y, f, cnt,failed);

  y := t_inv(nu, 1-ldexp(1,-20)); {0.99999904632568359375}
  f := 6.6189291385494687291;
  testrel(48, NE, y, f, cnt,failed);

  y := t_inv(100, 0.9921875);
  f := 2.4596086534573194710;
  testrel(49, NE, y, f, cnt,failed);

  {extreme cases}
  y := t_inv(10000, 0.9921875);
  f := 2.4179727636644442344;
  testrele(50, 1e-17, y, f, cnt,failed);

  y := t_inv(120000, 0.9921875);
  f := 2.4175934900459174302;
  testrele(51, 1e-15, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_t_pdf;
var
  y,f: double;
  cnt, failed: integer;
  nu: longint;
const
  NE = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','t_pdf');

  {statevalf[pdf,studentst[nu]](t);}

  nu := 1;

  y := t_pdf(nu, -1e5);
  f := 0.31830988615195968291e-10;
  testrel( 1, NE, y, f, cnt,failed);

  y := t_pdf(nu, -10);
  f := 0.31515830315226799161e-2;
  testrel( 2, NE, y, f, cnt,failed);

  y := t_pdf(nu, -2.0);
  f := 0.63661977236758134306e-1;
  testrel( 3, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.5);
  f := 0.97941503441166360472e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.0);
  f := 0.15915494309189533577;
  testrel( 5, NE, y, f, cnt,failed);

  y := t_pdf(nu, -0.5);
  f := 0.25464790894703253722;
  testrel( 6, NE, y, f, cnt,failed);

  y := t_pdf(nu, 0);
  f := 0.31830988618379067153;
  testrel( 7, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1.0);
  f := 0.15915494309189533577;
  testrel( 8, NE, y, f, cnt,failed);

  y := t_pdf(nu, 2.0);
  f := 0.63661977236758134306e-1;
  testrel( 9, NE, y, f, cnt,failed);

  y := t_pdf(nu, 4.0);
  f := 0.18724110951987686560e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := t_pdf(nu, 10.0);
  f := 0.31515830315226799161e-2;
  testrel(11, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1e5);
  f := 0.31830988615195968291e-10;
  testrel(12, NE, y, f, cnt,failed);

  nu := 2;

  y := t_pdf(nu, -1e5);
  f := 0.99999999970000000006e-15;
  testrel(13, NE, y, f, cnt,failed);

  y := t_pdf(nu, -10);
  f := 0.97073288527124932269e-3;
  testrel(14, NE, y, f, cnt,failed);

  y := t_pdf(nu, -2.0);
  f := 0.68041381743977169398e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.5);
  f := 0.11413441178180375225;
  testrel(16, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.0);
  f := 0.19245008972987525482;
  testrel(17, NE, y, f, cnt,failed);

  y := t_pdf(nu, -0.5);
  f := 0.29629629629629629630;
  testrel(18, NE, y, f, cnt,failed);

  y := t_pdf(nu, 0);
  f := 0.35355339059327376220;
  testrel(19, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1.0);
  f := 0.19245008972987525482;
  testrel(20, NE, y, f, cnt,failed);

  y := t_pdf(nu, 2.0);
  f := 0.68041381743977169398e-1;
  testrel(21, NE, y, f, cnt,failed);

  y := t_pdf(nu, 4.0);
  f := 0.13094570021973102304e-1;
  testrel(22, NE, y, f, cnt,failed);

  y := t_pdf(nu, 10.0);
  f := 0.97073288527124932269e-3;
  testrel(23, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1e5);
  f := 0.99999999970000000006e-15;
  testrel(24, NE, y, f, cnt,failed);

  nu := 5;

  y := t_pdf(nu, -1e5);
  f := 0.47450836156635549626e-28;
  testrel(25, NE, y, f, cnt,failed);

  y := t_pdf(nu, -10);
  f := 0.40989816415343314024e-4;
  testrel(26, NE, y, f, cnt,failed);

  y := t_pdf(nu, -2.0);
  f := 0.65090310326216466251e-1;
  testrel(27, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.5);
  f := 0.12451734464635513754;
  testrel(28, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.0);
  f := 0.21967979735098057360;
  testrel(29, NE, y, f, cnt,failed);

  y := t_pdf(nu, -0.5);
  f := 0.32791853132274651219;
  testrel(30, NE, y, f, cnt,failed);

  y := t_pdf(nu, 0);
  f := 0.37960668982249443118;
  testrel(31, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1.0);
  f := 0.21967979735098057360;
  testrel(32, NE, y, f, cnt,failed);

  y := t_pdf(nu, 2.0);
  f := 0.65090310326216466251e-1;
  testrel(33, NE, y, f, cnt,failed);

  y := t_pdf(nu, 4.0);
  f := 0.51237270519179142530e-2;
  testrel(34, NE, y, f, cnt,failed);

  y := t_pdf(nu, 10.0);
  f := 0.40989816415343314024e-4;
  testrel(35, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1e5);
  f := 0.47450836156635549626e-28;
  testrel(36, NE, y, f, cnt,failed);

  nu := 20;

  y := t_pdf(nu, -1e5);
  f := 0.18042577746105863731e-91;
  testrel(37, NE, y, f, cnt,failed);

  y := t_pdf(nu, -10);
  f := 0.26600849800550893824e-8;
  testrel(38, NE, y, f, cnt,failed);

  y := t_pdf(nu, -2.0);
  f := 0.58087215247356954399e-1;
  testrel(39, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.5);
  f := 0.12862738297214604673;
  testrel(40, NE, y, f, cnt,failed);

  y := t_pdf(nu, -1.0);
  f := 0.23604564912670098385;
  testrel(41, NE, y, f, cnt,failed);

  y := t_pdf(nu, -0.5);
  f := 0.34580861238374166972;
  testrel(42, NE, y, f, cnt,failed);

  y := t_pdf(nu, 0);
  f := 0.39398858571143259539;
  testrel(43, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1.0);
  f := 0.23604564912670098385;
  testrel(44, NE, y, f, cnt,failed);

  y := t_pdf(nu, 2.0);
  f := 0.58087215247356954399e-1;
  testrel(45, NE, y, f, cnt,failed);

  y := t_pdf(nu, 4.0);
  f := 0.82247430013313934155e-3;
  testrel(46, NE, y, f, cnt,failed);

  y := t_pdf(nu, 10.0);
  f := 0.26600849800550893824e-8;
  testrel(47, NE, y, f, cnt,failed);

  y := t_pdf(nu, 1e5);
  f := 0.18042577746105863731e-91;
  testrel(48, NE, y, f, cnt,failed);

  y := t_pdf(120, -1.0);
  f := 0.24096600519662408799;
  testrel(49, NE, y, f, cnt,failed);

  y := t_pdf(550, -1.0);
  f := 0.24175091768893257456;
  testrel(50, NE, y, f, cnt,failed);

  {extreme cases}
  y := t_pdf(120000, -1.0);
  f := 0.24196971631129191181;
  testrel(51, NE, y, f, cnt,failed);

  y := t_pdf(120000000, -1.0);
  f := 0.24197072351093200114;
  testrel(52, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_f_pdf;
var
  y,f: double;
  cnt, failed: integer;
  nu1, nu2: longint;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','f_pdf');

  {statevalf[pdf,fratio[nu1,nu2]](x);}

  nu1 := 10;
  nu2 := 15;

  y := f_pdf(nu1, nu2, 1/1024);
  f := 0.36198487346066596768e-9;
  testrel( 1, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.125);
  f := 0.36019453612611847297e-1;
  testrel( 2, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.5);
  f := 0.68796953428063353340;
  testrel( 3, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.75);
  f := 0.79896492175019497702;
  testrel( 4, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 1);
  f := 0.67657200953522187979;
  testrel( 5, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 2);
  f := 0.16137399972081001477;
  testrel( 6, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 10);
  f := 0.35143185961678413375e-4;
  testrel( 7, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 50);
  f := 0.15953419650749369941e-9;
  testrel( 8, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 100);
  f := 0.52934474606409138305e-12;
  testrel( 9, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 500);
  f := 0.70316530729462402384e-18;
  testrel(10, NE, y, f, cnt,failed);

  nu1 := 100;
  nu2 := 2;

  y := f_pdf(nu1, nu2, 1/1024);
  f := 0.12213827580441443014e-61;
  testrel(11, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.125);
  f := 0.33023456147108613310e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.5);
  f := 0.54120236666630659326;
  testrel(13, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.75);
  f := 0.46448845334507355336;
  testrel(14, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 1);
  f := 0.36424302169309984212;
  testrel(14, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 2);
  f := 0.15050465957647268865;
  testrel(16, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 10);
  f := 0.90312154274027085984e-2;
  testrel(17, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 50);
  f := 0.39192426751861425846e-3;
  testrel(18, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 100);
  f := 0.98985285309689047108e-4;
  testrel(19, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 500);
  f := 0.39918484804061809859e-5;
  testrel(20, NE, y, f, cnt,failed);

  nu1 := 6;
  nu2 := 4;

  y := f_pdf(nu1, nu2, 1/1024);
  f := 0.38342159511298883228e-4;
  testrel(21, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.125);
  f := 0.26798282298082588781;
  testrel(22, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.5);
  f := 0.61688582138394716487;
  testrel(23, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 0.75);
  f := 0.52575435413566295761;
  testrel(24, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 1);
  f := 0.41472000000000000000;
  testrel(25, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 2);
  f := 0.15820312500000000000;
  testrel(26, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 10);
  f := 0.38623809814453125000e-2;
  testrel(27, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 50);
  f := 0.39932552424196286174e-4;
  testrel(28, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 100);
  f := 0.51590565311656193475e-5;
  testrel(29, NE, y, f, cnt,failed);

  y := f_pdf(nu1, nu2, 500);
  f := 0.42383356469663650962e-7;
  testrel(30, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_f_cdf;
var
  y,f: double;
  cnt, failed: integer;
  nu1, nu2: longint;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','f_cdf');

  {statevalf[pdf,fratio[nu1,nu2]](x);}

  nu1 := 10;
  nu2 := 15;

  y := f_cdf(nu1, nu2, 1/1024);
  f := 0.70796121838458865266e-13;
  testrel( 1, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.125);
  f := 0.10696408932185717088e-2;
  testrel( 2, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.5);
  f := 0.13538507460415484644;
  testrel( 3, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.75);
  f := 0.32852562335634297149;
  testrel( 4, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 1);
  f := 0.51551998894148855316;
  testrel( 5, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 2);
  f := 0.89085985394213032552;
  testrel( 6, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 10);
  f := 0.99994212465657099760;
  testrel( 7, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 50);
  f := 0.99999999888891906392;
  testrel( 8, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 100);
  f := 0.99999999999278539289;
  testrel( 9, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 500);
  f := 0.99999999999999995292;
  testrel(10, NE, y, f, cnt,failed);

  nu1 := 100;
  nu2 := 2;

  y := f_cdf(nu1, nu2, 1/1024);
  f := 0.250199333598978229e-66;
  testrel(11, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.125);
  f := 0.598550142666343616e-3;
  testrel(12, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.5);
  f := 0.14071261533323971426;
  testrel(13, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.75);
  f := 0.26824208180677997706;
  testrel(14, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 1);
  f := 0.37152788212696183896;
  testrel(15, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 2);
  f := 0.60803882468894966212;
  testrel(16, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 10);
  f := 0.90492778582575140156;
  testrel(17, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 50);
  f := 0.98020259306405426042;
  testrel(18, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 100);
  f := 0.99005082366750984917;
  testrel(19, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 500);
  f := 0.99800203858634930828;
  testrel(20, NE, y, f, cnt,failed);

  nu1 := 6;
  nu2 := 4;

  y := f_cdf(nu1, nu2, 1/1024);
  f := 0.12504032119315787138e-7;
  testrel(21, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.125);
  f := 0.13881108954044244596e-1;
  testrel(22, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.5);
  f := 0.21366097459391920033;
  testrel(23, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 0.75);
  f := 0.35786209456304402484;
  testrel(24, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 1);
  f := 0.47520000000000000000;
  testrel(25, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 2);
  f := 0.73828125000000000000;
  testrel(26, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 10);
  f := 0.97846984863281250000;
  testrel(27, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 50);
  f := 0.99897935314531042579;
  testrel(28, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 100);
  f := 0.99973917146437760256;
  testrel(29, NE, y, f, cnt,failed);

  y := f_cdf(nu1, nu2, 500);
  f := 0.99998938059887170883;
  testrel(30, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_f_inv;
var
  y,f: double;
  cnt, failed: integer;
  nu1, nu2: longint;
const
  NE   = 3;
  NE2  = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','f_inv');

  {statevalf[icdf, fratio[nu1,nu2]](x);}

  nu1 := 10;
  nu2 := 15;

  y := f_inv(nu1, nu2, 1e-6);
  f := 0.27250601258372329764e-1;
  testrel( 1, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.001953125);
  f := 0.14449168758699289717;
  testrel( 2, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.03125);
  f := 0.30351965339193473178;
  testrel( 3, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.125);
  f := 0.48469389512563697769;
  testrel( 4, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.25);
  f := 0.65198482424610119623;
  testrel( 5, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.5);
  f := 0.97731605405255776796;
  testrel( 6, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.75);
  f := 1.4490744742928857851;
  testrel( 7, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.9);
  f := 2.0593194961568176081;
  testrel( 8, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.95);
  f := 2.5437185496928079783;
  testrel( 9, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.99);
  f := 3.8049397459502750704;
  testrel(10, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.9990234375);
  f := 6.1078009092810849030;
  testrel(11, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 1-ldexp(1,-20)); {0.99999904632568359375}
  f := 18.975080250131432137;
  testrel(12, NE, y, f, cnt,failed);

  nu1 := 100;
  nu2 := 2;

  y := f_inv(nu1, nu2, 0.001953125);
  f := 0.15050733921182494928;
  testrel(13, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.03125);
  f := 0.27865452345825932003;
  testrel(14, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.125);
  f := 0.47096765968297816859;
  testrel(15, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.25);
  f := 0.71139372966448604356;
  testrel(16, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.5);
  f := 1.4327181457209769326;
  testrel(17, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.75);
  f := 3.4660690861793311186;
  testrel(18, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.9);
  f := 9.4812250930468317113;
  testrel(19, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.95);
  f := 19.485727456000138944;
  testrel(20, NE2, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.99);
  f := 99.489162808433367622;
  testrel(21, NE2, y, f, cnt,failed);

  nu1 := 6;
  nu2 := 4;

  y := f_inv(nu1, nu2, 1e-6);
  f := 0.42330914491139243793e-2;
  testrel(22, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.001953125);
  f := 0.58281870550544859393e-1;
  testrel(23, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.03125);
  f := 0.17731846391979853619;
  testrel(24, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.125);
  f := 0.35652611565092988244;
  testrel(25, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.25);
  f := 0.55954869735653314273;
  testrel(26, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.5);
  f := 1.0616688782738316091;
  testrel(27, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.75);
  f := 2.0765682540568617845;
  testrel(28, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.9);
  f := 4.0097493126739445490;
  testrel(29, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.95);
  f := 6.1631322826886348580;
  testrel(20, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.99);
  f := 15.206864861157530011;
  testrel(31, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 0.9990234375);
  f := 51.141068844668912659;
  testrel(32, NE, y, f, cnt,failed);

  y := f_inv(nu1, nu2, 1-ldexp(1,-20)); {0.99999904632568359375}
  f := 1671.0737758294571665;
  testrele(33, 1e-16, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.
