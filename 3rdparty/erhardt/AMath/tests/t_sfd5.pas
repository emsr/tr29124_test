{Part 5 of regression test for SPECFUN unit  (c) 2010-2016  W.Ehrhardt}

unit t_sfd5;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_chebyshev_t;
procedure test_chebyshev_u;
procedure test_chebyshev_v;
procedure test_chebyshev_w;
procedure test_chebfun;
procedure test_gegenbauer_c;
procedure test_hermite_h;
procedure test_hermite_he;
procedure test_jacobi_p;
procedure test_laguerre;
procedure test_zernike_r;

procedure test_orthopoly;


implementation


uses
  amath, specfun, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_chebyshev_t;
var
  x,y,f: double;
  cnt,failed,i: integer;
const
  NE = 1;
  NE2 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chebyshev_t');

  i := 16384;
  x := 1 - ldexp(1,-30);
  y := chebyshev_t(i,x);
  f := 0.7602445970399789147;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-48);
  y := chebyshev_t(16384,x);
  f := 0.9999990463258351762;
  testrel( 2, NE, y, f, cnt,failed);

  i := 128;
  x := 0.25;
  y := chebyshev_t(i,x);
  f := 0.6001192699014801455;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-30);
  y := chebyshev_t(i,x);
  f := 0.9999847412497401993;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-48);
  y := chebyshev_t(i,x);
  f := 0.9999999999417923391;
  testrel( 5, NE, y, f, cnt,failed);

  i := 64;
  x := 0.25;
  y := chebyshev_t(i,x);
  f := -0.8944605273295966260;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.75;
  y := chebyshev_t(i,x);
  f := -0.6456684889424518930;
  testrel( 7, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-30);
  y := chebyshev_t(i,x);
  f := 0.9999961853051591020;
  testrel( 8, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-48);
  y := chebyshev_t(i,x);
  f := 0.9999999999854480848;
  testrel( 9, NE, y, f, cnt,failed);

  x := 1.25;
  y := chebyshev_t(70,x);
  f := 590295810358705651712.0;
  testrel(10, NE, y, f, cnt,failed);

  x := 1.25;
  y := chebyshev_t(50,x);
  f := 562949953421312.0;
  testrel(11, NE, y, f, cnt,failed);

  x := -0.375;
  y := chebyshev_t(15,x);
  f := -0.4944775826297700405120849609375;
  testrel(12, NE, y, f, cnt,failed);

  x := cosPi(5/20);
  y := chebyshev_t(10,x);
  f := 0.0;
  testrel(13, NE2, y, f, cnt,failed);

  x := cosPi(-15/42);
  y := chebyshev_t(21,x);
  f := 0.0;
  testrel(14, NE2, y, f, cnt,failed);

  x := cosPi(31/100);
  y := chebyshev_t(50,x);
  f := 0.0;
  testrel(15, NE2, y, f, cnt,failed);

  y := chebyshev_t(20,1e-7);
  f := 0.999999999998;
  testrel(16, NE, y, f, cnt,failed);

  y := chebyshev_t(21,1e-10);
  f := 0.2099999999999999998e-8;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_chebyshev_u;
var
  x,y,f: double;
  cnt,failed,i: integer;
const
  NE  = 1;
  NE2 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chebyshev_u');

  i := 16384;
  x := 1 - ldexp(1,-30);
  y := chebyshev_u(i,x);
  f := 15053.1566949272367874;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-48);
  y := chebyshev_u(16384,x);
  f := 16384.99479071345040;
  testrel( 2, NE, y, f, cnt,failed);

  i := 128;
  x := 0.25;
  y := chebyshev_u(i,x);
  f := 0.3935832582081333238;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-30);
  y := chebyshev_u(i,x);
  f := 128.9993336211039270;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-48);
  y := chebyshev_u(i,x);
  f := 128.9999999974579623;
  testrel( 5, NE, y, f, cnt,failed);

  i := 64;
  x := 0.25;
  y := chebyshev_u(i,x);
  f := -0.7790076899025733705;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.75;
  y := chebyshev_u(i,x);
  f := +0.2201927521250710782;
  testrel( 7, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-30);
  y := chebyshev_u(i,x);
  f := 64.9999147653914774;
  testrel( 8, NE, y, f, cnt,failed);

  x := 1 - ldexp(1,-48);
  y := chebyshev_u(i,x);
  f := 64.9999999996748556;
  testrel( 9, NE, y, f, cnt,failed);

  x := 1.25;
  y := chebyshev_u(70,x);
  f := 0.1574122160956548406e22;
  testrel(10, NE, y, f, cnt,failed);

  x := 1.25;
  y := chebyshev_u(50,x);
  f := 0.1501199875790165333e16;
  testrel(11, NE, y, f, cnt,failed);

  x := -0.375;
  y := chebyshev_u(15,x);
{  f := double(f12);}
  f := -1.4287276100367307663E-1;
  testrel(12, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-60);
  y := chebyshev_u(10,x);
  f := 10.99999999999999962;
  testrel(13, NE, y, f, cnt,failed);

  y := chebyshev_u(65,-1.25);
  f := -0.4919131752989213764e20;
  testrel(14, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-63);
  y := chebyshev_u(65,-x);
  f := -66.00000000000001039;
  testrel(15, NE, y, f, cnt,failed);

  x := -0.65486073394528506406;
  y := chebyshev_u(10,x);
  f := 0.0;
  testrel(16, NE2, y, f, cnt,failed);

  x := 0.3826834323650897717;
  y := chebyshev_u(15,x);
  f := 0.0;
  testrel(17, NE2, y, f, cnt,failed);

  x := 0.875;
  y := chebyshev_u(15,x);
  f := 2.01035215612500906;
  testrel(18, NE, y, f, cnt,failed);

  x := -0.375;
  y := chebyshev_u(48,x);
  f := 1.078612338977720964;
  testrel(19, NE, y, f, cnt,failed);

  x := 1.375;
  y := chebyshev_u(-48,x);
  f := -0.7778125523172374398e17;
  testrel(20, NE, y, f, cnt,failed);

  y := chebyshev_u(20,1e-7);
  f := 0.9999999999978;
  testrel(21, NE, y, f, cnt,failed);

  y := chebyshev_u(21,1e-10);
  f := 0.2199999999999999998e-8;
  testrel(22, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_chebyshev_v;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chebyshev_v');

  y := chebyshev_v(6,1);
  f := 1;
  testrel(1, 0, y, f, cnt,failed);

  y := chebyshev_v(7,1);
  f := 1;
  testrel(2, 0, y, f, cnt,failed);

  y := chebyshev_v(6,-1);
  f := 13;
  testrel(3, 0, y, f, cnt,failed);

  y := chebyshev_v(7,-1);
  f := -15;
  testrel(4, 0, y, f, cnt,failed);

  y := chebyshev_v(6,0);
  f := -1;
  testrel(5, 0, y, f, cnt,failed);

  y := chebyshev_v(7,0);
  f := 1;
  testrel(6, 0, y, f, cnt,failed);

  y := chebyshev_v(8,0);
  f := 1;
  testrel(7, 0, y, f, cnt,failed);

  y := chebyshev_v(9,0);
  f := -1;
  testrel(8, 0, y, f, cnt,failed);

  y := chebyshev_v(1,-0.75);
  f := -2.5;
  testrel(9, NE, y, f, cnt,failed);

  y := chebyshev_v(3,-1.75);
  f := -47.125;
  testrel(10, NE, y, f, cnt,failed);

  y := chebyshev_v(3,0.75);
  f := -0.875;
  testrel(11, NE, y, f, cnt,failed);

  y := chebyshev_v(20,1e-10);
  f := 1.000000001999999998;
  testrel(12, NE, y, f, cnt,failed);

  y := chebyshev_v(21,1e-10);
  f := -0.4999999988999999989*2;
  testrel(13, NE, y, f, cnt,failed);

  y := chebyshev_v(30,1.75);
  f := 0.953672986334881337e15;
  testrel(14, NE, y, f, cnt,failed);

  y := chebyshev_v(31,-0.75);
  f := 1.979055979754775763;
  testrel(15, NE, y, f, cnt,failed);

  y := chebyshev_v(65,-1.75);
  f := -0.7512911534276803238e33;
  testrel(16, NE, y, f, cnt,failed);

  y := chebyshev_v(66,0.75);
  f := -0.6322890154096164170;
  testrel(17, NE, y, f, cnt,failed);

  y := chebyshev_v(129,-1.75);
  f := -0.1215527424708477387e66;
  testrel(18, NE, y, f, cnt,failed);

  y := chebyshev_v(130,0.75);
  f := 1.066498186308875424;
  testrel(19, NE, y, f, cnt,failed);

  y := chebyshev_v(259,-1.75);
  f := -0.3230037956135866602e131;
  testrel(20, NE, y, f, cnt,failed);

  y := chebyshev_v(262,0.75);
  f := 0.3651824569614297354;
  testrel(21, NE, y, f, cnt,failed);

  y := chebyshev_v(519,-1.75);
  f := -0.2280832294625592145e262;
  testrel(22, NE, y, f, cnt,failed);

  y := chebyshev_v(520,-1.5);
  f := 0.3598599375293868940e218;
  testrel(23, NE, y, f, cnt,failed);

  y := chebyshev_v(520,0.75);
  f := 0.7386964251018360491;
  testrel(24, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-30);
  y := chebyshev_v(250,x);
  f := 1.000058441060756372;
  testrel(25, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-30);
  y := chebyshev_v(250,x);
  f := 0.999941560077637695;
  testrel(26, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-30);
  y := chebyshev_v(20,x);
  f := 1.000000391155506718;
  testrel(27, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-30);
  y := chebyshev_v(20,x);
  f := 0.999999608844544041;
  testrel(28, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-17);
  y := chebyshev_v(17,x);
  f := 0.997666307587320151;
  testrel(29, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-17);
  y := chebyshev_v(17,x);
  f := 1.002335497315851516;
  testrel(30, NE, y, f, cnt,failed);

  {recurrence}
  x := 1-ldexp(1,-20);
  y := chebyshev_v(15,x);
  f := 0.9997711268223232483;
  testrel(31, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-20);
  y := chebyshev_v(15,x);
  f := 1.0002288904944558755;
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_chebyshev_w;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chebyshev_w');

  y := chebyshev_w(6,1);
  f := 13;
  testrel(1, 0, y, f, cnt,failed);

  y := chebyshev_w(7,1);
  f := 15;
  testrel(2, 0, y, f, cnt,failed);

  y := chebyshev_w(6,-1);
  f := 1;
  testrel(3, 0, y, f, cnt,failed);

  y := chebyshev_w(7,-1);
  f := -1;
  testrel(4, 0, y, f, cnt,failed);

  y := chebyshev_w(6,0);
  f := -1;
  testrel(5, 0, y, f, cnt,failed);

  y := chebyshev_w(7,0);
  f := -1;
  testrel(6, 0, y, f, cnt,failed);

  y := chebyshev_w(8,0);
  f := 1;
  testrel(7, 0, y, f, cnt,failed);

  y := chebyshev_w(9,0);
  f := 1;
  testrel(8, 0, y, f, cnt,failed);

  y := chebyshev_w(1,-0.75);
  f := -0.5;
  testrel(9, 0, y, f, cnt,failed);

  y := chebyshev_w(3,-1.75);
  f := -24.625;
  testrel(10, 0, y, f, cnt,failed);

  y := chebyshev_w(3,0.75);
  f := 1.625;
  testrel(11, 0, y, f, cnt,failed);

  y := chebyshev_w(30,1.75);
  f := 0.1826144738103910391e16;
  testrel(12, NE, y, f, cnt,failed);

  y := chebyshev_w(31,-0.75);
  f := 0.7637629988603293896;
  testrel(13, NE, y, f, cnt,failed);

  y := chebyshev_w(65,-1.75);
  f := -0.3923490087868294021e33;
  testrel(14, NE, y, f, cnt,failed);

  y := chebyshev_w(66,0.75);
  f := -2.280674068547797772;
  testrel(15, NE, y, f, cnt,failed);

  y := chebyshev_w(129,-1.75);
  f := -0.6347884945293532023e65;
  testrel(16, NE, y, f, cnt,failed);

  y := chebyshev_w(128,0.75);
  f := -2.7751678941584253546;
  testrel(17, NE, y, f, cnt,failed);

  y := chebyshev_w(259,-1.75);
  f := -0.1686832308156194434e131;
  testrel(18, NE, y, f, cnt,failed);

  y := chebyshev_w(258,0.75);
  f := -2.814866999222654284;
  testrel(19, NE, y, f, cnt,failed);

  y := chebyshev_w(516,-1.75);
  f := 0.3682669992306201899e260;
  testrel(20, NE, y, f, cnt,failed);

  y := chebyshev_w(515,0.75);
  f := 2.709695098411944901;
  testrel(21, NE, y, f, cnt,failed);

  y := chebyshev_w(20,1e-10);
  f := 0.4999999989999999989*2;
  testrel(22, NE, y, f, cnt,failed);

  y := chebyshev_w(21,1e-10);
  f := 1.000000002199999998;
  testrel(23, NE, y, f, cnt,failed);

  {trig/hyp}
  x := 1+ldexp(1,-30);
  y := chebyshev_w(250,x);
  f := 501.0097596191237407;
  testrel(24, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-30);
  y := chebyshev_w(250,x);
  f := 500.9902404949433448;
  testrel(25, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-30);
  y := chebyshev_w(20,x);
  f := 41.00000534579178640;
  testrel(26, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-30);
  y := chebyshev_w(20,x);
  f := 40.99999465420862981;
  testrel(27, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-17);
  y := chebyshev_w(17,x);
  f := 34.97276937799598810;
  testrel(28, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-17);
  y := chebyshev_w(17,x);
  f := 35.02724325632614493;
  testrel(29, NE, y, f, cnt,failed);

  {Recurrence}
  x := 1-ldexp(1,-20);
  y := chebyshev_w(15,x);
  f := 30.99763494137675732;
  testrel(30, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-20);
  y := chebyshev_w(15,x);
  f := 31.00236516598727324;
  testrel(31, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_chebfun;
var
  v,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 10;
  NE2 = 20;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chebyshev_f1');

  {--------------}
  v := -Pi;

  y := chebyshev_f1(v, 2);
  f := 31.32614110218444419;
  testrel(1, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 1);
  f := 1;
  testrel(2, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.875);
  f := -0.1683974373270588857e-1;
  testrel(3, NE2, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.25);
  f := -0.5408145658400643331;
  testrel(4, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.0);
  f := 0.2205840407496980887;
  testrel(5, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, -0.125);
  f := 0.5778910756590656914;
  testrel(6, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, -0.5);
  f := 0.9563500657790059248;
  testrel(7, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, -0.75);
  f := 0.2521620302612723670;
  testrel(8, NE1, y, f, cnt,failed);

  y := chebyshev_f1(v, -1);
  f := -0.9026853619330710662;
  testrel(9, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, -2);
  f := -28.27764901879181877;
  testrel(10, NE, y, f, cnt,failed);

  {--------------}
  v := 7.5;
  y := chebyshev_f1(v, 2);
  f := 9740.395962177307707;
  testrel(11, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.875);
  f := -0.7969235294129488632;
  testrel(12, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.5);
  f := 0;
  testabs(13, 1, y, f, cnt,failed);

  y := chebyshev_f1(v, 0);
  f := 0.7071067811865475244;
  testrel(14, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, -0.75);
  f := 0.7595873626027365985;
  testrel(15, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, -1.5);
  f := 0;
  testabs(16, 1, y, f, cnt,failed);


  {--------------}
  v := 12.345;
  y := chebyshev_f1(v, 10);
  f := 0.5581248293906288603e16;
  testrel(17, NE1, y, f, cnt,failed);

  y := chebyshev_f1(v, 8);
  f := 0.348951097468463844515e15;
  testrel(18, NE1, y, f, cnt,failed);

  y := chebyshev_f1(v, 2);
  f := 5749928.95460841233665731957428;
  testrel(19, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.875);
  f := 0.9990096020328359641;
  testrel(20, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.25);
  f := -0.8450275517082062011;
  testrel(21, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, 0.0);
  f := 0.8567175188650496416;
  testrel(22, NE, y, f, cnt,failed);

  y := chebyshev_f1(v, -0.125);
  f := -0.4954045464085850867;
  testrel(23, NE1, y, f, cnt,failed);

  y := chebyshev_f1(v, -0.5);
  f := 0.7501110696304595415;
  testrel(24, NE1, y, f, cnt,failed);

  y := chebyshev_f1(v, -1);
  f := 0.4679298142605733772;
  testrel(25, NE2, y, f, cnt,failed);

  y := chebyshev_f1(v, -1.5);
  f := 33810.80876745418861;
  testrel(26, NE2, y, f, cnt,failed);

  y := chebyshev_f1(v, -2);
  f := 2690563.187741407234;
  testrel(27, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_gegenbauer_c;
var
  y,f: double;
  cnt,failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gegenbauer_c');

  y := gegenbauer_c(1, -0.2, 1.0);
  f := -0.4;
  testrel( 1, NE, y, f, cnt,failed);

  y := gegenbauer_c(1, 0.0, 1.0);
  f := 2.0;
  testrel( 2, NE, y, f, cnt,failed);

  y := gegenbauer_c(1, 1.0, 1.0);
  f := 2.0;
  testrel( 3, NE, y, f, cnt,failed);

  y := gegenbauer_c(1, 1.0, 0.5);
  f := 1.0;
  testrel( 4, NE, y, f, cnt,failed);

  y := gegenbauer_c(1, 5.0, 1.0);
  f := 10.0;
  testrel( 5, NE, y, f, cnt,failed);

  y := gegenbauer_c(1, 100.0, 0.5);
  f := 100.0;
  testrel( 6, NE, y, f, cnt,failed);

  y := gegenbauer_c(2, -0.2, 0.5);
  f := 0.12;
  testrel( 7, NE, y, f, cnt,failed);

  y := gegenbauer_c(2, 0.0, 1.0);
  f := 1.00;
  testrel( 8, NE, y, f, cnt,failed);

  y := gegenbauer_c(2, 1.0, 1.0);
  f := 3.00;
  testrel( 9, NE, y, f, cnt,failed);

  y := gegenbauer_c(2, 1.0, 0.1);
  f := -0.96;
  testrel(10, NE, y, f, cnt,failed);

  y := gegenbauer_c(2, 5.0, 1.0);
  f := 55.0;
  testrel(11, NE, y, f, cnt,failed);

  y := gegenbauer_c(2, 100.0, 0.5);
  f := 4950.0;
  testrel(12, NE, y, f, cnt,failed);

  y := gegenbauer_c(3, -0.2, 0.5);
  f := 0.112;
  testrel(13, NE, y, f, cnt,failed);

  y := gegenbauer_c(3, 0.0, 1.0);
  f := 2/3; {Fix311}
  testrel(14, NE, y, f, cnt,failed);

  y := gegenbauer_c(3, 1.0, 1.0);
  f := 4.000;
  testrel(15, NE, y, f, cnt,failed);

  y := gegenbauer_c(3, 1.0, 0.1);
  f := -0.392;
  testrel(16, NE, y, f, cnt,failed);

  y := gegenbauer_c(3, 5.0, 1.0);
  f := 220.0;
  testrel(17, NE, y, f, cnt,failed);

  y := gegenbauer_c(3, 100.0, 0.5);
  f := 161600.0;
  testrel(18, NE, y, f, cnt,failed);

  y := gegenbauer_c(10,-0.2, 0.5);
  f := 0.9934389248e-2;
  testrel(19, NE, y, f, cnt,failed);

  y := gegenbauer_c(10, 0.0, 1.0);
  f := 0.2;
  testrel(20, NE, y, f, cnt,failed);

  y := gegenbauer_c(10, 1.0, 1.0);
  f := 11.0;
  testrel(21, NE, y, f, cnt,failed);

  y := gegenbauer_c(10, 1.0, 0.1);
  f := -0.4542309376;
  testrel(22, NE, y, f, cnt,failed);

  y := gegenbauer_c(10, 5.0, 1.0);
  f := 9.23780e+4;
  testrel(23, NE, y, f, cnt,failed);

  y := gegenbauer_c(10, 100.0, 0.5);
  f := 1.5729338392690000e+13;
  testrel(24, NE, y, f, cnt,failed);

  y := gegenbauer_c(21, -0.2, 0.5);
  f := 0.9991386973302926030e-2;
  testrel(25, NE, y, f, cnt,failed);

  y := gegenbauer_c(21, 0.0, 1.0);
  f := 2/21;
  testrel(26, NE, y, f, cnt,failed);

  y := gegenbauer_c(21, 1.0, 1.0);
  f := 22.0;
  testrel(27, NE, y, f, cnt,failed);

  y := gegenbauer_c(21, 1.0, 0.1);
  f := 0.8103854265403986412;
  testrel(28, NE, y, f, cnt,failed);

  y := gegenbauer_c(21, 5.0, 1.0);
  f := 14307150.0;
  testrel(29, NE, y, f, cnt,failed);

  y := gegenbauer_c(21, 100.0, 0.5);
  f := 0.826932426083121018e20;
  testrel(30, NE, y, f, cnt,failed);

  y := gegenbauer_c(42, -0.2, 0.5);
  f := -0.4324011590989000524e-2;
  testrel(31, NE, y, f, cnt,failed);

  y := gegenbauer_c(42, 0.0, 1.0);
  f := 1/21;
  testrel(32, NE, y, f, cnt,failed);

  y := gegenbauer_c(42, 1.0, 1.0);
  f := 43.0;
  testrel(33, NE, y, f, cnt,failed);

  y := gegenbauer_c(42, 1.0, 0.1);
  f := 0.3961791364804087116;
  testrel(34, 3*NE, y, f, cnt,failed);    {!!!!!!!}

  y := gegenbauer_c(42, 5.0, 1.0);
  f := 3042312350.0;
  testrel(35, NE, y, f, cnt,failed);

  y := gegenbauer_c(42, 100.0, 0.5);
  f := -0.2728803040216282769e30;
  testrel(36, NE, y, f, cnt,failed);

  y := gegenbauer_c(1000, 100.0, 1.0);
  f := 3.335366613562732234e232;
  testrel(37, 60, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_hermite_h;
var
  y,f: double;
  cnt,failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hermite_h');

  y := hermite_h(11,-0.5);
  f := 107029.0;
  testrel(1, NE, y, f, cnt,failed);

  y := hermite_h(11,0);
  f := 0.0;
  testrel(2, NE, y, f, cnt,failed);

  y := hermite_h(11,0.1);
  f := -64328.09016829952;
  testrel(3, NE, y, f, cnt,failed);

  y := hermite_h(11,0.5);
  f := -107029.0;
  testrel(4, NE, y, f, cnt,failed);

  y := hermite_h(11,1);
  f := 230848.0;
  testrel(5, NE, y, f, cnt,failed);

  y := hermite_h(11,3);
  f := -10425024.0;
  testrel(6, NE, y, f, cnt,failed);

  y := hermite_h(11,5);
  f := 24329873600.0;
  testrel(7, NE, y, f, cnt,failed);

  y := hermite_h(11,10);
  f := 153373602947200.0;
  testrel(8, NE, y, f, cnt,failed);

  y := hermite_h(12,-0.5);
  f := -604031.0;
  testrel(9, NE, y, f, cnt,failed);

  y := hermite_h(12,0);
  f := 665280.0;
  testrel(10, NE, y, f, cnt,failed);

  y := hermite_h(12,0.1);
  f := 586769.878872887296;
  testrel(11, NE, y, f, cnt,failed);

  y := hermite_h(12,0.5);
  f := -604031.0;
  testrel(12, NE, y, f, cnt,failed);

  y := hermite_h(12,1);
  f := 280768.0;
  testrel(13, NE, y, f, cnt,failed);

  y := hermite_h(12,3);
  f := 5517504.0;
  testrel(14, NE, y, f, cnt,failed);

  y := hermite_h(12,5);
  f := 171237081280.0;
  testrel(15, NE, y, f, cnt,failed);

  y := hermite_h(12,10);
  f := 2889419938329280.0;
  testrel(16, NE, y, f, cnt,failed);

  y := hermite_h(101,0);
  f := 0;
  testrel(17, NE, y, f, cnt,failed);

  y := hermite_h(101,0.1);
  f := 0.4325734655712572854e95;
  testrel(18, NE, y, f, cnt,failed);

  y := hermite_h(101,1);
  f := 0.7146283559232673061e95;
  testrel(19, NE, y, f, cnt,failed);

  y := hermite_h(101,10);
  f := -0.2108624426594399435e117;
  testrel(20, NE, y, f, cnt,failed);

  y := hermite_h(101,101);
  f := 0.5400837669264194369e233;
  testrel(21, NE, y, f, cnt,failed);

  y := hermite_h(200,0);
  f := 0.8450550186924629496e217;
  testrel(22, NE, y, f, cnt,failed);

  y := hermite_h(200,0.1);
  f := -0.3553562424339073720e217;
  testrel(23, 2, y, f, cnt,failed);         {!!!!!!!!}

  y := hermite_h(200,1.5);
  f := 0.4262910729228930473e217;
  testrel(24, NE, y, f, cnt,failed);

  y := hermite_h(200,1.625);
  f := 0.1468781261584823922e218;
  testrel(25, NE, y, f, cnt,failed);

  y := hermite_h(200,1.75);
  f := -0.3539419442763632425e218;
  testrel(26, NE, y, f, cnt,failed);

  y := hermite_h(200,10.0);
  f := -0.4700538718627334916e239;
  testrel(27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_hermite_he;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
  NE1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hermite_he');

  y := hermite_he(11,-0.5);
  f := 3240.81005859375;
  testrel(1, NE, y, f, cnt,failed);

  y := hermite_he(11,0);
  f := 0.0;
  testrel(2, NE, y, f, cnt,failed);

  y := hermite_he(11,0.1);
  f := -1022.24420105499;
  testrel(3, NE, y, f, cnt,failed);

  y := hermite_he(11,0.5);
  f := -3240.81005859375;
  testrel(4, NE, y, f, cnt,failed);

  y := hermite_he(11,1);
  f := 936.0;
  testrel(5, NE, y, f, cnt,failed);

  y := hermite_he(11,3);
  f := 12312.0;
  testrel(6, NE, y, f, cnt,failed);

  y := hermite_he(11,5);
  f := -792600.0;
  testrel(7, NE, y, f, cnt,failed);

  y := hermite_he(11,10);
  f := 54224221050.0;
  testrel(8, NE, y, f, cnt,failed);

  y := hermite_he(12,-0.5);
  f := -2159.888427734375;
  testrel(9, NE, y, f, cnt,failed);

  y := hermite_he(12,0);
  f := 10395.0;
  testrel(10, NE, y, f, cnt,failed);

  y := hermite_he(12,0.1);
  f := 9776.483654843401;
  testrel(11, NE, y, f, cnt,failed);

  y := hermite_he(12,0.5);
  f := -2159.888427734375;
  testrel(12, NE, y, f, cnt,failed);

  y := hermite_he(12,1);
  f := -12440.0;
  testrel(13, NE, y, f, cnt,failed);

  y := hermite_he(12,3);
  f := -67608.0;
  testrel(14, NE, y, f, cnt,failed);

  y := hermite_he(12,5);
  f := -5939480.0;
  testrel(15, NE, y, f, cnt,failed);

  y := hermite_he(12,10);
  f := 475153523395.0;
  testrel(16, NE, y, f, cnt,failed);

  y := hermite_he(101,0);
  f := 0;
  testrel(17, NE, y, f, cnt,failed);

  y := hermite_he(101,0.1);
  f := 0.2315844485708444393e80;
  testrel(18, NE, y, f, cnt,failed);

  y := hermite_he(101,1);
  f := -0.2112875427125553680e80;
  testrel(19, NE, y, f, cnt,failed);

  y := hermite_he(101,10);
  f := 0.1711488399736831495e91;
  testrel(20, NE, y, f, cnt,failed);

  y := hermite_he(101,101);
  f := 0.1657027632229917093e203;
  testrel(21, NE, y, f, cnt,failed);

  y := hermite_he(200,0);
  f := 0.6666308670072953744e187;
  testrel(22, NE, y, f, cnt,failed);

  y := hermite_he(200,0.1);
  f := 0.1030515302661571644e187;
  testrel(23, NE1, y, f, cnt,failed);

  y := hermite_he(200,1.5);
  f := -0.8475471777212029840e187;
  testrel(24, NE, y, f, cnt,failed);

  y := hermite_he(200,1.75);
  f := 0.1338360150483440600e188;
  testrel(25, NE, y, f, cnt,failed);

  y := hermite_he(200,10.0);
  f := 0.4630390274429460851e198;
  testrel(26, NE, y, f, cnt,failed);

  y := hermite_he(16,1/3); {Fix311}
  f := 448614.8690102551597;
  testrel(27, NE1, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_jacobi_p;
var
  y,f: double;
  cnt,failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_p');

  y := jacobi_p(2,1,0.75,0.5);
  f := 0.291015625;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_p(2,1,0.75,0.1);
  f := -0.647109375;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_p(8,1,0.75,0.1);
  f := 0.2039343077673979243;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_p(8,1,0.75,2.5);
  f := 178708.9906612954856;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_p(8, 1.5, -0.5, -1);
  f := 0.196380615234375;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_p(8,1.5, -0.5, -2);
  f := 8899.118499755859377;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_p(8,1.5, -0.5, 0);
  f := 0.21820068359375e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_p(8,1.5, -0.5, 0.9375);
  f := 6.162780192089485354;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_p(8,1.5, -0.5, -0.9375);
  f := -0.2024995914052851731;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_p(8,1.5, -0.5, 0.25);
  f := -0.4237014055252075195;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_p(8,1.5, -0.5, -0.75);
  f := 0.2199906110763549806;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_p(25,1.5, -0.5, 0.75);
  f := -0.8993370671575232300;
  testrel(12, NE, y, f, cnt,failed);

  y := jacobi_p(100,1.5, -0.5, -0.25);
  f := 0.3442708728724186509e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_p(100,1.5, -0.5, -2);
  f := 0.1094511122764460740e57;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_p(200,1.5, -0.5, 0.25);
  f := -0.5283109520288860554e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := jacobi_p(200,1.5, -0.5, 2.5);
  f := 0.1562423735191907416e136;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_p(8, 0.5, 1.5, 2);
  f := 21852.069671630859375;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_p(8, 0.5, 1.5, -0.875);
  f := -0.8486263942904770374;
  testrel(18, NE, y, f, cnt,failed);

  y := jacobi_p(100, 9.5, -1.5, 2);
  f := 0.4661975827107720642e60;
  testrel(19, NE, y, f, cnt,failed);

  y := jacobi_p(100, 9.5, -1, 2);
  f := 0.5731010730126711304e60;
  testrel(20, NE, y, f, cnt,failed);

  y := jacobi_p(10, 9.5, -1., 0.5);
  f := -0.1811809968133711664e3;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_laguerre;
var
  y,f: double;
  cnt,failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','laguerre/_l/_ass');

  y := laguerre(1, -2.0, 1.0);
  f := -2.0;
  testrel( 1, NE, y, f, cnt,failed);

  y := laguerre(1, 0.5, 1.0);
  f := 0.5;
  testrel( 2, NE, y, f, cnt,failed);

  y := laguerre(2, -2.0, 1.0);
  f := 0.5;
  testrel( 3, NE, y, f, cnt,failed);

  y := laguerre(2, 0.5, -1.0);
  f := 4.875;
  testrel( 4, NE, y, f, cnt,failed);

  y := laguerre(2, 0.5,  1.0);
  f := -0.125;
  testrel( 5, NE, y, f, cnt,failed);

  y := laguerre_ass(2, 1,  1.0);
  f := 0.5;
  testrel( 6, NE, y, f, cnt,failed);

  y := laguerre(2,-1.0,  1.0);
  f := -0.5;
  testrel( 7, NE, y, f, cnt,failed);

  y := laguerre(2,-2.0,  1.0);
  f := 0.5;
  testrel( 8, NE, y, f, cnt,failed);

  y := laguerre(2,-3.0,  1.0);
  f := 2.5;
  testrel( 9, NE, y, f, cnt,failed);

  y := laguerre(3, -2.0, 1.0);
  f := 1/3;  {Fix311}
  testrel(10, NE, y, f, cnt,failed);

  y := laguerre(3, 0.5, -1.0);
  f := 8.479166666666666667;
  testrel(11, NE, y, f, cnt,failed);

  y := laguerre(3, 0.5,  1.0);
  f := -0.6041666666666666667;
  testrel(12, NE, y, f, cnt,failed);

  y := laguerre_ass(3, 1,  1.0);
  f := -1/6;
  testrel(13, NE, y, f, cnt,failed);

  y := laguerre_ass(3, 2,  1.0);
  f := 7/3;
  testrel(14, NE, y, f, cnt,failed);

  y := laguerre_ass(3, 2, 0);
  f := 10.0;
  testrel(15, NE, y, f, cnt,failed);

  y := laguerre_ass(3,2,0.0009765625);
  f := 64361610239.0/6442450944.0;
  testrel(16, NE, y, f, cnt,failed);

  y := laguerre_ass(3,2,0.125);
  f := 26999/3072;
  testrel(17, NE, y, f, cnt,failed);

  y := laguerre_ass(3,2,1.5);
  f := 1/16;
  testrel(18, NE, y, f, cnt,failed);

  y := laguerre_ass(3,2,100);
  f := -427970/3;
  testrel(19, NE, y, f, cnt,failed);

  y := laguerre_ass(3,2,-100);
  f := 578030/3;
  testrel(20, NE, y, f, cnt,failed);

  y := laguerre(3,-2.0,  1.0);
  f := 1/3; {Fix311}
  testrel(21, NE, y, f, cnt,failed);

  y := laguerre(3,-2,0.125);
  f := 0.7486979166666666667e-2;
  testrel(22, NE, y, f, cnt,failed);

  y := laguerre(3,-2,0.0009765625);
  f := 0.4766819377740224203e-6;
  testrel(23, NE, y, f, cnt,failed);

  y := laguerre(3,-3.0,  1.0);
  f := -1/6;
  testrel(24, NE, y, f, cnt,failed);

  y := laguerre(3,-4.0,  1.0);
  f := -8/3;
  testrel(25, NE, y, f, cnt,failed);

  y := laguerre_ass(4, 2, 0.5);
  f := 6.752604166666666667;
  testrel(26, NE, y, f, cnt,failed);

  y := laguerre(10,-0.5,0);
  f := 46189/262144;   {Fix311}
  testrel(27, NE, y, f, cnt,failed);

  y := laguerre(10,-0.5,0.125);
  f := -0.1198039588943549134;
  testrel(28, NE, y, f, cnt,failed);

  y := laguerre(10,-0.5,-0.125);
  f := 0.8053223525137110097;
  testrel(29, NE, y, f, cnt,failed);

  y := laguerre(10,-0.5,1.5);
  f := 0.2206821986607142857e-1;
  testrel(30, NE, y, f, cnt,failed);

  y := laguerre(10,-0.5,-1.5);
  f := 0.1102325383649553571e3;
  testrel(31, NE, y, f, cnt,failed);

  y := laguerre(10,-0.5,100);
  f := 0.9602106167023362142e13;
  testrel(32, NE, y, f, cnt,failed);

  y := laguerre(10,-0.5,-100);
  f := 0.6599656255927201471e14;
  testrel(33, NE, y, f, cnt,failed);

  y := laguerre(42,-1.5,0);
  f := -0.1045754671992811442e-2;
  testrel(34, NE, y, f, cnt,failed);

  y := laguerre(42,-1.5,1);
  f := -0.9140042224300855665e-2;
  testrel(35, NE, y, f, cnt,failed);

  y := laguerre(42,-1.5,10);
  f := -0.2285000442561404022e1;
  testrel(36, NE, y, f, cnt,failed);

  y := laguerre(42,-1.5,100);
  f := -0.870027424548133248e21;
  testrel(37, NE, y, f, cnt,failed);

  y := laguerre_ass(90, 2,  0.5);
  f := -48.79047157201507897;
  testrel(38, NE, y, f, cnt,failed);

  y := laguerre_ass(90, 2, -100.0);
  f := 2.5295879275042410902e+63;
  testrel(39, NE, y, f, cnt,failed);

  y := laguerre_ass(90, 2,  100.0);
  f := -2.0929042259546928670e+20;
  testrel(40, NE, y, f, cnt,failed);

  y := laguerre(99,1.5,0);
  f := 0.7550696187240360581e3;
  testrel(41, NE, y, f, cnt,failed);

  y := laguerre(99,1.5,1);
  f := -0.323309770429035454e1;
  testrel(42, NE, y, f, cnt,failed);

  y := laguerre(99,1.5,10);
  f := -0.818754021133035264e2;
  testrel(43, NE, y, f, cnt,failed);

  y := laguerre(99,1.5,100);
  f := 0.3142897935722162957e21;
  testrel(44, NE, y, f, cnt,failed);

  y := laguerre_l(100, 0.5);
  f := 0.1868226036769227880;
  testrel(45, NE, y, f, cnt,failed);

  y := laguerre_l(100, -0.5);
  f := 0.1199091342685978022e6;
  testrel(46, NE, y, f, cnt,failed);

  y := laguerre_l(100, 10.5);
  f := 0.9179690735405005987e1;
  testrel(47, NE, y, f, cnt,failed);

  y := laguerre_l(100, -10.5);
  f := 0.56329215744170606488e25;
  testrel(48, NE, y, f, cnt,failed);

  y := laguerre_l(100, 100.5);
  f := -3.9844782875811907525e20;
  testrel(49, NE, y, f, cnt,failed);

  y := laguerre_l(100, -100.5);
  f := 0.25531876914022721486e68;
  testrel(50, NE, y, f, cnt,failed);

  y := laguerre_l(100, 150);
  f := -1.4463204337261709595e31;
  testrel(51, NE, y, f, cnt,failed);

  y := laguerre_ass(100, 2, -0.5);
  f := 2.2521795545919391405e+07;
  testrel(52, NE, y, f, cnt,failed);

  y := laguerre_ass(100, 2,  0.5);
  f := -28.764832945909097418;
  testrel(53, NE, y, f, cnt,failed);

  y := laguerre_ass(1000, 2, -0.5);
  f := 2.4399915170947549589e+21;
  testrel(54, NE, y, f, cnt,failed);

  y := laguerre_l(15, 1.5);
  f := -0.445840471006807594304760529713;
  testrel(55, NE, y, f, cnt,failed);

  y := laguerre_l(15, -1.5);
  f := 0.981192886508793791401687389063e3;
  testrel(56, NE, y, f, cnt,failed);

  y := laguerre_l(15, 10.5);
  f := 0.1300848998942491891500823465e2;
  testrel(57, NE, y, f, cnt,failed);

  y := laguerre_l(15, -10.5);
  f := 0.976246944859354857723196069678e8;
  testrel(58, NE, y, f, cnt,failed);

  y := laguerre_l(15, 100.5);
  f := -0.56894945639281351968637746151e17;
  testrel(59, NE, y, f, cnt,failed);

  y := laguerre_l(15, -100.5);
  f := 0.594126093839947334749590608537e19;
  testrel(60, NE, y, f, cnt,failed);

  y := laguerre_l(15, 150);
  f := -0.62788150869447747429712145e20;
  testrel(61, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_orthopoly;
var
  x,y,f,t,c: double;
  failed: longint;
  i,n: integer;
begin
  writeln('Function: ','Orthogonal polynomials interrelations');

  {Check absolute values because recurrence relations do not always}
  {give high relative accuraca near the zeroes of the polynomials.}
  writeln(' - Jacobi / Legendre');
  failed := 0;
  for i:=1 to 200 do begin
    x := -1;
    while x<0.95 do begin
      {HMF[1], 22.5.35}
      y := jacobi_p(i,0,0,x);
      f := legendre_p(i,x);
      t := abs(y-f)/eps_d;
      if t > 1 then begin
        writeln('   i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;
      x := x + 0.1;
    end;
  end;
  inc(total_cnt);
  if failed > 0 then inc(total_failed);

  failed := 0;
  writeln(' - Jacobi / Chebyshev T');
  for i:=1 to 100 do begin
    x := -0.9;
    c := sqrtpi/Gamma(i+0.5)*fac(i);
    while x<0.95 do begin
      {HMF[1], 22.5.31}
      y := jacobi_p(i,-0.5,-0.5,x)*c;
      f := chebyshev_t(i,x);
      t := abs(y-f)/eps_d;
      if t > 1 then begin
        writeln('   i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;
      x := x + 0.1;
    end;
  end;
  inc(total_cnt);
  if failed > 0 then inc(total_failed);

  failed := 0;
  writeln(' - Jacobi / Chebyshev U');
  for i:=1 to 100 do begin
    x := -0.9;
    c := 0.5*sqrtpi/Gamma(i+1.5)*fac(i+1);
    while x<0.95 do begin
      {HMF[1], 22.5.32}
      y := jacobi_p(i,0.5,0.5,x)*c;
      f := chebyshev_u(i,x);
      t := abs(y-f)/eps_d;
      if t > 70 then begin
        writeln('   i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;
      x := x + 0.1;
    end;
  end;
  inc(total_cnt);
  if failed > 0 then inc(total_failed);

  failed := 0;
  writeln(' - Jacobi / Chebyshev V');
  for i:=1 to 100 do begin
    x := -0.9;
    c := fac(i)/pochhammer(0.5,i);
    while x<0.95 do begin
      y := jacobi_p(i,-0.5,0.5,x)*c;
      f := chebyshev_v(i,x);
      t := abs(y-f)/eps_d;
      if t > 10 then begin
        writeln('   i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;
      x := x + 0.1;
    end;
  end;
  inc(total_cnt);
  if failed > 0 then inc(total_failed);

  failed := 0;
  writeln(' - Jacobi / Chebyshev W');
  for i:=1 to 100 do begin
    x := -0.9;
    c := (2*i+1)*fac(i)/pochhammer(1.5,i);
    while x<0.95 do begin
      y := jacobi_p(i,0.5,-0.5,x)*c;
      f := chebyshev_w(i,x);
      t := abs(y-f)/eps_d;
      if t > 10 then begin
        writeln('   i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;
      x := x + 0.1;
    end;
  end;
  inc(total_cnt);
  if failed > 0 then inc(total_failed);


  failed := 0;
  writeln(' - Gegenbauer / Legendre / Chebyshev U');
  for i:=1 to 100 do begin
    x := -0.9;
    while x<0.95 do begin
      {HMF[1], 22.5.32}
      y := gegenbauer_c(i,0.5,x);
      f := legendre_p(i,x);
      t := abs(y-f)/eps_d;
      if t > 1 then begin
        writeln('  Pn: i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;

      y := gegenbauer_c(i,1,x);
      f := chebyshev_u(i,x);
      t := abs(y-f)/eps_d;
      if t > 70 then begin
        writeln('  Un: i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;
      x := x + 0.1;
    end;
  end;

  inc(total_cnt);
  if failed > 0 then inc(total_failed);

  writeln(' - Laguerre / Hermite');
  failed := 0;
  for i:=1 to 25 do begin
    x := 0.0;
    while x<=100 do begin
      n := i div 2;
      if odd(i) then begin
        y := x*laguerre(n,0.5,x*x)*fac(n);
        if odd(n) then y := -y;
        y := ldexp(y,i);
      end
      else begin
        y := laguerre(n,-0.5,x*x)*fac(n);
        if odd(n) then y := -y;
        y := ldexp(y,i);
      end;
      f := hermite_h(i,x);
      if f=0 then t := abs(y-f)/eps_d
      else t := abs(1-y/f)/eps_d;
      if t > 2 then begin {2 needed only for FPC240}
        writeln('   i/x/t = ',i:6, x:10:4, t:30:2);
        inc(failed);
      end;
      x := x + 0.875;
    end;
  end;
  inc(total_cnt);
  if failed > 0 then inc(total_failed);

  if failed>0 then writeln(' *** failed!')
  else writeln(' Tests OK');
end;


{---------------------------------------------------------------------------}
procedure test_zernike_r;
var
  y,f: double;
  cnt,failed: integer;
const
  NE = 1;
  NA = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','zernike_r');

  {Wolfram Alpha: ZernikeR[3,1,0.5] etc}

  y := zernike_r(3,1,0.5);
  f := -0.625;
  testrel( 1, NE, y, f, cnt,failed);

  y := zernike_r(6,2,0.75);
  f := -0.283447265625;
  testrel( 2, NE, y, f, cnt,failed);

  y := zernike_r(19,5,0.25);
  f := -0.2671170888206688687;
  testrel( 3, NE, y, f, cnt,failed);

  y := zernike_r(19,5,0.4287525723551651222);
  f := 0;
  testabs( 4, NA, y, f, cnt,failed);

  y := zernike_r(19,5,0.5);
  f := 0.28446197509765625;
  testrel( 5, NE, y, f, cnt,failed);

  y := zernike_r(19,5,2.5);
  f := 8.6222178509933948516845703125e11;
  testrel( 6, NE, y, f, cnt,failed);

  y := zernike_r(14,4,0.75);
  f := -0.298682145774364471435546875;
  testrel( 7, NE, y, f, cnt,failed);

  y := zernike_r(14,4,0.8246570394661101917);
  f := 0;
  testabs( 8, NA, y, f, cnt,failed);

  y := zernike_r(14,4,0.875);
  f := 0.3172262174016395875;
  testrel( 9, NE, y, f, cnt,failed);

  y := zernike_r(42,16,0.75);
  f := -0.1378478860818820635;
  testrel(10, NE, y, f, cnt,failed);

  y := zernike_r(42,16,0.8597234728144440095);
  f := 0;
  testabs(11, NA, y, f, cnt,failed);

  y := zernike_r(42,16,0.875);
  f := 0.1826304987525820983;
  testrel(12, NE, y, f, cnt,failed);

  y := zernike_r(42,16,0.984375);
  f := 0.3032452839343570153;
  testrel(13, NE, y, f, cnt,failed);

  y := zernike_r(42,16,2.5);
  f := 0.2842735123784586302e27;
  testrel(14, NE, y, f, cnt,failed);

  y := zernike_r(120,42,0.125);
  f := -0.7251933738244693132e-15;
  testrel(15, NE, y, f, cnt,failed);

  y := zernike_r(120,42,0.25);
  f := -0.6271346168560805383e-4;
  testrel(16, NE, y, f, cnt,failed);

  y := zernike_r(120,42,0.5);
  f := 0.9992364238854123563e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := zernike_r(120,42,0.75);
  f := -0.1509886728548490631e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := zernike_r(120,42,0.875);
  f := 0.1087275784879156651;
  testrel(19, NE, y, f, cnt,failed);

  y := zernike_r(120,42,7.5);
  f := 0.3591506386963877752e137;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.
