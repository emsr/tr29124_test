{Part 8 of regression test for SPECFUND unit  (c) 2012-2018  W.Ehrhardt}

unit t_sfd8b;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_fermi_dirac;
procedure test_fermi_dirac_r;
procedure test_fermi_dirac_half;
procedure test_bose_einstein;

implementation

uses
  specfund, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_fermi_dirac;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 4;
  NE1 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','fermi_dirac');

  {Maple: fd := (n,x) -> -polylog(n+1,-exp(x));}

  y := fermi_dirac(0,1);
  f := 1.313261687518222834;
  testrel(1, NE, y, f, cnt,failed);

  y := fermi_dirac(0,-1);
  f := 0.3132616875182228340;
  testrel(2, NE, y, f, cnt,failed);

  y := fermi_dirac(1,1);
  f := 1.806286070444774257;
  testrel(3, NE, y, f, cnt,failed);

  y := fermi_dirac(-1,1);
  f := 0.7310585786300048793;
  testrel(4, NE, y, f, cnt,failed);

  y := fermi_dirac(1,-1);
  f := 0.3386479964034521798;
  testrel(5, NE, y, f, cnt,failed);

  y := fermi_dirac(-1,-1);
  f := 0.2689414213699951207;
  testrel(6, NE, y, f, cnt,failed);

  y := fermi_dirac(-2,1);
  f := 0.1966119332414818525;
  testrel(7, NE, y, f, cnt,failed);

  y := fermi_dirac(-3,1);
  f := -0.9085774767294840944e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := fermi_dirac(-2,-1);
  f := 0.1966119332414818525;
  testrel(9, NE, y, f, cnt,failed);

  y := fermi_dirac(-3,-1);
  f := 0.9085774767294840944e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := fermi_dirac(-4,1e-5);
  f := -0.1249999999875000000;
  testrel(11, NE, y, f, cnt,failed);

  y := fermi_dirac(-4,-1e-5);
  f := -0.1249999999875000000;
  testrel(12, NE, y, f, cnt,failed);

  y := fermi_dirac(4,1e-5);
  f := 0.9721292408202815493;
  testrel(13, NE, y, f, cnt,failed);

  y := fermi_dirac(4,-1e-5);
  f := 0.9721103001636913303;
  testrel(14, NE, y, f, cnt,failed);

  y := fermi_dirac(8,10);
  f := 7946.319776573429079;
  testrel(15, NE, y, f, cnt,failed);

  y := fermi_dirac(8,-10);
  f := 0.4539992573679893686e-4;
  testrel(16, NE, y, f, cnt,failed);

  y := fermi_dirac(9,10);
  f := 10388.99016731291248;
  testrel(17, NE, y, f, cnt,failed);

  y := fermi_dirac(9,-10);
  f := 0.4539992774964110184e-4;
  testrel(18, NE, y, f, cnt,failed);

  y := fermi_dirac(-3,10);
  f := -0.4539168599011319565e-4;
  testrel(19, NE+1, y, f, cnt,failed);  {+1 for ARM}

  y := fermi_dirac(-3,-10);
  f := 0.4539168599011319565e-4;
  testrel(20, NE+1, y, f, cnt,failed);  {+1 for ARM}

  y := fermi_dirac(8,10000);
  f := 0.2755735186158236597e31;
  testrel(21, NE, y, f, cnt,failed);

  y := fermi_dirac(8,15000);
  f := 0.1059396483983069418e33;
  testrel(22, NE, y, f, cnt,failed);

  y := fermi_dirac(5,20000);
  f := 0.8888889985511638002e23;
  testrel(23, NE, y, f, cnt,failed);

  y := fermi_dirac(1,12000);
  f := 0.7200000164493406684823e8;
  testrel(24, NE, y, f, cnt,failed);

  y := fermi_dirac(0,12000);
  f := 12000;
  testrel(25, NE, y, f, cnt,failed);

  y := fermi_dirac(-1,12000);
  f := 1;
  testrel(26, NE, y, f, cnt,failed);

  y := fermi_dirac(-2,12000);
  f := 0;
  testabs(27, 0, y, f, cnt,failed);

  y := fermi_dirac(100,200);
  f := 0.4216103500655924161e73;
  testrel(28, NE1, y, f, cnt,failed);

  y := fermi_dirac(200,100);
  f := 0.2688117141816135437e44;
  testrel(29, NE1, y, f, cnt,failed);

  y := fermi_dirac(-100,200);
  f := 0.1383896526736737531e-86;
  testrel(30, NE, y, f, cnt,failed);

  y := fermi_dirac(-200,100);
  f := -0.1111917984507457867e-26;
  testrel(31, NE, y, f, cnt,failed);

  y := fermi_dirac(500,10);
  f := 22026.46579480671652;
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_fermi_dirac_half;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 4;
  NE2 = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','fermi_dirac_m05/p05/p15/p25');

  y := fermi_dirac_m05(-20);
  f := 0.2061153619434517730e-8;
  testrel(1, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(-10);
  f := 0.4539847236080549532e-4;
  testrel(2, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(-2);
  f := 0.1236656218012099427;
  testrel(3, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(-0.5);
  f := 0.4312314419263970869;
  testrel(4, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(0);
  f := 0.6048986434216303702;
  testrel(5, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(1);
  f := 1.027057125474350689;
  testrel(6, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(3.5);
  f := 2.026193699100460058;
  testrel(7, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(4);
  f := 2.185869690623811946;
  testrel(8, NE2, y, f, cnt,failed);

  y := fermi_dirac_m05(4.5);
  f := 2.334245766142244249;
  testrel(9, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(10);
  f := 3.552779239536617160;
  testrel(10, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(13.5);
  f := 4.136326785587186739;
  testrel(11, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(14);
  f := 4.212933653105178558;
  testrel(12, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(14.5);
  f := 4.288146102878802782;
  testrel(13, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(20);
  f := 5.041018507535328603;
  testrel(14, NE2, y, f, cnt,failed);

  y := fermi_dirac_m05(39.5);
  f := 7.089878704263408613;
  testrel(15, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(40);
  f := 7.134657233550764731;
  testrel(16, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(40.5);
  f := 7.179155892384824069;
  testrel(17, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(100);
  f := 11.28332744292768060;
  testrel(18, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(1000);
  f := 35.68246764915936962; {alpha}
  testrel(19, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(10000);
  f := 112.8379162455239043; {alpha}
  testrel(20, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(1000000);
  f := 1128.379167095048547; {alpha}
  testrel(21, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(-20);
  f := 0.2061153620936537778e-8;
  testrel(22, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(-10);
  f := 0.4539920105264132755e-4;
  testrel(23, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(-2);
  f := 0.1292985133200755911;
  testrel(24, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(-0.5);
  f := 0.5075371035546378361;
  testrel(25, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(0);
  f := 0.7651470246254079454;
  testrel(26, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(1);
  f := 1.575640776151300231;
  testrel(27, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(3.5);
  f := 5.458044388745459911;
  testrel(28, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(4);
  f := 6.511567592754791415;
  testrel(29, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(4.5);
  f := 7.642031275793455330;
  testrel(30, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(9);
  f := 20.62401821156795236;
  testrel(31, NE2, y, f, cnt,failed);

  y := fermi_dirac_p05(11.5);
  f := 29.61226879150170305;
  testrel(32, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(12);
  f := 31.54020328704424259;
  testrel(33, NE2, y, f, cnt,failed);

  y := fermi_dirac_p05(12.5);
  f := 33.50923983932172004;
  testrel(34, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(20);
  f := 67.49151222165892049;
  testrel(35, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(35.5);
  f := 159.2691126191554513;
  testrel(36, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(36);
  f := 162.6413796464579914;
  testrel(37, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(36.5);
  f := 166.0371774151569984;
  testrel(38, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(100);
  f := 752.3455915521961188;
  testrel(39, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(10000);
  f := 752252.7873442217908;  {alpha}
  testrel(40, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(1000000);
  f := 7.522527780646031039e8; {alpha}
  testrel(41, NE, y, f, cnt,failed);

  {---------}
  y := fermi_dirac_p15(-20);
  f := 0.2061153621687547803e-8;
  testrel(42, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(-10);
  f := 0.4539956540456176333e-4;
  testrel(43, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(-2);
  f := 0.1322467822517723668;
  testrel(44, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(-0.5);
  f := 0.5526495259473540854;
  testrel(45, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(0);
  f := 0.8671998890121841382;
  testrel(46, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(1);
  f := 2.002258148778464457;
  testrel(47, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(3.5);
  f := 10.27141115061732259;
  testrel(48, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(4);
  f := 13.26048817729106850;
  testrel(49, NE+1, y, f, cnt,failed);

  y := fermi_dirac_p15(4.5);
  f := 16.79579730928991532;
  testrel(50, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(9);
  f := 78.66626344736753771;
  testrel(51, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(11.5);
  f := 141.2287997372907052;
  testrel(52, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(12);
  f := 156.5151864279572804;
  testrel(53, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(12.5);
  f := 172.7758530334651074;
  testrel(54, NE2, y, f, cnt,failed);

  y := fermi_dirac_p15(20);
  f := 546.5630100657601959;
  testrel(55, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(35.5);
  f := 2270.464576509283294;
  testrel(56, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(36);
  f := 2350.941215700957481;
  testrel(57, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(36.5);
  f := 2433.109877926863218;
  testrel(58, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(100);
  f := 30108.67168135486936;
  testrel(59, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(10000);
  f := 3.009011297865632890e9; {alpha}
  testrel(60, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(20000);
  f := 1.702153755962129317e10; {alpha}
  testrel(61, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(1000000);
  f := 3.00901111227325e14;  {GSL/double}
  testrele(62, 5e-15, y, f, cnt,failed);

  {---------}
  y := fermi_dirac_p25(-20);
  f := 0.2061153622063052815e-8;
  testrel(63, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(-10);
  f := 0.4539974758252285430e-4;
  testrel(64, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(-2);
  f := 0.1337669290459732776;
  testrel(65, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(-0.5);
  f := 0.5779521605410086652;
  testrel(66, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(0);
  f := 0.9275535777739480351;
  testrel(67, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(1);
  f := 2.294832963121518129;
  testrel(68, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(3.5);
  f := 15.60768391014749828;
  testrel(69, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(4);
  f := 21.46870822801142616;
  testrel(70, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(4.5);
  f := 28.959226461407335057;
  testrel(71, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(9);
  f := 221.7903608693186891;
  testrel(72, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(11.5);
  f := 491.9765499165918583;
  testrel(73, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(12);
  f := 566.3723808356910136;
  testrel(74, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(12.5);
  f := 648.654118806843732657;
  testrel(75, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(20);
  f := 3186.735096831265885;
  testrel(76, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(35.5);
  f := 23178.76346773165267;
  testrel(77, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(36);
  f := 24334.04466016358338;
  testrel(78, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(36.5);
  f := 25529.98668772732444;
  testrel(79, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(100);
  f := 860954.9737352777100;
  testrel(80, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(1000);
  f := 0.2718704450106142809e10;
  testrel(81, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(5000);
  f := 0.7598904953966710865e12;
  testrel(82, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_fermi_dirac_r;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 2;
  NE2 = 5;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','fermi_dirac_r');

  {Note: more test in the wrapper functions for s=-1/2 .. 5/2}
  y := fermi_dirac_r(-0.25, 0);
  f := 0.6511156799649282540;
  testrel(1, NE, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, 1);
  f := 1.172413596752283625;
  testrel(2, NE2, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, 10);
  f := 6.099056345136382366;
  testrel(3, NE2, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, 100);
  f := 34.40658283378104481;
  testrel(4, NE, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, -1);
  f := 0.3043018436567979761;
  testrel(5, NE, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, -1);
  f := 0.3043018436567979761;
  testrel(6, NE, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, -10);
  f := 0.4539870423425797216e-4;
  testrel(7, NE, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, -100);
  f := 0.3720075976020835963e-43;
  testrel(8, NE, y, f, cnt,failed);

  y := fermi_dirac_r(1.75, 0);
  f := 0.8855362673934438180;
  testrel(9, NE, y, f, cnt,failed);

  y := fermi_dirac_r(1.75, 1);
  f := 2.087331304240103293;
  testrel(10, NE+1, y, f, cnt,failed);    {+1 CPUARM}

  y := fermi_dirac_r(1.75, 10);
  f := 137.18287504419966494;
  testrel(11, NE2, y, f, cnt,failed);

  y := fermi_dirac_r(1.75, 100);
  f := 71553.00102426850928;
  testrel(12, NE, y, f, cnt,failed);

  y := fermi_dirac_r(1.75, -1);
  f := 0.3498499005951593310;
  testrel(13, NE, y, f, cnt,failed);

  y := fermi_dirac_r(1.75, -10);
  f := 0.4539962337472687423e-4;
  testrel(14, NE, y, f, cnt,failed);

  y := fermi_dirac_r(1.75, -100);
  f := 0.3720075976020835963e-43;
  testrel(15, NE, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, 0);
  f := 0.9999999999999999902;
  testrel(16, NE, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, 1);
  f := 2.718281828459045163;
  testrel(17, NE, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, 10);
  f := 22026.46579480195556;
  testrel(18, NE2, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, 100);
  f := 0.3274009421164246242e38;
  testrel(19, NE, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, -1);
  f := 0.3678794411714423203;
  testrel(20, NE, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, -10);
  f := 0.4539992976248485154e-4;
  testrel(21, NE, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, -100);
  f := 0.3720075976020835963e-43;
  testrel(22, NE, y, f, cnt,failed);

  y := fermi_dirac_r(-0.25, -10000);
  f := 0.1135483865314736099e-4342;
  testrel(23, NE, y, f, cnt,failed);

  y := fermi_dirac_r(55.5, -10000);
  f := 0.1135483865314736099e-4342;
  testrel(24, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_bose_einstein;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 2;
  NE1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bose_einstein');

  {Maple: be := (n,x)->polylog(n+1,exp(x));}

  y := bose_einstein(0.5, -10);
  f := 0.4540065850834588316e-4;
  testrel(1, NE, y, f, cnt,failed);

  y := bose_einstein(0.5, -1);
  f := 0.4284407345998380098;
  testrel(2, NE, y, f, cnt,failed);

  y := bose_einstein(0.5, 0);
  f := 2.612375348685488343;
  testrel(3, NE, y, f, cnt,failed);

  y := bose_einstein(0.5, 1);
  f := 1.044217325937597869;
  testrel(4, NE, y, f, cnt,failed);

  y := bose_einstein(0.5, 20);
  f := -66.86797194901337887;
  testrel(5, NE, y, f, cnt,failed);

  y := bose_einstein(0.5, 100);
  f := -752.0671579633630026;
  testrel(6, NE, y, f, cnt,failed);

  y := bose_einstein(0.5, 1000);
  f := -23788.26285334388673;
  testrel(7, NE, y, f, cnt,failed);

  y := bose_einstein(0.5, 1e10);
  f := -752252778063675.0493;
  testrel(8, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, -10);
  f := 0.4540011194644879406e-4;
  testrel(9, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, -1);
  f := 0.3810793119677888313;
  testrel(10, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, 0);
  f := 1.126733867317056646;
  testrel(11, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, 1);
  f := 3.522154705061993875;
  testrel(12, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, 20);
  f := -2854.191162420754632;
  testrel(13, NE1, y, f, cnt,failed);

  y := bose_einstein(2.5, 100);
  f := -857242.5260654477744;
  testrel(14, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, 1000);
  f := -2.718587059372204487e9; {alpha}
  testrel(15, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, 10000);
  f := -8.597172131629545314e12; {alpha}
  testrel(16, NE, y, f, cnt,failed);

  y := bose_einstein(2.5, 1e10);
  f := -8.597174606442000561e33; {alpha}
  testrel(17, NE, y, f, cnt,failed);

  y := bose_einstein(5, -10);
  f := 0.4539996196813856563e-4;
  testrel(18, NE, y, f, cnt,failed);

  y := bose_einstein(5, -1);
  f := 0.3700673152590212488;
  testrel(19, NE, y, f, cnt,failed);

  y := bose_einstein(5, 0);
  f := 1.017343061984449140;
  testrel(20, NE, y, f, cnt,failed);

  y := bose_einstein(5, 1);
  f := 2.882630992436145203;
  testrel(21, NE, y, f, cnt,failed);

  y := bose_einstein(5, 20);
  f := -66521.47068463950671;
  testrel(22, NE, y, f, cnt,failed);

  y := bose_einstein(5, 50);
  f := -20841944.55296836914;
  testrel(23, NE, y, f, cnt,failed);

  y := bose_einstein(5, 100);
  f := -1375170279.731463767;
  testrel(24, NE, y, f, cnt,failed);

  y := bose_einstein(5, 1e9);
  f := -1.388888888888888752e51;  {alpha}
  testrel(25, NE, y, f, cnt,failed);

  y := bose_einstein(0, -1);
  f := 0.4586751453870818910;
  testrel(26, NE, y, f, cnt,failed);

  y := bose_einstein(0, 1.5);
  f := -1.247517541074546004;
  testrel(27, NE, y, f, cnt,failed);

  y := bose_einstein(0, 1000);
  f := -1000.0;
  testrel(28, NE, y, f, cnt,failed);

  y := bose_einstein(-0.5, 1.5);
  f := -1.79535421414243226645;
  testrel(29, NE, y, f, cnt,failed);

  y := bose_einstein(-0.5, -100);
  f := 0.3720075976020835963e-43;
  testrel(30, NE, y, f, cnt,failed);

  y := bose_einstein(-0.5, 5.5);
  f := -2.7258589382195530838;
  testrel(31, NE, y, f, cnt,failed);

  y := bose_einstein(-0.5, 1000);
  f := -35.682511670793217276; {alpha}
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


end.
