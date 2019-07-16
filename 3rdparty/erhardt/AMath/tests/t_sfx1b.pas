{Part 1b of regression test for SPECFUNX unit  (c) 2010-2018  W.Ehrhardt}

unit t_sfx1b;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_igammax;
procedure test_igammanx;
procedure test_igammalx;
procedure test_igammatx;
procedure test_incgammax;
procedure test_igamma_invx;
procedure test_igammap_derx;


implementation


uses
  amath, specfunx, t_sfx0;

{---------------------------------------------------------------------------}
procedure test_incgammax;
var
  a,x,p,q,cp,cq: extended;
  cnt, failed: integer;
const
  NE = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','incgammax');

  a := 0.125;
  x := 0.0009765625;
  q := 0.55359080783499591843;
  p := 0.44640919216500408157;
  incgammax(a,x,cp,cq);
  testrel( 1, NE, cp, p, cnt,failed);
  testrel( 2, NE, cq, q, cnt,failed);

  a := 0.125;
  x := 0.1240234375;
  q := 0.19290797934251178633;
  p := 0.80709202065748821367;
  incgammax(a,x,cp,cq);
  testrel( 3, NE, cp, p, cnt,failed);
  testrel( 4, NE, cq, q, cnt,failed);

  a := 0.125;
  x := 0.1259765625;
  q := 0.19149663495096560193;
  p := 0.80850336504903439807;
  incgammax(a,x,cp,cq);
  testrel( 5, NE, cp, p, cnt,failed);
  testrel( 6, NE, cq, q, cnt,failed);

  a := 0.125;
  x := 1.125;
  q := 0.025367388598510957966;
  p := 0.97463261140148904203;
  incgammax(a,x,cp,cq);
  testrel( 7, NE, cp, p, cnt,failed);
  testrel( 8, NE, cq, q, cnt,failed);

  a := 0.125;
  x := 0.25;
  q := 0.13033534544971364266;
  p := 0.86966465455028635734;
  incgammax(a,x,cp,cq);
  testrel( 9, NE, cp, p, cnt,failed);
  testrel(10, NE, cq, q, cnt,failed);

  a := 0.125;
  x := 100.0;
  q := 8.7052780108565328454E-47;
  p := 1.0;
  incgammax(a,x,cp,cq);
  testrel(11, NE, cp, p, cnt,failed);
  testrel(12, NE, cq, q, cnt,failed);

  a := 0.125;
  x := 0.001;
  q := 0.55226660028386113266;
  p := 0.44773339971613886734;
  incgammax(a,x,cp,cq);
  testrel(13, NE, cp, p, cnt,failed);
  testrel(14, NE, cq, q, cnt,failed);

  a := 0.5;
  x := 0.00390625;
  q := 0.92956802227761292195;
  p := 0.070431977722387078051;
  incgammax(a,x,cp,cq);
  testrel(15, NE, cp, p, cnt,failed);
  testrel(16, NE, cq, q, cnt,failed);

  a := 0.5;
  x := 0.49609375;
  q := 0.31920831748601688161;
  p := 0.68079168251398311839;
  incgammax(a,x,cp,cq);
  testrel(17, NE, cp, p, cnt,failed);
  testrel(18, NE, cq, q, cnt,failed);

  a := 0.5;
  x := 0.50390625;
  q := 0.31542746722370978261;
  p := 0.68457253277629021739;
  incgammax(a,x,cp,cq);
  testrel(19, NE, cp, p, cnt,failed);
  testrel(20, NE, cq, q, cnt,failed);

  a := 0.5;
  x := 1.5;
  q := 0.083264516663550401855;
  p := 0.91673548333644959815;
  incgammax(a,x,cp,cq);
  testrel(21, NE, cp, p, cnt,failed);
  testrel(22, NE, cq, q, cnt,failed);

  a := 0.5;
  x := 1.0;
  q := 0.15729920705028513066;
  p := 0.84270079294971486934;
  incgammax(a,x,cp,cq);
  testrel(23, NE, cp, p, cnt,failed);
  testrel(24, NE, cq, q, cnt,failed);

  a := 0.5;
  x := 100.0;
  q := 2.0884875837625447570E-45;
  p := 1.0;
  incgammax(a,x,cp,cq);
  testrel(25, NE, cp, p, cnt,failed);
  testrel(26, NE, cq, q, cnt,failed);

  a := 0.5;
  x := 0.001;
  q := 0.96432940827032011495;
  p := 0.035670591729679885046;
  incgammax(a,x,cp,cq);
  testrel(27, NE, cp, p, cnt,failed);
  testrel(28, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 0.009765625;
  q := 0.99730512061186844099;
  p := 0.0026948793881315590139;
  incgammax(a,x,cp,cq);
  testrel(29, NE, cp, p, cnt,failed);
  testrel(30, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 1.240234375;
  q := 0.38474791393713691285;
  p := 0.61525208606286308715;
  incgammax(a,x,cp,cq);
  testrel(31, NE, cp, p, cnt,failed);
  testrel(32, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 1.259765625;
  q := 0.37822004134435599166;
  p := 0.62177995865564400834;
  incgammax(a,x,cp,cq);
  testrel(33, NE, cp, p, cnt,failed);
  testrel(34, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 0.25;
  q := 0.86388362845858884258;
  p := 0.13611637154141115742;
  incgammax(a,x,cp,cq);
  testrel(35, NE, cp, p, cnt,failed);
  testrel(36, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 2.25;
  q := 0.15499954779476378341;
  p := 0.84500045220523621659;
  incgammax(a,x,cp,cq);
  testrel(37, NE, cp, p, cnt,failed);
  testrel(38, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 2.5;
  q := 0.12308857115265884524;
  p := 0.87691142884734115476;
  incgammax(a,x,cp,cq);
  testrel(39, NE, cp, p, cnt,failed);
  testrel(40, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 100.0;
  q := 1.301089352558546598E-43;
  p := 1.0;
  incgammax(a,x,cp,cq);
  testrel(41, NE, cp, p, cnt,failed);
  testrel(42, NE, cq, q, cnt,failed);

  a := 1.25;
  x := 0.001;
  q := 0.999843134425254734;
  p := 0.00015686557474526599535;
  incgammax(a,x,cp,cq);
  testrel(43, NE, cp, p, cnt,failed);
  testrel(44, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 0.0966796875;
  q := 1.0;
  p := 2.0443064129172721431E-22;
  incgammax(a,x,cp,cq);
  testrel(45, NE, cp, p, cnt,failed);
  testrel(46, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 12.2783203125;
  q := 0.47311590470245409997;
  p := 0.52688409529754590003;
  incgammax(a,x,cp,cq);
  testrel(47, NE, cp, p, cnt,failed);
  testrel(48, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 12.4716796875;
  q := 0.45133715451332652040;
  p := 0.54866284548667347960;
  incgammax(a,x,cp,cq);
  testrel(49, NE, cp, p, cnt,failed);
  testrel(50, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 11.375;
  q := 0.57797890730569987083;
  p := 0.42202109269430012917;
  incgammax(a,x,cp,cq);
  testrel(51, NE, cp, p, cnt,failed);
  testrel(52, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 13.375;
  q := 0.35521881768820689405;
  p := 0.64478118231179310595;
  incgammax(a,x,cp,cq);
  testrel(53, NE, cp, p, cnt,failed);
  testrel(54, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 24.75;
  q := 0.0022231460245010735082;
  p := 0.99777685397549892649;
  incgammax(a,x,cp,cq);
  testrel(55, NE, cp, p, cnt,failed);
  testrel(56, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 100.0;
  q := 2.3484886146232989951E-29;
  p := 1.0;
  incgammax(a,x,cp,cq);
  testrel(57, NE, cp, p, cnt,failed);
  testrel(58, NE, cq, q, cnt,failed);

  a := 12.375;
  x := 0.001;
  q := 1.0;
  p := 6.0319676270277627123E-47;
  incgammax(a,x,cp,cq);
  testrel(59, NE, cp, p, cnt,failed);
  testrel(60, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 0.193359375;
  q := 1.0;
  p := 2.6139638700232537181E-43;
  incgammax(a,x,cp,cq);
  testrel(61, NE, cp, p, cnt,failed);
  testrel(62, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 24.556640625;
  q := 0.48877461095124572049;
  p := 0.51122538904875427951;
  incgammax(a,x,cp,cq);
  testrel(63, NE, cp, p, cnt,failed);
  testrel(64, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 24.943359375;
  q := 0.45787485849439974988;
  p := 0.54212514150560025012;
  incgammax(a,x,cp,cq);
  testrel(65, NE, cp, p, cnt,failed);
  testrel(66, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 23.75;
  q := 0.55428166101076371554;
  p := 0.44571833898923628446;
  incgammax(a,x,cp,cq);
  testrel(67, NE, cp, p, cnt,failed);
  testrel(68, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 25.75;
  q := 0.39542497277452763030;
  p := 0.60457502722547236970;
  incgammax(a,x,cp,cq);
  testrel(69, NE, cp, p, cnt,failed);
  testrel(70, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 49.5;
  q := 0.000037467135608746412312;
  p := 0.99996253286439125359;
  incgammax(a,x,cp,cq);
  testrel(71, NE, cp, p, cnt,failed);
  testrel(72, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 100.0;
  q := 5.5034695332350412741E-20;
  p := 0.99999999999999999994;
  incgammax(a,x,cp,cq);
  testrel(73, NE, cp, p, cnt,failed);
  testrel(74, NE, cq, q, cnt,failed);

  a := 24.75;
  x := 0.001;
  q := 1.0;
  p := 8.1291472534766856537E-100;
  incgammax(a,x,cp,cq);
  testrel(75, NE, cp, p, cnt,failed);
  testrel(76, NE, cq, q, cnt,failed);

  a := 245.675;
  x := 1.9193359375;
  q := 1.0;
  p := 3.8186081058738398329E-414;
  incgammax(a,x,cp,cq);
  testrele(77, 1e-17, cp, p, cnt,failed);    {!!!!!!}
  testrel(78, NE, cq, q, cnt,failed);

  a := 245.675;
  x := 243.7556640625;
  q := 0.54041968599852003456;
  p := 0.45958031400147996544;
  incgammax(a,x,cp,cq);
  testrel(79, NE, cp, p, cnt,failed);
  testrel(80, NE, cq, q, cnt,failed);

  a := 245.675;
  x := 247.5943359375;
  q := 0.44299079551595765110;
  p := 0.55700920448404234890;
  incgammax(a,x,cp,cq);
  testrel(81, NE, cp, p, cnt,failed);
  testrel(82, NE, cq, q, cnt,failed);

  a := 245.675;
  x := 244.675;
  q := 0.51699406317227828678;
  p := 0.48300593682772171322;
  incgammax(a,x,cp,cq);
  testrel(83, NE, cp, p, cnt,failed);
  testrel(84, NE, cq, q, cnt,failed);

  a := 245.675;
  x := 246.675;
  q := 0.46614064932598652437;
  p := 0.53385935067401347563;
  incgammax(a,x,cp,cq);
  testrel(85, NE, cp, p, cnt,failed);
  testrel(86, NE, cq, q, cnt,failed);

  a := 245.675;
  x := 491.35;
  q := 4.5956807163396664910E-35;
  p := 1.0;
  incgammax(a,x,cp,cq);
  testrel(87, NE, cp, p, cnt,failed);
  testrele(88, 1e-17, cq, q, cnt,failed);        {!!!!!!!!}

  a := 245.675;
  x := 100.0;
  q := 1.0;
  p := 9.8734129137566074939E-35;
  incgammax(a,x,cp,cq);
  testrele(89, 4e-17, cp, p, cnt,failed);        {!!!!!!!! AMD 4e-17,  Intel 3e-17}
  testrel(90, NE, cq, q, cnt,failed);

  a := 245.675;
  x := 0.001;
  q := 1.0;
  p := 6.6606493267933284947E-1220;
  incgammax(a,x,cp,cq);
  testrele(91, 4e-17, cp, p, cnt,failed);        {!!!!!!!!}
  testrel(92, NE, cq, q, cnt,failed);

  a := 2456.75;
  x := 19.193359375;
  q := 1.0;
  p := 4.3493575976177286679E-4121;
  incgammax(a,x,cp,cq);
  testrele(93, 7e-17, cp, p, cnt,failed);        {!!!!!!!!}
  testrel(94, NE, cq, q, cnt,failed);

  a := 2456.75;
  x := 2437.556640625;
  q := 0.64858936942781918995;
  p := 0.35141063057218081005;
  incgammax(a,x,cp,cq);
  testrel(95, NE, cp, p, cnt,failed);
  testrel(96, NE, cq, q, cnt,failed);

  a := 2456.75;
  x := 2475.943359375;
  q := 0.34717883377068764690;
  p := 0.65282116622931235310;
  incgammax(a,x,cp,cq);
  testrel(97, NE, cp, p, cnt,failed);
  testrel(98, NE, cq, q, cnt,failed);

  a := 2456.75;
  x := 2455.75;
  q := 0.50536666063737615989;
  p := 0.49463333936262384011;
  incgammax(a,x,cp,cq);
  testrel( 99, NE, cp, p, cnt,failed);
  testrel(100, NE, cq, q, cnt,failed);

  a := 2456.75;
  x := 2457.75;
  q := 0.48927075527512612274;
  p := 0.51072924472487387726;
  incgammax(a,x,cp,cq);
  testrel(101, NE, cp, p, cnt,failed);
  testrel(102, NE, cq, q, cnt,failed);

  a := 2456.75;
  x := 4913.5;
  q := 3.2198475192983317356E-330;
  p := 1.0;
  incgammax(a,x,cp,cq);
  testrel(102, NE, cp, p, cnt,failed);
  testrele(104, 6e-17, cq, q, cnt,failed);       {!!!!!!!!}

  a := 2456.75;
  x := 100.0;
  q := 1.0;
  p := 4.7640315099546672134E-2395;
  incgammax(a,x,cp,cq);
  testrele(105, 4e-18, cp, p, cnt,failed);       {!!!!!!!!}
  testrel(106, NE, cq, q, cnt,failed);

  {test cases erf/erfc and igammapx/igammaqx}
  a := 0.5;

  x := ldexp(1,-20);
  q := 0.9988980675699281853;
  p := 0.0011019324300718147042;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(107, NE, cp, p, cnt,failed);
  testrel(108, NE, cq, q, cnt,failed);

  x := 0.0078125;
  q := 0.90052355033977421404;
  p := 0.099476449660225785959;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(109, NE, cp, p, cnt,failed);
  testrel(110, NE, cq, q, cnt,failed);

  x := 0.125;
  q := 0.61707507745197379272;
  p := 0.38292492254802620728;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(111, NE, cp, p, cnt,failed);
  testrel(112, NE, cq, q, cnt,failed);

  x := 0.5;
  q := 0.31731050786291410283;
  p := 0.68268949213708589717;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(113, NE, cp, p, cnt,failed);
  testrel(114, NE, cq, q, cnt,failed);

  x := 0.9;
  q := 0.17971249487899984212;
  p := 0.82028750512100015788;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(115, NE, cp, p, cnt,failed);
  testrel(116, NE, cq, q, cnt,failed);

  x := 1.0;
  q := 0.15729920705028513066;
  p := 0.84270079294971486934;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(117, NE, cp, p, cnt,failed);
  testrel(118, NE, cq, q, cnt,failed);

  x := 1.1;
  q := 0.13801073756865955964;
  p := 0.86198926243134044036;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(119, NE, cp, p, cnt,failed);
  testrel(120, NE, cq, q, cnt,failed);

  x := 2.0;
  q := 0.045500263896358414401;
  p := 0.9544997361036415856;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(121, NE, cp, p, cnt,failed);
  testrel(122, NE, cq, q, cnt,failed);

  x := 4.0;
  q := 0.0046777349810472658379;
  p := 0.99532226501895273416;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(123, NE, cp, p, cnt,failed);
  testrel(124, NE, cq, q, cnt,failed);

  x := 6.0;
  q := 0.00053200550513924969929;
  p := 0.9994679944948607503;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(125, NE, cp, p, cnt,failed);
  testrel(126, NE, cq, q, cnt,failed);

  x := 10.0;
  q := 0.0000077442164310440836377;
  p := 0.99999225578356895592;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(127, NE, cp, p, cnt,failed);
  testrel(128, NE, cq, q, cnt,failed);

  x := 40.0;
  q := 3.7440973842028987635e-19;
  p := 0.99999999999999999963;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(129, NE, cp, p, cnt,failed);
  testrel(130, NE, cq, q, cnt,failed);

  x := 100.0;
  q := 2.0884875837625447570E-45;
  p := 1.0;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(131, NE, cp, p, cnt,failed);
  testrel(132, NE, cq, q, cnt,failed);

  x := 10000.0;
  q := 6.4059614249217320390E-4346;
  p := 1.0;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(133, NE, cp, p, cnt,failed);
  testrel(134, NE, cq, q, cnt,failed);

  x := 20000.0;
  q := 0.0;
  p := 1.0;
  cp := igammapx(a,x);
  cq := igammaqx(a,x);
  testrel(135, NE, cp, p, cnt,failed);
  testrel(136, NE, cq, q, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_igammax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 10;
  NE1 = 25;
  NE2 = 75;
  NE3 = 750;   {extreme case}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','igammax');

  {Test values from GSL}
  y := igammax(-1.0/1048576.0, 1.0/1048576.0);
  f := 13.28581959629062427;
  testrel( 1, NE, y, f, cnt,failed);

  y := igammax(-0.001, 1.0/1048576.0);
  f := 13.38127512862532886;
  testrel( 2, NE, y, f, cnt,failed);

  y := igammax(-1.0, 1.0/1048576.0);
  f := 1.048561714271576866e+06;
  testrel( 3, NE, y, f, cnt,failed);

  y := igammax(-0.00001,0.001);
  f := 6.331768143456359214;
  testrel( 4, NE, y, f, cnt,failed);

  y := igammax(-0.0001,0.001);
  f := 6.333827643976718939;
  testrel( 5, NE, y, f, cnt,failed);

  y := igammax(-0.001, 0.001);
  f := 6.354470910251084379;
  testrel( 6, NE, y, f, cnt,failed);

  y := igammax(-0.5, 0.001);
  f := 59.763880515942196981;
  testrel( 7, NE, y, f, cnt,failed);

  y := igammax(-1.0, 0.001);
  f := 992.6689604692388423;
  testrel( 8, NE, y, f, cnt,failed);

  y := igammax(-3.5, 0.001);
  f := 9.0224404490639003711e+09;
  testrel( 9, NE1, y, f, cnt,failed);

  y := igammax(-10.5, 0.001);
  f := 3.008366155818481566e+30;
  testrel(10, NE1, y, f, cnt,failed);

  y := igammax(-0.001, 0.1);
  f := 1.824910960941862007;
  testrel(11, NE, y, f, cnt,failed);

  y := igammax(-0.5, 0.1);
  f := 3.401769336691615416;
  testrel(12, NE, y, f, cnt,failed);

  y := igammax(-10.0, 0.1);
  f := 8.949075748358698918e+08;
  testrel(13, NE, y, f, cnt,failed);

  y := igammax(-10.5, 0.1);
  f := 2.696740383422642177e+09;
  testrel(14, NE, y, f, cnt,failed);

  y := igammax(-0.001, 1.0);
  f := 0.2192861267907276634;
  testrel(15, NE, y, f, cnt,failed);

  y := igammax(-0.5, 1.0);
  f := 0.1781477117815606902;
  testrel(16, NE, y, f, cnt,failed);

  y := igammax(-1.0, 1.0);
  f := 0.1484955067759220479;
  testrel(17, NE, y, f, cnt,failed);

  y := igammax(-2.5, 1.0);
  f := 0.9655664863127516026e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := igammax(-1.0, 10.0);
  f := 3.830240465631608762e-07;
  testrel(19, NE, y, f, cnt,failed);

  y := igammax(-0.001, 10.0);
  f := 4.147056232480732096e-06;
  testrel(20, NE, y, f, cnt,failed);

  y := igammax(-0.5, 10.0);
  f := 1.260904261324157068e-06;
  testrel(21, NE, y, f, cnt,failed);

  y := igammax(-1.0, 10.0);
  f := 3.830240465631608762e-07;
  testrel(22, NE, y, f, cnt,failed);

  y := igammax(-10.5, 10.0);
  f := 6.840492732844156679e-17;
  testrel(23, NE, y, f, cnt,failed);

  y := igammax(-100.0, 10.0);
  f := 4.123832766985831400e-107;
  testrel(24, NE, y, f, cnt,failed);

  y := igammax(-200.0, 10.0);
  f := 2.161409183052934342e-207;
  testrel(25, NE, y, f, cnt,failed);

  y := igammax(0.0, 0.001);
  f := 6.331539364136149332;
  testrel(26, NE, y, f, cnt,failed);

  y := igammax(0.001, 0.001);
  f := 6.308715939486400726;
  testrel(27, NE, y, f, cnt,failed);

  y := igammax(1.0, 0.001);
  f := 0.99900049983337499167;
  testrel(28, NE, y, f, cnt,failed);

  y := igammax(10.0, 0.001);
  f := 362880.0;
  testrel(29, NE, y, f, cnt,failed);

  y := igammax(0.0, 1.0);
  f := 0.21938393439552027368;
  testrel(30, NE, y, f, cnt,failed);

  y := igammax(0.001, 1.0);
  f := 0.21948181320730279613;
  testrel(31, NE, y, f, cnt,failed);

  y := igammax(1.0, 1.0);
  f := 0.36787944117144232160;
  testrel(32, NE, y, f, cnt,failed);

  y := igammax(10.0, 1.0);
  f := 362879.9595659224205;
  testrel(33, NE, y, f, cnt,failed);

  y := igammax(100.0, 1.0);
  f := 9.332621544394415268e+155;
  testrel(34, NE, y, f, cnt,failed);

  y := igammax(0.0, 100.0);
  f := 3.683597761682032180e-46;
  testrel(35, NE, y, f, cnt,failed);

  y := igammax(0.001, 100.0);
  f := 3.700636767406355063e-46;
  testrel(36, NE, y, f, cnt,failed);

  y := igammax(1.0, 100.0);
  f := 3.720075976020835963e-44;
  testrel(37, NE, y, f, cnt,failed);

  y := igammax(10.0, 100.0);
  f := 4.083660630910611272e-26;
  testrel(38, NE, y, f, cnt,failed);

  y := igammax(100.0, 100.0);
  f := 4.542198120862669429e+155;
  testrel(39, NE, y, f, cnt,failed);

  {calculated with Pari};
  y := igammax(-10000,1);
  f := 3.678426532276930936E-5;
  testrel(40, NE, y, f, cnt,failed);

  y := igammax(-5000,0.998046875);
  f := 1.296585775691598779;
  testrel(41, NE, y, f, cnt,failed);

  y := igammax(-5.5,3);
  f := 1.334725498656031944e-5;
  testrel(42, NE, y, f, cnt,failed);

  y := igammax(-6,3);
  f := 7.309373616318302801e-6;
  testrel(43, NE, y, f, cnt,failed);

  y := igammax(0,4);
  f := 0.3779352409848906479e-2;
  testrel(44, NE, y, f, cnt,failed);

  y := igammax(100,50);
  f := 0.9332621541407915409e156;
  testrel(45, NE, y, f, cnt,failed);

  y := igammax(2,0);
  f := 1;
  testrel(46, NE, y, f, cnt,failed);

  y := igammax(2,-1);
  f := 0;
  testrel(47, NE, y, f, cnt,failed);

  y := igammax(-1,2.5);
  f := 0.7919081579289782572e-2;
  testrel(48, NE, y, f, cnt,failed);

  y := igammax(0,2.5);
  f := 0.24914917870269735496e-1;
  testrel(49, NE, y, f, cnt,failed);

  y := igammax(0,-2.5);
  f := -7.073765894578600712;
  testrel(50, NE, y, f, cnt,failed);

  y := igammax(0,-100);
  f := -0.2715552744853879822e42;
  testrel(51, NE, y, f, cnt,failed);

  y := igammax(5,2.5);
  f := 21.38827245393962982;
  testrel(52, NE, y, f, cnt,failed);

  y := igammax(5,-2.5);
  f := 189.5900622634478054;
  testrel(53, NE, y, f, cnt,failed);

  y := igammax(5,-10);
  f := 153832837.1109301082;
  testrel(54, NE, y, f, cnt,failed);

  y := igammax(5,-100);
  f := 0.2583754327050379842e52;
  testrel(55, NE, y, f, cnt,failed);

  y := igammax(100,-10);
{$ifdef BIT16}
  f := 0.18665243088788830536e157 * 0.5;
{$else}
  f := 0.9332621544394415268e156;
{$endif}
  testrel(56, NE, y, f, cnt,failed);

  {negative x}
  y := igammax(100,-100);
  f := -0.1347427096011818133e242;
  testrel(57, NE, y, f, cnt,failed);

  y := igammax(2.5, -10);
  f := 1.329340388179137020;
  testrel(58, NE, y, f, cnt,failed);

  y := igammax(-2.5, -10);
  f := -0.9453087204829418812;
  testrel(59, NE, y, f, cnt,failed);

  y := igammax(3.5, -4.5);
  f := 3.323350970447842551;
  testrel(60, NE, y, f, cnt,failed);

  y := igammax(-10.25, -20);
  f := -2.494160973094416239e-6;
  testrel(61, NE2, y, f, cnt,failed);

  y := igammax(-10.25, 20);
  f := 3.079892386669366206e-24;
  testrel(62, NE, y, f, cnt,failed);

  y := igammax(10.25, -20);
  f := -2.511354417362715712e20;
  testrel(63, NE, y, f, cnt,failed);

  y := igammax(-1/3, -10);
  f := -609.0884007167802168;
  testrel(64, NE, y, f, cnt,failed);

  y := igammax(1/3, 10);
  f := 9.216133861986501190e-6;
  testrel(65, NE, y, f, cnt,failed);

  y := igammax(-1/3, 10);
  f := 1.8763219551451487814e-6;
  testrel(66, NE, y, f, cnt,failed);

  y := igammax(10.25, -200);
  f := -9.401145720184358356e107;
  testrel(67, NE, y, f, cnt,failed);

  y := igammax(-0.75, -200);
  f := 4.846355047687434826e82;
  testrel(68, NE2, y, f, cnt,failed);


  {extreme cases}
  y := igammax(1200,11000);
  f := 2.758325061796967207E68;
  testrel(69, NE3, y, f, cnt,failed);

  y := igammax(5000,55000);
  f := 8.278720754962823876E-190;
  testrel(70, NE3, y, f, cnt,failed);

  y := igammax(1700,16500);
  f := 0.4936956875959810692;
  testrel(71, NE3, y, f, cnt,failed);


  {extended only}

  y := igammax(1700,0.1);
  f := 1.763737179740245938E4752;
  testrel(72, NE, y, f, cnt,failed);

  y := igammax(1800,2e4);
  f := 0.5061235229092104215e-948;
  testrel(73, NE3, y, f, cnt,failed);

  y := igammax(1755.5,1756);
  f := 4.079759411085472574E4931;
  testrel(74, NE3, y, f, cnt,failed);

  y := igammax(2,10000);
  f := 1.135597413701267572E-4339;
  testrel(75, NE3, y, f, cnt,failed);

  y := igammax(2,11356);
  f := 1.611117193805680611E-4928;
  testrel(76, NE3, y, f, cnt,failed);

  y := igammax(5,-1000);
  f := 0.1962214423179921981e447;
  testrel(77, NE, y, f, cnt,failed);

  {negative x}
  y := igammax(5,-10000);
  f := 8.803296554979500997e4358;
  testrel(78, NE, y, f, cnt,failed);

  y := igammax(200,-10);
  f := 0.3943289336823952518e373;
  testrel(79, NE, y, f, cnt,failed);

  y := igammax(500,-10);
  f := 2.440273651982220137e1131;
  testrel(80, NE, y, f, cnt,failed);

  y := igammax(1700, -475);
  f := -1.969477481521004251e4753;
  testrel(81, NE2, y, f, cnt,failed);

  y := igammax(1234.5, -200);
  f := 1.454091272346110182e3279;
  testrel(82, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_igammalx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE1 = 20;
  NE2 = 50;
  NE3 = 600;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','igammalx');

  y := igammalx(-1.0/1048576.0, 1.0/1048576.0);
  f := -1048589.86303620443;
  testrel( 1, NE, y, f, cnt,failed);

  y := igammalx(-0.001, 1.0/1048576.0);
  f := -1013.959480757983977;
  testrel( 2, NE, y, f, cnt,failed);

  y := igammalx(-0.00001,0.001);
  f := -100006.9089936990086;
  testrel( 3, NE, y, f, cnt,failed);

  y := igammalx(-0.001, 0.001);
  f := -1006.932676539609732;
  testrel( 4, NE, y, f, cnt,failed);

  y := igammalx(-0.5, 0.001);
  f := -63.30878821775322904;
  testrel( 5, NE, y, f, cnt,failed);

  y := igammalx(-0.001, 0.1);
  f := -1002.40311659030051;
  testrel( 6, NE, y, f, cnt,failed);

  y := igammalx(-0.5, 0.1);
  f := -6.946677038502647471;
  testrel( 7, NE, y, f, cnt,failed);

  y := igammalx(-0.001, 1.0);
  f := -1000.797491756149376;
  testrel( 8, NE, y, f, cnt,failed);

  y := igammalx(-0.5, 1.0);
  f := -3.723055413592592745;
  testrel( 9, NE, y, f, cnt,failed);

  y := igammalx(-2.5, 1.0);
  f := -1.041865369114217041;
  testrel(10, NE, y, f, cnt,failed);

  y := igammalx(-0.001, 10.0);
  f := -1000.578209776414880;
  testrel(11, NE, y, f, cnt,failed);

  y := igammalx(-0.5, 10.0);
  f := -3.544908962715293379;
  testrel(12, NE, y, f, cnt,failed);

  y := igammalx(-10.5, 10.0);
  f := -2.640121821231765590e-7;
  testrel(13, NE, y, f, cnt,failed);

  y := igammalx(0.001, 0.001);
  f := 993.1150565451090654;
  testrel(14, NE, y, f, cnt,failed);

  y := igammalx(1.0, 0.001);
  f := 9.99500166625008332e-4;
  testrel(15, NE, y, f, cnt,failed);

  y := igammalx(10.0, 0.001);
  f := 9.990913256294003857e-32;
  testrel(16, NE, y, f, cnt,failed);

  y := igammalx(0.001, 1.0);
  f := 999.2042906713881633;
  testrel(17, NE, y, f, cnt,failed);

  y := igammalx(1.0, 1.0);
  f := 0.6321205588285576784;
  testrel(18, NE, y, f, cnt,failed);

  y := igammalx(10.0, 1.0);
  f := 4.043407757955495940e-2;
  testrel(19, NE, y, f, cnt,failed);

  y := igammalx(100.0, 1.0);
  f := 3.715578714528098103e-3;
  testrel(20, NE1, y, f, cnt,failed);

  y := igammalx(0.001, 100.0);
  f := 999.4237724845954661;
  testrel(21, NE, y, f, cnt,failed);

  y := igammalx(1.0, 100.0);
  f := 1.0;
  testrel(22, NE, y, f, cnt,failed);

  y := igammalx(10.0, 100.0);
  f := 362880.0;
  testrel(23, NE, y, f, cnt,failed);

  y := igammalx(100.0, 100.0);
  f := 4.790423423531745839e155;
  testrel(24, NE, y, f, cnt,failed);

  y := igammalx(-5.5,3);
  f := 1.089930752692330267e-2;
  testrel(25, NE, y, f, cnt,failed);

  y := igammalx(100,50);
  f := 2.986499859169264780e146;
  testrel(26, NE1, y, f, cnt,failed);

  y := igammalx(-1.5,1);
  f := 2.236783981614100282;
  testrel(27, NE, y, f, cnt,failed);

  y := igammalx(-99.5,2);
  f := -1.487220346103121896e-33;
  testrel(28, NE1, y, f, cnt,failed);

  y := igammalx(-99.5,0.125);
  f := -6.379593929924726477e87;
  testrel(29, NE1, y, f, cnt,failed);

  y := igammalx(-1234.25,1.5);
  f := -8.241205460141793872e-222;
  testrel(30, NE1, y, f, cnt,failed);

  y := igammalx(-1 -1/8192, 0.5);
  f := 8190.924081011363356;
  testrel(31, NE, y, f, cnt,failed);

  y := igammalx(2,1);
  f := 0.2642411176571153568;
  testrel(32, NE, y, f, cnt,failed);

  y := igammalx(2,4);
  f := 0.9084218055563290985;
  testrel(33, NE, y, f, cnt,failed);

  y := igammalx(1e-10,3.5);
  f := 9999999999.415814195;
  testrel(34, NE, y, f, cnt,failed);

  y := igammalx(1e-10,4);
  f := 9999999999.419004983;
  testrel(35, NE, y, f, cnt,failed);

  y := igammalx(1/1024,4);
  f := 1023.419964213348909;
  testrel(36, NE, y, f, cnt,failed);

  y := igammalx(1e-100,3.5);
  f := 1e100;
  testrel(37, NE, y, f, cnt,failed);

  y := igammalx(10000,1+1/1024);
  f := 6.374004880783332249e-1;
  testrel(38, NE, y, f, cnt,failed);

  y := igammalx(10000,1);
  f := 3.679162291151916232e-5;
  testrel(39, NE, y, f, cnt,failed);

  y := igammalx(100,100);
  f := 4.790423423531745839e155;
  testrel(40, NE, y, f, cnt,failed);

  y := igammalx(100,99.5);
  f := 4.604031018526237350e155;
  testrel(41, NE, y, f, cnt,failed);

  {x < 0}
  y := igammalx(0.75, -1.25);
  f := -2.044071997352096675;
  testrel(42, NE, y, f, cnt,failed);

  y := igammalx(-0.75, -1.25);
  f := -2.676560824639220107;
  testrel(43, NE, y, f, cnt,failed);

  y := igammalx(0.5, -1);
  f := 0;
  testrel(44, NE, y, f, cnt,failed);

  y := igammalx(2, -1);
  f := 1;
  testrel(45, NE, y, f, cnt,failed);

  y := igammalx(3, -1);
  f := -0.7182818284590452354;
  testrel(46, NE, y, f, cnt,failed);

  y := igammalx(-1.5, -1);
  f := 0;
  testrel(47, NE, y, f, cnt,failed);

  y := igammalx(1.25, -10);
  f := -26938.36962856737398;
  testrel(48, NE, y, f, cnt,failed);

  y := igammalx(1.25, -100);
  f := -5.995670897521874460e43;
  testrel(49, NE, y, f, cnt,failed);

  y := igammalx(-1.25, -100);
  f := -6.150648625679935847e38;
  testrel(50, NE2, y, f, cnt,failed);

  y := igammalx(-32.125, -100);
  f := 2.104684115131538288e-23;
  testrel(51, NE2, y, f, cnt,failed);

  y := igammalx(32.125, -100);
  f := 3.361888520003071388e105;
  testrel(52, NE, y, f, cnt,failed);

  y := igammalx(32.125, 100);
  f := 1.265961294960876989e34;
  testrel(53, NE, y, f, cnt,failed);

  y := igammalx(-32.125, -10);
  f := -7.056503991519611836e-30;
  testrel(54, NE2, y, f, cnt,failed);

  {extended}
  y := igammalx(1700,0.125);
  f := 2.899467852314774796e-1539;
  testrel(55, NE3, y, f, cnt,failed);

  y := igammalx(1000,4);
  f := 2.111311496567870236e597;
  testrel(56, NE3, y, f, cnt,failed);

  y := igammalx(1000,1001);
  f := 2.079583445208785949e2564;
  testrel(57, NE, y, f, cnt,failed);

  y := igammalx(1000,999);
  f := 1.978081150765354015e2564;
  testrel(58, NE, y, f, cnt,failed);

  y := igammalx(2000,200);
  f := 8.826644856281366284e4511;
  testrel(59, 2500, y, f, cnt,failed); {FPC271}

  y := igammalx(1e-4930, 1000);
  f := 1e4930;
  testrel(60, NE, y, f, cnt,failed);

  y := igammalx(1755.5, 2000);
  f := 8.291075796613038194e4931;
  testrel(61, NE3, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_igammatx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 8;
{$ifdef FPC}
  NE2 = 44;
{$else}
  NE2 = 40;
{$endif}
  NE2A = 1100;          {Pure Double/AMD, NE2 is OK with I387 !!!!!!!!}
  NE3 = 400;
  NA  = 32;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','igammatx');

  {MAPLE: git := (a,x) -> (1-GAMMA(a,x)/GAMMA(a))/x^a;}
  {PARI : git(a,x) = incgamc(a,x)/gamma(a)/x^a }

  y := igammatx(-10.5,0);
  f := 360733.7150008375899;
  testrel( 1, NE, y, f, cnt,failed);

  y := igammatx(-10,2);
  f := 1024;
  testrel( 2, NE, y, f, cnt,failed);

  y := igammatx(-1.0/1048576.0, 1.0/1048576.0);
  f := 0.9999994495245583333;
  testrel( 3, NE, y, f, cnt,failed);

  y := igammatx(-0.001, 0.001);
  f := 0.9994231286718222221;
  testrel( 4, NE, y, f, cnt,failed);

  y := igammatx(-0.5, 0.001);
  f := 0.5647536791185097470;
  testrel( 5, NE, y, f, cnt,failed);

  y := igammatx(-0.001, 0.1);
  f := 0.9995197254661496346;
  testrel( 6, NE, y, f, cnt,failed);

  y := igammatx(-0.5, 0.1);
  f := 0.6196867015758642738;
  testrel( 7, NE, y, f, cnt,failed);

  y := igammatx(-0.001, 1.0);
  f := 1.000219159407587534;
  testrel( 8, NE, y, f, cnt,failed);

  y := igammatx(-0.5, 1.0);
  f := 1.050254541660012221;
  testrel( 9, NE, y, f, cnt,failed);

  y := igammatx(-2.5, 1.0);
  f := 1.102142978837586559;
  testrel(10, NE, y, f, cnt,failed);

  y := igammatx(-0.001, 10.0);
  f := 1.002305242232113866;
  testrel(11, NE, y, f, cnt,failed);

  y := igammatx(-0.5, 10.0);
  f := 3.162278784973229726;
  testrel(12, NE, y, f, cnt,failed);

  y := igammatx(-10.5, 10.0);
  f := 31622776609.87717939;
  testrel(13, NE, y, f, cnt,failed);

  y := igammatx(0.001, 0.001);
  f := 1.000575560417974683;
  testrel(14, NE, y, f, cnt,failed);

  y := igammatx(1.0, 0.001);
  f := 0.999500166625008332;
  testrel(15, NE, y, f, cnt,failed);

  y := igammatx(10.0, 0.001);
  f := 2.753227859428462262e-7;
  testrel(16, NE, y, f, cnt,failed);

  y := igammatx(0.001, 1.0);
  f := 0.9997803916424144436;
  testrel(17, NE, y, f, cnt,failed);

  y := igammatx(1.0, 1.0);
  f := 0.6321205588285576784;
  testrel(18, NE, y, f, cnt,failed);

  y := igammatx(10.0, 1.0);
  f := 0.1114254783387206774e-6;
  testrel(19, NE, y, f, cnt,failed);

  y := igammatx(100.0, 1.0);
  f := 3.981280818956854411E-159;
  testrel(20, NE, y, f, cnt,failed);

  y := igammatx(0.001, 100.0);
  f := 0.9954054173515269624;
  testrel(21, NE, y, f, cnt,failed);

  y := igammatx(1.0, 100.0);
  f := 0.01;
  testrel(22, NE, y, f, cnt,failed);

  y := igammatx(10.0, 100.0);
  f := 1e-20;
  testrel(23, NE, y, f, cnt,failed);

  y := igammatx(100.0, 100.0);
  f := 5.132987982791486649E-201;
  testrel(24, NE, y, f, cnt,failed);

  y := igammatx(-5.5,3);
  f := 420.3735582073240524;
  testrel(25, NE, y, f, cnt,failed);

  y := igammatx(100,50);
  f := 4.056564729479877823E-180;
  testrel(26, NE, y, f, cnt,failed);

  y := igammatx(-1.5,1);
  f := 0.9464776673048635452;
  testrel(27, NE, y, f, cnt,failed);

  y := igammatx(-99.5,2);
  f := -3.9552214374489233866e153;
  testrel(28, NE2, y, f, cnt,failed);

  y := igammatx(-99.5,0.125);
  f := -2.628150592774871258e154;
  testrel(29, NE, y, f, cnt,failed);

  y := igammatx(-1 -1/8192, 0.5);
  f := 0.4999178220658637586;
  testrel(30, NE, y, f, cnt,failed);

  y := igammatx(2,1);
  f := 0.2642411176571153568;
  testrel(31, NE, y, f, cnt,failed);

  y := igammatx(2,4);
  f := 0.5677636284727056866e-1;
  testrel(32, NE, y, f, cnt,failed);

  y := igammatx(1e-10,3.5);
  f := 0.9999999998740266892;
  testrel(33, NE, y, f, cnt,failed);

  y := igammatx(1e-10,4);
  f := 0.9999999998609926287;
  testrel(34, NE, y, f, cnt,failed);

  y := igammatx(1/1024,4);
  f := 0.998643419395128638;
  testrel(35, NE, y, f, cnt,failed);

  y := igammatx(0,3.5);
  f := 1;
  testrel(36, NE, y, f, cnt,failed);

  y := igammatx(160,1+1/1024);
  f := 7.843928041790720964e-286;
  testrel(37, NE, y, f, cnt,failed);

  y := igammatx(160,1);
  f := 7.851543950671356050e-286;
  testrel(38, NE, y, f, cnt,failed);

  y := igammatx(100,100);
  f := 5.132987982791486649e-201;
  testrel(39, NE, y, f, cnt,failed);

  y := igammatx(100,99.5);
  f := 8.143788976857973453e-201;
  testrel(40, NE, y, f, cnt,failed);

  y := igammatx(-10 - 1/1048576.0, 0.5);
  f := 0.1998981668022373637;
  testrel(41, NE, y, f, cnt,failed);

  y := igammatx(-10 - 1/1048576.0, 5);
  f := 9765639.99059019234;
  testrel(42, NE, y, f, cnt,failed);

  y := igammatx(-10 - 1/1048576.0, 20);
  f := 10240029255239.77433;
  testrel(43, NE, y, f, cnt,failed);

  y := igammatx(10, 12000);
  f := 1.615055828898457214e-41;
  testrel(44, NE, y, f, cnt,failed);

  y := igammatx(-99.5, 1e-15);
  f := -0.2981864024908420834e155;
  testrel(45, NE, y, f, cnt,failed);

  y := igammatx(-99.5, 1e-9);
  f := -0.2981864021896287091e155;
  testrel(46, NE, y, f, cnt,failed);

  y := igammatx(-99.5, 0.25);
  f := -0.2316399004959906513e155;
  testrel(47, NE, y, f, cnt,failed);

  y := igammatx(-99.5, 1);
  f := -0.1085942807021431442e155;
  testrel(48, NE, y, f, cnt,failed);

  y := igammatx(-1.5, 10);
  f := 31.62277504264856436;
  testrel(49, NE, y, f, cnt,failed);

  y := igammatx(-1.5, 10);
  f := 31.62277504264856436;
  testrel(50, NE, y, f, cnt,failed);

  y := igammatx(-100-1/1024,100);
  f := 1.004507364254462516e200;
  testrel(51, NE, y, f, cnt,failed);


  {-----------  zeroes -----------------}
  y := igammatx(-1.5,0.2920206138896944060);
  f := 5.314588259549414770e-20;
  testabs(52, NA, y, f, cnt,failed);

  y := igammatx(-Pi,0.54433923997351093287);
  f := -5.140208654063338670e-21;
  testabs(53, NA, y, f, cnt,failed);

  y := igammatx(-5.5,1.268373355325407569);
  f := -5.807062635698484658e-18;
  testabs(54, NA, y, f, cnt,failed);

  {-----------  x < 0 -----------------}
  y := igammatx(1.5,-2);
  f := 2.834828226167050019;
  testrel(55, NE, y, f, cnt,failed);

  y := igammatx(-1.5,-3);
  f := 3.216126974433846461;
  testrel(56, NE, y, f, cnt,failed);

  y := igammatx(1.5,-20);
  f := 0.266696136957315215745591863640e8;
  testrel(57, NE, y, f, cnt,failed);

  y := igammatx(-1.5,-20);
  f := 0.118445084627638054e8;
  testrel(58, NE, y, f, cnt,failed);

  y := igammatx(-2.5,-6.75);
  f := -284.1030440431725117;
  testrel(59, NE, y, f, cnt,failed);

  y := igammatx(-3,-6.75);
  f := -307.546875;
  testrel(60, NE, y, f, cnt,failed);

  y := igammatx(-4,-6.75);
  f := 2075.94140625;
  testrel(61, NE, y, f, cnt,failed);

  y := igammatx(-4.5,-6.75);
  f := -3683.869931833015079;
  testrel(62, NE, y, f, cnt,failed);

  y := igammatx(10,-234.5);
  f := 0.7865522244639469583e94;
  testrel(63, NE, y, f, cnt,failed);

  y := igammatx(-10,-234.5);
  f := 0.5028384401634493640e24;
  testrel(64, NE, y, f, cnt,failed);

  y := igammatx(-10.5,-234.5);
  f := -0.1180942036541155015e107;
  testrel(65, NE2, y, f, cnt,failed);

  {AMD is more inaccurate on exp(-414.905849945067989) and exp(585.094150054932011)}
  y := igammatx(300,-1000);
  f := 0.1486325482088438744e-180;
  testrel(66, NE2a, y, f, cnt,failed);

  y := igammatx(300,-2000);
  f := 0.1654696924976208728e254;
  testrel(67, NE2a, y, f, cnt,failed);

  {--------- extended ------------------}
  y := igammatx(-876.5,4321);
  f := 3.886954185910445972e3186;
  testrel(68, NE, y, f, cnt,failed);

  y := igammatx(-1234.25,1.5);
  f := 1.230663641811252244e3277;
  testrel(69, NE2, y, f, cnt,failed);

  y := igammatx(1000,4);
  f := 4.570005917775803474e-2570;
  testrel(70, NE, y, f, cnt,failed);

  y := igammatx(1000,1001);
  f := 1.902193310596208902e-3001;
  testrel(71, NE2, y, f, cnt,failed);

  y := igammatx(1000,998);
  f := 3.546213331125246836e-3000;
  testrel(72, NE, y, f, cnt,failed);

  y := igammatx(1700,0.125);
  f := 2.943488302318745557e-4756;
  testrel(73, NE, y, f, cnt,failed);

  y := igammatx(2000,-6000);
  f := 0.4407797138966629691e-3130;
  testrel(74, NE3, y, f, cnt,failed);

  y := igammatx(2000,-18000);
  f := 0.6026196015041197621e2081;
  testrel(75, NE3, y, f, cnt,failed);

  {----------------------------------------------}
  {extreme, must use log form}
  y := igammatx(-1755.5,487.9289845875);
  f := -0.874170972805294213426123048354e4710;
  testrele(76, 2e-6, y, f, cnt,failed);

  {debug warning "loss of accuracy"}
  y := igammatx(-1755.5,487.9289845876);
  f := 0.307889423407470407597774053272e4710;
  testrele(77, 5e-6, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_igamma_invx;
var
  a,x,p,q,z: extended;
  cnt, failed: integer;
const
  NE = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','igamma/p/q_invx');

  a := 0.125;
  x := 0.0009765625;
  q := 0.5535908078349959184;
  p := 0.4464091921650040816;
  z := igamma_invx(a,p,q);
  testrel( 1, 20, z, x, cnt,failed);  {!!!!!}

  a := 0.125;
  x := 0.1240234375;
  q := 0.1929079793425117863;
  p := 0.8070920206574882137;
  z := igamma_invx(a,p,q);
  testrel( 2, NE, z, x, cnt,failed);

  a := 0.125;
  x := 0.1259765625;
  q := 0.1914966349509656019;
  p := 0.8085033650490343981;
  z := igamma_invx(a,p,q);
  testrel( 3, NE, z, x, cnt,failed);

  a := 0.125;
  x := 1.125;
  q := 0.02536738859851095797;
  p := 0.9746326114014890420;
  z := igamma_invx(a,p,q);
  testrel( 4, NE, z, x, cnt,failed);

  a := 0.125;
  x := 0.25;
  q := 0.1303353454497136427;
  p := 0.8696646545502863573;
  z := igamma_invx(a,p,q);
  testrel( 5, NE, z, x, cnt,failed);

  a := 0.125;
  x := 100.0;
  q := 8.7052780108565328454E-47;
  p := 1.0;
  z := igamma_invx(a,p,q);
  testrel( 6, NE, z, x, cnt,failed);

  a := 0.125;
  x := 0.001;
  q := 0.5522666002838611327;
  p := 0.4477333997161388673;
  z := igamma_invx(a,p,q);
  testrel( 7, NE, z, x, cnt,failed);

  a := 0.5;
  x := 0.00390625;
  q := 0.9295680222776129220;
  p := 0.07043197772238707805;
  z := igamma_invx(a,p,q);
  testrel( 8, NE, z, x, cnt,failed);

  a := 0.5;
  x := 0.49609375;
  q := 0.3192083174860168816;
  p := 0.6807916825139831184;
  z := igamma_invx(a,p,q);
  testrel( 9, NE, z, x, cnt,failed);

  a := 0.5;
  x := 0.50390625;
  q := 0.3154274672237097826;
  p := 0.6845725327762902174;
  z := igamma_invx(a,p,q);
  testrel(10, NE, z, x, cnt,failed);

  a := 0.5;
  x := 1.5;
  q := 0.08326451666355040186;
  p := 0.9167354833364495982;
  z := igamma_invx(a,p,q);
  testrel(11, NE, z, x, cnt,failed);

  a := 0.5;
  x := 1.0;
  q := 0.1572992070502851307;
  p := 0.8427007929497148693;
  z := igamma_invx(a,p,q);
  testrel(12, NE, z, x, cnt,failed);

  a := 0.5;
  x := 100.0;
  q := 2.088487583762544757E-45;
  p := 1.0;
  z := igamma_invx(a,p,q);
  testrel(13, NE, z, x, cnt,failed);

  a := 0.5;
  x := 0.001;
  q := 0.9643294082703201150;
  p := 0.03567059172967988505;
  z := igamma_invx(a,p,q);
  testrel(14, NE, z, x, cnt,failed);

  a := 1.25;
  x := 0.009765625;
  q := 0.9973051206118684410;
  p := 0.002694879388131559014;
  z := igamma_invx(a,p,q);
  testrel(15, NE, z, x, cnt,failed);

  a := 1.25;
  x := 1.240234375;
  q := 0.3847479139371369129;
  p := 0.6152520860628630872;
  z := igamma_invx(a,p,q);
  testrel(16, NE, z, x, cnt,failed);

  a := 1.25;
  x := 1.259765625;
  q := 0.3782200413443559917;
  p := 0.6217799586556440083;
  z := igamma_invx(a,p,q);
  testrel(17, NE, z, x, cnt,failed);

  a := 1.25;
  x := 0.25;
  q := 0.8638836284585888426;
  p := 0.1361163715414111574;
  z := igamma_invx(a,p,q);
  testrel(18, NE, z, x, cnt,failed);

  a := 1.25;
  x := 2.25;
  q := 0.1549995477947637834;
  p := 0.8450004522052362166;
  z := igamma_invx(a,p,q);
  testrel(19, NE, z, x, cnt,failed);

  a := 1.25;
  x := 2.5;
  q := 0.1230885711526588452;
  p := 0.8769114288473411548;
  z := igamma_invx(a,p,q);
  testrel(20, NE, z, x, cnt,failed);

  a := 1.25;
  x := 100.0;
  q := 1.301089352558546598E-43;
  p := 1.0;
  z := igamma_invx(a,p,q);
  testrel(21, NE, z, x, cnt,failed);

  a := 1.25;
  x := 0.001;
  q := 0.999843134425254734;
  p := 0.00015686557474526599535;
  z := igamma_invx(a,p,q);
  testrel(22, NE, z, x, cnt,failed);

  a := 12.375;
  x := 0.0966796875;
  q := 1.0;
  p := 2.0443064129172721431E-22;
  z := igamma_invx(a,p,q);
  testrel(23, NE, z, x, cnt,failed);

  a := 12.375;
  x := 12.2783203125;
  q := 0.4731159047024541000;
  p := 0.5268840952975459000;
  z := igamma_invx(a,p,q);
  testrel(24, NE, z, x, cnt,failed);

  a := 12.375;
  x := 12.4716796875;
  q := 0.4513371545133265204;
  p := 0.5486628454866734796;
  z := igamma_invx(a,p,q);
  testrel(25, NE, z, x, cnt,failed);

  a := 12.375;
  x := 11.375;
  q := 0.5779789073056998708;
  p := 0.4220210926943001292;
  z := igamma_invx(a,p,q);
  testrel(26, NE, z, x, cnt,failed);

  a := 12.375;
  x := 13.375;
  q := 0.35521881768820689405;
  p := 0.64478118231179310595;
  z := igamma_invx(a,p,q);
  testrel(27, NE, z, x, cnt,failed);

  a := 12.375;
  x := 24.75;
  q := 0.0022231460245010735;
  p := 0.9977768539754989265;
{$ifdef BIT16}
  p := 1-q;
{$endif}
  z := igamma_invx(a,p,q);
  testrel(28, NE, z, x, cnt,failed);

  a := 12.375;
  x := 100.0;
  q := 2.3484886146232989951E-29;
  p := 1.0;
  z := igamma_invx(a,p,q);
  testrel(39, NE, z, x, cnt,failed);

  a := 12.375;
  x := 0.001;
  q := 1.0;
  p := 6.0319676270277627123E-47;
  z := igamma_invx(a,p,q);
  testrel(30, NE, z, x, cnt,failed);

  a := 24.75;
  x := 0.193359375;
  q := 1.0;
  p := 2.6139638700232537181E-43;
  z := igamma_invx(a,p,q);
  testrel(31, NE, z, x, cnt,failed);

  a := 24.75;
  x := 24.556640625;
  q := 0.4887746109512457205;
  p := 0.5112253890487542795;
  z := igamma_invx(a,p,q);
  testrel(32, NE, z, x, cnt,failed);

  a := 24.75;
  x := 24.943359375;
  q := 0.4578748584943997499;
  p := 0.5421251415056002501;
  z := igamma_invx(a,p,q);
  testrel(33, NE, z, x, cnt,failed);

  a := 24.75;
  x := 23.75;
  q := 0.5542816610107637155;
  p := 0.4457183389892362845;
  z := igamma_invx(a,p,q);
  testrel(34, NE, z, x, cnt,failed);

  a := 24.75;
  x := 25.75;
  q := 0.3954249727745276303;
  p := 0.6045750272254723697;
  z := igamma_invx(a,p,q);
  testrel(35, NE, z, x, cnt,failed);

  a := 24.75;
  x := 49.5;
  q := 0.00003746713560874641231;
  p := 0.9999625328643912536;
{$ifdef BIT16}
  p := 1-q;
{$endif}
  z := igamma_invx(a,p,q);
  testrel(36, NE, z, x, cnt,failed);

  a := 24.75;
  x := 100.0;
  q := 5.5034695332350412741E-20;
  p := 0.9999999999999999999;
{$ifdef BIT16}
  p := 1-q;
{$endif}
  z := igamma_invx(a,p,q);
  testrel(37, NE, z, x, cnt,failed);

  a := 24.75;
  x := 0.001;
  q := 1.0;
  p := 8.1291472534766856537E-100;
  z := igamma_invx(a,p,q);
  testrel(38, NE, z, x, cnt,failed);

  a := 245.675;
  x := 1.9193359375;
  q := 1.0;
  p := 3.8186081058738398329E-414;
  z := igamma_invx(a,p,q);
  testrel(39, NE, z, x, cnt,failed);

  a := 245.675;
  x := 243.7556640625;
  q := 0.5404196859985200346;
  p := 0.4595803140014799654;
  z := igamma_invx(a,p,q);
  testrel(40, NE, z, x, cnt,failed);

  a := 245.675;
  x := 247.5943359375;
  q := 0.4429907955159576511;
  p := 0.5570092044840423489;
  z := igamma_invx(a,p,q);
  testrel(41, NE, z, x, cnt,failed);

  a := 245.675;
  x := 244.675;
  q := 0.5169940631722782868;
  p := 0.4830059368277217132;
  z := igamma_invx(a,p,q);
  testrel(42, NE, z, x, cnt,failed);

  a := 245.675;
  x := 246.675;
  q := 0.4661406493259865244;
  p := 0.5338593506740134756;
  z := igamma_invx(a,p,q);
  testrel(43, NE, z, x, cnt,failed);

  a := 245.675;
  x := 491.35;
  q := 4.5956807163396664910E-35;
  p := 1.0;
  z := igamma_invx(a,p,q);
  testrel(44, NE, z, x, cnt,failed);

  a := 245.675;
  x := 100.0;
  q := 1.0;
  p := 9.8734129137566074939E-35;
  z := igamma_invx(a,p,q);
  testrel(45, NE, z, x, cnt,failed);

  a := 245.675;
  x := 0.001;
  q := 1.0;
  p := 6.6606493267933284947E-1220;
  z := igamma_invx(a,p,q);
  testrel(46, NE, z, x, cnt,failed);

  a := 2456.75;
  x := 19.193359375;
  q := 1.0;
  p := 4.3493575976177286679E-4121;
  z := igamma_invx(a,p,q);
  testrel(47, NE, z, x, cnt,failed);

  a := 2456.75;
  x := 2437.556640625;
  q := 0.64858936942781918995;
  p := 0.35141063057218081005;
  z := igamma_invx(a,p,q);
  testrel(48, NE, z, x, cnt,failed);

  a := 2456.75;
  x := 2475.943359375;
  q := 0.3471788337706876469;
  p := 0.6528211662293123531;
  z := igamma_invx(a,p,q);
  testrel(49, NE, z, x, cnt,failed);

  a := 2456.75;
  x := 2455.75;
  q := 0.5053666606373761599;
  p := 0.4946333393626238401;
  z := igamma_invx(a,p,q);
  testrel(50, NE, z, x, cnt,failed);

  a := 2456.75;
  x := 2457.75;
  q := 0.4892707552751261227;
  p := 0.5107292447248738773;
  z := igamma_invx(a,p,q);
  testrel(51, NE, z, x, cnt,failed);

  a := 2456.75;
  x := 4913.5;
  q := 3.2198475192983317356E-330;
  p := 1.0;
  z := igamma_invx(a,p,q);
  testrel(52, NE, z, x, cnt,failed);

  a := 2456.75;
  x := 100.0;
  q := 1.0;
  p := 4.7640315099546672134E-2395;
  z := igamma_invx(a,p,q);
  testrel(53, 12, z, x, cnt,failed);  {!!!!!}

  a := 1000.0;
  x := 1000.0+1/1024.0;
  q := 0.4957824368694063252;
  p := 0.5042175631305936748;
  z := igamma_invx(a,p,q);
  testrel(54, NE, z, x, cnt,failed);

  a := 0.125;
  x := 0.6766354068395209429;
  q := 0.053;
  p := 0.947;
  z := igamma_invx(a,p,q);
  testrel(55, NE, z, x, cnt,failed);

  a := 0.125;
  x := 1.4360412743024867748;
  q := 0.016;
  p := 0.984;
  z := igamma_invx(a,p,q);
  testrel(56, NE, z, x, cnt,failed);


  {test cases erf/erfc and igammap/q_invx}

  a := 0.5;
  x := ldexp(1,-20);
  p := 0.0011019324300718147042;
  z := igammap_invx(a,p);
  testrel(57, NE, z, x, cnt,failed);

  x := 0.0078125;
  p := 0.099476449660225785959;
  z := igammap_invx(a,p);
  testrel(58, NE, z, x, cnt,failed);

  x := 0.125;
  p := 0.38292492254802620728;
  z := igammap_invx(a,p);
  testrel(59, NE, z, x, cnt,failed);

  x := 0.5;
  q := 0.31731050786291410283;
  z := igammaq_invx(a,q);
  testrel(60, NE, z, x, cnt,failed);

  x := 0.9;
  q := 0.17971249487899984212;
  z := igammaq_invx(a,q);
  testrel(61, NE, z, x, cnt,failed);

  x := 1.0;
  q := 0.15729920705028513066;
  z := igammaq_invx(a,q);
  testrel(62, NE, z, x, cnt,failed);

  x := 1.1;
  q := 0.13801073756865955964;
  z := igammaq_invx(a,q);
  testrel(63, NE, z, x, cnt,failed);

  x := 2.0;
  q := 0.045500263896358414401;
  z := igammaq_invx(a,q);
  testrel(64, NE, z, x, cnt,failed);

  x := 4.0;
  q := 0.0046777349810472658379;
  z := igammaq_invx(a,q);
  testrel(65, NE, z, x, cnt,failed);

  x := 6.0;
  q := 0.00053200550513924969929;
  z := igammaq_invx(a,q);
  testrel(66, NE, z, x, cnt,failed);

  x := 10.0;
  q := 0.0000077442164310440836377;
  z := igammaq_invx(a,q);
  testrel(67, NE, z, x, cnt,failed);

  x := 40.0;
  q := 3.7440973842028987635e-19;
  z := igammaq_invx(a,q);
  testrel(68, NE, z, x, cnt,failed);

  x := 100.0;
  q := 2.0884875837625447570E-45;
  z := igammaq_invx(a,q);
  testrel(69, NE, z, x, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_igammap_derx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 8;
  NE2 = 32;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','igammap_derx');

  y := igammap_derx(-0.5, 0.25);
  f := -1.757565157870889588;
  testrel(1, NE, y, f, cnt,failed);

  y := igammap_derx(-0.5, 1);
  f := -0.1037768743551486758;
  testrel(2, NE1, y, f, cnt,failed);

  y := igammap_derx(-0.5, 10);
  f := -0.4049955478044558679e-6;
  testrel(3, NE, y, f, cnt,failed);

  y := igammap_derx(-0.5, 100);
  f := -0.1049414057838604223e-46;
  testrel(4, NE1, y, f, cnt,failed);

  y := igammap_derx(0.125, 0.25);
  f := 0.3477015467101832473;
  testrel(5, NE, y, f, cnt,failed);

  y := igammap_derx(0.125, 1);
  f := 0.4882961147855917348e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := igammap_derx(0.125, 10);
  f := 0.8035870541742101302e-6;
  testrel(7, NE, y, f, cnt,failed);

  y := igammap_derx(0.125, 100);
  f := 0.8780708511191650384e-46;
  testrel(8, NE, y, f, cnt,failed);

  y := igammap_derx(2.5, 0.25);
  f := 0.7323188157795373284e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := igammap_derx(2.5, 1);
  f := 0.2767383316137298022;
  testrel(10, NE, y, f, cnt,failed);

  y := igammap_derx(2.5, 10);
  f := 0.1079988127478548981e-2;
  testrel(11, NE1, y, f, cnt,failed);

  y := igammap_derx(2.5, 100);
  f := 0.2798437487569611260e-40;
  testrel(12, NE, y, f, cnt,failed);

  y := igammap_derx(2.5, 500);
  f := 0.5992083479155489790e-213;
  testrel(13, NE1, y, f, cnt,failed);  {FPC < 3}

  y := igammap_derx(100, 1);
  f := 0.3941866060050479210e-156;
  testrel(14, NE2, y, f, cnt,failed);

  y := igammap_derx(100, 10);
  f := 0.4864649182067610436e-61;
  testrel(15, NE1, y, f, cnt,failed);

  y := igammap_derx(100, 100);
  f := 0.3986099680914713523e-1;
  testrel(16, NE1, y, f, cnt,failed);

  y := igammap_derx(100, 1000);
  f := 0.5438942180826245858e-293;
  testrel(17, NE2, y, f, cnt,failed);

  y := igammap_derx(-1.5, 2);
  f := 0.1012330622122275974e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := igammap_derx(1, 0);
  f := 1.0;
  testrel(19, NE, y, f, cnt,failed);

  y := igammap_derx(1.5, 0);
  f := 0;
  testrel(20, NE, y, f, cnt,failed);

  y := igammap_derx(2, 0);
  f := 0;
  testrel(21, NE, y, f, cnt,failed);

  y := igammap_derx(0.5, 1e-10);
  f := 56418.95834913373286;
  testrel(22, NE, y, f, cnt,failed);

  y := igammap_derx(0.5, 1e-20);
  f := 5641895835.477562869;
  testrel(23, NE, y, f, cnt,failed);

  y := igammap_derx(-0.5, 1e-10);
  f := -0.2820947917456686643e15;
  testrel(24, NE, y, f, cnt,failed);

  y := igammap_derx(-0.75, 2);
  f := -0.8323169367186116874e-2;
  testrel(25, NE, y, f, cnt,failed);

  y := igammap_derx(1e-10, 2);
  f := 0.6766764162690259172e-11;
  testrel(26, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_igammanx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','igammanx');

  y := igammax(2,0);
  f := 1;
  testrel(1, NE, y, f, cnt,failed);

  y := igammax(2,-1);
  f := 0;
  testrel(2, NE, y, f, cnt,failed);

  y := igammax(-1,2.5);
  f := 0.7919081579289782572e-2;
  testrel(3, NE, y, f, cnt,failed);

  y := igammax(0,2.5);
  f := 0.24914917870269735496e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := igammax(0,-2.5);
  f := -7.073765894578600712;
  testrel(5, NE, y, f, cnt,failed);

  y := igammax(0,-100);
  f := -0.2715552744853879822e42;
  testrel(6, NE, y, f, cnt,failed);

  y := igammax(5,2.5);
  f := 21.38827245393962982;
  testrel(7, NE1, y, f, cnt,failed);

  y := igammax(5,-2.5);
  f := 189.5900622634478054;
  testrel(8, NE1, y, f, cnt,failed);

  y := igammax(5,-10);
  f := 153832837.1109301082;
  testrel(9, NE1, y, f, cnt,failed);

  y := igammax(5,-100);
  f := 0.2583754327050379842e52;
  testrel(10, NE1, y, f, cnt,failed);   {FPC < 3}

  y := igammax(100,-10);
{$ifdef BIT16}
  f := 0.18665243088788830536e157 * 0.5;
{$else}
  f := 0.9332621544394415268e156;
{$endif}
  testrel(11, NE1, y, f, cnt,failed);

  y := igammax(100,-100);
  f := -0.1347427096011818133e242;
  testrel(12, NE1, y, f, cnt,failed);   {FPC < 3}

  {--------- extended -------}
  y := igammax(5,-1000);
  f := 0.1962214423179921981e447;
  testrel(13, NE, y, f, cnt,failed);

  y := igammax(5,-10000);
  f := 8.803296554979500997e4358;
  testrel(14, NE, y, f, cnt,failed);

  y := igammax(200,-10);
  f := 0.3943289336823952518e373;
  testrel(15, NE, y, f, cnt,failed);

  y := igammax(500,-10);
  f := 2.440273651982220137e1131;
  testrel(16, NE1, y, f, cnt,failed);    {FPC < 3}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;

end.
