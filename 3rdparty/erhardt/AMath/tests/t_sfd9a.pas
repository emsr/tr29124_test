{Part 9a of regression test for SPECFUN unit  (c) 2013  W.Ehrhardt}

unit t_sfd9a;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_hyperg_u;
procedure test_whittaker;
procedure test_hyperg_0F1;
procedure test_hyperg_0F1r;

procedure test_cylinderd;
procedure test_cylinderu;
procedure test_cylinderv;
procedure test_hermiteh;

implementation

uses
  AMath,SpecFun,t_sfd0;

{---------------------------------------------------------------------------}
procedure test_hyperg_u;
const
  NE = 1;
  NE1 = 4;
  NE2 = 32;
  NE3 = 128;
var
  a,b,x,y,f: double;
  cnt, failed: integer;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_u');

  y := hyperg_u(2.5, 1000.5, 300);
  f := 0.1425867116494189890e218;
  testrel(1, NE, y, f, cnt,failed);

  y := hyperg_u(-2.5, 29000.5, 24000);
  f := -0.5886530457853419985e213;
  testrel(2, NE, y, f, cnt,failed);

  y := hyperg_u(1, -10, 10);
  f := 0.4656199948873187212e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := hyperg_u(12, 12, 10);
  f := 0.4656199948873187212e-12;
  testrel(4, NE, y, f, cnt,failed);

  y := hyperg_u(12, 12, 2);
  f := 0.3709781774548841033e-4;
  testrel(5, NE, y, f, cnt,failed);

  y := hyperg_u(12.5, 12.5, 10);
  f := 0.1439616407106279977e-12;
  testrel(6, NE, y, f, cnt,failed);

  y := hyperg_u(12.5, 12.5, 3);
  f := 0.2214855943221667552e-6;
  testrel(7, NE, y, f, cnt,failed);

  y := hyperg_u(-2.5, -1.5, 1.25);
  f := 1.746928107421710700;
  testrel(8, NE, y, f, cnt,failed);

  y := hyperg_u(-2.5, -0.5, 1.25);
  f := -1.746928107421710700;
  testrel(9, NE, y, f, cnt,failed);

  y := hyperg_u(-2.5, 0.5, 1.25);
  f := -1.048156864453026420;
  testrel(10, NE, y, f, cnt,failed);

  y := hyperg_u(-22.5, -10.25, 0.01);
  f := 0.1111201146250259309e15;
  testrel(11, NE, y, f, cnt,failed);

  y := hyperg_u(-22.5, -10.5, 0.01);
  f := -0.3883490058589971964e-9;
  testrel(12, NE1, y, f, cnt,failed);

  y := hyperg_u(-50,-100,0.01);
  f := 0.3083899381146898393e94;
  testrel(13, NE, y, f, cnt,failed);

  y := hyperg_u(-50,-100,200);
  f := 0.2849221879469261031e121;
  testrel(14, NE, y, f, cnt,failed);

  y := hyperg_u(1.5,1.25,4);
  f := 0.08877058989167655998;
  testrel(15, NE, y, f, cnt,failed);

  y := hyperg_u(1.5,1.25, 0.25);
  f := 1.347631561941240031;
  testrel(16, NE, y, f, cnt,failed);

  y := hyperg_u(1.5,1.25, 0.125);
  f := 2.156837399109530982;
  testrel(17, NE, y, f, cnt,failed);

  y := hyperg_u(15, 20, 4);
  f := 0.6033806130290031433e-6;
  testrel(18, NE, y, f, cnt,failed);

  y := hyperg_u(-1.5,-1.25,100);
  f := 1011.298996263445510;
  testrel(19, NE, y, f, cnt,failed);

  y := hyperg_u(-1.5,-1.25,1e10);
  f := 1000000000112500.000;
  testrel(20, NE, y, f, cnt,failed);

  y := hyperg_u(-1.5,-0.5,2);
  f := 2.828427124746190098;
  testrel(21, NE, y, f, cnt,failed);

  y := hyperg_u(150,200,120);
  f := 0.2796726136983788707e-293;
  testrel(22, NE, y, f, cnt,failed);

  y := hyperg_u(1e-16,1,10);
  f := 0.9999999999999997697;
  testrel(23, NE, y, f, cnt,failed);

  {Luke}
  y := hyperg_u(1.125,2.5,10);
  f := 0.7797136719091902962e-1;
  testrel(24, NE, y, f, cnt,failed);

  y := hyperg_u(1.125,-2,10);
  f := 0.5201024170621504198e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := hyperg_u(0.25,0.5,10);
  f := 0.5527654093095883258;
  testrel(26, NE, y, f, cnt,failed);

  y := hyperg_u(-0.25,-0.5,10);
  f := 1.829860857257179893;
  testrel(27, NE, y, f, cnt,failed);

  {b = 2a}
  y := hyperg_u(10,20,5);
  f := 0.9397105536532480000e-2;
  testrel(28, NE, y, f, cnt,failed);

  y := hyperg_u(5.125,10.25,15);
  f := 0.3433682497121066193e-5;
  testrel(29, NE, y, f, cnt,failed);

  y := hyperg_u(5,10,15);
  f := 0.4497655845145557080e-5;
  testrel(30, NE, y, f, cnt,failed);

  {both a+1-b and x small}
  a := 2.5;
  b := (a+1)-ldexp(1,-20);
  x := ldexp(1,-30);
  y := hyperg_u(a,b,x);
  f := 3.777815733924505203e22;
  testrel(31, NE, y, f, cnt,failed);

  a := -2.5;
  b := (a+1)-ldexp(1,-20);
  x := ldexp(1,-30);
  y := hyperg_u(a,b,x);
  f := 0.1267759335942076171e-5;
  testrel(32, NE, y, f, cnt,failed);

  a := 2.5;
  b := (a+1)-ldexp(1,-40);
  x := ldexp(1,-30);
  y := hyperg_u(a,b,x);
  f := 3.777893186221851076e22;
  testrel(33, NE, y, f, cnt,failed);

  a := -2.5;
  b := (a+1)-ldexp(1,-40);
  x := ldexp(1,-30);
  y := hyperg_u(a,b,x);
  f := 0.1209028041806229080e-11;
  testrel(34, NE, y, f, cnt,failed);

  {test for x=0}
  y := hyperg_u(1,0.5,0);
  f := 2;
  testrel(35, NE, y, f, cnt,failed);

  y := hyperg_u(1,0,0);
  f := 1;
  testrel(36, NE, y, f, cnt,failed);

  y := hyperg_u(1,-1.5,0);
  f := 0.4;
  testrel(37, NE, y, f, cnt,failed);

  y := hyperg_u(-2,-1.5,0);
  f := 0.75;
  testrel(38, NE, y, f, cnt,failed);

  y := hyperg_u(-2.5,-1.5,0);
  f := 0;
  testrel(39, NE, y, f, cnt,failed);

  {----------------------------}
  {Test cases borrowed from GSL}
  {----------------------------}
  y := hyperg_u(1, 1, 0.0001);
  f := 8.634088070212725330;
  testrel(40, NE, y, f, cnt,failed);

  y := hyperg_u(1, 1, 0.01);
  f := 4.078511443456425847;
  testrel(41, NE, y, f, cnt,failed);

  y := hyperg_u(1, 1, 0.5);
  f := 0.9229106324837304688;
  testrel(42, NE, y, f, cnt,failed);

  y := hyperg_u(1, 1, 2.0);
  f := 0.3613286168882225847;
  testrel(43, NE, y, f, cnt,failed);

  y := hyperg_u(1, 1, 1000);
  f := 0.0009990019940238807150;
  testrel(44, NE, y, f, cnt,failed);

  y := hyperg_u(1, 8, 0.01);
  f := 7.272361203006010000e+16;
  testrel(45, NE, y, f, cnt,failed);

  y := hyperg_u(1, 8, 1);
  f := 1957.0;
  testrel(46, NE, y, f, cnt,failed);

  y := hyperg_u(1, 8, 5);
  f := 1.042496;
  testrel(47, NE, y, f, cnt,failed);

  y := hyperg_u(1, 8, 50);
  f := 0.022660399001600000000;
  testrel(48, NE, y, f, cnt,failed);

  y := hyperg_u(1, 8, 1000);
  f := 0.0010060301203607207200;
  testrel(49, NE, y, f, cnt,failed);

  y := hyperg_u(1, 20, 1);
  f := 1.7403456103284421000e+16;
  testrel(50, NE, y, f, cnt,failed);

  y := hyperg_u(1, 20, 20);
  f := 0.22597813610531052969;
  testrel(51, NE, y, f, cnt,failed);

  y := hyperg_u(1, 50, 1);
  f := 3.374452117521520758e+61;
  testrel(52, NE, y, f, cnt,failed);

  y := hyperg_u(1, 50, 50);
  f := 0.15394136814987651785;
  testrel(52, NE, y, f, cnt,failed);

  y := hyperg_u(1, 100, 0.1);
  f := 1.04183251719908528577e+253;
  testrel(54, NE2, y, f, cnt,failed);

  y := hyperg_u(1, 100, 1);
  f := 2.5624945006073464385e+154;
  testrel(55, NE, y, f, cnt,failed);

  y := hyperg_u(1, 100, 50);
  f := 3.0978624160896431391e+07;
  testrel(56, NE, y, f, cnt,failed);

  y := hyperg_u(1, 100, 1000);
  f := 0.0011085142546061528661;
  testrel(57, NE, y, f, cnt,failed);

  y := hyperg_u(1, 1000, 2000);
  f := 0.0009970168547036318206;
  testrel(58, NE, y, f, cnt,failed);

  y := hyperg_u(1, -1, 1);
  f := 0.29817368116159703717;
  testrel(59, NE, y, f, cnt,failed);

  y := hyperg_u(1, -1, 10);
  f := 0.07816669698940409380;
  testrel(60, NE, y, f, cnt,failed);

  y := hyperg_u(1, -10, 1);
  f := 0.08271753756946041959;
  testrel(61, NE, y, f, cnt,failed);

  y := hyperg_u(1, -10, 5);
  f := 0.06127757419425055261;
  testrel(62, NE, y, f, cnt,failed);

  y := hyperg_u(1, -10, 20);
  f := 0.031606421847946077709;
  testrel(63, NE, y, f, cnt,failed);

  y := hyperg_u(1, -100, 0.01);
  f := 0.009900000099999796950;
  testrel(64, NE, y, f, cnt,failed);

  y := hyperg_u(1, -100, 1);
  f := 0.009802970197050404429;
  testrel(65, NE, y, f, cnt,failed);

  y := hyperg_u(1, -100, 10);
  f := 0.009001648897173103447;
  testrel(66, NE, y, f, cnt,failed);

  y := hyperg_u(1, -100, 20);
  f := 0.008253126487166557546;
  testrel(67, NE, y, f, cnt,failed);

  y := hyperg_u(1, -100, 90);
  f := 0.005222713769726871937;
  testrel(68, NE, y, f, cnt,failed);

  y := hyperg_u(1, -1000, 1);
  f := 0.0009980029970019970050;
  testrel(69, NE, y, f, cnt,failed);

  y := hyperg_u(1, -1000, 1010);
  f := 0.4971408839859245170e-3;
  testrel(70, NE, y, f, cnt,failed);

  y := hyperg_u(8, 1, 0.001);
  f := 0.0007505359326875706975;
  testrel(71, NE, y, f, cnt,failed);

  y := hyperg_u(8, 1, 0.5);
  f := 6.449509938973479986e-06;
  testrel(72, NE, y, f, cnt,failed);

  y := hyperg_u(8, 8, 1);
  f := 0.12289755012652317578;
  testrel(73, NE, y, f, cnt,failed);

  y := hyperg_u(8, 8, 10);
  f := 5.687710359507564272e-09;
  testrel(74, NE, y, f, cnt,failed);

  y := hyperg_u(8, 8, 20);
  f := 2.8175404594901039724e-11;
  testrel(75, NE, y, f, cnt,failed);

  y := hyperg_u(100, 100, 0.01);
  f := 1.0099979491941914867e+196;
  testrel(76, NE2, y, f, cnt,failed);

  y := hyperg_u(100, 100, 0.1);
  f := 1.0090713562719862834e+97;
  testrel(77, NE2, y, f, cnt,failed);

  y := hyperg_u(100, 100, 1);
  f := 0.009998990209084729106;
  testrel(78, NE, y, f, cnt,failed);

  y := hyperg_u(100, 100, 20);
  f := 1.3239363905866130603e-131;
  testrel(79, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 1, 0.01);
  f := 3.274012540759009536e+06;
  testrel(80, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 1, 1);
  f := 1.5202710000000000000e+06;
  testrel(81, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 1, 10);
  f := 1.0154880000000000000e+08;
  testrel(82, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 1, 100);
  f := 3.284529863685482880e+19;
  testrel(83, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 10, 1);
  f := 1.1043089864100000000e+11;
  testrel(84, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 100, 1);
  f := 1.3991152402448957897e+20;
  testrel(85, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 100, 10);
  f := 5.364469916567136000e+19;
  testrel(86, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 100, 100);
  f := 3.909797568000000000e+12;
  testrel(87, NE, y, f, cnt,failed);

  y := hyperg_u(-10, 100, 500);
  f := 8.082625576697984130e+25;
  testrel(88, NE, y, f, cnt,failed);

  y := hyperg_u(-50, 1, 0.01);
  f := 1.6973422555823855798e+64;
  testrel(89, NE, y, f, cnt,failed);

  y := hyperg_u(-50, 1, 1);
  f := 7.086160198304780325e+63;
  testrel(90, NE, y, f, cnt,failed);

  y := hyperg_u(-50, 1, 10);
  f := 5.332862895528712200e+65;
  testrel(91, NE, y, f, cnt,failed);

  y := hyperg_u(-50, 10, 1);
  f := -7.106713471565790573e+71;
  testrel(92, NE, y, f, cnt,failed);

  y := hyperg_u(-50, 100, 1);
  f := 2.4661377199407186476e+104;
  testrel(93, NE, y, f, cnt,failed);

  y := hyperg_u(-50, 10, 10);
  f := 5.687538583671241287e+68;
  testrel(94, NE, y, f, cnt,failed);

  y := hyperg_u(-50, 100, 10);
  f := 1.7880761664553373445e+102;
  testrel(95, NE, y, f, cnt,failed);

  y := hyperg_u(-90, 1, 0.01);
  f := 4.185245354032917715e+137;
  testrel(96, NE, y, f, cnt,failed);

  y := hyperg_u(-90, 1, 0.1);
  f := 2.423404340800784135e+137;
  testrel(97, NE, y, f, cnt,failed);

  y := hyperg_u(-90, 1, 10);
  f := -1.8987677149221888807e+139;
  testrel(98, NE, y, f, cnt,failed);

  y := hyperg_u(-90, 10, 10);
  f := -5.682999988842066677e+143;
  testrel(99, NE, y, f, cnt,failed);

  y := hyperg_u(-90, 100, 10);
  f := 2.3410029853990624280e+189;
  testrel(100, NE, y, f, cnt,failed);

  y := hyperg_u(-90, 1000, 10);
  f := 1.9799451517572225316e+271;
  testrel(101, NE, y, f, cnt,failed);

  y := hyperg_u(-50, -1, 10);
  f := -9.083195466262584149e+64;
  testrel(102, NE, y, f, cnt,failed);

  y := hyperg_u(-50, -10, 10);
  f := -1.4418257327071634407e+62;
  testrel(103, NE, y, f, cnt,failed);

  y := hyperg_u(-50, -100, 0.01);
  f := 3.0838993811468983931e+93;
  testrel(104, NE, y, f, cnt,failed);

  y := hyperg_u(-50, -100, 10);
  f := 4.014552630378340665e+95;
  testrel(105, NE, y, f, cnt,failed);

  y := hyperg_u(-100, -100, 10);
  f := 2.0556466922347982030e+162;
  testrel(106, NE, y, f, cnt,failed);

  y := hyperg_u(-100, -200, 10);
  f := 1.1778399522973555582e+219;
  testrel(107, NE, y, f, cnt,failed);

  y := hyperg_u(-100, -200, 100);
  f := 9.861313408898201873e+235;
  testrel(108, NE, y, f, cnt,failed);

  y := hyperg_u(-2.0, 4.0, 1.0);
  f := 11.0;
  testrel(109, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 0.0001, 0.0001);
  f := 1.0000576350699863577;
  testrel(110, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 0.0001, 1.0);
  f := 0.9999403679233247536;
  testrel(111, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 0.0001, 100.0);
  f := 0.9995385992657260887;
  testrel(112, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 1, 0.0001);
  f := 1.0009210608660065989;
  testrel(113, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 1.0, 1.0);
  f := 0.9999999925484179084;
  testrel(114, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 10, 1);
  f := 13.567851006281412726;
  testrel(115, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 100, 1);
  f := 2.5890615708804247881e+150;
  testrel(116, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 100, 10);
  f := 2.3127845417739661466e+55;
  testrel(117, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 100, 50);
  f := 6402.818715083582554;
  testrel(118, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 1000, 300);
  f := 2.5389557274938010716e+213;
  testrel(119, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 1000, 300);
  f := 1.1977955438214207486e+217;
  testrel(120, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 0.5, 1.0);
  f := 0.7578721561413121060;
  testrel(121, NE, y, f, cnt,failed);

  y := hyperg_u(1, 0.0001, 0.0001);
  f := 0.9992361337764090785;
  testrel(122, NE, y, f, cnt,failed);

  y := hyperg_u(1, 0.0001, 1);
  f := 0.4036664068111504538;
  testrel(123, NE, y, f, cnt,failed);

  y := hyperg_u(1, 0.0001, 100);
  f := 0.009805780851264329587;
  testrel(124, NE, y, f, cnt,failed);

  y := hyperg_u(1, 1.2, 2.0);
  f := 0.3835044780075602550;
  testrel(125, NE, y, f, cnt,failed);

  y := hyperg_u(1, -0.0001, 1);
  f := 0.4036388693605999482;
  testrel(126, NE, y, f, cnt,failed);

  y := hyperg_u(8, 10.5, 1);
  f := 27.981926466707438538;
  testrel(127, NE, y, f, cnt,failed);

  y := hyperg_u(8, 10.5, 100);
  f := 1.1226567526311488330e-16;
  testrel(128, NE, y, f, cnt,failed);

  y := hyperg_u(-2.0, 0.5, 3.14);
  f := 1.1896;
  testrel(129, NE1, y, f, cnt,failed);

  y := hyperg_u(-2.0, 0.5, 1.13);
  f := -1.3631;
  testrel(130, NE, y, f, cnt,failed);

  y := hyperg_u(-2.0, 0.5, 0.0);      {WE: via gamma_delta_ratio}
  f := 0.75;
  testrel(131, NE, y, f, cnt,failed);

  y := hyperg_u(-2.0, 0.5, 1e-20);
  f := 0.75;
  testrel(132, NE, y, f, cnt,failed);


  {--------------------------------------------}
  y := hyperg_u(-100, 0.5, 0.1);
  f := 0.5522118284507690095e157;
  testrel(133, NE, y, f, cnt,failed);

  y := hyperg_u(-100.1, 0.5, 0.1);
  f := 0.8463924522217037604e157;
  testrel(134, NE3, y, f, cnt,failed);   {!!!!}

  y := hyperg_u(-100.5, -110, 1);
  f := 0.3493093973937063178e173;
  testrel(135, NE, y, f, cnt,failed);

  y := hyperg_u(-100.5, -99.75, 1);
  f := 0.2227370160242474945e158;
  testrel(136, NE, y, f, cnt,failed);

  y := hyperg_u(-100.5, -99.75, 1);
  f := 0.2227370160242474945e158;
  testrel(137, NE, y, f, cnt,failed);

  {Special case a=-1}
  y := hyperg_u(-1, 2.25, 9.5);
  f := 7.25;
  testrel(138, NE, y, f, cnt,failed);

  {small x}
  y := hyperg_u(0.5, -0.5, 1e-9);
  f := 0.8862269245665732506;
  testrel(139, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 0, 1e-9);
  f := 1.128379155511377653;
  testrel(140, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 0.5, 1e-9);
  f := 1.772390607124724348;
  testrel(141, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 1, 1e-9);
  f := 12.14832450104719998;
  testrel(142, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 1.25, 1e-9);
  f := 362.4006708432691124;
  testrel(143, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 1.25, 1e-12);
  f := 2044.179387746884078;
  testrel(144, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 2, 1e-9);
  f := 0.5641895899040133278e9;
  testrel(145, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 2.5, 1e-9);
  f := 0.1581138833246467326e14;
  testrel(146, NE, y, f, cnt,failed);

  y := hyperg_u(0.5, 2.75, 1e-9);
  f := 0.2915883197062265245e16;
  testrel(147, NE, y, f, cnt,failed);

  {----------------------------------------}
  y := hyperg_u(0.1,0.1,0.1);
  f := 1.033514114257088719;
  testrel(148, NE, y, f, cnt,failed);

  y := hyperg_u(0.1,1.1,0.1);
  f := 1.258925411794167210;
  testrel(149, NE, y, f, cnt,failed);

  y := hyperg_u(0.1,1.2,0.1);
  f := 1.304388497882194119;
  testrel(150, NE, y, f, cnt,failed);

  y := hyperg_u(1.1,1.2,0.1);
  f := 2.584232565731920044;
  testrel(151, NE, y, f, cnt,failed);

  y := hyperg_u(1.0001,0.0001,0.0001);
  f := 0.9991938118206768393;
  testrel(152, NE, y, f, cnt,failed);

  y := hyperg_u(1.0001,1.0001,0.0001);
  f := 8.638232493095183934;
  testrel(153, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001,1.0001,0.0001);
  f := 1.000921458319295876;
  testrel(154, NE, y, f, cnt,failed);

  {some random 'cases'}
  y := hyperg_u(-10.5, 3.75, 42);
  f := 0.5385731532368236301e15;
  testrel(155, NE, y, f, cnt,failed);

  y := hyperg_u(4.25, -1.5, 0.875);
  f := 0.7385778980547234118e-3;
  testrel(156, NE, y, f, cnt,failed);

  y := hyperg_u( 5, 3.5, 17);
  f := 0.3852636454942120639e-6;
  testrel(157, NE, y, f, cnt,failed);

  y := hyperg_u( 0.25, -2, 8);
  f := 0.5480447978459433088;
  testrel(158, NE, y, f, cnt,failed);

  y := hyperg_u(-64.125, 12, 38);
  f := 0.2684017814775656072e98;
  testrel(159, NE, y, f, cnt,failed);

(*
  ----------------------------------------
  Alpha version problem test cases from GSL.
  Fixed with in Version 1.19.03.
*)
  y := hyperg_u(8, 1, 8);
  f := 6.190694573035761284e-10;
  testrel(160, NE, y, f, cnt,failed);

  y := hyperg_u(8, 1, 20);
  f := 3.647213845460374016e-12;
  testrel(161, NE, y, f, cnt,failed);

  y := hyperg_u(8, 1, 8);
  f := 6.190694573035761284e-10;
  testrel(162, NE, y, f, cnt,failed);

  y := hyperg_u(8, 1, 20);
  f := 3.647213845460374016e-12;
  testrel(163, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 10, 5);
  f := 1.0006265020064596364;
  testrel(164, NE, y, f, cnt,failed);

  y := hyperg_u(0.0001, 10, 10);
  f := 0.9999244381454633265;
  testrel(165, NE, y, f, cnt,failed);

  y := hyperg_u(8, 10.5, 10);
  f := 2.4370135607662056809e-8;
  testrel(166, NE, y, f, cnt,failed);

  y := hyperg_u(10, -2.5, 10);
  f := 6.734690720346560349e-14;
  testrel(167, NE, y, f, cnt,failed);

  y := hyperg_u(10, 2.5, 10);
  f := 6.787780794037971638e-13;
  testrel(168, NE, y, f, cnt,failed);

  y := hyperg_u(10, 2.5, 50);
  f := 2.409872007659608713e-18;
  testrel(169, NE, y, f, cnt,failed);

  y := hyperg_u(-10.5, 1.1, 99);
  f := 2.546328328942063558e20;
  testrel(170, NE, y, f, cnt,failed);

  y := hyperg_u(-50.5, 100.1, 10);
  f := -5.5740216266953697786e126;
  testrel(171, NE3, y, f, cnt,failed);  {!!!!}

  y := hyperg_u(-50.5, 100.1, 40);
  f := 5.937463786613929435e91;
  testrel(172, NE2, y, f, cnt,failed);

  {-----------------------------------}
  {NEW problem cases:                 }

  {Temme underflow and dchu very inaccurate with extended}
  y := hyperg_u(1, -2000, 1);
  f := 0.4995003748125624532e-3;
  testrel(173, NE1, y, f, cnt,failed);


  {-----------------------------------}
  {Negative x}
  y := hyperg_u(-3,2.5,-1.5);
  f := -144.0;
  testrel(174, NE, y, f, cnt,failed);

  y := hyperg_u(-4,-3,-4.5);
  f := 410.0625;
  testrel(175, NE, y, f, cnt,failed);

  y := hyperg_u(-4,-4,-15);
  f := 39489.0;
  testrel(176, NE, y, f, cnt,failed);

  y := hyperg_u(-4,-4.5,-15);
  f := 34709.0625;
  testrel(177, NE, y, f, cnt,failed);

  y := hyperg_u(-4, 4.5,-15);
  f := 234981.5625;
  testrel(178, NE, y, f, cnt,failed);

  {GSL}
  y := hyperg_u(3,4,-2);
  f := -0.125;
  testrel(179, NE, y, f, cnt,failed);

  y := hyperg_u(3,6,-0.5);
  f := -296.0;
  testrel(180, NE, y, f, cnt,failed);

  y := hyperg_u(1, 2, -0.1);
  f := -10.0;
  testrel(181, NE, y, f, cnt,failed);

  y := hyperg_u(-5, -2, -0.1);
  f := -0.02101;
  testrel(182, NE, y, f, cnt,failed);

  y := hyperg_u(-2, 0, -0.5);
  f := 1.25;
  testrel(183, NE, y, f, cnt,failed);

  y := hyperg_u(-8, 1.5, -0.5);
  f := 827341.0625;
  testrel(184, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_whittaker;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','WhittakerM/W');

  y := WhittakerM(1,2,0.5);
  f := 0.1606687379402938379;
  testrel(1, NE, y, f, cnt,failed);

  y := WhittakerM(0.25,-0.25, 5);
  f := 0.1227457026487926287;
  testrel(2, NE, y, f, cnt,failed);

  y := WhittakerM(0.25,-0.75, 5);
  f := 8.146924726654260291;
  testrel(3, NE, y, f, cnt,failed);

  y := WhittakerM(10.5, -12, 2);
  f := 0.8569214036105368264e-3;
  testrel(4, NE, y, f, cnt,failed);

  y := WhittakerM(1, -2.5, -3);
  f := -0.6224568153247312254e-1;
  testrel(5, NE, y, f, cnt,failed);

  y := WhittakerM(-1, -3.5, -4);
  f := -0.2309080030915828196e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := WhittakerM(-2, 1.5, -5);
  f := 2.052124965597469879;
  testrel(7, NE, y, f, cnt,failed);

  y := WhittakerM( 1.750, -2.125, 8.5);
  f := 2.829796856385884167;
  testrel(8, NE, y, f, cnt,failed);

  y := WhittakerM( 1.750, -0.875, 2.5);
  f := -5.549545003211265027;
  testrel(9, NE, y, f, cnt,failed);

  y := WhittakerM(-2.500, -2.375, 9.0);
  f := 3496.785910216931565;
  testrel(10, NE, y, f, cnt,failed);

  y := WhittakerM(-1.750, -0.250, 0.5);
  f := 2.830026884975734590;
  testrel(11, NE, y, f, cnt,failed);

  y := WhittakerM( 1.500,  1.375, 7.5);
  f := 6.424703378315435298;
  testrel(12, NE, y, f, cnt,failed);

  y := WhittakerM(-2.375, -0.750, 3.5);
  f := -704.3776507705635371;
  testrel(13, NE, y, f, cnt,failed);

  y := WhittakerM( 1.625,  1.625, 7.0);
  f := 10.15549651136315476;
  testrel(14, NE, y, f, cnt,failed);

  y := WhittakerM(-1.000, -1.875, 5.5);
  f := 12.79759428131255171;
  testrel(15, NE, y, f, cnt,failed);

  y := WhittakerM( 0.875, 2.500, -4.5);
  f := -234.7374288027617797;
  testrel(16, NE, y, f, cnt,failed);

  y := WhittakerM(-0.250,  0.500, 2.0);
  f := 2.955279926427901456;
  testrel(17, NE, y, f, cnt,failed);

  y := WhittakerW(1,2,0.5);
  f := 22.37777237182168845;
  testrel(18, NE, y, f, cnt,failed);

  y := WhittakerW(0.125,-0.375, 5);
  f := 0.1003771927930624768;
  testrel(19, NE, y, f, cnt,failed);

  y := WhittakerW(0.125,0.375, 5);
  f := 0.1003771927930624768;
  testrel(20, NE, y, f, cnt,failed);

  y := WhittakerW(0.125,0.375, 0.5);
  f := 0.7141634669274423374;
  testrel(21, NE, y, f, cnt,failed);

  y := WhittakerW(12.5, -12, 2);
  f := 2130.985349213752919;
  testrel(22, NE, y, f, cnt,failed);

  y := WhittakerW( 1.750, -2.125, 8.5);
  f := 0.8770608373367219502;
  testrel(23, NE, y, f, cnt,failed);

  y := WhittakerW( 1.750, -0.875, 2.5);
  f := 0.9048504560799111236;
  testrel(24, NE, y, f, cnt,failed);

  y := WhittakerW(-2.500, -2.375, 9.0);
  f := 0.3454177176066522701e-4;
  testrel(25, NE, y, f, cnt,failed);

  y := WhittakerW(-1.750, -0.500, 0.5);
  f := 0.185564017849227978;
  testrel(26, NE, y, f, cnt,failed);

  y := WhittakerW( 1.500,  1.375, 7.5);
  f := 0.5479055675519496603;
  testrel(27, NE, y, f, cnt,failed);

  y := WhittakerW(-2.375, -0.750, 3.5);
  f := 0.2408548565266728948e-2;
  testrel(28, NE, y, f, cnt,failed);

  y := WhittakerW( 1.625,  1.625, 7.0);
  f := 0.8819274542529732877;
  testrel(29, NE, y, f, cnt,failed);

  y := WhittakerW(-1.000, -1.875, 5.5);
  f := 0.1384268664814007373e-1;
  testrel(30, NE, y, f, cnt,failed);

  y := WhittakerW( 0.875, -2.000, 4.5);
  f := 0.8772280811512047670;
  testrel(31, NE, y, f, cnt,failed);

  y := WhittakerW(-0.250,  0.500, 2.0);
  f := 0.2780860689367105241;
  testrel(32, NE, y, f, cnt,failed);

  {--------------------------------------------}
  y := WhittakerM(1,2,1000);
  f := 0.3794422724775862325e216;
  testrel(33, NE, y, f, cnt,failed);

  y := WhittakerM(2,-0.25, 1000);
  f := 9.062112037420482495E210;
  testrel(34, NE, y, f, cnt,failed);

  y := WhittakerW(1,2, 1000);
  f := 0.7151343692149910402e-214;
  testrel(35, NE, y, f, cnt,failed);

  y := WhittakerW(2,-0.25, 1000);
  f := 7.108992856855041338E-212;
  testrel(36, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_hyperg_0F1;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_0F1');

  y := hyperg_0F1(0.0, 0.0);
  f := 1.0;
  testrel(1, NE, y, f, cnt,failed);

  y := hyperg_0F1(1, 0.5);
  f := 1.566082929756350537;
  testrel(2, NE, y, f, cnt,failed);

  y := hyperg_0F1(5, 0.5);
  f := 1.104267440482868457;
  testrel(3, NE, y, f, cnt,failed);

  y := hyperg_0F1(100, 30);
  f := 1.349259863948511018;
  testrel(4, NE, y, f, cnt,failed);

  y := hyperg_0F1(-0.5, 3);
  f := -39.29137997543434276;
  testrel(5, NE, y, f, cnt,failed);

  y := hyperg_0F1(-100.5, 50);   {CF better}
  f := 0.6087930289227538496;
  testrel(6, NE, y, f, cnt,failed);

  y := hyperg_0F1(-50.5, 50);    {CF better}
  f := 0.3751397961644844952;
  testrel(7, NE, y, f, cnt,failed);

  y := hyperg_0F1(-10.5, 50);
  f := -0.2333238666215689656e7;
  testrel(8, NE, y, f, cnt,failed);

  y := hyperg_0F1(-1.5, 50);
  f := 0.3710003667607782601e8;
  testrel(9, NE, y, f, cnt,failed);

  y := hyperg_0F1(1, -5.0);
  f := -0.3268752818235339109;
  testrel(10, NE, y, f, cnt,failed);

  y := hyperg_0F1(-0.5, -5.0);
  f := -4.581634759005381184;
  testrel(11, NE, y, f, cnt,failed);

  y := hyperg_0F1(-10.5, -555.0);
  f := 1.849189744827795158e8;
  testrel(12, NE, y, f, cnt,failed);

  y := hyperg_0F1(-10000.5, 50.0);
  f := 0.9950127291770990678;
  testrel(13, NE, y, f, cnt,failed);

  y := hyperg_0F1(10.5, 100000);
  f := 0.1398797368722340445e256;
  testrel(14, NE, y, f, cnt,failed);

  f := -0.3866144745970700341e102;
  y := hyperg_0F1(-10.5, 10000);
  testrel(15, NE, y, f, cnt,failed);

  y := hyperg_0F1(10000, 100);
  f := 1.010050162034428941;
  testrel(16, NE, y, f, cnt,failed);

  y := hyperg_0F1(-2000.5, -10000);      {CF better: < 2 eps}
  f := 149.1600197284238245;
  testrel(17, NE, y, f, cnt,failed);

  y := hyperg_0F1(-3000.5, 60);
  f := 0.9802020054319757207;
  testrel(18, NE, y, f, cnt,failed);

  y := hyperg_0F1(1000, 1e-15);
  f := 1.000000000000000001;
  testrel(19, NE, y, f, cnt,failed);

  y := hyperg_0F1(1, 1e-5);
  f := 1.000010000025000028;
  testrel(20, NE, y, f, cnt,failed);

  y := hyperg_0F1(1e-20, 1e-5);
  f := 0.1000005000008334340e16;
  testrel(21, NE, y, f, cnt,failed);

  y := hyperg_0F1(1e-20, 1e-4);
  f := 0.1000050000833340378e17;
  testrel(22, NE, y, f, cnt,failed);

  y := hyperg_0F1(1e-200, 1.5e-5);
  f := 0.1500011250028125035e196;
  testrel(23, NE, y, f, cnt,failed);

  y := hyperg_0F1(2000, 1);
  f := 1.000500124958335951;
  testrel(24, NE, y, f, cnt,failed);

  y := hyperg_0F1(-1000.5,0.25);
  f := 0.9997501561849311309;
  testrel(25, NE, y, f, cnt,failed);

  y := hyperg_0F1(-1000.5,-0.25);        {CF better}
  f := 1.000249906315084821;
  testrel(26, NE, y, f, cnt,failed);

  y := hyperg_0F1(1000.5,-0.25);
  f := 0.9997501561225403168;
  testrel(27, NE, y, f, cnt,failed);

  y := hyperg_0F1(1000.5,0.25);
  f := 1.000249906252662819;
  testrel(28, NE, y, f, cnt,failed);

  y := hyperg_0F1(-1000.5,0.125);
  f := 0.9998750702809415706;
  testrel(29, NE, y, f, cnt,failed);

  y := hyperg_0F1(-1000.5,-0.125);
  f := 1.00012494534406235607729185265;
  testrel(30, NE, y, f, cnt,failed);

  y := hyperg_0F1(1000.5,-0.125);
  f := 0.9998750702653419182;
  testrel(31, NE, y, f, cnt,failed);

  y := hyperg_0F1(1000.5,0.125);
  f := 1.000124945328458805;
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_hyperg_0F1r;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_0F1r');

  y := hyperg_0F1r(0, 0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := hyperg_0F1r(3, 0);
  f := 0.5;
  testrel(2, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-2.5, 0);
  f := -1.057855469152043038;
  testrel(3, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-2.5, -0.0625);
  f := -1.084875067753319410;
  testrel(4, NE, y, f, cnt,failed);

  y := hyperg_0F1r(0, 2);
  f := 4.789666198546809433;
  testrel(5, NE, y, f, cnt,failed);

  y := hyperg_0F1r(0, 0.1);
  f := 0.1050840312616016807;
  testrel(6, NE, y, f, cnt,failed);

  y := hyperg_0F1r(0, -0.1);
  f := -0.9508264234956454522e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-2, 10);
  f := 1339.263965463333012;
  testrel(8, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-10, 4);
  f := 0.1460392507048367632;
  testrel(9, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-10, -30);
  f := -0.2627191378120922238e8;
  testrel(10, NE, y, f, cnt,failed);

  y := hyperg_0F1r(170, 1);
  f := 0.2356251022619302925e-304;
  testrel(11, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-160.125, 1);
  f := -0.1076810946443946249e285;
  testrel(12, NE, y, f, cnt,failed);

  y := hyperg_0F1r(170, 2500);
  f := 0.3227168427177659957e-298;
  testrel(13, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-160, -120);        {DAMath: power overflow}
  f := -0.3511263902718893603e48;
  testrel(14, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-160.125, -0.25);   {DAMath: Jv overflow}
  f := -0.1085249736003132468e285;
  testrel(15, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-160.125, 0.1875);  {DAMath: Iv overflow}
  f := -0.1082288616245840387e285;
  testrel(16, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-1700, -700);
  f := -0.4222215366957806816e81;
  testrel(17, NE, y, f, cnt,failed);

  y := hyperg_0F1r(-1701, -800);        {power overflow}
  f := 0.8252189997101624829e179;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cylinderd;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','CylinderD');

  y := CylinderD(3.5, 1);
  f := -2.0368057121849972916;
  testrel(1, NE, y, f, cnt,failed);

  y := CylinderD(3.5, -1);
  f := 0.7877920070390601594;
  testrel(2, NE, y, f, cnt,failed);  {???, cancellation D(.,1) - 1F1 term}

  y := CylinderD(20, 20);
  f := 0.2368786440467208020e-17;
  testrel(3, NE, y, f, cnt,failed);

  y := CylinderD(20, -20);
  f := 0.2368786440467208020e-17;
  testrel(4, NE, y, f, cnt,failed);

  y := CylinderD(20.5, -20);   {N[ParabolicCylinderD[20.5, -20],30]}
  f := -0.4822736223072907922e35;
  testrel(5, NE, y, f, cnt,failed);

  y := CylinderD(4.5, 50);
  f := 1.621601576208132601e-264; {alpha}
  testrel(6, NE, y, f, cnt,failed);

  y := CylinderD(4.5, -30);  {Double: Overflow for x=-50}
  f := -0.1665643112041891717e92;
  testrel(7, NE, y, f, cnt,failed);

  y := CylinderD(-5.5, 10);
  f := 0.3714173811693538314e-16;
  testrel(8, NE, y, f, cnt,failed);

  y := CylinderD(-5.5, -10);
  f := 0.1177099087249377251e15;
  testrel(9, NE, y, f, cnt,failed);

  y := CylinderD(-5.5, 0.125);
  f := 0.7808084719547119203e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := CylinderD(-5.5, -0.125);
  f := 0.1367562241679692377;
  testrel(11, NE, y, f, cnt,failed);

  y := CylinderD(5.5, 1e-9);
  f := -6.841576188411541418;
  testrel(12, NE, y, f, cnt,failed);

  y := CylinderD(0, -10);
  f := 0.1388794386496402059e-10;
  testrel(13, NE, y, f, cnt,failed);

  y := CylinderD(-3, -8);
  f := 0.7239107161715334941e9;
  testrel(14, NE, y, f, cnt,failed);

  y := CylinderD(-3, 8);
  f := 0.2013035825596948112e-9;
  testrel(15, NE, y, f, cnt,failed);

  y := CylinderD(3, 8);
  f := 0.5491716526299844788e-4;
  testrel(16, NE, y, f, cnt,failed);

  y := CylinderD(3, -8);
  f := -0.5491716526299844788e-4;
  testrel(17, NE, y, f, cnt,failed);

  y := CylinderD(4, 8);
  f := 0.4180681740820476104e-3;
  testrel(18, NE, y, f, cnt,failed);

  y := CylinderD(4, -8);
  f := 0.4180681740820476104e-3;
  testrel(19, NE, y, f, cnt,failed);

  y := CylinderD(-4.5, 0);
  f := 0.2316724217633371968;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_hermiteh;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','HermiteH');

  y := HermiteH(-20.5,0);
  f := 0.5911985466323585402e-12;
  testrel(1, NE, y, f, cnt,failed);

  y := HermiteH(20.5,0);
  f := 0.1196056359093635307e13;
  testrel(2, NE, y, f, cnt,failed);

  y := HermiteH(0,-123);
  f := 1;
  testrel(3, NE, y, f, cnt,failed);

  y := HermiteH(3,1.5);
  f := 9;
  testrel(4, NE, y, f, cnt,failed);

  y := HermiteH(3,-1.5);
  f := -9;
  testrel(5, NE, y, f, cnt,failed);

  y := HermiteH(4,1.5);
  f := -15;
  testrel(6, NE, y, f, cnt,failed);

  y := HermiteH(4,-1.5);
  f := -15;
  testrel(7, NE, y, f, cnt,failed);

  y := HermiteH(0.5,2.0);
  f := 2.028395654496731894;
  testrel(8, NE, y, f, cnt,failed);

  y := HermiteH(5.5,55);
  f := 1.685663374942501320e11;
  testrel(9, NE, y, f, cnt,failed);

  y := HermiteH(0.5,-2.0);
  f := -13.597704371948080095;
  testrel(10, NE, y, f, cnt,failed);

  y := HermiteH(5.5,-25);
  f := 0.3686457341429633091e265;
  testrel(11, NE, y, f, cnt,failed);

  y := HermiteH(-2.5,10);
  f := 0.5471521794115101914e-3;  {Digits:=80;!!!}
  testrel(12, NE, y, f, cnt,failed);

  y := HermiteH(-2.5,-10);
  f := 0.1135536855006631709e46;
  testrel(13, NE, y, f, cnt,failed);

  y := HermiteH(5.5,0.125);
  f := -22.60805184613228380;
  testrel(14, NE, y, f, cnt,failed);

  y := HermiteH(5.5,1e-9);
  f := -46.02445508538673465;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cylinderu;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','CylinderU');

  y := CylinderU(3,10.5);
  f := 0.2666106211214159922e-15;
  testrel(1, NE, y, f, cnt,failed);

  y := CylinderU(3,20);
  f := 0.1019767537025606792e-47;
  testrel(2, NE, y, f, cnt,failed);

  y := CylinderU(100,0.5);
  f := 0.7807612330477519795e-81;
  testrel(3, NE, y, f, cnt,failed);

  y := CylinderU(5,5);
  f := 0.1552271294767621389e-6;
  testrel(4, NE, y, f, cnt,failed);

  y := CylinderU(5,-5);
  f := 45998.28922772747856;
  testrel(5, NE, y, f, cnt,failed);

  y := CylinderU(-5,5);
  f := 1.879976816273250088;
  testrel(6, NE, y, f, cnt,failed);

  y := CylinderU(5,-30);
  f := 0.1115127981951192897e104;
  testrel(7, NE, y, f, cnt,failed);

  y := CylinderU(13,17.5);
  f := 0.6851742325403592862e-50;
  testrel(8, NE, y, f, cnt,failed);

  y := CylinderU(13,-17.5);
  f := 0.1129967367045593180e41;
  testrel(9, NE, y, f, cnt,failed);

  y := CylinderU(20,20);
  f := 0.4701410003402728724e-70;
  testrel(10, NE, y, f, cnt,failed);

  y := CylinderU(-12,-8);
  f := 32436.55833588695948;
  testrel(11, NE, y, f, cnt,failed);

  y := CylinderU(4.5,50);
  f := 0.1170845970954047558e-279; {Digits:=150; !!}
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cylinderv;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','CylinderV');

  x := 1.5;
  y := CylinderV(0.0,x);
  f := 1.283927297107191901;
  testrel(1, NE, y, f, cnt,failed);

  y := CylinderV(0.5,x);
  f := 1.400331014153791334;
  testrel(2, NE, y, f, cnt,failed);

  y := CylinderV(1.0,x);
  f := 1.599899029475239644;
  testrel(3, NE, y, f, cnt,failed);

  y := CylinderV(1.5,x);
  f := 2.100496521230687000;
  testrel(4, NE, y, f, cnt,failed);

  y := CylinderV(2.0,x);
  f := 3.041812192766455417;
  testrel(5, NE, y, f, cnt,failed);

  y := CylinderV(2.5,x);
  f := 4.551075795999821834;
  testrel(6, NE, y, f, cnt,failed);

  y := CylinderV(3.0,x);
  f := 6.962566833362542592;
  testrel(7, NE, y, f, cnt,failed);

  y := CylinderV(3.5,x);
  f := 11.02760673646110675;
  testrel(8, NE, y, f, cnt,failed);

  y := CylinderV(4.0,x);
  f := 18.04838073195995243;
  testrel(9, NE, y, f, cnt,failed);

  y := CylinderV(4.5,x);
  f := 30.19463749269112563;
  testrel(10, NE, y, f, cnt,failed);

  y := CylinderV(20.0,x);
  f := 0.1040883971315052268e12;
  testrel(11, NE, y, f, cnt,failed);

  y := CylinderV(20.5,x);
  f := 0.2384630376131610783e12;
  testrel(12, NE, y, f, cnt,failed);

  y := CylinderV(200.0,x);
  f := 0.1165404358441727419e196;
  testrel(13, NE, y, f, cnt,failed);

  y := CylinderV(200.5,x);
  f := 0.4498864031728626495e196;
  testrel(14, NE, y, f, cnt,failed); {FPC}

  x := 12.5;
  y := CylinderV(0.0,x);
  f := 0.2085303487458286447e17;
  testrel(15, NE, y, f, cnt,failed);

  y := CylinderV(0.5,x);
  f := 0.73547558577480954217e17;
  testrel(16, NE, y, f, cnt,failed);

  y := CylinderV(1.0,x);
  f := 0.2598205934017800488e18;
  testrel(17, NE, y, f, cnt,failed);

  y := CylinderV(1.5,x);
  f := 0.9193444822185119277e18;
  testrel(18, NE, y, f, cnt,failed);

  y := CylinderV(2.0,x);
  f := 0.3258183934959542042e19;
  testrel(19, NE, y, f, cnt,failed);

  y := CylinderV(2.5,x);
  f := 0.1156535358630888005e20;
  testrel(20, NE, y, f, cnt,failed);

  y := CylinderV(3.0,x);
  f := 0.4111703007709694560e20;
  testrel(21, NE, y, f, cnt,failed);

  y := CylinderV(3.5,x);
  f := 0.1464056087932980245e21;
  testrel(22, NE, y, f, cnt,failed);

  y := CylinderV(4.0,x);
  f := 0.5221083358011106751e21;
  testrel(23, NE, y, f, cnt,failed);

  y := CylinderV(4.5,x);
  f := 0.1864766170675151946e22;
  testrel(24, NE, y, f, cnt,failed);

  y := CylinderV(20.0,x);
  f := 0.5113882356302268043e39;
  testrel(25, NE, y, f, cnt,failed);

  y := CylinderV(20.5,x);
  f := 0.1906065374402101901e40;
  testrel(26, NE, y, f, cnt,failed);

  y := CylinderV(200.0,x);
  f := 0.1081735411455045016e266;
  testrel(27, NE, y, f, cnt,failed);

  y := CylinderV(200.5,x);
  f := 0.5038731252253890156e266;
  testrel(28, NE, y, f, cnt,failed);

  y := CylinderV(1.5,1.5);
  f := 2.100496521230687000;
  testrel(29, NE, y, f, cnt,failed);

  y := CylinderV(0.5,0);
  f := 0.7978845608028653559;
  testrel(30, NE, y, f, cnt,failed);

  y := CylinderV(0.5,1e-5);
  f := 0.7978845608228124699;
  testrel(31, NE, y, f, cnt,failed);

  y := CylinderV(0.5,1.5);
  f := 1.400331014153791334;
  testrel(32, NE, y, f, cnt,failed);

  y := CylinderV(0.5,5);
  f := 413.3144351007517741;
  testrel(33, NE, y, f, cnt,failed);

  y := CylinderV(0.5,40);
  f := 0.4166130050162937620e174;
  testrel(34, NE, y, f, cnt,failed);

  y := CylinderV(1,0);
  f := 0.32800194866687646639;
  testrel(35, NE, y, f, cnt,failed);

  y := CylinderV(1,1e-10);
  f := 0.32800194870118709778;
  testrel(36, NE, y, f, cnt,failed);

  y := CylinderV(1,1.5);
  f := 1.599899029475239644;
  testrel(37, NE, y, f, cnt,failed);

  y := CylinderV(1,5);
  f := 919.3820780817696689;
  testrel(38, NE, y, f, cnt,failed);

  y := CylinderV(1,40);
  f := 0.2634686025645765032e175;
  testrel(39, NE, y, f, cnt,failed);

  y := CylinderV(0,1.5);
  f := 1.283927297107191901;
  testrel(40, NE, y, f, cnt,failed);

  y := CylinderV(0,1e-5);
  f := 0.68621590757881282595;
  testrel(41, NE, y, f, cnt,failed);

  y := CylinderV(0,1e-9);
  f := 0.68621262788732810586;
  testrel(42, NE, y, f, cnt,failed);

  y := CylinderV(0,0);
  f := 0.6862126275593261572;
  testrel(43, NE, y, f, cnt,failed);

  y := CylinderV(0,5);
  f := 187.9108321272246048;
  testrel(44, NE, y, f, cnt,failed);

  y := CylinderV(0,40);
  f := 0.6588775991761827045e173;
  testrel(45, NE, y, f, cnt,failed);

  {----- negative x -----}
  y := CylinderV(5/2, -2.0);
  f := 10.844375514192275467;
  testrel(46, NE, y, f, cnt,failed);

  y := CylinderV(3/2, -2.0);
  f := -4.337750205676910186;
  testrel(47, NE, y, f, cnt,failed);

  y := CylinderV(1/2, -2.0);
  f := 2.168875102838455093;
  testrel(48, NE, y, f, cnt,failed);

  y := CylinderV(2, -2.0);
  f := 0.1441906692704577264e-1;
  testrel(49, NE, y, f, cnt,failed);

  y := CylinderV(2, -5.0);
  f := 0.1253934231162904331e-4;
  testrel(50, NE, y, f, cnt,failed);

  y := CylinderV(2, -10.0);
  f := 1.781692680493916561e-14; {alpha}
  testrel(51, NE, y, f, cnt,failed);

  y := CylinderV(2, -20.0);
  f := 8.7048180138934015569e-48; {alpha}
  testrel(52, NE, y, f, cnt,failed);

  y := CylinderV(2, -50.0);
  f := 8.795313494255047943e-277; {alpha}
  testrel(53, NE, y, f, cnt,failed);

  y := CylinderV(19/2, -50.0);
  f := -4.2949336357374568326e286; {alpha}
  testrel(54, NE, y, f, cnt,failed);

  y := CylinderV(10, -50.0);
  f := 1.877207232604526005e-284; {alpha}
  testrel(55, NE, y, f, cnt,failed);

  y := CylinderV(21/2, -50.0);
  f := 2.155173106979384993e288;  {alpha}
  testrel(56, NE, y, f, cnt,failed);

  y := CylinderV(150, -2.0);
  f := 1.554325336318121333e119;  {alpha}
  testrel(57, NE, y, f, cnt,failed);

  y := CylinderV(200, -2.0);
  f := 3.590716405717152808e173;  {alpha}
  testrel(58, NE, y, f, cnt,failed);

  {a < 0}
  y := CylinderV(-2, -2.0);
  f := 0.7120211298549305625;
  testrel(59, NE, y, f, cnt,failed);

  y := CylinderV(-2, 0.125);
  f := -0.5319222596308388933;
  testrel(60, NE, y, f, cnt,failed);

  y := CylinderV(-2, 5.0);
  f := 9.206662277539298119;
  testrel(61, NE, y, f, cnt,failed);

  y := CylinderV(-2, 10.0);
  f := 0.1901626810882079821e9;
  testrel(62, NE, y, f, cnt,failed);

  y := CylinderV(-2.5, -2.0);
  f := 0.8679380128080086810e-1;
  testrel(63, NE, y, f, cnt,failed);

  y := CylinderV(-2.5, -0.0625);
  f := 0.4978666915953578549e-1;
  testrel(64, NE, y, f, cnt,failed);

  y := CylinderV(-2.5, 4.0);
  f := 1.218554345187601585;
  testrel(65, NE, y, f, cnt,failed);

  y := CylinderV(-2.5, -10.0);
  f := -0.6118449404325429709e8;
  testrel(66, NE, y, f, cnt,failed);

  y := CylinderV(-20, 0.25);
  f := 0.5430352848664535255e-9;
  testrel(67, NE, y, f, cnt,failed);

  y := CylinderV(-20, -0.5);
  f := -0.5706716204367091304e-9;
  testrel(68, NE, y, f, cnt,failed);

  y := CylinderV(-20, 4);
  f := -0.4303060953877838270e-9;
  testrel(69, NE, y, f, cnt,failed);

  y := CylinderV(-20, -5);
  f := -0.6292981467961312568e-9;
  testrel(70, NE, y, f, cnt,failed);

  {V(-12.5,x) is odd}
  y := CylinderV(-12.5, 1.5);
  f := -0.1867855278694266291e-4;
  testrel(71, NE, y, f, cnt,failed);

  y := CylinderV(-12.5, -2.5);
  f := -0.1560948544586090133e-4;
  testrel(72, NE, y, f, cnt,failed);

  {V(-9.5,x) is even}
  y := CylinderV(-9.5, 12);
  f := 84352.28865472105374;
  testrel(73, NE, y, f, cnt,failed);

  y := CylinderV(-9.5, -12.5);
  f := 1156953.0123654599948;
  testrel(74, NE, y, f, cnt,failed);

  y := CylinderV(-150, -2.0);
  f := -0.3717710476062361776e-131;
  testrel(75, NE, y, f, cnt,failed);

  y := CylinderV(-150.5, -4.0);
  f := -0.1060667074455632123e-131;
  testrel(76, NE, y, f, cnt,failed);

  {DAMath: calls pcfv_gen because |a| > MAXGAMD}
  y := CylinderV(-200, -2.0);
  f := -0.2283072372613928369e-187;
  testrel(77, NE, y, f, cnt,failed);

  y := CylinderV(-200.5, 2.0);
  f := -0.1841171876846235625e-189;
  testrel(78, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.


