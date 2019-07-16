{Part 6b of regression test for SPECFUNX unit  (c) 2018  W.Ehrhardt}

unit t_sfx6b;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface


procedure test_CoulombCLx;
procedure test_CoulombSLx;
procedure test_CoulombFFpx;
procedure test_CoulombGGpx;

procedure test_synchfx;
procedure test_synchgx;

implementation


uses sfbessl2,
  specfunx, t_sfx0;


{---------------------------------------------------------------------------}
procedure test_CoulombCLx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 5;
  NE1 = 12;
  NE2 = 100;
  NE3 = 500;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','CoulombCLx');

  {cl := (L,y) ->2^L*exp(-Pi/2*y)*abs(GAMMA(L+1+I*y))/(2*L+1)!;}

  y := CoulombCLx(0, 0);
  f := 1;
  testrel(1, NE, y, f, cnt,failed);

  y := CoulombCLx(1, 0);
  f := 1/3;
  testrel(2, NE, y, f, cnt,failed);

  y := CoulombCLx(100, 0);
  f := 0.7463087331140648667e-189;
  testrel(3, NE, y, f, cnt,failed);

  y := CoulombCLx(0, 3);
  f := 0.3503656340775789373e-3;
  testrel(4, NE, y, f, cnt,failed);

  y := CoulombCLx(0, -3);
  f := 4.3416075414867748016;
  testrel(5, NE, y, f, cnt,failed);

  y := CoulombCLx(0, 20);
  f := 0.5781996909122189450e-26;
  testrel(6, NE, y, f, cnt,failed);

  y := CoulombCLx(0, -1000);
  f := 79.26654595212022027;
  testrel(7, NE, y, f, cnt,failed);

  y := CoulombCLx(1, 10);
  f := 0.6030673659418024024e-12;
  testrel(8, NE, y, f, cnt,failed);

  y := CoulombCLx(1, -10);
  f := 26.553964257822392971;
  testrel(9, NE, y, f, cnt,failed);

  y := CoulombCLx(2, 10);
  f := 0.6150104533896502516e-12;
  testrel(10, NE, y, f, cnt,failed);

  y := CoulombCLx(2, -10);
  f := 27.07983638277634637;
  testrel(11, NE, y, f, cnt,failed);

  y := CoulombCLx(2, 0.125);
  f := 0.5461296580176828492e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := CoulombCLx(2, -0.125);
  f := 0.8088030980681890415e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := CoulombCLx(10, 1);
  f := 0.14417727229038187738e-10;
  testrel(14, NE, y, f, cnt,failed);

  y := CoulombCLx(10, 5);
  f := 0.8954369139191368151e-14;
  testrel(15, NE, y, f, cnt,failed);

  y := CoulombCLx(10, -5);
  f := 0.5941782675897785582e-7;
  testrel(16, NE, y, f, cnt,failed);

  y := CoulombCLx(100, 10);
  f := 0.6844263887409802494e-196;
  testrel(17, NE, y, f, cnt,failed);

  y := CoulombCLx(100, -10);
  f := 0.3013632454701968721e-182;
  testrel(18, NE, y, f, cnt,failed);

  y := CoulombCLx(100, 50);
  f := 0.3674392475607524167e-228;
  testrel(19, NE2, y, f, cnt,failed);

  y := CoulombCLx(100, 100);
  f := 0.4691031109905246504e-276;
  testrel(20, NE2, y, f, cnt,failed);

  y := CoulombCLx(100, -100);
  f := 0.1285001685407405358e-139;
  testrel(21, NE, y, f, cnt,failed);

  y := CoulombCLx(120,-20);
  f := 0.7652934282174957408e-223;
  testrel(22, NE, y, f, cnt,failed);

  y := CoulombCLx(120,20);
  f := 0.39473070212666459836e-250;
  testrel(23, NE1, y, f, cnt,failed);

  y := CoulombCLx(150,-50);
  f := 0.3257876285627941997e-278;
  testrel(24, NE, y, f, cnt,failed);

  {extended}
  y := CoulombCLx(150, 3);
  f := 0.7717453984666966428e-311;
  testrel(25, NE1, y, f, cnt,failed);

  y := CoulombCLx(400, 1e-100);
  f := 0.26772916875429901774e-990;
  testrel(26, NE1, y, f, cnt,failed);

  y := CoulombCLx(900, 1000);
  f := 0.6597023070360907704e-3431;
  testrel(27, NE3, y, f, cnt,failed);

  y := CoulombCLx(1000, 1);
  f := 0.1350010147148391328e-2870;
  testrel(28, NE1, y, f, cnt,failed);

  y := CoulombCLx(100, 1000);
  f := 0.3154286828558268651e-1409;
  testrel(29, NE3, y, f, cnt,failed);

  y := CoulombCLx(100, 200);
  f := 0.2468997056454290060e-386;
  testrel(30, NE2, y, f, cnt,failed);

  y := CoulombCLx(1500, 1);
  f := 0.2816689354341098661e-4568;
  testrel(31, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_CoulombSLx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 4;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','CoulombSLx');

  y := CoulombSLx(0,0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := CoulombSLx(0,1e-7);
  f := -0.5772156649015288538e-7;
  testrel(2, NE, y, f, cnt,failed);

  y := CoulombSLx(0,1e-5);
  f := -0.57721566486146429726e-5;
  testrel(3, NE, y, f, cnt,failed);

  y := CoulombSLx(0,0.25);
  f := -0.1382373401412416032;
  testrel(4, NE, y, f, cnt,failed);

  y := CoulombSLx(0,1);
  f := -0.3016403204675331979;
  testrel(5, NE2, y, f, cnt,failed);

  y := CoulombSLx(0,1.8055419921875);
  f := -0.3135957509996107496e-5;
  testrel(6, NE1, y, f, cnt,failed);

  y := CoulombSLx(0,2);
  f := 0.1296463163097883114;
  testrel(7, NE1, y, f, cnt,failed);

  y := CoulombSLx(0,-3);
  f := -1.053350771068613200;
  testrel(8, NE, y, f, cnt,failed);

  y := CoulombSLx(0,4);
  f := 2.309698056572534960;
  testrel(9, NE, y, f, cnt,failed);

  y := CoulombSLx(0,5);
  f := 3.8158985746149244778;
  testrel(10, NE, y, f, cnt,failed);

  y := CoulombSLx(0,10);
  f := 13.80291297422990069;
  testrel(11, NE1, y, f, cnt,failed);

  y := CoulombSLx(0,100);
  f := 361.3015834260953946;
  testrel(12, NE2, y, f, cnt,failed);

  y := CoulombSLx(0,1e10);
  f := 0.2202585093001899666e12;
  testrel(13, NE, y, f, cnt,failed);

  y := CoulombSLx(0,1e18);
  f := 0.4044653167389282231e20;
  testrel(14, NE, y, f, cnt,failed);

  y := CoulombSLx(0,-1e20);
  f := -0.4505170185988091368e22;
  testrel(15, NE, y, f, cnt,failed);

  y := CoulombSLx(0,-1e20);
  f := -0.4505170185988091368e22;
  testrel(15, NE, y, f, cnt,failed);

  y := CoulombSLx(0,1e300);
  f := 0.6897755278982137052e303;
  testrel(16, NE, y, f, cnt,failed);

  y := CoulombSLx(1,0);
  f := 0;
  testrel(17, NE, y, f, cnt,failed);

  y := CoulombSLx(1,1e-7);
  f := 0.4227843350984678129e-7;
  testrel(18, NE, y, f, cnt,failed);

  y := CoulombSLx(1,0.25);
  f := 0.1067413229856225509;
  testrel(19, NE2, y, f, cnt,failed);

  y := CoulombSLx(1,-1);
  f := -0.4837578429299151117;
  testrel(20, NE, y, f, cnt,failed);

  y := CoulombSLx(1,2);
  f := 1.236795034103878814;
  testrel(21, NE1, y, f, cnt,failed);

  y := CoulombSLx(10,0);
  f := 0;
  testrel(22, NE, y, f, cnt,failed);

  y := CoulombSLx(10,1e-7);
  f := 0.2351752589066721123e-6;
  testrel(23, NE, y, f, cnt,failed);

  y := CoulombSLx(10,-0.25);
  f := -0.5879617105592471912;
  testrel(24, NE1, y, f, cnt,failed);

  y := CoulombSLx(10,1);
  f := 2.353256829534718204;
  testrel(25, NE, y, f, cnt,failed);

  y := CoulombSLx(10,2);
  f := 4.715443169414506097;
  testrel(26, NE, y, f, cnt,failed);

  y := CoulombSLx(100,-0.5);
  f := -2.305082988957891971;
  testrel(27, NE, y, f, cnt,failed);

  y := CoulombSLx(1000,10);
  f := 69.08271845158168286;
  testrel(28, NE, y, f, cnt,failed);

  y := CoulombSLx(20000,100);
  f := 990.3516718754886751;
  testrel(29, NE1, y, f, cnt,failed);

  y := CoulombSLx(30000,10000);
  f := 103269.0066910751866;
  testrel(30, NE1, y, f, cnt,failed);


  y := CoulombSLx(0,1e4500);
  f := 0.1036063291847320558e4505;
  testrel(31, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_CoulombFFpx;
var
  f,fc,fcp: extended;
  cnt, failed, ierr: integer;

const
  NE  = 4;
  NE1 = 16;
  NE2 = 50;
  NE3 = 2000;

begin

  cnt := 0;
  failed := 0;
  writeln('Function: ','CoulombFFpx');

  {Some test cases from Barnett's COULTEST.REF, recomputed with Maple}
  {Barnett JWKB examples}

  CoulombFFpx(0, 100, 1, fc,fcp, ierr);
  f := 0.8999367996930662705e-125;
  testrel(1, NE2, fc, f, cnt,failed);
  f := 0.1292746745757069321e-123;
  testrel(2, NE2, fcp, f, cnt,failed);

  CoulombFFpx(30, 100, 1, fc,fcp, ierr);
  f := 0.7091353218216275425e-148;
  testrel(3, NE2, fc, f, cnt,failed);
  f := 0.2415384435355665713e-146;
  testrel(4, NE2, fcp, f, cnt,failed);

  CoulombFFpx(20, 100, 50, fc,fcp, ierr);
  f := 0.4535368874756233436e-55;
  testrel(5, NE2, fc, f, cnt,failed);
  f := 0.8103549101586879086e-55;
  testrel(6, NE2, fcp, f, cnt,failed);

  CoulombFFpx(10, -50, 5, fc,fcp, ierr);
  f := -0.3681143602184919517;
  testrel(7, NE1, fc, f, cnt,failed);
  f := 1.338467510317216956;
  testrel(8, NE1, fcp, f, cnt,failed);

  CoulombFFpx(20, 10, 5, fc,fcp, ierr);
  f := 0.4042838669803226134e-17;
  testrel(9, NE, fc, f, cnt,failed);
  f := 0.1837787576572863352e-16;
  testrel(10, NE, fcp, f, cnt,failed);

  CoulombFFpx(30, 10, 50, fc,fcp, ierr);
  f := -0.3526354364298439625;
  testrel(11, NE2, fc, f, cnt,failed);
  f := -0.6644367519319223979;
  testrel(12, NE, fcp, f, cnt,failed);

  {Next fails in COULTEST.REF}
  CoulombFFpx(0,-500, 0.001, fc,fcp, ierr);
  f := 0.3232536822545551523e-1;
  testrel(13, NE, fc, f, cnt,failed);
  f := 12.54904113796246721;
  testrel(14, NE, fcp, f, cnt,failed);


  {---- start GSL values for L, eta, x ----}
  CoulombFFpx(0, 1, 5, fc,fcp, ierr);
  f := 0.6849374120059439677;
  testrel(15, NE, fc, f, cnt,failed);
  f := -0.7236423862556063963;
  testrel(16, NE, fcp, f, cnt,failed);

  CoulombFFpx(10, 1, 5, fc,fcp, ierr);
  f := 0.6423773354915823698e-3;
  testrel(17, NE, fc, f, cnt,failed);
  f := 0.13299570958719702545e-2;
  testrel(18, NE, fcp, f, cnt,failed);

  CoulombFFpx(4, 50, 120, fc,fcp, ierr);
  f := 0.7351947118237984950e-1;
  testrel(19, NE3, fc, f, cnt,failed);
  f := 0.6368149124126783325;
  testrel(20, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, -1000, 1, fc,fcp, ierr);
  f := 0.9682225189913408376e-1;
  testrel(21, NE2, fc, f, cnt,failed);
  f := 5.120633962746318391;
  testrel(22, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, -50, 5, fc,fcp, ierr);
  f := 0.1522369757142367852;
  testrel(23, NE2, fc, f, cnt,failed);
  f := 2.030910411661365466;
  testrel(24, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, -50, 1000, fc,fcp, ierr);
  f := -0.2267212182760888524;
  testrel(25, NE3, fc, f, cnt,failed);
  f := -0.9961306810018401525;
  testrel(26, NE3, fcp, f, cnt,failed);

  CoulombFFpx(10, -50, 5, fc,fcp, ierr);
  f := -0.3681143602184919517;
  testrel(27, NE1, fc, f, cnt,failed);
  f := 1.338467510317216956;
  testrel(28, NE1, fcp, f, cnt,failed);

  CoulombFFpx(0, -4, 5, fc,fcp, ierr);
  f := 0.4078627230056171895;
  testrel(29, NE, fc, f, cnt,failed);
  f := 1.098212336357309464;
  testrel(30, NE, fcp, f, cnt,failed);

  CoulombFFpx(3, -4, 5, fc,fcp, ierr);
  f := -0.2568630935581316177;
  testrel(31, NE, fc, f, cnt,failed);
  f := 1.143229422014827760;
  testrel(32, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, 1, 2, fc,fcp, ierr);
  f := 0.6617816138326812982;
  testrel(33, NE, fc, f, cnt,failed);
  f := 0.4815574557099492028;
  testrel(34, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, 1, 0.5, fc,fcp, ierr);
  f := 0.8315404535022023302e-1;
  testrel(35, NE, fc, f, cnt,failed);
  f := 0.2269387461622278757;
  testrel(36, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, 8, 1.05, fc,fcp, ierr);
  f := 0.9882706082810274357e-8;
  testrel(37, NE1, fc, f, cnt,failed);
  f := 0.4005167028235547770e-7;
  testrel(38, NE1, fcp, f, cnt,failed);

  CoulombFFpx(0, 50, 0.1, fc,fcp, ierr);
  f := 0.2807788027954216071e-66;
  testrel(39, NE2, fc, f, cnt,failed);
  f := 0.9677600748751576606e-65;
  testrel(40, NE2, fcp, f, cnt,failed);

  CoulombFFpx(0, 10, 5, fc,fcp, ierr);
  f := 0.1720745409178793061e-5;
  testrel(41, NE, fc, f, cnt,failed);
  f := 0.3097599470640545805e-5;
  testrel(42, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, 25, 10, fc,fcp, ierr);
  f := 0.1545127450107611431e-15;
  testrel(43, NE2, fc, f, cnt,failed);
  f := 0.3139086939337863093e-15;
  testrel(44, NE2, fcp, f, cnt,failed);

  CoulombFFpx(0, 1, 9.2, fc,fcp, ierr);
  f := -0.2563201231975795565;
  testrel(45, NE2, fc, f, cnt,failed);
  f := 0.9151879228672422037;
  testrel(46, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, 10, 10, fc,fcp, ierr);
  f := 0.1626271125013587825e-2;
  testrel(47, NE, fc, f, cnt,failed);
  f := 0.1706047632079280601e-2;
  testrel(48, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, 100, 1, fc,fcp, ierr);
  f := 0.8999367996930662705e-125;
  testrel(49, NE2, fc, f, cnt,failed);
  f := 0.1292746745757069321e-123;
  testrel(50, NE2, fcp, f, cnt,failed);

  CoulombFFpx(1, 0, 3.25, fc,fcp, ierr);
  f := 0.9608388654558974880;
  testrel(51, NE, fc, f, cnt,failed);
  f := -0.403837862362692219;
  testrel(52, NE1, fcp, f, cnt,failed);

  CoulombFFpx(37, 0, 1.2693881947287221e-07, fc,fcp, ierr);
  f := 0.6589072427862341297e-317;
  testrel(53, NE1, fc, f, cnt,failed);
  f := 0.1972483699616239865e-308;
  testrel(54, NE1, fcp, f, cnt,failed);
  {---- end GSL ----}

  {CF2 does not converge}
  CoulombFFpx(0, -2.5, 1e-4, fc,fcp, ierr);
  f := 0.3962336840308557983e-3;
  testrel(55, NE, fc, f, cnt,failed);
  f := 3.961346160330035983;
  testrel(56, NE, fcp, f, cnt,failed);

  CoulombFFpx(0, -1, 1e-40, fc,fcp, ierr);
  f := 0.2508972050168545737e-39;
  testrel(57, NE, fc, f, cnt,failed);
  f := 2.508972050168545737;
  testrel(58, NE, fcp, f, cnt,failed);

  CoulombFFpx(5, 1, 1e-40, fc,fcp, ierr);
  f := 0.1827363156789148311e-244;
  testrel(59, NE, fc, f, cnt,failed);
  f := 0.1096417894073488987e-203;
  testrel(60, NE1, fcp, f, cnt,failed);

  CoulombFFpx(50, -10, 100, fc,fcp, ierr);
  f := -0.8047984148934464893;
  testrel(61, NE1, fc, f, cnt,failed);
  f := -0.5993651172838818544;
  testrel(62, NE2, fcp, f, cnt,failed);

  CoulombFFpx(4, 1, 100, fc,fcp, ierr);
  f := 0.88114358268583090838;
  testrel(63, NE1, fc, f, cnt,failed);
  f := -0.4792257866629417248;
  testrel(64, NE2, fcp, f, cnt,failed);

  CoulombFFpx(4, 10, 10, fc,fcp, ierr);
  f := 0.5812004954531329642e-3;
  testrel(65, NE, fc, f, cnt,failed);
  f := 0.6654184996305726707e-3;
  testrel(66, NE, fcp, f, cnt,failed);

  CoulombFFpx(4, 20, 10, fc,fcp, ierr);
  f := 0.3353643854604955045e-11;
  testrel(67, NE, fc, f, cnt,failed);
  f := 0.6115724294011260319e-11;
  testrel(68, NE, fcp, f, cnt,failed);

  CoulombFFpx(4, 20, 20, fc,fcp, ierr);
  f := 0.3274428895706980578e-5;
  testrel(69, NE, fc, f, cnt,failed);
  f := 0.3436373424890256195e-5;
  testrel(70, NE, fcp, f, cnt,failed);

  {Next uses JWKB for double}
  CoulombFFpx(4, 100, 100, fc,fcp, ierr);
  f := 0.7329253469930406898e-25;
  testrel(71, NE2, fc, f, cnt,failed);
  f := 0.7373136842030188765e-25;
  testrel(72, NE2, fcp, f, cnt,failed);

  CoulombFFpx(4, 0.05, 0.5, fc,fcp, ierr);
  f := 0.3036846891106875154e-4;
  testrel(73, NE, fc, f, cnt,failed);
  f := 0.3026065847565835231e-3;
  testrel(74, NE, fcp, f, cnt,failed);

  CoulombFFpx(2,10, 1e-10, fc,fcp, ierr);
  f := 0.6150104535946537361e-42;
  testrel(75, NE, fc, f, cnt,failed);
  f := 0.1845031360988964693e-31;
  testrel(76, NE, fcp, f, cnt,failed);

  CoulombFFpx(0,10, 1e-8, fc,fcp, ierr);
  f := 0.1800223551964554248e-20;
  testrel(77, NE, fc, f, cnt,failed);
  f := 0.1800223731986903383e-12;
  testrel(78, NE, fcp, f, cnt,failed);

  CoulombFFpx(0,10, 1, fc,fcp, ierr);
  f := 0.3696585553278718968e-10;
  testrel(79, NE, fc, f, cnt,failed);
  f := 0.1717251970258742264e-9;
  testrel(80, NE, fcp, f, cnt,failed);

  CoulombFFpx(10,1,20, fc,fcp, ierr);
  f := -0.4457416646219343703;
  testrel(81, NE, fc, f, cnt,failed);
  f :=  0.8226484280767422858;
  testrel(82, NE, fcp, f, cnt,failed);

  CoulombFFpx(5,0,20, fc,fcp, ierr);
  f := 0.3336781612619138553;
  testrel(83, NE, fc, f, cnt,failed);
  f := 0.9261034438714763083;
  testrel(84, NE, fcp, f, cnt,failed);

  CoulombFFpx(2,1,20, fc,fcp, ierr);
  f := 1.018006874993061383;
  testrel(85, NE1, fc, f, cnt,failed);
  f := -0.1550724192405322779;
  testrel(86, NE2, fcp, f, cnt,failed);

  CoulombFFpx(0,1,20, fc,fcp, ierr);
  f := -0.3292255362657539884;
  testrel(87, NE1, fc, f, cnt,failed);
  f := -0.9221468899352132127;
  testrel(88, NE1, fcp, f, cnt,failed);

  CoulombFFpx(0,1,2.5, fc,fcp, ierr);
  f := 0.8953375775746753376;
  testrel(89, NE, fc, f, cnt,failed);
  f := 0.4376968273858673597;
  testrel(90, NE, fcp, f, cnt,failed);

  CoulombFFpx(1,1,2.5, fc,fcp, ierr);
  f := 0.5768405868274806900;
  testrel(91, NE, fc, f, cnt,failed);
  f := 0.4586217235499059814;
  testrel(92, NE, fcp, f, cnt,failed);

  CoulombFFpx(2, 0.75, 1.75, fc,fcp, ierr);
  f := 0.1223277573186102112;
  testrel(93, NE, fc, f, cnt,failed);
  f := 0.2096540359862375426;
  testrel(94, NE, fcp, f, cnt,failed);

  {Single test for F function only}
  fc := CoulombFx(2, 0.75, 1.75);
  f := 0.1223277573186102112;
  testrel(95, NE, fc, f, cnt,failed);

  {extended only}
  CoulombFFpx(0, 500, 0.025, fc,fcp, ierr);
  f := 0.4853159323960411709e-679;
  testrel(96, NE3, fc, f, cnt,failed);
  f := 0.1023214917134743202e-676;
  testrel(97, NE3, fcp, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_CoulombGGpx;
var
  f,gc,gcp: extended;
  cnt, failed, ierr: integer;
const
  NE = 8;
  NE1 = 40;
  NE2 = 100;
  NE3 = 1800;
  ERel  = 0.008; {Direct relative error for JWKB cases}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','CoulombGGpx');

  CoulombGGpx(10, -50, 5, gc,gcp, ierr);
  f := 0.3315883246109352578;
  testrel(1, NE, gc, f, cnt,failed);
  f := 1.510888628136178449;
  testrel(2, NE, gcp, f, cnt,failed);

  CoulombGGpx(0, -10, 0.1, gc,gcp, ierr);
  f := -0.1533930267618161314;
  testrel(3, NE2, gc, f, cnt,failed);
  f := -3.392874637658605463;
  testrel(4, NE1, gcp, f, cnt,failed);

  CoulombGGpx(4, 5, 10, gc,gcp, ierr);
  f := 2.667598149114350864;
  testrel(5, NE, gc, f, cnt,failed);
  f := -0.8736758151710172844;
  testrel(5, NE, gcp, f, cnt,failed);

  CoulombGGpx(4, -5, 10, gc,gcp, ierr);
  f := 0.8214321599693482372;
  testrel(7, NE, gc, f, cnt,failed);
  f := 0.3636620155054198230;
  testrel(8, NE1, gcp, f, cnt,failed);

  CoulombGGpx(4, 0.05, 5, gc,gcp, ierr);
  f := 0.9719392861938234616;
  testrel(9, NE, gc, f, cnt,failed);
  f := -0.6576649675567845689;
  testrel(10, NE, gcp, f, cnt,failed);

  CoulombGGpx(0, 0.05, 1, gc,gcp, ierr);
  f := 0.6101942704496222569;
  testrel(11, NE, gc, f, cnt,failed);
  f := -0.7945505629724031144;
  testrel(12, NE, gcp, f, cnt,failed);

  CoulombGGpx(10,1,20, gc,gcp, ierr);
  f := 1.031338769952050517;
  testrel(13, NE, gc, f, cnt,failed);
  f := 0.3400417643544783185;
  testrel(14, NE, gcp, f, cnt,failed);

  CoulombGGpx(5,0,20, gc,gcp, ierr);
  f := 0.9634469551474556235;
  testrel(15, NE, gc, f, cnt,failed);
  f := -0.3229113240224271940;
  testrel(16, NE, gcp, f, cnt,failed);

  CoulombGGpx(2,1,20, gc,gcp, ierr);
  f := -0.1628453040171422116;
  testrel(17, NE1, gc, f, cnt,failed);
  f := -0.9575055028491291390;
  testrel(18, NE, gcp, f, cnt,failed);

  CoulombGGpx(0,1,20, gc,gcp, ierr);
  f := -0.9724283986971198997;
  testrel(19, NE, gc, f, cnt,failed);
  f := 0.3137003818968775765;
  testrel(20, NE1, gcp, f, cnt,failed);

  CoulombGGpx(0,1,2.5, gc,gcp, ierr);
  f := 0.9739143230656190246;
  testrel(21, NE, gc, f, cnt,failed);
  f := -0.6407871232241147580;
  testrel(22, NE, gcp, f, cnt,failed);

  CoulombGGpx(0, 0.5,2.5, gc,gcp, ierr);
  f := 0.1023855310825321012;
  testrel(23, NE1, gc, f, cnt,failed);
  f := -0.8946050863218798490;
  testrel(24, NE, gcp, f, cnt,failed);

  {GSL test cases for G, G': Note reduced accuracy, recomputed with Maple/MMA}
  CoulombGGpx(0, 1, 5, gc,gcp, ierr);
  f := -0.8984143590920205487;
  testrel(25, NE, gc, f, cnt,failed);
  f := -0.5108047585190350106;
  testrel(26, NE, gcp, f, cnt,failed);

  CoulombGGpx(0, -1000, 1, gc,gcp, ierr);
  f := 0.1139367843794720103;
  testrel(27, NE, gc, f, cnt,failed);
  f := -4.302434865224374720;
  testrel(28, NE1, gcp, f, cnt,failed);

  CoulombGGpx(0, -50, 5, gc,gcp, ierr);
  f := 0.4416806902362504456;
  testrel(29, NE, gc, f, cnt,failed);
  f := -0.6764853747668713568;
  testrel(30, NE1, gcp, f, cnt,failed);

  CoulombGGpx(0, -50, 1000, gc,gcp, ierr);
  f := -0.9497684438900352186;
  testrel(31, NE2, gc, f, cnt,failed);
  f := 0.2377656295411961399;
  testrel(32, NE3, gcp, f, cnt,failed);

  CoulombGGpx(10, -50, 5, gc,gcp, ierr);
  f := 0.3315883246109352578;
  testrel(33, NE, gc, f, cnt,failed);
  f := 1.510888628136178449;
  testrel(34, NE, gcp, f, cnt,failed);

  CoulombGGpx(0, -4, 5, gc,gcp, ierr);
  f := 0.6743270353832443793;
  testrel(35, NE, gc, f, cnt,failed);
  f := -0.6361104272804454129;
  testrel(36, NE, gcp, f, cnt,failed);

  CoulombGGpx(3, -4, 5, gc,gcp, ierr);
  f := 0.7879899223927996790;
  testrel(37, NE, gc, f, cnt,failed);
  f := 0.3859905878106711480;
  testrel(38, NE, gcp, f, cnt,failed);

  CoulombGGpx(0, 1, 2, gc,gcp, ierr);
  f := 1.275778784768276589;
  testrel(39, NE, gc, f, cnt,failed);
  f := -0.5827288130971847434;
  testrel(40, NE, gcp, f, cnt,failed);

  CoulombGGpx(0, 1, 0.5, gc,gcp, ierr);
  f := 3.106006927954887514;
  testrel(41, NE2, gc, f, cnt,failed);
  f := -3.549156038719924236;
  testrel(42, NE2, gcp, f, cnt,failed);

  CoulombGGpx(0, 1, 9.2, gc,gcp, ierr);
  f := 1.031205859189734661;
  testrel(43, NE, gc, f, cnt,failed);
  f := 0.2194632671749125019;
  testrel(44, NE1, gcp, f, cnt,failed);

  {GSL-JWKB: Large errors}
  CoulombGGpx(0, 8, 1.05, gc,gcp, ierr);
  f := 13331279.92006686320;
  testrele(45, ERel , gc, f, cnt,failed);
  f := -47159145.30842402330;
  testrele(46, ERel , gcp, f, cnt,failed);

  CoulombGGpx(0, 50, 0.1, gc,gcp, ierr);
  f := 0.5579810686998358766e65;
  testrele(47, ERel , gc, f, cnt,failed);
  f := -0.1638329512756321424e67;
  testrele(48, ERel , gcp, f, cnt,failed);

  CoulombGGpx(0, 10, 5, gc,gcp, ierr);
  f := 167637.5660945996762;
  testrele(49, ERel , gc, f, cnt,failed);
  f := -279370.7665536180307;
  testrele(50, ERel , gcp, f, cnt,failed);

  CoulombGGpx(0, 25, 10, gc,gcp, ierr);
  f := 1.617712900833631826e15;
  testrele(51, ERel , gc, f, cnt,failed);
  f := -3.185406201314974086e15;
  testrele(52, ERel , gcp, f, cnt,failed);

  CoulombGGpx(0, 10, 10, gc,gcp, ierr);
  f := 307.8732166109083799;
  testrele(53, ERel , gc, f, cnt,failed);
  f := -291.9277238082682287;
  testrele(54, ERel , gcp, f, cnt,failed);

  CoulombGGpx(0, 100, 1, gc,gcp, ierr);
  f := 0.3936654148133683610e124;
  testrele(55, ERel , gc, f, cnt,failed);
  f := -0.5456942268061526371e125;
  testrele(56, ERel , gcp, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_synchfx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','SynchFx');

  y := SynchFx(0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := SynchFx(1e-8);
  f := 0.4631000072771472311e-2;
  testrel(2, NE, y, f, cnt,failed);

  y := SynchFx(1e-5);
  f := 0.4629204411487711933e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := SynchFx(1/128);
  f := 0.4123549949580485690;
  testrel(4, NE, y, f, cnt,failed);

  y := SynchFx(1/64);
  f := 0.5090660116671645589;
  testrel(5, NE1, y, f, cnt,failed);

  y := SynchFx(0.025);
  f := 0.5832543586462289715;
  testrel(6, NE, y, f, cnt,failed);

  y := SynchFx(0.125);
  f := 0.8511257213236801121;
  testrel(7, NE, y, f, cnt,failed);

  y := SynchFx(0.25);
  f := 0.9157998878007424141;
  testrel(8, NE1, y, f, cnt,failed);

  y := SynchFx(0.5);
  f := 0.8708191468754688509;
  testrel(9, NE, y, f, cnt,failed);

  y := SynchFx(1.0);
  f := 0.6514228153553639698;
  testrel(10, NE, y, f, cnt,failed);

  y := SynchFx(2.0);
  f := 0.3016359028507394028;
  testrel(11, NE, y, f, cnt,failed);

  y := SynchFx(3.0);
  f := 0.1285657100090638130;
  testrel(12, NE, y, f, cnt,failed);

  y := SynchFx(4.0);
  f := 0.5282739669786681830e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := SynchFx(5.0);
  f := 0.21248129774981984268e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := SynchFx(10.0);
  f := 0.192238264300868974186e-3;
  testrel(15, NE, y, f, cnt,failed);

  y := SynchFx(25.0);
  f := 0.8956424677223712774e-10;
  testrel(16, NE, y, f, cnt,failed);

  y := SynchFx(50.0);
  f := 1.734785203976593787e-21;
  testrel(17, NE, y, f, cnt,failed);

  y := SynchFx(60.0);
  f := 0.8606939212110731423e-25;
  testrel(18, NE, y, f, cnt,failed);

  y := SynchFx(70.0);
  f := 0.4213333513123880200e-29;
  testrel(19, NE, y, f, cnt,failed);

  y := SynchFx(80.0);
  f := 0.2042253718758261554e-33;
  testrel(20, NE, y, f, cnt,failed);

  y := SynchFx(100.0);
  f := 0.4697593665922171859e-42;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_synchgx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','SynchGx');

  y := sfc_synchg(0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := sfc_synchg(1e-10);
  f := 0.4988613141720601489e-3;
  testrel(2, NE, y, f, cnt,failed);

  y := sfc_synchg(4e-10);
  f := 0.7918929749027351694e-3;
  testrel(3, NE1, y, f, cnt,failed);

  y := sfc_synchg(1e-9);
  f := 0.1074764120765973600e-2;
  testrel(4, NE1, y, f, cnt,failed);

  y := sfc_synchg(1.69e-8);
  f := 0.2758091869172511582e-2;
  testrel(5, NE, y, f, cnt,failed);

  y := sfc_synchg(1/1024);
  f := 0.1066180156185777283;
  testrel(6, NE1, y, f, cnt,failed);

  y := sfc_synchg(1/64);
  f := 0.2675041309969837707;
  testrel(7, NE, y, f, cnt,failed);

  y := sfc_synchg(0.125);
  f := 0.5040422411091107865;
  testrel(8, NE, y, f, cnt,failed);

  y := sfc_synchg(0.5);
  f := 0.6029652323601678511;
  testrel(9, NE, y, f, cnt,failed);

  y := sfc_synchg(1);
  f := 0.4944750621042082670;
  testrel(10, NE, y, f, cnt,failed);

  y := sfc_synchg(2);
  f := 0.2496778549762566211;
{$ifdef BIT16}
  testrel(11, 16, y, f, cnt,failed);  {????}
{$else}
  testrel(11, NE1, y, f, cnt,failed);
{$endif}

  y := sfc_synchg(3.75);
  f := 0.5841380816344249859e-1;
  testrel(12, NE1, y, f, cnt,failed);

  y := sfc_synchg(4);
  f := 0.4692320582610133071e-1;
  testrel(13, NE1, y, f, cnt,failed);

  y := sfc_synchg(5.0);
  f := 0.1922212317248410644e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := sfc_synchg(10.0);
  f := 0.1816118756953020428e-3;
  testrel(15, NE, y, f, cnt,failed);

  y := sfc_synchg(25);
  f := 0.87362343746221526073e-10;
  testrel(16, NE, y, f, cnt,failed);

  y := sfc_synchg(50);
  f := 0.1712604265071687323e-20;
  testrel(17, NE, y, f, cnt,failed);

  y := sfc_synchg(100);
  f := 0.4666936458728046656e-42;
  testrel(18, NE, y, f, cnt,failed);

  y := sfc_synchg(700);
  f := 0.3269880654608995997e-302;
  testrel(19, NE, y, f, cnt,failed);

  y := sfc_synchg(11000);
  f := 0.7576335246898575055e-4775;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.
