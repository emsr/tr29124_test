{Part 3a of regression test for SPECFUND unit  (c) 2010-2014  W.Ehrhardt}

unit t_sfd3a;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_cauchy_cdf;
procedure test_cauchy_inv;
procedure test_cauchy_pdf;

procedure test_chi2_cdf;
procedure test_chi2_inv;
procedure test_chi2_pdf;

procedure test_chi_dist;

procedure test_exp_cdf;
procedure test_exp_inv;
procedure test_exp_pdf;

procedure test_gamma_cdf;
procedure test_gamma_inv;
procedure test_gamma_pdf;

procedure test_invgamma_pdf;
procedure test_invgamma_cdf;
procedure test_invgamma_inv;


implementation


uses
  damath, specfund, t_sfd0;

{---------------------------------------------------------------------------}
procedure test_cauchy_cdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','cauchy_cdf');

  {statevalf[cdf,cauchy[a,b]](x);}
  a := 5;
  b := 2.0;
  y := cauchy_cdf(a, b, -70);
  f := 0.8486252456738483251e-2;
  testrel( 1, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, -1);
  f := 0.1024163823495667258;
  testrel( 2, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 0);
  f := 0.1211189415908433987;
  testrel( 3, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 2);
  f := 0.1871670418109988162;
  testrel( 4, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 5);
  f := 0.5;
  testrel( 5, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 10);
  f := 0.8788810584091566013;
  testrel( 6, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 100);
  f := 0.9932997290043335682;
  testrel( 7, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 1e8);
  f := 0.9999999936338019580;
  testrel( 8, NE, y, f, cnt,failed);


  a := -1.5;
  b := 0.5;
  y := cauchy_cdf(a, b, -100);
  f := 0.1615772346391954910e-2;
  testrel( 9, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, -3);
  f := 0.1024163823495667258;
  testrel(10, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, -1);
  f := 0.75;
  testrel(11, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 0);
  f := 0.8975836176504332742;
  testrel(12, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 1);
  f := 0.9371670418109988162;
  testrel(13, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 10);
  f := 0.9861691504333380231;
  testrel(14, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 1000);
  f := 0.9998410834449638694;
  testrel(15, NE, y, f, cnt,failed);

  y := cauchy_cdf(a, b, 1e8);
  f := 0.9999999984084505930;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cauchy_inv;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 4;
  NE1 = 8;
  NE2 = 16;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','cauchy_inv');

  a := 5;
  b := 2.0;

  y := cauchy_inv(a, b, ldexpd(1,-20));
  f := -667539.2144281116036;
  testrel( 1, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, ldexpd(1,-10));
  f := -646.8966015954026980;
  testrel( 2, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.125);
  f := 0.1715728752538099024;
  testrel( 3, NE2, y, f, cnt,failed);   {abs(err) ~ 2 eps}

  y := cauchy_inv(a, b, 0.5);
  f := 5.0;
  testrel( 4, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.75);
  f := 7.0;
  testrel( 5, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.875);
  f := 9.828427124746190098;
  testrel( 6, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.9);
  f := 11.15536707435050681;
  testrel( 7, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.99);          {!!!!! not exact}
  f := 68.64103190754791608;
  testrel( 8, NE1, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.9990234375);
  f := 656.8966015954026979;
  testrel( 9, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 1-ldexpd(1,-20));
  f := 667549.2144281116036;
  testrel(10, NE, y, f, cnt,failed);

  a := -1.5;
  b := 0.5;
  y := cauchy_inv(a, b, ldexpd(1,-20));
  f := -166887.55360702790089488646;
  testrel(11, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, ldexpd(1,-10));
  f := -164.47415039885067448789416;
  testrel(12, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.125);
  f := -2.7071067811865475244008445;
  testrel(13, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.5);
  f := -1.5;
  testrel(14, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.75);
  f := -1.0;
  testrel(15, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.875);
  f := -0.2928932188134524756;
  testrel(16, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.9375);
  f := 1.013669746062924052;
  testrel(17, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.99);             {!!!!! not exact}
  f := 14.41025797688697902;
  testrel(18, NE1, y, f, cnt,failed);

  y := cauchy_inv(a, b, 0.9990234375);
  f := 161.4741503988506745;
  testrel(19, NE, y, f, cnt,failed);

  y := cauchy_inv(a, b, 1-ldexpd(1,-20));
  f := 166884.5536070279009;
  testrel(20, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cauchy_pdf;
var
  y,f,a,b: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','cauchy_pdf');

  a := 5;
  b := 2.0;
  f := 0.1130964242969588458e-3;
  y := cauchy_pdf(a, b, -70);
  testrel( 1, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, -1);
  f := 0.1591549430918953358e-1;
  testrel( 2, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 0);
  f := 0.2195240594370970149e-1;
  testrel( 3, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 2);
  f := 0.4897075172058318024e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 5);
  f := 0.1591549430918953358;
  testrel( 5, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 10);
  f := 0.2195240594370970149e-1;
  testrel( 6, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 1000);
  f := 0.6430314388442978368e-6;
  testrel( 7, NE, y, f, cnt,failed);

  a := -1.5;
  b := 0.5;
  y := cauchy_pdf(a, b, -100);
  f := 0.1640349838617833917e-4;
  testrel( 8, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, -1);
  f := 0.3183098861837906715;
  testrel( 9, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 0);
  f := 0.6366197723675813431e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 2);
  f := 0.1273239544735162686e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 5);
  f := 0.3744822190397537312e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 100);
  f := 0.1544818666264453635e-4;
  testrel(13, NE, y, f, cnt,failed);

  y := cauchy_pdf(a, b, 10000);
  f := 0.1591072069521295231e-8;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_chi2_cdf;
var
  y,f: double;
  cnt, failed: integer;
  nu: longint;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chi2_cdf');

  {statevalf[cdf,chisquare[nu]](x);}

  nu := 1;

  y := chi2_cdf(nu, 1e-10);
  f := 0.7978845607895672798e-5;
  testrel( 1, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.001953125);
  f := 0.3525037386732282597e-1;
  testrel( 2, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.125);
  f := 0.27632639016823693296;
  testrel( 3, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.25);
  f := 0.38292492254802620727;
  testrel( 4, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.5);
  f := 0.52049987781304653768;
  testrel( 5, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.75);
  f := 0.61352376922876733507;
  testrel( 6, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 1.0);
  f := 0.68268949213708589717;
  testrel( 7, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 2.0);
  f := 0.84270079294971486934;
  testrel( 8, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 4.0);
  f := 0.95449973610364158560;
  testrel( 9, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 10.0);
  f := 0.99843459774199745032;
  testrel(10, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 50.0);
  f := 0.99999999999846254021;
  testrel(11, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 100.0);
  f := 1.0;
  testrel(12, NE, y, f, cnt,failed);

  nu := 2;

  y := chi2_cdf(nu, 0);
  f := 0.0;
  testrel(13, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.001953125);
  f := 0.97608581802433776529e-3;
  testrel(14, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.125);
  f := 0.6058693718652421388e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.25);
  f := 0.11750309741540459714;
  testrel(16, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.5);
  f := 0.22119921692859513175;
  testrel(17, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.75);
  f := 0.31271072120902780145;
  testrel(18, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 1.0);
  f := 0.39346934028736657640;
  testrel(19, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 2.0);
  f := 0.63212055882855767840;
  testrel(20, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 4.0);
  f := 0.86466471676338730811;
  testrel(21, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 10.0);
  f := 0.99326205300091453290;
  testrel(22, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 50.0);
  f := 0.99999999998611205614;
  testrel(23, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 100.0);
  f := 1.0;
  testrel(24, NE, y, f, cnt,failed);

  nu := 5;

  y := chi2_cdf(nu, 0);
  f := 0.0;
  testrel(25, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.001953125);
  f := 0.89612990307108155127e-8;
  testrel(26, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.125);
  f := 0.28104397654051540762e-3;
  testrel(27, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.25);
  f := 0.15208185533684396863e-2;
  testrel(28, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.5);
  f := 0.78767067673704077948e-2;
  testrel(29, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.75);
  f := 0.19887707187131059486e-1;
  testrel(30, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 1.0);
  f := 0.37434226752703631043e-1;
  testrel(31, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 2.0);
  f := 0.15085496391539036377;
  testrel(32, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 4.0);
  f := 0.45058404864721976739;
  testrel(33, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 10.0);
  f := 0.92476475385348782128;
  testrel(34, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 50.0);
  f := 0.99999999861420266330;
  testrel(35, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 100.0);
  f := 0.99999999999999999995;
  testrel(36, NE, y, f, cnt,failed);

  nu := 20;

  y := chi2_cdf(nu, 0);
  f := 0.0;
  testrel(37, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.001953125);
  f := 0.21719600981878520805e-36;
  testrel(38, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.125);
  f := 0.23679208135235717612e-18;
  testrel(39, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.25);
  f := 0.22909148211916296278e-15;
  testrel(40, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.5);
  f := 0.20942485399973611160e-12;
  testrel(41, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 0.75);
  f := 0.10782224479191469117e-10;
  testrel(42, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 1.0);
  f := 0.17096700293489033565e-9;
  testrel(43, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 2.0);
  f := 0.11142547833872067735e-6;
  testrel(44, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 4.0);
  f := 0.46498075017263808251e-4;
  testrel(45, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 10.0);
  f := 0.31828057306204811737e-1;
  testrel(46, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 50.0);
  f := 0.99977852336175121642;
  testrel(47, NE, y, f, cnt,failed);

  y := chi2_cdf(nu, 100.0);
  f := 0.99999999999874039154;
  testrel(48, NE, y, f, cnt,failed);

  y := chi2_cdf(100, 200);
  f := 0.99999998821549927902;
  testrel(49, NE, y, f, cnt,failed);

  y := chi2_cdf(1000, 1234.5);
  f := 0.99999951255538488480;
  testrel(50, NE, y, f, cnt,failed);

  y := chi2_cdf(100000, 101010);
  f := 0.98785020993089822507;
  testrel(51, NE, y, f, cnt,failed);

  y := chi2_cdf(1200000, 1201234);
  f := 0.78718764724057492688;
  testrel(52, NE, y, f, cnt,failed);

  y := chi2_cdf(2,MinDouble);
  f := 1.1125369292536006915e-308;
  testrel(53, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_chi2_inv;
var
  y,f: double;
  cnt, failed: integer;
  nu: longint;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chi2_inv');

  {statevalf[icdf,chisquare[nu]](p);}

  nu := 1;

  y := chi2_inv(nu, 1e-6);
  f := 0.15707963267957190874e-11;
  testrel( 1, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.001953125);
  f := 0.59921244211799117736e-5;
  testrel( 2, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.03125);
  f := 0.15347656753460162514e-2;
  testrel( 3, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.125);
  f := 0.24746651492520595311e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.25);
  f := 0.10153104426762154521;
  testrel( 5, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.5);
  f := 0.45493642311957275194;
  testrel( 6, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.75);
  f := 1.3233036969314659497;
  testrel( 7, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9);
  f := 2.7055434540954145671;
  testrel( 8, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.95);
  f := 3.8414588206941259584;
  testrel( 9, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.99);
  f := 6.6348966010212151384;
  testrel(10, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9990234375);
  f := 10.871483958875362877;
  testrel(11, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 1-ldexpd(1,-20)); {0.99999904632568359375}
  f := 24.019450167736287801;
  testrel(12, NE, y, f, cnt,failed);

  nu := 2;

  y := chi2_inv(nu, 1e-6);
  f := 0.2000001000000666667e-5;
  testrel(13, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.001953125);
  f := 0.39100696716067011153e-2;
  testrel(14, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.03125);
  f := 0.63497396629160602314e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.125);
  f := 0.26706278524904524629;
  testrel(16, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.25);
  f := 0.57536414490356185488;
  testrel(17, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.5);
  f := 1.3862943611198906188;
  testrel(18, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.75);
  f := 2.7725887222397812377;
  testrel(19, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9);
  f := 4.6051701859880913680;
  testrel(20, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.95);
  f := 5.9914645471079819869;
  testrel(21, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.99);
  f := 9.2103403719761827361;
  testrel(22, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9990234375);
  f := 13.862943611198906188;
  testrel(23, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 1-ldexpd(1,-20));  {0.99999904632568359375}
  f := 27.725887222397812377;
  testrel(24, NE, y, f, cnt,failed);

  nu := 5;

  y := chi2_inv(nu, 1e-6);
  f := 0.128961602064970965978034e-1; {Wolfram Alpha: Quantile[ChiSquareDistribution[5], 1e-6]}
  testrel(25, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.001953125);
  f := 0.27738682467162264626;
  testrel(26, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.03125);
  f := 0.92009253480946400288;
  testrel(27, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.125);
  f := 1.8081716938818009134;
  testrel(28, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.25);
  f := 2.6746028094321631155;
  testrel(29, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.5);
  f := 4.3514601910955273172;
  testrel(30, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.75);
  f := 6.6256797638292507690;
  testrel(31, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9);
  f := 9.2363568997811184514;
  testrel(32, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.95);
  f := 11.070497693516354178;
  testrel(33, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.99);
  f := 15.086272469388990113;
  testrel(34, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9990234375);
  f := 20.569688573352666502;
  testrel(35, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 1-ldexpd(1,-20)); {0.99999904632568359375}
  f := 35.991186423468273709;
  testrel(36, NE, y, f, cnt,failed);

  nu := 20;

  y := chi2_inv(nu, 1e-6);
  f := 2.5536375757288159624;
  testrel(37, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.001953125);
  f := 6.4927915055156692220;
  testrel(38, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.03125);
  f := 9.9679153842225177677;
  testrel(39, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.125);
  f := 13.055271826703046031;
  testrel(40, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.25);
  f := 15.451773539047727154;
  testrel(41, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.5);
  f := 19.337429229428262304;
  testrel(42, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.75);
  f := 23.827692043030858503;
  testrel(43, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9);
  f := 28.411980584305633252;
  testrel(44, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.95);
  f := 31.410432844230926553;
  testrel(45, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.99);
  f := 37.566234786625051400;
  testrel(46, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 0.9990234375);
  f := 45.390383389813432931;
  testrel(47, NE, y, f, cnt,failed);

  y := chi2_inv(nu, 1-ldexpd(1,-20)); {0.99999904632568359375}
  f := 65.549649187779351815;
  testrel(48, NE, y, f, cnt,failed);

  y := chi2_inv(100, 0.9921875);
  f := 137.39138405740706884;
  testrel(49, NE, y, f, cnt,failed);

  {extreme cases}
  y := chi2_inv(10000, 0.9921875);
  f := 10345.121945562315395;
{$ifdef CPUARM}
  testrele(50, 2e-16, y, f, cnt,failed);
{$else}
  testrele(50, 1e-17, y, f, cnt,failed);
{$endif}

  y := chi2_inv(120000, 0.9921875);
  f := 121187.58629154332684;
  testrele(51, 1e-15, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_chi2_pdf;
var
  y,f,x: double;
  cnt, failed: integer;
  nu: longint;
const
  NE = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chi2_pdf');

  {statevalf[pdf,chisquare[nu]](x);}

  nu := 1;

  y := chi2_pdf(nu, 1e-10);
  f := 39894.228038148556392;
  testrel( 1, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.001953125);
  f := 9.0182221775452522375;
  testrel( 2, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.125);
  f := 1.0600141293761142435;
  testrel( 3, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.25);
  f := 0.70413065352859895555;
  testrel( 4, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.5);
  f := 0.43939128946772239705;
  testrel( 5, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.75);
  f := 0.31660589975553934697;
  testrel( 6, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1.0);
  f := 0.24197072451914334980;
  testrel( 7, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 2.0);
  f := 0.10377687435514867584;
  testrel( 8, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 4.0);
  f := 0.26995483256594025975e-1;
  testrel( 9, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 10.0);
  f := 0.85003666025203418132e-3;
  testrel(10, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 100.0);
  f := 0.76945986267064193463e-23;
  testrel(11, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1400.0);
  f := 0.10512565523214444743e-305;
  testrel(12, NE, y, f, cnt,failed);

  nu := 2;

  y := chi2_pdf(nu, 0);
  f := 0.5;
  testrel(13, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.001953125);
  f := 0.49951195709098783112;
  testrel(14, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.125);
  f := 0.46970653140673789306;
  testrel(15, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.25);
  f := 0.44124845129229770143;
  testrel(16, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.5);
  f := 0.38940039153570243413;
  testrel(17, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.75);
  f := 0.34364463939548609928;
  testrel(18, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1.0);
  f := 0.30326532985631671180;
  testrel(19, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 2.0);
  f := 0.18393972058572116080;
  testrel(20, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 4.0);
  f := 0.67667641618306345945e-1;
  testrel(21, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 10.0);
  f := 0.33689734995427335483e-2;
  testrel(22, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 100.0);
  f := 0.96437492398195889150e-22;
  testrel(23, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1400.0);
  f := 0.49298382718798854284e-304;
  testrel(24, NE, y, f, cnt,failed);

  nu := 5;

  y := chi2_pdf(nu, 0);
  f := 0.0;
  testrel(25, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.001953125);
  f := 0.11467262493826868994e-4;
  testrel(26, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.125);
  f := 0.55209069238339283517e-2;
  testrel(27, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.25);
  f := 0.14669388615179144908e-1;
  testrel(28, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.5);
  f := 0.36615940788976866421e-1;
  testrel(29, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.75);
  f := 0.59363606204163627561e-1;
  testrel(30, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1.0);
  f := 0.80656908173047783271e-1;
  testrel(31, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 2.0);
  f := 0.13836916580686490112;
  testrel(32, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 4.0);
  f := 0.14397591070183480521;
  testrel(33, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 10.0);
  f := 0.28334555341734472710e-1;
  testrel(34, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 100.0);
  f := 0.25648662089021397821e-19;
  testrel(35, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1400);
  f := 0.68682094751667705661e-300;
  testrel(36, NE, y, f, cnt,failed);

  nu := 20;

  y := chi2_pdf(nu, 0);
  f := 0.0;
  testrel(37, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.001953125);
  f := 0.11119448455436073856e-32;
  testrel(38, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.125);
  f := 0.18835784907520618654e-16;
  testrel(39, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.25);
  f := 0.90596261839205305002e-14;
  testrel(40, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.5);
  f := 0.40934871274926930941e-11;
  testrel(41, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 0.75);
  f := 0.13887624019009048206e-9;
  testrel(42, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1.0);
  f := 0.16322616219566208601e-8;
  testrel(43, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 2.0);
  f := 0.50688855981514870149e-6;
  testrel(44, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 4.0);
  f := 0.95474626621948988987e-4;
  testrel(45, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 10.0);
  f := 0.18132788707821873516e-1;
  testrel(46, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 100.0);
  f := 0.51905444593316343142e-12;
  testrel(47, NE, y, f, cnt,failed);

  y := chi2_pdf(nu, 1400);
  f := 0.54821636959049836799e-284;
  testrel(48, NE, y, f, cnt,failed);

  y := chi2_pdf(100, 99);
  f := 0.28375467778026354860e-1;
  testrel(49, NE, y, f, cnt,failed);

  y := chi2_pdf(1000, 1001);
  f := 0.89079979170919584589e-2;
  testrel(50, 2*NE, y, f, cnt,failed);      {!!64-bit}

  y := chi2_pdf(100000, 100000);
  f := 0.89206057130752775749e-3;
  testrel(51, NE, y, f, cnt,failed);

  y := chi2_pdf(1200000, 1200000);
  f := 0.25751609891599905370e-3;
  testrel(52, NE, y, f, cnt,failed);

  x := 2.01e-17;
  y := chi2_pdf(1,x);
  f := 88984023.13185840146;
  testrel(53, NE, y, f, cnt,failed);

  y := chi2_pdf(2,x);
  f := 0.4999999999999999950;
  testrel(54, NE, y, f, cnt,failed);

  y := chi2_pdf(3,x);
  f := 0.1788578864950353869e-8;
  testrel(55, NE, y, f, cnt,failed);

  y := chi2_pdf(36,x);
  f := 0.1530116643226734334e-303;
  testrel(56, NE, y, f, cnt,failed);

  y := chi2_pdf(38,x);
  f := 0.8543151258015933366e-322;
  testrel(57, NE, y, f, cnt,failed);

  y := chi2_pdf(39,x);
  f := 0;
  testrel(58, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_chi_dist;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 2;
  NE1 = 4;
  NE2 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','chi_pdf/cdf/inv');

  y := chi_pdf(1,0);
  f := 0.7978845608028653559;
  testrel( 1, NE, y, f, cnt,failed);

  y := chi_pdf(1,1e-10);
  f := 0.7978845608028653559;
  testrel( 2, NE, y, f, cnt,failed);

  y := chi_pdf(2,1e-10);
  f := 1e-10;
  testrel( 3, NE, y, f, cnt,failed);

  y := chi_pdf(2,5);
  f := 0.1863326586039335496e-4;
  testrel( 4, NE, y, f, cnt,failed);

  y := chi_cdf(2,1);
  f := 0.3934693402873665764;
  testrel( 5, NE, y, f, cnt,failed);

  y := chi_inv(2,0.25);
  f := 0.7585276164409321326;
  testrel( 6, NE, y, f, cnt,failed);

  y := chi_pdf(3,4);
  f := 0.4282567224476331257e-2;
  testrel( 7, NE, y, f, cnt,failed);

  y := chi_cdf(3,2);
  f := 0.7385358700508893778;
  testrel( 8, NE, y, f, cnt,failed);

  y := chi_inv(3,0.75);
  f := 2.026905260645478921;
  testrel( 9, NE, y, f, cnt,failed);

  y := chi_pdf(3,1e-20);
  f := 7.978845608028653559e-41;
  testrel(10, NE, y, f, cnt,failed);

  y := chi_inv(200,1/8);
  f := 13.31333421845250538;
  testrel(11, NE, y, f, cnt,failed);

  y := chi_inv(10,1/2);
  f := 3.056438739054320951;
  testrel(12, NE, y, f, cnt,failed);

  y := chi_inv(10,1e-7);
  f := 0.4594606793283617908;
  testrel(13, NE, y, f, cnt,failed);

  y := chi_pdf(200,40);
  f := 0.4002944107993220147e-214;
  testrel(14, NE1, y, f, cnt,failed);

  y := chi_cdf(200,15);
  f := 0.8914767979691302269;
  testrel(15, NE, y, f, cnt,failed);

  y := chi_cdf(200,10);
  f := 3.200065324585125294e-10;
  testrel(16, NE2, y, f, cnt,failed);

  y := chi_pdf(200,10);
  f := 3.260638704295460052e-9;
  testrel(17, NE1, y, f, cnt,failed);

  y := chi_inv(200,1/100);
  f := 12.50727652639021383;
  testrel(18, NE, y, f, cnt,failed);

  y := chi_pdf(5,3);
  f := 0.2393198142446523875;
  testrel(19, NE1, y, f, cnt,failed);   {ARM}

  y := chi_cdf(5,8);
  f := 0.9999999999981934109;
  testrel(20, NE, y, f, cnt,failed);

  y := chi_inv(5,9/10);
  f := 3.039137525644589591;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;




{---------------------------------------------------------------------------}
procedure test_exp_cdf;
var
  y,f,a,alpha: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','exp_cdf');

  {statevalf[cdf,exponential[alpha, a]](x);}
  a := 5;
  alpha := 2.0;
  y := exp_cdf(a, alpha, 5);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 5.015625);
  f := 0.3076676552365591815e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 5.25);
  f := 0.3934693402873665764;
  testrel(3, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 5.5);
  f := 0.6321205588285576784;
  testrel(4, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 6);
  f := 0.8646647167633873081;
  testrel(5, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 8);
  f := 0.9975212478233336416;
  testrel(6, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 10);
  f := 0.9999546000702375152;
  testrel(7, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 20);
  f := 0.9999999999999064238;
  testrel(8, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 28);
  f := 1.0;
  testrel(9, NE, y, f, cnt,failed);

  a := -1.5;
  alpha := 0.5;
  y := exp_cdf(a, alpha, -2);
  f := 0.0;
  testrel(10, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, -1.5);
  f := 0.0;
  testrel(11, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, -1);
  f := 0.2211992169285951318;
  testrel(12, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, -0.5);
  f := 0.3934693402873665764;
  testrel(13, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 0);
  f := 0.5276334472589852929;
  testrel(14, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 1);
  f := 0.7134952031398098997;
  testrel(15, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 2);
  f := 0.8262260565495548733;
  testrel(16, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 20);
  {$ifdef BIT16}
  f := 0.999978554591683411;
  {$else}
  f := 0.9999785545916834108;
  {$endif}
  testrel(17, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 50);
  f := 0.9999999999934397998;
  testrel(18, NE, y, f, cnt,failed);

  y := exp_cdf(a, alpha, 90);
  f := 1.0;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_exp_pdf;
var
  y,f,a,alpha: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','exp_pdf');

  {statevalf[pdf,exponential[alpha,a]](x);}
  a := 5;
  alpha := 2.0;
  y := exp_pdf(a, alpha, 2.5);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 5);
  f := 2.0;
  testrel(2, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 7.5);
  f := 0.1347589399817093419e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 10);
  f := 0.9079985952496970307e-4;
  testrel(4, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 100);
  f := 0.6096469901943713570e-82;
  testrel(5, NE+2, y, f, cnt,failed);  {+2 for ARM}


  a := -1.5;
  alpha := 0.5;
  y := exp_pdf(a, alpha, -2);
  f := 0;
  testrel(6, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, -1.5);
  f := 0.5;
  testrel(7, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, -1);
  f := 0.3894003915357024341;
  testrel(8, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 0);
  f := 0.2361832763705073536;
  testrel(9, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 2.5);
  f := 0.6766764161830634595e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 5);
  f := 0.1938710391586100494e-1;
  testrel(11, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 10);
  f := 0.15913903982548335340e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 25);
  f := 0.88017315607808464929e-6;
  testrel(13, NE, y, f, cnt,failed);

  y := exp_pdf(a, alpha, 1000);
  f := 0.16827057984961732511e-217;
  testrel(14, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;

{---------------------------------------------------------------------------}
procedure test_exp_inv;
var
  y,f,a,alpha: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','exp_inv');
  a := 5;
  alpha := 2.0;
  y := exp_inv(a, alpha, ldexpd(1,-20));;
  f := 5.000000476837385577;
  testrel(1, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.125);
  f := 5.066765696312261312;
  testrel(2, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.5);
  f := 5.346573590279972655;
  testrel(3, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.75);
  f := 5.693147180559945309;
  testrel(4, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.875);
  f := 6.039720770839917964;
  testrel(5, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.9);
  f := 6.151292546497022842;
  testrel(6, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.99);
  f := 7.302585092994045684;
  testrel(7, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.9990234375);
  f := 8.465735902799726547;
  testrel(8, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 1-ldexpd(1,-20));
  f := 11.93147180559945309;
  testrel(9, NE, y, f, cnt,failed);

  a := -1.5;
  alpha := 0.5;
  y := exp_inv(a, alpha, ldexpd(1,-20));;
  f := -1.499998092650457692;
  testrel(10, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.125);
  f := -1.232937214750954754;
  testrel(11, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.5);
  f := -0.113705638880109381167;
  testrel(12, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.75);
  f := 1.272588722239781238;
  testrel(13, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.875);
  f := 2.658883083359671857;
  testrel(14, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.9);
  f := 3.105170185988091368;
  testrel(15, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.99);  {!!! not exact}
  f := 7.710340371976182736;
  testrel(16, 2, y, f, cnt,failed);

  y := exp_inv(a, alpha, 0.9990234375);
  f := 12.36294361119890619;
  testrel(17, NE, y, f, cnt,failed);

  y := exp_inv(a, alpha, 1-ldexpd(1,-20));
  f := 26.22588722239781238;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_gamma_pdf;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gamma_pdf');

  {statevalf[pdf,gamma[a,b](x);}

  a := 1.0; b := 2.0;
  y := gamma_pdf(a,b, 1e-10);
  f := 0.499999999975;
  testrel( 1, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.001953125);
  f := 0.49951195709098783112;
  testrel( 2, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.125);
  f := 0.46970653140673789306;
  testrel( 3, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.5);
  f := 0.38940039153570243413;
  testrel( 4, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1.0);
  f := 0.30326532985631671180;
  testrel( 5, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 2.0);
  f := 0.18393972058572116080;
  testrel( 6, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 4.0);
  f := 0.67667641618306345945e-1;
  testrel( 7, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 10.0);
  f := 0.33689734995427335483e-2;
  testrel( 8, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 100.0);
  f := 0.96437492398195889150e-22;
  testrel(09, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1400.0);
  f := 0.49298382718798854284e-304;
  testrel(10, NE, y, f, cnt,failed);

  a := 2.0; b := 2.0;
  y := gamma_pdf(a,b, 1e-10);
  f := 0.2499999999875e-10;
  testrel(11, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.001953125);
  f := 0.48780464559666780383e-3;
  testrel(12, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.125);
  f := 0.29356658212921118318e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.5);
  f := 0.97350097883925608533e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1.0);
  f := 0.15163266492815835590;
  testrel(15, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 2.0);
  f := 0.18393972058572116080;
  testrel(16, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 4.0);
  f := 0.13533528323661269189;
  testrel(17, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 10.0);
  f := 0.16844867497713667742e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 100.0);
  f := 0.48218746199097944575e-20;
  testrel(19, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1400);
  f := 0.34508867903159197998e-301;
  testrel(20, NE, y, f, cnt,failed);


  a := 3.0; b := 2.0;
  y := gamma_pdf(a,b, 1e-10);
  f := 0.62499999996875e-21;
  testrel(21, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.001953125);
  f := 0.23818586210774795109e-6;
  testrel(22, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.125);
  f := 0.91739556915378494738e-3;
  testrel(23, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.5);
  f := 0.12168762235490701066e-1;
  testrel(24, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1.0);
  f := 0.37908166232039588975e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 2.0);
  f := 0.91969860292860580400e-1;
  testrel(26, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 4.0);
  f := 0.13533528323661269189;
  testrel(27, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 10.0);
  f := 0.42112168744284169354e-1;
  testrel(28, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 100.0);
  f := 0.12054686549774486144e-18;
  testrel(29, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1400);
  f := 0.12078103766105719299e-298;
  testrel(30, NE, y, f, cnt,failed);

  a := 9.0; b := 0.5;
  y := gamma_pdf(a,b, 1e-10);
  f := 0.12698412695873015873e-81;
  testrel(31, 2*NE, y, f, cnt,failed);       {!!!!!!!!!}

  y := gamma_pdf(a,b, 0.001953125);
  f := 0.26785100912581535303e-23;
  testrel(32, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.125);
  f := 0.58946214635894780554e-9;
  testrel(33, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.5);
  f := 0.18247988153345353254e-4;
  testrel(34, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1.0);
  f := 0.17185432791950818018e-2;
  testrel(35, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 2.0);
  f := 0.59540362609726351176e-1;
  testrel(36, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 4.0);
  f := 0.27917306390119385230;
  testrel(37, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 10.0);
  f := 0.26173379332553115276e-2;
  testrel(38, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 100.0);
  f := 0.17573289228403016261e-72;
  testrel(39, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 350);
  f := 0.28193984721797058926e-285;
  testrel(40, NE, y, f, cnt,failed);


  a := 5.0; b := 1.0;
  y := gamma_pdf(a,b, 1e-10);
  f := 0.41666666662500000001e-41;
  testrel(41, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.001953125);
  f := 0.60514671901878529313e-12;
  testrel(42, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.125);
  f := 0.89772227232319682097e-5;
  testrel(43, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 0.5);
  f := 0.15795069263349828740e-2;
  testrel(44, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 1.0);
  f := 0.15328310048810096733e-1;
  testrel(45, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 2.0);
  f := 0.90223522157741794592e-1;
  testrel(46, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 4.0);
  f := 0.19536681481316458981;
  testrel(47, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 10.0);
  f := 0.18916637401035354807e-1;
  testrel(48, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 100.0);
  f := 0.15500316566753483179e-36;
  testrel(49, NE, y, f, cnt,failed);

  y := gamma_pdf(a,b, 700);
  f := 0.98637847423196707612e-294;
  testrel(50, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_gamma_cdf;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gamma_cdf');

  {statevalf[cdf,gamma[a,b](x);}

  a := 1.0; b := 2.0;
  y := gamma_cdf(a,b, 1e-10);
  f := 0.49999999998750000000e-10;
  testrel( 1, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.001953125);
  f := 0.97608581802433776529e-3;
  testrel( 2, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.125);
  f := 0.60586937186524213880e-1;
  testrel( 3, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.5);
  f := 0.22119921692859513175;
  testrel( 4, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 1.0);
  f := 0.39346934028736657640;
  testrel( 5, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 2.0);
  f := 0.63212055882855767840;
  testrel( 6, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 4.0);
  f := 0.86466471676338730811;
  testrel( 7, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 10.0);
  f := 0.99326205300091453290;
  testrel( 8, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 100.0);
  f := 1.0;
  testrel(09, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.25);
  f := 0.11750309741540459714;
  testrel(10, NE, y, f, cnt,failed);


  a := 2.0; b := 2.0;
  y := gamma_cdf(a,b, 1e-10);
  f := 0.12499999999583333333e-20;
  testrel(11, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.001953125);
  f := 0.47652683100215763712e-6;
  testrel(12, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.125);
  f := 0.18736207606819772478e-2;
  testrel(13, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.5);
  f := 0.26499021160743914694e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 1.0);
  f := 0.90204010431049864594e-1;
  testrel(15, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 2.0);
  f := 0.26424111765711535681;
  testrel(16, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 4.0);
  f := 0.59399415029016192432;
  testrel(17, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 10.0);
  f := 0.95957231800548719742;
  testrel(18, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 100.0);
  f := 0.99999999999999999999;
  testrel(19, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.125);
  f := 0.18736207606819772478e-2;
  testrel(20, NE, y, f, cnt,failed);


  a := 3.0; b := 2.0;
  y := gamma_cdf(a,b, 1e-10);
  f := 0.20833333332552083333e-31;
  testrel(21, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.001953125);
  f := 0.15510678666173495128e-9;
  testrel(22, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.125);
  f := 0.38829622374407353042e-4;
  testrel(23, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.5);
  f := 0.21614966897625125609e-2;
  testrel(24, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 1.0);
  f := 0.14387677966970686644e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 2.0);
  f := 0.80301397071394196011e-1;
  testrel(26, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 4.0);
  f := 0.32332358381693654053;
  testrel(27, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 10.0);
  f := 0.87534798051691885871;
  testrel(28, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 100.0);
  f := 0.99999999999999999975;
  testrel(29, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.25);
  f := 0.29647754088802019211e-3;
  testrel(30, NE, y, f, cnt,failed);


  a := 9.0; b := 0.5;
  y := gamma_cdf(a,b, 1e-10);
  f := 0.14109347440141093475e-92;
  testrel(31, 2*NE, y, f, cnt,failed);             {!!!!!!!!!}

  y := gamma_cdf(a,b, 0.001953125);
  f := 0.58150103210426417249e-27;
  testrel(32, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.125);
  f := 0.83963991089851222999e-11;
  testrel(33, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.5);
  f := 0.11252025979690180803e-5;
  testrel(34, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 1.0);
  f := 0.23744732826116178623e-3;
  testrel(35, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 2.0);
  f := 0.21363434487984163417e-1;
  testrel(36, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 4.0);
  f := 0.40745265856240858839;
  testrel(37, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 10.0);
  f := 0.99791274095086498120;
  testrel(38, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 100.0);
  f := 1.0;
  testrel(39, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.25);
  f := 0.34354902468481320560e-8;
  testrel(40, NE, y, f, cnt,failed);


  a := 5.0; b := 1.0;
  y := gamma_cdf(a,b, 1e-10);
  f := 0.83333333326388888889e-52;
  testrel(41, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.001953125);
  f := 0.23646240697640727994e-15;
  testrel(42, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.125);
  f := 0.22919102136524120147e-6;
  testrel(43, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.5);
  f := 0.17211562995584077811e-3;
  testrel(44, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 1.0);
  f := 0.36598468273437123455e-2;
  testrel(45, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 2.0);
  f := 0.52653017343711156742e-1;
  testrel(46, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 4.0);
  f := 0.37116306482012647658;
  testrel(47, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 10.0);
  f := 0.97074731192303892733;
  testrel(48, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 100.0);
  f := 1.0;
  testrel(49, NE, y, f, cnt,failed);

  y := gamma_cdf(a,b, 0.25);
  f := 0.66117105610342470461853e-5;
  testrel(50, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_gamma_inv;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','gamma_inv');

  {statevalf[icdf,gamma[a,b](p);}

  a := 1.0; b := 2.0;

  y := gamma_inv(a, b, 1e-6);
  f := 0.20000010000006666672e-5;
  testrel( 1, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.125);
  f := 0.26706278524904524629;
  testrel( 2, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.25);
  f := 0.57536414490356185488;
  testrel( 3, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.5);
  f := 1.3862943611198906188;
  testrel( 4, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.75);
  f := 2.7725887222397812377;
  testrel( 5, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9);
  f := 4.6051701859880913680;
  testrel( 6, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.99);
  f := 9.2103403719761827361;
  testrel( 7, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9990234375);
  f := 13.862943611198906188;
  testrel( 8, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 1-ldexpd(1,-20));
  f := 27.725887222397812377;
  testrel( 9, NE, y, f, cnt,failed);


  a := 2.0; b := 2.0;

  y := gamma_inv(a, b, 1e-6);
  f := 0.28297613229586858205e-2;
  testrel(11, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.125);
  f := 1.2187621354461568880;
  testrel(12, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.25);
  f := 1.9225575262295541917;
  testrel(13, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.5);
  f := 3.3566939800333213068;
  testrel(14, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.75);
  f := 5.3852690577793915415;
  testrel(15, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9);
  f := 7.7794403397348581158;
  testrel(16, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.99);
  f := 13.276704135987624539;
  testrel(17, 2*NE, y, f, cnt,failed);       {!!!!!!!!!!!!}

  y := gamma_inv(a, b, 0.9990234375);
  f := 18.519389817856053111;
  testrel(18, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 1-ldexpd(1,-20));
  f := 33.477383689881828147;
  testrel(19, NE, y, f, cnt,failed);


  a := 3.0; b := 2.0;

  y := gamma_inv(a, b, 1e-6);
  f := 0.36508565926558585216e-1;
  testrel(21, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.125);
  f := 2.4411036298316021684;
  testrel(22, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.25);
  f := 3.4545988357210387862;
  testrel(23, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.5);
  f := 5.3481206274471206358;
  testrel(24, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.75);
  f := 7.8408041205851200986;
  testrel(25, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9);
  f := 10.644640675668419809;
  testrel(26, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.99);
  f := 16.811893829770931055;
  testrel(27, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9990234375);
  f := 22.514365804408070772;
  testrel(28, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 1-ldexpd(1,-20));
  f := 38.363624414049344648;
  testrel(29, NE, y, f, cnt,failed);

  {Maple V R4 failed on most of the following test cases}
  {Wolfram Alpha: InverseGammaRegularized[a,1-p]*b to xx digits}
  a := 9.0; b := 0.5;

  y := gamma_inv(a, b, 1e-6);
  f := 0.4927844783647370284;
  testrel(31, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.125);
  f := 2.8587200022579131718;
  testrel(32, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.25);
  f := 3.4188225875995730573;
  testrel(33, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.5);
  f := 4.3344755921851860761;
  testrel(34, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.75);
  f := 5.4012224489320413344;
  testrel(35, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9);
  f := 6.4973557706593028874;
  testrel(36, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.99);
  f := 8.7013264336762682160;
  testrel(37, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9990234375);
  f := 10.596459041791414033;
  testrel(38, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 1-ldexpd(1,-20));
  f := 15.510096022825544320;
  testrel(39, NE, y, f, cnt,failed);


  a := 5.0; b := 1.0;

  y := gamma_inv(a, b, 1e-6);
  f := 0.16906300162147725805;
  testrel(41, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.125);
  f := 2.6170585017706683721;
  testrel(42, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.25);
  f := 3.3686003859773210579;
  testrel(43, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.5);
  f := 4.6709088827959837203;
  testrel(44, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.75);
  f := 6.2744306984446885035;
  testrel(45, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9);
  f := 7.9935895860526304375;
  testrel(46, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.99);
  f := 11.604625579477179839;
  testrel(47, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 0.9990234375);
  f := 14.825755928410414354;
  testrel(18, NE, y, f, cnt,failed);

  y := gamma_inv(a, b, 1-ldexpd(1,-20));
  f := 23.488170156554407308;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_invgamma_pdf;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','invgamma_pdf');

  {Maple:  ig_pdf := (a,b,x) -> b^a*x^(-a-1)*exp(-b/x)/GAMMA(a);}
  a := 1.0;  b := 2.0;
  y := invgamma_pdf(a,b,0.03125);
  f := 0.3284604703843610323e-24;
  testrel(1, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,0.125);
  f := 0.1440450236406516666e-4;
  testrel(2, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,1);
  f := 0.2706705664732253838;
  testrel(3, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,10.0);
  f := 0.1637461506155963717e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,100.0);
  f := 0.1960397346613510604e-3;
  testrel(5, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,1000);
  f := 0.1996003997334666134e-5;
  testrel(6, NE, y, f, cnt,failed);

  a := 2.0; b := 1.0;
  y := invgamma_pdf(a,b,0.03125);
  f := 0.4149793767127179501e-9;
  testrel(7, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,0.125);
  f := 0.1717568654860860615;
  testrel(8, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,0.5);
  f := 1.082682265892901535;
  testrel(9, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,1);
  f := 0.3678794411714423216;
  testrel(10, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,10.0);
  f := 0.9048374180359595732e-3;
  testrel(11, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,100.0);
  f := 0.9900498337491680536e-6;
  testrel(12, NE, y, f, cnt,failed);

  a := 0.5; b := 4.0;
  y := invgamma_pdf(a,b,0.03125);
  f := 0.5253954932694545881e-53;
  testrel(13, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,0.125);
  f := 0.3233453493463611064e-12;
  testrel(14, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,0.5);
  f := 0.1070641806119082814e-2;
  testrel(15, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,5.0);
  f := 0.4534866089799165245e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,50.0);
  f := 0.2946161122426586462e-2;
  testrel(17, NE, y, f, cnt,failed);

  y := invgamma_pdf(a,b,5000.0);
  f := 0.3188986033636839945e-5;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_invgamma_cdf;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 4;
  NE1 = 20;
  NE2 = 80;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','invgamma_cdf');

  {Maple: ig_cdf := (a,b,x) -> GAMMA(a,b/x)/GAMMA(a); }
  a := 1.0; b := 2.0;
  y := invgamma_cdf(a,b,0.03125);
  f := 0.1603810890548637853e-27;
  testrel(1, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,0.125);
  f := 0.1125351747192591145e-6;
  testrel(2, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,1);
  f := 0.1353352832366126919;
  testrel(3, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,10.0);
  f := 0.8187307530779818587;
  testrel(4, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,100.0);
  f := 0.9801986733067553022;
  testrel(5, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,1000);
  f := 0.9980019986673330668;
  testrel(6, NE, y, f, cnt,failed);

  a := 2.0; b := 1.0;
  y := invgamma_cdf(a,b,0.03125);
  f := 0.4179174631201077989e-12;
  testrel(7, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,0.125);
  f := 0.3019163651122606549e-2;
  testrel(8, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,0.5);
  f := 0.4060058497098380757;
  testrel(9, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,1);
  f := 0.7357588823428846432;
  testrel(10, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,10.0);
  f := 0.9953211598395555305;
  testrel(11, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,100.0);
  f := 0.9999503320866597341;
  testrel(12, NE, y, f, cnt,failed);

  a := 0.5; b := 4.0;
  y := invgamma_cdf(a,b,0.03125);
  f := 0.1277750880107617456e-56;
  testrel(13, NE2, y, f, cnt,failed);

  y := invgamma_cdf(a,b,0.125);
  f := 0.1244192114854356825e-14;
  testrel(14, NE1, y, f, cnt,failed);

  y := invgamma_cdf(a,b,1);
  f := 0.4677734981047265838e-2;
  testrel(15, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,5.0);
  f := 0.2059032107320683089;
  testrel(16, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,50.0);
  f := 0.6891565167793516665;
  testrel(17, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,250.0);
  f := 0.8580276569875211913;
  testrel(18, NE, y, f, cnt,failed);

  y := invgamma_cdf(a,b,5000.0);
  f := 0.9680931262943384772;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_invgamma_inv;
var
  a,b,y,f: double;
  cnt, failed: integer;
const
  NE = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','invgamma_inv');

  {Maple: ig_inv := (a,b,y) -> b/stats[statevalf,icdf,gamma[a,1]](1-y);}
  a := 1.0; b := 2.0;
  y := invgamma_inv(a,b,1e-6);
  f := 0.1447648273010839426;
  testrel(1, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.125);
  f := 0.961796693925975605;
  testrel(2, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.5);
  f := 2.885390081777926815;
  testrel(3, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.75);
  f := 6.952118993564413821;
  testrel(4, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.9);
  f := 18.98244316205980605;
  testrel(5, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.9990234375);
  f := 2046.999837160061286;
  testrel(6, NE, y, f, cnt,failed);

  a := 2.0; b := 1.0;
  y := invgamma_inv(a,b,1e-6);
  f := 0.5992178723991007878e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.125);
  f := 0.2772368934910238455;
  testrel(8, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.5);
  f := 0.5958243473776976101;
  testrel(9, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.75);
  f := 1.040280965700060485;
  testrel(10, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.9);
  f := 1.880365122205807787;
  testrel(11, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.9990234375);
  f := 22.29221271800046843;
  testrel(12, NE, y, f, cnt,failed);

  a := 0.5; b := 4.0;
  y := invgamma_inv(a,b,1e-10);
  f := 0.1912893690316530870;
  testrel(13, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,1e-6);
  f := 0.3343345681720714721;
  testrel(14, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.125);
  f := 3.399155364424810784;
  testrel(15, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.5);
  f := 17.58487470654185923;
  testrel(16, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.75);
  f := 78.79363457459499543;
  testrel(17, NE, y, f, cnt,failed);

  y := invgamma_inv(a,b,0.9990234375);
  f := 5340351.048773805654;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;

end.
